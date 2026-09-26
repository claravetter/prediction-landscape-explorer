#!/usr/bin/env Rscript
file<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]);root<-normalizePath(file.path(dirname(file),".."))
suppressPackageStartupMessages({library(dplyr);library(tibble);library(purrr);library(ggplot2);library(pROC);library(patchwork)})
for (f in c("engine","demo_data","results")) source(file.path(root,"R",paste0(f,".R")))
n<-0L;check<-function(x,label){stopifnot(isTRUE(x));n<<-n+1L;cat("PASS:",label,"\n")};fails<-function(x)inherits(tryCatch(force(x),error=function(e)e),"error")
raw<-make_demo_data(); map<-default_col_map_for(raw); d<-prep_data_for_app(raw,map)
check(nrow(d)==1200 && identical(make_demo_data(),raw),"deterministic in-memory synthetic example")
check(identical(binary_mapping(factor(c("no","yes",NA)),"no","yes"),c(0L,1L,NA_integer_)),"explicit factor outcome mapping")
check(fails(binary_mapping(c(0,.5,1))),"multiclass/fractional labels rejected")
raw$EXP_LABEL<-10+seq_len(nrow(raw));bad<-map;bad$moderators<-c("age","EXP_LABEL")
check(fails(prep_data_for_app(raw,bad)),"reserved moderator cannot overwrite outcome")
raw$Mean_Score<-raw$calibrated_prob;bad<-map;bad$prob<-"Mean_Score"
z<-prep_data_for_app(raw,bad)
check(identical(z$Mean_Score,raw$risk_score)&&identical(z$prob_work,raw$Mean_Score),"probability source name cannot overwrite score")
raw$EXP_LABEL<-raw$calibrated_prob;bad$prob<-"EXP_LABEL";z<-prep_data_for_app(raw,bad)
check(identical(z$EXP_LABEL,raw$event)&&identical(z$prob_work,raw$calibrated_prob),"probability source name cannot overwrite outcome")
raw$age_z<-seq_len(nrow(raw));bad<-map;bad$moderators<-c("age_z","age");z<-prep_data_for_app(raw,bad)
check(identical(z$age_z,as.numeric(raw$age_z)),"supplied z moderator survives generated sibling collision")
missing_prob<-raw;missing_prob$calibrated_prob<-NA_real_
check(all(is.na(prep_data_for_app(missing_prob,map)$prob_work)),"explicit missing probabilities never silently become transformed scores")
unmapped<-raw;unmapped$Probs_platt<-seq(2,3,length.out=nrow(raw))
score_only<-map;score_only$prob<-NULL;score_only$moderators<-c("age","Probs_platt")
z<-prep_data_for_app(unmapped,score_only)
check(identical(z$prob_work,plogis(z$decision_z))&&grepl("logistic_standardized_score",attr(z,"probability_method"),fixed=TRUE),"unselected probability-like moderator uses the documented score transform")
check(fails(make_prob_work(d,prob_col="missing_probability")),"explicit unknown probability cannot silently use a score transform")
raw$calibrated_prob[1]<-1.1
check(fails(prep_data_for_app(raw,map)),"out-of-range probability rejected")
# Non-finite original scores cannot manufacture finite substitute probabilities.
finite_raw<-make_demo_data();score_map<-default_col_map_for(finite_raw);score_map$prob<-NULL
infinite<-finite_raw;infinite$risk_score[1]<-Inf
z<-prep_data_for_app(infinite,score_map);finite_reference<-prep_data_for_app(finite_raw[-1,],score_map)
check(is.na(z$Mean_Score[1])&&is.na(z$prob_work[1])&&isTRUE(all.equal(z$prob_work[-1],finite_reference$prob_work)),"one infinite score is excluded without changing finite-score standardization")
infinite$risk_score<-rep(c(-Inf,Inf),nrow(infinite)/2);z<-prep_data_for_app(infinite,score_map)
check(all(is.na(z$Mean_Score))&&all(is.na(z$prob_work)),"all-infinite scores yield no derived probabilities")
constant<-finite_raw;constant$risk_score<-3;constant$risk_score[1:4]<-c(Inf,-Inf,NaN,NA_real_);z<-prep_data_for_app(constant,score_map)
check(all(is.na(z$prob_work[1:4]))&&all(z$prob_work[-(1:4)]==.5),"constant finite scores use rank fallback without including nonfinite values")
constant$risk_score[]<-NA_real_;z<-prep_data_for_app(constant,score_map)
check(all(is.na(z$prob_work)),"all-missing scores remain missing through rank fallback")
z<-prep_data_for_app(infinite,default_col_map_for(infinite))
check(all(is.na(z$Mean_Score))&&identical(z$prob_work,infinite$calibrated_prob),"supplied valid probabilities remain independent of infinite original scores")
z<-make_prob_work(tibble(Mean_Score=c(-Inf,1,2,Inf,NA)))
check(all(is.na(z$prob_work[c(1,4,5)]))&&identical(z$prob_work[2:3],c(.25,.75)),"direct rank fallback only ranks finite scores")
check(tail(tail(make_windows(103,30,20),1)[[1]],1)==103,"fixed-count final window reaches last row")
grid<-expand.grid(x_center=0:3,y_center=0:3);grid$auc<-.5;grid$auc[grid$x_center==1&grid$y_center==1]<-NA
nd<-data.frame(x_center=c(.5,2.5,4,1),y_center=c(.5,2.5,2,1))
check(identical(supported_grid_points(grid,nd,"auc"),c(FALSE,TRUE,FALSE,FALSE)),"GAM support cannot bridge a hole or extrapolate")
for(dim in c("1d","2d")) for(method in if(dim=="1d")c("fixedrange","fixedn") else c("fixedrange","knn")) {
 s<-result_defaults(dim);s$method<-method;s$moderator<-"age";s$x_var<-"age";s$y_var<-"cognitive_score";s$surface_mode<-"raw"
 request<-list(data=d,settings=s,source="synthetic_demo",mapping=list(outcome="event",positive="1"))
 r<-compute_result(request,dim,FALSE); e<-result_export(r)
 check(nrow(e)>0&&all(e$source=="synthetic_demo")&&all(e$software_version==APP_VERSION),paste(dim,method,"computed export has provenance"))
 meta<-jsonlite::fromJSON(e$applied_settings_json[1]);check(meta$settings$method==method&&meta$mapping$positive=="1",paste(dim,method,"export settings match applied request"))
 check(grepl(method,result_filename(r,"csv"),fixed=TRUE),paste(dim,method,"filename identifies actual method"))
 file<-tempfile(fileext=".png");ggsave(file,plot_result(r),width=11,height=9,dpi=70);check(file.info(file)$size>1000,paste(dim,method,"PNG render"));unlink(file)
 if(dim=="2d") {
  check(result_summary(r)$event_memberships>r$eligible_events,"overlap counts are memberships, not events")
  check(hist_top_from_grid(r$table)$labels$y=="Event memberships","2D marginal label states repeated memberships")
 }
}
uneven<-expand.grid(x_center=c(0,1,4),y_center=c(0,2,3));uneven$auc<-seq(.5,.9,length.out=9);uneven$n_transition<-10
layer<-ggplot_build(plot_grid_raw(uneven,overlay_points=FALSE))$data[[1]]
check(identical(sort(unique(layer$xmin)),c(-.5,.5,2.5))&&identical(sort(unique(layer$xmax)),c(.5,2.5,5.5)),"uneven grid centres retain their display boundaries")
decimal_grid<-expand.grid(x_center=c(.0123456789012345,1.23456789012345,4.01234567890123),y_center=c(.123456789,2,3));decimal_grid$auc<-.6
probes<-expand.grid(x_center=c(.5,2),y_center=c(.5,2.5))
check(all(supported_grid_points(decimal_grid,probes,"auc")),"decimal grid lookup does not blank supported cells")
decimal_grid$auc[decimal_grid$y_center==2]<-NA_real_
check(!any(supported_grid_points(decimal_grid,probes,"auc")),"entire unsupported grid row remains a hole")
s<-result_defaults("2d");s$x_var<-"age";s$y_var<-"cognitive_score";s$method<-"knn"
r<-compute_result(list(data=d,settings=s,source="synthetic_demo"),"2d",FALSE)
smooth<-gam_surface_from_grid(r$table,metric="bacc")
check(all(is.finite(smooth$surface$pred)),"fully supported KNN grid retains its complete GAM surface")
hole_grid<-expand.grid(x_center=0:4,y_center=0:4);hole_grid$auc<-.6+.01*sin(hole_grid$x_center+hole_grid$y_center);hole_grid$n_obs<-100;hole_grid$n_transition<-20
hole_grid$auc[hole_grid$y_center==2]<-NA_real_
hole_surface<-gam_surface_from_grid(hole_grid,k=10,family="gaussian")$surface
check(all(is.na(hole_surface$pred[hole_surface$y_center>1 & hole_surface$y_center<3])),"GAM integration retains an entirely unsupported row")
cat("TOTAL",n,"engine assertions passed\n")
