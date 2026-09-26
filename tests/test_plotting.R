#!/usr/bin/env Rscript
file<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]);root<-normalizePath(file.path(dirname(file),".."))
suppressPackageStartupMessages({library(dplyr);library(tibble);library(purrr);library(ggplot2);library(pROC);library(patchwork)})
for(f in c("engine","demo_data","results"))source(file.path(root,"R",paste0(f,".R")))
n<-0L;check<-function(x,label){stopifnot(isTRUE(x));n<<-n+1L;cat("PASS:",label,"\n")}
raw<-make_demo_data();raw$calibrated_prob[raw$age<33]<-NA_real_
d<-prep_data_for_app(raw,default_col_map_for(raw))
s<-result_defaults("2d");s$x_var<-"age";s$y_var<-"cognitive_score";s$metric<-"bacc"
r<-compute_result(list(data=d,settings=s,source="generated_test"),"2d",FALSE)
p<-prepare_result_plot(r)
check(sum(is.finite(r$table$bacc))==7&&!p$ok&&is.null(p$plot),"seven supported cells cannot enable a GAM plot")
check(grepl("at least 20",p$message)&&nrow(result_export(r))>0,"sparse GAM has an actionable message and exportable numerical results")
r$request$settings$surface_mode<-"raw";p<-prepare_result_plot(r)
check(p$ok&&!is.null(p$plot),"raw display works for the same sparse grid")
file<-tempfile(fileext=".png");ggsave(file,p$plot,width=11,height=9,dpi=70)
check(file.info(file)$size>1000,"prepared sparse raw plot renders to PNG");unlink(file)
r$table$bacc<-NA_real_;p<-prepare_result_plot(r)
check(!p$ok&&grepl("No supported values",p$message),"zero finite values cannot enable a raw plot")
r$table<-r$table[FALSE,];p<-prepare_result_plot(r)
check(!p$ok&&grepl("No windows",p$message),"empty result cannot enable a plot")

# A spatially varied grid tests the exact display minimum independently of
# the number of windows produced by a particular participant dataset.
g<-expand.grid(x_center=0:4,y_center=0:3)
g$bacc<-.55+.03*sin(g$x_center+g$y_center);g$n_obs<-100;g$n_transition<-30
g$threshold_n_obs<-80;g$threshold_n_events<-20
r$table<-g;r$request$settings$surface_mode<-"gam";r$request$settings$gam_family<-"gaussian";r$request$settings$gam_k<-10
p<-prepare_result_plot(r)
check(p$ok,"exactly twenty usable cells can fit a GAM")
file<-tempfile(fileext=".png");ggsave(file,p$plot,width=11,height=9,dpi=70)
check(file.info(file)$size>1000,"prepared twenty-cell GAM renders to PNG");unlink(file)
r$table<-g[-1,];p<-prepare_result_plot(r)
check(!p$ok&&grepl("at least 20",p$message),"nineteen cells are insufficient for GAM display")
r$table<-g;r$table$threshold_n_obs[1]<-0;p<-prepare_result_plot(r)
check(!p$ok,"zero-weight cells do not count towards the GAM display minimum")
r$table<-g;r$request$settings$gam_k<-200;p<-prepare_result_plot(r)
check(!p$ok&&grepl("Unable to fit",p$message)&&is.null(p$plot),"real fit failure above the count minimum becomes an unavailable plot")
r$request$settings$gam_k<-10
check(prepare_result_plot(r)$ok,"valid settings recover after a fit failure")
cat("TOTAL",n,"plotting assertions passed\n")
