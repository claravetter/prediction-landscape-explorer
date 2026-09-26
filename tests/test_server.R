#!/usr/bin/env Rscript
file<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]);setwd(file.path(dirname(normalizePath(file)),".."))
app<-source("app.R")$value
n<-0L;check<-function(x,label){stopifnot(isTRUE(x));n<<-n+1L;cat("PASS:",label,"\n")};fails<-function(x)inherits(tryCatch(force(x),error=function(e)e),"error")
shiny::testServer(server, {
 session$setInputs(data_source="example")
 check(nrow(raw_rv())==1200,"server generates demo without data file")
 session$setInputs(col_outcome="event",col_score="risk_score",col_prob="calibrated_prob",col_cohort="group",col_moderators=c("age","cognitive_score"),col_outcome_negative="0",col_outcome_positive="1",cohorts=unique(make_demo_data()$group))
 s1<-result_defaults("1d");s1$moderator<-"age"
 s2<-result_defaults("2d");s2$x_var<-"age";s2$y_var<-"cognitive_score";s2$surface_mode<-"raw"
 do.call(session$setInputs,c(setNames(s1,paste0(names(s1),"_1d")),setNames(s2,paste0(names(s2),"_2d"))))
 session$setInputs(update_1d=1,update_2d=1)
 check(nrow(result_1d()$table)>0&&nrow(result_2d()$table)>0,"both dimensions compute from applied settings")
 old<-cached_1d()
 session$setInputs(width_pct_1d=.3)
 check(fails(result_1d())&&identical(cached_1d(),old),"changed width invalidates cached 1D output")
 check(grepl("Click Update",output$result_status_1d),"visible outdated status")
 session$setInputs(update_1d=2)
 check(result_1d()$request$settings$width_pct==.3,"Update binds new settings")
 session$setInputs(metric_1d="bacc")
 check(fails(result_1d()),"changed display metric invalidates cached labels and plot")
 session$setInputs(cohorts="Cohort_A")
 check(fails(result_2d()),"cohort filtering invalidates 2D results")
 session$setInputs(col_outcome_positive="invalid")
 check(is.null(dat_rv())&&fails(result_1d())&&fails(result_2d()),"invalid mapping clears previous prepared data and outputs")
 session$setInputs(col_outcome_positive="1",cohorts=unique(make_demo_data()$group),metric_1d="auc",method_2d="knn")
 session$setInputs(update_2d=2)
 check(result_2d()$request$settings$method=="knn","2D alternative recalculates and labels KNN")
 # CSV and RDS use the actual upload/reactive paths; generated fixtures are temporary.
 for(ext in c("csv","rds")) {
  f<-tempfile(fileext=paste0(".",ext));df<-make_demo_data();df$risk_score<- -df$risk_score
  if(ext=="csv")readr::write_csv(df,f) else saveRDS(df,f)
  session$setInputs(data_source="upload",data_upload=data.frame(name=paste0("example.",ext),size=file.info(f)$size,type="",datapath=f))
  check(!is_example_rv()&&nrow(raw_rv())==1200&&fails(result_2d()),paste(ext,"upload invalidates earlier data"))
  session$setInputs(update_1d=3+match(ext,c("csv","rds")))
  r<-result_1d();check(r$request$source=="uploaded_data"&&mean(r$table$auc)<.5,paste(ext,"upload actually computes reversed scores"))
  check(all(result_export(r)$source=="uploaded_data"),paste(ext,"download metadata identifies uploaded source"))
  unlink(f)
 }
 f<-tempfile(fileext=".csv");df<-make_demo_data()
 df$calibrated_prob<-NA_real_
 df$calibrated_prob[which(df$event==1)[1]]<-.8
 df$calibrated_prob[which(df$event==0)[1]]<-.2
 readr::write_csv(df,f)
 session$setInputs(data_source="upload",data_upload=data.frame(name="sparse.csv",size=file.info(f)$size,type="",datapath=f),metric_1d="bacc",metric_2d="bacc")
 session$setInputs(update_1d=10,update_2d=10)
 for(dim in c("1d","2d")) {
  r<-if(dim=="1d")result_1d() else result_2d()
  check(nrow(r$table)>0&&all(is.na(r$table$bacc))&&all(r$table$threshold_n_obs<=2),paste(dim,"uploaded sparse probabilities cannot pass support minima"))
  check(grepl("No supported values",output[[paste0("result_status_",dim)]]),paste(dim,"shows missing metric support status"))
  check(all(!result_export(r)$threshold_supported)&&all(result_export(r)$source=="uploaded_data"),paste(dim,"export retains unsupported counts and upload provenance"))
 }
 unlink(f)
 # Sparse GAM support and a real fit failure must share plot/status gating.
 f<-tempfile(fileext=".csv");df<-make_demo_data();df$calibrated_prob[df$age<33]<-NA_real_;readr::write_csv(df,f)
 session$setInputs(data_upload=data.frame(name="sparse_surface.csv",size=file.info(f)$size,type="",datapath=f),metric_2d="bacc",method_2d="fixedrange",surface_mode_2d="gam",gam_k_2d=30)
 session$setInputs(update_2d=11)
 check(sum(is.finite(result_2d()$table$bacc))==7&&!prepared_plot_2d()$ok,"seven-cell GAM is unavailable while metric values are retained")
 check(grepl("at least 20",output$result_status_2d)&&nrow(result_export(result_2d()))>0,"GAM status explains cell minimum and CSV stays usable")
 session$setInputs(surface_mode_2d="raw");check(fails(prepared_plot_2d()),"changing surface mode invalidates prepared plot")
 session$setInputs(update_2d=12)
 check(prepared_plot_2d()$ok&&grepl("Current results",output$result_status_2d),"raw display remains available for sparse metric cells")
 session$setInputs(metric_2d="auc",surface_mode_2d="gam",gam_k_2d=200);session$setInputs(update_2d=13)
 check(sum(is.finite(result_2d()$table$auc))>=20&&!prepared_plot_2d()$ok&&grepl("Unable to fit",output$result_status_2d),"actual GAM fit failure is caught after passing the cell minimum")
 session$setInputs(gam_k_2d=10);session$setInputs(update_2d=14)
 check(prepared_plot_2d()$ok,"a valid GAM setting recovers after a failed fit")
 unlink(f)
 # Metadata follows the effective filter through uploads and remapping.
 f<-tempfile(fileext=".csv")
 df<-tibble(event=rep(0:1,50),risk_score=rep(c(-1,1),50),calibrated_prob=rep(c(.2,.8),50),age=seq_len(100))
 readr::write_csv(df,f)
 session$setInputs(data_upload=data.frame(name="no_groups.csv",size=file.info(f)$size,type="",datapath=f),col_cohort="",col_moderators="age",width_pct_1d=1,step_pct_1d=1,metric_1d="auc")
 session$setInputs(update_1d=20)
 r<-result_1d();meta<-jsonlite::fromJSON(result_export(r)$applied_settings_json[1])
 check(r$windowing_rows==100&&is.null(meta$mapping$group)&&is.null(meta$mapping$selected_groups)&&!meta$mapping$group_filter_applied,"upload without groups clears inactive selection in export metadata")
 df$group<-rep(c("New_A","New_B"),each=50);readr::write_csv(df,f)
 session$setInputs(data_upload=data.frame(name="new_groups.csv",size=file.info(f)$size,type="",datapath=f),col_cohort="group",cohorts=c("Cohort_A","New_B"))
 session$setInputs(update_1d=21)
 r<-result_1d();meta<-jsonlite::fromJSON(result_export(r)$applied_settings_json[1])
 check(r$windowing_rows==50&&all(r$table$n_obs==50&r$table$auc==1)&&identical(meta$mapping$selected_groups,"New_B")&&meta$mapping$group_filter_applied,"new upload exports only applied groups and computes the same subset")
 session$setInputs(cohorts=NULL);session$setInputs(update_1d=22)
 check(nrow(data_filtered())==0&&result_1d()$windowing_rows==0&&!prepared_plot_1d()$ok&&grepl("No groups selected",output$result_status_1d),"clearing an active group selection excludes all rows and explains the cause")
 session$setInputs(col_cohort="",cohorts="New_B");session$setInputs(update_1d=23)
 check(result_1d()$windowing_rows==100&&is.null(result_1d()$request$mapping$selected_groups),"removing the group mapping restores all rows without stale filter metadata")
 unlink(f)
 session$setInputs(data_upload=data.frame(name="bad.rds",size=1,type="",datapath="absent-file.rds"))
 check(is.null(raw_rv())&&is.null(dat_rv())&&fails(result_1d()),"failed upload cannot reuse earlier data")
})
shiny::testServer(server, {
 session$setInputs(data_source="example",col_outcome="event",col_score="risk_score",col_prob="calibrated_prob",col_cohort="group",col_moderators=c("age","cognitive_score"),col_outcome_negative="0",col_outcome_positive="1")
 s1<-result_defaults("1d");s1$moderator<-"age"
 s2<-result_defaults("2d");s2$x_var<-"age";s2$y_var<-"cognitive_score"
 do.call(session$setInputs,c(setNames(s1,paste0(names(s1),"_1d")),setNames(s2,paste0(names(s2),"_2d"))))
 session$setInputs(update_1d=1,update_2d=1)
 for(dim in c("1d","2d")) {
  r<-if(dim=="1d")result_1d() else result_2d()
  check(is.null(input$cohorts)&&r$windowing_rows==0&&grepl("No groups selected",output[[paste0("result_status_",dim)]]),paste(dim,"Update before group selector initialisation explains the empty selection"))
 }
})
cat("TOTAL",n,"server assertions passed\n")
