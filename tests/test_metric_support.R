#!/usr/bin/env Rscript
file<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);root<-normalizePath(file.path(dirname(file),'..'))
suppressPackageStartupMessages({library(dplyr);library(tibble);library(purrr);library(ggplot2);library(pROC);library(patchwork)})
source(file.path(root,'R/engine.R'));source(file.path(root,'R/results.R'))
n<-0L;check<-function(x,label){stopifnot(isTRUE(x));n<<-n+1L;cat('PASS:',label,'\n')}
d<-tibble(EXP_LABEL=rep(c(0L,1L),50),Mean_Score=rep(c(-1,1),50),prob_work=rep(c(.1,.9),50),x=seq_len(100),y=rev(seq_len(100)))
methods<-list(
 fixedn=function(x,k=7)calc_metrics_trajectory_1d(x,'x',window_size=100,min_ev=k,min_nonev=k),
 fixedrange1d=function(x,k=7)calc_metrics_fixedrange_1d(x,'x',width_pct=1,step_pct=1,min_ev=k,min_nonev=k),
 knn=function(x,k=7)calc_grid_2d_knn(x,'x','y',grid_res=2,cell_n_mode='fixed',cell_n=100,min_ev=k,min_nonev=k,show_progress=FALSE),
 fixedrange2d=function(x,k=7)calc_grid_2d_fixedrange(x,'x','y',width_pct_x=1,step_pct_x=1,width_pct_y=1,step_pct_y=1,min_ev=k,min_nonev=k,show_progress=FALSE))
for(name in names(methods)) {
 fun<-methods[[name]]
 sparse<-d;sparse$prob_work[3:100]<-NA_real_;z<-fun(sparse)
 check(nrow(z)>0&&all(is.na(as.matrix(z[c('bacc','sens','spec')]))),paste(name,'one finite probability per class is unsupported'))
 check(all(z$threshold_n_obs==2&z$threshold_n_events==1&z$threshold_n_nonevents==1&!z$threshold_supported),paste(name,'exports actual probability denominators'))
 check(all(z$auc==1&z$auc_n_obs==100&z$auc_supported),paste(name,'missing probabilities do not suppress supported AUC'))
 sparse<-d;sparse$prob_work[15:100]<-NA_real_;z<-fun(sparse)
 check(all(z$bacc==1&z$sens==1&z$spec==1&z$threshold_n_obs==14&z$threshold_supported),paste(name,'exact seven-per-class threshold is supported'))
 check(all(is.na(fun(sparse,8)$bacc)),paste(name,'raising the support threshold masks the same observations'))
 sparse<-d;sparse$Mean_Score[3:100]<-NA_real_;z<-fun(sparse)
 check(all(is.na(z$auc)&z$auc_n_obs==2&!z$auc_supported),paste(name,'score support uses finite metric cases'))
 check(all(z$bacc==1&z$threshold_n_obs==100&z$threshold_supported),paste(name,'supported probability metric survives missing scores'))
 sparse<-d;sparse$Mean_Score[15:100]<-NA_real_;sparse$prob_work[c(1:14,29:100)]<-NA_real_;z<-fun(sparse)
 check(all(z$auc==1&z$bacc==1&z$auc_n_obs==14&z$threshold_n_obs==14),paste(name,'independent complete-case masks preserve supported metrics'))
 sparse<-d;sparse$prob_work[sparse$EXP_LABEL==0]<-NA_real_;z<-fun(sparse)
 check(all(is.na(z$bacc)&z$threshold_n_nonevents==0&!z$threshold_supported),paste(name,'one-class probability inputs remain unsupported'))
}
for(dim in c('1d','2d')) {
 sparse<-d;sparse$prob_work[3:100]<-NA_real_
 s<-result_defaults(dim);s$method<-'fixedrange';s$moderator<-'x';s$x_var<-'x';s$y_var<-'y';s$width_pct<-1;s$step_pct<-1;s$width_pct_x<-1;s$step_pct_x<-1;s$width_pct_y<-1;s$step_pct_y<-1;s$metric<-'bacc'
 r<-compute_result(list(data=sparse,settings=s,source='generated_test'),dim,FALSE)
 exported<-result_export(r);meta<-jsonlite::fromJSON(exported$applied_settings_json[1])
 check(r$eligible_rows==2&&r$eligible_events==1&&r$windowing_rows==100,paste(dim,'result separates metric eligibility from window geometry'))
 check(meta$eligible_rows_before_windowing==2&&meta$windowing_rows_before_metric_exclusions==100,paste(dim,'CSV metadata records the metric eligibility rule'))
 check(all(exported$threshold_n_obs==2&is.na(exported$bacc)),paste(dim,'unsupported metrics and denominators survive CSV export'))
 check(result_summary(r)$metric_row_memberships==sum(exported$threshold_n_obs),paste(dim,'summary uses metric-specific memberships'))
 check(inherits(tryCatch(plot_result(r),error=function(e)e),'error'),paste(dim,'unsupported selected metric has an explicit plotting error'))
}
bad<-d;bad$x[1]<-Inf;bad$y[2]<-NA;bad$EXP_LABEL[3]<-NA
check(nrow(window_input(bad,c('x','y')))==97,'nonfinite moderators and missing outcomes cannot enter window geometry')
empty<-d;empty$Mean_Score<-Inf;empty$prob_work<-NA_real_;z<-window_metrics(empty,7,7,.5)
check(z$auc_n_obs==0&&z$threshold_n_obs==0&&!z$auc_supported&&!z$threshold_supported,'nonfinite metric values never contribute support')
grid<-expand.grid(x_center=0:4,y_center=0:4);grid$bacc<-.55+.02*sin(grid$x_center+grid$y_center)
grid$n_obs<-100;grid$n_transition<-50;grid$threshold_n_obs<-20:44;grid$threshold_n_events<-10:34
fit<-gam_surface_from_grid(grid,metric='bacc',family='gaussian',k=10)$model
check(identical(as.numeric(fit$prior.weights),as.numeric(grid$threshold_n_obs)),'GAM display weights use metric-specific eligible counts')
cat('TOTAL',n,'metric-support assertions passed\n')
