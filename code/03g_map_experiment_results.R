# Audit nested folds and summarize the exploratory MAP model comparison.
source(if(file.exists("code/setup.R")) "code/setup.R" else "setup.R")
out <- "tables/map_model_experiment"
paths <- file.path("models/map_model_experiment",sprintf("outer_%02d.rds",1:10))
stopifnot(all(file.exists(paths)))
x <- lapply(paths,readRDS)
d <- do.call(rbind,lapply(x,`[[`,"predictions"))
stopifnot(nrow(d)==92L,!anyDuplicated(d$site),all(is.finite(as.matrix(d[,-1]))))
truth <- read.csv("data/dat_site.csv");rownames(truth)<-truth$Site
stopifnot(setequal(d$site,truth$Site),max(abs(d$observed-log(truth[d$site,"map"])))<1e-8)
old <- read.csv("tables/pip_cv_site_uncertainty.csv");old<-old[old$target=="log_map",]
stopifnot(max(abs(d$pip-old$estimate[match(d$site,old$site)]))<1e-7)
old_lm <- read.csv("tables/loso_cv_site_predictions.csv")
old_lm <- old_lm[match(d$site,old_lm$site),]
stopifnot(max(abs(d$sp_lm-old_lm$lm_sp_site_impute_log_map))<1e-7,
          max(abs(d$site_lm-old_lm$lm_site_site_specimen_impute_log_map))<1e-7)
for(i in seq_along(x)){
 z <- x[[i]]
 stopifnot(all(z$predictions$fold==i),!length(intersect(z$train_sites,z$held_sites)),
           setequal(z$predictions$site,z$held_sites),setequal(c(z$train_sites,z$held_sites),d$site))
 ip <- z$inner_predictions
 for(j in 1:5){
  m <- z$inner_membership[[j]]
  stopifnot(!length(intersect(m$train_sites,m$held_sites)),
            !length(intersect(c(m$train_sites,m$held_sites),z$held_sites)),
            setequal(c(m$train_sites,m$held_sites),z$train_sites),
            setequal(ip$site[ip$inner==j],m$held_sites))
 }
 for(model in unique(ip$model)){
  v<-ip[ip$model==model,]
  stopifnot(nrow(v)==length(z$train_sites),!anyDuplicated(v$site),setequal(v$site,z$train_sites))
 }
 fit<-lm(observed~prediction,data=ip[ip$model=="pip",])
 stopifnot(max(abs(coef(fit)-z$calibration))<1e-10,
           max(abs(as.numeric(predict(fit,data.frame(prediction=z$predictions$pip)))-z$predictions$pip_recalibrated))<1e-10)
 for(k in seq_len(nrow(z$selections))){
  s<-z$selections[k,];eligible<-paste(s$level,z$configs$id[z$configs$method==s$method],sep="_")
  losses<-tapply((ip$prediction-ip$observed)^2,ip$model,mean)
  stopifnot(s$selected %in% eligible,abs(losses[s$selected]-min(losses[eligible]))<1e-12)
 }
}
# A train-only offset control separates correcting mean bias from stretching.
d$pip_offset_only <- d$pip
for(i in seq_along(x)) {
 ip <- x[[i]]$inner_predictions
 ip <- ip[ip$model=="pip",]
 idx <- d$fold==i
 d$pip_offset_only[idx] <- d$pip[idx]+mean(ip$observed-ip$prediction)
}
models<-setdiff(names(d),c("site","fold","observed"))
labels<-c(pip_offset_only="PIP + offset only",null="Mean only",pip="PIP",pip_recalibrated="PIP + recalibration",sp_lm="Species LM",sp_gam="Species GAM",sp_rf="Species RF",site_lm="Site LM",site_gam="Site GAM",site_rf="Site RF")
base_mse<-mean((d$pip-d$observed)^2)
metrics<-do.call(rbind,lapply(models,function(m){
 p<-d[[m]];e<-p-d$observed
 data.frame(model=m,label=labels[[m]],rmse=sqrt(mean(e^2)),mae=mean(abs(e)),bias=mean(e),
  correlation=cor(p,d$observed),prediction_sd=sd(p),observed_sd=sd(d$observed),spread_ratio=sd(p)/sd(d$observed),
  mse_improvement_vs_pip=1-mean(e^2)/base_mse)
}))
metrics<-metrics[order(metrics$rmse),];rownames(metrics)<-NULL
fold_metrics<-do.call(rbind,lapply(1:10,function(i)do.call(rbind,lapply(models,function(m){
 a<-d[d$fold==i,];data.frame(fold=i,model=m,n_sites=nrow(a),rmse=sqrt(mean((a[[m]]-a$observed)^2)))
}))))
calibration<-do.call(rbind,lapply(1:10,function(i)data.frame(fold=i,intercept=unname(x[[i]]$calibration[1]),slope=unname(x[[i]]$calibration[2]))))
selections<-do.call(rbind,lapply(x,`[[`,"selections"))
write.csv(d,file.path(out,"site_predictions.csv"),row.names=FALSE)
write.csv(metrics,file.path(out,"model_comparison.csv"),row.names=FALSE)
write.csv(fold_metrics,file.path(out,"fold_metrics.csv"),row.names=FALSE)
write.csv(calibration,file.path(out,"recalibration_coefficients.csv"),row.names=FALSE)
write.csv(selections,file.path(out,"selected_settings.csv"),row.names=FALSE)

# Every point is a site excluded from model fitting, tuning, and calibration.
plot_models<-function(selected,path,width,height){
 png(path,width=width,height=height,res=160)
 par(mfrow=if(length(selected)==4)c(2,2) else c(3,3),mar=c(4,4,3,1),oma=c(2,0,1,0))
 limits<-range(d[,c("observed",selected)])
 for(m in selected){
  score<-metrics[metrics$model==m,]
  plot(d$observed,d[[m]],xlim=limits,ylim=limits,pch=16,col="#48788baa",cex=.85,
       xlab="Observed ln(MAP in cm)",ylab="Predicted ln(MAP in cm)",
       main=sprintf("%s | RMSE %.3f",labels[[m]],score$rmse))
  abline(0,1,col="#999999",lty=2)
 }
 mtext("Nested site-grouped CV: all 92 sites held out; dashed line = perfect prediction",side=1,outer=TRUE,cex=.8,line=.4)
 dev.off()
}
bestgam<-metrics$model[metrics$model %in% c("sp_gam","site_gam")][1]
bestrf<-metrics$model[metrics$model %in% c("sp_rf","site_rf")][1]
plot_models(c("pip","pip_recalibrated",bestgam,bestrf),"plots/map_model_experiment_comparison.png",1500,1400)
plot_models(setdiff(models,"pip_offset_only"),"plots/map_model_experiment_all_models.png",1800,1800)
print(metrics,row.names=FALSE)
cat("All nested membership, prediction parity, recalibration, and selection checks passed.\n")
