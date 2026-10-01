# Exploratory application of the extant site recalibration to fossil MAP.
# Point estimates only: this does not validate fossil transfer or update intervals.
source(if(file.exists("code/setup.R")) "code/setup.R" else "setup.R")
cv <- read.csv("tables/pip_cv_site_uncertainty.csv")
cv <- cv[cv$target=="log_map",]
stopifnot(nrow(cv)==92L,!anyDuplicated(cv$site),all(is.finite(cv$estimate)),all(is.finite(cv$observed)))
# Final correction uses all extant out-of-fold site predictions. The preceding
# nested experiment assessed this fitting procedure; do not average its slopes.
correction <- lm(observed~estimate,data=cv)
a <- unname(coef(correction)[1]);b <- unname(coef(correction)[2])
rows <- lapply(c("formal_only","include_informal"),function(scenario){
 d <- read.csv(file.path("tables",paste0("fossil_site_uncertainty_",scenario,".csv")))
 d <- d[d$target=="log_map",]
 stopifnot(nrow(d)==10L,!anyDuplicated(d$site),all(is.finite(d$estimate)),
           all(abs(exp(d$estimate)-d$estimate_response)<1e-7))
 corrected <- a+b*d$estimate
 data.frame(site=d$site,age_ma=d$age_ma,taxonomy_scenario=scenario,n_species=d$n_species,
  original_log_map=d$estimate,recalibrated_log_map=corrected,
  original_map_cm=exp(d$estimate),recalibrated_map_cm=exp(corrected),
  change_cm=exp(corrected)-exp(d$estimate),percent_change=100*(exp(corrected-d$estimate)-1),
  within_extant_prediction_range=d$estimate>=min(cv$estimate)&d$estimate<=max(cv$estimate))
})
result <- do.call(rbind,rows)
result <- result[order(result$taxonomy_scenario,-result$age_ma),]
write.csv(result,"tables/fossil_map_recalibration_trial.csv",row.names=FALSE)
write.csv(data.frame(intercept=a,slope=b,unchanged_map_cm=exp(-a/(b-1)),
                    n_calibration_sites=nrow(cv),source="ten-fold out-of-fold extant PIP predictions"),
          "tables/fossil_map_recalibration_coefficients.csv",row.names=FALSE)
saveRDS(correction,"models/fossil_map_recalibration_trial.rds")
d <- result[result$taxonomy_scenario=="formal_only",];y<-rev(seq_len(nrow(d)))
png("plots/fossil_map_recalibration_trial.png",width=1550,height=1050,res=160)
par(mar=c(5,13,4,2),oma=c(2,0,0,0))
plot(d$original_map_cm,y,type="n",yaxt="n",ylab="",xlab="Mean annual precipitation (cm/year)",
     xlim=range(c(d$original_map_cm,d$recalibrated_map_cm))*c(.93,1.08),
     main="Fossil MAP: applying the extant recalibration")
abline(h=y,col="#eeeeee")
segments(d$original_map_cm,y,d$recalibrated_map_cm,y,col="#a0a0a0",lwd=2)
points(d$original_map_cm,y,pch=1,cex=1.2,lwd=2,col="#555555")
points(d$recalibrated_map_cm,y,pch=16,cex=1.1,col="#397d91")
axis(2,at=y,labels=d$site,las=1,tick=FALSE,cex.axis=.85)
legend("bottomright",c("Current PIP","Recalibrated PIP"),pch=c(1,16),col=c("#555555","#397d91"),bty="n")
mtext("Formal taxonomy; point estimates only. Transfer of this extant correction to fossils is unvalidated.",side=1,outer=TRUE,cex=.8,line=.5)
dev.off()
print(coef(correction));cat("Unchanged at",exp(-a/(b-1)),"cm/year\n")
print(d[,c("site","original_map_cm","recalibrated_map_cm","percent_change")],row.names=FALSE)
