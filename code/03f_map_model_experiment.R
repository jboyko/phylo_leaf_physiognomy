# Exploratory, nested site-grouped comparison of MAP models. Does not replace
# production models or fossil outputs. All tuning/recalibration is outer-train only.
if (nzchar(Sys.getenv("PIP_EXPERIMENT_R_LIB"))) .libPaths(c(Sys.getenv("PIP_EXPERIMENT_R_LIB"), .libPaths()))
if (nzchar(Sys.getenv("DILP_SOURCE"))) pkgload::load_all(Sys.getenv("DILP_SOURCE"), quiet=TRUE)
stopifnot(requireNamespace("mgcv",quietly=TRUE),requireNamespace("ranger",quietly=TRUE))
Sys.setenv(PIP_CV_SCHEME="ten_fold")
# Reuse exact production preprocessing and fold definitions, stopping before CV.
found_stop <- FALSE
for (expr in parse("code/03_loso_cv.R")) {
 if (is.call(expr) && identical(expr[[1]],as.name("<-")) && identical(expr[[2]],as.name("dilp_cv"))) {found_stop <- TRUE;break}
 eval(expr,envir=globalenv())
}
stopifnot(found_stop, length(all_sites)==92L)
EXP_DIR <- "tables/map_model_experiment"
EXP_MODELS <- "models/map_model_experiment"
dir.create(EXP_DIR,recursive=TRUE,showWarnings=FALSE)
dir.create(EXP_MODELS,recursive=TRUE,showWarnings=FALSE)

prepare_split <- function(train_sites,test_sites,seed) {
 stopifnot(!length(intersect(train_sites,test_sites)))
 filled <- fill_tooth_traits(raw_dat[raw_dat$Site %in% train_sites,])
 sp <- agg_species(filled);sp <- sp[intersect(rownames(sp),rownames(full_vcv)),]
 pn <- active_pred_names(sp,fossil_traits)
 site <- agg_site_direct(filled)
 held <- fill_tooth_traits(raw_dat[raw_dat$Site %in% test_sites,])
 heldsp <- aggregate(held[,pn,drop=FALSE],list(Site=held$Site,species=held$genusSpecies),mean,na.rm=TRUE)
 heldsp <- nan_to_na(heldsp)
 heldsite <- agg_site_direct(held)
 # Same training-only bagged trait imputation as the existing benchmark.
 set.seed(seed); impsp <- preProcess(sp[,pn,drop=FALSE],method="bagImpute")
 set.seed(seed); impsite <- preProcess(site[,pn,drop=FALSE],method="bagImpute")
 trainsp <- predict(impsp,sp[,pn,drop=FALSE]); testsp <- predict(impsp,heldsp[,pn,drop=FALSE])
 trainsite <- predict(impsite,site[,pn,drop=FALSE]); testsite <- predict(impsite,heldsite[,pn,drop=FALSE])
 for(z in list(trainsp,testsp,trainsite,testsite))stopifnot(all(is.finite(as.matrix(z))))
 # Explicit train-only scaling and removal of constant predictors for smooths.
 make_level <- function(x,y,newx,sites){
  center <- colMeans(x); scale <- vapply(x,sd,numeric(1));keep <- is.finite(scale)&scale>1e-10
  transform <- function(a)as.data.frame(sweep(sweep(as.matrix(a[,keep,drop=FALSE]),2,center[keep]),2,scale[keep],"/"))
  list(x=transform(x),y=y,newx=transform(newx),sites=sites)
 }
 list(sp=make_level(trainsp,sp$log_map,testsp,heldsp$Site),
      site=make_level(trainsite,site$log_map,testsite,heldsite$Site),
      rawsp=sp,pn=pn,heldsp=heldsp,testsp=testsp,seed=seed,
      train_sites=train_sites,test_sites=test_sites)
}

fit_pip <- function(d) {
 # Match production PGLS imputation, which includes the training response.
 dr <- d$rawsp[,c("log_map",d$pn),drop=FALSE]
 set.seed(d$seed); imp <- preProcess(dr,method="bagImpute");df <- predict(imp,dr)
 stopifnot(all(is.finite(as.matrix(df))))
 ids <- rownames(df);V <- full_vcv[ids,ids,drop=FALSE]
 form <- reformulate(d$pn,response="log_map")
 fit <- pglmEstLambda(formula=form,data=df,phylomat=V)
 X <- model.matrix(form,df);beta <- t(as.matrix(coef(fit)));vars <- intersect(colnames(X),rownames(beta))
 X <- X[,vars,drop=FALSE];beta <- beta[vars,,drop=FALSE]
 Vl <- V*fit$lambda;diag(Vl)<-diag(V)
 list(beta=beta,lambda=fit$lambda,K_train=solve(Vl),
      epsilon=setNames(df$log_map-as.numeric(X%*%beta),ids),sp_fit=ids,
      formula=form,common_vars=vars)
}
predict_pip <- function(pc,d){
 new <- d$testsp;new$log_map <- 0
 X <- model.matrix(pc$formula,new)[,rownames(pc$beta),drop=FALSE]
 pred <- as.numeric(X%*%pc$beta)
 ids <- d$heldsp$species;represented <- ids %in% rownames(full_vcv)
 C <- full_vcv[ids[represented],pc$sp_fit,drop=FALSE]*pc$lambda
 pred[represented] <- pred[represented]+as.numeric(C%*%pc$K_train%*%pc$epsilon[pc$sp_fit])
 tapply(pred,d$heldsp$Site,mean)
}
configs <- data.frame(id=c("gam_k3","gam_k4","rf_small","rf_wide","rf_smooth_small","rf_smooth_wide"),
 method=c("gam","gam",rep("rf",4)),k=c(3,4,rep(NA,4)),
 mtry=c(NA,NA,3,8,3,8),leaf=c(NA,NA,3,3,10,10))
fit_candidate <- function(d,cfg,seed,trees=400L){
 train <- d$x;train$y <- d$y
 if(cfg$method=="gam"){
  terms <- vapply(names(d$x),function(nm){
   k <- min(cfg$k,length(unique(d$x[[nm]])))
   if(k<3)nm else sprintf("s(%s,k=%d,bs='ts')",nm,k)
  },character(1))
  form <- reformulate(terms,response="y",env=asNamespace("mgcv"))
  fit <- mgcv::gam(form,data=train,method="REML",gamma=1.4)
  pred <- as.numeric(predict(fit,newdata=d$newx))
  details <- list(edf=sum(fit$edf))
 }else{
  fit <- ranger::ranger(y~.,data=train,num.trees=trees,mtry=min(cfg$mtry,ncol(d$x)),
                       min.node.size=cfg$leaf,seed=seed,num.threads=1)
  pred <- predict(fit,data=d$newx)$predictions
  details <- list(edf=NA_real_)
 }
 stopifnot(all(is.finite(pred)))
 list(pred=tapply(pred,d$sites,mean),details=details)
}

fold_arg <- Sys.getenv("PIP_MAP_FOLDS","")
outer_folds <- if(nzchar(fold_arg))as.integer(strsplit(fold_arg,",",fixed=TRUE)[[1]]) else 1:10
stopifnot(!anyNA(outer_folds),!anyDuplicated(outer_folds),all(outer_folds %in% 1:10))
for(outer in outer_folds){
 started <- proc.time()["elapsed"]
 archive <- readRDS(sprintf("models/loso_cv_fold_%02d.rds",outer))
 train_sites <- archive$train_sites;held_sites <- archive$held_sites
 stopifnot(setequal(held_sites,names(fold_assignment)[fold_assignment==outer]))
 cat("\nOUTER",outer,"training sites",length(train_sites),"held",length(held_sites),"\n")
 set.seed(7400+outer);inner_assignment <- setNames(sample(rep(1:5,length.out=length(train_sites))),train_sites)
 inner_predictions <- list();membership <- list()
 for(inner in 1:5){
  itest <- names(inner_assignment)[inner_assignment==inner];itrain <- setdiff(train_sites,itest)
  stopifnot(!length(intersect(c(itrain,itest),held_sites)))
  seed <- 17000+outer*10+inner
  d <- prepare_split(itrain,itest,seed)
  cat(" outer",outer,"inner",inner,"PIP\n")
  pip <- predict_pip(fit_pip(d),d)
  rows <- list(data.frame(site=names(pip),model="pip",prediction=as.numeric(pip)))
  for(level in c("sp","site"))for(i in seq_len(nrow(configs))){
   ans <- fit_candidate(d[[level]],configs[i,],seed)
   rows[[paste(level,i)]] <- data.frame(site=names(ans$pred),model=paste(level,configs$id[i],sep="_"),prediction=as.numeric(ans$pred))
  }
  ip <- do.call(rbind,rows);ip$observed <- dat_site_obs[ip$site,"log_map"];ip$inner<-inner
  inner_predictions[[inner]] <- ip
  membership[[inner]] <- list(train_sites=itrain,held_sites=itest,seed=seed)
  cat(" outer",outer,"inner",inner,"complete\n")
 }
 ip <- do.call(rbind,inner_predictions)
 inner_pip <- ip[ip$model=="pip",]
 stopifnot(nrow(inner_pip)==length(train_sites),setequal(inner_pip$site,train_sites),!anyDuplicated(inner_pip$site))
 calibration <- lm(observed~prediction,data=inner_pip)
 scores <- aggregate((ip$prediction-ip$observed)^2,list(model=ip$model),mean);names(scores)[2]<-"mse"
 outer_data <- prepare_split(train_sites,held_sites,42+outer)
 # Archive avoids refitting the already validated outer PIP model.
 pip <- predict_pip(archive$pgls_res$impute$log_map,outer_data)
 saved <- read.csv("tables/pip_cv_site_uncertainty.csv");saved<-saved[saved$fold==outer & saved$target=="log_map",]
 stopifnot(max(abs(pip[saved$site]-saved$estimate))<1e-7)
 predictions <- data.frame(site=names(pip),fold=outer,observed=dat_site_obs[names(pip),"log_map"],
                           null=mean(dat_site_obs[train_sites,"log_map"]),pip=as.numeric(pip))
 predictions$pip_recalibrated <- as.numeric(predict(calibration,newdata=data.frame(prediction=predictions$pip)))
 selections <- list()
 for(level in c("sp","site")){
  level_data <- outer_data[[level]]
  lmdata <- level_data$x;lmdata$y <- level_data$y
  lmpred <- tapply(predict(lm(y~.,data=lmdata),newdata=level_data$newx),level_data$sites,mean)
  predictions[[paste0(level,"_lm")]] <- lmpred[predictions$site]
  for(method in c("gam","rf")){
   options <- paste(level,configs$id[configs$method==method],sep="_")
   score <- scores[scores$model %in% options,];selected <- score$model[which.min(score$mse)]
   cfg <- configs[paste(level,configs$id,sep="_")==selected,]
   ans <- fit_candidate(level_data,cfg,42+outer,trees=800L)
   predictions[[paste(level,method,sep="_")]] <- ans$pred[predictions$site]
   selections[[paste(level,method)]] <- data.frame(fold=outer,level,method,selected,inner_rmse=sqrt(min(score$mse)),edf=ans$details$edf)
  }
 }
 stopifnot(all(is.finite(as.matrix(predictions[,-1]))))
 checkpoint <- list(predictions=predictions,inner_predictions=ip,selections=do.call(rbind,selections),
  calibration=coef(calibration),train_sites=train_sites,held_sites=held_sites,inner_membership=membership,
  inner_scores=scores,configs=configs,session=sessionInfo())
 saveRDS(checkpoint,file.path(EXP_MODELS,sprintf("outer_%02d.rds",outer)))
 cat("OUTER",outer,"DONE",round(proc.time()["elapsed"]-started,1),"seconds\n")
}
