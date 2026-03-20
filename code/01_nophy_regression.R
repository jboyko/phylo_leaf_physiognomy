setwd("/Users/jboyko/University\ of\ Michigan\ Dropbox/James\ Boyko/James\ Boyko’s\ files/Home/phylo_leaf_physiognomy")

# imports
library(ape)
library(caret)
library(glmnet)
library(ranger)
library(pdp)

phy <- read.tree("data/tre_pruned.tre")
dat <- read.csv("data/data_species.csv")
dat_site <- read.csv("data/dat_site.csv")
summary(dat)

# initial models
mat_model <- lm(mat ~ margin.score, data = dat)
summary(mat_model)

map_model <- lm(log(map) ~ ln.leaf.area.mm2, data = dat)
summary(map_model)

joint_model <- lm(cbind(mat = mat, map = log(map)) ~ margin.score + ln.leaf.area.mm2, data = dat)
summary(joint_model)

#inital predictions
predict(mat_model, newdata = data.frame(margin.score = 0.7))
predict(map_model, newdata = data.frame(margin.score = 0.7, ln.leaf.area.mm2 = 7.5))
predict(joint_model, newdata = data.frame(margin.score = 0.7, ln.leaf.area.mm2 = 7.5))

# more complex (non-phylogenetic) regressions
ctrl <- trainControl(method = "cv", 
  number = 10, 
  preProcOptions = list(thresh = 0.95),
  savePredictions = "final")

predictors <- dat[, 5:50]


exclude_vars <- c("Group.1", "Group.2", "Group.3", "Group.4", "site.info",
  "latitude", "longitude", "mat", "map", "gsmt", "gsp", 
  "gdd", "gsdd", "gsl", "cmmt", "wmmt", "mart")

vars_0 <- c(names(dat)[grep("teeth", names(dat))][-6], names(dat)[grep("tooth", names(dat))][-3])
vars_1 <- names(dat)[grep("ratio", names(dat))][5:6]

dat[,vars_0][is.na(dat[,vars_0])] <- 0
dat[,vars_1][is.na(dat[,vars_1])] <- 1

na_pct <- colSums(is.na(dat)) / nrow(dat)

clean_predictor_names <- names(dat)[na_pct < 0.40 & !(names(dat) %in% exclude_vars)]
predictors_clean <- dat[, clean_predictor_names]

dat$log_map <- log(dat$map)
target_vars <- c("mat", "log_map")

methods <- c(LM = "lm", ENet = "glmnet", RF = "ranger")
all_results <- list()

for (target in target_vars) {
  cat("\n--- Processing Target:", target, "---\n")
  all_results[[target]] <- list()
  for (m_label in names(methods)) {
    cat("Training", m_label, "...\n")
    importance_val <- if(methods[m_label] == "ranger") "permutation" else NULL
    tune_len <- if(methods[m_label] == "glmnet") 10 else 1
    all_results[[target]][[m_label]] <- train(
      x = predictors_clean, 
      y = dat[[target]],
      method = methods[m_label],
      trControl = ctrl,
      importance = importance_val,
      tuneLength = tune_len,
      preProcess = c("center", "scale", "bagImpute")
    )
  }
}
saveRDS(all_results, file = "models/nophy_models.rds")

# Compare MAT models
mat_comparison <- resamples(all_results$mat)
summary(mat_comparison)
dotplot(mat_comparison, main = "MAT Model Comparison")
mat_importance <- varImp(all_results$mat[["RF"]], scale = FALSE)
plot(mat_importance, top = 10, main = "Top 10 Predictors for MAT")

# Compare log(MAP) models
map_comparison <- resamples(all_results$log_map)
summary(map_comparison)
dotplot(map_comparison, main = "log(MAP) Model Comparison")
map_importance <- varImp(all_results$log_map[["RF"]], scale = FALSE)
plot(map_importance, top = 10, main = "Top 10 Predictors for MAP")

# Visualize the effect of margin.score on MAT while holding other variables constant
# partial(all_results$mat[["RF"]], pred.var = "evergreen", plot = TRUE)
# partial(all_results$mat[["RF"]], pred.var = "major.length.cm", plot = TRUE)
# partial(all_results$log_map[["RF"]], pred.var = "perim.area.cm2", plot = TRUE)
# partial(all_results$log_map[["RF"]], pred.var = "ln.leaf.area.mm2", plot = TRUE)

# Define a helper function to export resamples summaries
export_metrics <- function(comparison, target_name) {
  all_stats <- summary(comparison)$statistics
  for (metric in names(all_stats)) {
    file_path <- paste0("tables/", target_name, "_", metric, "_summary.csv")
    write.csv(all_stats[[metric]], file = file_path)
  }
  imp_data <- varImp(all_results[[target_name]][["RF"]], scale = FALSE)$importance
  imp_data <- imp_data[order(imp_data[,1], decreasing = TRUE), , drop = FALSE]
  write.csv(imp_data, file = paste0("tables/", target_name, "_RF_variable_importance.csv"))
}

export_metrics(mat_comparison, "mat")
export_metrics(map_comparison, "log_map")

library(gridExtra)
library(ggplot2)

export_pdp_grid <- function(target_name) {
  imp <- varImp(all_results[[target_name]][["RF"]], scale = FALSE)
  top_vars <- rownames(imp$importance)[order(imp$importance$Overall, decreasing = TRUE)][1:4]
  p_list <- lapply(top_vars, function(v) {
    pd_data <- partial(all_results[[target_name]][["RF"]], pred.var = v)
    ggplot(pd_data, aes(x = .data[[v]], y = .data$yhat)) +
      geom_line(linewidth = 0.8) +
      theme_minimal() +
      labs(title = v, 
        y = "Predicted Response (yhat)", 
        x = v)
  })
  g <- grid.arrange(grobs = p_list, ncol = 2, nrow = 2, 
    top = paste("Top 4 Drivers for", target_name))
  ggsave(paste0("plots/", target_name, "_top4_PDP_grid.png"), g, width = 10, height = 8)
}
export_pdp_grid("mat")
export_pdp_grid("log_map")

