# Module 5: Predictive Model — Corrected (features selected inside training fold)
suppressMessages({
  library(randomForest)
  library(pROC)
  library(glmnet)
  library(ggplot2)
})

workdir <- "~/ASD_multiomics"
outdir  <- file.path(workdir, "module5")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("===== Module 5: Predictive Model (corrected) =====\n\n")

# ====================================================================
# 1. LOAD DATA
# ====================================================================
mod1 <- readRDS(file.path(workdir, "module1/module1_results.rds"))
expr <- mod1$expr; meta <- mod1$meta

# Only blood samples
blood_idx <- meta$Tissue != "Brain"
blood_expr <- expr[, blood_idx]; blood_meta <- meta[blood_idx, ]

# Train = Prenatal (maternal) + Birth (cord); Test = Childhood (child)
train_idx <- blood_meta$Stage %in% c("Prenatal", "Birth")
test_idx  <- blood_meta$Stage == "Childhood"
train_x <- t(blood_expr[, train_idx, drop=FALSE])
train_y <- as.numeric(ifelse(blood_meta$Diagnosis[train_idx] == "ASD", 1, 0))
test_x  <- t(blood_expr[, test_idx, drop=FALSE])
test_y  <- ifelse(blood_meta$Diagnosis[test_idx] == "ASD", 1, 0)
cat(sprintf("Train: %d (Prenatal+Birth), Test: %d (Childhood)\n",
            nrow(train_x), nrow(test_x)))

# ====================================================================
# 2. FEATURE SELECTION (INSIDE TRAINING FOLD ONLY)
# ====================================================================
cat("\n2. Feature selection inside training fold...\n")

# Step A: Pre-filter top 500 by variance (in training data only)
train_var <- apply(train_x, 2, var)
top500 <- names(sort(train_var, decreasing = TRUE))[1:min(500, length(train_var))]
cat(sprintf("  Variance filter: %d genes\n", length(top500)))

# Step B: LASSO (10-fold CV, lambda.1se)
set.seed(123)
cv_lasso <- cv.glmnet(x = train_x[, top500, drop=FALSE], y = train_y,
                       family = "binomial", alpha = 1, nfolds = 10)
lasso_coef <- coef(cv_lasso, s = "lambda.1se")
lasso_genes <- setdiff(rownames(lasso_coef)[lasso_coef[, 1] != 0], "(Intercept)")
cat(sprintf("  LASSO (lambda.1se): %d genes\n", length(lasso_genes)))

if (length(lasso_genes) < 10) {
  lasso_coef <- coef(cv_lasso, s = "lambda.min")
  lasso_genes <- setdiff(rownames(lasso_coef)[lasso_coef[, 1] != 0], "(Intercept)")
  cat(sprintf("  LASSO (lambda.min): %d genes\n", length(lasso_genes)))
}

# Step C: Random Forest importance on training data
rf_select <- randomForest(x = train_x[, lasso_genes, drop=FALSE],
                           y = factor(train_y), ntree = 1000)
rf_imp <- importance(rf_select)[, "MeanDecreaseGini"]
rf_imp <- sort(rf_imp, decreasing = TRUE)
features <- names(rf_imp)[1:min(20, length(rf_imp))]
cat(sprintf("  RF top features: %d\n", length(features)))
cat(sprintf("  Genes: %s\n", paste(features, collapse = ", ")))

# Save feature importance
write.csv(data.frame(Gene = names(rf_imp), Importance = rf_imp),
          file.path(outdir, "feature_importance.csv"), row.names = FALSE)

# ====================================================================
# 3. TRAIN FINAL MODEL
# ====================================================================
cat("\n3. Training final RF model...\n")
set.seed(789)
rf_final <- randomForest(
  x = train_x[, features, drop = FALSE],
  y = factor(train_y),
  ntree = 2000,
  mtry = floor(sqrt(length(features)))
)
oob_err <- rf_final$err.rate[2000, 1]
cat(sprintf("  OOB error rate: %.4f\n", oob_err))

# ====================================================================
# 4. TEST ON CHILDHOOD BLOOD (independent)
# ====================================================================
cat("\n4. Testing on childhood blood...\n")
test_pred <- predict(rf_final, test_x[, features, drop = FALSE], type = "prob")[, 2]
test_roc <- roc(test_y, test_pred, quiet = FALSE)
test_auc <- auc(test_roc)
cat(sprintf("  Test AUC: %.4f (95%% CI: %.4f-%.4f)\n",
            test_auc, ci.auc(test_roc)[1], ci.auc(test_roc)[3]))

pred_class <- ifelse(test_pred > 0.5, 1, 0)
cm <- table(Predicted = pred_class, Actual = test_y)
sensitivity <- cm[2, 2] / sum(cm[, 2])
specificity <- cm[1, 1] / sum(cm[, 1])
cat(sprintf("  Sensitivity: %.3f, Specificity: %.3f\n", sensitivity, specificity))
cat("  Confusion matrix:\n"); print(cm)

# ====================================================================
# 5. 10-FOLD CV (WITHIN TRAINING DATA)
# ====================================================================
cat("\n5. 10-fold CV (within training)...\n")
set.seed(101)
nfolds <- 10
folds <- sample(rep(1:nfolds, length.out = nrow(train_x)))
cv_aucs <- numeric(nfolds)
for (k in 1:nfolds) {
  fold_test <- which(folds == k)
  fold_train <- which(folds != k)

  # Feature selection inside each fold (lightweight: just RF on pre-filtered genes)
  fold_var <- apply(train_x[fold_train, , drop=FALSE], 2, var)
  fold_top <- names(sort(fold_var, decreasing = TRUE))[1:min(500, length(fold_var))]

  set.seed(101 + k)
  fold_cv <- cv.glmnet(x = train_x[fold_train, fold_top, drop=FALSE],
                        y = train_y[fold_train], family = "binomial", nfolds = 5)
  fold_lg <- setdiff(rownames(coef(fold_cv, s = "lambda.1se"))[coef(fold_cv, s = "lambda.1se")[, 1] != 0], "(Intercept)")
  if (length(fold_lg) < 5) {
    fold_lg <- setdiff(rownames(coef(fold_cv, s = "lambda.min"))[coef(fold_cv, s = "lambda.min")[, 1] != 0], "(Intercept)")
  }

  rf_cv <- randomForest(x = train_x[fold_train, fold_lg, drop=FALSE],
                         y = factor(train_y[fold_train]), ntree = 500)
  fold_imp <- importance(rf_cv)[, "MeanDecreaseGini"]
  fold_feat <- names(sort(fold_imp, decreasing = TRUE))[1:min(20, length(fold_imp))]

  rf_cv2 <- randomForest(x = train_x[fold_train, fold_feat, drop=FALSE],
                          y = factor(train_y[fold_train]), ntree = 500)
  cv_pred <- predict(rf_cv2, train_x[fold_test, fold_feat, drop=FALSE], type = "prob")[, 2]
  cv_aucs[k] <- auc(roc(train_y[fold_test], cv_pred, quiet = TRUE))
}
cat(sprintf("  Nested CV AUC: %.4f (SD: %.4f)\n", mean(cv_aucs), sd(cv_aucs)))

# ====================================================================
# 6. SUBTYPE-SPECIFIC MODELS
# ====================================================================
cat("\n6. Subtype-specific models...\n")
subtype_results <- list()

mod4 <- tryCatch(readRDS(file.path(workdir, "module4/module4_results.rds")), error = function(e) NULL)
if (!is.null(mod4)) {
  st_df <- mod4$subtypes
  st_names <- st_df$Subtype_Name_K2
  names(st_names) <- st_df$Sample

  for (st in unique(st_names)) {
    st_samples <- names(st_names)[st_names == st]
    st_train <- intersect(st_samples, colnames(blood_expr)[train_idx])
    ctrl_train <- colnames(blood_expr)[train_idx][blood_meta$Diagnosis[train_idx] == "Control"]

    if (length(st_train) >= 5) {
      # Binary: subtype vs Control, features from inside this training
      sub_samples <- c(st_train, ctrl_train)
      sub_x <- t(blood_expr[features, sub_samples, drop=FALSE])
      sub_y <- factor(c(rep(1, length(st_train)), rep(0, length(ctrl_train))))

      sub_rf <- randomForest(x = sub_x, y = sub_y, ntree = 500)

      # Test
      st_test <- intersect(st_samples, colnames(blood_expr)[test_idx])
      ctrl_test <- colnames(blood_expr)[test_idx][blood_meta$Diagnosis[test_idx] == "Control"]
      sub_test_x <- t(blood_expr[features, c(st_test, ctrl_test), drop=FALSE])
      sub_test_y <- c(rep(1, length(st_test)), rep(0, length(ctrl_test)))
      sub_pred <- predict(sub_rf, sub_test_x, type = "prob")[, 2]
      sub_auc <- auc(roc(sub_test_y, sub_pred, quiet = TRUE))
      cat(sprintf("  %s: AUC=%.4f (train=%d, test=%d)\n",
                  st, sub_auc, length(st_train), length(st_test)))
      subtype_results[[st]] <- list(auc = sub_auc)
    }
  }
}

# ====================================================================
# 7. SAVE
# ====================================================================
cat("\n7. Saving...\n")
saveRDS(list(
  features = features,
  model = rf_final,
  test_auc = test_auc,
  test_roc = test_roc,
  sensitivity = sensitivity,
  specificity = specificity,
  cv_aucs = cv_aucs,
  oob_error = oob_err,
  subtype_results = subtype_results,
  confusion_matrix = cm
), file.path(outdir, "module5_results.rds"))

write.csv(data.frame(Gene = features), file.path(outdir, "predictive_model_genes.csv"), row.names = FALSE)

sink(file.path(outdir, "module5_summary.txt"))
cat(sprintf("Module 5: Predictive Model (corrected)\nDate: %s\n\n", Sys.Date()))
cat(sprintf("Training: %d samples (Prenatal+Birth)\n", nrow(train_x)))
cat(sprintf("Testing: %d samples (Childhood)\n", nrow(test_x)))
cat(sprintf("Feature selection: variance filter -> LASSO(10-fold CV) -> RF importance\n"))
cat(sprintf("Features: %d genes\n", length(features)))
cat(sprintf("Feature genes: %s\n\n", paste(features, collapse = ", ")))
cat(sprintf("OOB error: %.4f\n", oob_err))
cat(sprintf("Nested CV AUC: %.4f +/- %.4f\n", mean(cv_aucs), sd(cv_aucs)))
cat(sprintf("Test AUC: %.4f (95%% CI: %.4f-%.4f)\n",
            test_auc, ci.auc(test_roc)[1], ci.auc(test_roc)[3]))
cat(sprintf("Sensitivity: %.3f\n", sensitivity))
cat(sprintf("Specificity: %.3f\n", specificity))
cat(sprintf("\nConfusion matrix:\n"))
print(cm)
if (length(subtype_results) > 0) {
  cat("\nSubtype-specific:\n")
  for (nm in names(subtype_results)) {
    cat(sprintf("  %s: AUC=%.4f\n", nm, subtype_results[[nm]]$auc))
  }
}
sink()

cat(sprintf("\n===== Module 5 DONE: Test AUC=%.4f (nested CV=%.4f) =====\n",
            test_auc, mean(cv_aucs)))
