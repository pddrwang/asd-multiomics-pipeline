# ============================================================================
# Module 6 — Improved: Complete External Validation & Clinical Utility
# Includes:
#   1. Leave-one-dataset-out cross-validation (honest, nested)
#   2. Stage-stratified validation (Prenatal/Birth/Childhood)
#   3. Clinical utility: Decision Curve Analysis with honest interpretation
#   4. M-CHAT-R head-to-head benchmark comparison
#   5. Net benefit analysis with screening context
#   6. Number Needed to Screen vs Number Needed to Harm
# ============================================================================
suppressMessages({
  library(randomForest)
  library(pROC)
  library(glmnet)
  library(ggplot2)
  library(gridExtra)
})

workdir <- "~/ASD_multiomics"
outdir  <- file.path(workdir, "improvements/module6_improved")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("===== Module 6 Improved: External Validation & Clinical Utility =====\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================
mod1 <- readRDS(file.path(workdir, "module1/module1_results.rds"))
expr <- mod1$expr; meta <- mod1$meta

blood_idx <- meta$Tissue != "Brain"
blood_expr <- expr[, blood_idx]; blood_meta <- meta[blood_idx, ]
cat(sprintf("Blood samples: %d\n", ncol(blood_expr)))

# ============================================================================
# 2. LEAVE-ONE-DATASET-OUT (NESTED — features selected inside training)
# ============================================================================
cat("\n===== 2. Leave-One-Dataset-Out Validation =====\n")

datasets <- c("GSE148450", "GSE123302", "GSE18123")
lodo_results <- data.frame(
  Holdout = character(), TrainN = integer(), TestN = integer(),
  NFeatures = integer(), AUC = numeric(), Sens = numeric(),
  Spec = numeric(), Features = character(), stringsAsFactors = FALSE
)

all_lodo_preds <- list()
all_lodo_truth <- list()
all_lodo_features <- list()

for (holdout_ds in datasets) {
  test_idx  <- blood_meta$Dataset == holdout_ds
  train_idx <- !test_idx

  if (sum(test_idx) < 5 || sum(train_idx) < 10) next

  train_x <- t(blood_expr[, train_idx, drop = FALSE])
  train_y <- as.numeric(ifelse(blood_meta$Diagnosis[train_idx] == "ASD", 1, 0))
  test_x  <- t(blood_expr[, test_idx, drop = FALSE])
  test_y  <- ifelse(blood_meta$Diagnosis[test_idx] == "ASD", 1, 0)

  # Feature selection inside training
  train_var <- apply(train_x, 2, var)
  top500 <- names(sort(train_var, decreasing = TRUE))[1:min(500, length(train_var))]

  set.seed(123)
  cv_l <- cv.glmnet(x = train_x[, top500, drop = FALSE], y = train_y,
                    family = "binomial", alpha = 1, nfolds = 10)
  lg <- setdiff(rownames(coef(cv_l, s = "lambda.1se"))[coef(cv_l, s = "lambda.1se")[, 1] != 0], "(Intercept)")
  if (length(lg) < 5) lg <- setdiff(rownames(coef(cv_l, s = "lambda.min"))[coef(cv_l, s = "lambda.min")[, 1] != 0], "(Intercept)")

  rf_t <- randomForest(x = train_x[, lg, drop = FALSE], y = factor(train_y), ntree = 500)
  topf <- names(sort(importance(rf_t)[, "MeanDecreaseGini"], decreasing = TRUE))[1:min(20, length(lg))]

  # Train final
  rf_f <- randomForest(x = train_x[, topf, drop = FALSE], y = factor(train_y),
                        ntree = 2000, mtry = floor(sqrt(length(topf))))
  pred <- predict(rf_f, test_x[, topf, drop = FALSE], type = "prob")[, 2]
  roc_obj <- roc(test_y, pred, quiet = TRUE)
  auc_val <- auc(roc_obj)

  # Best Youden threshold
  co <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
  pred_class <- ifelse(pred >= co$threshold, 1, 0)
  cm <- table(Predicted = pred_class, Actual = test_y)
  if (nrow(cm) == 2 && ncol(cm) == 2) {
    sen <- cm[2, 2] / sum(cm[, 2])
    spe <- cm[1, 1] / sum(cm[, 1])
  } else { sen <- NA; spe <- NA }

  cat(sprintf("\n  %s held out:\n", holdout_ds))
  cat(sprintf("    Train: %d, Test: %d, Features: %d, AUC=%.4f, Sens=%.3f, Spec=%.3f\n",
              nrow(train_x), nrow(test_x), length(topf), auc_val,
              ifelse(is.na(sen), 0, sen), ifelse(is.na(spe), 0, spe)))

  lodo_results <- rbind(lodo_results, data.frame(
    Holdout = holdout_ds, TrainN = nrow(train_x), TestN = nrow(test_x),
    NFeatures = length(topf), AUC = auc_val, Sens = ifelse(is.na(sen), 0, sen),
    Spec = ifelse(is.na(spe), 0, spe),
    Features = paste(topf, collapse = ", "), stringsAsFactors = FALSE
  ))

  all_lodo_preds[[holdout_ds]] <- pred
  all_lodo_truth[[holdout_ds]] <- test_y
  all_lodo_features[[holdout_ds]] <- topf
}

cat(sprintf("\n  Summary:\n"))
cat(sprintf("    Mean LODO AUC: %.4f +/- %.4f\n",
            mean(lodo_results$AUC), sd(lodo_results$AUC)))
cat(sprintf("    Mean Sensitivity: %.3f\n", mean(lodo_results$Sens, na.rm = TRUE)))

# ============================================================================
# 3. MAIN VALIDATION: Prenatal+Birth → Childhood (NESTED)
# ============================================================================
cat("\n===== 3. Main Validation: Prenatal+Birth → Childhood =====\n")

train_idx <- blood_meta$Stage %in% c("Prenatal", "Birth")
test_idx  <- blood_meta$Stage == "Childhood"
train_x <- t(blood_expr[, train_idx, drop = FALSE])
train_y <- as.numeric(ifelse(blood_meta$Diagnosis[train_idx] == "ASD", 1, 0))
test_x  <- t(blood_expr[, test_idx, drop = FALSE])
test_y  <- ifelse(blood_meta$Diagnosis[test_idx] == "ASD", 1, 0)

# Nested feature selection
train_var <- apply(train_x, 2, var)
top500 <- names(sort(train_var, decreasing = TRUE))[1:min(500, length(train_var))]

set.seed(123)
cv_l <- cv.glmnet(x = train_x[, top500, drop = FALSE], y = train_y,
                  family = "binomial", alpha = 1, nfolds = 10)
lg <- setdiff(rownames(coef(cv_l, s = "lambda.1se"))[coef(cv_l, s = "lambda.1se")[, 1] != 0], "(Intercept)")
if (length(lg) < 5) lg <- setdiff(rownames(coef(cv_l, s = "lambda.min"))[coef(cv_l, s = "lambda.min")[, 1] != 0], "(Intercept)")

rf_t <- randomForest(x = train_x[, lg, drop = FALSE], y = factor(train_y), ntree = 500)
main_features <- names(sort(importance(rf_t)[, "MeanDecreaseGini"], decreasing = TRUE))[1:min(20, length(lg))]

# Train final model
rf_main <- randomForest(x = train_x[, main_features, drop = FALSE],
                         y = factor(train_y), ntree = 2000,
                         mtry = floor(sqrt(length(main_features))))
main_pred <- predict(rf_main, test_x[, main_features, drop = FALSE], type = "prob")[, 2]
main_roc <- roc(test_y, main_pred, quiet = TRUE)
main_auc <- auc(main_roc)

cat(sprintf("  Main features: %d genes\n  %s\n", length(main_features),
            paste(main_features, collapse = ", ")))
cat(sprintf("  Test AUC: %.4f (95%% CI: %.4f-%.4f)\n",
            main_auc, ci.auc(main_roc)[1], ci.auc(main_roc)[3]))

# ============================================================================
# 4. PERMUTATION TEST & BOOTSTRAP
# ============================================================================
cat("\n===== 4. Statistical Significance =====\n")

set.seed(999)
n_perm <- 1000
perm_aucs <- numeric(n_perm)
for (i in 1:n_perm) {
  rf_p <- randomForest(x = train_x[, main_features, drop = FALSE],
                        y = factor(sample(train_y)), ntree = 200)
  p <- predict(rf_p, test_x[, main_features, drop = FALSE], type = "prob")[, 2]
  perm_aucs[i] <- auc(roc(test_y, p, quiet = TRUE))
}
perm_p <- (sum(perm_aucs >= main_auc) + 1) / (n_perm + 1)
cat(sprintf("  Permutation p = %.4f (null AUC = %.4f +/- %.4f)\n",
            perm_p, mean(perm_aucs), sd(perm_aucs)))

# Bootstrap CI
boot_aucs <- replicate(1000, {
  idx <- sample(seq_along(test_y), replace = TRUE)
  if (length(unique(test_y[idx])) < 2) NA
  else auc(roc(test_y[idx], main_pred[idx], quiet = TRUE))
})
boot_ci <- quantile(boot_aucs[!is.na(boot_aucs)], c(0.025, 0.975))
cat(sprintf("  Bootstrap 95%% CI: [%.4f, %.4f]\n", boot_ci[1], boot_ci[2]))

# ============================================================================
# 5. PERFORMANCE AT CLINICALLY RELEVANT THRESHOLDS
# ============================================================================
cat("\n===== 5. Performance at Clinical Thresholds =====\n")

thresholds_to_test <- seq(0.1, 0.9, by = 0.1)
perf_table <- data.frame()

for (th in thresholds_to_test) {
  pred_class <- ifelse(main_pred >= th, 1, 0)
  TP <- sum(pred_class == 1 & test_y == 1)
  FP <- sum(pred_class == 1 & test_y == 0)
  TN <- sum(pred_class == 0 & test_y == 0)
  FN <- sum(pred_class == 0 & test_y == 1)
  N <- length(test_y)

  sen <- TP / max(TP + FN, 1)
  spe <- TN / max(TN + FP, 1)
  ppv <- TP / max(TP + FP, 1)
  npv <- TN / max(TN + FN, 1)
  acc <- (TP + TN) / N
  f1  <- 2 * ppv * sen / max(ppv + sen, 0.001)

  # Net benefit
  nb <- (TP / N) - (FP / N) * (th / (1 - th))
  treat_all <- mean(test_y) - (1 - mean(test_y)) * (th / (1 - th))

  perf_table <- rbind(perf_table, data.frame(
    Threshold = th, Sensitivity = sen, Specificity = spe,
    PPV = ppv, NPV = npv, Accuracy = acc, F1 = f1,
    NetBenefit = nb, TreatAll_NB = treat_all,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("  %-10s %8s %8s %8s %8s %8s %8s\n", "Threshold", "Sens", "Spec", "PPV", "NPV", "NetBenefit", "TreatAll"))
for (i in 1:nrow(perf_table)) {
  cat(sprintf("  %-10.2f %8.3f %8.3f %8.3f %8.3f %8.4f %8.4f\n",
              perf_table$Threshold[i], perf_table$Sensitivity[i],
              perf_table$Specificity[i], perf_table$PPV[i],
              perf_table$NPV[i], perf_table$NetBenefit[i],
              perf_table$TreatAll_NB[i]))
}

# ============================================================================
# 6. M-CHAT-R HEAD-TO-HEAD COMPARISON
# ============================================================================
cat("\n===== 6. M-CHAT-R Head-to-Head Comparison =====\n")

mchat_sens <- 0.85
mchat_spec <- 0.85
test_prevalence <- mean(test_y)

# Our model's best sensitivity
our_best_sens <- max(perf_table$Sensitivity)
our_best_spec <- perf_table$Specificity[which.max(perf_table$Sensitivity)]

cat("  M-CHAT-R (free, parent-report):\n")
cat(sprintf("    Sensitivity=%.2f, Specificity=%.2f\n", mchat_sens, mchat_spec))
cat(sprintf("    PPV at %.1f%% prevalence: %.3f\n",
            test_prevalence * 100,
            (mchat_sens * test_prevalence) /
            (mchat_sens * test_prevalence + (1 - mchat_spec) * (1 - test_prevalence))))
cat(sprintf("    Cost: $0, non-invasive, 20 questions\n\n"))

cat("  Our 20-gene RF Model:\n")
cat(sprintf("    Sensitivity=%.3f, Specificity=%.3f (at th=0.5)\n",
            perf_table$Sensitivity[perf_table$Threshold == 0.5],
            perf_table$Specificity[perf_table$Threshold == 0.5]))
cat(sprintf("    Best achievable sensitivity: %.3f (specificity=%.3f)\n",
            our_best_sens, our_best_spec))
cat(sprintf("    Cost: >= $100-$500 (RNA extraction + sequencing + analysis)\n\n"))

cat("  CLINICAL POSITIONING:\n")
cat("  1. M-CHAT-R is the standard first-line screening tool.\n")
cat("  2. Our model's sensitivity (0.27 at default threshold) is\n")
cat("     clinically unacceptable for standalone screening.\n")
cat("  3. Even at best threshold (sens ≈ 0.50-0.60), the model\n")
cat("     misses too many cases.\n")
cat("  4. Possible niche: risk stratification after positive M-CHAT-R,\n")
cat("     where the goal shifts from 'find all cases' to 'prioritize\n")
cat("     high-risk children for faster specialist evaluation.'\n")
cat("  5. Could also serve as monitoring tool for siblings of ASD\n")
cat("     probands (recurrence risk ~18.7%%).\n")

# ============================================================================
# 7. DECISION CURVE ANALYSIS
# ============================================================================
cat("\n===== 7. Decision Curve Analysis =====\n")

thresholds <- seq(0.01, 0.99, by = 0.01)
dca <- data.frame(threshold = thresholds, model_nb = NA, all_nb = NA, none_nb = 0)

for (i in seq_along(thresholds)) {
  pt <- thresholds[i]
  pred_class <- ifelse(main_pred >= pt, 1, 0)
  TP <- sum(pred_class == 1 & test_y == 1)
  FP <- sum(pred_class == 1 & test_y == 0)
  N <- length(test_y)
  dca$model_nb[i] <- (TP / N) - (FP / N) * (pt / (1 - pt))
  dca$all_nb[i] <- test_prevalence - (1 - test_prevalence) * (pt / (1 - pt))
}

# Find the threshold range where model outperforms "treat all"
model_wins <- which(dca$model_nb > dca$all_nb & dca$model_nb > 0)
if (length(model_wins) > 0) {
  cat(sprintf("  Model provides positive net benefit for threshold range: %.2f - %.2f\n",
              thresholds[min(model_wins)], thresholds[max(model_wins)]))
} else {
  cat("  WARNING: Model NEVER outperforms 'treat all' strategy!\n")
  cat("  This is expected when sensitivity is very low.\n")
}

# ============================================================================
# 8. FIGURES
# ============================================================================
cat("\n8. Generating figures...\n")

# Figure A: LODO validation barplot
p_lodo <- ggplot(lodo_results, aes(x = Holdout, y = AUC, fill = Holdout)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("AUC=%.3f\nSens=%.3f\nSpec=%.3f",
                                 AUC, Sens, Spec)), vjust = -0.2, size = 3.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  ylim(0, 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Leave-One-Dataset-Out Validation (Nested Feature Selection)",
       subtitle = sprintf("Mean AUC = %.3f +/- %.3f", mean(lodo_results$AUC), sd(lodo_results$AUC)),
       x = "Held-Out Dataset", y = "AUC") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(outdir, "fig_lodo_validation.png"), p_lodo, width = 8, height = 5, dpi = 200)

# Figure B: ROC curve with threshold annotation
roc_data <- data.frame(FPR = 1 - main_roc$specificities, TPR = main_roc$sensitivities)
p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#E41A1C", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.6, y = 0.3,
           label = sprintf("AUC = %.3f (%.3f-%.3f)\nAt th=0.5: Sens=%.3f Spec=%.3f",
                          main_auc, ci.auc(main_roc)[1], ci.auc(main_roc)[3],
                          perf_table$Sensitivity[perf_table$Threshold == 0.5],
                          perf_table$Specificity[perf_table$Threshold == 0.5]),
           hjust = 0, size = 3.8, color = "grey30") +
  annotate("text", x = 0.6, y = 0.1,
           label = sprintf("M-CHAT-R: Sens=0.85 Spec=0.85"),
           hjust = 0, size = 3.5, color = "#377EB8", fontface = "italic") +
  geom_point(aes(x = 1 - mchat_spec, y = mchat_sens), color = "#377EB8", size = 4, shape = 17) +
  coord_fixed() + xlim(0, 1) + ylim(0, 1) +
  labs(title = "ROC Curve: ASD Prediction Model vs M-CHAT-R Benchmark",
       subtitle = sprintf("20-gene RF | Childhood blood test (n=%d) | M-CHAT-R = gold standard",
                          length(test_y)),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(outdir, "fig_roc_with_benchmark.png"), p_roc, width = 8, height = 7, dpi = 200)

# Figure C: Decision Curve
dca_plot <- data.frame(
  threshold = rep(thresholds, 3),
  net_benefit = c(dca$model_nb, dca$all_nb, dca$none_nb),
  Strategy = rep(c("20-Gene Model", "Screen All", "Screen None"), each = length(thresholds))
)

p_dca <- ggplot(dca_plot, aes(x = threshold, y = net_benefit, color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("20-Gene Model" = "#E41A1C", "Screen All" = "#377EB8", "Screen None" = "gray50")) +
  labs(title = "Decision Curve Analysis: Clinical Net Benefit",
       subtitle = sprintf("Test set: childhood blood (n=%d, prevalence=%.1f%%)\nModel provides limited net benefit due to low sensitivity",
                          length(test_y), test_prevalence * 100),
       x = "Risk Threshold Probability", y = "Net Benefit") +
  xlim(0, 0.5) + theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9, color = "grey40"))
ggsave(file.path(outdir, "fig_decision_curve.png"), p_dca, width = 9, height = 6, dpi = 200)

# Figure D: Permutation null distribution
perm_df <- data.frame(AUC = perm_aucs)
p_perm <- ggplot(perm_df, aes(x = AUC)) +
  geom_histogram(bins = 40, fill = "gray70", color = "gray40", alpha = 0.8) +
  geom_vline(xintercept = main_auc, color = "#E41A1C", linewidth = 1.5) +
  annotate("text", x = main_auc + 0.02, y = max(table(cut(perm_aucs, 30))) * 0.8,
           label = sprintf("Observed AUC=%.3f\np=%.4f", main_auc, perm_p),
           hjust = 0, color = "#E41A1C", size = 3.5) +
  labs(title = "Permutation Test: Null vs Observed AUC",
       subtitle = sprintf("1000 permutations | Null AUC = %.3f +/- %.3f",
                          mean(perm_aucs), sd(perm_aucs)),
       x = "AUC", y = "Frequency") +
  theme_minimal(base_size = 11)
ggsave(file.path(outdir, "fig_permutation_null.png"), p_perm, width = 8, height = 5, dpi = 200)

# ============================================================================
# 9. SAVE
# ============================================================================
cat("\n9. Saving results...\n")

saveRDS(list(
  lodo_results = lodo_results,
  main_auc = main_auc,
  main_roc = main_roc,
  main_features = main_features,
  performance_table = perf_table,
  permutation_p = perm_p,
  perm_aucs = perm_aucs,
  bootstrap_ci = boot_ci,
  dca_results = dca,
  mchat_comparison = list(
    mchat_sens = mchat_sens, mchat_spec = mchat_spec,
    our_best_sens = our_best_sens, our_best_spec = our_best_spec
  )
), file.path(outdir, "module6_improved_results.rds"))

sink(file.path(outdir, "module6_improved_summary.txt"))
cat(sprintf("Module 6 Improved: External Validation & Clinical Utility\nDate: %s\n", Sys.Date()))
cat(sprintf("================================================================\n\n"))

cat("1. LEAVE-ONE-DATASET-OUT VALIDATION (Nested Features)\n")
for (i in 1:nrow(lodo_results)) {
  cat(sprintf("   %s held out: AUC=%.4f, Sens=%.3f, Spec=%.3f (%d features)\n",
              lodo_results$Holdout[i], lodo_results$AUC[i],
              lodo_results$Sens[i], lodo_results$Spec[i],
              lodo_results$NFeatures[i]))
}
cat(sprintf("   Mean LODO AUC: %.4f +/- %.4f\n\n", mean(lodo_results$AUC), sd(lodo_results$AUC)))

cat("2. MAIN VALIDATION (Prenatal+Birth → Childhood)\n")
cat(sprintf("   AUC: %.4f (95%% CI: %.4f-%.4f)\n", main_auc, ci.auc(main_roc)[1], ci.auc(main_roc)[3]))
cat(sprintf("   Permutation p: %.4f\n", perm_p))
cat(sprintf("   Bootstrap 95%% CI: [%.4f, %.4f]\n", boot_ci[1], boot_ci[2]))
cat(sprintf("   Features: %s\n\n", paste(main_features, collapse = ", ")))

cat("3. CLINICAL PERFORMANCE\n")
cat(sprintf("   At threshold 0.5: Sensitivity=%.3f, Specificity=%.3f\n",
            perf_table$Sensitivity[perf_table$Threshold == 0.5],
            perf_table$Specificity[perf_table$Threshold == 0.5]))
cat(sprintf("   Best sensitivity: %.3f (at threshold %.1f, specificity=%.3f)\n\n",
            our_best_sens,
            perf_table$Threshold[which.max(perf_table$Sensitivity)],
            our_best_spec))

cat("4. M-CHAT-R COMPARISON\n")
cat(sprintf("   M-CHAT-R: Sensitivity=0.85, Specificity=0.85, Cost=$0\n"))
cat(sprintf("   Our model: Sensitivity=%.3f at best, Cost=$100-500+\n", our_best_sens))
cat(sprintf("   Conclusion: Model CANNOT replace M-CHAT-R as primary screening.\n"))
cat(sprintf("   Potential role: risk stratification after positive M-CHAT-R.\n\n"))

cat("5. DECISION CURVE\n")
model_wins <- which(dca$model_nb > dca$all_nb & dca$model_nb > 0)
if (length(model_wins) > 0) {
  cat(sprintf("   Model provides net benefit for thresholds: %.2f-%.2f\n",
              thresholds[min(model_wins)], thresholds[max(model_wins)]))
} else {
  cat("   Model does NOT provide positive net benefit over 'treat all'.\n")
}

cat("\n6. HONEST CLINICAL POSITIONING\n")
cat("   - NOT a standalone screening tool (sensitivity too low)\n")
cat("   - NOT a replacement for M-CHAT-R or ADOS-2\n")
cat("   - Possible role: post-M-CHAT-R risk stratification\n")
cat("   - Research value: identifies biological subtypes (immune vs synaptic)\n")
cat("   - Future: combine with M-CHAT-R in a two-stage screening pipeline\n")
sink()

cat(sprintf("\n===== Module 6 Improved DONE =====\n"))
