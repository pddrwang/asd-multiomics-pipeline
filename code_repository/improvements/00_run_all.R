# ============================================================================
# Master Script: Run All Improved Modules
# Execute: Rscript improvements/00_run_all_improvements.R
# ============================================================================
cat("===================================================================\n")
cat("  ASD Multi-Omics Study — Comprehensive Improvement Suite\n")
cat("  Date:", as.character(Sys.Date()), "\n")
cat("===================================================================\n\n")

workdir <- "~/ASD_multiomics"
imp_dir <- file.path(workdir, "improvements")

cat("This script contains all improvements addressing the audit findings.\n")
cat("Run each section individually, or source this file to run all.\n\n")

cat("Modules to run:\n")
cat("  1. module1_improved.R   — ComBat diagnostics & PCA audit\n")
cat("  2. module2_improved.R   — Permutation test enhancement & power analysis\n")
cat("  3. module3_improved.R   — MR F-statistics & sensitivity analyses\n")
cat("  4. module5_improved.R   — Prediction model audit & clinical utility\n")
cat("  5. module6_improved.R   — External validation & M-CHAT-R benchmark\n")
cat("  6. module8_improved.R   — BH correction audit (TLN2/CCK)\n\n")

# Uncomment to run:
# source(file.path(imp_dir, "module1_improved.R"))
# source(file.path(imp_dir, "module2_improved.R"))
# source(file.path(imp_dir, "module3_improved.R"))
# source(file.path(imp_dir, "module5_improved.R"))
# source(file.path(imp_dir, "module6_improved.R"))
# source(file.path(imp_dir, "module8_improved.R"))

cat("Done. See improvements/*/ for each module's output.\n")
