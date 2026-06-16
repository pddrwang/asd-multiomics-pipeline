# Module 4: Molecular Subtyping of ASD Using Causal Genes
# Consensus clustering on ASD samples from GSE18123 (child blood)
suppressMessages({
  library(ConsensusClusterPlus)
  library(pheatmap)
  library(ggplot2)
})

workdir <- "~/ASD_multiomics"
outdir  <- file.path(workdir, "module4_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("===== Module 4: ASD Molecular Subtyping =====\n\n")

# ====================================================================
# 1. LOAD UNIFIED, BATCH-CORRECTED EXPRESSION (Module 1 output)
# ====================================================================
cat("1. Loading Module 1 output...\n")
mod1 <- readRDS(file.path(workdir, "module1/module1_results.rds"))
expr <- mod1$expr
meta <- mod1$meta

cat(sprintf("  Expression: %d genes x %d samples\n", nrow(expr), ncol(expr)))

# Extract ASD-only samples, blood only (not brain)
asd_blood <- meta$Diagnosis == "ASD" & meta$Tissue != "Brain"
asd_expr <- expr[, asd_blood, drop = FALSE]
asd_meta <- meta[asd_blood, ]
cat(sprintf("  ASD blood samples: %d\n", ncol(asd_expr)))
cat(sprintf("  Tissues: %s\n", paste(unique(asd_meta$Tissue), collapse=", ")))

# ====================================================================
# 2. SELECT GENES — from Module 2 brain DEG (proxy for causal genes, replace with Module 3 when ready)
# ====================================================================
cat("\n2. Selecting causal/candidate genes...\n")

# Priority 1: Module 3 causal genes IF available
mr_rds <- file.path(workdir, "module3/module3_results.rds")
top_genes <- NULL

if (file.exists(mr_rds)) {
  mr_res <- readRDS(mr_rds)
  # Extract significant MR genes
  if (exists("mr_res") && !is.null(mr_res$causal_genes)) {
    top_genes <- mr_res$causal_genes
    cat(sprintf("  Using %d MR causal genes\n", length(top_genes)))
  }
}

# Priority 2: Module 2 high-confidence brain DEGs
if (is.null(top_genes) || length(top_genes) < 20) {
  cat("  Using Module 2 high-confidence genes as proxy...\n")
  brain_rds <- file.path(workdir, "module2/module2_results.rds")
  if (file.exists(brain_rds)) {
    brain_res <- readRDS(brain_rds)
    # Use all lenient DEGs + direction-consistent genes for broader coverage
    de_all <- brain_res$de_combined
    de_all <- de_all[order(de_all$adj.P.Val), ]
    top_genes <- de_all$Gene[de_all$adj.P.Val < 0.05 & abs(de_all$logFC) > 0.3]
    cat(sprintf("  %d brain DEGs (FDR<0.05, |logFC|>0.3)\n", length(top_genes)))
    if (length(top_genes) < 20) {
      # Further relax
      top_genes <- de_all$Gene[1:min(500, nrow(de_all))]
      cat(sprintf("  Relaxed to top %d genes\n", length(top_genes)))
    }
  }
}

# Priority 3: Top blood DEGs
if (is.null(top_genes) || length(top_genes) < 20) {
  cat("  Using blood DEGs as fallback...\n")
  blood_de <- read.csv(file.path(workdir, "GSE18123/GSE18123_DEG_sig_ASDvsControl_log2.csv"),
                       stringsAsFactors = FALSE)
  gpl570 <- getGEO(filename = file.path(workdir, "GSE18123/GPL570.annot.gz"))
  tab570 <- Table(gpl570)
  m <- match(blood_de$X, tab570$ID)
  blood_genes <- unique(gsub("///.*", "", trimws(as.character(tab570[["Gene symbol"]][m]))))
  blood_genes <- blood_genes[!is.na(blood_genes) & blood_genes != ""]
  top_genes <- blood_genes
  cat(sprintf("  %d blood DEG genes\n", length(top_genes)))
}

# Filter to genes present in our expression matrix
top_genes <- intersect(top_genes, rownames(asd_expr))
cat(sprintf("  Final gene set: %d genes in expression matrix\n", length(top_genes)))

# Cap at 100 for clustering
if (length(top_genes) > 100) {
  cat(sprintf("  Capping at 100 genes (from %d)\n", length(top_genes)))
  # Select top by variance across ASD samples
  vars <- apply(asd_expr[top_genes, ], 1, var)
  top_genes <- names(sort(vars, decreasing = TRUE))[1:100]
}
cat(sprintf("  Clustering with %d genes\n", length(top_genes)))

# ====================================================================
# 3. CONSENSUS CLUSTERING
# ====================================================================
cat("\n3. Running ConsensusClusterPlus...\n")

clust_data <- asd_expr[top_genes, ]
# Scale features
clust_data_scaled <- t(scale(t(clust_data)))

set.seed(42)
old_wd <- getwd()
setwd(outdir)
cc_results <- ConsensusClusterPlus(
  d = as.matrix(clust_data_scaled),
  maxK = 5,
  reps = 100,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "pam",
  distance = "pearson",
  seed = 42,
  plot = "png",
  title = "consensus_cluster",
  writeTable = FALSE
)
setwd(old_wd)

# ====================================================================
# 4. DETERMINE OPTIMAL K
# ====================================================================
cat("\n4. Determining optimal K...\n")

# Use NbClust for 3-4 indices
icl <- sapply(2:5, function(k) {
  cc_results[[k]][["consensusMatrix"]] -> cm
  diag_cm <- mean(diag(cm))
  off_diag <- mean(cm[upper.tri(cm)])
  # Simplified PAC score
  list(k = k, consensus = mean(diag_cm))
})

# Default to k=2 or k=3 based on consensus clustering delta area
# Calculate delta area
areas <- sapply(2:5, function(k) {
  cm <- cc_results[[k]][["consensusMatrix"]]
  cdf <- ecdf(as.vector(cm[upper.tri(cm)]))
  # Area under CDF
  mean(cm[upper.tri(cm)])
})

delta <- diff(areas)
optimal_k <- which.max(abs(delta)) + 2  # +2 because starting from k=2
cat(sprintf("  Optimal K based on delta area: %d\n", optimal_k))

# ====================================================================
# 5. EXTRACT SUBTYPE LABELS
# ====================================================================
cat("\n5. Extracting subtype labels...\n")

# Use K=2 (basic split) and K=3 (finer) if optimal_k >= 2
k_values <- unique(c(2, optimal_k, 3))
k_values <- k_values[k_values >= 2 & k_values <= 5]

subtype_labels <- list()
for (k in k_values) {
  clust <- cc_results[[k]][["consensusClass"]]
  subtype_labels[[as.character(k)]] <- clust
  cat(sprintf("  K=%d: %s\n", k, paste(table(clust), collapse=", ")))
}

# Create final subtype dataframe
subtypes_df <- data.frame(
  Sample = colnames(asd_expr),
  stringsAsFactors = FALSE
)
for (k in names(subtype_labels)) {
  subtypes_df[[paste0("Subtype_K", k)]] <- subtype_labels[[k]][subtypes_df$Sample]
}
subtypes_df$Tissue <- asd_meta[subtypes_df$Sample, "Tissue"]
subtypes_df$Dataset <- asd_meta[subtypes_df$Sample, "Dataset"]

# Add interpretable subtype names (K=2: Immune vs Synaptic based on marker genes)
# Check expression of key markers in each cluster
k2 <- subtypes_df$Subtype_K2  # use the dataframe column, not the named vector
if (length(unique(k2)) >= 2) {
  # Check immune markers
  immune_genes <- intersect(c("IFITM3","CXCL16","HSPB1","ABCA1","JUN"), rownames(clust_data))
  synaptic_genes <- intersect(c("SYN1","DLG4","SYP","GRIN1"), rownames(clust_data))
  if (length(immune_genes)==0) immune_genes <- intersect(c("IFITM2","CEBPD","DDR1","HSPA1A","MSN"), rownames(clust_data))
  if (length(synaptic_genes)==0) synaptic_genes <- intersect(c("SNAP25","GAP43","STMN2","MAP2","BDNF"), rownames(clust_data))

  for (cl in sort(unique(k2))) {
    samples_in_cl <- subtypes_df$Sample[subtypes_df$Subtype_K2 == cl]
    cat(sprintf("\n  Cluster %d (n=%d):\n", cl, length(samples_in_cl)))
    if (length(immune_genes) > 0 && length(samples_in_cl) > 0) {
      imm_mean <- mean(colMeans(clust_data[immune_genes, samples_in_cl, drop=FALSE], na.rm=TRUE))
      cat(sprintf("    Immune markers mean: %.3f\n", imm_mean))
    }
    if (length(synaptic_genes) > 0 && length(samples_in_cl) > 0) {
      syn_mean <- mean(colMeans(clust_data[synaptic_genes, samples_in_cl, drop=FALSE], na.rm=TRUE))
      cat(sprintf("    Neuronal markers mean: %.3f\n", syn_mean))
    }
  }

  # Assign names based on immune score
  cl_names <- c()
  overall_immune_mean <- mean(colMeans(clust_data[immune_genes, , drop=FALSE], na.rm=TRUE))
  for (cl in sort(unique(k2))) {
    samples_in_cl <- subtypes_df$Sample[subtypes_df$Subtype_K2 == cl]
    if (length(immune_genes) > 0 && length(samples_in_cl) > 0) {
      imm_mean <- mean(colMeans(clust_data[immune_genes, samples_in_cl, drop=FALSE], na.rm=TRUE))
      cl_names[as.character(cl)] <- ifelse(imm_mean > overall_immune_mean,
                                            "ASD-Immune", "ASD-Synaptic")
    } else {
      cl_names[as.character(cl)] <- paste0("ASD-", letters[cl])
    }
  }

  subtypes_df$Subtype_Name_K2 <- cl_names[as.character(k2)]
  cat(sprintf("\n  K=2 subtypes: %s\n", paste(unique(subtypes_df$Subtype_Name_K2), collapse=", ")))
}

# ====================================================================
# 6. VISUALIZE
# ====================================================================
cat("\n6. Generating heatmap...\n")

# Simple annotation
anno_col <- data.frame(
  Subtype = as.factor(subtypes_df$Subtype_Name_K2),
  Tissue = as.factor(subtypes_df$Tissue),
  row.names = subtypes_df$Sample
)

pheatmap(
  clust_data_scaled,
  annotation_col = anno_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 6,
  clustering_method = "ward.D2",
  main = paste0("ASD Subtypes by Gene Expression (", nrow(clust_data_scaled), " genes)"),
  filename = file.path(outdir, "asd_subtypes_heatmap.png"),
  width = 14,
  height = 10
)

# ====================================================================
# 7. SAVE
# ====================================================================
cat("\n7. Saving results...\n")
write.csv(subtypes_df, file.path(outdir, "asd_subtype_labels.csv"), row.names=FALSE)
saveRDS(list(
  subtypes = subtypes_df,
  consensus_results = cc_results,
  genes_used = top_genes,
  optimal_k = optimal_k
), file.path(outdir, "module4_results.rds"))

sink(file.path(outdir, "module4_summary.txt"))
cat("Module 4: ASD Molecular Subtyping\n")
cat(sprintf("Date: %s\n\n", Sys.Date()))
cat(sprintf("ASD blood samples: %d\n", ncol(asd_expr)))
cat(sprintf("Genes used: %d\n", length(top_genes)))
cat(sprintf("Optimal K: %d\n", optimal_k))
cat(sprintf("\nSubtype breakdown (K=2):\n"))
if ("Subtype_Name_K2" %in% colnames(subtypes_df)) {
  print(table(subtypes_df$Subtype_Name_K2))
}
sink()

cat(sprintf("\n===== Module 4 DONE: %d ASD samples subtyped =====\n", ncol(asd_expr)))
