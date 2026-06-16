# Module 1: Multi-tissue Integration (incremental, robust version)
# Load pre-computed expression CSVs, map probes -> genes, intersect, batch-correct
suppressMessages({
  library(limma)
  library(sva)
  library(org.Hs.eg.db)
  library(GEOquery)
})

workdir <- "~/ASD_multiomics"
outdir  <- file.path(workdir, "module1_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("===== Step 1: Load expression matrices =====\n")

# ---- GSE18123 (Child Blood, GPL570, log2) ----
cat("\n1. GSE18123...\n")
e1 <- read.csv(file.path(workdir, "raw_data/GEO/GSE18123/GSE18123_expression_matrix_log2.csv"),
               row.names = 1, check.names = FALSE)
e1 <- as.matrix(e1)
# Map probes to gene symbols via GPL570 annotation
gpl570 <- getGEO(filename = file.path(workdir, "raw_data/GEO/GSE18123/GPL570.annot.gz"))
tab570 <- Table(gpl570)
m1 <- match(rownames(e1), tab570$ID)
g1 <- as.character(tab570[["Gene symbol"]])[m1]
g1 <- gsub("///.*", "", g1); g1 <- trimws(g1)
keep1 <- !is.na(g1) & g1 != ""
e1 <- e1[keep1, ]; g1 <- g1[keep1]
# Collapse to gene-level (take max-var probe per gene)
gene_list1 <- split(seq_len(nrow(e1)), g1)
e1_gene <- t(sapply(gene_list1, function(idx) {
  if (length(idx) == 1) e1[idx, ] else { v <- apply(e1[idx, , drop=FALSE], 1, var, na.rm=TRUE); e1[idx[which.max(v)], ] }
}))
cat(sprintf("  %d probes -> %d genes, %d samples\n", length(g1), nrow(e1_gene), ncol(e1_gene)))

# Get metadata
pheno1 <- read.csv(file.path(workdir, "raw_data/GEO/GSE18123/GSE18123_group_info.csv"), stringsAsFactors=FALSE)
rownames(pheno1) <- pheno1$Sample
pheno1 <- pheno1[colnames(e1), ]
# Determine ASD/Control
diag_clean1 <- trimws(gsub("^diagnosis:\\s*", "", pheno1$Diagnosis, ignore.case=FALSE))
grp1 <- ifelse(diag_clean1 == "CONTROL", "Control",
        ifelse(diag_clean1 %in% c("AUTISM", "ASPERGER'S DISORDER", "PDD-NOS"), "ASD", "Other"))
keep_samp1 <- grp1 %in% c("ASD", "Control")
e1_gene <- e1_gene[, keep_samp1]; grp1 <- grp1[keep_samp1]
cat(sprintf("  After filtering: %d ASD + %d Control\n", sum(grp1=="ASD"), sum(grp1=="Control")))

# ---- GSE123302 (Cord Blood, GPL16686) ----
cat("\n2. GSE123302...\n")
e2 <- read.csv(file.path(workdir, "raw_data/GEO/GSE123302/GSE123302_expression_matrix_ASDvsTD.csv"),
               row.names = 1, check.names = FALSE)
e2 <- as.matrix(e2)
# Map probes -> gene symbols via GPL16686 annotation
gpl16686 <- getGEO(filename = file.path(workdir, "raw_data/GEO/GSE123302/GPL16686.soft.gz"))
tab16686 <- Table(gpl16686)
m2 <- match(rownames(e2), tab16686$ID)
acc2 <- as.character(tab16686$GB_ACC)[m2]
# Map RefSeq -> Gene Symbol via select()
cat("  Mapping RefSeq to Gene Symbols...\n")
valid2 <- acc2 != "" & !is.na(acc2)
sym2 <- rep(NA, length(acc2))
if (any(valid2)) {
  mapped <- suppressMessages(select(org.Hs.eg.db, keys = acc2[valid2],
                                    columns = "SYMBOL", keytype = "REFSEQ"))
  mapped <- mapped[!duplicated(mapped$REFSEQ), ]
  mm <- match(acc2, mapped$REFSEQ)
  sym2 <- mapped$SYMBOL[mm]
}
keep2 <- !is.na(sym2) & sym2 != ""
e2 <- e2[keep2, ]; sym2 <- sym2[keep2]
gene_list2 <- split(seq_len(nrow(e2)), sym2)
e2_gene <- t(sapply(gene_list2, function(idx) {
  if (length(idx) == 1) e2[idx, ] else { v <- apply(e2[idx, , drop=FALSE], 1, var, na.rm=TRUE); e2[idx[which.max(v)], ] }
}))
cat(sprintf("  %d probes -> %d genes, %d samples\n", length(sym2), nrow(e2_gene), ncol(e2_gene)))

pheno2 <- read.csv(file.path(workdir, "raw_data/GEO/GSE123302/GSE123302_group_info.csv"), stringsAsFactors=FALSE)
grp2 <- pheno2$Group
grp2[grp2 == "TD"] <- "Control"
cat(sprintf("  %d ASD + %d Control\n", sum(grp2=="ASD"), sum(grp2=="Control")))

# ---- GSE148450 (Maternal Blood, GPL25483, ALREADY gene symbols) ----
cat("\n3. GSE148450...\n")
e3 <- read.csv(file.path(workdir, "raw_data/GEO/GSE148450/GSE148450_expression_matrix_ASDvsTD.csv"),
               row.names = 1, check.names = FALSE)
e3 <- as.matrix(e3)
g3 <- rownames(e3)
keep3 <- !is.na(g3) & g3 != "" & !grepl("^\\d", g3)
e3 <- e3[keep3, ]; g3 <- g3[keep3]
gene_list3 <- split(seq_len(nrow(e3)), g3)
e3_gene <- t(sapply(gene_list3, function(idx) {
  if (length(idx) == 1) e3[idx, ] else { v <- apply(e3[idx, , drop=FALSE], 1, var, na.rm=TRUE); e3[idx[which.max(v)], ] }
}))
cat(sprintf("  %d probes -> %d genes, %d samples\n", length(g3), nrow(e3_gene), ncol(e3_gene)))

pheno3 <- read.csv(file.path(workdir, "raw_data/GEO/GSE148450/GSE148450_group_info.csv"), stringsAsFactors=FALSE)
grp3 <- pheno3$Group
grp3[grp3 == "TD"] <- "Control"
cat(sprintf("  %d ASD + %d Control\n", sum(grp3=="ASD"), sum(grp3=="Control")))

# ---- GSE38322 (Brain, GPL10558) ----
cat("\n4. GSE38322...\n")
e4 <- read.csv(file.path(workdir, "raw_data/GEO/GSE38322/GSE38322_expression_matrix.csv"),
               row.names = 1, check.names = FALSE)
e4 <- as.matrix(e4)
gpl10558 <- getGEO(filename = file.path(workdir, "raw_data/GEO/GSE38322/GPL10558.annot.gz"))
tab10558 <- Table(gpl10558)
m4 <- match(rownames(e4), tab10558$ID)
g4 <- as.character(tab10558[["Gene symbol"]])[m4]
g4 <- trimws(g4)
keep4 <- !is.na(g4) & g4 != ""
e4 <- e4[keep4, ]; g4 <- g4[keep4]
gene_list4 <- split(seq_len(nrow(e4)), g4)
e4_gene <- t(sapply(gene_list4, function(idx) {
  if (length(idx) == 1) e4[idx, ] else { v <- apply(e4[idx, , drop=FALSE], 1, var, na.rm=TRUE); e4[idx[which.max(v)], ] }
}))
cat(sprintf("  %d probes -> %d genes, %d samples\n", length(g4), nrow(e4_gene), ncol(e4_gene)))

# Parse diagnosis from saved file (pre-extracted from source_name_ch1)
diag4_rds <- readRDS(file.path(workdir, "raw_data/GEO/GSE38322/GSE38322_diagnosis.rds"))
diag4 <- diag4_rds[colnames(e4_gene)]
cat(sprintf("  %d ASD + %d Control\n", sum(diag4=="ASD"), sum(diag4=="Control")))

# ====================================================
# Step 2: Intersect genes
# ====================================================
cat("\n===== Step 2: Gene intersection =====\n")
common <- Reduce(intersect, list(rownames(e1_gene), rownames(e2_gene),
                                  rownames(e3_gene), rownames(e4_gene)))
cat(sprintf("Common genes: %d\n", length(common)))

m1f <- e1_gene[common, ]
m2f <- e2_gene[common, ]
m3f <- e3_gene[common, ]
m4f <- e4_gene[common, ]
colnames(m1f) <- paste0("GSE18123_", colnames(m1f))
colnames(m2f) <- paste0("GSE123302_", colnames(m2f))
colnames(m3f) <- paste0("GSE148450_", colnames(m3f))
colnames(m4f) <- paste0("GSE38322_", colnames(m4f))

unified <- cbind(m1f, m2f, m3f, m4f)
cat(sprintf("Unified matrix: %d genes x %d samples\n", nrow(unified), ncol(unified)))

# ====================================================
# Step 3: Low-expression filter
# ====================================================
cat("\n===== Step 3: Low-expression filter =====\n")
gene_rate <- rowMeans(unified > 2)
keep_genes <- gene_rate >= 0.10
unified <- unified[keep_genes, ]
cat(sprintf("After filter: %d genes (removed %d low-expr)\n",
            nrow(unified), sum(!keep_genes)))

# ====================================================
# Step 4: Quantile normalize within each dataset
# ====================================================
cat("\n===== Step 4: Quantile normalization =====\n")
batch <- rep(c("GSE18123","GSE123302","GSE148450","GSE38322"),
             c(ncol(m1f), ncol(m2f), ncol(m3f), ncol(m4f)))
for (b in unique(batch)) {
  idx <- which(batch == b)
  unified[, idx] <- normalizeBetweenArrays(unified[, idx], method = "quantile")
  cat(sprintf("  QN %s: %d samples\n", b, length(idx)))
}

# ====================================================
# Step 5: Build metadata
# ====================================================
cat("\n===== Step 5: Metadata =====\n")
diagnosis <- c(grp1, grp2, grp3, diag4)
names(diagnosis) <- colnames(unified)

# Tissue (confounded with batch)
tissue <- c(rep("Peripheral_Blood", ncol(m1f)),
            rep("Cord_Blood", ncol(m2f)),
            rep("Maternal_Blood", ncol(m3f)),
            rep("Brain", ncol(m4f)))

# Stage
stage <- c(rep("Childhood", ncol(m1f)),
           rep("Birth", ncol(m2f)),
           rep("Prenatal", ncol(m3f)),
           rep("Postmortem", ncol(m4f)))

meta <- data.frame(
  Sample    = colnames(unified),
  Dataset   = batch,
  Tissue    = tissue,
  Stage     = stage,
  Diagnosis = diagnosis,
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$Sample

cat("Sample breakdown:\n")
print(table(meta$Tissue, meta$Diagnosis))

# ====================================================
# Step 6: ComBat batch correction
# ====================================================
cat("\n===== Step 6: ComBat batch correction =====\n")
mod <- model.matrix(~ Diagnosis, data = meta)
unified_corrected <- ComBat(dat = unified, batch = batch, mod = mod,
                             par.prior = TRUE, prior.plots = FALSE)
cat("ComBat completed.\n")

# ====================================================
# Step 7: Save
# ====================================================
cat("\n===== Step 7: Save =====\n")
write.csv(unified_corrected, file.path(outdir, "unified_expression_matrix_batch_corrected.csv"))
write.csv(meta, file.path(outdir, "unified_metadata.csv"), row.names = FALSE)
saveRDS(list(expr = unified_corrected, meta = meta, common_genes = common),
        file.path(outdir, "module1_results.rds"))

sink(file.path(outdir, "module1_summary.txt"))
cat(sprintf("Module 1: Multi-tissue Integration\n%s\n", Sys.Date()))
cat(sprintf("Unified matrix: %d genes x %d samples\n", nrow(unified_corrected), ncol(unified_corrected)))
cat(sprintf("Common genes (before low-expr filter): %d\n", length(common)))
cat("Batch correction: ComBat (batch=Dataset, protect=Diagnosis)\n")
cat("\nSamples per tissue x diagnosis:\n")
print(table(meta$Tissue, meta$Diagnosis))
sink()

cat("\n===== Module 1 DONE =====\n")
cat(sprintf("Output: %d genes x %d samples in %s\n",
            nrow(unified_corrected), ncol(unified_corrected), outdir))
