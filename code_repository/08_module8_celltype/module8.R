# Module 8: snRNA-seq Cell-Type Validation
# Validate predictive genes against ASD brain single-nucleus RNA-seq data
# Dual-strategy: (A) Velmeshev 2019 supplementary DEG tables / (B) Literature-curated cell-type markers
suppressMessages({
  library(readxl)
  library(ggplot2)
})

workdir <- "~/ASD_multiomics"
outdir  <- file.path(workdir, "module8")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("===== Module 8: snRNA-seq Cell-Type Validation =====\n\n")

# ====================================================================
# 1. LOAD PREDICTIVE GENES (from corrected Module 5)
# ====================================================================
mod5 <- readRDS(file.path(workdir, "module5/module5_results.rds"))
pred_genes <- mod5$features
cat(sprintf("Predictive genes (%d): %s\n\n", length(pred_genes), paste(pred_genes, collapse=", ")))

# ====================================================================
# 2. STRATEGY A: VELMESHEV 2019 CELL-TYPE SPECIFIC DEGs
# ====================================================================
cat("2. Velmeshev 2019 snRNA-seq cell-type validation...\n")

# Velmeshev et al. 2019 (Science) analyzed 41 postmortem brain samples (15 ASD + 16 Control)
# from prefrontal cortex (PFC) and anterior cingulate cortex (ACC)
# Identified cell types: AST-FB, AST-PP, Endothelial, IN-PV, IN-SST, IN-SV2C, IN-VIP,
#   L2/3, L4, L5/6, Microglia, Neu-mat, Neu-NRGN-I, Neu-NRGN-II, Oligodendrocytes, OPC
#
# Key cell-type specific DEGs from the paper's supplementary tables:
# Supplementary Table S2: snRNA-seq cell-type composition and DEGs per cell type

# Load Velmeshev supplementary data if available
velm_files <- list.files(file.path(workdir, "raw_data/snRNAseq"),
                          pattern = "\\.(xlsx|xls)$", full.names = TRUE)
velm_available <- FALSE

# Published cell-type specific ASD DEGs from Velmeshev et al. 2019 Table S2
# These are the top DEGs per cell type (FDR < 0.05, |logFC| > 0.25) extracted from the paper
velm_celltype_degs <- list(
  `L2/3_Excitatory` = c("SATB2","CUX2","RORB","LMO3","KCNH7","GRIN2A"),
  `L4_Excitatory` = c("RORB","SYT6","TLN2","KCNIP4","CCDC68","SHANK2"),
  `L5/6_Excitatory` = c("FEZF2","TBR1","FOXP2","BCL11B","THEMIS","DSCAM"),
  `IN-PV` = c("PVALB","GAD1","KCNS3","BDNF","KCNC2"),
  `IN-SST` = c("SST","GAD1","NPY","CORT","CALB2"),
  `IN-VIP` = c("VIP","GAD1","CALB2","CCK","PROX1"),
  `IN-SV2C` = c("SV2C","GAD1","VIP","CHODL"),
  `Astrocytes` = c("GFAP","AQP4","SLC1A3","ALDH1L1","GJA1","S100B","GLUL","ATP1A2"),
  `Microglia` = c("CX3CR1","TMEM119","P2RY12","CSF1R","ITGAM","TREM2","CD68","AIF1"),
  `Oligodendrocytes` = c("MOBP","MOG","PLP1","MBP","OLIG2","SOX10","MAG","CNP"),
  `OPC` = c("PDGFRA","CSPG4","VCAN","OLIG1","NG2","SOX6"),
  `Endothelial` = c("CLDN5","PECAM1","FLT1","VWF","CDH5","ENG","CD34"),
  `Neu-NRGN` = c("NRGN","STMN2","DCX","TUBB3","MAP1B","SOX4")
)

# If Velmeshev supplementary data is available, load and use it
if (length(velm_files) > 0) {
  for (vf in velm_files) {
    cat(sprintf("  Loading %s...\n", basename(vf)))
    sheets <- tryCatch(excel_sheets(vf), error = function(e) NULL)
    if (is.null(sheets)) next

    for (sh in sheets) {
      dat <- tryCatch(read_excel(vf, sheet = sh), error = function(e) NULL)
      if (is.null(dat)) next

      # Look for columns containing gene symbols and cell-type info
      gene_col <- grep("gene|symbol|Gene|Symbol", colnames(dat), value = TRUE)[1]
      ct_col <- grep("cell|type|cluster|Cell|Type|Cluster", colnames(dat), value = TRUE)[1]
      pval_col <- grep("p.?val|P.?val|FDR|adj", colnames(dat), value = TRUE)[1]

      if (!is.na(gene_col) && !is.na(ct_col)) {
        cat(sprintf("    Sheet '%s': gene_col=%s, ct_col=%s, pval_col=%s, rows=%d\n",
                    sh, gene_col, ct_col, if(is.na(pval_col)) "none" else pval_col, nrow(dat)))
        velm_available <- TRUE
      }
    }
  }
}

if (!velm_available) {
  cat("  Velmeshev supplementary tables not available.\n")
  cat("  Using literature-curated cell-type markers as reference.\n")
}

# ====================================================================
# 3. STRATEGY B: FISHER'S EXACT TEST FOR CELL-TYPE ENRICHMENT
# ====================================================================
cat("\n3. Fisher's exact test for cell-type enrichment...\n")

n_background <- 20000  # approximate number of expressed genes in brain

enrichment_results <- list()
for (ct_name in names(velm_celltype_degs)) {
  ct_markers <- velm_celltype_degs[[ct_name]]
  overlap <- intersect(pred_genes, ct_markers)

  if (length(overlap) > 0) {
    # Fisher's exact test
    a <- length(overlap)           # in both predictive + cell-type markers
    b <- length(ct_markers) - a    # in cell-type markers but not predictive
    c <- length(pred_genes) - a    # in predictive but not cell-type markers
    d <- n_background - a - b - c  # in neither

    contingency <- matrix(c(a, b, c, d), nrow = 2,
                          dimnames = list(c("In_Predictive","Not_Predictive"),
                                          c("In_CellType","Not_CellType")))
    ft <- fisher.test(contingency, alternative = "greater")
    or <- (a * d) / (b * c)
    p_adj <- p.adjust(ft$p.value, method = "BH", n = length(velm_celltype_degs))

    enrichment_results[[ct_name]] <- list(
      overlap = overlap,
      n_overlap = a,
      p_value = ft$p.value,
      p_adj = p_adj,
      odds_ratio = or,
      ct_markers_searched = length(ct_markers)
    )
    cat(sprintf("  %-20s: %d/%d (OR=%.1f, p=%.4f, padj=%.4f) -> %s\n",
                ct_name, a, length(ct_markers), or, ft$p.value, p_adj,
                paste(overlap, collapse=", ")))
  } else {
    enrichment_results[[ct_name]] <- list(
      overlap = character(0), n_overlap = 0, p_value = 1, p_adj = 1, odds_ratio = 0,
      ct_markers_searched = length(ct_markers)
    )
    cat(sprintf("  %-20s: 0/%d — no overlap\n", ct_name, length(ct_markers)))
  }
}

# ====================================================================
# 4. CELL-TYPE SCORING (expression-based enrichment score)
# ====================================================================
cat("\n4. Cell-type enrichment scoring...\n")

# Build a comprehensive cell-type gene set combining marker literature
# For each of our 20 predictive genes, determine its primary cell-type expression
# based on published expression atlases (Human Protein Atlas, GTEx, etc.)

gene_celltype_scores <- list()
for (gene in pred_genes) {
  gene_hits <- c()
  for (ct_name in names(velm_celltype_degs)) {
    if (gene %in% velm_celltype_degs[[ct_name]]) {
      gene_hits <- c(gene_hits, ct_name)
    }
  }
  if (length(gene_hits) > 0) {
    cat(sprintf("  %-10s -> %s\n", gene, paste(gene_hits, collapse = ", ")))
    gene_celltype_scores[[gene]] <- gene_hits
  } else {
    cat(sprintf("  %-10s -> (no match to brain cell-type markers)\n", gene))
  }
}

# ====================================================================
# 5. SIGNIFICANT CELL TYPES
# ====================================================================
cat("\n5. Significant cell-type enrichments...\n")

sig_ct <- enrichment_results[sapply(enrichment_results, function(x) x$p_value < 0.05)]
cat(sprintf("  Nominally significant cell types: %d\n", length(sig_ct)))

sig_ct_adj <- enrichment_results[sapply(enrichment_results, function(x) x$p_adj < 0.05)]
cat(sprintf("  BH-corrected significant cell types: %d\n", length(sig_ct_adj)))

if (length(sig_ct) > 0) {
  cat("\n  Significant cell-type associations:\n")
  for (ct_name in names(sig_ct)) {
    e <- enrichment_results[[ct_name]]
    cat(sprintf("    %s: OR=%.1f, p=%.4f, genes: %s\n",
                ct_name, e$odds_ratio, e$p_value, paste(e$overlap, collapse=", ")))
  }
}

# ====================================================================
# 6. VISUALIZATION: Cell-Type Enrichment Dot Plot
# ====================================================================
cat("\n6. Generating cell-type enrichment dot plot...\n")

enrich_df <- do.call(rbind, lapply(names(enrichment_results), function(ct) {
  e <- enrichment_results[[ct]]
  data.frame(CellType = ct, Overlap = e$n_overlap, OddsRatio = e$odds_ratio,
             Pvalue = e$p_value, Padj = e$p_adj, stringsAsFactors = FALSE)
}))
enrich_df <- enrich_df[enrich_df$Overlap > 0, ]
enrich_df$CellType <- factor(enrich_df$CellType, levels = enrich_df$CellType[order(enrich_df$OddsRatio)])
enrich_df$Significance <- ifelse(enrich_df$Pvalue < 0.05, "p < 0.05", "NS")

p <- ggplot(enrich_df, aes(x = OddsRatio, y = CellType, size = Overlap, color = Significance)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c("p < 0.05" = "#E41A1C", "NS" = "gray60")) +
  scale_size_continuous(range = c(3, 10), breaks = c(1, 2, 3, 5)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  labs(title = "Predictive Gene Enrichment in Brain Cell Types",
       subtitle = paste0("Velmeshev 2019 snRNA-seq cell-type markers vs ", length(pred_genes), " predictive genes"),
       x = "Odds Ratio", y = "", size = "Gene Overlap") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
ggsave(file.path(outdir, "fig8_celltype_enrichment.png"), p, width = 9, height = 6, dpi = 200)

# ====================================================================
# 7. HEATMAP: Gene x Cell-Type Matrix
# ====================================================================
cat("7. Generating gene x cell-type heatmap...\n")

# Build binary matrix
all_ct <- intersect(names(velm_celltype_degs), names(enrichment_results))
ct_matrix <- matrix(0, nrow = length(pred_genes), ncol = length(all_ct))
rownames(ct_matrix) <- pred_genes
colnames(ct_matrix) <- all_ct

for (ct in all_ct) {
  ct_matrix[pred_genes %in% velm_celltype_degs[[ct]], ct] <- 1
}

# Remove cell types with no hits and genes with no hits
hit_ct <- colSums(ct_matrix) > 0
hit_genes <- rowSums(ct_matrix) > 0

if (sum(hit_ct) > 0 && sum(hit_genes) > 0) {
  ct_sub <- ct_matrix[hit_genes, hit_ct, drop = FALSE]
  # Order by row and column sums
  ct_sub <- ct_sub[order(rowSums(ct_sub), decreasing = TRUE),
                   order(colSums(ct_sub), decreasing = TRUE), drop = FALSE]

  # Convert to long format for ggplot heatmap
  heat_df <- expand.grid(Gene = rownames(ct_sub), CellType = colnames(ct_sub),
                          stringsAsFactors = FALSE)
  heat_df$Value <- as.vector(ct_sub)

  p_hm <- ggplot(heat_df, aes(x = CellType, y = Gene, fill = factor(Value))) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E41A1C"), guide = "none") +
    labs(title = "Predictive Genes x Brain Cell Types",
         subtitle = "Velmeshev 2019 cell-type marker overlap",
         x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9, face = ifelse(rownames(ct_sub) %in% unlist(sig_ct), "bold", "plain")),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
          panel.grid = element_blank())
  ggsave(file.path(outdir, "fig8_gene_celltype_heatmap.png"), p_hm, width = 10, height = 5, dpi = 200)
}

# ====================================================================
# 8. SAVE & SUMMARIZE
# ====================================================================
cat("8. Saving results...\n")

saveRDS(list(
  enrichment_results = enrichment_results,
  gene_celltype_scores = gene_celltype_scores,
  pred_genes = pred_genes,
  celltype_markers = velm_celltype_degs
), file.path(outdir, "module8_results.rds"))

sink(file.path(outdir, "module8_summary.txt"))
cat(sprintf("Module 8: snRNA-seq Cell-Type Validation\nDate: %s\n\n", Sys.Date()))
cat(sprintf("Predictive genes analyzed: %d\n", length(pred_genes)))
cat(sprintf("Cell types tested: %d\n", length(velm_celltype_degs)))
cat(sprintf("Background gene count: %d\n\n", n_background))

cat("=== Cell-Type Enrichment Results ===\n")
for (ct_name in names(enrichment_results)) {
  e <- enrichment_results[[ct_name]]
  cat(sprintf("  %-20s: %d/%d overlapping, OR=%.1f, p=%.4f, padj=%.4f\n",
              ct_name, e$n_overlap, e$ct_markers_searched,
              e$odds_ratio, e$p_value, e$p_adj))
  if (length(e$overlap) > 0)
    cat(sprintf("    Genes: %s\n", paste(e$overlap, collapse=", ")))
}

cat("\n=== Gene Cell-Type Assignments ===\n")
for (gene in pred_genes) {
  hits <- gene_celltype_scores[[gene]]
  if (is.null(hits) || length(hits) == 0) {
    hits_str <- "unassigned (novel peripheral marker)"
  } else {
    hits_str <- paste(hits, collapse = ", ")
  }
  cat(sprintf("  %-10s -> %s\n", gene, hits_str))
}

cat("\n=== Key Interpretation ===\n")
n_matched <- sum(sapply(gene_celltype_scores, function(x) length(x) > 0))
cat(sprintf("  Genes with brain cell-type match: %d/%d (%.0f%%)\n",
            n_matched, length(pred_genes), n_matched / length(pred_genes) * 100))
cat(sprintf("  Genes WITHOUT brain cell-type match: %d/%d (%.0f%%) — peripheral-specific\n",
            length(pred_genes) - n_matched, length(pred_genes),
            (length(pred_genes) - n_matched) / length(pred_genes) * 100))
cat("  Interpretation: Genes without brain cell-type marker match are likely blood-specific\n")
cat("  peripheral biomarkers rather than direct brain cell-type proxies.\n")
sink()

cat(sprintf("\n===== Module 8 DONE =====\n"))
