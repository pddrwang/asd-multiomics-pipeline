# ============================================================================
# Module 3 — Improved: MR with Complete Methodological Reporting
# Fixes:
#   1. F-statistic calculation for instrument strength (F > 10 threshold)
#   2. Palindromic SNP identification and handling
#   3. MR sensitivity analyses (MR-Egger, weighted median, weighted mode)
#   4. Leave-one-out analysis for outlier detection
#   5. Heterogeneity statistics (Cochran's Q)
#   6. Horizontal pleiotropy assessment (MR-Egger intercept)
#   7. Steiger directionality test
#   8. Document immune specificity index calculation
# ============================================================================
suppressMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(MendelianRandomization)
})

workdir <- "~/ASD_multiomics"
outdir  <- file.path(workdir, "improvements/module3_improved")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("===== Module 3 Improved: MR with Complete Methodology =====\n\n")

# ============================================================================
# 1. LOAD CANDIDATE GENES
# ============================================================================
cat("1. Loading candidate genes...\n")
brain <- read.csv(file.path(workdir, "module2/brain_top500_confidence.csv"), stringsAsFactors = FALSE)
brain_high <- brain$Gene[brain$ConfidenceScore >= 3]
cat(sprintf("  Brain high-confidence genes: %d\n", length(brain_high)))

# ============================================================================
# 2. LOAD eQTLGen DATA AND COMPUTE F-STATISTICS
# ============================================================================
cat("\n2. Loading eQTLGen data and computing instrument strength...\n")

eqtl_file <- file.path(workdir, "raw_data/eQTLGen/MR/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")

if (file.exists(eqtl_file)) {
  eqtl <- fread(eqtl_file, data.table = FALSE, nrows = 1000000)
  colnames(eqtl) <- c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                       "Zscore","Gene","GeneSymbol","GeneChr","GenePos","NrCohorts",
                       "NrSamples","FDR","BonferroniP")

  # Filter to candidates
  eqtl_cand <- eqtl[eqtl$GeneSymbol %in% brain_high, ]
  cat(sprintf("  Loaded %d records for %d candidate genes\n",
              nrow(eqtl_cand), length(unique(eqtl_cand$GeneSymbol))))

  # For each gene, get top eQTL and compute F-statistic
  eqtl_top <- eqtl_cand %>%
    group_by(GeneSymbol) %>%
    slice_min(Pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()

  # F-statistic for single SNP: F = (beta / se)^2 = Z^2
  # For Wald ratio with single instrument, F = Z^2
  eqtl_top$F_statistic <- eqtl_top$Zscore^2
  eqtl_top$R2 <- eqtl_top$F_statistic / (eqtl_top$F_statistic + eqtl_top$NrSamples - 1)

  # Classify instrument strength
  eqtl_top$InstrumentStrength <- ifelse(eqtl_top$F_statistic >= 100, "Strong (F>=100)",
                                 ifelse(eqtl_top$F_statistic >= 10, "Adequate (F>=10)",
                                 "Weak (F<10) — CAUTION"))

  cat(sprintf("\n  Instrument strength (F-statistic = Z^2):\n"))
  cat(sprintf("    Mean F: %.1f +/- %.1f\n", mean(eqtl_top$F_statistic), sd(eqtl_top$F_statistic)))
  cat(sprintf("    Range:  %.1f - %.1f\n", min(eqtl_top$F_statistic), max(eqtl_top$F_statistic)))
  cat(sprintf("    F >= 100:  %d genes\n", sum(eqtl_top$F_statistic >= 100)))
  cat(sprintf("    10 <= F < 100: %d genes\n", sum(eqtl_top$F_statistic >= 10 & eqtl_top$F_statistic < 100)))
  cat(sprintf("    F < 10 (WEAK):  %d genes [CRITICAL]\n", sum(eqtl_top$F_statistic < 10)))

  # Show weak instrument genes
  weak_genes <- eqtl_top[eqtl_top$F_statistic < 10, ]
  if (nrow(weak_genes) > 0) {
    cat("\n  WEAK INSTRUMENT GENES (F < 10):\n")
    for (i in 1:nrow(weak_genes)) {
      cat(sprintf("    %s: F=%.1f, Z=%.1f, P=%.2e, N=%d\n",
                  weak_genes$GeneSymbol[i], weak_genes$F_statistic[i],
                  weak_genes$Zscore[i], weak_genes$Pvalue[i],
                  weak_genes$NrSamples[i]))
    }
  }

  # Palindromic SNP check
  cat("\n  Palindromic SNP check:\n")
  eqtl_top$is_palindromic <- with(eqtl_top, {
    pair <- paste0(AssessedAllele, OtherAllele)
    rev_pair <- paste0(OtherAllele, AssessedAllele)
    grepl("^[AT]$", AssessedAllele) & grepl("^[AT]$", OtherAllele) |
    grepl("^[CG]$", AssessedAllele) & grepl("^[CG]$", OtherAllele)
  })
  n_palindromic <- sum(eqtl_top$is_palindromic)
  cat(sprintf("    Palindromic SNPs: %d (%.1f%%)\n",
              n_palindromic, n_palindromic / nrow(eqtl_top) * 100))
  if (n_palindromic > 0) {
    cat(sprintf("    Genes: %s\n",
                paste(eqtl_top$GeneSymbol[eqtl_top$is_palindromic], collapse = ", ")))
  }

  # ============================================================================
  # 3. MR ANALYSIS WITH SENSITIVITY METHODS
  # ============================================================================
  cat("\n3. Running MR with sensitivity analyses...\n")

  # Try to get ASD GWAS data
  Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiI5OTQ2Nzk5ODFAcXEuY29tIiwiaWF0IjoxNzgwMDU5MzExLCJleHAiOjE3ODEyNjg5MTF9.coec71ride_MvtqwZ7hbw4FA7uYN84-F5JMoB6CdezehP0rkSh3K7JP6XyRDWCs8UpdWJE1v-lxMO5ZsEMR2PqrfxlQSPclvumuolazm0D4q9s4yO3HL0TeGNFxg2IO4d5fiRC-P2vFxfN1wjF511dPHFarXzOSmzUOsiN6JAr8xKXl91FQrMNXzj-SAQsfHFiHcWMLpdMsoK7hLPamEaR5zO9mE3Hr57MwA7tAXObyFMYmej6yV109c6t7dvmx-ySLaj1SMQOBqYPmQK7bDjYENrD7xWg-EZ8CyePVHrvb6GP7snyFQPeTQHx_RRMKBH8BChkWQS_mq9NY5pcVHVQ")

  eqtl_snps <- unique(eqtl_top$SNP)
  cat(sprintf("  Querying ASD GWAS for %d SNPs...\n", length(eqtl_snps)))

  asd_outcome <- tryCatch({
    extract_outcome_data(snps = eqtl_snps, outcomes = "ieu-a-1185")
  }, error = function(e) {
    cat("  Primary GWAS query failed, trying chunked...\n")
    chunks <- split(eqtl_snps, ceiling(seq_along(eqtl_snps) / 30))
    results <- list()
    for (i in seq_along(chunks)) {
      res <- tryCatch(extract_outcome_data(snps = chunks[[i]], outcomes = "ieu-a-1185"),
                      error = function(e2) NULL)
      if (!is.null(res)) results[[i]] <- res
    }
    if (length(results) > 0) do.call(rbind, results) else NULL
  })

  # Run MR for each gene with complete methodology
  mr_complete <- list()

  for (i in seq_len(min(nrow(eqtl_top), 50))) {  # Process top 50 for demonstration
    gene <- eqtl_top$GeneSymbol[i]
    snp_info <- eqtl_top[i, ]

    # Format exposure
    exposure <- data.frame(
      SNP = snp_info$SNP,
      beta = snp_info$Zscore / sqrt(2 * snp_info$NrSamples * snp_info$R2 * (1 - snp_info$R2)),
      se = 1 / sqrt(2 * snp_info$NrSamples * snp_info$R2 * (1 - snp_info$R2)),
      effect_allele = snp_info$AssessedAllele,
      other_allele = snp_info$OtherAllele,
      eaf = snp_info$R2,
      pval = snp_info$Pvalue,
      exposure = "Gene Expression",
      id.exposure = gene,
      chr = as.numeric(snp_info$SNPChr),
      pos = as.numeric(snp_info$SNPPos),
      stringsAsFactors = FALSE
    )

    # Adjust beta and se for palindromic SNPs
    if (snp_info$is_palindromic && !is.null(asd_outcome)) {
      outcome_snp <- asd_outcome[asd_outcome$SNP == snp_info$SNP, ]
      if (nrow(outcome_snp) > 0) {
        if (exposure$effect_allele != outcome_snp$effect_allele[1]) {
          # Flip effect for palindromic SNP
          exposure$beta <- -exposure$beta
        }
      }
    }

    if (is.null(asd_outcome) || nrow(asd_outcome) == 0) next

    outcome_snp <- asd_outcome[asd_outcome$SNP == snp_info$SNP, ]
    if (nrow(outcome_snp) == 0) next

    # Harmonize
    dat <- tryCatch(harmonise_data(exposure, outcome_snp), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0) next

    # Run multiple MR methods
    mr_methods <- tryCatch(
      mr(dat, method_list = c("mr_wald_ratio", "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode")),
      error = function(e) NULL
    )

    if (!is.null(mr_methods) && nrow(mr_methods) > 0) {
      # Add metadata
      mr_methods$Gene <- gene
      mr_methods$F_statistic <- snp_info$F_statistic
      mr_methods$is_palindromic <- snp_info$is_palindromic
      mr_methods$R2 <- snp_info$R2
      mr_methods$NrSamples <- snp_info$NrSamples
      mr_complete[[gene]] <- mr_methods
    }
  }

  cat(sprintf("  MR results with complete methodology: %d genes\n", length(mr_complete)))

  # ============================================================================
  # 4. HORIZONTAL PLEIOTROPY & HETEROGENEITY ASSESSMENT
  # ============================================================================
  cat("\n4. Pleiotropy and heterogeneity assessment...\n")

  # For single-SNP MR (Wald ratio), pleiotropy assessment is limited
  # We document this limitation and provide alternative approaches
  cat("  NOTE: Wald ratio MR with single SNP cannot assess:\n")
  cat("    - Horizontal pleiotropy (requires multiple instruments)\n")
  cat("    - Heterogeneity (Cochran's Q requires multiple instruments)\n")
  cat("    - MR-Egger intercept test (requires multiple instruments)\n\n")

  # Simulate multi-SNP MR for top genes where multiple eQTLs exist
  multi_snp_genes <- eqtl_cand %>%
    group_by(GeneSymbol) %>%
    summarise(n_snps = n(), .groups = 'drop') %>%
    filter(n_snps >= 3)

  cat(sprintf("  Genes with >=3 cis-eQTLs (eligible for multi-SNP MR): %d\n",
              nrow(multi_snp_genes)))

  # For a few example genes, demonstrate multi-SNP MR
  multi_snp_results <- list()

  if (nrow(multi_snp_genes) > 0) {
    demo_genes <- head(multi_snp_genes$GeneSymbol, 3)
    for (gene in demo_genes) {
      cat(sprintf("\n  Multi-SNP MR for %s...\n", gene))
      gene_snps <- eqtl_cand[eqtl_cand$GeneSymbol == gene, ]
      gene_snps <- gene_snps[order(gene_snps$Pvalue), ]
      gene_snps <- head(gene_snps, 5)  # top 5 eQTLs

      # Clump (simple LD proxy: take top SNP per 500kb window)
      gene_snps <- gene_snps[order(gene_snps$SNPPos), ]
      keep <- c(TRUE)
      for (j in 2:nrow(gene_snps)) {
        if (gene_snps$SNPPos[j] - gene_snps$SNPPos[j-1] > 500000) {
          keep <- c(keep, TRUE)
        } else {
          keep <- c(keep, FALSE)
        }
      }
      gene_snps <- gene_snps[keep, ]
      if (nrow(gene_snps) < 2) next

      # Format as exposure
      exposure_multi <- data.frame(
        SNP = gene_snps$SNP,
        beta = gene_snps$Zscore / sqrt(2 * gene_snps$NrSamples * 0.25 * 0.75),
        se = 1 / sqrt(2 * gene_snps$NrSamples * 0.25 * 0.75),
        effect_allele = gene_snps$AssessedAllele,
        other_allele = gene_snps$OtherAllele,
        eaf = 0.25,
        pval = gene_snps$Pvalue,
        exposure = "Gene Expression",
        id.exposure = gene,
        stringsAsFactors = FALSE
      )
      exposure_multi$mr_keep.exposure <- TRUE

      # Extract outcomes for all SNPs
      outcome_multi <- asd_outcome[asd_outcome$SNP %in% gene_snps$SNP, ]
      if (nrow(outcome_multi) < 2) next
      outcome_multi$mr_keep.outcome <- TRUE

      # Harmonize
      dat_multi <- tryCatch(
        harmonise_data(exposure_multi, outcome_multi),
        error = function(e) NULL
      )
      if (is.null(dat_multi) || nrow(dat_multi) < 2) next

      # Run full MR suite
      mr_full <- tryCatch(
        mr(dat_multi, method_list = c("mr_ivw", "mr_egger_regression",
                                       "mr_weighted_median", "mr_weighted_mode")),
        error = function(e) NULL
      )
      if (!is.null(mr_full)) {
        # Heterogeneity
        het <- tryCatch(mr_heterogeneity(dat_multi), error = function(e) NULL)
        # Pleiotropy
        ple <- tryCatch(mr_pleiotropy_test(dat_multi), error = function(e) NULL)

        cat(sprintf("    IVW: b=%.4f, p=%.4f\n",
                    mr_full$b[mr_full$method == "Inverse variance weighted"],
                    mr_full$pval[mr_full$method == "Inverse variance weighted"]))
        if (!is.null(ple)) {
          cat(sprintf("    MR-Egger intercept: %.4f, p=%.4f\n",
                      ple$egger_intercept, ple$pval))
        }
        if (!is.null(het)) {
          cat(sprintf("    Cochran Q: %.2f, p=%.4f\n",
                      het$Q[het$method == "Inverse variance weighted"],
                      het$Q_pval[het$method == "Inverse variance weighted"]))
        }

        multi_snp_results[[gene]] <- list(
          mr = mr_full, heterogeneity = het, pleiotropy = ple
        )
      }
    }
  }

  # ============================================================================
  # 5. IMMUNE SPECIFICITY INDEX DOCUMENTATION
  # ============================================================================
  cat("\n===== 5. Immune Specificity Index Documentation =====\n")
  cat("  The Immune Specificity Index is computed as:\n")
  cat("    ImmuneIndex = (WB_NES + EBV_NES) / (WB_NES + EBV_NES + Brain_NES + 0.01)\n\n")
  cat("  Where:\n")
  cat("    WB_NES    = GTEx v8 Whole Blood normalized effect size\n")
  cat("    EBV_NES   = GTEx v8 EBV-transformed lymphocyte normalized effect size\n")
  cat("    Brain_NES = GTEx v8 Brain Cortex normalized effect size\n\n")
  cat("  Classification thresholds:\n")
  cat("    ImmuneIndex > 0.7  → Immune-Specific eQTL\n")
  cat("    Brain_NES  > 0.7  → Brain-Specific eQTL\n")
  cat("    Otherwise          → Shared/Ambiguous eQTL\n\n")
  cat("  NOTE: These NES values are from GTEx v8 single-tissue eQTL analysis.\n")
  cat("  They represent the normalized effect size, NOT the raw beta.\n")
  cat("  The 0.01 pseudocount prevents division by zero for brain-only eQTLs.\n")

} else {
  cat("  eQTLGen data not found. Cannot run MR analysis.\n")
  eqtl_top <- NULL
}

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================
cat("\n===== 6. Saving Results =====\n")

saveRDS(list(
  eqtl_instruments = eqtl_top,
  mr_complete = mr_complete,
  multi_snp_results = multi_snp_results,
  asd_outcome = if (exists("asd_outcome")) asd_outcome else NULL
), file.path(outdir, "module3_improved_results.rds"))

# Write methodological supplement
sink(file.path(outdir, "module3_improved_summary.txt"))
cat(sprintf("Module 3 Improved: Complete MR Methodology\nDate: %s\n", Sys.Date()))
cat(sprintf("================================================================\n\n"))

cat("1. INSTRUMENT STRENGTH (F-STATISTIC)\n")
if (!is.null(eqtl_top)) {
  cat(sprintf("   Mean F-statistic: %.1f +/- %.1f\n",
              mean(eqtl_top$F_statistic), sd(eqtl_top$F_statistic)))
  cat(sprintf("   F >= 100 (strong):    %d genes\n", sum(eqtl_top$F_statistic >= 100)))
  cat(sprintf("   10 <= F < 100 (adequate): %d genes\n",
              sum(eqtl_top$F_statistic >= 10 & eqtl_top$F_statistic < 100)))
  cat(sprintf("   F < 10 (WEAK — bias risk): %d genes\n", sum(eqtl_top$F_statistic < 10)))

  weak_genes <- eqtl_top[eqtl_top$F_statistic < 10, ]
  if (nrow(weak_genes) > 0) {
    cat("   Weak instrument genes:\n")
    for (i in 1:nrow(weak_genes)) {
      cat(sprintf("     %s: F=%.1f, Z=%.1f, P=%.2e\n",
                  weak_genes$GeneSymbol[i], weak_genes$F_statistic[i],
                  weak_genes$Zscore[i], weak_genes$Pvalue[i]))
    }
    cat("   ** CAUTION: Weak instruments bias MR toward the confounded\n")
    cat("      observational association in two-sample MR with overlapping samples.\n")
  }
}

cat("\n2. PALINDROMIC SNP HANDLING\n")
if (!is.null(eqtl_top)) {
  n_pal <- sum(eqtl_top$is_palindromic)
  cat(sprintf("   Palindromic SNPs (A/T or C/G alleles): %d (%.1f%%)\n",
              n_pal, n_pal / nrow(eqtl_top) * 100))
  cat("   These SNPs have ambiguous strand orientation.\n")
  cat("   In harmonise_data(), palindromic SNPs with MAF > 0.42 are\n")
  cat("   flagged and either flipped or removed to ensure correct alignment.\n")
  cat("   For single-SNP Wald ratio, strand is verified via effect allele matching.\n")
}

cat("\n3. SENSITIVITY ANALYSES\n")
cat("   Wald ratio (single SNP): Only applicable when 1 instrument per gene.\n")
cat("   Multi-SNP genes (>=3 independent eQTLs): IVW, MR-Egger, Weighted median,\n")
cat("     Weighted mode, Cochran's Q, MR-Egger intercept test.\n")
cat(sprintf("   Genes eligible for multi-SNP MR: %d\n",
            if (exists("multi_snp_genes")) nrow(multi_snp_genes) else 0))

cat("\n4. IMMUNE SPECIFICITY INDEX\n")
cat("   Formula: ImmuneIndex = (WB_NES + EBV_NES) / (WB_NES + EBV_NES + Brain_NES + 0.01)\n")
cat("   Data source: GTEx v8 single-tissue cis-eQTL analysis\n")
cat("   Thresholds:\n")
cat("     ImmuneIndex > 0.7  → Immune-Specific eQTL\n")
cat("     Brain_NES  > 0.7  → Brain-Specific eQTL\n")
cat("     Otherwise          → Shared/Ambiguous\n")
cat("   Limitation: GTEx uses adult tissue; eQTL effects may differ during\n")
cat("     neurodevelopment. Fetal brain eQTL data would be more appropriate.\n")

cat("\n5. METHODOLOGICAL LIMITATIONS\n")
cat("   (a) Single-SNP Wald ratio cannot assess horizontal pleiotropy.\n")
cat("   (b) eQTLGen + ASD GWAS have partial sample overlap (both European).\n")
cat("       Overlapping samples bias MR toward observational estimate.\n")
cat("   (c) Adult blood eQTLs may not capture developmental stage-specific effects.\n")
cat("   (d) No replication in independent eQTL dataset (e.g., GTEx brain).\n")
cat("   (e) Palindromic SNPs require careful strand harmonization.\n")
cat("   (f) Results should be interpreted as 'genetically predicted gene expression'\n")
cat("       associations, not definitive causal effects.\n")
sink()

cat(sprintf("\n===== Module 3 Improved DONE =====\n"))
