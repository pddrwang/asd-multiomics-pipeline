# Module 7: Functional Analysis — Corrected Features
suppressMessages({ library(ggplot2) })
workdir <- "~/ASD_multiomics"
outdir  <- file.path(workdir, "module7")

cat("===== Module 7: Functional Analysis =====\n\n")

mod5 <- readRDS(file.path(workdir, "module5/module5_results.rds"))
core_genes <- mod5$features
cat(sprintf("Core genes (%d):\n  %s\n\n", length(core_genes), paste(core_genes, collapse=", ")))

# ---- Pathway enrichment ----
asd_pathways <- list(
  Immune_Inflammatory = c("CXCL8","CXCL10","CCL2","IL1B","IL6","TNF","NFKB1","NFKBIA","TLR2","TLR4","IFNG","IFNGR1","STAT1","STAT3","JAK2","C1QA","C1QB","C3","C5AR1","CX3CR1","ITGAM","CD14","FCGR3A"),
  Synaptic_Transmission = c("DLG4","SYP","SYN1","SNAP25","GRIN1","GRIN2A","GRIN2B","GABRA1","GABRB2","GAD1","GAD2","SHANK3","NLGN3","NLGN4X","NRXN1","CNTNAP2","CACNA1C","SCN2A","KCNMA1","BDNF"),
  Mitochondrial = c("NDUFS1","NDUFS4","SDHA","SDHB","COX5A","COX6C","ATP5A1","ATP5B","TFAM","POLG","SOD2","SIRT1","PPARGC1A","MFN1","MFN2","DNM1L","OPA1"),
  Autophagy_Proteostasis = c("SQSTM1","LAMP1","LAMP2","ATG5","ATG7","BECN1","MAP1LC3B","UBB","UBC","UBE3A","PARK2","HSPA1A","HSPA1B","HSP90AA1","DNAJB1","HSPB1"),
  Chromatin_Epigenetics = c("MECP2","HDAC1","HDAC2","EHMT1","EHMT2","DNMT1","DNMT3A","SETDB1","CHD8","ARID1B","SMARCA4","KMT2A","KMT2C","KDM5B"),
  Cell_Cycle_Development = c("CDKN1A","CDKN1B","CCND1","CCNE1","TP53","RB1","MYC","FOS","JUN","EGR1","CTNNB1","NOTCH1","HES1","SOX2","PAX6"))

enrichment <- list()
for (pw_name in names(asd_pathways)) {
  pw_genes <- asd_pathways[[pw_name]]
  overlap <- intersect(core_genes, pw_genes)
  if (length(overlap) > 0) {
    fold_enrich <- (length(overlap) / length(core_genes)) / (length(pw_genes) / 20000)
    enrichment[[pw_name]] <- list(overlap=overlap, pct=length(overlap)/length(core_genes)*100, fe=fold_enrich)
    cat(sprintf("  %s: %d genes (%.1f%%, %.1fx)\n    %s\n", pw_name, length(overlap), length(overlap)/length(core_genes)*100, fold_enrich, paste(overlap, collapse=", ")))
  } else { cat(sprintf("  %s: none\n", pw_name)) }
}

# ---- Cell-type mapping ----
ct_markers <- list(
  Microglia = c("CX3CR1","TMEM119","P2RY12","CSF1R","ITGAM","TREM2","CD68","AIF1"),
  Astrocytes = c("GFAP","AQP4","SLC1A3","ALDH1L1","GJA1","S100B","GLUL"),
  Excitatory_Neurons = c("SLC17A7","SATB2","CUX2","RORB","FEZF2","TBR1","CAMK2A"),
  Inhibitory_Neurons = c("GAD1","GAD2","SST","PVALB","VIP","LHX6","CALB2"),
  Oligodendrocytes = c("MOBP","MOG","PLP1","MBP","OLIG2","SOX10","MAG"),
  Endothelial = c("CLDN5","PECAM1","FLT1","VWF","CDH5","ENG"),
  OPC = c("PDGFRA","CSPG4","VCAN","NG2","OLIG1"))
ct_hits <- list()
cat("\nCell-type mapping:\n")
for (ct in names(ct_markers)) {
  m <- intersect(ct_markers[[ct]], core_genes)
  ct_hits[[ct]] <- m
  cat(sprintf("  %s: %s\n", ct, if(length(m)>0) paste(m, collapse=", ") else "-"))
}

# ASD GWAS overlap
asd_gwas <- c("CHD8","SCN2A","SYNGAP1","SHANK3","GRIN2B","ADNP","ANK2","ARID1B","ASH1L","CHD2","CTNNB1","DDX3X","DSCAM","DYRK1A","FOXP1","KMT2A","KMT5B","MECP2","MED13L","NAA15","POGZ","PTEN","SETD5","SUV420H1","TBR1","TCF4","TRIO","WDFY3","CACNA1C","CUL3","DEAF1","KATNAL2","KDM5B","KDM6B","NRXN1","PHF2","RELN","RORB","SLC6A1","SPAST","STXBP1","UBE3A","NCKAP1","FMR1","NLGN3","NLGN4X","CNTNAP2","GABRB3","PTCHD1")
gwas_ol <- intersect(core_genes, asd_gwas)
cat(sprintf("\nGWAS overlap: %d/%d: %s\n", length(gwas_ol), length(core_genes), paste(gwas_ol, collapse=", ")))

# Gene annotations
ann <- c(
  PSPH="磷酸丝氨酸磷酸酶/丝氨酸合成", HOXD1="同源盒转录因子/发育模式",
  MCL1="抗凋亡BCL2家族/细胞存活", PEG3_AS1="lncRNA/父源表达印记",
  DDX6="RNA解旋酶/转录后调控", TFAP2A="AP-2转录因子/神经嵴发育",
  BRD3="溴结构域蛋白/染色质乙酰化阅读器", STK17B="丝苏氨酸激酶/DNA损伤应答",
  UQCRB="线粒体复合物III/氧化磷酸化", CCK="胆囊收缩素/神经肽与突触传递",
  MALAT1="lncRNA/核散斑与剪接调控", MSN="Moesin/ERM细胞骨架/免疫突触",
  OLFML3="Olfactomedin样/分泌糖蛋白", RBPMS="RNA结合蛋白/神经元mRNA转运",
  IL22="白介素22/黏膜免疫", TLN2="Talin 2/突触后骨架锚定",
  CALD1="Caldesmon/细胞骨架调控", TBX3="T-box转录因子/发育多能性",
  NRP1="Neuropilin-1/轴突导向与血管生成", DSP="Desmoplakin/桥粒与细胞黏附")
cat("\nAnnotations:\n")
for (g in core_genes) cat(sprintf("  %-10s %s\n", g, if(g %in% names(ann)) ann[[g]] else "?"))

saveRDS(list(core_genes=core_genes, enrichment=enrichment, cell_types=ct_hits, gwas_overlap=gwas_ol, annotations=ann), file.path(outdir, "module7_results.rds"))

sink(file.path(outdir, "module7_summary.txt"))
cat(sprintf("Module 7: Functional Analysis\nDate: %s\n\nCore genes (%d): %s\n\n", Sys.Date(), length(core_genes), paste(core_genes, collapse=", ")))
cat("Pathway enrichment:\n")
for (nm in names(enrichment)) { e <- enrichment[[nm]]; cat(sprintf("  %s: %d genes (%.1f%%, %.1fx): %s\n", nm, length(e$overlap), e$pct, e$fe, paste(e$overlap, collapse=", "))) }
cat("\nCell-type expression:\n")
for (ct in names(ct_hits)) cat(sprintf("  %s: %s\n", ct, if(length(ct_hits[[ct]])>0) paste(ct_hits[[ct]], collapse=", ") else "-"))
cat(sprintf("\nASD GWAS overlap: %d/%d: %s\n", length(gwas_ol), length(core_genes), paste(gwas_ol, collapse=", ")))
cat("\nGene annotations:\n")
for (g in core_genes) cat(sprintf("  %-10s %s\n", g, if(g %in% names(ann)) ann[[g]] else "?"))
sink()
cat("\n===== Module 7 DONE =====\n")
