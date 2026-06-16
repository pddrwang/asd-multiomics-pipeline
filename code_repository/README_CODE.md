# ASD Multi-Omics Integration — Analysis Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20715220.svg)](https://doi.org/10.5281/zenodo.20715220)

**Data (Zenodo):** https://doi.org/10.5281/zenodo.20715220

**Manuscript:** Multi-Tissue Transcriptomic Integration and Mendelian
Randomization Identify a Myeloid-Mediated Immune Subtype in Autism
Spectrum Disorder

## Overview

This repository contains the complete R analysis pipeline for an
11-module multi-omics integration study of autism spectrum disorder
(ASD). The pipeline integrates five GEO microarray datasets (n=438)
spanning four tissue types across three developmental stages.

## Pipeline Modules

| Module | Description | Key Method |
|--------|-------------|------------|
| 01 | Data Integration | ComBat batch correction |
| 02 | Brain DEG Analysis | Platform-conditioned permutation test |
| 03 | Causal Inference | Two-sample Mendelian randomization (eQTLGen → ASD GWAS) |
| 04 | Molecular Subtyping | ConsensusClusterPlus (K=2) |
| 05 | Predictive Modeling | Nested Random Forest (3-stage feature selection) |
| 06 | Model Validation | Leave-one-dataset-out CV + Decision curve analysis |
| 07 | Functional Annotation | Pathway enrichment (KEGG/GO/Reactome) |
| 08 | Brain Cell-Type Overlap | snRNA-seq marker enrichment (Velmeshev 2019) |
| 09 | Immune Deconvolution | ssGSEA + LM22 immune signatures |
| 10 | Tissue-Stratified MR | GTEx v8 cis-eQTL tissue specificity |
| 11 | Regulatory Integration | Multi-omic TLN2/CCK annotation |

## Quick Start

### 设置工作目录

所有脚本使用变量 `workdir <- "~/ASD_multiomics"` 指定项目根目录。
**首次使用前**，将每个脚本中的 `workdir` 修改为你的实际项目路径。

```r
# 每个模块开头类似：
workdir <- "~/ASD_multiomics"   # ← 改为你的路径
```

脚本期望的工作目录结构：
```
ASD_multiomics/
├── raw_data/GEO/          ← 下载原始 GEO 数据到此
├── module1/ → module11/   ← 各模块运行后产生输出
└── improvements/          ← 改进模块输出
```

### Prerequisites

```r
install.packages(c("BiocManager"))
BiocManager::install(c("limma", "sva", "GSVA", "GEOquery"))
install.packages(c("randomForest", "glmnet", "pROC",
                   "ConsensusClusterPlus", "TwoSampleMR"))
```

### Execution Order

Run modules sequentially (01 → 11). Each module reads outputs from
prior modules. See `00_RUN_ORDER.txt` for details.

```r
source("01_module1_integration/module1.R")
source("02_module2_brain_DEG/module2.R")
# ... continue through module 11
```

### Input Data

Due to file size constraints, raw GEO datasets are not included in this
repository. The pipeline expects raw data in a `00_raw_data/GEO/`
directory. See the accompanying data repository (Zenodo DOI above) for
the complete processed data package including the ComBat-corrected
unified expression matrix.

## Improvement Modules

The `improvements/` directory contains supplementary analyses added
during the peer review process, addressing:
- ComBat batch correction diagnostics
- Permutation test calibration and statistical power analysis
- MR methodological completeness
- Predictive model audit (AUC gap decomposition)
- External validation framework
- Benjamini-Hochberg correction verification

## Repository Structure

```
├── 01_module1_integration/
├── 02_module2_brain_DEG/
├── 03_module3_MR/
├── 04_module4_subtyping/
├── 05_module5_prediction/
├── 06_module6_validation/
├── 07_module7_functional/
├── 08_module8_celltype/
├── 09_module9_immune/
├── 10_module10_eQTL_MR/
├── 11_module11_regulation/
├── improvements/
│   ├── 00_run_all.R
│   ├── module1_comdat/
│   ├── module2_permutation/
│   ├── module3_MR/
│   ├── module5_prediction/
│   ├── module6_validation/
│   └── module8_celltype/
└── 00_RUN_ORDER.txt
```

## Citation

If you use this pipeline, please cite:
[Paper citation — to be added after publication]

## License

[To be specified]
