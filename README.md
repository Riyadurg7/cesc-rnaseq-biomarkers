# Identification of Prognostic and Diagnostic Biomarkers in Cervical Cancer (TCGA-CESC) by RNA-seq Analysis

MSc Bioinformatics mini-project. The pipeline follows the TCGA RNA-seq workflow of Wadapurkar et al. (2023, ovarian cancer) and applies it to cervical cancer (TCGA-CESC), with added cross-method validation (DESeq2), survival modelling (Cox regression, ROC), and external validation (GEO, GEPIA2, KM-plotter).

## Aim

Find genes whose expression differs between clinically meaningful patient groups in cervical cancer, and test whether any of them work as biomarkers.

Three comparisons on 304 TCGA-CESC tumours:

| Comparison | Groups | Significant genes (edgeR, \|logFC\| > 2, FDR < 0.05) |
|---|---|---|
| Stage | Early (I-II) vs Late (III-IV) | 453 (188 up in late, 265 down) |
| Histology | Squamous vs Adenocarcinoma | 3,415 (1,849 up in adeno, 1,566 up in squamous) |
| Age | Younger vs Older at diagnosis | 772 (549 up in older, 223 down) |

## Key findings

- **EGFR** is the headline gene. It is a **prognostic** biomarker (multivariable Cox: HR = 1.40 per doubling of expression, 95% CI 1.15-1.71, p = 0.00098, adjusted for stage and histology, concordance 0.695) and a **diagnostic typing** biomarker (ROC AUC = 0.89 for squamous vs adenocarcinoma; ~4-fold higher in squamous; DESeq2 padj = 3.8e-26). It is **not** a detection biomarker - GEPIA2 shows tumour = normal expression. EGFR is druggable (cetuximab), which gives the finding clinical relevance.
- **30 unique hub genes** were identified from STRING/Cytoscape networks (top-20 by Degree ∩ Betweenness per comparison). **14 of 30 were confirmed at full stringency by an independent method (DESeq2 run from raw counts)**: KRT20, CDH17, EGFR, MUC2, LTF, CHGA, PAX3, FGF3, MUC5AC, CFTR, CEACAM5, EPCAM, FOXJ1, KRT5. Main claims rest only on these.
- Method agreement (edgeR vs DESeq2) tracks signal strength: histology 61% > stage 16% > age 10%; direction agreement 82%. One candidate (ALDH8A1) was rejected because its signal came from a single extreme patient.
- A label-scrambling bug (patient labels paired by row position instead of barcode) was caught because the Cox result contradicted the KM curve. Fix: always match labels by barcode with `match()`. Kept in the project as a worked example of sanity-checking.

## Pipeline (9 R scripts, run in order)

| Script | What it does |
|---|---|
| `01_download.R` | Download TCGA-CESC RNA-seq counts + clinical data (TCGAbiolinks) |
| `02_preprocess.R` | Normalisation (gcContent) + filtering (quantile 0.25) → 45,169 genes × 304 samples |
| `03_dea_R.R` | Differential expression with edgeR for the three comparisons; volcano plots |
| `04_overlap.R` | Venn of the three gene lists; core set of 111 shared genes |
| `05_networkhubs.R` | Hub genes from STRING (confidence 0.40) + Cytoscape networks (Degree ∩ Betweenness) |
| `06_enrichment.R` | GO and KEGG enrichment (clusterProfiler) per comparison + core-111 |
| `07_survival_R.R` | Kaplan-Meier curves, multivariable Cox (with PH checks), ROC/AUC (pROC) |
| `08_deseq2_R.R` | Independent re-analysis from raw counts with DESeq2; hub confirmation check |
| `09_Geo_validation_R.R` | External validation on GEO microarray GSE63514 (normal vs cancer, squamous) |

## Folder structure

```
scripts/   the 9 analysis scripts (R)
figures/   all result plots (volcano, Venn, networks, GO/KEGG dotplots, KM, ROC, GEO)
results/   key tables: significant DEG lists, hub lists, dual-confirmed gene lists,
           Cox and ROC results, GO/KEGG enrichment tables
```

Raw data (GDCdata download, ~GBs) and intermediate `.rds` objects are not in the repo - `01_download.R` and `02_preprocess.R` recreate them.

## Requirements

R ≥ 4.5 with Bioconductor packages: TCGAbiolinks, edgeR, DESeq2, clusterProfiler, org.Hs.eg.db, GEOquery, hgu133plus2.db; CRAN packages: survival, survminer, pROC, ggplot2, VennDiagram. Network analysis uses STRING (string-db.org) and Cytoscape 3.x with MCODE.

## Data sources

TCGA-CESC (GDC portal, via TCGAbiolinks) · GEO GSE63514 · GEPIA2 · The Human Protein Atlas · KM-plotter.

## Reference

Wadapurkar R. et al. (2023). Analysis of next generation sequencing data for comprehensive characterization of high-grade serous ovarian carcinoma. (Methodology template applied here to TCGA-CESC.)
