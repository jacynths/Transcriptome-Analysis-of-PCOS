# Transcriptome-Analysis-of-PCOS
## Differential Gene Expression and Functional Enrichment Analysis of GEO Dataset GSE277906

Polycystic Ovary Syndrome (PCOS) is the most prevalent endocrine disorder in women of reproductive age, affecting approximately 5–15% of the global female population, yet its molecular underpinnings remain incompletely understood. This project undertakes a comprehensive transcriptome-level investigation of PCOS using publicly available **RNA-seq-derived FPKM expression data** from GEO accession GSE277906.

### About GEO GSE277906
These are mRNA expression profiles of human cumulus cells from patients with PCOS and non-PCOS women. It contains **40 RNA expression profiles** (23 PCOS and 17 control) pooled seperately via Illumina Novaseq 6000. 

This dataset is cited as: _Chen Y, Xie M, Wu S, Deng Z et al. Multi-omics approach to reveal follicular metabolic changes and their effects on oocyte competence in PCOS patients. Front Endocrinol (Lausanne) 2024;15:1426517. PMID: 39464191._ 

### Workflow
The workflow encompasses three analytical phases:
1. Preprocessing and QC in Python (including low-expression gene filtering, log2transformation, PCA, and sample correlation analysis)
2. Differential Gene Expression Analysis in R using the limma linear modelling framework (thresholds: |log FC| > 0.3 and raw p-value < 0.05)
   *applied in lieu of FDR due to the attenuated dynamic range characteristic of FPKM data
4. Functional Enrichment analysis using clusterProfiler against GO, KEGG, and Reactome databases

### Results
After filtering, **15,834** of 28,001 genes were retained. Differential analysis identified **157 significant DEGs** — 141 upregulated and 16 downregulated in PCOS relative to controls. PCA revealed a partial but discernible separation of the two cohorts along PC1. Enrichment analysis implicated **inflammatory signalling, steroid hormone biosynthesis, and insulin response pathways** among upregulated genes, with downregulated genes clustering in DNA repair and cell cycle regulation. 

#### This repository contains Python and R scripts used for data preprocessing, DEG and functional enrinchment analysis respectively. 
