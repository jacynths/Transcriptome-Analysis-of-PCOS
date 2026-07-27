# Transcriptome-Analysis-of-PCOS
(Differential Gene Expression and Functional Enrichment Analysis of GEO Dataset GSE277906)

Polycystic Ovary Syndrome (PCOS) is the most prevalent endocrine disorder in women of reproductive age, affecting approximately 5–15% of the global female population, yet its molecular underpinnings remain incompletely understood. This project undertakes a comprehensive transcriptome-level investigation of PCOS using publicly available RNA-seq-derived FPKM expression data from GEO accession GSE277906, comprising 17 control and 23 PCOS patient samples. 

The workflow encompasses three analytical phases: 
(i) preprocessing and quality control in Python, including low-expression gene filtering, log2transformation, PCA, and sample correlation analysis
(ii) differential expression analysis in R using the limma linear modelling framework, with thresholds of |log FC| > 0.3 and raw p-value < 0.05 applied in lieu of FDR due to the attenuated dynamic range characteristic of FPKM data
(iii) functional enrichment analysis using clusterProfiler against GO, KEGG, and Reactome databases. 

After filtering, 15,834 of 28,001 genes were retained. Differential analysis identified 157 significant DEGs — 141 upregulated and 16 downregulated in PCOS relative to controls. 
PCA revealed a partial but discernible separation of the two cohorts along PC1. 
Enrichment analysis implicated inflammatory signalling, steroid hormone biosynthesis, and insulin response pathways among upregulated genes, with downregulated genes clustering in DNA repair and cell cycle regulation. 
