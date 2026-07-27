#Mini Project: DEG Analysis on PCOS dataset 
#Part III: Enrichment Analysis
#by Ananya V

#install libraries
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db","ReactomePA", "enrichplot"))
install.packages(c("ggplot2", "ggupset"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

#A. Loading Dataset
degs <- read.table("outputs/DEGs_significant.tsv", sep="\t", header=TRUE,
                   stringsAsFactors=FALSE)
cat("Rows loaded:", nrow(degs), "\n")
cat("Columns:", paste(colnames(degs), collapse=", "), "\n")
output_dir <- "outputs"
dir.create(output_dir, showWarnings=FALSE)

#B. Gene Mapping
id_map <- bitr(degs$gene, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat(sprintf("Mapped: %d / %d genes\n", nrow(id_map), nrow(degs)))

degs_id     <- merge(degs, id_map, by.x = "gene", by.y = "SYMBOL")
all_entrez  <- unique(degs_id$ENTREZID)
up_entrez   <- unique(degs_id$ENTREZID[degs_id$sig == "Up in PCOS"])
down_entrez <- unique(degs_id$ENTREZID[degs_id$sig == "Down in PCOS"])

#C. Helper Function Set Up
 #run enrichGo and then filter manually 
 run_enrichGO <- function(entrez_ids, ont_type, label) {
  cat(sprintf("\n-- GO %s: %s (%d genes) --\n", ont_type, label, length(entrez_ids)))
  
  res <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ont_type,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,       # retrieve ALL terms
                  qvalueCutoff  = 1,
                  minGSSize     = 3,
                  maxGSSize     = 500,
                  readable      = TRUE)
  
  df <- as.data.frame(res)
  cat(sprintf("  Total terms returned: %d\n", nrow(df)))
  
  # Manually filter on nominal p-value < 0.05
  df_sig <- df[df$pvalue < 0.05, ]
  cat(sprintf("  Significant (nominal p < 0.05): %d\n", nrow(df_sig)))
  cat(sprintf("  Of which p.adjust < 0.05:       %d\n", sum(df_sig$p.adjust < 0.05)))
  
  return(list(result = res, df = df, df_sig = df_sig))
}

 #save enrichment tsv
 save_tsv <- function(df_sig, filename) {
   if (nrow(df_sig) > 0) {
     write.table(df_sig, file.path(output_dir, filename),
                 sep = "\t", quote = FALSE, row.names = FALSE)
     cat(sprintf("  Saved → outputs/%s\n", filename))
   }
 }
  #dotplot from filtered result 
plot_enrichment <- function(res_obj, df_sig, title_text, filename,
                            w = 9, h = 8) {
  if (nrow(df_sig) == 0) return(invisible(NULL))
  
  # Subset the enrichResult object to significant terms only
  res_obj@result <- res_obj@result[res_obj@result$pvalue < 0.05, ]
  
  p <- dotplot(res_obj, showCategory = min(20, nrow(df_sig)), font.size = 9) +
    labs(title   = title_text,
         caption = "Nominal p < 0.05 | BH-adjusted p shown as colour scale") +
    theme(plot.title   = element_text(face = "bold", size = 11),
          plot.caption = element_text(size = 7, colour = "grey50"))
  
  ggsave(file.path(output_dir, filename), p, width = w, height = h, dpi = 150)
  cat(sprintf("  Plot → outputs/%s\n", filename))
}

#D. GO: Biological Process
bp      <- run_enrichGO(all_entrez, "BP", "All DEGs")
save_tsv(bp$df_sig, "GO_BP_enrichment.tsv")
plot_enrichment(bp$result, bp$df_sig,
                "GO Biological Process — PCOS vs Control",
                "GO_BP_dotplot.png")

 #enrichment map 
if (nrow(bp$df_sig) > 5) {
  bp$result@result <- bp$result@result[bp$result@result$pvalue < 0.05, ]
  bp2 <- pairwise_termsim(bp$result)
  p_emap <- emapplot(bp2, showCategory = 30) +
    labs(title = "GO BP Enrichment Map") +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(output_dir, "GO_BP_enrichment_map.png"),
         p_emap, width = 10, height = 9, dpi = 150)
  cat("  Plot → outputs/GO_BP_enrichment_map.png\n")
}
 
#map up and down genes seperately
bp_up <- run_enrichGO(up_entrez,   "BP", "Up in PCOS")
bp_dn <- run_enrichGO(down_entrez, "BP", "Down in PCOS")
save_tsv(bp_up$df_sig, "GO_BP_upregulated.tsv")
save_tsv(bp_dn$df_sig, "GO_BP_downregulated.tsv")

if (nrow(bp_up$df_sig) > 0 && nrow(bp_dn$df_sig) > 0) {
  top_up <- head(bp_up$df_sig[order(bp_up$df_sig$pvalue), ], 10)
  top_dn <- head(bp_dn$df_sig[order(bp_dn$df_sig$pvalue), ], 10)
  top_up$direction <- "Up in PCOS";   top_up$log10p <- -log10(top_up$pvalue)
  top_dn$direction <- "Down in PCOS"; top_dn$log10p <- -log10(top_dn$pvalue)
  
  comb <- rbind(top_up[, c("Description","direction","log10p")],
                top_dn[, c("Description","direction","log10p")])
  comb$Description <- factor(comb$Description, levels = rev(unique(comb$Description)))
  
  p_ud <- ggplot(comb, aes(x = log10p, y = Description, fill = direction)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = c("Up in PCOS"="hotpink","Down in PCOS"="orange")) +
    facet_wrap(~ direction, scales = "free_y", ncol = 1) +
    labs(title = "GO BP — Up vs Down DEGs",
         x = "-log10(Nominal P-value)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "none")
  ggsave(file.path(output_dir, "GO_BP_UpDown_barplot.png"),
         p_ud, width = 9, height = 8, dpi = 150)
  cat("  Plot → outputs/GO_BP_UpDown_barplot.png\n")
}

#E. GO: Molecular Function 
mf <- run_enrichGO(all_entrez, "MF", "All DEGs")
save_tsv(mf$df_sig, "GO_MF_enrichment.tsv")
plot_enrichment(mf$result, mf$df_sig,
                "GO Molecular Function — PCOS vs Control",
                "GO_MF_dotplot.png", w = 8, h = 7)

#F. KEGG
cat("\n-- KEGG Pathways --\n")
kegg_res <- enrichKEGG(gene          = all_entrez,
                       organism      = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,       # retrieve all, filter manually
                       qvalueCutoff  = 1,
                       minGSSize     = 3)
kegg_res  <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df   <- as.data.frame(kegg_res)
kegg_sig  <- kegg_df[kegg_df$pvalue < 0.05, ]
cat(sprintf("  KEGG terms nominal p<0.05: %d\n", nrow(kegg_sig)))
save_tsv(kegg_sig, "KEGG_enrichment.tsv")

if (nrow(kegg_sig) > 0) {
  kegg_res@result <- kegg_res@result[kegg_res@result$pvalue < 0.05, ]
  p_kegg <- dotplot(kegg_res, showCategory = min(20, nrow(kegg_sig)), font.size = 9) +
    labs(title   = "KEGG Pathways — PCOS vs Control",
         caption = "Nominal p < 0.05") +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.caption = element_text(size = 7, colour = "grey50"))
  ggsave(file.path(output_dir, "KEGG_dotplot.png"), p_kegg,
         width = 8, height = 8, dpi = 150)
  cat("  Plot → outputs/KEGG_dotplot.png\n")
}
#G. Reactome
cat("\n-- Reactome Pathways --\n")
react_res <- enrichPathway(gene          = all_entrez,
                           organism      = "human",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,    # retrieve all, filter manually
                           qvalueCutoff  = 1,
                           minGSSize     = 3,
                           readable      = TRUE)
react_df  <- as.data.frame(react_res)
react_sig <- react_df[react_df$pvalue < 0.05, ]
cat(sprintf("  Reactome terms nominal p<0.05: %d\n", nrow(react_sig)))
save_tsv(react_sig, "Reactome_enrichment.tsv")

if (nrow(react_sig) > 0) {
  react_res@result <- react_res@result[react_res@result$pvalue < 0.05, ]
  p_react <- dotplot(react_res, showCategory = min(20, nrow(react_sig)), font.size = 9) +
    labs(title   = "Reactome — PCOS vs Control",
         caption = "Nominal p < 0.05") +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.caption = element_text(size = 7, colour = "grey50"))
  ggsave(file.path(output_dir, "Reactome_dotplot.png"), p_react,
         width = 8, height = 8, dpi = 150)
  cat("  Plot → outputs/Reactome_dotplot.png\n")
}

cat("\n Enrichment Analysis complete.\n")