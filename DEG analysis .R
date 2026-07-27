#Part II: Differential Gene Esxpression Analysis

#install libraries
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"))
install.packages(c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer"))

library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

#paths 
input_file  <- "log2_FPKM_filtered.tsv"
output_dir  <- "outputs"
dir.create(output_dir, showWarnings = FALSE)

#A. Load expression data 
cat("Loading log2 FPKM matrix...\n")
expr <- read.table(input_file, sep = "\t", header = TRUE,
                   row.names = 1, check.names = FALSE)
cat(sprintf("  Matrix: %d genes × %d samples\n", nrow(expr), ncol(expr)))

#B. Build sample metadata
samples <- colnames(expr)
group   <- factor(ifelse(grepl("^Control", samples), "Control", "PCOS"),
                  levels = c("Control", "PCOS"))   # Control = reference
cat(sprintf("  Controls: %d  |  PCOS: %d\n",
            sum(group == "Control"), sum(group == "PCOS")))

#C. Design matrix and limma model
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "PCOS_vs_Control")

fit  <- lmFit(as.matrix(expr), design)
fit2 <- eBayes(fit)

   #extract results
all_results <- topTable(fit2, coef = "PCOS_vs_Control",
                        number = Inf, adjust.method = "BH", sort.by = "P")
all_results$gene <- rownames(all_results)
cat(sprintf("\nTotal genes tested: %d\n", nrow(all_results)))

#D. Define DEGs
   #thresholds are |logFc| > 0.3 and raw p value + logFC threshold instead of FDR as pAdj. is high and logFC is lower than 1
FC_THRESH  <- 0.3
P_THRESH  <- 0.05

   #diagnostics of threshold
cat(sprintf("  logFC range    : %.3f to %.3f\n", min(all_results$logFC), max(all_results$logFC)))
cat(sprintf("  AveExp range    : %.3f to %.3f\n", min(all_results$AveExpr), max(all_results$AveExpr)))
cat(sprintf("  Raw P range    : %.4f to %.4f\n", min(all_results$P.Value), max(all_results$P.Value)))
cat(sprintf("  Adj P range    : %.4f to %.4f\n", min(all_results$adj.P.Val), max(all_results$adj.P.Val)))
cat(sprintf("  |logFC| > 1  : %d genes\n",     sum(abs(all_results$logFC) > 0.3)))
cat(sprintf("  raw P < 0.05   : %d genes\n",     sum(all_results$P.Value < 0.05)))
cat(sprintf("  Both combined  : %d genes\n", sum(abs(all_results$logFC) > FC_THRESH & all_results$P.Value < P_THRESH)))

all_results$sig <- "Not significant"
all_results$sig[all_results$logFC >  FC_THRESH & all_results$P.Value < P_THRESH] <- "Up in PCOS"
all_results$sig[all_results$logFC < -FC_THRESH & all_results$P.Value < P_THRESH] <- "Down in PCOS"

deg_table <- all_results[all_results$sig != "Not significant", ]
cat(sprintf("\nDEGs (|logFC|>0.3, raw P<0.05): %d total\n", nrow(deg_table)))
cat(sprintf("  Up in PCOS:   %d\n", sum(deg_table$sig == "Up in PCOS")))
cat(sprintf("  Down in PCOS: %d\n", sum(deg_table$sig == "Down in PCOS")))

deg_path <- file.path(output_dir, "DEG_results.tsv")
write.table(all_results, deg_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Full results saved → %s\n", deg_path))

sig_path <- file.path(output_dir, "DEGs_significant.tsv")
write.table(deg_table, sig_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Significant DEGs saved → %s\n", sig_path))

#E. Volcano Plot
cat("\nGenerating Volcano plot...\n")

 #label top 15 significant genes
top_labels <- head(deg_table[order(deg_table$P.Value), ], 15)

sig_colours <- c(
  "Up in PCOS"      = "magenta",
  "Down in PCOS"    = "blue",
  "Not significant" = "grey"
)

volcano_plot <- ggplot(all_results, aes(x = logFC, y = -log10(P.Value),
                                        colour = sig, size = sig)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(P_THRESH), linetype = "dashed",
             colour = "grey35", linewidth = 0.5) +
  geom_vline(xintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed",
             colour = "grey35", linewidth = 0.5) +
  geom_text_repel(data = top_labels,
                  aes(label = gene), colour = "black",
                  size = 2.8, max.overlaps = 20,
                  box.padding = 0.35, segment.size = 0.3) +
  scale_colour_manual(values = sig_colours, name = "Regulation") +
  scale_size_manual(values = c("Up in PCOS" = 1.8, "Down in PCOS" = 1.8,
                               "Not significant" = 0.8), guide = "none") +
  labs(title    = "Volcano Plot: PCOS vs Control",
       subtitle = sprintf("DEGs: %d up, %d down  (|logFC|>0.3, raw P<0.05)",
                          sum(deg_table$sig == "Up in PCOS"),
                          sum(deg_table$sig == "Down in PCOS")),
       x = "log2 Fold Change (PCOS / Control)",
       y = "-log10(Raw P-value)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(output_dir, "Volcano_plot.png"),
       volcano_plot, width = 7, height = 6, dpi = 150)
cat("  Saved → outputs/Volcano_plot.png\n")

#F. MA Plot
cat("Generating MA plot...\n")

 #AveExpr from limma is the average log2 expression across all samples
ma_plot <- ggplot(all_results, aes(x = AveExpr, y = logFC, colour = sig)) +
  geom_point(alpha = 0.5, size = 0.9) +
  geom_hline(yintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed",
             colour = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  scale_colour_manual(values = sig_colours, name = "Regulation") +
  labs(title = "MA Plot: PCOS vs Control",
       x = "Average log2 Expression (A)",
       y = "log2 Fold Change (M)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(output_dir, "MA_plot.png"),
       ma_plot, width = 7, height = 5, dpi = 150)
cat("  Saved → outputs/MA_plot.png\n")

#G. Heatmap
cat("Generating DEG heatmap...\n")

  #pick top 50 DEGs ranked by |logFC|
top50 <- head(deg_table[order(abs(deg_table$logFC), decreasing = TRUE), ], 50)
heat_mat <- as.matrix(expr[top50$gene, ])

  #Z-score each gene across samples for visualisation
heat_z <- t(scale(t(heat_mat)))

  #annotation bar for samples
anno_col <- data.frame(Group = group)
rownames(anno_col) <- samples
anno_colours <- list(Group = c(Control = "darkorange", PCOS = "darkred"))

png(file.path(output_dir, "DEG_Heatmap_top50.png"),
    width = 1600, height = 2000, res = 150)
pheatmap(heat_z,
         annotation_col  = anno_col,
         annotation_colors = anno_colours,
         show_colnames   = TRUE,
         show_rownames   = TRUE,
         fontsize_row    = 7,
         fontsize_col    = 7,
         color           = colorRampPalette(rev(brewer.pal(9, "PuOr")))(100),
         clustering_method = "ward.D2",
         main            = "Top 50 DEGs — Z-scored log2 FPKM\n(PCOS vs Control)")
dev.off()
cat("  Saved → outputs/DEG_Heatmap_top50.png\n")

cat("\n DEG analysis complete.\n")
