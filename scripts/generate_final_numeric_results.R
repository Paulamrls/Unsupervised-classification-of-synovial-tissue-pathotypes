# Generación final de resultados descriptivos para la solución jerárquica k=5.
# No recalcula ni modifica las asignaciones de cluster.

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
local_lib <- file.path(project_root, ".r_libs", "4.5")
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

source(file.path(project_root, "pipeline", "functions", "01_preprocessing.R"))

old_results <- file.path(project_root, "analysis_output", "results")
out <- file.path(project_root, "results")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)

assignments <- read.csv(file.path(old_results, "clusters", "cluster_assignments.csv"),
                        check.names = FALSE)
assignments <- assignments[, intersect(c("sample_id", "cluster_id",
                                          "pathotype_histological", "PC1", "PC2"),
                                        names(assignments)), drop = FALSE]
stopifnot(nrow(assignments) == 210L,
          !anyDuplicated(assignments$sample_id),
          setequal(sort(unique(assignments$cluster_id)), 1:5))
assignments$cluster_label <- factor(paste("Cluster", assignments$cluster_id),
                                    levels = paste("Cluster", 1:5))

palette <- c(
  "Cluster 1" = "#009E73",
  "Cluster 2" = "#E69F00",
  "Cluster 3" = "#0072B2",
  "Cluster 4" = "#CC79A7",
  "Cluster 5" = "#D55E00"
)

# PCA definitivo: mismas coordenadas y asignaciones jerárquicas, etiquetas neutras.
p_pca <- ggplot(assignments, aes(PC1, PC2, colour = cluster_label)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = palette, name = "Cluster") +
  labs(title = "PCA - hierarchical clustering (k = 5)",
       x = "PC1 (25.9%)", y = "PC2 (17.8%)") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(out, "figures", "pca_hierarchical_k5_numeric.png"),
       p_pca, width = 10, height = 7, dpi = 300, bg = "white")

# Composición histológica posterior al clustering, con IDs numéricos.
composition_counts <- as.data.frame(table(
  cluster_id = assignments$cluster_id,
  pathotype_histological = assignments$pathotype_histological
))
composition_counts <- composition_counts[composition_counts$Freq > 0, ]
composition_counts$pct <- ave(
  composition_counts$Freq,
  composition_counts$cluster_id,
  FUN = function(x) 100 * x / sum(x)
)
composition_counts$cluster_label <- factor(
  paste("Cluster", composition_counts$cluster_id),
  levels = rev(paste("Cluster", 1:5))
)

p_composition <- ggplot(
  composition_counts,
  aes(pathotype_histological, cluster_label, fill = pct)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(pct, 1), "%")), size = 4) +
  scale_fill_gradient(low = "white", high = "#2166AC",
                      name = "% de muestras\npor cluster") +
  labs(title = "Histological composition of transcriptomic clusters",
       x = "Histological pathotype", y = "Transcriptomic cluster") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

ggsave(file.path(out, "figures", "histological_composition_clusters.png"),
       p_composition, width = 10, height = 7, dpi = 300, bg = "white")
write.csv(composition_counts[, c("cluster_id", "pathotype_histological",
                                 "Freq", "pct")],
          file.path(out, "tables", "histological_composition_clusters.csv"),
          row.names = FALSE)

# Marcadores positivos: reutiliza exactamente los contrastes DESeq2 existentes.
old_file_labels <- c("Endothelial", "Fibroid", "Myeloid", "Lymphoid", "IFN_high")
names(old_file_labels) <- as.character(1:5)
lfc_threshold <- 0.5
fdr_threshold <- 0.05
all_deg <- list()
top20 <- list()
for (cluster_id in names(old_file_labels)) {
  df <- read.csv(file.path(old_results, "DEGs",
                           paste0("DEGs_ALL_", old_file_labels[[cluster_id]], ".csv")),
                 check.names = FALSE)
  if ("X" %in% names(df)) df$X <- NULL
  df$cluster_id <- as.integer(cluster_id)
  df$foldChange <- 2^df$log2FoldChange
  df$status <- "Not significant"
  df$status[df$padj < fdr_threshold & df$log2FoldChange >= lfc_threshold] <- "Up"
  df$status[df$padj < fdr_threshold & df$log2FoldChange <= -lfc_threshold] <- "Down"
  all_deg[[cluster_id]] <- df
  up <- df[df$status == "Up", , drop = FALSE]
  up <- up[order(up$padj, -up$log2FoldChange), , drop = FALSE]
  top20[[cluster_id]] <- head(up, 20)
}
top_table <- do.call(rbind, lapply(names(top20), function(cluster_id) {
  df <- top20[[cluster_id]]
  data.frame(cluster_id = as.integer(cluster_id), gene = df$gene,
             log2FoldChange = df$log2FoldChange, foldChange = df$foldChange,
             pvalue = df$pvalue, FDR = df$padj, baseMean = df$baseMean,
             row.names = NULL)
}))
write.csv(top_table, file.path(out, "tables", "top20_up_markers_by_cluster.csv"),
          row.names = FALSE)

# Volcano: cinco contrastes en una sola imagen.
volcano_df <- do.call(rbind, all_deg)
volcano_df$cluster_label <- factor(paste("Cluster", volcano_df$cluster_id),
                                   levels = paste("Cluster", 1:5))
volcano_df$minus_log10_fdr <- -log10(pmax(volcano_df$padj, .Machine$double.xmin))
label_df <- do.call(rbind, lapply(top20, head, 10))
label_df$cluster_label <- factor(paste("Cluster", label_df$cluster_id),
                                 levels = paste("Cluster", 1:5))
label_df$minus_log10_fdr <- -log10(pmax(label_df$padj, .Machine$double.xmin))
p_volcano <- ggplot(volcano_df,
                    aes(log2FoldChange, minus_log10_fdr, colour = status)) +
  geom_point(alpha = 0.55, size = 0.55) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
  geom_text_repel(data = label_df, aes(label = gene), colour = "black",
                  size = 2.2, max.overlaps = Inf, min.segment.length = 0) +
  facet_wrap(~cluster_label, ncol = 2, scales = "free") +
  scale_colour_manual(values = c("Down" = "#377EB8",
                                 "Not significant" = "grey78", "Up" = "#E41A1C")) +
  labs(title = "Differential expression: each cluster versus the rest",
       x = "log2 fold change", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom",
                                   strip.text = element_text(face = "bold"))
ggsave(file.path(out, "figures", "volcano_all_clusters.png"),
       p_volcano, width = 13, height = 14, dpi = 300, bg = "white")

# Reconstrucción exacta del preprocesamiento previo para obtener la matriz VST.
load(file.path(project_root, "data", "strap_counts.RData"))
stopifnot(exists("strap_counts"), ncol(strap_counts) == 210L,
          !anyDuplicated(colnames(strap_counts)),
          setequal(colnames(strap_counts), assignments$sample_id))
meta <- data.frame(pathotype_clean = assignments$pathotype_histological,
                   row.names = assignments$sample_id)
expr_scaled <- preprocess_expression(strap_counts, meta = meta,
                                     var_quantile_cutoff = 0.25,
                                     n_top = 1000, kw_pval = 0.01,
                                     blind = TRUE, seed = 42)
expr_scaled <- expr_scaled[assignments$sample_id, , drop = FALSE]

# Heatmap de resultados: diez marcadores positivos por cluster para legibilidad.
# Los DEGs se calcularon sobre todos los genes expresados, no solo sobre las
# variables usadas como entrada del clustering; por ello se genera una VST
# descriptiva con el mismo filtro empleado por DESeq2 en los contrastes.
heatmap_genes <- unique(unlist(lapply(top20, function(df) head(df$gene, 10))))
counts_for_deg <- strap_counts[rowSums(strap_counts >= 5) >= 2, , drop = FALSE]
dds_heatmap <- DESeqDataSetFromMatrix(
  countData = round(counts_for_deg),
  colData = data.frame(condition = factor(rep("all", ncol(counts_for_deg))),
                       row.names = colnames(counts_for_deg)),
  design = ~1
)
vsd_heatmap <- vst(dds_heatmap, blind = TRUE)
expr_vst_all <- t(assay(vsd_heatmap))
heatmap_genes <- intersect(heatmap_genes, colnames(expr_vst_all))
stopifnot(length(heatmap_genes) >= 45)
sample_order <- order(assignments$cluster_id, assignments$sample_id)
mat <- t(expr_vst_all[assignments$sample_id[sample_order], heatmap_genes, drop = FALSE])
mat <- t(scale(t(mat)))
mat[!is.finite(mat)] <- 0
ann_col <- data.frame(Cluster = assignments$cluster_label[sample_order],
                      row.names = assignments$sample_id[sample_order])
gaps <- cumsum(as.integer(table(ann_col$Cluster)))
gaps <- gaps[-length(gaps)]
heatmap_plot <- pheatmap(
  mat, annotation_col = ann_col,
  annotation_colors = list(Cluster = palette),
  show_rownames = TRUE, show_colnames = FALSE,
  cluster_rows = TRUE, cluster_cols = FALSE, gaps_col = gaps,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main = "Top differentially expressed genes by cluster",
  fontsize = 17, fontsize_row = 16, fontsize_main = 19,
  silent = TRUE
)
heatmap_plot$gtable <- gtable::gtable_add_cols(
  heatmap_plot$gtable, grid::unit(2.5, "cm")
)
png(file.path(out, "figures", "heatmap_final_top_markers.png"),
    width = 16, height = 14, units = "in", res = 300, bg = "white")
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)
dev.off()

cluster_sizes <- as.data.frame(table(cluster_id = assignments$cluster_id))
cluster_sizes$cluster_id <- as.integer(as.character(cluster_sizes$cluster_id))
names(cluster_sizes)[2] <- "n_samples"
write.csv(cluster_sizes, file.path(out, "tables", "cluster_sizes.csv"),
          row.names = FALSE)
write.csv(assignments[, c("sample_id", "cluster_id")],
          file.path(out, "tables", "sample_assignments_hierarchical_k5.csv"),
          row.names = FALSE)

audit <- data.frame(
  check = c("total_samples", "unique_samples", "clusters", "cluster_size_sum",
            "top_markers", "all_markers_pass_FDR", "all_markers_pass_log2FC"),
  value = c(nrow(assignments), length(unique(assignments$sample_id)),
            length(unique(assignments$cluster_id)), sum(cluster_sizes$n_samples),
            nrow(top_table), all(top_table$FDR < fdr_threshold),
            all(top_table$log2FoldChange >= lfc_threshold))
)
write.csv(audit, file.path(out, "tables", "audit_checks.csv"), row.names = FALSE)
saveRDS(expr_scaled, file.path(out, "tables", "expression_matrix_scaled.rds"))
cat("Final numeric results created in:", out, "\n")
