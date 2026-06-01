# ==============================================================================
# tutor_followup_analysis.R
# Análisis de seguimiento — Revisión del tutor (junio 2026)
#
# Tareas:
#   1. Tabla de marcadores top por cluster (DEGs UP)
#   2. Enrichment analysis del cluster Myeloid (enrichR)
#   3. Por qué las muestras contaminadas van al cluster Myeloid
#   4. PCA con marcadores canónicos y propios (feature plots)
# ==============================================================================

setwd("C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes")

# ── Paquetes ──────────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("enrichR",     quietly = TRUE)) install.packages("enrichR")
if (!requireNamespace("ggplot2",     quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",       quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("patchwork",   quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("pheatmap",    quietly = TRUE)) install.packages("pheatmap")

library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)
library(enrichR)

# ── Parámetros ────────────────────────────────────────────────────────────────
RESULTS_DIR <- "results_k5_v3/results"
DATA_FILE   <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes/strap_counts.RData"
META_FILE   <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes/E-MTAB-13733.sdrf.txt"

PATHOTYPE_COLORS <- c(
  "Fibroid"  = "#E41A1C",
  "Vascular" = "#FF7F00",
  "Lymphoid" = "#4DAF4A",
  "Myeloid"  = "#377EB8",
  "IFN-high" = "#984EA3"
)

N_TOP <- 15  # marcadores por cluster a mostrar


# ==============================================================================
# CARGAR DATOS Y RESULTADOS PREVIOS
# ==============================================================================

cat("\n=== Cargando datos y resultados previos ===\n")

# Cargar datos de expresión
source("functions/00_setup.R")
source("functions/01_preprocessing.R")
source("functions/05_annotation.R")

setup_environment()

expr_raw <- load_counts(DATA_FILE)

META_COLS <- c(
  "Characteristics.sampleid.",
  "Characteristics.pathotype.",
  "Characteristics.tender_joints_counts_tjc.",
  "Characteristics.swollen_joints_counts_sjc.",
  "Characteristics.arthritis.activity.",
  "Characteristics.esr.",
  "Characteristics.crp.",
  "Characteristics.physicians.global.assessment."
)

meta <- load_metadata(META_FILE, META_COLS)

expr_vst    <- preprocess_expression(expr_raw, n_top = 1000, blind = TRUE, seed = 42)
aligned     <- align_samples(expr_vst, meta)
expr_scaled <- aligned$expr
meta_clean  <- aligned$meta

# Cargar asignaciones de cluster
clusters_df <- read.csv(file.path(RESULTS_DIR, "clusters/cluster_assignments.csv"),
                        row.names = 1)
clusters    <- setNames(clusters_df$cluster_id, clusters_df$sample_id)
cluster_bio <- setNames(clusters_df$cluster_bio, clusters_df$sample_id)

# Alinear muestras
common_samples <- intersect(rownames(expr_scaled), names(clusters))
expr_scaled    <- expr_scaled[common_samples, ]
clusters       <- clusters[common_samples]
cluster_bio    <- cluster_bio[common_samples]
meta_clean     <- meta_clean[meta_clean$sample_id %in% common_samples, ]
rownames(meta_clean) <- meta_clean$sample_id

cat(sprintf("  Muestras alineadas: %d\n", length(common_samples)))

# Cargar DEGs
deg_files <- list(
  Fibroid  = read.csv(file.path(RESULTS_DIR, "DEGs/DEGs_UP_Fibroid.csv"),  row.names = 1),
  Vascular = read.csv(file.path(RESULTS_DIR, "DEGs/DEGs_UP_Vascular.csv"), row.names = 1),
  Lymphoid = read.csv(file.path(RESULTS_DIR, "DEGs/DEGs_UP_Lymphoid.csv"), row.names = 1),
  Myeloid  = read.csv(file.path(RESULTS_DIR, "DEGs/DEGs_UP_Myeloid.csv"),  row.names = 1),
  `IFN-high` = read.csv(file.path(RESULTS_DIR, "DEGs/DEGs_UP_IFN_high.csv"), row.names = 1)
)


# ==============================================================================
# TAREA 1: TABLA DE MARCADORES TOP POR CLUSTER
# ==============================================================================

cat("\n=== [1/4] Tabla de marcadores top por cluster ===\n")

marker_table <- do.call(rbind, lapply(names(deg_files), function(cl) {
  df <- deg_files[[cl]]
  df <- df[order(df$padj), ]
  top <- head(df, N_TOP)
  data.frame(
    cluster      = cl,
    gene         = top$gene,
    log2FC       = round(top$log2FoldChange, 3),
    FDR          = formatC(top$padj, format = "e", digits = 2),
    baseMean     = round(top$baseMean, 1),
    stringsAsFactors = FALSE
  )
}))

write.csv(marker_table,
          file.path(RESULTS_DIR, "tables/top_markers_by_cluster.csv"),
          row.names = FALSE)
cat(sprintf("  → Tabla guardada: %d genes × %d clusters\n",
            N_TOP, length(deg_files)))

# Heatmap de top marcadores
top_genes <- unique(marker_table$gene)
genes_in_expr <- intersect(top_genes, colnames(expr_scaled))

mat_markers <- t(expr_scaled[, genes_in_expr])

# Ordenar muestras por cluster
sample_order <- order(cluster_bio[common_samples])
mat_markers  <- mat_markers[, sample_order]

anno_col <- data.frame(
  Cluster = cluster_bio[colnames(mat_markers)],
  row.names = colnames(mat_markers)
)

ann_colors <- list(Cluster = PATHOTYPE_COLORS)

png(file.path(RESULTS_DIR, "figures/heatmap_top_markers.png"),
    width = 3600, height = 3000, res = 300)
pheatmap(mat_markers,
         annotation_col  = anno_col,
         annotation_colors = ann_colors,
         cluster_cols    = FALSE,
         cluster_rows    = TRUE,
         show_colnames   = FALSE,
         fontsize_row    = 7,
         color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main            = "Top marcadores DEG por cluster",
         border_color    = NA)
dev.off()
cat("  → Heatmap guardado: heatmap_top_markers.png\n")


# ==============================================================================
# TAREA 2: ENRICHMENT ANALYSIS DEL CLUSTER MYELOID
# ==============================================================================

cat("\n=== [2/4] Enrichment analysis — cluster Myeloid ===\n")

myeloid_genes <- deg_files$Myeloid$gene

# Bases de datos a consultar
dbs <- c("GO_Biological_Process_2023",
         "KEGG_2021_Human",
         "MSigDB_Hallmark_2020")

cat("  Consultando enrichR (requiere conexión a internet)...\n")
enrich_results <- tryCatch({
  enrichr(myeloid_genes, dbs)
}, error = function(e) {
  warning("enrichR falló: ", e$message,
          "\nVerifica conexión a internet o ejecuta manualmente.")
  NULL
})

if (!is.null(enrich_results)) {

  # Guardar tablas completas
  for (db in dbs) {
    db_clean <- gsub("[^a-zA-Z0-9]", "_", db)
    res_df   <- enrich_results[[db]]
    res_df   <- res_df[order(res_df$Adjusted.P.value), ]
    write.csv(res_df,
              file.path(RESULTS_DIR, paste0("enrichment/enrichR_Myeloid_", db_clean, ".csv")),
              row.names = FALSE)
  }

  # Plot top 15 términos por base de datos
  plots_enrich <- lapply(dbs, function(db) {
    res_df <- enrich_results[[db]]
    res_df <- res_df[res_df$Adjusted.P.value < 0.05, ]
    if (nrow(res_df) == 0) return(NULL)
    res_df <- head(res_df[order(res_df$Adjusted.P.value), ], 15)
    res_df$Term <- factor(res_df$Term,
                          levels = rev(res_df$Term[order(res_df$Adjusted.P.value)]))
    res_df$neg_log10_padj <- -log10(res_df$Adjusted.P.value)

    ggplot(res_df, aes(x = neg_log10_padj, y = Term, fill = neg_log10_padj)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
      labs(title = gsub("_", " ", db),
           x = "-log10(FDR)", y = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 10))
  })

  plots_enrich <- Filter(Negate(is.null), plots_enrich)

  if (length(plots_enrich) > 0) {
    combined_enrich <- wrap_plots(plots_enrich, ncol = 1)
    ggsave(file.path(RESULTS_DIR, "figures/enrichR_Myeloid.png"),
           combined_enrich, width = 10, height = 4 * length(plots_enrich),
           dpi = 300, bg = "white")
    cat("  → Figura guardada: enrichR_Myeloid.png\n")
  }

} else {
  cat("  ⚠ enrichR no disponible. Tablas no generadas.\n")
}


# ==============================================================================
# TAREA 3: ¿POR QUÉ LOS CONTAMINADOS VAN AL CLUSTER MYELOID?
# ==============================================================================

cat("\n=== [3/4] Investigación cluster Myeloid (artefacto) ===\n")

# Identificar muestras Myeloid
myeloid_samples <- names(cluster_bio)[cluster_bio == "Myeloid"]
cat(sprintf("  Muestras Myeloid: %s\n", paste(myeloid_samples, collapse = ", ")))

# ── 3a. Scores de firma génica para cada muestra ─────────────────────────────
RA_SIGS <- list(
  Fibroid  = c("CXCL12", "PDPN", "PRG4", "THY1", "CDH11", "FAP",
               "ACTA2", "COL1A1", "COL3A1", "FN1"),
  Myeloid  = c("CD68", "MRC1", "CSF1R", "IL1B", "TNF", "IL6",
               "ITGAM", "CD14", "FCGR3A", "CCL2"),
  Lymphoid = c("CD3D", "CD3E", "CD19", "MS4A1", "SELL", "PTPRC",
               "CD4", "CD8A", "FOXP3", "CXCR5"),
  IFN_high = c("MX1", "OAS1", "IFIT1", "IFIT3", "ISG15", "RSAD2",
               "IFI44L", "CXCL10", "IRF7", "STAT1"),
  Vascular = c("EMCN", "JAM2", "AQP1", "SPARCL1", "PKN3",
               "TRPC1", "TNFRSF11B", "DIO2", "CLEC3A", "SNTB2"),
  Muscle   = c("TNNI1", "NEB", "TNNI2", "MYH1", "MYH2", "ACTN2",
               "TNNC2", "MYL1", "TPM1", "TTN")
)

score_sig <- function(expr, genes) {
  avail <- intersect(genes, colnames(expr))
  if (length(avail) == 0) return(rep(0, nrow(expr)))
  rowMeans(expr[, avail, drop = FALSE])
}

sig_scores <- as.data.frame(sapply(RA_SIGS, score_sig, expr = expr_scaled))
sig_scores$sample_id   <- rownames(sig_scores)
sig_scores$cluster_bio <- cluster_bio[rownames(sig_scores)]

# Guardar scores
write.csv(sig_scores,
          file.path(RESULTS_DIR, "tables/signature_scores_all_samples.csv"),
          row.names = FALSE)

# ── 3b. Comparar scores entre clusters ───────────────────────────────────────
sig_long <- tidyr::pivot_longer(sig_scores,
                                 cols = names(RA_SIGS),
                                 names_to = "signature",
                                 values_to = "score")

p_scores <- ggplot(sig_long, aes(x = cluster_bio, y = score,
                                  fill = cluster_bio)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~signature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(PATHOTYPE_COLORS, Muscle = "#8B4513")) +
  labs(title = "Scores de firma génica por cluster",
       subtitle = "Incluye firma muscular para detectar contaminación",
       x = NULL, y = "Score medio (z-score)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

ggsave(file.path(RESULTS_DIR, "figures/signature_scores_by_cluster.png"),
       p_scores, width = 14, height = 10, dpi = 300, bg = "white")
cat("  → Figura guardada: signature_scores_by_cluster.png\n")

# ── 3c. Score medio por cluster (tabla resumen) ───────────────────────────────
score_summary <- sig_scores %>%
  group_by(cluster_bio) %>%
  summarise(across(all_of(names(RA_SIGS)), mean, .names = "{.col}")) %>%
  arrange(match(cluster_bio, c("Fibroid","Vascular","Lymphoid","Myeloid","IFN-high")))

write.csv(score_summary,
          file.path(RESULTS_DIR, "tables/signature_score_means_by_cluster.csv"),
          row.names = FALSE)

cat("\n  Score medio de CADA FIRMA por cluster:\n")
print(as.data.frame(score_summary), digits = 3)

# ── 3d. Radar/heatmap de scores medios ───────────────────────────────────────
score_mat <- as.matrix(score_summary[, names(RA_SIGS)])
rownames(score_mat) <- score_summary$cluster_bio

# Normalizar por firma para visualización
score_mat_norm <- scale(score_mat)

png(file.path(RESULTS_DIR, "figures/signature_heatmap_clusters.png"),
    width = 2400, height = 1600, res = 300)
pheatmap(score_mat_norm,
         color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_rows    = FALSE,
         cluster_cols    = FALSE,
         display_numbers = round(score_mat, 3),
         number_format   = "%.3f",
         fontsize         = 11,
         main            = "Score medio de firma génica por cluster (incluyendo Muscular)",
         border_color    = NA)
dev.off()
cat("  → Heatmap guardado: signature_heatmap_clusters.png\n")

# ── 3e. Expresión de genes musculares en muestras Myeloid ────────────────────
muscle_genes <- c("TNNI1", "NEB", "TNNI2", "MYH1", "MYH2", "ACTN2", "TNNC2")
muscle_avail <- intersect(muscle_genes, colnames(expr_scaled))

if (length(muscle_avail) > 0) {
  muscle_expr <- expr_scaled[, muscle_avail, drop = FALSE]
  muscle_df   <- as.data.frame(muscle_expr)
  muscle_df$sample_id   <- rownames(muscle_df)
  muscle_df$cluster_bio <- cluster_bio[rownames(muscle_df)]

  muscle_long <- tidyr::pivot_longer(muscle_df,
                                      cols = all_of(muscle_avail),
                                      names_to = "gene",
                                      values_to = "expression")

  p_muscle <- ggplot(muscle_long, aes(x = cluster_bio, y = expression,
                                       fill = cluster_bio)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
    facet_wrap(~gene, scales = "free_y") +
    scale_fill_manual(values = PATHOTYPE_COLORS) +
    labs(title = "Expresión de genes musculares por cluster",
         subtitle = "El cluster Myeloid muestra expresión extrema de marcadores musculares",
         x = NULL, y = "Expresión escalada (z-score)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(RESULTS_DIR, "figures/muscle_genes_by_cluster.png"),
         p_muscle, width = 14, height = 8, dpi = 300, bg = "white")
  cat("  → Figura guardada: muscle_genes_by_cluster.png\n")
}

# ── 3f. Metadatos clínicos del cluster Myeloid ───────────────────────────────
if ("crp" %in% colnames(meta_clean)) {
  meta_clean$crp <- suppressWarnings(as.numeric(meta_clean$crp))
  meta_clean$cluster_bio <- cluster_bio[meta_clean$sample_id]

  meta_myeloid <- meta_clean[meta_clean$cluster_bio == "Myeloid", ]
  cat("\n  Metadatos clínicos del cluster Myeloid (artefacto):\n")
  print(meta_myeloid[, intersect(c("sample_id","pathotype_clean","crp",
                                    "sex","cluster_bio"), colnames(meta_myeloid))])
}


# ==============================================================================
# TAREA 4: PCA CON MARCADORES CANÓNICOS Y PROPIOS (FEATURE PLOTS)
# ==============================================================================

cat("\n=== [4/4] PCA con marcadores canónicos — Feature plots ===\n")

# Calcular PCA
pca_res   <- prcomp(expr_scaled, center = FALSE, scale. = FALSE)
pca_df    <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample_id   <- rownames(pca_df)
pca_df$cluster_bio <- cluster_bio[rownames(pca_df)]

var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

# ── 4a. PCA base coloreado por cluster ───────────────────────────────────────
p_base <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_bio)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = PATHOTYPE_COLORS) +
  labs(title = "PCA — Clusters K-means (k=5)",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       color = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# ── 4b. Feature plots — marcadores canónicos ─────────────────────────────────
# Marcadores canónicos por patotipo
canonical_markers <- list(
  Fibroid  = c("FAP", "PDPN", "CDH11", "COL1A1"),
  Lymphoid = c("CD3D", "CD19", "MS4A1", "FOXP3"),
  `IFN-high` = c("MX1", "ISG15", "IFIT1", "OAS1"),
  Vascular = c("EMCN", "AQP1", "JAM2", "SPARCL1"),
  Myeloid  = c("CD68", "CSF1R", "TNNI1", "NEB")
)

# Marcadores propios (top 3 DEG por cluster)
own_markers <- lapply(deg_files, function(df) {
  head(df$gene[order(df$padj)], 3)
})

all_markers <- unique(c(
  unlist(canonical_markers),
  unlist(own_markers)
))
all_markers <- intersect(all_markers, colnames(expr_scaled))

feature_plot <- function(gene, pca_df, expr_scaled) {
  if (!gene %in% colnames(expr_scaled)) return(NULL)

  expr_gene <- expr_scaled[rownames(pca_df), gene]

  df_plot <- pca_df
  df_plot$expr <- expr_gene

  ggplot(df_plot, aes(x = PC1, y = PC2, color = expr)) +
    geom_point(size = 1.8, alpha = 0.85) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0) +
    labs(title = gene, x = NULL, y = NULL, color = "z-score") +
    theme_minimal(base_size = 9) +
    theme(plot.title  = element_text(face = "bold", size = 10),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 7))
}

# Feature plots por grupo de marcadores
for (group_name in names(canonical_markers)) {
  genes  <- intersect(canonical_markers[[group_name]], colnames(expr_scaled))
  if (length(genes) == 0) next

  plots  <- lapply(genes, feature_plot, pca_df = pca_df, expr_scaled = expr_scaled)
  plots  <- Filter(Negate(is.null), plots)

  combined <- wrap_plots(c(list(p_base + theme(legend.position = "none")),
                            plots),
                          ncol = min(3, length(plots) + 1))

  fname <- gsub("[^a-zA-Z0-9]", "_", group_name)
  ggsave(file.path(RESULTS_DIR, sprintf("figures/pca_markers_%s.png", fname)),
         combined, width = 14, height = 5, dpi = 300, bg = "white")
  cat(sprintf("  → Figura guardada: pca_markers_%s.png\n", fname))
}

# ── 4c. Feature plots — top 3 marcadores propios por cluster ─────────────────
for (cl in names(own_markers)) {
  genes  <- intersect(own_markers[[cl]], colnames(expr_scaled))
  if (length(genes) == 0) next

  plots  <- lapply(genes, feature_plot, pca_df = pca_df, expr_scaled = expr_scaled)
  plots  <- Filter(Negate(is.null), plots)

  combined <- wrap_plots(c(list(p_base + theme(legend.position = "none")),
                            plots),
                          ncol = min(4, length(plots) + 1))

  fname <- gsub("[^a-zA-Z0-9]", "_", cl)
  ggsave(file.path(RESULTS_DIR, sprintf("figures/pca_own_markers_%s.png", fname)),
         combined, width = 14, height = 5, dpi = 300, bg = "white")
  cat(sprintf("  → Figura guardada: pca_own_markers_%s.png\n", fname))
}


# ==============================================================================
# RESUMEN FINAL
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("ANÁLISIS COMPLETADO\n")
cat(strrep("=", 60), "\n")
cat("Archivos generados en", file.path(RESULTS_DIR), ":\n\n")
cat("  TABLAS:\n")
cat("    tables/top_markers_by_cluster.csv\n")
cat("    tables/signature_scores_all_samples.csv\n")
cat("    tables/signature_score_means_by_cluster.csv\n")
if (!is.null(enrich_results)) {
  cat("    enrichment/enrichR_Myeloid_*.csv  (3 bases de datos)\n")
}
cat("\n  FIGURAS:\n")
cat("    figures/heatmap_top_markers.png\n")
if (!is.null(enrich_results)) cat("    figures/enrichR_Myeloid.png\n")
cat("    figures/signature_scores_by_cluster.png\n")
cat("    figures/signature_heatmap_clusters.png\n")
cat("    figures/muscle_genes_by_cluster.png\n")
cat("    figures/pca_markers_[Fibroid|Lymphoid|IFN_high|Vascular|Myeloid].png\n")
cat("    figures/pca_own_markers_[cluster].png\n")
cat(strrep("=", 60), "\n")
