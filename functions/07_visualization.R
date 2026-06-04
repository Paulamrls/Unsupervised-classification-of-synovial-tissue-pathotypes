# ==============================================================================
# functions/06_visualization.R
# Todas las visualizaciones — colores consistentes en TODOS los gráficos
# ==============================================================================

# ── Tema ggplot2 corporativo ──────────────────────────────────────────────────
theme_strap <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(colour = "grey40"),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold")
    )
}

# ── Helpers internos ──────────────────────────────────────────────────────────

# Convierte cluster_ids a factor con nombres biológicos.
# Clusters sin nombre (métodos que generan más de k clusters) → "Unresolved"
cl_to_factor <- function(cl, cluster_labels) {
  bio <- cluster_labels[as.character(cl)]
  bio[is.na(bio)] <- "Unresolved"
  lvls <- c(unique(cluster_labels[order(as.numeric(names(cluster_labels)))]),
            "Unresolved")
  factor(bio, levels = lvls)
}

# Añade columna de label al data.frame de coordenadas
add_cluster_col <- function(coords_df, cl, cluster_labels) {
  coords_df$cluster_id  <- cl[rownames(coords_df)]
  coords_df$cluster_bio <- cl_to_factor(coords_df$cluster_id, cluster_labels)
  coords_df
}

# ── PCA coloreado por patotipo biológico ──────────────────────────────────────
plot_pca_clusters <- function(pca_res, clusters, cluster_labels,
                               cluster_colors, var1 = 1, var2 = 2,
                               title = "PCA — Clusters") {
  df <- as.data.frame(pca_res$coords[, c(var1, var2)])
  colnames(df) <- c("PC_x", "PC_y")
  df <- add_cluster_col(df, clusters, cluster_labels)

  pct <- round(pca_res$var_explained[c(var1, var2)], 1)

  ggplot(df, aes(x = PC_x, y = PC_y, colour = cluster_bio)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_manual(values = cluster_colors, name = "Patotipo") +
    labs(
      title    = title,
      x        = paste0("PC", var1, " (", pct[1], "%)"),
      y        = paste0("PC", var2, " (", pct[2], "%)")
    ) +
    theme_strap()
}

# ── PCA coloreado por patotipo histológico (referencia externa) ───────────────
plot_pca_histological <- function(pca_res, meta, pathotype_col,
                                   pathotype_order, pathotype_colors,
                                   var1 = 1, var2 = 2) {
  common <- intersect(rownames(pca_res$coords), rownames(meta))
  df     <- as.data.frame(pca_res$coords[common, c(var1, var2)])
  colnames(df) <- c("PC_x", "PC_y")
  df$pathotype <- factor(meta[common, pathotype_col],
                          levels = pathotype_order)

  pct <- round(pca_res$var_explained[c(var1, var2)], 1)

  # Usar los colores de pathotype_order
  pt_cols <- pathotype_colors[pathotype_order]
  pt_cols <- pt_cols[!is.na(pt_cols)]

  ggplot(df, aes(x = PC_x, y = PC_y, colour = pathotype)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_manual(values = pt_cols, name = "Patotipo\nhistológico") +
    labs(
      title = "PCA — Patotipos histológicos (referencia)",
      x     = paste0("PC", var1, " (", pct[1], "%)"),
      y     = paste0("PC", var2, " (", pct[2], "%)")
    ) +
    theme_strap()
}

# ── Diffusion Maps coloreado ──────────────────────────────────────────────────
# NOTA: recorta los ejes al percentil 2.5–97.5 para evitar que outliers
# colapsen la visualización (problema frecuente con diffusionMap).
plot_dm_clusters <- function(dm_res, clusters, cluster_labels,
                              cluster_colors, dc1 = 1, dc2 = 2,
                              title = "Diffusion Map — Clusters") {
  df <- as.data.frame(dm_res$coords[, c(dc1, dc2)])
  colnames(df) <- c("DC_x", "DC_y")
  df <- add_cluster_col(df, clusters, cluster_labels)

  # Límites robustos (percentil 2.5–97.5) para excluir outliers extremos
  xlim <- quantile(df$DC_x, probs = c(0.025, 0.975), na.rm = TRUE)
  ylim <- quantile(df$DC_y, probs = c(0.025, 0.975), na.rm = TRUE)

  # Paleta extendida con "Unresolved"
  cols_ext <- c(cluster_colors, "Unresolved" = "#999999")

  ggplot(df, aes(x = DC_x, y = DC_y, colour = cluster_bio)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_manual(values = cols_ext, name = "Patotipo") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title,
         subtitle = "Ejes recortados al percentil 2.5–97.5",
         x = paste0("DC", dc1), y = paste0("DC", dc2)) +
    theme_strap()
}

# ── Comparación PCA vs Diffusion Maps (panel) ─────────────────────────────────
plot_dimred_comparison <- function(pca_res, dm_res, clusters,
                                    cluster_labels, cluster_colors) {
  p1 <- plot_pca_clusters(pca_res, clusters, cluster_labels,
                           cluster_colors, title = "PCA")
  p2 <- plot_dm_clusters(dm_res,  clusters, cluster_labels,
                          cluster_colors, title = "Diffusion Map")
  p1 + p2 + plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# ── Silhouette plot ───────────────────────────────────────────────────────────
plot_silhouette_comparison <- function(validation, cluster_colors) {
  df <- validation$ranking
  # Mantener solo filas con silhouette no-NA
  df <- df[!is.na(df$silhouette), ]
  df$method <- factor(df$method, levels = df$method)   # orden por ranking

  ggplot(df, aes(x = method, y = silhouette)) +
    geom_col(fill = "#2166AC", width = 0.6) +
    geom_errorbar(aes(ymin = silhouette - 0.02,
                      ymax = silhouette + 0.02), width = 0.2) +
    geom_text(aes(label = round(silhouette, 3)),
              vjust = -0.5, size = 3.5) +
    labs(title = "Silhouette score por método de clustering",
         x = "Método", y = "Silhouette medio") +
    theme_strap() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ── Bootstrap stability ───────────────────────────────────────────────────────
plot_bootstrap_stability <- function(validation) {
  df <- validation$ranking
  df <- df[!is.na(df$boot_ari), ]
  df$method <- factor(df$method, levels = df$method)

  ggplot(df, aes(x = method, y = boot_ari)) +
    geom_col(fill = "#4DAF4A", width = 0.6) +
    geom_errorbar(aes(ymin = boot_ari - boot_sd,
                      ymax = boot_ari + boot_sd), width = 0.25) +
    geom_text(aes(label = round(boot_ari, 3)), vjust = -0.5, size = 3.5) +
    labs(title = "Estabilidad bootstrap (ARI medio ± SD)",
         x = "Método", y = "ARI medio") +
    theme_strap() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ── Heatmap de expresión ──────────────────────────────────────────────────────
plot_heatmap <- function(expr, clusters, cluster_labels, cluster_colors,
                          n_genes = 50, base_dir = ".") {
  # Seleccionar genes más variables para el heatmap
  var_g     <- apply(expr, 2, var)
  top_genes <- names(sort(var_g, decreasing = TRUE))[seq_len(min(n_genes, ncol(expr)))]

  mat <- t(expr[, top_genes])   # genes × muestras (para pheatmap)

  # Anotación de muestras
  ann_col <- data.frame(
    Patotipo = cluster_labels[as.character(clusters)],
    row.names = names(clusters)
  )

  ann_colors <- list(Patotipo = cluster_colors)

  path <- file.path(base_dir, "results", "figures", "heatmap_clusters.png")
  pheatmap(mat,
           annotation_col  = ann_col,
           annotation_colors = ann_colors,
           show_rownames   = FALSE,
           show_colnames   = FALSE,
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           color           = colorRampPalette(rev(RColorBrewer::brewer.pal(
             11, "RdBu")))(100),
           main            = "Heatmap — Top genes variables por patotipo",
           filename        = path,
           width           = 12, height = 10)
  cat("  Figura guardada: heatmap_clusters.png\n")
}

# ── Tabla de contingencia normalizada (%) ─────────────────────────────────────
plot_confusion_pct <- function(validation, best_method, cluster_labels,
                                base_dir = ".") {
  bio   <- validation$bio_overlap[[best_method]]
  ct    <- bio$table_pct
  if (is.null(ct)) return(invisible(NULL))

  # Convertir a data.frame largo
  df <- as.data.frame(ct)
  colnames(df) <- c("Cluster_id", "Patotipo_histologico", "Pct")

  # Sustituir IDs numéricos por nombres biológicos en el eje Y
  df$Cluster <- cluster_labels[as.character(df$Cluster_id)]
  df$Cluster[is.na(df$Cluster)] <- paste0("Cluster_", df$Cluster_id[is.na(df$Cluster)])

  # Orden canónico en eje Y
  df$Cluster <- factor(df$Cluster,
                        levels = rev(unique(cluster_labels)))

  path <- file.path(base_dir, "results", "figures", "confusion_pct.png")
  png(path, width = 2400, height = 1800, res = 300, bg = "white")
  p <- ggplot(df, aes(x = Patotipo_histologico, y = Cluster, fill = Pct)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(round(Pct, 1), "%")), size = 4) +
    scale_fill_gradient(low = "white", high = "#2166AC",
                        name = "% muestras\npor cluster") +
    labs(title = paste("Clusters vs Patotipos histológicos (%) —", best_method),
         x = "Patotipo histológico", y = "Cluster") +
    theme_strap()
  print(p)
  dev.off()
  cat("  Figura guardada: confusion_pct.png\n")
}

# ── Sankey — clusters → patotipos histológicos ────────────────────────────────
# Usa ggalluvial para mantener la paleta de colores consistente
plot_sankey <- function(clusters, meta, cluster_labels, cluster_colors,
                         pathotype_col = "pathotype_clean",
                         pathotype_order, base_dir = ".") {
  common <- intersect(names(clusters), rownames(meta))
  df <- data.frame(
    cluster_bio = cluster_labels[as.character(clusters[common])],
    pathotype   = meta[common, pathotype_col],
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$pathotype), ]

  # Orden canónico en ejes
  df$cluster_bio <- factor(df$cluster_bio, levels = names(cluster_colors))
  df$pathotype   <- factor(df$pathotype,   levels = pathotype_order)

  # Frecuencias para ggalluvial
  df_freq <- df %>%
    dplyr::count(cluster_bio, pathotype, name = "Freq") %>%
    dplyr::filter(Freq > 0)

  # Paleta unificada: cubre nombres de clusters Y nombres de patotipos histológicos.
  # Como "Fibroid","Myeloid","Lymphoid" aparecen en ambos ejes, reciben el mismo
  # color → consistencia visual total en ambos lados del Sankey.
  all_labels    <- unique(c(levels(df_freq$cluster_bio),
                             levels(df_freq$pathotype)))
  unified_colors <- cluster_colors[all_labels]
  # Para etiquetas sin color asignado (raro), usar gris
  unified_colors[is.na(unified_colors)] <- "#999999"

  p <- ggplot(df_freq,
              aes(axis1 = cluster_bio, axis2 = pathotype, y = Freq)) +
    geom_alluvium(aes(fill = cluster_bio), alpha = 0.75, width = 1/4) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/4,
                 colour = "white", linewidth = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = unified_colors, guide = "none") +
    scale_x_discrete(limits = c("Cluster predicho", "Patotipo histológico"),
                     expand = c(0.05, 0.05)) +
    labs(title = "Sankey — Clusters predichos vs Patotipos histológicos",
         y = "Número de muestras") +
    theme_strap() +
    theme(axis.text.x = element_text(face = "bold", size = 12))

  path <- file.path(base_dir, "results", "figures", "sankey.png")
  ggsave(path, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  cat("  Figura guardada: sankey.png\n")
}

# ── Boxplots de genes marcadores ──────────────────────────────────────────────
plot_marker_boxplots <- function(deg_results, expr, clusters, cluster_labels,
                                  cluster_colors, n_genes = 10,
                                  base_dir = ".") {
  bio_labels <- cluster_labels[as.character(clusters)]

  for (cl_name in names(deg_results$up)) {
    # Top n genes up y down
    up_genes   <- head(rownames(deg_results$up[[cl_name]]),   n_genes)
    down_genes <- head(rownames(deg_results$down[[cl_name]]), n_genes)
    genes      <- c(up_genes, down_genes)
    genes      <- intersect(genes, colnames(expr))
    if (length(genes) == 0) next

    # Long format
    df_long <- reshape2::melt(
      cbind(cluster = bio_labels,
            as.data.frame(expr[, genes, drop = FALSE])),
      id.vars      = "cluster",
      variable.name = "gene",
      value.name    = "expression"
    )
    df_long$direction <- ifelse(df_long$gene %in% up_genes,
                                "Upregulated", "Downregulated")
    df_long$gene      <- factor(df_long$gene,
                                levels = c(up_genes, down_genes))

    p <- ggplot(df_long, aes(x = cluster, y = expression, fill = cluster)) +
      geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
      geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
      facet_wrap(~ gene, scales = "free_y", ncol = 5) +
      scale_fill_manual(values = cluster_colors, guide = "none") +
      labs(title = paste("Genes marcadores —", cl_name),
           x = NULL, y = "Expresión (VST escalada)") +
      theme_strap() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
            strip.text  = element_text(size = 7))

    fname <- paste0("boxplots_markers_",
                    gsub("[^A-Za-z0-9]", "_", cl_name), ".png")
    path  <- file.path(base_dir, "results", "figures", fname)
    ggsave(path, plot = p,
           width  = max(10, length(genes) * 1.5),
           height = ceiling(length(genes) / 5) * 3 + 2,
           dpi    = 300, bg = "white",
           limitsize = FALSE)
    cat("  Figura guardada:", fname, "\n")
  }
}

# ── Varianza explicada PCA (scree plot) ───────────────────────────────────────
plot_scree <- function(pca_res, n_pcs = 15, base_dir = ".") {
  df <- data.frame(
    PC   = seq_len(n_pcs),
    var  = pca_res$var_explained[seq_len(n_pcs)],
    cum  = cumsum(pca_res$var_explained[seq_len(n_pcs)])
  )
  p <- ggplot(df, aes(x = PC)) +
    geom_col(aes(y = var), fill = "#2166AC", alpha = 0.8) +
    geom_line(aes(y = cum), colour = "#E41A1C", linewidth = 1) +
    geom_point(aes(y = cum), colour = "#E41A1C", size = 2.5) +
    geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = seq_len(n_pcs)) +
    labs(title = "Scree plot — Varianza explicada por PC",
         x = "Componente principal",
         y = "% varianza explicada / acumulada") +
    theme_strap()

  fname <- "scree_plot.png"
  ggsave(file.path(base_dir, "results", "figures", fname),
         plot = p, width = 9, height = 5, dpi = 300, bg = "white")
  cat("  Figura guardada:", fname, "\n")
}

# ── Orquestador de visualizaciones ────────────────────────────────────────────
generate_all_plots <- function(pca_res, dm_res, all_clusters, final_clusters,
                                cluster_labels, cluster_colors,
                                validation, deg_results, meta, expr, k,
                                pathotype_order, pathotype_col,
                                base_dir = ".") {
  best_method <- validation$ranking$method[1]

  # 1. Scree plot
  plot_scree(pca_res, base_dir = base_dir)

  # 2. PCA por método
  for (m in names(all_clusters)) {
    p <- plot_pca_clusters(pca_res, all_clusters[[m]],
                            cluster_labels, cluster_colors,
                            title = paste("PCA —", m))
    fname <- paste0("pca_", m, ".png")
    ggsave(file.path(base_dir, "results", "figures", fname),
           plot = p, width = 8, height = 6, dpi = 300, bg = "white")
    cat("  Figura guardada:", fname, "\n")
  }

  # 3. PCA con patotipos histológicos
  if (pathotype_col %in% colnames(meta)) {
    p <- plot_pca_histological(pca_res, meta, pathotype_col,
                                pathotype_order, cluster_colors)
    ggsave(file.path(base_dir, "results", "figures", "pca_histological.png"),
           plot = p, width = 8, height = 6, dpi = 300, bg = "white")
    cat("  Figura guardada: pca_histological.png\n")
  }

  # 4. Diffusion Map
  p_dm <- plot_dm_clusters(dm_res, final_clusters,
                            cluster_labels, cluster_colors,
                            title = paste("Diffusion Map —", best_method))
  ggsave(file.path(base_dir, "results", "figures", "diffusion_map.png"),
         plot = p_dm, width = 8, height = 6, dpi = 300, bg = "white")
  cat("  Figura guardada: diffusion_map.png\n")

  # 5. Comparación PCA vs Diffusion Maps
  p_comp <- plot_dimred_comparison(pca_res, dm_res, final_clusters,
                                    cluster_labels, cluster_colors)
  ggsave(file.path(base_dir, "results", "figures", "dimred_comparison.png"),
         plot = p_comp, width = 14, height = 6, dpi = 300, bg = "white")
  cat("  Figura guardada: dimred_comparison.png\n")

  # 6. Silhouette comparativo
  p_sil <- plot_silhouette_comparison(validation, cluster_colors)
  ggsave(file.path(base_dir, "results", "figures", "silhouette_comparison.png"),
         plot = p_sil, width = 9, height = 6, dpi = 300, bg = "white")
  cat("  Figura guardada: silhouette_comparison.png\n")

  # 7. Bootstrap stability
  p_boot <- plot_bootstrap_stability(validation)
  ggsave(file.path(base_dir, "results", "figures", "bootstrap_stability.png"),
         plot = p_boot, width = 9, height = 6, dpi = 300, bg = "white")
  cat("  Figura guardada: bootstrap_stability.png\n")

  # 8. Heatmap
  plot_heatmap(expr, final_clusters, cluster_labels, cluster_colors,
               base_dir = base_dir)

  # 9. Confusion matrix (%) — con nombres biológicos en eje Y
  plot_confusion_pct(validation, best_method,
                     cluster_labels = cluster_labels,
                     base_dir = base_dir)

  # 10. Sankey
  if (pathotype_col %in% colnames(meta)) {
    plot_sankey(final_clusters, meta, cluster_labels, cluster_colors,
                pathotype_col, pathotype_order, base_dir = base_dir)
  }

  # 11. Boxplots de genes marcadores
  if (!is.null(deg_results)) {
    plot_marker_boxplots(deg_results, expr, final_clusters,
                          cluster_labels, cluster_colors,
                          base_dir = base_dir)
  }

  cat("  Todas las visualizaciones generadas.\n")
}
