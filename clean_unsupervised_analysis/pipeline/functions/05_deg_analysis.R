################################################################################
# functions/07_deg_analysis.R
################################################################################

# 1. DEGs uno vs. el resto
run_deseq2_one_vs_rest <- function(expr_raw, clusters, cluster_id, seed = 42) {
  set.seed(seed)

  samples <- names(clusters)

  # Alineación de muestras
  common   <- intersect(samples, colnames(expr_raw))
  counts   <- round(expr_raw[, common])
  cl_sub   <- clusters[common]

  # Factor binario
  group <- ifelse(cl_sub == cluster_id, "target", "rest")
  col_data <- data.frame(group = factor(group, levels = c("rest", "target")),
                          row.names = common)

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData   = col_data,
                                 design    = ~group)

  # Filtrar genes de baja expresión para ir más rápidos
  dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]

  dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)

  res <- results(dds, contrast = c("group", "target", "rest"),
                 alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df <- res_df[order(res_df$padj), ]

  return(res_df)
}

# 2. Separar genes UP y DOWN, ordenándolos por relevancia

split_deg_table <- function(deg_df, lfc_threshold = 1.0, padj_threshold = 0.05) {
  sig <- deg_df[deg_df$padj < padj_threshold &
                  abs(deg_df$log2FoldChange) >= lfc_threshold, ]

  up   <- sig[sig$log2FoldChange > 0, ]
  down <- sig[sig$log2FoldChange < 0, ]

  # Ordenar primero  p-adj primero, y luego por |LFC|
  up   <- up[order(up$padj, -abs(up$log2FoldChange)), ]
  down <- down[order(down$padj, -abs(down$log2FoldChange)), ]

  list(up = up, down = down, all_sig = sig)
}

# 3.DEG para todos los clusters

run_deg_analysis <- function(expr, expr_raw, clusters, cluster_labels,
                              n_markers = 20, seed = 42,
                              lfc_threshold = 0.5, padj_threshold = 0.05) {
  cluster_ids <- sort(unique(clusters))
  up_tables   <- list()
  down_tables <- list()
  all_tables  <- list()

  for (cl_id in cluster_ids) {
    cl_name <- cluster_labels[as.character(cl_id)]
    cat(sprintf("  DEG: %s (cluster %d) vs. resto...\n", cl_name, cl_id))

    deg_df <- tryCatch(
      run_deseq2_one_vs_rest(expr_raw, clusters, cl_id, seed),
      error = function(e) {
        warning("DESeq2 falló para cluster ", cl_id, ": ", e$message)
        NULL
      }
    )

    if (is.null(deg_df)) next

    split <- split_deg_table(deg_df, lfc_threshold, padj_threshold)

    up_tables[[cl_name]]  <- split$up
    down_tables[[cl_name]] <- split$down
    all_tables[[cl_name]] <- deg_df

    cat(sprintf("    → %d up, %d down (padj<%.2f, |LFC|>=%.1f)\n",
                nrow(split$up), nrow(split$down),
                padj_threshold, lfc_threshold))
  }

  #  Top x marcadores por cada cluster (para una mejor visualización de los resultados)
  top_markers <- lapply(names(all_tables), function(cl_name) {
    df <- all_tables[[cl_name]]
    df <- df[!is.na(df$padj) & df$padj < padj_threshold, ]
    df <- df[order(df$padj), ]
    head(df, n_markers)
  })
  names(top_markers) <- names(all_tables)

  return(list(
    up          = up_tables,
    down        = down_tables,
    all         = all_tables,
    top_markers = top_markers
  ))
}
