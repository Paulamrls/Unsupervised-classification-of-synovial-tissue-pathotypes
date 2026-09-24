# ==============================================================================
# functions/07_deg_analysis.R
# Análisis de expresión diferencial por patotipo/cluster
# Usa DESeq2 (counts crudos) para estadística robusta
# Exporta tablas separadas de genes UP y DOWN por cluster
# ==============================================================================

# ── DEGs con DESeq2: un cluster vs. el resto ──────────────────────────────────
run_deseq2_one_vs_rest <- function(expr_raw, clusters, cluster_id, seed = 42) {
  set.seed(seed)

  samples <- names(clusters)

  # Alinear muestras
  common   <- intersect(samples, colnames(expr_raw))
  counts   <- round(expr_raw[, common])
  cl_sub   <- clusters[common]

  # Factor binario: este cluster vs. el resto
  group <- ifelse(cl_sub == cluster_id, "target", "rest")
  col_data <- data.frame(group = factor(group, levels = c("rest", "target")),
                          row.names = common)

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData   = col_data,
                                 design    = ~group)

  # Filtrar genes de baja expresión para acelerar
  dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]

  dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)

  res <- results(dds, contrast = c("group", "target", "rest"),
                 alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df <- res_df[order(res_df$padj), ]

  return(res_df)
}

# ── Separar UP y DOWN, ordenados por relevancia ───────────────────────────────
split_deg_table <- function(deg_df, lfc_threshold = 1.0, padj_threshold = 0.05) {
  sig <- deg_df[deg_df$padj < padj_threshold &
                  abs(deg_df$log2FoldChange) >= lfc_threshold, ]

  up   <- sig[sig$log2FoldChange > 0, ]
  down <- sig[sig$log2FoldChange < 0, ]

  # Ordenar: p-adj primero, luego por |LFC|
  up   <- up[order(up$padj, -abs(up$log2FoldChange)), ]
  down <- down[order(down$padj, -abs(down$log2FoldChange)), ]

  list(up = up, down = down, all_sig = sig)
}

# ── Orquestador DEG para todos los clusters ───────────────────────────────────
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

  # Top n marcadores positivos por cluster. Se toman exclusivamente de la
  # tabla UP ya filtrada por FDR y log2FC, no de todos los genes significativos.
  top_markers <- lapply(names(up_tables), function(cl_name) {
    df <- up_tables[[cl_name]]
    df <- df[order(df$padj, -df$log2FoldChange), , drop = FALSE]
    head(df, n_markers)
  })
  names(top_markers) <- names(up_tables)

  return(list(
    up          = up_tables,
    down        = down_tables,
    all         = all_tables,
    top_markers = top_markers
  ))
}
