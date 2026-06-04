# ==============================================================================
# functions/08_export.R
# Exportación organizada de todos los resultados
# ==============================================================================

export_all_results <- function(all_clusters, final_clusters, cluster_labels,
                                validation, deg_results, meta, k,
                                base_dir = ".") {

  best_method <- validation$ranking$method[1]

  # ── Asignaciones de cluster ─────────────────────────────────────────────────
  cat("  Exportando asignaciones de cluster...\n")

  # Tabla maestra de muestras + cluster
  sample_table <- data.frame(
    sample_id   = names(final_clusters),
    cluster_id  = final_clusters,
    cluster_bio = cluster_labels[as.character(final_clusters)],
    row.names   = NULL
  )
  # Añadir patotipo histológico si está disponible
  if ("pathotype_clean" %in% colnames(meta)) {
    common <- intersect(sample_table$sample_id, rownames(meta))
    idx    <- match(sample_table$sample_id, rownames(meta))
    sample_table$pathotype_histological <- meta$pathotype_clean[idx]
  }
  save_table(sample_table, "cluster_assignments.csv",
             base_dir = base_dir, subdir = "clusters")

  # Tabla por cada método
  for (m in names(all_clusters)) {
    cl_df <- data.frame(
      sample_id  = names(all_clusters[[m]]),
      cluster_id = all_clusters[[m]],
      row.names  = NULL
    )
    save_table(cl_df,
               paste0("clusters_", m, ".csv"),
               base_dir = base_dir, subdir = "clusters")
  }

  # ── Ranking y métricas de validación ────────────────────────────────────────
  cat("  Exportando validación...\n")
  save_table(validation$ranking,
             "ranking_methods.csv",
             base_dir = base_dir, subdir = "validation")

  # Tabla de solapamiento biológico (% por cluster) para el mejor método
  bio_table <- validation$bio_overlap[[best_method]]$table_pct
  if (!is.null(bio_table)) {
    df_bio <- as.data.frame.matrix(bio_table)
    df_bio$cluster_bio <- cluster_labels[rownames(df_bio)]
    save_table(df_bio, "confusion_pct_best_method.csv",
               base_dir = base_dir, subdir = "validation")
  }

  # ── DEGs separados por patotipo (up y down) ──────────────────────────────────
  cat("  Exportando genes diferenciales...\n")
  for (cl_name in names(deg_results$up)) {
    cl_safe <- gsub("[^A-Za-z0-9]", "_", cl_name)

    # Upregulated
    df_up <- deg_results$up[[cl_name]]
    if (nrow(df_up) > 0) {
      df_up$gene <- rownames(df_up)
      save_table(df_up[, c("gene", "log2FoldChange", "padj",
                            "pvalue", "baseMean", "stat")],
                 paste0("DEGs_UP_", cl_safe, ".csv"),
                 base_dir = base_dir, subdir = "DEGs")
    }

    # Downregulated
    df_down <- deg_results$down[[cl_name]]
    if (nrow(df_down) > 0) {
      df_down$gene <- rownames(df_down)
      save_table(df_down[, c("gene", "log2FoldChange", "padj",
                              "pvalue", "baseMean", "stat")],
                 paste0("DEGs_DOWN_", cl_safe, ".csv"),
                 base_dir = base_dir, subdir = "DEGs")
    }

    # Todos los DEGs significativos juntos (ordenados por padj)
    df_all <- deg_results$all[[cl_name]]
    if (!is.null(df_all) && nrow(df_all) > 0) {
      df_all$gene <- rownames(df_all)
      save_table(df_all[, c("gene", "log2FoldChange", "padj",
                             "pvalue", "baseMean", "stat")],
                 paste0("DEGs_ALL_", cl_safe, ".csv"),
                 base_dir = base_dir, subdir = "DEGs")
    }
  }

  # ── Tabla resumen de la anotación ───────────────────────────────────────────
  cat("  Exportando tabla de anotación...\n")
  ann_summary <- data.frame(
    cluster_id  = as.integer(names(cluster_labels)),
    cluster_bio = unname(cluster_labels),
    n_samples   = as.integer(table(final_clusters)[names(cluster_labels)]),
    row.names   = NULL
  )
  save_table(ann_summary, "cluster_annotation.csv",
             base_dir = base_dir, subdir = "tables")

  # ── Parámetros del análisis (reproducibilidad) ─────────────────────────────
  cat("  Exportando parámetros...\n")
  params <- data.frame(
    parameter = c("k", "global_seed", "var_quantile_cutoff", "best_method",
                  "n_methods_run", "date"),
    value     = c(k, GLOBAL_SEED, VAR_CUTOFF, best_method,
                  length(all_clusters), as.character(Sys.Date()))
  )
  save_table(params, "analysis_parameters.csv",
             base_dir = base_dir, subdir = "tables")

  cat("  Exportación completada.\n")
}
