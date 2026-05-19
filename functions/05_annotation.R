# ==============================================================================
# functions/05_annotation.R
# Anotación biológica automática de clusters
#
# Estrategia en dos pasos:
#   1. Para k <= 3: asignar fibroide/mieloide/linfoide basándose en el
#      solapamiento máximo con patotipos histológicos conocidos.
#   2. Para k > 3 (nuevos subtipos): el cluster más parecido a cada patotipo
#      histórico recibe su nombre; los clusters sin correspondencia clara
#      reciben nombres basados en sus genes marcadores principales.
#
# Firmas de referencia (genes canónicos de AR):
#   Fibroid:  genes de fibroblastos sinoviales (CXCL12, PDPN, PRG4, THY1)
#   Myeloid:  genes de macrófagos/monocitos   (CD68, MRC1, CSF1R, IL1B)
#   Lymphoid: genes de linfocitos              (CD3D, CD19, MS4A1, SELL)
#   IFN-high: interferón                       (MX1, OAS1, IFIT1, ISG15)
# ==============================================================================

# ── Firmas génicas de referencia para AR sinovial ─────────────────────────────
RA_SIGNATURES <- list(
  Fibroid  = c("CXCL12", "PDPN", "PRG4", "THY1", "CDH11", "FAP",
               "ACTA2", "COL1A1", "COL3A1", "FN1"),
  Myeloid  = c("CD68", "MRC1", "CSF1R", "IL1B", "TNF", "IL6",
               "ITGAM", "CD14", "FCGR3A", "CCL2"),
  Lymphoid = c("CD3D", "CD3E", "CD19", "MS4A1", "SELL", "PTPRC",
               "CD4", "CD8A", "FOXP3", "CXCR5"),
  IFN_high = c("MX1", "OAS1", "IFIT1", "IFIT3", "ISG15", "RSAD2",
               "IFI44L", "CXCL10", "IRF7", "STAT1")
)

# ── Puntuación de firma para cada muestra ────────────────────────────────────
score_signature <- function(expr_scaled, sig_genes) {
  # expr_scaled: muestras × genes
  available <- intersect(sig_genes, colnames(expr_scaled))
  if (length(available) == 0) return(rep(0, nrow(expr_scaled)))
  rowMeans(expr_scaled[, available, drop = FALSE])
}

# ── Asignar nombre biológico a cada cluster ───────────────────────────────────
annotate_clusters <- function(clusters, expr, meta, k,
                               pathotype_order, pathotype_colors,
                               n_markers = 20, seed = 42) {
  set.seed(seed)

  n_clusters  <- length(unique(clusters))
  cluster_ids <- sort(unique(clusters))

  # Paso 1: Puntuar cada muestra en cada firma
  sig_scores <- lapply(RA_SIGNATURES, score_signature, expr_scaled = expr)
  sig_mat    <- do.call(cbind, sig_scores)   # muestras × firmas

  # Paso 2: Puntuación media de firma por cluster
  cluster_sig <- sapply(cluster_ids, function(cl) {
    idx <- which(clusters == cl)
    colMeans(sig_mat[idx, , drop = FALSE])
  })
  colnames(cluster_sig) <- cluster_ids

  # Paso 3: Usar solapamiento histológico para anclar nombres clásicos
  if ("pathotype_clean" %in% colnames(meta)) {
    overlap <- compute_biological_overlap(clusters, meta, "pathotype_clean")
    ct_pct  <- overlap$table_pct   # cluster × patotipo, % por cluster

    # Para cada cluster, patotipo histológico dominante
    histological_match <- if (!is.null(ct_pct)) {
      apply(ct_pct, 1, function(row) {
        colnames(ct_pct)[which.max(row)]
      })
    } else NULL
  } else {
    histological_match <- NULL
  }

  # Paso 4: Construir labels combinando firma génica + solapamiento histológico
  labels <- character(length(cluster_ids))
  names(labels) <- cluster_ids
  used_names <- character(0)

  # Orden de prioridad de asignación
  priority <- c("Fibroid", "Myeloid", "Lymphoid", "IFN_high")

  # Asignar primero los clusters con correspondencia histológica clara (>40%)
  if (!is.null(histological_match) && !is.null(ct_pct)) {
    for (cl in cluster_ids) {
      cl_char   <- as.character(cl)
      hist_name <- histological_match[cl_char]
      max_pct   <- max(ct_pct[cl_char, ], na.rm = TRUE)

      if (!is.na(hist_name) && max_pct >= 40 && !hist_name %in% used_names) {
        labels[cl_char] <- hist_name
        used_names       <- c(used_names, hist_name)
      }
    }
  }

  # Para clusters sin asignación histológica clara: usar firma génica
  sig_name_map <- c(Fibroid  = "Fibroid",
                    Myeloid  = "Myeloid",
                    Lymphoid = "Lymphoid",
                    IFN_high = "IFN-high")

  for (cl in cluster_ids) {
    cl_char <- as.character(cl)
    if (labels[cl_char] == "") {
      # Firma con mayor score para este cluster
      best_sig  <- names(which.max(cluster_sig[, cl_char]))
      candidate <- sig_name_map[best_sig]

      if (!candidate %in% used_names) {
        labels[cl_char] <- candidate
        used_names       <- c(used_names, candidate)
      } else {
        # Si el nombre ya está usado, asignar nombre genérico numérico
        fallback         <- paste0("Subtype", cl)
        labels[cl_char] <- fallback
        used_names       <- c(used_names, fallback)
      }
    }
  }

  # Paso 5: Ordenar labels según orden biológico canónico
  # (fibroide, mieloide, linfoide, luego nuevos por orden de cluster_id)
  canonical_order <- c(pathotype_order, "IFN-high",
                       paste0("Subtype", 1:10))
  ordered_labels  <- labels[order(match(labels, canonical_order,
                                        nomatch = 999))]

  # Paso 6: Asignar colores consistentes en el orden canónico
  colors_ordered <- pathotype_colors[seq_len(n_clusters)]
  if (length(colors_ordered) < n_clusters) {
    extra <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(8, "Set1"))(n_clusters - length(colors_ordered))
    colors_ordered <- c(colors_ordered, extra)
  }
  names(colors_ordered) <- unname(ordered_labels)

  cat("  Anotación biológica completada:\n")
  for (cl in names(ordered_labels)) {
    cat(sprintf("    Cluster %s → %s\n", cl, ordered_labels[cl]))
  }

  return(list(
    labels       = labels,           # cluster_id → nombre
    ordered      = ordered_labels,   # en orden canónico
    colors       = colors_ordered,   # nombre → color (hex)
    sig_scores   = sig_mat,
    cluster_sig  = cluster_sig
  ))
}

# ── Reetiqueta un vector de clusters con los nombres biológicos ───────────────
relabel_clusters <- function(cl, cluster_labels) {
  # cl: vector integer (cluster_id)
  # cluster_labels: vector nombrado cluster_id → nombre biológico
  factor(cluster_labels[as.character(cl)],
         levels = unname(cluster_labels))
}
