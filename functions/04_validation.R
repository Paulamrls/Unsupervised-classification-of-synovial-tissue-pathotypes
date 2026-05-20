# ==============================================================================
# functions/04_validation.R
# Validación matemática y biológica de clusters
#
# VALIDACIÓN MATEMÁTICA:
#   - Silhouette score (por método y global)
#   - Bootstrap stability (Adjusted Rand Index entre submuestras)
#   - Comparación entre métodos (ranking)
#
# VALIDACIÓN BIOLÓGICA:
#   - Solapamiento con patotipos histológicos (ARI, tabla de contingencia)
#   - Enriquecimiento en firmas de AR conocidas (preparado)
# ==============================================================================

# ── Silhouette score para un vector de clusters ───────────────────────────────
compute_silhouette <- function(cl, expr) {
  # Eliminar muestras sin asignación de cluster válida
  valid <- !is.na(cl) & cl > 0
  if (sum(valid) < 3 || length(unique(cl[valid])) < 2) return(NA_real_)

  dist_mat <- dist(expr[valid, ])
  sil      <- silhouette(cl[valid], dist_mat)
  mean(sil[, "sil_width"])
}

# ── Estabilidad bootstrap ─────────────────────────────────────────────────────
# Remuestrea el 80% de muestras B veces, reclustra, y mide ARI vs. referencia.
compute_bootstrap_stability <- function(cl_ref, expr, k, method_fn,
                                        boot_n = 100, frac = 0.8,
                                        seed = 42) {
  set.seed(seed)
  n       <- nrow(expr)
  ari_vec <- numeric(boot_n)

  for (b in seq_len(boot_n)) {
    idx    <- sample(n, size = floor(n * frac), replace = FALSE)
    expr_b <- expr[idx, ]

    cl_b   <- tryCatch(
      method_fn(expr_b, k, seed = seed + b),
      error = function(e) NULL
    )

    if (is.null(cl_b)) {
      ari_vec[b] <- NA_real_
      next
    }

    # Alinear muestras entre referencia y bootstrap
    common   <- names(cl_b)
    ari_vec[b] <- adjustedRandIndex(cl_ref[common], cl_b)
  }

  list(
    mean_ari = mean(ari_vec, na.rm = TRUE),
    sd_ari   = sd(ari_vec,   na.rm = TRUE),
    ari_vec  = ari_vec
  )
}

# ── Validación biológica: solapamiento con patotipos histológicos ─────────────
compute_biological_overlap <- function(cl, meta, pathotype_col = "pathotype_clean") {
  if (!pathotype_col %in% colnames(meta))
    return(list(ari = NA, table_pct = NULL))

  common    <- intersect(names(cl), rownames(meta))
  cl_common <- cl[common]
  pt_common <- meta[common, pathotype_col]

  ari <- adjustedRandIndex(cl_common, pt_common)

  # Tabla de contingencia normalizada a porcentajes por cluster
  ct      <- table(cluster = cl_common, pathotype = pt_common)
  ct_pct  <- prop.table(ct, margin = 1) * 100   # % por fila (cluster)

  list(ari = ari, table_abs = ct, table_pct = round(ct_pct, 1))
}

# ── Orquestador de validación completa ────────────────────────────────────────
run_validation <- function(expr, clusters, meta, k,
                            boot_n = 100, seed = 42,
                            pathotype_col = "pathotype_clean") {
  cat("  Calculando silhouette por método...\n")

  # Mapeo método → función de clustering (para bootstrap)
  method_fns <- list(
    kmeans       = function(e, k, seed) cluster_kmeans(e, k, seed),
    hierarchical = function(e, k, seed) cluster_hierarchical(e, k, seed),
    spectral     = function(e, k, seed) cluster_spectral(e, k, seed),
    leiden       = function(e, k, seed) cluster_leiden(e, k, seed),
    mcl          = function(e, k, seed) cluster_mcl(e, k, seed),
    nmf          = function(e, k, seed) cluster_nmf(e, k, seed),
    gmm          = function(e, k, seed) cluster_gmm(e, k, seed)
  )

  methods_avail <- names(clusters)

  # ── Silhouette ────────────────────────────────────────────────────────────
  sil_scores <- sapply(methods_avail, function(m) {
    compute_silhouette(clusters[[m]], expr)
  })

  # ── Bootstrap stability ───────────────────────────────────────────────────
  cat("  Calculando estabilidad bootstrap (", boot_n, "iters)...\n")
  boot_results <- lapply(methods_avail, function(m) {
    if (!m %in% names(method_fns)) {
      return(list(mean_ari = NA, sd_ari = NA, ari_vec = NULL))
    }
    compute_bootstrap_stability(
      cl_ref    = clusters[[m]],
      expr      = expr,
      k         = k,
      method_fn = method_fns[[m]],
      boot_n    = boot_n,
      seed      = seed
    )
  })
  names(boot_results) <- methods_avail

  boot_ari <- sapply(boot_results, `[[`, "mean_ari")
  boot_sd  <- sapply(boot_results, `[[`, "sd_ari")

  # ── Validación biológica ──────────────────────────────────────────────────
  cat("  Calculando solapamiento biológico con patotipos histológicos...\n")
  bio_overlap <- lapply(methods_avail, function(m) {
    compute_biological_overlap(clusters[[m]], meta, pathotype_col)
  })
  names(bio_overlap) <- methods_avail

  bio_ari <- sapply(bio_overlap, `[[`, "ari")

  # ── Ranking combinado ─────────────────────────────────────────────────────
  ranking <- data.frame(
    method      = methods_avail,
    silhouette  = round(sil_scores,  4),
    boot_ari    = round(boot_ari,    4),
    boot_sd     = round(boot_sd,     4),
    bio_ari     = round(bio_ari,     4),
    stringsAsFactors = FALSE
  )

  # Ordenar por silhouette (mayor → mejor)
  ranking <- ranking[order(ranking$silhouette, decreasing = TRUE), ]
  rownames(ranking) <- NULL

  cat("  Ranking de métodos:\n")
  print(ranking)

  return(list(
    sil_scores   = sil_scores,
    boot_results = boot_results,
    bio_overlap  = bio_overlap,
    ranking      = ranking
  ))
}
