# ==============================================================================
# functions/03_clustering.R
# Todos los métodos de clustering — arquitectura modular para añadir más fácilmente
# ==============================================================================
# MÉTODOS IMPLEMENTADOS:
#   kmeans       — K-means clásico
#   hierarchical — Ward D2 jerárquico
#   spectral     — Spectral clustering
#   consensus    — ConsensusClusterPlus (reproducible con seed)
#   leiden       — Leiden (detección de comunidades en grafo)
#   mcl          — Markov Clustering (MCL)
#   nmf          — NMF: descompone en k programas transcriptómicos (metagenes)
#                  → especialmente adecuado para bulk RNA-seq (señales mixtas)
#   gmm          — GMM: mezcla de k gaussianas multivariantes sobre PCA
#                  → clustering probabilístico con base estadística formal
#
# ESTRUCTURA MODULAR: Para añadir un nuevo método, añade una función
#   cluster_METODO(expr, k, seed, ...) que devuelva un vector integer con
#   los clusters (nombres = rownames(expr)), y regístralo en run_all_clustering().
# ==============================================================================


# ── Constructor de grafo k-NN compartido por métodos de grafo ─────────────────
build_knn_graph <- function(expr, k_nn = 10, seed = 42) {
  set.seed(seed)

  # Matriz de distancias euclídeas
  dist_mat <- as.matrix(dist(expr))

  n <- nrow(dist_mat)
  adj <- matrix(0, n, n, dimnames = list(rownames(expr), rownames(expr)))

  for (i in seq_len(n)) {
    # Índices de los k_nn vecinos más cercanos (excluir la propia muestra)
    nn_idx <- order(dist_mat[i, ])[-1][seq_len(k_nn)]
    adj[i, nn_idx] <- 1
    adj[nn_idx, i] <- 1   # grafo no dirigido
  }

  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = NULL)
  return(g)
}

# ── Normalizar clusters a enteros consecutivos 1..k ──────────────────────────
normalize_clusters <- function(cl) {
  cl_fac <- factor(cl)
  as.integer(cl_fac)
}


# ── K-means ───────────────────────────────────────────────────────────────────
cluster_kmeans <- function(expr, k, seed = 42, nstart = 50) {
  set.seed(seed)
  km <- kmeans(expr, centers = k, nstart = nstart, iter.max = 300)
  cl <- setNames(km$cluster, rownames(expr))
  cat(sprintf("    K-means: tabla = %s\n",
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# ── Jerárquico (Ward D2) ──────────────────────────────────────────────────────
cluster_hierarchical <- function(expr, k, seed = 42) {
  dist_mat <- dist(expr, method = "euclidean")
  hc       <- hclust(dist_mat, method = "ward.D2")
  cl       <- cutree(hc, k = k)
  cl       <- setNames(as.integer(cl), rownames(expr))
  cat(sprintf("    Jerárquico: tabla = %s\n",
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# ── Spectral ──────────────────────────────────────────────────────────────────
cluster_spectral <- function(expr, k, seed = 42) {
  set.seed(seed)
  sp <- kernlab::specc(expr, centers = k)
  cl <- setNames(as.integer(sp), rownames(expr))
  cat(sprintf("    Spectral: tabla = %s\n",
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# ── Consensus clustering (reproducible) ───────────────────────────────────────
# ConsensusClusterPlus tiene randomness interno; se controla con seed.
cluster_consensus <- function(expr, k, seed = 42,
                               max_k = 6, reps = 100,
                               p_item = 0.8, p_feature = 1.0,
                               alg = "hc", dist_metric = "pearson",
                               base_dir = ".") {
  set.seed(seed)

  out_dir <- file.path(base_dir, "results", "clusters", "consensus_plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ConsensusClusterPlus espera genes × muestras (transponer)
  cc_res <- ConsensusClusterPlus(
    d           = t(expr),
    maxK        = max_k,
    reps        = reps,
    pItem       = p_item,
    pFeature    = p_feature,
    clusterAlg  = alg,
    distance    = dist_metric,
    seed        = seed,     # seed interno para reproducibilidad
    plot        = "png",
    title       = out_dir
  )

  cl <- setNames(as.integer(cc_res[[k]]$consensusClass), rownames(expr))
  cat(sprintf("    Consensus (k=%d): tabla = %s\n", k,
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# ── Leiden ────────────────────────────────────────────────────────────────────
# igraph::cluster_leiden usa "resolution" en igraph >= 2.0
# (antes era "resolution_parameter"); probamos ambos para compatibilidad.
cluster_leiden <- function(expr, k, seed = 42, k_nn = 10, resolution = 1.0) {
  set.seed(seed)
  g       <- build_knn_graph(expr, k_nn = k_nn, seed = seed)
  res_val <- resolution
  cl_vec  <- NULL

  for (attempt in 1:20) {
    cl_raw <- tryCatch(
      igraph::cluster_leiden(g, resolution = res_val, n_iterations = 10),
      error = function(e)
        igraph::cluster_leiden(g, resolution_parameter = res_val,
                               objective_function = "modularity",
                               n_iterations = 10)
    )
    cl_vec  <- as.integer(igraph::membership(cl_raw))
    n_found <- length(unique(cl_vec))
    if (n_found == k) break
    res_val <- if (n_found < k) res_val * 1.2 else res_val * 0.85
  }

  cl <- setNames(normalize_clusters(cl_vec), rownames(expr))
  cat(sprintf("    Leiden (res=%.3f): %d clusters -> tabla = %s\n",
              res_val, length(unique(cl)),
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# ── NMF (Non-negative Matrix Factorization) ───────────────────────────────────
# Descompone la matriz de expresión en k programas transcriptómicos (metagenes).
# Cada muestra queda descrita por su peso en cada metagén; el cluster asignado
# es el metagén dominante.
# Especialmente apropiado para bulk RNA-seq: captura señales mixtas de tipos
# celulares que coexisten en la misma biopsia.
# NOTA: nrun iteraciones para seleccionar la mejor solución — puede ser lento.
cluster_nmf <- function(expr, k, seed = 42, nrun = 20) {
  set.seed(seed)

  # NMF requiere valores no negativos
  # Desplazar la matriz z-score: sumar |min| para que el mínimo sea 0
  expr_nn <- expr - min(expr)

  # NMF espera genes × muestras → transponer
  res <- NMF::nmf(t(expr_nn), rank = k, nrun = nrun, seed = seed,
                   .options = "-v")   # sin verbose

  # Asignación: metagén de mayor peso para cada muestra
  cl_vec <- as.integer(NMF::predict(res, what = "samples"))
  cl     <- setNames(normalize_clusters(cl_vec), rownames(expr))

  cat(sprintf("    NMF (k=%d, nrun=%d): tabla = %s\n", k, nrun,
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# ── GMM (Gaussian Mixture Models) ─────────────────────────────────────────────
# Ajusta una mezcla de k gaussianas multivariantes sobre los primeros PCs.
# La reducción PCA previa es necesaria: mclust no escala bien en alta dimensión
# (1000 genes). Con 30 PCs capturamos la mayor parte de la varianza estructural.
cluster_gmm <- function(expr, k, seed = 42, n_pcs = 30) {
  set.seed(seed)

  # Reducción dimensional previa (expr ya está centrado/escalado)
  pca    <- prcomp(expr, center = FALSE, scale. = FALSE)
  n_use  <- min(n_pcs, ncol(pca$x), nrow(expr) - 1)
  expr_pc <- pca$x[, seq_len(n_use), drop = FALSE]

  res <- mclust::Mclust(expr_pc, G = k, verbose = FALSE)

  if (is.null(res)) {
    warning("GMM no convergió para k=", k, ". Devolviendo NULL.")
    return(NULL)
  }

  cl <- setNames(as.integer(res$classification), rownames(expr))
  cat(sprintf("    GMM (k=%d, modelo=%s, PCs=%d): tabla = %s\n",
              k, res$modelName, n_use,
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# ── MCL (Markov Clustering) ───────────────────────────────────────────────────
cluster_mcl <- function(expr, k, seed = 42, k_nn = 10, expansion = 2,
                         inflation = 2) {
  set.seed(seed)
  g       <- build_knn_graph(expr, k_nn = k_nn, seed = seed)
  adj_mat <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))

  mcl_res <- mcl(adj_mat, addLoops = TRUE, expansion = expansion,
                 inflation = inflation, allow1 = TRUE, max.iter = 100,
                 ESM = FALSE)

  cl_vec <- mcl_res$Cluster
  cl     <- setNames(normalize_clusters(cl_vec), rownames(expr))
  cat(sprintf("    MCL: %d clusters → tabla = %s\n",
              length(unique(cl)),
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}


# ── Orquestador: ejecuta todos los métodos ────────────────────────────────────
# Para añadir un nuevo método: añade su función cluster_X() arriba,
# y súmala a la lista `methods` de esta función.
run_all_clustering <- function(expr, k, seed = 42,
                                cc_max_k = 6, cc_reps = 100,
                                cc_p_item = 0.8, cc_p_feature = 1.0,
                                cc_alg = "hc", cc_dist = "pearson",
                                graph_k_nn = 10, graph_res = 1.0,
                                base_dir = ".") {
  results <- list()

  # ── Métodos clásicos ─────────────────────────────────────────────────────
  cat("  → K-means\n")
  results$kmeans       <- tryCatch(
    cluster_kmeans(expr, k, seed), error = function(e) {
      warning("K-means falló: ", e$message); NULL })

  cat("  → Jerárquico\n")
  results$hierarchical <- tryCatch(
    cluster_hierarchical(expr, k, seed), error = function(e) {
      warning("Jerárquico falló: ", e$message); NULL })

  cat("  → Spectral\n")
  results$spectral     <- tryCatch(
    cluster_spectral(expr, k, seed), error = function(e) {
      warning("Spectral falló: ", e$message); NULL })

  cat("  → Consensus\n")
  results$consensus    <- tryCatch(
    cluster_consensus(expr, k, seed, cc_max_k, cc_reps,
                      cc_p_item, cc_p_feature, cc_alg, cc_dist, base_dir),
    error = function(e) { warning("Consensus falló: ", e$message); NULL })

  # ── Métodos de grafo ─────────────────────────────────────────────────────
  cat("  → Leiden\n")
  results$leiden       <- tryCatch(
    cluster_leiden(expr, k, seed, graph_k_nn, graph_res),
    error = function(e) { warning("Leiden falló: ", e$message); NULL })

  cat("  → MCL\n")
  results$mcl          <- tryCatch(
    cluster_mcl(expr, k, seed, graph_k_nn), error = function(e) {
      warning("MCL falló: ", e$message); NULL })

  # ── Métodos específicos para bulk RNA-seq ─────────────────────────────────
  cat("  → NMF (puede tardar unos minutos)\n")
  results$nmf          <- tryCatch(
    cluster_nmf(expr, k, seed), error = function(e) {
      warning("NMF falló: ", e$message); NULL })

  cat("  → GMM\n")
  results$gmm          <- tryCatch(
    cluster_gmm(expr, k, seed), error = function(e) {
      warning("GMM falló: ", e$message); NULL })

  # Eliminar métodos que fallaron
  results <- Filter(Negate(is.null), results)

  cat(sprintf("  Métodos ejecutados con éxito: %d/%d\n",
              length(results), 8))
  return(results)
}


# ==============================================================================
# ARQUITECTURA PREPARADA PARA MÉTODOS AVANZADOS (implementación futura)
# ==============================================================================
# Los siguientes módulos tienen la estructura definida y están listos para
# implementación completa cuando se requiera:
#
# cluster_nmf()      — NMF: detecta programas transcriptómicos superpuestos
#                      (rank = k; usar NMF::nmf con nrun >= 30)
# cluster_biclustering() — biclúster simultáneo de muestras y genes
#                          (usar biclust o QUBIC)
# cluster_leiden_nmf() — Leiden + NMF combinado:
#                         1. Leiden para partición de grafo
#                         2. NMF dentro de cada comunidad para metagenes
#                         → Identifica estados híbridos y programas superpuestos
#
# Para implementar, copia la firma:
#   cluster_X <- function(expr, k, seed = 42, ...) { ... }
# y añade el tryCatch correspondiente en run_all_clustering().
# ==============================================================================
