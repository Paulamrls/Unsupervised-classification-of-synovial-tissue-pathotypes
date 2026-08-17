###################################################################################
# functions/03_clustering.R
# Todos los métodos de clustering — arquitectura modular para añadir más fácilmente
###################################################################################


# 1. Constructor de grafo k-NN compartido por métodos basados en grafo
build_knn_graph <- function(expr, k_nn = 10, seed = 42) {
  set.seed(seed)

  # 1.1. Matriz de distancias euclídeas
  dist_mat <- as.matrix(dist(expr))

  n <- nrow(dist_mat)
  adj <- matrix(0, n, n, dimnames = list(rownames(expr), rownames(expr)))

  for (i in seq_len(n)) {
    # 1.2.Índices de los k_nn vecinos más cercanos 
    nn_idx <- order(dist_mat[i, ])[-1][seq_len(k_nn)]
    adj[i, nn_idx] <- 1
    adj[nn_idx, i] <- 1  
  }

  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = NULL)
  return(g)
}

# 2. Normalización de clusters a enteros consecutivos 1..k

normalize_clusters <- function(cl) {
  cl_fac <- factor(cl)
  as.integer(cl_fac)
}

# 3. Clustering no supervisado

#  3.1. K-means
cluster_kmeans <- function(expr, k, seed = 42, nstart = 50) {
  set.seed(seed)
  km <- kmeans(expr, centers = k, nstart = nstart, iter.max = 300)
  cl <- setNames(km$cluster, rownames(expr))
  cat(sprintf("    K-means: tabla = %s\n",
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

#  3.2. Jerárquico (Ward D2)
cluster_hierarchical <- function(expr, k, seed = 42) {
  dist_mat <- dist(expr, method = "euclidean")
  hc       <- hclust(dist_mat, method = "ward.D2")
  cl       <- cutree(hc, k = k)
  cl       <- setNames(as.integer(cl), rownames(expr))
  cat(sprintf("    Jerárquico: tabla = %s\n",
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# 3.3. Spectral
cluster_spectral <- function(expr, k, seed = 42) {
  set.seed(seed)
  sp <- kernlab::specc(expr, centers = k)
  cl <- setNames(as.integer(sp), rownames(expr))
  cat(sprintf("    Spectral: tabla = %s\n",
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

#  3.4. Consensus clustering 
cluster_consensus <- function(expr, k, seed = 42,
                               max_k = 6, reps = 100,
                               p_item = 0.8, p_feature = 1.0,
                               alg = "hc", dist_metric = "pearson",
                               base_dir = ".") {
  set.seed(seed) # necesitamos seed

  out_dir <- file.path(base_dir, "results", "clusters", "consensus_plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ConsensusClusterPlus espera genes × muestras 
  cc_res <- ConsensusClusterPlus(
    d           = t(expr),
    maxK        = max_k,
    reps        = reps,
    pItem       = p_item,
    pFeature    = p_feature,
    clusterAlg  = alg,
    distance    = dist_metric,
    seed        = seed,     
    plot        = "png",
    title       = out_dir
  )

  cl <- setNames(as.integer(cc_res[[k]]$consensusClass), rownames(expr))
  cat(sprintf("    Consensus (k=%d): tabla = %s\n", k,
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# 3.5. Leiden

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

# 3.6.  NMF (Non-negative Matrix Factorization)

cluster_nmf <- function(expr, k, seed = 42, nrun = 5, n_genes = 200) { # He elegido 5 run para que el tiempo de ejecución no sean muy largo
  set.seed(seed)

  var_genes <- apply(expr, 2, var)
  top_genes <- names(sort(var_genes, decreasing = TRUE))[seq_len(min(n_genes, ncol(expr)))]
  expr_sub  <- expr[, top_genes, drop = FALSE]

  # NMF necesita valores no negativos
  expr_nn <- expr_sub - min(expr_sub)

  # NMF espera genes × muestras, por lo que debemos trasponer transponer
  # He usado lee porque es más rápido que brunet.
  
  res <- NMF::nmf(t(expr_nn), rank = k, nrun = nrun, seed = seed,
                   method = "lee", .options = "-v")

  # Asignación de un metagén de mayor peso para cada muestra
  cl_vec <- as.integer(NMF::predict(res, what = "samples"))
  cl     <- setNames(normalize_clusters(cl_vec), rownames(expr))

  cat(sprintf("    NMF (k=%d, nrun=%d, genes=%d): tabla = %s\n", k, nrun, n_genes,
              paste(names(table(cl)), table(cl), sep = ":", collapse = " | ")))
  return(cl)
}

# 3.7. GMM (Gaussian Mixture Models)

cluster_gmm <- function(expr, k, seed = 42, n_pcs = 30) {
  set.seed(seed)

  # Reducción dimensional previa (igual que antes, expr ya está centrado/escalado)
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

# 3.8. MCL (Markov Clustering)
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


# 4. Ejecución de todos los métodos

run_all_clustering <- function(expr, k, seed = 42,
                                cc_max_k = 6, cc_reps = 100,
                                cc_p_item = 0.8, cc_p_feature = 1.0,
                                cc_alg = "hc", cc_dist = "pearson",
                                graph_k_nn = 10, graph_res = 1.0,
                                base_dir = ".") {
  results <- list()

  
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

  cat("  → Leiden\n")
  results$leiden       <- tryCatch(
    cluster_leiden(expr, k, seed, graph_k_nn, graph_res),
    error = function(e) { warning("Leiden falló: ", e$message); NULL })

  cat("  → MCL\n")
  results$mcl          <- tryCatch(
    cluster_mcl(expr, k, seed, graph_k_nn), error = function(e) {
      warning("MCL falló: ", e$message); NULL })

  cat("  → NMF (puede tardar unos minutos)\n")
  results$nmf          <- tryCatch(
    cluster_nmf(expr, k, seed), error = function(e) {
      warning("NMF falló: ", e$message); NULL })

  cat("  → GMM\n")
  results$gmm          <- tryCatch(
    cluster_gmm(expr, k, seed), error = function(e) {
      warning("GMM falló: ", e$message); NULL })

  # Eliminación de métodos fallidos
  results <- Filter(Negate(is.null), results)

  cat(sprintf("  Métodos ejecutados con éxito: %d/%d\n",
              length(results), 8))
  return(results)
}


