################################################################################
# functions/02_dimreduction.R
# PCA y Diffusion Maps (sin UMAP)
################################################################################


# 1. PCA

run_pca <- function(expr_scaled, seed = 42) {
  set.seed(seed)
  pca <- prcomp(expr_scaled, center = FALSE, scale. = FALSE) # ya escalado en el módulo de preprocesamiento
 
  var_explained <- summary(pca)$importance[2, ] * 100  # varianza en %
  cat(sprintf("  PCA: PC1=%.1f%% | PC2=%.1f%% | PC3=%.1f%%\n",
              var_explained[1], var_explained[2], var_explained[3]))

  return(list(
    pca          = pca,
    coords       = as.data.frame(pca$x),
    var_explained = var_explained
  ))
}

# 2. Diffusion Maps

run_diffusion_maps <- function(expr_scaled, n_dims = 5, seed = 42) {
  set.seed(seed)
  dist_mat <- as.matrix(dist(expr_scaled))

  eps_val <- diffusionMap::epsilonCompute(dist_mat)

  dm <- diffusionMap::diffuse(dist_mat,
                               eps.val = eps_val,
                               t       = 0,
                               maxdim  = n_dims,
                               delta   = 1e-10)

  # Extracción de coordenadas
  dm_coords <- as.data.frame(dm$X[, seq_len(n_dims), drop = FALSE])
  colnames(dm_coords) <- paste0("DC", seq_len(n_dims))
  rownames(dm_coords) <- rownames(expr_scaled)

  cat(sprintf("  Diffusion Maps: %d componentes calculados.\n", n_dims))

  return(list(
    dm     = dm,
    coords = dm_coords
  ))
}

# 3.  Varianza explicada por PCA 

pca_elbow_data <- function(pca_res, n_pcs = 20) {
  var_exp <- pca_res$var_explained[seq_len(min(n_pcs, length(pca_res$var_explained)))]
  data.frame(
    PC      = seq_along(var_exp),
    var_pct = var_exp,
    cumvar  = cumsum(var_exp)
  )
}
