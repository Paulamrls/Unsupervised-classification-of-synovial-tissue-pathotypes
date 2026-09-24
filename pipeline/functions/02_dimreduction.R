# ==============================================================================
# functions/02_dimreduction.R
# PCA y Diffusion Maps (sin UMAP)
# ==============================================================================

# ── PCA ───────────────────────────────────────────────────────────────────────
run_pca <- function(expr_scaled, seed = 42) {
  set.seed(seed)
  pca <- prcomp(expr_scaled, center = FALSE, scale. = FALSE)
  # (expr_scaled ya está centrado/escalado en preprocesado)

  var_explained <- summary(pca)$importance[2, ] * 100  # % varianza

  cat(sprintf("  PCA: PC1=%.1f%% | PC2=%.1f%% | PC3=%.1f%%\n",
              var_explained[1], var_explained[2], var_explained[3]))

  return(list(
    pca          = pca,
    coords       = as.data.frame(pca$x),
    var_explained = var_explained
  ))
}

# ── Diffusion Maps ────────────────────────────────────────────────────────────
# Captura trayectorias biológicas continuas y heterogeneidad no lineal.
# Especialmente útil para detectar estados intermedios/híbridos en AR.
# Usa el paquete diffusionMap (compatible con R 4.5, a diferencia de destiny).
run_diffusion_maps <- function(expr_scaled, n_dims = 5, seed = 42) {
  set.seed(seed)

  # Calcular matriz de distancias euclídeas entre muestras
  dist_mat <- as.matrix(dist(expr_scaled))

  # Estimar epsilon óptimo para el kernel de difusión
  eps_val <- diffusionMap::epsilonCompute(dist_mat)

  # Calcular Diffusion Map
  dm <- diffusionMap::diffuse(dist_mat,
                               eps.val = eps_val,
                               t       = 0,
                               maxdim  = n_dims,
                               delta   = 1e-10)

  # Extraer coordenadas (columnas = componentes de difusión)
  dm_coords <- as.data.frame(dm$X[, seq_len(n_dims), drop = FALSE])
  colnames(dm_coords) <- paste0("DC", seq_len(n_dims))
  rownames(dm_coords) <- rownames(expr_scaled)

  cat(sprintf("  Diffusion Maps: %d componentes calculados.\n", n_dims))

  return(list(
    dm     = dm,
    coords = dm_coords
  ))
}

# ── Varianza explicada por PCA (para elegir n. componentes óptimo) ────────────
pca_elbow_data <- function(pca_res, n_pcs = 20) {
  var_exp <- pca_res$var_explained[seq_len(min(n_pcs, length(pca_res$var_explained)))]
  data.frame(
    PC      = seq_along(var_exp),
    var_pct = var_exp,
    cumvar  = cumsum(var_exp)
  )
}
