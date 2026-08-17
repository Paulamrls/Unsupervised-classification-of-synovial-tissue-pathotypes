################################################################################
# PIPELINE DE CLUSTERING TRANSCRIPTOMICO NO SUPERVISADO - COHORTE STRAP
################################################################################

# Configuracion comun. La unica diferencia entre las dos ejecuciones es k.
K_VALUES <- c(4L, 5L)
GLOBAL_SEED <- 42L
GENE_SELECTION_METHOD <- "top_variable"
N_TOP_GENES <- 3000L
MIN_COUNT <- 10L
MIN_SAMPLE_PROP <- 0.10
N_CLUSTERING_PCS <- 30L
FINAL_CLUSTER_METHOD <- "kmeans"
VAR_CUTOFF <- 0.25
KW_PVAL <- 0.01
VST_BLIND <- TRUE

CC_REPS <- 20L
CC_P_ITEM <- 0.80
CC_P_FEATURE <- 1.00
CC_ALG <- "hc"
CC_DIST <- "pearson"
DM_N_DIMS <- 5L
GRAPH_K_NN <- 10L
GRAPH_RES <- 1.0
BOOT_N <- 20L
N_TOP_MARKERS <- 20L

PATHOTYPE_ORDER <- c("Fibroid", "Myeloid", "Lymphoid")
PATHOTYPE_COLORS <- c(
  Fibroid = "#E41A1C", Myeloid = "#377EB8", Lymphoid = "#4DAF4A",
  `IFN-high` = "#984EA3", Endothelial = "#1B9E77", Subtype1 = "#FF7F00",
  Subtype2 = "#A65628", Subtype3 = "#F781BF", Subtype4 = "#FFFF33",
  Subtype5 = "#66C2A5", Unresolved = "#999999"
)

CONTAMINATED_SAMPLES <- c(
  "STRAPPAT00113-baseline", "STRAPPAT00210-baseline",
  "STRAPPAT00199-baseline", "STRAPPAT00007-baseline",
  "STRAPPAT00051-baseline"
)

PROJECT_DIR <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes"
DATA_FILE <- file.path(PROJECT_DIR, "strap_counts.RData")
META_FILE <- file.path(PROJECT_DIR, "E-MTAB-13733.sdrf.txt")
OUTPUT_ROOT <- Sys.getenv(
  "STRAP_OUTPUT_ROOT",
  unset = file.path(PROJECT_DIR, "results_clean_unsupervised_relaxed_p10_g3000_pc30")
)

META_COLS <- c(
  "Characteristics.sampleid.", "Characteristics.pathotype.",
  "Characteristics.tender_joints_counts_tjc.",
  "Characteristics.swollen_joints_counts_sjc.",
  "Characteristics.arthritis.activity.", "Characteristics.esr.",
  "Characteristics.crp.",
  "Characteristics.physicians.global.assessment."
)

source("functions/00_setup.R")
source("functions/01_preprocessing.R")
source("functions/02_dimreduction.R")
source("functions/03_clustering.R")
source("functions/04_validation.R")
source("functions/05_deg_analysis.R")
source("functions/06_annotation.R")
source("functions/07_visualization.R")
source("functions/08_export.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop("COMPROBACION FALLIDA: ", message, call. = FALSE)
}

assert_no_contaminated <- function(sample_ids, stage) {
  found <- intersect(sample_ids, CONTAMINATED_SAMPLES)
  assert_true(length(found) == 0L,
              paste0(stage, ": siguen presentes: ", paste(found, collapse = ", ")))
}

write_csv_plain <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
}

run_one_k <- function(k, cluster_input, expr_scaled, expr_raw, meta_clean,
                      pca_res, dm_res) {
  set.seed(GLOBAL_SEED)
  base_dir <- file.path(OUTPUT_ROOT, paste0("k", k))
  create_output_dirs(base_dir)
  cc_max_k <- max(k + 2L, 6L)

  cat(sprintf("\n=== ANALISIS k = %d ===\n", k))
  all_clusters <- run_all_clustering(
    expr = cluster_input, k = k, seed = GLOBAL_SEED,
    cc_max_k = cc_max_k, cc_reps = CC_REPS,
    cc_p_item = CC_P_ITEM, cc_p_feature = CC_P_FEATURE,
    cc_alg = CC_ALG, cc_dist = CC_DIST,
    graph_k_nn = GRAPH_K_NN, graph_res = GRAPH_RES,
    base_dir = base_dir
  )

  validation <- run_validation(
    expr = cluster_input, clusters = all_clusters, meta = meta_clean,
    k = k, boot_n = BOOT_N, seed = GLOBAL_SEED,
    pathotype_col = "pathotype_clean"
  )
  assert_true(FINAL_CLUSTER_METHOD %in% names(all_clusters),
              paste0("el metodo final no esta disponible: ", FINAL_CLUSTER_METHOD))
  validation$selected_method <- FINAL_CLUSTER_METHOD
  best_method <- FINAL_CLUSTER_METHOD
  final_clusters <- all_clusters[[FINAL_CLUSTER_METHOD]]

  assert_true(!anyDuplicated(names(final_clusters)),
              paste0("k=", k, ": hay muestras duplicadas"))
  assert_no_contaminated(names(final_clusters), paste0("clustering k=", k))
  assert_true(length(final_clusters) == sum(table(final_clusters)),
              paste0("k=", k, ": n muestras no coincide con suma de clusters"))
  assert_true(setequal(names(final_clusters), rownames(expr_scaled)),
              paste0("k=", k, ": las muestras difieren de la matriz comun"))
  assert_true(identical(rownames(cluster_input), rownames(expr_scaled)),
              paste0("k=", k, ": los PCs no coinciden con las muestras comunes"))

  numeric_labels <- setNames(as.character(sort(unique(final_clusters))),
                             as.character(sort(unique(final_clusters))))
  deg_raw <- run_deg_analysis(
    expr = expr_scaled, expr_raw = expr_raw, clusters = final_clusters,
    cluster_labels = numeric_labels, n_markers = N_TOP_MARKERS,
    seed = GLOBAL_SEED
  )
  annotation <- annotate_clusters(
    clusters = final_clusters, expr = expr_scaled, meta = meta_clean, k = k,
    pathotype_order = PATHOTYPE_ORDER, pathotype_colors = PATHOTYPE_COLORS,
    n_markers = N_TOP_MARKERS, seed = GLOBAL_SEED
  )
  labels <- annotation$labels
  deg_results <- deg_raw
  names(deg_results$up) <- labels[names(deg_results$up)]
  names(deg_results$down) <- labels[names(deg_results$down)]
  names(deg_results$all) <- labels[names(deg_results$all)]
  names(deg_results$top_markers) <- labels[names(deg_results$top_markers)]

  generate_all_plots(
    pca_res = pca_res, dm_res = dm_res, all_clusters = all_clusters,
    final_clusters = final_clusters, cluster_labels = labels,
    cluster_colors = annotation$colors, validation = validation,
    deg_results = deg_results, meta = meta_clean, expr = expr_scaled, k = k,
    pathotype_order = PATHOTYPE_ORDER, pathotype_col = "pathotype_clean",
    base_dir = base_dir
  )
  export_all_results(
    all_clusters = all_clusters, final_clusters = final_clusters,
    cluster_labels = labels, validation = validation, deg_results = deg_results,
    meta = meta_clean, k = k, pca_res = pca_res, base_dir = base_dir
  )

  write_csv_plain(
    data.frame(cluster_id = names(table(final_clusters)),
               n_samples = as.integer(table(final_clusters))),
    file.path(base_dir, "cluster_sizes.csv")
  )
  write_csv_plain(validation$ranking[, c("method", "silhouette")],
                  file.path(base_dir, "silhouette.csv"))
  write_csv_plain(validation$ranking[, c("method", "boot_ari", "boot_sd")],
                  file.path(base_dir, "bootstrap_ari.csv"))
  write_csv_plain(
    data.frame(sample_id = names(final_clusters), cluster_id = unname(final_clusters)),
    file.path(base_dir, "cluster_assignments.csv")
  )

  list(k = k, samples = names(final_clusters), genes = colnames(expr_scaled),
       clustering_features = colnames(cluster_input),
       clusters = final_clusters, best_method = best_method,
       validation = validation)
}

set.seed(GLOBAL_SEED)
setup_environment()
dir.create(file.path(OUTPUT_ROOT, "common"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_ROOT, "comparison_k4_k5"), recursive = TRUE,
           showWarnings = FALSE)

cat("\n=== PREPROCESADO COMUN NO SUPERVISADO ===\n")
expr_raw <- load_counts(DATA_FILE)
assert_true(!anyDuplicated(colnames(expr_raw)), "counts contiene muestras duplicadas")

# Exclusion obligatoria antes de cualquier transformacion o analisis.
present_before <- intersect(CONTAMINATED_SAMPLES, colnames(expr_raw))
assert_true(setequal(present_before, CONTAMINATED_SAMPLES),
            "no se localizaron exactamente las cinco muestras contaminadas en counts")
expr_raw <- remove_samples(expr_raw, CONTAMINATED_SAMPLES)
assert_no_contaminated(colnames(expr_raw), "counts tras exclusion")

# No se filtra ninguna muestra por su histologia. Los metadatos se usan despues.
meta <- load_metadata(META_FILE, META_COLS, filter_pathotype = FALSE)

expr_scaled <- preprocess_expression(
  expr = expr_raw, meta = meta, method = GENE_SELECTION_METHOD,
  min_count = MIN_COUNT,
  min_samples = ceiling(MIN_SAMPLE_PROP * ncol(expr_raw)),
  var_quantile_cutoff = VAR_CUTOFF, n_top = N_TOP_GENES,
  kw_pval = KW_PVAL, blind = VST_BLIND, seed = GLOBAL_SEED
)
selected_genes <- attr(expr_scaled, "selected_genes")
assert_true(GENE_SELECTION_METHOD == "top_variable",
            "el analisis principal debe usar top_variable")
assert_true(length(selected_genes) == N_TOP_GENES,
            "no se seleccionaron exactamente 3.000 genes")
assert_true(!anyDuplicated(selected_genes), "hay genes seleccionados duplicados")
assert_no_contaminated(rownames(expr_scaled), "matriz VST")

aligned <- align_samples(expr_scaled, meta)
expr_scaled <- aligned$expr
meta_clean <- aligned$meta
assert_true(!anyDuplicated(rownames(expr_scaled)), "hay muestras duplicadas tras alinear")
assert_no_contaminated(rownames(expr_scaled), "datos comunes")

write_csv_plain(data.frame(gene = selected_genes),
                file.path(OUTPUT_ROOT, "common", "selected_genes_3000.csv"))
write_csv_plain(data.frame(sample_id = rownames(expr_scaled)),
                file.path(OUTPUT_ROOT, "common", "samples_used.csv"))

pca_res <- run_pca(expr_scaled, seed = GLOBAL_SEED)
dm_res <- run_diffusion_maps(expr_scaled, n_dims = DM_N_DIMS, seed = GLOBAL_SEED)
n_pcs_use <- min(N_CLUSTERING_PCS, ncol(pca_res$coords))
cluster_input <- as.matrix(pca_res$coords[, seq_len(n_pcs_use), drop = FALSE])
assert_true(n_pcs_use == N_CLUSTERING_PCS,
            "no fue posible obtener exactamente 30 PCs para clustering")
assert_true(identical(rownames(cluster_input), rownames(expr_scaled)),
            "orden de muestras distinto entre VST y PCs")
write_csv_plain(data.frame(
  parameter = c("gene_selection_method", "min_count", "min_sample_prop",
                "min_samples", "n_top_genes", "n_clustering_pcs",
                "final_cluster_method", "seed"),
  value = c(GENE_SELECTION_METHOD, MIN_COUNT, MIN_SAMPLE_PROP,
            ceiling(MIN_SAMPLE_PROP * ncol(expr_raw)), N_TOP_GENES,
            N_CLUSTERING_PCS, FINAL_CLUSTER_METHOD, GLOBAL_SEED)
), file.path(OUTPUT_ROOT, "common", "common_parameters.csv"))

runs <- lapply(K_VALUES, run_one_k, cluster_input = cluster_input,
               expr_scaled = expr_scaled, expr_raw = expr_raw,
               meta_clean = meta_clean, pca_res = pca_res, dm_res = dm_res)
names(runs) <- paste0("k", K_VALUES)

assert_true(identical(runs$k4$samples, runs$k5$samples),
            "k=4 y k=5 no usan las mismas muestras y el mismo orden")
assert_true(identical(runs$k4$genes, runs$k5$genes),
            "k=4 y k=5 no usan los mismos genes y el mismo orden")
assert_true(identical(runs$k4$genes, selected_genes),
            "los genes usados no coinciden con la lista exportada")
assert_true(identical(runs$k4$clustering_features,
                      runs$k5$clustering_features),
            "k=4 y k=5 no usan exactamente los mismos PCs")

comparison <- data.frame(
  sample_id = runs$k4$samples,
  cluster_k4 = unname(runs$k4$clusters[runs$k4$samples]),
  cluster_k5 = unname(runs$k5$clusters[runs$k4$samples])
)
write_csv_plain(comparison,
                file.path(OUTPUT_ROOT, "comparison_k4_k5", "sample_comparison_k4_k5.csv"))
contingency <- as.data.frame.matrix(table(k4 = comparison$cluster_k4,
                                         k5 = comparison$cluster_k5))
contingency$cluster_k4 <- rownames(contingency)
contingency <- contingency[, c("cluster_k4", setdiff(names(contingency), "cluster_k4"))]
write_csv_plain(contingency,
                file.path(OUTPUT_ROOT, "comparison_k4_k5", "contingency_k4_k5.csv"))

checks <- data.frame(
  check = c("contaminated_removed", "same_samples_k4_k5", "same_3000_genes_k4_k5",
            "same_30_pcs_k4_k5", "histology_not_used_for_gene_selection",
            "no_duplicate_samples",
            "cluster_sizes_sum_to_n_k4", "cluster_sizes_sum_to_n_k5"),
  passed = TRUE
)
write_csv_plain(checks, file.path(OUTPUT_ROOT, "common", "pipeline_checks.csv"))
saveRDS(runs, file.path(OUTPUT_ROOT, "comparison_k4_k5", "runs_k4_k5.rds"))

cat("\nPipeline completado: ", OUTPUT_ROOT, "\n", sep = "")
