# ==============================================================================
# PIPELINE DE CLUSTERING TRANSCRIPTÓMICO EN AR — COHORTE STRAP
# ==============================================================================
# Autor  : Paula Morales Sandoica
# Fecha  : 2026-05-08
# Script : STRAP_pipeline_main.R
#
# Objetivo: Descubrimiento no supervisado de patotipos transcriptómicos en
#           tejido sinovial de artritis reumatoide (cohorte STRAP).
#
# Uso: Modificar los parámetros de la sección "CONFIGURACIÓN GLOBAL" y
#      ejecutar este script completo. Todo el pipeline se adapta solo.
#
# Estructura de salida:
#   results/
#     figures/      → todos los gráficos
#     tables/       → todas las tablas CSV
#     clusters/     → asignaciones de cluster
#     validation/   → métricas de validación
#     enrichment/   → resultados de enriquecimiento
#     DEGs/         → genes diferenciales (up y down separados)
# ==============================================================================


# ==============================================================================
# 0. CONFIGURACIÓN GLOBAL ← ÚNICO LUGAR DONDE MODIFICAR PARÁMETROS
# ==============================================================================

# ── Parámetro central ─────────────────────────────────────────────────────────
k <- 5        # <── CAMBIA SOLO ESTE VALOR. Todo lo demás se adapta solo.

# ── Reproducibilidad ──────────────────────────────────────────────────────────
GLOBAL_SEED    <- 42
set.seed(GLOBAL_SEED)

# ── Preprocesado ──────────────────────────────────────────────────────────────
VAR_CUTOFF     <- 0.25   # eliminar el 25% de genes de menor varianza (recomendacion del tutor)
VST_BLIND      <- TRUE   # VST ciega (recomendado para clustering)

# ── Consensus clustering ──────────────────────────────────────────────────────
CC_MAX_K       <- max(k + 2, 6)   # explora hasta k+2 (mínimo 6)
CC_REPS        <- 100             # iteraciones de remuestreo
CC_P_ITEM      <- 0.80            # fracción de muestras por iteración
CC_P_FEATURE   <- 1.00            # fracción de genes por iteración
CC_ALG         <- "hc"            # algoritmo interno: "hc", "km" o "pam"
CC_DIST        <- "pearson"       # métrica de distancia

# ── Diffusion Maps ────────────────────────────────────────────────────────────
DM_N_DIMS      <- 5     # número de componentes de difusión a calcular

# ── Clustering basado en grafos ───────────────────────────────────────────────
GRAPH_K_NN     <- 10    # vecinos más cercanos para construir el grafo
GRAPH_RES      <- 1.0   # resolución Leiden/Louvain (>1 → más clusters)

# ── Validación bootstrap ──────────────────────────────────────────────────────
BOOT_N         <- 100   # iteraciones bootstrap para estabilidad

# ── Visualización y anotación ─────────────────────────────────────────────────
N_TOP_MARKERS  <- 20    # marcadores a mostrar por cluster

# ── Orden canónico de patotipos (fijoado biológicamente) ─────────────────────
# Los nuevos subtipos descubiertos se añaden al final automáticamente
PATHOTYPE_ORDER <- c("Fibroid", "Myeloid", "Lymphoid")

# ── Paleta de colores consistente en TODOS los gráficos ─────────────────────
# Orden alineado con canonical_order en 05_annotation.R:
#   Fibroid → Myeloid → Lymphoid → IFN-high → Vascular → Subtype1…
# Si k > 3 se usan los colores adicionales en orden
PATHOTYPE_COLORS <- c(
  "Fibroid"   = "#E41A1C",   # rojo        — fibroblastos sinoviales
  "Myeloid"   = "#377EB8",   # azul        — macrófagos/monocitos
  "Lymphoid"  = "#4DAF4A",   # verde       — linfocitos
  "IFN-high"  = "#984EA3",   # morado      — interferón tipo I
  "Vascular"  = "#1B9E77",   # verde azulado — endotelio / angiogénesis
  "Subtype1"  = "#FF7F00",   # naranja
  "Subtype2"  = "#A65628",   # marrón
  "Subtype3"  = "#F781BF",   # rosa
  "Subtype4"  = "#FFFF33",   # amarillo
  "Subtype5"  = "#66C2A5",   # turquesa
  "Unresolved"= "#999999"    # gris
)

# ── Rutas de datos ────────────────────────────────────────────────────────────
BASE_DIR  <- paste0("results_k", k)   # cada k guarda en su propia carpeta
DATA_FILE <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes/strap_counts.RData"
META_FILE <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes/E-MTAB-13733.sdrf.txt"

# ── Columnas de metadatos de interés ─────────────────────────────────────────
META_COLS <- c(
  "Characteristics.sampleid.",
  "Characteristics.pathotype.",
  "Characteristics.tender_joints_counts_tjc.",
  "Characteristics.swollen_joints_counts_sjc.",
  "Characteristics.arthritis.activity.",
  "Characteristics.esr.",
  "Characteristics.crp.",
  "Characteristics.physicians.global.assessment."
)


# ==============================================================================
# 1. CARGA DE FUNCIONES Y SETUP
# ==============================================================================

source("functions/00_setup.R")
source("functions/01_preprocessing.R")
source("functions/02_dimreduction.R")
source("functions/03_clustering.R")
source("functions/04_validation.R")
source("functions/05_deg_analysis.R")
source("functions/06_annotation.R")
source("functions/07_visualization.R")
source("functions/08_export.R")

# Cargar librerías y crear estructura de carpetas
setup_environment()
create_output_dirs(BASE_DIR)


# ==============================================================================
# 2. PREPROCESADO
# ==============================================================================

cat("\n=== [1/8] PREPROCESADO ===\n")

# Carga de datos brutos
expr_raw <- load_counts(DATA_FILE)

# Carga y limpieza de metadatos
meta <- load_metadata(META_FILE, META_COLS)

# Filtrado por varianza y normalización VST
expr_vst <- preprocess_expression(expr_raw, var_quantile_cutoff = VAR_CUTOFF,
                                  blind = VST_BLIND, seed = GLOBAL_SEED)

# Alinear muestras entre expresión y metadatos
aligned <- align_samples(expr_vst, meta)
expr_scaled  <- aligned$expr    # muestras × genes, escalado
meta_clean   <- aligned$meta    # metadatos limpios, alineados

cat(sprintf("  Muestras finales: %d | Genes: %d\n",
            nrow(expr_scaled), ncol(expr_scaled)))


# ==============================================================================
# 3. REDUCCIÓN DIMENSIONAL
# ==============================================================================

cat("\n=== [2/8] REDUCCIÓN DIMENSIONAL ===\n")

# PCA
pca_res <- run_pca(expr_scaled, seed = GLOBAL_SEED)

# Diffusion Maps (captura trayectorias biológicas continuas)
dm_res  <- run_diffusion_maps(expr_scaled, n_dims = DM_N_DIMS,
                               seed = GLOBAL_SEED)

cat("  PCA y Diffusion Maps calculados.\n")


# ==============================================================================
# 4. CLUSTERING
# ==============================================================================

cat("\n=== [3/8] CLUSTERING (k =", k, ") ===\n")

all_clusters <- run_all_clustering(
  expr         = expr_scaled,
  k            = k,
  seed         = GLOBAL_SEED,
  cc_max_k     = CC_MAX_K,
  cc_reps      = CC_REPS,
  cc_p_item    = CC_P_ITEM,
  cc_p_feature = CC_P_FEATURE,
  cc_alg       = CC_ALG,
  cc_dist      = CC_DIST,
  graph_k_nn   = GRAPH_K_NN,
  graph_res    = GRAPH_RES,
  base_dir     = BASE_DIR
)

# all_clusters es una lista con: kmeans, hierarchical, spectral,
#   consensus, leiden, louvain, infomap, mcl
cat("  Métodos ejecutados:", paste(names(all_clusters), collapse = ", "), "\n")


# ==============================================================================
# 5. VALIDACIÓN
# ==============================================================================

cat("\n=== [4/8] VALIDACIÓN ===\n")

validation_results <- run_validation(
  expr         = expr_scaled,
  clusters     = all_clusters,
  meta         = meta_clean,
  k            = k,
  boot_n       = BOOT_N,
  seed         = GLOBAL_SEED,
  pathotype_col = "pathotype_clean"
)

# Identificar el mejor método según silhouette
best_method <- validation_results$ranking$method[1]
cat(sprintf("  Mejor método (silhouette): %s\n", best_method))

# Usar el mejor método como clusters definitivos
final_clusters <- all_clusters[[best_method]]


# ==============================================================================
# 6. ANÁLISIS DE EXPRESIÓN DIFERENCIAL
# ==============================================================================
# Los DEGs se calculan ANTES de la anotación: los genes marcadores de cada
# cluster informan directamente qué nombre biológico asignarle.
# En este paso se usan IDs numéricos; se reetiquetan tras la anotación.
# ==============================================================================

cat("\n=== [5/8] GENES DIFERENCIALES ===\n")

# Etiquetas numéricas temporales (la anotación biológica aún no existe)
numeric_labels <- setNames(
  as.character(sort(unique(final_clusters))),
  as.character(sort(unique(final_clusters)))
)

deg_results_raw <- run_deg_analysis(
  expr           = expr_scaled,
  expr_raw       = expr_raw,
  clusters       = final_clusters,
  cluster_labels = numeric_labels,
  n_markers      = N_TOP_MARKERS,
  seed           = GLOBAL_SEED
)


# ==============================================================================
# 7. ANOTACIÓN BIOLÓGICA
# ==============================================================================
# La anotación usa las firmas de referencia + el solapamiento con patotipos
# histológicos. Los DEGs calculados en el paso anterior respaldan y validan
# el nombre asignado a cada cluster.
# ==============================================================================

cat("\n=== [6/8] ANOTACIÓN BIOLÓGICA ===\n")

annotation_res <- annotate_clusters(
  clusters         = final_clusters,
  expr             = expr_scaled,
  meta             = meta_clean,
  k                = k,
  pathotype_order  = PATHOTYPE_ORDER,
  pathotype_colors = PATHOTYPE_COLORS,
  n_markers        = N_TOP_MARKERS,
  seed             = GLOBAL_SEED
)

# cluster_labels: vector nombrado cluster_id → nombre biológico
cluster_labels <- annotation_res$labels
cluster_colors <- annotation_res$colors
cat("  Anotación:", paste(names(cluster_labels), "→", cluster_labels,
                           collapse = " | "), "\n")

# Reetiqueta los DEGs con los nombres biológicos definitivos
deg_results <- deg_results_raw
names(deg_results$up)          <- cluster_labels[names(deg_results$up)]
names(deg_results$down)        <- cluster_labels[names(deg_results$down)]
names(deg_results$all)         <- cluster_labels[names(deg_results$all)]
names(deg_results$top_markers) <- cluster_labels[names(deg_results$top_markers)]


# ==============================================================================
# 8. VISUALIZACIÓN
# ==============================================================================

cat("\n=== [7/8] VISUALIZACIÓN ===\n")  # 7: con clusters anotados + DEGs disponibles

generate_all_plots(
  pca_res        = pca_res,
  dm_res         = dm_res,
  all_clusters   = all_clusters,
  final_clusters = final_clusters,
  cluster_labels = cluster_labels,
  cluster_colors = cluster_colors,
  validation     = validation_results,
  deg_results    = deg_results,
  meta           = meta_clean,
  expr           = expr_scaled,
  k              = k,
  pathotype_order = PATHOTYPE_ORDER,
  pathotype_col  = "pathotype_clean",
  base_dir       = BASE_DIR
)


# ==============================================================================
# 9. EXPORTACIÓN
# ==============================================================================

cat("\n=== [8/8] EXPORTACIÓN ===\n")  # 8: exporta todo una vez el pipeline está completo

export_all_results(
  all_clusters   = all_clusters,
  final_clusters = final_clusters,
  cluster_labels = cluster_labels,
  validation     = validation_results,
  deg_results    = deg_results,
  meta           = meta_clean,
  k              = k,
  base_dir       = BASE_DIR
)

cat("\n✓ Pipeline completado. Resultados en:", file.path(BASE_DIR, "results"), "\n")
cat("  k =", k, "| Mejor método:", best_method, "\n")
cat("  Patotipos:", paste(cluster_labels, collapse = " | "), "\n\n")
