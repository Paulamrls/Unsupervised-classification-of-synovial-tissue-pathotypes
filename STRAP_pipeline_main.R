################################################################################
# PIPELINE DE CLUSTERING TRANSCRIPTÓMICO EN AR — COHORTE STRAP
################################################################################
# Autor  : Paula Morales Sandoica
# Tutor:Cankut ÇUBUK
# Fecha  : 2026-05-08
# Script : STRAP_pipeline_main.R
## MÁSTER EN BIOINFORMÁTICA--UNIVERSIDAD EUROPEA DE MADRID
################################################################################
#################################################################################


# 0. CONFIGURACIÓN GLOBAL. Se trata del único script donde se pueden modificar
  # parámetros. El tutor sugirió que así sería más cómodo ir probando distintas k 
  # sin necesidad de cambiar todos los módulos y que el pipeline fuera óptimo y 
  # totalmente reproducible.


## Parámetro central

k <- 5

## Modo de ejecución:
#   "all"         → todas las muestras, sin filtro de patotipo
#   "pathotypes"  → solo muestras con patotipo válido, contaminadas incluidas
#   "clean"       → solo muestras con patotipo válido, contaminadas eliminadas
RUN_MODE <- "pathotypes"

## Reproducibilidad

GLOBAL_SEED    <- 42
set.seed(GLOBAL_SEED)

# 1. Preprocesado
VAR_CUTOFF     <- 0.25   # eliminar el 25% de genes de menor varianza (también recomendado por el tutor)
KW_PVAL        <- 0.01   # umbral p-valor Kruskal-Wallis (~5000 genes objetivo)
VST_BLIND      <- TRUE   # VST ciega

# 2.  Parámetros de Consensus clustering
CC_MAX_K       <- max(k + 2, 6)   
CC_REPS        <- 20
CC_P_ITEM      <- 0.80            
CC_P_FEATURE   <- 1.00            
CC_ALG         <- "hc"            
CC_DIST        <- "pearson"       

# 3. Número de componentes de Diffusion Maps a calcular

DM_N_DIMS      <- 5     

# 4. Clustering basado en grafos
GRAPH_K_NN     <- 10    
GRAPH_RES      <- 1.0   

# 5. Validación bootstrap 
BOOT_N         <- 20

# 6. Visualización y anotación

N_TOP_MARKERS  <- 20    

# Fijación del orden canónico de patotipos

 ## Los nuevos subtipos descubiertos se añadirán al final de forma automática

PATHOTYPE_ORDER <- c("Fibroid", "Myeloid", "Lymphoid")

# Fijación de una paleta de colores consistente en todos los gráficos

# Si k > 3 se usan los colores adicionales en orden
PATHOTYPE_COLORS <- c(
  "Fibroid"   = "#E41A1C",   # fibroblastos sinoviales
  "Myeloid"   = "#377EB8",   # macrófagos/monocitos
  "Lymphoid"  = "#4DAF4A",   # linfocitos
  "IFN-high"  = "#984EA3",   # interferón tipo I
  "Endothelial" = "#1B9E77",   # endotelio
  "Subtype1"  = "#FF7F00",   # otros
  "Subtype2"  = "#A65628",   
  "Subtype3"  = "#F781BF",   
  "Subtype4"  = "#FFFF33",  
  "Subtype5"  = "#66C2A5",   
  "Unresolved"= "#999999"    
)

# Muestras con contaminación muscular confirmada (ver TUTOR_FOLLOWUP_REPORT.md)
CONTAMINATED_SAMPLES <- c(
  "STRAPPAT00113-baseline",
  "STRAPPAT00210-baseline",
  "STRAPPAT00199-baseline",
  "STRAPPAT00007-baseline",
  "STRAPPAT00051-baseline"
)

# RUTAS
BASE_DIR  <- paste0("results_k", k, "_", RUN_MODE)  # carpeta separada por modo y k
DATA_FILE <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes/strap_counts.RData"
META_FILE <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes/E-MTAB-13733.sdrf.txt"

  ## Columnas de los metadatos de que son de interés

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


################################################################################
# 1. CARGA DE FUNCIONES Y SETUP

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


################################################################################
# 2. PREPROCESADO

cat("\n=== [1/8] PREPROCESADO ===\n")
cat(sprintf("  Modo de ejecucion: %s\n", RUN_MODE))

# Carga de datos brutos
expr_raw <- load_counts(DATA_FILE)

# Modo "clean": eliminar muestras contaminadas antes del análisis
if (RUN_MODE == "clean") {
  expr_raw <- remove_samples(expr_raw, CONTAMINATED_SAMPLES)
}

# Carga y limpieza de metadatos
meta <- load_metadata(META_FILE, META_COLS,
                      filter_pathotype = (RUN_MODE != "all"))

# Filtrado por varianza, Kruskal-Wallis y normalización VST
expr_vst <- preprocess_expression(expr_raw, meta = meta,
                                  var_quantile_cutoff = VAR_CUTOFF,
                                  n_top = 1000,
                                  kw_pval = KW_PVAL,
                                  blind = VST_BLIND, seed = GLOBAL_SEED)

# Alineación de muestras entre expresión y metadatos
aligned <- align_samples(expr_vst, meta)
expr_scaled  <- aligned$expr
meta_clean   <- aligned$meta

cat(sprintf("  Muestras finales: %d | Genes: %d\n",
            nrow(expr_scaled), ncol(expr_scaled)))


################################################################################
# 3. REDUCCIÓN DIMENSIONAL


cat("\n=== [2/8] REDUCCIÓN DIMENSIONAL ===\n")

# 3.1.PCA
pca_res <- run_pca(expr_scaled, seed = GLOBAL_SEED)

# 3.2.Diffusion Maps 
dm_res  <- run_diffusion_maps(expr_scaled, n_dims = DM_N_DIMS,
                               seed = GLOBAL_SEED)

cat("  PCA y Diffusion Maps calculados.\n")


################################################################################
# 4. CLUSTERING


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


cat("  Métodos ejecutados:", paste(names(all_clusters), collapse = ", "), "\n")


################################################################################
# 5. VALIDACIÓN


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

best_method <- validation_results$ranking$method[1]
cat(sprintf("  Mejor método (silhouette): %s\n", best_method))

final_clusters <- all_clusters[[best_method]]


################################################################################
# 6. ANÁLISIS DE EXPRESIÓN DIFERENCIAL


cat("\n=== [5/8] GENES DIFERENCIALES ===\n")

# Etiquetas numéricas temporales porque la anotación biológica no existe todavía
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

################################################################################
# 7. ANOTACIÓN BIOLÓGICA


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

# cluster_labels
cluster_labels <- annotation_res$labels
cluster_colors <- annotation_res$colors
cat("  Anotación:", paste(names(cluster_labels), "→", cluster_labels,
                           collapse = " | "), "\n")

# Reetiquetado de los DEGs con los nombres biológicos definitivos
deg_results <- deg_results_raw
names(deg_results$up)          <- cluster_labels[names(deg_results$up)]
names(deg_results$down)        <- cluster_labels[names(deg_results$down)]
names(deg_results$all)         <- cluster_labels[names(deg_results$all)]
names(deg_results$top_markers) <- cluster_labels[names(deg_results$top_markers)]


################################################################################
# 8. VISUALIZACIÓN


cat("\n=== [7/8] VISUALIZACIÓN ===\n")  
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


################################################################################
# 9. EXPORTACIÓN


cat("\n=== [8/8] EXPORTACIÓN ===\n")  
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
