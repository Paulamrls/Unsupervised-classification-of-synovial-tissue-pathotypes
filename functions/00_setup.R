# ==============================================================================
# functions/00_setup.R
# Instalación de paquetes, carga de librerías y creación de directorios
# ==============================================================================

# Lista completa de paquetes requeridos
REQUIRED_CRAN <- c(
  "cluster",       # silhouette, pam
  "factoextra",    # fviz_nbclust, fviz_pca
  "ggplot2",       # gráficos base
  "ggrepel",       # etiquetas sin solapamiento
  "pheatmap",      # heatmaps
  "networkD3",     # sankey diagram
  "igraph",        # grafos + Louvain, Infomap, Leiden
  "MCL",           # Markov Clustering
  "dendextend",    # dendrogramas coloreados
  "RColorBrewer",  # paletas de color
  "dplyr",         # manipulación de datos
  "tidyr",         # pivotado
  "tibble",        # tablas tidy
  "scales",        # formato de ejes
  "patchwork",     # combinar ggplots
  "ggalluvial",    # sankey en ggplot2
  "viridis",       # paletas perceptualmente uniformes
  "boot",          # bootstrap
  "mclust",        # adjustedRandIndex + GMM clustering
  "kernlab",       # Spectral clustering (specc)
  "diffusionMap",  # Diffusion Maps (alternativa a destiny, compatible con R 4.5)
  "reshape2",      # melt para boxplots de marcadores
  "NMF"            # Non-negative Matrix Factorization (bulk RNA-seq)
)

REQUIRED_BIOC <- c(
  "DESeq2",              # normalización VST + DEG
  "ConsensusClusterPlus",# consensus clustering reproducible
  "limma",               # DEG moderado (análisis de expresión diferencial)
  "clusterProfiler",     # enriquecimiento funcional GO/KEGG
  "org.Hs.eg.db",        # base de datos anotación humana
  "fgsea",               # GSEA rápido
  "BiocParallel"         # paralelización Bioc
)

# ── Función principal de setup ─────────────────────────────────────────────────
setup_environment <- function(install_missing = TRUE) {
  cat("Configurando entorno...\n")

  # Instalar BiocManager si no existe
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  if (install_missing) {
    # CRAN
    missing_cran <- REQUIRED_CRAN[!sapply(REQUIRED_CRAN, requireNamespace,
                                           quietly = TRUE)]
    if (length(missing_cran) > 0) {
      cat("  Instalando paquetes CRAN:", paste(missing_cran, collapse = ", "), "\n")
      install.packages(missing_cran, quiet = TRUE)
    }

    # Bioconductor
    missing_bioc <- REQUIRED_BIOC[!sapply(REQUIRED_BIOC, requireNamespace,
                                           quietly = TRUE)]
    if (length(missing_bioc) > 0) {
      cat("  Instalando paquetes Bioc:", paste(missing_bioc, collapse = ", "), "\n")
      BiocManager::install(missing_bioc, ask = FALSE, quiet = TRUE)
    }
  }

  # Cargar todos los paquetes
  pkgs_all <- c(REQUIRED_CRAN, REQUIRED_BIOC)
  invisible(lapply(pkgs_all, function(p) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }))

  cat("  Todos los paquetes cargados.\n")
}

# ── Crear estructura de carpetas de salida ─────────────────────────────────────
create_output_dirs <- function(base_dir = ".") {
  dirs <- c(
    file.path(base_dir, "results", "figures"),
    file.path(base_dir, "results", "tables"),
    file.path(base_dir, "results", "clusters"),
    file.path(base_dir, "results", "validation"),
    file.path(base_dir, "results", "enrichment"),
    file.path(base_dir, "results", "DEGs")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  cat("  Estructura de carpetas creada en:", file.path(base_dir, "results"), "\n")
}

# ── Guardar figura automáticamente ────────────────────────────────────────────
save_figure <- function(plot_obj = NULL, filename, base_dir = ".",
                         width = 10, height = 8, dpi = 300) {
  path <- file.path(base_dir, "results", "figures", filename)
  if (inherits(plot_obj, "ggplot")) {
    ggsave(path, plot = plot_obj, width = width, height = height,
           dpi = dpi, bg = "white")
  } else {
    # Para gráficos base R
    png(path, width = width * dpi, height = height * dpi,
        res = dpi, bg = "white")
    if (!is.null(plot_obj)) print(plot_obj)
    dev.off()
  }
  cat("  Figura guardada:", filename, "\n")
}

# ── Guardar tabla CSV con encabezado informativo ──────────────────────────────
save_table <- function(df, filename, base_dir = ".", subdir = "tables") {
  path <- file.path(base_dir, "results", subdir, filename)
  write.csv(df, path, row.names = TRUE)
  cat("  Tabla guardada:", file.path(subdir, filename), "\n")
}

# ── Construir paleta de colores para k clusters ──────────────────────────────
#   Usa los colores de PATHOTYPE_COLORS en orden; si k > longitud, extiende
build_cluster_palette <- function(cluster_labels, pathotype_colors) {
  # cluster_labels: vector nombrado id → nombre biológico
  n <- length(cluster_labels)
  base_palette <- pathotype_colors
  all_colors   <- unname(base_palette)

  # Asignar colores en el orden en que aparecen los labels
  palette <- setNames(all_colors[seq_len(n)], cluster_labels)
  return(palette)
}

# ── Mensaje de progreso con timestamp ────────────────────────────────────────
log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}
