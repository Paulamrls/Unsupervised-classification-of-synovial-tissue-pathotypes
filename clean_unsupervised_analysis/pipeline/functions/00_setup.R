# ##############################################################################
# functions/00_setup.R
# Instalación de paquetes, carga de librerías y creación de directorios
################################################################################

# 0. Paquetes necesarios

REQUIRED_CRAN <- c(
  "cluster",      
  "factoextra",    
  "ggplot2",       
  "ggrepel",       
  "pheatmap",      
  "networkD3",     
  "igraph",       
  "MCL",           
  "dendextend",    
  "RColorBrewer",  
  "dplyr",        
  "tidyr",         
  "tibble",        
  "scales",        
  "patchwork",     
  "ggalluvial",    
  "viridis",       
  "boot",          
  "mclust",        
  "kernlab",       
  "diffusionMap", 
  "reshape2",      
  "NMF"          
)

REQUIRED_BIOC <- c(
  "DESeq2",              
  "ConsensusClusterPlus",
  "limma",               
  "clusterProfiler",     
  "org.Hs.eg.db",        
  "fgsea",               
  "BiocParallel"         
)

# 1.Función principal de setup 
setup_environment <- function(install_missing = TRUE) {
  cat("Configurando entorno...\n")

  # Evita fallos no interactivos de install.packages() cuando R no tiene mirror.
  if (is.null(getOption("repos")) || identical(unname(getOption("repos")["CRAN"]),
                                                "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }

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

  # Cargar paquetes
  pkgs_all <- c(REQUIRED_CRAN, REQUIRED_BIOC)
  invisible(lapply(pkgs_all, function(p) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }))

  cat("  Todos los paquetes cargados.\n")
}

# 2. Crear estructura de carpetas de salida 
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

# 3. Guardar figura automáticamente
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

# 4. Guardar tabla CSV con encabezado informativo
save_table <- function(df, filename, base_dir = ".", subdir = "tables") {
  path <- file.path(base_dir, "results", subdir, filename)
  write.csv(df, path, row.names = TRUE)
  cat("  Tabla guardada:", file.path(subdir, filename), "\n")
}

# 5. Construir paleta fija de colores para k clusters 

#  colores de PATHOTYPE_COLORS en orden; si k > longitud, extiende
build_cluster_palette <- function(cluster_labels, pathotype_colors) {
  # cluster_labels: vector nombrado id → nombre biológico
  n <- length(cluster_labels)
  base_palette <- pathotype_colors
  all_colors   <- unname(base_palette)

  # colores en el orden en que aparecen las etiquetas
  palette <- setNames(all_colors[seq_len(n)], cluster_labels)
  return(palette)
}

# 6. Mensaje de progreso con timestamp
log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}
