################################################################################
# functions/01_preprocessing.R
# Carga de datos, filtrado por varianza, normalización VST, alineación
################################################################################

# 1. Carga de counts crudos
load_counts <- function(data_file) {
  if (!file.exists(data_file))
    stop("No se encuentra el archivo de datos: ", data_file)

  env <- new.env()
  load(data_file, envir = env)

  # Búsqueda del primer objeto matrix/data.frame en el entorno cargado
  obj_name <- ls(env)[1]
  expr <- get(obj_name, envir = env)
  expr <- as.matrix(expr)

  # Aseguramos que existan nombres de fila y columna
  if (is.null(rownames(expr)))
    rownames(expr) <- paste0("Gene_", seq_len(nrow(expr)))
  if (is.null(colnames(expr)))
    colnames(expr) <- paste0("Sample_", seq_len(ncol(expr)))

  cat(sprintf("  Datos cargados: %d genes x %d muestras\n",
              nrow(expr), ncol(expr)))
  return(expr)
}

# 2. Carga y limpieza de metadatos
load_metadata <- function(meta_file, cols_keep, filter_pathotype = TRUE) {
  if (!file.exists(meta_file))
    stop("No se encuentra el archivo de metadatos: ", meta_file)

  meta <- read.delim(meta_file, stringsAsFactors = FALSE, check.names = TRUE)

  # Columnas disponibles para el diagnóstico
  cat("  Columnas disponibles en metadatos:\n")
  cat("   ", paste(colnames(meta)[1:min(10, ncol(meta))], collapse = "\n    "), "\n")

  # Normalización de los nombres de columna 
  actual_cols <- colnames(meta)

  # Detección columna de sampleid por coincidencia parcial
  id_col <- detect_column(actual_cols, c("sampleid", "sample.id", "SampleID",
                                          "Characteristics.sampleid."))
  if (is.null(id_col))
    stop("No se encontró columna de ID de muestra en metadatos. ",
         "Columnas disponibles: ", paste(actual_cols, collapse = ", "))

  # Detección de columna de patotipo por coincidencia parcial
  path_col <- detect_column(actual_cols, c("pathotype", "Pathotype",
                                            "Characteristics.pathotype."))

  # Selección de columnas existentes de cols_keep + las detectadas automáticamente
  cols_found <- intersect(cols_keep, actual_cols)
  extra_cols <- unique(c(id_col, path_col))
  extra_cols <- extra_cols[!is.na(extra_cols)]
  cols_final <- unique(c(cols_found, extra_cols))

  meta_sub <- meta[, cols_final, drop = FALSE]

  # Eliminación de duplicados por sampleid
  meta_unique <- meta_sub[!duplicated(meta_sub[[id_col]]), ]
  rownames(meta_unique) <- as.character(meta_unique[[id_col]])

  # Limpieza de la columna de patotipo
  if (!is.null(path_col) && path_col %in% colnames(meta_unique)) {
    meta_unique$pathotype_clean <- as.character(meta_unique[[path_col]])

    # Unificación de variantes de nombre
    meta_unique$pathotype_clean[
      meta_unique$pathotype_clean == "Fibrous"] <- "Fibroid"

    # Eliminación de muestras sin clasificación válida (solo si filter_pathotype = TRUE)
    if (filter_pathotype) {
      meta_unique <- meta_unique[
        !is.na(meta_unique$pathotype_clean) &
          meta_unique$pathotype_clean != "Ungraded", ]
    }

    cat(sprintf("  Patotipos encontrados: %s\n",
                paste(names(table(meta_unique$pathotype_clean)),
                      table(meta_unique$pathotype_clean),
                      sep = "=", collapse = " | ")))
  } else {
    warning("No se encontró columna de patotipo. La validacion biologica ",
            "no estara disponible.")
    meta_unique$pathotype_clean <- NA
  }

  cat(sprintf("  Metadatos cargados: %d muestras unicas\n", nrow(meta_unique)))
  return(meta_unique)
}

# 3. Detección de columna por coincidencia parcial

detect_column <- function(actual_cols, candidates) {
  # Coincidencia exacta
  exact <- intersect(candidates, actual_cols)
  if (length(exact) > 0) return(exact[1])

  # Coincidencia parcial 
  for (cand in candidates) {
    hits <- grep(tolower(gsub("[^a-z0-9]", "", cand)),
                 tolower(gsub("[^a-z0-9]", "", actual_cols)),
                 value = FALSE)
    if (length(hits) > 0) return(actual_cols[hits[1]])
  }
  return(NULL)
}

# 4. Eliminación de muestras contaminadas
remove_samples <- function(expr, samples_to_remove) {
  to_drop <- intersect(colnames(expr), samples_to_remove)
  if (length(to_drop) == 0) {
    cat("  No se encontraron muestras a eliminar en la matriz.\n")
    return(expr)
  }
  expr_clean <- expr[, !colnames(expr) %in% to_drop, drop = FALSE]
  cat(sprintf("  Muestras eliminadas (%d): %s\n",
              length(to_drop), paste(to_drop, collapse = ", ")))
  cat(sprintf("  Muestras restantes: %d\n", ncol(expr_clean)))
  return(expr_clean)
}

# 5. Filtrado, Kruskal-Wallis, top1000 y normalización VST
preprocess_expression <- function(expr, meta = NULL,
                                   var_quantile_cutoff = 0.25,
                                   n_top = 1000,
                                   kw_pval = 0.01,
                                   blind = TRUE, seed = 42) {
  set.seed(seed)

  # 5.1. Filtrar genes de baja expresión (al menos 10 counts en TODAS las muestras)
  keep      <- rowSums(expr >= 10) >= ncol(expr)
  expr_keep <- expr[keep, , drop = FALSE]
  cat(sprintf("  Genes tras filtro de expresion minima (todas las muestras): %d\n", nrow(expr_keep)))

  # 5.2. Eliminar el 25% de genes de menor varianza
  var_genes  <- apply(expr_keep, 1, var)
  var_cutoff <- quantile(var_genes, probs = var_quantile_cutoff, na.rm = TRUE)
  expr_var   <- expr_keep[var_genes >= var_cutoff, , drop = FALSE]
  cat(sprintf("  Genes tras filtro de varianza (75%% superior): %d\n", nrow(expr_var)))

  # 5.3. Top n_top genes mas variables (dentro del 75% filtrado)
  var_genes2 <- apply(expr_var, 1, var)
  top_n      <- names(sort(var_genes2, decreasing = TRUE))[1:min(n_top, nrow(expr_var))]

  # 5.4. Kruskal-Wallis: genes que difieren significativamente entre patotipos
  kw_genes <- character(0)
  if (!is.null(meta) && "pathotype_clean" %in% colnames(meta)) {
    common_s <- intersect(colnames(expr_var), rownames(meta))
    if (length(common_s) >= 10) {
      expr_kw          <- expr_var[, common_s, drop = FALSE]
      pathotype_factor <- factor(meta[common_s, "pathotype_clean"])
      kw_pvals <- apply(expr_kw, 1, function(g) {
        tryCatch(kruskal.test(g ~ pathotype_factor)$p.value, error = function(e) 1)
      })
      kw_genes <- names(kw_pvals[kw_pvals < kw_pval])
      cat(sprintf("  Genes significativos Kruskal-Wallis (p < %.3f): %d\n",
                  kw_pval, length(kw_genes)))
    } else {
      cat("  Kruskal-Wallis omitido: muestras en comun insuficientes.\n")
    }
  } else {
    cat("  Kruskal-Wallis omitido: metadatos de patotipo no disponibles.\n")
  }

  # 5.5. Union: top n_top + genes KW significativos
  final_genes <- union(top_n, kw_genes)
  expr_filt   <- expr_var[final_genes, , drop = FALSE]
  cat(sprintf("  Genes finales (union top%d + KW): %d\n", n_top, nrow(expr_filt)))

  # Asegurar nombres de columna antes de DESeq2
  if (is.null(colnames(expr_filt)) || any(is.na(colnames(expr_filt))))
    colnames(expr_filt) <- paste0("Sample_", seq_len(ncol(expr_filt)))

  # 5.6. Normalización VST con DESeq2
  col_data <- data.frame(
    condition = rep("all", ncol(expr_filt)),
    row.names = colnames(expr_filt)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = round(expr_filt),
    colData   = col_data,
    design    = ~1
  )
  vsd      <- vst(dds, blind = blind)
  expr_vst <- assay(vsd)
  cat("  Normalizacion VST aplicada.\n")

  # 5.7. Transponer (muestras x genes) y escalar
  expr_scaled <- scale(t(expr_vst))

  return(expr_scaled)
}

# 5. Alineación de muestras entre expresión y metadatos

align_samples <- function(expr_scaled, meta) {
  common <- intersect(rownames(expr_scaled), rownames(meta))

  if (length(common) == 0) {
    cat("  IDs en expresion (primeros 5):",
        paste(head(rownames(expr_scaled), 5), collapse = ", "), "\n")
    cat("  IDs en metadatos (primeros 5):",
        paste(head(rownames(meta), 5), collapse = ", "), "\n")
    stop("No hay muestras en comun entre expresion y metadatos. ",
         "Comprueba que los IDs de muestra coinciden (ver arriba).")
  }

  cat(sprintf("  Muestras en comun: %d\n", length(common)))

  expr_aligned <- expr_scaled[common, , drop = FALSE]
  meta_aligned <- meta[common, , drop = FALSE]

  return(list(expr = expr_aligned, meta = meta_aligned))
}
