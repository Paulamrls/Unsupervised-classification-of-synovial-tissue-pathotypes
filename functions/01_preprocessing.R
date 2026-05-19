# ==============================================================================
# functions/01_preprocessing.R
# Carga de datos, filtrado por varianza, normalización VST, alineación
# ==============================================================================

# ── Carga de counts crudos ────────────────────────────────────────────────────
load_counts <- function(data_file) {
  if (!file.exists(data_file))
    stop("No se encuentra el archivo de datos: ", data_file)

  env <- new.env()
  load(data_file, envir = env)

  # Buscar el primer objeto matrix/data.frame en el entorno cargado
  obj_name <- ls(env)[1]
  expr <- get(obj_name, envir = env)
  expr <- as.matrix(expr)

  # Garantizar que existan nombres de fila y columna
  if (is.null(rownames(expr)))
    rownames(expr) <- paste0("Gene_", seq_len(nrow(expr)))
  if (is.null(colnames(expr)))
    colnames(expr) <- paste0("Sample_", seq_len(ncol(expr)))

  cat(sprintf("  Datos cargados: %d genes x %d muestras\n",
              nrow(expr), ncol(expr)))
  return(expr)
}

# ── Carga y limpieza de metadatos ─────────────────────────────────────────────
load_metadata <- function(meta_file, cols_keep) {
  if (!file.exists(meta_file))
    stop("No se encuentra el archivo de metadatos: ", meta_file)

  meta <- read.delim(meta_file, stringsAsFactors = FALSE, check.names = TRUE)

  # Mostrar columnas disponibles para diagnóstico
  cat("  Columnas disponibles en metadatos:\n")
  cat("   ", paste(colnames(meta)[1:min(10, ncol(meta))], collapse = "\n    "), "\n")

  # Normalizar nombres de columna (check.names=TRUE ya reemplaza espacios y
  # caracteres especiales por puntos, igual que cols_keep definidos en main)
  actual_cols <- colnames(meta)

  # Detectar columna de sampleid por coincidencia parcial
  id_col <- detect_column(actual_cols, c("sampleid", "sample.id", "SampleID",
                                          "Characteristics.sampleid."))
  if (is.null(id_col))
    stop("No se encontró columna de ID de muestra en metadatos. ",
         "Columnas disponibles: ", paste(actual_cols, collapse = ", "))

  # Detectar columna de patotipo por coincidencia parcial
  path_col <- detect_column(actual_cols, c("pathotype", "Pathotype",
                                            "Characteristics.pathotype."))

  # Seleccionar columnas existentes de cols_keep + las detectadas automáticamente
  cols_found <- intersect(cols_keep, actual_cols)
  extra_cols <- unique(c(id_col, path_col))
  extra_cols <- extra_cols[!is.na(extra_cols)]
  cols_final <- unique(c(cols_found, extra_cols))

  meta_sub <- meta[, cols_final, drop = FALSE]

  # Eliminar duplicados por sampleid
  meta_unique <- meta_sub[!duplicated(meta_sub[[id_col]]), ]
  rownames(meta_unique) <- as.character(meta_unique[[id_col]])

  # Limpiar columna de patotipo
  if (!is.null(path_col) && path_col %in% colnames(meta_unique)) {
    meta_unique$pathotype_clean <- as.character(meta_unique[[path_col]])

    # Unificar variantes de nombre
    meta_unique$pathotype_clean[
      meta_unique$pathotype_clean == "Fibrous"] <- "Fibroid"

    # Eliminar muestras sin clasificación válida
    meta_unique <- meta_unique[
      !is.na(meta_unique$pathotype_clean) &
        meta_unique$pathotype_clean != "Ungraded", ]

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

# ── Helper: detectar columna por coincidencia parcial ─────────────────────────
detect_column <- function(actual_cols, candidates) {
  # Primero busca coincidencia exacta
  exact <- intersect(candidates, actual_cols)
  if (length(exact) > 0) return(exact[1])

  # Luego busca coincidencia parcial (insensible a mayúsculas)
  for (cand in candidates) {
    hits <- grep(tolower(gsub("[^a-z0-9]", "", cand)),
                 tolower(gsub("[^a-z0-9]", "", actual_cols)),
                 value = FALSE)
    if (length(hits) > 0) return(actual_cols[hits[1]])
  }
  return(NULL)
}

# ── Filtrado y normalización VST ──────────────────────────────────────────────
preprocess_expression <- function(expr, n_top = 1000, blind = TRUE,
                                   seed = 42) {
  set.seed(seed)

  # 1. Filtrar genes de baja expresión (al menos 10 counts en >= 2 muestras)
  keep      <- rowSums(expr >= 10) >= 2
  expr_keep <- expr[keep, , drop = FALSE]
  cat(sprintf("  Genes tras filtro de expresion minima: %d\n", nrow(expr_keep)))

  # 2. Seleccionar los n_top genes más variables (por varianza)
  var_genes <- apply(expr_keep, 1, var)
  n_select  <- min(n_top, length(var_genes))
  top_genes <- names(sort(var_genes, decreasing = TRUE))[seq_len(n_select)]
  expr_filt <- expr_keep[top_genes, , drop = FALSE]
  cat(sprintf("  Top %d genes mas variables seleccionados.\n", n_select))

  # Garantizar nombres de columna antes de DESeq2
  if (is.null(colnames(expr_filt)) || any(is.na(colnames(expr_filt))))
    colnames(expr_filt) <- paste0("Sample_", seq_len(ncol(expr_filt)))

  # 3. Normalización VST con DESeq2
  col_data <- data.frame(
    condition = rep("all", ncol(expr_filt)),
    row.names = colnames(expr_filt)
  )

  dds <- DESeqDataSetFromMatrix(
    countData = round(expr_filt),   # VST requiere enteros
    colData   = col_data,
    design    = ~1
  )
  vsd      <- vst(dds, blind = blind)
  expr_vst <- assay(vsd)
  cat("  Normalizacion VST aplicada.\n")

  # 4. Transponer (muestras x genes) y escalar genes
  expr_t      <- t(expr_vst)
  expr_scaled <- scale(expr_t)

  return(expr_scaled)
}

# ── Alinear muestras entre expresión y metadatos ──────────────────────────────
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
