################################################################################
# functions/01_preprocessing.R
# Carga, filtro minimo, VST, seleccion de genes y alineacion.
################################################################################

load_counts <- function(data_file) {
  if (!file.exists(data_file)) stop("No se encuentra el archivo: ", data_file)
  env <- new.env()
  load(data_file, envir = env)
  candidates <- ls(env)[vapply(ls(env), function(x) {
    is.matrix(get(x, env)) || is.data.frame(get(x, env))
  }, logical(1))]
  if (length(candidates) == 0L) stop("El RData no contiene una matriz/data.frame.")
  expr <- as.matrix(get(candidates[1], envir = env))
  if (is.null(rownames(expr))) rownames(expr) <- paste0("Gene_", seq_len(nrow(expr)))
  if (is.null(colnames(expr))) colnames(expr) <- paste0("Sample_", seq_len(ncol(expr)))
  cat(sprintf("  Datos cargados: %d genes x %d muestras\n", nrow(expr), ncol(expr)))
  expr
}

detect_column <- function(actual_cols, candidates) {
  exact <- intersect(candidates, actual_cols)
  if (length(exact) > 0L) return(exact[1])
  actual_clean <- tolower(gsub("[^a-z0-9]", "", actual_cols))
  for (cand in candidates) {
    hits <- grep(tolower(gsub("[^a-z0-9]", "", cand)), actual_clean)
    if (length(hits) > 0L) return(actual_cols[hits[1]])
  }
  NULL
}

load_metadata <- function(meta_file, cols_keep, filter_pathotype = FALSE) {
  if (!file.exists(meta_file)) stop("No se encuentra el archivo: ", meta_file)
  meta <- read.delim(meta_file, stringsAsFactors = FALSE, check.names = TRUE)
  actual_cols <- colnames(meta)
  id_col <- detect_column(actual_cols, c(
    "sampleid", "sample.id", "SampleID", "Characteristics.sampleid."
  ))
  if (is.null(id_col)) stop("No se encontro la columna de ID de muestra.")
  path_col <- detect_column(actual_cols, c(
    "pathotype", "Pathotype", "Characteristics.pathotype."
  ))
  cols_final <- unique(c(intersect(cols_keep, actual_cols), id_col, path_col))
  cols_final <- cols_final[!is.na(cols_final)]
  meta_sub <- meta[, cols_final, drop = FALSE]
  meta_unique <- meta_sub[!duplicated(meta_sub[[id_col]]), , drop = FALSE]
  rownames(meta_unique) <- as.character(meta_unique[[id_col]])
  if (!is.null(path_col) && path_col %in% colnames(meta_unique)) {
    meta_unique$pathotype_clean <- as.character(meta_unique[[path_col]])
    meta_unique$pathotype_clean[meta_unique$pathotype_clean == "Fibrous"] <- "Fibroid"
    if (filter_pathotype) {
      meta_unique <- meta_unique[!is.na(meta_unique$pathotype_clean) &
                                   meta_unique$pathotype_clean != "Ungraded", , drop = FALSE]
    }
  } else {
    meta_unique$pathotype_clean <- NA_character_
  }
  cat(sprintf("  Metadatos cargados: %d muestras unicas (sin filtrar por histologia)\n",
              nrow(meta_unique)))
  meta_unique
}

remove_samples <- function(expr, samples_to_remove) {
  to_drop <- intersect(colnames(expr), samples_to_remove)
  expr_clean <- expr[, !colnames(expr) %in% samples_to_remove, drop = FALSE]
  cat(sprintf("  Muestras eliminadas (%d): %s\n", length(to_drop),
              paste(to_drop, collapse = ", ")))
  cat(sprintf("  Muestras restantes: %d\n", ncol(expr_clean)))
  expr_clean
}

preprocess_expression <- function(expr, meta = NULL,
                                  method = c("top_variable", "top_variable_plus_kw"),
                                  min_count = 10, min_samples = NULL,
                                  var_quantile_cutoff = 0.25,
                                  n_top = 1000, kw_pval = 0.01,
                                  blind = TRUE, seed = 42) {
  set.seed(seed)
  method <- match.arg(method)

  # El metodo principal no consulta meta ni ninguna etiqueta histologica.
  if (is.null(min_samples)) min_samples <- ncol(expr)
  keep <- rowSums(expr >= min_count) >= min_samples
  expr_keep <- expr[keep, , drop = FALSE]
  cat(sprintf("  Genes tras filtro de expresion minima: %d\n", nrow(expr_keep)))
  if (nrow(expr_keep) < n_top) stop("El filtro deja menos de ", n_top, " genes.")

  col_data <- data.frame(condition = rep("all", ncol(expr_keep)),
                         row.names = colnames(expr_keep))
  dds <- DESeqDataSetFromMatrix(countData = round(expr_keep),
                                colData = col_data, design = ~1)
  vsd <- vst(dds, blind = blind)
  expr_vst_all <- assay(vsd)
  cat("  Normalizacion VST aplicada antes de seleccionar genes.\n")

  vst_variances <- apply(expr_vst_all, 1, var)
  top_n <- names(sort(vst_variances, decreasing = TRUE))[seq_len(n_top)]
  final_genes <- top_n

  # Compatibilidad reproducible con el enfoque anterior. Esta es la unica rama
  # de seleccion que puede usar pathotype_clean y esta desactivada por defecto.
  if (identical(method, "top_variable_plus_kw")) {
    raw_variances <- apply(expr_keep, 1, var)
    cutoff <- quantile(raw_variances, probs = var_quantile_cutoff, na.rm = TRUE)
    expr_var <- expr_keep[raw_variances >= cutoff, , drop = FALSE]
    kw_genes <- character(0)
    if (!is.null(meta) && "pathotype_clean" %in% colnames(meta)) {
      common <- intersect(colnames(expr_var), rownames(meta))
      valid <- common[!is.na(meta[common, "pathotype_clean"]) &
                        meta[common, "pathotype_clean"] != "Ungraded"]
      if (length(valid) >= 10L &&
          length(unique(meta[valid, "pathotype_clean"])) >= 2L) {
        groups <- factor(meta[valid, "pathotype_clean"])
        pvals <- apply(expr_var[, valid, drop = FALSE], 1, function(g) {
          tryCatch(kruskal.test(g ~ groups)$p.value, error = function(e) 1)
        })
        kw_genes <- names(pvals[pvals < kw_pval])
      }
    }
    final_genes <- union(top_n, kw_genes)
    cat(sprintf("  top_variable_plus_kw: %d genes (top%d + KW).\n",
                length(final_genes), n_top))
  } else {
    cat(sprintf("  top_variable: exactamente %d genes; histologia no consultada.\n",
                length(final_genes)))
  }

  expr_scaled <- scale(t(expr_vst_all[final_genes, , drop = FALSE]))
  attr(expr_scaled, "selected_genes") <- final_genes
  attr(expr_scaled, "gene_selection_method") <- method
  expr_scaled
}

align_samples <- function(expr_scaled, meta) {
  common <- intersect(rownames(expr_scaled), rownames(meta))
  if (length(common) == 0L) stop("No hay muestras comunes entre expresion y metadatos.")
  cat(sprintf("  Muestras en comun: %d\n", length(common)))
  list(expr = expr_scaled[common, , drop = FALSE],
       meta = meta[common, , drop = FALSE])
}
