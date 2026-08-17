suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(cluster)
  library(mclust)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(clusterProfiler)
  library(fgsea)
})

set.seed(42)
project_dir <- "C:/Users/Paula/Desktop/bioinformatica/TFM/codigo/Unsupervised-classification-of-synovial-tissue-pathotypes"
task_dir <- "C:/Users/Paula/Documents/Codex/2026-08-05/necesito-que-revises-mi-pipeline-de"
root <- file.path(task_dir, "outputs", "results_clean_unsupervised_relaxed_p10_g3000_pc30")
out <- file.path(root, "comparison_k4_k5", "functional_validation_ID3_ID4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

contaminated <- c("STRAPPAT00113-baseline", "STRAPPAT00210-baseline",
                  "STRAPPAT00199-baseline", "STRAPPAT00007-baseline",
                  "STRAPPAT00051-baseline")
ifn_genes <- c("MX1", "OAS1", "IFIT1", "IFIT3", "ISG15", "RSAD2",
               "IFI44L", "CXCL10", "IRF7", "STAT1")

load_counts <- function() {
  env <- new.env()
  load(file.path(project_dir, "strap_counts.RData"), envir = env)
  candidates <- ls(env)[vapply(ls(env), function(x) {
    is.matrix(get(x, env)) || is.data.frame(get(x, env))
  }, logical(1))]
  x <- as.matrix(get(candidates[1], env))
  x[, !colnames(x) %in% contaminated, drop = FALSE]
}

counts_raw <- load_counts()
k4 <- read.csv(file.path(root, "k4", "cluster_assignments.csv"),
               stringsAsFactors = FALSE)
k5 <- read.csv(file.path(root, "k5", "cluster_assignments.csv"),
               stringsAsFactors = FALSE)
stopifnot(identical(k4$sample_id, k5$sample_id))
samples <- k5$sample_id
counts <- counts_raw[, samples, drop = FALSE]

# Reconstruir exactamente la matriz comun usada para PCA/clustering.
keep <- rowSums(counts >= 10) >= ceiling(0.10 * ncol(counts))
counts_filt <- counts[keep, , drop = FALSE]
cd_all <- data.frame(condition = factor(rep("all", ncol(counts_filt))),
                     row.names = colnames(counts_filt))
dds_all <- DESeqDataSetFromMatrix(round(counts_filt), cd_all, ~1)
vst_all <- assay(vst(dds_all, blind = TRUE))
vars <- apply(vst_all, 1, var)
genes3000 <- names(sort(vars, decreasing = TRUE))[seq_len(3000)]
expr_scaled <- scale(t(vst_all[genes3000, , drop = FALSE]))
pca <- prcomp(expr_scaled, center = FALSE, scale. = FALSE)
pcs30 <- pca$x[, seq_len(30), drop = FALSE]

# Verificacion de reproducibilidad exacta de k-means.
set.seed(42)
km <- kmeans(pcs30, centers = 5, nstart = 50, iter.max = 300)
cl_recomputed <- setNames(km$cluster, rownames(pcs30))
cl_exported <- setNames(k5$cluster_id, k5$sample_id)
ari_recomputed <- adjustedRandIndex(cl_exported[names(cl_recomputed)], cl_recomputed)
if (ari_recomputed != 1) stop("La reconstruccion no reproduce exactamente k5")

# 1. DESeq2 directo: log2FC positivo = ID4 mayor que ID3.
sub_ids <- k5$sample_id[k5$cluster_id %in% c(3, 4)]
groups <- factor(paste0("ID", k5$cluster_id[match(sub_ids, k5$sample_id)]),
                 levels = c("ID3", "ID4"))
cd <- data.frame(group = groups, row.names = sub_ids)
dds <- DESeqDataSetFromMatrix(round(counts_raw[, sub_ids, drop = FALSE]), cd, ~group)
dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("group", "ID4", "ID3")))
res$gene <- rownames(res)
res <- res[order(res$padj, na.last = TRUE), ]
write.csv(res, file.path(out, "DESeq2_ID4_vs_ID3_all.csv"), row.names = FALSE)
sig <- subset(res, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
write.csv(sig, file.path(out, "DESeq2_ID4_vs_ID3_significant.csv"), row.names = FALSE)

# 2. ORA GO-BP para cada direccion.
run_ego <- function(symbols, label) {
  if (length(symbols) < 10) return(NULL)
  eg <- enrichGO(gene = unique(symbols), OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.20, readable = TRUE)
  df <- as.data.frame(eg)
  write.csv(df, file.path(out, paste0("GO_BP_ORA_", label, ".csv")),
            row.names = FALSE)
  df
}
ego4 <- run_ego(sig$gene[sig$log2FoldChange > 0], "higher_ID4")
ego3 <- run_ego(sig$gene[sig$log2FoldChange < 0], "higher_ID3")

# 3. GSEA GO-BP sobre todos los genes con estadistico DESeq2.
ranked <- res$stat
names(ranked) <- res$gene
ranked <- sort(ranked[is.finite(ranked) & !duplicated(names(ranked))],
               decreasing = TRUE)
go_map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(names(ranked)),
                                columns = c("GO", "ONTOLOGY"),
                                keytype = "SYMBOL")
go_map <- unique(go_map[!is.na(go_map$GO) & go_map$ONTOLOGY == "BP",
                        c("GO", "SYMBOL")])
pathways <- split(go_map$SYMBOL, go_map$GO)
pathways <- pathways[lengths(pathways) >= 15 & lengths(pathways) <= 500]
fg <- fgseaMultilevel(pathways = pathways, stats = ranked,
                      minSize = 15, maxSize = 500, eps = 0)
fg <- as.data.frame(fg)
go_terms <- AnnotationDbi::select(GO.db, keys = unique(fg$pathway),
                                  columns = "TERM", keytype = "GOID")
term_map <- setNames(go_terms$TERM, go_terms$GOID)
fg$description <- unname(term_map[fg$pathway])
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", character(1))
fg <- fg[order(fg$padj, -abs(fg$NES)), ]
write.csv(fg, file.path(out, "GSEA_GO_BP_ID4_vs_ID3.csv"), row.names = FALSE)

# 4. Firma IFN usando genes disponibles en toda la VST relajada.
ifn_avail <- intersect(ifn_genes, rownames(vst_all))
ifn_z <- t(scale(t(vst_all[ifn_avail, , drop = FALSE])))
ifn_score <- colMeans(ifn_z, na.rm = TRUE)
ifn_df <- data.frame(sample_id = samples, cluster_k5 = k5$cluster_id,
                     ifn_score = unname(ifn_score[samples]))
write.csv(ifn_df, file.path(out, "IFN_score_all_k5_samples.csv"), row.names = FALSE)
ifn_sub <- subset(ifn_df, cluster_k5 %in% c(3, 4))
ifn_test <- wilcox.test(ifn_score ~ factor(cluster_k5), data = ifn_sub,
                        exact = FALSE, conf.int = TRUE)
ifn_summary <- do.call(rbind, lapply(split(ifn_sub$ifn_score, ifn_sub$cluster_k5),
                                     function(x) data.frame(n = length(x),
                                       mean = mean(x), sd = sd(x), median = median(x),
                                       q1 = quantile(x, .25), q3 = quantile(x, .75))))
ifn_summary$cluster_k5 <- as.integer(rownames(ifn_summary))
ifn_summary$p_value_ID3_vs_ID4 <- ifn_test$p.value
write.csv(ifn_summary, file.path(out, "IFN_summary_ID3_ID4.csv"), row.names = FALSE)

# 5. Comparaciones clinicas no parametricas y sexo por Fisher.
meta <- read.delim(file.path(project_dir, "E-MTAB-13733.sdrf.txt"),
                   stringsAsFactors = FALSE, check.names = TRUE)
detect <- function(pattern) grep(pattern, names(meta), ignore.case = TRUE, value = TRUE)[1]
id_col <- detect("sampleid")
clinical_patterns <- c(tender = "tender.*joints", swollen = "swollen.*joints",
                       activity = "arthritis.*activity", esr = "(^|\\.)esr($|\\.)",
                       crp = "(^|\\.)crp($|\\.)", physician_global = "physicians.*global")
meta_ids <- as.character(meta[[id_col]])
clinical_rows <- list()
for (nm in names(clinical_patterns)) {
  col <- detect(clinical_patterns[[nm]])
  if (is.na(col)) next
  vals <- suppressWarnings(as.numeric(meta[[col]][match(sub_ids, meta_ids)]))
  g <- groups
  valid <- is.finite(vals)
  x3 <- vals[valid & g == "ID3"]
  x4 <- vals[valid & g == "ID4"]
  if (length(x3) >= 3 && length(x4) >= 3) {
    wt <- wilcox.test(x3, x4, exact = FALSE)
    clinical_rows[[nm]] <- data.frame(
      variable = nm, n_ID3 = length(x3), median_ID3 = median(x3),
      IQR_ID3 = IQR(x3), n_ID4 = length(x4), median_ID4 = median(x4),
      IQR_ID4 = IQR(x4), p_value = wt$p.value)
  }
}
clinical <- do.call(rbind, clinical_rows)
clinical$padj_BH <- p.adjust(clinical$p_value, method = "BH")
write.csv(clinical, file.path(out, "clinical_comparison_ID3_ID4.csv"), row.names = FALSE)
sex_col <- detect("(^|\\.)sex($|\\.)")
if (!is.na(sex_col)) {
  sex <- meta[[sex_col]][match(sub_ids, meta_ids)]
  sex_tab <- table(cluster = groups, sex = sex, useNA = "no")
  write.csv(as.data.frame.matrix(sex_tab), file.path(out, "sex_contingency_ID3_ID4.csv"))
  write.csv(data.frame(p_value = fisher.test(sex_tab)$p.value),
            file.path(out, "sex_fisher_test_ID3_ID4.csv"), row.names = FALSE)
}

# 6. Silhouette individual y Jaccard bootstrap por cluster.
sil <- silhouette(k5$cluster_id, dist(pcs30))
sil_df <- data.frame(sample_id = rownames(pcs30), cluster_k5 = sil[, "cluster"],
                     silhouette = sil[, "sil_width"])
write.csv(sil_df, file.path(out, "silhouette_by_sample_k5.csv"), row.names = FALSE)
sil_summary <- do.call(rbind, lapply(split(sil_df$silhouette, sil_df$cluster_k5),
                                     function(x) data.frame(
                                       n = length(x), mean = mean(x),
                                       median = median(x), min = min(x),
                                       negative_pct = 100 * mean(x < 0))))
sil_summary$cluster_k5 <- as.integer(rownames(sil_summary))
rownames(sil_summary) <- NULL
sil_summary <- sil_summary[, c("cluster_k5", "n", "mean", "median", "min",
                               "negative_pct")]
write.csv(sil_summary, file.path(out, "silhouette_by_cluster_k5.csv"), row.names = FALSE)

B <- 100
jaccard <- matrix(NA_real_, nrow = B, ncol = 5,
                  dimnames = list(NULL, paste0("ID", 1:5)))
for (b in seq_len(B)) {
  set.seed(42 + b)
  idx <- sample(seq_len(nrow(pcs30)), floor(0.8 * nrow(pcs30)), replace = FALSE)
  xb <- pcs30[idx, , drop = FALSE]
  set.seed(42 + b)
  kb <- kmeans(xb, centers = 5, nstart = 50, iter.max = 300)$cluster
  ref <- k5$cluster_id[idx]
  for (r in 1:5) {
    a <- ref == r
    scores <- vapply(1:5, function(c) {
      z <- kb == c
      sum(a & z) / sum(a | z)
    }, numeric(1))
    jaccard[b, r] <- max(scores)
  }
}
write.csv(data.frame(iteration = seq_len(B), jaccard),
          file.path(out, "bootstrap_cluster_jaccard_iterations.csv"), row.names = FALSE)
j_summary <- data.frame(
  cluster_k5 = 1:5, mean_jaccard = colMeans(jaccard),
  sd = apply(jaccard, 2, sd), median = apply(jaccard, 2, median),
  q05 = apply(jaccard, 2, quantile, .05),
  q95 = apply(jaccard, 2, quantile, .95),
  pct_ge_075 = 100 * colMeans(jaccard >= .75)
)
write.csv(j_summary, file.path(out, "bootstrap_cluster_jaccard_summary.csv"), row.names = FALSE)

# Figuras compactas.
p_ifn <- ggplot(ifn_sub, aes(factor(cluster_k5), ifn_score,
                             fill = factor(cluster_k5))) +
  geom_boxplot(outlier.shape = NA, alpha = .75) +
  geom_jitter(width = .12, alpha = .6) + theme_bw() +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Firma IFN: k5-ID3 frente a k5-ID4", x = "ID k5", y = "IFN score")
ggsave(file.path(out, "IFN_ID3_vs_ID4.png"), p_ifn, width = 6, height = 5, dpi = 300)

p_sil <- ggplot(subset(sil_df, cluster_k5 %in% c(3,4)),
                aes(factor(cluster_k5), silhouette, fill = factor(cluster_k5))) +
  geom_boxplot(outlier.shape = NA, alpha = .75) + geom_jitter(width = .12, alpha = .6) +
  theme_bw() + scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Silhouette individual: k5-ID3 frente a k5-ID4",
       x = "ID k5", y = "Silhouette")
ggsave(file.path(out, "silhouette_ID3_vs_ID4.png"), p_sil,
       width = 6, height = 5, dpi = 300)

summary <- data.frame(
  item = c("n_ID3", "n_ID4", "recomputed_ARI", "DEG_absLFC1_padj005",
           "IFN_genes_used", "IFN_p_value"),
  value = c(sum(groups == "ID3"), sum(groups == "ID4"), ari_recomputed,
            nrow(sig), length(ifn_avail), ifn_test$p.value)
)
write.csv(summary, file.path(out, "validation_summary.csv"), row.names = FALSE)

cat("Validation complete\n")
print(summary)
print(ifn_summary)
print(clinical)
print(subset(sil_summary, cluster_k5 %in% c(3,4)))
print(subset(j_summary, cluster_k5 %in% c(3,4)))
