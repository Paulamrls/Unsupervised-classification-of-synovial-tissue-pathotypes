# =============================================================================
# Subtype4 clinical characterisation
# Checks whether Subtype4 samples are clinically unusual compared to the other
# four STRAP pathotype clusters (Fibroid, Myeloid, Lymphoid, IFN-high).
# =============================================================================

library(dplyr)

# ---- 0. File paths (adjust if needed) ---------------------------------------

cluster_file <- "results_k5/results/clusters/cluster_assignments.csv"
sdrf_file    <- "E-MTAB-13733.sdrf.txt"

# ---- 1. Load cluster assignments --------------------------------------------

clusters <- tryCatch(
  read.csv(cluster_file, stringsAsFactors = FALSE),
  error = function(e) stop("Cannot read cluster file: ", cluster_file, "\n", e)
)

# Expected columns: sample_id, cluster_bio
# Keep only what we need and rename for clarity
clusters <- clusters[, c("sample_id", "cluster_bio")]
colnames(clusters) <- c("sampleid", "cluster")

cat("Cluster counts:\n")
print(table(clusters$cluster))

# ---- 2. Load SDRF metadata --------------------------------------------------

sdrf <- tryCatch(
  read.delim(sdrf_file, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop("Cannot read SDRF file: ", sdrf_file, "\n", e)
)

# Select and rename the columns of interest
col_map <- c(
  sampleid  = "Characteristics[sampleid]",
  pathotype = "Characteristics[pathotype]",
  age       = "Characteristics[age]",
  sex       = "Characteristics[sex]",
  tjc       = "Characteristics[tender_joints_counts_tjc]",
  sjc       = "Characteristics[swollen_joints_counts_sjc]",
  crp       = "Characteristics[crp]",
  esr       = "Characteristics[esr]",
  rf_status = "Characteristics[rf_status]",
  bcell     = "Characteristics[cell.type (bcell)]"
)

# Keep only columns that actually exist in the file
present <- col_map[col_map %in% colnames(sdrf)]
missing <- setdiff(names(col_map), names(present))
if (length(missing) > 0)
  warning("These metadata columns were not found and will be skipped: ",
          paste(missing, collapse = ", "))

meta <- sdrf[, present, drop = FALSE]
colnames(meta) <- names(present)

# Deduplicate by sampleid (SDRF often has one row per file, not per sample)
meta <- meta[!duplicated(meta$sampleid), ]

# ---- 3. Merge ---------------------------------------------------------------

df <- merge(clusters, meta, by = "sampleid", all.x = TRUE)

# Coerce numeric columns
for (col in intersect(c("age", "tjc", "sjc", "crp", "esr"), colnames(df))) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}

cat("\nMerged dataset: ", nrow(df), "samples,", ncol(df), "columns\n")
cat("Samples without metadata match:", sum(is.na(df$age)), "\n\n")

# ---- 4. Subtype4 summary vs other clusters ----------------------------------

df$group <- ifelse(df$cluster == "Subtype4", "Subtype4", "Other clusters")

cat("=================================================================\n")
cat("SUBTYPE4 CHARACTERISATION\n")
cat("=================================================================\n\n")

## 4a. Sample count
cat("-- Sample counts per cluster --\n")
print(table(df$cluster))
cat("\n")

## 4b. Pathotype distribution
if ("pathotype" %in% colnames(df)) {
  cat("-- Pathotype distribution (%) --\n")
  pt <- prop.table(table(df$cluster, df$pathotype), margin = 1) * 100
  print(round(pt, 1))
  cat("\n")
}

## 4c. Sex distribution
if ("sex" %in% colnames(df)) {
  cat("-- Sex distribution (%) --\n")
  sx <- prop.table(table(df$cluster, df$sex), margin = 1) * 100
  print(round(sx, 1))
  cat("\n")
}

## 4d. RF status distribution
if ("rf_status" %in% colnames(df)) {
  cat("-- RF status distribution (%) --\n")
  rf <- prop.table(table(df$cluster, df$rf_status), margin = 1) * 100
  print(round(rf, 1))
  cat("\n")
}

## 4e. Continuous variable means per cluster
num_vars <- intersect(c("age", "crp", "esr", "tjc", "sjc"), colnames(df))
if (length(num_vars) > 0) {
  cat("-- Mean continuous variables per cluster --\n")
  summ <- df %>%
    group_by(cluster) %>%
    summarise(
      n        = n(),
      across(all_of(num_vars),
             ~ round(mean(.x, na.rm = TRUE), 1),
             .names = "mean_{.col}"),
      .groups = "drop"
    ) %>%
    arrange(cluster != "Subtype4", cluster)   # Subtype4 first
  print(as.data.frame(summ))
  cat("\n")
}

## 4f. Quick Kruskal-Wallis tests (non-parametric, no normality assumption)
cat("-- Kruskal-Wallis p-values (cluster differences) --\n")
for (v in num_vars) {
  kt <- kruskal.test(df[[v]] ~ df$cluster)
  cat(sprintf("  %-6s  p = %.4f\n", v, kt$p.value))
}
cat("\n")

# ---- 5. Boxplots: CRP and TJC across all 5 clusters -------------------------

cluster_order <- c("Fibroid", "Myeloid", "Lymphoid", "IFN-high", "Subtype4")
df$cluster_f  <- factor(df$cluster, levels = cluster_order)

# Colour palette: Subtype4 highlighted in orange
pal <- c("Fibroid"   = "#4DAF4A",
         "Myeloid"   = "#377EB8",
         "Lymphoid"  = "#984EA3",
         "IFN-high"  = "#E41A1C",
         "Subtype4"  = "#FF7F00")

old_par <- par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))

if ("crp" %in% colnames(df)) {
  boxplot(crp ~ cluster_f, data = df,
          col     = pal[levels(df$cluster_f)],
          ylab    = "CRP (mg/L)",
          xlab    = "",
          main    = "CRP by cluster",
          las     = 2,
          outline = FALSE)
  stripchart(crp ~ cluster_f, data = df,
             method = "jitter", add = TRUE, pch = 20,
             col = adjustcolor("black", alpha.f = 0.3), cex = 0.6)
}

if ("tjc" %in% colnames(df)) {
  boxplot(tjc ~ cluster_f, data = df,
          col     = pal[levels(df$cluster_f)],
          ylab    = "Tender joint count (TJC)",
          xlab    = "",
          main    = "TJC by cluster",
          las     = 2,
          outline = FALSE)
  stripchart(tjc ~ cluster_f, data = df,
             method = "jitter", add = TRUE, pch = 20,
             col = adjustcolor("black", alpha.f = 0.3), cex = 0.6)
}

par(old_par)

cat("=================================================================\n")
cat("Interpretation hint:\n")
cat("  If Subtype4 shows SIMILAR CRP/TJC/pathotype to other clusters,\n")
cat("  the muscle-contamination hypothesis is less supported (the\n")
cat("  cluster could be a real biological variant).\n")
cat("  If Subtype4 is clinically unremarkable but transcriptomically\n")
cat("  distinct, contamination or batch effects become more likely.\n")
cat("=================================================================\n")
