suppressPackageStartupMessages({
  library(cluster)
  library(ggplot2)
})

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(z) z * (z - 1) / 2
  nij <- sum(choose2(tab))
  ai <- sum(choose2(rowSums(tab)))
  bj <- sum(choose2(colSums(tab)))
  n2 <- choose2(sum(tab))
  expected <- ai * bj / n2
  maximum <- (ai + bj) / 2
  if (maximum == expected) return(1)
  (nij - expected) / (maximum - expected)
}

out <- "results"
dir.create(file.path(out, "figures"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "tables"), recursive=TRUE, showWarnings=FALSE)

x <- readRDS("results/tables/expression_matrix_scaled.rds")
d <- dist(x, method="euclidean")
hc <- hclust(d, method="ward.D2")
ks <- 2:10
set.seed(42)

one_k <- function(k, B=100, frac=.8) {
  ref <- cutree(hc, k=k)
  sil <- silhouette(ref, d)[, "sil_width"]
  ari <- replicate(B, {
    idx <- sample(seq_len(nrow(x)), floor(nrow(x)*frac), replace=FALSE)
    db <- dist(x[idx,,drop=FALSE], method="euclidean")
    cb <- cutree(hclust(db, method="ward.D2"), k=k)
    adjusted_rand(ref[idx], cb)
  })
  data.frame(k=k, silhouette_mean=mean(sil), silhouette_sd=sd(sil),
             bootstrap_ARI_mean=mean(ari), bootstrap_ARI_sd=sd(ari),
             min_cluster_size=min(table(ref)), max_cluster_size=max(table(ref)))
}

res <- do.call(rbind, lapply(ks, one_k))
write.csv(res, file.path(out,"tables","k_selection_metrics.csv"), row.names=FALSE)

long <- rbind(
  data.frame(k=res$k, metric="Mean silhouette", value=res$silhouette_mean, sd=res$silhouette_sd),
  data.frame(k=res$k, metric="Bootstrap stability (mean ARI)", value=res$bootstrap_ARI_mean, sd=res$bootstrap_ARI_sd)
)
p <- ggplot(long, aes(k, value, color=metric, group=metric)) +
  geom_ribbon(aes(ymin=pmax(0,value-sd), ymax=pmin(1,value+sd), fill=metric), alpha=.10, color=NA) +
  geom_line(linewidth=1) + geom_point(size=2.8) +
  geom_vline(xintercept=5, linetype="dashed", color="black", linewidth=.6) +
  annotate("text", x=5.15, y=.98, label="Selected k = 5", hjust=0, vjust=1, size=4, color="black") +
  scale_x_continuous(breaks=ks) + scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values=c("Mean silhouette"="#2C7FB8", "Bootstrap stability (mean ARI)"="#D95F02")) +
  scale_fill_manual(values=c("Mean silhouette"="#2C7FB8", "Bootstrap stability (mean ARI)"="#D95F02")) +
  labs(title="Selection of the number of hierarchical clusters",
       subtitle="Ward.D2 clustering; Euclidean distance; 100 subsamples (80%) per k",
       x="Number of clusters (k)", y="Metric value", color=NULL, fill=NULL) +
  theme_bw(base_size=13) + theme(text=element_text(color="black"), axis.text=element_text(color="black"),
                                 legend.position="bottom", plot.title=element_text(face="bold"))
ggsave(file.path(out,"figures","k_selection_silhouette_bootstrap.png"), p, width=8.2, height=5.6, dpi=450, bg="white")
print(res)
