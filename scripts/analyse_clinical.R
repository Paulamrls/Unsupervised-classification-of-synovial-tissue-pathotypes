suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

out <- "results"
dir.create(file.path(out,"figures"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,"tables"), recursive=TRUE, showWarnings=FALSE)

s <- read.delim("data/E-MTAB-13733.sdrf.txt", check.names=FALSE, stringsAsFactors=FALSE)
a <- read.csv("results/tables/sample_assignments_hierarchical_k5.csv", stringsAsFactors=FALSE)

take <- c(
  sample_id="Characteristics[sampleid]", age="Characteristics[age]", sex="Characteristics[sex]",
  rf_status="Characteristics[rf_status]", ccp_status="Characteristics[ccp_status]",
  TJC="Characteristics[tender_joints_counts_tjc]", SJC="Characteristics[swollen_joints_counts_sjc]",
  arthritis_activity="Characteristics[arthritis.activity]", ESR="Characteristics[esr]",
  CRP="Characteristics[crp]", physician_global="Characteristics[physicians.global.assessment]"
)
meta <- s[, unname(take), drop=FALSE]
names(meta) <- names(take)
meta <- meta[!duplicated(meta$sample_id),]
dat <- merge(a, meta, by="sample_id", all.x=TRUE, sort=FALSE)
dat$cluster <- factor(dat$cluster_id, levels=1:5, labels=paste("Cluster",1:5))

continuous <- c("age","TJC","SJC","arthritis_activity","ESR","CRP","physician_global")
categorical <- c("sex","rf_status","ccp_status")
for (v in continuous) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
for (v in categorical) dat[[v]][dat[[v]] %in% c("", "NA", "not available")] <- NA

cont_tests <- do.call(rbind, lapply(continuous, function(v) {
  z <- dat[,c("cluster",v)]; names(z)[2] <- "value"; z <- z[!is.na(z$value),]
  kt <- kruskal.test(value ~ cluster, data=z)
  data.frame(variable=v, n=nrow(z), test="Kruskal-Wallis", statistic=unname(kt$statistic), df=unname(kt$parameter), p_value=kt$p.value)
}))
cont_tests$p_BH <- p.adjust(cont_tests$p_value, method="BH")

cat_tests <- do.call(rbind, lapply(categorical, function(v) {
  z <- dat[,c("cluster",v)]; names(z)[2] <- "value"; z <- z[!is.na(z$value),]
  tab <- table(z$cluster, z$value)
  ft <- fisher.test(tab, simulate.p.value=FALSE)
  data.frame(variable=v, n=sum(tab), test="Fisher exact", statistic=NA_real_, df=NA_real_, p_value=ft$p.value)
}))
cat_tests$p_BH <- p.adjust(cat_tests$p_value, method="BH")
tests <- rbind(cont_tests, cat_tests)
write.csv(tests, file.path(out,"tables","clinical_tests_by_cluster.csv"), row.names=FALSE)
write.csv(dat, file.path(out,"tables","clinical_metadata_with_clusters.csv"), row.names=FALSE)

summ <- do.call(rbind, lapply(continuous, function(v) do.call(rbind, lapply(levels(dat$cluster), function(cl) {
  x <- dat[dat$cluster==cl,v]; x <- x[!is.na(x)]
  data.frame(variable=v, cluster=cl, n=length(x), median=median(x), Q1=unname(quantile(x,.25)), Q3=unname(quantile(x,.75)), min=min(x), max=max(x))
}))))
write.csv(summ, file.path(out,"tables","clinical_continuous_summary.csv"), row.names=FALSE)

catsumm <- dat %>% select(cluster, all_of(categorical)) %>%
  pivot_longer(-cluster, names_to="variable", values_to="category") %>%
  filter(!is.na(category)) %>% count(variable, cluster, category, name="n") %>%
  group_by(variable,cluster) %>% mutate(percent=100*n/sum(n)) %>% ungroup()
write.csv(catsumm, file.path(out,"tables","clinical_categorical_summary.csv"), row.names=FALSE)

# Pairwise exploratory comparisons for continuous variables retaining BH significance.
sig_cont <- cont_tests$variable[cont_tests$p_BH < .05]
pairwise <- do.call(rbind, lapply(sig_cont, function(v) {
  z <- dat[,c("cluster",v)]; names(z)[2] <- "value"; z <- z[!is.na(z$value),]
  pw <- pairwise.wilcox.test(z$value, z$cluster, p.adjust.method="BH", exact=FALSE)
  idx <- which(!is.na(pw$p.value), arr.ind=TRUE)
  data.frame(variable=v, cluster_1=rownames(pw$p.value)[idx[,1]],
             cluster_2=colnames(pw$p.value)[idx[,2]], p_BH=pw$p.value[idx])
}))
write.csv(pairwise, file.path(out,"tables","clinical_pairwise_wilcoxon.csv"), row.names=FALSE)

pretty <- c(age="Age (years)", TJC="Tender joint count", SJC="Swollen joint count",
            arthritis_activity="Arthritis activity", ESR="ESR", CRP="CRP",
            physician_global="Physician global assessment")
labels <- setNames(sprintf("%s\nKruskal-Wallis p = %.3g", pretty[cont_tests$variable], cont_tests$p_value), cont_tests$variable)
long <- dat %>% select(cluster, all_of(continuous)) %>% pivot_longer(-cluster, names_to="variable", values_to="value")
long$panel <- factor(labels[long$variable], levels=labels[continuous])
cluster_palette <- c("Cluster 1"="#009E73", "Cluster 2"="#E69F00", "Cluster 3"="#0072B2",
                     "Cluster 4"="#CC79A7", "Cluster 5"="#D55E00")
p1 <- ggplot(long, aes(cluster,value,fill=cluster)) + geom_boxplot(width=.62, outlier.shape=NA, alpha=.78) +
  geom_jitter(width=.14, size=.75, alpha=.45, color="black") + facet_wrap(~panel, scales="free_y", ncol=2) +
  scale_fill_manual(values=cluster_palette) + labs(title="Clinical variables across transcriptomic clusters", x=NULL, y=NULL) +
  theme_bw(base_size=12) + theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"),
      axis.text.x=element_text(angle=25,hjust=1), strip.text=element_text(face="bold",color="black"), plot.title=element_text(face="bold"))
ggsave(file.path(out,"figures","clinical_continuous_by_cluster.png"), p1, width=9, height=11, dpi=450, bg="white")

cat_pretty <- c(sex="Sex", rf_status="Rheumatoid factor status", ccp_status="Anti-CCP status")
clab <- setNames(sprintf("%s\nFisher p = %.3g", cat_pretty[cat_tests$variable], cat_tests$p_value), cat_tests$variable)
catsumm$panel <- factor(clab[catsumm$variable], levels=clab[categorical])
p2 <- ggplot(catsumm, aes(cluster, percent, fill=category)) + geom_col(color="white", linewidth=.25) +
  facet_wrap(~panel, ncol=1) + scale_y_continuous(labels=function(x) paste0(x,"%"), limits=c(0,100)) +
  labs(title="Categorical clinical variables across transcriptomic clusters", x=NULL, y="Samples", fill=NULL) +
  theme_bw(base_size=12) + theme(text=element_text(color="black"), axis.text=element_text(color="black"),
      strip.text=element_text(face="bold",color="black"), plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(out,"figures","clinical_categorical_by_cluster.png"), p2, width=8.2, height=8.5, dpi=450, bg="white")

cat("Matched samples:", nrow(dat), "\n")
print(tests)
