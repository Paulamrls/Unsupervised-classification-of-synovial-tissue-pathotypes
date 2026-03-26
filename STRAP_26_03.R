
# Carga y exploración de datos 

load("strap_counts.RData")

dim(strap_counts)
colnames(strap_counts)
head(strap_counts)

# Renombrar strap_counts
expr <- strap_counts

#################################################################################

# Quitar genes con baja varianza y seleccionar los 1000 más significativos 
var_genes <- apply(expr, 1, var)

top_genes <- names(sort(var_genes, decreasing = TRUE))[1:1000]

expr_filt <- expr[top_genes, ]


# Log-transform para estabilizar varianza
expr_log <- log2(expr_filt + 1)

# Transponer samples en filas
expr_t <- t(expr_log)


# Visualización previa de los datos 

  # pca

pca <- prcomp(expr_t, scale. = TRUE)

plot(pca$x[,1], pca$x[,2],
     xlab="PC1", ylab="PC2",
     main="PCA RNA-seq (STRAP)",
     pch=19,
     col="red",
     cex=1.5)

# cluster no supervisado

  # 1). k means

set.seed(123)

k <- 3 
km <- kmeans(expr_t, centers = k)

cluster_colors <- rainbow(k)[km$cluster]

plot(pca$x[,1], pca$x[,2],
     col=cluster_colors, pch=19, cex=1.5,
     main="PCA con clusters RNA-seq",
     xlab="PC1", ylab="PC2")

legend("bottomright", legend=paste("Cluster", 1:k),
       col=rainbow(k), pch=19)


  # 2. Jerárquico 

# Usando expr_filt (los 1000 genes más variables)
expr_scaled <- scale(expr_filt)
expr_scaled_t <- t(expr_scaled)

dist_mat <- dist(expr_scaled_t, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")

plot(hc,
     main = "Hierarchical Clustering - Ward",
     xlab = "",
     sub = "",
     cex = 0.7)
k_hc <- 3  # Ajusta según tu criterio
hc_clusters <- cutree(hc, k = k_hc)
table(hc_clusters)

library(dendextend)

dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = k_hc)

plot(dend,
     main = "Hierarchical Clustering coloreado",
     ylab = "Altura")

pca <- prcomp(expr_scaled_t, scale. = TRUE)

cluster_colors_hc <- rainbow(k_hc)[hc_clusters]

plot(pca$x[,1], pca$x[,2],
     col = cluster_colors_hc,
     pch = 19,
     cex = 1.5,
     main = "PCA coloreado por clusters jerárquicos",
     xlab = "PC1", ylab = "PC2")

legend("bottomleft",
       legend = paste("Cluster", 1:k_hc),
       col = rainbow(k_hc), pch = 19)

  # 3. umap

library(umap)

umap_res <- umap(expr_t)

plot(umap_res$layout,
     col=cluster_colors, pch=19, cex=1.5,
     main="UMAP RNA-seq clusters",
     xlab="UMAP1", ylab="UMAP2")


    #  Visualización de los resultados: heatmap

library(pheatmap)

annotation <- data.frame(Cluster = as.factor(km$cluster))
rownames(annotation) <- colnames(expr_log)

pheatmap(expr_log,
         scale = "row",
         annotation_col = annotation,
         show_rownames = FALSE,
         main = "Heatmap RNA-seq")


###############################################################################

# COMPARACIÓN CON PATOTIPOS

# Quitar duplicados del fasq
list.files()
meta <- read.delim("E-MTAB-13733.sdrf.txt")
colnames(meta)
SampleID <- "Characteristics.sampleid."

# Elegimos columnas de interés
cols_keep <- c("Characteristics.sampleid.",
               "Characteristics.pathotype.",
               "Characteristics.tender_joints_counts_tjc.",
               "Characteristics.swollen_joints_counts_sjc.",
               "Characteristics.arthritis.activity.",
               "Characteristics.esr.",
               "Characteristics.crp.",
               "Characteristics.physicians.global.assessment.")

meta_sub <- meta[, cols_keep]

# Colapsar duplicados por sample
meta_unique <- meta_sub[!duplicated(meta_sub$`Characteristics.sampleid.`), ]
rownames(meta_unique) <- meta_unique$`Characteristics.sampleid.`

# Ahora meta_unique tiene **una fila por muestra**
dim(meta_unique)
head(meta_unique)

###############################################################################

# Alinear con expr t 
# rownames de strap_counts deben coincidir con `Characteristics.sampleid.`
common_samples <- intersect(rownames(expr_t), rownames(meta_unique))

expr_t_filt <- expr_t[common_samples, ]
meta_filt <- meta_unique[common_samples, ]

# Clusters que ya calculé
clusters <- km$cluster[match(common_samples, rownames(expr_t_filt))]


###############################################################################

# Evaluación biológica

# 1.pca con patotipos
plot(pca$x[,1], pca$x[,2],
     col=as.factor(meta_filt$`Characteristics.pathotype.`),
     pch=19,
     main="PCA coloreado por patotipo")

# Matriz de confusión con patotipo
table_clusters <- table(clusters, meta_filt$`Characteristics.pathotype.`)

# Opcional: visualización en heatmap
library(pheatmap)
pheatmap(as.matrix(table_clusters),
         color = colorRampPalette(c("white", "blue"))(50),
         main = "Clusters vs Patotipos",
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE)

# 2. Sankey 

install.packages("networkD3")
library(networkD3)

df <- data.frame(
  cluster = as.factor(clusters),
  pathotype = as.factor(meta_filt$`Characteristics.pathotype.`)
)

links <- as.data.frame(table(df))
nodes <- data.frame(name = unique(c(as.character(links$cluster),
                                    as.character(links$pathotype))))

links$IDsource <- match(links$cluster, nodes$name) - 1
links$IDtarget <- match(links$pathotype, nodes$name) - 1

sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "Freq", NodeID = "name")

################################################################################

# Evaluación de métodos 
library(cluster)
# recalcular todo limpio
expr_t <- t(log2(expr_filt + 1))

km <- kmeans(expr_t, centers = k)

dist_mat <- dist(expr_t)

sil <- silhouette(km$cluster, dist_mat)

plot(sil)
