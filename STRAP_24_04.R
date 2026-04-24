
# 1). Carga y exploración de datos 

load("strap_counts.RData")

dim(strap_counts)
colnames(strap_counts)
head(strap_counts)

# Renombrar strap_counts
expr <- strap_counts

#################################################################################

# 2). Procesado de datos 

  # Quitar genes con baja varianza y seleccionar los 1000 más           significativos 

var_genes <- apply(expr, 1, var)

top_genes <- names(sort(var_genes, decreasing = TRUE))[1:1000]

expr_filt <- expr[top_genes, ]


# Normalización (VST) ES mejor para RNA-seq

library(DESeq2)

   # Crear objeto DESeq (sin diseño porque es un modelo no supervisado)
dds <- DESeqDataSetFromMatrix(countData = expr_filt,
                              colData = data.frame(row.names =   colnames(expr_filt)),
                              design = ~1)

   # Aplicar VST
vsd <- vst(dds, blind = TRUE)

expr_vst <- assay(vsd)


 # Transponer samples en filas

expr_t <- t(expr_vst)

 # Escalado 

expr_scaled <- scale(expr_t)

############################################################################

# 3). Reducción de dimensionalidad y Visualización previa de los datos 

  # PCA ÚNICO 

pca <- prcomp(expr_scaled)

plot(pca$x[,1], pca$x[,2],
     xlab="PC1", ylab="PC2",
     main="PCA RNA-seq (VST)",
     pch=19, col="grey")
# No consigo visualizar ninguna separación natural de patotipos

##############################################################

# 4). Selección de k antes del clustering 

library(factoextra)

fviz_nbclust(expr_scaled, kmeans, method = "gap_stat")

k <- 3 

  # El método Gap sugiere k = 3, aunque el incremento progresivo del estadístico indica una estructura de clustering débil y posiblemente continua, lo cual es coherente con los datos que están siendo tratados, ya que es una enfermedad heterogénea contínua. 

##############################################################################################################################################

# 5). Clustering 

  # 5.1). k-means

set.seed(123)
km <- kmeans(expr_scaled, centers = k)

cluster_colors <- rainbow(k)[km$cluster]

  # 5.2). Jerárquico 

dist_mat <- dist(expr_scaled)
hc <- hclust(dist_mat, method = "ward.D2")

hc_clusters <- cutree(hc, k = k)

  # 5.3). Densidad
library(dbscan)

hdb <- hdbscan(expr_scaled, minPts = 10)

hdb_clusters <- hdb$cluster

table(hdb_clusters)

  # 5.4. Espectral 

library(kernlab)

set.seed(123)

spec <- specc(expr_scaled, centers = k)
spec_clusters <- as.numeric(spec)
table(spec_clusters)

######################################################################################################################################

# 6). Visualización de cluster + PCA

  # 6.1). K-means + PCA

plot(pca$x[,1], pca$x[,2],
     col = cluster_colors,
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "Clusters (K-means) en PCA")

legend("bottomright",
       legend = paste("Cluster", 1:k),
       col = rainbow(k),
       pch = 19)

  # 6.2). Jerárquico + PCA
cluster_colors_hc <- rainbow(k)[hc_clusters]

plot(pca$x[,1], pca$x[,2],
     col=cluster_colors_hc, pch=19,
     main="PCA con clusters jerárquicos")

  # 6.3). Densidad + pca

plot(pca$x[,1], pca$x[,2],
     col = hdb_clusters + 1,  # +1 para evitar 0 (ruido)
     pch = 19,
     main = "HDBSCAN en PCA",
     xlab = "PC1", ylab = "PC2")

  # 6.4). Espectral + pca

plot(pca$x[,1], pca$x[,2],
     col = spec_clusters,
     pch = 19,
     main = "Spectral Clustering en PCA",
     xlab = "PC1", ylab = "PC2")

##########################################################################################################################

# 7). UMAP 

library(umap)

umap_res <- umap(expr_scaled)

plot(umap_res$layout,
     col=cluster_colors, pch=19, cex=1.5,
     main="UMAP clusters",
     xlab="UMAP1", ylab="UMAP2")

##########################################################################################################################

# 8). Heatmap

library(pheatmap)

annotation <- data.frame(Cluster = as.factor(km$cluster))
rownames(annotation) <- rownames(expr_scaled)

pheatmap(t(expr_scaled),
         annotation_col = annotation,
         show_rownames = FALSE,
         main = "Heatmap RNA-seq (VST)")


###############################################################################

# 9). Evaluación sin recalcular nada 

library(cluster)

sil <- silhouette(km$cluster, dist(expr_scaled))
plot(sil)

######################################################################################################

# 10). COMPARACIÓN CON PATOTIPOS

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

# Alinear con expr scaled 
# rownames de strap_counts deben coincidir con `Characteristics.sampleid.`
common_samples <- intersect(rownames(expr_scaled), rownames(meta_unique))

expr_scaled_filt <- expr_scaled[common_samples, ]
meta_filt <- meta_unique[common_samples, ]

# Clusters que ya calculé
clusters <- km$cluster[match(common_samples, rownames(expr_scaled_filt))]

   # Limpiar patotipos 
meta_filt$pathotype_clean <- meta_filt$`Characteristics.pathotype.`

   # Unir Fibrous → Fibroid
meta_filt$pathotype_clean[
  meta_filt$pathotype_clean == "Fibrous"
] <- "Fibroid"

   # Eliminar Ungraded
keep <- meta_filt$pathotype_clean != "Ungraded"

meta_filt <- meta_filt[keep, ]
expr_scaled <- expr_scaled[keep, ]
km$cluster <- km$cluster[keep]

###############################################################################

# Evaluación biológica

  # Unión
# 1.pca con patotipos
pca <- prcomp(expr_scaled)

plot(pca$x[,1], pca$x[,2],
     col = as.factor(meta_filt$pathotype_clean),
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA coloreado por patotipo")

legend("topright",
       legend = levels(as.factor(meta_filt$pathotype_clean)),
       col = 1:length(levels(as.factor(meta_filt$pathotype_clean))),
       pch = 19)


# 2.Matriz de confusión con patotipo

clusters <- km$cluster 

table_clusters <- table(clusters, meta_filt$pathotype_clean)
table_clusters


# Opcional: visualización en heatmap
library(pheatmap)

pheatmap(as.matrix(table_clusters),
         color = colorRampPalette(c("white", "blue"))(50),
         main = "Clusters vs Patotipos",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE)

# 2. Sankey 
library(networkD3)

df <- data.frame(
  cluster = as.factor(clusters),
  pathotype = as.factor(meta_filt$pathotype_clean)
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

sil <- silhouette(km$cluster, dist(expr_scaled))
plot(sil)
