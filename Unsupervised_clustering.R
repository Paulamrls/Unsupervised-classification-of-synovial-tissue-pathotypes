###############################################################################
###############################################
# clustering no supervisado
# Prueba 
# Thu Mar 12 18:16:08 2026
# TFM AR
###############################################
##############################################################################
 
# Directorio 
setwd("C:/Users/Paula/Desktop/bioinformatica/TFM/codigo")

# Archivos 

list.files()

# Cargar matriz de datos normalizados 

expr <- read.table("logCPM_matrix.tsv",
                   header = TRUE,
                   sep = "\t",
                   row.names = 1,
                   check.names = FALSE)



# 1). Preparar la matriz antes del clustering no supervisado

 # Selecciono los 1000 genes más variables para reducir ruido y mejorar el resultado de los clusters

var_genes <- apply(expr, 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:1000]
expr_filt <- expr[top_genes, ]

#=============================================================================

# 2).VISUAL INICIAL 

#    2.1. PCA para ver patrones generales (LINEAL)

pca <- prcomp(t(expr_filt), scale. = TRUE)

# Gráfico de PC1 vs PC2

plot(pca$x[,1], pca$x[,2],
     xlab="PC1", ylab="PC2",
     main="PCA - exploración inicial",
     pch=19,           # punto sólido
     col="red",  # color visible
     cex=1.5)

# Muy útil para ver si hay separación natural entre muestras
# Detectar outliers
# No hace clustering, solo visualización

# A pesar de las muy buenas propiedades que tiene el PCA, sufre de algunas limitaciones, por ejemplo, solo tiene en cuenta combinaciones lineales de las variables originales. En determinados escenarios, el no poder considerar otro tipo de combinaciones supone perder mucha información.


  #    2.2. t-SNE para reducir variables (NO LINEAL)
install.packages("Rtsne")
library(Rtsne)

# También se limita el número de iteraciones (epoch) a 100, aunque los
# resultados podrían mejorar si se aumentara

expr_t <- t(expr_filt)
set.seed(123)  # reproducibilidad
tsne_res <- Rtsne(expr_t, 
                  dims = 2,      # dos dimensiones
                  perplexity = 30,  # ajustar según número de samples
                  verbose = TRUE, 
                  max_iter = 1000)

# Hice primero el k-means, después descubrí Rt-SNE

cluster_colors <- rainbow(k)[km$cluster]
plot(tsne_res$Y, pch=19, col=cluster_colors, cex=1.5,
      xlab="t-SNE 1", ylab="t-SNE 2", main="t-SNE con clusters")
 legend("topright", legend=paste("Cluster", 1:k), col=rainbow(k), pch=19)
 
 
 
#       2.3. UMAP (SUPERA A T-SNE)
 
 install.packages("umap")
 library(umap)
 
 umap_res <- umap(t(expr_filt))
 
 plot(umap_res$layout,
      col=cluster_colors, pch=19, cex=1.5,
      main="UMAP con k-means clusters",
      xlab="UMAP1", ylab="UMAP2")
 legend("topleft", legend=paste("Cluster", 1:k),
        col=rainbow(k), pch=19)
 
 # Captura estructuras no lineales
 # Muy útil para visualizar clusters complejos
 # Complementa PCA y clustering jerárquico
 
 
 
 #     2.4. Deep learning
 
 #     2.5. CIDR
 
 #     2.6. seurat (single cell)
#==============================================================================

# CLUSTERING

  # 3).Hierarchical clustering (dendrograma)

dist_mat <- dist(t(expr_filt))
hc <- hclust(dist_mat, method = "ward.D2")

plot(hc, main="Hierarchical Clustering - Ward", xlab="", sub="")

# dendrograma en k clusters

k_hc <- 3
hc_clusters <- cutree(hc, k = k_hc)
cluster_colors_hc <- rainbow(k_hc)[hc_clusters]

# Resaltar dendrograma con colores si quieres

install.packages("dendextend")
library(dendextend)
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = k_hc)
plot(dend, main="Hierarchical Clustering coloreado")

# Mostrar estructura jerárquica de los datos
# Fácil de visualizar con dendrogramas
# No requiere decidir número de clusters de inicio (puedes “cortar” el dendrograma)

#=================================================================================

  # 4).k-means clustering

#  tiende a identificar grupos globulares, lo que resulta en fallas en la detección de tipos de células raros.

set.seed(123)
k <- 4
km <- kmeans(t(expr_filt), centers = k)

cluster_colors <- rainbow(k)[km$cluster]

# PCA coloreado por clusters de k-means

plot(pca$x[,1], pca$x[,2],
     col=cluster_colors, pch=19, cex=1.5,
     xlab="PC1", ylab="PC2",
     main="PCA con k-means clusters")


legend("bottomright", legend=paste("Cluster", 1:k),
       col=rainbow(k), pch=19)

# Simple y rápido
# Funciona bien si los clusters son globulares y bien separados
# Necesita decidir k de antemano

#=================================================================================

  # 5).Consensus clustering (robustez)

BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

results <- ConsensusClusterPlus(
  as.matrix(expr_filt),
  maxK = 6,        # probamos hasta 6 clusters
  reps = 100,      # re-muestreo
  pItem = 0.8,     # 80% de samples
  pFeature = 1,    # todos los genes
  clusterAlg = "hc",
  distance = "pearson",
  seed = 123,
  plot = "png"
)

# Evalúa la robustez de los clusters
# Muy usado en subtipos de enfermedad (cáncer, artritis)
# Te ayuda a decidir cuántos clusters son confiables

#=================================================================================
 
   # 6). Densidad DBSCAN (Density-Based Clustering)


install.packages("dbscan")
library(dbscan)

set.seed(123)
# eps = radio de vecindad, minPts = mínimo número de puntos para formar cluster
db <- dbscan(expr_t, eps = 3, minPts = 5)  

# Cluster asignados
db$cluster  # 0 = outlier, 1,2,... = clusters

# Colores
cluster_colors <- rainbow(length(unique(db$cluster)))[db$cluster + 1]

# PCA para visualizar
pca <- prcomp(expr_t, scale. = TRUE)
plot(pca$x[,1], pca$x[,2], col=cluster_colors, pch=19, cex=1.5,
     main="DBSCAN Clusters en PCA", xlab="PC1", ylab="PC2")
legend("topright", legend=paste("Cluster", sort(unique(db$cluster))),
       col=rainbow(length(unique(db$cluster))), pch=19)

# Detecta clusters de forma arbitraria y outliers

#================================================================================


  # 7). NMF (Non-negative Matrix Factorization)

 BiocManager::install("NMF")
library(NMF)

# NMF requiere datos no negativos; si tus datos tienen negativos, usa log-transform + shift
expr_shift <- expr_t - min(expr_t) + 1  # asegura >0

set.seed(123)
# rank = número de clusters que queremos probar
nmf_res <- nmf(expr_shift, rank=3, nrun=30, seed=123)

# Cluster asignados a cada sample
nmf_clusters <- predict(nmf_res)

# Colores
cluster_colors <- rainbow(3)[nmf_clusters]

# Visualización PCA
plot(pca$x[,1], pca$x[,2], col=cluster_colors, pch=19, cex=1.5,
     main="NMF Clusters en PCA", xlab="PC1", ylab="PC2")
legend("topright", legend=paste("Cluster", 1:3), col=rainbow(3), pch=19)

# Muy usado para detectar patotipos
# Patrones de expresión, metagenes, útil para patotipos biológicos

#===============================================================================

#  8). Spectral Clustering

install.packages("kernlab")
library(kernlab)

set.seed(123)
# k = número de clusters
k <- 4
sc <- specc(expr_t, centers = k)

spectral_clusters <- sc@.Data

# Colores
cluster_colors <- rainbow(k)[spectral_clusters]

# Visualización PCA
plot(pca$x[,1], pca$x[,2], col=cluster_colors, pch=19, cex=1.5,
     main="Spectral Clustering en PCA", xlab="PC1", ylab="PC2")
legend("topright", legend=paste("Cluster", 1:k), col=rainbow(k), pch=19)


# Detecta clusters complejos que no sean globulares


# Pienso que el espectral, el NMF, la densidad, k means son iguales. Igual porque todos colorean el PCA.

#===============================================================================

# Una forma de visualizar los resultados, pero creo que no es cluster. Volcano también se puede usar.

#  Heatmap de expresión

install.packages("pheatmap")
library(pheatmap)

pheatmap(expr_filt,
         scale = "row",           # estandariza genes
         cluster_cols = TRUE,     # agrupa muestras
         cluster_rows = TRUE,     # agrupa genes
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap de genes variables")

# Visualizar patrones de expresión por cluster
# Muy útil para comunicar resultados a tu tutor/paper



# When clustering the TCGA-BRCA matrix, three natural groups consistently emerge that reflect intrinsic breast cancer subtypes. This helps me understand how the most variable genes can generate robust clusters and prepares me to see if the same occurs in RA synovial pathotypes.
