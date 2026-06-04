# Tabla de librerías utilizadas

**Entorno:** R 4.5.2 (2025-04-11) — "How About a Twenty-Six"  
**Sistema:** Windows 11 x86_64  
**Script principal:** `STRAP_pipeline_main.R` + `tutor_followup_analysis.R`

---

## Paquetes CRAN

| Paquete | Versión | Fuente | Módulo | Función en el análisis |
|---------|---------|--------|--------|------------------------|
| **boot** | 1.3.32 | CRAN | 04_validation | Infraestructura de remuestreo para validación de estabilidad bootstrap |
| **cluster** | 2.1.8.2 | CRAN | 04_validation | Cálculo del coeficiente de silueta (`silhouette()`) para evaluar la cohesión de los clusters |
| **dendextend** | 1.19.1 | CRAN | 03_clustering | Visualización y coloración de dendrogramas jerárquicos |
| **diffusionMap** | 1.2.0 | CRAN | 02_dimreduction | Reducción dimensional no lineal mediante Diffusion Maps (`diffuse()`, `epsilonCompute()`); alternativa a UMAP compatible con R 4.5 |
| **dplyr** | 1.2.1 | CRAN | 06_visualization, follow-up | Manipulación de tablas de datos (`count()`, `filter()`, `mutate()`) |
| **enrichR** | 3.4 | CRAN | tutor_followup | Análisis de enriquecimiento funcional vía API de Enrichr (GO, KEGG, MSigDB Hallmark) |
| **factoextra** | 2.0.0 | CRAN | 02_dimreduction | Visualización de resultados de PCA y clustering |
| **ggalluvial** | 0.12.6 | CRAN | 06_visualization | Diagrama de Sankey (`geom_alluvium()`, `geom_stratum()`) para comparar clusters con patotipos histológicos |
| **ggplot2** | 4.0.3 | CRAN | 06_visualization, follow-up | Sistema principal de visualización; todos los gráficos del pipeline |
| **ggrepel** | 0.9.8 | CRAN | 06_visualization | Etiquetas de texto sin solapamiento en gráficos de dispersión |
| **igraph** | 2.3.1 | CRAN | 03_clustering | Construcción del grafo k-NN (`graph_from_adjacency_matrix()`), algoritmo de Leiden (`cluster_leiden()`) y Markov Clustering |
| **kernlab** | 0.9.33 | CRAN | 03_clustering | Spectral Clustering (`specc()`) basado en kernel RBF |
| **MCL** | 1.0 | CRAN | 03_clustering | Markov Clustering Algorithm (`mcl()`): propagación de flujo sobre grafo de similitud |
| **mclust** | 6.1.2 | CRAN | 03_clustering, 04_validation | Gaussian Mixture Models (`Mclust()`); cálculo del Adjusted Rand Index (`adjustedRandIndex()`) |
| **networkD3** | 0.4.1 | CRAN | — | Diagramas de red interactivos (incluido en setup; en el pipeline final se usa `ggalluvial`) |
| **NMF** | 0.28 | CRAN | 03_clustering | Non-negative Matrix Factorization (`nmf()`): descomposición en k programas transcriptómicos (metagenes) |
| **patchwork** | 1.3.2 | CRAN | 06_visualization, follow-up | Composición de múltiples gráficos ggplot2 en un único panel (`plot_layout()`) |
| **pheatmap** | 1.0.13 | CRAN | 06_visualization, follow-up | Heatmaps con clustering jerárquico y anotación de columnas |
| **RColorBrewer** | 1.1.3 | CRAN | 05_annotation, 06_visualization | Paletas de color cualitativamente distinguibles (`brewer.pal()`) |
| **reshape2** | 1.4.5 | CRAN | 06_visualization | Conversión de matrices a formato largo (`melt()`) para boxplots de marcadores |
| **scales** | 1.4.0 | CRAN | 06_visualization | Formato de ejes y escalas en ggplot2 |
| **tibble** | 3.3.1 | CRAN | General | Estructura de datos tabular moderna (reemplaza `data.frame` en contextos tidyverse) |
| **tidyr** | 1.3.2 | CRAN | General | Pivotado y reestructuración de tablas (`pivot_longer()`, `pivot_wider()`) |
| **viridis** | 0.6.5 | CRAN | 06_visualization | Paletas de color perceptualmente uniformes y accesibles para daltónicos |

---

## Paquetes Bioconductor

| Paquete | Versión | Fuente | Módulo | Función en el análisis |
|---------|---------|--------|--------|------------------------|
| **BiocParallel** | 1.44.0 | Bioconductor | General | Gestión de paralelización para paquetes Bioconductor |
| **clusterProfiler** | 4.18.4 | Bioconductor | Follow-up | Enriquecimiento funcional GO y KEGG (`enrichGO()`, `enrichKEGG()`) |
| **ConsensusClusterPlus** | 1.74.0 | Bioconductor | 03_clustering | Consensus Clustering: remuestrea iterativamente submuestras y subconjuntos de genes para obtener una partición robusta y reproducible (`ConsensusClusterPlus()`) |
| **DESeq2** | 1.50.2 | Bioconductor | 01_preprocessing, 07_deg | Normalización VST de counts crudos (`vst()`); análisis de expresión diferencial one-vs-rest (`DESeq()`, `results()`) |
| **fgsea** | 1.36.2 | Bioconductor | Follow-up | Gene Set Enrichment Analysis rápido (`fgsea()`) sobre rankings de genes diferenciales |
| **limma** | 3.66.0 | Bioconductor | Follow-up | Análisis de expresión diferencial moderado con modelos lineales (`lmFit()`, `eBayes()`) |
| **org.Hs.eg.db** | 3.22.0 | Bioconductor | Follow-up | Base de datos de anotación genómica humana; conversión de IDs entre Entrez, Ensembl y símbolos HGNC |

---

## Entorno base de R

| Componente | Versión | Función |
|------------|---------|---------|
| **R** | 4.5.2 | Lenguaje y entorno de computación estadística |
| **stats** (base) | 4.5.2 | `prcomp()` (PCA), `kmeans()`, `hclust()`, `dist()`, `cutree()` |
| **base** (base) | 4.5.2 | Operaciones matriciales, lectura de datos, gestión de entorno |
| **grDevices** (base) | 4.5.2 | Dispositivos gráficos, `colorRampPalette()` |

---

## Resumen por módulo

| Módulo | Archivo | Paquetes principales |
|--------|---------|----------------------|
| Setup | `00_setup.R` | Todos los anteriores (instalación y carga) |
| Preprocesado | `01_preprocessing.R` | **DESeq2** |
| Reducción dimensional | `02_dimreduction.R` | **diffusionMap**, stats (prcomp) |
| Clustering | `03_clustering.R` | **kernlab**, **igraph**, **MCL**, **ConsensusClusterPlus**, **NMF**, **mclust** |
| Validación | `04_validation.R` | **cluster**, **mclust**, **boot** |
| Anotación | `05_annotation.R` | **RColorBrewer** |
| Visualización | `06_visualization.R` | **ggplot2**, **pheatmap**, **ggalluvial**, **patchwork**, **reshape2**, **RColorBrewer**, **viridis**, **scales** |
| Expresión diferencial | `07_deg_analysis.R` | **DESeq2** |
| Exportación | `08_export.R` | Base R |
| Seguimiento tutor | `tutor_followup_analysis.R` | **enrichR**, **ggplot2**, **pheatmap**, **patchwork**, **dplyr**, **DESeq2**, **clusterProfiler**, **fgsea** |
