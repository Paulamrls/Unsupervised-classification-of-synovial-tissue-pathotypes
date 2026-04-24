# Análisis de Clustering No Supervisado de Patotipos en Artritis Reumatoide --- STRAP RNA-seq

## 1. Introducción

La artritis reumatoide (AR) es una enfermedad autoinmune crónica
caracterizada por una elevada heterogeneidad biológica, tanto a nivel
clínico como molecular. En el tejido sinovial, esta heterogeneidad se
traduce en distintos patotipos histológicos que reflejan mecanismos
inmunológicos diferenciados.

Los principales patotipos descritos en la literatura son:

-   **Lymphoid**: caracterizado por infiltración de células B y T
    organizadas, con expresión de genes como *CXCL13*, *CD79A* o
    *MS4A1*.
-   **Myeloid**: dominado por la respuesta inflamatoria innata, con alta
    expresión de *IL1B*, *TNF* y *S100A8/A9*.
-   **Fibroid**: asociado a fibroblastos y remodelación de la matriz
    extracelular (*COL1A1*, *FN1*, *THY1*).

El objetivo de este trabajo es evaluar si estos patotipos pueden ser
identificados mediante técnicas de clustering no supervisado sobre datos
de RNA-seq.

------------------------------------------------------------------------

## 2. Datos y preprocesamiento

### 2.1 Datos

-   Dataset: STRAP
-   Tipo: RNA-seq bulk de tejido sinovial
-   Muestras: \~210

### 2.2 Preprocesamiento

-   Transformación: VST (Variance Stabilizing Transformation)
-   Selección de genes: alta varianza
-   Escalado: aplicado para métodos basados en distancia

### Interpretación biológica

El uso de RNA-seq bulk implica que cada muestra representa una mezcla de
múltiples tipos celulares. Esto introduce una limitación fundamental:
los patotipos no se observan de forma pura, sino como combinaciones.

------------------------------------------------------------------------

## 3. Metodología de clustering

Se aplicaron los siguientes métodos:

-   K-means
-   Clustering jerárquico (Ward)
-   Spectral clustering
-   HDBSCAN

### Limitaciones

-   K-means asume clusters esféricos
-   Distancia euclídea no captura bien la coexpresión génica
-   No se incorpora conocimiento biológico

------------------------------------------------------------------------

## 4. Reducción dimensional

### PCA

El PCA muestra una distribución continua sin separaciones claras.

Interpretación biológica: la variabilidad capturada refleja gradientes
biológicos, no subgrupos discretos.

### UMAP

UMAP revela separaciones locales más claras, pero puede introducir
artefactos visuales.

------------------------------------------------------------------------

## 5. Interpretación biológica de resultados

### Clustering

Se observa:

-   Un cluster bien definido asociado a Lymphoid
-   Dos clusters con fuerte solapamiento (Myeloid y Fibroid)

Esto sugiere que:

-   El patotipo Lymphoid tiene una firma transcriptómica más clara
-   Myeloid y Fibroid comparten rutas inflamatorias

------------------------------------------------------------------------

## 6. Evaluación

Silhouette score ≈ 0.19

Interpretación:

-   Clusters débiles
-   Alta ambigüedad en asignación
-   Evidencia de estructura continua

------------------------------------------------------------------------

## 7. Comparación con patotipos

-   Lymphoid: bien capturado
-   Myeloid: disperso
-   Fibroid: parcialmente identificado

------------------------------------------------------------------------

## 8. Discusión

### Continuum biológico

Los patotipos no son discretos, sino extremos de un espectro.

### RNA-seq bulk

La señal se diluye por mezcla celular.

### Alta dimensionalidad

Reduce la capacidad de clustering.

------------------------------------------------------------------------

## 9. Propuestas de mejora

-   Firmas génicas específicas
-   Deconvolución celular
-   Modelos probabilísticos
-   Clustering de consenso
-   Integración multimodal

------------------------------------------------------------------------

## 10. Conclusión

El clustering no supervisado permite capturar parcialmente la estructura
biológica, especialmente el patotipo Lymphoid, pero falla en separar
completamente los subtipos debido a limitaciones biológicas y técnicas.
