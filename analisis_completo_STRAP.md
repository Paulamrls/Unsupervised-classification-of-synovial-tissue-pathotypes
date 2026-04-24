# Análisis de Clustering No Supervisado de Patotipos en Artritis Reumatoide 24/04

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
