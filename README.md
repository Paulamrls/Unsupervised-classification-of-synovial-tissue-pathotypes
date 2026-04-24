# Análisis de Clustering No Supervisado de Patotipos en Artritis Reumatoide — Datos STRAP RNA-seq

## 1. Descripción del análisis

### 1.1 Tipo de datos

El análisis parte de datos de expresión génica obtenidos mediante RNA-seq del estudio STRAP (*Stratification of Biologic Therapies for Rheumatoid Arthritis by Pathobiology*), asociado al repositorio ArrayExpress **E-MTAB-13733**. La matriz de conteos utilizada contiene muestras de tejido sinovial de pacientes con artritis reumatoide (AR). El objetivo del análisis es explorar si los perfiles transcriptómicos permiten identificar, mediante técnicas de aprendizaje no supervisado, subgrupos de pacientes compatibles con los patotipos histopatológicos descritos en la enfermedad: **Lymphoid**, **Myeloid** y **Fibroid/Fibrous**.

A diferencia de un enfoque supervisado, en este análisis los algoritmos de clustering no reciben como entrada la etiqueta de patotipo. Las etiquetas histopatológicas se emplean posteriormente como referencia externa para evaluar si los grupos generados de forma automática tienen coherencia biológica.

### 1.2 Preprocesamiento de los datos

El pipeline implementado en el archivo `STRAP_24_04.R` sigue una estrategia más robusta que una simple transformación logarítmica, ya que incorpora normalización mediante **VST** de DESeq2. Los pasos principales son los siguientes:

1. **Carga de la matriz de conteos**: se carga el objeto `strap_counts.RData` y se renombra como `expr` para trabajar con una matriz de genes por muestras.
2. **Selección de genes por varianza**: se calcula la varianza de cada gen a lo largo de todas las muestras y se seleccionan los **1.000 genes con mayor varianza**, con el objetivo de conservar los genes con mayor capacidad discriminativa.
3. **Normalización VST**: se crea un objeto `DESeqDataSetFromMatrix` con diseño `~1`, al tratarse de un análisis no supervisado, y se aplica la transformación `vst(dds, blind = TRUE)`. Esta transformación estabiliza la varianza y es más adecuada para RNA-seq que una transformación log₂ simple.
4. **Transposición de la matriz**: la matriz se transpone para que las filas representen muestras y las columnas genes.
5. **Escalado z-score**: se aplica `scale(expr_t)` para centrar y escalar cada variable. Esto es importante porque métodos como K-means, clustering jerárquico, PCA y clustering espectral son sensibles a la escala de las variables.

Este preprocesamiento es metodológicamente más consistente que aplicar clustering directamente sobre conteos crudos o sobre datos log-transformados sin normalización específica para RNA-seq.

### 1.3 Métodos de reducción dimensional y visualización

Se utilizaron dos enfoques principales de visualización:

- **PCA (Análisis de Componentes Principales)**: aplicado sobre la matriz normalizada y escalada (`expr_scaled`) mediante `prcomp(expr_scaled)`. El PCA permite observar la estructura global de la variabilidad transcriptómica en las dos primeras componentes principales.
- **UMAP**: aplicado también sobre `expr_scaled` mediante la librería `umap`. UMAP permite representar relaciones locales y no lineales entre muestras, por lo que puede revelar estructuras que no son evidentes en PCA.

Estas visualizaciones no deben interpretarse como prueba definitiva de existencia de clusters, sino como herramientas exploratorias para observar tendencias, solapamientos y posibles outliers.

### 1.4 Selección del número de clusters

Antes de ejecutar los métodos principales de clustering se utilizó el estadístico de **Gap** mediante `fviz_nbclust(expr_scaled, kmeans, method = "gap_stat")`. El gráfico sugiere **k = 3** como solución razonable, lo cual además coincide con el conocimiento biológico previo de los tres patotipos principales de AR.

No obstante, el estadístico Gap muestra un incremento progresivo para valores superiores de k, lo que indica que la estructura de agrupamiento no es completamente nítida. Por tanto, k = 3 debe interpretarse como una elección biológicamente razonable y operativa, no como una prueba de que existan tres clusters perfectamente separados.

### 1.5 Métodos de clustering aplicados

El script evalúa cuatro métodos de clustering:

- **K-means** con `centers = 3`, aplicado sobre la matriz VST escalada. Es el método principal utilizado para la comparación con patotipos.
- **Clustering jerárquico Ward.D2**, basado en distancia euclídea y corte del dendrograma en k = 3.
- **HDBSCAN**, método basado en densidad, con `minPts = 10`, utilizado para comprobar si existen agrupaciones densas naturales sin imponer k.
- **Spectral Clustering**, mediante `specc()` de `kernlab`, configurado con k = 3.

### 1.6 Evaluación

La evaluación se realiza desde dos perspectivas:

1. **Evaluación interna**, mediante silhouette score para la solución K-means.
2. **Evaluación biológica externa**, mediante comparación de los clusters K-means con los patotipos anotados en los metadatos STRAP.

Para la comparación con patotipos, el script realiza una limpieza de etiquetas: une **Fibrous** dentro de **Fibroid** y elimina las muestras **Ungraded**, ya que no representan un patotipo histológico definido.

---

## 2. Interpretación biológica de cada gráfico

### Figura 1 — PCA previo sobre datos RNA-seq normalizados con VST

![PCA RNA-seq VST](img/pca_previo_vst.png)

Esta figura muestra la proyección inicial de las muestras sobre las dos primeras componentes principales, sin colorear por cluster ni por patotipo. El objetivo es comprobar si existe una separación natural evidente en los datos antes de aplicar algoritmos de agrupamiento.

La nube de puntos aparece distribuida de forma continua, sin fronteras claras entre grupos. La mayoría de las muestras se concentran en una región central-superior del espacio PCA, mientras que algunas muestras se alejan hacia valores negativos de PC2, actuando como observaciones más extremas. Esta distribución indica que la variabilidad transcriptómica dominante no forma grupos discretos fácilmente separables, sino un gradiente continuo de expresión.

Este resultado es coherente con la naturaleza de la artritis reumatoide, una enfermedad inflamatoria heterogénea en la que los perfiles sinoviales pueden presentar características mixtas. Desde el punto de vista metodológico, esta figura anticipa que cualquier algoritmo de clustering tendrá dificultades para encontrar grupos completamente definidos.

### Figura 2 — Selección del número óptimo de clusters mediante Gap Statistic

![Gap Statistic](img/gap_statistic.png)

El gráfico del estadístico Gap evalúa diferentes valores de k para estimar cuántos clusters podrían representar mejor la estructura de los datos. La línea discontinua vertical marca la selección de **k = 3**.

La elección de k = 3 resulta razonable porque coincide con los tres patotipos principales de referencia: Lymphoid, Myeloid y Fibroid. Sin embargo, el gráfico también muestra que el valor del Gap no alcanza una meseta clara en k = 3, sino que continúa aumentando para valores superiores, especialmente hacia k = 8–10. Esto sugiere que los datos podrían contener una estructura más continua o subestructuras internas, en lugar de tres grupos perfectamente separados.

Por tanto, el uso de k = 3 se justifica principalmente por coherencia biológica y comparabilidad con los patotipos conocidos, pero la propia métrica advierte de que la estructura de clusters es débil.

### Figura 3 — PCA con clusters K-means

![K-means PCA](img/kmeans_pca_v1.png)

Esta figura muestra la proyección PCA coloreada según la asignación de K-means con k = 3. Los clusters aparecen parcialmente diferenciados, pero con solapamiento visible entre ellos.

El cluster rojo se sitúa preferentemente en la zona superior-derecha del gráfico, el cluster azul ocupa regiones más centrales e izquierdas, y el cluster verde se concentra sobre todo en valores positivos de PC1 y negativos de PC2. Aunque existe cierta organización espacial, las fronteras entre clusters no son nítidas. Esto indica que K-means está particionando una distribución continua en tres regiones, más que detectando tres grupos completamente separados.

La solución K-means parece capturar gradientes principales de variación transcriptómica, pero no permite afirmar por sí sola que existan tres subtipos moleculares discretos. En datos bulk RNA-seq, esta situación es esperable, ya que cada muestra representa una mezcla de tipos celulares y estados inflamatorios.

### Figura 4 — PCA con clusters K-means tras filtrado de metadatos

![K-means PCA segunda visualización](img/kmeans_pca_v2.png)

Esta segunda visualización de K-means en PCA presenta una distribución similar, aunque con un rango más limitado en PC2 debido al filtrado posterior de muestras para la comparación con patotipos. El patrón general se mantiene: los clusters muestran una cierta orientación en el espacio PCA, pero no una separación completamente clara.

La persistencia del solapamiento refuerza la interpretación de que los clusters representan tendencias transcriptómicas parciales. El cluster rojo agrupa una región superior, el azul ocupa una zona amplia hacia PC1 negativo y el verde se concentra en una región inferior-derecha. Esta distribución puede ser útil para describir subgrupos, pero debe interpretarse con cautela porque el PCA solo resume una parte de la variabilidad total de los 1.000 genes seleccionados.

### Figura 5 — UMAP coloreado por clusters K-means

![UMAP clusters](img/umap_clusters.png)

El UMAP muestra una separación visual más clara que el PCA. El cluster rojo aparece principalmente en la zona superior y derecha, el cluster azul se localiza en regiones inferiores y centrales, y el cluster verde se concentra especialmente en la zona izquierda.

Esta visualización sugiere que K-means ha capturado cierta estructura local de los datos. Sin embargo, UMAP puede amplificar separaciones locales y generar agrupamientos visualmente más definidos de lo que realmente existe en el espacio original de alta dimensión. Por ello, esta figura debe leerse junto con el silhouette score y la comparación con patotipos.

Biológicamente, la distribución en UMAP sugiere que algunas muestras comparten perfiles transcriptómicos más próximos, especialmente el grupo verde de la zona izquierda y el grupo rojo de la zona superior-derecha. Aun así, la presencia de zonas de transición y puntos intermedios indica que la separación entre perfiles no es absoluta.

### Figura 6 — Silhouette plot de K-means

![Silhouette plot](img/silhouette_kmeans.png)

El gráfico de silueta evalúa la calidad de la asignación de cada muestra a su cluster. En este análisis se obtiene una **silhouette media global de 0.19**, con los siguientes valores por cluster:

| Cluster | n | Silhouette media |
|---|---:|---:|
| 1 | 109 | 0.26 |
| 2 | 31 | 0.18 |
| 3 | 70 | 0.07 |

El resultado global de 0.19 indica una estructura de clustering débil. Aunque el cluster 1 presenta una cohesión moderadamente superior al resto, el cluster 3 tiene un valor muy bajo, próximo a 0, lo que significa que muchas de sus muestras se encuentran cerca de la frontera con otros clusters.

Desde el punto de vista estadístico, estos valores no permiten considerar la partición como robusta. En términos biológicos, el resultado sugiere que los pacientes no se distribuyen en tres grupos claramente separados, sino en un continuo de estados transcriptómicos con regiones parcialmente diferenciadas.

### Figura 7 — HDBSCAN en PCA

![HDBSCAN PCA](img/hdbscan_pca.png)

La visualización de HDBSCAN en PCA aparece prácticamente como una única nube sin separación cromática clara, lo que indica que el método basado en densidad no ha identificado clusters densos bien definidos bajo los parámetros utilizados.

Este resultado es importante porque HDBSCAN no fuerza un número fijo de clusters. A diferencia de K-means, que siempre divide las muestras en k grupos, HDBSCAN solo genera clusters si detecta regiones de densidad suficientemente separadas. El hecho de que no aparezcan agrupaciones claras apoya la interpretación de que la estructura de los datos no contiene grupos densos naturales fácilmente detectables.

Biológicamente, esto refuerza la hipótesis de que los patotipos no se comportan como entidades discretas en el espacio transcriptómico bulk, sino como estados parcialmente solapados.

### Figura 8 — Spectral Clustering en PCA

![Spectral Clustering PCA](img/spectral_pca.png)

El clustering espectral genera una partición distinta a K-means. La mayoría de las muestras quedan agrupadas en un gran cluster, mientras que otro conjunto se localiza hacia valores altos de PC1. También se observa una muestra aislada en color negro, probablemente correspondiente a un caso extremo o a una asignación minoritaria.

Spectral Clustering es útil cuando los clusters no tienen forma esférica, pero en este caso tampoco produce una solución claramente alineada con tres grupos equilibrados y bien separados. El método parece capturar principalmente una separación entre el núcleo central de muestras y un subconjunto desplazado hacia PC1 positivo, más que una correspondencia clara con los tres patotipos histológicos.

Esto sugiere que la estructura principal de los datos puede estar dominada por gradientes transcriptómicos o por muestras atípicas, no por tres clases biológicas claramente separadas.

### Figura 9 — PCA coloreado por patotipo antes de la limpieza final

![PCA por patotipo primera versión](img/pca_patotipo_v1.png)

Esta figura colorea directamente las muestras por patotipo de referencia. Se observa un solapamiento considerable entre los grupos, especialmente entre Lymphoid y Myeloid. Algunas muestras Fibroid se concentran en regiones más negativas de PC1, pero no forman un grupo completamente aislado.

El resultado muestra que las etiquetas histopatológicas no se traducen de forma directa en grupos claramente separados en las dos primeras componentes principales. Esto no invalida los patotipos, sino que indica que la señal histológica puede no estar representada de forma simple por las principales fuentes de variación transcriptómica del tejido completo.

La AR sinovial combina infiltración inmune, activación estromal, angiogénesis, inflamación local y composición celular variable. Por ello, muestras con distinto patotipo pueden compartir una parte importante de su expresión génica global.

### Figura 10 — PCA coloreado por patotipo tras unificar Fibrous y Fibroid

![PCA por patotipo limpio](img/pca_patotipo_v2.png)

Tras limpiar las etiquetas y mantener los tres patotipos principales, el PCA sigue mostrando un fuerte solapamiento. El grupo Fibroid aparece con mayor presencia en la región izquierda-superior, mientras que Lymphoid y Myeloid se distribuyen de forma más amplia y mezclada.

Esta figura confirma que los patotipos no son linealmente separables en el espacio de PC1 y PC2. Aunque existen tendencias globales, como cierta acumulación de Fibroid en una zona concreta, no hay una frontera clara que permita distinguir los tres grupos únicamente con las dos primeras componentes principales.

El resultado es coherente con la baja silueta de K-means: si los patotipos de referencia ya se solapan en PCA, es esperable que los clusters no supervisados también presenten una concordancia parcial y no perfecta.

### Figura 11 — Heatmap de contingencia completo: Clusters vs Patotipos

![Heatmap completo](img/heatmap_clusters_vs_patotipos_completo.png)

Este heatmap cruza los clusters K-means con todos los patotipos presentes en los metadatos originales, incluyendo **Fibrous** y **Ungraded**. Los valores observados son:

| Cluster | Fibroid | Fibrous | Lymphoid | Myeloid | Ungraded |
|---|---:|---:|---:|---:|---:|
| 1 | 4 | 0 | 78 | 27 | 0 |
| 2 | 0 | 0 | 21 | 10 | 0 |
| 3 | 29 | 1 | 14 | 22 | 4 |

El resultado muestra una asociación clara entre el **Cluster 1** y el patotipo **Lymphoid**, con 78 muestras Lymphoid dentro de este cluster. Sin embargo, también contiene 27 muestras Myeloid, lo que indica que no es un cluster exclusivamente linfocitario.

El **Cluster 3** está enriquecido en Fibroid/Fibrous, con 30 muestras si se combinan ambas etiquetas, pero también contiene muestras Lymphoid y Myeloid. Por tanto, este cluster captura parcialmente un eje fibroide, aunque no de forma pura.

El **Cluster 2** es el más pequeño y no muestra una asociación biológica clara, ya que contiene principalmente Lymphoid y Myeloid. Podría representar una zona intermedia del espacio transcriptómico.

### Figura 12 — Heatmap de contingencia limpio: Fibrous unido a Fibroid y Ungraded eliminado

![Heatmap limpio](img/heatmap_clusters_vs_patotipos_limpio.png)

Tras unificar **Fibrous** dentro de **Fibroid** y eliminar **Ungraded**, la tabla queda reducida a los tres patotipos principales:

| Cluster | Fibroid | Lymphoid | Myeloid |
|---|---:|---:|---:|
| 1 | 4 | 78 | 27 |
| 2 | 0 | 21 | 10 |
| 3 | 30 | 14 | 22 |

Esta versión es más apropiada para la interpretación biológica porque compara los clusters con las tres categorías principales de referencia. El patrón confirma tres ideas clave:

- **Cluster 1** está claramente enriquecido en **Lymphoid**, por lo que K-means captura parcialmente el eje linfocitario.
- **Cluster 3** está enriquecido en **Fibroid**, aunque mantiene una mezcla importante con Myeloid y Lymphoid.
- **Myeloid** no queda aislado en un único cluster, sino repartido entre los tres grupos.

La principal conclusión biológica es que el clustering recupera parcialmente los extremos Lymphoid y Fibroid, pero no consigue separar adecuadamente el patotipo Myeloid. Esto puede deberse a que Myeloid comparte señales inflamatorias tanto con Lymphoid como con Fibroid, y a que RNA-seq bulk mezcla la expresión de múltiples poblaciones celulares.

### Figura 13 — Diagrama de Sankey: Clusters hacia patotipos

![Sankey](img/sankey_clusters_patotipos.png)

El diagrama de Sankey representa de forma visual el flujo entre los clusters generados por K-means y los patotipos de referencia. La anchura de cada flujo indica el número de muestras que conectan un cluster con un patotipo.

El flujo más relevante es el que conecta el **Cluster 1** con **Lymphoid**, confirmando la fuerte asociación observada en el heatmap. También destaca el flujo desde el **Cluster 3** hacia **Fibroid**, aunque este cluster se reparte también hacia Myeloid y Lymphoid. El **Cluster 2** presenta flujos hacia Lymphoid y Myeloid, sin una identidad clara.

Este gráfico resume bien el resultado global del análisis: existe una concordancia parcial entre clusters y patotipos, pero no una correspondencia uno-a-uno. Los clusters generados por expresión génica no reproducen perfectamente las etiquetas histopatológicas, sino que capturan algunos ejes dominantes de variación molecular.

---

## 3. Comparación entre métodos de clustering

### 3.1 Tabla resumen

| Criterio | K-means | Jerárquico Ward.D2 | HDBSCAN | Spectral Clustering |
|---|---|---|---|---|
| Entrada | VST + z-score | VST + z-score | VST + z-score | VST + z-score |
| Número de clusters | k = 3 impuesto | k = 3 por corte | No impuesto | k = 3 impuesto |
| Separación en PCA | Parcial, con solapamiento | No mostrada en las figuras finales aportadas | No detecta grupos claros | Partición dominada por un gran grupo y un subconjunto desplazado |
| Evaluación interna | Silhouette medio = 0.19 | No evaluado formalmente | Exploratorio | No evaluado formalmente |
| Concordancia con patotipos | Parcial: Lymphoid y Fibroid mejor capturados que Myeloid | No cruzado formalmente | No útil para patotipos | No cruzado formalmente |
| Interpretabilidad biológica | Moderada | Limitada sin validación externa | Baja en estos parámetros | Limitada |
| Limitación principal | Fuerza k = 3 aunque la estructura sea continua | Sensible a distancia euclídea y corte elegido | No identifica densidades separadas | Sensible al kernel y a outliers |

### 3.2 Interpretación comparativa

El método más útil para este trabajo es **K-means**, no porque produzca una partición muy robusta, sino porque permite comparar directamente los clusters con los tres patotipos principales. Además, su solución tiene una interpretación parcial: el Cluster 1 se asocia fuertemente con Lymphoid y el Cluster 3 muestra enriquecimiento en Fibroid.

El método **HDBSCAN** aporta una conclusión importante: al no forzar k, no encuentra clusters densos claramente separados. Esto respalda la idea de que la estructura real de los datos es débil o continua.

El **Spectral Clustering** detecta una partición diferente, pero no muestra una correspondencia evidente con los tres patotipos. Parece separar principalmente un subconjunto desplazado en PC1, lo que puede reflejar un gradiente transcriptómico más que una clasificación histopatológica.

El **clustering jerárquico Ward.D2** está implementado en el script, pero en la información visual aportada no se dispone de una evaluación equivalente mediante silhouette o cruce con patotipos. Por tanto, no puede considerarse superior a K-means en este análisis concreto.

---

## 4. Interpretación de por qué los resultados son parcialmente limitados

### 4.1 Los patotipos de AR no son clases completamente discretas

La principal explicación biológica es que los patotipos sinoviales de artritis reumatoide no se comportan como categorías completamente separadas. Un paciente puede presentar infiltrado linfocitario, activación mieloide y componente fibroblástico en diferentes proporciones. Por ello, las muestras pueden situarse en posiciones intermedias entre patotipos.

Este comportamiento es incompatible con la idea de clusters perfectamente separados. K-means obliga a cada muestra a pertenecer a un único grupo, pero la biología subyacente puede ser más gradual.

### 4.2 RNA-seq bulk diluye la señal celular

Los datos proceden de tejido sinovial completo. Esto significa que cada muestra contiene una mezcla de fibroblastos, macrófagos, linfocitos B, linfocitos T, células endoteliales y otros tipos celulares. La expresión medida es un promedio de todas esas poblaciones.

Los patotipos histológicos se definen en gran medida por composición celular, pero en RNA-seq bulk esa composición se traduce en una señal mezclada. Por este motivo, un patotipo Myeloid puede compartir expresión inflamatoria con Lymphoid y también componentes estromales con Fibroid.

### 4.3 El patotipo Myeloid es el más difícil de aislar

La matriz de contingencia muestra que Myeloid está repartido entre los tres clusters. Esto indica que no existe un cluster puramente mieloide en la solución K-means.

Una posible explicación es que Myeloid representa un eje inflamatorio transversal, compartido parcialmente con muestras Lymphoid y Fibroid. También puede reflejar que las firmas mieloides en bulk RNA-seq no son suficientemente dominantes frente a otras fuentes de variación.

### 4.4 La estructura interna es débil según silhouette

El silhouette medio de 0.19 indica que la separación interna entre clusters es baja. Muchas muestras se encuentran cerca de fronteras entre grupos, especialmente en el Cluster 3, cuya silueta media es 0.07.

Esto significa que, aunque los clusters puedan tener cierto significado exploratorio, no deben utilizarse como clasificación definitiva de pacientes sin validación adicional.

### 4.5 La selección de los 1.000 genes más variables puede introducir ruido

Seleccionar genes por varianza es una estrategia habitual, pero no necesariamente selecciona los genes más relevantes para patotipos. Algunos genes muy variables pueden reflejar efectos técnicos, diferencias de calidad de muestra, profundidad de secuenciación residual o procesos biológicos no relacionados directamente con la clasificación sinovial.

Una selección guiada por firmas inmunológicas, fibroblásticas o mieloides podría mejorar la interpretabilidad.

### 4.6 El estadístico Gap sugiere k = 3, pero no una solución fuerte

Aunque el método Gap marca k = 3, la curva continúa aumentando para valores mayores. Esto indica que los datos podrían admitir particiones adicionales o subestructuras, pero no necesariamente clusters biológicos claros.

En enfermedades heterogéneas, un número fijo de grupos puede ser una simplificación excesiva. Una alternativa sería modelar los datos como gradientes o puntuaciones de firma, en lugar de como clases cerradas.

---

## 5. Propuestas de mejora

### 5.1 Utilizar firmas génicas específicas de patotipo

En lugar de seleccionar únicamente los 1.000 genes más variables, se recomienda combinar la selección por varianza con firmas biológicas conocidas:

- **Lymphoid**: genes asociados a linfocitos B/T, organización linfoide, quimiocinas y señalización inmune adaptativa.
- **Myeloid**: genes asociados a macrófagos, monocitos, inflamación innata, IL-1, TNF y metaloproteinasas.
- **Fibroid/Fibroblast**: genes asociados a matriz extracelular, colágeno, fibroblastos sinoviales y remodelación tisular.

Esto reduciría ruido y aumentaría la relación señal/ruido del análisis.

### 5.2 Aplicar puntuaciones de firma en lugar de clustering directo

Una alternativa más interpretable sería calcular scores por muestra para firmas Lymphoid, Myeloid y Fibroid, utilizando métodos como ssGSEA, GSVA o medias estandarizadas de genes marcador.

Después, las muestras podrían representarse en un espacio de tres puntuaciones biológicas, más directamente relacionado con los patotipos que el espacio de 1.000 genes.

### 5.3 Evaluar estabilidad mediante clustering de consenso

Se recomienda aplicar métodos como `ConsensusClusterPlus`, repitiendo el clustering sobre subconjuntos de muestras y genes. Esto permitiría comprobar si los grupos son estables o si dependen excesivamente de pequeñas variaciones en los datos.

Dado el bajo silhouette, esta evaluación es especialmente importante antes de extraer conclusiones clínicas.

### 5.4 Usar modelos probabilísticos o soft clustering

Los métodos de asignación rígida obligan a cada muestra a pertenecer a un único cluster. Para una enfermedad continua como la AR, puede ser más adecuado utilizar modelos de mezcla gaussianos o enfoques probabilísticos, en los que cada muestra tenga una probabilidad de pertenencia a cada grupo.

Esto permitiría identificar muestras intermedias, por ejemplo pacientes con perfil parcialmente Lymphoid y parcialmente Myeloid.

### 5.5 Incorporar deconvolución celular

Dado que los patotipos están muy relacionados con la composición celular del tejido, sería recomendable aplicar herramientas de deconvolución celular para estimar proporciones de tipos celulares a partir del RNA-seq bulk.

Posteriormente, el clustering podría realizarse sobre proporciones celulares estimadas o combinarlas con expresión génica. Este enfoque estaría más alineado con la definición histopatológica de los patotipos.

### 5.6 Validar los resultados con variables clínicas

El script ya selecciona variables clínicas como tender joint count, swollen joint count, actividad de la artritis, ESR, CRP y evaluación global del médico. Una mejora natural sería comparar los clusters con estas variables para comprobar si los grupos tienen relevancia clínica.

Por ejemplo, se podría evaluar si el cluster enriquecido en Lymphoid presenta mayor inflamación sistémica, mayor CRP o mayor actividad clínica que los clusters enriquecidos en Fibroid.

---

## 6. Conclusión general

El análisis no supervisado de los datos STRAP RNA-seq mediante VST, escalado, PCA, UMAP y varios métodos de clustering muestra que la estructura transcriptómica de las muestras sinoviales de artritis reumatoide es **heterogénea, continua y parcialmente solapada**.

La solución K-means con k = 3 es útil como aproximación exploratoria porque recupera parcialmente algunos ejes biológicos: el **Cluster 1** está claramente enriquecido en muestras **Lymphoid**, mientras que el **Cluster 3** muestra enriquecimiento en **Fibroid**. Sin embargo, el patotipo **Myeloid** no queda separado en un único grupo, sino distribuido entre los tres clusters.

La evaluación interna confirma esta limitación: el silhouette medio global es **0.19**, un valor bajo que indica baja separación entre clusters. Además, HDBSCAN no identifica agrupaciones densas claras, lo que refuerza la idea de que los datos no contienen clusters discretos fuertes.

En conjunto, los resultados sugieren que los patotipos de AR no deben interpretarse como clases moleculares perfectamente separadas en RNA-seq bulk, sino como estados biológicos parcialmente superpuestos. El análisis es válido como exploración inicial, pero para mejorar su robustez sería recomendable incorporar firmas génicas, deconvolución celular, clustering de consenso, modelos probabilísticos y validación clínica.

---

## 7. Nota sobre coherencia entre código, gráficos e interpretación

- El script utiliza **VST de DESeq2**, lo que mejora la calidad del preprocesamiento respecto a una transformación log₂ simple.
- El análisis se realiza sobre los **1.000 genes más variables**, seleccionados antes de la normalización VST.
- K-means, HDBSCAN, Spectral Clustering, PCA y UMAP se aplican sobre `expr_scaled`, es decir, datos VST escalados.
- La comparación con patotipos se realiza después de alinear las muestras de expresión con los metadatos de `E-MTAB-13733.sdrf.txt`.
- En el análisis final se unifica **Fibrous** dentro de **Fibroid** y se elimina **Ungraded**, lo cual es adecuado para comparar con los tres patotipos principales.
- La silueta mostrada corresponde a la solución K-means con k = 3, con media global **0.19**.
- La matriz de contingencia limpia indica una concordancia parcial: Lymphoid se asocia principalmente al Cluster 1, Fibroid al Cluster 3, y Myeloid queda repartido entre los tres clusters.
