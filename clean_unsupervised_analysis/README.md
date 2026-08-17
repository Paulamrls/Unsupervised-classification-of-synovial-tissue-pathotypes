# Análisis transcriptómico no supervisado limpio - cohorte STRAP

## Objetivo de esta versión

Esta carpeta documenta una revisión **provisional** del pipeline realizada para avanzar antes de la reunión con el tutor. No pretende imponer una solución final, sino corregir posibles fuentes de supervisión y comparar de forma justa `k = 4` y `k = 5`.

## Mensaje principal para la reunión

- El análisis principal es ahora realmente no supervisado: la histología no participa en la selección de genes ni en el clustering.
- Las cinco biopsias con contaminación muscular se excluyen antes de normalizar o analizar los datos.
- `k = 4` y `k = 5` se ejecutan con exactamente las mismas 205 muestras, 3.000 genes, 30 PCs, semilla y parámetros.
- `k = 4` es ligeramente más parsimonioso, cohesivo y estable que `k = 5`.
- `k = 5` subdivide principalmente el cluster 1 de `k = 4` y se conserva como hipótesis exploratoria.
- Esa subdivisión presenta señal biológica, pero no evidencia suficiente para denominarla formalmente `IFN-high`.

## Cambios realizados y justificación

1. **Exclusión temprana de cinco muestras contaminadas.** Se evita que la señal muscular afecte al VST, la varianza, PCA, clustering, validación, DEGs y figuras.
2. **Kruskal-Wallis desactivado en el análisis principal.** Seleccionar genes usando el patotipo histológico introducía supervisión y circularidad.
3. **Método anterior conservado como opción.** `top_variable_plus_kw` permite reproducirlo, pero la configuración principal usa `top_variable`.
4. **Filtro mínimo de expresión explícito.** Se conserva un gen si tiene al menos 10 cuentas en 21 muestras, equivalentes al 10 % de las 205 muestras.
5. **VST anterior a la selección por varianza.** La varianza se calcula sobre valores estabilizados mediante DESeq2.
6. **Selección de los 3.000 genes VST más variables.** Se amplió desde 1.000 para no perder señales biológicas menos dominantes sin usar etiquetas.
7. **Clustering sobre 30 PCs.** Reduce ruido y mantiene una representación común para ambas soluciones.
8. **Comparación controlada de `k = 4` y `k = 5`.** La única diferencia de la solución final es el valor de `k`.
9. **k-means fijado como método final.** Evita que cada valor de `k` termine comparándose con un algoritmo distinto.
10. **Semilla global 42.** Hace reproducibles normalización, reducción dimensional, clustering y validación.
11. **Histología reservada para interpretación posterior.** Solo se usa después del clustering en matrices de confusión, Sankey, composición y anotación.
12. **Comprobaciones automáticas.** Verifican contaminación, duplicados, identidad de muestras/genes/PCs y suma de tamaños de clusters.
13. **Salidas separadas.** Los resultados nuevos no sobrescriben los análisis anteriores.
14. **Validación funcional de ID3 e ID4 de `k = 5`.** Incluye DESeq2, GO/GSEA, puntuación IFN, variables clínicas, silhouette y bootstrap Jaccard.

## Configuración exacta

| Parámetro | Valor |
|---|---:|
| Selección génica | `top_variable` |
| Conteo mínimo | 10 |
| Proporción mínima | 0,10 |
| Muestras mínimas | 21 de 205 |
| Genes seleccionados | 3.000 |
| PCs para clustering | 30 |
| Método final | k-means |
| Semilla | 42 |
| Valores de k | 4 y 5 |

## Muestras excluidas

- `STRAPPAT00113-baseline`
- `STRAPPAT00210-baseline`
- `STRAPPAT00199-baseline`
- `STRAPPAT00007-baseline`
- `STRAPPAT00051-baseline`

## Resultados comparativos

| Solución | Silhouette k-means | Bootstrap ARI medio | DE bootstrap | Tamaños |
|---|---:|---:|---:|---|
| `k = 4` | 0,1581 | 0,8676 | 0,0828 | 63, 65, 15, 62 |
| `k = 5` | 0,1522 | 0,8400 | 0,0773 | 60, 62, 25, 43, 15 |

Interpretación: ambas soluciones son posibles, pero `k = 4` obtiene una cohesión y estabilidad ligeramente superiores y es más parsimoniosa. Las silhouettes bajas indican que los estados transcriptómicos no forman fronteras completamente discretas.

## Cómo se transforma k = 4 en k = 5

La tabla de contingencia demuestra que tres grupos permanecen casi intactos. La principal diferencia es:

- `k4-ID1` (63 muestras) se divide principalmente en `k5-ID3` (22) y `k5-ID4` (41).
- `k4-ID2` conserva 60 de 65 muestras en `k5-ID1`.
- `k4-ID3` conserva sus 15 muestras en `k5-ID5`.
- `k4-ID4` conserva 61 de 62 muestras en `k5-ID2`.

Los IDs son numéricos y no implican un nombre biológico previo.

## Validación de la subdivisión ID3/ID4 de k = 5

- Se compararon los clusters completos `ID3` (`n = 25`) e `ID4` (`n = 43`).
- Se detectaron 2.747 DEGs con `|log2FC| >= 1` y FDR < 0,05: 1.917 superiores en ID4 y 830 en ID3.
- ID4 muestra enriquecimiento inmunitario, inflamatorio, quimiotáctico, angiogénico y de defensa antiviral.
- ID3 muestra señales relacionadas con cartílago/tejido conectivo, cilios y homeostasis de metales.
- La PCR fue superior en ID4 y permaneció significativa tras corrección múltiple (`FDR = 0,032`).
- La puntuación IFN fue mayor como tendencia en ID4, pero no significativa (`p = 0,124`).
- La vía específica de interferón tipo I tampoco alcanzó FDR < 0,05 (`FDR = 0,075`).
- Por tanto, ID4 puede describirse como **más inflamatorio/inmunitario**, pero no debe llamarse todavía `IFN-high`.
- ID3: silhouette media 0,192 y Jaccard bootstrap medio 0,638.
- ID4: silhouette media 0,095 y Jaccard bootstrap medio 0,736.
- La división posee señal biológica, aunque su separación y estabilidad son moderadas.

## Conclusión provisional

Se propone `k = 4` como solución principal por estabilidad, cohesión y parsimonia. `k = 5` se presenta como análisis secundario que identifica una posible subdivisión inflamatoria que requeriría validación externa. Este resultado no es negativo: indica que la heterogeneidad sinovial puede organizarse como estados parcialmente continuos y no necesariamente como cinco patotipos discretos.

## Cambios específicos del código

### `pipeline/STRAP_pipeline_main.R`

- Centraliza los parámetros comunes y ejecuta `K_VALUES <- c(4L, 5L)`.
- Define `GENE_SELECTION_METHOD <- "top_variable"`, 3.000 genes, filtro 10/10 %, 30 PCs, k-means y semilla 42.
- Define y elimina explícitamente las cinco muestras contaminadas antes del preprocesamiento.
- Ejecuta una única preparación común y reutiliza la misma matriz en ambos valores de `k`.
- Añade comprobaciones con parada inmediata si muestras, genes, PCs o tamaños no coinciden.
- Exporta parámetros, genes, muestras, asignaciones y contingencia común.

### `pipeline/functions/01_preprocessing.R`

- Añade el filtro mínimo de expresión.
- Aplica VST antes de ordenar los genes por varianza.
- Implementa `top_variable` sin consultar metadatos histológicos.
- Conserva `top_variable_plus_kw` como rama opcional y desactivada.
- Evita filtrar muestras según su patotipo durante la carga de metadatos.

### `pipeline/functions/02_dimreduction.R`

- Hace explícita la semilla en PCA y diffusion maps.
- Devuelve coordenadas reutilizables para una entrada común de clustering.

### `pipeline/functions/03_clustering.R`

- Propaga la misma semilla y los mismos parámetros a los algoritmos.
- Mantiene k-means como solución final común para la comparación.

### `pipeline/functions/04_validation.R`

- Calcula silhouette y bootstrap ARI de forma reproducible.
- Permite validar cada `k` sobre la misma matriz de 30 PCs.

### `pipeline/functions/05_deg_analysis.R`

- Mantiene la expresión diferencial posterior al clustering.
- Alinea conteos crudos y asignaciones antes de ejecutar DESeq2.

### `pipeline/functions/06_annotation.R`

- Separa la creación de clusters de su interpretación biológica posterior.
- Usa histología y firmas únicamente después de obtener las asignaciones.

### `pipeline/functions/07_visualization.R`

- Genera PCA, diffusion map, heatmap, Sankey, confusión, silhouette, bootstrap y marcadores por cada valor de `k`.

### `pipeline/functions/08_export.R`

- Organiza las salidas por ejecución y exporta parámetros, tablas, figuras y asignaciones sin sobrescribir análisis anteriores.

### `additional_validation/validate_k5_subgroups_ID3_ID4.R`

- Reconstruye la solución y verifica que coincide con la original mediante ARI = 1.
- Compara ID3 e ID4 mediante DESeq2, GO, GSEA, puntuación IFN, clínica, silhouette y bootstrap Jaccard.

## Figuras recomendadas para enseñar al tutor

1. `figures/k4/pca_kmeans.png`: estructura PCA de la solución principal.
2. `figures/k5/pca_kmeans.png`: estructura PCA de la solución alternativa.
3. `figures/k4/silhouette_comparison.png`: cohesión interna en `k = 4`.
4. `figures/k5/silhouette_comparison.png`: cohesión interna en `k = 5`.
5. `figures/k4/bootstrap_stability.png`: estabilidad de `k = 4`.
6. `figures/k5/bootstrap_stability.png`: estabilidad de `k = 5`.
7. `figures/k4/sankey.png`: relación posterior con los patotipos histológicos.
8. `figures/k5/sankey.png`: relación posterior de los cinco clusters.
9. `figures/functional_validation_ID3_ID4/IFN_ID3_vs_ID4.png`: puntuación IFN de la subdivisión.
10. `figures/functional_validation_ID3_ID4/silhouette_ID3_vs_ID4.png`: separación de ID3 e ID4.

También se incluyen todas las demás figuras de ambos análisis para consulta.

## Estructura de esta carpeta

```text
clean_unsupervised_analysis/
├── README.md
├── pipeline/                         # pipeline final reproducible
├── additional_validation/            # validación funcional ID3/ID4
├── figures/
│   ├── k4/
│   ├── k5/
│   └── functional_validation_ID3_ID4/
└── tables/
    ├── common/
    ├── k4/
    ├── k5/
    └── comparison_k4_k5/
```

## Estado

Versión de trabajo para revisión con el tutor. Las decisiones biológicas y la elección definitiva de `k` deben confirmarse antes de congelar la memoria del TFM.
