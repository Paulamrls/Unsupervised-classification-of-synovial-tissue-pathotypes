# Unsupervised Classification of Synovial Tissue Pathotypes in Rheumatoid Arthritis
### Transcriptomic Analysis of the STRAP Cohort (ArrayExpress E-MTAB-13733)

---

## Overview

Rheumatoid arthritis (RA) is not a single disease — it is a spectrum of synovial pathologies with distinct cellular compositions and molecular programs. Histopathological classification defines three canonical pathotypes: **Fibroid**, **Myeloid**, and **Lymphoid**. However, whether these categories are truly discrete at the transcriptomic level, or whether they blend into a continuous landscape, remains an open question.

This project applies unsupervised clustering to bulk RNA-seq data from the STRAP cohort (n ≈ 300 synovial biopsies) to ask: *can gene expression alone reconstruct biologically meaningful subgroups that align with — or go beyond — the canonical pathotype classification?*

The pipeline benchmarks **8 clustering methods** (K-means, Hierarchical, Spectral, Consensus, Leiden, Louvain, Infomap, MCL), selects the most robust solution, and characterises each cluster through differential expression, known RA gene signatures, and comparison with histological ground truth.

---

## Key Findings

> **K-means with k = 4 identifies four transcriptomic subgroups.** Three align broadly with the canonical Fibroid, Myeloid and Lymphoid pathotypes. A fourth, designated *Subtype4*, co-segregates with histological Lymphoid samples but carries a distinct molecular identity, suggesting a previously unresolved subdivision of the lymphoid-rich synovial state.

| Pathotype | Colour | n (approx.) | Dominant biology |
|-----------|--------|------------|-----------------|
| **Fibroid** | 🔴 | ~55 | Stromal activation, antigen presentation (HLA genes), ECM remodelling |
| **Myeloid** | 🔵 | ~10 | Innate inflammatory programme; small, transcriptomically concentrated cluster |
| **Lymphoid** | 🟢 | ~40 | Ribosomal/translational activity (RPL28, RPS19); T/B cell infiltration |
| **Subtype4** | 🟣 | ~110 | CD163⁺ macrophage-like signature, immune cell trafficking (PLXNB2), stress response |

---

## Dimensionality Reduction

Two complementary approaches were used to visualise the transcriptomic landscape:

- **PCA** — captures linear variance structure; PC1 (22.2%) and PC2 (21.8%) together explain ~44% of total variance.
- **Diffusion Maps** — captures non-linear geometry and biological trajectories; axes clipped to the 2.5–97.5 percentile range to prevent outlier collapse.

<table>
<tr>
<th align="center">PCA — K-means solution</th>
<th align="center">Diffusion Map — K-means solution</th>
</tr>
<tr>
<td><img src="results/figures/pca_kmeans.png" width="480"/></td>
<td><img src="results/figures/diffusion_map.png" width="480"/></td>
</tr>
</table>

**PCA reading:** Fibroid (red) clusters in the upper-left quadrant; Lymphoid (green) extends toward positive PC1; Subtype4 (purple) forms a dense mass in the lower portion; Myeloid (blue) is a small, spatially distinct group. The separation is real but not sharp — consistent with the known biological continuum of RA pathotypes.

**Diffusion Map reading:** The non-linear geometry reveals a striking **diagonal trajectory for Subtype4**, suggesting it represents a continuous state rather than a discrete subtype. Fibroid samples spread along a separate arm. Myeloid appears almost isolated, indicating strong transcriptomic distinctiveness despite its small size.

---

## The Four Pathotypes — Cluster-Level Interpretation

### 🔴 Fibroid

The Fibroid cluster is the second-largest group and shows the most complex composition in the confusion matrix: ~48% of samples derive from histological Fibroid, but ~33% come from histological Myeloid. This mixture indicates that at the transcriptomic level, stromal and innate inflammatory programs co-exist in a substantial fraction of samples.

**Top DEG markers (upregulated):** *PRKAR2B, MTURN, GSN, MAOA, CD74, HLA-DRA, HLA-DMB, HLA-DPB1, HLA-DPA1*

The prevalence of **HLA class II genes** (HLA-DRA, HLA-DMB, HLA-DPB1, HLA-DPA1) is biologically informative: these are the major histocompatibility complex molecules responsible for antigen presentation to CD4⁺ T cells. Their upregulation in the Fibroid cluster may reflect that synovial fibroblasts (FLS) in this state are actively participating in local antigen presentation, blurring the boundary between stromal and immune cell functions. *GSN* (gelsolin, an actin-remodelling protein) and *MAOA* (monoamine oxidase A) point to cytoskeletal dynamics and neuroinflammatory components.

<img src="results/figures/boxplots_markers_Fibroid.png" width="900"/>

---

### 🔵 Myeloid

The Myeloid cluster is the smallest by sample count (~10 samples), which is the primary reason only one statistically significant DEG was recovered after DESeq2 one-vs-rest testing. Despite its size, **60% of its samples come from histological Myeloid** tissue — the highest pathotype purity of all four clusters.

**Top DEG marker:** *ZFP36L1* — an RNA-binding protein that destabilises pro-inflammatory mRNAs (TNF, IL-6). Its downregulation in Myeloid samples suggests impaired post-transcriptional dampening of inflammatory signals, consistent with the hyperactivated macrophage state characteristic of the Myeloid pathotype.

**Notable feature — TNNC1 and CRYAB:** In the Subtype4 marker boxplot, *TNNC1* (troponin C1) and *CRYAB* (αB-crystallin) appear dramatically elevated specifically in Myeloid samples. Both are stress-response proteins expressed in fibroblasts under mechanical or thermal stress, and CRYAB has documented roles in macrophage survival. Their high expression here may reflect either synovial fibroblast co-activation within Myeloid biopsies or a stress-adapted subpopulation.

<img src="results/figures/boxplots_markers_Myeloid.png" width="480"/>

---

### 🟢 Lymphoid

The Lymphoid cluster achieves the strongest single-pathotype overlap: **70% of its samples come from histological Lymphoid** tissue. Its DEG profile centres on ribosomal proteins (*RPL28, RPS19*) and extracellular matrix components (*LAMB2*).

The upregulation of **ribosomal proteins** in the Lymphoid cluster is not trivial. Lymphocytes that have been activated and expanded — as occurs in germinal centre-like structures of RA synovium — undergo dramatic upregulation of the translational machinery to sustain rapid proliferation and antibody/cytokine production. This transcriptional signature aligns with the known biology of lymphoid aggregates (LAs) in RA synovium.

*LAMB2* (laminin subunit beta-2) is a structural component of basement membranes. Its differential expression between Lymphoid and other clusters may reflect the vascular remodelling required to sustain lymphoid aggregate organisation.

<img src="results/figures/boxplots_markers_Lymphoid.png" width="700"/>

---

### 🟣 Subtype4 — A Candidate Lymphoid Subdivision

Subtype4 is the largest cluster and the most biologically unexpected result of this analysis. Its confusion matrix profile mirrors Lymphoid: **70.9% of its samples carry histological Lymphoid labels**, almost identical to the Lymphoid cluster itself. Yet its gene expression profile is fundamentally different.

**Top DEG markers (upregulated):** *CD163, PLXNB2, HSP90B1, ADIRF, TPP1*  
**Top DEG markers (downregulated):** *TNNC1, CRYAB*

*CD163* is a scavenger receptor and a definitive marker of alternatively activated (M2-like) macrophages. Its upregulation positions Subtype4 as enriched for **immunosuppressive or anti-inflammatory macrophage activity**, in contrast to the pro-inflammatory bias expected in canonical Myeloid samples.

*PLXNB2* (Plexin B2) mediates immune cell chemotaxis and is expressed on dendritic cells and macrophages during tissue surveillance. *HSP90B1* (GRP94) is an endoplasmic reticulum chaperone critical for secretory pathway proteins, including MHC class I/II complexes — hinting at active antigen processing machinery.

**Interpretation:** Subtype4 may represent a **lymphoid-rich synovial state dominated by regulatory or alternatively activated macrophages** rather than the classically cytotoxic or effector immune programme of Lymphoid. If validated, this distinction would have therapeutic implications: conventional anti-inflammatory strategies may behave differently in these two histologically indistinguishable but transcriptomically separable subgroups.

<img src="results/figures/boxplots_markers_Subtype4.png" width="800"/>

---

## Validation Against Histological Ground Truth

### Confusion Matrix

The confusion matrix quantifies, for each predicted cluster, what fraction of its samples comes from each histological pathotype.

<img src="results/figures/confusion_pct.png" width="600"/>

Key observations:
- **Lymphoid** cluster: 70% Lymphoid histological → highest purity, strongest biological recovery
- **Myeloid** cluster: 60% Myeloid histological → high purity, but very small cluster
- **Subtype4**: 70.9% Lymphoid histological → transcriptomically distinct subgroup within the Lymphoid category
- **Fibroid**: 48.3% Fibroid histological → most heterogeneous cluster; mixed stromal-immune state

### Sankey Diagram — Sample Flow Between Prediction and Ground Truth

<img src="results/figures/sankey.png" width="800"/>

The Sankey diagram makes the cross-cluster flow immediately readable. The most notable feature is the **double flow from histological Lymphoid** into both the green (Lymphoid) and purple (Subtype4) clusters — direct visual evidence of the proposed Lymphoid subdivision. The Myeloid band on the left is thin, reflecting its small size, but its flow is concentrated toward histological Myeloid on the right.

---

## Heatmap — Global Expression Structure

<img src="results/figures/heatmap_clusters.png" width="900"/>

The heatmap of the top variable genes per pathotype confirms a block structure: columns (samples) cluster by pathotype with visible colour bar separations. The most prominent pattern is a **warm stripe in Myeloid samples** (right side of the heatmap), corresponding to a small set of genes highly specific to this cluster. The Fibroid and Subtype4 blocks are large but relatively uniform, consistent with their dominant presence in the dataset.

---

## Clustering Method Benchmarking

Eight algorithms were evaluated on the same VST-normalised, z-scored expression matrix. Two independent metrics were used:

### Internal validity: Silhouette score

<img src="results/figures/silhouette_comparison.png" width="750"/>

### External stability: Bootstrap ARI (Adjusted Rand Index)

<img src="results/figures/bootstrap_stability.png" width="750"/>

| Method | Silhouette | Bootstrap ARI | Verdict |
|--------|-----------|--------------|---------|
| **K-means** | **0.19** | **0.829** | ✅ Best overall — high stability, highest silhouette |
| Hierarchical | 0.168 | 0.655 | ✅ Strong second |
| Spectral | 0.138 | 0.528 | Moderate |
| Leiden | 0.094 | 0.633 | Graph method, decent stability |
| Infomap | 0.084 | 0.491 | Borderline |
| Consensus | 0.073 | — | Lower than expected; small clusters affect convergence |
| Louvain | 0.054 | 0.528 | Produces extra clusters beyond k=4 |
| MCL | 0.015 | 0.582 | Near-random separation for this data |

**Why K-means wins:** The transcriptomic data, preprocessed through PCA space after VST normalisation, has an approximately **Euclidean geometry** that favours centroid-based methods. Graph-based methods (Leiden, Louvain, MCL) are designed for manifold-structured data like single-cell RNA-seq, and underperform on bulk RNA-seq matrices where local neighbourhood graphs are less informative.

The **bootstrap ARI of 0.829** for K-means means that when the same algorithm is re-run on 80% subsamples of the data, the resulting partitions agree with the full solution 83% of the time (by chance-corrected concordance). This is a strong result and justifies using K-means as the reference solution.

---

## A Note on Biological Continuity

A silhouette score of 0.19, while the best among all methods, is still low in absolute terms. This is expected and biologically meaningful: **RA synovial pathotypes are not discrete molecular entities**. Each biopsy contains a heterogeneous mixture of cell types, and bulk RNA-seq measures a population average. The clustering reveals dominant transcriptomic axes, but patients exist along a continuum rather than in sharply bounded categories.

The discovery of Subtype4 as a possible Lymphoid subdivision illustrates this: histology collapses it into the Lymphoid category, but the transcriptome separates it. This finding motivates future work with single-cell resolution or spatial transcriptomics, where cell-type-level resolution could resolve the macrophage-lymphocyte interactions hypothesised in Subtype4.

---

## Pipeline Architecture

The analysis is implemented as a modular R pipeline. A single global parameter (`k <- 4` in `STRAP_pipeline_main.R`) controls the number of pathotypes across all steps.

```
STRAP_pipeline_main.R          ← Entry point; all global parameters here
functions/
  00_setup.R                   ← Package installation and output directory creation
  01_preprocessing.R           ← Count loading, VST normalisation, z-score scaling
  02_dimreduction.R            ← PCA (prcomp) + Diffusion Maps (diffusionMap)
  03_clustering.R              ← 8 clustering algorithms (modular, easily extensible)
  04_validation.R              ← Silhouette + bootstrap ARI
  05_annotation.R              ← RA gene signatures → automatic pathotype naming
  06_visualization.R           ← All figures (ggplot2, pheatmap, ggalluvial)
  07_deg_analysis.R            ← DESeq2 one-vs-rest; separate UP/DOWN tables per cluster
  08_export.R                  ← Organised export to results/
results/
  figures/                     ← All 21 publication-ready figures (PNG, 300 dpi)
  tables/                      ← Cluster assignments, validation ranking
  DEGs/                        ← Up- and down-regulated genes per pathotype
  clusters/                    ← Raw cluster vectors + consensus plots
  validation/                  ← Bootstrap and silhouette results
```

To reproduce the full analysis, open `STRAP_pipeline_main.R` in RStudio and press **Source**. The working directory is set automatically.

---

## Data

- **Cohort:** STRAP (*Stratification of Biologic Therapies for RA by Pathobiology*)
- **Accession:** ArrayExpress [E-MTAB-13733](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13733)
- **Data type:** Bulk RNA-seq, synovial tissue biopsies
- **Histological reference:** Fibroid, Myeloid, Lymphoid pathotype labels from matched SDRF metadata

---

## Dependencies

R ≥ 4.4. Key packages: `DESeq2`, `ConsensusClusterPlus`, `clusterProfiler`, `igraph`, `kernlab`, `diffusionMap`, `MCL`, `ggplot2`, `ggalluvial`, `pheatmap`, `mclust`. Full list managed automatically by `functions/00_setup.R`.

---

*Master's Thesis in Bioinformatics — Analysis of synovial tissue pathotypes in rheumatoid arthritis using unsupervised transcriptomic classification.*
