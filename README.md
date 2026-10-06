# Single-cell RNA-seq analysis of conditional intestinal adenoma

This repository contains the R scripts and supporting files used for the single-cell RNA sequencing (scRNA-seq) analysis of conditional intestinal adenomas described in:

> **Surakhy M, Matheson J, Barnes DJ, et al.**
> *Smad4 and TGFβ1 dependent gene expression signatures in conditional intestinal adenoma, organoids and colorectal cancer.*
> **Scientific Reports.** 2025;15:16330.
> https://doi.org/10.1038/s41598-025-00908-4

The study investigated the consequences of **Smad4 loss in an Apc-dependent murine intestinal adenoma model**, with particular emphasis on Smad4- and TGF-β1-dependent gene expression and changes in epithelial progenitor populations.

## Study overview

Conditional activation of the Lgr5-CreER^T2^ system in **Apc^fl/fl^Smad4^fl/fl^** mice generates intestinal adenomas with conditional Smad4 loss. The study compared adenomas with intact Smad4 signalling (**Smad4^+/+**) and adenomas with Smad4 loss (**Smad4^Δ/Δ^**).

Single-cell RNA sequencing of primary caecal adenomas was used to:

* characterise the cellular composition of the adenomas;
* identify transcriptionally distinct epithelial and non-epithelial cell populations;
* identify marker genes defining individual cell populations;
* compare cell-type proportions between Smad4 genotypes;
* identify genotype-associated differential gene expression;
* investigate progenitor populations associated with Smad4 loss.

The published study identified expansion of **Lgr5^low^, Pak3^high^ and Id1^low^ progenitor populations** in Smad4-deficient caecal adenomas.

## scRNA-seq data

The raw scRNA-seq data generated for the study are publicly available through the European Nucleotide Archive:

**ENA accession: PRJEB85413**

The bulk RNA-seq data generated in the study are available under:

**ENA accession: PRJEB85309**

The complete experimental and biological context is described in the associated publication.

## Analysis workflow

The repository implements the following analysis workflow:

```text
scRNA-seq data
      │
      ↓
01. Seurat object creation and quality control
      │
      ↓
02. SCTransform and dataset integration
      │
      ↓
03. Clustering and marker-gene identification
      │
      ↓
04. Automated cell-type annotation using clustifyr
      │
      ↓
05. Manual biological cell-type annotation
    and cell-type proportion analysis
      │
      ├───────────────┐
      ↓               ↓
Cell proportions     Differential expression
scProportionTest     MAST
                      │
                      ↓
             07. Pseudo-bulk analysis
                 and heatmap generation
```

## Repository structure

```text
scRNAseq/
│
├── README.md
├── LICENSE
├── Sotware versions
│
├── scripts/
│   ├── 01_CreatingSeuratObj_and_QC.r
│   ├── 02_SCtransform_And_integration.R
│   ├── 03-findAllMarkers.r
│   ├── 04-DataLabelling_usingClustifyr.r
│   ├── 05_AssignCellTypeAndScproportiopn.r
│   ├── 06_DGE.r
│   └── 07-pseudoBulk_DGE_and_heatmap.r
│
├── figures/
│
└── shinyApp_SCRNA/
    ├── README.txt
    ├── app.R
    └── RDS files
```

## Analysis scripts

### 01 — Creating Seurat object and quality control

`01_CreatingSeuratObj_and_QC.r`

Creates the Seurat object from the processed scRNA-seq data and performs initial quality-control procedures prior to downstream analysis.

### 02 — SCTransform and integration

`02_SCtransform_And_integration.R`

Performs SCTransform-based normalisation and integration of the scRNA-seq datasets. Dimensionality reduction and neighbourhood-based analyses are performed to prepare the data for clustering.

### 03 — Cluster identification and marker genes

`03-findAllMarkers.r`

Identifies transcriptionally distinct cell clusters and determines cluster-associated marker genes using Seurat's `FindAllMarkers`.

Multiple clustering resolutions are evaluated. Resolution **1.2** is used to define the principal cell populations, while resolution **2.2** provides additional subdivision of epithelial populations, including enterocyte precursor populations.

### 04 — Automated cell-type annotation

`04-DataLabelling_usingClustifyr.r`

Uses **clustifyr** and reference expression profiles to provide automated cell-type annotations.

These annotations are used as a reference to support biological interpretation of the clusters.

### 05 — Cell-type annotation and cell proportions

`05_AssignCellTypeAndScproportiopn.r`

Assigns final biological cell-type identities based on cluster marker expression and published intestinal cell-type markers.

The script also evaluates differences in cell-type proportions between experimental groups using **scProportionTest**.

The resulting cell populations include epithelial progenitor populations as well as immune and stromal populations.

Representative marker-expression and cell-composition figures are generated in this step.

### 06 — Differential gene expression

`06_DGE.r`

Performs cell-type-specific differential expression analysis between **WT and Floxed** conditions.

Differential expression is performed using **MAST** within Seurat. Analyses are performed separately for the annotated cell populations to identify genotype-associated transcriptional changes within individual cell types.

### 07 — Pseudo-bulk differential expression and heatmap

`07-pseudoBulk_DGE_and_heatmap.r`

Performs the condition-level expression analysis described in the associated publication using Seurat's MAST implementation of `FindMarkers`.

The analysis compares the WT and Floxed conditions at resolution 0 and generates heatmaps of the most strongly differentially expressed genes.

The script is retained under the original project terminology of **pseudo-bulk analysis**, consistent with the associated publication.

## Cell-type annotation

Cell identities were determined using a combination of:

1. transcriptional cluster markers;
2. published intestinal scRNA-seq reference datasets;
3. automated reference-based annotation using **clustifyr**; and
4. manual biological interpretation of marker-gene expression.

Automated annotation was therefore used as a supporting annotation method rather than as the sole determinant of the final biological cell-type labels.

## Differential expression

Differential expression analyses were performed using the **MAST** statistical framework implemented in Seurat.

For the cell-type-specific analysis, WT and Floxed cells were compared within each annotated cell population.

The resolution-0 analysis provides a condition-level comparison across the full cellular population and was used to generate the corresponding differential-expression heatmap.

## Biological findings supported by the scRNA-seq analysis

The associated study identified several cellular and transcriptional features associated with Smad4 loss.

In particular, Smad4-deficient caecal adenomas showed expansion of:

* **Lgr5^low^ progenitor cells**
* **Pak3^high^ progenitor cells**
* **Id1^low^ progenitor cells**

The single-cell analysis therefore provided cellular context for the Smad4-dependent transcriptional signatures identified in the complementary bulk RNA-seq and organoid experiments.

The broader study identified **Id1^low^** and **Spp1^high^** as prominent Smad4/TGF-β1-dependent expression features and subsequently evaluated these signatures in human colorectal cancer.

## Figures

The `figures/` directory contains representative figures generated during the scRNA-seq analysis, including:

* cluster and cell-type UMAP visualisations;
* cell-type marker expression;
* cell-type proportion analysis;
* feature plots;
* differential-expression heatmaps; and
* genotype-associated expression plots.

## Shiny application

The `shinyApp_SCRNA/` directory contains an interactive Shiny application for exploring the scRNA-seq data.

The application includes the RDS objects required for interactive exploration and is provided separately from the main analysis scripts.

The application can be launched from the `shinyApp_SCRNA` directory using:

```r
shiny::runApp("app.R")
```

## Software and reproducibility

The analysis was performed using R and the packages documented in the repository's software-version record.

Key packages used in the workflow include:

* **Seurat**
* **SCTransform**
* **clustifyr**
* **clustifyrdata**
* **MAST**
* **scProportionTest**
* **dittoSeq**
* **tidyverse**
* **scCustomize**

Software versions used for the analysis are documented in:

`Sotware versions`

The analysis scripts are provided in sequential order and are intended to document the computational workflow used to generate the reported scRNA-seq results.

## Data availability

The raw scRNA-seq data generated in this study are available from the European Nucleotide Archive under:

**PRJEB85413**

Bulk RNA-seq data generated in the study are available under:

**PRJEB85309**

Supplementary data associated with the publication are available through the journal article.

## Citation

If you use this repository or the associated scRNA-seq analysis, please cite:

Surakhy M, Matheson J, Barnes DJ, et al. *Smad4 and TGFβ1 dependent gene expression signatures in conditional intestinal adenoma, organoids and colorectal cancer.* Scientific Reports. 2025;15:16330. doi:10.1038/s41598-025-00908-4.

