# ifpan-michkor-gr/dextis-analysis-paper

Please note that all code is offered "as-is", so use at your own peril!

## Overview of Scripts and Analysis for GR-Dependent Genes

This repository contains scripts and documentation related to the preprocessing, clustering, and downstream analysis of GR-dependent genes. Below is a detailed description of each component:

### 1. Clustering GR-Dependent Genes Across Tissues

The code for clustering GR-dependent genes based on similar expression signatures across different tissues can be found [here](https://github.com/ippas/ifpan-michkor-gr/tree/dextis-analysis-paper/preprocessing/gr-dependent-gene-clusters-definition).

### 2. Preliminary Processing of Results from Individual Studies

Scripts responsible for extracting and preprocessing data from individual studies are available [here](https://github.com/ippas/ifpan-michkor-gr/tree/dextis-analysis-paper/preprocessing/extract-genes-from-papers).

### 3. Aggregating Results into a Unified Database of GR-Dependent Genes

The scripts used to aggregate results from various studies and organize them into a unified database of GR-dependent genes are located [here](https://github.com/ippas/ifpan-michkor-gr/tree/dextis-analysis-paper/preprocessing/prepare-gr-gene-database). A more detailed description of the database preparation process is available in the [README](https://github.com/ippas/ifpan-michkor-gr/blob/dextis-analysis-paper/preprocessing/prepare-gr-gene-database/README.md).

### 4. Downstream Analysis of GR-Dependent Genes

The code used for downstream analysis is located in the [gr-gene-downstream-analysis](https://github.com/ippas/ifpan-michkor-gr/tree/dextis-analysis-paper/preprocessing/gr-gene-downstream-analysis) folder. This analysis utilized both the GR-dependent gene clusters and the unified GR-dependent gene database. Several comparisons and analyses were performed:

#### Overlap Analysis

The following overlap analyses were performed to compare GR-dependent gene signatures with various external datasets:

- **Overlap with GR-dependent gene lists from studies**: Comparison of GR-dependent gene signatures with gene lists derived from individual studies. Script available [here](https://github.com/ippas/ifpan-michkor-gr/blob/dextis-analysis-paper/preprocessing/gr-gene-downstream-analysis/overlapping-gr-genes/papers-figure-paper.R).
- **Overlap with phenotype-associated gene lists from [Pan-UK Biobank](https://pan.ukbb.broadinstitute.org/)**: Comparison of GR-dependent gene signatures with gene lists associated with phenotypes. Script available [here](https://github.com/ippas/ifpan-michkor-gr/blob/dextis-analysis-paper/preprocessing/gr-gene-downstream-analysis/overlapping-gr-genes/phenotypes-figure-paper.R).
- **Overlap with metabolic traits gene lists from [Metabolon](https://www.omicspred.org/Scores/Metabolon/INTERVAL) and [Nightingale](https://www.omicspred.org/Scores/Nightingale/INTERVAL) databases**: Comparison of GR-dependent gene signatures with gene lists related to metabolic traits. Script available [here](https://github.com/ippas/ifpan-michkor-gr/blob/dextis-analysis-paper/preprocessing/gr-gene-downstream-analysis/overlapping-gr-genes/metabolism-figure-paper.R#L16).

#### Enrichr Analysis

- **Analysis using target genes of transcription factors from the [ChEA Transcription Factor Targets 2022 dataset](https://maayanlab.cloud/Harmonizome/dataset/CHEA+Transcription+Factor+Targets+2022)**: Comparison of GR-dependent gene signatures with transcription factor target gene sets. Script available [here](https://github.com/ippas/ifpan-michkor-gr/blob/dextis-analysis-paper/preprocessing/gr-gene-downstream-analysis/enrichr/chea-enrichr-figure-clean.R).

A more detailed description of the downstream analysis is provided [here](https://github.com/ippas/ifpan-michkor-gr/blob/dextis-analysis-paper/preprocessing/gr-gene-downstream-analysis/overlapping-gr-genes.md).
