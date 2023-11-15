# Overlapping GR genes documentation
## 1. Gene Database preprocessing for glucocorticoid-dependent Genes
Before the analysis was initialed, a data cleaning process for the glucocorticoid-dependent genes database was conducted using the [gene-database-preprocessing.R](https://github.com/ippas/ifpan-michkor-gr/blob/overlapping-gr-genes/preprocessing/overlapping-gr-genes/gene-database-preprocessing.R) script. Entries without an assigned HGNC (HUGO Gene Nomenclature Committee) symbol were removed. A vector con# 1. Gene Database preprocessing for glucocorticoid-dependent Genes
Before the analysis was initialed, a data cleaning process for the glucocorticoid-dependent genes database was conducted using the [gene-database-preprocessing.R](https://github.com/ippas/ifpan-michkor-gr/blob/overlapping-gr-genes/preprocessing/overlapping-gr-genes/gene-database-preprocessing.R) script. Entries without an assigned HGNC (HUGO Gene Nomenclature Committee) symbol were removed. A vector containing HGNC symbols for protein-coding genes was then generated, using data sourced from BioMart version 110 (n = 19437). This allowed for the removal of genes from the database that are not responsible for protein encoding. 
Based on microarray data, genes were subsequently clustered. After this classification 15 HGNC symbols that were associated with more than one gene and fell into different clusters were eliminated. This step was taken to avoid situations where identical HGNC symbols appear across multiple clusters. When several genes or transcripts shared the same HGNC but were categorized within the same cluster, the one with the lowest p-value ANOVA was chosen

The following HGNC symbols were removed from the database as they were described to more than one cluster:
- ARRDC2
- C12orf75
- CDC37L1
- DOC2B
- FAM117B
- HPX
- HR
- KIAA0513
- MBD1
- NDRG1
- NFKBIA
- RIN2
- SETDB2
- KLK2
- SERPINA3


# 2. Analysis of Genes Database

### 2.1. Overlap data from microarray vs papers list

Analysis done using [overlap-dextis-papers.R](https://github.com/ippas/ifpan-michkor-gr/blob/overlapping-gr-genes/preprocessing/overlapping-gr-genes/overlap-dextis-papers.R) script. 
Based on prepared database of glucocorticoid-dependent genes, a list of genes that appeared in at last six datasets from articles was compiled and added  to the gene lists. In the first comparison, gene lists from papers were compared to groups of genes divided into clusters and genes expressed in tissues obtained from microarrays. For this purpose, chi-square tests were conducted in which the overlap of paired gene lists was compared. As a background in the tests, HGNC symbols prepared from Biomart v110 were used. The use of a comprehensive list of all genes ensures that the statistical analysis takes into account the overall genomic context and the relative sizes of gene sets. In the tests conducted, to minimize obtaining significant results by chance, those gene lists that had fewer than 5 genes were removed. On the obtained results, FDR correction was applied, and a significance threshold of 0.05 was set. Additionally, to reduce the impact of chance, results that had at least two overlapping genes between two lists were considered significant

### 2.2. Overlap dextis vs metabolism list

In the analysis of gene overlap associated with metabolic processes, gene lists from two databases, [Metabolon](https://www.omicspred.org/Scores/Metabolon/INTERVAL) and [Nightingale](https://www.omicspred.org/Scores/Nightingale/INTERVAL), were utilized. A total of 841 gene lists were examined (701 from Metabolon and 138 from Nightingale). The analyses were conducted using the [overlap-dextis-metabolism.R](https://github.com/ippas/ifpan-michkor-gr/blob/overlapping-gr-genes/preprocessing/overlapping-gr-genes/overlap-dextis-metabolism.R) script. Before executes the main analysis, gene lists containing fewer than three genes related to a metabolic process were eliminated. This step was necessary because too short lists could lead to statistical significance by chance, generating falsely positive associations. After this selection, 654 metabolic gene lists remained (517 from Metabolon and 137 from Nihtingale).

On the selected data, chi-square tests were performed, using HGNC symbols (n = 19437) from Biomart v110 as a background. Two type of comparisons were made using the metabolic gene lists. The first compared lists of genes with similar expression profiles among the tissues studied. The second compared lists of genes expressed in the studied tissues. On the obtained results, a correction for multiple testing (FDR) was applied. An FDR threshold of < 0.05 was used for significance. From the statistically significant results, four metabolic processes without specific assigned names (marked as "X-11564", "X-11261", "X-21470", "X-21467") were excluded. In cases where several biochemical processes had the same list of overlapping genes, only one list with the lowest FDR correction value was chosen for visualization. 

### 2.3. Overlap dextis vs phenotypes

Analysis done using [overlap-dextis-phenotypes-biobank.R](https://github.com/ippas/ifpan-michkor-gr/blob/overlapping-gr-genes/preprocessing/overlapping-gr-genes/overlap-dextis-phenotypes-biobank.R) script.
In the phenotypic analysis, a genetic model based on phenotypes from the [PAN-UK Biobank](https://pan.ukbb.broadinstitute.org/) was used. These models were developed as part of an earlier project and are available on [GitHub](https://github.com/ippas/imdik-zekanowski-sportwgs). For comparisons, only those genetic models were selected that were based on at least fifteen genetic variants, with a GWAS significance threshold of 10e-8 (n = 729). The analysis included performing chi-squere tests in two arrangements: (I) comparison with lists of genes divided according to expression in tissue types. HGNC symbols from Biomart, version v110, were used as a reference background in the tests. The results were subjected to correction for multiple tests using the FDR adjustment. A significance threshold of FDR < 0.05 was adopted. Results where the number of overlapping genes was less than two were removed. In cases where several phenotypes had an identical list of overlapping genes, only one phenotype was chosen for visualization purposes. The criterion of selection was the lowest FDR correction value for the given phenotype.
