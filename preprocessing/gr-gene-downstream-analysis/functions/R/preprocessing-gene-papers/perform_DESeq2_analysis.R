perform_DESeq2_analysis <- function(expression_data, pvalue_threshold = 0.05, log2fc_threshold = 2, condition_vector) {
  library(DESeq2)
  library(dplyr)
  
  # Create sample information data frame
  sample_info <- data.frame(
    sample = colnames(expression_data),
    condition = condition_vector
  )
  
  print(sample_info)
  
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = expression_data,
                                colData = sample_info,
                                design = ~ condition)
  
  # Perform DESeq2 analysis
  dds <- DESeq(dds)
  
  # Get differential expression results and filter
  significant_genes <- results(dds) 
  
  
  significant_genes %>% as.data.frame %>% filter(padj < pvalue_threshold, abs(log2FoldChange) > log2(log2fc_threshold)) -> significant_genes
  
  significant_genes %>% 
    select(-c(baseMean, lfcSE, stat)) %>% 
    set_colnames(c("log2ratio", "pvalue", "fdr")) %>% 
    as.data.frame() %>% 
    rownames_to_column(., var = "ensembl_id") %>% 
    mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> significant_genes
  
  return(significant_genes)
}