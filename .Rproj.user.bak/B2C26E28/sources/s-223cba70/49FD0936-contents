enrichr_download_gene_lists <- function(gene_list, database, grep_term = NULL, filter_gene = NULL) {
  # 1. Split the gene list into two parts
  n <- length(gene_list)
  part1 <- gene_list[1:ceiling(n / 2)]
  part2 <- gene_list[(ceiling(n / 2) + 1):n]
  
  # 2. Query Enrichr for each part separately
  enr1 <- enrichr(part1, databases = database)[[database]]
  enr2 <- enrichr(part2, databases = database)[[database]]
  
  # 3. Select only relevant columns
  df1 <- enr1 %>% select(Term, Genes, Overlap)
  df2 <- enr2 %>% select(Term, Genes, Overlap)
  
  # 4. Combine results and construct combined gene lists
  combined <- bind_rows(df1, df2) %>%
    group_by(Term) %>%
    summarise(
      Combined_Genes = paste(unique(unlist(strsplit(Genes, ";"))), collapse = "; "),
      Overlap = unique(Overlap)[1]
    ) %>%
    ungroup()
  
  # 5. Add new columns: number of overlapping genes, reformatted overlap, database name
  combined <- combined %>%
    mutate(
      n_genes_overlap = sapply(strsplit(Combined_Genes, ";\\s*"), function(x) length(unique(trimws(x)))),
      n_all_genes = as.numeric(sub(".*/", "", Overlap)),
      Overlap = paste0(n_genes_overlap, "/", n_all_genes),
      database = database
    ) %>%
    select(Term, Overlap, n_genes_overlap, n_all_genes, Combined_Genes, database)
  
  # 6. Optional filtering by term name
  if (!is.null(grep_term)) {
    combined <- combined %>%
      filter(grepl(grep_term, Term, ignore.case = TRUE))
  }
  
  # 7. Optional filtering by specific gene
  if (!is.null(filter_gene)) {
    combined <- combined %>%
      filter(grepl(filter_gene, Combined_Genes, ignore.case = TRUE))
  }
  
  # 8. Convert column names to lowercase
  names(combined) <- tolower(names(as.data.frame(combined)))
  
  return(combined)
}
