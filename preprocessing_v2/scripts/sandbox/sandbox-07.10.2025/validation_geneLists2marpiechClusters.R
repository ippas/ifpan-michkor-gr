# ##############################################################################
# ---- uses data ----
# ##############################################################################

AllBrainGeneLists %>% length()

AllBrain2BrainSignatures2GlobalMarpiechClusters[c(str_detect(names(AllBrain2BrainSignatures2GlobalMarpiechClusters), "^cluster.*"))]

AllBrain2BrainSignatures2GlobalMarpiechClusters[c(str_detect(names(AllBrain2BrainSignatures2GlobalMarpiechClusters), "^cluster.*"))]

AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P %>% length()

AllBrainGeneLists$`michkor-cells_NA_astrocyte_dexamethasone_NA_NA_NA_in-vitro_NA_up`

intersect(AllBrainGeneLists$`michkor-cells_NA_astrocyte_dexamethasone_NA_NA_NA_in-vitro_NA_up`, AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P)


# ##############################################################################
# ---- functions ----
# ##############################################################################
intersect_with_all_lists <- function(gene_lists, cluster_genes) {
  result <- lapply(names(gene_lists), function(list_name) {
    genes_in_list <- gene_lists[[list_name]]
    common_genes <- intersect(genes_in_list, cluster_genes)
    data.frame(
      list_name = list_name,
      n_overlap = length(common_genes),
      n_list = length(genes_in_list),
      percent_overlap = round(100 * length(common_genes) / length(genes_in_list), 2),
      genes = paste(common_genes, collapse = "|")
    )
  })
  dplyr::bind_rows(result)
}


res <- intersect_with_all_lists(
  gene_lists = AllBrainGeneLists,
  cluster_genes = AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P
)

res %>% class

res %>% 
  filter(genes != "") %>% 
  .$genes %>% 
  str_split(pattern = "\\|") %>% 
  unlist %>% table

res %>%   filter(!grepl("FKBP5", genes)) %>% 
  filter(grepl("_up", list_name))

  filter(grepl("FKBP5|DDIT4|TSC22D3", genes)) %>% 
  filter(!grepl("FKBP5", genes))
