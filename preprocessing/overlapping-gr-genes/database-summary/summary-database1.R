
clusters_papers_data$gene_list_sizes

intersect(names(papers_mapping), names(clusters_papers_data$gene_list_sizes)) %>%
  as.data.frame() %>%
  set_colnames("original_name") %>% 
  mutate(name = papers_mapping[original_name],
         size = clusters_papers_data$gene_list_sizes[original_name]) %>% 
  mutate(new_name = paste0(name, " ", "(", size, ")")) %>% 
  select(c(original_name, new_name)) %>% 
  mutate(new_name = sapply(new_name, function(s) {
    parts <- strsplit(s, "PMID", fixed = TRUE)[[1]]
    result <- paste0("PMID", parts[2], ", ", parts[1])
    # Remove trailing comma
    result <- gsub(",\\s*$", "", result)
    # Add four spaces at the beginning
    result <- paste0("    ", result)
    result
  }))-> row_mapping_df




intersect(names(papers_mapping), names(clusters_papers_data$gene_list_sizes)) %>%
  as.data.frame() %>%
  set_colnames("original_name") %>% 
  mutate(name = papers_mapping[original_name],
         size = clusters_papers_data$gene_list_sizes[original_name]) %>% 
  mutate(new_name = paste0(name, "|", size)) %>% 
  select(c(original_name, new_name)) %>%
  mutate(tmp = str_replace_all(new_name, " ", ""))  %>% 
  mutate(tmp = str_replace(tmp, ",(?=[^,]*$)", "|"))  %>% 
  mutate(tmp = str_replace(tmp, ",(?=[^,]*$)", "|")) %>% 
  separate(tmp, into = c("tissue", "regulation", "pmid", "size"), sep = "\\|", remove = TRUE)

