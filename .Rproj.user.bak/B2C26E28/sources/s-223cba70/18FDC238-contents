genes_juszczak_df %>% 
  rename(category = "signature_name") %>% 
  select(hgnc_symbol, signature_name)

top50_genes_sum_rs %>% 
  mutate(signature_name = paste0(simple_tissue, "_", regulation)) %>% 
  ungroup %>% 
  select(hgnc_symbol, signature_name)

cell_top50_genes_sum_rs %>% 
  mutate(signature_name = paste0(simple_tissue, "_", regulation)) %>% 
  ungroup %>% 
  select(hgnc_symbol, signature_name)


rbind(top50_per_tissue, top50_per_celltype) %>%
  unnest(data) %>%
  group_by(simple_tissue, regulation, hgnc_symbol) %>%
  summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(simple_tissue, regulation) %>%
  arrange(desc(sum_rs)) %>%
  slice_head(n = 50) %>% 
  mutate(signature_name = paste(simple_tissue, regulation, sep = "_")) %>% 
  ungroup %>% 
  select(hgnc_symbol, signature_name) %>% 
  mutate(signature_derivation = "tissue_cell")

# universal
top50_genes_sum_rs %>% 
  unnest(data) %>%
  group_by(simple_tissue, regulation, hgnc_symbol) %>%
  summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(simple_tissue, regulation) %>%
  arrange(desc(sum_rs)) %>%
  nest %>% 
  mutate(data = map(data, ~ .x %>% 
                      arrange(desc(sum_rs)) %>% 
                      mutate(second_rs = row_number(sum_rs)))) %>%
  unnest(data) %>% 
  ungroup() %>% 
  group_by(hgnc_symbol, regulation) %>% 
  nest() %>% 
  mutate(second_sum_rs = map(data, ~ .x %>% .$second_rs %>% sum)) %>% 
  unnest(second_sum_rs) %>% 
  arrange(desc(second_sum_rs)) %>% 
  ungroup() %>% 
  group_by(regulation) %>% 
  slice_max(second_sum_rs, n = 50) %>% 
  select(regulation, hgnc_symbol) %>% 
  ungroup %>% 
  mutate(signature_name = paste0("universal_", regulation)) %>% 
  select(hgnc_symbol, signature_name)



hgnc_new_regulation %>%
  mutate(
    new_regulation_group = case_when(
      new_regulation_perc %in% c("up", "up_weak") ~ "frequently_up",
      new_regulation_perc %in% c("down", "down_weak") ~ "frequently_down",
      new_regulation_perc == "both" ~ "frequently_bidirectional",
      TRUE ~ NA_character_ # warto dodać domyślną opcję na wszelki wypadek
    )
  ) %>% 
  ungroup %>% 
  select(hgnc_symbol, new_regulation_group) %>% 
  rename(new_regulation_group = "signature_name")



# Combine multiple gene signature datasets into one dataframe
gr_genes_signatures_multi_approach_df <- bind_rows(
  # Rename the category column to signature_name and select relevant columns
  genes_juszczak_df %>% 
    rename(category = "signature_name") %>%
    select(hgnc_symbol, signature_name) %>% 
    mutate(signature_derivation = "juszczak_PMID:35405299"),
  
  # # Create a signature_name by combining simple_tissue and regulation columns
  # top50_genes_sum_rs %>%
  #   mutate(signature_name = paste0(simple_tissue, "_", regulation)) %>%
  #   ungroup() %>%
  #   select(hgnc_symbol, signature_name) %>% 
  #   mutate(signature_derivation = "tissue_cell"),
  # 
  # # Create a signature_name by combining simple_tissue and regulation columns
  # cell_top50_genes_sum_rs %>%
  #   mutate(signature_name = paste0(simple_tissue, "_", regulation)) %>%
  #   ungroup() %>%
  #   select(hgnc_symbol, signature_name) %>% 
  #   mutate(signature_derivation = "tissue_cell"),
  
  rbind(top50_per_tissue, top50_per_celltype) %>%
    unnest(data) %>%
    group_by(simple_tissue, regulation, hgnc_symbol) %>%
    summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
    group_by(simple_tissue, regulation) %>%
    arrange(desc(sum_rs)) %>%
    slice_head(n = 50) %>% 
    mutate(signature_name = paste(simple_tissue, regulation, sep = "_")) %>% 
    ungroup %>% 
    select(hgnc_symbol, signature_name) %>% 
    mutate(signature_derivation = "tissue_cell"),
  
  # add universal signatures gr genes
  
  # universal
  top50_genes_sum_rs %>% 
    unnest(data) %>%
    group_by(simple_tissue, regulation, hgnc_symbol) %>%
    summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
    group_by(simple_tissue, regulation) %>%
    arrange(desc(sum_rs)) %>%
    nest %>% 
    mutate(data = map(data, ~ .x %>% 
                        arrange(desc(sum_rs)) %>% 
                        mutate(second_rs = row_number(sum_rs)))) %>%
    unnest(data) %>% 
    ungroup() %>% 
    group_by(hgnc_symbol, regulation) %>% 
    nest() %>% 
    mutate(second_sum_rs = map(data, ~ .x %>% .$second_rs %>% sum)) %>% 
    unnest(second_sum_rs) %>% 
    arrange(desc(second_sum_rs)) %>% 
    ungroup() %>% 
    group_by(regulation) %>% 
    slice_max(second_sum_rs, n = 50) %>% 
    select(regulation, hgnc_symbol) %>% 
    ungroup %>% 
    mutate(signature_name = paste0("universal_", regulation)) %>% 
    select(hgnc_symbol, signature_name) %>% 
    mutate(signature_derivation = "universal"),
  
  
  # Define new regulation groups and map them to signature_name
  hgnc_new_regulation %>%
    mutate(
      signature_name = case_when(
        new_regulation_perc %in% c("up", "up_weak") ~ "frequently_up",
        new_regulation_perc %in% c("down", "down_weak") ~ "frequently_down",
        new_regulation_perc == "both" ~ "frequently_bidirectional",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    select(hgnc_symbol, signature_name) %>% 
    mutate(signature_derivation = "frequent")
) %>%
  distinct()  # Remove duplicate entries

# Display the final combined dataframe
gr_genes_signatures_multi_approach_df 


# Convert the dataframe into a list of gene vectors, grouped by signature_name
gr_genes_signatures_multi_approach_list <- gr_genes_signatures_multi_approach_df %>%
  group_by(signature_name) %>%
  summarise(genes = list(hgnc_symbol), .groups = "drop") %>%
  deframe()  # Convert to a named list where names are signature_name

# Display the structure of the list
str(gr_genes_signatures_multi_approach_list)

gr_genes_signatures_multi_approach_list 


gr_database_blocked_gene_lists$marpiech_cluster_dex_letters %>%
  enframe(name = "signature_name", value = "hgnc_symbol") %>%
  unnest(hgnc_symbol) %>% 
  as.data.frame() %>% 
  mutate(signature_derivation = "marpiech_cluster") %>% 
  select(c("hgnc_symbol", "signature_name", "signature_derivation")) %>% 
  rbind(., gr_genes_signatures_multi_approach_df) -> gr_genes_signatures_multi_approach_df
  

nested_gr_genes_signatures_multi_approach_list <- gr_genes_signatures_multi_approach_df %>%
  select(signature_derivation, signature_name, hgnc_symbol) %>%
  split(.$signature_derivation) %>%
  map(~ .x %>%
        select(signature_name, hgnc_symbol) %>%
        split(.$signature_name) %>%
        map(~ .x %>% .$hgnc_symbol))

gr_genes_signatures_multi_approach_df %>%
  bind_rows(
    map_dfr(names(slezak_signatures_list), function(signature) {
      tibble(
        hgnc_symbol = slezak_signatures_list[[signature]],
        signature_name = signature
      )
    }) %>%
      mutate(signature_derivation = "slezak")
  )  %>% 
  write.table("results/gr-signatures/gr-signatures-multi-approach-24.04.2025.tsv", 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
            )

