compute_rank_score_log2ratioPlusFDR <- function(data, n_top = 100){
  data %>% names() %>% str_remove_all(., pattern = "_fdr|_log2ratio") %>% unique() -> pattern_list
  
  results_rs_log2ratioPlusFDR <- lapply(pattern_list, function(name){
    tryCatch({
      name_fdr <- paste0(name, "_fdr")
      name_log2ratio <- paste0(name, "_log2ratio")
      
      intersect(
        data[[name_log2ratio]]$metadata$label,
        data[[name_fdr]]$metadata$label
      ) -> lists_names
      
      data[[name_fdr]]$genes_sum_score_df %>% 
        ungroup %>% select(n_papers_raw, n_lists_raw, n_papers_filt, n_lists_filt) %>% unique() -> summary_fdr
      
      data[[name_log2ratio]]$genes_sum_score_df %>% 
        ungroup %>% select(n_papers_raw, n_lists_raw, n_papers_filt, n_lists_filt) %>% unique() -> summary_log2ratio
      
      rbind(summary_log2ratio, summary_fdr) %>% 
        summarise(across(everything(), min, na.rm = TRUE)) -> summary_minLog2ratioFDR
      
      print(lists_names)
      result_lists <- lapply(lists_names, function(list_name){
        
        rs_log2ratio <- data[[name_log2ratio]]$ranked_scores_data[[list_name]] %>% 
          select(-c(log2ratio, rank))
        rs_fdr <- data[[name_fdr]]$ranked_scores_data[[list_name]] %>% 
          select(-c(fdr, rank))
        
        rbind(rs_log2ratio, rs_fdr) %>% 
          group_by(hgnc_symbol) %>% 
          nest() %>% 
          mutate(score = map(data, ~sum(.x$score))) %>% 
          select(-data) %>% unnest()
      }) %>% bind_rows() %>% 
        group_by(hgnc_symbol) %>% 
        nest() %>% 
        mutate(score = map(data, ~sum(.x$score))) %>% 
        select(-data) %>% 
        unnest(score) %>% arrange(desc(score)) %>% head(n_top) %>% 
        mutate(
          n_papers_raw = summary_minLog2ratioFDR$n_papers_raw,
          n_lists_raw = summary_minLog2ratioFDR$n_lists_raw,
          n_papers_filt = summary_minLog2ratioFDR$n_papers_filt,
          n_lists_filt = summary_minLog2ratioFDR$n_lists_filt
        )
      
      return(result_lists)
    }, error = function(e){
      cat("Error in processing pattern:", name, "- returning NULL.\n")
      return(NULL)
    })
  })
  
  names(results_rs_log2ratioPlusFDR) <- pattern_list %>% paste0(., "_log2ratioPlusFDR")
  
  # Remove NULL elements from the list
  results_rs_log2ratioPlusFDR <- Filter(Negate(is.null), results_rs_log2ratioPlusFDR)
  
  return(results_rs_log2ratioPlusFDR)
}