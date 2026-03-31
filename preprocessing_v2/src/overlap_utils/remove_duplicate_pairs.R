remove_duplicate_pairs <- function(df,
                                   col_a = "Var1",
                                   col_b = "Var2",
                                   remove_self_pairs = TRUE) {
  # Zamiana nazw kolumn na symbole dla dplyr
  col_a_sym <- rlang::sym(col_a)
  col_b_sym <- rlang::sym(col_b)
  
  df_filtered <- df %>%
    dplyr::mutate(
      pair_id = ifelse(
        !!col_a_sym < !!col_b_sym,
        paste(!!col_a_sym, !!col_b_sym, sep = "_"),
        paste(!!col_b_sym, !!col_a_sym, sep = "_")
      )
    ) %>%
    dplyr::distinct(pair_id, .keep_all = TRUE)
  
  # Opcjonalne usuwanie par, gdzie A == B
  if (remove_self_pairs) {
    df_filtered <- df_filtered %>%
      dplyr::filter(!!col_a_sym != !!col_b_sym)
  }
  
  df_filtered %>%
    dplyr::select(-pair_id)
}


remove_duplicate_pairs <- function(df,
                                   col_a = "Var1",
                                   col_b = "Var2",
                                   remove_self_pairs = TRUE) {
  col_a_sym <- rlang::sym(col_a)
  col_b_sym <- rlang::sym(col_b)
  
  df_filtered <- df %>%
    dplyr::mutate(
      a_chr = as.character(!!col_a_sym),
      b_chr = as.character(!!col_b_sym),
      pair_id = ifelse(
        a_chr < b_chr,
        paste(a_chr, b_chr, sep = "_"),
        paste(b_chr, a_chr, sep = "_")
      )
    ) %>%
    dplyr::distinct(pair_id, .keep_all = TRUE)
  
  if (remove_self_pairs) {
    df_filtered <- df_filtered %>%
      dplyr::filter(a_chr != b_chr)
  }
  
  df_filtered %>%
    dplyr::select(-pair_id, -a_chr, -b_chr)
}
