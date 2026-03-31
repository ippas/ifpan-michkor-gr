# Dane
associated_tissue <- 33
total_tissue <- 1012
associated_global <- 0
total_global <- 327

# Tabela 2x2
tbl <- matrix(
  c(associated_tissue,
    total_tissue - associated_tissue,
    associated_global,
    total_global - associated_global),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("Tissue", "Global"),
    c("Associated", "NotAssociated")
  )
)

# Test Fishera
fisher_result <- fisher.test(tbl)

# Test chi² (bez poprawki Yatesa)
chi_result <- chisq.test(tbl, correct = FALSE)

fisher_result
chi_result
chi_result$expected
