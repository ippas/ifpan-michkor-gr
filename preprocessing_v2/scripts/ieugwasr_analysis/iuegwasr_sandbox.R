install.packages("remotes")
remotes::install_github("MRCIEU/ieugwasr")
# 
Sys.setenv(IEU_OPEN_GWAS_TOKEN = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJtYXRldXN6emllYmE5N0BnbWFpbC5jb20iLCJpYXQiOjE3NzE4ODA3MzMsImV4cCI6MTc3MzA5MDMzM30.n7fROUBG300xc2OIOVPSpj3738k_vuVZwu6HCvIzXHj7icoFvPqbcvu6cH8hpYbMQ1DQ1n3gI_R1oEKHC1-cu6YmDUvkNhPaFN9xlGhBp-VjdHBhh8InfpI44JcBfWJvQXSlObDekZJEtEi3OF2U5x3JkMjlfv6162s_K1swUtNbW9RUIhUdxGvkJE-kXdkKPOs_v7lCgj1BdoPRlTz0EvCKqa44VExY6ZEzVJ87zOSbl3zZXUsAjQ0CZ9pJfFL0Dm0vxqjIprfpuaiJMIUOMs6HDbfLszhPQ_iTGm9MEny18PaljAF9dECcsAgaBTvVzKIyOm_tOS5pgy3KbFNt7w")
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJtYXRldXN6emllYmE5N0BnbWFpbC5jb20iLCJpYXQiOjE3NzE4ODA3MzMsImV4cCI6MTc3MzA5MDMzM30.n7fROUBG300xc2OIOVPSpj3738k_vuVZwu6HCvIzXHj7icoFvPqbcvu6cH8hpYbMQ1DQ1n3gI_R1oEKHC1-cu6YmDUvkNhPaFN9xlGhBp-VjdHBhh8InfpI44JcBfWJvQXSlObDekZJEtEi3OF2U5x3JkMjlfv6162s_K1swUtNbW9RUIhUdxGvkJE-kXdkKPOs_v7lCgj1BdoPRlTz0EvCKqa44VExY6ZEzVJ87zOSbl3zZXUsAjQ0CZ9pJfFL0Dm0vxqjIprfpuaiJMIUOMs6HDbfLszhPQ_iTGm9MEny18PaljAF9dECcsAgaBTvVzKIyOm_tOS5pgy3KbFNt7w")
library(ieugwasr)
library(data.table)

# 1) Pobierz metadane WSZYSTKICH datasetów (to może chwilę potrwać)
info2 <- gwasinfo()

dim(info)
names(info)[1:30]
head(info[, c("id","trait")])

info2 %>% head %>% as.data.frame()

subset(info, grepl("bipolar", trait, ignore.case = TRUE)) %>% 
  as.data.frame() %>% 
  filter(trait == "Bipolar disorder")

info %>% 
  as.data.frame() %>% 
  filter(subcategory == "Lipid")
  .$subcategory %>% table


  
  
subset(info, grepl("depress", trait, ignore.case = TRUE))[, c("id","trait","population","sample_size")]
subset(info, grepl("insulin", trait, ignore.case = TRUE))[, c("id","trait","population","sample_size")]




info$subcategory %>% table %>% 
  as.data.frame() %>% 
  set_colnames(c("subcategory", "freq")) %>% 
  filter(freq > 10)


devtools::install_github('MRCIEU/gpmapr')
#to search for a trait, gene, or variant:
# gpmapr::search_gpmapr('Haemoglobin')
#get data for a gene:
gpmapr::gene('APOD') -> tmp


c(
  "ADARB1",
  "ADM",
  "ADRB2",
  "APOD",
  "BCR",
  "DUSP6",
  "FKBP5",
  "FZD4",
  "IL6",
  "MAOA",
  "PER2",
  "PTGDS",
  "RPS6KA2",
  "S100A10",
  "S1PR1",
  "SAT1",
  "SOD2",
  "ST3GAL1",
  "STAB1",
  "TGM2"
)


tmp$rare_results %>% 
  arrange(min_p) %>% 
  filter(situated_gene == "FKBP5") %>% 
  # head(10) %>% 
  filter(!grepl("Whole blood", trait_name, ignore.case = TRUE)) %>% 
  select(min_p, trait_name) %>% head(1)




tmp$variants %>% head

tmp$coloc_groups %>% 
  arrange(min_p) %>% 
  filter(gene == "APOD") %>% 
  # head(10) %>% 
  filter(!grepl("Whole blood", trait_name, ignore.case = TRUE)) %>% 
  head(10)


tmp$rare_results %>% 
  arrange(min_p) %>% 
  filter(situated_gene == "FKBP5") %>% 
  # head(10) %>% 
  filter(!grepl("Whole blood", trait_name, ignore.case = TRUE)) %>% 
  select(min_p, trait_name) %>% head(1)


tmp$study_extractions %>% 
  arrange(min_p) %>% 
  filter(gene == "FKBP5") %>% 
  # head(10) %>% 
  filter(!grepl("Whole blood", trait_name, ignore.case = TRUE)) 




genes <- c(
  "ADARB1","ADM","ADRB2","APOD","BCR","DUSP6","FKBP5","FZD4","IL6",
  "MAOA","PER2","PTGDS","RPS6KA2","S100A10","S1PR1","SAT1",
  "SOD2","ST3GAL1","STAB1","TGM2"
)

genes <- c("ADARB1", "ADM", "ALDH1A1", "APOD", "BTG1", "CD163", "CKB", "CP", "DHCR24", "GLUL", "IL6", "IL6R", "KDR", "LCN2", "LRRFIP1", "MAOA", "PAG1", "PTGS2", "RGS2", "SOD2", "ST3GAL1", "TGM2")

genes <- sample(hgnc_symbols_vector_v110, 22)

results_df <- map_dfr(genes, function(g) {
  print(g)
  
  tmp <- gpmapr::gene(g)
  
  # przypadek: gen nie znaleziony
  if (!is.null(tmp$detail) && grepl("not found", tmp$detail, ignore.case = TRUE)) {
    return(tibble(
      gene = g,
      min_p = NA_real_,
      trait_name = NA_character_
    ))
  }
  
  # przypadek: brak rare_results albo pusty obiekt
  if (is.null(tmp$rare_results) || nrow(tmp$rare_results) == 0) {
    return(tibble(
      gene = g,
      min_p = NA_real_,
      trait_name = NA_character_
    ))
  }
  
  res <- tmp$rare_results %>%
    arrange(min_p) %>%
    filter(situated_gene == g) %>%
    filter(!grepl("Whole blood", trait_name, ignore.case = TRUE)) %>%
    select(min_p, trait_name) %>%
    slice(1)
  
  # przypadek: po filtrach nic nie zostało
  if (nrow(res) == 0) {
    tibble(
      gene = g,
      min_p = NA_real_,
      trait_name = NA_character_
    )
  } else {
    tibble(
      gene = g,
      min_p = res$min_p,
      trait_name = res$trait_name
    )
  }
})

