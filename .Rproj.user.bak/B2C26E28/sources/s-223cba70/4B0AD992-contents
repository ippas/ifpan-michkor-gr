install.packages("gwasrapidd")
install.packages("stringdist")


# Wczytaj listę nazw z UK Biobank (do przetłumaczenia na EFO / GWAS Catalog)
input_traits <- c("Ever taken cannabis")

# Funkcja dopasowująca nazwę tekstową do najbliższego EFO trait
find_matching_efo <- function(query) {
  traits <- gwasrapidd::get_traits()@traits
  traits$dist <- stringdist::stringdist(
    a = tolower(traits$trait),
    b = tolower(query),
    method = "jw"
  )
  best_match <- traits[order(traits$dist), ][1, ]
  return(best_match[, c("efo_id", "trait", "dist")])
}

lapply(input_traits, find_matching_efo)


# ⏳ Pobieramy dane tylko raz (na początku skryptu)
traits <- gwasrapidd::get_traits()@traits

search_gwas_trait <- function(keyword, max_results = 5) {
  matches <- traits[grepl(keyword, traits$trait, ignore.case = TRUE), ]
  
  if (nrow(matches) > 0) {
    return(head(matches[order(matches$trait), c("efo_id", "trait")], max_results))
  } else {
    message("❌ No matches found.")
    return(NULL)
  }
}

search_gwas_trait("Ever taken cannabis")




GWASCatalog_traits@traits %>% class


setNames(
  lapply(excel_sheets("data/metaphenotypes/PMID38965376_FA_STable5_41562_2024_1909_MOESM5_ESM.xlsx")[-1],
         \(s) read_excel("data/metaphenotypes/PMID38965376_FA_STable5_41562_2024_1909_MOESM5_ESM.xlsx", sheet = s)),
  excel_sheets("data/metaphenotypes/PMID38965376_FA_STable5_41562_2024_1909_MOESM5_ESM.xlsx")[-1]
)
