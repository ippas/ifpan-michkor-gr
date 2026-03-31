# 📦 wymagane pakiety
if (!requireNamespace("ontologyIndex", quietly = TRUE))
  install.packages("ontologyIndex")
library(ontologyIndex)



url <- "https://raw.githubusercontent.com/MaayanLab/Enrichr-Viz-Appyter/master/Enrichr-Processed-Library-Storage/Clustered_Scatterplots/MGI_Mammalian_Phenotype_Level_4_2024.csv"
mgi2024 <- read.csv(url, stringsAsFactors = FALSE)


mp_get_descendants_by_name <- function(category_name, mp = NULL, file = "mp.obo") {
  # Jeśli ontologia nie została przekazana, wczytaj ją
  if (is.null(mp)) {
    if (!file.exists(file)) {
      message("Pobieranie pliku mp.obo...")
      download.file("http://purl.obolibrary.org/obo/mp.obo", file, quiet = TRUE)
    }
    mp <- get_ontology(file, extract_tags = "everything")
  }
  
  # Znajdź ID kategorii po nazwie (case-insensitive)
  category_id <- mp$id[tolower(mp$name) == tolower(category_name)]
  if (length(category_id) == 0) {
    stop("Nie znaleziono kategorii o nazwie: ", category_name)
  }
  
  # Pobierz potomne terminy
  descendants <- get_descendants(mp, category_id)
  
  # Utwórz ramkę danych bez nazw wierszy
  df <- data.frame(
    term_id = descendants,
    name = mp$name[descendants],
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL
  
  message("✅ Znaleziono ", nrow(df), " fenotypów dla kategorii: ", category_name)
  return(df)
}

behaviorNeurological_mgi <- mp_get_descendants_by_name(category_name = "behavior/neurological phenotype")

mgi2024 %>% 
  select(c(term, genes)) %>% 
  mutate(mp_id = str_extract(term, "MP:\\d{7}")) %>% 
  select(c(term, mp_id, genes)) %>% 
  filter(mp_id %in% behaviorNeurological_mgi$term_id)
