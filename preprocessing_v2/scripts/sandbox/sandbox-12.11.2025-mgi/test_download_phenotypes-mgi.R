install.packages("ontologyIndex")

library(ontologyIndex)

# 1. pobierz plik
url <- "http://purl.obolibrary.org/obo/mp.obo"
download.file(url, "mp.obo")

# 2. wczytaj
mp <- get_ontology("mp.obo", extract_tags = "everything")

# 3. struktura
df <- data.frame(term_id = mp$id,
                 name = mp$name,
                 parents = sapply(mp$parents, paste, collapse=";"))

# 4. wybór top-level
top_levels <- subset(df, parents == "MP:0000001")

top_levels


# 2️⃣ Zdefiniuj kategorię, dla której chcesz pobrać wszystkie fenotypy
target_category <- "MP:0005386"  # nervous system phenotype

# 3️⃣ Pobierz wszystkie potomne terminy (rekurencyjnie)
descendants <- get_descendants(mp, target_category)

# 4️⃣ Utwórz ramkę danych z nazwami
mp_df <- data.frame(
  term_id = descendants,
  name = mp$name[descendants],
  stringsAsFactors = FALSE
)

# 5️⃣ Dodaj nazwę kategorii nadrzędnej
mp_df$category <- mp$name[[target_category]]

# 6️⃣ Posortuj alfabetycznie
mp_df <- mp_df %>% arrange(name)

