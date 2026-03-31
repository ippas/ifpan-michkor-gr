library(httr)
library(jsonlite)
library(readr)

# 🧠 Funkcja: wyszukuje gen w HPA i pobiera jego stronę HTML (np. /brain)
get_hpa_html <- function(gene_symbol, section = "brain", save_file = TRUE) {
  # 1️⃣ znajdź Ensembl ID przez API
  search_url <- sprintf("https://www.proteinatlas.org/search/%s?format=json",
                        URLencode(gene_symbol, reserved = TRUE))
  res <- GET(search_url)
  stop_for_status(res)
  
  js <- fromJSON(rawToChar(res$content))
  if (length(js) == 0) stop("Nie znaleziono genu: ", gene_symbol)
  df <- as.data.frame(js)
  ensg_col <- intersect(names(df), c("Ensembl","ensembl"))[1]
  ensg <- df[[ensg_col]][1]
  
  # 2️⃣ zbuduj adres URL np. https://www.proteinatlas.org/ENSG00000157514-TSC22D3/brain
  url <- sprintf("https://www.proteinatlas.org/%s-%s/%s",
                 ensg, gene_symbol, section)
  message("➡️ Pobieram: ", url)
  
  # 3️⃣ pobierz stronę HTML
  html_resp <- GET(url)
  stop_for_status(html_resp)
  
  html_text <- content(html_resp, "text", encoding = "UTF-8")
  
  # 4️⃣ zapisz lokalnie (jeśli chcesz)
  if (save_file) {
    file_name <- sprintf("HPA_%s_%s.html", gene_symbol, section)
    write_file(html_text, file_name)
    message("💾 Zapisano: ", file_name)
  }
  
  return(html_text)
}

# 🔹 przykład użycia:
html_page <- get_hpa_html("TSC22D3")   # domyślnie pobiera sekcję "brain"


# Wczytaj HTML do obiektu
page <- read_html(html_page)


# pobierz wszystkie tabele
tables <- page %>% html_nodes("table")



# pobierz tabelę
tbl3 <- tables[[3]] %>%
  html_table(header = TRUE, fill = TRUE) %>%
  as_tibble(.name_repair = "unique")  # 👉 automatycznie unikalne nazwy kolumn

# podgląd nazw
names(tbl3)

# czyszczenie
clean_tbl3 <- tbl3 %>%
  rename(col1 = 1, col2 = 2) %>%
  mutate(
    # usuń taby i nowe linie
    across(everything(), ~ str_replace_all(., "[\\t\\n]+", " ")),
    # usuń tekst po 'i' tylko jeśli jest to 'i' od "information tooltip"
    col1 = str_replace(col1, "i\\s.*", ""),
    # przytnij spacje
    col1 = str_trim(col1),
    col2 = str_trim(col2)
  ) %>%
  rename(attribute = col1, value = col2)


clean_tbl3


tbl4 <- tables[[4]] %>%
  rvest::html_table(header = TRUE, fill = TRUE) %>%
  as_tibble(.name_repair = "unique")

clean_tbl4 <- tbl4 %>%
  # automatycznie nadaj unikalne nazwy
  as_tibble(.name_repair = "unique") %>%
  # wybierz tylko pierwsze cztery kolumny (reszta to puste/duplikaty)
  select(1:4) %>%
  rename_with(~ c("col1", "col2", "col3", "col4")) %>%
  # usuń znaki sterujące i tooltipy
  mutate(across(everything(), ~ str_replace_all(., "[\\t\\n]+", " "))) %>%
  mutate(col1 = str_replace(col1, "i\\s.*", "")) %>%
  mutate(across(everything(), str_trim)) %>%
  # usuń pusty wiersz z nazwami podnagłówków
  filter(!col1 %in% c("")) %>%
  # nadaj ostateczne nazwy
  rename(
    attribute = col1,
    human = col2,
    pig = col3,
    mouse = col4
  )

clean_tbl4



# ##############################################################################
# ---- brain expression ----
# ##############################################################################

# ---- FUNKCJA ----
extract_hpa_human_brain_dataset <- function(page, verbose = TRUE) {
  # Znajdź tytuł sekcji
  span_node <- page %>%
    html_elements(xpath = "//span[contains(., 'HPA Human brain dataset')]")
  
  if (length(span_node) == 0) {
    message("❌ Nie znaleziono sekcji 'HPA Human brain dataset' na stronie.")
    return(NULL)
  }
  
  # Przejdź do nadrzędnego <table> (tam siedzą linki <a class='brainrna'>)
  parent_node <- span_node %>%
    html_element(xpath = "./ancestor::table[1]")
  
  # Wyciągnij linki <a class='brainrna'>
  nodes <- parent_node %>% html_nodes("a.brainrna")
  
  if (length(nodes) == 0) {
    message("⚠️ Nie znaleziono żadnych linków <a class='brainrna'> w sekcji.")
    return(NULL)
  }
  
  df <- tibble(
    region = nodes %>% html_text(trim = TRUE),
    brainregion_id = nodes %>% html_attr("brainregion"),
    nTPM = as.numeric(nodes %>% html_attr("nx")),
    color = nodes %>% html_attr("color")
  ) %>%
    filter(!is.na(nTPM)) %>%
    arrange(desc(nTPM))
  
  if (verbose) {
    message("✅ Znaleziono ", nrow(df), " regionów mózgu.")
    print(df)
  }
  
  return(df)
}

human_brain_df <- extract_hpa_human_brain_dataset(page)

# ##############################################################################
# ---- pig expression brain ----
# ##############################################################################

extract_hpa_pig_brain_dataset <- function(page, verbose = TRUE) {
  # znajdź <script> zawierający JSON dla Pig dataset
  script_text <- page %>%
    html_elements("script") %>%
    html_text2() %>%
    str_subset("RNAChart66") %>%
    first()
  
  if (is.na(script_text)) {
    message("❌ Nie znaleziono sekcji 'HPA Pig brain RNA-Seq dataset' na stronie.")
    return(NULL)
  }
  
  # wyodrębnij zawartość JSON z barChart([...])
  json_match <- str_match(script_text, "barChart\\((\\[\\{.*?\\}\\])")[, 2]
  if (is.na(json_match)) {
    message("⚠️ Nie udało się wyodrębnić JSON-a z sekcji Pig.")
    return(NULL)
  }
  
  # parsowanie JSON → tibble
  df <- jsonlite::fromJSON(json_match) %>%
    as_tibble() %>%
    transmute(
      region = label,
      brainregion_id = str_extract(url, "(?<=tissue\\/)[^#]+"),  # poprawione
      nTPM = as.numeric(value),
      color
    ) %>%
    arrange(desc(nTPM))
  
  if (verbose) {
    message("✅ Znaleziono ", nrow(df), " regionów mózgu (Pig).")
    print(df)
  }
  
  return(df)
}


pig_brain_df <- extract_hpa_pig_brain_dataset(page)



# ##############################################################################
# ---- mouse ----
# ##############################################################################
extract_hpa_mouse_brain_dataset <- function(page, verbose = TRUE) {
  # znajdź skrypt zawierający JSON dla Mouse dataset
  script_text <- page %>%
    html_elements("script") %>%
    html_text2() %>%
    str_subset("RNAChart46") %>%
    first()
  
  if (is.na(script_text)) {
    message("❌ Nie znaleziono sekcji 'HPA Mouse brain RNA-Seq dataset' na stronie.")
    return(NULL)
  }
  
  # wyodrębnij zawartość JSON z barChart([...])
  json_match <- str_match(script_text, "barChart\\((\\[\\{.*?\\}\\])")[, 2]
  if (is.na(json_match)) {
    message("⚠️ Nie udało się wyodrębnić JSON-a z sekcji Mouse.")
    return(NULL)
  }
  
  # parsowanie JSON → tibble
  df <- jsonlite::fromJSON(json_match) %>%
    as_tibble() %>%
    transmute(
      region = label,
      brainregion_id = URLdecode(str_extract(url, "(?<=tissue\\/)[^#]+")),
      nTPM = as.numeric(value),
      color
    ) %>%
    arrange(desc(nTPM))
  
  if (verbose) {
    message("✅ Znaleziono ", nrow(df), " regionów mózgu (Mouse).")
    print(df)
  }
  
  return(df)
}


mouse_brain_df <- extract_hpa_mouse_brain_dataset(page)


# ##############################################################################
# ---- HPA stereo-seq cerebral cortex ----
# ##############################################################################
extract_hpa_stereoseq_cerebral_cortex <- function(page, verbose = TRUE) {
  # znajdź <script> z JSON-em dla HPA Stereo-seq
  script_text <- page %>%
    html_elements("script") %>%
    html_text2() %>%
    str_subset("RNAChart105") %>%
    first()
  
  if (is.na(script_text)) {
    message("❌ Nie znaleziono sekcji 'HPA stereo-seq cerebral cortex' na stronie.")
    return(NULL)
  }
  
  # wyodrębnij zawartość JSON-a
  json_match <- str_match(script_text, "barChart\\((\\[\\{.*?\\}\\])")[, 2]
  if (is.na(json_match)) {
    message("⚠️ Nie udało się wyodrębnić JSON-a z sekcji Stereo-seq.")
    return(NULL)
  }
  
  # parsowanie JSON → tibble
  df <- jsonlite::fromJSON(json_match) %>%
    as_tibble() %>%
    transmute(
      cell_type = label,
      enrichment_change = as.numeric(value),
      color
    ) %>%
    arrange(desc(enrichment_change))
  
  if (verbose) {
    message("✅ Znaleziono ", nrow(df), " typów komórek (Stereo-seq cerebral cortex).")
    print(df)
  }
  
  return(df)
}


stereo_df <- extract_hpa_stereoseq_cerebral_cortex(page)

# ##############################################################################
# ---- GTEX ----
# ##############################################################################
extract_hpa_gtex_brain_dataset <- function(page, verbose = TRUE) {
  # znajdź <script> z JSON-em dla GTEx dataset
  script_text <- page %>%
    html_elements("script") %>%
    html_text2() %>%
    str_subset("RNAChart68") %>%
    first()
  
  if (is.na(script_text)) {
    message("❌ Nie znaleziono sekcji 'GTEx Human brain RNA-Seq dataset' na stronie.")
    return(NULL)
  }
  
  # wyodrębnij zawartość JSON z barChart([...])
  json_match <- str_match(script_text, "barChart\\((\\[\\{.*?\\}\\])")[, 2]
  if (is.na(json_match)) {
    message("⚠️ Nie udało się wyodrębnić JSON-a z sekcji GTEx.")
    return(NULL)
  }
  
  # parsowanie JSON → tibble
  df <- jsonlite::fromJSON(json_match) %>%
    as_tibble() %>%
    transmute(
      region = label,
      brainregion_id = URLdecode(str_extract(url, "(?<=tissue\\/)[^#]+")),
      nTPM = as.numeric(value),
      color
    ) %>%
    arrange(desc(nTPM))
  
  if (verbose) {
    message("✅ Znaleziono ", nrow(df), " regionów mózgu (GTEx).")
    print(df)
  }
  
  return(df)
}

gtex_brain_df <- extract_hpa_gtex_brain_dataset(page)

# ##############################################################################
# ---- fantom5 ----
# ##############################################################################
extract_hpa_fantom5_brain_dataset <- function(page, verbose = TRUE) {
  # znajdź <script> z JSON-em dla FANTOM5
  script_text <- page %>%
    html_elements("script") %>%
    html_text2() %>%
    str_subset("RNAChart69") %>%
    first()
  
  if (is.na(script_text)) {
    message("❌ Nie znaleziono sekcji 'FANTOM5 Human brain CAGE dataset' na stronie.")
    return(NULL)
  }
  
  # wyodrębnij zawartość JSON z barChart([...])
  json_match <- str_match(script_text, "barChart\\((\\[\\{.*?\\}\\])")[, 2]
  if (is.na(json_match)) {
    message("⚠️ Nie udało się wyodrębnić JSON-a z sekcji FANTOM5.")
    return(NULL)
  }
  
  # parsowanie JSON → tibble
  df <- jsonlite::fromJSON(json_match) %>%
    as_tibble() %>%
    transmute(
      region = label,
      brainregion_id = URLdecode(str_extract(url, "(?<=tissue\\/)[^#]+")),
      scaled_TPM = as.numeric(value),
      color
    ) %>%
    arrange(desc(scaled_TPM))
  
  if (verbose) {
    message("✅ Znaleziono ", nrow(df), " regionów mózgu (FANTOM5).")
    print(df)
  }
  
  return(df)
}


fantom5_brain_df <- extract_hpa_fantom5_brain_dataset(page)


# ##############################################################################
# ---- mouse brain ish ----
# ##############################################################################

extract_allen_mouse_brain_ish_dataset <- function(page, verbose = TRUE) {
  # Znajdź <script> zawierający dane dla Allen ISH
  script_text <- page %>%
    html_elements("script") %>%
    html_text2() %>%
    str_subset("RNAChart67") %>%
    first()
  
  if (is.na(script_text)) {
    message("❌ Nie znaleziono sekcji 'Allen Mouse brain ISH dataset' na stronie.")
    return(NULL)
  }
  
  # Wyodrębnij zawartość JSON z barChart([...])
  json_match <- str_match(script_text, "barChart\\((\\[\\{.*?\\}\\])")[, 2]
  if (is.na(json_match)) {
    message("⚠️ Nie udało się wyodrębnić JSON-a z sekcji Allen ISH.")
    return(NULL)
  }
  
  # Parsowanie JSON → tibble
  df <- jsonlite::fromJSON(json_match) %>%
    as_tibble() %>%
    transmute(
      region = label,
      brainregion_id = URLdecode(str_extract(url, "(?<=tissue\\/)[^#]+")),
      expression_energy = as.numeric(value),
      color
    ) %>%
    arrange(desc(expression_energy))
  
  if (verbose) {
    message("✅ Znaleziono ", nrow(df), " regionów mózgu (Allen ISH).")
    print(df)
  }
  
  return(df)
}


allen_mouse_ish_df <- extract_allen_mouse_brain_ish_dataset(page)

# ##############################################################################
# ---- get correlated genes ----
# ##############################################################################

clean_hpa_brain_table <- function(tbl_raw, verbose = TRUE) {
  # 1️⃣ nadaj bezpieczne nazwy
  names(tbl_raw) <- make.names(names(tbl_raw), unique = TRUE)
  
  # 2️⃣ usuń kolumny całkowicie puste
  tbl <- tbl_raw %>%
    select(where(~ !all(is.na(.x) | .x == "")))
  
  # 3️⃣ znajdź wiersz z opisem klastra ("is part of")
  cluster_row <- tbl %>%
    filter(if_any(everything(), ~ str_detect(.x, "is part of"))) %>%
    slice_head(n = 1)
  
  if (nrow(cluster_row) == 0) {
    message("❌ Nie znaleziono sekcji 'is part of'.")
    return(NULL)
  }
  
  # 4️⃣ ekstrakcja tekstu i informacji o klastrze
  txt <- paste(cluster_row, collapse = " ")
  cluster_info <- tibble(
    gene = str_extract(txt, "^[A-Z0-9]+"),
    cluster_number = str_extract(txt, "(?<=cluster )\\d+"),
    cluster_description = str_match(txt, "cluster \\d+ (.*?) with")[,2],
    confidence = str_extract(txt, "\\d+\\.\\d+") %>% as.numeric(),
    n_genes_in_cluster = str_extract(txt, "\\d+(?= genes in cluster)") %>% as.numeric()
  )
  
  # 5️⃣ --- prosty i niezawodny sposób na wyciągnięcie genów skorelowanych ---
  correlated_genes <- tbl_raw %>%
    .[-c(1:3), c(1:4)] %>%
    as.data.frame() %>%
    setNames(c("hgnc_symbol", "Description", "Correlation", "Cluster")) %>%
    mutate(
      Correlation = suppressWarnings(as.numeric(Correlation)),
      Cluster = suppressWarnings(as.numeric(Cluster))
    ) %>%
    filter(!is.na(hgnc_symbol) & hgnc_symbol != "")
  
  if (verbose) {
    message(glue::glue("✅ Klaster {cluster_info$cluster_number} ({cluster_info$cluster_description}), {nrow(correlated_genes)} genów skorelowanych."))
  }
  
  list(
    cluster_info = cluster_info,
    correlated_genes = correlated_genes
  )
}
tbl_raw <- tables[[12]] %>%
  rvest::html_table(header = TRUE, fill = TRUE)

clean_hpa_brain_table(tbl_raw)


