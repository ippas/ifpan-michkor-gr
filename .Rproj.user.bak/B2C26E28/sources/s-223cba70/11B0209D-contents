# Przygotowanie danych
df <- filt_1_12_gr_database %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(hgnc_occurence > 35) %>%  
  filter(!is.na(hgnc_symbol)) %>% 
  select(label, hgnc_symbol)

# Tworzymy macierz: wiersze = label, kolumny = hgnc_symbol
gene_label_matrix <- table(df$label, df$hgnc_symbol)

# Zamieniamy na macierz binarną: jeśli gen występuje (czyli liczba > 0), ustawiamy 1, w przeciwnym razie 0
gene_label_matrix_binary <- (gene_label_matrix > 0) * 1

# Sprawdzamy wynik
print(gene_label_matrix_binary)

# Funkcja obliczająca współczynnik Jaccarda dla dwóch wektorów binarnych
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  if (union == 0) return(NA) else return(intersection / union)
}

# Uzyskujemy nazwy genów (kolumn) z macierzy
genes <- colnames(gene_label_matrix_binary)

# Inicjujemy macierz do przechowania wyników
jaccard_matrix <- matrix(NA, nrow = length(genes), ncol = length(genes), 
                         dimnames = list(genes, genes))

# Obliczamy współczynnik Jaccarda dla każdej pary genów
for (i in seq_along(genes)) {
  for (j in seq_along(genes)) {
    jaccard_matrix[i, j] <- jaccard_similarity(gene_label_matrix_binary[, i],
                                               gene_label_matrix_binary[, j])
  }
}

# Wyświetlamy macierz
print(jaccard_matrix)