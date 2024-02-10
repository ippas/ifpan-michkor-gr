# Define the number of empty character vectors you want in the list
N <- 127

# Define the sizes for each vector
sizes <- c(3, 5, 2, 4, 6)
sizes <- papers_gene_list %>% lapply(., length) %>% unname() %>% unlist

# Create an empty list to store the character vectors
empty_list <- list()

# Generate and name the empty character vectors with specified sizes
for (i in 1:N) {
  empty_list[[paste0("random", i)]] <- vector("character", length = sizes[i])
}

# Print the list
print(empty_list)


# prepare random gene list
papers_gene_list %>% unname() %>% unlist() %>% unique() %>% length() -> number_genes

papers_gene_list %>% unname() %>% unlist() %>% table() %>% as.data.frame() %>% set_colnames(c("gene", "n")) %>% 
  arrange(desc(n)) %>% .$n -> frequency_genes

sample(hgnc_symbols_vector_v110, number_genes) -> random_genes

random_genes[1:10]


# check which vector in the list have empty fields
# Example list of character vectors with empty elements
# Example list of character vectors with empty elements
empty_list <- list(random1 = c("a", "d", "b", "d", "c", "d"),
                   random2 = c("", "", "b", "d", "", "f"),
                   random3 = c("", "", "", "", "", ""))

# Function to check if a vector contains at least one empty element
has_empty <- function(x) any(x == "")

# Find names of vectors with at least one empty element
names_with_empty <- names(Filter(has_empty, empty_list))

# Print names of vectors with at least one empty element
print(names_with_empty)
lapply(empty_list, function(x) which(x == ""))

# Function to replace the first empty element with a specified value in a vector
replace_empty_with_value <- function(x, value) {
  empty_index <- which(x == "")
  if (length(empty_index) > 0) {
    x[empty_index[1]] <- value
  }
  return(x)
}

for(i in 1:number_genes){
  # print(i)
  # print(random_genes[i])
  # print(frequency_genes[i])
  # print("")
  # Function to check if a vector contains at least one empty element
  has_empty <- function(x) any(x == "")
  
  # Find names of vectors with at least one empty element
  names_with_empty <- names(Filter(has_empty, empty_list))
  
  print(names_with_empty)
  
  random_empty_vectors <- sample(names_with_empty, frequency_genes[i])
  print(random_empty_vectors)
  
  for(random_vector in random_empty_vectors) {
    empty_list[[random_vector]] <- replace_empty_with_value(empty_list[[random_vector]], random_genes[i])
  }
}


empty_list %>% 

replace_empty_with_value(empty_list[["random"]], random_genes[i])

create_gene_pie_chart(gene_list = empty_list, threshold_percent = 99)
