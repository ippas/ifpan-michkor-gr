install.packages("umap")
install.packages("Rtsne")

source("preprocessing/functions/R/install-load-packages.R")


read_genes_from_phenotype_models(directory_path = "data/prs-models-pan-biobank-uk/") -> genes_phenotypes_PanUkBiobank

genes_phenotypes_PanUkBiobank %>% 
  unname %>% unlist %>% unique() %>% length()

create_gene_matrix <- function(gene_list) {
  # Identify all unique genes
  all_genes <- unique(unlist(gene_list))
  
  # Initialize matrix with zeros
  gene_matrix <- matrix(0, nrow = length(gene_list), ncol = length(all_genes),
                        dimnames = list(names(gene_list), all_genes))
  
  # Fill the matrix
  for (i in seq_along(gene_list)) {
    # Find which genes are present for this phenotype
    present_genes <- gene_list[[i]]
    # Ensure present_genes is a character vector to avoid indexing issues
    if (!is.character(present_genes)) {
      present_genes <- as.character(present_genes)
    }
    # Set matrix cells to 1 for present genes
    gene_matrix[i, present_genes] <- 1
  }
  
  # Return the matrix
  return(gene_matrix)
}

# Example usage
phenotypes <- list(
  phenotype1 = c("gene1", "gene3"),
  phenotype2 = c("gene2", "gene3", "gene4"),
  phenotype3 = c("gene1", "gene4", "gene5")
)


create_gene_matrix(genes_phenotypes_PanUkBiobank) -> genes_phenotypes_PanUkBiobank_matrix

genes_phenotypes_PanUkBiobank_matrix %>% dim

genes_phenotypes_PanUkBiobank_matrix[1:10, 1:10]

################################################################################
# perform PCA
pca_result <-prcomp(genes_phenotypes_PanUkBiobank_matrix, center = TRUE, scale. = TRUE)
 
# View summary of PCA results
summary(pca_result)

plot(pca_result, type = "l")

# Plot the first two principal components
plot(pca_result$x[,1], pca_result$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA Scatter Plot")
text(pca_result$x[,1], pca_result$x[,2], labels = rownames(pca_result$x), pos = 4, cex = 0.6)

# Convert PCA results to a data frame for ggplot2
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], RowNames = rownames(pca_result$x))

# Create the scatter plot
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point() +
  # geom_text(vjust = "inward", hjust = "inward", angle = 45, check_overlap = TRUE, size = 3) +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") +
  ggtitle("PCA Scatter Plot") +
  theme_minimal()

################################################################################
# UMAP
library(umap)

# Assuming `gene_matrix` is your matrix
umap_result <- umap(genes_phenotypes_PanUkBiobank_matrix)

# Plot UMAP results
plot(umap_result$layout[,1], umap_result$layout[,2], xlab = "UMAP 1", ylab = "UMAP 2", main = "UMAP Projection")

# Convert UMAP results to a data frame
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

umap_df %>% .$UMAP2 %>% min

umap_df %>% filter(UMAP2 > -30) -> umap_df

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.8) + # Optional: color points by label if available
  theme_minimal() +
  labs(title = "UMAP Projection",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2",
       color = "Label") + # Adjust the label as per your data
  scale_color_viridis_d() # Optional: Use a different color scale


################################################################################
library(Rtsne)

# Run t-SNE
tsne_result <- Rtsne(tmp, dims = 2, perplexity = 80, verbose = TRUE)

# Plot t-SNE results
plot(tsne_result$Y[,1], tsne_result$Y[,2], xlab = "t-SNE 1", ylab = "t-SNE 2", main = "t-SNE Projection", asp = 1, col = rainbow(length(unique(tsne_result$Y[,1]))))


################################################################################
# Compute the distance matrix
genes_phenotypes_PanUkBiobank_matrix -> genes_phenotypes_PanUkBiobank_matrix2

rownames(genes_phenotypes_PanUkBiobank_matrix2) <- gsub("biobankuk-", "", rownames(genes_phenotypes_PanUkBiobank_matrix2))



dist_matrix <- dist(genes_phenotypes_PanUkBiobank_matrix2, method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "complete")

# Plot the dendrogram
plot(hc, hang = -1)



# Open PDF device
pdf("results/figures-phenotypes/phenotypes-hclust.pdf", width = 1000, height = 50)

# Plotting (replace with your actual plot code)
plot(hc, hang = -1) # Make sure 'hc' is defined and contains your hierarchical clustering
rect.hclust(hc, k = 10, border = 2:10)
# Close PDF device
dev.off()


sub_grp <- cutree(hc, k = 10)

table(sub_grp)
