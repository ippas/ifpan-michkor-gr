# Process data: filter for relevant regulations, select necessary columns,
# remove rows with missing values, and eliminate duplicate entries
venn_data <- filt_1_12_gr_database %>% 
  filter(regulation %in% c("up", "down")) %>%
  select(hgnc_symbol, treatment_type) %>% 
  drop_na() %>% 
  distinct()

# Helper function to extract genes for a specific treatment type
get_genes_by_treatment <- function(data, treatment) {
  data %>% 
    filter(treatment_type == treatment) %>% 
    pull(hgnc_symbol)
}

# Create a list of gene vectors for each treatment type
venn_list <- list(
  Chronic = get_genes_by_treatment(venn_data, "chronic"),
  Acute = get_genes_by_treatment(venn_data, "acute"),
  Stress_Induction = get_genes_by_treatment(venn_data, "stress_induction")
)

# Draw the Venn diagram
venn_plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = rep("white", length(venn_list)),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("blue", "red", "green"),
  margin = 0.1
)

# Display the Venn diagram
grid.newpage()

svg("results/figures/gene-database-summary/filtered-1-12-gene-database-summary/venn-treatment-type.svg", width = 8, height = 8)
grid.draw(venn_plot)
dev.off()
