# Load necessary libraries
install.packages("rentrez")
install.packages("httr")
install.packages("XML")
library(rentrez)
library(httr)
library(XML)

# Function to download a paper given its PMID
download_paper <- function(pmid, file_path) {
  # Fetch the paper's details using Entrez
  paper_xml <- entrez_fetch(db = "pubmed", id = pmid, rettype = "xml")
  
  # Parse the XML content
  paper_details <- xmlParse(paper_xml)
  
  # Extract the URL of the full text (if available)
  # Note: The method of extracting the URL will depend on the structure of the XML and the availability of the full text
  # This is a placeholder and may need to be adjusted
  full_text_nodes <- getNodeSet(paper_details, "//FullTextURL")
  
  if (length(full_text_nodes) > 0) {
    full_text_url <- xmlValue(full_text_nodes[[1]])
    
    # Check if the URL is available
    if (!is.null(full_text_url) && full_text_url != "") {
      # Download the PDF
      GET(url = full_text_url, write_disk(file_path, overwrite = TRUE))
      cat("Downloaded paper:", pmid, "to", file_path, "\n")
    } else {
      cat("Full text URL not found for PMID:", pmid, "\n")
    }
  } else {
    cat("No full text URL available for PMID:", pmid, "\n")
  }
}

# Example usage
pmid <- "30478444" # Replace with the actual PMID
file_path <- "data/papers/psychiatric-genomic-consortium/paper.pdf" # Replace with your desired file path
download_paper(pmid, file_path)
