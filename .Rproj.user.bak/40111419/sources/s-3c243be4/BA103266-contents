save_data_to_tsv <- function(data_list, path) {
  # Create the directory if it doesn't exist
  if (!dir.exists(path)) {
    dir.create(path)
  }
  
  # Loop through each data frame in the list
  for (i in seq_along(data_list)) {
    # Get the current data frame
    current_data <- data_list[[i]]
    current_name <- names(data_list)[i]
    
    print(names(current_data))
    # Get the name for the file from the 'name' column
    filename <- paste0(current_name, ".tsv")
    print(filename)
    
    # Create the full path for the file
    full_path <- file.path(path, filename)
    
    # Write the data frame to a TSV file
    write.table(current_data, file = full_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = T)
  }
}
