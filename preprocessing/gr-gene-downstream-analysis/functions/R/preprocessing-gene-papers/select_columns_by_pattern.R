select_columns_by_pattern <- function(data, pattern, sort_columns = FALSE) {
  # Check if the data is a matrix and convert it to a data frame
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  
  # Get the names of columns that match the pattern
  matching_columns <- grep(pattern, names(data), value = TRUE)
  
  # Sort the matching column names alphabetically if sort_columns is TRUE
  if (sort_columns) {
    matching_columns <- sort(matching_columns)
  }
  
  # Select and return columns that match the pattern
  return(data[, matching_columns])
}
