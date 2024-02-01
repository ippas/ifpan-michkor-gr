# Function to read Excel sheets with an option to skip rows
read_excel_sheets <- function(file_path, skip_rows = 0) {
  # Read the names of the Excel sheets from the given file
  sheets <- excel_sheets(file_path)
  
  # Load all sheets from the Excel file into a list of data frames
  data <- lapply(sheets, function(sheet) {
    read_excel(file_path, sheet = sheet, skip = skip_rows)
  })
  
  return(data)
}
