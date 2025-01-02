download_and_read_geo_file <- function(file_url, output_dir, separator = "\t") {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract the file name from the URL
  file_name <- basename(file_url)
  
  # Full path for the downloaded file
  dest_file <- file.path(output_dir, file_name)
  
  # Download the file
  download.file(file_url, destfile = dest_file, mode = "wb")
  
  data <- read.table(gzfile(dest_file), header = TRUE, sep = separator,  row.names = 1)
  
  return(data)
}
