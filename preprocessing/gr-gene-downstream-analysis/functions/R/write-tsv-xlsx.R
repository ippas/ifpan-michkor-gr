write_tsv_xlsx <- function(data, tsv_file){
  tsv_file
  xlsx_file <- sub("\\.tsv$", ".xlsx", tsv_file)
  
  
  write.table(data,
              file = tsv_file, 
              sep = "\t", 
              row.names = F,
              col.names = T,
              quote = F)
  
  write_xlsx(data, xlsx_file)
}