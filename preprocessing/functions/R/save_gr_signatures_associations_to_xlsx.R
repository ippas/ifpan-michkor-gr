save_gr_signature_associations_to_xlsx <- function(data_list, output_file) {
  #' Save GR gene signature associations to an Excel file (fixed 4-sheet version)
  #'
  #' This function writes GR-related gene association data to an Excel file,
  #' using a fixed set of four gene signatures: metasignature_up, metasignature_down,
  #' brain_up, and brain_down. Each signature is saved to a separate worksheet.
  #'
  #' @param data_list A named list of data frames with association results,
  #'        e.g., randomization_results$original_genebass_association
  #' @param output_file A character string specifying the output .xlsx file path
  #'
  #' @return Invisible NULL (writes file to disk)
  #'
  #' @examples
  #' save_gr_signature_associations_to_xlsx(
  #'   data_list = randomization_results$original_genebass_association,
  #'   output_file = "results/gr_signature_associations.xlsx"
  #' )
  
  stopifnot(is.list(data_list), is.character(output_file))
  
  # Define the 4 GR signature sheet names
  sheet_names <- c("metasignature_up", "metasignature_down", "brain_up", "brain_down")
  
  # Create a new workbook
  wb <- openxlsx::createWorkbook()
  
  # Write each signature to a separate worksheet
  for (sheet in sheet_names) {
    df <- data_list[[sheet]]
    openxlsx::addWorksheet(wb, sheetName = sheet)
    openxlsx::writeData(wb, sheet = sheet, x = df, withFilter = TRUE)
    openxlsx::freezePane(wb, sheet = sheet, firstRow = TRUE)
  }
  
  # Save the workbook to file
  openxlsx::saveWorkbook(wb, file = output_file, overwrite = TRUE)
  invisible(NULL)
}
