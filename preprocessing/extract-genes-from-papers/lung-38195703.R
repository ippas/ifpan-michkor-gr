source("preprocessing/functions/R/install-load-packages.R")

data_path <- "data/supplement-genes-papers/lung-38195703/41598_2024_51301_MOESM1_ESM.xlsx"

# Read sheet names and load each sheet into a list of data frames
sheets <- excel_sheets(data_path)
data <- lapply(sheets, function(sheet) read_excel(data_path, sheet = sheet))

data
