# Function to install and load packages
install_and_load_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)
  sapply(packages, require, character.only = TRUE)
}

# List of packages to be installed and loaded
required_packages <- c("readxl", "writexl", "magrittr", "tidyverse")

# Install and load the packages
install_and_load_packages(required_packages)
