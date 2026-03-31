# ============================================================
# Required packages for DisGeNET category → gene extraction
# ============================================================
message("📦 Loading required packages...")
library(dplyr)
library(purrr)
library(rlang)
library(magrittr)
library(stringr)
library(disgenet2r)

# ------------------------------------------------------------
# Load local functions
# ------------------------------------------------------------
message("🔧 Loading utility functions (source files)...")
source("preprocessing_v2/src/disgenet_utils/disgenet_category2genes.R")
source("preprocessing_v2/src/disgenet_utils/disgenet_multiDisease2genes.R")

# ##############################################################################
# ---- Load API key ----
# ##############################################################################
message("🔑 Setting DisGeNET API key...")
api_key <- "70b3d3dc-6ddb-4167-a052-1306c238c3da"
Sys.setenv(DISGENET_API_KEY = api_key)

# ##############################################################################
# ---- Read input data ----
# ##############################################################################
message("📥 Loading input data files...")

disgenet_metadataDiseases <- readRDS(
  "data/databases/disgenet/disgenet_metadataDiseases.rds"
)
message("   ✔ Loaded: disgenet_metadataDiseases")

hgnc_symbols_vector_v110 <- readRDS(
  file = "data/databases/biomart/hgnc_symbols_vector_v110.rds"
)
message("   ✔ Loaded: hgnc_symbols_vector_v110")


# # ##############################################################################
# # ---- Mental Disorders (F03) ----
# # ##############################################################################
# message("\n🧠 Running DisGeNET pipeline for category: Mental Disorders (F03)")
# message("──────────────────────────────────────────────")
# 
# output_path <- "data/databases/disgenet/disgenet_mentalDisordersF03.rds"
# 
# disgenet_mentalDisordersF03 <- disgenet_category2genes(
#   msh_categories = "Mental Disorders (F03)",
#   metadata_df    = disgenet_metadataDiseases,
#   filter_genes   = hgnc_symbols_vector_v110,
#   save_to_rds    = output_path
# )
# 
# message("\n🎉 DONE! DisGeNET extraction for Mental Disorders (F03) completed.")
# message("💾 Results saved to: ", output_path)
# message("──────────────────────────────────────────────\n")


# ##############################################################################
# ---- Behavior and Behavior Mechanisms (F01) ----
# ##############################################################################
message("\n🧠 Running DisGeNET pipeline for category: Behavior and Behavior Mechanisms (F01)")
message("──────────────────────────────────────────────")

output_path <- "data/databases/disgenet/disgenet_MSH_BehaviorAndBehaviorMechanismsF01.rds"

disgenet_mentalDisordersF03 <- disgenet_category2genes(
  msh_categories = "Behavior and Behavior Mechanisms (F01)",
  metadata_df    = disgenet_metadataDiseases,
  filter_genes   = hgnc_symbols_vector_v110,
  save_to_rds    = output_path
)

message("\n🎉 DONE! DisGeNET extraction for Behavior and Behavior Mechanisms (F01) completed.")
message("💾 Results saved to: ", output_path)
message("──────────────────────────────────────────────\n")


# ##############################################################################
# ---- Cardiovascular Diseases (C14) ----
# ##############################################################################
message("\n🧠 Running DisGeNET pipeline for category: Cardiovascular Diseases (C14)")
message("──────────────────────────────────────────────")

output_path <- "data/databases/disgenet/disgenet_MSH_cardiovascularDiseasesC14.rds"

disgenet_mentalDisordersF03 <- disgenet_category2genes(
  msh_categories = "Cardiovascular Diseases (C14)",
  metadata_df    = disgenet_metadataDiseases,
  filter_genes   = hgnc_symbols_vector_v110,
  save_to_rds    = output_path
)

message("\n🎉 DONE! DisGeNET extraction for Cardiovascular Diseases (C14) completed.")
message("💾 Results saved to: ", output_path)
message("──────────────────────────────────────────────\n")


# ##############################################################################
# ---- Nutritional and Metabolic Diseases (C18) ----
# ##############################################################################
message("\n🧠 Running DisGeNET pipeline for category: Nutritional and Metabolic Diseases (C18)")
message("──────────────────────────────────────────────")

output_path <- "data/databases/disgenet/disgenet_MSH_NutritionalAndMetabolicDiseasesC18.rds"

disgenet_mentalDisordersF03 <- disgenet_category2genes(
  msh_categories = "Nutritional and Metabolic Diseases (C18)",
  metadata_df    = disgenet_metadataDiseases,
  filter_genes   = hgnc_symbols_vector_v110,
  save_to_rds    = output_path
)

message("\n🎉 DONE! DisGeNET extraction for Nutritional and Metabolic Diseases (C18) completed.")
message("💾 Results saved to: ", output_path)
message("──────────────────────────────────────────────\n")


# ##############################################################################
# ---- Immune System Diseases (C20) ----
# ##############################################################################
message("\n🧠 Running DisGeNET pipeline for category: Immune System Diseases (C20)")
message("──────────────────────────────────────────────")

output_path <- "data/databases/disgenet/disgenet_MSH_ImmuneSystemDiseasesC20.rds"

disgenet_mentalDisordersF03 <- disgenet_category2genes(
  msh_categories = "Immune System Diseases (C20)",
  metadata_df    = disgenet_metadataDiseases,
  filter_genes   = hgnc_symbols_vector_v110,
  save_to_rds    = output_path
)

message("\n🎉 DONE! DisGeNET extraction for Immune System Diseases (C20) completed.")
message("💾 Results saved to: ", output_path)
message("──────────────────────────────────────────────\n")


# ##############################################################################
# ---- Immune System Diseases (C20) ----
# ##############################################################################
message("\n🧠 Running DisGeNET pipeline for category: Immune System Diseases (C20)")
message("──────────────────────────────────────────────")

output_path <- "data/databases/disgenet/disgenet_MSH_BehaviorAndBehaviorMechanisms.rds"

disgenet_BehaviorAndBehaviorMechanismsF01 <- disgenet_category2genes(
  msh_categories = "Behavior and Behavior Mechanisms (F01)",
  metadata_df    = disgenet_metadataDiseases,
  filter_genes   = hgnc_symbols_vector_v110,
  save_to_rds    = output_path
)

message("\n🎉 DONE! DisGeNET extraction for Immune System Diseases (C20) completed.")
message("💾 Results saved to: ", output_path)
message("──────────────────────────────────────────────\n")


# ##############################################################################
# ---- Respiratory Tract Diseases (C08) ----
# ##############################################################################
message("\n🧠 Running DisGeNET pipeline for category: Respiratory Tract Diseases (C08)")
message("──────────────────────────────────────────────")

output_path <- "data/databases/disgenet/disgenet_MSH_RespiratoryTractDisease.rds"

disgenet_RespiratoryTractDisease <- disgenet_category2genes(
  msh_categories = "Respiratory Tract Diseases (C08)",
  metadata_df    = disgenet_metadataDiseases,
  filter_genes   = hgnc_symbols_vector_v110,
  save_to_rds    = output_path
)

message("\n🎉 DONE! DisGeNET extraction for Immune System Diseases (C20) completed.")
message("💾 Results saved to: ", output_path)
message("──────────────────────────────────────────────\n")

# disgenet_metadataDiseases %>%
#   filter(map_lgl(diseaseClasses_MSH, ~ "Psychological Phenomena (F02)" %in% .x))

