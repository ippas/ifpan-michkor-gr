install.packages("remotes")
remotes::install_github("ropensci/gtexr")

library(gtexr)
?get_gene_expression



g <- gtexr::get_gene_search("NR3C1")
gencode_id <- g$gencodeId[1]

med <- gtexr::get_median_gene_expression(gencodeIds = gencode_id)

expr <- gtexr::get_gene_expression(gencodeId = gencode_id)

med %>% 
  mutate(
    tissue = tissue_mapper_vec[tissueSiteDetailId]
  )


tissue_detailed_mapper <- list(
  
  ## BRAIN
  Brain = c(
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra"
  ),
  
  ## Adipose_Tissue
  Adipose_Tissue = c(
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum"
  ),
  
  Blood_Vessel = c(
    "Artery_Aorta", 
    "Artery_Coronary", 
    "Artery_Tibial",
    "Cells_EBV-transformed_lymphocytes"
  ),
  
  Cervix_Uteri = c(
    "Cervix_Ectocervix",
    "Cervix_Endocervix"
  ),
  
  Colon = c(
    "Colon_Sigmoid",
    "Colon_Transverse"
  ),
  
  Esophagus = c(
    "Esophagus_Gastroesophageal_Junction",
    "Esophagus_Mucosa",
    "Esophagus_Muscularis" 
  ),
  
  Heart = c(
    "Heart_Atrial_Appendage",
    "Heart_Left_Ventricle" 
  ),
  
  Kidney = c(
    "Kidney_Cortex",
    "Kidney_Medulla"
  ),
  
  Skin = c(
    "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg"
  ),
  
  Adrenal_Gland        = c("Adrenal_Gland"),
  Bladder              = c("Bladder"),
  Fallopian_Tube       = c("Fallopian_Tube"),
  Liver                = c("Liver"),
  Lung                 = c("Lung"),
  Minor_Salivary_Gland = c("Minor_Salivary_Gland"),
  Muscle_Skeletal      = c("Muscle_Skeletal"),
  Nerve_Tibial         = c("Nerve_Tibial"),
  Ovary                = c("Ovary"),
  Pancreas             = c("Pancreas"),
  Pituitary            = c("Pituitary"),
  Prostate             = c("Prostate"),
  Spleen               = c("Spleen"),
  Stomach              = c("Stomach"),
  Testis               = c("Testis"),
  Thyroid              = c("Thyroid"),
  Uterus               = c("Uterus"),
  Vagina               = c("Vagina"),
  Whole_Blood          = c("Whole_Blood")
)

tissue_mapper_vec <- unlist(
  lapply(names(tissue_detailed_mapper), function(tissue) {
    setNames(
      rep(tissue, length(tissue_detailed_mapper[[tissue]])),
      tissue_detailed_mapper[[tissue]]
    )
  }),
  use.names = TRUE
)
