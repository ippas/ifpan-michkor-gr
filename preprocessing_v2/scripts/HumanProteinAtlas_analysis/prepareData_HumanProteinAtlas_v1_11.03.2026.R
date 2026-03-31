# ##############################################################################
# ---- read data ----
# ##############################################################################
hpa_consensus_organ_mapper <- c(
  
  # brain
  "cerebral cortex" = "brain",
  "cerebellum" = "brain",
  "basal ganglia" = "brain",
  "hypothalamus" = "brain",
  "midbrain" = "brain",
  "amygdala" = "brain",
  "choroid plexus" = "brain",
  "hippocampal formation" = "brain",
  "spinal cord" = "brain",
  
  # eye
  "retina" = "eye",
  
  # endocrine
  "thyroid gland" = "endocrine tissues",
  "parathyroid gland" = "endocrine tissues",
  "adrenal gland" = "endocrine tissues",
  "pituitary gland" = "endocrine tissues",
  
  # respiratory
  "lung" = "respiratory system",
  
  # proximal digestive
  "salivary gland" = "proximal digestive tract",
  "esophagus" = "proximal digestive tract",
  "tongue" = "proximal digestive tract",
  
  # gastrointestinal
  "stomach" = "gastrointestinal tract",
  "duodenum" = "gastrointestinal tract",
  "small intestine" = "gastrointestinal tract",
  "colon" = "gastrointestinal tract",
  "rectum" = "gastrointestinal tract",
  
  # liver
  "liver" = "liver & gallbladder",
  "gallbladder" = "liver & gallbladder",
  
  # pancreas
  "pancreas" = "pancreas",
  
  # kidney
  "kidney" = "kidney & urinary bladder",
  "urinary bladder" = "kidney & urinary bladder",
  
  # male
  "testis" = "male tissues",
  "epididymis" = "male tissues",
  "seminal vesicle" = "male tissues",
  "prostate" = "male tissues",
  
  # female
  "vagina" = "female tissues",
  "ovary" = "female tissues",
  "fallopian tube" = "female tissues",
  "endometrium" = "female tissues",
  "cervix" = "female tissues",
  "placenta" = "female tissues",
  
  # muscle / vascular
  "heart muscle" = "muscle & vascular tissue",
  "smooth muscle" = "muscle & vascular tissue",
  "skeletal muscle" = "muscle & vascular tissue",
  "blood vessel" = "muscle & vascular tissue",
  
  # connective
  "adipose tissue" = "connective & soft tissue",
  "breast" = "connective & soft tissue",
  
  # skin
  "skin" = "skin",
  
  # immune
  "appendix" = "bone marrow & lymphoid tissues",
  "spleen" = "bone marrow & lymphoid tissues",
  "lymph node" = "bone marrow & lymphoid tissues",
  "tonsil" = "bone marrow & lymphoid tissues",
  "bone marrow" = "bone marrow & lymphoid tissues",
  "thymus" = "bone marrow & lymphoid tissues"
)


hpa_single_nuclei_brain_mapper <- c(
  
  # neuronal cells
  "amygdala excitatory" = "neuronal cells",
  "cerebellar inhibitory" = "neuronal cells",
  "c-f neuron" = "neuronal cells",
  "deep-layer corticothalamic and 6b" = "neuronal cells",
  "deep-layer intratelencephalic" = "neuronal cells",
  "deep-layer near-projecting" = "neuronal cells",
  "eccentric medium spiny neuron" = "neuronal cells",
  "hippocampal CA1-3" = "neuronal cells",
  "hippocampal CA4" = "neuronal cells",
  "lamp5-lhx6 and chandelier" = "neuronal cells",
  "lower rhombic lip" = "neuronal cells",
  "mammillary body" = "neuronal cells",
  "medium spiny neuron" = "neuronal cells",
  "mixed neuron" = "neuronal cells",
  "midbrain-derived inhibitory" = "neuronal cells",
  "miscellaneous" = "neuronal cells",
  "obn-satellite" = "neuronal cells",
  "upper-layer intratelencephalic" = "neuronal cells",
  "CGE interneuron" = "neuronal cells",
  "hippocampal dentate gyrus" = "neuronal cells",
  "LAMP5-LHX6 and Chandelier" = "neuronal cells",
  "MGE interneuron" = "neuronal cells",
  "splatter" = "neuronal cells",
  "thalamic excitatory" = "neuronal cells",
  "upper rhombic lip" = "neuronal cells", 
  
  
  # glial cells
  "astrocyte" = "glial cells",
  "Bergmann glia" = "glial cells",
  "committed oligodendrocyte precursor" = "glial cells",
  "oligodendrocyte" = "glial cells",
  "oligodendrocyte precursor cell" = "glial cells",
  "central nervous system macrophage" = "glial cells",
  
  # related cells
  "ependymal cell" = "ciliated cells",
  "choroid plexus epithelial cell" = "ciliated cells",
  
  # endothelial and mural cells
  "endothelial cell" = "endothelial and mural cells",
  "pericyte" = "endothelial and mural cells",
  "vascular associated smooth muscle cell" = "endothelial and mural cells",
  
  # mesenchymal cells
  "fibroblast" = "mesenchymal cells",
  
  # blood and immune cells
  "leukocyte" = "blood and immune cells"
)

# ##############################################################################

HumanProteinAtlas_proteinatlas <- read_tsv("/home/mateusz/projects/ifpan-michkor-gr/data/HumanProteinAtlas_data/proteinatlas.tsv")

HumanProteinAtlas_rnaTissueConsensus <- read_tsv("/home/mateusz/projects/ifpan-michkor-gr/data/HumanProteinAtlas_data/rna_tissue_consensus.tsv") %>% 
  mutate(
    Organ = recode(Tissue, !!!hpa_consensus_organ_mapper)
  ) %>% 
  dplyr::rename(gene_symbol = "Gene name")

HumanProteinAtlas_rnaTissueConsensus %>% 
  # head(100) %>% 
  group_by(Gene, gene_symbol, Organ) %>% 
  nest() %>% 
  mutate(
    organ_nTPM_sum    = map_dbl(data, ~ sum(.x$nTPM, na.rm = TRUE)),
    organ_nTPM_mean   = map_dbl(data, ~ mean(.x$nTPM, na.rm = TRUE)),
    organ_nTPM_median = map_dbl(data, ~ median(.x$nTPM, na.rm = TRUE)),
    organ_nTPM_max    = map_dbl(data, ~ max(.x$nTPM, na.rm = TRUE))
  ) %>% unnest %>% ungroup -> HumanProteinAtlas_rnaTissueConsensus

# write_tsv(
#   HumanProteinAtlas_rnaTissueConsensus,
#   "/home/mateusz/projects/ifpan-michkor-gr/data/HumanProteinAtlas_data/rna_tissue_consensus_withOrganSummary_11.03.2026.tsv"
# )
# 
# # ##############################################################################
# HumanProteinAtlas_rnaTissueConsensus <- read_tsv(
#   "/home/mateusz/projects/ifpan-michkor-gr/data/HumanProteinAtlas_data/rna_tissue_consensus_withOrganSummary_11.03.2026.tsv"
# )

HumanProteinAtlas_rnaTissueConsensus %>% 
  filter(Organ == "brain")


# ##############################################################################
HumanProteinAtlas_13rnaBrainRegion <- read_tsv("/home/mateusz/projects/ifpan-michkor-gr/data/HumanProteinAtlas_data/rna_brain_region_hpa.tsv") %>% 
  dplyr::rename(gene_symbol = "Gene name") %>% 
  dplyr::rename(ensembl_id = "Gene") %>% 
  dplyr::rename(brain_region = "Brain region")

HumanProteinAtlas_13rnaBrainRegion %>% 
  filter(gene_symbol == "FKBP5")
  
  
# ##############################################################################
HumanProteinAtlas_34rnaSingleNucleiClusterType <- read_tsv("/home/mateusz/projects/ifpan-michkor-gr/data/HumanProteinAtlas_data/rna_single_nuclei_cluster_type.tsv") %>% 
  dplyr::rename(gene_symbol = "Gene name") %>% 
  dplyr::rename(ensembl_id = "Gene") %>% 
  dplyr::rename(cluster_type = "Cluster type") %>% 
  mutate(
    cellType_group = recode(cluster_type, !!!hpa_single_nuclei_brain_mapper)
  ) 

HumanProteinAtlas_34rnaSingleNucleiClusterType %>% 
  group_by(ensembl_id, gene_symbol, cluster_type) %>% 
  nest() %>% 
  mutate(
    cellType_nCPM_sum    = map_dbl(data, ~ sum(.x$nCPM, na.rm = TRUE)),
    cellType_nCPM_mean   = map_dbl(data, ~ mean(.x$nCPM, na.rm = TRUE)),
    cellType_nCPM_median = map_dbl(data, ~ median(.x$nCPM, na.rm = TRUE)),
    cellType_nCPM_max    = map_dbl(data, ~ max(.x$nCPM, na.rm = TRUE))
  ) %>% unnest %>% ungroup -> HumanProteinAtlas_34rnaSingleNucleiClusterType
  
  
