# load data


factors_rsidGenes_P1e4tss100kb <- read.table(file = "data/factors-pmid38965376/factorID-pmid38965376-rsidP1e4-proteinCodingv110-tss100kb.tsv", 
           sep = "\t", 
           header = TRUE) 

factors_rsidGenes_P1e4tss100kb %>% head
           
lite_grSignatures


metasignatureUp_GWASCatalog <- c(
  "CDKN1A", "CRISPLD2", "IL1R2", "FAM107A", "CXCL5",
  "DCN", "ART3", "IL17RA", "PNMT", "ANGPTL4",
  "CD163", "HP", "MAP3K6", "ZBTB16"
)

metasignatureDown_GWASCatalog <- c(
  "IL11", "CHST3", "MME", "IFITM2", "ZBTB16", "KLRC1",
  "CELA2A", "TRPV6", "TOP2A", "CCR7", "FOSL1", "CSF2",
  "BCL3", "TNFRSF11B", "KLRG1", "GDF15", "IL12B"
)

brainUp_GWASCatalog <- c(
  "F13A1", "CYP1B1", "IL1R1", "TNFRSF11A", "S100A11",
  "ERRFI1", "ZBTB16", "HTRA1", "TNS1", "FKBP5",
  "SPSB1", "ANKDD1B", "RAMP2", "IGF2", "PIM3",
  "RHOU", "PER1", "PABPC1L", "FAM107A", "LMOD1",
  "OMG"
)

brainDown_GWASCatalog <- c(
  "DSP", "CSF1", "LPL", "CES2", "RMI2", "TNXB", "EHD4", "IL1R1", "TNFRSF11B",
  "ACAD10", "CTSC", "FAM171A1", "SLCO4A1", "PIEZO2", "FMO2", "LRRC8E", "DENND2C",
  "LAPTM5", "PTPN7", "LAMA1", "RAC2", "PLAU", "CXCL10", "HSD17B14", "CCL18",
  "CORO1A", "ADAP2", "ZNF385C", "WNT7A", "HPRT1", "FCGR3A", "TDRD12"
)

################################################################################
# preprocessing data
################################################################################
# grSignatureLite_association_factorsP1e4Tss100kb

grSignatureLite_association_factorsP1e4Tss100kb <- factors_rsidGenes_P1e4tss100kb %>% 
  inner_join(lite_grSignatures, by = "gene_symbol") %>%
  dplyr::select(signature_name, everything()) %>% 
  group_by(signature_name, factor_id, gene_symbol) %>%
  slice_min(pvalue, with_ties = FALSE) %>% 
  ungroup %>% 
  mutate(
    in_GWASCatalog = case_when(
      signature_name == "metasignature_up"   & gene_symbol %in% metasignatureUp_GWASCatalog   ~ 1,
      signature_name == "metasignature_down" & gene_symbol %in% metasignatureDown_GWASCatalog ~ 1,
      signature_name == "brain_up"           & gene_symbol %in% brainUp_GWASCatalog           ~ 1,
      signature_name == "brain_down"         & gene_symbol %in% brainDown_GWASCatalog         ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  select(signature_name, gene_symbol, in_GWASCatalog, everything()) %>%
  arrange(pvalue) %>%   group_split(signature_name)

names(grSignatureLite_association_factorsP1e4Tss100kb) <- c("brain_down", "brain_up", "metasignature_down", "metasignature_up")

grSignatureLite_association_factorsP1e4Tss100kb

wb <- createWorkbook()

for (sheetname in c("metasignature_up", "metasignature_down", "brain_up", "brain_down")) {
  df <- grSignatureLite_association_factorsP1e4Tss100kb[[sheetname]]
  
  addWorksheet(wb, sheetName = sheetname)
  writeData(wb, sheet = sheetname, x = df, withFilter = TRUE)
  freezePane(wb, sheet = sheetname, firstRow = TRUE)
}

saveWorkbook(
  wb,
  file = "data/factors-pmid38965376/grSignatureLite_association_factorsP1e4Tss100kb_topRsidByGene.xlsx",
  overwrite = TRUE
)
