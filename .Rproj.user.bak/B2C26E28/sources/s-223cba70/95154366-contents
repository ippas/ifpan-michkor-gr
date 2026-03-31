acute_stress <- c(
  "Cldn5", "Sox2", "Fabp7", "Col1a2", "Sema3b",
  "Homer1", "Sult1a1", "Nr4a3", "Tnfrsf25",
  "Rasl11b", "Mertk", "Midn", "Zfp189", "Nr4a1",
  "Fam107a", "Tsc22d3", "Il1Rap", "Plin4",
  "Tenm4", "Ddit4", "Arl4d", "Npas4", "Adipor2", "Arc",
  "Junb", "Gjb6", "Egr1", "Mfsd2a",
  "Arrdc2", "Gadd45g", "Dusp1", "Gadd45b", "Dio2",
  "Ier2", "Slc2a1", "Nfkbia", "Dusp6", "Errfi1", "Btg2",
  "Per1", "Fosb", "Egr2", "Fos", "Sgk1", "Cdkn1a", "Ccn1",
  "Htra1", "Csrnp1", "Tiparp", "Cacna2d1", "Irs2",
  "Kcna1", "Mt1", "Phf13", "Pim3", "Sdc4", "Tob2"
)


prolonged_stress <- c(
  "Adamts4", "Gjb1", "Fcrls", "Emp2", "Car2", "Ttyh2",
  "Adcy5", "Ido1", "Plxnb3", "Pllp", "Foxn3",
  "Tmem125", "Atp2b2", "Bcas1", "Cachd1", "Cryab",
  "Plp1", "Ppp1r14a", "Rgs9", "Bfsp2", "Gjc2", "C1qb",
  "Aldh1a3", "Cwc22", "Fabp7", "Mbp", "Ptn", "Vamp1",
  "Map1b", "Ttr", "Mog", "Igf2", "Hba-a2", "Tcf7l2",
  "Il6ra", "Sgk1", "Ptgs2", "Bmp4", "Fmo2", "Per2",
  "Tsc22d3", "Hbb-bs", "Ccn2", "Prg4", "S100a8",
  "Hba-a1", "Pten", "Sult1a1", "Hbb-bt", "Alas2",
  "Sparc", "Galnt15", "S100a9", "Fkbp5", "Mgp",
  "Mpzl2", "Cyp1b1", "Lcn2", "Lrg1"
)

vulnerable <- c(
  "Fgf18", "Ocln", "4930556M19Rik",
  "A330015K06Rik", "Adcy5", "Atp2a3", "Bche",
  "Bfsp2", "Casc1", "Fabp5", "Glrp1", "Id1", "Klf12",
  "Lama2", "Nkx6–2", "Npas3", "Pde4c", "Plekhg3",
  "Rps4l-ps", "Tgds", "Galnt6", "Dcdc2c", "Klf10",
  "Tom1", "Cemip", "Slc25a42", "Ctla2a", "Ifi30",
  "Krt1", "Nabp1", "Slc4a1", "Eif3a", "Fam234a",
  "Klra2", "Ly6c1", "Ncmap", "Pmaip1", "Prl",
  "Usp27x", "Zfyve26", "Pparg", "S100a10",
  "Clcn5", "Fdxr", "Nlrp10", "Trpm1", "Abcc12",
  "4921536K21Rik", "Paqr5", "Slc25a47",
  "Wrap73"
)

resistant <- c(
  "Cacng7", "Ntrk2", "Prepl", "Luzp2", "Syt1", "Ccr5",
  "Col6A1", "Stx1A", "Fancd2", "Gsx2", "Kcnk12",
  "Mmp23", "Pih1D1", "Tdrd12", "Zbtb18", "Pparg",
  "S100a10", "Clcn5", "Fdxr", "Nlrp10", "Trpm1",
  "Abcc12", "4921536K21Rik", "Paqr5",
  "Slc25a47", "Wrap73"
)
  
ptsd1 <- c(
  "LSR", "NKX2–2", "NOTCH4", "PCBP4",
  "SEMA3B", "SLC25A48", "STARD3",
  "ZCCHC24", "ARHGEF10", "CAPN3",
  "CDC42SE1", "CDK5RAP2", "CLMN", "CNKSR3",
  "DNAH1", "DTX2", "DYSF", "ERMN", "FA2H",
  "FBXO9", "GAB2", "GAB3", "GAL3ST1", "GMEB2",
  "LPAR1", "MTUS1", "OPALIN", "PLLP", "POU3F4",
  "ST18", "TTYH2", "ZBED3", "ZNF77",
  "HSP90AB1", "APP", "CACNA2D3", "CCT6A",
  "MRPL55", "SLC9A6", "SSTR1"
)

ptsd2 <- c(
  "A2M", "AQP1", "CGNL1", "EGR2", "LAMA4",
  "MBP", "MOBP", "NR4A3", "OPALIN", "PENK",
  "PER2", "PLLP", "SEMA3B", "SGK1", "TSC22D4",
  "WNK1", "ZEB2", "ADAMTS1", "ADGRL3",
  "ANK2", "APP", "BHLHE22", "CADPS2",
  "CARTPT", "CCBE1", "CDKN1A", "COL1A1",
  "COL1A2", "EGR4", "FMOD", "FOSB", "FUBP1",
  "GADD45B", "GNAS", "GPR88", "HBA1", "HBA2",
  "HSPA1A", "IGF2", "IGFBP2", "MAP1B", "MEF2C",
  "NFIL3", "NPAS4", "PRKCD", "PTGS2", "PTPRK",
  "RASD1", "RXFP1", "SLC13A4", "SLC17A7",
  "SOSTDC1", "THBD", "TSHZ2", "VAT1L"
)


genes_juszczak_df <- data.frame(
  category = rep(c("acute_stress", "prolonged_stress", "vulnerable", "resistant", "ptsd1", "ptsd2"),
                 times = c(length(acute_stress), length(prolonged_stress), length(vulnerable), 
                           length(resistant), length(ptsd1), length(ptsd2))),
  gene_symbol = c(acute_stress, prolonged_stress, vulnerable, resistant, ptsd1, ptsd2)
)

install.packages("gprofiler2")

library(gprofiler2)

result <- gorth(query = genes_juszczak_df$gene_symbol,
                source_organism = "mmusculus",
                target_organism = "hsapiens")

mapping_df <- result[, c("input", "ortholog_name")]
colnames(mapping_df) <- c("gene_name", "hgnc_symbol")

print(mapping_df)

# removing 21 genes,
genes_juszczak_df %>%
  left_join(., mapping_df, by = c("gene_symbol" = "gene_name")) %>% 
  filter(!is.na(hgnc_symbol)) %>% 
  mutate(category = paste0(category, "_juszczak")) -> genes_juszczak_df 




  



