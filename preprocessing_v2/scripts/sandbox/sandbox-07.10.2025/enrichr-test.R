# ##############################################################################
# ---- uses data ----
# ##############################################################################

AllBrain2BrainSignatures2GlobalMarpiechClusters

nuclear_receptors_go0004879 <- c(
  # Glucocorticoid & Mineralocorticoid
  "NR3C1", "NR3C2",
  # Sex hormone receptors
  "AR", "ESR1", "ESR2", "PGR",
  # Thyroid hormone
  "THRA", "THRB",
  # Retinoic acid and retinoid X
  "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG",
  # Peroxisome proliferator-activated
  "PPARA", "PPARD", "PPARG",
  # Liver X and Farnesoid X
  "NR1H2", "NR1H3", "NR1H4",
  # Vitamin D receptor
  "VDR",
  # Orphan ROR family
  "RORA", "RORB", "RORC",
  # Pregnane X and Constitutive androstane
  "NR1I2", "NR1I3",
  # Hepatocyte nuclear factors
  "HNF4A", "HNF4G",
  # COUP-TFs
  "NR2F1", "NR2F2", "NR2F6",
  # Testicular receptors
  "NR2C1", "NR2C2",
  # NR4A family (orphan)
  "NR4A1", "NR4A2", "NR4A3",
  # DAX/SHP and other orphan receptors
  "NR0B1", "NR0B2",
  # Germ cell nuclear factor
  "NR6A1"
)

# ##############################################################################
# ---- run test enrichr ----
# ##############################################################################
AllBrain2BrainSignatures2Global[!c("metasignature_up", "metasignature_down", "brain_up", "brain_down")]


enrichr_test <-run_enrichr(
  gene_list = AllBrain2BrainSignatures2GlobalMarpiechClusters$`michkor-cells_NA_astrocyte_dexamethasone_NA_NA_NA_in-vitro_NA_up`,
  database = "ChEA_2016"
)

enrichr_test %>% 
  filter(grepl("NR3C1|NR3C2", Term))


# ##############################################################################
# ---- run test enrichr for all gene lists ----
# ##############################################################################
nms <- setdiff(
  names(AllBrain2BrainSignatures2Global),
  c("metasignature_up","metasignature_down","brain_up","brain_down")
)

res_list <- setNames(lapply(nms, function(nm) {
  gl <- AllBrain2BrainSignatures2GlobalMarpiechClusters[[nm]]
  if (is.null(gl) || length(gl) == 0) return(NULL)
  out <- tryCatch(run_enrichr(gene_list = gl, database = "ChEA_2016"), error = function(e) NULL)
  if (is.null(out) || !nrow(out)) return(NULL)
  out <- out[grepl("NR3C1|NR3C2", out$Term, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out$list_name <- nm
  out
}), nms)
  
rbindlist(res_list, fill = TRUE, use.names = TRUE, idcol = "list_name") %>% 
  rownames_to_column(var = "rowname") %>% select(-rowname) %>% 
  as.data.frame() %>% 
  mutate(
    n_overlap = as.integer(str_extract(Overlap, "\\d+(?=/)"))
  ) %>% 
  select(-list_name.1) %>% 
  filter(Term != "NR3C1 23031785 ChIP-Seq PC12 Mouse") %>% filter(n_overlap > 5) %>% 
  filter(grepl("TSC22", Genes))
  

rm(enrichr_test)