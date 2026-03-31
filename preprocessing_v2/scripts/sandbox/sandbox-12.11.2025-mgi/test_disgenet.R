library(devtools)
install_gitlab("medbio/disgenet2r")

library(disgenet2r)


api_key <- "70b3d3dc-6ddb-4167-a052-1306c238c3da"
Sys.setenv(DISGENET_API_KEY= api_key)

results <- gene2disease( gene = 3953, vocabulary = "ENTREZ",
                         database = "CURATED")

results@qresult


results <- gene2disease(
  gene = "BCL2L1",
  vocabulary = "HGNC",
  database = "CURATED"
)

mental_disorders <- diseaseClass2disease(
  diseaseClass = "Mental Disorders (F03)",
  database = "CURATED"  # dostępne w profilu akademickim
)


flat_allGrSignatures_31.10.2025$global_GR_genes_globalUp5TissuesDerivedCells


results <- gene2disease(
  gene = flat_allGrSignatures_31.10.2025$global_GR_genes_globalUp6TissuesDerivedCells,
  vocabulary = "HGNC",
  database = "ALL"
)

myListOfGenes <- c("KCNE1", "KCNE2", "KCNH1", "KCNH2", "KCNG1")




hgnc_symbols_vector_v110

results <- gene2disease(
  gene      = myListOfGenes,
  database  = "CURATED",
  vocabulary = "HGNC",
  verbose   = TRUE
)

results@qresult %>% 
  dplyr::select(!c("gene_symbol", 
           "geneid", 
           "ensemblid", 
           "geneNcbiType", 
           "geneDSI", 
           "geneDPI",
           "genepLI",
           "uniprotids",
           "protein_classid",
           "protein_class_name",
           "evidence_index", 
           "yearInitial", 
           "yearFinal", 
           "score",
           "numberPmidsWithChemsFiltered",
           "numNCTSWithChemsIncludedInEvidences",
           "evidence_level",
           "chemsIncludedInEvidenceBySource", 
           "numChemsIncludedInEvidences", 
           "numPMIDSWithChemsIncludedInEvidences",
           "numPMIDs",
           "disease_prevalence_class",
           "disease_prevalence_geo_area",
           "disease_prevalence_type",
           "disease_inheritance",
           "numCTsupportingAssociation",
           "numPMIDs",
           "diseaseType"
           ))

