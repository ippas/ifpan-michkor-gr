setwd('/home/dzesikah/project/dextis/')

#~~~~~~~~~~~~~~~~~~ read data
 
data_non_normalised = read.csv("matrix_non_normalized.csv", header=T)

rownames(data_non_normalised) = data_non_normalised[, 1]
data_non_normalised = data_non_normalised[, 2:97] 
colnames(data_non_normalised) = substr(colnames(data_non_normalised), 2, 20)


data_normalised = read.csv('matrix_normalized.csv', header=T)

rownames(data_normalised) = data_normalised[, 1]
data_normalised = data_normalised[, 2:97] 
colnames(data_normalised) = substr(colnames(data_normalised), 2, 20)


#~~~~~~~~~~~~~~~~~~ sample info

sample_info = read.csv('sample_nfo_dextis_fixed_sort.csv', colClasses = 'character')

names2 = vector(length=dim(sample_info)[1])
for (i in 1:length(names2)){
  names2[i] = paste(sample_info$Tissue_code[i], sample_info$Treat[i], sep='_')
}

rm(i)

sample_info = cbind(sample_info, names2)

rm(names2)


lista_sample = match(colnames(data_non_normalised), sample_info$sample)
sample_info = sample_info[lista_sample, ]

rm(lista_sample)

#~~~~~~~~~~~~~~~~ annotacje
setwd('/home/dzesikah/Annotacje/')
annotacje = read.csv('annotacje.csv', header=T, sep='\t', colClasses='character')

annotacje = annotacje[ , c('Entrez_Gene_ID', 
                           'Symbol', 'Probe_Id', 
                           'Probe_Sequence', 'Chromosome', 
                           'Probe_Coordinates')]

lista_probeId = match(rownames(data_normalised), annotacje[,"Probe_Id"])
annotacje = annotacje[lista_probeId, ]

rm(lista_probeId)


# uzupelnienie bazy danych

x <- illuminaMousev2SYMBOL
# Get the probe identifiers that are mapped to a gene symbol

# Convert to a list
xx = as.list(x[rownames(data_normalised)])


for (i in 1:length(rownames(data_normalised))){
  if (sum(is.na(xx[[i]])) == 0)
    annotacje$Symbol[i] = xx[[i]]
}

rm(i)
rm(x)
rm(xx)


#~~~~~~~~~~~~~~ out out

matrix_normalised = as.matrix(data_normalised)
colnames(matrix_normalised) = sample_info$names2
rownames(matrix_normalised) = annotacje$Symbol

matrix_normalised_out = matrix_normalised[, -c(which(sample_info$descr == "out"))]


sample_info_out = sample_info[-c(which(sample_info$descr == "out")), ]



