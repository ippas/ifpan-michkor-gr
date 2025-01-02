setwd("/home/dzesikah/GEO/")

probki = read.csv('Probki_ILMN.csv', header=F)
probki = as.character(unlist(probki))

GEO_non_normalised = data_non_normalised[match(probki, rownames(data_non_normalised)), -(which(sample_info$descr=="out"))]
GEO_normalised = data_normalised[match(probki, rownames(data_normalised)), -(which(sample_info$descr=="out"))]

write.csv(GEO_normalised, file='GEO_normalised.csv', quote=F)
write.csv(GEO_non_normalised, file='GEO_non_normalised.csv', quote=F)


write.csv(sample_info_out, file = "sample_info.csv", quote = F, row.names = F)
