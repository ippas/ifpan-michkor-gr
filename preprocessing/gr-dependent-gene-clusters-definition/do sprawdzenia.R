#~~~~~~~~~~ sample_info
setwd('/home/dzesikah/project/dextis')

sample_info = read.csv('sample_nfo_dextis_fixed_sort.csv', header=T, colClasses='character')

names2 = vector(length=length(colnames(data)))
for (i in 1:length(names2)){
  names2[i] = paste(sample_info$Tissue_code[i], sample_info$Treat[i], sep='_')
}

rm(i)
sample_info = cbind(sample_info, names2)

rm(names2)

sample_info = sample_info[order(sample_info$names2), ]

#~~~~~~~~~~~~~~~~~~ read data


data = read.csv('matrix_non_normalized.csv', header=T)

row_names = data[, 1]
data = data[, 2:97]
rownames(data) = row_names
colnames(data) = substr(colnames(data), 2, 30)
rm(row_names)

sample_id = match(sample_info$sample, colnames(data))
data = data[,sample_id]

data_normalised = read.csv('matrix_normalized.csv', header=T)

row_names = data_normalised[, 1]
data_normalised = data_normalised[, 2:97]
rownames(data_normalised) = row_names
colnames(data_normalised) = substr(colnames(data_normalised), 2, 30)
rm(row_names)

sample_id = match(sample_info$sample, colnames(data_normalised))
data_normalised = data_normalised[, sample_id]
rm(sample_id)

#~~~~~~~~~~ annotacje
setwd('/home/dzesikah/Annotacje')
annotacje = read.csv('annotacje.csv', header=T, sep='\t', colClasses='character')

annotacje = annotacje[ , c('Entrez_Gene_ID', 'Symbol', 'Probe_Id', 'Probe_Sequence', 'Chromosome', 'Probe_Coordinates')]


list_probeId = (match(rownames(data_normalised), annotacje[,"Probe_Id"]))
annotacje = annotacje[list_probeId, ]

rm(list_probeId)


# uzupelnienie bazy danych

x <- illuminaMousev2SYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx = as.list(x[rownames(data_normalised)])


for (i in 1:length(rownames(data_normalised))){
  if (sum(is.na(xx[[i]])) == 0)
    annotacje$Symbol[i] = xx[[i]]
}

rm(i)
rm(x)
rm(xx)

#~~~~~~~~~~~~~~~~~~~~~~~~ data -> matrix

matrix_normalised = as.matrix(data_normalised)

colnames(matrix_normalised) = sample_info$names2
rownames(matrix_normalised) = annotacje$Symbol


#~~~~~~~~~~~~~~~~~~ standarize

install.packages("http://darwin.cremag.org/resources/rscripts_0.1.tar.gz", repo=NULL, type="source")
require('rscripts')

matrix_normalised_st = apply(matrix_normalised, 1, function(x) standardize(x, factor = list(sample_info$Tissue_code),
                                                                                           control.sample = (sample_info$Treat == 'Sal')))

matrix_normalised_st = t(matrix_normalised_st)
#~~~~~~~~~~~~~~~~ two way anova

cells = as.factor(sample_info$Tissue_code)
treat = as.factor(sample_info$Treat)

anova = apply(matrix_normalised, 1, function(x) anova(aov(x ~ cells*treat))[1:3, 5])
anova = t(anova)
colnames(anova) = c('cells', 'treat', 'cells:treat')

rm(cells)
rm(treat)

#~~~~~~~~~~~~~~~~~~~~~~~ FDR treat

pvalue_fdr_treat = p.adjust(anova[, 2], method='fdr')

#~~~~~~~~~~~~~~~~~~~~~~~~ FDR int

pvalue_fdr_int = p.adjust(anova[, 3], method='fdr')

#~~~~~~~~~~~~~~~~~~~~~~~~ bonferroni treat

pvalue_bonferroni_treat = p.adjust(anova[, 2], method='bonferroni')

#~~~~~~~~~~~~~~~~~~~~~~~~ bonferroni int

pvalue_bonferroni_int = p.adjust(anova[, 3], method='bonferroni')


#~~~~~~~~~~~~~~~~~~~~~~~~ fold 

fold = matrix(nrow=dim(matrix_normalised)[1], ncol=length(levels(as.factor(sample_info$Tissue_code)))*2)
colnames(fold) = c(paste(levels(as.factor(sample_info$Tissue_code)), "fold-mean", sep='_'),
                   paste(levels(as.factor(sample_info$Tissue_code)), "fold-median", sep='_'))
rownames(fold) = rownames(matrix_normalised)

counter = 1
for (tissue in levels(as.factor(sample_info$Tissue_code))){
  fold[, counter] =apply(matrix_normalised, 1, function(x) {mean(x[sample_info$Tissue_code==tissue & sample_info$Treat=='Dex'])-
                                                            mean(x[sample_info$Tissue_code==tissue & sample_info$Treat=='Sal'])}) 
  fold[, counter+9] =apply(matrix_normalised, 1, function(x) {median(x[sample_info$Tissue_code==tissue & sample_info$Treat=='Dex'])-
                                                              median(x[sample_info$Tissue_code==tissue & sample_info$Treat=='Sal'])}) 
  counter = counter+1  
}

rm(counter)
rm(tissue)

fold_mean = fold[ , 1:9]

max_fold = apply(fold_mean, 1, function(x) mean(abs(x)))

fold_mean = cbind(fold_mean, max_fold)

#~~~~~~~~~~~~~~~~~~~~~~~~


percent=c(0.05, 0.01, 0.001)
n_pvalue_maxfold = matrix(ncol =4, nrow=length(percent))

n=0.5
counter=1
for (i in percent){
  n_pvalue_maxfold[counter, ] = c(sum(pvalue_fdr_treat<i & max_fold > n), sum(pvalue_fdr_int<i & max_fold > n), 
                            sum(pvalue_bonferroni_treat<i & max_fold > n), sum(pvalue_bonferroni_int<i & max_fold > n))
  counter = counter + 1
}
colnames(n_pvalue_maxfold)=paste(c('fdr_treat, fold >', 'fdr_int, fold >', 'bonferroni_treat, fold >', 'bonferroni_int, fold >'), n, sep=' ')
rownames(n_pvalue_maxfold)=as.character(percent)

rm(counter)
rm(percent)
rm(i)


#~~~~~~~~~~~~~~~~~~~~~~~~ heaemap bonferoni

heatmap.2(matrix_normalised[which(pvalue_bonferroni_treat<0.05)[1:50], ], distfun = dist.pear, col = bluered(1000), trace = "none", Colv=NA, scale="row")
heatmap.2(matrix_normalised_st[union(which(pvalue_bonferroni_treat<0.05), which(pvalue_bonferroni_int<0.05)), ], distfun = dist.pear, col = bluered(1000), 
          trace = "none", Colv=NA, scale="row")

heatmap.2(matrix_normalised_st[union(which(pvalue_fdr_treat<i & max_fold > n), which(pvalue_fdr_int<i & max_fold > n)), ], 
          distfun = dist.pear, col = bluered(1000), trace = "none", Colv=NA, scale="row")

#~~~~~~~~~~~~~~~~~~~~~~~~~ usuniecie wartosci odstajacych z probki mus_sal

matrix_normalised_out = matrix_normalised[, -which(sample_info$descr == "out")]
sample_info_out = sample_info[-which(sample_info$descr == 'out'), ]

#~~~~~~~~~~~~~~~~ two way anova

cells = as.factor(sample_info_out$Tissue_code)
treat = as.factor(sample_info_out$Treat)

anova = apply(matrix_normalised_out, 1, function(x) anova(aov(x ~ cells*treat))[1:3, 5])
anova = t(anova)
colnames(anova) = c('cells', 'treat', 'cells:treat')

rm(cells)
rm(treat)

#~~~~~~~~~~~~~~~~~~~~~~~ FDR treat

pvalue_fdr_treat = p.adjust(anova[, 2], method='fdr')

#~~~~~~~~~~~~~~~~~~~~~~~~ FDR int

pvalue_fdr_int = p.adjust(anova[, 3], method='fdr')

#~~~~~~~~~~~~~~~~~~~~~~~~ bonferroni treat

pvalue_bonferroni_treat = p.adjust(anova[, 2], method='bonferroni')

#~~~~~~~~~~~~~~~~~~~~~~~~ bonferroni int

pvalue_bonferroni_int = p.adjust(anova[, 3], method='bonferroni')


#~~~~~~~~~~~~~~~~~~~~~~~~ fold 

fold = matrix(nrow=dim(matrix_normalised_out)[1], ncol=length(levels(as.factor(sample_info_out$Tissue_code))))
colnames(fold) = paste(levels(as.factor(sample_info_out$Tissue_code)), "fold-mean", sep='_')
rownames(fold) = rownames(matrix_normalised_out)

counter = 1
for (tissue in levels(as.factor(sample_info_out$Tissue_code))){
  fold[, counter] =apply(matrix_normalised_out, 1, function(x) {mean(x[sample_info_out$Tissue_code==tissue & sample_info_out$Treat=='Dex'])-
                                                              mean(x[sample_info_out$Tissue_code==tissue & sample_info_out$Treat=='Sal'])}) 

  counter = counter+1  
}

rm(counter)
rm(tissue)


max_fold = apply(fold, 1, function(x) mean(abs(x)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

percent=c(0.05, 0.01, 0.001)
n_pvalue_maxfold = matrix(ncol =4, nrow=length(percent))

n=0.5
counter=1
for (i in percent){
  n_pvalue_maxfold[counter, ] = c(sum(pvalue_fdr_treat<i & max_fold > n), sum(pvalue_fdr_int<i & max_fold > n), 
                                  sum(pvalue_bonferroni_treat<i & max_fold > n), sum(pvalue_bonferroni_int<i & max_fold > n))
  counter = counter + 1
}
colnames(n_pvalue_maxfold)=paste(c('fdr_treat, fold >', 'fdr_int, fold >', 'bonferroni_treat, fold >', 'bonferroni_int, fold >'), n, sep=' ')
rownames(n_pvalue_maxfold)=as.character(percent)

rm(counter)
rm(percent)
rm(i)

#~~~~~~~~~~~~~~~~~~ standarize

install.packages("http://darwin.cremag.org/resources/rscripts_0.1.tar.gz", repo=NULL, type="source")
require('rscripts')

matrix_normalised_st = apply(matrix_normalised_out, 1, function(x) standardize(x, factor = list(sample_info_out$Tissue_code),
                                                                           control.sample = (sample_info_out$Treat == 'Sal')))

matrix_normalised_st = t(matrix_normalised_st)
#~~~~~~~~~~~~~~~~~~~~ heatmap

dist.pear <- function(x) as.dist(1-cor(t(x)))

n=0.5
i=0.05

heatmap.2(matrix_normalised_st[union(which(pvalue_fdr_treat<i & max_fold > n), which(pvalue_fdr_int<i & max_fold > n)), ], 
          distfun = dist.pear, col = bluered(1000), trace = "none", Colv=NA, scale="row")

#~~~~~~~~~~~~~~~~~~~~~~

dist = dist.pear(matrix_normalised_st)
hclust = hclust(dist.pear, method = 'complete')

