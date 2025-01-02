#~~~~~~~~~~~~~~~~~~~~` standaryzacja przez srednia

install.packages("http://darwin.cremag.org/resources/rscripts_0.1.tar.gz", repo=NULL, type="source")
require('rscripts')

matrix_normalised_st_mean = apply(
  matrix_normalised_out, 1, function(x) standardize(x, factor = list(sample_info_out$Tissue_code),
                                                    control.sample = (sample_info_out$Treat == 'Sal')))

matrix_normalised_st_mean = t(matrix_normalised_st_mean)


#~~~~~~~~~~ standaryzacja na mediane

matrix_normalised_st_median = apply(
  matrix_normalised_out, 1, function(x) standardize_median(x, factor = list(sample_info_out$Tissue_code),
                                                           control.sample = (sample_info_out$Treat == 'Sal')))

matrix_normalised_st_median = t(matrix_normalised_st_median)


#~~~~~~~~~~~~~~~ anova

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


#~~~~~~~~~~~~~~ fold

fold = matrix(ncol = 9, nrow = 45281)
sal = substr(colnames(matrix_normalised_out), 5, 8) == "Sal" 
dex = substr(colnames(matrix_normalised_out), 5, 8) == "Dex"
names = levels(as.factor(substr(levels(factor(colnames(matrix_normalised_out))), 1, 3)))
for (i in 1:9){
  tissue = substr(colnames(matrix_normalised_out), 1, 3) == names[i] 
  fold[,i] = apply(matrix_normalised_out, 1, 
                       function(x) {mean(x[dex&tissue])-mean(x[sal&tissue])})
}

colnames(fold)=names

rm(i)
rm(tissue)
rm(names)
rm(sal)
rm(dex)

max_fold = apply(fold, 1, function(x) max(abs(x)))
max_fold = t(max_fold)


# Pierwsza analiza
# 1. fdr interakacja 10%
# 2. foldmax > 0.5
# 3. 30 klastrow

# Druga analiza
# 1. fdr treat 10% bez interakcji
# 2. foldmax >0.05
# 3. 2 klasty



require(magrittr)
cut.above.thresh <- function(x, thresh) {
  x[x > thresh] <- thresh; x[x < -thresh] <- -thresh; x
}

#~~~~~~~~~~~~~~~ odleglosc - korelacja

dist.pear = function(x){ 
  as.dist(1 - cor(t(x)))
}


#~~~~~~~~~~~~~~

drow_heatmap = function(row, col, k, dist, matrix){
  matrix_heatmap = matrix[row, col]
  matrix_heatmap %>% apply(1, scale) %>% t() %>% apply(1, cut.above.thresh, 4) %>% t() -> matrix_heatmap_out
  colnames(matrix_heatmap_out) = colnames(matrix_heatmap)
  #matrix_heatmap_out = matrix_heatmap_out[, 45:93]
  dist_matrix = dist(matrix_heatmap_out)
  hclust_matrix = hclust(dist_matrix, method = "complete")
  dev.off()
  heatmap.2(
    (matrix_heatmap_out[hclust_matrix$order, ]),
    #keysize=0.3,    
    distfun = dist, 
    col = bluered(500), 
    scale = 'none', 
    trace = "none", 
    Rowv=NA, 
    Colv = NA ,
    RowSideColors = as.character(cutree(hclust_matrix, k))[hclust_matrix$order],
    margins = c(5, 5))
}
drow_heatmap(row, col, 2, dist.pear, matrix_normalised_st_median)

fold.thr = 1
fdr.thr = 0.01

col = sal_dex_n
row = which(pvalue_fdr_treat< fdr.thr & pvalue_fdr_int> fdr.thr & max_fold > fold.thr)
row = which(pvalue_fdr_int< fdr.thr & max_fold > fold.thr)

which(pvalue_fdr_treat< fdr.thr & pvalue_fdr_int> fdr.thr & max_fold > fold.thr) %>% length
which(pvalue_fdr_int< fdr.thr & max_fold > fold.thr) %>% length

#~~~~~~~~~ posportowanie wierszy 

sal= paste(c('HTH', 'LIV', 'UNG', 'ADR', 'SPL', 'MUS', 'FAT', 'PIT', 'KID'), '_Sal', sep='')
dex = paste(c('HTH', 'LIV', 'UNG', 'ADR', 'SPL', 'MUS', 'FAT', 'PIT', 'KID'), '_Dex', sep='')
sal_dex = c(sal, dex)

factor_sal_dex = factor(colnames(matrix_normalised_out), levels = sal_dex)

sal_dex_n = order(factor_sal_dex, sort(factor_sal_dex))

rm(sal)
rm(dex)
rm(sal_dex) 
rm(factor_sal_dex)


#~~~~~~~~~~~~ wyciągnicie nazw genów

row = which(pvalue_fdr_int< fdr.thr & max_fold > fold.thr)

matrix_tree = matrix_normalised_st_median[row, ]
matrix_tree %>% apply(1, scale) %>% t() %>% apply(1, cut.above.thresh, 4) %>% t() -> matrix_tree_out
colnames(matrix_tree_out) = colnames(matrix_tree)
dist_matrix = dist.pear(matrix_tree_out)
hclust_matrix = hclust(dist_matrix, method = "complete")

trees = cutree(hclust_matrix, 2)

genes=list()
for (i in 3:22){
  genes[[i]] = names(trees[trees==i-2])
}

#~~~~~~~~~~~ upreg downreg

row = which(pvalue_fdr_treat< fdr.thr & pvalue_fdr_int> fdr.thr & max_fold > fold.thr)

matrix_tree = matrix_normalised_st_median[row, ]
matrix_tree %>% apply(1, scale) %>% t() %>% apply(1, cut.above.thresh, 4) %>% t() -> matrix_tree_out
colnames(matrix_tree_out) = colnames(matrix_tree)
#matrix_heatmap_out = matrix_heatmap_out[, 45:93]
dist_matrix = dist.pear(matrix_tree_out)
hclust_matrix = hclust(dist_matrix, method = "complete")


trees = cutree(hclust_matrix, 2)

for (i in 1:2){
  genes[[i]] = names(trees[trees==i])
}

#~~~~~~~~~~~~~

names_cl = c("Down", "Up", 1:20)

matrix_cl = matrix(ncol = 22, nrow = max(unlist(lapply(genes, length))) )
colnames(matrix_cl) = names_cl

for (k in 1:22){
  for (w in 1:133){
    print(genes[[k]][w], quote=FALSE) -> matrix_cl[w, k]
  }
}

rm(trees)
rm(matrix_tree)
rm(matrix_tree_out)
rm(names_cl)
rm(hclust_matrix)
rm(k)
rm(w) 



#~~~~~~~~~~~~~~

setwd('/home/dzesikah/html/dextis')
write.csv(matrix_cl, 'klastry.csv')

#~~~~~~~~~~~~~~~~~ barplot

order_fold = apply(fold, 2, function(x) order(x, decreasing=T))
mean_barplot = vector(length=9)

mean_barplot = apply(fold, 2, function(x) mean(x[order(x, decreasing=T)[1:20]]))

barplot(mean_barplot)

#~~~~~~ podpisanie klastrów

clast = as.integer(trees[hclust_matrix$order][sort(match(1:20, trees[hclust_matrix$order]))])

#~~~~~~~ klastry tkanek

tissue = c('HTH', 'LIV', 'UNG', 'ADR', 'SPL', 'MUS', 'FAT', 'PIT', 'KID')

tissue_sd = apply(matrix_normalised_out, 1, sd)


tissue_sd_1000 = tissue_sd[order(tissue_sd, decreasing=T)[1:1000]]

n_tissue = order(tissue_sd, decreasing=T)[1:1000]


baz_mean = matrix(ncol = 9, nrow = 1000)
colnames(baz_mean) = tissue

where_sal = substr(colnames(matrix_normalised_out), 5, 8) == 'Sal'

for (i in 1:9){
  where_tissue = substr(colnames(matrix_normalised_out), 1, 3) == tissue[i]
  baz_mean[, i] = apply(matrix_normalised_out[n_tissue, which(where_sal&where_tissue)], 1, mean)
}

hc_dist = dist.pear(t(baz_mean))
hc = hclust(hc_dist)
plot(as.dendrogram(hc))


#~~~~~~~~~~~~~~~~ klastry odpowiedzi


row = which(pvalue_fdr_int< fdr.thr & max_fold > fold.thr)
row = which(pvalue_fdr_treat< fdr.thr & pvalue_fdr_int> fdr.thr & max_fold > fold.thr)

ind_mean = matrix(ncol = 9, nrow = length(row))
colnames(ind_mean) = tissue

where_dex = substr(colnames(matrix_normalised_st_median), 5, 8) == 'Dex'

matrix_normalised_st_median[row, ] %>% apply(1, scale) %>% t() %>% apply(1, cut.above.thresh, 4) %>% t() -> tmp_stand

for (i in 1:9){
  where_tissue = substr(colnames(matrix_normalised_st_median), 1, 3) == tissue[i]
  ind_mean[, i] = apply(tmp_stand[, which(where_dex&where_tissue)], 1, mean)
}

hc_dist = dist.pear(t(ind_mean))
hc = hclust(hc_dist)
plot(as.dendrogram(hc))


