# ----------------------------------------------------------------------------
# Author: Shamim Mollah  
# Created: 12-10-2016
# 
# Generate functional modules using NMF
#-----------------------------------------------------------------------------

#set directory
setwd("/Users/smollah")
rm(list = ls())
histon_data<-read.table("60histone_24hr.txt", sep="\t", header=T, row.names=1)
YNormZ<-scale(as.matrix(histon_data[,]))

library(NMF)

#Normalized Data to [0,1]
fun <- function(x){ 
     a <- min(x) 
     b <- max(x) 
     (x - a)/(b - a) 
} 

# centered histone data
hi=as.data.frame(YNormZ)
#hi=as.data.frame(histon_data)
#2^ log fc value
YNormZ_anti=as.data.frame(sapply(hi, function(x) 2^x))

#[0-1] range
mat_hist <- apply(YNormZ_anti, 2, fun) 
#mat_hist <- YNormZ2


#to identify optimal factorization rank
V.random <- randomize(mat_hist)
estim.r.random <- nmf(V.random, 2:10, nrun = 1000, seed = 123456)
estim.r <- nmf(mat_hist, 2:10, nrun = 1000, seed = 123456)
plot(estim.r, estim.r.random)

#with nndsvd seed, k=4
res_hist <- nmf(mat_hist, 4, "brunet", nrun=1000, seed = "nndsvd")
#consensusmap(res_hist)
#to visually plot the basis and loading matrices
layout(cbind(1, 2))
basismap(res_hist, labRow= row.names(YNormZ))
coefmap(res_hist)

# to get matrix of basis and loadings 
w_hist<-basis(res_hist)
h_hist<-coef(res_hist)
write.table(w_hist, "centered_histone_drug_result_w_basis_4clust.txt", sep="\t")
write.table(h_hist, "centered_histone_drug_result_h_loading_4clust.txt", sep="\t")