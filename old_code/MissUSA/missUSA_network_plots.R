set.seed(1400)  
library(igraph)
library(Matrix)
library(dplyr)
#library(R.utils)
printf <- function(...) invisible(cat(sprintf(...)))

# set working directory to current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../posterior_analysis.R", chdir = T)
source("../hsbm_infer_module.R", chdir = T)
source("../network_models.R", chdir = T)

n <- 50
Tn <- 5
eta <- matrix(c(0.9,0.75,0.5,0.75,0.60,0.25,0.5,0.25,0.10),3,3)

zb <- lapply(1:Tn, function(i) sort( sample(c(1,2,3),n,replace=TRUE,prob=c(0.45,0.35,0.25)) ) )
G <- lapply(1:Tn, function(t) sample_sbm2(n, eta, zb[[t]]))
A <- lapply(1:Tn, function(t) as_adj(G[[t]]) )

for (j in 1:Tn) {
  pdf(paste('missUSA_G',j,'.pdf',sep=""))
  plotG( A[[j]], zb[[j]], deg_func = function(deg) log(deg+3)*3 )
  dev.off()
}

for (j in 1:Tn) {
  pdf(paste('missUSA_A',j,'.pdf',sep=""))
  image( A[[j]], xlab="", sub="", ylab="",useRaster=TRUE)
  dev.off()
}

