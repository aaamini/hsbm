#set.seed(1000)  
set.seed(1400)  
library(igraph)
library(Matrix)
library(dplyr)
library(tidyr)
library(parallel)
#library(R.utils)
printf <- function(...) invisible(cat(sprintf(...)))

# set working directory to current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../posterior_analysis.R", chdir = T)
source("../hsbm_infer_module.R", chdir = T)
source("../network_models.R", chdir = T)

n = 50
Tn = 5
eta = matrix(c(0.9,0.75,0.5,0.75,0.60,0.25,0.5,0.25,0.10),3,3)

zb <- lapply(1:Tn, function(i) sort( sample(c(1,2,3),n,replace=TRUE,prob=c(0.45,0.35,0.25)) ) )
G <- lapply(1:Tn, function(t) sample_sbm2(n, eta, zb[[t]]))
A <- lapply(1:Tn, function(t) as_adj(G[[t]]) )

ITRmax <- 5000
burnin <- ceiling(ITRmax/2)
fname <- paste("MissUSA_postplot_itr",ITRmax,'_burnin',burnin,'.RData',sep='')
RE_ESTIMATE <- F
if (RE_ESTIMATE) {
  zh_raw <- hsbm_infer(A, beta0=1, gam0=1, ITRmax=ITRmax, Kcap=20, Gcap=20)$zb 
  zhi_raw <- dpsbm_slice_infer(A, gam0=1, ITRmax=ITRmax, Zcap=20)

  save(zh_raw,zhi_raw,ITRmax,burnin,file=fname)
}  else {
  load(fname)
}

# fix the zero-base label
for (layer in 1:Tn) {
  hsbm_nmi <- sapply((burnin+1):ITRmax, function(t) compute_nmi(zh_raw[[t]][[layer]]+1,zb[[layer]]))
  dpsbm_nmi <- sapply((burnin+1):ITRmax, function(t) compute_nmi(zhi_raw[[t]][[layer]]+1,zb[[layer]]))
  breaks = seq(0,1,by=0.025)
  # h1 <- hist(dpsbm_nmi,breaks=breaks)
  # h2 <- hist(hsbm_nmi,breaks=breaks)
  # cbind(h1$counts,h2$counts)
  
  # p <- plot_ly(alpha = 0.6) %>%
  #   add_histogram(x = ~dpsbm_nmi) %>%
  #   add_histogram(x = ~hsbm_nmi) %>%
  #   layout(barmode = "overlay")
  
  pdf(sprintf('missUSA_post%d.pdf',layer),width=5, height=4)
  h1 <- hist(hsbm_nmi, breaks=breaks, probability=F, col=rgb(1,0,0,0.5),
       ylab=" ", main="",xlim=c(0,1),cex.lab=1.5, cex.axis = 1.5,cex.main=2,xlab="",lty="blank")
  h2 <- hist(dpsbm_nmi, breaks=breaks, probability=F, col=rgb(0,0,1,0.5), add=T,lty="blank")
  legend(0.75, max(h1$counts)/1.1, legend=c("HSBM","DP-SBM"), col=c(rgb(1,0,0,0.5),
                                                                               rgb(0,0,1,0.5)), pt.cex=3, pch=15, box.lty=0 )
  dev.off()
}


hsbm_nmi <- sapply((burnin+1):ITRmax, function(t) compute_nmi(unlist(zh_raw[[t]])+1,zb))
dpsbm_nmi <- sapply((burnin+1):ITRmax, function(t) compute_nmi(unlist(zhi_raw[[t]])+1,zb))

pdf(sprintf('missUSA_post_agg.pdf',layer),width=5, height=4)
h1 <- hist(hsbm_nmi, breaks=breaks, probability=F, col=rgb(1,0,0,0.5),
     ylab="", main="",xlim=c(0,1),cex.lab=1.5, cex.axis = 1.5,cex.main=2,xlab="",lty="blank")
h2 <- hist(dpsbm_nmi, breaks=breaks, probability=F, col=rgb(0,0,1,0.5), add=T,lty="blank")
legend(0.75, max(h1$counts)/1.1, legend=c("HSBM","DP-SBM"), col=c(rgb(1,0,0,0.5),                                                                  rgb(0,0,1,0.5)), pt.cex=3, pch=15, box.lty=0 )
dev.off()




# library(plotly)
# plot_ly(x= h1$mids, y=h1$counts, alpha=0.6) %>% add_bars() %>%
#   add_bars(x= h2$mids, y=h2$counts, mode="bars") %>%
#   layout(barmode = "overlay")
# p2 <- plot_ly(x= h2$mids, y=h2$counts, alpha=0.6) %>% add_bars() 
# plot(density(hsbm_nmi,bw=0.01),col="blue")
# lines(density(dpsbm_nmi,bw = .01),col="red")
