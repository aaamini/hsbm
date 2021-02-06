#set.seed(1000)  
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

################## Data-set example simulation###############################
sample_sbm2 <- function(n, eta, z) {
  tab <- table(z)
  idx <- as.integer(names(tab))
  
  sample_sbm(n, pref.matrix = eta[idx,idx], block.sizes = tab)
}
#############################################################

# All layers jointy
#Enter with eta and z to sample (and also Tn and n)
# n is the number of nodes; Tn the number of layers 

n <- 50
Tn <- 5
eta <- matrix(c(0.9,0.75,0.5,0.75,0.60,0.25,0.5,0.25,0.10),3,3)


ITRmax <- 5000
burnin <- ceiling(0.6*ITRmax)
Nrep <- 20
fname <- paste("MissUSA_itr",ITRmax,'_rep',Nrep,'.RData',sep='')
RE_ESTIMATE <- F


run_sims <- function(r) {
  nmi <- dt <- data.frame(hsbm=0, dpsbm=0)
  
  zb <- lapply(1:Tn, function(i) sort( sample(c(1,2,3),n,replace=TRUE,prob=c(0.45,0.35,0.25)) ) )
  G <- lapply(1:Tn, function(t) sample_sbm2(n, eta, zb[[t]]))
  A <- lapply(1:Tn, function(t) as_adj(G[[t]]) )
  
  
  dt["hsbm"] <- system.time( zh_raw <- hsbm_infer(A, beta0=1, gam0=1, ITRmax=ITRmax, Kcap=20, Gcap=20)$zb )["elapsed"]
  
  dt["dpsbm"] <- system.time( zhi_raw <- dpsbm_slice_infer(A, gam0=1, ITRmax=ITRmax, Zcap=20) )["elapsed"]
  
  zh <- compMAPLabel(zh_raw, burnin = burnin, consecutive = T)$MAPlabs
  nmi["hsbm"] <- compute_nmi(zb, zh)                                                                                                                               
  hsbm_slice_nmi <- compute_slice_nmi(zb,zh)
  
  zhi <- compMAPLabel(zhi_raw, burnin = burnin, consecutive = T)$MAPlabs
  nmi["dpsbm"] <- compute_nmi(zb, zhi)
  dpsbm_slice_nmi <- compute_slice_nmi(zb,zhi)
  
  list(dt=dt,nmi=nmi,hsbm_slice_nmi=hsbm_slice_nmi, dpsbm_slice_nmi=dpsbm_slice_nmi)
}

library(parallel)
if (RE_ESTIMATE) {
  out <- mclapply(1:Nrep, run_sims, mc.cores = 1)
  
  # unpack mclapply output
  fields <- names(out[[1]])
  for (j in 1:length(fields)) {
    var_name <- fields[[j]]
    assign(var_name, do.call(rbind, lapply(1:Nrep, function(r) out[[r]][[var_name]] )))
  }
  
  save(dt,nmi,hsbm_slice_nmi,dpsbm_slice_nmi,ITRmax,file=fname)
} else {
  load(fname)  
}

library(reshape2)
library(ggplot2)


#nmi <- data.frame(hsbm=hsbm_nmi, dpsbm=dpsbm_nmi)
colnames(nmi) <- c("HSBM","DP-SBM")
nmi_m <- melt(nmi)
pdf('missUSA_agg_nmi.pdf')
ggplot(nmi_m, aes(x=variable,y=value,fill=variable)) + 
  geom_boxplot() + theme_bw() + xlab("") + ylab("Aggregate NMI") +
  scale_y_continuous(limits=c(0, 1), expand = c(0, 0.01)) +
  theme(#axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.9),
        legend.text = element_text(size=16),
        text = element_text(size=16)) +
  guides(fill=guide_legend(keywidth=0.5,keyheight=0.5,default.unit="inch")) +
  scale_fill_discrete(labels=c("HSBM","DP-SBM"))
dev.off()
 
slice_nmi_m <- rbind( cbind(melt(hsbm_slice_nmi), type="hsbm"), cbind(melt(dpsbm_slice_nmi), type="dpsbm") ) 
years <- 2014:2018
slice_nmi_m <- slice_nmi_m %>% mutate(Var2 = years[Var2])
pdf('missUSA_slice_nmi.pdf')
ggplot(slice_nmi_m, aes(x=factor(Var2),y=value,fill=type)) + 
  geom_boxplot() + theme_bw() + 
  xlab("Year") + ylab("Slicewise NMI") + 
  scale_y_continuous(limits=c(0, 1), expand = c(0, 0.01)) +
  theme(legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.1),
        legend.text = element_text(size=16),
        text = element_text(size=16)) +
    guides(fill=guide_legend(keywidth=0.5,keyheight=0.5,default.unit="inch"))+
  scale_fill_discrete(labels=c("HSBM","DP-SBM"))
dev.off()
#dpsbm_slice_nmi <- compute_slice_nmi(zb,zhi)

colnames(nmi) <- c("hsbm","dpsbm")
printf("\n--- Aggregate NMI, eplapsed time (s) ---\n%10s = %1.5f, %3.2f\n%10s = %1.5f, %3.2f\n",
       "HSBM", mean(nmi$hsbm), mean(dt$hsbm),
       "DP-SBM", mean(nmi$dpsbm), mean(dt$dpsbm))


