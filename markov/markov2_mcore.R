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

# 
# write_progress <- function(msg) {

#     writeLines(paste0(i,'/',nrow(loc),' ',(i/nrow(loc)*100)), fileConn)
#     close(fileConn)
#   }

# All layers jointy
# Enter with eta and z to sample (and also Tn and n)
# n is the number of nodes; Tn the number of layers 
n = 50
Tn_vec <- c(2,4,8,12)
 eta = matrix(c(0.8,0.1,0.3,0.1,0.9,0.2,0.3,0.2,0.7),3,3)
# eta = matrix(c(0.6,0.4,0.3,0.4,0.5,0.2,0.3,0.2,0.1),3,3)
p <- 0.9
K <- 3
CPU_CORES_TO_USE <- detectCores()
#plotG( A[[2]], zb[[2]], deg_func = function(deg) log(deg+3)*2 );

ITRmax <- 2000
burnin <- ceiling(ITRmax/2)
Nrep <- 8
fname <- paste("markov_itr",ITRmax,'_rep',Nrep,'.RData',sep='')
RE_ESTIMATE <- T


run_sims <- function(r) {
  out <- do.call(rbind, lapply( Tn_vec, run_fixed_layer_cound ))
  out$rep <- r
  #setTxtProgressBar(pb, r)
  out
}
run_fixed_layer_cound <- function(Tn) {
  #time1 = proc.time()[3]
  
  zb <- markov_label_process(zinit=sample(1:K,n,replace=TRUE), Tn = Tn, p = p, K = K)
  zb <- lapply(zb, sort) 
  G <- lapply(1:Tn, function(t) sample_sbm2(n, eta, zb[[t]]))
  A <- lapply(1:Tn, function(t) as_adj(G[[t]]) )
  
  result <-  data.frame(type=factor("hsbm", levels=c("hsbm","dpsbm")), agg_nmi=0, slice_nmi= 0, dt = 0, Tn=Tn)
  result[1,"dt"] <- system.time( zh_raw <- hsbm_infer(A, beta0=1, gam0=1, ITRmax=ITRmax, Kcap=20, Gcap=20)$zb )["elapsed"]
  zh <- compMAPLabel(zh_raw, burnin = burnin, consecutive = T)$MAPlabs
  result[1,"agg_nmi"] <- compute_nmi(zb, zh)                                                                                                  
  result[1,"slice_nmi"] <- mean(compute_slice_nmi(zb,zh))
  
  result[2,"Tn"] <- Tn
  result[2,"type"] <- "dpsbm"
  result[2,"dt"] <- system.time( zhi_raw <- dpsbm_slice_infer(A, gam0=1, ITRmax=ITRmax, Zcap=20) )["elapsed"]
  zhi <- compMAPLabel(zhi_raw, burnin = burnin, consecutive = T)$MAPlabs
  result[2,"agg_nmi"] <- compute_nmi(zb, zhi)                        
  result[2,"slice_nmi"] <- mean(compute_slice_nmi(zb,zhi))

  result
}

library(parallel)
if (RE_ESTIMATE) {
  
  elapsed_time <- system.time( 
    result <- do.call(rbind, mclapply(1:Nrep, run_sims, mc.cores = CPU_CORES_TO_USE))
  )["elapsed"]
  printf('Total time = %3.1f', elapsed_time)
  
  save(result, ITRmax, file=fname)
  #save(dt,nmi,hsbm_slice_nmi,dpsbm_slice_nmi,ITRmax,file=fname)
} else {
  load(fname)  
}


library(reshape2)
library(ggplot2)

result2 <- result %>% mutate(Tn = as.factor(Tn))

ggplot(result2, aes(x=Tn,y=agg_nmi,fill=type)) +
  geom_boxplot() + theme_bw() + xlab("Number of layers") + ylab("Aggregate NMI") +
  scale_y_continuous(limits=c(0, 1), expand = c(0, 0.01)) +
  theme(#axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.background=element_blank(), 
    legend.title=element_blank(), 
    legend.position = c(0.85, 0.9),
    legend.text = element_text(size=14),
    text = element_text(size=16)) +
  guides(fill=guide_legend(keywidth=0.25,keyheight=0.25,default.unit="inch")) +
  scale_fill_discrete(labels=c("HSBM","DP-SBM"))
ggsave('markov_agg_nmi.pdf',width = 6,height = 5)


ggplot(result2, aes(x=Tn,y=slice_nmi,fill=type)) +
  geom_boxplot() + theme_bw() + 
  xlab("Number of layers") + ylab("Average Slicewise NMI") + 
  scale_y_continuous(limits=c(0, 1), expand = c(0, 0.01)) +
  theme(legend.background=element_blank(), 
        legend.title=element_blank(), 
        legend.position = c(0.85, 0.1),
        legend.text = element_text(size=14),
        text = element_text(size=16)) +
  guides(fill=guide_legend(keywidth=0.25,keyheight=0.25,default.unit="inch"))+
  scale_fill_discrete(labels=c("HSBM","DP-SBM"))
ggsave('markov_slice_nmi.pdf',width = 6,height = 5)



nmi <- result %>% 
  group_by(type) %>%
  summarise(agg = mean(agg_nmi, na.rm=T), nmi = mean(slice_nmi, na.rm=T), dt= mean(dt, na.rm=T))

printf("\n--- Aggregate NMI, eplapsed time (s) ---\n")
for (j in 1:nrow(nmi)) {
  printf("%10s = %1.5f, %3.2f\n", as.character(nmi[j,]$type), nmi[j,"agg"],nmi[j,"dt"])
}
printf("\n--- Total simulation time (h) = %3.2f\n", result2 %>% summarise(sum(dt)) /3600 )

       


