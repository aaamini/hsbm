library(Matrix)
library(hsbm)
library(nett)
library(ggplot2)
library(dplyr)

# Load data
faonet_path = "data/FAO_Multiplex_Trade/Dataset/"
el = read.table(paste0(faonet_path,"fao_trade_multiplex.edges"))
colnames(el) = c("lay","e1","e2","wei")
node_names = read.table(paste0(faonet_path,"fao_trade_nodes_simplified.txt"), header = T)$nodeLabel %>% 
  gsub("_"," ", . )
layer_names = read.table(paste0(faonet_path,"fao_trade_layers.txt"), header = T)$layerLabel %>% 
  gsub("_"," ", . )


n = 214
nlayers = 364
A = vector("list", nlayers)
for (nl in 1:nlayers) {
  els = el[el$lay == nl,]
  ii = c(els$e1, els$e2)
  jj = c(els$e2, els$e1)
  A[[nl]] = sparseMatrix(i = ii, j = jj, dims = c(n,n)) * 1
  diag(A[[nl]]) = 0
}

# Parameters
nlayers = 20
deg_thresh = 20
grp_thresh = 8
seed = 1400
niter = 2500
burnin = round(niter/2)
Kcap = Gcap = 15
tag = sprintf("faonet")
flags = list(
  load_estimate = T,
  save_estimate = F,
  wclouds = F,
  seq_nmi = T,
  plot_edge_counts = T,
  insamp_dist_plots = T,
  outsamp_dist_plots = F,
  net_plots = F
)

# Select layers
nedges = sapply(A, sum)
idx = order(nedges, decreasing = T)
if (flags$plot_edge_counts) {
  data.frame(layer_idx = seq_along(nedges[idx]), edge_count = nedges[idx]/2) %>% 
    ggplot(aes(layer_idx, edge_count)) + 
    geom_point() +
    scale_x_continuous(trans="log10", breaks=c(1, 2, 5, 10, 20, 50, 100, 200)) +
    xlab("Ordered Layer Index") +
    ylab("Edge Count")+
    theme_minimal()
    ggsave(sprintf("%s_edge_count.pdf", tag), width = 5, height=4)
}

As = A[idx[1:nlayers]]
layer_names = layer_names[idx[1:nlayers]]

# Select nodes
Asum = Reduce(`+`, As)
# image(Asum)
degs = rowSums(Asum)
node_idx = order(degs, decreasing = T)
high_deg_idx = degs[node_idx] > deg_thresh
node_names = node_names[high_deg_idx]
# image(Asum[high_deg_idx, high_deg_idx])

# New adjacency matrix list
Anew = lapply(1:nlayers, function(t) As[[t]][high_deg_idx, high_deg_idx])
nnew = nrow(Anew[[1]])

# Fit HSBM
set.seed(seed)
if (!flags$load_estimate) {
  zh_all = fit_hsbm(Anew, beta0=.1, gam0=.5, niter=niter, Kcap=Kcap, Gcap=Gcap,
                    seq_g_update = F, verb = T)$zb
  if (flags$save_estimate) {
    saveRDS(list(nlayers = nlayers, deg_thresh = deg_thresh, grp_thresh = grp_thresh,
                 seed = seed, tag = tag,
                 niter = niter, burnin = burnin, zh_all = zh_all, 
                 KCap = Kcap, Gcap = Gcap), 
            file = sprintf("results/%s.rds", tag))
    
  }
} else {
  res = readRDS(sprintf("results/%s.rds", tag))
  zh_all = res$zh_all
}

out = get_map_labels(zh_all, burnin = burnin, consecutive = T)
zh = out$labels
conf <- out$conf

# Plot diagnostic sequential NMI plot (for the sampler-chain output)
if (flags$seq_nmi) {
  # zh_all_compact =  compactify_chain_multi_labels(zh_all) # sapply(1:niter, function(t) as.vector(do.call(rbind, zh_all[[t]]))) 
  seq_nmi_plot(zh_all)
  ggsave(sprintf("%s_seq_nmi.pdf", tag), width = 5, height=4)
}
  

# Compute most groups based on most frequent label 
Kest <- max(unlist(zh))
# nett::printf("Estimated number of comm. = %d\n", Kest)
# table(unlist(zh))
xx = lapply(1:nrow(Anew[[1]]), function(j) sapply(1:nlayers, function(t) zh[[t]][j]))
Z = t(sapply(seq_along(xx), function(j) tabulate(xx[[j]], Kest)))
grps = lapply(1:Kest, function(j) which(Z[,j] > grp_thresh))
grps_freq = lapply(1:Kest, function(j) Z[grps[[j]], j])
(frequent_grps = lapply(1:Kest, function(j) node_names[grps[[j]]]))

# Create word clouds
if (flags$wclouds) {
  library(wordcloud)
  for (gid in 1:nlayers) {
    pdf(sprintf("wcloud_g%s.pdf",gid), width = 4, height=4)
    wordcloud(words = frequent_grps[[gid]], freq = grps_freq[[gid]], rot.per = .35, 
              colors=brewer.pal(8, "Dark2"), scale=c(3,.25), random.order = T)
    dev.off()
  }
}

# Plot within/random pair H-Hamming distributions
if (flags$insamp_dist_plots) {
  cat("Calculating in-sample pair nHamming ....")
  out = pair_hamming_test(Anew, grps)
  cat("done.\n")
  out$plot
  out$df %>% filter(type == "Within") %>% pull() %>% median()
  out$df %>% filter(type == "Random") %>% pull() %>% median()
  ggsave(sprintf("%s_pair_hamming.pdf", tag), width = 5, height = 4)
  
}

# This takes time
if (flags$outsamp_dist_plots) {
  cat("Calculating out-sample pair nHamming ....")
  Aout = lapply((nlayers+1):364, function(t) A[[idx[t]]][high_deg_idx, high_deg_idx])
  out2 = pair_hamming_test(Aout, grps)
  cat("done.\n")
  out2$plot
  out2$df %>% filter(type == "Within") %>% pull() %>% median()
  out2$df %>% filter(type == "Random") %>% pull() %>% median()
  ggsave(sprintf("%s_pair_hamming_outsamp.pdf", tag), width = 5, height = 4)
}

# Plot networks
if (flags$net_plots) {
  for (j in 1:nlayers) {
    pdf(paste(tag, layer_names[j],'_gr.pdf',sep=''))
    par(mar=rep(0,4))
    netp = nett::plot_net(igraph::graph_from_adjacency_matrix(Anew[[j]], mode = "undirected"), 
                   zh[[j]], vsize_func = function(deg) log(deg+3)^1.2, 
                   vertex_alpha = 0.75, extract_lcc = T, niter = 2000)
    dev.off()
    
  }
}


# Other tests
reverse_lookup = function(name) {  
  which(grepl(name, node_names))
}
my_corr = function(i,j) {
  # mean(sapply(seq_along(Anew), function(t) sum(Anew[[t]][i,] != Anew[[t]][j,]) / ncol(Anew[[t]])))
  mean(sapply(seq_along(Anew), function(t) hsbm::nhamming(Anew[[t]],i-1,j-1)))
}

my_corr(reverse_lookup("Iraq"),reverse_lookup("^Guinea"))
my_corr(reverse_lookup("Guam"),reverse_lookup("^Grenada"))
my_corr(reverse_lookup("Mauri"),reverse_lookup("Bahr"))
my_corr(reverse_lookup("Macao"),reverse_lookup("Botsw"))
my_corr(reverse_lookup("Macao"),reverse_lookup("Rwa"))
my_corr(reverse_lookup("Fiji"),reverse_lookup("Uganda"))
my_corr(reverse_lookup("Fiji"),reverse_lookup("Vanuatu"))
my_corr(reverse_lookup("Mauri"),reverse_lookup("Bahrain"))
my_corr(reverse_lookup("Germany"),reverse_lookup("USA"))
my_corr(reverse_lookup("Nauru"),reverse_lookup("Falkland Isl."))
my_corr(reverse_lookup("Mont"),reverse_lookup("Alb"))
my_corr(reverse_lookup("Eritrea"),reverse_lookup("American Samoa"))
degs[c(reverse_lookup("Iraq"),reverse_lookup("^Guinea"))]
degs[c(reverse_lookup("Mauri"),reverse_lookup("Bahrain"))]
Anew[[8]][c(reverse_lookup("Mauri"),reverse_lookup("Bahrain")),]
nhamming(Anew[[4]], reverse_lookup("Mauri"), reverse_lookup("Bahrain"))
