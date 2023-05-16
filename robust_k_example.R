rm(list=ls())
source('helper_functions.R')
random_rb <- readRDS('data/example_for_robust_k_8.RDS')
n = 8

set.seed(1)
seq_length <- 500000
seq <- mvrnorm(n = seq_length, mu = rep(0, n), Sigma = random_rb$B)
k = Nnode(random_rb$tree)
colnames(seq) <- random_rb$tree$tip.label
dist <- GaussianDist(seq)
step1 <- reconstruct_median_distances(dist)
X = step1$dist
comb = step1$comb
weight = step1$weight
h = 500
Z_known <- getReluZ(X, k = k, h)
ECR_adj_known <- postprocessing(X, Z_known, dist, comb)

Z_unknown <- getReluZ(X, k = n-2, h)
ECR_adj_unknown <- postprocessing(X, Z_unknown, dist, comb)

plot_igraph <- function(igraph, tip_label){
  vertex_name <- names(igraph[1])
  vertex_name <- sapply(vertex_name, function(x){ifelse(x %in% tip_label, x, '')})
  plot(igraph, vertex.label = vertex_name, vertex.size = 25)
}


RF_adj_tree(list(ECR_adj_known, ECR_adj_unknown),random_rb$tree)

pdf(file = "plots/robust_k.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 2)
par(mfrow = c(1, 3), mar = rep(0, 4))
igraph <- as.igraph(random_rb$tree)
plot_igraph(igraph, random_rb$tree$tip.label)
igraph <- graph_from_adjacency_matrix(ECR_adj_known, "undirected", diag = F, weighted = T)
plot_igraph(igraph, random_rb$tree$tip.label)
igraph <- graph_from_adjacency_matrix(ECR_adj_unknown, "undirected", diag = F, weighted = T)
plot_igraph(igraph, random_rb$tree$tip.label)
dev.off()

