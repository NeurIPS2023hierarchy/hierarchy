library(MASS)

Get_result_binary_tree_gaussian <- function(n_rep, n, h = 500, seq_length, type){
  # theta_dist: 'constant', 'rexp', 'rexp square'
  # n_rep: replicate time
  # n: number of leaf nodes
  # h: expansion dimension
  
  # k: number of internal nodes
  nmethods = 2
  nseq = length(seq_length)
  RF_dist <-matrix(0, n_rep, nseq*nmethods)
  
  for(run in 1:n_rep){
    rb = Get_rb(n, type)
    k = rb$tree$Nnode
    print(paste('run:', run, '; k:', k))
    RF_dist_temp <-c()
    for(j in 1:nseq){
      l = seq_length[j]
      print(paste('sequence length:', l))
      seq <- mvrnorm(n = l, mu = rep(0, n), Sigma = rb$B^0.6)
      colnames(seq) <- rb$tree$tip.label
      dist <- GaussianDist(seq)
      rownames(dist) <- colnames(dist) <- rb$tree$tip.label
      
      # now know the number of internal nodes
      nj_tree <- nj(dist)
      igraph <- as.igraph(nj_tree, directed = F)
      adj <- as.matrix(as_adjacency_matrix(igraph))
      diag(adj) <- 1
      nj_adj <- adjustOrderAdj(adj, rb$tree$tip.label)
      
      # do step1 just once as its time consuming
      step1 <- reconstruct_median_distances(dist)
      X = step1$dist
      comb = step1$comb
      Z <- getReluZ(X, k = n-2, h)
      ECR_adj <- postprocessing(X, Z, dist, comb)

      dist_rf <- RF_adj_tree(list(nj_adj,
                                  ECR_adj),
                             rb$tree)
      RF_dist_temp = c(RF_dist_temp, dist_rf)
    }
    RF_dist[run, ] = RF_dist_temp
    if(run==1){
      print(t(matrix(RF_dist[1:run, ], ncol = nseq)))
    }else{
      print(t(matrix(colMeans(RF_dist[1:run, ]), ncol = nseq)))
    }
  }
  return(RF_dist)
}



Get_result_non_binary_tree_gaussian <- function(n_rep, n, h = 500, seq_length){
  
  nmethods = 4
  nseq = length(seq_length)
  RF_dist <-matrix(0, n_rep, nseq*nmethods)
  
  for(run in 1:n_rep){
    rb = Get_rb(n, type = 'random')
    k = rb$tree$Nnode
    print(paste('run:', run, '; k:', k))
    RF_dist_temp <-c()
    for(j in 1:nseq){
      l = seq_length[j]
      print(paste('sequence length:', l))
      seq <- mvrnorm(n = l, mu = rep(0, n), Sigma = rb$B^0.6)
      colnames(seq) <- rb$tree$tip.label
      dist <- GaussianDist(seq)
      rownames(dist) <- colnames(dist) <- rb$tree$tip.label
      
      nj_shrink <- nj_shrink(dist, k)
      igraph <- as.igraph(nj_shrink$tree_shrink, directed = F)
      adj <- as.matrix(as_adjacency_matrix(igraph))
      diag(adj) <- 1
      split <- as.splits(rb$tree)
      split <- as.matrix(split)
      nj_adj <- adjustOrderAdj(adj, colnames(split))
      igraph <- as.igraph(nj_shrink$tree, directed = F)
      adj <- as.matrix(as_adjacency_matrix(igraph))
      diag(adj) <- 1
      split <- as.splits(rb$tree)
      split <- as.matrix(split)
      nj_adj_no_shrink <- adjustOrderAdj(adj, colnames(split))
      
      # do step1 just once as its time consuming
      step1 <- reconstruct_median_distances(dist)
      X = step1$dist
      comb = step1$comb
      
      Z_known <- getReluZ(X, k = k, h)
      ECR_adj_known <- postprocessing(X, Z_known, dist, comb)
      
      Z_unknown <- getReluZ(X, k = n-2, h)
      ECR_adj_unknown <- postprocessing(X, Z_unknown, dist, comb)
      
      dist_rf <- RF_adj_tree(list(nj_adj,
                                  nj_adj_no_shrink,
                                  ECR_adj_known,
                                  ECR_adj_unknown),
                             rb$tree)
      RF_dist_temp = c(RF_dist_temp, dist_rf)
    }
    RF_dist[run, ] = RF_dist_temp
    if(run==1){
      print(t(matrix(RF_dist[1:run, ], ncol = nseq)))
    }else{
      print(t(matrix(colMeans(RF_dist[1:run, ]), ncol = nseq)))
    }
  }
  return(RF_dist)
}
