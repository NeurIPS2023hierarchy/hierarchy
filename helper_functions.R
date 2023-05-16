library(ape)
library(igraph)


Get_rb = function(n, type){
  
  # this function generate random trees, 
  # together with the corresponding correlation matrix to generate Gaussian vectors
  
  # @param n: number of leaf nodes
  # @param type: one of 'random', 'caterpillar', or 'balanced'
  # @return tree: phylo object
  # @return B: n*n matrix, pairwise correlations
  
  if(!(type %in% c('random', 'caterpillar', 'balanced'))){
    stop("invalid type of tree")
  }
  if(type == 'random'){
    tree <- rtree(n) 
    tree$edge.length = tree$edge.length
    tree <- unroot(tree)
    tree <- di2multi(tree, tol = 0.5)
    }
  if(type == 'caterpillar'){
    tree <- stree(n, type = 'left')
    tree <- unroot(tree)
    tree$edge.length = rep(0.5, dim(tree$edge)[1])
    }
  if(type == 'balanced'){
    tree <- stree(n, type = 'balanced')
    tree$edge.length <- rep(0.5, dim(tree$edge)[1]) 
    tree <- unroot(tree)
    }
  dist <- cophenetic(tree)
  B <- exp(-dist)
  return(list(tree = tree, B = B))
}


reconstruct_median_distances = function(dist){
  
  #  each triad (set of three leaves) uniquely determines a "median" node
  #  each latent nodes have multiple triads which have that node as a median. 
  #  this code takes all possible triads,
  #  estimates the distance from each of those three leaves to their median.
  #  then, extends this to compute the distance from this median to all leaves.
  
  # @param dist: pariwise distance matrix of distances among leaves in a tree.
  # @return X: a matrix with dimension m*n, where m = choose(n, 3)
  # @return comb: combn(n, 3), this can be created as long as we know n, but we still save it for later usage
  
  n = nrow(dist)
  comb <- combn(n, 3)
  dist_to_median = matrix(0, nrow = ncol(comb),ncol = n)

  for(i in 1:ncol(comb)){
      triad = comb[,i]
      # get pairwise distances within the triad
      d3 = dist[triad,triad]
      # make it a vector of d12, d13, d23
      d3 = d3[upper.tri(d3)]
      # this estimates the edge length from each triad-node to the median
      w3 = matrix(c(1,1,-1,1,-1,1,-1,1,1), byrow=T,nrow = 3) %*%d3/2
      
      # this sample code computes the distance to the median for one additional node.
      # it is easier to read than the chunk that vectorize it and does it for all nodes:
      # other_nodes = setdiff(1:k, triad)
      # node4 = sample(other_nodes,1)
      # d4 = full_dist[triad,node4]
      # deltas = (d4-w3)/2
      # sum(deltas)-min(deltas)
      # again, this code is a check if you have the full_distance matrix:
      # full_dist[median_node,node4]
      
      #  here is vectorized code that computes distance from all leaves to the median:
      d4 = dist[triad,1:n, drop = F]
      dist_to_median[i,] = apply(d4, 2, function(x,w3){
        deltas = (x-w3)/2
        deltas = sort(deltas)
        deltas[deltas<0] = 0
        sum(deltas)-min(deltas)
      }, w3=w3)
    }
  return(list(dist = dist_to_median, comb = comb))
}


reconstruct_median_distances_to_internal = function(leaf_leaf_dist,leaf_internal_dist, comb, Z){
  
  # for each col in combination, take the median of these three leaf nodes
  # calculate the distance from median to other internal nodes
  # then calcualte the pairwise distance matrix between internal nodes
  
  
  # @param leaf_leaf_dist: n*n matrix, distance between leaf nodes 
  # @param leaf_internal_dist: n*h matrix, distance between leaf nodes and internal nodes
  # @param comb: 3*m matrix, each col being a combination (m = choose(n, 3))
  # @param Z: m*k matrix, class membership matrix estimated by the clustering step
  # @return internal_internal: k*k pairwise distance matrix between latent (internal) nodes
  
  k = ncol(leaf_internal_dist)
  dist_to_median = matrix(0, nrow = ncol(comb),ncol = k)
  for(i in 1:ncol(comb)){
    triad = comb[,i]
    # get pairwise distances within the triad
    d3 = leaf_leaf_dist[triad,triad]
    # make it a vector of d11, d12, d23
    d3 = d3[upper.tri(d3)]
    # this estimates the edge length from each triad-node to the median
    w3 = matrix(c(1,1,-1,1,-1,1,-1,1,1), byrow=T,nrow = 3) %*%d3/2
    #  here is vectorized code that computes distance from all internals to the median:
    d4 = leaf_internal_dist[triad,,drop = F]
    dist_to_median[i,] = apply(d4, 2, function(x,w3){
      deltas = (x-w3)/2
      deltas[deltas<0] = 0
      sum(deltas)-min(deltas)
    }, w3=w3)
  }
  internal_internal <- t(dist_to_median)%*%Z%*%diag(1/colSums(Z), ncol = ncol(Z))
  return(internal_internal)
}


getReluZ <- function(X, k, h){
  n = dim(X)[2]
  trans <- matrix(rnorm(n*h), nrow = n)
  trans_r <- apply(trans, 2, function(x) 1/sqrt(sum(x^2)))
  trans <- trans%*%diag(trans_r)
  trans2 <- X%*%trans
  trans2[trans2<0] <- 0
  trans2 <- trans2[,colSums(trans2)>0]
  expansion2 <- trans2
  fa <- svd(expansion2, nu = k)
  va <- varimax(fa$u, normalize = F)
  Z <- va$loadings
  sign <- ifelse(colSums(Z^3)>0, 1, -1)
  Z <- Z%*%diag(sign)
  Z[Z<0] = 0
  return(Z)
}


postprocessing <- function(X, Z, dist, comb){
  # reconstruct the tree via mst
  
  # @param X: m*n matrix, distance embedding 
  # @param Z: m*k matrix, class membership estimation
  # @param dist: n*n matrix, pairwise distance between leaf nodes
  # @return tree_adj: adj matrix of the recovered tree
  
  leaf_internal_dist <- t(X)%*%Z%*%diag(1/colSums(Z), ncol = ncol(Z))
  internal_internal_dist <- reconstruct_median_distances_to_internal(
    dist, leaf_internal_dist, comb, Z)
  full_dist_est <- rbind(cbind(dist, leaf_internal_dist), 
                         cbind(t(leaf_internal_dist), internal_internal_dist))
  n <- ncol(dist)
  k <- ncol(leaf_internal_dist)
  colnames(full_dist_est)[(n+1):(n+k)] <- paste('m', 1:k, sep = "")
  rownames(full_dist_est) <- colnames(full_dist_est)
  full_dist_est[full_dist_est<0] = 0
  tree_adj <- MaST(exp(-full_dist_est))
  return(tree_adj = tree_adj)
}


GaussianDist <- function(seq){
  rho = cor(seq)
  rho[rho<0] = 0
  e = 0.0001*mean(rho)
  rho = (rho + e)/(1 + e)
  dist = -log(rho)
  return(dist)
}



adjustOrderAdj <- function(adj, tip_labels){
  
  # adjust col/row order of adj matrix. 
  # make sure the first k are tip nodes, following tip_lanels order, 
  # then internal nodes, the order of internal nodes doesn't matter
  
  if(!all(tip_labels %in% colnames(adj))){
    print("invalid tip label")
    return(NA)
  }
  k = length(tip_labels)
  n = ncol(adj)
  tip <- rep(0, k)
  non_tip <- c()
  for(i in 1:n){
    label = colnames(adj)[i]
    if(label %in% tip_labels){
      tip[which(tip_labels==label)] = i
    }else{
      non_tip <- c(non_tip, i)
    }
  }
  order = c(tip, non_tip)
  return(adj[order, order])
}

adjtoSplit <- function(adj, k){
  
  # get split matrix from adj matrix, the first k in adj needs to be leaf nodes
  # the split are of shape something * k
  # the column order is the same as the first k in adj
  
  n = dim(adj)[1]
  split <- c()
  for(i in 1:n){
    for(j in i:n){
      if(adj[i, j]>0){
        curr_adj = adj
        curr_adj[i,j] = 0
        curr_adj[j,i] = 0
        g=graph.adjacency(adjmatrix=curr_adj,
                          mode="undirected",
                          weighted=TRUE,
                          diag=FALSE)
        comp <- components(g)
        split <- c(split, comp$membership[1:k]-1)
      }
    }
  }
  split <- matrix(split, ncol = k, byrow = T)
  colnames(split) <- colnames(adj)[1:k]
  split <- normalizeSplit(split)
  return(split)
}


normalizeSplit <- function(split){
  
  # make sure each row has the first element being 1
  # remove all zero and all 1 split
  # remove duplicated split
  
  n = nrow(split)
  for(i in 1:n){
    if(split[i, 1]==0){
      split[i,] = 1-split[i,]
    }
  }
  select <- (rowSums(split)<ncol(split)-1 & rowSums(split)>1)
  split <- split[select,,drop=FALSE]
  return(unique(split))
}


RF_adj_tree <- function(adj_list, tree){
  
  # calculate the RF distance
  # @param tree: phylo object
  # @param adj_list: list of adj matrix
  
  treesplit = as.matrix(as.splits(tree))
  k = ncol(treesplit)
  treesplit <- normalizeSplit(treesplit)
  dist <- c()
  for(adj in adj_list){
    adjsplit <- adjtoSplit(adj, k)
    n = nrow(treesplit)
    m = nrow(adjsplit)
    mn = nrow(unique(rbind(treesplit,adjsplit)))
    dist <- c(dist, ifelse(m+n>0, (2*mn-m-n)/(m+n), 0))
  }
  return(dist) 
}


nj_shrink <- function(dist, m){
  
  # NJ with post processing that shrink edges
  
  tree <- nj(dist)
  tree <- unroot(tree)
  sort_value = sort(tree$edge.length)
  delete = tree$Nnode - m
  threshold = (sort_value[delete] + sort_value[delete+1])/2
  tree_shrink <- di2multi(tree, tol = threshold) 
  return(list(tree = tree, tree_shrink = tree_shrink))
}








