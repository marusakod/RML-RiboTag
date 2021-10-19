# R implementation of the DIAMOnD algorithm proposed by Ghiassian et.al (https://doi.org/10.1371/journal.pcbi.1004120)

library(igraph)
library(tidyverse)

#===============================================================================
#                calculate connectivity significance p-values 
#===============================================================================

# ARGUMENTS
# ks: number of connections to all seeds (initial & diamond), ks0 (number of connections to inital seeds),
# N: number of all the genes in the network
# s: numer of all seed genes (inital seeds & diamond genes that become seeds), s0: number of initial seed genes
# k: number of all connections (degree)
# alpha.: seed weight

get_CS <- function(ks, N, s, s0, k, alpha){ 
  m <- s + (alpha-1)*s0
  x <- ks
  n <- N - s
  K <- k 
  
  p_val <- phyper(q = x-1, m = m, n = n, k = K, lower.tail = FALSE)
  p_val
}



#=======================================================================================
#   DIAMOnD algorithm function: find disease modules based on connectivity significance
#=======================================================================================

# ARGUMENTS
# seeds: vector of seed gene IDs (same ID type as vertex names in PPI), 
# PPI: igraph object
# iterations: number of diamond genes to add to a disease module

diamond <- function(seeds, PPI, iterations = 100, alpha = 1) {
  
  
  initial_seeds <- seeds # a vector of inital seed genes
  
  s0 <- length(initial_seeds) # number of inital seeds
  N <- gorder(PPI) # number of all the vertices in the network
  
  all_seeds <- initial_seeds  # after every iteration one diamond gene is added to all seeds
  
  diamond_genes <- c()  #after every iteration one diamond gene is added to diamond_genes
  pvalues <- c() # after every iteration connectivity significance p-value of the added diamond gene is added to pvalues
  number_of_links_to_seeds <- c() # after every iteration ks of added diamond gene is added to number_of_links_to_seeds
  number_of_all_links <- c() # after every iterationk(degree) of added diamond gene is added to number_of_all_links
  
  
  
  while(length(diamond_genes) < iterations){
    
    s <- length(all_seeds) # number of all seeds
    
    # find neighbors of all seeds
    neighbors_all <- igraph::adjacent_vertices(graph = PPI, v = all_seeds)
    neighbors_all <- unlist(sapply(X = neighbors_all, FUN = as_ids, simplify = TRUE, USE.NAMES = TRUE))
    neighbors_all <- unique(neighbors_all)
    
    # exclude from neighbors_all all the genes that are already in the module:
    neighbors  <- neighbors_all[!neighbors_all %in% all_seeds]
    
    # find neighbors of neighbors
    
    neighbors_of_neighbors<- igraph::adjacent_vertices(graph = PPI, v = neighbors)
    neighbors_of_neighbors<- sapply(X = neighbors_of_neighbors, FUN = as_ids, simplify = TRUE, USE.NAMES = TRUE)
    
    # count how many neighbors of each neighbor are seeds (inital seeds and all seeds) 
    
    ks_of_neighbors <- sapply(neighbors_of_neighbors, function(x){
      kb <- length(dplyr::intersect(x, all_seeds))
      return(kb)
    }
    , simplify = TRUE, USE.NAMES =TRUE)
    
    ks0_of_neighbors <- sapply(neighbors_of_neighbors, function(x){
      kb <- length(dplyr::intersect(x, initial_seeds))
      return(kb)
    }
    , simplify = TRUE, USE.NAMES =TRUE)
    
    
    # get degree of all neighbors:
    
    neighbors_degrees <- degree(graph = PPI, v = neighbors)
    
    #add weights to number of connections to seeds & total number of connections
    
    ks_final <- ks_of_neighbors + (alpha-1)*ks0_of_neighbors
    k_final <- neighbors_degrees +(alpha-1)*ks0_of_neighbors
    
    #reduce the number of genes for which the p value calculation will be performed
    ks_k_df <- tibble(gene = names(ks_final), ks = ks_final, k = k_final) %>%
      # classify the nodes based on their ks and rank the node with lowest k highest within that class
      group_by(ks) %>% 
      filter(k == min(k)) %>% 
      # classify top ranks of each class by their degree k and choose the ones with highest ks
      ungroup() %>% 
      group_by(k) %>%
      filter(ks == max(ks))
    
    # calculate connectivity significance p-value for genes that are in ks_k_df
    
    remaining_genes_ks <- ks_k_df$ks
    names(remaining_genes_ks) <- ks_k_df$gene
    
    remaining_genes_k <- ks_k_df$k
    names(remaining_genes_k) <- ks_k_df$gene
    
    
    CS_pvals <- mapply(ks = remaining_genes_ks, k = remaining_genes_k, N = N,
                       s = s, s0 = s0, alpha = alpha, FUN = get_CS, SIMPLIFY = TRUE, USE.NAMES = TRUE)
    
    
    
    # find the gene with the lowest connectivity significance p-value
    
    new_gene <- names(CS_pvals)[which(CS_pvals == min(CS_pvals))]
    new_gene_ks <- remaining_genes_ks[new_gene]
    new_gene_k <- remaining_genes_k[new_gene]
    p_value <- CS_pvals[new_gene]
    
    # add new gene to the list of DIAMOnD genes 
    
    diamond_genes <- c(diamond_genes, new_gene)
    pvalues <- c(pvalues, p_value)
    number_of_links_to_seeds <- c(number_of_links_to_seeds, new_gene_ks) 
    number_of_all_links <- c(number_of_all_links, new_gene_k) 
    
    # add new gene to the list of seed genes
    
    all_seeds <- c(all_seeds, new_gene)
    
    
  }
  
  # return a list object containing a vector with all module gene IDs and a  dataframe with information on all predicted module genes
  
  added_genes_df <- data.frame(gene = diamond_genes,
                               degree = number_of_all_links,
                               connectivity = number_of_links_to_seeds,
                               pvalue = pvalues)
  
  module_genes_vec <- all_seeds
  
  results_list <- list(module_genes = module_genes_vec,
                       added_genes = added_genes_df)
  
  
  return(results_list)
  
  
}

