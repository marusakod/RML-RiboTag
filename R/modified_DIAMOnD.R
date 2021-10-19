# R implementation of the modified DIAMOnD algorithm for identification of disease modules based on connectivity significance and differential expression data.
# The algorithm consists of the following steps:
# 1.) for all genes connected to any of the seed nodes Wald test p-values and connectivity significance p-values are determined.
# 2.) candidate disease module genes are ranked based on their respective Wald test p-values and connectivity signficance
# 3.) individual rankings are combined into a single score given as the sum of the inverted ranks
# 4.) the gene with the highest combined score is included into a disease module
# 5.) steps 1-4 can be repeated until all the nodes in the network are included in the disease module



library(tidyverse)
library(igraph)


#===============================================================================
#                calculate connectivity significance p values 
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



#================================================
#      modified DIAMOnD algorithm function
#================================================

# ARGUMENTS
# seeds: a vector of seed gene IDs (ENTREZ)
# PPI: igraph object
# iterations: the number of genes in the final disease module
# alpha: seed weight
# gene_expression_data: a dataframe with a column named ENTREZ (ENTREZ gene IDs) and a column named P_value (Wald test p-value)


modified_diamond <- function(seeds, PPI, iterations = 100, alpha = 1, gene_expression_data) {
  
  
  initial_seeds <- seeds # vector of inital seed genes
  
  s0 <- length(initial_seeds) # number of initial seeds
  N <- gorder(PPI) # number of all the vertices in the network
  
  all_seeds <- vector("character", length = s0 + iterations) #after every iteration one diamond gene is added to all_seeds
  all_seeds[1:s0] <- initial_seeds
  
  diamond_genes <- vector("character", length = iterations)  #after every iteration one diamond gene is added to diamond_genes.
  CS_pvals_allGenes <- vector("numeric", length = iterations) # after every iteration connectivity significance p-value of a new gene is added to CS_pvals_allGenes
  ks_all_genes <- vector("integer", length = iterations) # after every iteration ks (nummber of links to seeds) of a new gene is added to ks_all_genes
  k_all_genes <- vector("integer", length = iterations) # after every iteration k (degree) of a new gene is added to k_all_genes
  DE_pval_allGenes <- vector("numeric", length = iterations)   # after every iteration Wald test p-value of a new gene is added to DE_pval_allGenes
  
  
  
  
  i <- 1 #initialze a variable counting iterations steps
  
  while(i <= iterations){
    
    s <- length(which(all_seeds != "")) # number of all seeds
    
    all_seeds_defined <- all_seeds[which(all_seeds != "")]
    
    # find neighbors of all seeds
    neighbors_all <- igraph::adjacent_vertices(graph = PPI, v = all_seeds_defined)
    neighbors_all <- unlist(sapply(X = neighbors_all, FUN = as_ids, simplify = TRUE, USE.NAMES = TRUE))
    neighbors_all <- unique(neighbors_all)
    
    # exclude from neighbors_all all the genes that are already in the module
    neighbors  <- neighbors_all[!neighbors_all %in% all_seeds_defined]
    
    # find neighbors of neighbors
    
    neighbors_of_neighbors<- igraph::adjacent_vertices(graph = PPI, v = neighbors)
    neighbors_of_neighbors<- sapply(X = neighbors_of_neighbors, FUN = as_ids, simplify = TRUE, USE.NAMES = TRUE)
    
    # count how many neighbors of each neighbor are seeds (inital seeds and all seeds)
    
    ks_of_neighbors <- sapply(neighbors_of_neighbors, function(x){
      kb <- length(dplyr::intersect(x, all_seeds_defined))
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
    
    
    # calculate connectivity significance p value for all neighbors
    
    CS_pvals <- mapply(ks = ks_final, k = k_final, N = N,
                       s = s, s0 = s0, alpha = alpha, FUN = get_CS, SIMPLIFY = TRUE, USE.NAMES = TRUE)
    
    # rank CS_pvals
    
    CS_pvals_ranks <- rank(CS_pvals, na.last = TRUE, ties.method = "average" )
    
    #first df with ranks (just CS ranks):
    
    rank_df_1 <- data.frame(ENTREZ = neighbors,
                            CS_Pvals = CS_pvals,
                            CS_Pval_ranks = CS_pvals_ranks)
    
    
    # get pvalues of all neighbors
    
    diff_expr_data_of_neigbors <- gene_expression_data %>% 
      filter(ENTREZ %in% neighbors) %>%
      group_by(ENTREZ) %>%
      top_n(1, -log10(P_value)) %>% # if two genes have the same ENTREZ ID keep the gene with lower pval
      select(ENTREZ, P_value) %>%
      as.data.frame()
    
    # make 2nd df with ranks (both ranks) and a combined score
    
    rank_df_2 <- merge(x = rank_df_1, y = diff_expr_data_of_neigbors, by = "ENTREZ")
    rank_df_2$diff_exp_pval_rank <- rank(rank_df_2$P_value, na.last = TRUE, ties.method = "average")
    
    rank_df_2  <- rank_df_2 %>%
      rowwise() %>%
      mutate(combined_score = 1/CS_Pval_ranks + 1/diff_exp_pval_rank) %>% 
      as.data.frame()
    
    
    # gene with the maximum combined score
    
    max_score <-  max(rank_df_2[, "combined_score"])
    new_gene <- rank_df_2[rank_df_2$combined_score == max_score, "ENTREZ"]
    
    ks_new_gene <- ks_final[new_gene]
    k_new_gene <- k_final[new_gene]
    CS_new_gene <- CS_pvals[new_gene]
    DE_pval_new_gene <- rank_df_2[rank_df_2$ENTREZ == new_gene, "P_value"]
    
    
    
    # add new gene to the list of diamond genes 
    diamond_genes[i] <- new_gene
    CS_pvals_allGenes[i] <- CS_new_gene
    ks_all_genes[i] <- ks_new_gene
    k_all_genes[i] <- k_new_gene
    DE_pval_allGenes[i] <- DE_pval_new_gene
    
    
    
    # add new gene to the list of seed genes
    
    all_seeds[s0 +i] <- new_gene
    
    # increase i by 1
    i <- i +1
    
  }
  
  
  # return a list object containing a vector with all module gene IDs and a dataframe with information on all predicted module genes
  
  added_genes_df <- data.frame(gene = diamond_genes,
                               degree = k_all_genes,
                               connectivity = ks_all_genes,
                               CS_pvalue = CS_pvals_allGenes,
                               DE_pvalue = DE_pval_allGenes)
  
  module_genes_vec <- all_seeds
  
  results_list <- list(module_genes = module_genes_vec,
                       added_genes = added_genes_df)
  
  
  return(results_list)
  
  
}

