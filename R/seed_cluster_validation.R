#===========================================================================
#              DETERMINE LCC SIZE SIGNIFICANCE
#===========================================================================
library(tidyverse)
library(igraph)



allCells_trimmed_networks <- readRDS("output/allCells_trimmed_networks.rds")


all_cells_degree_dfs <- vector("list", 3) # df with vertex name, vertex degree and degree interval

for(i in 1:length(Comparison_DE_vec_10W)){
  
  degrees <- degree(allCells_trimmed_networks[[i]], v = V(allCells_trimmed_networks[[i]]), loops= TRUE)
  
  # arrange vertex degrees into intervals such that each interval has at leas 100 degrees
  breaks <- c()
  j <- max(degrees)
  
  while (j >= min(degrees)){
    x <- sort(degrees, decreasing = TRUE)
    k <- j
    if(length(x[x %in% k:j]) < 100){
      while (length(x[x %in% k:j]) < 100){
        k <- k-1
      }
      breaks[length(breaks)+1] <- k              
      j <- k-1
    }else{
      breaks[length(breaks)+1] <- j
      j <- j -1
    }
  }
  
  breaks <- c(sort(breaks), max(degrees))
  
  degree_interval <- cut(degrees, breaks = breaks, include.lowest = TRUE, right = FALSE)
  
  # map each vertex to appropriate degree interval
  
  degree_df <- data.frame(vertex = names(degrees), degree = degrees, interval = degree_interval)
  
  
  all_cells_degree_dfs[[i]] <- degree_df
  names(all_cells_degree_dfs)[i] <- gsub("10Weeks__", "", Comparison_DE_vec_10W[i])
  
}




######################################
# DEGREE PRESERVING RANDOMIZATION:
#=====================================
# Select a random node set (n_randomizations-times) from PPI of size equal to ENTREZ_vector (vector with ENTREZ IDs of DEGs). Perform node sampling from degree bins (defined in degree_df) containing seed genes


DP_LCC_size_randomization <- function(degree_df, n_randomizations, ENTREZ_vector, PPI){
  
  
  sample_random <- function(interval_chr, n_cat){
    
    filtered_df <- degree_df
    filtered_df$interval <- as.character(filtered_df$interval)
    filtered_df_2 <- filtered_df %>% 
      dplyr::filter(interval == interval_chr) %>%
      as.data.frame()
    
    filtered_df_vertices <- as.vector(as.character(filtered_df_2$vertex))
    random_vertices <- sample(filtered_df_vertices, size = n_cat, replace = FALSE)
    return(random_vertices)
  }
  
  
  similar_vertices_LCCs <- function(ENTREZ_vector){
    x <- degree_df %>% filter(vertex %in% ENTREZ_vector)
    x$interval <- factor(x$interval, levels = unique(x$interval))
    y <- as.vector(table(x$interval))
    names(y) <- as.character(levels(x$interval))
    
    
    random_vertices_vec <- mapply(interval_chr = names(y), n_cat = y, FUN = sample_random, SIMPLIFY = FALSE, USE.NAMES = TRUE)
    random_vertices_vec <- unlist(random_vertices_vec)
    
    random_subgraph <- induced_subgraph(PPI, vids = random_vertices_vec, impl = "auto")
    
    C_random_subgraph <- components(random_subgraph)
    LCC_size <- max(C_random_subgraph$csize)
    
    return(LCC_size)
    
    
  }
  
  
  LCC_vec<- replicate(n = n_randomizations, similar_vertices_LCCs(ENTREZ_vector), simplify = TRUE)
  
  
  return(LCC_vec)
  
}

# define vectors with seed genes (for PV and SST get genes with p unadjusted <0.01)


treshold_value_vec <- c(rep(0.1, 3), rep(0.01, 2))
treshold_vec <- c(rep("FDR", 3), rep("P_value", 2))
names_vec <- gsub("10Weeks__", "", Comparison_DE_vec_10W)

get_seeds <- function(treshold, treshold_value, cell_time, PPI){
  if(treshold == "FDR"){
    df <- RML_data_26_May_ENTREZ_NAs_filtered %>% filter(ENTREZ %in% as_ids(V(PPI)), Comparison_DE == cell_time, FDR < treshold_value) %>% as.data.frame()
    vec <- c(as.vector(as.character(df$ENTREZ)))
    return(vec)
  } else{
    df <- RML_data_26_May_ENTREZ_NAs_filtered %>% filter(ENTREZ %in% as_ids(V(PPI)), Comparison_DE == cell_time, P_value < treshold_value) %>% as.data.frame()
    vec <- c(as.vector(as.character(df$ENTREZ)))
    return(vec)
  }
}

seeds_allCells<- mapply(cell_time = Comparison_DE_vec_10W, treshold = treshold_vec, treshold_value = treshold_value_vec, PPI = allCells_trimmed_networks,
                        FUN = get_seeds, USE.NAMES = TRUE, SIMPLIFY = FALSE)


names(seeds_allCells) <- names_vec


# DP randomization results

random_LCCs_all_Cells_10000x <- mapply(degree_df = all_cells_degree_dfs,
                                       ENTREZ_vector = seeds_allCells,
                                       PPI = allCells_unweighted_networks,
                                       n_randomizations = 100000,
                                       FUN = DP_LCC_size_randomization, SIMPLIFY = FALSE, USE.NAMES = TRUE)



saveRDS(random_LCCs_all_Cells_10000x, "random_LCCs_all_Cells_10000x.rds")



# CALCULATE LCC SIZE SIGNIFICANCE


LCC_pval_zscore <- function(ENTREZ, random_LCCs, PPI){
  
  #get LCC of seeds:
  
  ENTREZ_graph <- induced_subgraph(PPI, vids = ENTREZ, impl = "auto")
  Cs_ENTREZ_graph <- components(ENTREZ_graph)
  observed_LCC <- max(Cs_ENTREZ_graph$csize)
  
  # calculate empirical p-value
  pval <- sum(random_LCCs >= observed_LCC)/10000
  
  # calculate z-score
  z_score <- (observed_LCC - mean(random_LCCs))/sd(random_LCCs)
  
  p_z <- c(pval, z_score)
  names(p_z) <- c("p_value", "z_score")
  p_z
}




LCCs_pvals_zscores_allCells <- mapply(random_LCCs = random_LCCs_all_Cells_10000x,
                                       ENTREZ = seeds_allCells, PPI = allCells_trimmed_networks, FUN = LCC_pval_zscore, SIMPLIFY = TRUE, USE.NAMES = TRUE)

