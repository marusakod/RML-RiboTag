#===========================================================================================
#  CLUSTERING DISEASE MODULE GENES, OVERREPRESENTATION ANALYSIS OF DISEASE MODULE CLUSTERS
#===========================================================================================

library(igraph)
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(simplifyEnrichment)


# Import diamond results 

modified_diamond_results <- list(alpha_1 = readRDS("output/modified_diamond_results/low_pval_diamond_alpha_1_iterations_500.rds"),
                                         alpha_2 = readRDS("output/modified_diamond_results/low_pval_diamond_alpha_2_iterations_500.rds"),
                                         alpha_5 = readRDS("output/modified_diamond_results/low_pval_diamond_alpha_5_iterations_500.rds"))

# import RML data

RML_data_26_May_ENTREZ_NAs_filtered <- readRDS("data/RML_data_26_May_ENTREZ_NAs_filtered.rds")
RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE <- gsub(" ", "", RML_data_26_May_ENTREZ_NAs_filtered$Comparison_DE)


# import PPI networks

allCells_trimmed_networks <- readRDS("output/allCells_trimmed_networks.rds")


# get list of vectors with seed gene IDs

Comparison_DE_vec_10W <- c("10Weeks__vGluT2", "10Weeks__Gad2","10Weeks__Cx43", "10Weeks__PV", "10Weeks__SST")

DEGs_list <- vector("list", 3)

for(i in 1:length(Comparison_DE_vec_10W[1:3])){
  
  vertices <- as_ids(V(allCells_trimmed_networks[[i]]))
  df <- RML_data_26_May_ENTREZ_NAs_filtered %>%
        dplyr::filter(ENTREZ %in% vertices) %>%
        filter(Comparison_DE == Comparison_DE_vec_10W[i] & FDR< 0.1) %>% as.data.frame()
 
   vec <- as.character(df$ENTREZ)
  DEGs_list[[i]] <- vec
  names(DEGs_list)[i] <- gsub("10Weeks__", "", Comparison_DE_vec_10W[i])
  
}

#==============================================================================
# assign cluster membership to module genes and make a df with all node data
#==============================================================================
# ARGS:
# seeds: vector of ENTREZ gene IDS
# diamond_results
# threshold: the # of diamond genes to be included in the module
# PPI: graph object with ENTREZ gene IDs as vertex names (different for all cell types)
# celltype: string (e.g., "10Weeks__Gad2")

get_cluster_membership_in_modules <- function(seeds, diamond_results, threshold, PPI, celltype){
  
  # get diamond genes in the module
  diamond_genes <- diamond_results$added_genes$gene[1:threshold]
  
  # merge diamonds and seeds to get all genes in the module
  all_module_genes <- c(seeds, diamond_genes)
  
  # make df with module genes and iteration #
  all_module_genes_df <- data.frame(Id = all_module_genes,
                                    iteration = c(rep(0, length(seeds)), c(1:threshold)))
  
  # make a graph object with module vertices from the filtered PPI
  module_subgraph <- induced_subgraph(graph = PPI, vids = all_module_genes, impl = "auto")
  
  # decompose module igraph object to remove all single genes from module
  decomposed_graph <- decompose(module_subgraph, min.vertices = 2)
  
  # merge components with more than 2 vertices in a single igraph
  decomposed_graph_vertices <- unlist(sapply(X = decomposed_graph, FUN = function(x){as_ids(V(x))}, simplify = TRUE))
  
  recomposed_graph <- igraph::simplify(induced_subgraph(graph = PPI, vids = decomposed_graph_vertices, impl = "auto")) # use igraph::simplify to remove loops
  
  # do clustering on modules with no singletons
  cluster_fg_result <- cluster_fast_greedy(graph = recomposed_graph, weights = NULL)
  cluster_lp_result <- cluster_label_prop(graph = recomposed_graph, weights = NULL)
  cluster_Wt_results <- cluster_walktrap(graph = recomposed_graph, weights = NULL, steps = 8)
  
  #make dataframe with vertices which can be assigned cluster membership
  df1 <- data.frame(ENTREZ = cluster_fg_result$names,
                    clust_membership_fg = cluster_fg_result$membership,
                    clust_membership_lp = cluster_lp_result$membership,
                    clust_membership_wt = cluster_Wt_results$membership)
  
  #make dataframe with remaining vertices)
  remaining_genes <- setdiff(all_module_genes, decomposed_graph_vertices)
  
  df2 <- data.frame(ENTREZ = remaining_genes, 
                    clust_membership_fg = rep(0, length(remaining_genes)),
                    clust_membership_lp = rep(0, length(remaining_genes)),
                    clust_membership_wt = rep(0, length(remaining_genes)))
  
  
  # add gene origin (seed or diamond) to df
  
  df3 <- bind_rows(df1, df2) %>%
    mutate(origin = ifelse(ENTREZ %in% seeds, "seed", "diamond")) %>%
    mutate(cell_type = gsub("10Weeks__", "", celltype))
  
  # add gene symbol to df
  
  df3$gene_symbol <- mapIds(org.Mm.eg.db, df3$ENTREZ, "SYMBOL", "ENTREZID")
  
  
  # add DE P-value to df
  
  RML_data_filtered <- RML_data_26_May_ENTREZ_NAs_filtered %>% 
    filter(Comparison_DE == celltype) %>%
    group_by(ENTREZ) %>%
    top_n(1, -log10(P_value)) %>% # if more P-values correspond to the same ENTREZ Id, keep only the min p-value
    ungroup() %>%
    as.data.frame()
  
  df3 <- merge(df3, RML_data_filtered, by = "ENTREZ") 
  
  df3 <- df3 %>%
    dplyr::select(ENTREZ, clust_membership_fg, clust_membership_lp, clust_membership_wt, origin, cell_type, gene_symbol, P_value, FDR, L2FC) %>%
    dplyr::rename(Id = ENTREZ) %>%
    dplyr::rename(Label = gene_symbol) %>%
    dplyr::rename(DE_pvalue = P_value) %>%
    dplyr::mutate(direction = ifelse(L2FC < 0, "downregulated", "upregulated"))
  
  # merge df3 with df with df that has iteration #
  
  
  df4 <- merge(df3, all_module_genes_df, by = "Id")
  
  df4 <- df4 %>% dplyr::select(Id, Label, clust_membership_fg, clust_membership_lp, clust_membership_wt, origin, cell_type, DE_pvalue, FDR, iteration, L2FC, direction)
  
  df4
  
}


all_cells_cluster_memberships <- mapply(seeds = DEGs_list,
                                        diamond_results = modified_diamond_results$alpha_1,
                                        threshold = c(120, 100, 215),
                                        PPI = allCells_trimmed_networks[c(1:3)],
                                        celltype = Comparison_DE_vec_10W[1:3],
                                        FUN = get_cluster_membership_in_modules,
                                        SIMPLIFY = FALSE,
                                        USE.NAMES = TRUE
                                        
                                    
)

#=================================================================================
#       Overrepresentation analysis on module clusters
#=================================================================================


source("R/thresholds_validation.R")


# ORA_for_clusters ARGUMENTS:
# clustering_results: output of get_cluster_membership_in_modules()
# cluster_method: one of "clust_membership_fg", "clust_membership_wt", "clust_membership_lp" (fg = fast greedy, lp = label propagation, wt = walktrap)
# PPI: graph object
# gene_set = "GO" / "KEGG" / "REACTOME"/ "BIOCARTA"


ORA_for_clusters <- function(clustering_results, cluster_method = c("clust_membership_fg", "clust_membership_wt", "clust_membership_lp"), PPI, gene_set){ 
 
  #select method used for clustering
   clust_method <- match.arg(cluster_method)
  
   # remove 0 from cluster indices 
   
    cluster_indices <- unique(clustering_results[, clust_method])
    cluster_indices <- cluster_indices[which(cluster_indices != 0)]
    
   #get all genes from clusters that have index different that 0
  
   vecs_with_cluster_genes <- vector("list", length(cluster_indices))
  
  for(i in 1:length(cluster_indices)){
    
    keep <- which(clustering_results[, clust_method] == cluster_indices[i])
    
    keep_genes <- clustering_results[keep, "Id"]   
    
    vecs_with_cluster_genes[[i]] <- keep_genes
    
  }
  
  # do ORA on all clusters
  
  ORA_results <- sapply(X = vecs_with_cluster_genes,
                        FUN = supersets_enricher,
                        simplify = FALSE,
                        USE.NAMES = TRUE,
                        allgenes = as_ids(V(PPI)), # gene universe = all annotated vertices in the PPI network
                        padjust = 0.05,
                        set = gene_set,
                        maxSS = 500)
  
  
  # make list with cluster genes and corresponding ORA results 
  
  final <- mapply(x = vecs_with_cluster_genes,
                  y = ORA_results,
                  FUN = function(x,y){list(cluster_genes = x,
                                           ORA_results = y)},
                  SIMPLIFY = FALSE,
                  USE.NAMES = FALSE)
  names(final) <- paste("cluster", cluster_indices, sep = "_")
  
  final
  
}


allCells_FG_clustering_GO_ORA <- mapply(clustering_results = all_cells_cluster_memberships,
                                        cluster_method = "clust_membership_fg",
                                        PPI = allCells_trimmed_networks[c(1:3)],
                                        gene_set = "GO",
                                        FUN = ORA_for_clusters, 
                                        SIMPLIFY = FALSE,
                                        USE.NAMES = TRUE)




###################################################################################
#         remove redundat GO terms from ORA results 
#==================================================================================


reduce_GO_enrichment_results <- function(ORA_result){ # oRA_result = output of supersets_enricher()
  
  
  # Add GO ID to ORA_results
  full_term_df <- MSigDb_GO_KEGG_REACTOME_BIOCARTA[, c(1,3)] %>% distinct() %>% dplyr::rename(ID = gs_name)
  
  ORA_result_full <- merge(ORA_result, full_term_df, by = "ID")
  
  
  # separate ORA_results db by GO category

  ORA_results_filtered <- sapply(X = c("GO:BP", "GO:CC", "GO:MF"),
                               FUN = function(x){
                                 ORA_result_full %>% filter(str_detect(ID, paste("^", x, sep = "")))
                               },
                               simplify = FALSE, 
                               USE.NAMES = TRUE)
  
  # make vectors with significant GO IDs
  
  significant_GO_IDS <- sapply(X = ORA_results_filtered,
                               FUN = function(x){
                                 x %>% select(gs_exact_source) %>% flatten_chr()
                               },
                               simplify = FALSE,
                               USE.NAMES = TRUE)
  
  # for each GO category get the # of significant IDs
  significant_GO_IDs_count <- sapply(X = significant_GO_IDS, FUN = length, simplify = TRUE, USE.NAMES = TRUE)
  
  # kepp only the vecs with > 1 GO ID
  
  significant_GO_IDs_final <- significant_GO_IDS[which(significant_GO_IDs_count > 1)] 
  
  # kepp only dfs with > 1 GO ID
  
  ORA_results_filtered_final <- ORA_results_filtered[which(significant_GO_IDs_count > 1)]
  
  # get categories to which GO IDs in significant_GO_IDs_final belong to 
  
  ontology_vec <- c("BP", "CC", "MF")[which(significant_GO_IDs_count > 1)]
  
  # make similarity matrices 

  similarity_matrix_all_cats <- mapply(go_id = significant_GO_IDs_final, ontology_vec,
                                       MoreArgs = list(db = org.Mm.eg.db), FUN = GO_similarity, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  #cluster GO terms
  clustering_res <- sapply(X = similarity_matrix_all_cats, FUN = simplifyGO, plot = FALSE, simplify = FALSE, USE.NAMES = TRUE)
  
  # in reduced_ORA_results keep only the GO term with the lowest p.adjust from each cluster
  
  reduced_ORA_results <- mapply(x = clustering_res,
                                    y = ORA_results_filtered_final,
                                    FUN = function(x, y){
                                      x$p.adjust <- y$p.adjust
                                      x <- x %>% group_by(cluster) %>% dplyr::slice(which.min(p.adjust))
                                      y_min <- y %>% filter(gs_exact_source %in% x$id)
                                      y_min
                                    },
                                    SIMPLIFY = FALSE,
                                    USE.NAMES = TRUE)
  
  # if there is a df with 1 GO ID add it to the result
  if(length(which(significant_GO_IDs_count <= 1)) != 0){
    reduced_ORA_results_complete <- bind_rows(bind_rows(reduced_ORA_results), bind_rows(ORA_results_filtered[which(significant_GO_IDs_count <=1)]))
    return(reduced_ORA_results_complete)
  }else{
    return(bind_rows(reduced_ORA_results))
  }

}



allCells_FG_clustering_GO_ORA_no_geneID_vecs <- vector("list", 3)

for(i in 1:3){
  ORA_results_merged <- unlist(allCells_FG_clustering_GO_ORA[[i]], recursive = FALSE)
  just_ORA <- ORA_results_merged[str_which(names(ORA_results_merged), "ORA_results")]
  allCells_FG_clustering_GO_ORA_no_geneID_vecs[[i]] <- just_ORA
  names(allCells_FG_clustering_GO_ORA_no_geneID_vecs)[i] <- names(allCells_FG_clustering_GO_ORA)[i]
}


allCells_FG_clustering_GO_ORA_nonredundant.results <-sapply(X = allCells_FG_clustering_GO_ORA_no_geneID_vecs, 
                                                            FUN = function(x){
                                                            sapply(X = x, FUN = reduce_GO_enrichment_results, simplify = FALSE, USE.NAMES = TRUE)
                                                              }, simplify = FALSE, USE.NAMES = TRUE)



# edit ORA results

# convert specific words in gene_set_name to caps
case_convert <- function(x){

  x <- tolower(substring(x, 1))
  
  caps_words <- c("atp", "gtp", "rna", "snap", "5s", "dna", "amp", "nadh")
  
  position <- na.omit(as.data.frame(str_locate(x, pattern  = caps_words)))
  

  if(nrow(position) == 0){
    a <- toupper(substring(x, 1,1))
    b <- tolower(substring(x, 2))
    return(paste(a, b, sep = ""))
      
    
  } else{
    a <- toupper(substring(x, 1,1))
    b <- tolower(substring(x, 2))
    ab <- paste(a, b, sep = "")
      i <- 1
      while(i <= nrow(position)){
      ab <- str_replace_all(ab, pattern = substr(ab, position[i,1], position[i,2]), replacement = toupper(substr(ab, position[i,1], position[i,2])))
       i <- i+1 
      }
      
      return(ab)
      }
    }
      

# covert GeneRatio to numeric

make_GeneRatio_numeric <- function(x){
  y <- parse_number(str_split_fixed(x, pattern = "/", n =2))[1]/parse_number(str_split_fixed(x, pattern = "/", n =2))[2]
  y
}



# arrange gene_set_name into multiple lines

line_splitter <- function(x){
  
  y <- x
  x <- unlist(strsplit(x, split = " "))
  
  if (length(x) < 5){
    return(y)
  }else if(length(x) == 5) {
    paste0(x[1], " ", x[2], "\n",
           paste0(x[3:length(x)], sep = " ", collapse = "")
    )
  }else if(length(x) == 6) {
    paste0(x[1], " ", x[2], " ", x[3],"\n",
           paste0(x[4:length(x)], sep = " ", collapse = "")
    )
  }else if(length(x) == 7) {
    paste0(x[1], " ", x[2], " ", x[3],"\n",
           paste0(x[4:length(x)], sep = " ", collapse = "")
    )
  }else if(length(x) == 8) {
    paste0(x[1], " ", x[2], " ", x[3], " ", x[4],"\n",
           paste0(x[5:length(x)], sep = " ", collapse = "")
    )
  }else{
    paste0(x[1], " ", x[2], " ", x[3], " ", x[4], " ", x[5],"\n",
           paste0(x[6:length(x)], sep = " ", collapse = "")
    )
  }
}




allCells_FG_clustering_GO_ORA_nonredundant_results_final <- sapply(X = unlist(allCells_FG_clustering_GO_ORA_nonredundant_results, recursive = FALSE),
                                                               FUN = function(x){
                                                                 y <- x %>% 
                                                                      mutate(category = substr(Description, 4, 5)) %>%
                                                                      separate(Description, into = c(NA, NA, "gene_set_name"), sep = "_", remove = TRUE, extra = "merge")
                                                                 
                                                                 y$gene_set_name <- str_replace_all(y$gene_set_name, pattern = "_", replacement = " ")
                                                                 y$gene_set_name <- sapply(y$gene_set_name, FUN = case_convert)
                                                                 y$gene_set_name <- sapply(y$gene_set_name, FUN = line_splitter)
                                                                 y$GeneRatio_numeric <- sapply(y$GeneRatio, FUN = make_GeneRatio_numeric)
                                                                 y$category <- factor(y$category, levels = c("BP", "CC", "MF"))
                                                                 y <- y %>% group_by(category) %>% arrange(desc(GeneRatio_numeric), .by_group = TRUE )
                                                                 y
                                                               },
                                                        simplify = FALSE,
                                                        USE.NAMES = TRUE)


saveRDS(allCells_FG_clustering_GO_ORA_nonredundant_results_final, "output/allCells_FG_clustering_GO_ORA_nonredundant_results_final.rds", compress = TRUE)


###########################################################################
#       make edgelists and node tables for Gephi
##########################################################################


# ARGUMENTS
# clustering_results: output of get_cluster_membership_in_modules
# PPI: graph object
# min_vertices: the minimum # of vertices a connected component has to contain in order to be considered a part of a disease module
# min_cluster_size: the clusters with size < min_cluster_size are assigned cluster membership 0

STRING_mouse_700_removed_duplicates <- readRDS("./data/STRING_mouse_700_removed_duplicates_extended.rds")

get_final_node_and_edge_data <- function(clustering_results, PPI, min_vertices, min_cluster_size){
  
  all_module_genes <- clustering_results$Id
  
  #make a graph object with module vertices from the filtered PPI
  module_subgraph <- induced_subgraph(graph = PPI, vids = all_module_genes, impl = "auto")
  
  # decompose module igraph object to remove all single genes from module
  decomposed_graph <- decompose(module_subgraph, min.vertices = min_vertices)
  
  #merge components with more than min_vertices vertices in a single igraph
  decomposed_graph_vertices <- unlist(sapply(X = decomposed_graph, FUN = function(x){as_ids(V(x))}, simplify = TRUE))
  
  clustering_results_filtered<- clustering_results %>% dplyr::filter(Id %in% decomposed_graph_vertices)

  
  clustering_results_filtered$clust_membership_fg[which(clustering_results_filtered$clust_membership_fg %in% as.double(names(which(table(clustering_results_filtered$clust_membership_fg) < min_cluster_size))))] <- 0
  clustering_results_filtered$clust_membership_wt[which(clustering_results_filtered$clust_membership_wt %in% as.double(names(which(table(clustering_results_filtered$clust_membership_wt) < min_cluster_size))))] <- 0
  clustering_results_filtered$clust_membership_lp[which(clustering_results_filtered$clust_membership_lp %in% as.double(names(which(table(clustering_results_filtered$clust_membership_lp) < min_cluster_size))))] <- 0
  
  
  
  nodes <- clustering_results_filtered$Id
  
  edges <- STRING_mouse_700_removed_duplicates %>% filter(protein1 %in% nodes, protein2 %in% nodes) %>%
    rowwise() %>%
    dplyr::mutate(fg_clustering = ifelse(clustering_results_filtered[which(clustering_results_filtered$Id == protein1), "clust_membership_fg"] == clustering_results_filtered[which(clustering_results_filtered$Id == protein2), "clust_membership_fg"], 
                                         clustering_results_filtered[which(clustering_results_filtered$Id == protein1), "clust_membership_fg"], 0)) %>%
    
    dplyr::mutate(wt_clustering = ifelse(clustering_results_filtered[which(clustering_results_filtered$Id == protein1), "clust_membership_wt"] == clustering_results_filtered[which(clustering_results_filtered$Id == protein2), "clust_membership_wt"], 
                                         clustering_results_filtered[which(clustering_results_filtered$Id == protein1), "clust_membership_wt"], 0)) %>%
    
    dplyr::mutate(lp_clustering = ifelse(clustering_results_filtered[which(clustering_results_filtered$Id == protein1), "clust_membership_lp"] == clustering_results_filtered[which(clustering_results_filtered$Id == protein2), "clust_membership_lp"], 
                                         clustering_results_filtered[which(clustering_results_filtered$Id == protein1), "clust_membership_lp"], 0)) %>%
    dplyr::rename(Source = protein1) %>%
    dplyr::rename(Target = protein2)
  
  
  edges
  
  list(nodes = clustering_results_filtered,
       edges = edges)
  
}


allCells_all_node_and_edge_data <- mapply(clustering_results = all_cells_cluster_memberships,
                                          PPI = allCells_trimmed_networks[c(1:3)],
                                          min_vertices = c(132, 16, 286),
                                          min_cluster_size = rep(1,3),
                                          FUN = get_final_node_and_edge_data,
                                          SIMPLIFY = FALSE,
                                          USE.NAMES = TRUE)


#=======================================================================================
#                         visualize ORA results
#=======================================================================================


# select ORA results for clusters that are in the final module
final_modules_cluster_numbers <- vector("list", 3)

for(i in 1:3){
  x <- unique(allCells_all_node_and_edge_data[[i]]$nodes$clust_membership_fg)
  final_modules_cluster_numbers[[i]] <- paste(names(allCells_all_node_and_edge_data)[i], ".cluster_", x, ".ORA_results", sep = "")
  names(final_modules_cluster_numbers)[i] <- names(allCells_all_node_and_edge_data)[i]
}

vGluT2_ORA_results_for_plots <- allCells_FG_clustering_GO_ORA_nonredundant_results_final[final_modules_cluster_numbers$vGluT2] 
Gad2_ORA_results_for_plots <- allCells_FG_clustering_GO_ORA_nonredundant_results_final[final_modules_cluster_numbers$Gad2] 
Cx43_ORA_results_for_plots <- allCells_FG_clustering_GO_ORA_nonredundant_results_final[final_modules_cluster_numbers$Cx43] 


# select top 3 categories for each cluster:

get_top_GO <- function(x){
  y <- x %>% 
       group_by(category) %>% 
       top_n(3, -log10(p.adjust)) %>% 
       select(gene_set_name, category, GeneRatio_numeric, p.adjust)
  y
  
}

Cx43_top_ORA <- sapply(Cx43_ORA_results_for_plots, FUN = get_top_GO, simplify = FALSE, USE.NAMES = TRUE)
Gad2_top_ORA <- sapply(Gad2_ORA_results_for_plots, FUN = get_top_GO, simplify = FALSE, USE.NAMES = TRUE)
vGluT2_top_ORA <- sapply(vGluT2_ORA_results_for_plots, FUN = get_top_GO, simplify = FALSE, USE.NAMES = TRUE)

# make ORA plots for all clusters in disease modules

make_ORA_bar_plots <- function(ORA_df, cluster_color){


p <- ggplot(ORA_df, aes(x = GeneRatio_numeric, y = reorder(gene_set_name, GeneRatio_numeric))) +
  theme_gray() +

  geom_bar(stat = "identity",fill = cluster_color, color = "black", width = 0.98) +
  
  labs(y = NULL, x = "Gene ratio") +
  
  facet_grid(rows = vars(category), scales = "free_y", space = "free_y") +
  

  theme(strip.text = element_text(size = 8.9),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.x = element_text(size = 8.9, face = "bold"),
        axis.text.x  = element_text(size = 6.2),
        axis.text.y = element_text(size = 6.5),
        panel.background = element_rect(fill = "#E5E7E9", color = "#E5E7E9"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.12, "lines")
       
        )

p

}


vGluT2_cluster_colors <- c("#E9F28E", "#69C404", "#00C7FF", "#FFDF2A", "#FA9200", "#D3B3B0", "#EDB8FF", "#CDD9EF", "#FF5D85", "#54D6B1")
Gad2_cluster_colors <- c("#FF5D85", "#69C404", "#FA9200", "#00C7FF", "#EDB8FF", "#FFDF2A", "#D3B3B0")
Cx43_cluster_colors <- c("#69C404", "#FFDF2A", "#FF5D85", "#FA9200", "#D3B3B0", "#EDB8FF", "#00C7FF")


Cx43_clusters_ORA_plots <- mapply(ORA_df = Cx43_top_ORA, cluster_color = Cx43_cluster_colors, FUN = make_ORA_bar_plots, SIMPLIFY = FALSE, USE.NAMES = TRUE)
Gad2_clusters_ORA_plots <- mapply(ORA_df = Gad2_top_ORA, cluster_color = Gad2_cluster_colors, FUN = make_ORA_bar_plots, SIMPLIFY = FALSE, USE.NAMES = TRUE)
vGluT2_clusters_ORA_plots <- mapply(ORA_df = vGluT2_top_ORA, cluster_color = vGluT2_cluster_colors, FUN = make_ORA_bar_plots, SIMPLIFY = FALSE, USE.NAMES = TRUE)


#===========================================
# merge final node data with ORA results 
#===========================================

merge_node_table_with_ORA <- function(network_data, ORA_results){
  
dfs_list <- vector("list", length(ORA_results))

for(i in 1:length(ORA_results)){
  df <- data.frame(Id = unlist(str_split(ORA_results[[i]][1, 8], pattern = "/")))
  GO_term  <- as.character(ORA_results[[i]][1,2])
  df$GO_term <- rep(GO_term, nrow(df))
  
  dfs_list[[i]] <- df
  
}

one_df <- bind_rows(dfs_list)

complete_df <- merge(network_data$nodes, one_df, by = "Id", all.x = TRUE)

complete_df$GO_term <- replace_na(complete_df$GO_term, "not_in_top_GO")

list(nodes = complete_df, edges = network_data$edges)

}


allCells_all_node_and_edge_data_extended <- mapply(network_data = allCells_all_node_and_edge_data, 
                                                   ORA_results = list(vGluT2_ORA_results_for_plots, Gad2_ORA_results_for_plots, Cx43_ORA_results_for_plots),
                                                   FUN = merge_node_table_with_ORA,
                                                   SIMPLIFY = FALSE,
                                                   USE.NAMES = TRUE)


library(writexl)

for(i in 1:3){
  writexl::write_xlsx(allCells_all_node_and_edge_data_extended[[i]]$edges, paste(names(allCells_all_node_and_edge_data_extended)[i], "modules_edges_main", ".xlsx", sep = "" ))
}

for(i in 1:3){
  writexl::write_xlsx(allCells_all_node_and_edge_data_extended[[i]]$nodes, paste(names(allCells_all_node_and_edge_data_extended)[i], "modules_nodes_main", ".xlsx", sep = "" ))
}




