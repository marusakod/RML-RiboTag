########################################################################
#  Determination of disease module thresholds and threshold validation
########################################################################

library(msigdbr)
library(tidyverse)
library(igraph)
library(clusterProfiler)


# function that counts the number of genes with p-value <= threshold at each iteration 

count_low_pval_genes <- function(diamond_results, threshold, celltype){
  
  df <- diamond_results$added_genes %>%
    mutate(low_pval = ifelse(DE_pvalue <= threshold, TRUE, FALSE)) %>% 
    mutate(low_pval_count = replace(low_pval, which(low_pval == TRUE), c(1 : length(which(low_pval == TRUE))))) 
  
  df$true_low_pval_count <- zeros_converter(added_genes = df, logical_column = "low_pval", count_column = "low_pval_count")
  
  df$iteration_no <- 1:nrow(diamond_results$added_genes)
  
  df$cell <- celltype
  df
}


zeros_converter <- function(added_genes, logical_column, count_column){
  
  final_vec <- c() 
  for(i in 1:nrow(added_genes)){
    if(added_genes[i, logical_column] == FALSE){
      vec <- c(1:i)
      vec_2 <- added_genes[1:i, count_column]
      final <- max(vec_2)
      final_vec[i] <- final
      
    }else{
      final <- added_genes[i, count_column]
      final_vec[i] <- final
      
    } 
  }
  
  return(final_vec)
}



modified_diamond_results_alpha1 <- readRDS("output/modified_diamond_results/low_pval_diamond_alpha_1_iterations_500.rds")


allCells_low_pvalue_count_dfs <- mapply(diamond_results = modified_diamond_results_alpha1,
                                        celltype = c("vGluT2", "Gad2", "Cx43"),
                                        threshold = 0.01,
                                        FUN = count_low_pval_genes,
                                        SIMPLIFY = FALSE,
                                        USE.NAMES = TRUE)

allCells_low_pvalue_count_df_combined <- bind_rows(allCells_low_pvalue_count_dfs)
allCells_low_pvalue_count_df_combined$cell <- factor(allCells_low_pvalue_count_df_combined$cell, levels = c("Cx43", "vGluT2", "Gad2"))




thresholds_df <- data.frame(cell = c("vGluT2", "Gad2", "Cx43"),
                            threshold = c(120, 100, 215))


threshold_text_df <- data.frame(cell = c("vGluT2", "Gad2", "Cx43"),  
                                         x_text = c(275,275,350),
                                         y_text = c(70,40,115),
                                         # p-values from all_cells_final_modules_enrich_signif
                                         lab = c("threshold: 120\np-value: 2.12 x 10-6", "threshold: 100\np-value: 1.71 x 10-2", "threshold: 215\np-value: 2.26 x 10-21"))




low_pval_count_plot  <- ggplot(allCells_low_pvalue_count_df_combined, aes(x = iteration_no, y = true_low_pval_count)) +
  theme_bw() +
  geom_point(size = 0.5)  +
  facet_wrap(~ cell, scale = "free_y", ncol = 1, strip.position = "right") +
  geom_vline(data = thresholds_df, aes(xintercept = threshold), linetype ="dashed", size = 0.5, color = "#808B96") +
  
  labs(x= "Iteration", y = "Number of genes with Wald p-value < 0.01") +
  
  geom_text(threshold_text_df, inherit.aes = FALSE, 
  mapping = aes(x = x_text, y = y_text, label = lab), size = 4) +
  
  
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "#D7DBDD", color = "#D7DBDD"),
        panel.background = element_rect(fill = "#F2F3F4", color = "#F2F3F4"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.4, "lines"),
        strip.text = element_text(size =14)) 






#########################################################################################################
# threshold validation: determine the significance of the number of true positives in disease modules
#########################################################################################################


# source all KEGG, REACTOME. BIOCARTA and GO terms from MsigDb

MSigDb_GO_KEGG_REACTOME_BIOCARTA<- bind_rows(msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% unite(gs_name, gs_subcat, gs_name, sep = "_", remove = TRUE),
                                             msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:CC") %>% unite(gs_name, gs_subcat, gs_name, sep = "_", remove = TRUE),
                                             msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF") %>% unite(gs_name, gs_subcat, gs_name, sep = "_", remove = TRUE),
                                             msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(-gs_subcat),
                                             msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:REACTOME") %>% dplyr::select(-gs_subcat),
                                             msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:BIOCARTA") %>% dplyr::select(-gs_subcat)) %>%
  dplyr::select(gs_name, entrez_gene, gs_exact_source)



# function for hypergeometric test

supersets_enricher <- function(genes, allgenes, padjust, set, maxSS) {
  
  db <- MSigDb_GO_KEGG_REACTOME_BIOCARTA %>% dplyr::select(gs_name, entrez_gene) %>%  filter(str_detect(gs_name, paste("^", set, sep = "")))
  results <- enricher(gene = genes, 
                      pvalueCutoff = padjust,
                      pAdjustMethod = "BH",
                      universe = allgenes,
                      qvalueCutoff = 0.2,
                      minGSSize = 10,
                      maxGSSize = maxSS,
                      TERM2GENE = db)@result %>%
    filter(p.adjust <= padjust)
  return(results)
  
}



# Get list of vectors with seed genes ---> in module_genes_clustering_and_ORA.R:

source("R/module_genes_clustering_and_ORA.R")


# Get list of vectors with background genes (all network genes) :
universe_list <- vector("list", 3)


for(i in 1:length(Comparison_DE_vec_10W[1:3])){
  
  vertices <- as_ids(V(allCells_trimmed_networks[[i]]))
 
  universe_list[[i]] <- vertices
  names(universe_list)[i] <- gsub("10Weeks__", "", Comparison_DE_vec_10W[i])
  
}

# Identify GO terms and PATHWAYS significantly enriched within the given set of seed genes:

seeds_GO_PATHWAYS_enrichr_results_p0.05_maxSS_500 <- mapply(genes = DEGs_list, allgenes = universe_list,
                                                            padjust = 0.05, set = "GO|KEGG|REACTOME|BIOCARTA", maxSS = 500,  FUN = supersets_enricher, SIMPLIFY = FALSE, USE.NAMES = TRUE)

# For each DIAMOnD gene check whether it is a true positive (annotated with any of terms that are enriched within the seed gene set). 

diamondGenes_function_check <- function(diamond_result, seed_enrich_result, celltype, set){
  iterations <- nrow(diamond_result$added_genes)
  
  db <- MSigDb_GO_KEGG_REACTOME_BIOCARTA %>% filter(str_detect(gs_name, paste("^", set, sep = "")))
  
  
  sig_terms <- as.vector(seed_enrich_result$ID)
  filtered_supersets <- db %>% filter(gs_name %in% sig_terms)
  all_sig_genes <- as.vector(filtered_supersets$entrez_gene)
  
  
  added_genes_df <- data.frame(gene = diamond_result$added_genes$gene) %>%
    mutate(cell = gsub("10Weeks__", "", celltype)) %>%
    mutate(iteration_no = c(1: iterations)) %>%
    mutate(True_positive  = ifelse(gene %in% all_sig_genes, TRUE, FALSE)) %>%
    mutate(TP_count = replace(True_positive, which(True_positive == TRUE), c(1 : length(which(True_positive == TRUE))))) 
  
  
  added_genes_df$true_TP_count <- zeros_converter(added_genes = added_genes_df, "True_positive", "TP_count")
  
  return(added_genes_df)
  
}


allCells_TP_count_dfs <- mapply(diamond_result = modified_diamond_results_alpha1,
                                seed_enrich_result = seeds_GO_PATHWAYS_enrichr_results_p0.05_maxSS_500,
                                celltype = Comparison_DE_vec_10W[1:3],
                                set = "GO|KEGG|REACTOME|BIOCARTA",
                                FUN = diamondGenes_function_check, 
                                SIMPLIFY = FALSE,
                                USE.NAMES = TRUE)


# Determine the significance of the number of true positives in disease modules


ckeck_functional_similarity_significance <- function(TP_result, threshold,  seed_enrich_result, set, PPI){
  
  db <- MSigDb_GO_KEGG_REACTOME_BIOCARTA %>% filter(str_detect(gs_name, paste("^", set, sep = "")))
  sig_terms <- as.vector(seed_enrich_result$ID)
  
  all_vertices <- as_ids(V(PPI))
  
  
# get the number of genes annotated to significant gene sets (enriched in seed genes) and in the network --> to be used as the number of all possible successes in Fischer's exact test
  
  filtered_supersets <- db %>% filter(gs_name %in% sig_terms) # filter db so that just the significant sets are left
  all_sig_genes <- as.vector(filtered_supersets$entrez_gene) # extract all genes annotated to significant gene sets 
  
  # all genes annotated with terms significantly enriched within the given set of seed genes and in the STRING network
  all_successes <- length(dplyr::intersect(all_sig_genes, all_vertices))
  
  # get the number of all ANNOTATED genes in the STRING network --> to be used as population size in Fischer's exact test
  population_size <- length(dplyr::intersect(all_vertices, unique(db$entrez_gene)))
  
  # get number of all annotated disease module genes:

  sample_size <- length(intersect(TP_result$gene[1:threshold], db$entrez_gene))
  
  
  successes <- TP_result[threshold, "true_TP_count"]
  
  p_val <- phyper(successes-1, all_successes, population_size - all_successes, sample_size, lower.tail = FALSE)
  
  p_val
  
}



all_cells_final_modules_enrich_signif <- mapply(TP_result = allCells_TP_count_dfs,
                                                threshold = c(120, 100, 215),
                                                seed_enrich_result = seeds_GO_PATHWAYS_enrichr_results_p0.05_maxSS_500,
                                                set = "GO|KEGG|REACTOME|BIOCARTA" ,
                                                PPI = allCells_trimmed_networks[c(1:3)],
                                                FUN = ckeck_functional_similarity_significance,
                                                SIMPLIFY = TRUE,
                                                USE.NAMES = TRUE)













