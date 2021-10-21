# FIND TRUE POSITIVES  AND DEFINE THRESHOLDS


library(msigdbr)
library(tidyverse)


MSigDb_GO_KEGG_REACTOME_BIOCARTA<- bind_rows(msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% unite(gs_name, gs_subcat, gs_name, sep = "_", remove = TRUE),
                                             msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:CC") %>% unite(gs_name, gs_subcat, gs_name, sep = "_", remove = TRUE),
                                             msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF") %>% unite(gs_name, gs_subcat, gs_name, sep = "_", remove = TRUE),
                                             msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(-gs_subcat),
                                             msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:REACTOME") %>% dplyr::select(-gs_subcat),
                                             msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:BIOCARTA") %>% dplyr::select(-gs_subcat)) %>%
  dplyr::select(gs_name, entrez_gene, gs_exact_source)



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



