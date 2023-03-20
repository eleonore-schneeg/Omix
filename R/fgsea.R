#' Performs functional gene set enrichment analysis
#'
#' @param gene_list
#' @param pval
#'
#' @return
#' @export
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr filter
#'
#' @examples
GSEA = function(gene_list, pval) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)

  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }

  cgp_gene_sets =msigdbr::msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
  msigdbr_list = split(x = cgp_gene_sets$gene_symbol, f = cgp_gene_sets$gs_name)
  names(msigdbr_list)=sub('.*GOBP_', '', names(msigdbr_list))
  names(msigdbr_list)=gsub("_", " ", names(msigdbr_list))
  names(msigdbr_list)=tolower(names(msigdbr_list))
  myGO=msigdbr_list

  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        maxSize=400, ## maximum gene set size
                        nperm=10000) %>%
    as.data.frame() %>%
    dplyr::filter(padj < !!pval) %>%
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))

  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))

  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")

  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))

  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header)

  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

#' Plots FGSEA results
#'
#' @param gs_results
#' @param myGO
#' @param min.sz
#' @param main
#'
#' @return
#' @export
#'
#' @examples
plot_geneset_clusters = function( gs_results, myGO, min.sz = 4, main="GSEA clusters"){
  library(ggplot2)
  library(ggrepel)
  library(stringr)


  df = matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
  rownames(df) = colnames(df) = gs_results$pathway

  for ( i in 1:nrow(gs_results)) {
    genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
    for (j in 1:nrow(gs_results)) {
      genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
      ## Jaccards distance  1 - (intersection / union )
      overlap = sum(!is.na(match(genesI, genesJ )))
      jaccards = overlap / length(unique(c(genesI, genesJ) ))
      df[i,j] = 1-jaccards
    }
  }

  ## Cluster nodes using dynamic tree cut, for colors
  distMat = as.dist(df)
  dendro = hclust(distMat, method = "average" )
  clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = min.sz )
  ## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it

  gs_results$Cluster = clust
  gs_results = gs_results[gs_results$Cluster != 0, ]

  ## select gene sets to label for each clusters
  bests = gs_results %>%
    group_by( Cluster ) %>%
    top_n(wt = abs(size), n = 10) %>%
    .$pathway
  ## determine cluster order for plotting
  clust_ords = gs_results %>%
    group_by( Cluster ) %>%
    summarise("Average" = NES ) %>%
    arrange(desc(Average)) %>%
    .$Cluster %>%
    unique

  gs_results$Cluster = factor(gs_results$Cluster, levels = clust_ords)

  gs_results$Label = ""
  gs_results$Label[gs_results$pathway %in% bests ] = gs_results$pathway[gs_results$pathway %in% bests ]
  gs_results$Label = str_remove(gs_results$Label, "GO_")
  gs_results$Label = tolower(gs_results$Label)

  g1 = ggplot(gs_results, aes(x = Cluster, y = NES, label = Label )) +
    geom_jitter( aes(color = Cluster,  size = size), alpha = 0.8, height = 0, width = 0.2 ) +
    scale_size_continuous(range = c(0.5,5)) +
    geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5)+
    ggtitle(main)

  return(g1)
}
