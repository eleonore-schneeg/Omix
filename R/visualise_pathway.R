#' visualise_pathway
#'
#' @param pathways dataframe of enrichment results
#' @param gene_data LogFC
#' @param num_path number of pathways to plot
#' @family Plotting
#' @import RColorBrewer
#' @import visNetwork
#' @import dplyr
#' @return viridis
#' @export
#'
visualise_pathway <- function(pathways,gene_data=NULL,num_path=10){
  
  pathways <- pathways %>%
    arrange(FDR) %>%       # Sort the dataframe by the FDR column in ascending order
    slice_head(n = num_path) 
  
  color_palette <- colorRampPalette(rev(brewer.pal(9, "Blues")))(100)
  fdr_colors <- rev(viridis(100, option = "D")[as.numeric(cut(pathways$FDR, breaks = 100))])
  
  # Initialize vectors for edges
  gene_names <- c()
  pathway_ids <- c()
  
  # Loop through each row and extract genes
  for (i in 1:nrow(pathways)) {
    genes <- unlist(strsplit(as.character(pathways$genes[i]), ";"))
    gene_names <- c(gene_names, genes)
    pathway_ids <- c(pathway_ids, rep(pathways$description[i], length(genes)))
  }
  
  # Create an edge dataframe
  edges <- data.frame(from = pathway_ids, to = gene_names)
  
  # Create a nodes dataframe with labels and types
  pathway_labels <- setNames(pathways$description, pathways$geneset)
  names(pathway_labels) <-pathways$description
  gene_labels <- unique(gene_names)
  node_labels <- c(pathway_labels, setNames(gene_labels, gene_labels))
  nodes <- data.frame(id = names(node_labels), label = node_labels)
  
  # Distinguish between pathway and gene nodes
  nodes$type <- ifelse(nodes$id %in% names(pathway_labels), "Pathway", "Gene")
  
  # Add FDR colors to pathway nodes
  pathway_fdr_colors <- setNames(fdr_colors, names(pathway_labels))
  
  nodes$type <- ifelse(nodes$id %in% names(pathway_labels), "Pathway", "Gene")
  
  nodes$color <- ifelse(nodes$type == "Pathway", pathway_fdr_colors[nodes$id], NA)
  edges$color <- 'lightgrey'
  # Customize node sizes and label font sizes
  nodes$size <- ifelse(nodes$type == "Pathway", 30, 10)
  nodes$font.size <- ifelse(nodes$type == "Pathway", 30, 20)
  nodes$font.color <- ifelse(nodes$type == "Pathway", "black", "black")
  nodes$font.style <- ifelse(nodes$type == "Pathway", "bold", "normal")
  
  insertLineBreaks <- function(label, words_per_line = 2) {
    # Split the label into words
    words <- unlist(strsplit(label, " "))
    
    # Initialize an empty string for the new label
    new_label <- ""
    
    # Iterate over the words and insert line breaks
    for (i in seq_along(words)) {
      new_label <- paste0(new_label, words[i], 
                          ifelse(i %% words_per_line == 0 && i < length(words), "\n", " "))
    }
    
    return(trimws(new_label))
  }
  
  # Example usage
  new_label <- insertLineBreaks("This is an example label with more than three words", 3)
  
  # Apply the function to the label column
  nodes$label <- sapply(nodes$label, insertLineBreaks)
  
  # Merge the gene_data with the nodes dataframe
  if(!is.null(gene_data)){
    nodes <- merge(nodes, gene_data, by.x = "id", by.y = "gene", all.x = TRUE)
  }
  
  # Scale the size of gene nodes based on avg_log2FC
  # You may need to adjust the scaling factor to get a desirable range of sizes
  max_size <- 40  # Maximum size for any node
  min_size <- 5  # Minimum size for gene nodes
  scale_factor <- 25  # Adjust this factor to scale sizes
  
  # Handling missing or NA avg_log2FC values
  if(!is.null(gene_data)){
    nodes$avg_log2FC[is.na(nodes$avg_log2FC)] <- min(nodes$avg_log2FC, na.rm = TRUE)
    
    
    # Scale node sizes
    nodes$size <- ifelse(nodes$type == "Pathway", max_size, 
                         min_size + scale_factor * abs(nodes$avg_log2FC))
    
  }
  # Generate Viridis color scale
  fdr_scale <- seq(0, 1, length.out = 100)
  fdr_colors <- rev(viridis(100, option = "D"))
  
  # Create a dataframe for the legend
  legend_df <- data.frame(FDR = fdr_scale, Color = fdr_colors)
  
  # Create a legend plot using ggplot2
  legend=ggplot(legend_df, aes(x = FDR, y = 1, fill = Color)) +
    geom_tile() +
    scale_fill_identity() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10)) +
    labs(x = "FDR", y = "")
  
  
  ##
  library(igraph)
  network_igraph <- graph_from_data_frame(d = edges, vertices = nodes)
  # Convert to undirected graph
  network_igraph_undirected <- as.undirected(network_igraph, mode = "collapse")
  
  communities <- cluster_fast_greedy(network_igraph_undirected)
  
  nodes$group <- communities$membership
  
  # Plot the network with customized node properties
  visNetwork(nodes, edges, width = "100%",height = "800px") %>% 
    visNodes(size = "size", font = list(size = "font.size" )) %>%  visIgraphLayout()%>% 
    visPhysics(
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(gravitationalConstant = 0, centralGravity = 0.01,
                              springLength = 20, springConstant = 0.08)
    )%>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
               selectedBy = "group")
  
  
  
}