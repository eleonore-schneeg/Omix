#' Creates a circos plot of identified TF and gene targets in modules 
#'
#' @param TF output of `Transcription_Factor_enrichment()`
#'
#' @return circos plot
#' @family Multi-omic integration downstream analysis
#' @import RColorBrewer
#' @import circlize
#' @importFrom purrr is_empty

#' @export
#'

circos_TF <- function(TF){
  

suppressWarnings({
TF_list=names(TF@go.nested.list)
list_genes=list()
list_genes_OR=list()
modules=names(TF@go.nested.list[[1]])

for(j in TF_list){
for(i in modules){
  list_genes[[j]][[i]]=TF@go.nested.list[[paste(j)]][[paste(i)]]@intersection
  
}
}

for(i in names(list_genes)){
  plot=(lapply(list_genes[[i]], function(x){purrr::is_empty(x)}))
}

if(all(plot)==FALSE){
for(j in TF_list){
for(i in modules){
  list_genes_OR[[j]][[i]]=TF@go.nested.list[[paste(j)]][[paste(i)]]@odds.ratio
  
}
}


list_enrichment=list_genes_OR

# Get the unique genes and transcription factors
unique_genes <- unique(unlist(lapply(list_genes, unlist)))
unique_tfs <-  names(list_enrichment)

# Create an empty matrix with appropriate dimensions
link_matrix <- NULL
link_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(unique_tfs))
colnames(link_matrix) <- unique_tfs
# Populate the matrix with enrichment scores

for (i in 1:length(list_genes)) {
  for (j in 1:length(list_genes[[i]])) {
    genes <- list_genes[[i]][[j]]
    enrichment <- list_enrichment[[i]][[j]]
    
    gene_indices <- match(genes, unique_genes)
    tf_indices <- match(unique_tfs[i], colnames(link_matrix))
    
    link_matrix[gene_indices, tf_indices] <- enrichment
  }
}

# Assign row and column names to the link_matrix
row.names(link_matrix) <- unique_genes

unique_genes <- unique(unlist(unlist(list_genes, recursive = FALSE), use.names = FALSE))
gene_module <- sapply(unique_genes, function(gene) {
  module <- sapply(list_genes, function(tf_genes) {
    modules <- names(tf_genes)[sapply(tf_genes, function(module_genes) gene %in% module_genes)]
    if (length(modules) > 0) modules else NA
  })
  module <- na.omit(module)
  if (length(module) > 1) {
    warning(paste("Gene", gene, "belongs to multiple modules:", paste(module, collapse = ", ")))
  }
  setNames(module[1], sub(".*\\.", "", gene))
})

# Create the vector of modules named by corresponding unique gene
gene_module_vector <- unlist(gene_module)
names(gene_module_vector)=sub(".*\\.", "", names(gene_module_vector))

tf_vector <- setNames(paste0("TF", seq_along(unique_tfs)), unique_tfs)




num_modules <- length(unique( c(gene_module_vector,tf_vector)))
color_palette <- colorRampPalette(brewer.pal(9, "Spectral"))
module_colors <- color_palette(num_modules)

#module_colors <- rainbow(length(unique( c(gene_module_vector,tf_vector))))
module_color_vector <- setNames(module_colors, unique(c(gene_module_vector,tf_vector)))

vector_all=c(gene_module_vector,tf_vector)
for(i in names(module_color_vector)){
  vector_all[vector_all==i ]=module_color_vector[i]
}

# Create a chord diagram with gene grouping
par(cex = 0.8)
chordDiagram(
  t(link_matrix),
  group = c(gene_module_vector,tf_vector),
  grid.col =vector_all,
  preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(link_matrix))))),
  annotationTrack = "grid")

# Customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

for(i in unique(gene_module_vector)){
  highlight.sector(names(gene_module_vector)[gene_module_vector==i], track.index = 1, col = "lightgrey",text = i, cex = 0.8, text.col = "black", niceFacing = TRUE,padding = c(-0.9, 0, 0.27, 0),font=2) 
}
}else{cat('No target genes found in modules!')}
})

}
