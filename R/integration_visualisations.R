#' Get multiomics networks from multimodal object
#'
#' @param multimodal_object
#' @param list
#' @param correlation_threshold
#'
#' @return
#' @export
#'
#' @examples
.multiomics_network_matrix <- function(multimodal_object,
                                       list,
                                       correlation_threshold = 0.5) {

  features_interest <- c(list[[1]], list[[2]])
  pk <- c(length(list[[1]]), length(list[[2]]))

  matrix1 <- t(multimodal_object[[1]])
  matrix1 <- matrix1[, list[[1]]]


  matrix2 <- t(multimodal_object[[2]])
  matrix2 <- matrix2[, list[[2]]]

  matrix <- scale(cbind(matrix1, matrix2))

  colnames(matrix)[1:pk[1]] <- paste0(colnames(matrix)[1:pk[1]], "_rna")
  colnames(matrix)[(pk[1] + 1):(pk[1] + pk[2])] <- paste0(colnames(matrix)[(pk[1] + 1):(pk[1] + pk[2])], "_protein")

  rho <- stats::cor(matrix)
  thr_cor <- correlation_threshold
  rho <- ifelse(abs(rho) >= thr_cor, yes = rho, no = 0)
  #A <- ifelse(rho != 0, yes = 1, no = 0)
  #diag(A) <- 0

  rho[rho  < correlation_threshold ] <- 0
  diag(rho) <- 0
  mygraph <- igraph::graph.adjacency(rho, weighted=TRUE, mode="lower")


  a <- apply(rho, 1, function(x) {
    .color.gradient(x)
  })
  rownames(a) <- colnames(a)


  mynode_colours <- c(rep("#B39DDB", pk[1]), rep("#00796B", pk[2]))
  # mynode_labels <- sub("*\\.[0-9]", "", features_interest)
  mynode_labels <- colnames(matrix)
  #mygraph <- .GetGraph(adjacency = A, node_color = mynode_colours, node_label = mynode_labels)

  # data <- igraph::as_edgelist(mygraph)
  # V1 <- data[, 1]
  # V2 <- data[, 2]
  # b <- numeric()
  #
  # for (i in 1:length(V1)) {
  #   b[i] <- a[which(rownames(a) == V1[i]), which(colnames(a) == V2[i])]
  # }

  igraph::V(mygraph)$label <- mynode_labels
  igraph::V(mygraph)$color <-  mynode_colours
  igraph::V(mygraph)$label.family <- "sans"
  igraph::V(mygraph)$label.color <- "grey20"
  igraph::E(mygraph)$weight
  igraph::V(mygraph)$label.cex=0.5
  #igraph::E(mygraph)$color <- b
  set.seed(1)
  l <- layout.fruchterman.reingold(mygraph)
  #plot=plot(mygraph,layout=l,vertex.size=(hub_score(mygraph)$vector*20)+1)

  list <- list(
    graph = mygraph,
    matrix = matrix,
    #plot=plot,
    hubs=hub_score(mygraph)$vector
  )
  return(list)
}


#' Get multiomics networks from multiassay object

#'
#' @param multiassay
#' @param list
#' @param correlation_threshold
#'
#' @return
#' @export
#'
#' @examples
.multiomics_network <- function(multiassay,
                                list,
                                correlation_threshold = 0.5) {
  library(igraph)
  features_interest <- c(list[[1]], list[[2]])
  pk <- c(length(list[[1]]), length(list[[2]]))

  multimodal_object <- multiassay@metadata$multimodal_object


  matrix1 <- t(multimodal_object[[1]])
  matrix1 <- matrix1[, list[[1]]]


  matrix2 <- t(multimodal_object[[2]])
  matrix2 <- matrix2[, list[[2]]]

  matrix <- scale(cbind(matrix1, matrix2))

  colnames(matrix)[1:pk[1]] <- paste0(colnames(matrix)[1:pk[1]], "_rna")
  colnames(matrix)[(pk[1] + 1):(pk[1] + pk[2])] <- paste0(colnames(matrix)[(pk[1] + 1):(pk[1] + pk[2])], "_protein")

  rho <- stats::cor(matrix)
  thr_cor <- correlation_threshold
  rho <- ifelse(abs(rho) >= thr_cor, yes = rho, no = 0)
  #A <- ifelse(rho != 0, yes = 1, no = 0)
  #diag(A) <- 0

  rho[rho  < correlation_threshold ] <- 0
  diag(rho) <- 0
  mygraph <- igraph::graph.adjacency(rho, weighted=TRUE, mode="lower")


  a <- apply(rho, 1, function(x) {
    .color.gradient(x)
  })
  rownames(a) <- colnames(a)


  mynode_colours <- c(rep("#B39DDB", pk[1]), rep("#00796B", pk[2]))
  # mynode_labels <- sub("*\\.[0-9]", "", features_interest)
  mynode_labels <- colnames(matrix)
  #mygraph <- .GetGraph(adjacency = A, node_color = mynode_colours, node_label = mynode_labels)

  # data <- igraph::as_edgelist(mygraph)
  # V1 <- data[, 1]
  # V2 <- data[, 2]
  # b <- numeric()
  #
  # for (i in 1:length(V1)) {
  #   b[i] <- a[which(rownames(a) == V1[i]), which(colnames(a) == V2[i])]
  # }

  igraph::V(mygraph)$label <- mynode_labels
  igraph::V(mygraph)$color <-  mynode_colours
  igraph::V(mygraph)$label.family <- "sans"
  igraph::V(mygraph)$label.color <- "grey20"
  igraph::E(mygraph)$weight
  igraph::V(mygraph)$label.cex=0.5
  #igraph::E(mygraph)$color <- b
  set.seed(1)
  l <- layout.fruchterman.reingold(mygraph)
  #plot=plot(mygraph,layout=l,vertex.size=(hub_score(mygraph)$vector*20)+1)

  list <- list(
    graph = mygraph,
    matrix = matrix,
    #plot=plot,
    hubs=hub_score(mygraph)$vector
  )
  return(list)
}

.color.gradient <- function(x, colors = c("blue", "white", "red"), colsteps = 100) {
  return(colorRampPalette(colors)(colsteps)[findInterval(x, seq(-1, 1, length.out = colsteps))])
}

.GetGraph <- function(calib_object = NULL, adjacency = NULL,
                      node_label = NULL, node_color = NULL, node_shape = NULL,
                      weighted = NULL, satellites = FALSE) {
  # either out or adjacency have to be provided

  if (is.null(adjacency)) {
    if (is.null(calib_object)) {
      stop("Either 'calib_object' or 'adjacency' needs to be provided.")
    }
    adjacency <- CalibratedAdjacency(calib_object)
  }

  if (is.null(node_color)) {
    node_color <- rep("skyblue", ncol(adjacency))
  }

  if (is.null(node_shape)) {
    node_shape <- rep("circle", ncol(adjacency))
  }

  if (is.null(node_label)) {
    node_label <- colnames(adjacency)
  }

  names(node_color) <- colnames(adjacency)
  names(node_label) <- colnames(adjacency)
  names(node_shape) <- colnames(adjacency)

  mygraph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = weighted)
  V(mygraph)$label <- node_label[V(mygraph)$name]
  V(mygraph)$color <- node_color[V(mygraph)$name]
  V(mygraph)$shape <- node_shape[V(mygraph)$name]
  V(mygraph)$frame.color <- V(mygraph)$color
  V(mygraph)$label.family <- "sans"
  E(mygraph)$color <- "grey60"
  V(mygraph)$label.color <- "grey20"
  E(mygraph)$width <- 0.5

  return(mygraph)
}

#' Get interactive network from igraph object

#'
#' @param igraph
#'
#' @return
#' @export
#'
#' @examples
.interactive_network <- function(igraph,
                                 communities,
                                 cluster) {



  data1 <- visNetwork::toVisNetworkData(igraph)
  data1$nodes$font.color <- "black"
  nodes1 <- data1$nodes
  edges1 <- data1$edges
  nodes1$size <- (hub_score(igraph)$vector*30)+1

  #Create group column
  if(communities==TRUE){
    cluster=communities$community_object
    cluster_df <- data.frame(as.list(membership(cluster)))
    cluster_df <- as.data.frame(t(cluster_df))
    cluster_df$label <- rownames(cluster_df)


  nodes1 <- merge(x = nodes1, y = cluster_df, by = "label", all.x = TRUE)
  colnames(nodes1)[8] <- "group"
  }

  visNetwork::visNetwork(nodes1,
    edges1,
    width = "100%", height = 800
  ) %>%
    visNetwork::visIgraphLayout(layout = "layout.fruchterman.reingold") %>%
    visNetwork::visNodes(
      color = list(
        background = "#0085AF",
        border = "#013848",
        highlight = "#FF8000",
        labelHighlightBold = T
      ),
      scaling = list(
        max = 1,
        min = 1,
        label = list(
          enabled = F,
          max = 1,
          min = 1
        )
      ),
      shadow = list(enabled = TRUE, size = 10)
    ) %>%
    visNetwork::visEdges(
      shadow = FALSE,
      color = list(color = "grey", highlight = "#C62F4B")
    ) %>%
    visNetwork::visLayout(randomSeed = 11)
}


.communities_network <- function(igraph,
                                 community_detection='louvain') {

  if(community_detection=='edge_betweeness'){
  lc1 <- igraph::cluster_edge_betweenness(igraph, weights = NULL)
  }
  if(community_detection=='leading_eigen'){
    lc1 <- igraph::cluster_leading_eigen(igraph, weights = NULL)
  }

  if(community_detection=='walktrap'){
    lc1 <- igraph::cluster_walktrap(igraph, weights = NULL)
  }

  if(community_detection=='fastgreedy'){
    lc1 <- igraph::cluster_fast_greedy(igraph, weights = NULL)
  }
  if(community_detection=='louvain'){
    lc1 <- igraph::cluster_louvain(igraph, weights = NULL)
  }
  communities <- igraph::communities(lc1)
  names <- attributes(communities)$dimnames[[1]]
  attributes(communities) <- NULL
  names(communities) <- names

  layout <-layout.fruchterman.reingold(igraph)
  plot(lc1, igraph,layout=layout,
       vertex.label=NA, vertex.size=5, edge.arrow.size=.2)

  return(list(community_object=lc1,
              communities=communities))
}


.multiomics_network_cluster <- function(multiassay,
                                        integration = "iCluster",
                                        cluster = 1,
                                        list,
                                        correlation_threshold = 0.5) {
  library(igraph)
  features_interest <- c(list[[1]], list[[2]])
  multimodal_object <- multiassay@metadata$multimodal_object

  clusters <- multiassay@metadata$integration[[paste(integration)]]$clust.res$clust
  keep <- clusters == cluster

  matrix1 <- t(multimodal_object[[1]])
  matrix1 <- matrix1[keep, ]
  matrix1 <- matrix1[, colnames(matrix1) %in% list[[1]]]

  matrix2 <- t(multimodal_object[[2]])
  matrix2 <- matrix2[keep, ]
  matrix2 <- matrix2[, colnames(matrix2) %in% list[[2]]]

  matrix <- cbind(matrix1, matrix2)

  rho <- cor(matrix)
  thr_cor <- correlation_threshold
  rho <- ifelse(abs(rho) >= thr_cor, yes = rho, no = 0)
  A <- ifelse(rho != 0, yes = 1, no = 0)
  diag(A) <- 0


  a <- apply(rho, 1, function(x) {
    .color.gradient(x)
  })
  rownames(a) <- colnames(a)

  pk <- c(length(list[[1]]), length(list[[2]]))
  mynode_colours <- c(rep("#B39DDB", pk[1]), rep("#00796B", pk[2]))
  mynode_labels <- sub("*\\.[0-9]", "", features_interest)

  mygraph <- .GetGraph(adjacency = A, node_color = mynode_colours, node_label = mynode_labels)

  data <- as_edgelist(mygraph)
  V1 <- data[, 1]
  V2 <- data[, 2]
  b <- numeric()

  for (i in 1:length(V1)) {
    b[i] <- a[which(rownames(a) == V1[i]), which(colnames(a) == V2[i])]
  }

  V(mygraph)$name <- mynode_labels
  E(mygraph)$color <- b
  E(mygraph)$width <- 3
  set.seed(1)

  return(mygraph)
}

.color.gradient <- function(x, colors = c("blue", "white", "red"), colsteps = 100) {
  return(colorRampPalette(colors)(colsteps)[findInterval(x, seq(-1, 1, length.out = colsteps))])
}

.community_graph<- function(igraph,
                            community_object,
                            community=1){

  keep = which(community_object$membership == community)

  ## Which nodes should be kept?
  Keep = V(igraph)[keep]
  igraph_sub = igraph::induced_subgraph(igraph, Keep,impl ='create_from_scratch')
  E(igraph_sub)$width <- 0.5
  layout=igraph::layout.fruchterman.reingold(igraph_sub)
  #plot=plot(igraph_sub, vertex.size=(hub_score(igraph_sub)$vector*20)+1, layout=layout)

   res=list(
    graph = igraph_sub,
    #plot = plot,
    hubs=hub_score(igraph_sub)$vector)

return(res)
}
