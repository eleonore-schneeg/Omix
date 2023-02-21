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
  library(igraph)
  features_interest <- c(list[[1]], list[[2]])
  pk <- c(length(list[[1]]), length(list[[2]]))

  matrix1 <- t(multimodal_object[[1]])
  matrix1 <- matrix1[, colnames(matrix1) %in% list[[1]]]

  matrix2 <- t(multimodal_object[[2]])
  matrix2 <- matrix2[, colnames(matrix2) %in% list[[2]]]

  matrix <- scale(cbind(matrix1, matrix2))

  colnames(matrix)[1:pk[1]] <- paste0(colnames(matrix)[1:pk[1]], "_rna")
  colnames(matrix)[(pk[1] + 1):(pk[1] + pk[2])] <- paste0(colnames(matrix)[(pk[1] + 1):(pk[1] + pk[2])], "_protein")

  rho <- cor(matrix)
  thr_cor <- correlation_threshold
  rho <- ifelse(abs(rho) >= thr_cor, yes = rho, no = 0)
  A <- ifelse(rho != 0, yes = 1, no = 0)
  diag(A) <- 0


  a <- apply(rho, 1, function(x) {
    .color.gradient(x)
  })
  rownames(a) <- colnames(a)


  mynode_colours <- c(rep("#B39DDB", pk[1]), rep("#00796B", pk[2]))
  # mynode_labels <- sub("*\\.[0-9]", "", features_interest)
  mynode_labels <- colnames(matrix)
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

  list <- list(
    graph = mygraph,
    matrix = matrix
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
  multimodal_object <- multiassay@metadata$multimodal_object

  matrix1 <- t(multimodal_object[[1]])
  matrix1 <- matrix1[, colnames(matrix1) %in% list[[1]]]

  matrix2 <- t(multimodal_object[[2]])
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
.interactive_network <- function(igraph) {
  data1 <- visNetwork::toVisNetworkData(igraph)
  data1$nodes$font.color <- "black"
  nodes1 <- data1$nodes
  edges1 <- data1$edges

  visNetwork::visNetwork(nodes1,
    edges1,
    width = "100%", height = 800
  ) %>%
    visNetwork::visIgraphLayout(layout = "layout_with_fr") %>%
    visNetwork::visNodes(
      size = 1000,
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
      color = list(color = "#0085AF", highlight = "#C62F4B")
    ) %>%
    visNetwork::visLayout(randomSeed = 11)
}


.communities_network <- function(igraph) {
  lc1 <- igraph::cluster_louvain(igraph, weights = NULL, resolution = 1)
  communities <- igraph::communities(lc1)
  names <- attributes(communities)$dimnames[[1]]
  attributes(communities) <- NULL
  names(communities) <- names
  return(communities)
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
