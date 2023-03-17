#' Vertical integration with DIABLO including model tuning
#'
#' @param multimodal_omics
#' @param ncomp
#' @param design Design matrix (design = "full"): The strength of all
#' relationships between dataframes is maximised (= 1) – a “fully connected” design.
#' If design is set on cor, the correlation between PC1 of each dataset will be
#' set for the design matrix

#' @param range List of the range of numbers of features to keep in the tuning phase.
#' First element must be for rna, second for proteins
#'
#' @return
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment listToMap colData
#'  getWithColData sampleMap
#' @importFrom mixOmics block.splsda tune.block.splsda
#' @export
#'
#' @examples
integrate_with_DIABLO <- function(multimodal_omics,
                                  Y,
                                  ncomp,
                                  design = c("cor", "full"),
                                  range = list(
                                    mRNA = seq(5, 10, by = 10),
                                    proteins = seq(5, 10, by = 10)
                                  )) {
  multimodal_omics <- lapply(multimodal_omics, t)
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )


  A <- X$mRNA
  B <- X$proteins
  slot1 <- "mRNA"
  slot2 <- "proteins"

  if (design == "cor") {
    Apca <- prcomp(A, rank. = 1)
    Bpca <- prcomp(B, rank. = 1)
    cor <- cor(Apca$x, Bpca$x)

    design <- matrix(cor,
      ncol = length(X), nrow = length(X),
      dimnames = list(names(X), names(X))
    )
    diag(design) <- 0
  }

  if (design == "full") {
    design <- "full"
  }


  Y <- factor(Y)
  cli::cli_h2("MODEL TUNING")
  tune <- suppressWarnings({mixOmics::tune.block.splsda(
    X = X,
    Y = Y,
    ncomp = ncomp,
    design = design,
    test.keepX = list(
      "mRNA" = range[[1]],
      "proteins" = range[[2]]
    ),
    validation = "Mfold",
    folds = 5,
    nrepeat = 1,
    progressBar = TRUE
  )
  })

  list.keepX <- tune$choice.keepX
  tuned.diablo <- suppressWarnings({mixOmics::block.splsda(
    X = X,
    Y = Y,
    keepX = list.keepX,
    ncomp = ncomp,
    design = design,
    scale = T
  )})

  model <- tuned.diablo
  X <- lapply(X, t)

  return(list(
    multimodal_object = X,
    model = model
  ))
}

#' Vertical integration with SMBPLS
#'
#' @param multimodal_omics
#' @param Y
#' @param design
#' @param ncomp
#' @param list.keepX
#'
#' @return
#' @importFrom mixOmics block.spls
#' @export
#'
#' @examples
integrate_with_sMBPLS <- function(multimodal_omics,
                                  Y,
                                  design = c("cor", "full", "avg"),
                                  ncomp,
                                  list.keepX = list(mRNA = c(50), proteins = c(50))) {
  multimodal_omics <- lapply(multimodal_omics, t)
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )

  A <- multimodal_omics[[1]]
  B <- multimodal_omics[[2]]
  if (design == "cor") {
    Apca <- prcomp(A, rank. = 1)
    Bpca <- prcomp(B, rank. = 1)
    cor <- cor(Apca$x, Bpca$x)

    design <- matrix(cor,
      ncol = length(X), nrow = length(X),
      dimnames = list(names(X), names(X))
    )
    diag(design) <- 0

    design <- matrix(cor,
      ncol = length(X), nrow = length(X),
      dimnames = list(names(X), names(X))
    )
    diag(design) <- 0

    if (design == "avg") {
      design <- matrix(0.3,
        ncol = length(X), nrow = length(X),
        dimnames = list(names(X), names(X))
      )
      diag(design) <- 0
    }
  } else {
    design <- "full"
  }
  tuned.diablo <- mixOmics::block.spls(
    X = X,
    Y = Y,
    keepX = list.keepX,
    ncomp = ncomp,
    design = design
  )
  model <- tuned.diablo
  X <- lapply(X, t)


  return(list(
    multimodal_object = X,
    model = model
  ))
}


#' Vertical integration with MBPLS
#'
#' @param multimodal_omics
#' @param Y
#' @param design
#' @param ncomp
#'
#' @return
#' @importFrom mixOmics block.pls
#' @export
#'
#' @examples
integrate_with_MBPLS <- function(multimodal_omics,
                                 Y,
                                 design = c("cor", "full", "avg"),
                                 ncomp) {
  multimodal_omics <- lapply(multimodal_omics, t)
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )

  A <- multimodal_omics[[1]]
  B <- multimodal_omics[[2]]
  if (design == "cor") {
    Apca <- prcomp(A, rank. = 1)
    Bpca <- prcomp(B, rank. = 1)
    cor <- cor(Apca$x, Bpca$x)

    design <- matrix(cor,
      ncol = length(X), nrow = length(X),
      dimnames = list(names(X), names(X))
    )
    diag(design) <- 0

    design <- matrix(cor,
      ncol = length(X), nrow = length(X),
      dimnames = list(names(X), names(X))
    )
    diag(design) <- 0

    if (design == "avg") {
      design <- matrix(0.3,
        ncol = length(X), nrow = length(X),
        dimnames = list(names(X), names(X))
      )
      diag(design) <- 0
    }
  } else {
    design <- "full"
  }
  tuned.diablo <- mixOmics::block.pls(
    X = X,
    Y = Y,
    ncomp = ncomp,
    design = design
  )
  model <- tuned.diablo
  X <- lapply(X, t)


  return(list(
    multimodal_object = X,
    model = model
  ))
}


#' Vertical integration with MOFA
#'
#' @param multimodal_omics
#' @param num_factors
#' @param scale_views
#'
#' @return
#' @export
#' @importFrom MOFA2 create_mofa get_default_data_options
#' get_default_model_options get_default_training_options prepare_mofa run_mofa
#'
#' @examples
integrate_with_MOFA <- function(multimodal_omics,
                                num_factors = 5,
                                scale_views = T,
                                metadata) {
  # cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )

  MOFAobject <- MOFA2::create_mofa(X)
  data_opts <- MOFA2::get_default_data_options(MOFAobject)
  data_opts$scale_views <- scale_views

  model_opts <- MOFA2::get_default_model_options(MOFAobject)
  model_opts$num_factors <- num_factors

  train_opts <- MOFA2::get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- "medium"
  train_opts$seed <- 42

  MOFAobject <- MOFA2::prepare_mofa(MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )

  MOFAobject <- MOFA2::run_mofa(MOFAobject)
  metadata$sample <- rownames(metadata)
  MOFA2::samples_metadata(MOFAobject) <- metadata
  model <- MOFAobject


  return(list(
    multimodal_object = X,
    model = model
  ))
}


#' Vertical integration with MEIFESTO
#'
#' @param multimodal_omics
#' @param num_factors
#' @param scale_views
#' @param metadata
#' @param time
#'
#' @return
#' @export
#'
#' @examples
integrate_with_MEIFESTO <- function(multimodal_omics,
                                    num_factors = 5,
                                    scale_views = T,
                                    metadata,
                                    time = "pseudotime") {
  time <- metadata[, time]
  names(time) <- rownames(metadata)
  time <- data.frame(time)
  time <- t(time)

  # cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )

  MOFAobject <- MOFA2::create_mofa(X)
  MOFAobject <- MOFA2::set_covariates(MOFAobject, covariates = time)
  data_opts <- MOFA2::get_default_data_options(MOFAobject)
  data_opts$scale_views <- scale_views

  model_opts <- MOFA2::get_default_model_options(MOFAobject)
  model_opts$num_factors <- num_factors


  train_opts <- MOFA2::get_default_training_options(MOFAobject)
  train_opts$seed <- 2020
  train_opts$maxiter <- 1000
  train_opts$convergence_mode <- "medium"

  mefisto_opts <- MOFA2::get_default_mefisto_options(MOFAobject)

  MOFAobject <- MOFA2::prepare_mofa(MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts,
    mefisto_options = mefisto_opts
  )

  MOFAobject <- MOFA2::run_mofa(MOFAobject)
  metadata$sample <- rownames(metadata)
  MOFA2::samples_metadata(MOFAobject) <- metadata
  model <- MOFAobject


  return(list(
    multimodal_object = X,
    model = model
  ))
}


#'  Vertical integration with iCluster
#'
#' @param multimodal_omics
#' @param try.N.clust
#'
#' @return

#' @importFrom MOVICS getClustNum getiClusterBayes
#' @export
#'
#' @examples
integrate_with_iCluster <- function(multimodal_omics,
                                    try.N.clust = 2:4) {
  cli::cli_alert_success("Optimising the number of cluster (this make take a while")

  optk.i <- MOVICS::getClustNum(
    data = multimodal_omics,
    is.binary = c(F, F),
    try.N.clust = try.N.clust,
    fig.name = "Cluster number tuning"
  )


  optimal_n <- optk.i$N.clust

  cli::cli_alert_success("Integration in progress (this make take a while)")
  model <- MOVICS::getiClusterBayes(
    data = multimodal_omics,
    N.clust = optimal_n,
    type = c("gaussian", "gaussian"),
    n.burnin = 1800,
    n.draw = 1200,
    prior.gamma = c(0.5, 0.5),
    sdev = 0.05,
    thin = 3
  )

  return(list(
    multimodal_object = multimodal_omics,
    model = model
  ))
}
