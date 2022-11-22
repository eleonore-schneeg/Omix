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
#' @export
#'
#' @examples
integrate_with_DIABLO <- function(multimodal_omics,
                                  Y,
                                  ncomp,
                                  design=c('cor','full'),
                                  range = list(mRNA=seq(5,100,by=10),
                                               proteins=seq(5,100,by=10))) {
  multimodal_omics <- lapply(multimodal_omics, t)
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )


  A <- X$mRNA
  B <- X$proteins
  slot1 <- "mRNA"
  slot2 <- "proteins"

  if(design=='cor'){
  Apca <- prcomp(A, rank. = 1)
  Bpca <- prcomp(B, rank. = 1)
  cor <- cor(Apca$x, Bpca$x)

  design <- matrix(cor,
    ncol = length(X), nrow = length(X),
    dimnames = list(names(X), names(X))
  )
  diag(design) <- 0
  }

  if(design=='full'){
  design='full'
  }

  Y <- factor(Y)
  cli::cli_h2("MODEL TUNING")
  tune <- mixOmics::tune.block.splsda(
    X = X,
    Y = Y,
    ncomp = ncomp,
    design = design,
    test.keepX = list(
      "mRNA" = range[[1]],
      "proteins" = range[[2]]
    ),
    progressBar = TRUE
  )

  list.keepX <- tune$choice.keepX
  #cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
  tuned.diablo <- mixOmics::block.splsda(
    X = X,
    Y = Y,
    keepX = list.keepX,
    ncomp = ncomp,
    design = design,
    scale = T
  )

  model <- tuned.diablo
  X <- lapply(X, t)

  return(list(
    multimodal_object = X,
    model = model
  ))
}

integrate_with_sMBPLS <- function(multimodal_omics,
                                  Y,
                                  design=c('cor','full'),
                                  ncomp,
                                  list.keepX=list(mRNA = c(50), proteins = c(50))) {

  multimodal_omics <- lapply(multimodal_omics, t)
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )

  A <- multimodal_omics[[1]]
  B <- multimodal_omics[[2]]
  if(design=='cor'){
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
  }else{
    design='full'
  }
  #cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
  tuned.diablo <-  mixOmics::block.spls(
    X = X,
    Y = Y,
    keepX = list.keepX,
    ncomp = 1,
    design = design
  )
  model <- tuned.diablo
  X <- lapply(X, t)


  return(list(
    multimodal_object = X,
    model = model
  ))
}

integrate_with_MOFA <- function(multimodal_omics) {
  #cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )
  MOFAobject <- MOFA2::create_mofa(X)
  data_opts <- MOFA2::get_default_data_options(MOFAobject)
  data_opts$scale_views <- FALSE

  model_opts <- MOFA2::get_default_model_options(MOFAobject)
  model_opts$num_factors <- 10

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


integrate_with_iCluster <- function(){

  optk.brca <- getClustNum(data        = mo.data,
                           is.binary   = c(F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                           try.N.clust = 2:8, # try cluster number from 2 to 8
                           fig.name    = "CLUSTER NUMBER OF TCGA-BRCA")

  # perform iClusterBayes (may take a while)
  iClusterBayes.res <- getiClusterBayes(data        = mo.data,
                                        N.clust     = 5,
                                        type        = c("gaussian","gaussian","gaussian","binomial"),
                                        n.burnin    = 1800,
                                        n.draw      = 1200,
                                        prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                                        sdev        = 0.05,
                                        thin        = 3)
}
