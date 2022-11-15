integrate_with_DIABLO <- function(multimodal_omics,
                                  ncomp,
                                  range) {
  multimodal_omics <- lapply(multimodal_omics, t)
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )


  A <- X$mRNA
  B <- X$proteins
  slot1 <- "mRNA"
  slot2 <- "proteins"

  Apca <- prcomp(A, rank. = 1)
  Bpca <- prcomp(B, rank. = 1)
  cor <- cor(Apca$x, Bpca$x)

  design <- matrix(cor,
    ncol = length(X), nrow = length(X),
    dimnames = list(names(X), names(X))
  )
  diag(design) <- 0
  Y <- factor(Y)
  cli::cli_h2("MODEL TUNING")
  tune <- mixOmics::tune.block.splsda(
    X = X,
    Y = Y,
    ncomp = ncomp,
    design = design,
    test.keepX = list(
      "mRNA" = range,
      "proteins" = range
    ),
    progressBar = TRUE
  )

  list.keepX <- tune$choice.keepX
  cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
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
                                  ncomp,
                                  range) {
  multimodal_omics <- lapply(multimodal_omics, t)
  X <- list(
    mRNA = multimodal_omics[[1]],
    proteins = multimodal_omics[[2]]
  )

  A <- multimodal_omics[[1]]
  B <- multimodal_omics[[2]]

  Apca <- prcomp(A, rank. = 1)
  Bpca <- prcomp(B, rank. = 1)
  cor <- cor(Apca$x, Bpca$x)

  design <- matrix(cor,
    ncol = length(X), nrow = length(X),
    dimnames = list(names(X), names(X))
  )
  diag(design) <- 0
  Y <- factor(Y)
  list.keepX <- list(mRNA = c(50), proteins = c(50))

  design <- matrix(cor,
    ncol = length(X), nrow = length(X),
    dimnames = list(names(X), names(X))
  )
  diag(design) <- 0
  cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
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
  cli::cli_h2("VERTICAL INTEGRATION IN PROCESS")
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
  metadata$sample <- metadata[[paste(map_by_column)]]
  MOFA2::samples_metadata(MOFAobject) <- metadata
  model <- MOFAobject


  return(list(
    multimodal_object = X,
    model = model
  ))
}
