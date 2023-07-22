#' @name getStdiz
#' @title Get standardized omics data
#' @description This function prepare standardized data for generating heatmap. Omics data, especially for expression, should be centered or scaled or z-scored (both centered and scaled). Generally, DNA methylation beta matrix and somatic mutation (0 and 1 binary matrix) should not be normalized. This function also provides an argument of `halfwidth` 
#' for continuous omics data; such argument is used to truncate the 'extremum' after normalization; specifically, normalized values that exceed the halfwidth boundaries will be replaced by the halfwidth, which is vary useful to map colors in heatmap. Function adapted from package `MOVICS`  (see ref)
#' @param data A list of data.frame or matrix storing raw multiple omics data with rows for features and columns for samples.
#' @param halfwidth A numeric vector to assign marginal cutoff for truncating values in data; 1 by default.
#' @param centerFlag A logical vector to indicate if each subdata should be centered; TRUE by default.
#' @param scaleFlag A logical vector to indicate if each subdata should be scaled; TRUE by default.
#' @export
#' @family Plotting
#' @return A standardized data.frame containing multi-omics data.
#' @references Lu, X., et al. (2020). MOVICS: an R package for multi-omics integration and visualization in cancer subtyping. Bioinformatics, 36(22-23), 5539–5541. 

getStdiz <- function(data       = NULL,
                     halfwidth  = rep(1, length(data)),
                     centerFlag = rep(TRUE, length(data)),
                     scaleFlag  = rep(TRUE, length(data))) {

  # check data
  if(is.null(names(data))){
    names(data) <- sprintf("dat%s", 1:length(data))
  }

  n_dat <- length(data)

  outdata <- list()
  for (i in 1:n_dat) {
    tmp <- t(scale(t(data[[i]]), center = centerFlag[i], scale = scaleFlag[i]))
    if (!is.na(halfwidth[i])) {
      tmp[tmp > halfwidth[i]] = halfwidth[i]
      tmp[tmp < (-halfwidth[i])] = -halfwidth[i]
    }
    outdata[[names(data)[i]]] <- tmp
  }

  return(outdata)
}


#' Multi-omics heatmap
#'
#' @param multimodal_omics multimodal_omics list
#' @param metadata metdata dataframe
#' @param integration_model Possible integration methods are `MOFA`,
#' `DIABLO`,`sMBPLS`,`iCluster`,`MEIFESTO`
#' @param covariates clinical covariates
#' @param omic `rna`,`protein` or `both`
#' @param weights weights from integrative results
#' @param de_features If integration model is clustering, use differentially expressed features
#' @param direction `positive` or `negative`
#' @param cluster default to NULL
#' @param num_features default to NULL will plot all the selected features. Otherwise numerical variable to define the maximum number of features
#'
#' @return Heatmap
#' @family Plotting
#' @export

multiomics_heatmap <- function(multimodal_omics,
                               metadata,
                               integration_model = "iCluster",
                               covariates=c('PHF1','diagnosis'),
                               omic=c('rna','protein'),
                               weights,
                               de_features,
                               direction=c("both",'positive','negative'),
                               cluster=NULL,
                               num_features=NULL){

  data =  multimodal_omics

  plot_data <- getStdiz(
    data = data,
    halfwidth = c(2, 2),
    centerFlag = c(T, T),
    scaleFlag = c(T, T)
  )

  if(integration_model!='iCluster'){
  pos=weights$weights$ranked_weights_positive
  neg=weights$weights$ranked_weights_negative
  }

  if(integration_model=='iCluster'){
    pos=de_features[[cluster]]$Up
    neg=de_features[[cluster]]$Down
  }


  if(!is.null( num_features)){
    pos=lapply(pos, function(x){ x[1:num_features]})
    neg=lapply(neg, function(x){ x[1:num_features]})
  }

  if(!is.null( num_features)){
    pos=lapply(pos, function(x){ x[!is.na(x)]})
    neg=lapply(neg, function(x){ x[!is.na(x)]})
  }

  if(omic=='rna'){

  data=plot_data$rna_processed

    if(direction=='both'){
      annotation_row = data.frame(direction=c(rep('pos',length(pos$rna)),rep('neg',length(neg$rna))))
      rownames(annotation_row )= rownames(data[c(pos$rna,neg$rna),])
      plot= pheatmap::pheatmap(data[c(pos$rna,neg$rna),],
           annotation_col = metadata[,covariates],
            annotation_row =annotation_row)
    }
    if(direction=='positive'){
      plot=pheatmap::pheatmap(data[c(pos$rna),],
               annotation_col = metadata[,covariates])
    }
    if(direction=='negative'){
      plot=pheatmap::pheatmap(data[c(neg$rna),],
               annotation_col = metadata[,covariates])
    }
  }

  if(omic=='protein'){

  data=plot_data$protein_processed

    if(direction=='both'){
      annotation_row = data.frame(direction=c(rep('pos',length(pos$protein)),rep('neg',length(neg$protein))))
      rownames(annotation_row )= rownames(data[c(pos$protein,neg$protein),])
      plot=pheatmap::pheatmap(data[c(pos$protein,neg$protein),],
               annotation_col = metadata[,covariates],
              annotation_row =annotation_row
              )
    }
    if(direction=='positive'){
      plot=pheatmap::pheatmap(data[c(pos$protein),],
               annotation_col = metadata[,covariates])
    }
    if(direction=='negative'){
      plot=pheatmap::pheatmap(data[c(neg$protein),],
               annotation_col = metadata[,covariates])
    }
  }

  if(omic=='both'){

    data=plot_data$rna_processed
    data2=plot_data$protein_processed

    if(direction=='both'){
      plot_d=rbind(data[c(pos$rna,neg$rna),],data2[c(pos$protein,neg$protein),])
      annotation_row = data.frame(omic=
                                  c(rep('rna',length(c(pos$rna,neg$rna))),
                                  rep('protein',length(c(pos$protein,neg$protein)))))
      rownames(annotation_row )= rownames(plot_d)
      plot=pheatmap::pheatmap( plot_d,
                               annotation_col = metadata[,covariates],
                               annotation_row =annotation_row, cluster_rows = F,
                               fontsize_row = 5)


    }
    if(direction=='positive'){

      plot_d=rbind(data[c(pos$rna),],data2[c(pos$protein),])
      annotation_row = data.frame(omic=c(rep('rna',length(pos$rna)),
                                         rep('protein',length(pos$protein))))
      rownames(annotation_row )= rownames(plot_d)
      plot=pheatmap::pheatmap( plot_d,
                              annotation_col = metadata[,covariates],
                              annotation_row =annotation_row, cluster_rows = F,
                              fontsize_row = 5)

    }
    if(direction=='negative'){
      plot_d=rbind(data[c(neg$rna),],data2[c(neg$protein),])
      annotation_row = data.frame(omic=c(rep('rna',length(neg$rna)),
                                         rep('protein',length(neg$protein))))
      rownames(annotation_row )= rownames(plot_d)
      plot=pheatmap::pheatmap( plot_d,
                               annotation_col = metadata[,covariates],
                               annotation_row =annotation_row, cluster_rows = F,
                               fontsize_row = 5)
    }

  }
plot
}
