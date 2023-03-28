utils::globalVariables(c("nimmo", "umap.defaults"))

#' The function to obtain the scaling factor and dimension of a data matrix
#'
#' @param dat a data matrix whose rows stand for samples and columns stand for PCs.
#' @param sample.idx the chosen sample indices which are used for the calculation of scaling factors
#' @return a vector of length 2. The first number is the scaling factor and the second one is the number of PCs.
#' @importFrom dplyr %>%

dist.fcn <- function(dat, sample.idx){
  dat.norm <- dat / sqrt(rowSums(dat * dat))
  dat.sub <- dat.norm[sample.idx, ]
  dat.cos <- dat.sub %*% t(dat.norm)

  scaling.factor <- mean(1 - dat.cos)
  dims <- ncol(dat)

  res <- c(scaling.factor, dims)
  names(res) <- c("sca.fac", "dims")
  return(res)
}

#' Perform batch integration of multimodal single-cell data using NIMMO.
#'
#' @param dat.seu a seurat object with multiple assays. Data preprocessing steps like normalization, high-variable feature selection and
#' dimension reduction are assumed to be performed.
#' @param assay.ls the names of assays used for integration. This information can be found using \code{Assays(dat.seu)}.
#' @param red.name.ls the names of dimension reductions. This information can be found using \code{names(dat.seu@reductions)}. Note that \code{assay.ls}
#' and \code{red.name.ls} should match.
#' @param q the parameter which decides the Lq norm. The default is \code{q = 2} which means using L2 norm to integrate data.
#' @param sampling.rate the proportion of cells that are used to calculate scaling factors. The default is \code{sampling.rate = 0.1}.
#' @param n.core the number of cores used for the calculation. The default is \code{n.core = 4}
#' @return a seurat object \code{dat} with one new object called \code{nimmo} under \code{dat@reductions}. The new object is a two-dimensional data matrix
#' obtained by UMAP dimension reduction.
#' @importFrom dplyr %>%
#' @export

NIMMO <- function(dat.seu, assay.ls, red.name.ls,
                  q = 2, sampling.rate = 0.1, n.core = 4) {
  if (!length(assay.ls) == length(red.name.ls)) {
    stop("Error: the length of assays is not equal to the length of the reduction names")
  }

  dat.ls <- c()
  for (i in 1:length(red.name.ls)) {
    assay.i <- assay.ls[i]
    red.name.i <- red.name.ls[i]
    dat.ls[[assay.i]] <- dat.seu[[red.name.i]]@cell.embeddings
  }

  dat.ls.1 <- dat.ls[[1]]
  n <- nrow(dat.ls.1)

  sample.idx <- sample(n, round(n * sampling.rate))

  fx <- function(dat.idx){
    dat <- dat.ls[[dat.idx]]
    dat.res <- dist.fcn(dat, sample.idx)
    return(dat.res)
  }

  dat.idx.ls <- c(1:length(dat.ls))
  scale.res.ls <- parallel::mclapply(dat.idx.ls, fx, mc.cores = n.core)
  names(scale.res.ls) <- names(dat.ls)

  # summarize scaling factors and dimensions
  scales.mat <- matrix(rep(NA, length(scale.res.ls)), 1, length(scale.res.ls))
  dims.mat <- matrix(rep(NA, length(scale.res.ls)), 1, length(scale.res.ls))

  for (i in 1:length(scale.res.ls)) {
    scales.mat[i] <- scale.res.ls[[i]]["sca.fac"]
    dims.mat[i] <- scale.res.ls[[i]]["dims"]
  }

  # Put all data together
  dat.all <- c()
  for (i in 1:length(dat.ls)) {
    dat.all <- cbind(dat.all, dat.ls[[i]])
  }

  # run python code
  set.ev <- paste0("os.environ['NIMMO_EV'] = ", "'", q, "'")
  reticulate::py_run_string("import os")
  reticulate::py_run_string(set.ev)
  reticulate::source_python(system.file("python", "nimmo_dist.py", package = "NIMMO"))
  # reticulate::source_python("./inst/python/nimmo_dist.py")

  config = umap.defaults
  config$n_neighbors = 30
  config$min_dist = 0.3
  config$metric = nimmo # the function in the .py file
  config$metric_kwds = list(dims = dims.mat, scales = scales.mat)
  for(x in c("n_epochs", "n_neighbors", "n_components", "random_state", "negative_sample_rate", "transform_state")) {
    config[[x]] = as.integer(config[[x]])
  }

  res <- umap:::umap.learn(dat.all, config = config)
  res.umap <- res$layout
  colnames(res.umap) <- c("UMAP.1", "UMAP.2")
  res.umap <- Seurat::CreateDimReducObject(res.umap, key = "nimmo_", assay = "RNA")
  dat.seu[['nimmo']] <- res.umap
  return(dat.seu)
}
