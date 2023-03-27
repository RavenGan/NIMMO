


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


NIMMO <- function(dat.seu, assay.ls, red.name.ls,
                  q = 2, sampling.rate = 0.1, n.core) {
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
  reticulate::source_python("./R.fcn.0327.2023/nimmo_dist.py")

  config = umap.defaults
  config$n_neighbors = 30
  config$min_dist = 0.3
  config$metric = nimmo # the function in the .py file
  config$metric_kwds = list(dims = dims.mat, scales = scales.mat)
  for(x in c("n_epochs", "n_neighbors", "n_components", "random_state", "negative_sample_rate", "transform_state")) {
    config[[x]] = as.integer(config[[x]])
  }

  res <- umap::umap.learn(dat.all, config = config)
  res.umap <- res$layout
  colnames(res.umap) <- c("UMAP.1", "UMAP.2")
  res.umap <- Seurat::CreateDimReducObject(res.umap, key = "nimmo_", assay = "RNA")
  dat.seu[['nimmo']] <- res.umap
  return(dat.seu)
}
