

get.dat.list <- function(dat.seu, assay.ls, red.name.ls) {
  if (!length(assay.ls) == length(red.name.ls)) {
    stop("Error: the length of assays is not equal to the length of the reduction names")
  }
  
  dat.ls <- c()
  for (i in 1:length(red.name.ls)) {
    assay.i <- assay.ls[i]
    red.name.i <- red.name.ls[i]
    dat.ls[[assay.i]] <- dat.seu[[red.name.i]]@cell.embeddings %>% t()
  }
  
  return(dat.ls)
}

# cell.type is a vector of true cell labels
# dnorm is a vector that stores the normalization factors, which can be calculated using the cal.dnorm function
cal.modality.weights <- function(dat.list, cell.type, q, dnorm, subsample=TRUE, n.subsample=200, seed=NULL) {
  # set seed if needed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # get constants
  ct <- unique(cell.type)
  n.ct <- length(ct)
  M <- length(dat.list)
  
  # initialize the weigth matrices
  mw <- list()
  for (m in 1 : M) {
    mw[[m]] <- matrix(NA, nrow=n.ct, ncol=n.ct)
    rownames(mw[[m]]) <- ct
    colnames(mw[[m]]) <- ct
  }
  
  names(mw) <- names(dat.list)
  
  for (i in 1 : (n.ct - 1)) {
    di <- pick.with.subsample(dat.list, cell.type, ct[i], subsample, n.subsample)
    for (j in (i + 1) : n.ct) {
      dj <- pick.with.subsample(dat.list, cell.type, ct[j], subsample, n.subsample)
      
      # distance matrix for individual modalities
      dmat <- list()
      for (m in 1 : M) {
        dmat[[m]] <- (1 - t(di[[m]]) %*% dj[[m]]) / dnorm[m] # dnorm is used for scaling!
      }
      
      # integrated distance
      dinter <- dmat[[1]] ^ q
      for (m in 2 : M) {
        dinter <- dinter + dmat[[m]] ^ q
      }
      dinter <- dinter ^ (1 / q)
      
      # weights
      v <- list()
      vsum <- matrix(0, nrow=ncol(di[[m]]), ncol=ncol(dj[[m]]))
      for (m in 1 : M) {
        v[[m]] <- (dmat[[m]] / dinter) ^ (q - 1)
        vsum <- vsum + v[[m]]
      }
      for (m in 1 : M) {
        mw[[m]][i, j] <- mean(v[[m]] / vsum)
        mw[[m]][j, i] <- mw[[m]][i, j]
      }
    }
  }
  
  return(mw)
}

# this is for internal use only
pick.cell.type.with.subsample <- function(dat, cell.type, cti, subsample, n.subsample) {
  # pick up data from cell type cti
  dm <- dat[, cell.type == cti]
  
  # subsample if needed
  if (subsample & ncol(dm) > n.subsample) {
    dm <- dm[, sample(1 : ncol(dm), n.subsample)]
  }
  
  # normalize to prepare for computing the cosine distance
  dm <- scale(dm, center=FALSE, scale=TRUE) / sqrt(nrow(dm) - 1)
  
  return(dm) 
}

# this is for internal use only
pick.with.subsample <- function(dat.list, cell.type, cti, subsample, n.subsample) {
  d <- list()
  for (m in 1 : length(dat.list)) {
    d[[m]] <- pick.cell.type.with.subsample(dat.list[[m]], cell.type, cti, subsample, n.subsample)
  }
  
  return(d)
}

# This function calculates the normalization factor using subsampling if necessary
# in dat.list, the rows are genes and columns are cells
cal.dnorm <- function(dat.list, subsample=TRUE, n.subsample=1000, seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  M <- length(dat.list)
  dnorm <- rep(0, M)
  
  # subsample if needed
  ind <- 1 : ncol(dat.list[[1]])
  if (subsample & length(ind) > n.subsample) {
    ind <- sample(ind, n.subsample)
  }
  
  for (m in 1 : M) {
    dat <- dat.list[[m]]
    dat <- dat[, ind]
    dat <- scale(dat, center=FALSE, scale=TRUE) / sqrt(nrow(dat) - 1)
    dnorm[m] <- mean(1 - t(dat) %*% dat)
  }
  
  return(dnorm)
}