
### These are helper functions for soft unknown interventions
# this one combines covariance matrices

combinecovs <- function(sigmas, mus, Ns, iis, jjs) {
  sigma <- matrix(0, length(jjs), length(jjs))
  mu <- rep(0, length(jjs))
  N <- 0
  for (ii in iis) {
    sigma <- sigma + sigmas[[ii]][jjs, jjs] + Ns[[ii]]*mus[[ii]][jjs]%*%t(mus[[ii]][jjs])
    mu <- mu + Ns[[ii]]*mus[[ii]][jjs]
    N <- N + Ns[[ii]]
  }
  sigma <- sigma - mu%*%t(mu)/N # mu is not rescaled yet!
  out <- list()
  out$sigma <- sigma
  out$mu <- mu/N
  out$N <- N
  return(out)
}

# this one turns them into statistics needed for the BGe score

BGeaugment <- function(sigma, mu, N, n, am, aw, logedgepmat) {
  mu0 <- numeric(length(mu))
  T0scale <- am * (aw - n - 1)/(am + 1)
  T0 <- diag(T0scale, length(mu), length(mu))
  TN <- T0 + sigma + ((am * N)/(am + N)) * (mu0 - mu) %*% t(mu0 - mu)
  awpN <- aw + N
  constscorefact <- -(N/2) * log(pi) + (1/2) * log(am/(am + N))
  muN <- (N * mu + am * mu0)/(N + am)
  SigmaN <- TN/(awpN - n - 1)
  scoreconstvec <- numeric(length(mu))
  for (j in (1:length(mu))) {
    awp <- aw - n + j
    scoreconstvec[j] <- constscorefact - lgamma(awp/2) + 
      lgamma((awp + N)/2) + ((awp + j - 1)/2) * log(T0scale)
  }
  localparam <- list()
  localparam$type <- "bge"
  localparam$TN <- TN
  localparam$awpN <- awpN
  localparam$scoreconstvec <- scoreconstvec
  localparam$DBN <- FALSE
  localparam$MDAG <- FALSE
  localparam$logedgepmat <- logedgepmat
  return(localparam)
}


### This function gives edges weights between the bounds
# with both positive and negative signs

wFUN <- function(m, lb, ub){ # function for edge weights
  runif(m, lb, ub)*sample(c(-1, 1), m, replace = TRUE)
}


### This function generates Gaussian data from a DAG
# following the topological order
### Intervened data points are set to random values
# without worrying about the parents

irmvDAG <- function(trueDAGedges, Tmat, shift = NULL, damp = NULL, scale = TRUE) {
  trueDAG <- 1*(trueDAGedges != 0) # the edge presence in the DAG
  n <- ncol(trueDAG) # number of variables
  if (n != ncol(Tmat)) {
    stop("The DAG and intervention data don't have the same number of variables.")
  }
  N <- nrow(Tmat) # number of observations
  data <- matrix(0, nrow = N, ncol = n) # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy) # go down order
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1) # find parents
    lp <- length(parents) # number of parents
    if (lp == 0) { # no parents
      data[, jj] <- 0
    } else if (lp == 1) { # one parent
      data[, jj] <- data[, parents]*trueDAGedges[parents, jj]
    } else { # more than one parent
      data[, jj] <- colSums(t(data[, parents])*trueDAGedges[parents, jj])
    }
    int_obs <- which(Tmat[, jj] == 1) # which observations were intervened upon
    if (length(int_obs) > 0) { 
      if (is.null(shift) && is.null(damp)) { # intervene and set to 0
        data[int_obs, jj] <- 0
      } else {
        if (!is.null(damp)) { # dampen contribution from parents
          data[int_obs, jj] <- data[int_obs, jj]*damp_vec[int_obs]
        }
        if (!is.null(shift)) { # shift means
          data[int_obs, jj] <- data[int_obs, jj] + shift_vec[int_obs]
        }
      }
    }
    # add random noise
    data[, jj] <- data[, jj] + rnorm(N)
  }
  if (scale) {
    scale(data)
  } else {
    data
  }
}


### This function turns a DAG into an EG
DAG2EG <- function(incidence, targets) {
  as(dag2essgraph(as(incidence, "graphNEL"), targets), "matrix")
}

### This function extracts directed edges from an EG
EGdedges <- function(incidence) {
  incidence*(1 - t(incidence))
}

### This function extracts the skeleton from an EG
EGskel <- function(incidence) {
  1*(incidence|t(incidence))
}


### This function compares an estimated EG to the true one

compareEGs <- function (estEG, trueEG) {
  estSkel <- EGskel(estEG) # estimated skeleton
  trueSkel <- EGskel(trueEG) # true skeleton
  P <- sum(trueSkel)/2 # number of positives
  diffSkel <- estSkel - trueSkel
  extra_edges <- which(diffSkel > 0) # edges in estimated but not true EG
  FP <- length(extra_edges)/2 # count to FPs
  estEG[extra_edges] <- 0 # remove them from further comparisons
  missing_edges <- which(diffSkel < 0) # edges in true but not estimated EG
  FN <- length(missing_edges)/2 # count to FNs
  trueEG[missing_edges] <- 0 # remove them from further comparisons
  # modified graphs have the same skeletons, so now just need to count mismatches
  mismatches <- 1*(estEG != trueEG)
  wrong_order <- sum(EGskel(mismatches))/2 # number of wrongly oriented edges
  FP <- FP + wrong_order/2 # include half in FP
  FN <- FN + wrong_order/2 # and half in FN
  SHD <- FP + FN # shd is the sum of errors
  TP <- P - FN # true positives are without false negatives
  # TPR, FPR_P, FPR_N
  if (P == 0) { # true graph is empty
    if (FP >= 0) {
      TPR <- 0
      FPR_P <- 1
    } else {
      TPR <- 1
      FPR_P <- 0
    }
  } else { # true graph is non-empty
    TPR <- TP/P
    FPR_P <- FP/P
  }
  compEGs <- c(TP, FP, SHD, TPR, FPR_P)
  names(compEGs) <- c("TP","FP", "SHD", "TPR", "FPR_P")
  return(compEGs)
}


### This function takes in a chain of DAGs,
# converts then to EGs and computes the posterior of the 
# skeleton and directed edges 

EGedgePosterior <- function (MCMCchain, burnin = 0.2, targets) {
  endstep <- length(MCMCchain)
  startstep <- max(as.integer(burnin*endstep), 1)
  EGs <- lapply(MCMCchain[startstep:endstep], DAG2EG, targets)
  skels <- lapply(EGs, EGskel)
  dedges <- lapply(EGs, EGdedges)
  EG.post <- Reduce("+", EGs)/(endstep - startstep + 1)
  skel.post <- Reduce("+", skels)/(endstep - startstep + 1)
  dedges.post <- Reduce("+", dedges)/(endstep - startstep + 1)
  return(list(EG = EG.post, skel = skel.post, dedges = dedges.post))
}


### This function takes in the skeleton and edge posteriors and returns a 
# consensus above a threshold

consensusEG <- function (posteriors, p_threshold) {
  #skel <- 1*(posteriors$skel > p_threshold)
  #dedges <- 1*(posteriors$dedges > p_threshold)
  #skel - t(dedges) # remove the backward version of the edges from the skeleton
  EG <- 1*(posteriors$EG > p_threshold)
  EG
}


### This function takes the string output of UT-IGSP and turns it into 
# an adjacency matrix

SP_to_DAG <- function(x) {
  if(x == "error"){
    NULL
  } else {
    y <- strsplit(x, "\\[|\\]")[[1]]
    y <- y[which(y != "")]
    n <- length(y)
    DAG <- matrix(0, n, n)
    for (ii in 1:n) {
      z <- as.numeric(strsplit(y[ii], ",|\\|")[[1]]) + 1 # remember offset
      if(length(z) > 1) {
        DAG[z[-1], z[1]] <- 1
      }
    }
    DAG
  }
}
