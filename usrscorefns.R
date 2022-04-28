
### These user defined score functions are for soft unknown interventions

### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam, usrpar = list(Imat = NULL, pctesttype = "bge", am = 1, weightvector = NULL, edgepmat = NULL, bgremove = TRUE)){

  n <- initparam$n
  Imat <- usrpar$Imat
  bgn <- ncol(Imat)
  colnames(Imat) <- paste0("i", 1:bgn)
  data <- initparam$data
  weightvector <- usrpar$weightvector
  
  exps <- mgcv::uniquecombs(Imat) # experimental conditions
  expsrows <- attr(exps, "index") # rows with each condition
 
  if (bgn > 0 && nrow(exps) > 1) {
    initparam <- scoreparameters(scoretype = "bge", data = cbind(Imat,data), weightvector = weightvector, bgnodes = 1:ncol(Imat), bgepar = list(am = am_value),
                                 edgepmat = usrpar$edgepmat)
    initparam$type <- "usr" # make sure it knows that we have redefined the score
    initparam$pctesttype <- "bge"
    initparam$bgremove <- usrpar$bgremove
    sigmas <- vector("list", nrow(exps))
    mus <- vector("list", nrow(exps))
    Ns <- vector("list", nrow(exps))
    
    for (ii in 1:nrow(exps)) {
      datalocal <- data[which(expsrows == ii), , drop = FALSE]
      if (is.null(weightvector)) {
        N <- nrow(datalocal)
        covmat <- cov(datalocal) * (N - 1)
        means <- colMeans(datalocal)
      }
      else {
        weightvectorlocal <- weightvector[which(expsrows == ii)]
        N <- sum(weightvectorlocal)
        forcov <- cov.wt(datalocal, wt = weightvectorlocal, cor = TRUE, 
                         method = "ML")
        covmat <- forcov$cov * N
        means <- forcov$center
      }
      sigmas[[ii]] <- covmat
      mus[[ii]] <- means
      Ns[[ii]] <- N
    }
    initparam$sigmas <- sigmas
    initparam$mus <- mus
    initparam$Ns <- Ns
  } else {
    initparam <- scoreparameters(scoretype = "bge", data = data, bgepar = list(am = am_value),
                                 edgepmat = usrpar$edgepmat)
  }
  
  initparam$exps <- exps
  initparam$expsrows <- expsrows
  initparam
}


### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j, parentnodes, n, param) {

  iparents <- parentnodes[which(parentnodes %in% param$bgnodes)]
  
  if (length(iparents) == 0 || nrow(param$exps) < 2) {# use standard BGe score
    localparam <- param
    localparam$type <- "bge"
    outscore <- BiDAG:::DAGcorescore(j, parentnodes, n, localparam)
  } else {
    parents <- setdiff(parentnodes, iparents)
    # find the different exp conditions for these parents
    local_exps <- attr(mgcv::uniquecombs(param$exps[, iparents, drop = FALSE]), "index")
    outscore <- 0
    for (ii in 1:max(local_exps)) {
      local_stats <- combinecovs(param$sigmas, param$mus, param$Ns, which(local_exps == ii), c(j, parents) - param$bgn)
      localparam <- BGeaugment(local_stats$sigma, local_stats$mu, local_stats$N, param$n, param$am, param$aw, param$logedgepmat)
      if (length(parents) > 0) {
        outscore <- outscore + BiDAG:::DAGcorescore(1, 1:length(parents) + 1, param$n, localparam)
      }  else {
        outscore <- outscore + BiDAG:::DAGcorescore(1, parents, param$n, localparam)
      }
    }
  }
  outscore
}
