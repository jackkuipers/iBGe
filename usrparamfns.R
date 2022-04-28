
### This user defined parameter function is estimating hard causal effects
# for data with soft unknown interventions

usrDAGparametersCore <- function(j, parentNodes, param) {
  # check which parents are interventional
  iparents <- parentNodes[which(parentNodes %in% param$bgnodes)]

  if (length(iparents) == 0 || nrow(param$exps) < 2) {# use standard BGe score
    localparam <- param
    localparam$type <- "bge"
    outscore <- Bestie:::DAGparametersCore(j, parentNodes, localparam)
  } else {
    parents <- setdiff(parentNodes, iparents)
    # find the different exp conditions for these parents
    exp_conds <- mgcv::uniquecombs(param$exps[, iparents, drop = FALSE])
    local_exps <- attr(exp_conds, "index")
    ii <- which(rowSums(exp_conds) == 0) # this defines the observational state
    local_stats <- combinecovs(param$sigmas, param$mus, param$Ns, which(local_exps == ii), c(j, parents) - param$bgn)
    localparam <- BGeaugment(local_stats$sigma, local_stats$mu, local_stats$N, param$n, param$am, param$aw, param$logedgepmat)
    if (length(parents) > 0) {
      outscore <- Bestie:::DAGparametersCore(1, 1:length(parents) + 1, localparam)
    }  else {
      outscore <- Bestie:::DAGparametersCore(1, parents, localparam)
    }
  }
  outscore
}