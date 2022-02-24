
### These user defined score functions are for known perfect interventions


### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam, usrpar = list(Tmat = NULL, pctesttype = "bge", am = 1, chi = 1, edgepf = 1, edgepmat = NULL)){

  n <- initparam$n
  Tmat <- usrpar$Tmat
  nodeparams <- vector("list", n)
  for (jj in 1:n){
    nint_obs <- which(Tmat[, jj] == 0)
    if (length(nint_obs) < 2) {
      stop("Not enough observational data.")
    } else {
      nodeparams[[jj]] <- scoreparameters(scoretype = usrpar$pctesttype, data = initparam$data[nint_obs, ],
                                          weightvector = initparam$weightvector[nint_obs],
                                          bgepar = list(am = usrpar$am), bdepar = list(chi = usrpar$chi, edgepf = usrpar$edgepf),
                                          bdecatpar = list(chi = usrpar$chi, edgepf = usrpar$edgepf),
                                          edgepmat = usrpar$edgepmat)
    }

  }
  initparam$nodeparams <- nodeparams
  initparam
}



### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j, parentnodes, n, param) {
  BiDAG:::DAGcorescore(j, parentnodes, n, param$nodeparams[[j]])
}
