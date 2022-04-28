source("./sachssetup.R")

# if we want UT-IGSP evaluated 
# (if needs to be run before in python, after making the data)

eval_SP <- TRUE

# load libraries
library(BiDAG)

# Use BiDAG with intervention scoring
insertSource("../usrscorefns.R", package = "BiDAG")

seed_number <- 42

method_vec <- rep(NA, 4)
names(method_vec) <- c("method", "parameter", "value", "gt")

if (!file.exists("./sachsresults/sachs_CandP.Rdata")) {
  
CandPedges <- list(c(2,1), #Raf, Mek,
                c(3,4), c(3,5), c(5,4), #Plcg and PIPs 
                c(6,8), c(7,6), c(7, 8), # Akt, Erk and PKA
                c(9,10), c(9,11), c(11,10))
CandPDAG <- matrix(0, 11, 11)
for (t in CandPedges) {
  CandPDAG[t[1], t[2]] <- 1
}
CandPEG <- DAG2EG(CandPDAG, targets)

method_vec[1] <- "C+P"
method_vec[2] <- "median"

result_df <- NULL
method_vec[4] <- "bn"
CandPresult <- compareEGs(CandPEG, trueEGbn)
result_df <- rbind(result_df, data.frame(t(c(method_vec, CandPresult))))
method_vec[4] <- "sachs"
CandPresult <- compareEGs(CandPEG, trueEGsachs)
result_df <- rbind(result_df, data.frame(t(c(method_vec, CandPresult))))
method_vec[4] <- "wang"
CandPresult <- compareEGs(CandPEG, trueEGwang)
result_df <- rbind(result_df, data.frame(t(c(method_vec, CandPresult))))

save(result_df, file = "./sachsnets/sachs_CandP.Rdata")
}

### MAP and MCMC search
# fill out missing simulations
if (!file.exists("./sachsnets/sachs_BiDAG.Rdata")) {
  
result_df <- NULL

MAPtime <- 0
MCMCtime <- 0

### MAP search

MAP_ams <- signif(exp(c(-4:0*3/2, 1:4/2))/10, 3)

am_samp <- MAP_ams[5]
am_value <- am_samp

#method_vec[2] <- "am"

method_vec[2] <- "pf"

for (pf in c(3*(3:-2),-8:-10,-10-1*(1:2)/4)){
  method_vec[1] <- "MAP"  
  method_vec[3] <- pf
  
  scoreObject <- scoreparameters(scoretype = "usr", data = data, usrpar = list(pctesttype = "bge", Imat = Imat, am = am_value,
                                                                               edgepmat = matrix(2^pf, n+6, n+6)))
  
  set.seed(seed_number) # set seed

  bestDAGs <- iterativeMCMC(scoreObject, scoreout = TRUE, verbose = TRUE)

  bestDAG <- bestDAGs$DAG
#(bestDAGscore <- bestDAGs$score)
#DAGscore(scoreObject, bestDAG)
  
  if (ncol(Imat) > 0) {
    bestDAG <- bestDAG[-c(1:ncol(Imat)), -c(1:ncol(Imat))] 
  }

## Turn DAG into essential graph
  bestEG <- DAG2EG(bestDAG, targets)

  method_vec[4] <- "bn"
  MAPresult <- compareEGs(bestEG, trueEGbn)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MAPresult))))
  method_vec[4] <- "sachs"
  MAPresult <- compareEGs(bestEG, trueEGsachs)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MAPresult))))
  method_vec[4] <- "wang"
  MAPresult <- compareEGs(bestEG, trueEGwang)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MAPresult))))
  
  print(result_df)
  
## Order sample

set.seed(seed_number) # set seed

orderresult <- orderMCMC(scoreObject, startspace = bestDAGs$endspace, 
                         MAP = FALSE, plus1 = TRUE, chainout = TRUE, 
                         startorder = bestDAGs$maxorder,
                         iterations = 20*round(n*n*log(n)),
                         verbose = TRUE)

orderchain <- orderresult$traceadd$incidence

if (ncol(Imat) > 0) {
  for (ii in 1:length(orderchain)) {
    orderchain[[ii]] <- orderchain[[ii]][-c(1:ncol(Imat)), -c(1:ncol(Imat))]
  }
}

  posteriors <- EGedgePosterior(orderchain, targets = targets)

# Consensus graphs
  method_vec[1] <- "MCMC"

  MCMCEG <- consensusEG(posteriors, 0.5)

  method_vec[4] <- "bn"
  MCMCresult <- compareEGs(MCMCEG, trueEGbn)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MCMCresult))))
  method_vec[4] <- "sachs"
  MCMCresult <- compareEGs(MCMCEG, trueEGsachs)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MCMCresult))))
  method_vec[4] <- "wang"
  MCMCresult <- compareEGs(MCMCEG, trueEGwang)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MCMCresult))))

}


save(result_df, file = "./sachsnets/sachs_BiDAG.Rdata")

}

if (eval_SP) {

### SP

# fill out missing simulations
if (!file.exists("./sachsnets/sachs_UT-IGSP.Rdata")) {
  
  result_df <- NULL
  SP_DAGs <- NULL

  if (file.exists("./sachsnets/sachs_UT-IGSP_DAGs.csv")) {
    SP_DAGs <- read.csv("./sachsnets/sachs_UT-IGSP_DAGs.csv", header = FALSE)
  }

  method_vec[1:2] <- c("UT-IGSP", "alpha")
  SP_alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)  
  
  if (!is.null(SP_DAGs)) {
  
    for (ii in 1:length(SP_alphas)) {
  
      alpha <- SP_alphas[ii]
      method_vec[3] <- alpha

      SP_fit <- SP_to_DAG(as.character(SP_DAGs[ii, ]))
  
      SPEG <- 1*as(dag2essgraph(as(SP_fit, "graphNEL"), targets), "matrix")
 
      method_vec[4] <- "bn"
      SPresult <- compareEGs(SPEG, trueEGbn)
      result_df <- rbind(result_df, data.frame(t(c(method_vec, SPresult))))
      method_vec[4] <- "sachs"
      SPresult <- compareEGs(SPEG, trueEGsachs)
      result_df <- rbind(result_df, data.frame(t(c(method_vec, SPresult)))) 
      method_vec[4] <- "wang"
      SPresult <- compareEGs(SPEG, trueEGwang)
      result_df <- rbind(result_df, data.frame(t(c(method_vec, SPresult)))) 
    }

    save(result_df, file = "./sachsnets/sachs_UT-IGSP.Rdata")
  }

}

}
  



