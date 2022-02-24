# if we want UT-IGSP evaluated 
# (if needs to be run before in python, after making the data)

eval_SP <- TRUE

# load libraries
library(pcalg)
library(graph)
library(BiDAG)

# Use BiDAG with intervention scoring
insertSource("./usrscorefns.R", package = "BiDAG")
source("./Rfns/intfns.R") # load other functions

seed_number <- 42

## Sachs data
sachs_data <- read.csv("./data/sachs.csv", sep = "\t")
# we take the first 7 experiments following Wang et al
sachs_data_sel <- sachs_data[which(sachs_data[, 12] < 8), ]
Ts <- sachs_data_sel[, 12] - 1
n_exps <- table(Ts)
# and treat the first two as observational
t.list <- list(integer(0), integer(0), c(7), c(9), c(4), c(2), c(5))

# generate intervention matrix
Tmat <- matrix(0, ncol = ncol(sachs_data_sel) - 1, nrow = nrow(sachs_data_sel))
for (ii in 1:length(n_exps)) {
  rows_sel <- which(sachs_data_sel[, 12] == ii)
  Tmat[rows_sel, t.list[[ii]]] <- 1
}

# turn intervention matrix into format for pcalg
targetsmat <- mgcv::uniquecombs(Tmat)
if (is.null(nrow(targetsmat))) {
  targets <- vector("list", 1)
  targets[[1]] <- integer(0)
} else {
  targets <- vector("list", nrow(targetsmat))
  for (ii in 1:nrow(targetsmat)) {
    targets[[ii]] <- which(targetsmat[ii,] == 1)
  }
}

# generate intervention matrix
Imat <- matrix(0, ncol = max(Ts), nrow = nrow(sachs_data_sel))
for (ii in 1:nrow(sachs_data_sel)){
  if (Ts[ii] > 0) {
    Imat[ii, Ts[ii]] <- 1
  }
}

# remove interventions that do not appear
# this should not happen in the current setting!
Imat <- Imat[, which(colSums(Imat) > 0)]


data <- as.matrix(log(sachs_data_sel[, -12] + 0.5)) # log transform
n <- ncol(data)

#ground truth dag from bnlearn
# this is quite different from the one from Wang et al
# but if we correct PLC to PIP3 the only differences are the missing edges
# in the original figure
ground.truth.bnlearn <- list(c(1,2), c(2,6), c(6,7), #Raf, Mek, Erk
                             c(3,4), c(3,5), c(5,4), #Plcg and PIPs 
                             c(8,1), c(8,2), c(8,6), c(8,7), c(8,10), c(8,11), #PKA
                             c(9,1), c(9,2), c(9,8), c(9,10), c(9,11)) #PKCtrueDAG <- matrix(0, nrow=11, ncol=11)
trueDAGbn <- matrix(0, 11, 11)
for (t in ground.truth.bnlearn) {
  trueDAGbn[t[1], t[2]] <- 1
}
colnames(trueDAGbn) <- colnames(data)
rownames(trueDAGbn) <- colnames(data)

# ground truth bn learn with missing edges included
ground.truth.sachs <- list(c(1,2), c(2,6), c(3,4), c(3,9), c(4,9), 
                           c(3,5), # this edge is wrong way in WANG et al!!!
                           c(5,4), c(5,7), c(8,1), c(8,2), c(8,6), c(8,7), c(8,10), c(8,11), c(9,1), c(9,2), c(9,10), c(9,11))
trueDAGsachs <- matrix(0, 11, 11)
for (t in ground.truth.sachs) {
  trueDAGsachs[t[1], t[2]] <- 1
}
colnames(trueDAGsachs) <- colnames(data)
rownames(trueDAGsachs) <- colnames(data)


# ground truth dag from Wang et al
# they misinterpreted the direction of PLC to PIP3!
ground.truth.Wang <- list(c(1,2), c(2,6), c(3,4), c(3,9), c(4,9), 
                          c(5, 3), # this edge!!!
                          c(5,4), c(5,7), c(8,1), c(8,2), c(8,6), c(8,7), c(8,10), c(8,11), c(9,1), c(9,2), c(9,10), c(9,11))
trueDAGwang <- matrix(0, 11, 11)
for (t in ground.truth.Wang) {
  trueDAGwang[t[1], t[2]] <- 1
}
colnames(trueDAGwang) <- colnames(data)
rownames(trueDAGwang) <- colnames(data)


trueEGbn <- DAG2EG(trueDAGbn, targets)
trueEGsachs <- DAG2EG(trueDAGsachs, targets)
trueEGwang <- DAG2EG(trueDAGwang, targets)


# write files
if (!file.exists("./data/sachs_data_sel.csv")) {
  write.csv(data, "./data/sachs_data_sel.csv", 
            row.names = FALSE)
  write.csv(Ts, "./data/sachs_Ts.csv", 
            row.names = FALSE)
}


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

save(result_df, file = "./sachsresults/sachs_CandP.Rdata")
}

### MAP and MCMC search
# fill out missing simulations
if (!file.exists("./sachsresults/sachs_BiDAG.Rdata")) {
  
result_df <- NULL

MAPtime <- 0
MCMCtime <- 0

### MAP search

MAP_ams <- signif(exp(c(-4:0*3/2, 1:4/2))/10, 3)

am_samp <- MAP_ams[5]
am_value <- am_samp

#method_vec[2] <- "am"

method_vec[2] <- "pf"

for (pf in c(14:0*-1)) {
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
  
  print(result_df)
}


save(result_df, file = "./sachsresults/sachs_BiDAG.Rdata")

}

if (eval_SP) {

### SP

# fill out missing simulations
if (!file.exists("./sachsresults/sachs_UT-IGSP.Rdata")) {
  
  result_df <- NULL
  SP_DAGs <- NULL

  if (file.exists("./sachsresults/sachs_UT-IGSP_DAGs.csv")) {
    SP_DAGs <- read.csv("./sachsresults/sachs_UT-IGSP_DAGs.csv", header = FALSE)
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

    save(result_df, file = "./sachsresults/sachs_UT-IGSP.Rdata")
  }

}

}
  



