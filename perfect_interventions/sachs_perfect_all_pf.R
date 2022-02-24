
norm_flag <- TRUE

# load libraries
library(pcalg)
library(graph)
library(BiDAG)
# Use BiDAG with intervention scoring
insertSource("./usrscorefns.R", package = "BiDAG")
source("../Rfns/intfns.R") # load other functions
source("../Rfns/sp_gauss.R") # load SP functions

seed_number <- 42

## Sachs data
sachs_data <- read.csv("../data/sachs.csv", sep = "\t")
# we take the first 7 experiments following Yang et al
sachs_data_sel <- sachs_data[which(!sachs_data[, 12] %in% 8:9), ]
n_exps <- table(sachs_data_sel[, 12])
# and treat the first two as observational
t.list <- list(integer(0), integer(0), c(7), c(9), c(4), c(2), c(5), 
               c(7), c(9), c(4), c(2), c(5))

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

data <- as.matrix(log(sachs_data_sel[, -12] + 0.5)) # log transform
n <- ncol(data)

if(norm_flag) { # correct each component of the data
for (ii in 1:length(n_exps)) { # batch correction
  rows_sel <- which(sachs_data_sel[, 12] == ii)
  data[rows_sel, ] <- scale(data[rows_sel, ], scale = FALSE)
}
}

# ground truth dag from Wang et al
ground.truth.Wang <- list(c(1,2), c(2,6), c(3,4), c(3,9), c(4,9), c(5,3), c(5,4), c(5,7), c(8,1), c(8,2), c(8,6), c(8,7), c(8,10), c(8,11), c(9,1), c(9,2), c(9,10), c(9,11))
trueDAGwang <- matrix(0, 11, 11)
for (t in ground.truth.Wang) {
  trueDAGwang[t[1], t[2]] <- 1
}

#ground truth dag from bnlearn
# this is very different from the one from Wang et al!!!
ground.truth.bnlearn <- list(c(1,2), c(2,6), c(6,7), #Raf, Mek, Erk
                     c(3,4), c(3,5), c(5,4), #Plcg and PIPs 
                     c(8,1), c(8,2), c(8,6), c(8,7), c(8,10), c(8,11), #PKA
                     c(9,1), c(9,2), c(9,8), c(9,10), c(9,11)) #PKCtrueDAG <- matrix(0, nrow=11, ncol=11)
trueDAGbn <- matrix(0, 11, 11)
for (t in ground.truth.bnlearn) {
  trueDAGbn[t[1], t[2]] <- 1
}

trueEGbn <- DAG2EG(trueDAGbn, targets)
trueEGwang <- DAG2EG(trueDAGwang, targets)


method_vec <- rep(NA, 4)
names(method_vec) <- c("method", "parameter", "value", "gt")


### MAP and MCMC search
# fill out missing results
if (!file.exists(paste0("./sachsresults/sachs_BiDAG", ifelse(norm_flag,"_norm",""), "_all.Rdata"))) {
  
result_df <- NULL

### MAP search

MAP_ams <- signif(exp(c(-4:0*3/2, 1:4/2))/10, 3)

am_samp <- MAP_ams[5]
am_value <- am_samp

#method_vec[2] <- "am"

method_vec[2] <- "pf"

#for (am_value in MAP_ams) {
for (pf in c(5:1*-3, 0:5*100)) {
  method_vec[1] <- "MAP"  
  method_vec[3] <- pf
  
  scoreObject <- scoreparameters(scoretype = "usr", data = data, 
                                 usrpar = list(pctesttype = "bge", Tmat = Tmat, am = am_value,
                                               edgepmat = matrix(2^pf, n, n)))

  set.seed(seed_number) # set seed

  bestDAGs <- iterativeMCMC(scoreObject, scoreout = TRUE, verbose = TRUE)

  bestDAG <- bestDAGs$DAG
#(bestDAGscore <- bestDAGs$score)
#DAGscore(scoreObject, bestDAG)

## Turn DAG into essential graph
  bestEG <- DAG2EG(bestDAG, targets)

  method_vec[4] <- "bn"
  MAPresult <- compareEGs(bestEG, trueEGbn)
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

  posteriors <- EGedgePosterior(orderresult$traceadd$incidence, targets = targets)

# Consensus graphs
  method_vec[1] <- "MCMC"

  MCMCEG <- consensusEG(posteriors, 0.5)

  method_vec[4] <- "bn"
  MCMCresult <- compareEGs(MCMCEG, trueEGbn)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MCMCresult))))
  method_vec[4] <- "wang"
  MCMCresult <- compareEGs(MCMCEG, trueEGwang)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, MCMCresult))))
  
  print(result_df)
}

save(result_df, file = paste0("./sachsresults/sachs_BiDAG", ifelse(norm_flag,"_norm",""), "_all.Rdata"))

}


## GIES
# fill out missing simulations
if (!file.exists(paste0("./sachsresults/sachs_GIES", ifelse(norm_flag,"_norm",""), "_all.Rdata"))) {

  result_df <- NULL  

method_vec[1:2] <- c("GIES", "lambda")
gies_lambdas <- signif(0.5*exp(-20:10/3)*log(n), 3) # default is 0.5 
# we need more penalisation with larger n

for (lambda_scale in gies_lambdas) {

  method_vec[3] <- lambda_scale
  gies_score <- new("GaussL0penIntScore", data, targets, attr(targetsmat,"index"), 
                  lambda = lambda_scale*log(nrow(data)))

  set.seed(seed_number) # set seed
  
  gies_fit <- gies(gies_score)

  giesEG <- 1*as(gies_fit$essgraph, "matrix")

  method_vec[4] <- "bn"
  GIESresult <- compareEGs(giesEG, trueEGbn)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, GIESresult))))
  method_vec[4] <- "wang"
  GIESresult <- compareEGs(giesEG, trueEGwang)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, GIESresult))))
}
  
save(result_df, file = paste0("./sachsresults/sachs_GIES", ifelse(norm_flag,"_norm",""), "_all.Rdata"))

}

### SP

# fill out missing simulations
if (!file.exists(paste0("./sachsresults/sachs_IGSP", ifelse(norm_flag,"_norm",""), "_all.Rdata"))) {
  
  result_df <- NULL  

# turn data into format for SP 

SP_data.list <- lapply(1:length(targets), function(t) data[which(attr(targetsmat,"index")==t), ])

suffstat <- list(C = SP_cor(SP_data.list[[1]]), n = SP_nrow(SP_data.list[[1]]))
if (length(targets) > 1) {
  intsuffstat <- lapply(2:length(targets), function(t) list(C = SP_cor(SP_data.list[[t]]), n = SP_nrow(SP_data.list[[t]])))
  inttargets <- targets[2:length(targets)]
} else {
  intsuffstat <- list(NULL)
  inttargets <- list(NULL)
}

method_vec[1:2] <- c("IGSP", "alpha")
SP_alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)  

for (alpha in SP_alphas) {
  
  method_vec[3] <- alpha
  set.seed(seed_number) # set seed
  
  SP_fit <- sp.restart.alg(suffstat, intsuffstat, inttargets, alpha)
  SPEG <- 1*as(dag2essgraph(as(SP_fit, "graphNEL"), targets), "matrix")
 
  method_vec[4] <- "bn"
  SPresult <- compareEGs(SPEG, trueEGbn)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, SPresult))))
  method_vec[4] <- "wang"
  SPresult <- compareEGs(SPEG, trueEGwang)
  result_df <- rbind(result_df, data.frame(t(c(method_vec, SPresult))))
}

save(result_df, file = paste0("./sachsresults/sachs_IGSP", ifelse(norm_flag,"_norm",""), "_all.Rdata"))

}
