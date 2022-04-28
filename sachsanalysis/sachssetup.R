# load libraries
library(pcalg)
library(graph)
source("../Rfns/intfns.R") # load other functions

## Sachs data
sachs_data <- read.csv("./sachsdata/sachs.csv", sep = "\t")
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
if (!file.exists("./sachsdata/sachs_data_sel.csv")) {
  write.csv(data, "./sachsdata/sachs_data_sel.csv", 
            row.names = FALSE)
  write.csv(Ts, "./sachsdata/sachs_Ts.csv", 
            row.names = FALSE)
}
