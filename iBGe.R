
# when we don't know the interventions and they just lead to a mean shift

local_run <- FALSE#TRUE

standardize <- TRUE # whether to standardise the data

# if we just want the data

just_data <- FALSE#TRUE

# if we want UT-IGSP evaluated 
# (if needs to be run before in python, after making the data)

eval_SP <- TRUE

# if we only want a mean shift
only_shift <- FALSE


# to run locally we need to set the seed and interventions
# this is done externally on the cluster runs

if (local_run) { 
  kk <- 1
  seed_number <- 101 # the seed
}

# load libraries
library(pcalg)
library(graph)
library(BiDAG)

# Use BiDAG with intervention scoring
insertSource("./usrscorefns.R", package = "BiDAG")
source("./Rfns/intfns.R") # load other functions

### Settings
# for 100 nodes takes up to an hour for MAP and MCMC
# with all settings for ROC curve

#scaleN <- 10 # normally 1
scaleN <- 1

n <- 100 # number of nodes
N <- 4*n*scaleN # number of observations

exp_ints <- scaleN*n*c(5, 10, 20)[kk]/100 # number of observations per intervention 
exp_parents <- 2 # expected number of parents

ni <- 10 # number of possible interventions

# store values for later dataframe
setup_vec <- c(n, N, exp_ints, exp_parents, ni, seed_number)
names(setup_vec) <- c("n", "N", "ints", "parents", "ni", "seed")
# create a name for the directory to store stuff
f_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
subdir_name <- paste(paste(names(setup_vec[-6]), setup_vec[-6], sep = "_"), collapse = "_")
if (only_shift) {
  if (standardize) {
    top_dir_name <- "./simresults_s"
  } else {
    top_dir_name <- "./simresults_s_unscaled"
  }
} else {
  if (standardize) {
    top_dir_name <- "./simresults_g"
  } else {
    top_dir_name <- "./simresults_g_unscaled"
  }
}
dir_name <- paste(top_dir_name, subdir_name, sep = "/")

if (!dir.exists(top_dir_name)) { # create directory if none exists
  dir.create(top_dir_name)
}
if (!dir.exists(dir_name)) { # create directory if none exists
  dir.create(dir_name)
}


### Generate data

set.seed(seed_number) # set seed


# generate which intervention each node gets
#Ts <- c(rep(0, N*(1-exp_ints)), 
#        sort(sample.int(ni, N*exp_ints, replace = TRUE)))
Ts <- c(rep(0, N - exp_ints*ni), 
        rep(1:ni, each = exp_ints))

# can remove non-existent interventions,
# but we remove later anyway
#Ts <- rep(1:length(unique(Ts))-1, table(Ts))

# which nodes are targeted
Tlist <- vector("list", ni)
# how much mean shift they create
Tshifts <- rnorm(ni, 1, 0.1)
# how much they damp the contribution from parents
Tdamps <- runif(ni, min = 0.1) #rnorm(ni, 1, 0.1)

for (ii in 1:ni){
  Tlist[[ii]] <- sample.int(n, 1+rpois(1, 1))
}

# generate intervention matrix
Tmat <- matrix(0, ncol = n, nrow = N)
Imat <- matrix(0, ncol = ni, nrow = N)
shift_vec <- rep(0, N)
damp_vec <- rep(0, N) 
for (ii in 1:N){
  if (Ts[ii] > 0) {
    Tmat[ii, Tlist[[Ts[ii]]]] <- 1
    Imat[ii, Ts[ii]] <- 1
    shift_vec[ii] <- Tshifts[Ts[ii]]
    damp_vec[ii] <- Tdamps[Ts[ii]]
  }
}

# remove interventions that do not appear
# this should not happen in the current setting!
Imat <- Imat[, which(colSums(Imat) > 0)]

# turn intervention matrix into format for pcalg
# should be same as used elements of Tlist 
# with a padded entry at the start for observational data
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

set.seed(seed_number) # set seed

# generate random DAG
trueDAGedges <- as(pcalg::randDAG(n = n, d = 2*exp_parents, wFUN = list(wFUN, 0.25, 1)), "matrix")
trueDAG <- 1*(trueDAGedges != 0)

trueEG <- DAG2EG(trueDAG, targets)

# generate simulated data
if (only_shift) {
  data <- irmvDAG(trueDAGedges, Tmat, shift = shift_vec, scale = standardize)
} else {
  data <- irmvDAG(trueDAGedges, Tmat, shift = shift_vec, damp = damp_vec, scale = standardize)
}

method_vec <- rep(NA, 3)
names(method_vec) <- c("method", "parameter", "value")

if(just_data) { # write files
if (!file.exists(paste0(dir_name, "/", f_name, "_data.csv"))) {
  write.csv(data, paste0(dir_name, "/", f_name, "_data.csv"), 
                         row.names = FALSE)
  write.csv(Ts, paste0(dir_name, "/", f_name, "_Ts.csv"), 
            row.names = FALSE)
}
}

if (!just_data) {

### MAP and MCMC search
# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_BiDAG.Rdata"))) {
  
result_df <- NULL

MAPtime <- 0
MCMCtime <- 0

### MAP search

MAP_ams <- signif(exp(c(-4:0*3/2, 1:4/2))/10, 3)

am_samp <- MAP_ams[5]

method_vec[2] <- "am"

# Start the clock!
ptm <- proc.time()

for (am_value in MAP_ams) {
  method_vec[1] <- "MAP"  
  method_vec[3] <- am_value
  
  scoreObject <- scoreparameters(scoretype = "usr", data = data, usrpar = list(pctesttype = "bge", Imat = Imat, am = am_value))
  
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

  MAPresult <- compareEGs(bestEG, trueEG)

  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, MAPresult))))

# Stop the clock
MAPtime <- MAPtime + (proc.time() - ptm)[1]

## Order sample

# Start the clock!
ptm <- proc.time()

set.seed(seed_number) # set seed

orderresult <- orderMCMC(scoreObject, startspace = bestDAGs$endspace, 
                         MAP = FALSE, plus1 = TRUE, chainout = TRUE, 
                         startorder = bestDAGs$maxorder,
                         iterations = 20*round(n*n*log(n)),
                         verbose = TRUE)

# Stop the clock
MCMCtime <- MCMCtime + (proc.time() - ptm)[1]

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

  MCMCresult <- compareEGs(MCMCEG, trueEG)
  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, MCMCresult))))
}

time_df <- data.frame(n = n, N = N, ints = exp_ints, 
                      parents = exp_parents, ni = ni, seed = seed_number,
                      method = c("MAP", "MCMC"), 
                      time = c(MAPtime, MCMCtime))

save(result_df, time_df, file = paste0(dir_name, "/", f_name, "_BiDAG.Rdata"))

}

if (eval_SP) {

### SP

# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_UT-IGSP.Rdata"))) {
  
  result_df <- NULL
  SP_DAGs <- NULL

  if (file.exists(paste0(dir_name, "/", f_name, "_UT-IGSP_DAGs.csv"))) {
    SP_DAGs <- read.csv(paste0(dir_name, "/", f_name, "_UT-IGSP_DAGs.csv"), header = FALSE)
  }

  method_vec[1:2] <- c("UT-IGSP", "alpha")
  SP_alphas <- signif(exp(c(-4:0*3/2, 1:4/2))/1e5, 3)  
  
  if (!is.null(SP_DAGs)) {
  
    for (ii in 1:length(SP_alphas)) {
  
      alpha <- SP_alphas[ii]
      method_vec[3] <- alpha

      SP_fit <- SP_to_DAG(as.character(SP_DAGs[ii, ]))
  
      if (is.null(SP_fit)) {
        SPresult <- rep(NA, 5)
        names(SPresult) <- c("TP", "FP", "SHD", "TPR", "FPR_P")
      } else {
        SPEG <- 1*as(dag2essgraph(as(SP_fit, "graphNEL"), targets), "matrix")
        SPresult <- compareEGs(SPEG, trueEG)
      }
  
      result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, SPresult))))
    }

    save(result_df, file = paste0(dir_name, "/", f_name, "_UT-IGSP.Rdata"))
  }

}

}
  
}


