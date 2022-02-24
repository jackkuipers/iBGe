
# when we know the interventions and they are perfect

local_run <- FALSE

standardize <- TRUE # whether to standardise the data

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
source("../Rfns/intfns.R") # load other functions
source("../Rfns/sp_gauss.R") # load SP functions

### Settings
# for 100 nodes takes up to an hour for all the thresholds for MAP and MCMC
# a few seconds for GIES and can take several hours for IGSP

n <- 100 # number of nodes
N <- 4*n # number of observations

exp_ints <- c(0, 0.01, 0.03, 0.1, 0.3, 1)[kk] # expected number of interventions per observation 
exp_parents <- 2 # expected number of parents

inode_frac <- 0.1 # fraction of nodes to intervene upon
ni <- inode_frac*n # number of possible intervention nodes

# store values for later dataframe
setup_vec <- c(n, N, exp_ints, exp_parents, ni, seed_number)
names(setup_vec) <- c("n", "N", "ints", "parents", "ni", "seed")
# create a name for the directory to store stuff
f_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
subdir_name <- paste(paste(names(setup_vec[-6]), setup_vec[-6], sep = "_"), collapse = "_")
if (standardize) {
  top_dir_name <- "./simresults_p"
} else {
  top_dir_name <- "./simresults_p_unscaled"
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

# generate intervention matrix
Tmat <- matrix(0, ncol = n, nrow = N)
inodes <- sample(1:n, size = ni)
# old way with random interventions
#int_prob <- exp_ints/ni # probability of intervening per data point
#Tmat[, inodes] <- matrix(rbinom(ni*N, 1, int_prob), ncol = ni)
# new way with fixed set of interventions
int_rows <- exp_ints*N
if (int_rows > 0) {
  ints <- sample(inodes, int_rows, replace = TRUE)
  for (ii in 1:int_rows) {
    Tmat[ii, ints[ii]] <- 1
  }
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

set.seed(seed_number) # set seed

# generate random DAG
trueDAGedges <- as(pcalg::randDAG(n = n, d = 2*exp_parents, wFUN = list(wFUN, 0.25, 1)), "matrix")
trueDAG <- 1*(trueDAGedges != 0)

trueEG <- DAG2EG(trueDAG, targets)

# generate simulated data
data <- irmvDAG(trueDAGedges, Tmat, scale = standardize)


method_vec <- rep(NA, 3)
names(method_vec) <- c("method", "parameter", "value")


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
  scoreObject <- scoreparameters(scoretype = "usr", data = data, usrpar = list(pctesttype = "bge", Tmat = Tmat, am = am_value))

  set.seed(seed_number) # set seed

  bestDAGs <- iterativeMCMC(scoreObject, scoreout = TRUE, verbose = TRUE)

  bestDAG <- bestDAGs$DAG
#(bestDAGscore <- bestDAGs$score)
#DAGscore(scoreObject, bestDAG)

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

  posteriors <- EGedgePosterior(orderresult$traceadd$incidence, targets = targets)

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


## GIES
# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_GIES.Rdata"))) {

  result_df <- NULL  
    
# Start the clock!
ptm <- proc.time()

method_vec[1:2] <- c("GIES", "lambda")
gies_lambdas <- signif(0.5*exp(-4:4/3)*log(n), 3) # default is 0.5 
# we need more penalisation with larger n

for (lambda_scale in gies_lambdas) {

  method_vec[3] <- lambda_scale
  gies_score <- new("GaussL0penIntScore", data, targets, attr(targetsmat,"index"), 
                  lambda = lambda_scale*log(nrow(data)))

  set.seed(seed_number) # set seed
  
  gies_fit <- gies(gies_score)

  giesEG <- 1*as(gies_fit$essgraph, "matrix")

  GIESresult <- compareEGs(giesEG, trueEG)
  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, GIESresult))))
}
  
# Stop the clock
GIEStime <- (proc.time() - ptm)[1]

time_df <- data.frame(n = n, N = N, ints = exp_ints, 
                      parents = exp_parents, ni = ni, seed = seed_number,
                      method = c("GIES"), 
                      time = c(GIEStime))

save(result_df, time_df, file = paste0(dir_name, "/", f_name, "_GIES.Rdata"))

}

### SP

# fill out missing simulations
if (!file.exists(paste0(dir_name, "/", f_name, "_IGSP.Rdata"))) {
  
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

# Start the clock!
ptm <- proc.time()

method_vec[1:2] <- c("IGSP", "alpha")
SP_alphas <- signif(exp(c(-4:0*3/2, 1:4/2))/1e5, 3)  

for (alpha in SP_alphas) {
  
  method_vec[3] <- alpha
  set.seed(seed_number) # set seed
  
  SP_fit <- sp.restart.alg(suffstat, intsuffstat, inttargets, alpha)
  SPEG <- 1*as(dag2essgraph(as(SP_fit, "graphNEL"), targets), "matrix")
  
  SPresult <- compareEGs(SPEG, trueEG)
  result_df <- rbind(result_df, data.frame(t(c(setup_vec, method_vec, SPresult))))
}

# Stop the clock
SPtime <- (proc.time() - ptm)[1]


time_df <- data.frame(n = n, N = N, ints = exp_ints, 
                      parents = exp_parents, ni = ni, seed = seed_number,
                      method = c("IGSP"), 
                      time = c(SPtime))

save(result_df, time_df, file = paste0(dir_name, "/", f_name, "_IGSP.Rdata"))

}
