
### range of seeds
seed_numbers <- 100 + 1:100 # the seeds

### Settings
n <- 100 # number of nodes
N <- 4*n # number of observations
exp_ints_vec <- c(0, 0.01, 0.03, 0.1, 0.3, 1) # expected number of interventions per observation 
exp_parents <- 2 # expected number of parents

inode_frac <- 0.1 # fraction of nodes to intervene upon

ni <- inode_frac*n # number of possible intervention nodes

for (exp_ints in exp_ints_vec) {

# store values for later dataframe
setup_vec <- c(n, N, exp_ints, exp_parents, ni)
names(setup_vec) <- c("n", "N", "ints", "parents", "ni")
# create a name for the directory to store stuff
subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
dir_name <- paste("./simresults_p", subdir_name, sep = "/")
dir_name2 <- paste("./mergedresults_p", subdir_name, sep = "/")

method_files <- c("BiDAG", "GIES", "IGSP")

#if (!file.exists(paste0(dir_name, ".Rdata"))) {

results_df <- NULL
times_df <- NULL

for (seed_number in seed_numbers) {
  seed_part <- paste("", "seed", seed_number, sep = "_")
  for (ii in 1:length(method_files)) {
    if (file.exists(paste0(dir_name, "/", subdir_name, seed_part, "_", method_files[ii], ".Rdata"))) {
      load(paste0(dir_name, "/", subdir_name, seed_part, "_", method_files[ii], ".Rdata"))
      results_df <- rbind(results_df, result_df)
      times_df <- rbind(times_df, time_df)
    }
  }
}

save(results_df, times_df, file = paste0(dir_name2, ".Rdata"))

}

