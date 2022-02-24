
# if we only want a mean shift
only_shift <- FALSE#TRUE

### range of seeds
seed_numbers <- 100 + 1:100 # the seeds

for (scaleN in c(1, 10)) {

### Settings
n <- 100 # number of nodes
N <- 4*n*scaleN # number of observations
exp_ints_vec <- scaleN*n*c(5, 10, 20)/100 
exp_parents <- 2 # expected number of parents
ni <- 10 # number of possible interventions

for (exp_ints in exp_ints_vec) {

# store values for later dataframe
setup_vec <- c(n, N, exp_ints, exp_parents, ni)
names(setup_vec) <- c("n", "N", "ints", "parents", "ni")
# create a name for the directory to store stuff
subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
dir_name <- paste(paste0("./simresults_", ifelse(only_shift, "s", "g")), subdir_name, sep = "/")
dir_name2 <- paste(paste0("./mergedresults_", ifelse(only_shift, "s", "g")), subdir_name, sep = "/")

method_files <- c("BiDAG", "UT-IGSP")

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

}
