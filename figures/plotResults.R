for (standardize in c(TRUE, FALSE)) {
  
  std_string <- ""
  if (!standardize) {
    std_string <- "_unscaled"
  }
  

# if we only have a mean shift
only_shift <- FALSE#TRUE

### Settings
n <- 100 # number of nodes
exp_parents <- 2 # expected number of parents
ni <- 10

all_results_df <- NULL

for (scaleN in c(1, 10)) {

  exp_ints_vec <- scaleN*c(5, 10, 20) # expected number of interventions per observation 
  N <- 4*n*scaleN # number of observations
  
for (exp_ints in exp_ints_vec) { 

  # store values for later dataframe
  setup_vec <- c(n, N, exp_ints, exp_parents, ni)
  names(setup_vec) <- c("n", "N", "ints", "parents", "ni")
  # create a name for the directory to store stuff
  subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
  dir_name <- paste(paste0("../mergedresults_", ifelse(only_shift, "s", "g"), std_string), subdir_name, sep = "/")

  load(paste0(dir_name, ".Rdata"))

  all_results_df <- rbind(all_results_df, results_df)
}
  
}

results_df <- all_results_df

results_df$TPR <- as.numeric(as.character(results_df$TPR))
results_df$FPRp <- as.numeric(as.character(results_df$FPR_P))
results_df$value <- as.numeric(as.character(results_df$value))
results_df$ints <- as.numeric(as.character(results_df$ints))
results_df$N <- as.numeric(as.character(results_df$N))

results_df$method <- factor(results_df$method, levels = c("UT-IGSP", "MAP", "MCMC"))

library(tidyverse)

results_df %>% group_by(method, value, N, ints) %>% 
  summarise(n = n(), missing = mean(is.na(TPR)*is.na(FPRp)), 
            TPR = mean(TPR, na.rm = TRUE), 
            FPRp = mean(FPRp, na.rm = TRUE)) -> summarised_results_df

# check which runs led to an empty graph
results_df %>% group_by(n, N, ints, parents, ni, seed, method) %>% 
  summarise(empty = prod((TPR == 0)*(FPRp == 0))) -> empty_results_df


### Make plots 

for (plot_window_size in round(c(1.2,10)/10, 2)) {

p100 <- ggplot() + 
  geom_point(data = results_df, aes(x = FPRp, y = TPR, fill = method, color=method), 
              alpha = 0.5, size = 1, shape = 21, stroke = 0) + 
  geom_path(data = summarised_results_df, aes(x = FPRp, y = TPR, colour = method),
             size = 1.5) + 
  geom_point(data = summarised_results_df %>% filter(value == 2.23e-06 | value == 0.1),
             aes(x = FPRp, y = TPR, color=method), 
             alpha = 0.75, size = 3, shape = 19) +
  coord_fixed(ratio = 1, xlim = c(0, plot_window_size), ylim = c(1 - plot_window_size, 1), expand = FALSE) + 
  annotate("segment", x = 0, xend = 1, y = 1, yend = 1) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 1) +
  scale_color_manual(values = c("darkkhaki", "dodgerblue", "darkorchid4")) + 
  facet_wrap(~ N + ints, labeller = label_bquote(N[I] == .(ints)*"," ~ N == .(N)), nrow = 2) + 
  scale_fill_manual(values = c(paste0("#", paste(as.hexmode(round(rowSums(cbind(col2rgb("darkkhaki"), rep(255, 3)))/2)), collapse="")),
paste0("#", paste(as.hexmode(round(rowSums(cbind(col2rgb("dodgerblue"), rep(255, 3)))/2)), collapse="")),
paste0("#", paste(as.hexmode(round(rowSums(cbind(col2rgb("darkorchid4"), rep(255, 3)))/2)), collapse=""))))

plot_name <- paste(paste(names(setup_vec[-c(2:3)]), setup_vec[-c(2:3)], sep = "_"), collapse = "_")
plot_name <- paste0(ifelse(only_shift, "mean_shift_", "general_"), plot_name)

p100 + theme(panel.spacing.x = unit(1.4, "lines"))
ggsave(paste0("./", plot_name, "_zoom_", plot_window_size, std_string, ".pdf"), width=8, height=5)  

}

}
