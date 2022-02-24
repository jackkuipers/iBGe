for (standardize in c(TRUE, FALSE)) {

std_string <- ""
if (!standardize) {
  std_string <- "_unscaled"
}
  
### Settings
n <- 100 # number of nodes
N <- 4*n # number of observations
exp_ints_vec <- c(0, 0.01, 0.03, 0.1, 0.3, 1) # expected number of interventions per observation 
exp_parents <- 2 # expected number of parents

times_combined_df <- NULL

inode_frac <- 0.1 
ni <- inode_frac*n # number of possible intervention nodes

all_results_df <- NULL
all_times_df <- NULL

for (exp_ints in exp_ints_vec) { 

  # store values for later dataframe
  setup_vec <- c(n, N, exp_ints, exp_parents, ni)
  names(setup_vec) <- c("n", "N", "ints", "parents", "ni")
  # create a name for the directory to store stuff
  subdir_name <- paste(paste(names(setup_vec), setup_vec, sep = "_"), collapse = "_")
  dir_name <- paste(paste0("../perfect_interventions/mergedresults_p", std_string), subdir_name, sep = "/")

  load(paste0(dir_name, ".Rdata"))

  all_results_df <- rbind(all_results_df, results_df)
  all_times_df <- rbind(all_times_df, times_df)
}

results_df <- all_results_df

results_df$TPR <- as.numeric(as.character(results_df$TPR))
results_df$FPRp <- as.numeric(as.character(results_df$FPR_P))
results_df$value <- as.numeric(as.character(results_df$value))
results_df$ints <- as.numeric(as.character(results_df$ints))

results_df$method <- factor(results_df$method, levels = c("GIES", "IGSP", "MAP", "MCMC"))

times_df <- all_times_df

times_df$time <- as.numeric(as.character(times_df$time))
times_df$ints <- as.numeric(as.character(times_df$ints))

times_df$method <- factor(times_df$method, levels = c("GIES", "IGSP", "MAP", "MCMC"))



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
  geom_point(data = summarised_results_df %>% filter(value == 2.23e-06 | value == 0.1 | value == 1.650),
             aes(x = FPRp, y = TPR, color=method), 
             alpha = 0.75, size = 3, shape = 19) +
  coord_fixed(ratio = 1, xlim = c(0, plot_window_size), ylim = c(1 - plot_window_size, 1), expand = FALSE) + 
  annotate("segment", x = 0, xend = 1, y = 1, yend = 1) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 1) +
  scale_color_manual(values = c("firebrick3", "forestgreen", "dodgerblue", "darkorchid4")) + 
  facet_wrap(~ ints, labeller = label_bquote(rho == .(ints)), nrow = 2) #+ 
  scale_fill_manual(values = c(paste0("#", paste(as.hexmode(round(rowSums(cbind(col2rgb("firebrick3"), rep(255, 3)))/2)), collapse="")),
paste0("#", paste(as.hexmode(round(rowSums(cbind(col2rgb("forestgreen"), rep(255, 3)))/2)), collapse="")),
paste0("#", paste(as.hexmode(round(rowSums(cbind(col2rgb("dodgerblue"), rep(255, 3)))/2)), collapse="")),
paste0("#", paste(as.hexmode(round(rowSums(cbind(col2rgb("darkorchid4"), rep(255, 3)))/2)), collapse=""))))

plot_name <- paste(paste(names(setup_vec[-3]), setup_vec[-3], sep = "_"), collapse = "_")
plot_name <- paste0("perfect_", plot_name)

p100 + theme(panel.spacing.x = unit(1.4, "lines"))
ggsave(paste0("./", plot_name, "_zoom_", plot_window_size, std_string, ".pdf"), width=8, height=5)  

}
# remove times of empty networks of IGSP
times_df %>% group_by(n, N, ints, parents, ni, seed, method) %>% 
  summarise(time = time) %>% head
times_df <- merge(times_df, empty_results_df)

times_combined_df <- rbind(times_combined_df, times_df)


colrysall <- c("#cd2626",
               "#228b22",
               "#1e90ff",
               "#68228b")

labels <- times_combined_df$method %>% levels()

ggplot(times_combined_df %>% filter(empty == 0), 
       aes(x = method, y = time, colour = method, fill = method)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(shape = 16, alpha = 0.25, position=position_jitterdodge()) + scale_y_log10() +
  scale_color_manual(values = colrysall, name= " ", labels = labels) + 
  scale_fill_manual(values = colrysall, name= " ", labels = labels) +
  theme(legend.key.height = unit(0.9, "cm")) + 
  facet_wrap(~ ints, labeller = label_bquote(rho == .(ints)), nrow = 2) -> p_t

cairo_pdf(paste0("./", plot_name, "_times", std_string, ".pdf"), width=8, height=5)
print(p_t)
dev.off()
#ggsave(paste0("./figures/", plot_name, "_times.pdf"), width=8, height=5)  

}

