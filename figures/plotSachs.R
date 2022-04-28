
all_results_df <- NULL

for (method in c("GIES", "IGSP")) { 

  dir_name <- paste0("../perfect_interventions/sachsresults/sachs_", method)

  load(paste0(dir_name, ".Rdata"))
  
  all_results_df <- rbind(all_results_df, result_df)
}

for (method in c("BiDAG", "CandP", "UT-IGSP")) { 
  
  dir_name <- paste0("../sachsanalysis/sachsnets/sachs_", method)
  
  load(paste0(dir_name, ".Rdata"))
  
  all_results_df <- rbind(all_results_df, result_df)
}


results_df <- all_results_df

results_df$TP <- as.numeric(as.character(results_df$TP))
results_df$FP <- as.numeric(as.character(results_df$FP))
results_df$value <- as.numeric(as.character(results_df$value))

results_df$SHD <- as.numeric(as.character(results_df$SHD))

results_df$method <- factor(results_df$method, levels = c("GIES", "IGSP", "UT-IGSP", "C+P", "MAP", "MCMC"))


library(tidyverse)

results_df <- filter(results_df, FP < 21, gt != "wang")

### Make plots 

p100 <- ggplot() + 
  geom_point(data = results_df, aes(x = FP, y = TP, color=method, shape = method),
             size = 4, alpha = 0.67) + 
  geom_line(data = results_df, aes(x = FP, y = TP, color=method), linetype = 2) +
  scale_shape_manual(values = c(rep(20, 3), 15, rep(20, 2))) + 
  scale_color_manual(values = c("firebrick3", "forestgreen", "darkkhaki", "darkorange", "dodgerblue", "darkorchid4")) +
  facet_wrap(~ gt, ncol = 2, labeller = labeller(gt = c("bn" = "vs Sachs reduced model", "sachs"="vs Sachs full model"))) + xlim(0, 20) + ylim(0, 12.5)

p100 
ggsave(paste0("./sachs.pdf"), width=8, height=4)  

