# ==============================================================================
# SCRIPT 2: Calibration Curve (Target Fst = 0.039)
# ==============================================================================

library(slendr)
library(tidyverse)

init_env()

res_dir <- "/cfs/klemming/projects/snic/naiss2024-6-260/mafalda/simulations/results"
fig_dir <- "/cfs/klemming/projects/snic/naiss2024-6-260/mafalda/simulations/results/figures"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
recomb_rate <- 2.54e-8 
mut_rate <- 2e-9       
seq_length <- 50000  
maf_threshold <- 0.05
n_fragments_calib <- 500   # Fast run just to find the mean/median Fst

# Target Parameters
target_mean_fst <- 0.039
target_median_fst <- 0.022

cat("\n--- Running Calibration Loop ---\n")
test_times <- c(50000, 100000, 200000, 300000, 400000, 500000)

calibration_results <- data.frame(Years = numeric(), Mean_Fst = numeric(), Median_Fst = numeric())
raw_fst_list <- list() # List to store raw output

for (t_years in test_times) {
  cat(sprintf("Testing split time: %d years...\n", t_years))
  atl <- population("Atlantic", time = t_years + 5000, N = 1000000)
  bal <- population("Baltic", parent = atl, time = t_years, N = 1000000)
  model <- compile_model(populations = list(atl, bal), generation_time = 3, direction = "backward")
  sampling_plan <- schedule_sampling(model, times = 0, list(atl, 50), list(bal, 50))
  
  calib_fst <- rep(NA, n_fragments_calib)
  for (i in 1:n_fragments_calib) {
    ts <- msprime(model, sequence_length = seq_length, recombination_rate = recomb_rate, samples = sampling_plan) %>% ts_mutate(mutation_rate = mut_rate)
    if (ts$num_mutations == 0) next
    
    all_nodes <- ts_nodes(ts)
    atl_nodes <- all_nodes %>% filter(sampled == TRUE, pop == "Atlantic") %>% pull(node_id)
    bal_nodes <- all_nodes %>% filter(sampled == TRUE, pop == "Baltic") %>% pull(node_id)
    
    geno_mat <- ts$genotype_matrix(samples = as.integer(c(atl_nodes, bal_nodes)))
    daf <- rowSums(geno_mat) / length(c(atl_nodes, bal_nodes))
    maf <- pmin(daf, 1 - daf)
    
    site_fst <- as.numeric(ts$Fst(sample_sets = list(as.integer(atl_nodes), as.integer(bal_nodes)), windows = "sites"))
    valid_indices <- which(maf >= maf_threshold & !is.nan(site_fst))
    
    if (length(valid_indices) > 0) calib_fst[i] <- site_fst[sample(valid_indices, 1)]
  }
  
  # Clean up NA values to calculate metrics and store raw data
  valid_calib_fst <- na.omit(calib_fst)
  
  # Calculate both metrics and bind to results frame
  mean_val <- mean(valid_calib_fst)
  median_val <- median(valid_calib_fst)
  calibration_results <- rbind(calibration_results, data.frame(Years = t_years, Mean_Fst = mean_val, Median_Fst = median_val))
  
  # Store raw values for this time iteration
  raw_fst_list[[as.character(t_years)]] <- data.frame(Years = t_years, Fst = as.numeric(valid_calib_fst))
}

# Combine raw data and save both files
raw_fst_df <- bind_rows(raw_fst_list)
write.csv(calibration_results, file.path(res_dir, "Data_Calibration_Results.csv"), row.names = FALSE)
write.csv(raw_fst_df, file.path(res_dir, "Data_Calibration_Results_Raw.csv"), row.names = FALSE)

# Plotting
p_calib <- ggplot(calibration_results, aes(x = Years)) +
  geom_line(aes(y = Mean_Fst, color = "Mean Fst"), linewidth = 1.2) + 
  geom_point(aes(y = Mean_Fst, color = "Mean Fst"), size = 3) +
  geom_line(aes(y = Median_Fst, color = "Median Fst"), linewidth = 1.2, linetype = "dotted") + 
  geom_point(aes(y = Median_Fst, color = "Median Fst"), size = 3) +
  # Target lines
  geom_hline(yintercept = target_mean_fst, color = "#C44E52", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = target_median_fst, color = "#C44E52", linetype = "dotted", linewidth = 1) +
  scale_color_manual(values = c("Mean Fst" = "#4C72B0", "Median Fst" = "#55A868")) +
  labs(title = "Calibration Curve", 
       x = "Divergence Time (Years)", 
       y = "Simulated per-SNP Fst",
       color = "Metric") +
  theme_minimal() + 
  scale_x_continuous(labels = scales::comma) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Fig_Calibration_Curve.png"), plot = p_calib, width = 8, height = 6)
cat("\n*** Script 2 Complete! Check 'Fig_Calibration_Curve.png' to set your target for Script 3. ***\n")