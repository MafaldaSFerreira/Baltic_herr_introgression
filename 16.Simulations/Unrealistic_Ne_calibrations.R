# ==============================================================================
# SCRIPT 4: Calibration for Unrealistic Ne (8,000 years fixed) - HPC VERSION
# ==============================================================================

library(slendr)
library(dplyr)
library(ggplot2)

init_env()

cat("\n*** Running Calibration Locally (Sequential) ***\n")

# Working directory:
res_dir <- "/cfs/klemming/projects/snic/naiss2024-6-260/mafalda/simulations/results"
fig_dir <- "/cfs/klemming/projects/snic/naiss2024-6-260/mafalda/simulations/results/figures"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# Parameters
recomb_rate <- 2.54e-8 
mut_rate <- 2e-9       
seq_length <- 50000  
maf_threshold <- 0.05
n_fragments_calib <- 5000  # Smaller run just to find the Fst sweet spot

# Fixed Parameters
fixed_split_time <- 8000
fixed_gen_time <- 3
target_mean_fst <- 0.039
target_median_fst <- 0.022 # Added target median

# Expanded Unrealistic Ne Grid to test (from 100 to 30,000)
test_Ne_unrealistic <- c(100, 500, 1000, 2500, 5000, 10000, 15000, 20000, 30000)

# Create an empty data frame to store the summary statistics
calib_results <- data.frame(Ne = integer(), Mean_Fst = numeric(), Median_Fst = numeric(), Max_Fst = numeric())

# Create a list to store the raw Fst values
raw_fst_list <- list()

cat(sprintf("\n--- Running Calibration for Fixed Time: %d years ---\n", fixed_split_time))

for (N in test_Ne_unrealistic) {
  cat(sprintf("\nTesting Unrealistic Ne = %d ...\n", N))
  
  atl <- population("Atlantic", time = fixed_split_time + 5000, N = N)
  bal <- population("Baltic", parent = atl, time = fixed_split_time, N = N)
  model <- compile_model(populations = list(atl, bal), generation_time = fixed_gen_time, direction = "backward")
  sampling_plan <- schedule_sampling(model, times = 0, list(atl, 50), list(bal, 50))
  
  simulate_fragment <- function(i) {
    ts <- msprime(model, sequence_length = seq_length, recombination_rate = recomb_rate, samples = sampling_plan) %>% 
      ts_mutate(mutation_rate = mut_rate)
    if (ts$num_mutations == 0) return(NA)
    
    all_nodes <- ts_nodes(ts)
    atl_nodes <- all_nodes %>% dplyr::filter(sampled == TRUE, pop == "Atlantic") %>% dplyr::pull(node_id)
    bal_nodes <- all_nodes %>% dplyr::filter(sampled == TRUE, pop == "Baltic") %>% dplyr::pull(node_id)
    
    geno_mat <- ts$genotype_matrix(samples = as.integer(c(atl_nodes, bal_nodes)))
    daf <- rowSums(geno_mat) / length(c(atl_nodes, bal_nodes))
    maf <- pmin(daf, 1 - daf)
    
    site_fst <- as.numeric(ts$Fst(sample_sets = list(as.integer(atl_nodes), as.integer(bal_nodes)), windows = "sites"))
    valid_indices <- which(maf >= maf_threshold & !is.nan(site_fst))
    
    if (length(valid_indices) > 0) return(site_fst[sample(valid_indices, 1)]) else return(NA)
  }
  
  # Run sequentially on a single core
  results_list <- lapply(1:n_fragments_calib, simulate_fragment)
  collected_fst <- as.numeric(na.omit(unlist(results_list)))
  
  mean_fst <- mean(collected_fst)
  median_fst <- median(collected_fst) # Calculate median
  max_fst <- max(collected_fst)
  
  # Save the summary results for this Ne
  calib_results <- rbind(calib_results, data.frame(Ne = N, Mean_Fst = mean_fst, Median_Fst = median_fst, Max_Fst = max_fst))
  
  # Save the raw Fst values for this Ne
  raw_fst_list[[as.character(N)]] <- data.frame(Ne = N, Fst = collected_fst)
  
  cat(sprintf("-> Results for Ne=%d | Mean Fst: %.4f | Median Fst: %.4f | Max Fst: %.4f\n", N, mean_fst, median_fst, max_fst))
}

# Combine raw data and save both files
raw_fst_df <- bind_rows(raw_fst_list)
saveRDS(calib_results, file="Data_Calibration_UnrealisticNe.rds")
saveRDS(raw_fst_df, file="Data_Calibration_UnrealisticNe_Raw.rds")

# ==============================================================================
# Plot the Results
# ==============================================================================
cat("\n--- Generating Calibration Plot ---\n")

p_calib <- ggplot(calib_results, aes(x = Ne)) +
  geom_line(aes(y = Max_Fst, color = "Max Fst"), linewidth = 1) +
  geom_point(aes(y = Max_Fst, color = "Max Fst"), size = 2) +
  geom_line(aes(y = Mean_Fst, color = "Mean Fst"), linewidth = 1) +
  geom_point(aes(y = Mean_Fst, color = "Mean Fst"), size = 2) +
  geom_line(aes(y = Median_Fst, color = "Median Fst"), linewidth = 1, linetype = "dotted") +
  geom_point(aes(y = Median_Fst, color = "Median Fst"), size = 2) +
  # Add the target Fst lines
  geom_hline(yintercept = target_mean_fst, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = target_median_fst, linetype = "dotted", color = "black", linewidth = 0.8) +
  annotate("text", x = max(test_Ne_unrealistic) * 0.8, y = target_mean_fst + 0.002, 
           label = sprintf("Target Mean = %.3f", target_mean_fst), color = "black", fontface = "bold") +
  annotate("text", x = max(test_Ne_unrealistic) * 0.8, y = target_median_fst - 0.002, 
           label = sprintf("Target Median = %.3f", target_median_fst), color = "black", fontface = "bold") +
  scale_color_manual(values = c("Max Fst" = "#C44E52", "Mean Fst" = "#4C72B0", "Median Fst" = "#55A868")) +
  labs(title = "Calibration: Fst vs. Bottleneck Size (Ne)",
       subtitle = "Fixed Divergence Time: 8,000 years",
       x = "Effective Population Size (Ne)",
       y = "Simulated Fst",
       color = "Metric") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot
ggsave(file.path(fig_dir, "Fig_UnrealisticNe_Calibration.png"), plot = p_calib, width = 8, height = 6, bg = "white")

cat("\n*** Calibration Complete! Plot saved as 'Fig_UnrealisticNe_Calibration.png' ***\n")