# ==============================================================================
# SCRIPT 5: Final Unrealistic Ne Simulations (8,000 years fixed)
# ==============================================================================

library(slendr)
library(dplyr)
library(ggplot2)
library(parallel)

init_env()

num_cores <- 8
cat(sprintf("\n*** Running strictly with %d cores ***\n", num_cores))

res_dir <- "/cfs/klemming/projects/snic/naiss2024-6-260/mafalda/simulations/results"
fig_dir <- "/cfs/klemming/projects/snic/naiss2024-6-260/mafalda/simulations/results/figures"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
recomb_rate <- 2.54e-8 
mut_rate <- 2e-9       
seq_length <- 50000  
maf_threshold <- 0.05
n_fragments_main <- 50000  

# Load Empirical Data
cat("\n--- Loading Empirical Data ---\n")
load("/cfs/klemming/projects/snic/naiss2024-6-260/mafalda/simulations/input/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05.RData")
df_spring_filtered <- df_spring_maf0.05 %>% filter(FST >= 0)
names(df_spring_filtered)[3] <- "Fst"
df_spring_filtered$name <- "Observed"
set.seed(42) 
df_spring_filtered_target <- df_spring_filtered[sample(nrow(df_spring_filtered), n_fragments_main), ]
observed_fst_vector <- df_spring_filtered_target$Fst

# ==============================================================================
# CALIBRATED TARGETS
test_Ne_unrealistic <- c(15000, 17500, 20000) 
fixed_split_time <- 8000
fixed_gen_time <- 3
# ==============================================================================

for (target_Ne in test_Ne_unrealistic) {
  cat(sprintf("\n--- Running Final Unrealistic Ne (%d) Simulation ---\n", target_Ne))

  atl <- population("Atlantic", time = fixed_split_time + 5000, N = target_Ne)
  bal <- population("Baltic", parent = atl, time = fixed_split_time, N = target_Ne)
  model <- compile_model(populations = list(atl, bal), generation_time = fixed_gen_time, direction = "backward")
  # Updated to match empirical pool sizes: Atlantic = 482, Baltic = 992
  #sampling_plan <- schedule_sampling(model, times = 0, list(atl, 50), list(bal, 50))
  sampling_plan <- schedule_sampling(model, times = 0, list(atl, 482), list(bal, 992))

  simulate_fragment <- function(i) {
    ts <- msprime(model, sequence_length = seq_length, recombination_rate = recomb_rate, samples = sampling_plan) %>% 
          ts_mutate(mutation_rate = mut_rate)
    if (ts$num_mutations == 0) return(NA)
    
    all_nodes <- ts_nodes(ts)
    
    # Keeping the explicit dplyr:: just to be absolutely safe
    atl_nodes <- all_nodes %>% dplyr::filter(sampled == TRUE, pop == "Atlantic") %>% dplyr::pull(node_id)
    bal_nodes <- all_nodes %>% dplyr::filter(sampled == TRUE, pop == "Baltic") %>% dplyr::pull(node_id)
    
    geno_mat <- ts$genotype_matrix(samples = as.integer(c(atl_nodes, bal_nodes)))
    daf <- rowSums(geno_mat) / length(c(atl_nodes, bal_nodes))
    maf <- pmin(daf, 1 - daf)
    
    site_fst <- as.numeric(ts$Fst(sample_sets = list(as.integer(atl_nodes), as.integer(bal_nodes)), windows = "sites"))
    valid_indices <- which(maf >= maf_threshold & !is.nan(site_fst))
    
    if (length(valid_indices) > 0) return(site_fst[sample(valid_indices, 1)]) else return(NA)
  }

  cat(sprintf("Simulating 50,000 fragments for Ne = %d...\n", target_Ne))
  results_list <- mclapply(1:n_fragments_main, simulate_fragment, mc.cores = num_cores)
  collected_fst <- unlist(results_list)

  final_clean_fst <- as.numeric(na.omit(collected_fst))
  sim_name_label <- sprintf("Unrealistic Ne (%d)", target_Ne)
  final_sim_df <- data.frame(Fst = final_clean_fst, name = sim_name_label)

  # _CorrectN suffix marks files produced with matched empirical sample sizes
  # Old filename commented out to avoid overwriting previous results
  #saveRDS(final_sim_df, file = file.path(res_dir, sprintf("Data_UnrealisticNe_%d.rds", target_Ne)))
  saveRDS(final_sim_df, file = file.path(res_dir, sprintf("Data_UnrealisticNe_%d_CorrectN.rds", target_Ne)))

  # Histogram — _CorrectN suffix added to avoid overwriting old figures
  to_plot_final <- rbind(df_spring_filtered_target[, c("Fst", "name")], final_sim_df)
  p_final_hist <- ggplot(to_plot_final, aes(x = Fst, fill = name)) + 
    geom_histogram(position = "dodge", bins = 50, color = "black", linewidth = 0.5, alpha = 0.8) +
    scale_fill_manual(values = setNames(c("#C44E52", "#4C72B0"), c("Observed", sim_name_label))) +
    labs(title = sprintf("Empirical vs. Unrealistic Ne (Ne: %d, 8k yrs)", target_Ne), x = "Per-SNP Fst", y = "Count") +
    theme_minimal() + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(xlim = c(0, 0.5)) + theme(legend.position = "bottom", legend.title = element_blank())
  
  # _CorrectN suffix added to avoid overwriting old figures
  ggsave(file.path(fig_dir, sprintf("Fig_UnrealisticNe_Hist_%d_CorrectN.png", target_Ne)), plot = p_final_hist, width = 8, height = 6, bg = "white")

  # QQ-Plot — already had _CorrectN suffix, kept as-is
  qq_data <- qqplot(x = final_sim_df$Fst, y = observed_fst_vector, plot.it = FALSE)
  qq_df <- data.frame(Simulated_Quantiles = qq_data$x, Observed_Quantiles = qq_data$y)

  p_final_qq <- ggplot(qq_df, aes(x = Simulated_Quantiles, y = Observed_Quantiles)) +
    geom_point(alpha = 0.5, color = "#4C72B0", size = 1.5) +
    geom_abline(intercept = 0, slope = 1, color = "#C44E52", linetype = "dashed", linewidth = 1) +
    labs(title = sprintf("QQ-Plot: Unrealistic Ne = %d (8k yrs)", target_Ne), x = "Simulated Neutral Fst (Expected)", y = "Observed Fst (Empirical)") +
    theme_minimal() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, sprintf("Fig_UnrealisticNe_QQ_%d_CorrectN.png", target_Ne)), plot = p_final_qq, width = 8, height = 6, bg = "white")
  
  cat(sprintf("-> Completed Ne = %d\n", target_Ne))
}

cat("\n*** Script 5 Complete! All unrealistic Ne targets simulated and plotted. ***\n")