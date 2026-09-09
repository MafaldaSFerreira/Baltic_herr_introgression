# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
# Resolved to an absolute path before we setwd() below, so it stays valid
# regardless of the working directory.
FIGSHARE_ROOT <- normalizePath("../Figshare")

# Raw simulation output from the Dardel cluster run -- not part of the
# Figshare deposit. Set this to wherever your own simulation results live.
setwd("<path_to_your_simulation_results>")

# ==============================================================================
# Plot the Results
# ==============================================================================
cat("\n--- Generating Calibration Plot ---\n")

p_calib <- ggplot(calib_results, aes(x = Ne)) +
  geom_line(aes(y = Max_Fst, color = "Max Fst"), linewidth = 1) +
  geom_point(aes(y = Max_Fst, color = "Max Fst"), size = 2) +
  geom_line(aes(y = Mean_Fst, color = "Mean Fst"), linewidth = 1) +
  geom_point(aes(y = Mean_Fst, color = "Mean Fst"), size = 2) +
  # Add the target Fst line
  geom_hline(yintercept = target_fst, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = max(test_Ne_unrealistic) * 0.8, y = target_fst + 0.002, 
           label = sprintf("Target Fst = %.3f", target_fst), color = "black", fontface = "bold") +
  scale_color_manual(values = c("Max Fst" = "#C44E52", "Mean Fst" = "#4C72B0")) +
  labs(title = "Calibration: Fst vs. Bottleneck Size (Ne)",
       subtitle = "Fixed Divergence Time: 8,000 years",
       x = "Effective Population Size (Ne)",
       y = "Simulated Fst",
       color = "Metric") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot
ggsave(file.path("./2026-02-22_dardel/figures/Fig_UnrealisticNe_Calibration.png"), plot = p_calib, width = 8, height = 6, bg = "white")

cat("\n*** Calibration Complete! Plot saved as 'Fig_UnrealisticNe_Calibration.png' ***\n")


#### Extended Data Figure ####
# ==============================================================================
# Extended Data Figure: Neutral drift cannot explain observed FST between
# Atlantic and Baltic spring-spawning herring
#
# Produces:
#   1. Main Extended Data figure (4 panels):
#      a: Calibration curve — divergence time vs mean + median FST
#      b: Calibration curve — Ne (bottleneck) vs mean + median FST
#      c: QQ-plot — most conservative Ne (Ne = 20,000, 8k years)
#      d: QQ-plot — most conservative divergence time (475k years)
#
#   2. Supplementary figure — FST frequency distributions:
#      Overlaid density plots comparing simulated neutral FST vs observed
#      (genome-wide and introgressed regions) for each scenario
#
#   3. Supplementary QQ-plots — remaining Ne and divergence time scenarios
#
# Input files:
#   Data_Calibration_Results.csv
#   Data_Calibration_UnrealisticNe.rds
#   Data_ImpossibleTime_425k.rds
#   Data_ImpossibleTime_475k.rds
#   Data_UnrealisticNe_15000.rds
#   Data_UnrealisticNe_17500.rds
#   Data_UnrealisticNe_20000.rds
#
# ==============================================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
library(tidyr)

# ------------------------------------------------------------------------------
# 0. PATHS AND EMPIRICAL VALUES
# ------------------------------------------------------------------------------

DATA_DIR   <- "."
OUTPUT_DIR <- "./figures"

EMPIRICAL_DATA_PATH <- file.path(FIGSHARE_ROOT, "5.processing_pool_seq_data/fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05.RData")
EMPIRICAL_DATA_PATH_INTROGRESSION <- file.path(FIGSHARE_ROOT, "5.processing_pool_seq_data/fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05_introgression_regions.RData")

# Empirical FST summary statistics
FST_GENOME_MEAN   <- 0.0391534
FST_GENOME_MEDIAN <- 0.02199862
FST_INTRO_MEAN    <- 0.2125821
FST_INTRO_MEDIAN  <- 0.1676966

# Realistic Ne range for herring
NE_REALISTIC_LOW  <- 100000
NE_REALISTIC_HIGH <- 1000000

# Baltic Sea opening
BALTIC_OPEN_YRS <- 8000

# Number of windows to sample from empirical data for QQ quantile matching
N_FRAGMENTS_MAIN <- 50000

# ------------------------------------------------------------------------------
# 0.5. LOAD EMPIRICAL FST DATA
# ------------------------------------------------------------------------------

message("Loading empirical FST data...")
load(EMPIRICAL_DATA_PATH)
load(EMPIRICAL_DATA_PATH_INTROGRESSION)

# Genome-wide FST (all windows, MAF-filtered)
df_genome <- df_spring_maf0.05 %>%
  filter(FST >= 0) %>%
  rename(Fst = FST)

# Introgressed-region FST
df_intro <- overlaps_spring_df %>%
  filter(FST >= 0) %>%
  rename(Fst = FST)

# Sample N_FRAGMENTS_MAIN windows for QQ-plot quantile matching
set.seed(42)
observed_fst_vector <- df_genome %>%
  slice_sample(n = N_FRAGMENTS_MAIN) %>%
  pull(Fst)

# ------------------------------------------------------------------------------
# 1. SHARED THEME AND COLOUR PALETTE
# ------------------------------------------------------------------------------

theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0,
                                      colour = "grey40"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}

COL_MEAN       <- "#2166ac"  # blue   — mean FST simulation line
COL_MEDIAN     <- "#5e9c6a"  # green  — median FST simulation line
COL_MAX        <- "#d6604d"  # red    — max FST simulation line
COL_GWIDE_REF  <- "#d6604d"  # red    — genome-wide FST reference line
COL_INTRO_REF  <- "#4dac26"  # green  — introgressed FST reference line
COL_SIM_DENS   <- "#2166ac"  # blue   — simulated density fill
COL_GWIDE_DENS <- "#969696"  # grey   — genome-wide observed density fill
COL_INTRO_DENS <- "#fc8d59"  # orange — introgressed region density fill

# ------------------------------------------------------------------------------
# 2. PANEL A: Divergence time calibration — mean AND median
# ------------------------------------------------------------------------------

cal_time <- read_csv(file.path(DATA_DIR, "Data_Calibration_Results.csv"),
                     show_col_types = FALSE)
# Expected columns: Years, Mean_Fst, Median_Fst
# If Median_Fst is absent, see NOTE at bottom of script for how to add it

cal_time_long <- cal_time %>%
  pivot_longer(cols = c(Mean_Fst, Median_Fst),
               names_to  = "metric",
               values_to = "fst") %>%
  mutate(metric = recode(metric,
                         "Mean_Fst"   = "Mean FST",
                         "Median_Fst" = "Median FST"))

# Y-axis range for positioning annotations
y_max_a <- max(cal_time_long$fst, na.rm = TRUE)
y_step  <- y_max_a * 0.025

panel_A <- ggplot(cal_time_long, aes(x = Years, y = fst,
                                     colour = metric, group = metric)) +
  # Genome-wide reference lines (mean = dashed, median = dotdash)
  geom_hline(yintercept = FST_GENOME_MEAN,
             linetype = "dashed",  colour = COL_GWIDE_REF, linewidth = 0.6) +
  geom_hline(yintercept = FST_GENOME_MEDIAN,
             linetype = "dotdash", colour = COL_GWIDE_REF, linewidth = 0.6) +
  # Introgressed-region reference lines
  geom_hline(yintercept = FST_INTRO_MEAN,
             linetype = "dashed",  colour = COL_INTRO_REF, linewidth = 0.6) +
  geom_hline(yintercept = FST_INTRO_MEDIAN,
             linetype = "dotdash", colour = COL_INTRO_REF, linewidth = 0.6) +
  # Annotations
  annotate("text", x = max(cal_time$Years) * 0.97,
           y = FST_GENOME_MEAN + y_step,
           label = sprintf("Genome-wide mean = %.3f", FST_GENOME_MEAN),
           hjust = 1, size = 2.6, colour = COL_GWIDE_REF) +
  annotate("text", x = max(cal_time$Years) * 0.97,
           y = FST_GENOME_MEDIAN - y_step,
           label = sprintf("Genome-wide median = %.3f", FST_GENOME_MEDIAN),
           hjust = 1, size = 2.6, colour = COL_GWIDE_REF) +
  annotate("text", x = max(cal_time$Years) * 0.97,
           y = FST_INTRO_MEAN + y_step,
           label = sprintf("Introgressed mean = %.3f", FST_INTRO_MEAN),
           hjust = 1, size = 2.6, colour = COL_INTRO_REF) +
  annotate("text", x = max(cal_time$Years) * 0.97,
           y = FST_INTRO_MEDIAN - y_step,
           label = sprintf("Introgressed median = %.3f", FST_INTRO_MEDIAN),
           hjust = 1, size = 2.6, colour = COL_INTRO_REF) +
  # Baltic opening
  geom_vline(xintercept = BALTIC_OPEN_YRS,
             linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  annotate("text",
           x = BALTIC_OPEN_YRS + max(cal_time$Years) * 0.01,
           y = y_max_a * 0.90,
           label = "Baltic\nopening",
           hjust = 0, size = 2.6, colour = "grey40") +
  # Data
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Mean FST"   = COL_MEAN,
                                 "Median FST" = COL_MEDIAN)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x      = "Divergence time (years)",
    y      = expression(bold("Simulated F")[bold("ST")]),
    colour = NULL,
    title  = "a"
  ) +
  theme_publication() +
  theme(legend.position = c(0.18, 0.88))

# ------------------------------------------------------------------------------
# 3. PANEL B: Ne calibration — mean, median, AND max
# ------------------------------------------------------------------------------

cal_ne <- readRDS(file.path(DATA_DIR, "Data_Calibration_UnrealisticNe.rds"))
# Expected columns: Ne, Mean_Fst, Max_Fst
# Median_Fst will be plotted if present; see NOTE at bottom for how to add it

pivot_cols <- intersect(c("Mean_Fst", "Median_Fst", "Max_Fst"), colnames(cal_ne))

cal_ne_long <- cal_ne %>%
  pivot_longer(cols = all_of(pivot_cols),
               names_to  = "metric",
               values_to = "fst") %>%
  mutate(metric = recode(metric,
                         "Mean_Fst"   = "Mean FST",
                         "Median_Fst" = "Median FST",
                         "Max_Fst"    = "Max FST"))

colour_vals_ne <- c("Mean FST"   = COL_MEAN,
                    "Median FST" = COL_MEDIAN,
                    "Max FST"    = COL_MAX)

y_max_b <- max(cal_ne_long$fst, na.rm = TRUE)

panel_B <- ggplot(cal_ne_long, aes(x = Ne, y = fst,
                                   colour = metric, group = metric)) +
  # Reference lines
  geom_hline(yintercept = FST_GENOME_MEAN,
             linetype = "dashed",  colour = COL_GWIDE_REF, linewidth = 0.6) +
  geom_hline(yintercept = FST_GENOME_MEDIAN,
             linetype = "dotdash", colour = COL_GWIDE_REF, linewidth = 0.6) +
  annotate("text", x = max(cal_ne$Ne) * 0.97,
           y = FST_GENOME_MEAN + y_max_b * 0.025,
           label = sprintf("Genome-wide mean = %.3f", FST_GENOME_MEAN),
           hjust = 1, size = 2.6, colour = COL_GWIDE_REF) +
  annotate("text", x = max(cal_ne$Ne) * 0.97,
           y = FST_GENOME_MEDIAN - y_max_b * 0.03,
           label = sprintf("Genome-wide median = %.3f", FST_GENOME_MEDIAN),
           hjust = 1, size = 2.6, colour = COL_GWIDE_REF) +
  # Realistic Ne shaded region
  annotate("rect",
           xmin = NE_REALISTIC_LOW, xmax = NE_REALISTIC_HIGH,
           ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  annotate("text",
           x = sqrt(NE_REALISTIC_LOW * NE_REALISTIC_HIGH),
           y = y_max_b * 0.94,
           label = "Realistic Ne\nfor herring",
           size = 2.6, colour = "grey40") +
  # Data
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = colour_vals_ne) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x      = "Effective population size (Ne)",
    y      = expression(bold("Simulated F")[bold("ST")]),
    colour = NULL,
    title  = "b"
  ) +
  theme_publication() +
  theme(legend.position = c(0.75, 0.85))

# ------------------------------------------------------------------------------
# 4. QQ-PLOT HELPER
# ------------------------------------------------------------------------------

make_qq_plot <- function(rds_path, panel_label, obs_fst_vector,
                         plot_subtitle = NULL) {
  
  dat     <- readRDS(rds_path)
  sim_fst <- sort(dat$Fst)
  obs_fst <- sort(obs_fst_vector)
  
  probs <- seq(0, 1, length.out = min(length(sim_fst), length(obs_fst)))
  qq_df <- data.frame(
    simulated = quantile(sim_fst, probs = probs, na.rm = TRUE),
    observed  = quantile(obs_fst, probs = probs, na.rm = TRUE)
  )
  
  ggplot(qq_df, aes(x = simulated, y = observed)) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = COL_GWIDE_REF, linewidth = 0.7) +
    geom_point(colour = COL_MEAN, size = 0.8, alpha = 0.6) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x        = expression(bold("Simulated neutral F")[bold("ST")]),
      y        = expression(bold("Observed F")[bold("ST")]),
      title    = panel_label,
      subtitle = plot_subtitle
    ) +
    theme_publication() +
    theme(plot.subtitle = element_text(size = 9, colour = "grey30"))
}

# ------------------------------------------------------------------------------
# 5. PANELS C AND D: QQ-plots (most conservative scenarios)
# ------------------------------------------------------------------------------

panel_C <- make_qq_plot(
  rds_path       = file.path(DATA_DIR, "Data_UnrealisticNe_20000_CorrectN.rds"),
  panel_label    = "c",
  obs_fst_vector = observed_fst_vector,
  plot_subtitle  = "Most conservative Ne (Ne = 20,000 | 8,000 years)"
)

panel_D <- make_qq_plot(
  rds_path       = file.path(DATA_DIR, "Data_ImpossibleTime_475k_CorrectN.rds"),
  panel_label    = "d",
  obs_fst_vector = observed_fst_vector,
  plot_subtitle  = "Most conservative divergence time (475,000 years)"
)

# ------------------------------------------------------------------------------
# 6. ASSEMBLE AND SAVE MAIN EXTENDED DATA FIGURE
# ------------------------------------------------------------------------------

figure_extended_data <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_annotation(
    caption = paste0(
      "Simulations performed with msprime via slendr. ",
      "Red dashed/dotdash lines: observed genome-wide mean/median FST ",
      "(Atlantic vs Baltic spring-spawning herring). ",
      "Green dashed/dotdash lines: mean/median FST at introgressed regions. ",
      "Grey shaded region in b: realistic Ne range for Atlantic herring ",
      "(Ne = ", scales::comma(NE_REALISTIC_LOW), "\u2013",
      scales::comma(NE_REALISTIC_HIGH), "). ",
      "QQ-plots compare the empirical genome-wide FST distribution to neutral ",
      "expectations; points above the 1:1 line indicate excess differentiation."
    ),
    theme = theme(plot.caption = element_text(size = 7.5, colour = "grey40",
                                              hjust = 0, lineheight = 1.3))
  )

ggsave(file.path(OUTPUT_DIR, "ExtendedData_SimulationFST_CorrectNe.pdf"),
       plot = figure_extended_data,
       width = 20, height = 26, units = "cm", device = "pdf", bg = "white")

ggsave(file.path(OUTPUT_DIR, "ExtendedData_SimulationFST_CorrectNe.png"),
       plot = figure_extended_data,
       width = 20, height = 26, units = "cm", dpi = 300, bg = "white")

message("Main Extended Data figure saved.")

# ------------------------------------------------------------------------------
# 7. DENSITY OVERLAY HELPER
#    PI suggestion: show how the simulated neutral FST distribution compares to
#    both the genome-wide and introgressed-region observed distributions.
#    Dashed vertical lines mark group medians (less sensitive to outliers than
#    means, and directly comparable to the median-based calibration in panels
#    A and B).
# ------------------------------------------------------------------------------

make_density_plot <- function(rds_path, panel_label, plot_subtitle,
                              obs_genome_df, obs_intro_df) {
  
  dat <- readRDS(rds_path)
  
  plot_df <- bind_rows(
    data.frame(Fst = dat$Fst,           group = "Simulated neutral"),
    data.frame(Fst = obs_genome_df$Fst, group = "Observed genome-wide"),
    data.frame(Fst = obs_intro_df$Fst,  group = "Observed introgressed")
  ) %>%
    mutate(group = factor(group,
                          levels = c("Simulated neutral",
                                     "Observed genome-wide",
                                     "Observed introgressed")))
  
  medians_df <- plot_df %>%
    group_by(group) %>%
    summarise(med = median(Fst, na.rm = TRUE), .groups = "drop")
  
  ggplot(plot_df, aes(x = Fst, fill = group, colour = group)) +
    geom_density(alpha = 0.35, linewidth = 0.5, adjust = 1.5) +
    geom_vline(data = medians_df,
               aes(xintercept = med, colour = group),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = c(
      "Simulated neutral"     = COL_SIM_DENS,
      "Observed genome-wide"  = COL_GWIDE_DENS,
      "Observed introgressed" = COL_INTRO_DENS
    )) +
    scale_colour_manual(values = c(
      "Simulated neutral"     = COL_SIM_DENS,
      "Observed genome-wide"  = COL_GWIDE_DENS,
      "Observed introgressed" = COL_INTRO_DENS
    )) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x        = expression(bold("F")[bold("ST")]),
      y        = "Density",
      fill     = NULL,
      colour   = NULL,
      title    = panel_label,
      subtitle = plot_subtitle
    ) +
    theme_publication() +
    theme(
      legend.position  = c(0.65, 0.80),
      legend.key.size  = unit(0.35, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    )
}

# ------------------------------------------------------------------------------
# 8. BUILD SUPPLEMENTARY DISTRIBUTION FIGURE (4 panels)
# ------------------------------------------------------------------------------

dist_scenarios <- list(
  list(path  = "Data_UnrealisticNe_20000_CorrectN.rds",   label = "a",
       sub   = "Ne = 20,000 | 8,000 years (most conservative Ne)"),
  list(path  = "Data_ImpossibleTime_475k_CorrectN.rds",   label = "b",
       sub   = "475,000 years divergence (most conservative time)"),
  list(path  = "Data_UnrealisticNe_15000_CorrectN.rds",   label = "c",
       sub   = "Ne = 15,000 | 8,000 years"),
  list(path  = "Data_ImpossibleTime_425k_CorrectN.rds",   label = "d",
       sub   = "425,000 years divergence")
)

dist_plots <- lapply(dist_scenarios, function(x) {
  make_density_plot(
    rds_path      = file.path(DATA_DIR, x$path),
    panel_label   = x$label,
    plot_subtitle = x$sub,
    obs_genome_df = df_genome,
    obs_intro_df  = df_intro
  )
})

figure_distributions <- (dist_plots[[1]] | dist_plots[[2]]) /
  (dist_plots[[3]] | dist_plots[[4]]) +
  plot_annotation(
    title   = "Supplementary Figure XX \u2014 FST frequency distributions: neutral simulations vs observed",
    caption = paste0(
      "Density plots comparing FST distributions under each neutral simulation ",
      "scenario (blue) with the observed genome-wide FST distribution (grey) ",
      "and FST at introgressed regions (orange). ",
      "Dashed vertical lines indicate group medians. ",
      "Even under the most conservative (least realistic) demographic parameters, ",
      "the simulated neutral distribution cannot account for the elevated FST ",
      "values observed at introgressed regions."
    ),
    theme = theme(
      plot.title   = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 7.5, colour = "grey40",
                                  hjust = 0, lineheight = 1.3)
    )
  )

ggsave(file.path(OUTPUT_DIR, "Supplementary_FstDistributions_CorrectN.pdf"),
       plot = figure_distributions,
       width = 18, height = 16, units = "cm", device = "pdf", bg = "white")

ggsave(file.path(OUTPUT_DIR, "Supplementary_FstDistributions_CorrectN.png"),
       plot = figure_distributions,
       width = 18, height = 16, units = "cm", dpi = 300, bg = "white")

message("Supplementary FST distribution figure saved.")

# ------------------------------------------------------------------------------
# 9. SUPPLEMENTARY QQ-PLOTS (remaining scenarios)
# ------------------------------------------------------------------------------

supp_qq_files <- list(
  list(path = "Data_UnrealisticNe_15000_CorrectN.rds", label = "a",
       sub  = "Ne = 15,000 | 8,000 years"),
  list(path = "Data_UnrealisticNe_17500_CorrectN.rds", label = "b",
       sub  = "Ne = 17,500 | 8,000 years"),
  list(path = "Data_ImpossibleTime_425k_CorrectN.rds", label = "c",
       sub  = "425,000 years divergence")
)

supp_qq_plots <- lapply(supp_qq_files, function(x) {
  make_qq_plot(
    rds_path       = file.path(DATA_DIR, x$path),
    panel_label    = x$label,
    obs_fst_vector = observed_fst_vector,
    plot_subtitle  = x$sub
  )
})

figure_supp_qq <- (supp_qq_plots[[1]] | supp_qq_plots[[2]] | supp_qq_plots[[3]]) +
  plot_annotation(
    title   = "Supplementary Figure XX \u2014 Additional QQ-plots for neutral drift simulations",
    caption = paste0(
      "QQ-plots for Ne = 15,000, Ne = 17,500 (8,000 years divergence), and ",
      "divergence time = 425,000 years. ",
      "See Extended Data Figure XX for the most conservative scenarios."
    ),
    theme = theme(
      plot.title   = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 7.5, colour = "grey40", hjust = 0)
    )
  )

ggsave(file.path(OUTPUT_DIR, "Supplementary_QQplots_additional_CorrectN.pdf"),
       plot = figure_supp_qq,
       width = 24, height = 9, units = "cm", device = "pdf", bg = "white")

ggsave(file.path(OUTPUT_DIR, "Supplementary_QQplots_additional_CorrectN.png"),
       plot = figure_supp_qq,
       width = 24, height = 9, units = "cm", dpi = 300, bg = "white")

message("Supplementary QQ-plot figure saved.")
message("All outputs written to: ", OUTPUT_DIR)

# ==============================================================================
# NOTE: Adding Median_Fst to simulation output
# ==============================================================================
# Panels A and B plot Median_Fst alongside Mean_Fst. Panel B uses
# pivot_longer with intersect() so it degrades gracefully if Median_Fst is
# absent — but you should add it by updating the simulation calibration loops.
#
# In the Ne calibration loop (Script 4), replace:
#   calib_results <- rbind(calib_results,
#     data.frame(Ne = N, Mean_Fst = mean_fst, Max_Fst = max_fst))
# with:
#   calib_results <- rbind(calib_results,
#     data.frame(Ne = N,
#                Mean_Fst   = mean(collected_fst),
#                Median_Fst = median(collected_fst),
#                Max_Fst    = max(collected_fst)))
#
# Apply the same change to the divergence-time calibration loop, and add a
# Median_Fst column to Data_Calibration_Results.csv.
# ==============================================================================


# ==============================================================================
# SUPPLEMENTARY FIGURE: ECDF plots — FST distributions
# Copy-paste this block at the end of plot_simulations_extended_data_v2.R
#
# Replaces the density overlay figure with ECDF plots, which avoid the
# sample-size distortion of kernel density estimates. The vertical gap between
# curves at any FST value reads directly as the difference in the fraction of
# windows exceeding that value — i.e. how many introgressed windows have FST
# values that are essentially impossible under the neutral expectation.
#
# Requires objects already created earlier in the script:
#   df_genome          — data frame with column Fst (genome-wide, MAF-filtered)
#   df_intro           — data frame with column Fst (introgressed regions)
#   DATA_DIR, OUTPUT_DIR, theme_publication(), COL_* colour objects
# ==============================================================================

# ------------------------------------------------------------------------------
# ECDF HELPER FUNCTION
# ------------------------------------------------------------------------------

make_ecdf_plot <- function(rds_path, panel_label, plot_subtitle,
                           obs_genome_df, obs_intro_df) {
  
  dat <- readRDS(rds_path)
  
  # Build combined data frame
  # We downsample genome-wide to 50k for visual clarity (ECDF is stable at
  # this n); introgressed windows are kept in full because there are few of them
  set.seed(42)
  n_genome_sample <- min(nrow(obs_genome_df), 50000)
  
  plot_df <- bind_rows(
    data.frame(Fst   = dat$Fst,
               group = "Simulated neutral"),
    data.frame(Fst   = obs_genome_df %>%
                 slice_sample(n = n_genome_sample) %>%
                 pull(Fst),
               group = "Observed genome-wide"),
    data.frame(Fst   = obs_intro_df$Fst,
               group = "Observed introgressed")
  ) %>%
    mutate(group = factor(group,
                          levels = c("Simulated neutral",
                                     "Observed genome-wide",
                                     "Observed introgressed")))
  
  # Colours and linetypes per group
  ecdf_colours   <- c("Simulated neutral"     = COL_SIM_DENS,
                      "Observed genome-wide"  = COL_GWIDE_DENS,
                      "Observed introgressed" = COL_INTRO_DENS)
  ecdf_linetypes <- c("Simulated neutral"     = "solid",
                      "Observed genome-wide"  = "solid",
                      "Observed introgressed" = "solid")
  ecdf_linewidth  <- c("Simulated neutral"     = 0.7,
                       "Observed genome-wide"  = 0.7,
                       "Observed introgressed" = 1.0)   # slightly thicker
  
  ggplot(plot_df, aes(x = Fst, colour = group, linewidth = group)) +
    stat_ecdf(geom = "step", pad = FALSE) +
    # Rug along the top for introgressed raw values — shows individual windows
    # geom_rug(data = filter(plot_df, group == "Observed introgressed"),
    #          aes(x = Fst),
    #          colour = COL_INTRO_DENS,
    #          alpha  = 0.6,
    #          length = unit(0.04, "npc"),
    #          sides  = "t",
    #          inherit.aes = FALSE) +
    scale_colour_manual(values = ecdf_colours) +
    scale_linewidth_manual(values = ecdf_linewidth) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0),
                       breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0),
                       labels = scales::percent_format(accuracy = 1)) +
    # Vertical reference lines for empirical medians
    geom_vline(xintercept = FST_GENOME_MEDIAN,
               linetype = "dashed", colour = COL_GWIDE_DENS,
               linewidth = 0.45, alpha = 0.8) +
    geom_vline(xintercept = FST_INTRO_MEDIAN,
               linetype = "dashed", colour = COL_INTRO_DENS,
               linewidth = 0.45, alpha = 0.8) +
    labs(
      x        = expression(bold("F")[bold("ST")]),
      y        = "Cumulative proportion of windows",
      colour   = NULL,
      linewidth = NULL,
      title    = panel_label,
      subtitle = plot_subtitle
    ) +
    theme_publication() +
    theme(
      legend.position  = c(0.62, 0.25),
      legend.key.size  = unit(0.4,  "cm"),
      legend.spacing.y = unit(0.08, "cm")
    ) +
    guides(linewidth = "none")   # suppress redundant linewidth legend
}

# ------------------------------------------------------------------------------
# BUILD THE 4-PANEL ECDF FIGURE
# Same scenario layout as the density figure it replaces
# ------------------------------------------------------------------------------

ecdf_scenarios <- list(
  list(path  = "Data_UnrealisticNe_20000_CorrectN.rds",  label = "a",
       sub   = "Most conservative Ne (Ne = 20,000 | 8,000 years)"),
  list(path  = "Data_ImpossibleTime_475k_CorrectN.rds",  label = "b",
       sub   = "Most conservative divergence time (475,000 years)"),
  list(path  = "Data_UnrealisticNe_15000_CorrectN.rds",  label = "c",
       sub   = "Ne = 15,000 | 8,000 years"),
  list(path  = "Data_ImpossibleTime_425k_CorrectN.rds",  label = "d",
       sub   = "425,000 years divergence"),
  list(path  = "Data_UnrealisticNe_17500_CorrectN.rds",  label = "e",
       sub   = "Ne = 17,500 | 8,000 years")
)

ecdf_plots <- lapply(ecdf_scenarios, function(x) {
  make_ecdf_plot(
    rds_path      = file.path(DATA_DIR, x$path),
    panel_label   = x$label,
    plot_subtitle = x$sub,
    obs_genome_df = df_genome,
    obs_intro_df  = df_intro
  )
})

figure_ecdf <- (ecdf_plots[[1]] | ecdf_plots[[2]]) /
  (ecdf_plots[[3]] | ecdf_plots[[4]]) /
  (ecdf_plots[[5]] | plot_spacer()) +
  plot_annotation(
    title   = "Supplementary Figure XX \u2014 FST cumulative distributions: neutral simulations vs observed",
    caption = paste0(
      "ECDF plots comparing the FST distribution under each neutral simulation scenario (blue) ",
      "with the observed genome-wide FST distribution (grey) and FST at introgressed regions (orange). ",
      "Tick marks along the top of each panel show the raw FST values of individual introgressed windows. ",
      "Dashed vertical lines indicate observed genome-wide and introgressed-region medians. ",
      "The rightward shift of the introgressed ECDF relative to the neutral simulation ",
      "shows that introgressed-region FST values occupy the extreme tail of, or lie entirely ",
      "outside, the neutral expectation under any biologically plausible demographic scenario."
    ),
    theme = theme(
      plot.title   = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 7.5, colour = "grey40",
                                  hjust = 0, lineheight = 1.3)
    )
  )

ggsave(file.path(OUTPUT_DIR, "Supplementary_FstECDF_CorrectN.pdf"),
       plot = figure_ecdf,
       width = 200, height = 260, units = "mm", device = "pdf", bg = "white")

ggsave(file.path(OUTPUT_DIR, "Supplementary_FstECDF_CorrectN.png"),
       plot = figure_ecdf,
       width = 200, height = 260, units = "mm", dpi = 300, bg = "white")

message("Supplementary ECDF figure saved.")



# ==============================================================================
# Histogram comparison: observed genome-wide FST vs simulated neutral FST
# Five scenarios, two versions:
#   1. Full range (FST 0 to 1)
#   2. Zoom-in (FST 0.5 to 1)
#
# Requires objects already created in plot_simulations_extended_data_v2.R:
#   df_genome          — data frame with column Fst (genome-wide, MAF-filtered)
#   DATA_DIR, OUTPUT_DIR, theme_publication(), COL_* colour objects
#
# NOTE: introgressed-region distribution is intentionally excluded here.
# ==============================================================================

# ------------------------------------------------------------------------------
# HISTOGRAM HELPER FUNCTION
# ------------------------------------------------------------------------------

make_histogram_plot <- function(rds_path, panel_label, plot_subtitle,
                                obs_genome_df,
                                x_min = 0, x_max = 1,
                                n_bins = 100) {
  
  dat <- readRDS(rds_path)
  
  # Downsample genome-wide to match simulation size for visual comparability
  set.seed(42)
  n_sample <- min(nrow(obs_genome_df), length(dat$Fst))
  
  plot_df <- bind_rows(
    data.frame(Fst   = dat$Fst,
               group = "Simulated neutral"),
    data.frame(Fst   = obs_genome_df %>%
                 slice_sample(n = n_sample) %>%
                 pull(Fst),
               group = "Observed genome-wide")
  ) %>%
    mutate(group = factor(group,
                          levels = c("Simulated neutral",
                                     "Observed genome-wide")))
  
  # Bin width computed over the full 0-1 range so bins are consistent
  # across full and zoomed versions
  binwidth <- 1 / n_bins
  
  ggplot(plot_df, aes(x = Fst, fill = group, colour = group)) +
    geom_histogram(binwidth  = binwidth,
                   boundary  = 0,
                   alpha     = 0.5,
                   position  = "identity",  # overlay, not stack
                   linewidth = 0.2) +
    scale_fill_manual(values = c("Simulated neutral"    = COL_SIM_DENS,
                                 "Observed genome-wide" = COL_GWIDE_DENS)) +
    scale_colour_manual(values = c("Simulated neutral"    = COL_SIM_DENS,
                                   "Observed genome-wide" = COL_GWIDE_DENS)) +
    scale_x_continuous(limits  = c(x_min, x_max),
                       expand  = c(0.01, 0),
                       breaks  = seq(x_min, x_max, by = 0.1)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x        = expression(bold("F")[bold("ST")]),
      y        = "Number of windows",
      fill     = NULL,
      colour   = NULL,
      title    = panel_label,
      subtitle = plot_subtitle
    ) +
    theme_publication() +
    theme(
      legend.position  = c(0.70, 0.85),
      legend.key.size  = unit(0.4,  "cm"),
      legend.spacing.y = unit(0.08, "cm")
    )
}

# ------------------------------------------------------------------------------
# SCENARIO LIST (all five)
# ------------------------------------------------------------------------------

hist_scenarios <- list(
  list(path  = "Data_UnrealisticNe_20000_CorrectN.rds",  label = "a",
       sub   = "Most conservative Ne (Ne = 20,000 | 8,000 years)"),
  list(path  = "Data_ImpossibleTime_475k_CorrectN.rds",  label = "b",
       sub   = "Most conservative divergence time (475,000 years)"),
  list(path  = "Data_UnrealisticNe_15000_CorrectN.rds",  label = "c",
       sub   = "Ne = 15,000 | 8,000 years"),
  list(path  = "Data_ImpossibleTime_425k_CorrectN.rds",  label = "d",
       sub   = "425,000 years divergence"),
  list(path  = "Data_UnrealisticNe_17500_CorrectN.rds",  label = "e",
       sub   = "Ne = 17,500 | 8,000 years")
)

# ------------------------------------------------------------------------------
# FIGURE 1: Full range (FST 0 to 1)
# ------------------------------------------------------------------------------

hist_plots_full <- lapply(hist_scenarios, function(x) {
  make_histogram_plot(
    rds_path      = file.path(DATA_DIR, x$path),
    panel_label   = x$label,
    plot_subtitle = x$sub,
    obs_genome_df = df_genome,
    x_min = 0, x_max = 1
  )
})

figure_hist_full <- (hist_plots_full[[1]] | hist_plots_full[[2]]) /
  (hist_plots_full[[3]] | hist_plots_full[[4]]) /
  (hist_plots_full[[5]] | plot_spacer()) +
  plot_annotation(
    title   = "FST distributions: neutral simulations vs observed genome-wide (full range)",
    caption = paste0(
      "Histograms (100 bins, bin width = 0.01) comparing the FST distribution ",
      "under each neutral simulation scenario (blue) with the observed ",
      "genome-wide FST distribution (grey; downsampled to match simulation n). ",
      "Distributions are overlaid (not stacked). ",
      "Introgressed-region windows are excluded."
    ),
    theme = theme(
      plot.title   = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 7.5, colour = "grey40",
                                  hjust = 0, lineheight = 1.3)
    )
  )

ggsave(file.path(OUTPUT_DIR, "Histogram_FstFull_CorrectN.pdf"),
       plot = figure_hist_full,
       width = 18, height = 24, units = "cm", device = "pdf", bg = "white")

ggsave(file.path(OUTPUT_DIR, "Histogram_FstFull_CorrectN.png"),
       plot = figure_hist_full,
       width = 18, height = 24, units = "cm", dpi = 300, bg = "white")

message("Full-range histogram figure saved.")

# ------------------------------------------------------------------------------
# FIGURE 2: Zoom-in (FST 0.5 to 1)
# ------------------------------------------------------------------------------

hist_plots_zoom <- lapply(hist_scenarios, function(x) {
  make_histogram_plot(
    rds_path      = file.path(DATA_DIR, x$path),
    panel_label   = x$label,
    plot_subtitle = x$sub,
    obs_genome_df = df_genome,
    x_min = 0.5, x_max = 1
  )
})

figure_hist_zoom <- (hist_plots_zoom[[1]] | hist_plots_zoom[[2]]) /
  (hist_plots_zoom[[3]] | hist_plots_zoom[[4]]) /
  (hist_plots_zoom[[5]] | plot_spacer()) +
  plot_annotation(
    title   = "FST distributions: neutral simulations vs observed genome-wide (FST 0.5\u20131.0)",
    caption = paste0(
      "Same as above but restricted to FST \u2265 0.5 to show the high-differentiation tail. ",
      "The observed genome-wide distribution (grey) retains windows reaching FST = 1, ",
      "while the neutral simulations produce far fewer such extreme values. ",
      "Introgressed-region windows are excluded."
    ),
    theme = theme(
      plot.title   = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 7.5, colour = "grey40",
                                  hjust = 0, lineheight = 1.3)
    )
  )

ggsave(file.path(OUTPUT_DIR, "Histogram_FstZoom_CorrectN.pdf"),
       plot = figure_hist_zoom,
       width = 18, height = 24, units = "cm", device = "pdf", bg = "white")

ggsave(file.path(OUTPUT_DIR, "Histogram_FstZoom_CorrectN.png"),
       plot = figure_hist_zoom,
       width = 18, height = 24, units = "cm", dpi = 300, bg = "white")

message("Zoom-in histogram figure saved.")
message("Both histogram figures written to: ", OUTPUT_DIR)