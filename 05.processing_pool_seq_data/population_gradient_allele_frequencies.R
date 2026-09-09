library(stringr)
library(ggplot2)
library(viridis)
library(tidyverse)
library(cowplot)

# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
FIGSHARE_ROOT <- "../Figshare"

# ==========================================
# 1. READ & ANNOTATE TARGET SNPS (TABLE 1)
# ==========================================

table1_snps <- read.table(file.path(FIGSHARE_ROOT, "5.processing_pool_seq_data/results/results_bs_as_with_fst.txt"),
                          header=TRUE, sep="\t") %>%
  mutate(CHROM = str_split_fixed(intro_reg, "_", 2)[,1],
         POS   = Position_Spring) %>%
  select(CHROM, POS, intro_reg)

# Create the explicit "Chr:Coordinate" label you requested
table1_snps <- table1_snps %>%
  mutate(chr_pos = paste0(CHROM, ":", POS))

# correct Snps positions
table1_snps<-table1_snps[c(c(24,25),c(1:23)),]

# Define Y-axis factor level order matching your dataset rows
region_order <- unique(table1_snps$chr_pos)
n_regions    <- length(region_order)

# ==========================================
# 2. LOAD FILES & ESTABLISH METADATA
# ==========================================

load(file.path(FIGSHARE_ROOT, "5.processing_pool_seq_data/inputs/60.Neff.AF.2024-08-14.Rdata"))
names <- names(pops_freq_df)

REPO <- "."

# The live 5.processing_pool_seq_data/plotting_files/pool_order is STILL the old
# version -- Orkney in Ireland_Britain, N_NorthSea in NEAtlantic, no
# BrackishSemilandlocked group. Read the corrected one from the repo.
pool_order_raw <- read.table(
  file.path(REPO, "00.revisions_relabeling/pool_order_20260813.txt"),
  header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
colnames(pool_order_raw) <- c("populations", "raw_group")

# Same lookup as Figures 3 and 4, so every figure says the same thing.
# No hardcoded per-population overrides: Landvik and Ringkobing are now
# BrackishSemilandlocked in the file itself, which is a habitat description
# rather than a geographic claim (Ringkobing is Danish).
pool_group_display <- c(
  "PacificOceanPacificherring"   = "Pacific herring",
  "ArcticPacificherring"         = "Arctic Pacific herring",
  "BalsfjordHybridPopulation"    = "Balsfjord",
  "BalticSummer"                 = "Baltic Summer",
  "BalticSpring"                 = "Baltic Spring",
  "BalticAutum"                  = "Baltic Autumn",      # raw code really is "Autum"
  "BrackishSemilandlocked"       = "Semi-landlocked brackish",
  "NEAtlanticTransition" = "NE Atlantic transition",
  "NEAtlanticFjords"             = "NE Atlantic fjords",
  "NEAtlantic"                   = "NE Atlantic",
  "NorthSea"                     = "North Sea",
  "Ireland_Britain"              = "Britain and Ireland",
  "NWAtlantic"                   = "NW Atlantic"
)
stopifnot(all(pool_order_raw$raw_group %in% names(pool_group_display)))

pool_order <- pool_order_raw %>%
  mutate(clean_group = factor(unname(pool_group_display[raw_group]),
                              levels = unname(pool_group_display[unique(raw_group)])))
# ==========================================
# 3. MAIN HEATMAP (SPLIT BY GROUP)
# ==========================================
snp_freqs_all_pops <- pops_freq_df %>%
  inner_join(table1_snps %>% select(CHROM, POS, intro_reg),
             by = c("CHROM", "POS")) %>%
  pivot_longer(cols = -c(CHROM, POS, intro_reg),
               names_to = "populations", values_to = "frequency") %>%
  select(intro_reg, populations, frequency)

# sanity check: every Table 1 SNP should survive the join, once per population
message(nrow(table1_snps), " target SNPs x ",
        ncol(pops_freq_df) - 2, " populations = ",
        nrow(table1_snps) * (ncol(pops_freq_df) - 2), " expected rows; got ",
        nrow(snp_freqs_all_pops))
stopifnot(dplyr::n_distinct(snp_freqs_all_pops$intro_reg) == nrow(table1_snps))


heatmap_data <- snp_freqs_all_pops %>%
  left_join(table1_snps, by = "intro_reg") %>%
  left_join(pool_order, by = "populations") %>%
  filter(!is.na(chr_pos)) %>%
  mutate(chr_pos     = factor(chr_pos, levels = rev(region_order)),
         populations = factor(populations, levels = pool_order$populations))

# polarize frequency:
heatmap_data_polarized <- heatmap_data %>%
  # 1. Group by each individual SNP/locus
  group_by(chr_pos) %>%
  mutate(
    # 2. Get the baseline frequencies for your two reference groups at this locus
    p_pacific     = mean(frequency[clean_group == "Arctic Pacific herring"], na.rm = TRUE),
    p_nw_atlantic = mean(frequency[clean_group == "NW Atlantic"], na.rm = TRUE),
    # 3. If Baltic Spring approaches 0 and NW Atlantic approaches 1, invert the frequency scale
    frequency = case_when(
      p_pacific < 0.8 & p_nw_atlantic > 0.5 ~ 1 - frequency,
      TRUE ~ frequency # Keep it exactly as it is otherwise
    )
  ) %>%
  # 4. Clean up the temporary reference columns and ungroup
  select(-p_pacific, -p_nw_atlantic) %>%
  ungroup()

main_heatmap <- ggplot(heatmap_data_polarized, aes(x = populations, y = chr_pos, fill = frequency)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis(option = "inferno", direction = -1, limits = c(0, 1),
                     na.value = "grey90", name = "Allele\nfrequency") +
  # This magic line splits the columns by group with a small physical space
  facet_grid(. ~ clean_group, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(
    axis.text.x       = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
    axis.text.y       = element_text(size = 7, colour = "black"),
    axis.title        = element_blank(),
    strip.background  = element_blank(), # Hide default facet boxes completely
    strip.text        = element_blank(), # Hide facet titles since legend is on top
    panel.spacing.x   = unit(2, "mm"),   # Adjust width of white gap between groups
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 7),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.2, "cm")
  )

# ==========================================
# 4. TOP TRACKING LEGEND BAR
# ==========================================

top_legend_bar <- ggplot(pool_order, aes(x = populations, y = 1, fill = clean_group)) +
  geom_tile() +
  scale_fill_viridis_d(option = "turbo", guide = "none") + 
  facet_grid(. ~ clean_group, scales = "free_x", space = "free_x") +
  theme_minimal() + 
  theme(
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    strip.background = element_blank(),
    
    # CHANGED: 45-degree angle with custom justifications so the labels 
    # sit cleanly directly above their colored tracking blocks
    strip.text       = element_text(size = 7, angle = 90, hjust = 0,  colour = "black"),
    
    panel.spacing.x  = unit(2, "mm")
    
    # CHANGED: Generous top margin (35mm) to safely containerize angled text strings
    #plot.margin      = margin(t = 35, r = 10, b = 2, l = 5, unit = "mm") 
  )

# ==========================================
# 5. ALIGN & COMBINE WITH COWPLOT
# ==========================================

# Align both plots precisely by their vertical coordinate grids
aligned_plots <- align_plots(top_legend_bar, main_heatmap, align = 'v', axis = 'lr')

# CHANGED: Increase rel_heights so the top section doesn't get squished
final_composite_plot <- plot_grid(
  aligned_plots[[1]], 
  aligned_plots[[2]], 
  ncol = 1, 
  rel_heights = c(1.5,4.5) # Gives 18% vertical space to the labels/bars, 82% to the tiles
)

# Save the final masterpiece
ggsave(plot = final_composite_plot,
       filename = "supplementary_allele_freq_heatmap_segmented.pdf"  # saved to your R working directory; set to wherever you want to save it,
       height = 180, width = 240, units = "mm")





gradient_distribution_plot <- ggplot(heatmap_data_polarized, aes(x = clean_group, y = frequency, fill = clean_group)) +
  # 1. The Violin plot shows the shape/density of the distribution
  geom_violin(scale = "width", alpha = 0.7, color = "white", linewidth = 0.3) +
  
  # 2. An internal narrow boxplot shows the median and quartiles clearly
  geom_boxplot(width = 0.15, fill = "white", color = "black", 
               outlier.size = 0.5, outlier.alpha = 0.5, alpha = 0.9) +
  
  # Use the same vibrant turbo palette to match your heatmap tracking bar
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  
  # Labels and limits
  labs(
    x = "Population Groups",
    y = "Polarized Allele Frequency\n(Arctic Pacific-associated allele)",
    title = "Introgression Gradient Across All Target SNPs"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  # Clean, publication-ready theme
  theme_classic() +
  theme(
    plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x     = element_text(size = 7, angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y     = element_text(size = 7, colour = "black"),
    axis.title.x    = element_text(size = 9, face = "bold", margin = margin(t = 10)),
    axis.title.y    = element_text(size = 9, face = "bold", margin = margin(r = 10)),
    plot.margin     = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm")
  )


#### Combine into a single figure ####

# ==========================================
# 2. MERGE INTO A SINGLE MULTI-PANEL FIGURE
# ==========================================

final_composite_figure <- plot_grid(
  final_composite_plot,
  gradient_distribution_plot,
  ncol = 1,
  labels = c("a", "b"),          # Adds the "A" and "B" panel letters
  label_size = 12,               # Font size for the letters
  label_fontface = "bold",
  rel_heights = c(0.65, 0.35)    # Allocates 65% of vertical space to heatmap, 35% to violins
)

# ==========================================
# 3. SAVE WITH YOUR EXACT DIMENSIONS
# ==========================================

ggsave(
  plot = final_composite_figure,
  filename = "Figure_Introgression_Heatmap_and_Gradient.pdf"  # saved to your R working directory; set to wherever you want to save it,
  width = 180,                   # Locked to your requested 180mm width
  height = 250,                  # Expanded height to give both panels breathing room
  units = "mm",
  bg = "white"                   # Ensures background doesn't default to transparent
)

ggsave(
  plot = final_composite_figure,
  filename = "Figure_Introgression_Heatmap_and_Gradient.png"  # saved to your R working directory; set to wherever you want to save it,
  width = 180,                   # Locked to your requested 180mm width
  height = 250,                  # Expanded height to give both panels breathing room
  units = "mm",
  bg = "white"                   # Ensures background doesn't default to transparent
)
