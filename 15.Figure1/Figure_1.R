# Plot figure1 · R
# =============================================================================
# Figure 1 -- unified script
# 2026-08-17
#
# Panels:
#   a) Sampling map + species distributions        (left column)
#   b) Admixture, K = 6, Top50angela_woInv         (right column, top)
#   c) D-statistic heatmap, P1 x P2, faceted by P3 (right column, middle)
#   d) f4-ratio heatmap,   P1 x P2, faceted by P3  (right column, bottom)
#
# Replaces the previous hand-assembly of separate panels in Affinity.
# =============================================================================

library(sf)
library(rnaturalearth)
library(tidyverse)
library(cowplot)
library(ggsci)

# ---- 0. Paths ---------------------------------------------------------------
# Everything is referenced from these four roots so the script can be run from
# anywhere -- no setwd() juggling between the map folder and Figshare.
# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
REPO     <- "."
FIGSHARE <- "../Figshare"

SAMPLE_INFO <- file.path(REPO, "00.revisions_relabeling/sampleinfo125_consolidated.txt")

# comment.char = "" is essential: the `color` column holds hex codes like
# "#10345C", and read.table()'s default comment.char = "#" would silently
# truncate every row at that field.
indTable <- read.table(SAMPLE_INFO, header = TRUE, sep = "\t", comment.char = "")


# ---- 0b. Shared label lookup ------------------------------------------------
# One source of truth for display text, used by every panel. Keeps the map,
# the admixture strip and the heatmap axes saying the same thing.
# Covers BOTH spellings on purpose: `bars_plot` carries the Figshare-side
# labels (Norway_Coastal, Baltic_Sea_Spring, ...) unless 00.population_labels.R
# Step 7 has since overwritten them with correct_label (Norway_stationary,
# Baltic_Spring, ...). Mapping both means the figure reads the same either way,
# and matches the wording used in the manuscript text.
label_overrides <- c(
  "Canada_Spring_Summer"     = "Canada Spring\nand Summer",
  "Britain_Ireland_NorthSea" = "Britain, Ireland\nand North Sea",
  "Norway_Coastal"           = "Norway stationary",
  "Norway_Costal"            = "Norway stationary",
  "Norway_stationary"        = "Norway stationary",
  "Norway_Pelagic"           = "Norway migratory",
  "Norway_migratory"         = "Norway migratory",
  "Baltic_Sea_Spring"        = "Baltic Spring",
  "Baltic_Sea_Autumn"        = "Baltic Autumn",
  "Ireland_Britain"          = "Britain and Ireland",
  "Barents_Sea"              = "Barents Sea",
  "White_Sea"                = "White Sea"
)
clean_label <- function(x) {
  ifelse(x %in% names(label_overrides), label_overrides[x], gsub("_", " ", x))
}

# Raw values in the `species` column are "C.harengus" / "C.pallasii"
species_overrides <- c("C.harengus" = "Atlantic herring",
                       "C.pallasii" = "Pacific herring")
clean_species <- function(x) {
  ifelse(x %in% names(species_overrides), species_overrides[x], gsub("_", " ", x))
}

dsuite_label_lookup <- c(
  "Baltic_Autumn"  = "Baltic Autumn",
  "Baltic_Spring"  = "Baltic Spring",
  "NorthSea"       = "North Sea",
  "Norway_Costal"  = "Norway stationary",
  "Norway_Pelagic" = "Norway migratory",
  "Canada_Autumn"  = "Canada Autumn",
  "Canada_Spring"  = "Canada Spring",
  "Canada_Summer"  = "Canada Summer",
  "UK"             = "Britain and Ireland",
  "BarentSea"      = "Barents Sea",
  "WhiteSea"       = "White Sea",
  "Balsfjord"      = "Balsfjord",
  "Japan"          = "Japan",
  "Vancouver"      = "Vancouver"
)


# =============================================================================
# PANEL A -- sampling map
# =============================================================================
shp_1 <- st_read(file.path(FIGSHARE, "15.Figure_1/redlist_species_data_840fe3a6-564c-4062-a02e-bd3882b89829/data_0.shp"))
shp_2 <- st_read(file.path(FIGSHARE, "15.Figure_1/redlist_species_data_6eb4b35e-3281-4f47-958f-af0d5ba45835/data_0.shp"))
land  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

points_sf <- st_as_sf(indTable, coords = c("longitude", "latitude"), crs = 4326)

# Lambert Azimuthal Equal-Area near the pole: the samples span Pacific AND
# Atlantic, so a plain lat/lon map would split the Pacific across the antimeridian.
polar_crs <- "+proj=laea +lat_0=75 +lon_0=20 +datum=WGS84 +units=m +no_defs"

# Norwegian Spring-Spawning herring distribution (Institute of Marine Research).
# "Validity area" is the WFS metadata bounding box, not a real range -- excluded.
# st_make_valid() is required: the raw union throws
#   "Loop 46 is not valid: Edge 7 is degenerate (duplicate vertex)"
# from duplicated vertices in the WFS tessellation.
nvg       <- st_read(file.path(FIGSHARE, "15.Figure_1/utbredelse_NVG_Sild/utbredelse_fiskPolygon.shp"))
nvg_range <- nvg %>% filter(map_type_e != "Validity area") %>% st_union() %>% st_sf()
atlantic_herring_sh <- st_union(st_make_valid(nvg_range), shp_1)

# Crop from the actual sample points rather than hardcoded limits, so nothing
# can be clipped by accident (an earlier hand-picked window cut off Vancouver).
points_proj <- st_transform(points_sf, polar_crs)
bb          <- st_bbox(points_proj)
margin_x    <- 0.15 * (bb["xmax"] - bb["xmin"])
margin_y    <- 0.15 * (bb["ymax"] - bb["ymin"])
xlim_range  <- c(bb["xmin"] - margin_x, bb["xmax"] + margin_x)
ylim_range  <- c(bb["ymin"] - margin_y, bb["ymax"] + margin_y)

panel_map <- ggplot() +
  geom_sf(data = land, fill = "gray80", colour = "gray60", linewidth = 0.2) +
  geom_sf(data = shp_2, colour = "#5189C7", fill = "#8AA7C7", alpha = 0.5) +
  geom_sf(data = atlantic_herring_sh, colour = "#FF7F00", fill = "#FDBF6F", alpha = 0.5) +
  geom_sf(data = points_sf, aes(shape = species, fill = species),
          colour = "gray20", size = 1.2, stroke = 0.3) +
  scale_fill_manual(values = c("C.harengus" = "#B25800",   # darker #FF7F00
                               "C.pallasii" = "#385F8B"),  # darker #5189C7
                    labels = c("Atlantic herring", "Pacific herring")) +
  scale_shape_manual(values = c("C.harengus" = 21, "C.pallasii" = 24),
                     labels = c("Atlantic herring", "Pacific herring")) +
  coord_sf(crs = polar_crs, xlim = xlim_range, ylim = ylim_range) +
  theme_minimal(base_size = 7) +
  theme(legend.position = "none",
        axis.title      = element_blank(),
        axis.text       = element_text(size = 5),
        panel.grid      = element_line(colour = "gray90"))



# =============================================================================
# PANEL B -- admixture, K = 6
# =============================================================================
K       <- 6
SEED    <- 10   # the ".10." in the file name is the sNMF seed, not a K value
BASE    <- "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv"
Q_FILE  <- file.path(FIGSHARE, "2.population_structure/admixture/results",
                     paste0(BASE, ".geno.", K, ".", SEED, ".Q"))

tbl <- read.table(Q_FILE, header = FALSE)

# join on `id`; order by order_in_admixture_bars (NOT order_in_admixture --
# those are two different numbering schemes from two different source tables,
# and only the "_bars" one matches bars_plot/color)
merged  <- merge(tbl, indTable, by.x = "V1", by.y = "id")
ordered <- merged[order(merged$order_in_admixture_bars), ]
ordered$V1        <- factor(ordered$V1, levels = ordered$V1)
ordered$bars_plot <- factor(ordered$bars_plot, levels = unique(ordered$bars_plot))

# Population colours: taken from indTable$color and matched BY NAME, not by
# list position. Same 12 hex codes as the hardcoded population_colors vector,
# but positional matching breaks now that Ireland_Britain/North_Sea and
# Canada_Spring/Canada_Summer deliberately share a colour.
population_colors <- ordered %>%
  distinct(bars_plot, color) %>%
  tibble::deframe()

# --- Ancestry-component colours -----------------------------------------------
# sNMF's component order (V2...V7) is arbitrary and changes between runs, so
# hardcoding a colour to a position would silently break on any re-run. Instead
# each colour is anchored to the POPULATION that component dominates, and the
# component is looked up from the data.
comp_cols <- paste0("V", 2:(K + 1))

# mean ancestry of every component within every population group
comp_by_pop <- ordered %>%
  group_by(bars_plot) %>%
  summarise(across(all_of(comp_cols), mean), .groups = "drop")
print(as.data.frame(comp_by_pop))   # rows = populations, cols = components

# anchor population -> colour (edit the names if your bars_plot labels differ;
# the print above lists exactly what's available)
anchors <- c(
  "Vancouver"                = "#10345C",   # Pacific herring, Vancouver
  "White_Sea"                = "#6DB2FF",   # White Sea
  "Canada_Autumn"            = "#D95F02",   # Canada Autumn
  "Canada_Spring_Summer"     = "#1B9E77",   # Canada Spring+Summer / Norway
  "Britain_Ireland_NorthSea" = "#72FEFF",   # Britain + Ireland
  "Baltic_Sea_Spring"        = "#7570B3"    # Baltic Spring
)

missing <- setdiff(names(anchors), comp_by_pop$bars_plot)
if (length(missing))
  stop("anchor population(s) not found in bars_plot: ", paste(missing, collapse = ", "),
       "\navailable: ", paste(comp_by_pop$bars_plot, collapse = ", "))

# for each anchor, which component peaks there?
top_comp <- sapply(names(anchors), function(p) {
  row <- comp_by_pop %>% filter(bars_plot == p) %>% select(all_of(comp_cols))
  comp_cols[which.max(as.numeric(row))]
})

if (anyDuplicated(top_comp)) {
  print(data.frame(anchor = names(anchors), component = top_comp))
  stop("two anchors resolved to the same component -- see table above, ",
       "pick a different anchor population for one of them")
}

admix_colors <- setNames(unname(anchors[names(top_comp)]), top_comp)[comp_cols]
print(admix_colors)

pop_strip <- ggplot(ordered, aes(x = V1, y = 1, fill = bars_plot)) +
  geom_tile() +
  scale_fill_manual(values = population_colors) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none")

admix_bars <- ordered %>%
  pivot_longer(cols = all_of(comp_cols),
               names_to = "component", values_to = "admix_proportion") %>%
  ggplot(aes(x = V1, y = admix_proportion, fill = component)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab(paste0("Ancestry (K", K, ")")) +
  scale_fill_manual(values = admix_colors) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 7) +
  theme(axis.title.y = element_text(size = 6),
        axis.text.y  = element_text(size = 5),
        axis.text.x  = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank())

# --- Headers: population names (45 deg) and species ---------------------------
# CRITICAL: the x-domain of these continuous scales must match the discrete
# panels below EXACTLY, or cowplot's alignment silently drifts. pop_strip and
# admix_bars use scale_x_discrete(expand = c(0,0)), so positions 1..n map to a
# panel running 0.5 to n+0.5 -- hence the limits below. (With ggplot's DEFAULT
# discrete expansion it would be 0.4 to n+0.6 instead.)
n_ind  <- nrow(ordered)
x_lims <- c(0.5, n_ind + 0.5)

label_positions <- ordered %>%
  mutate(x_num       = as.numeric(V1),
         label_clean = vapply(as.character(bars_plot), clean_label, character(1))) %>%
  group_by(bars_plot) %>%
  summarize(x_mid = mean(x_num), label_clean = dplyr::first(label_clean), .groups = "drop")

top_label_grob <- ggplot(label_positions, aes(x = x_mid, y = 0, label = label_clean)) +
  geom_text(angle = 45, hjust = 0, vjust = 0, size = 1.8) +
  scale_x_continuous(limits = x_lims, expand = c(0, 0)) +
  ylim(0, 1) +
  coord_cartesian(clip = "off") +   # let the diagonal labels overflow the panel
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 10, b = 1, l = 1, unit = "mm"))

# Species header, driven by the real `species` column rather than a hardcoded
# "first 5 blocks are Pacific", so it stays correct if the ordering changes.
species_positions <- ordered %>%
  mutate(x_num         = as.numeric(V1),
         species_clean = vapply(as.character(species), clean_species, character(1))) %>%
  group_by(species) %>%
  summarize(x_mid = mean(range(x_num)), species_clean = dplyr::first(species_clean),
            .groups = "drop")

species_label_grob <- ggplot(species_positions, aes(x = x_mid, y = 0.5, label = species_clean)) +
  geom_text(size = 2.2) +
  scale_x_continuous(limits = x_lims, expand = c(0, 0)) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 10, b = 1, l = 1, unit = "mm"))

panel_admix <- cowplot::plot_grid(
  species_label_grob, top_label_grob, pop_strip, admix_bars,
  ncol        = 1,
  align       = "v",
  axis        = "lr",
  rel_heights = c(0.5, 0.75, 0.125, 2)
)


# =============================================================================
# PANELS C & D -- D-statistic and f4-ratio heatmaps
# =============================================================================
atlantic_pops <- c("Baltic_Autumn", "Baltic_Spring", "NorthSea", "Norway_Costal",
                   "Norway_Pelagic", "Canada_Autumn", "Canada_Spring",
                   "Canada_Summer", "UK")
pacific_pops  <- c("Vancouver", "Japan", "WhiteSea", "BarentSea", "Balsfjord")

D_BBAA <- read.table(
  file.path(FIGSHARE, "3.dsuite/results_dtrios",
            "all_vs_all_clusters_v03_maf5_JK1000_run_2209_2023-11-24_BBAA.txt"),
  header = TRUE, as.is = TRUE)

# Restrict to the biologically relevant tests FIRST (P3 Pacific, both P1 and P2
# Atlantic), then BH-correct once within that set -- 180 tests = 36 pairs x 5 P3.
pops_F <- D_BBAA %>%
  filter(P3 %in% pacific_pops) %>%
  filter(P1 %in% atlantic_pops & P2 %in% atlantic_pops)
pops_F$p.adjust <- p.adjust(pops_F$p.value, method = "BH")

plot_order <- c("Baltic_Autumn", "Baltic_Spring", "NorthSea", "Norway_Costal",
                "Norway_Pelagic", "Canada_Autumn", "Canada_Spring",
                "Canada_Summer", "UK")
p3_order   <- c("WhiteSea", "BarentSea", "Balsfjord", "Japan", "Vancouver")
pops_F$P1 <- factor(pops_F$P1, levels = plot_order)
pops_F$P2 <- factor(pops_F$P2, levels = plot_order)
pops_F$P3 <- factor(pops_F$P3, levels = p3_order)

heatmap_panel <- function(fill_var, legend_name, show_x) {
  ggplot(pops_F, aes(x = P2, y = P1, fill = .data[[fill_var]])) +
    geom_tile(color = "white") +
    geom_point(data = subset(pops_F, p.adjust < 0.05),
               aes(x = P2, y = P1), colour = "white", shape = 8, size = 0.6,
               inherit.aes = FALSE) +
    scale_fill_viridis_c(name = legend_name) +
    scale_x_discrete(labels = dsuite_label_lookup) +
    scale_y_discrete(limits = rev(plot_order), labels = dsuite_label_lookup) +
    facet_wrap(~ P3, nrow = 1, labeller = as_labeller(dsuite_label_lookup)) +
    theme_bw(base_size = 7) +
    theme(axis.text.x  = if (show_x)
      element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
      else element_blank(),
      axis.text.y     = element_text(size = 5),
      strip.text      = element_text(size = 6),
      legend.title    = element_text(size = 6),
      legend.text     = element_text(size = 5),
      legend.key.width  = unit(2, "mm"),
      legend.key.height = unit(3, "mm"),
      panel.spacing   = unit(0.5, "mm")) +
    labs(x = NULL, y = NULL)
}

panel_dstat <- heatmap_panel("Dstatistic", "D-stats", show_x = FALSE)
panel_f4    <- heatmap_panel("f4.ratio",   "f4-ratio",    show_x = TRUE)

cowplot::plot_grid(panel_dstat,panel_f4, align="h", nrow=2)

# =============================================================================
# ASSEMBLY -- two columns: map on the left, the other three stacked on the right
# =============================================================================
heatmaps <- cowplot::plot_grid(
  panel_dstat, panel_f4,
  ncol        = 1,
  align       = "v",
  axis        = "lr",
  labels      = c("c", "d"),
  label_size  = 9,
  rel_heights = c(0.6, 1)
)

right_column <- cowplot::plot_grid(
  panel_admix, heatmaps,
  ncol        = 1,
  labels      = c("b", ""),
  label_size  = 9,
  rel_heights = c(1, 1.6)
)

figure1 <- cowplot::plot_grid(
  panel_map, right_column,
  ncol       = 2,
  labels     = c("a", ""),
  label_size = 9,
  rel_widths = c(0.9, 1.6)
)

ggsave(figure1,
       filename = file.path(FIGSHARE, "15.Figure_1/Figure1_20260817.pdf"),
       units = "mm", width = 180, height = 100)

ggsave(figure1,
       filename = "Figure1_20260817.pdf",  # saved to your R working directory; set to wherever you want to save it
       units = "mm", width = 180, height = 100)



figure1


