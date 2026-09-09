# 00.population_labels.R
#
# Single source of truth for population labels across the Baltic herring
# introgression figures. This script:
#   1. Consolidates the two sampleinfo125 files (GitHub curated labels +
#      Figshare admixture-plotting metadata) into one table.
#   2. Flags rows where the label actually changed vs. the old figure, so you
#      can eyeball the diff before trusting it.
#   3. Defines the display-label lookup for the pool-seq pipeline (pool_order).
#
# Run this with the working directory set to the Figshare root, same
# convention as the other plotting scripts (see e.g. plot_dsuite.R's setwd()).

library(dplyr)

# ---- 0. Paths -----------------------------------------------------------
# Assumes: (1) your R working directory is the root of this repo (e.g. open
# Baltic_herr_introgression as an RStudio project before running), and
# (2) the Figshare data folder is downloaded as a sibling of this repo, e.g.
# .../some_folder/Baltic_herr_introgression and .../some_folder/Figshare.
# Adjust both paths below if your local layout differs.
GITHUB_REPO   <- "."
FIGSHARE_ROOT <- "../Figshare/"

# ---- 1. Read both source files -------------------------------------------
github_labels <- read.table(
  file.path(GITHUB_REPO, "00.revisions_relabeling/sampleinfo125_20260813.txt"),
  header = TRUE, sep = "\t", as.is = TRUE, comment.char = ""
)

figshare_admix <- read.table(
  file.path(FIGSHARE_ROOT, "2.population_structure/plotting_files/sampleinfo125.txt.ordered.mafalda"),
  header = TRUE, sep = "\t", as.is = TRUE, comment.char = ""
)

# ---- 2. Sanity check before merging ---------------------------------------
# Both files should describe exactly the same 125 individuals. If either
# setdiff below is non-empty, STOP and resolve it before going further --
# don't silently proceed with a partial join.
missing_from_figshare <- setdiff(github_labels$name_in_vcf, figshare_admix$name_in_vcf)
missing_from_github    <- setdiff(figshare_admix$name_in_vcf, github_labels$name_in_vcf)

if (length(missing_from_figshare) > 0) {
  stop("These individuals are in the GitHub table but not in Figshare's ordered.mafalda: ",
       paste(missing_from_figshare, collapse = ", "))
}
if (length(missing_from_github) > 0) {
  stop("These individuals are in Figshare's ordered.mafalda but not in the GitHub table: ",
       paste(missing_from_github, collapse = ", "))
}

# ---- 3. Join on name_in_vcf + order_in_vcf (double key) -------------------
# order_in_vcf is included as a belt-and-braces check: a handful of
# individuals (the Vancouver samples) have name_in_vcf != id in BOTH files,
# consistently -- this looks like a real historical sample re-ID, not a
# copy-paste error (confirmed: both files agree on the pairing when matched
# by order_in_vcf). Joining on both columns protects against silently
# mismatching those rows if one file's ordering ever shifts.
combined <- github_labels %>%
  select(name_in_vcf, id, order_in_vcf, sampling_location, region, longitude, latitude,
         collection_season, species, order_in_admixture, correct_label, correct_label_2) %>%
  left_join(
    figshare_admix %>%
      select(name_in_vcf, order_in_vcf, bars_plot, color,
             order_in_admixture_bars = order_in_admixture),
    by = c("name_in_vcf", "order_in_vcf")
  )

# ---- 4. Confirm the join actually worked -----------------------------------
missing_join <- combined %>% filter(is.na(bars_plot))
if (nrow(missing_join) > 0) {
  warning(nrow(missing_join), " individuals failed to join -- check these by hand:")
  print(missing_join %>% select(name_in_vcf, id, order_in_vcf))
} else {
  message("All ", nrow(combined), " individuals matched cleanly.")
}

# ---- 5. Changelog: which individuals get a NEW label vs. the old figure? --
# Review this before trusting the consolidated file -- it's your check that
# nothing moved groups by accident.
relabel_changes <- combined %>%
  filter(correct_label != bars_plot) %>%
  select(name_in_vcf, sampling_location, old_label = bars_plot,
         new_label = correct_label, new_label_coarse = correct_label_2)

message(nrow(relabel_changes), " of ", nrow(combined), " individuals have a label that differs from the old figure.")
print(relabel_changes)

# Modify bars_plot such that Norwegian populations have the correct name:

combined <- combined %>%
  mutate(
    bars_plot = if_else(
      (grepl("Norway", region, ignore.case = TRUE) |
         grepl("Norway", sampling_location, ignore.case = TRUE)) &
        bars_plot != correct_label,
      correct_label,
      bars_plot
    )
  )

# ---- 6. Write the consolidated file ----------------------------------------
write.table(combined,
            file.path(GITHUB_REPO, "00.revisions_relabeling/sampleinfo125_consolidated.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ============================================================================
# Pool-seq (pool_order) display-label lookup
# ============================================================================
# pool_order gives each pool a raw_group code (2nd column). This is the
# single place that maps those codes to the text that actually appears on a
# figure -- update THIS, not the case_when() blocks scattered across scripts.

pool_group_display <- c(
  "PacificOceanPacificherring"    = "Pacific herring",
  "ArcticPacificherring"          = "Arctic Pacific herring",
  "BalsfjordHybridPopulation"     = "Balsfjord",
  "BalticSummer"                  = "Baltic Summer",
  "BalticSpring"                  = "Baltic Spring",
  "BalticAutum"                   = "Baltic Autumn",
  "BrackishSemilandlocked"        = "Semi-landlocked brackish",  # <- fix the pool_order typo to match this
  "BalticAtlanticTransitionZone"  = "NE Atlantic Transition",
  "NEAtlanticFjords"              = "NE Atlantic Fjords",
  "NEAtlantic"                    = "NE Atlantic",
  "NorthSea"                      = "North Sea",
  "Ireland_Britain"               = "Ireland/Britain",
  "NWAtlantic"                    = "NW Atlantic"
)



# Example usage in a plotting script, replacing the per-script case_when():
#   pool_order <- read.table("5.processing_pool_seq_data/plotting_files/pool_order",
#                             header = FALSE, stringsAsFactors = FALSE) %>%
#     rename(populations = V1, raw_group = V2) %>%
#     mutate(clean_group = pool_group_display[raw_group])
#
#   # sanity check -- catches typos immediately instead of silent NAs:
#   stopifnot(!any(is.na(pool_order$clean_group)))