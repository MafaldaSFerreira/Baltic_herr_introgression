# =============================================================================
# Figure 4 (THRB, chr19) and Supplementary Figure 25
# Baltic herring introgression paper
#
# Rewritten 2026-08-26. Same figures as before; the difference is that nothing
# here writes to a global that another section depends on.
#
# WHY THE REWRITE
#   The previous version drove every panel off shared globals (c, s, e, chr,
#   padding, padding_mb, candidate_matrix, tmp_*) and then REDEFINED them
#   further down for the supplementary figures. That made the script
#   order-dependent: re-running any block out of sequence silently produced a
#   different figure -- e.g. re-running the main panels after the supplement
#   gave an FST scan with no FST colouring, over the wrong coordinates, with no
#   error. It also only worked at all because a stale `c` was left in the
#   environment by a previous run.
#
#   Here, a genomic region is a value (make_region()), every panel is a
#   function of (region, data), and each figure is assembled in one place.
#   Blocks can be run in any order, any number of times, in a fresh session.
#
# LAYOUT -- Figure 4
#   gene track   THRB + GMPR, unlettered annotation strip   6.00-6.75 Mb
#   a            dAF coloured by FST                        6.00-6.75 Mb
#   b            allele-frequency heatmap                   6.25-6.50 Mb
#   c            XPEHH                                      6.25-6.50 Mb
#   d            TWISST                                     6.25-6.50 Mb
#   e            THRB H40Q vs RHO scatterplots + legend
#
# LAYOUT -- Supplementary Figure 25
#   a            ENSEMBL 114 / 109 gene models + chi-square over the THRB span
#   b            THRB H40Q allele frequency, pool-seq populations
#   c            THRB H40Q genotypes, 125 individuals + European sprat
# =============================================================================

library(tidyverse)
library(cowplot)
library(viridis)
library(ggrastr)


# =============================================================================
# 1. CONFIGURATION -- the only paths in the script
# =============================================================================
# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
FIGSHARE <- "../Figshare"
REPO     <- "."
OUTDIR   <- "."  # set to wherever you want to save these figures locally
STAMP    <- "20260903"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
fig_path <- function(name) file.path(OUTDIR, paste0(name, "_", STAMP, ".pdf"))


# =============================================================================
# 2. DATA -- loaded once, never reassigned
# =============================================================================

## TWISST ----------------------------------------------------------------
twisst_weights <- read.table(file.path(FIGSHARE, "6.twisst/outputs/output_outgroup_sprat/output.wg.phyml.w100.Mi10.weights.csv.gz"),
                             skip = 3, header = TRUE)
twisst_windows <- read.table(file.path(FIGSHARE, "6.twisst/outputs/output_outgroup_sprat/output.wg.phyml.w100.Mi10.data.tsv"),
                             header = TRUE)

twisst_weights <- twisst_weights / apply(twisst_weights, 1, sum)
keep_rows      <- which(!is.na(apply(twisst_weights, 1, sum)))

twisst_df <- cbind(twisst_windows[keep_rows, ], twisst_weights[keep_rows, ]) %>%
  mutate(topo1_frac = topo1 / (topo1 + topo2 + topo3),
         topo2_frac = topo2 / (topo1 + topo2 + topo3),
         topo3_frac = topo3 / (topo1 + topo2 + topo3))

## Pool-seq allele frequencies -------------------------------------------
load(file.path(FIGSHARE, "5.processing_pool_seq_data/inputs/60.Neff.AF.2024-08-14.Rdata"))  # -> pops_freq_df

## Chi-square + FST -------------------------------------------------------
load(file.path(FIGSHARE, "5.processing_pool_seq_data/chisquare_results/baltic_spring_vs_baltic_autumn_han_pops_chiseq.output.RData"))
load(file.path(FIGSHARE, "5.processing_pool_seq_data/chisquare_results/baltic_spring_vs_atlantic_spring_han_pops_chiseq.output.RData"))
load(file.path(FIGSHARE, "5.processing_pool_seq_data/fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05.RData"))  # -> df_spring_maf0.05

## Mutations, introgression regions, annotation ---------------------------
moderate_mutations <- read.table(file.path(FIGSHARE, "5.processing_pool_seq_data/results/moderate_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt"),
                                 sep = "\t", header = TRUE)

gene_annot <- read.table(file.path(FIGSHARE, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.maxgap20K.modified.txt"),
                         sep = "\t", header = TRUE, row.names = NULL)

## XPEHH ------------------------------------------------------------------
xpehh_raw <- read.table(file.path(FIGSHARE, "8.xpehh/results_xpehh/xpehh_group1+2.out"), skip = 1)
xpehh_raw$CHROM <- str_split_fixed(xpehh_raw$V2, "_", 2)[, 1]
xpehh_raw$POS   <- as.numeric(str_split_fixed(xpehh_raw$V2, "_", 2)[, 2])
xpehh_raw$V2    <- NULL
colnames(xpehh_raw) <- c("Index", "Freq", "iHH_A1", "iHH_B1", "iHH_P1",
                         "XPEHH", "std_XPEHH", "CHROM", "POS")

xpehh_results <- xpehh_raw %>%
  arrange(desc(abs(std_XPEHH))) %>%
  mutate(rank      = row_number(),
         rank.perc = rank / n(),
         rank.log  = -log10(rank.perc))

## ENSEMBL transcript models ----------------------------------------------
# Collapse exons to one row per transcript to get the full transcript span,
# then join back so both the span and the individual exons can be drawn.
read_ensembl <- function(path) {
  ann <- read.table(path, header = TRUE, sep = "\t")
  spans <- ann %>%
    group_by(Transcript.stable.ID) %>%
    summarise(min = min(Exon.region.start..bp.),
              max = max(Exon.region.end..bp.), .groups = "drop") %>%
    mutate(n = row_number())
  left_join(ann, spans, by = "Transcript.stable.ID")
}
ensembl114 <- read_ensembl(file.path(FIGSHARE, "13.THRB_figures/annotations/Ch_v2.0.2v2_ensembl114_chr19-5.88Mb-6.88Mb_annotation_genestart.txt"))
ensembl109 <- read_ensembl(file.path(FIGSHARE, "13.THRB_figures/annotations/Ch_v2.0.2_ensembl109_chr19-5.88Mb-6.88Mb_annotation_genestart.txt"))

## Genotypes (125 herring + European sprat) --------------------------------
genotype_header <- read.table(file.path(FIGSHARE, "1.mapping_variant_calling/genotypes/header_vcf_20231116.txt"))
individual_order <- read.table(file.path(FIGSHARE, "1.mapping_variant_calling/genotypes/individuals_v01_20231116.txt"))

genotypes_raw <- read.table(file.path(FIGSHARE, "1.mapping_variant_calling/genotypes/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.intro_reg.genotypes.2025-08-06.txt"))
colnames(genotypes_raw) <- genotype_header[, 1]

sprat_raw <- read.table(file.path(FIGSHARE, "1.mapping_variant_calling/genotypes/all_mutations_combined.sprat.genotypes.txt"),
                        sep = "\t", header = FALSE)
sprat_raw$V4 <- NULL
colnames(sprat_raw) <- c("CHROM", "POS", "SPRAT")

# Recode only the genotype columns, leaving CHROM/POS alone. The old version
# ran the substitution over the whole data frame, which happened to be safe but
# would corrupt CHROM/POS if either ever contained one of these strings.
recode_gt <- function(df) {
  gt_cols <- setdiff(names(df), c("CHROM", "POS"))
  df[gt_cols] <- lapply(df[gt_cols], function(x) {
    x <- as.character(x)
    x[x == "0/0"]            <- "0"
    x[x %in% c("0/1", "1/0")] <- "1"
    x[x == "1/1"]            <- "2"
    x[x == "./."]            <- NA
    x
  })
  df
}

genotypes_all <- left_join(recode_gt(genotypes_raw), recode_gt(sprat_raw),
                           by = c("CHROM", "POS"))

## Population metadata -----------------------------------------------------
# comment.char = "" is essential: the `color` column holds hex codes like
# "#10345C" and read.table()'s default comment.char = "#" truncates every row
# at that field, silently.
sample_info <- read.table(file.path(REPO, "00.revisions_relabeling/sampleinfo125_consolidated.txt"),
                          header = TRUE, sep = "\t", comment.char = "")

# The live 5.processing_pool_seq_data/plotting_files/pool_order is STILL the old
# version -- Orkney in Ireland_Britain, N_NorthSea in NEAtlantic, no
# BrackishSemilandlocked group. Read the corrected one from the repo instead.
pool_order <- read.table(file.path(REPO, "00.revisions_relabeling/pool_order_20260813.txt"),
                         header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         comment.char = "")
colnames(pool_order) <- c("population", "raw_group")

## THRB/RHO allele frequencies (Goodall et al. 2025) ------------------------
rho_thrb <- read.table(file.path(FIGSHARE, "13.THRB_figures/goodall_et_al_2025_data/jake_goodall-THRB/merged_refalt_freq_long_with_metadata.tsv"),
                       header = TRUE, sep = "\t")


# =============================================================================
# 3. DISPLAY LABELS -- one source of truth for every figure
# =============================================================================

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
stopifnot(all(pool_order$raw_group %in% names(pool_group_display)))

# Both spellings on purpose: sample_info$bars_plot carries the Figshare-side
# labels unless 00.population_labels.R Step 7 has overwritten them with
# correct_label. Mapping both means the figure reads the same either way.
individual_group_display <- c(
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
clean_individual_label <- function(x) {
  ifelse(x %in% names(individual_group_display),
         individual_group_display[x], gsub("_", " ", x))
}

# Population labels for the genotype heatmap. Deliberately coarser than
# bars_plot on the Pacific side, to match the grouping the pool-seq allele
# frequency heatmaps use: Vancouver + Japan read as one group, Barents + White
# Sea as another, Balsfjord kept separate. Hardcoded on purpose -- these are a
# presentation choice for this figure, not a property of the sample table.
#
# Both spellings are covered because bars_plot carries the Figshare-side labels
# unless 00.population_labels.R Step 7 has overwritten them with correct_label.
GENOTYPE_GROUP_DISPLAY <- c(
  "Vancouver"                = "Pacific Pacific herring",
  "Japan"                    = "Pacific Pacific herring",
  "Barents_Sea"              = "Arctic Pacific herring",
  "White_Sea"                = "Arctic Pacific herring",
  "Pechora_Sea"              = "Arctic Pacific herring",
  "Balsfjord"                = "Balsfjord hybrid\npopulation",
  "Canada_Autumn"            = "Canada Autumn",
  "Canada_Spring_Summer"     = "Canada Spring\nand Summer",
  "Canada_Spring"            = "Canada Spring\nand Summer",
  "Canada_Summer"            = "Canada Spring\nand Summer",
  "Britain_Ireland_NorthSea" = "Britain, Ireland\nand North Sea",
  "Ireland_Britain"          = "Britain, Ireland\nand North Sea",
  "North_Sea"                = "Britain, Ireland\nand North Sea",
  "Norway_Coastal"           = "Norway stationary",
  "Norway_Costal"            = "Norway stationary",
  "Norway_stationary"        = "Norway stationary",
  "Norway_Pelagic"           = "Norway migratory",
  "Norway_migratory"         = "Norway migratory",
  "Baltic_Sea_Autumn"        = "Baltic Autumn",
  "Baltic_Autumn"            = "Baltic Autumn",
  "Baltic_Sea_Spring"        = "Baltic Spring",
  "Baltic_Spring"            = "Baltic Spring",
  "SPRAT"                    = "European sprat"
)

species_display <- c("C.harengus" = "Atlantic herring",
                     "C.pallasii" = "Pacific herring",
                     "SPRAT"      = "European sprat")

# =============================================================================
# 4. HELPERS
# =============================================================================

# A genomic region is a value, not a set of globals. Pass it to a panel builder;
# nothing is left behind in the environment.
make_region <- function(chrom, start, end, padding_mb = 0) {
  list(chrom      = chrom,
       chr_num    = sub("^chr", "", chrom),
       start      = start,
       end        = end,
       padding_mb = padding_mb,
       xr_zoom    = c(start, end) / 1e6,
       xr_wide    = c(start - padding_mb, end + padding_mb) / 1e6)
}

# Introgression block drawn as a grey rectangle. Returned in Mb so every panel
# uses it the same way -- the old script had one version in bp and one in Mb,
# which is exactly the kind of thing that produces a silently misplaced box.
candidate_region_mb <- function(start_bp, end_bp) {
  data.frame(start = start_bp / 1e6, end = end_bp / 1e6)
}
thrb_candidate <- candidate_region_mb(6360001, 6380000)

# Consecutive heatmap rows sharing a group form a "run"; the group name goes on
# the axis at that run's midpoint. Replaces hand-added labels and the side
# colour-strip panels, and cannot drift out of register with the rows.
group_runs <- function(row_levels, row_groups) {
  r      <- rle(as.character(row_groups))
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  mids   <- floor((starts + ends) / 2)
  
  if (anyDuplicated(r$values)) {
    dup <- unique(r$values[duplicated(r$values)])
    warning("Groups not contiguous in the plotting order, so each is labelled ",
            "more than once: ", paste(dup, collapse = ", "), call. = FALSE)
  }
  
  data.frame(group = r$values, start_i = starts, end_i = ends,
             mid_row = row_levels[mids], boundary = ends + 0.5,
             stringsAsFactors = FALSE)
}

# Order pools so each group is one contiguous block (moves Orkney beside
# N_NorthSea), preserving the order in which groups first appear.
pool_order_blocked <- pool_order %>%
  mutate(raw_group = factor(raw_group, levels = unique(raw_group))) %>%
  arrange(raw_group) %>%
  mutate(raw_group = as.character(raw_group))

# levels[1] is drawn at the BOTTOM of a discrete axis, so rev() puts pool_order
# row 1 (Pacific herring) at the TOP.
af_levels  <- rev(pool_order_blocked$population)
af_runs    <- group_runs(af_levels, rev(pool_order_blocked$raw_group))
af_labels  <- unname(pool_group_display[af_runs$group])

# Long ticks act as leader lines: they point at the exact row a label belongs
# to, so text can be nudged by hand afterwards without losing the reference.
theme_group_axis <- function(size = 5) {
  theme(axis.text.y         = element_text(size = size, colour = "black", lineheight = 0.8),
        axis.ticks.y        = element_line(colour = "grey40", linewidth = 0.2),
        axis.ticks.length.y = unit(2.5, "mm"),
        axis.line.y         = element_blank())
}

theme_track <- function() {
  theme_classic() +
    theme(axis.text.y  = element_text(colour = "black", size = 6),
          axis.text.x  = element_text(colour = "black", size = 6),
          axis.title.y = element_text(colour = "black", size = 6),
          axis.title.x = element_blank(),
          legend.position = "none")
}


# =============================================================================
# 5. PANEL BUILDERS -- each is a pure function of (region, data)
# =============================================================================
# coord_cartesian() is used throughout rather than xlim(): it crops the view
# without deleting rows, so nothing feeding a line or a smooth is silently lost.

panel_gene_track <- function(region, genes_wanted, candidate,
                             mark_zoom = NULL, label_size = 2.2) {
  genes <- gene_annot %>%
    filter(seqnames == region$chr_num,
           tolower(external_gene_name) %in% tolower(genes_wanted)) %>%
    mutate(start_mb = as.numeric(start) / 1e6,
           end_mb   = as.numeric(end) / 1e6,
           mid_mb   = (start_mb + end_mb) / 2,
           label    = toupper(external_gene_name))
  
  missing <- setdiff(tolower(genes_wanted), tolower(genes$external_gene_name))
  if (length(missing)) {
    in_window <- gene_annot %>%
      filter(seqnames == region$chr_num,
             as.numeric(end) / 1e6   >= region$xr_wide[1],
             as.numeric(start) / 1e6 <= region$xr_wide[2])
    warning("gene(s) not found: ", paste(missing, collapse = ", "),
            "\navailable in this window: ",
            paste(sort(unique(in_window$external_gene_name)), collapse = ", "),
            call. = FALSE)
  }
  
  p <- ggplot(genes) +
    geom_rect(data = candidate, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
              fill = "azure3", colour = NA, alpha = 0.5) +
    geom_rect(aes(xmin = start_mb, xmax = end_mb, ymin = 0.15, ymax = 0.45),
              colour = "black", fill = "lightblue", linewidth = 0.25) +
    geom_text(aes(x = mid_mb, y = 0.55, label = label),
              size = label_size, hjust = 0.5, vjust = 0, fontface = "italic")
  
  if (!is.null(mark_zoom)) {
    p <- p + geom_vline(xintercept = mark_zoom, linetype = "dashed",
                        colour = "grey45", linewidth = 0.2)
  }
  
  p + coord_cartesian(xlim = region$xr_wide, ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank(),
          text = element_text(size = 6), legend.position = "none")
}


panel_fst_scan <- function(region, candidate, xr, mark_zoom = NULL) {
  dat <- baltic_spring_vs_atlantic_spring_chiseq %>%
    filter(CHROM == region$chrom,
           POS >= region$start - region$padding_mb,
           POS <= region$end   + region$padding_mb) %>%
    left_join(df_spring_maf0.05 %>% filter(FST >= 0) %>% select(Chr, Pos, FST),
              by = c("CHROM" = "Chr", "POS" = "Pos"))
  
  p <- ggplot(dat) +
    geom_rect(data = candidate, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0, ymax = 45),
              fill = "azure3", colour = NA, alpha = 0.5) +
    # SNPs with no FST value drawn first so they sit behind
    geom_point(data = dat %>% filter(is.na(FST)),
               aes(x = POS / 1e6, y = -log10(chisq_p)),
               colour = "darkgreen", size = 0.5) +
    geom_point(data = dat %>% filter(!is.na(FST)),
               aes(x = POS / 1e6, y = -log10(chisq_p), colour = FST), size = 0.5) +
    scale_colour_viridis(option = "magma", direction = -1,
                         breaks = c(0, 0.1, 0.2, 0.25),
                         labels = c("0", "0.1", "0.2", ">0.25"),
                         limits = c(0, 0.25), oob = scales::squish,
                         name = expression(F[ST]))
  
  if (!is.null(mark_zoom)) {
    p <- p + geom_vline(xintercept = mark_zoom, linetype = "dashed",
                        colour = "grey45", linewidth = 0.2)
  }
  
  p + coord_cartesian(xlim = xr) +
    ylab(expression(paste(-log[10], " (P value)"))) +
    theme_track() +
    theme(legend.position   = "right",
          legend.key.size   = unit(0.3, "cm"),
          legend.title      = element_text(size = 6),
          legend.text       = element_text(size = 5))
}


panel_af_heatmap <- function(region, xr) {
  dat <- pops_freq_df %>%
    filter(CHROM == region$chrom, POS >= region$start, POS <= region$end) %>%
    # by name, not by column index: pivot_longer(cols = 3:62) breaks silently
    # if a column is ever added upstream
    pivot_longer(cols = -c(CHROM, POS),
                 names_to = "populations", values_to = "frequencies")
  
  ggplot(dat) +
    rasterise(geom_tile(aes(y = factor(populations, levels = af_levels),
                            x = POS / 1e6, colour = frequencies)), dpi = 300) +
    geom_hline(yintercept = head(af_runs$boundary, -1), colour = "white",
               linewidth = 0.3) +
    # the legend sits OUTSIDE the panel area, and cowplot aligns on the panel,
    # so sibling panels are simply padded -- the x-axes stay aligned
    scale_colour_viridis(direction = -1, name = "Allele\nfrequency",
                         breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    scale_y_discrete(breaks = af_runs$mid_row, labels = af_labels,
                     expand = c(0, 0)) +
    coord_cartesian(xlim = xr) +
    ylab(NULL) +
    theme_classic() +
    theme(axis.title.x    = element_blank(),
          axis.text.x     = element_text(size = 6),
          legend.position = "right",
          legend.title    = element_text(size = 5),
          legend.text     = element_text(size = 4),
          legend.key.width  = unit(2, "mm"),
          legend.key.height = unit(3, "mm")) +
    theme_group_axis()
}


panel_xpehh <- function(region, candidate, xr) {
  dat <- xpehh_results %>%
    filter(CHROM == region$chrom, POS >= region$start, POS <= region$end)
  
  ggplot(dat) +
    geom_rect(data = candidate, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -6, ymax = 6),
              fill = "azure3", colour = NA, alpha = 0.5) +
    geom_point(aes(x = POS / 1e6, y = std_XPEHH,
                   colour = factor(abs(std_XPEHH) >= 2)), size = 0.5) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    coord_cartesian(xlim = xr, ylim = c(-6, 6)) +
    ylab("XPEHH score") +
    theme_track()
}


panel_twisst <- function(region, candidate, xr) {
  twisst_df %>%
    filter(scaffold == region$chrom) %>%
    ggplot() +
    geom_rect(data = candidate, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
              fill = "azure3", colour = NA, alpha = 0.5) +
    geom_line(aes(x = mid / 1e6, y = topo1_frac), colour = "#882E72", linewidth = 0.5) +
    geom_line(aes(x = mid / 1e6, y = topo2_frac), colour = "#E8601C", linewidth = 0.5) +
    geom_line(aes(x = mid / 1e6, y = topo3_frac), colour = "#5289C7", linewidth = 0.5) +
    coord_cartesian(xlim = xr) +
    ylab("topology \nsupport (%)") +
    xlab(paste0("Chr ", region$chr_num, " position (Mb)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          axis.line        = element_line(colour = "black"),
          axis.text.y      = element_text(colour = "black", size = 6),
          axis.text.x      = element_text(colour = "black", size = 6),
          axis.title.y     = element_text(colour = "black", size = 6),
          axis.title.x     = element_text(colour = "black", size = 6),
          legend.position  = "none")
}


# --- Panel e: THRB vs RHO ------------------------------------------------
# Stats are computed here rather than pasted in, so the numbers on the figure
# cannot drift from the data. Regression direction does not affect either
# number: for a simple linear regression R-squared and the slope P value are
# identical whether you fit y ~ x or x ~ y.
fit_label <- function(dat, yvar, xvar) {
  fit <- summary(lm(reformulate(xvar, yvar), data = dat))
  p   <- stats::pf(fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3],
                   lower.tail = FALSE)
  sprintf("adj. R² = %.2f, %s", fit$adj.r.squared,
          if (p < 2.2e-16) "P < 2.2e-16" else sprintf("P = %.3g", p))
}

# Jake's REF/ALT are swapped relative to our reference FASTA, so to plot the
# Baltic allele at both loci we take ALT for THRB and REF for RHO.
rho_thrb_pair <- function(rho_snp) {
  rho_thrb %>%
    filter(SNP %in% c("19_6359267", rho_snp)) %>%
    pivot_wider(names_from = "SNP",
                values_from = c("REF_FREQ", "ALT_FREQ", "CHR", "REF", "ALT")) %>%
    filter(sample_n >= 10)
}

panel_rho_scatter <- function(dat, yvar, ylab_txt, season_colors) {
  ggplot(dat, aes(x = ALT_FREQ_19_6359267, y = .data[[yvar]])) +
    geom_point(aes(colour = SEASON), size = 1) +
    geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.6) +
    scale_colour_manual(values = season_colors) +
    labs(x = "THRB Q40H frequency", y = ylab_txt,
         title = fit_label(dat, yvar, "ALT_FREQ_19_6359267")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text  = element_text(size = 6),
          axis.title = element_text(size = 6),
          plot.title = element_text(size = 5.5, hjust = 0.5, face = "plain"))
}


# =============================================================================
# 6. FIGURE 4
# =============================================================================
region_fig4 <- make_region("chr19", start = 6.25e6, end = 6.5e6, padding_mb = 2.5e5)

fig4_gene_track <- panel_gene_track(region_fig4, c("thrb", "gmpr"), thrb_candidate,
                                    mark_zoom = region_fig4$xr_zoom)
fig4_a <- panel_fst_scan(region_fig4, thrb_candidate, region_fig4$xr_wide,
                         mark_zoom = region_fig4$xr_zoom)
fig4_b <- panel_af_heatmap(region_fig4, region_fig4$xr_zoom)
fig4_c <- panel_xpehh(region_fig4, thrb_candidate, region_fig4$xr_zoom)
fig4_d <- panel_twisst(region_fig4, thrb_candidate, region_fig4$xr_zoom)

season_colors <- c("#F2DF0D", "#00BA38", "#545454")
dat_F261Y <- rho_thrb_pair("4_11217516")
dat_T213I <- rho_thrb_pair("4_11217660")

fig4_e1 <- panel_rho_scatter(dat_F261Y, "REF_FREQ_4_11217516",
                             "RHO F261Y frequency", season_colors)
fig4_e2 <- panel_rho_scatter(dat_T213I, "REF_FREQ_4_11217660",
                             "RHO T213I frequency", season_colors)

fig4_e_legend <- get_legend(
  fig4_e1 + theme(legend.position = "right",
                  legend.title    = element_text(size = 6),
                  legend.text     = element_text(size = 5),
                  legend.key.size = unit(3, "mm")) +
    guides(colour = guide_legend(override.aes = list(size = 1.5)))
)

# Three blocks, because they do not share an x axis: the wide pair, the zoomed
# trio, and the scatterplots. Nesting keeps each alignment local -- aligning
# across blocks would line up panels whose axes mean different things.
#
# PANEL LETTERS: the gene track is an unlettered annotation strip above a.
# To letter it, change the first labels vector to c("a", "b") and shift the rest.
fig4_block_wide <- plot_grid(fig4_gene_track, fig4_a,
                             ncol = 1, align = "v", axis = "lr",
                             labels = c("", "a"), label_size = 9,
                             rel_heights = c(0.7, 2))

fig4_block_zoom <- plot_grid(fig4_b, fig4_c, fig4_d,
                             ncol = 1, align = "v", axis = "lr",
                             labels = c("b", "c", "d"), label_size = 9,
                             rel_heights = c(2.4, 1, 1))

fig4_block_scatter <- plot_grid(fig4_e1, fig4_e2, fig4_e_legend,
                                ncol = 3, labels = c("e", "", ""), label_size = 9,
                                rel_widths = c(1, 1, 0.5))

figure4 <- plot_grid(fig4_block_wide, fig4_block_zoom, fig4_block_scatter,
                     ncol = 1, rel_heights = c(1.2, 4.0, 2.0))

ggsave(figure4, filename = fig_path("thrb_panel"),
       height = 185, width = 110, units = "mm")


# =============================================================================
# 7. SUPPLEMENTARY FIGURE 25
# =============================================================================
# Its own region -- the full THRB transcript span, no padding. Because this is a
# value rather than a reassignment of s/e/padding, building it cannot disturb
# Figure 4 above, and the two sections can be run in either order.
region_supp <- make_region("chr19", start = 6344165, end = 6477279, padding_mb = 0)
THRB_MISSENSE_POS <- 6359267

## 25a -- transcript models + chi-square over the THRB span ---------------
panel_transcripts <- function(ensembl_tbl, gene_name, region, fill_colour) {
  dat <- ensembl_tbl %>% filter(Gene.name == gene_name)
  if (nrow(dat) == 0) {
    warning("no transcripts for '", gene_name, "'; Gene.name values present: ",
            paste(sort(unique(ensembl_tbl$Gene.name)), collapse = ", "), call. = FALSE)
  }
  ggplot(dat) +
    geom_vline(xintercept = THRB_MISSENSE_POS / 1e6, colour = "pink") +
    # full transcript span, then exons drawn on top with outlines
    geom_rect(aes(xmin = min / 1e6, xmax = max / 1e6, ymin = n - 0.1, ymax = n + 0.1),
              fill = fill_colour) +
    geom_rect(aes(xmin = as.numeric(Exon.region.start..bp.) / 1e6,
                  xmax = as.numeric(Exon.region.end..bp.) / 1e6,
                  ymin = n - 0.1, ymax = n + 0.1),
              fill = fill_colour, colour = "black") +
    coord_cartesian(xlim = region$xr_wide) +
    theme_classic() +
    theme(legend.position = "none", axis.text = element_blank(),
          axis.title = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank())
}

supp_chisq <- function(region, candidate) {
  dat  <- baltic_spring_vs_atlantic_spring_chiseq %>%
    filter(CHROM == region$chrom,
           POS >= region$start - region$padding_mb,
           POS <= region$end   + region$padding_mb) %>%
    mutate(class = ifelse(POS == THRB_MISSENSE_POS, "missense",
                          as.character(-log10(chisq_p) > 9.1)))
  

  # Named colour vector: the old version relied on the three levels ("FALSE",
  # "TRUE", "red") happening to sort into the same order as an unnamed
  # values = c("gray","black","red"). Naming them removes that dependence.
  ggplot(dat) +
    geom_rect(data = candidate, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0, ymax = 45),
              fill = "azure3", colour = NA, alpha = 0.5) +
    geom_vline(xintercept = THRB_MISSENSE_POS / 1e6, colour = "pink") +
    geom_point(data = dat %>% filter(class != "missense"),
               aes(x = POS / 1e6, y = -log10(chisq_p), colour = class), size = 0.5) +
    geom_point(data = dat %>% filter(class == "missense"),
               aes(x = POS / 1e6, y = -log10(chisq_p), colour = class), size = 0.8) +
    geom_point(data = dat %>% filter(class == "missense"),
               aes(x = POS / 1e6, y = -log10(chisq_p), colour = class), size = 0.8) +
    geom_text(data = dat %>% filter(class == "missense"),
              aes(x = POS / 1e6, y = -log10(chisq_p),
                  label = paste0(region$chrom, ":", POS)),
              size = 1.8, hjust = 1, vjust = -0.4, colour = "red") +
    scale_colour_manual(values = c("FALSE" = "gray", "TRUE" = "black",
                                   "missense" = "red")) +
    coord_cartesian(xlim = region$xr_wide, clip = "off") +
    ylab(expression(paste(-log[10], " (P value)"))) +
    theme_track()
}

supp25a <- plot_grid(
  panel_transcripts(ensembl114, "THRB", region_supp, "lightgreen"),
  panel_transcripts(ensembl109, "thrb", region_supp, "lightblue"),
  supp_chisq(region_supp, thrb_candidate),
  nrow = 3, align = "v", axis = "lr",
  rel_heights = c(1, 1, 1.4)
)

ggsave(supp25a, filename = fig_path("supp25a_gene_models_chisq_20260903"),
       units = "mm", width = 90, height = 90)


## 25b -- THRB missense allele frequency, pool-seq populations ------------
thrb_snps <- moderate_mutations %>%
  filter(seqnames == "chr19") %>%
  pull(start) %>% sort() %>% unique()

message("Supplementary 25b/c: ", length(thrb_snps), " chr19 missense SNP(s): ",
        paste(thrb_snps, collapse = ", "))

thrb_af <- pops_freq_df %>%
  filter(CHROM == "chr19", POS %in% thrb_snps) %>%
  pivot_longer(cols = -c(CHROM, POS),
               names_to = "populations", values_to = "frequencies") %>%
  mutate(POS_f = factor(as.character(POS), levels = as.character(thrb_snps)))

# =============================================================================
# 7b/7c -- horizontal (transposed) heatmaps
# =============================================================================
# Populations/individuals run along x, SNPs along y. With a single missense SNP
# each heatmap is one row tall -- a thin block. Note the level order is NOT
# reversed here: on an x axis levels[1] sits at the LEFT, so pool_order row 1
# (Pacific herring) leads naturally. (On the y axis it was the opposite, hence
# the rev() used elsewhere.)
af_levels_x <- pool_order_blocked$population
af_runs_x   <- group_runs(af_levels_x, pool_order_blocked$raw_group)
af_runs_x$label  <- unname(pool_group_display[af_runs_x$group])
af_runs_x$mid_i  <- (af_runs_x$start_i + af_runs_x$end_i) / 2
n_pools     <- length(af_levels_x)

GROUP_HEADER_ANGLE <- 45     # 0 = horizontal, 90 = vertical

# Group names sit in a header strip above the block, bracketed over the columns
# they cover, because the x axis itself is taken up by the individual pool names.
# x-domain must match the heatmap's: scale_x_discrete(expand = c(0,0)) puts n
# columns in a panel running 0.5 to n+0.5.
pool_group_header <- ggplot(af_runs_x) +
  geom_segment(aes(x = start_i - 0.4, xend = end_i + 0.4, y = 0, yend = 0),
               linewidth = 0.4) +
  geom_text(aes(x = mid_i, y = 0.12, label = label),
            angle = GROUP_HEADER_ANGLE, hjust = 0, vjust = 0, size = 2.1,  # was 1.8
            lineheight = 0.8) +
  scale_x_continuous(limits = c(0.5, n_pools + 0.5), expand = c(0, 0)) +
  ylim(-0.05, 0.95) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 12, b = 0, l = 1, unit = "mm"))

supp25b <- ggplot(thrb_af) +
  geom_tile(aes(x = factor(populations, levels = af_levels_x),
                y = POS_f, fill = frequencies)) +
  geom_vline(xintercept = head(af_runs_x$boundary, -1), colour = "white",
             linewidth = 0.3) +
  scale_fill_viridis(direction = -1, name = "Allele\nfrequency",
                     breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x     = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5,
                                       colour = "black"),          # was 3.5
        axis.text.y     = element_text(size = 6, colour = "black"), # was 5
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        legend.position = "right",
        legend.title    = element_text(size = 6),                   # was 5
        legend.text     = element_text(size = 6),                   # was 4
        legend.key.width  = unit(2.5, "mm"),
        legend.key.height = unit(3.5, "mm"))

panel_b <- plot_grid(pool_group_header, supp25b, ncol = 1, align = "v",
                     axis = "lr", rel_heights = c(0.3, 1))


## 25c -- THRB missense genotypes, individuals + sprat --------------------
# Individual names in the genotype matrix are `id` values -- the VCF was
# re-headered, note "newID" in its filename -- NOT `name_in_vcf`. 49 of the 125
# differ between those two columns, so joining on the wrong one loses them.
#
# Set drop_balsfjord = TRUE to exclude the four Balsfjord hybrids, matching the
# Figure 3 genotype heatmap. The previous version computed a filtered object but
# then plotted the unfiltered one, so the two figures disagreed.

# --- genotypes, transposed, with the Figure 3 legend ------------------------
# 0 = yellow, 1 = teal, 2 = dark blue, NA = white: viridis' own three-level
# values with 0 and 2 swapped, so the palette matches the previous
# scale_fill_viridis(discrete = TRUE) figure. NA is made an explicit factor
# level -- that is what gives it a legend key at all.
GENOTYPE_LEVELS <- c("0", "1", "2", "NA")
GENOTYPE_COLORS <- c("0"  = "#FDE725",
                     "1"  = "#21918C",
                     "2"  = "#440154",
                     "NA" = "white")

GENOTYPE_LABELS <- c("Homozygote AA", "Heterozygote", "Homozygote BB", "Missing")


build_genotype_panel_h <- function(positions, drop_balsfjord = FALSE,
                                   reverse = TRUE, header_size = 2.1) {
  balsfjord <- c("HWS61_Balsfjord_Atlantic", "HWS62_Balsfjord_Atlantic",
                 "HWS63_Balsfjord_Atlantic", "HWS64_Balsfjord_Atlantic")
  
  ind_levels <- as.character(individual_order[, 1])
  if (drop_balsfjord) ind_levels <- setdiff(ind_levels, balsfjord)
  if (reverse)        ind_levels <- rev(ind_levels)
  n_ind <- length(ind_levels)
  
  # population AND species come from the same join, so the two strips can never
  # disagree about which individual is where
  ind_meta <- data.frame(id = ind_levels, stringsAsFactors = FALSE) %>%
    left_join(sample_info %>% select(id, bars_plot, species), by = "id") %>%
    mutate(bars_plot = case_when(id == "SPRAT" ~ "SPRAT", TRUE ~ bars_plot),
           species   = case_when(id == "SPRAT" ~ "SPRAT", TRUE ~ species),
           group_display = unname(GENOTYPE_GROUP_DISPLAY[bars_plot]))
  
  unmapped <- unique(ind_meta$bars_plot[is.na(ind_meta$group_display)])
  if (length(unmapped)) {
    stop("bars_plot values with no entry in GENOTYPE_GROUP_DISPLAY: ",
         paste(unmapped, collapse = ", "))
  }
  
  # run-length encode the COLLAPSED label, not bars_plot -- otherwise Vancouver
  # and Japan come out as two adjacent runs both labelled "Pacific herring",
  # with a separator line drawn between them
  runs_pop <- group_runs(ind_levels, ind_meta$group_display)
  runs_sp  <- group_runs(ind_levels, ind_meta$species)
  
  # Alternating grey blocks make the group boundaries unmistakable without
  # introducing a colour scheme that needs explaining.
  runs_pop <- runs_pop %>% mutate(fill_alt = rep_len(c("a", "b"), n()))
  
  pop_strip <- ggplot(runs_pop) +
    geom_rect(aes(xmin = start_i - 0.5, xmax = end_i + 0.5,
                  ymin = 0, ymax = 1, fill = fill_alt),
              colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = c(a = "grey85", b = "grey60")) +
    scale_x_continuous(limits = c(0.5, n_ind + 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none")
  
  
  # alternate the species labels between two rows: "European sprat" sits over a
  # single column and would otherwise collide with "Pacific herring" next to it
  runs_sp <- runs_sp %>%
    mutate(label = unname(species_display[group]),
           mid_i = (start_i + end_i) / 2,
           row   = rep_len(c(0, 1), n()),
           lab_y = ifelse(row == 0, 0.30, 0.72))
  
  # x-domain must match the heatmap's: scale_x_discrete(expand = c(0,0)) puts
  # n columns in a panel running 0.5 to n+0.5
  sp_header <- ggplot(runs_sp) +
    geom_segment(aes(x = start_i - 0.4, xend = end_i + 0.4, y = 0, yend = 0),
                 linewidth = 0.8) +           # was 0.4
    geom_segment(aes(x = mid_i, xend = mid_i, y = 0.03, yend = lab_y - 0.05),
                 linewidth = 0.15, colour = "grey50") +
    geom_text(aes(x = mid_i, y = lab_y, label = label),
              size = header_size, hjust = 0.5, vjust = 0) +
    scale_x_continuous(limits = c(0.5, n_ind + 0.5), expand = c(0, 0)) +
    ylim(-0.05, 1.1) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 1, r = 6, b = 0, l = 1, unit = "mm"))
  
  dat <- genotypes_all %>%
    filter(CHROM == "chr19", POS %in% positions) %>%
    pivot_longer(cols = -c(CHROM, POS),
                 names_to = "Individuals", values_to = "Genotype") %>%
    filter(Individuals %in% ind_levels) %>%
    mutate(Genotype_f = factor(ifelse(is.na(Genotype), "NA", as.character(Genotype)),
                               levels = GENOTYPE_LEVELS),
           POS_f      = factor(as.character(POS), levels = as.character(positions)))
  
  heatmap <- ggplot(dat) +
    geom_tile(aes(x = factor(Individuals, levels = ind_levels),
                  y = POS_f, fill = Genotype_f)) +
    geom_vline(xintercept = head(runs_pop$boundary, -1), colour = "black",
               linewidth = 0.3) +
    scale_fill_manual(values = GENOTYPE_COLORS, labels = GENOTYPE_LABELS,
                      drop = FALSE, name = NULL) +
    scale_x_discrete(breaks = runs_pop$mid_row,
                     labels = runs_pop$group,      # already display text
                     expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x  = element_text(size = 6, angle = 45, hjust = 1, vjust = 1,
                                      colour = "black", lineheight = 0.8),
          axis.text.y  = element_text(size = 6, colour = "black"),
          axis.ticks.x = element_line(colour = "grey40", linewidth = 0.2),
          axis.ticks.length.x = unit(1.5, "mm"),
          axis.ticks.y = element_blank(),
          axis.line    = element_blank(),
          legend.position = "bottom",
          legend.text     = element_text(size = 6),
          legend.key.size = unit(3, "mm"),
          legend.key      = element_rect(colour = "grey40", linewidth = 0.2)) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(colour = "grey40")))
  
  plot_grid(sp_header, pop_strip, heatmap, ncol = 1, align = "v", axis = "lr",
            rel_heights = c(1, 0.25, 3))
}

panel_c <- build_genotype_panel_h(thrb_snps, drop_balsfjord = FALSE, reverse = TRUE)

ggsave(panel_b, filename = fig_path("tmp_b"), width = 180, height = 70, units = "mm")
ggsave(panel_c, filename = fig_path("tmp_c"), width = 180, height = 35, units = "mm")

# =============================================================================
# 8. SUPPLEMENTARY FIGURE 25 -- assembled
# =============================================================================
# b and c are stacked but NOT aligned to each other: b has 59 pools, c has 126
# individuals, and they are different sample sets (Han et al. pool-seq vs the
# WGS panel), so there is no column correspondence to line up.
supp25 <- plot_grid(
  supp25a, panel_b, panel_c,
  ncol        = 1,
  labels      = c("a", "b", "c"),
  label_size  = 9,
  rel_heights = c(90, 70, 55)      # same numbers, as proportions
)

ggsave(supp25, filename = fig_path("supp25_THRB_missense_combined"),
       width = 180, height = 210, units = "mm")


