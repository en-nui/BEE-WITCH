# --- setup / config --------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(glue)
library(scales)
library(viridisLite)
library(tibble)
library(tidyr)
library(patchwork)
library(MASS)
library(reshape2)
library(stats)
library(purrr)
library(rstatix)
library(tidyr)
library(stringr)
library(ggpubr)
library(minpack.lm)
library(parallel)
library(forcats)
# Filepaths (edit if needed)
SWEEP_FILE <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/sweep_gene_merged.csv"
GENE_SUM_FILE <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/final_analysis_all_colonies_aggregated_polymorphic_trajectory_SNVs/master_per_gene_summary.csv"
ANN_RDS <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/gene_level_annotations_all_colonies.rds"
TAX_FILE <- "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv"
OUTDIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# Season map (consistent)
season_map <- c("May"="Spring","June"="Spring",
                "July"="Summer","August"="Summer","September"="Summer",
                "October"="Fall","November"="Fall",
                "January"="Winter","February"="Winter")

# Minimal helper to safely unlist list-columns (ann)
safe_unlist <- function(x){
  if (is.null(x)) return(NA_character_)
  y <- unlist(x)
  if (length(y)==0) return(NA_character_)
  paste(unique(as.character(y)), collapse=";;")
}

# --- load data -------------------------------------------------------------
message("Loading data (may take a bit)...")
sweep_dt <- fread(SWEEP_FILE)               # sweep -> gene mapping (should include Month_t1/Month_t2 or t1/t2)
gene_dt  <- fread(GENE_SUM_FILE)            # per-gene summary (pN/pS, pi, etc)
tax_dt   <- fread(TAX_FILE)                 # MAG -> Genus
ann_dt   <- readRDS(ANN_RDS)                # annotations (list-columns)

# canonicalize types
for (c in c("MAG_id","colony_id","corresponding_gene_call","Month")) {
  if (c %in% names(gene_dt)) gene_dt[[c]] <- as.character(gene_dt[[c]])
  if (c %in% names(sweep_dt)) sweep_dt[[c]] <- as.character(sweep_dt[[c]])
}
tax_dt[, MAG_id := as.character(MAG_id)]

# collapse annotation list-columns into compact strings (for export)
if (!is.data.table(ann_dt)) ann_dt <- as.data.table(ann_dt)
if ("corresponding_gene_call" %in% names(ann_dt)) ann_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]
if (!"GO_str" %in% names(ann_dt)) {
  ann_dt[, GO_str := sapply(GO_slim_name_list, safe_unlist)]
  ann_dt[, pfam_str := sapply(pfam_label_list, safe_unlist)]
}

# attach Genus to gene_dt and sweep_dt (minimal columns)
tax_dt_unique <- unique(tax_dt[, .(MAG_id, Genus)])
gene_dt  <- merge(gene_dt,  tax_dt_unique, by = "MAG_id", all.x = TRUE)
sweep_dt <- merge(sweep_dt, tax_dt_unique, by = "MAG_id", all.x = TRUE)

# add Season fields (for sweep t1/t2 we use Month_t1/Month_t2 if available)
if ("Month" %in% names(gene_dt)) gene_dt[, Season := season_map[Month]]
if ("Month_t1" %in% names(sweep_dt)) sweep_dt[, Season_t1 := season_map[Month_t1]]
if ("Month_t2" %in% names(sweep_dt)) sweep_dt[, Season_t2 := season_map[Month_t2]]

# attach annotation strings to gene_dt / sweep_dt (left joins)
ann_sub <- ann_dt[, .(corresponding_gene_call, GO_str, pfam_str)]
gene_dt  <- merge(gene_dt,  ann_sub, by = "corresponding_gene_call", all.x = TRUE)
sweep_dt <- merge(sweep_dt, ann_sub, by = "corresponding_gene_call", all.x = TRUE)

# --- Build parallelism lists -----------------------------------------------
message("Building parallelism lists...")

# 1) pN/pS > 1 events per gene -> collect which MAGs a gene shows diversifying signal in
pnps_events <- gene_dt[!is.na(pN_pS_ratio) & pN_pS_ratio > 1,
                       .(MAGs = list(unique(MAG_id)), Genus_list = list(unique(Genus)), Colonies = list(unique(colony_id))),
                       by = corresponding_gene_call]

# counts
pnps_events[, n_MAG := lengths(MAGs)]
pnps_events[, n_Genus := lengths(Genus_list)]
pnps_events[, n_Col := lengths(Colonies)]

# parallel lists (you wanted two ways - here are three useful cutoffs)
parallel_pnps_A <- pnps_events[n_MAG >= 2 & n_Genus >= 2 & n_Col >= 2]  # strict: 2+ MAGs, 2+ genera, 2+ colonies
parallel_pnps_B <- pnps_events[n_MAG >= 2 & n_Genus >= 2]               # 2+ MAGs + 2+ genera
parallel_pnps_C <- pnps_events[n_MAG >= 2 & n_Col >= 2]                 # 2+ MAGs + 2+ colonies
parallel_pnps_D <- pnps_events[n_MAG >= 2]               # 2+ MAGs


fwrite(parallel_pnps_A, file = file.path(OUTDIR, "parallel_pnps_A_genes.csv"))
fwrite(parallel_pnps_B, file = file.path(OUTDIR, "parallel_pnps_B_genes.csv"))
fwrite(parallel_pnps_C, file = file.path(OUTDIR, "parallel_pnps_C_genes.csv"))

# 2) 0D sweep events: collapse sweep_dt to gene-level events (we expect degeneracy column)
zeroD_dt <- sweep_dt[degeneracy == "0D"]

# which MAGs each gene has a 0D sweep in
zeroD_events <- zeroD_dt[, .(MAGs = list(unique(MAG_id)), Genus_list = list(unique(Genus)), Colonies = list(unique(colony_id))),
                         by = corresponding_gene_call]
zeroD_events[, n_MAG := lengths(MAGs)]
zeroD_events[, n_Genus := lengths(Genus_list)]
zeroD_events[, n_Col := lengths(Colonies)]

parallel_zeroD_A <- zeroD_events[n_MAG >= 2 & n_Genus >= 2 & n_Col >= 2]
parallel_zeroD_B <- zeroD_events[n_MAG >= 2 & n_Genus >= 2]
parallel_zeroD_C <- zeroD_events[n_MAG >= 2 & n_Col >= 2]
parallel_zeroD_D <- zeroD_events[n_MAG >= 2]               # 2+ MAGs



fwrite(parallel_zeroD_A, file = file.path(OUTDIR, "parallel_zeroD_A_genes.csv"))
fwrite(parallel_zeroD_B, file = file.path(OUTDIR, "parallel_zeroD_B_genes.csv"))
fwrite(parallel_zeroD_C, file = file.path(OUTDIR, "parallel_zeroD_C_genes.csv"))

# Also write counts summary (how many parallel genes in each category)
summary_dt <- data.table(
  test = c("pnps_A","pnps_B","pnps_C","zeroD_A","zeroD_B","zeroD_C"),
  n_genes = c(nrow(parallel_pnps_A), nrow(parallel_pnps_B), nrow(parallel_pnps_C),
              nrow(parallel_zeroD_A), nrow(parallel_zeroD_B), nrow(parallel_zeroD_C))
)
fwrite(summary_dt, file = file.path(OUTDIR,"parallel_summary_counts.csv"))

# --- Prepare matrix for the proportional-season tile plot -------------------
message("Preparing proportional-season tile data...")

# pick sweep (0D) or pnps-based matrix: user-specified choice (we'll produce for 0D sweeps as requested)
# We'll create a table: per (MAG_id, GO_slim_name) count of genes with 0D sweeps in that GO,
# and the season breakdown (proportion of swept genes in each season).
# note: our ann strings use 'GO_str' which can contain multiple GO_slm separated by ';;'
zeroD_dt
# explode zeroD_dt gene->GO rows (one GO per row)
zeroD_dt2 <- zeroD_dt[!is.na(GO_str), .(corresponding_gene_call, MAG_id, Season = Season_t1, GO = unlist(strsplit(GO_str, ";;")))]
# drop empties
zeroD_dt2 <- zeroD_dt2[!is.na(GO) & GO != ""]
mag_go
# compute counts by MAG x GO x Season
mag_go_season <- zeroD_dt2[, .(n_genes = uniqueN(corresponding_gene_call)), by = .(MAG_id, GO, Season)]
# total per MAG x GO
mag_go_total <- mag_go_season[, .(total = sum(n_genes)), by = .(MAG_id, GO)]
mag_go <- merge(mag_go_season, mag_go_total, by = c("MAG_id","GO"))
mag_go[, prop := n_genes / total]  # proportion within this MAG_id x GO cell that occur in this Season

# filter to GO categories that appear in >= 2 MAGs (parallelizable across genera later)
go_mag_counts <- mag_go_total[, .(n_MAGs = uniqueN(MAG_id)), by = GO]
parallel_GO_candidates <- go_mag_counts[n_MAGs >= 2]$GO

plot_dt <- mag_go[GO %in% parallel_GO_candidates]

# for plotting, we will create one tile per MAG_id x GO at y = GO, x = MAG_id.
# inside that tile we draw 4 sub-rectangles (stacked vertically) whose heights equal the season proportions.
# to do this we need cumulative positions per cell.

# set month/season order consistently
season_levels <- c("Spring","Summer","Fall","Winter")
plot_dt[, Season := factor(Season, levels = season_levels, ordered = TRUE)]

# build coordinates: for each MAG x GO produce cumulative ymin/ymax for each season slice
plot_dt <- plot_dt[order(MAG_id, GO, Season)]
plot_dt[, ymin := cumsum(shift(prop, fill = 0)), by = .(MAG_id, GO)]
plot_dt[, ymax := ymin + prop, by = .(MAG_id, GO)]

# create small tile coordinate system:
# create x integer index for MAGs (so we can control tile width)
mag_order <- unique(plot_dt$MAG_id)
mag_index <- data.table(MAG_id = mag_order, x = seq_along(mag_order))
plot_dt <- merge(plot_dt, mag_index, by = "MAG_id", all.x = TRUE)

# y will be factor GO; we map GO to integer index
go_order <- unique(plot_dt$GO)
go_index <- data.table(GO = go_order, y = seq_along(go_order))
plot_dt <- merge(plot_dt, go_index, by = "GO", all.x = TRUE)

# tile size controls (these control how big the cell is)
cell_width <- 0.9
cell_height <- 0.9

# for each slice compute rectangle coordinates in plot coordinates
plot_dt[, xmin := x - cell_width/2]
plot_dt[, xmax := x + cell_width/2]
# convert ymin/ymax (0..1) into plot y-range centered at y
plot_dt[, slice_ymin := (y - cell_height/2) + (ymax - 1) * cell_height]  # placeholder, we'll compute exact below
# Instead we compute slice_ymin/ymax as: y - cell_height/2 + (y - y)?? simpler to map 0..1 to [y-cellh/2, y+cellh/2]
plot_dt[, slice_ymin := (y - cell_height/2) + (y - y) * 0 + (ymin * cell_height)]
plot_dt[, slice_ymax := (y - cell_height/2) + (ymin * cell_height) + (prop * cell_height)]

# Add genus info for ordering/grouping by genus
plot_dt <- merge(plot_dt, tax_dt_unique, by = "MAG_id", all.x = TRUE)
plot_dt[, Genus := ifelse(is.na(Genus), "no_match", Genus)]

# For color mapping of seasons
season_colors <- c("Spring"="#66c2a5", "Summer"="#fc8d62", "Fall"="#8da0cb", "Winter"="#e78ac3")

# --- Plot: proportional-season tiles using geom_rect ------------------------
message("Plotting proportional-season tiles (0D sweeps) ...")
p <- ggplot() +
  geom_rect(data = plot_dt,
            aes(xmin = xmin, xmax = xmax, ymin = slice_ymin, ymax = slice_ymax, fill = Season),
            color = "white", size = 0.09) +
  scale_fill_manual(values = season_colors, na.value = "grey80") +
  scale_x_continuous(breaks = mag_index$x, labels = mag_index$MAG_id, expand = c(0,0)) +
  scale_y_continuous(breaks = go_index$y, labels = go_index$GO, expand = c(0,0)) +
  coord_fixed(ratio = 0.6) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(title = "0D-sweep genes: per MAG × GO, bar = seasonal proportion",
       x = "MAGs (not labeled individually)", y = "GO_slim pathways")





p
# thresholds
ALPHA_FDR <- 0.05
MIN_GENE_SAMPLES <- 2   # filter genes seen at least this many samples (optional)

# -------- LOAD --------
dt <- fread(GENE_SUM_FILE)

# expecting columns (per-sample rows): corresponding_gene_call, MAG_id, colony_id, Month,
# SNVs_4D_sites (count of 4D SNVs in that gene/sample), gene_length_4D_sites (4D site count),
# pN_pS_ratio (per-sample or NA)
# If names differ, adjust accordingly.

# Quick checks:
stopifnot("corresponding_gene_call" %in% names(dt))
if (!("SNVs_4D_sites" %in% names(dt) || "SNVs_4D" %in% names(dt))) {
  stop("Expected column named SNVs_4D_sites (or similar). Adjust script to match your file.")
}
# normalize column names
if ("SNVs_4D" %in% names(dt) && !("SNVs_4D_sites" %in% names(dt))) setnames(dt, "SNVs_4D", "SNVs_4D_sites")
if (!"effective_length_4D_sites" %in% names(dt)) {
  # fallback: maybe gene_length_all_sites or gene_length_4D_sites; check
  if ("full_gene_length_all_sites" %in% names(dt)) {
    stop("No full_gene_length_4D_sites column found. Please provide 4D site counts per gene.")
  } else {
    stop("need gene_length_4D_sites column")
  }
}

# restrict to genes observed (optionally filter low-sample genes)
per_gene_nobs <- dt[, .N, by=.(corresponding_gene_call)]
keep_genes <- per_gene_nobs[N >= MIN_GENE_SAMPLES, corresponding_gene_call]

dt <- dt[corresponding_gene_call %in% keep_genes]

# ---------- Estimate neutral per-site mutation rate ----------
# We'll estimate mu_site = total_4D_SNV_observations / total_4D_site_observations
# where:
# - total_4D_SNV_observations = sum over all sample-rows of SNVs_4D_sites
# - total_4D_site_observations = sum over all sample-rows of gene_length_4D_sites
# This gives the per-site *per-sample* mutation rate (observed).
# Expected count for a gene across its observed samples = mu_site * gene_length_4D_sites * n_samples_gene

total_4D_snvs <- dt[, sum(as.numeric(SNVs_4D_sites), na.rm=TRUE)]
total_4D_site_obs <- dt[, sum(as.numeric(effective_length_4D_sites), na.rm=TRUE)]

mu_site_per_sample <- total_4D_snvs / total_4D_site_obs

cat("Total 4D SNVs:", total_4D_snvs, "\n")
cat("Total 4D site observations:", total_4D_site_obs, "\n")
cat("Estimated mu_site_per_sample:", mu_site_per_sample, "\n")

# ---------- Per-gene observed & expected ----------
# observed_total_4D_snvs per gene (summed across samples)
gene_obs <- dt[, .(
  obs_4D = sum(as.numeric(SNVs_4D_sites), na.rm=TRUE),
  n_samples = .N,
  gene_length_4D = as.numeric(first(effective_length_4D_sites)),  # assume constant per gene
  mean_pnps = mean(as.numeric(pN_pS_ratio), na.rm=TRUE),     # per gene mean pN/pS
  median_pnps = median(as.numeric(pN_pS_ratio), na.rm=TRUE)
), by = corresponding_gene_call]

# expected lambda per gene (Poisson) across observed samples:
# lambda = mu_site_per_sample * gene_length_4D * n_samples
gene_obs[, expected_lambda := mu_site_per_sample * gene_length_4D * n_samples]

# Guard against zeros or NA gene_length
gene_obs <- gene_obs[!is.na(gene_length_4D) & gene_length_4D > 0]

# ---------- Poisson test (upper tail) ----------
# p-value = P(Poisson(lambda) >= obs_4D) = 1 - ppois(obs_4D - 1, lambda)
gene_obs[, pval_poisson := 1 - ppois(pmax(0, obs_4D - 1), lambda = expected_lambda)]
# handle small lambdas numerically stable
gene_obs[expected_lambda == 0 & obs_4D > 0, pval_poisson := 0]

# FDR correction
gene_obs[, qval_poisson := p.adjust(pval_poisson, method="BH")]

# Flag mutation-prone genes (excess 4D SNVs)
gene_obs[, high_mutation := qval_poisson <= ALPHA_FDR]

# Save table
fwrite(gene_obs, file=file.path(OUTDIR, "gene_neutral_model_poisson_tests.csv"))

# Summary counts
cat("Genes tested:", nrow(gene_obs), "\n")
cat("High-mutation genes (FDR <= ", ALPHA_FDR, "): ", sum(gene_obs$high_mutation, na.rm=TRUE), "\n", sep="")

# ---------- Diagnostics & Overdispersion test ----------
# check distribution of obs_4D vs expected_lambda
# Basic overdispersion: var(obs)/mean(obs)
obs_mean <- mean(gene_obs$obs_4D, na.rm=TRUE)
obs_var  <- var(gene_obs$obs_4D, na.rm=TRUE)
cat("Observed mean obs_4D:", obs_mean, " variance:", obs_var, " var/mean:", obs_var/obs_mean, "\n")

# If var >> mean, Poisson is underfitting; consider NB model or permutation null.

# ---------- Compare pN/pS distributions ----------
# Build two groups: high_mutation vs rest
grpA <- gene_obs[high_mutation == TRUE]
grpB <- gene_obs[high_mutation == FALSE]

# drop NA pN/pS
grpA_pnps <- na.omit(grpA$mean_pnps)
grpB_pnps <- na.omit(grpB$mean_pnps)

# Wilcoxon rank-sum (one-sided if directional)
wil <- wilcox.test(grpA_pnps, grpB_pnps, alternative = "greater", exact = FALSE)
cat("Wilcoxon test result:\n"); print(wil)

# Permutation test (empirical): permute labels many times to get p
perm_test <- function(x, y, nperm = 5000){
  obs_stat <- median(x, na.rm=TRUE) - median(y, na.rm=TRUE)
  combined <- c(x, y)
  n1 <- length(x)
  perm_stats <- replicate(nperm, {
    s <- sample(combined, length(combined), replace = FALSE)
    median(s[1:n1], na.rm=TRUE) - median(s[(n1+1):length(combined)], na.rm=TRUE)
  })
  p_emp <- mean(perm_stats >= obs_stat)
  list(obs_stat = obs_stat, p_emp = p_emp)
}
pt <- perm_test(grpA_pnps, grpB_pnps, nperm = 5000)
cat("Permutation median diff (A - B):", pt$obs_stat, " p_emp:", pt$p_emp, "\n")

# Cliff's delta (effect size)
cliff_delta <- function(x, y){
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  n_x <- length(x); n_y <- length(y)
  gt <- sum(outer(x, y, FUN = ">"))
  lt <- sum(outer(x, y, FUN = "<"))
  delta <- (gt - lt) / (n_x * n_y)
  return(list(delta = delta, n_x = n_x, n_y = n_y, gt = gt, lt = lt))
}
cd <- cliff_delta(grpA_pnps, grpB_pnps)
cat("Cliff's delta:", cd$delta, " (nA=", cd$n_x, ", nB=", cd$n_y, ")\n", sep="")

# Optional rstatix wilcox effsize
if ("wilcox_effsize" %in% ls("package:rstatix")) {
  # Build a tidy table
  tdf <- data.table(
    mean_pnps = c(grpA_pnps, grpB_pnps),
    group = c(rep("high_mut", length(grpA_pnps)), rep("rest", length(grpB_pnps)))
  )
  eff <- wilcox_effsize(tdf, mean_pnps ~ group, alternative = "greater", ci = TRUE)
  print(eff)
}

# ---------- Plots ----------
# 1) histogram of expected vs observed 4D counts
p1 <- ggplot(gene_obs, aes(x = obs_4D)) +
  geom_histogram(bins = 80) + labs(title = "Observed total 4D SNVs per gene")

p1
ggsave(file.path(OUTDIR, "hist_obs_4D_per_gene.png"), p1, width = 7, height = 4)

p2 <- ggplot(gene_obs, aes(x = expected_lambda)) +
  geom_histogram(bins = 80) + labs(title = "Expected 4D SNVs per gene (Poisson lambda)")
p2

ggsave(file.path(OUTDIR, "hist_expected_lambda_per_gene.png"), p2, width = 7, height = 4)

# 2) ECDF / violin of mean_pnps for high vs rest
gene_obs[, group := ifelse(high_mutation, "high_mut", "rest")]
p3 <- ggplot(gene_obs[!is.na(mean_pnps)], aes(x = group, y = mean_pnps)) +
  geom_violin(trim = TRUE) + geom_boxplot(width = 0.1) + labs(title = "pN/pS by mutation-prone genes")
p3

ggsave(file.path(OUTDIR, "violin_pnps_high_vs_rest.png"), p3, width = 6, height = 4)

# 3) ECDF plot
p4 <- ggplot(gene_obs[!is.na(mean_pnps)], aes(x = mean_pnps, color = group)) +
  stat_ecdf() + labs(title = "ECDF of mean pN/pS: high_mut vs rest")
p4

ggsave(file.path(OUTDIR, "ecdf_pnps_high_vs_rest.png"), p4, width = 6, height = 4)

# Save summary stats
summary_dt <- data.table(
  group = c("high_mut", "rest"),
  n = c(nrow(grpA), nrow(grpB)),
  median_pnps = c(median(grpA_pnps, na.rm=TRUE), median(grpB_pnps, na.rm=TRUE)),
  mean_pnps = c(mean(grpA_pnps, na.rm=TRUE), mean(grpB_pnps, na.rm=TRUE))
)
fwrite(summary_dt, file = file.path(OUTDIR, "pnps_summary_high_vs_rest.csv"))

# Save high-mutation gene list for parallelism analysis
fwrite(gene_obs[high_mutation == TRUE], file = file.path(OUTDIR, "high_mutation_genes_poissonFDR.csv"))

# ---------- OPTIONAL: Empirical null by permutation across genes ----------
# If Poisson seems underpowered (var >> mean), do:
# - Shuffle gene labels across SNV counts / lengths and recalc expected counts distribution,
# - Or fit negative-binomial (MASS::glm.nb) across genes with length * n_samples as offset.

# End
cat("Done. Outputs in:", OUTDIR, "\n")




# --- 1. Functional Enrichment of Parallel Genes ----------------------------
message("Running Fisher's Exact Test for Enrichment in Parallel Genes...")

# 1. Define Universe: All genes analyzed that have a GO annotation
# Explode the universe so one row per gene-GO pair
universe_dt <- gene_dt[!is.na(GO_str) & GO_str != "", .(corresponding_gene_call, GO_str)]
universe_dt <- universe_dt[, .(GO = unlist(strsplit(GO_str, ";;"))), by = corresponding_gene_call]
universe_unique <- unique(universe_dt$corresponding_gene_call)

# 2. Define Target Set: The Parallel Genes (e.g., from your strict list A)
# (Ensure we match IDs)
target_genes <- parallel_pnps_A$corresponding_gene_call # Or parallel_zeroD_A
target_genes <- target_genes[target_genes %in% universe_unique]

# 3. Contingency Table Calculation per GO term
# We only test GO terms that appear in the target set to save time
test_GOs <- universe_dt[corresponding_gene_call %in% target_genes, unique(GO)]

enrichment_res <- lapply(test_GOs, function(go_term) {
  # Genes with this GO
  genes_with_go <- universe_dt[GO == go_term, unique(corresponding_gene_call)]
  
  # A: Target & GO+
  a <- length(intersect(target_genes, genes_with_go))
  # B: Universe & GO+ (excluding target)
  b <- length(setdiff(genes_with_go, target_genes))
  # C: Target & GO-
  c <- length(setdiff(target_genes, genes_with_go))
  # D: Universe & GO- (excluding target)
  d <- length(setdiff(setdiff(universe_unique, genes_with_go), target_genes))
  
  mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  
  data.table(
    GO = go_term,
    p_val = ft$p.value,
    odds_ratio = ft$estimate,
    n_parallel_with_GO = a,
    n_universe_with_GO = a + b
  )
})

enrich_dt <- rbindlist(enrichment_res)
enrich_dt[, padj := p.adjust(p_val, method = "BH")]

# Filter for significant results for plotting
sig_enrich <- enrich_dt[padj < 0.05 & n_parallel_with_GO >= 2] # Filter as needed
fwrite(sig_enrich, file = file.path(OUTDIR, "parallel_GO_enrichment.csv"))

# --- Plotting the Enrichment ----------------------------------------------
p_enrich <- ggplot(sig_enrich[order(odds_ratio)], aes(x = odds_ratio, y = reorder(GO, odds_ratio))) +
  geom_segment(aes(x = 0, xend = odds_ratio, yend = GO), color = "grey50") +
  geom_point(aes(size = n_parallel_with_GO, color = padj)) +
  scale_color_viridis_c(direction = -1, option = "plasma", name = "FDR") +
  labs(title = "Functional Functions driving Parallel Evolution",
       x = "Odds Ratio (Enrichment)", y = NULL, size = "Count") +
  theme_minimal()
p_enrich
ggsave(file.path(OUTDIR, "parallel_enrichment_lollipop.png"), p_enrich, width = 8, height = 6)


# --- 3. Seasonality of Parallel Sweeps -------------------------------------
# Use zeroD_dt (sweeps) and filter for only the genes identified as parallel in parallel_zeroD_A

target_genes_Z <- parallel_zeroD_A$corresponding_gene_call
subset_sweeps <- sweep_dt[corresponding_gene_call %in% target_genes_Z & degeneracy == "0D"]

# Ensure Month is ordered factor
month_levs <- c("January","February","March","April","May","June",
                "July","August","September","October","November")
subset_sweeps[, Month_t1 := factor(Month_t1, levels = month_levs)]

# Aggregate
season_counts <- subset_sweeps[, .N, by = .(Month_t1, Genus)]

p_polar <- ggplot(season_counts, aes(x = Month_t1, y = N, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_polar(start = 0) +
  theme_minimal() +
  labs(title = "Timing of Parallel Sweep Events", x = NULL, y = "Number of Sweeps")
p_polar
ggsave(file.path(OUTDIR, "parallel_sweep_polar_plot.png"), p_polar, width = 6, height = 6)



library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
library(viridis)

# --- 0. Config -------------------------------------------------------------
OUTDIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"
# GENE_SUM_FILE defined in previous steps... assuming loaded as gene_dt

# Target Genera
target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium", 
                   "commensalibacter", "frischella", "gilliamella", 
                   "lactobacillus", "snodgrassella")
month_levels <- c("May", "June", "July", "August", "September", 
                  "October", "November", "December", "January", "February")

# --- 1. Prepare Data (Long Format) -----------------------------------------
message("Preparing Data for Combined Plot...")

# Filter and Aggregate to MAG level first
div_dt <- gene_dt[Genus %in% target_genera]
mag_summary <- div_dt[, .(
  pi_1D = mean(pi_1D, na.rm = TRUE),
  pi_4D = mean(pi_4D, na.rm = TRUE)
), by = .(MAG_id, Genus, Month)]

# Reshape to LONG format
# Explicitly use data.table::melt to avoid confusion, or setDT afterwards
mag_long <- melt(mag_summary, 
                 id.vars = c("MAG_id", "Genus", "Month"), 
                 measure.vars = c("pi_1D", "pi_4D"),
                 variable.name = "Metric", 
                 value.name = "Pi_Value")

# --- THE FIX: Force it to be a data.table ---
setDT(mag_long) 
# --------------------------------------------

# Clean up Metric names for the legend
mag_long[, Metric_Label := fcase(
  Metric == "pi_4D", "Synonymous (4D)",
  Metric == "pi_1D", "Non-synonymous (1D)"
)]

# Set Month Order and Numeric Value
mag_long[, Month := factor(Month, levels = month_levels)]
mag_long[, Month_Num := as.numeric(Month)]
mag_long <- mag_long[!is.na(Month_Num)]

# Filter Zeros/NAs for Log Plot
plot_data <- mag_long[Pi_Value > 0 & !is.na(Pi_Value)]

# --- 2. Plotting -----------------------------------------------------------
message("Generating Combined Faceted Plot...")

p_combined <- ggplot(plot_data, aes(x = Month_Num, y = Pi_Value, group = Metric_Label)) +
  
  # A. The Points (Raw MAG Data)
  # Color by Metric (Red for 1D, Blue for 4D style)
  geom_jitter(aes(color = Metric_Label), alpha = 0.2, size = 0.8, width = 0.2) +
  
  # B. The Trend Lines (LOESS)
  # Linetype: 1D = Dotted, 4D = Solid
  geom_smooth(aes(color = Metric_Label, fill = Metric_Label, linetype = Metric_Label), 
              method = "loess", size = 1.2, alpha = 0.15) +
  
  # C. Faceting
  facet_wrap(~Genus, ncol = 4, scales = "fixed") +
  
  # D. Scales and Colors
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = 1:length(month_levels), labels = month_levels) +
  
  # Manual Colors: High contrast (e.g., Orange/Red for Functional, Blue/Grey for Neutral)
  scale_color_manual(values = c("Non-synonymous (1D)" = "#D55E00", "Synonymous (4D)" = "#0072B2")) +
  scale_fill_manual(values = c("Non-synonymous (1D)" = "#D55E00", "Synonymous (4D)" = "#0072B2")) +
  scale_linetype_manual(values = c("Non-synonymous (1D)" = "dotted", "Synonymous (4D)" = "solid")) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(title = "Comparative Evolutionary Trajectories: Function vs. Drift",
       subtitle = "Solid Line = Synonymous (Neutral Drift) | Dotted Line = Non-synonymous (Selection)",
       y = expression(log[10](Mean ~ pi)),
       x = NULL)
p_combined
ggsave(file.path(OUTDIR, "trend_combined_1D_4D_faceted.pdf"), p_combined, width = 12, height = 8)
message("Plot saved to: ", file.path(OUTDIR, "trend_combined_1D_4D_faceted.pdf"))


# Filter to target genera
stat_dt <- gene_dt[Genus %in% target_genera]

# --- 1. Define Statistical Function ----------------------------------------
# This function runs Kruskal-Wallis + Pairwise Wilcoxon and returns a plotting table
run_pairwise_stats <- function(data, metric_col, metric_name) {
  
  # A. Global Kruskal-Wallis Test
  # "Are there ANY differences among the 8 groups?"
  kw_res <- kruskal.test(formula(paste(metric_col, "~ Genus")), data = data)
  message(paste("Kruskal-Wallis for", metric_name, "- p-value:", format.pval(kw_res$p.value)))
  
  # B. Pairwise Wilcoxon Test (Non-parametric t-test equivalent)
  # Corrects for multiple testing using "BH" (Benjamini-Hochberg FDR)
  pwc <- pairwise.wilcox.test(data[[metric_col]], data$Genus, 
                              p.adjust.method = "BH", 
                              paired = FALSE)
  
  # C. Convert 'htest' object to a Tidy Data Frame for Plotting
  p_mat <- pwc$p.value
  p_df <- as.data.frame(as.table(p_mat))
  colnames(p_df) <- c("Genus_A", "Genus_B", "p_adj")
  p_df <- p_df[!is.na(p_df$p_adj), ] # Remove NA diagonal/upper triangle if any
  
  # Add significance stars for plotting
  p_df$significance <- cut(p_df$p_adj, 
                           breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                           labels = c("***", "**", "*", "ns"))
  
  p_df$metric <- metric_name
  return(p_df)
}

# --- 2. Run Tests for 4D and 1D --------------------------------------------
message("Running stats for Pi 4D (Synonymous)...")
stats_4d <- run_pairwise_stats(stat_dt, "pi_4D", "Synonymous (4D)")

message("Running stats for Pi 1D (Non-synonymous)...")
stats_1d <- run_pairwise_stats(stat_dt, "pi_1D", "Functional (1D)")

# Combine results
all_stats <- rbind(stats_4d, stats_1d)

# --- 3. Visualization: The Pairwise Significance Heatmap -------------------
# This visualizes the p-values. Darker = More Significant.

# Create a mirrored dataset so the heatmap is full (not just a triangle)
all_stats_mirror <- all_stats %>% rename(Genus_A = Genus_B, Genus_B = Genus_A)
all_stats_full <- rbind(all_stats, all_stats_mirror)

p_heatmap <- ggplot(all_stats_full, aes(x = Genus_A, y = Genus_B)) +
  geom_tile(aes(fill = p_adj), color = "white") +
  
  # Overlay significance stars
  geom_text(aes(label = significance), color = "white", size = 3, vjust = 0.8) +
  
  # Facet by Metric (1D vs 4D)
  facet_wrap(~metric, ncol = 2) +
  
  # Scales
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt", 
                       name = "FDR p-value", limits = c(0, 0.05), oob = scales::squish) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(title = "Statistical Differences in Diversity Between Genera",
       subtitle = "Pairwise Wilcoxon Tests (FDR adjusted) | *** p<0.001, ** p<0.01, * p<0.05",
       x = NULL, y = NULL)
p_heatmap
ggsave(file.path(OUTDIR, "stats_pairwise_heatmap.pdf"), p_heatmap, width = 12, height = 6)


# --- 4. Alternative: Boxplot with Letters (Compact Letter Display) ---------
# This is often preferred if you just want to show "Who is highest/lowest"
# We simply calculate median Pi per genus and plot it.

# Simple violin plot to accompany the stats
p_violin <- ggplot(stat_dt, aes(x = Genus, y = pi_1D, fill = Genus)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_fill_viridis_d(option = "turbo") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(title = "Distribution of Non-Synonymous Diversity (pi 1D)",
       y = expression(log[10](pi[1][D])), x = NULL)
p_violin
ggsave(file.path(OUTDIR, "stats_distribution_violin.pdf"), p_violin, width = 8, height = 6)




head(parallel_pnps_A)



# We use D because it is the superset (Union) of A, B, and C.
# Condition: n_MAG >= 2.
# Ensure 'parallel_pnps_D' is loaded from your previous step.
# If not, recreate it:
# parallel_pnps_D <- pnps_events[n_MAG >= 2] 

# --- 1. Expand Data (Gene -> Colony -> Function) ---------------------------
message("Processing Parallel Events with Pfam+GO Labels...")

# A. Unroll Colony List & Handle Genera List-Column
# We convert the list of genera to a string for robust grouping
expanded_colony <- parallel_pnps_D[, .(
  Colony = unlist(Colonies), 
  Genera_Str = sapply(Genus_list, function(x) paste(sort(unique(x)), collapse = ";")) 
), by = corresponding_gene_call]

# B. Get Annotations (Pfam + GO)
# We select relevant columns from gene_dt
ann_map <- unique(gene_dt[, .(corresponding_gene_call, pfam_str, GO_str)])
annotated_dt <- merge(expanded_colony, ann_map, by = "corresponding_gene_call")

# Filter out rows with missing annotations to keep the plot clean
annotated_dt <- annotated_dt[!is.na(pfam_str) & pfam_str != "" & !is.na(GO_str) & GO_str != ""]

# C. Create the Composite Label (Pfam, GO)
# We use a code block {} inside the data.table call to handle the splitting
heatmap_dt <- annotated_dt[, {
  # 1. Split the strings into vectors
  p_list <- unlist(strsplit(pfam_str, ";;"))
  g_list <- unlist(strsplit(GO_str, ";;"))
  
  # 2. Use CJ (Cross Join) to create all unique pairwise combinations
  # This handles the length mismatch (e.g., 2 Pfams vs 1 GO) automatically
  CJ(Pfam_Single = p_list, GO_Single = g_list)
  
}, by = .(corresponding_gene_call, Colony, Genera_Str)]

# Now this column creation will work because the columns exist
heatmap_dt[, Y_Label := paste0("(", Pfam_Single, ", ", str_remove(GO_Single, "^GO:\\d+~?"), ")")]

# --- 2. Aggregate per Cell (Colony x Y_Label) ------------------------------
cell_stats <- heatmap_dt[, .(
  n_genes = uniqueN(corresponding_gene_call),
  
  # Extract unique genera involved (splitting the string back)
  unique_genera = list(unique(unlist(strsplit(Genera_Str, ";"))))
), by = .(Colony, Y_Label)]

# Calculate Genera Counts
cell_stats[, n_genera_count := lengths(unique_genera)]

# Create Guide Label for Manual Icons
# "G" for single genus, "2 gen" for multi-genus
cell_stats[, Guide_Label := ifelse(n_genera_count > 1, 
                                   paste0(n_genera_count, " gen"), 
                                   substring(sapply(unique_genera, `[`, 1), 1, 1))]

# --- 3. Filtering (The "2 Colonies" Rule) ----------------------------------
# We ONLY keep Y-axis rows that appear in at least 2 different colonies.
label_counts <- cell_stats[, .(n_cols = uniqueN(Colony)), by = Y_Label]
valid_labels <- label_counts[n_cols >= 2]$Y_Label

plot_dt <- cell_stats[Y_Label %in% valid_labels]

# --- 4. Clustering & Ordering ----------------------------------------------
# Cluster the Y-axis based on occurrence pattern so similar functions group together
if(length(unique(plot_dt$Y_Label)) > 1) {
  mat <- dcast(plot_dt, Y_Label ~ Colony, value.var = "n_genes", fill = 0)
  mat_num <- as.matrix(mat[, -1])
  rownames(mat_num) <- mat$Y_Label
  
  dd.row <- as.dendrogram(hclust(dist(mat_num)))
  ordered_labels <- rownames(mat_num)[order.dendrogram(dd.row)]
  
  plot_dt[, Y_Label := factor(Y_Label, levels = ordered_labels)]
}

# --- 5. Plotting The Blueprint ---------------------------------------------
message("Generating Blueprint Heatmap...")

p_blueprint <- ggplot(plot_dt, aes(x = as.factor(Colony), y = Y_Label)) +
  
  # 1. The Tile (Intensity)
  geom_tile(aes(fill = n_genes), color = "white", size = 0.2) +
  
  # 2. The Highlight (Black Border = Multi-Genus Convergence)
  # This highlights cells where >1 Genus is driving the selection
  geom_tile(data = plot_dt[n_genera_count > 1], 
            aes(color = "Multi-Genus"), fill = NA, size = 1) +
  
  # 3. The Text Guide (For Icon Placement)
  geom_text(aes(label = Guide_Label), size = 2.5, color = "grey20", fontface = "bold") +
  
  # Scales
  scale_fill_viridis_c(option = "plasma", name = "Parallel Genes", direction = -1) +
  scale_color_manual(values = c("Multi-Genus" = "black"), labels = "Cross-Genera Event") +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8), # Smaller font for long labels
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(title = "Functional Parallelism Blueprint (Pfam + GO)",
       subtitle = "Rows = Functions in >1 Colony. Black Border = Multi-Genus Convergence.",
       x = "Colony", y = NULL)
p_blueprint





# Ensure inputs are loaded:
# 1. parallel_pnps_D (The Union of pN/pS events)
# 2. parallel_zeroD_D (The Union of 0D events - needed for overlap check)
# 3. gene_dt (For annotations)

# ---------------------------------------------------------------------------
# PART 1: The pN/pS Blueprint (Diversifying Selection)
# ---------------------------------------------------------------------------
message("Processing pN/pS Parallel Events (Single Best Annotation)...")

# A. Unroll Colony & Genera
# Convert list of genera to sorted string
expanded_pnps <- parallel_pnps_D[, .(
  Colony = unlist(Colonies), 
  Genera_Str = sapply(Genus_list, function(x) paste(sort(unique(x)), collapse = ";")) 
), by = corresponding_gene_call]

# B. Get Annotations
ann_map <- unique(gene_dt[!is.na(pfam_str) & pfam_str != "" & !is.na(GO_str) & GO_str != "", 
                          .(corresponding_gene_call, pfam_str, GO_str)])
annotated_pnps <- merge(expanded_pnps, ann_map, by = "corresponding_gene_call")

# C. Select Single Best Annotation
heatmap_dt_pnps <- annotated_pnps[, {
  p_list <- unlist(strsplit(pfam_str, ";;"))
  g_list <- unlist(strsplit(GO_str, ";;"))
  
  p_list <- p_list[!is.na(p_list) & p_list != ""]
  g_list <- g_list[!is.na(g_list) & g_list != ""]
  
  p_pick <- if(length(p_list) > 0) p_list[1] else "No_Pfam"
  g_pick <- if(length(g_list) > 0) g_list[1] else "Uncharacterized"
  
  list(Pfam_Single = p_pick, GO_Single = g_pick)
}, by = .(corresponding_gene_call, Colony, Genera_Str)]

heatmap_dt_pnps[, Y_Label := paste0("(", Pfam_Single, ", ", str_remove(GO_Single, "^GO:\\d+~?"), ")")]

# D. Aggregate per Cell
cell_stats_pnps <- heatmap_dt_pnps[, .(
  n_genes = uniqueN(corresponding_gene_call),
  unique_genera = list(unique(unlist(strsplit(Genera_Str, ";"))))
), by = .(Colony, Y_Label)]

cell_stats_pnps[, n_genera_count := lengths(unique_genera)]
cell_stats_pnps[, Guide_Label := ifelse(n_genera_count > 1, 
                                        paste0(n_genera_count, " gen"), 
                                        substring(sapply(unique_genera, `[`, 1), 1, 1))]

# E. Filter (>= 2 Colonies)
label_counts_pnps <- cell_stats_pnps[, .(n_cols = uniqueN(Colony)), by = Y_Label]
valid_labels_pnps <- label_counts_pnps[n_cols >= 2]$Y_Label
plot_dt_pnps <- cell_stats_pnps[Y_Label %in% valid_labels_pnps]

# F. Plot
message("Generating pN/pS Blueprint Heatmap...")

# Order Y-axis
if(length(unique(plot_dt_pnps$Y_Label)) > 1) {
  mat <- dcast(plot_dt_pnps, Y_Label ~ Colony, value.var = "n_genes", fill = 0)
  mat_num <- as.matrix(mat[, -1])
  rownames(mat_num) <- mat$Y_Label
  dd.row <- as.dendrogram(hclust(dist(mat_num)))
  plot_dt_pnps[, Y_Label := factor(Y_Label, levels = rownames(mat_num)[order.dendrogram(dd.row)])]
}

p_blueprint_pnps <- ggplot(plot_dt_pnps, aes(x = as.factor(Colony), y = Y_Label)) +
  geom_tile(aes(fill = n_genes), color = "white", size = 0.2) +
  geom_tile(data = plot_dt_pnps[n_genera_count > 1], aes(color = "Multi-Genus"), fill = NA, size = 1) +
  geom_text(aes(label = Guide_Label), size = 2.5, color = "grey20", fontface = "bold") +
  
  # Using 'magma' or 'inferno' to distinguish from the 0D 'mako' plot
  scale_fill_viridis_c(option = "inferno", name = "pN/pS > 1 Genes\n(Count)", direction = 1) +
  scale_color_manual(values = c("Multi-Genus" = "black"), labels = "Cross-Genera Event") +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(title = "Diversifying Selection Blueprint (pN/pS > 1)",
       subtitle = "Rows = Functions under parallel selection in >1 Colony.",
       x = "Colony", y = NULL)
p_blueprint_pnps
ggsave(file.path(OUTDIR, "parallelism_blueprint_pnps_simplified.pdf"), p_blueprint_pnps, width = 12, height = 16)


# ---------------------------------------------------------------------------
# PART 2: The Overlap Analysis
# ---------------------------------------------------------------------------
message("Analyzing Functional Overlap between 0D Sweeps and pN/pS...")

# We need the 'Pfam_Single' or 'Y_Label' from the *filtered* datasets 
# to see what made it to the final plots.

# 1. Extract sets from pN/pS (From PART 1 above)
pnps_pfams <- unique(str_extract(as.character(plot_dt_pnps$Y_Label), "\\(.*?,")) # Extract "(Pfam," part
pnps_labels <- unique(as.character(plot_dt_pnps$Y_Label))

# 2. Extract sets from 0D (From previous user block - assuming plot_dt_0D or similar exists)
# If the previous block named it 'plot_dt', we use that. 
# Ideally, re-run the extraction logic for 0D here if variables are lost, 
# but assuming 'plot_dt' from the previous code block is the 0D data:
if(exists("plot_dt")) {
  # Rename for safety if it came from the 0D block
  plot_dt_0D <- plot_dt 
  
  zeroD_pfams <- unique(str_extract(as.character(plot_dt_0D$Y_Label), "\\(.*?,"))
  zeroD_labels <- unique(as.character(plot_dt_0D$Y_Label))
  
  # 3. Find Overlaps
  overlap_labels <- intersect(pnps_labels, zeroD_labels)
  
  # Clean up strings for display
  overlap_pfams_clean <- str_remove_all(intersect(pnps_pfams, zeroD_pfams), "[\\(, ]")
  
  message(paste0("Total Unique Functions in pN/pS Plot: ", length(pnps_labels)))
  message(paste0("Total Unique Functions in 0D Plot:   ", length(zeroD_labels)))
  message(paste0("Functions found in BOTH:             ", length(overlap_labels)))
  
  # 4. Save/Display the Overlap List
  if(length(overlap_labels) > 0) {
    overlap_dt <- data.table(
      Overlapping_Function = overlap_labels,
      Type = "Present in both pN/pS and 0D Parallel sets"
    )
    
    print(overlap_dt)
    fwrite(overlap_dt, file.path(OUTDIR, "overlapping_functions_pnps_vs_0D.csv"))
    
    # Optional: Venn Diagram Logic
    # Just printing the count is usually enough, but here is the list:
    cat("\nTop overlapping functions:\n")
    print(head(overlap_labels, 10))
  } else {
    message("No exact overlap found between the filtered sets.")
  }
} else {
  warning("The 0D plot data (plot_dt) was not found in environment. Run the 0D block first!")
}







########## PLOTTING #############
# --- 0. Setup --------------------------------------------------------------
OUTDIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"

# --- 1. Process 0D Sweeps & Annotate ---------------------------------------
message("Processing 0D Sweeps...")

# A. Unroll Colony, Genus, and MAGs
# FIX: Convert MAG list to a String ("mmag_1;mmag_2") to allow grouping later
expanded_0D <- parallel_zeroD_D[, .(
  Colony = unlist(Colonies), 
  Genera_Str = sapply(Genus_list, function(x) paste(sort(unique(x)), collapse = ";")),
  MAG_Str = sapply(MAGs, function(x) paste(sort(unique(x)), collapse = ";")) # <--- THE FIX
), by = corresponding_gene_call]

# B. Annotate (First Valid Logic)
ann_map <- unique(gene_dt[!is.na(pfam_str) & pfam_str != "" & !is.na(GO_str) & GO_str != "", 
                          .(corresponding_gene_call, pfam_str, GO_str)])
annotated_0D <- merge(expanded_0D, ann_map, by = "corresponding_gene_call")

# Group by MAG_Str (Character) instead of MAG_list (List)
heatmap_dt_0D <- annotated_0D[, {
  p_list <- unlist(strsplit(pfam_str, ";;")); g_list <- unlist(strsplit(GO_str, ";;"))
  p_list <- p_list[!is.na(p_list) & p_list != ""]; g_list <- g_list[!is.na(g_list) & g_list != ""]
  
  p_pick <- if(length(p_list) > 0) p_list[1] else "No_Pfam"
  g_pick <- if(length(g_list) > 0) g_list[1] else "Uncharacterized"
  
  list(Pfam_Single = p_pick, GO_Single = g_pick)
}, by = .(corresponding_gene_call, Colony, Genera_Str, MAG_Str)] # Grouping works now

heatmap_dt_0D[, Y_Label := paste0("(", Pfam_Single, ", ", str_remove(GO_Single, "^GO:\\d+~?"), ")")]

# --- 2. Aggregate per Cell (Colony x Function) -----------------------------
cell_stats <- heatmap_dt_0D[, .(
  n_genes = uniqueN(corresponding_gene_call),
  
  # Local Genera Count (For "2 gen" Label)
  # Split the string back into a vector to count
  cell_genera = list(unique(unlist(strsplit(Genera_Str, ";")))),
  
  # Collect ALL MAGs involved in this cell (For Global Check later)
  # Split the MAG string back into a vector
  cell_mags = list(unique(unlist(strsplit(MAG_Str, ";"))))
), by = .(Colony, Y_Label)]

cell_stats[, n_genera_local := lengths(cell_genera)]

# --- 3. Define Logic -------------------------------------------------------

# A. Cell Label: "2 gen" if multiple genera involved LOCALLY
cell_stats[, Text_Label := ifelse(n_genera_local >= 2, "2 gen", "")]

# B. Row Logic: Global Recurrence (The Black Border)
# Group by Y_Label (Function) and check the "Different Colony" rule
global_stats <- cell_stats[, .(
  n_unique_colonies = uniqueN(Colony),
  # Check if total unique MAGs involved globally > 1
  n_unique_mags_global = uniqueN(unlist(cell_mags)) 
), by = Y_Label]

# Join back to cell_stats
cell_stats <- merge(cell_stats, global_stats, by = "Y_Label")

# Apply Border Logic: 
# Highlight IF: (Occurs in >= 2 Colonies) AND (Involves >= 2 MAGs globally)
cell_stats[, Border_Highlight := ifelse(n_unique_colonies >= 2 & n_unique_mags_global >= 2, "Yes", "No")]

# --- 4. Filtering & pN/pS Overlap ------------------------------------------

# Filter: Show rows that have EITHER recurrence (Border) OR local complexity (2 gen)
plot_dt <- cell_stats[Border_Highlight == "Yes" | Text_Label == "2 gen"]

# pN/pS Overlap Check
expanded_pnps <- parallel_pnps_D[, .(Colony = unlist(Colonies)), by = corresponding_gene_call]
annotated_pnps <- merge(expanded_pnps, ann_map, by = "corresponding_gene_call")
heatmap_dt_pnps <- annotated_pnps[, {
  p_list <- unlist(strsplit(pfam_str, ";;")); g_list <- unlist(strsplit(GO_str, ";;"))
  p_list <- p_list[!is.na(p_list) & p_list != ""]; g_list <- g_list[!is.na(g_list) & g_list != ""]
  p_pick <- if(length(p_list) > 0) p_list[1] else "No_Pfam"
  g_pick <- if(length(g_list) > 0) g_list[1] else "Uncharacterized"
  list(Pfam_Single = p_pick, GO_Single = g_pick)
}, by = .(corresponding_gene_call, Colony)]

heatmap_dt_pnps[, Y_Label := paste0("(", Pfam_Single, ", ", str_remove(GO_Single, "^GO:\\d+~?"), ")")]
pnps_overlap_labels <- unique(heatmap_dt_pnps$Y_Label)

# Mark Overlap
plot_dt[, Is_Overlap := Y_Label %in% pnps_overlap_labels]

# --- 5. Plotting -----------------------------------------------------------
message("Generating Recurrence Map...")

# Cluster Y-axis
if(length(unique(plot_dt$Y_Label)) > 1) {
  mat <- dcast(plot_dt, Y_Label ~ Colony, value.var = "n_genes", fill = 0)
  mat_num <- as.matrix(mat[, -1])
  rownames(mat_num) <- mat$Y_Label
  dd.row <- as.dendrogram(hclust(dist(mat_num)))
  ordered_labels <- rownames(mat_num)[order.dendrogram(dd.row)]
  plot_dt[, Y_Label := factor(Y_Label, levels = ordered_labels)]
}

# Axis Colors
axis_colors <- ifelse(levels(plot_dt$Y_Label) %in% pnps_overlap_labels, "#D55E00", "black")
axis_faces  <- ifelse(levels(plot_dt$Y_Label) %in% pnps_overlap_labels, "bold", "plain")

p_recurrence <- ggplot(plot_dt, aes(x = as.factor(Colony), y = Y_Label)) +
  
  # 1. BASE: Intensity
  geom_tile(aes(fill = n_genes), color = "white", size = 0.2) +
  
  # 2. BORDER: Colony Recurrence (Global Signal)
  geom_tile(data = plot_dt[Border_Highlight == "Yes"], 
            aes(color = "Recurrent"), fill = NA, size = 1) +
  
  # 3. TEXT: Local Diversity (Local Signal)
  geom_text(aes(label = Text_Label), size = 3, color = "white", fontface = "bold") +
  
  # Scales
  scale_fill_viridis_c(option = "mako", name = "0D Genes", direction = -1) +
  scale_color_manual(values = c("Recurrent" = "black"), labels = "Found in ≥2 Colonies") +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8, 
                               color = axis_colors,
                               face = axis_faces),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Consolidated Functional Sweep Blueprint",
       subtitle = "Black Border = Function sweeps in ≥2 unique Colonies (Recurrence).\n'2 gen' = Sweep involves ≥2 Genera in that Colony.\nOrange Axis = Function also under Diversifying Selection.",
       x = "Colony", y = NULL)

p_recurrence





# --- 0. Setup --------------------------------------------------------------
OUTDIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"

# Ensure Gene Length exists (Use 4D sites as proxy if full length missing)
if(!"gene_length" %in% names(gene_dt)) {
  if("gene_length_4D_sites" %in% names(gene_dt)) {
    gene_dt[, gene_length := gene_length_4D_sites * 4]
  } else {
    stop("Column 'gene_length' is missing from gene_dt.")
  }
}

# 1. Define "Hits"
# A gene is a hit if it has a 0D sweep OR pN/pS > 1
hits_dt <- gene_dt[, .(
  MAG_id, 
  corresponding_gene_call, 
  gene_length,
  Is_Hit = ( (!is.na(pN_pS_ratio) & pN_pS_ratio > 1) | 
               (corresponding_gene_call %in% parallel_zeroD_D$corresponding_gene_call) ) 
)]

hits_dt[is.na(Is_Hit), Is_Hit := FALSE]
hits_dt <- unique(hits_dt) # Ensure unique gene entries

# --- 1. The Statistical Test Function --------------------------------------

test_likelihoods <- function(mag_A, mag_B, n_sims = 2000) {
  
  # A. Define the Shared Universe (Orthologs)
  genes_A <- hits_dt[MAG_id == mag_A]
  genes_B <- hits_dt[MAG_id == mag_B]
  
  common_genes <- intersect(genes_A$corresponding_gene_call, genes_B$corresponding_gene_call)
  
  # If fewer than 10 shared genes, we can't efficiently test stats
  if(length(common_genes) < 10) return(NULL)
  
  # Filter to this shared universe
  set_A <- genes_A[corresponding_gene_call %in% common_genes]
  set_B <- genes_B[corresponding_gene_call %in% common_genes]
  
  # B. Calculate Independent Likelihoods (Question 1 & 2)
  # "What is the likelihood I observe a sweep in MAG A?"
  # (Count of Hits / Total Shared Genes)
  n_hits_A <- sum(set_A$Is_Hit)
  prob_hit_A <- n_hits_A / nrow(set_A)
  
  n_hits_B <- sum(set_B$Is_Hit)
  prob_hit_B <- n_hits_B / nrow(set_B)
  
  # C. Calculate Observed Joint Probability (Overlap)
  hits_A_names <- set_A[Is_Hit == TRUE]$corresponding_gene_call
  hits_B_names <- set_B[Is_Hit == TRUE]$corresponding_gene_call
  
  obs_overlap <- length(intersect(hits_A_names, hits_B_names))
  
  # D. Simulation (Question 3: Conditional Likelihood)
  # We fix n_hits_A and n_hits_B. We randomly shuffle WHO gets hit based on LENGTH.
  
  # Weights
  w_A <- set_A$gene_length
  w_B <- set_B$gene_length
  
  # Safety check for weights
  if(any(is.na(w_A)) || sum(w_A) == 0) w_A <- rep(1, length(w_A))
  if(any(is.na(w_B)) || sum(w_B) == 0) w_B <- rep(1, length(w_B))
  
  sim_overlaps <- replicate(n_sims, {
    # Randomly pick 'n_hits' genes, weighted by length
    sim_A <- sample(set_A$corresponding_gene_call, n_hits_A, prob = w_A)
    sim_B <- sample(set_B$corresponding_gene_call, n_hits_B, prob = w_B)
    length(intersect(sim_A, sim_B))
  })
  
  # E. Stats
  exp_overlap <- mean(sim_overlaps)
  z_score <- (obs_overlap - exp_overlap) / sd(sim_overlaps)
  
  # P-value: Probability of observing this much overlap (or more) by chance
  p_val <- (sum(sim_overlaps >= obs_overlap) + 1) / (n_sims + 1)
  
  return(data.table(
    MAG_A = mag_A,
    MAG_B = mag_B,
    Total_Shared_Genes = length(common_genes),
    
    # Q1: Likelihood in A
    Hits_A = n_hits_A,
    Prob_Hit_A = round(prob_hit_A, 4),
    
    # Q2: Likelihood in B
    Hits_B = n_hits_B,
    Prob_Hit_B = round(prob_hit_B, 4),
    
    # Q3: Conditional Joint Likelihood
    Observed_Overlap = obs_overlap,
    Expected_Overlap = round(exp_overlap, 2),
    Enrichment = round(obs_overlap / exp_overlap, 2),
    P_Value = p_val
  ))
}

# --- 2. Run on Identified Pairs --------------------------------------------
message("Identifying pairs to test...")

# Extract pairs from parallel_zeroD_D where n_MAG >= 2
events_to_test <- rbind(parallel_pnps_D, parallel_zeroD_D, fill=TRUE)
events_to_test <- events_to_test[n_MAG >= 2]

# Build Pair List
pair_list <- list()
for(i in 1:nrow(events_to_test)) {
  mags <- unlist(events_to_test$MAGs[i])
  # Clean names
  mags <- mags[!is.na(mags) & mags != ""]
  if(length(mags) >= 2) {
    mags <- sort(unique(mags))
    combos <- t(combn(mags, 2))
    for(r in 1:nrow(combos)) {
      pair_id <- paste(combos[r,1], combos[r,2], sep="_")
      if(is.null(pair_list[[pair_id]])) {
        pair_list[[pair_id]] <- data.table(MAG_A = combos[r,1], MAG_B = combos[r,2])
      }
    }
  }
}
pairs_dt <- rbindlist(pair_list)

message(paste("Running likelihood tests on", nrow(pairs_dt), "pairs..."))

# Serial execution with error handling
results_list <- lapply(1:nrow(pairs_dt), function(i) {
  tryCatch({
    test_likelihoods(pairs_dt$MAG_A[i], pairs_dt$MAG_B[i])
  }, error = function(e) { return(NULL) })
})

final_stats <- rbindlist(results_list)
final_stats[, FDR := p.adjust(P_Value, method = "BH")]

# Add Taxonomy Labels
tax_map <- unique(tax_dt[, .(MAG_id, Genus)])
final_stats <- merge(final_stats, tax_map, by.x = "MAG_A", by.y = "MAG_id"); setnames(final_stats, "Genus", "Genus_A")
final_stats <- merge(final_stats, tax_map, by.x = "MAG_B", by.y = "MAG_id"); setnames(final_stats, "Genus", "Genus_B")




# --- 0. Setup --------------------------------------------------------------
OUTDIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"
# Load the stats table generated in the previous step
final_stats <- fread(file.path(OUTDIR, "parallelism_likelihood_stats.csv"))

# Create a Label for plotting (Genus A vs Genus B)
final_stats[, Label := paste0(Genus_A, " vs ", Genus_B)]

# Classify the Interaction
final_stats[, Type := ifelse(Genus_A == Genus_B, "Within-Genus", "Cross-Genus")]

# Filter for plotting (Only show pairs with at least 1 overlap to keep plot clean)
plot_dt <- final_stats[Observed_Overlap > 0]

# --- PLOT 1: Observed vs. Expected Scatter ---------------------------------
# This is the "Proof of Signal" plot

p_oe <- ggplot(plot_dt, aes(x = Expected_Overlap, y = Observed_Overlap)) +
  
  # 1. The Null Line (Random Chance)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  
  # 2. The Points (MAG Pairs)
  geom_point(aes(color = Type, size = -log10(FDR)), alpha = 0.7) +
  

  
  # Scales
  scale_color_manual(values = c("Cross-Genus" = "#D55E00", "Within-Genus" = "#0072B2")) +
  scale_size_continuous(name = "Significance\n(-log10 FDR)", range = c(2, 6)) +
  
  theme_bw() +
  labs(title = "Excess of Parallel Evolution",
       subtitle = "Points above the dashed line indicate convergent evolution exceeding random chance.",
       x = "Expected Shared Genes (Null Model)",
       y = "Observed Shared Genes (Real Data)")

p_oe
# --- PLOT 2: The Significance Matrix ---------------------------------------
# This shows WHO is converging with WHO.

# We need to aggregate because we might have multiple MAG pairs for the same Genus pair
genus_matrix <- final_stats[, .(
  mean_enrichment = mean(Enrichment),
  max_significance = max(-log10(FDR)),
  n_pairs = .N
), by = .(Genus_A, Genus_B)]

p_matrix <- ggplot(genus_matrix, aes(x = Genus_A, y = Genus_B)) +
  
  # Tile Color = Significance
  geom_tile(aes(fill = max_significance), color = "white") +
  
  # Tile Label = Enrichment Score
  geom_text(aes(label = round(mean_enrichment, 1)), size = 3, color = "white") +
  
  scale_fill_viridis_c(option = "magma", name = "Max Significance\n(-log10 FDR)", direction = 1) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "Genus-Level Convergence Map",
       subtitle = "Color = Statistical Significance | Number = Fold Enrichment (Obs/Exp)",
       x = NULL, y = NULL)
p_matrix
ggsave(file.path(OUTDIR, "stats_observed_vs_expected.pdf"), p_oe, width = 7, height = 6)


# --- PLOT 2: The Significance Matrix ---------------------------------------
# This shows WHO is converging with WHO.

# We need to aggregate because we might have multiple MAG pairs for the same Genus pair
genus_matrix <- final_stats[, .(
  mean_enrichment = mean(Enrichment),
  max_significance = max(-log10(FDR)),
  n_pairs = .N
), by = .(Genus_A, Genus_B)]

p_matrix <- ggplot(genus_matrix, aes(x = Genus_A, y = Genus_B)) +
  
  # Tile Color = Significance
  geom_tile(aes(fill = max_significance), color = "white") +
  
  # Tile Label = Enrichment Score
  geom_text(aes(label = round(mean_enrichment, 1)), size = 3, color = "white") +
  
  scale_fill_viridis_c(option = "magma", name = "Max Significance\n(-log10 FDR)", direction = 1) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "Genus-Level Convergence Map",
       subtitle = "Color = Statistical Significance | Number = Fold Enrichment (Obs/Exp)",
       x = NULL, y = NULL)

ggsave(file.path(OUTDIR, "stats_genus_matrix.pdf"), p_matrix, width = 8, height = 7)



# --- 0. Setup --------------------------------------------------------------
OUTDIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"

# Load the stats from the previous step
final_stats <- fread(file.path(OUTDIR, "parallelism_likelihood_stats.csv"))

# 1. Select the "Top Pairs" to visualize
# We don't want to plot 100 histograms. Let's pick the top 12 significant pairs.
# Priority: FDR < 0.05, then sort by Enrichment magnitude
top_pairs <- final_stats[FDR < 0.05][order(-Enrichment)][1:min(.N, 12)]

if(nrow(top_pairs) == 0) stop("No significant pairs found (FDR < 0.05). cannot generate plots.")

message(paste("Generating Null Histograms for top", nrow(top_pairs), "pairs..."))

# --- 1. Generate Distribution Data -----------------------------------------
# We need to re-run the sim for these specific pairs to get the raw counts for the histogram

# (Ensure hits_dt is loaded from previous steps!)
# If not, strictly define it again:
if(!exists("hits_dt")) stop("Please run the 'Setup' block from the previous turn to define 'hits_dt'")

get_sim_distribution <- function(mag_A, mag_B, obs_val, n_sims = 2000) {
  
  # A. Setup Universe
  genes_A <- hits_dt[MAG_id == mag_A]
  genes_B <- hits_dt[MAG_id == mag_B]
  common <- intersect(genes_A$corresponding_gene_call, genes_B$corresponding_gene_call)
  
  set_A <- genes_A[corresponding_gene_call %in% common]
  set_B <- genes_B[corresponding_gene_call %in% common]
  
  # B. Simulation Parameters
  n_hits_A <- sum(set_A$Is_Hit)
  n_hits_B <- sum(set_B$Is_Hit)
  w_A <- set_A$gene_length; if(sum(w_A)==0) w_A <- rep(1, length(w_A))
  w_B <- set_B$gene_length; if(sum(w_B)==0) w_B <- rep(1, length(w_B))
  
  # C. Run Sim
  sim_counts <- replicate(n_sims, {
    sim_A <- sample(set_A$corresponding_gene_call, n_hits_A, prob = w_A)
    sim_B <- sample(set_B$corresponding_gene_call, n_hits_B, prob = w_B)
    length(intersect(sim_A, sim_B))
  })
  
  return(data.table(
    Pair_Label = paste0(unique(genes_A$Genus), " vs\n", unique(genes_B$Genus)),
    Simulated_Overlap = sim_counts,
    Observed_Overlap = obs_val
  ))
}

plot_data_list <- list()

for(i in 1:nrow(top_pairs)) {
  # Add TryCatch to prevent crash on single bad pair
  try({
    res <- get_sim_distribution(top_pairs$MAG_A[i], top_pairs$MAG_B[i], top_pairs$Observed_Overlap[i])
    plot_data_list[[i]] <- res
  })
}

plot_df <- rbindlist(plot_data_list)

# --- 2. Plotting -----------------------------------------------------------

p_hist <- ggplot(plot_df, aes(x = Simulated_Overlap)) +
  
  # A. The Null Distribution (Grey Bars)
  geom_histogram(binwidth = 1, fill = "grey70", color = "grey90") +
  
  # B. The Observed Value (Red Dashed Line)
  geom_vline(aes(xintercept = Observed_Overlap), 
             color = "#D55E00", linetype = "dashed", size = 1) +
  
  # C. Faceting
  facet_wrap(~Pair_Label, scales = "free") +
  
  # D. Theme & Labels
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(size = 8)
  ) +
  labs(title = "Significance of Parallel Evolution",
       subtitle = "Grey = Null Distribution (Random Chance) | Red Line = Observed Parallel Genes",
       x = "Number of Parallel Genes (Simulated)",
       y = "Frequency (Count)")

ggsave(file.path(OUTDIR, "parallelism_null_histograms.pdf"), p_hist, width = 10, height = 8)

# Print to screen
print(p_hist)





# ==========================================
# 1. AGGREGATION: GENE -> MAG (THE WEIGHTED WAY)
# ==========================================

# We assume 'gene_dt' has the columns from your Python output:
# dN, dS, nN_gene_reference, nS_gene_reference

df_mag <- gene_dt %>%
  filter(!is.na(pS)) %>% # Remove failed calculations
  group_by(MAG_id, colony_id, Month) %>% 
  summarise(
    # 1. Calculate MAG-wide pS (Sum of dS / Sum of Sites)
    total_dS = sum(dS, na.rm = TRUE),
    total_S_sites = sum(nS_gene_reference, na.rm = TRUE),
    pS = total_dS / total_S_sites,
    
    # 2. Calculate MAG-wide pN (Sum of dN / Sum of Sites)
    total_dN = sum(dN, na.rm = TRUE),
    total_N_sites = sum(nN_gene_reference, na.rm = TRUE),
    pN = total_dN / total_N_sites,
    
    # 3. Gene Count for record keeping
    n_genes = n() 
  ) %>%
  # Now calculate the robust ratio from the sums
  mutate(pN_pS_ratio = pN / pS) %>%
  ungroup()

# ==========================================
# 2. FILTERING
# ==========================================

# We filter out pS = 0 here. 
# This removes the "infinite ratio" problems and the "discrete" noise at 0.
MAX_PS_THRESHOLD <- 0.1 

df_analysis <- df_mag %>%
  filter(is.finite(pN_pS_ratio)) %>% # Remove the "99" or Infs if they snuck in
  filter(pS > 0) %>%
  filter(pS <= MAX_PS_THRESHOLD)

# ==========================================
# 3. LOG-BINNING THE DATA (Weighted by Site Count)
# ==========================================

min_val <- 10e-5
max_val <- max(df_analysis$pS)
breaks <- 10^seq(log10(min_val), log10(max_val), length.out = 25)

df_binned <- df_analysis %>%
  mutate(bin = cut(pS, breaks = breaks, include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    # We weight the bin average by the amount of data (S_sites) in that MAG
    # This prevents a tiny, noisy MAG from skewing the bin
    mean_pS = weighted.mean(pS, w = total_S_sites, na.rm = TRUE),
    mean_pN_pS = weighted.mean(pN_pS_ratio, w = total_S_sites, na.rm = TRUE),
    n_mags = n()
  ) %>%
  filter(!is.na(mean_pS))

# ==========================================
# 4. MODEL FITTING (Garud/Rocha Model)
# ==========================================

purifying_model <- function(x, f_del, K) {
  (1 - f_del) + f_del * ((1 - exp(-K * x)) / (K * x))
}

fit <- nlsLM(mean_pN_pS ~ purifying_model(mean_pS, f_del, K),
             data = df_binned,
             start = list(f_del = 0.85, K = 1000), 
             control = nls.lm.control(maxiter = 100))

params <- coef(fit)
f_est <- params["f_del"]
K_est <- params["K"]

# ==========================================
# 5. VISUALIZATION
# ==========================================

x_pred <- 10^seq(log10(min_val), log10(max_val), length.out = 200)
y_pred <- purifying_model(x_pred, f_est, K_est)
df_pred <- data.frame(pS = x_pred, pN_pS = y_pred)

ggplot() +
  # Raw MAG Data
  geom_point(data = df_analysis, aes(x = pS, y = pN_pS_ratio), 
             alpha = 0.2, color = "grey50", size = 1) +
  # Binned Averages
  geom_point(data = df_binned, aes(x = mean_pS, y = mean_pN_pS), 
             color = "black", size = 3, shape = 21, fill = "white", stroke = 1.2) +
  # Model Curve
  geom_line(data = df_pred, aes(x = pS, y = pN_pS), 
            color = "#D55E00", size = 1.5) +
  scale_x_continuous(trans="log10")+
  scale_y_continuous(trans="log10")+
  annotation_logticks() +
  labs(
    title = "Decay of Selection Efficacy (MAG Level)",
    subtitle = paste0("Model Fit: K = ", round(K_est, 0), " (Efficiency), f_del = ", round(f_est, 3)),
    x = "Synonymous Diversity (pS)",
    y = "pN/pS Ratio"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())







# 1. Define Hits (The Observed Reality)
# Ensure gene_length is present (use 4D proxy if needed)
if(!"gene_length" %in% names(gene_dt)) gene_dt[, gene_length := gene_length_4D_sites * 4]

# A gene is a hit if 0D sweep OR pN/pS > 1
hits_dt <- gene_dt[, .(
  MAG_id, 
  corresponding_gene_call, 
  gene_length,
  Is_Hit = ( (!is.na(pN_pS_ratio) & pN_pS_ratio > 1) | 
               (corresponding_gene_call %in% parallel_zeroD_D$corresponding_gene_call) ) 
)]
hits_dt[is.na(Is_Hit), Is_Hit := FALSE]
hits_dt <- unique(hits_dt)

# --- 1. Define the Pairs to Test -------------------------------------------
# We want to test the pairs that exhibited parallelism in your Heatmap/Blueprint.
# Let's extract them from the parallel dataset.
events <- rbind(parallel_pnps_D, parallel_zeroD_D, fill=TRUE)
events <- events[n_MAG >= 2]

# Extract unique pairs found in the data
pair_list <- list()
for(i in 1:min(nrow(events), 500)) { # Check top 500 events to build pair list
  mags <- unlist(events$MAGs[i])
  mags <- mags[!is.na(mags) & mags != ""]
  if(length(mags) >= 2) {
    combos <- t(combn(sort(unique(mags)), 2))
    for(r in 1:nrow(combos)) {
      pid <- paste(combos[r,1], combos[r,2], sep="_")
      if(is.null(pair_list[[pid]])) pair_list[[pid]] <- data.table(MAG_A=combos[r,1], MAG_B=combos[r,2])
    }
  }
}
pairs_to_test <- rbindlist(pair_list)
pairs_to_test <- unique(pairs_to_test)

# Optional: To keep the plot readable, let's just pick the top 6 distinct pairs 
# (e.g. involving different genera) or just run all and filter later.
message(paste("Identified", nrow(pairs_to_test), "pairs to simulate..."))

# --- 2. The Simulation Function (The Null Model) ---------------------------
run_permutation_test <- function(mag_A, mag_B, n_sims = 1000) {
  
  # A. Define Universe (Shared Orthologs)
  genes_A <- hits_dt[MAG_id == mag_A]
  genes_B <- hits_dt[MAG_id == mag_B]
  common <- intersect(genes_A$corresponding_gene_call, genes_B$corresponding_gene_call)
  
  if(length(common) < 10) return(NULL)
  
  set_A <- genes_A[corresponding_gene_call %in% common]
  set_B <- genes_B[corresponding_gene_call %in% common]
  
  # B. Get Observed Data
  obs_hits_A <- set_A[Is_Hit == TRUE]$corresponding_gene_call
  obs_hits_B <- set_B[Is_Hit == TRUE]$corresponding_gene_call
  obs_overlap <- length(intersect(obs_hits_A, obs_hits_B))
  
  # If 0 overlap observed, not interesting for this plot
  if(obs_overlap == 0) return(NULL)
  
  # C. Simulation Loop (Shuffling Mutations)
  # We preserve the Number of Hits, but shuffle location weighted by Length
  n_darts_A <- length(obs_hits_A)
  n_darts_B <- length(obs_hits_B)
  w_A <- set_A$gene_length; if(sum(w_A)==0) w_A <- rep(1, length(w_A))
  w_B <- set_B$gene_length; if(sum(w_B)==0) w_B <- rep(1, length(w_B))
  
  sim_results <- replicate(n_sims, {
    # Randomly assign hits
    sim_A <- sample(set_A$corresponding_gene_call, n_darts_A, prob = w_A)
    sim_B <- sample(set_B$corresponding_gene_call, n_darts_B, prob = w_B)
    length(intersect(sim_A, sim_B))
  })
  
  # Calculate P-value
  p_val <- (sum(sim_results >= obs_overlap) + 1) / (n_sims + 1)
  
  # Return data structure for plotting
  return(data.table(
    MAG_A = mag_A, MAG_B = mag_B,
    Simulated_Overlap = sim_results,
    Observed_Overlap = obs_overlap,
    P_Val = p_val
  ))
}

# --- 3. Run Simulations ----------------------------------------------------
message("Running simulations...")
sim_list <- lapply(1:nrow(pairs_to_test), function(i) {
  run_permutation_test(pairs_to_test$MAG_A[i], pairs_to_test$MAG_B[i])
})
plot_data <- rbindlist(sim_list)

# Add Genus Names for nice labels
tax_map <- unique(tax_dt[, .(MAG_id, Genus)])
plot_data <- merge(plot_data, tax_map, by.x="MAG_A", by.y="MAG_id"); setnames(plot_data, "Genus", "Genus_A")
plot_data <- merge(plot_data, tax_map, by.x="MAG_B", by.y="MAG_id"); setnames(plot_data, "Genus", "Genus_B")
plot_data[, Pair_Label := paste0(Genus_A, " x ", Genus_B, "\n(", MAG_A, " vs ", MAG_B, ")")]

# Filter: Only show top 9 most significant/enriched pairs for the figure
# (Sort by P-value then Observed Overlap)
top_pairs <- unique(plot_data[, .(Pair_Label, P_Val, Observed_Overlap)])[order(P_Val, -Observed_Overlap)][1:9]
plot_subset <- plot_data[Pair_Label %in% top_pairs$Pair_Label]

# --- 4. Plotting The Null Histogram ----------------------------------------
message("Generating Null Distribution Histogram...")

p_null_test <- ggplot(plot_subset, aes(x = Simulated_Overlap)) +
  
  # A. The Null Distribution (Grey Histogram)
  # This represents the "Neutral Model" expectation
  geom_histogram(binwidth = 1, fill = "grey70", color = "grey50", alpha = 0.8) +
  
  # B. The Observed Value (Red Line)
  geom_vline(aes(xintercept = Observed_Overlap), 
             color = "#D55E00", linetype = "dashed", size = 1) +
  
  # C. P-value Annotation
  geom_text(data = unique(plot_subset[, .(Pair_Label, P_Val, Observed_Overlap)]),
            aes(x = Inf, y = Inf, label = paste0("p < ", round(P_Val, 4))),
            hjust = 1.1, vjust = 1.5, size = 3, fontface = "bold", color = "black") +
  
  # Facet
  facet_wrap(~Pair_Label, scales = "free") +
  
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Significance of Parallel Evolution (Null Model)",
       subtitle = "Grey = Expected overlap by chance (1000 simulations) | Red = Observed Overlap",
       x = "Number of Parallel Genes (Simulated)",
       y = "Frequency")

ggsave(file.path(OUTDIR, "parallelism_null_model_histogram.pdf"), p_null_test, width = 12, height = 10)
print(p_null_test)


head(gene_dt)























# ==========================================
# 1. AGGREGATION: GENE -> MAG (Weighted)
# ==========================================
# We use the raw counts (SNVs and Site lengths) to calculate the 
# true weighted mean for the MAG, rather than averaging the pre-calculated 'pi' columns.
gene_dt
df_mag <- gene_dt %>%
  # Filter out rows where lengths are 0 to avoid division by zero
  filter(effective_length_4D_sites > 0, effective_length_0D_sites > 0) %>%
  
  group_by(MAG_id, colony_id, Month) %>% 
  summarise(
    # 1. Calculate MAG-wide pi_4D (The Clock)
    # Sum of all 4D SNVs / Sum of all 4D Sites
    total_SNVs_4D = sum(SNVs_4D_sites, na.rm = TRUE),
    total_sites_4D = sum(effective_length_4D_sites, na.rm = TRUE),
    mag_pi_4D = total_SNVs_4D / total_sites_4D,
    
    # 2. Calculate MAG-wide pi_1D (The Target)
    # Sum of all 1D SNVs / Sum of all 1D Sites
    total_SNVs_1D = sum(SNVs_0D_sites, na.rm = TRUE),
    total_sites_1D = sum(effective_length_0D_sites, na.rm = TRUE),
    mag_pi_1D = total_SNVs_1D / total_sites_1D,
    
    n_genes = n()
  ) %>%
  # 3. Calculate the Ratio
  mutate(pi_ratio = mag_pi_1D / mag_pi_4D) %>%
  ungroup()

# ==========================================
# 2. CONFIGURATION & FILTERING
# ==========================================

# Set the cutoff for "Recent Evolution"
# 4D sites evolve faster than "pS" (which mixes 2D sites), so this threshold 
# might represent a slightly "younger" absolute time than 0.1 pS.
MAX_PI4D_THRESHOLD <- 0.1 

df_analysis <- df_mag %>%
  filter(is.finite(pi_ratio)) %>%   # Remove Infs
  filter(mag_pi_4D > 0) %>%         # Need > 0 for log plot
  filter(mag_pi_4D <= MAX_PI4D_THRESHOLD)

# ==========================================
# 3. BINNING (Weighted by 4D Sites)
# ==========================================

min_val <- min(df_analysis$mag_pi_4D)
max_val <- max(df_analysis$mag_pi_4D)
breaks <- 10^seq(log10(min_val), log10(max_val), length.out = 25)

df_binned <- df_analysis %>%
  mutate(bin = cut(mag_pi_4D, breaks = breaks, include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    # Weight the bin average by total_sites_4D (information content)
    mean_pi_4D = weighted.mean(mag_pi_4D, w = total_sites_4D, na.rm = TRUE),
    mean_ratio = weighted.mean(pi_ratio, w = total_sites_4D, na.rm = TRUE),
    n_mags = n()
  ) %>%
  filter(!is.na(mean_pi_4D))

# ==========================================
# 4. MODEL FITTING
# ==========================================

# The same equation applies, but the parameters might shift slightly.
# f_del might be higher (1D sites are strictly constrained).
purifying_model <- function(x, f_del, K) {
  (1 - f_del) + f_del * ((1 - exp(-K * x)) / (K * x))
}

fit <- nlsLM(mean_ratio ~ purifying_model(mean_pi_4D, f_del, K),
             data = df_binned,
             start = list(f_del = 0.90, K = 1000), 
             control = nls.lm.control(maxiter = 100))

params <- coef(fit)
f_est <- params["f_del"]
K_est <- params["K"]

# ==========================================
# 5. VISUALIZATION
# ==========================================

x_pred <- 10^seq(log10(min_val), log10(max_val), length.out = 200)
y_pred <- purifying_model(x_pred, f_est, K_est)
df_pred <- data.frame(mag_pi_4D = x_pred, pi_ratio = y_pred)

ggplot() +
  # Raw Data
  geom_point(data = df_analysis, aes(x = mag_pi_4D, y = pi_ratio), 
             alpha = 0.15, color = "grey60", size = 1.5) +
  
  # Binned Means
  geom_point(data = df_binned, aes(x = mean_pi_4D, y = mean_ratio), 
             color = "black", size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  # Fitted Curve
  geom_line(data = df_pred, aes(x = mag_pi_4D, y = pi_ratio), 
            color = "#0072B2", size = 1.5) + # Blue for 1D/4D distinction
  
  scale_x_log10(labels = scales::label_number()) +
  scale_y_log10() +
  annotation_logticks() +

  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())









# ==========================================
# 1. SETUP & GENE SPLITTING
# ==========================================

clean_genes <- gene_dt %>%
  filter(!is.na(dS), !is.na(dN), !is.na(nS_gene_reference), !is.na(nN_gene_reference)) %>%
  filter(nS_gene_reference > 0, nN_gene_reference > 0)

set.seed(123) 

genes_split <- clean_genes %>%
  mutate(split_group = sample(c("A", "B"), n(), replace = TRUE))

# ==========================================
# 2. AGGREGATION (INDEPENDENT AXES)
# ==========================================

# GROUP A: X-AXIS (pS only)
mag_stats_A <- genes_split %>%
  filter(split_group == "A") %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(
    total_dS_A = sum(dS, na.rm=TRUE),
    total_S_sites_A = sum(nS_gene_reference, na.rm=TRUE),
    pS_A = total_dS_A / total_S_sites_A,
    .groups = "drop"
  )

# GROUP B: Y-AXIS (pN/pS ratio)
mag_stats_B <- genes_split %>%
  filter(split_group == "B") %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(
    total_dS_B = sum(dS, na.rm=TRUE),
    total_S_sites_B = sum(nS_gene_reference, na.rm=TRUE),
    total_dN_B = sum(dN, na.rm=TRUE),
    total_N_sites_B = sum(nN_gene_reference, na.rm=TRUE),
    
    pS_B = total_dS_B / total_S_sites_B,
    pN_B = total_dN_B / total_N_sites_B,
    .groups = "drop"
  ) %>%
  mutate(ratio_B = pN_B / pS_B)

# JOIN & FILTER
df_split_analysis <- inner_join(mag_stats_A, mag_stats_B, 
                                by = c("MAG_id", "colony_id", "Month", "Genus")) %>%
  mutate(
    X_axis = pS_A,       
    Y_axis = ratio_B     
  ) %>%
  filter(X_axis > 0, X_axis <= 0.1) %>%
  filter(is.finite(Y_axis), Y_axis >= 0.01) # Filter Y < 10^-2

# ==========================================
# 3. CALCULATE MEDIANS PER GENUS (NEW)
# ==========================================

genus_summary <- df_split_analysis %>%
  group_by(Genus) %>%
  summarise(
    median_X = median(X_axis, na.rm = TRUE),
    median_Y = median(Y_axis, na.rm = TRUE),
    n_count = n()
  )

# ==========================================
# 4. WEIGHTED BINNING & FITTING
# ==========================================

min_val <- 10e-5 
max_val <- max(df_split_analysis$X_axis)
breaks <- 10^seq(log10(min_val), log10(max_val), length.out = 25)

df_binned <- df_split_analysis %>%
  mutate(bin = cut(X_axis, breaks = breaks, include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    mean_X = weighted.mean(X_axis, w = total_S_sites_A, na.rm = TRUE),
    mean_Y = weighted.mean(Y_axis, w = total_S_sites_B, na.rm = TRUE),
    n_mags = n()
  ) %>%
  filter(!is.na(mean_X))

purifying_model <- function(x, f_del, K) {
  (1 - f_del) + f_del * ((1 - exp(-K * x)) / (K * x))
}

fit <- nlsLM(mean_Y ~ purifying_model(mean_X, f_del, K),
             data = df_binned,
             start = list(f_del = 0.85, K = 1000),
             lower = c(f_del = 0, K = 0),
             upper = c(f_del = 1, K = 50000),
             control = nls.lm.control(maxiter = 100))

params <- coef(fit)
f_est <- params["f_del"]
K_est <- params["K"]

# ==========================================
# 5. VISUALIZATION
# ==========================================

x_seq <- 10^seq(log10(min_val), log10(max_val), length.out = 200)
y_seq <- purifying_model(x_seq, f_est, K_est)
df_pred <- data.frame(X_axis = x_seq, Y_axis = y_seq)

p <- ggplot() +
  # A. Raw Data (Faded Background)
  geom_point(data = df_split_analysis, aes(x = X_axis, y = Y_axis, color = Genus), 
             alpha = 0.3, size = 1.5) +
  
  # B. Genus Medians (Big "X" Markers)
  geom_point(data = genus_summary, aes(x = median_X, y = median_Y, color = Genus), 
             shape = 4, size = 6, stroke = 2) + # shape 4 is "X", stroke controls thickness
  
  # C. The Model Curve (On top of everything)
  geom_line(data = df_pred, aes(x = X_axis, y = Y_axis), 
            color = "#D55E00", size = 1.5) +
  
  # Neutral Line
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  
  # Scales
  scale_x_log10(labels = scales::label_number()) +
  scale_y_log10(limits = c(0.01, 10), breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
  annotation_logticks() +
  
  labs(
    title = "Decay of Selection Efficacy (Split-Half Independent)",
    subtitle = paste0("Model Fit: K = ", round(K_est, 0), ", f_del = ", round(f_est, 3)),
    x = "Synonymous Diversity (pS) [Set A]",
    y = "pN/pS Ratio [Set B]"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

# Show plot in R window
print(p)

























# ==========================================
# 1. SETUP
# ==========================================
clean_genes <- gene_dt %>%
  filter(!is.na(dS), !is.na(dN), !is.na(nS_gene_reference), !is.na(nN_gene_reference)) %>%
  filter(nS_gene_reference > 0, nN_gene_reference > 0)

# CONFIG
N_BOOTSTRAPS <- 100
MAX_PS <- 0.1
MIN_VAL <- 10e-5

# Storage
boot_params <- data.frame()
boot_preds <- data.frame()

# Prediction sequence for the ribbon (Smooth x-axis)
x_seq <- 10^seq(log10(MIN_VAL), log10(MAX_PS), length.out = 200)

message(paste("Starting", N_BOOTSTRAPS, "bootstrap iterations..."))

# ==========================================
# 2. THE BOOTSTRAP LOOP
# ==========================================
set.seed(42)

for (i in 1:N_BOOTSTRAPS) {
  if (i %% 10 == 0) message(paste("Iteration:", i))
  
  # A. POISSON RESAMPLING (Resample Genes with Replacement)
  resampled <- clean_genes %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    mutate(split = sample(c("A", "B"), n(), replace = TRUE)) # Partition pS1 / pS2
  
  # B. AGGREGATE (Split-Half)
  # Set A (X-axis)
  stats_A <- resampled %>% filter(split == "A") %>%
    group_by(MAG_id, colony_id, Month) %>%
    summarise(pS_A = sum(dS, na.rm=TRUE)/sum(nS_gene_reference, na.rm=TRUE), .groups="drop")
  
  # Set B (Y-axis)
  stats_B <- resampled %>% filter(split == "B") %>%
    group_by(MAG_id, colony_id, Month) %>%
    summarise(
      pS_B = sum(dS, na.rm=TRUE)/sum(nS_gene_reference, na.rm=TRUE),
      pN_B = sum(dN, na.rm=TRUE)/sum(nN_gene_reference, na.rm=TRUE),
      .groups="drop"
    ) %>% mutate(ratio_B = pN_B / pS_B)
  
  # Join & Filter
  df_boot <- inner_join(stats_A, stats_B, by=c("MAG_id","colony_id","Month")) %>%
    mutate(X = pS_A, Y = ratio_B) %>%
    filter(X > 0, X <= MAX_PS, is.finite(Y), Y >= 0.01)
  
  # C. BINNING (Weighted)
  # We require enough points to fit
  if(nrow(df_boot) > 20) {
    breaks <- 10^seq(log10(MIN_VAL), log10(max(df_boot$X)), length.out = 25)
    
    # We join weights back in (using Set A site counts as approximation for bin weight)
    # Ideally we'd pull sum(nS) from A and B, but for bootstrapping speed, unweighted binning 
    # of the MAG averages is often stable enough. Let's stick to the MAG means here.
    df_bin <- df_boot %>%
      mutate(bin = cut(X, breaks = breaks, include.lowest = TRUE)) %>%
      group_by(bin) %>%
      summarise(mean_X = mean(X), mean_Y = mean(Y), .groups="drop") %>%
      filter(!is.na(mean_X))
    
    # D. FIT (With Bounds!)
    tryCatch({
      fit <- nlsLM(mean_Y ~ (1-f) + f*((1-exp(-K*mean_X))/(K*mean_X)),
                   data = df_bin,
                   start = list(f = 0.85, K = 1000),
                   lower = c(f = 0, K = 0),     # Lower Bounds
                   upper = c(f = 1, K = 50000), # Upper Bounds (Prevents Infinity)
                   control = nls.lm.control(maxiter=50))
      
      # Store Params
      p <- coef(fit)
      boot_params <- rbind(boot_params, data.frame(iter=i, f=p["f"], K=p["K"]))
      
      # Store Line
      y_pred <- (1-p["f"]) + p["f"]*((1-exp(-p["K"]*x_seq))/(p["K"]*x_seq))
      boot_preds <- rbind(boot_preds, data.frame(iter=i, X=x_seq, Y=y_pred))
      
    }, error = function(e) {})
  }
}

# ==========================================
# 3. GENERATE PLOT DATA (SINGLE REALIZATION)
# ==========================================
# We need one "Clean" dataset to plot the points and X markers
set.seed(123)
real_split <- clean_genes %>%
  mutate(split = sample(c("A", "B"), n(), replace = TRUE))

# Aggregation for points
real_A <- real_split %>% filter(split=="A") %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(pS_A = sum(dS)/sum(nS_gene_reference), .groups="drop")

real_B <- real_split %>% filter(split=="B") %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(ratio_B = (sum(dN)/sum(nN_gene_reference)) / (sum(dS)/sum(nS_gene_reference)), .groups="drop")

df_points <- inner_join(real_A, real_B, by=c("MAG_id","colony_id","Month","Genus")) %>%
  mutate(X = pS_A, Y = ratio_B) %>%
  filter(X > 0, X <= MAX_PS, is.finite(Y), Y >= 0.01)

# Genus Medians (The X markers)
genus_summary <- df_points %>%
  group_by(Genus) %>%
  summarise(med_X = median(X), med_Y = median(Y), .groups="drop")

# ==========================================
# 4. CALCULATE RIBBON (95% CI)
# ==========================================
ribbon_data <- boot_preds %>%
  group_by(X) %>%
  summarise(
    ymin = quantile(Y, 0.025),
    ymax = quantile(Y, 0.975),
    ymean = median(Y), # The central line
    .groups="drop"
  )

# Robust Param Stats
med_K <- median(boot_params$K)
ci_K <- quantile(boot_params$K, c(0.025, 0.975))
med_f <- median(boot_params$f)

message(sprintf("Median K: %.0f [%.0f - %.0f]", med_K, ci_K[1], ci_K[2]))

# ==========================================
# 5. PLOT
# ==========================================
p <- ggplot() +
  # A. The Bootstrap Ribbon (Underneath)
  geom_ribbon(data = ribbon_data, aes(x = X, ymin = ymin, ymax = ymax), 
              fill = "#D55E00", alpha = 0.2) +
  
  # B. The Raw Data Cloud (Single Realization)
  geom_point(data = df_points, aes(x = X, y = Y, color = Genus), 
             alpha = 0.3, size = 1.5) +
  
  # C. The Genus Medians (Big X)
  geom_point(data = genus_summary, aes(x = med_X, y = med_Y, color = Genus), 
             shape = 4, size = 6, stroke = 2) +
  
  # D. The Fit Line (Median of Bootstraps)
  geom_line(data = ribbon_data, aes(x = X, y = ymean), 
            color = "#D55E00", size = 1.5) +
  
  # E. Neutral Line
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  
  # F. Scales
  scale_x_log10(labels = scales::label_number()) +
  scale_y_log10(limits = c(0.01, 10), breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
  annotation_logticks() +
  
  labs(
    title = "Decay of Selection Efficacy (Bootstrapped Split-Half)",
    subtitle = paste0("Median Fit: K = ", round(med_K, 0), " (95% CI: ", round(ci_K[1],0), "-", round(ci_K[2],0), ")",
                      "\nPoints: Split-half realization | Ribbon: 95% Bootstrap CI"),
    x = "Synonymous Diversity (pS) [Set A]",
    y = "pN/pS Ratio [Set B]"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

# Print
print(p)


























# Background: one row per gene × colony × Month where that gene is observed
bg_genes <- unique(
  gene_dt[, .(colony_id, MAG_id, corresponding_gene_call,
              Month, Season, Genus)],
  by = c("colony_id", "MAG_id", "corresponding_gene_call", "Month")
)

# sanity
bg_genes[, colony_id := as.character(colony_id)]
sweep_dt[, colony_id := as.character(colony_id)]

# helper to build a "universe" for a given sweep subset
make_sweep_universe <- function(
    sweep_table,
    degeneracy_filter = NULL,
    label = "any"
) {
  sweeps_sub <- copy(sweep_table)
  
  if (!is.null(degeneracy_filter)) {
    sweeps_sub <- sweeps_sub[degeneracy %in% degeneracy_filter]
  }
  
  # Assign sweeps to their endpoint month/season (t2)
  sweeps_gene <- unique(
    sweeps_sub[, .(
      colony_id,
      MAG_id,
      corresponding_gene_call,
      Month  = Month_t2,
      Season = Season_t2,
      Genus
    )],
    by = c("colony_id", "MAG_id", "corresponding_gene_call", "Month")
  )
  
  # Merge with background and flag sweeps
  universe <- copy(bg_genes)
  universe[, is_sweep := 0L]
  
  universe[sweeps_gene,
           is_sweep := 1L,
           on = .(colony_id, MAG_id, corresponding_gene_call, Month, Season, Genus)]
  
  universe[, sweep_type := label]
  universe
}

# Build three universes: all sweeps, 0D, 4D
univ_any  <- make_sweep_universe(sweep_dt, degeneracy_filter = NULL,   label = "any")
univ_0D   <- make_sweep_universe(sweep_dt, degeneracy_filter = "0D",   label = "0D")
univ_4D   <- make_sweep_universe(sweep_dt, degeneracy_filter = "4D",   label = "4D")

sweep_universe <- rbindlist(list(univ_any, univ_0D, univ_4D), use.names = TRUE, fill = TRUE)



# Function to summarise sweep rates by a grouping variable
summarise_sweep_rates <- function(universe, group_var) {
  universe %>%
    as.data.table() %>%
    .[, .(
      n_bg     = .N,
      n_sweep  = sum(is_sweep == 1L),
      sweep_rate = sum(is_sweep == 1L) / .N
    ), by = c(group_var, "sweep_type")] %>%
    .[order(sweep_type, get(group_var))]
}

month_rates <- summarise_sweep_rates(sweep_universe, "Month")
season_rates <- summarise_sweep_rates(sweep_universe, "Season")
genus_rates <- summarise_sweep_rates(sweep_universe, "Genus")



# any sweeps, by Month
tab_month_any <- sweep_universe[sweep_type == "any"] %>%
  count(Month, is_sweep) %>%
  tidyr::pivot_wider(names_from = is_sweep, values_from = n, values_fill = 0)

tab_month_any

# χ²: is is_sweep independent of Month?
chisq_month_any <- chisq.test(as.matrix(tab_month_any[, c("0","1")]))
chisq_month_any


for (st in c("any", "0D", "4D")) {
  cat("\n=== Month × is_sweep, type:", st, "===\n")
  tab <- sweep_universe[sweep_type == st] %>%
    count(Month, is_sweep) %>%
    tidyr::pivot_wider(names_from = is_sweep, values_from = n, values_fill = 0)
  
  print(tab)
  print(chisq.test(as.matrix(tab[, c("0","1")])))
}



for (st in c("any", "0D", "4D")) {
  cat("\n=== Season × is_sweep, type:", st, "===\n")
  tab <- sweep_universe[sweep_type == st] %>%
    count(Season, is_sweep) %>%
    tidyr::pivot_wider(names_from = is_sweep, values_from = n, values_fill = 0)
  
  print(tab)
  print(chisq.test(as.matrix(tab[, c("0","1")])))
}

min_bg <- 50
for (st in c("any", "0D", "4D")) {
  cat("\n=== Genus × is_sweep, type:", st, "===\n")
  tab <- sweep_universe[sweep_type == st] %>%
    as.data.table() %>%
    .[, .(n_bg = .N, n_sweep = sum(is_sweep == 1L)), by = .(Genus)] %>%
    .[n_bg >= min_bg] %>%
    .[, .(Genus,
          `0` = n_bg - n_sweep,
          `1` = n_sweep)]
  
  print(tab)
  print(chisq.test(as.matrix(tab[, c("0","1")])))
}



univ_any_only <- sweep_universe[sweep_type == "any" & !is.na(Season)]

univ_any_only[, Season := factor(Season,
                                 levels = c("Spring", "Summer", "Fall", "Winter"))]

fit_season <- glm(is_sweep ~ Season, data = univ_any_only, family = binomial)

summary(fit_season)
# exp(coef(fit_season)) gives odds ratios relative to the baseline season
exp(coef(fit_season))

gene_sites <- unique(
  gene_dt[, .(MAG_id,
              colony_id,
              corresponding_gene_call,
              Season,
              nN_gene_reference,
              nS_gene_reference)],
  by = c("MAG_id", "colony_id", "corresponding_gene_call", "Season")
)

# sum site "opportunity" per Season
sites_by_season <- gene_sites[, .(
  N_sites_nonsyn = sum(nN_gene_reference, na.rm = TRUE),
  N_sites_syn    = sum(nS_gene_reference, na.rm = TRUE),
  N_genes        = .N
), by = Season][order(Season)]

sites_by_season




# Sweep counts by Season of t2
sweeps_by_season <- sweep_dt[, .(
  sweeps_any = .N,
  sweeps_0D  = sum(degeneracy == "0D", na.rm = TRUE),
  sweeps_4D  = sum(degeneracy == "4D", na.rm = TRUE)
), by = Season_t2]

setnames(sweeps_by_season, "Season_t2", "Season")
sweeps_by_season <- sweeps_by_season[order(Season)]

sweeps_by_season

season_summary <- merge(
  sites_by_season,
  sweeps_by_season,
  by = "Season",
  all.x = TRUE
)

# Replace NA sweep counts with 0
for (col in c("sweeps_any", "sweeps_0D", "sweeps_4D")) {
  season_summary[is.na(get(col)), (col) := 0L]
}

# Rates per million sites
season_summary[, `:=`(
  rate_0D_per_1e6 = sweeps_0D / N_sites_nonsyn * 1e6,
  rate_4D_per_1e6 = sweeps_4D / N_sites_syn    * 1e6,
  rate_any_per_1e6_nonsyn = sweeps_any / N_sites_nonsyn * 1e6  # optional
)]

season_summary


# Only Seasons with nonzero exposure
season_summary_0D <- season_summary[N_sites_nonsyn > 0]
season_summary_4D <- season_summary[N_sites_syn > 0]

# 0D sweeps per nonsyn site
fit_0D <- glm(
  sweeps_0D ~ Season,
  offset = log(N_sites_nonsyn),
  family = poisson,
  data = season_summary_0D
)
summary(fit_0D)
exp(coef(fit_0D))  # incidence rate ratios vs baseline Season

# 4D sweeps per syn site
fit_4D <- glm(
  sweeps_4D ~ Season,
  offset = log(N_sites_syn),
  family = poisson,
  data = season_summary_4D
)
summary(fit_4D)
exp(coef(fit_4D))







gene_sites <- unique(
  gene_dt[, .(
    MAG_id,
    colony_id,
    corresponding_gene_call,
    Season,
    nN_gene_reference,
    nS_gene_reference
  )],
  by = c("MAG_id", "colony_id", "corresponding_gene_call", "Season")
)

## 0D sweeps per gene × Season (using Season_t2)
sweeps_gene_0D <- sweep_dt[degeneracy == "0D",
                           .N,
                           by = .(
                             Season = Season_t2,
                             MAG_id,
                             colony_id,
                             corresponding_gene_call
                           )
]

setnames(sweeps_gene_0D, "N", "n_sweeps_0D")

## 4D sweeps per gene × Season
sweeps_gene_4D <- sweep_dt[degeneracy == "4D",
                           .N,
                           by = .(
                             Season = Season_t2,
                             MAG_id,
                             colony_id,
                             corresponding_gene_call
                           )
]

setnames(sweeps_gene_4D, "N", "n_sweeps_4D")

## Merge sweeps onto gene_sites
gene_rates <- merge(
  gene_sites,
  sweeps_gene_0D,
  by = c("Season", "MAG_id", "colony_id", "corresponding_gene_call"),
  all.x = TRUE
)
gene_rates <- merge(
  gene_rates,
  sweeps_gene_4D,
  by = c("Season", "MAG_id", "colony_id", "corresponding_gene_call"),
  all.x = TRUE
)

# Replace NA sweep counts with 0
gene_rates[is.na(n_sweeps_0D), n_sweeps_0D := 0L]
gene_rates[is.na(n_sweeps_4D), n_sweeps_4D := 0L]

# Per-site gene-level sweep rates (per 1e6 sites for readability)
gene_rates[, rate_0D_per_1e6 := ifelse(
  !is.na(nN_gene_reference) & nN_gene_reference > 0,
  n_sweeps_0D / nN_gene_reference * 1e6,
  NA_real_
)]

gene_rates[, rate_4D_per_1e6 := ifelse(
  !is.na(nS_gene_reference) & nS_gene_reference > 0,
  n_sweeps_4D / nS_gene_reference * 1e6,
  NA_real_
)]






# Option A: restrict to genes that actually had at least one sweep of that type
gene_rates_long <- gene_rates %>%
  as.data.table() %>%
  mutate(
    has_0D = n_sweeps_0D > 0,
    has_4D = n_sweeps_4D > 0
  ) %>%
  pivot_longer(
    cols = c(rate_0D_per_1e6, rate_4D_per_1e6),
    names_to = "sweep_class",
    values_to = "rate_per_1e6"
  )

# Filter to non-NA and at least one sweep of that class
gene_rates_long <- gene_rates_long %>%
  filter(
    sweep_class == "rate_0D_per_1e6" & has_0D |
      sweep_class == "rate_4D_per_1e6" & has_4D
  )

gene_rates_long$sweep_class <- factor(
  gene_rates_long$sweep_class,
  levels = c("rate_0D_per_1e6", "rate_4D_per_1e6"),
  labels = c("0D (nonsyn)", "4D (syn)")
)





gene_rates_long_test <- gene_rates_long[gene_rates_long$sweep_class == "0D (nonsyn)",]


ggplot(gene_rates_long_test,
       aes(x = Season, y = rate_per_1e6)) +
  geom_boxplot(outlier.alpha = 0.3, position = position_dodge(width = 0.8))+
  scale_y_log10()

irr_0D  <- exp(coef(fit_0D))
irr_4D  <- exp(coef(fit_4D))

sub_txt <- glue::glue(
  "Poisson model per-site IRR vs Fall: 0D – Spring {round(irr_0D['SeasonSpring'],2)}×, Summer {round(irr_0D['SeasonSummer'],2)}×, Winter {round(irr_0D['SeasonWinter'],2)}×; ",
  "4D – Spring {round(irr_4D['SeasonSpring'],2)}×, Summer {round(irr_4D['SeasonSummer'],2)}×, Winter {round(irr_4D['SeasonWinter'],2)}×"
)

ggplot(gene_rates_long,
       aes(x = Season, y = rate_per_1e6, fill = sweep_class)) +
  geom_boxplot(outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
  scale_y_log10() +
  scale_fill_manual(
    values = c("0D (nonsyn)" = "#d95f02", "4D (syn)" = "#1b9e77"),
    name   = "Sweep class"
  ) +
  labs(
    x = "Season",
    y = "Gene-level sweep rate (per 10^6 sites)",
    title = "Seasonal bias in sweep rates normalized by nonsyn/syn sites",
    subtitle = sub_txt
  ) +
  theme_minimal(base_size = 12)





### starting from gene_rates

gene_0D_sweeps <- gene_rates %>%
  as.data.table() %>%
  filter(
    !is.na(nN_gene_reference),
    nN_gene_reference > 0,
    n_sweeps_0D > 0        # only genes that actually sweep at 0D
  ) %>%
  mutate(
    rate_0D_per_site = n_sweeps_0D / nN_gene_reference
  )

# Season order
gene_0D_sweeps[, Season := factor(
  Season,
  levels = c("Spring", "Summer", "Fall", "Winter")
)]



stats_0D <- gene_0D_sweeps %>%
  rstatix::pairwise_wilcox_test(
    rate_0D_per_site ~ Season,
    p.adjust.method = "BH"
  ) %>%
  rstatix::add_xy_position(x = "Season")

irr_0D <- exp(coef(fit_0D))
sub_txt_0D <- glue(
  "Poisson per-site IRR vs Fall (0D): ",
  "Spring {round(irr_0D['SeasonSpring'], 2)}×, ",
  "Summer {round(irr_0D['SeasonSummer'], 2)}×, ",
  "Winter {round(irr_0D['SeasonWinter'], 2)}×"
)



p_0D_sweeps <- ggplot(gene_0D_sweeps,
                      aes(x = Season, y = rate_0D_per_site)) +
  geom_boxplot(
    fill = "#d95f02",
    alpha = 0.7,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.1,
    alpha = 0.2,
    size  = 1
  ) +
  scale_y_log10(
    labels = label_scientific(digits = 2)
  ) +
  labs(
    x = "Season",
    y = "0D sweep rate per gene (per site, log₁₀ scale)",
    title = "Seasonal bias in 0D sweep rates among sweep genes",
    subtitle = sub_txt_0D
  ) +
  theme_minimal(base_size = 13)
p_0D_sweeps

p_0D_sweeps +
  stat_pvalue_manual(
    stats_0D,
    label = "p.adj.signif",  # or "p.adj" if you want numbers
    tip.length = 0.01,
    step.increase = 0.06
  )


# 1. Make sure pi_ratio is defined
gene_dt <- gene_dt %>%
  mutate(pi_ratio = pi_nonsyn / pi_syn)

# 2. Filter to usable genes for this analysis
gene_clean <- gene_dt %>%
  filter(
    !is.na(pi_syn),
    !is.na(pi_ratio),
    pi_syn > 0,
    pi_ratio > 0
  )

# 3. Collapse to MAG × colony × Month
magcoltime <- gene_clean %>%
  group_by(MAG_id, colony_id, Month) %>%
  summarise(
    mean_pi_syn   = mean(pi_syn, na.rm = TRUE),
    mean_pi_ratio = mean(pi_ratio, na.rm = TRUE),
    n_genes       = n(),
    .groups = "drop"
  )


# Optional: require some minimum number of genes per MAG×colony×Month
magcoltime <- magcoltime %>%
  filter(n_genes >= 10)  # tweak this threshold as you like

summary(magcoltime$mean_pi_syn)
summary(magcoltime$mean_pi_ratio)


magcoltime <- magcoltime %>%
  left_join(
    gene_dt %>% distinct(MAG_id, Genus),
    by = "MAG_id"
  ) %>%
  mutate(
    Season = season_map[Month]  # reuse your season map if in scope
  )

gene_dt

ggplot(magcoltime,
       aes(x = mean_pi_syn, y = mean_pi_ratio)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(trans="log10")+
  scale_y_continuous(trans="log10")+
  labs(
    x = "Mean synonymous diversity per MAG × colony × month (π_syn)",
    y = "Mean nonsyn/syn diversity ratio (π_nonsyn / π_syn)",
    title = "Purifying-selection pattern across MAG × colony × time snapshots"
  ) +
  theme_minimal(base_size = 13)

fit_pur <- lm(pi_ratio ~ pi_syn, data = gene_dt)
fit_pur




fit_data <- magcoltime %>%
  filter(
    is.finite(mean_pi_syn),
    is.finite(mean_pi_ratio),
    mean_pi_syn > 0,
    mean_pi_ratio > 0
  )

fit_pur_mag <- lm(
  log10(mean_pi_ratio) ~ log10(mean_pi_syn),
  data = fit_data
)

summary(fit_pur_mag)



coef_est <- coef(fit_pur_mag)["log10(mean_pi_syn)"]
coef_se  <- summary(fit_pur_mag)$coefficients["log10(mean_pi_syn)", "Std. Error"]
t_stat   <- coef_est / coef_se
df       <- df.residual(fit_pur_mag)

# one-sided p-value for β1 < 0
p_one_sided <- pt(t_stat, df = df, lower.tail = TRUE)
coef_est
p_one_sided


summ_dt <- gene_dt %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(
    pi_syn   = sum(pi_syn * effective_length_4D_sites, na.rm=TRUE) / sum(effective_length_4D_sites, na.rm=TRUE),
    pi_nonsyn= sum(pi_nonsyn * effective_length_0D_sites, na.rm=TRUE) / sum(effective_length_0D_sites, na.rm=TRUE),
    .groups  = "drop"
  ) %>%
  filter(pi_syn > 0, pi_nonsyn > 0) %>%
  mutate(pi_ratio = pi_nonsyn / pi_syn)

fit <- lm(
  log10(pi_ratio) ~ log10(pi_syn),
  data = summ_dt
)
fit
summary(fit)
with(summ_dt, cor.test(log10(pi_syn), log10(pi_ratio), method = "spearman"))

ggplot(summ_dt, aes(x = pi_syn, y = pi_ratio)) +
  geom_hex(bins = 40) +  # or geom_point(alpha=0.1)
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "red"
  ) +
  labs(
    x = expression(pi[S] ~ "(synonymous nucleotide diversity)"),
    y = expression(pi[N]/pi[S] ~ "(nonsyn/syn diversity ratio)"),
    title = "Selection efficacy increases with neutral diversity",
    subtitle = glue::glue(
      "log10(piN/piS) ~ log10(piS); slope = {round(coef(fit)[2], 3)}, R² = {round(summary(fit)$r.squared, 3)}"
    )
  ) +
  theme_minimal(base_size = 12)

binned <- summ_dt %>%
  mutate(bin = cut(log10(pi_syn), breaks = pretty(log10(pi_syn), n = 15))) %>%
  group_by(bin) %>%
  summarize(
    pi_syn_mid = 10^mean(log10(pi_syn), na.rm = TRUE),
    pi_ratio_mean = mean(pi_ratio, na.rm = TRUE)
  )

# choose target genera
target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium",
                   "commensalibacter", "frischella", "gilliamella",
                   "lactobacillus", "snodgrassella")

# make a plotting column: keep targets as-is, collapse others to "other"
summ_dt <- summ_dt %>%
  mutate(
    Genus_plot = ifelse(Genus %in% target_genera, Genus, "other"),
    # optional alpha column so 'other' is more transparent
    alpha_plot = ifelse(Genus_plot == "other", 0.35, 0.85)
  )
library(RColorBrewer)
# color palette for targets + grey for 'other'
pal_targets <- brewer.pal(n = max(3, length(target_genera)), name = "Set1")
pal_targets <- pal_targets[1:length(target_genera)]         # trim if brewer returns more colors
cols <- c(setNames(pal_targets, target_genera), other = "grey70")

# control legend order: targets then "other"
legend_order <- c(target_genera, "other")

# binned line as before
binned <- summ_dt %>%
  mutate(bin = cut(log10(pi_syn), breaks = pretty(log10(pi_syn), n = 15))) %>%
  group_by(bin) %>%
  summarize(
    pi_syn_mid = 10^mean(log10(pi_syn), na.rm = TRUE),
    pi_ratio_mean = mean(pi_ratio, na.rm = TRUE),
    .groups = "drop"
  )
magcoltime_split
# plot


ggplot(magcoltime_split, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  # MAG points: colored only for target genera, others grey + lower alpha
  geom_point(
    aes(color = Genus, alpha = 0.85),
    size = 1.6,
    position = position_jitter(width = 0, height = 0)
  ) +

  # linear fit in log-log space (same as before)
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "red"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = cols,
    breaks = legend_order,
    labels = legend_order,
    name = "Genus (highlighted)"
  ) +
  scale_alpha_identity() +       # use alpha values as provided (not shown in legend)
  guides(alpha = "none") +       # hide alpha legend
  labs(
    x = expression(pi[S] ~ "(synonymous nucleotide diversity)"),
    y = expression(pi[N]/pi[S] ~ "(nonsyn/syn diversity ratio)"),
    title = "Selection efficacy increases with neutral diversity",
    subtitle = glue::glue(
      "log10(piN/piS) ~ log10(piS); slope = {round(coef(fit)[2], 3)}, R² = {round(summary(fit)$r.squared, 3)}"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    axis.text = element_text(color = "black")
  )




mag_plot <- magcoltime_split %>%
  filter(
    is.finite(pi_syn_2),
    is.finite(pi_ratio_decoupled),
    pi_syn_2 > 0,
    pi_ratio_decoupled > 0
  ) %>%
  mutate(
    Genus_plot = ifelse(Genus %in% target_genera, Genus, "other"),
    alpha_plot = ifelse(Genus_plot == "other", 0.35, 0.85)
  )

## 2. Color palette: target genera + grey for 'other'

pal_targets <- brewer.pal(
  n = max(3, length(target_genera)),
  name = "Set1"
)[seq_along(target_genera)]

cols <- c(
  setNames(pal_targets, target_genera),
  other = "grey70"
)

legend_order <- c(target_genera, "other")

## 3. Good–Garud model predictions on a grid over pi_syn_2

x_grid <- seq(
  min(mag_plot$pi_syn_2, na.rm = TRUE),
  max(mag_plot$pi_syn_2, na.rm = TRUE),
  length.out = 300
)

pred_good <- predict(
  fit_model,
  newdata = data.frame(pi_syn_2 = x_grid)
)

pred_df <- data.frame(
  pi_syn_2   = x_grid,
  pred_good  = pred_good
)

## 4. Final plot: points by Genus, LM smooth, and Good-style curve

p_good_overlay <- ggplot(mag_plot, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  # points (all genera; targets colored, others grey)
  geom_point(
    aes(color = Genus_plot, alpha = alpha_plot),
    size = 1.6,
    position = position_jitter(width = 0, height = 0)
  ) +
  
  # linear fit in log-log space (your original red line)
  geom_smooth(
    method  = "lm",
    formula = y ~ x,
    se      = TRUE,
    color   = "red"
  ) +
  
  # Good-style model (blue dashed curve)
  geom_line(
    data        = pred_df,
    aes(x = pi_syn_2, y = pred_good),
    inherit.aes = FALSE,
    color       = "blue",
    linewidth   = 1.0,
    linetype    = "dashed"
  ) +
  
  # log10 axes
  scale_x_log10() +
  scale_y_log10() +
  
  # manual colors for target genera + 'other'
  scale_color_manual(
    values = cols,
    breaks = legend_order,
    labels = legend_order,
    name   = "Genus"
  ) +
  
  scale_alpha_identity() +  # use alpha as given (no alpha legend)
  guides(alpha = "none") +
  
  labs(
    x = expression(pi[S] ~ "(synonymous diversity, split set 2)"),
    y = expression(pi[N] / pi[S] ~ "(split set 1)"),
    title = "Selection efficacy increases with neutral diversity",
    subtitle = glue(
      "Good-style fit: fd = {round(coef(fit_model)['fd'], 2)}, ",
      "k = {round(coef(fit_model)['k'], 0)}; ",
      "global log10(piN/piS) ~ log10(piS) slope = {round(coef(fit)[2], 3)}, ",
      "R² = {round(summary(fit)$r.squared, 3)}"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    axis.text       = element_text(color = "black")
  )

p_good_overlay



ggplot(magcoltime_split, aes(x = pi_syn_2, y = pi_ratio_decoupled, color = Genus)) +
  geom_point(alpha = 0.5, size = 1.4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(pi[S] ~ "(synonymous diversity, split set 2)"),
    y = expression(pi[N]/pi[S] ~ "(split set 1)"),
    title = "Decoupled purifying selection cloud (all MAG × colony × month)"
  ) +
  theme_minimal()













# re-use pred_df from previous answer
ggplot(magcoltime_split, aes(x = pi_syn_2, y = pi_ratio_decoupled, color = Genus)) +
  geom_point(alpha = 0.7, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  geom_line(
    data = pred_df,
    aes(x = pi_syn_2, y = pred_good),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_hline(yintercept = 1)+
  scale_x_continuous(limits = c(0.001, 0.1), trans="log10") +
  scale_y_continuous(limits = c(0.1, 1), trans="log10") +
  labs(
    x = expression(pi[S] ~ "(synonymous diversity, split set 2)"),
    y = expression(pi[N]/pi[S] ~ "(split set 1)"),
    title = "Selection efficacy vs neutral diversity",
    subtitle = glue(
      "Good-style fd = {round(coef(fit_model)['fd'], 2)}, ",
      "k = {round(coef(fit_model)['k'], 0)}; ",
      "LM slope = {round(coef(fit)[2], 3)}, R² = {round(summary(fit)$r.squared, 3)}"
    )
  ) +
  theme_minimal()

gene_dt

ggplot(data = magcoltime,aes(x=mean_pi_syn, y = mean_pi_ratio))+
  geom_point()+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")

max(magcoltime$mean_pi_syn)



## Filter to valid log-scale points
mag_fit <- magcoltime %>%
  filter(
    is.finite(mean_pi_syn),
    is.finite(mean_pi_ratio),
    mean_pi_syn > 0,
    mean_pi_ratio > 0
  )

summary(mag_fit$mean_pi_syn)
summary(mag_fit$mean_pi_ratio)


fit_good_mag <- tryCatch({
  nls(
    mean_pi_ratio ~ (1 - fd) +
      fd * ((1 - exp(-k * mean_pi_syn)) / (k * mean_pi_syn)),
    data      = mag_fit,
    start     = list(fd = 0.5, k = 1000),
    algorithm = "port",
    lower     = c(fd = 0,   k = 1),
    upper     = c(fd = 1,   k = 1e6)
  )
}, error = function(e) {
  message("Good-style nls failed on magcoltime: ", e$message)
  NULL
})

if (!is.null(fit_good_mag)) coef(fit_good_mag)

# total synonymous sites per MAG×colony×Month (and Genus to match magcoltime)
syn_weight <- gene_dt %>%
  as_tibble() %>%
  filter(
    is.finite(effective_length_4D_sites),
    effective_length_4D_sites > 0
  ) %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(
    syn_sites_total = sum(effective_length_4D_sites, na.rm = TRUE),
    .groups = "drop"
  )
mag_fit
mag_fit_w <- mag_fit %>%
  left_join(syn_weight,
            by = c("MAG_id", "colony_id", "Month", "Genus")) %>%
  mutate(
    w = syn_sites_total,
    w = ifelse(is.na(w) | w <= 0, NA_real_, w)
  ) %>%
  filter(is.finite(w))

fit_good_mag_w <- tryCatch({
  nls(
    sqrt(w) * mean_pi_ratio ~ sqrt(w) * (
      (1 - fd) + fd * ((1 - exp(-k * mean_pi_syn)) / (k * mean_pi_syn))
    ),
    data      = mag_fit_w,
    start     = list(fd = 0.5, k = 1000),
    algorithm = "port",
    lower     = c(fd = 0,   k = 1),
    upper     = c(fd = 1,   k = 1e6)
  )
}, error = function(e) {
  message("Weighted Good-style nls failed: ", e$message)
  NULL
})

if (!is.null(fit_good_mag_w)) coef(fit_good_mag_w)

fit_logistic_mag <- tryCatch({
  nls(
    mean_pi_ratio ~ c + (1 - c) / (1 + (mean_pi_syn / x0)^h),
    data      = mag_fit,
    start     = list(
      c  = 0.2,
      x0 = median(mag_fit$mean_pi_syn, na.rm = TRUE),
      h  = 1
    ),
    algorithm = "port",
    lower     = c(
      c  = 0,
      x0 = min(mag_fit$mean_pi_syn[mag_fit$mean_pi_syn > 0]),
      h  = 0.1
    ),
    upper     = c(
      c  = 1,
      x0 = max(mag_fit$mean_pi_syn),
      h  = 10
    )
  )
}, error = function(e) {
  message("Logistic model failed on magcoltime: ", e$message)
  NULL
})

if (!is.null(fit_logistic_mag)) coef(fit_logistic_mag)



gam_mag <- gam(mean_pi_ratio ~ s(mean_pi_syn, k = 5), data = mag_fit)
summary(gam_mag)


fit_lin_mag <- lm(
  log10(mean_pi_ratio) ~ log10(mean_pi_syn),
  data = mag_fit
)
summary(fit_lin_mag)



x_grid <- seq(
  from = min(mag_fit$mean_pi_syn, na.rm = TRUE),
  to   = max(mag_fit$mean_pi_syn, na.rm = TRUE),
  length.out = 400
)

pred_df <- data.frame(mean_pi_syn = x_grid)

# Good
if (!is.null(fit_good_mag)) {
  pred_df$Good <- predict(fit_good_mag, newdata = pred_df)
}

# Weighted Good
if (!is.null(fit_good_mag_w)) {
  pred_df$Good_w <- predict(fit_good_mag_w, newdata = pred_df)
}

# Logistic
if (!is.null(fit_logistic_mag)) {
  pred_df$Logistic <- predict(fit_logistic_mag, newdata = pred_df)
}

# GAM
pred_df$GAM <- predict(gam_mag, newdata = pred_df)

# Linear log–log (back-transform)
pred_df$Lin_log <- 10^predict(
  fit_lin_mag,
  newdata = data.frame(
    `log10(mean_pi_syn)` = log10(x_grid)
  )
)
min(magcoltime$mean_pi_syn)
# Long format for supplemental multi-line plot
pred_long <- pred_df %>%
  pivot_longer(
    cols      = -mean_pi_syn,
    names_to  = "Model",
    values_to = "fit"
  ) %>%
  filter(!is.na(fit))

min(mag_fit$mean_pi_syn)
p_main <- ggplot(mag_fit, aes(x = mean_pi_syn, y = mean_pi_ratio)) +
  # points by Genus
  geom_point(aes(color = Genus), alpha = 0.6, size = 1.4) +
  
  # Good-style (unweighted) as the main “featured” curve (blue, solid)
  geom_line(
    data = subset(pred_long, Model == "Good"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1.1
  ) +
  # weighted Good (blue dashed)
  geom_line(
    data = subset(pred_long, Model == "Good_w"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  # Logistic (green)
  geom_line(
    data = subset(pred_long, Model == "Logistic"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "darkgreen",
    linewidth = 0.9,
    linetype = "dotdash"
  ) +
  # GAM (black solid)
  geom_line(
    data = subset(pred_long, Model == "GAM"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.0
  ) +
  # Linear log–log (red dotted)
  geom_line(
    data = subset(pred_long, Model == "Lin_log"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 0.8,
    linetype = "dotted"
  ) +
  

  
  labs(
    x = expression(pi[S] ~ "(synonymous diversity, MAG × colony × month)"),
    y = expression(pi[N] / pi[S] ~ "(nonsyn / syn diversity)"),
    title = "Purifying selection: Good-style fit vs alternative models",
    subtitle = {
      if (!is.null(fit_good_mag)) {
        coefs <- coef(fit_good_mag)
        glue(
          "Good-style: fd = {round(coefs['fd'], 2)}, k = {round(coefs['k'], 0)}; ",
          "log–log LM slope = {round(coef(fit_lin_mag)[2], 3)}, R² = {round(summary(fit_lin_mag)$r.squared, 3)}"
        )
      } else {
        "Good-style fit failed; showing other models only"
      }
    },
    color = "Genus"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text  = element_text(color = "black"),
    legend.position = "right"
  )

p_main
max(mag_fit$mean_pi_ratio)




y_min <- min(
  mag_fit$mean_pi_ratio,
  pred_df$Good, pred_df$Good_w,
  pred_df$Logistic, pred_df$GAM, pred_df$Lin_log,
  na.rm = TRUE
)

y_max <- max(
  mag_fit$mean_pi_ratio,
  pred_df$Good, pred_df$Good_w,
  pred_df$Logistic, pred_df$GAM, pred_df$Lin_log,
  na.rm = TRUE
)

## give a tiny pad to avoid points on the frame
y_limits <- c(y_min * 0.95, y_max * 1.05)
x_limits <- c(x_min * 0.95, x_max * 1.05)

y_limits
x_limits




p_main <- ggplot(mag_fit, aes(x = mean_pi_syn, y = mean_pi_ratio)) +
  geom_point(aes(color = Genus), alpha = 0.6, size = 1.4) +
  
  geom_line(
    data = subset(pred_long, Model == "Good"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1.1
  ) +
  geom_line(
    data = subset(pred_long, Model == "Good_w"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  geom_line(
    data = subset(pred_long, Model == "Logistic"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "darkgreen",
    linewidth = 0.9,
    linetype = "dotdash"
  ) +
  geom_line(
    data = subset(pred_long, Model == "GAM"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.0
  ) +
  geom_line(
    data = subset(pred_long, Model == "Lin_log"),
    aes(x = mean_pi_syn, y = fit),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 0.8,
    linetype = "dotted"
  ) +
  scale_x_log10(limits = x_limits, expand = expansion(mult = 0)) + 
  scale_y_log10(limits = y_limits, expand = expansion(mult = 0)) +
  geom_vline(xintercept = 0.001)+
  labs(
    x = expression(pi[S] ~ "(synonymous diversity, MAG × colony × month)"),
    y = expression(pi[N] / pi[S] ~ "(nonsyn / syn diversity)"),
    title = "Purifying selection: Good-style fit vs alternative models",
    subtitle = {
      if (!is.null(fit_good_mag)) {
        coefs <- coef(fit_good_mag)
        glue::glue(
          "Good-style: fd = {round(coefs['fd'], 2)}, k = {round(coefs['k'], 0)}; ",
          "log–log LM slope = {round(coef(fit_lin_mag)[2], 3)}, R² = {round(summary(fit_lin_mag)$r.squared, 3)}"
        )
      } else {
        "Good-style fit failed; showing other models only"
      }
    },
    color = "Genus"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text  = element_text(color = "black"),
    legend.position = "right"
  )

p_main



ggplot(data=mag_fit,aes(x=mean_pi_syn))+
  geom_histogram()








# nice labels for the models
model_labels <- c(
  Good     = "Good-style (unweighted)",
  Good_w   = "Good-style (weighted)",
  Logistic = "Hill / logistic",
  GAM      = "GAM spline",
  Lin_log  = "Linear in log–log"
)

p_supp <- ggplot() +
  # grey background points
  geom_point(
    data = mag_split_fit,
    aes(x = mean_pi_syn, y = mean_pi_ratio),
    color = "grey70",
    alpha = 0.4,
    size  = 1.2
  ) +
  # all model curves, color + linetype by Model
  geom_line(
    data = pred_long,
    aes(x = mean_pi_syn, y = fit, color = Model, linetype = Model),
    linewidth = 1
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      Good     = "blue",
      Good_w   = "blue4",
      Logistic = "darkgreen",
      GAM      = "black",
      Lin_log  = "red3"
    ),
    labels = model_labels,
    name   = "Model"
  ) +
  scale_linetype_manual(
    values = c(
      Good     = "solid",
      Good_w   = "dashed",
      Logistic = "dotdash",
      GAM      = "solid",
      Lin_log  = "dotted"
    ),
    labels = model_labels,
    name   = "Model"
  ) +
  labs(
    x = expression(pi[S]),
    y = expression(pi[N] / pi[S]),
    title = "Comparison of functional forms for purifying selection",
    subtitle = "Grey points: MAG × colony × month; colored lines: fitted models"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text  = element_text(color = "black"),
    legend.position = "right"
  )

p_supp
# ggsave("PurifyingSelection_model_comparison_noGenus.png",
#        p_supp, width = 7, height = 5.5, dpi = 300)


# ggsave("PurifyingSelection_Good_vs_others_colored_Genus.png",
#        p_main, width = 7, height = 5.5, dpi = 300)



model_comp <- tibble(
  Model = character(),
  AIC   = numeric()
)

if (!is.null(fit_good_mag)) {
  model_comp <- add_row(model_comp,
                        Model = "Good",
                        AIC   = AIC(fit_good_mag))
}
if (!is.null(fit_good_mag_w)) {
  model_comp <- add_row(model_comp,
                        Model = "Good_w",
                        AIC   = AIC(fit_good_mag_w))
}
if (!is.null(fit_logistic_mag)) {
  model_comp <- add_row(model_comp,
                        Model = "Logistic",
                        AIC   = AIC(fit_logistic_mag))
}

# For GAM, use extractAIC
gam_aic <- extractAIC(gam_mag)[2]
model_comp <- add_row(model_comp,
                      Model = "GAM",
                      AIC   = gam_aic)

# For linear (log–log), compare on original scale via a simple lm
fit_lin_orig <- lm(mean_pi_ratio ~ log10(mean_pi_syn), data = mag_fit)
model_comp <- add_row(model_comp,
                      Model = "LogLog_lm",
                      AIC   = AIC(fit_lin_orig))

model_comp <- model_comp %>% arrange(AIC)
model_comp




mag_split_fit <- magcoltime_split %>%
  filter(
    is.finite(pi_syn_2),
    is.finite(pi_ratio_decoupled),
    pi_syn_2 > 0,
    pi_ratio_decoupled > 0
  )



## Good-style (unweighted) on split data
fit_good_split <- fit_good_model(mag_split_fit)   # uses pi_syn_2, pi_ratio_decoupled
if (!is.null(fit_good_split)) coef(fit_good_split)

## Optional: weighted Good using total 4D length as weights (like before)

# 4D length per sample (MAG × colony × Month × Genus)
syn_weight_split <- gene_clean %>%
  as_tibble() %>%
  filter(
    is.finite(effective_length_4D_sites),
    effective_length_4D_sites > 0
  ) %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(
    syn_sites_total = sum(effective_length_4D_sites, na.rm = TRUE),
    .groups = "drop"
  )

mag_split_w <- mag_split_fit %>%
  left_join(
    syn_weight_split,
    by = c("MAG_id", "colony_id", "Month", "Genus")
  ) %>%
  mutate(
    w = syn_sites_total,
    w = ifelse(is.na(w) | w <= 0, NA_real_, w)
  ) %>%
  filter(is.finite(w))

fit_good_split_w <- tryCatch({
  nls(
    sqrt(w) * pi_ratio_decoupled ~ sqrt(w) * (
      (1 - fd) + fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2))
    ),
    data      = mag_split_w,
    start     = list(fd = 0.5, k = 1000),
    algorithm = "port",
    lower     = c(fd = 0,   k = 1),
    upper     = c(fd = 1,   k = 1e6)
  )
}, error = function(e) NULL)

if (!is.null(fit_good_split_w)) coef(fit_good_split_w)

## Logistic alternative on split data
fit_logistic_split <- tryCatch({
  nls(
    pi_ratio_decoupled ~ c + (1 - c) / (1 + (pi_syn_2 / x0)^h),
    data      = mag_split_fit,
    start     = list(
      c  = 0.2,
      x0 = median(mag_split_fit$pi_syn_2, na.rm = TRUE),
      h  = 1
    ),
    algorithm = "port",
    lower     = c(
      c  = 0,
      x0 = min(mag_split_fit$pi_syn_2[mag_split_fit$pi_syn_2 > 0]),
      h  = 0.1
    ),
    upper     = c(
      c  = 1,
      x0 = max(mag_split_fit$pi_syn_2),
      h  = 10
    )
  )
}, error = function(e) {
  message("Logistic model failed on split data: ", e$message)
  NULL
})
if (!is.null(fit_logistic_split)) coef(fit_logistic_split)

## GAM on split data
gam_split <- gam(pi_ratio_decoupled ~ s(pi_syn_2, k = 5), data = mag_split_fit)
summary(gam_split)

## Log–log linear model on split data
fit_lin_split <- lm(
  log10(pi_ratio_decoupled) ~ log10(pi_syn_2),
  data = mag_split_fit
)
summary(fit_lin_split)




## Prediction grid over pi_syn_2 on log scale for nicer coverage
x_grid_split <- exp(seq(
  log(min(mag_split_fit$pi_syn_2)),
  log(max(mag_split_fit$pi_syn_2)),
  length.out = 200
))

pred_df_split <- data.frame(pi_syn_2 = x_grid_split)

## Helper: Good-style function
good_fun <- function(x, fd, k) {
  (1 - fd) + fd * ((1 - exp(-k * x)) / (k * x))
}

## Good-style (unweighted) prediction
if (!is.null(fit_good_split)) {
  coefs_good <- coef(fit_good_split)
  pred_df_split$Good <- good_fun(
    x  = pred_df_split$pi_syn_2,
    fd = coefs_good["fd"],
    k  = coefs_good["k"]
  )
} else {
  pred_df_split$Good <- NA_real_
}

## Weighted Good-style prediction (use its fd, k in the same Good function)
if (!is.null(fit_good_split_w)) {
  coefs_good_w <- coef(fit_good_split_w)
  pred_df_split$Good_w <- good_fun(
    x  = pred_df_split$pi_syn_2,
    fd = coefs_good_w["fd"],
    k  = coefs_good_w["k"]
  )
} else {
  pred_df_split$Good_w <- NA_real_
}

## Logistic prediction
pred_df_split$Logistic <- if (!is.null(fit_logistic_split)) {
  predict(fit_logistic_split, newdata = pred_df_split)
} else NA_real_

## GAM prediction
pred_df_split$GAM <- predict(gam_split, newdata = pred_df_split)

## Log–log linear: back-transform to original scale
pred_df_split$Lin_log <- 10^predict(fit_lin_split, newdata = pred_df_split)

## Long format for ggplot
pred_long_split <- pred_df_split %>%
  tidyr::pivot_longer(
    cols      = c("Good", "Good_w", "Logistic", "GAM", "Lin_log"),
    names_to  = "Model",
    values_to = "fit"
  )




p_main_split <- ggplot(mag_split_fit, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  geom_point(
    aes(color = Genus),
    alpha = 0.6,
    size  = 1
  ) +
  geom_line(
    data = subset(pred_long_split, Model == "Good"),
    aes(x = pi_syn_2, y = fit),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 1.1
  ) +
  geom_line(
    data = subset(pred_long_split, Model == "Good_w"),
    aes(x = pi_syn_2, y = fit),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  geom_line(
    data = subset(pred_long_split, Model == "Logistic"),
    aes(x = pi_syn_2, y = fit),
    inherit.aes = FALSE,
    color = "darkgreen",
    linewidth = 0.9,
    linetype = "dotdash"
  ) +
  geom_line(
    data = subset(pred_long_split, Model == "GAM"),
    aes(x = pi_syn_2, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.0
  ) +
  geom_line(
    data = subset(pred_long_split, Model == "Lin_log"),
    aes(x = pi_syn_2, y = fit),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 0.8,
    linetype = "dotted"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(pi[S][2] ~ "(synonymous diversity, independent gene set)"),
    y = expression(pi[N][1] / pi[S][1] ~ "(nonsyn / syn diversity, gene set 1)"),
    title = "Purifying selection with decoupled πS: Good-style vs alternative models",
    subtitle = {
      if (!is.null(fit_good_split)) {
        coefs <- coef(fit_good_split)
        glue::glue(
          "Good-style (split): fd = {round(coefs['fd'], 2)}, k = {round(coefs['k'], 0)}; ",
          "log–log LM slope = {round(coef(fit_lin_split)[2], 3)}, R² = {round(summary(fit_lin_split)$r.squared, 3)}"
        )
      } else {
        "Good-style fit failed on split data; showing other models only"
      }
    },
    color = "Genus"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text       = element_text(color = "black"),
    legend.position = "right"
  )

p_main_split




model_comp




























fit_good_model <- function(dat) {
  # Optionally drop rows with missing or zero pi_syn_2
  dat <- subset(dat, is.finite(pi_ratio_decoupled) & is.finite(pi_syn_2) & pi_syn_2 > 0)
  
  if (nrow(dat) < 5) return(NULL)
  
  tryCatch({
    nls(
      pi_ratio_decoupled ~ (1 - fd) +
        fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2)),
      data      = dat,
      start     = list(fd = 0.5, k = 1000),
      algorithm = "port",               # bounds
      lower     = c(fd = 0,   k = 1),
      upper     = c(fd = 1,   k = 1e6)
    )
  }, error = function(e) NULL)
}

# Check it still works on your original data:
fit_model <- fit_good_model(magcoltime_split)
if (!is.null(fit_model)) coef(fit_model)



########### RANDOM SPLITTING OF PS ###########

make_magcoltime_split_random <- function(gene_dt,
                                         min_genes_per_split = 10,
                                         split_prob = 0.5) {
  gene_dt %>%
    as_tibble() %>%
    # keep genes with usable 0D/4D π and lengths
    filter(
      is.finite(pi_nonsyn),
      is.finite(pi_syn),
      is.finite(effective_length_0D_sites),
      is.finite(effective_length_4D_sites),
      effective_length_0D_sites > 0,
      effective_length_4D_sites > 0
    ) %>%
    group_by(MAG_id, colony_id, Month, Genus) %>%
    # random assignment of each gene to group 1 or 2 == Poisson thinning at gene level
    mutate(split_group = sample(c(1L, 2L), size = n(), replace = TRUE)) %>%
    group_by(MAG_id, colony_id, Month, Genus, split_group) %>%
    summarise(
      L0D = sum(effective_length_0D_sites, na.rm = TRUE),
      L4D = sum(effective_length_4D_sites, na.rm = TRUE),
      pi_nonsyn = ifelse(
        L0D > 0,
        sum(pi_nonsyn * effective_length_0D_sites, na.rm = TRUE) / L0D,
        NA_real_
      ),
      pi_syn = ifelse(
        L4D > 0,
        sum(pi_syn * effective_length_4D_sites, na.rm = TRUE) / L4D,
        NA_real_
      ),
      n_genes = n(),
      .groups = "drop"
    ) %>%
    # wide format: one row per sample with split 1 and split 2 side by side
    pivot_wider(
      names_from  = split_group,
      values_from = c(L0D, L4D, pi_nonsyn, pi_syn, n_genes),
      names_glue  = "{.value}_{split_group}"
    ) %>%
    # keep only samples with reasonable support in both splits
    filter(
      n_genes_1 >= min_genes_per_split,
      n_genes_2 >= min_genes_per_split,
      is.finite(pi_nonsyn_1),
      is.finite(pi_syn_1),
      is.finite(pi_syn_2),
      pi_syn_1 > 0,
      pi_syn_2 > 0
    ) %>%
    transmute(
      MAG_id, colony_id, Month, Genus,
      pi_nonsyn_1 = pi_nonsyn_1,
      pi_syn_1    = pi_syn_1,
      pi_syn_2    = pi_syn_2,
      n_genes_1   = n_genes_1,
      n_genes_2   = n_genes_2,
      pi_ratio_decoupled = pi_nonsyn_1 / pi_syn_1,
      # weight = total 4D length used for both splits
      w = L4D_1 + L4D_2
    ) %>%
    filter(
      is.finite(pi_ratio_decoupled),
      is.finite(pi_syn_2),
      w > 0
    )
}

set.seed(123)
mag_split_test <- make_magcoltime_split_random(gene_dt)
head(mag_split_test)



fit_good_model <- function(dat) {
  dat <- subset(dat, is.finite(pi_ratio_decoupled) & is.finite(pi_syn_2) & pi_syn_2 > 0)
  if (nrow(dat) < 5) return(NULL)
  
  tryCatch({
    nls(
      pi_ratio_decoupled ~ (1 - fd) +
        fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2)),
      data      = dat,
      start     = list(fd = 0.5, k = 1000),
      algorithm = "port",
      lower     = c(fd = 0,   k = 1),
      upper     = c(fd = 1,   k = 1e6)
    )
  }, error = function(e) NULL)
}



set.seed(42)

R_splits <- 100  # or 200 if it's still fast
split_estimates <- replicate(R_splits, {
  mag_split <- make_magcoltime_split_random(gene_dt)
  fit_s     <- fit_good_model(mag_split)
  
  if (is.null(fit_s)) return(c(fd = NA, k = NA))
  coef(fit_s)[c("fd", "k")]
})

split_estimates <- as.data.frame(t(split_estimates))
split_estimates <- subset(split_estimates, is.finite(fd) & is.finite(k))

# Summaries across random splits
fd_split_ci <- quantile(split_estimates$fd, c(0.025, 0.5, 0.975))
k_split_ci  <- quantile(split_estimates$k,  c(0.025, 0.5, 0.975))

fd_split_ci
k_split_ci

# Optional diagnostics
hist(split_estimates$fd, breaks = 30, main = "fd across random splits")
hist(split_estimates$k,  breaks = 30, main = "k across random splits")



syn_weight <- gene_dt %>%
  as_tibble() %>%
  filter(
    is.finite(effective_length_4D_sites),
    effective_length_4D_sites > 0
  ) %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(
    syn_sites_total = sum(effective_length_4D_sites, na.rm = TRUE),
    .groups = "drop"
  )

magcoltime_split_w <- magcoltime_split %>%
  left_join(syn_weight,
            by = c("MAG_id", "colony_id", "Month", "Genus")) %>%
  mutate(
    w = syn_sites_total,
    w = ifelse(is.na(w) | w <= 0, NA_real_, w)
  ) %>%
  filter(
    is.finite(pi_ratio_decoupled),
    is.finite(pi_syn_2),
    is.finite(w)
  )

# Unweighted fit on the subset used for weighting
fit_unweighted_subset <- fit_good_model(magcoltime_split_w)

fit_good_weighted <- tryCatch({
  nls(
    sqrt(w) * pi_ratio_decoupled ~ sqrt(w) * (
      (1 - fd) + fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2))
    ),
    data      = magcoltime_split_w,
    start     = list(fd = 0.5, k = 1000),
    algorithm = "port",
    lower     = c(fd = 0,   k = 1),
    upper     = c(fd = 1,   k = 1e6)
  )
}, error = function(e) NULL)

coef(fit_unweighted_subset)
coef(fit_good_weighted)













fit_simple <- tryCatch({
  nls(
    pi_ratio_decoupled ~ c + (1 - c) / (1 + (pi_syn_2 / x0)^h),
    data      = magcoltime_split,
    start     = list(c = 0.2, x0 = median(magcoltime_split$pi_syn_2, na.rm = TRUE), h = 1),
    algorithm = "port",
    lower     = c(c = 0,    x0 = min(magcoltime_split$pi_syn_2[magcoltime_split$pi_syn_2 > 0]), h = 0.1),
    upper     = c(c = 1,    x0 = max(magcoltime_split$pi_syn_2),                                h = 10)
  )
}, error = function(e) NULL)

if (!is.null(fit_simple)) coef(fit_simple)


if (!is.null(fit_model) && !is.null(fit_simple)) {
  AIC(fit_model, fit_simple)
}


library(mgcv)
gam_fit <- gam(pi_ratio_decoupled ~ s(pi_syn_2, k = 5), data = magcoltime_split)

summary(gam_fit)

# Predict on a grid for overlay plots
x_grid <- seq(min(magcoltime_split$pi_syn_2, na.rm = TRUE),
              max(magcoltime_split$pi_syn_2, na.rm = TRUE),
              length.out = 200)

pred_good <- if (!is.null(fit_model)) {
  predict(fit_model, newdata = data.frame(pi_syn_2 = x_grid))
} else NA

pred_simple <- if (!is.null(fit_simple)) {
  predict(fit_simple, newdata = data.frame(pi_syn_2 = x_grid))
} else NA

pred_gam <- predict(gam_fit, newdata = data.frame(pi_syn_2 = x_grid))

# Basic plot to compare curves
plot(pi_ratio_decoupled ~ pi_syn_2, data = magcoltime_split,
     pch = 16, cex = 0.5)

lines(x_grid, pred_gam,    lwd = 2)                # GAM smoothing
lines(x_grid, pred_good,   lwd = 2, lty = 2)       # Good-style model
lines(x_grid, pred_simple, lwd = 2, lty = 3)       # Hill/logistic model
legend("topright",
       legend = c("GAM", "Good-style", "Hill/logistic"),
       lty    = c(1, 2, 3),
       lwd    = 2,
       bty    = "n")






## 1. Grid over x range
x_grid <- seq(
  min(magcoltime_split$pi_syn_2, na.rm = TRUE),
  max(magcoltime_split$pi_syn_2, na.rm = TRUE),
  length.out = 300
)

grid_df <- data.frame(pi_syn_2 = x_grid)

## 2. Predictions from the three models
grid_df$GAM   <- predict(gam_fit,   newdata = grid_df)
grid_df$Good  <- predict(fit_model, newdata = grid_df)
grid_df$Hill  <- predict(fit_simple, newdata = grid_df)

## 3. Long format for ggplot
pred_long <- grid_df |>
  pivot_longer(
    cols = c(GAM, Good, Hill),
    names_to = "model",
    values_to = "pi_ratio_pred"
  )

## 4. Plot: points colored by Genus, lines by model (linetype)
p_compare_genus <- ggplot(
  magcoltime_split,
  aes(x = pi_syn_2, y = pi_ratio_decoupled)
) +
  # points colored by Genus
  geom_point(
    aes(color = Genus),
    shape = 16,
    size  = 1.5,
    alpha = 0.6
  ) +
  
  # model fits (one line per model, linetype-coded)
  geom_line(
    data = pred_long,
    aes(x = pi_syn_2, y = pi_ratio_pred, linetype = model),
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  
  # linetype legend matching your base-R example
  scale_linetype_manual(
    values = c(
      "GAM"  = "solid",   # lty = 1
      "Good" = "dashed",  # lty = 2
      "Hill" = "dotted"   # lty = 3
    ),
    name   = NULL,
    breaks = c("GAM", "Good", "Hill"),
    labels = c("GAM", "Good-style", "Hill/logistic")
  ) +
  
  # log axes; tweak limits if you want tighter framing
  scale_x_log10() +
  scale_y_log10() +
  
  labs(
    x = expression(pi[S] ~ "(independent gene set 2)"),
    y = expression(pi[N] / pi[S] ~ "(gene set 1)"),
    title    = "Purifying selection across MAG × colony × month",
    subtitle = "Points colored by Genus; lines = GAM, Good-style, Hill/logistic"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text       = element_text(color = "black"),
    legend.position = "right"
  )

p_compare_genus

















fit_good_model <- function(dat) {
  dat <- subset(dat,
                is.finite(pi_ratio_decoupled),
                is.finite(pi_syn_2),
                pi_syn_2 > 0)
  if (nrow(dat) < 5) return(NULL)
  
  tryCatch({
    nls(
      pi_ratio_decoupled ~ (1 - fd) +
        fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2)),
      data      = dat,
      start     = list(fd = 0.5, k = 1000),
      algorithm = "port",
      lower     = c(fd = 0,   k = 1),
      upper     = c(fd = 1,   k = 1e6)
    )
  }, error = function(e) NULL)
}



ggplot(magcoltime_split,aes(x=pi_syn_2,y=pi_ratio_decoupled,color=Genus))+
  geom_point()
  
if (!is.null(fit_model)) {
  coefs <- coef(fit_model)
  print(coefs)
  
  ## grid of x values over the observed range
  x_grid <- seq(
    min(magcoltime_split$pi_syn_2, na.rm = TRUE),
    max(magcoltime_split$pi_syn_2, na.rm = TRUE),
    length.out = 300
  )
  
  pred_df <- data.frame(pi_syn_2 = x_grid)
  pred_df$pi_ratio_pred <- predict(fit_model, newdata = pred_df)
  
  ggplot(magcoltime_split,
         aes(x = pi_syn_2, y = pi_ratio_decoupled, color = Genus)) +
    geom_point(alpha = 0.5, size = 1.4) +
    ## Good–Garud fit
    geom_line(
      data = pred_df,
      aes(x = pi_syn_2, y = pi_ratio_pred),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1.1
    ) +
    labs(
      x = expression(pi[S] ~ "(independent gene set 2)"),
      y = expression(pi[N] / pi[S] ~ "(gene set 1)"),
      title = "Decoupled purifying selection with Good–Garud fit",
      subtitle = glue::glue(
        "fd ≈ {round(coefs['fd'], 2)}, k ≈ {round(coefs['k'], 0)}"
      )
    ) +
    theme_minimal(base_size = 13)
}+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")



# Extract coefficients if fit worked
if(!is.null(fit_model)) {
  coefs <- coef(fit_model)
  print(coefs)

  # Add the fitted line to the plot
  magcoltime_split$predicted <- predict(fit_model, newdata = magcoltime_split)
  
  ggplot(magcoltime_split, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
    geom_point(
      aes(color = Genus, alpha = 0.85),
      size = 1.6,
      position = position_jitter(width = 0, height = 0)
    ) +    geom_line(aes(y = predicted), color = "blue", size = 1.2) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      color = "red"
    ) +
    scale_x_log10() + 
    scale_y_log10() +
    # Limit Y view if noise is high
    labs(
      title = "Fit of Time-Dependent Purifying Selection (Eq S8)",
      subtitle = glue::glue("Est. deleterious fraction (fd): {round(coefs['fd'], 2)}; Selection param (k): {round(coefs['k'], 0)}")
    ) +
    theme_minimal()
} else {
  message("Model failed to converge - data might be too noisy or not fit the decay shape.")
}





if (!is.null(fit_model)) {
  coefs <- coef(fit_model)
  print(coefs)
  
  # Filter to positive, finite values for log scale
  mag_plot <- magcoltime_split %>%
    filter(
      is.finite(pi_syn_2),
      is.finite(pi_ratio_decoupled),
      pi_syn_2 > 0,
      pi_ratio_decoupled > 0
    )
  
  # Add Good-style model predictions
  mag_plot$predicted <- predict(fit_model, newdata = mag_plot)
  
  p_good <- ggplot(mag_plot, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
    # points
    geom_point(
      aes(color = Genus),
      alpha = 0.5,
      size  = 1.4
    ) +
    # Good-style model (main line)
    geom_line(
      aes(y = predicted),
      color = "blue",
      linewidth = 1.1
    ) +
 
    # log10 axes, visually restricted to [0.01, 1]
    scale_y_log10(
      limits = c(1e-1, 1),
      breaks = c(0.01, 0.1, 1)
    ) +
    scale_x_log10(
      limits = c(1e-3, 1e-1),
      breaks = c(0.001, 0.01, 0.1)
    ) +
    labs(
      x = expression(pi[S] ~ "(independent gene set 2)"),
      y = expression(pi[N] / pi[S] ~ "(gene set 1)"),
      title = "Time-dependent purifying selection (Good-style model)",
      subtitle = glue(
        "Est. deleterious fraction (fd): {round(coefs['fd'], 2)}; ",
        "selection parameter (k): {round(coefs['k'], 0)}"
      ),
      color = "Genus"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text  = element_text(color = "black"),
      legend.position = "right"
    )
  
  print(p_good)
  
  ## optional save
  # ggsave("Good_purifying_selection_fit_log_axes_0.01_1.png",
  #        p_good, width = 6.5, height = 5, dpi = 300)
  
} else {
  message("Model failed to converge - data might be too noisy or not fit the decay shape.")
}






# --- 1. Prepare Data & Split Genes ---

# Ensure we are working with the clean gene table
# We need to assign a random "Split Group" (1 or 2) to every gene
set.seed(12345) # Critical for reproducibility
gene_clean_split <- gene_clean %>%
  mutate(
    # Randomly assign every gene to Group 1 or Group 2
    split_group = sample(c(1, 2), size = n(), replace = TRUE)
  )

# --- 2. Aggregate independently ---

magcoltime_split <- gene_clean_split %>%
  group_by(MAG_id, colony_id, Month, Genus) %>%
  summarise(
    # --- Set 1 (Y-axis): piN and piS ---
    # We use weighted means based on effective length, just like your original code
    pi_nonsyn_1 = sum(pi_nonsyn[split_group == 1] * effective_length_0D_sites[split_group == 1], na.rm = TRUE) / 
      sum(effective_length_0D_sites[split_group == 1], na.rm = TRUE),
    
    pi_syn_1    = sum(pi_syn[split_group == 1] * effective_length_4D_sites[split_group == 1], na.rm = TRUE) / 
      sum(effective_length_4D_sites[split_group == 1], na.rm = TRUE),
    
    # --- Set 2 (X-axis): Independent Clock ---
    pi_syn_2    = sum(pi_syn[split_group == 2] * effective_length_4D_sites[split_group == 2], na.rm = TRUE) / 
      sum(effective_length_4D_sites[split_group == 2], na.rm = TRUE),
    
    # Keep track of gene counts to filter out sparse groups
    n_genes_1 = sum(split_group == 1),
    n_genes_2 = sum(split_group == 2),
    .groups = "drop"
  ) %>%
  # Filter: Need valid data in BOTH splits
  filter(
    n_genes_1 >= 5, 
    n_genes_2 >= 5,
    pi_syn_1 > 0, 
    pi_syn_2 > 0,
    pi_nonsyn_1 > 0
  ) %>%
  mutate(
    # The decoupled ratio
    pi_ratio_decoupled = pi_nonsyn_1 / pi_syn_1
  )

# --- 3. The Corrected Plot ---

# Note: X-axis is pi_syn_2, Y-axis uses pi_syn_1
ggplot(magcoltime_split, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  geom_point(alpha = 0.4, color = "grey30") +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "loess", color = "#d95f02", fill = "#d95f02", alpha = 0.2) +
  labs(
    x = expression(pi[S] ~ "(Independent Gene Set 2)"),
    y = expression(pi[N] / pi[S] ~ "(Gene Set 1)"),
    title = "Decoupled Purifying Selection Analysis",
    subtitle = "Using split-gene method to eliminate spurious correlation"
  ) +
  theme_minimal(base_size = 14)


magcoltime_split












# gene-level table with pi_nonsyn, pi_syn, MAG_id, Genus
summ_dt_clean <- summ_dt %>%
  mutate(
    pi_ratio = pi_nonsyn / pi_syn
  ) %>%
  filter(
    !is.na(MAG_id),
    !is.na(Genus),
    is.finite(pi_syn),  pi_syn > 0,
    is.finite(pi_ratio), pi_ratio > 0
  )

# per-MAG summarized values, carrying Genus
mag_obs <- summ_dt_clean %>%
  group_by(MAG_id, Genus) %>%
  summarize(
    mean_pi_ratio = mean(pi_ratio, na.rm = TRUE),
    n_genes       = sum(!is.na(pi_ratio)),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 5)   # or your preferred cutoff

set.seed(1)
B <- 2000

all_ratios  <- summ_dt_clean$pi_ratio
global_mean <- mean(all_ratios, na.rm = TRUE)

permute_one_mag <- function(mean_obs, k, all_ratios, B, center = global_mean) {
  null_means <- replicate(
    B,
    mean(sample(all_ratios, size = k, replace = FALSE))
  )
  
  ci_low  <- quantile(null_means, 0.025, na.rm = TRUE)
  ci_high <- quantile(null_means, 0.975, na.rm = TRUE)
  
  obs_dev   <- abs(mean_obs - center)
  null_devs <- abs(null_means - center)
  
  p_two  <- (sum(null_devs >= obs_dev) + 1) / (B + 1)
  p_high <- (sum(null_means >= mean_obs) + 1) / (B + 1)
  p_low  <- (sum(null_means <= mean_obs) + 1) / (B + 1)
  
  tibble::tibble(
    null_ci_low  = as.numeric(ci_low),
    null_ci_high = as.numeric(ci_high),
    p_two        = p_two,
    p_high       = p_high,
    p_low        = p_low
  )
}

perm_res <- purrr::map2_dfr(
  mag_obs$mean_pi_ratio,
  mag_obs$n_genes,
  ~ permute_one_mag(.x, .y, all_ratios, B)
)

mag_perm <- dplyr::bind_cols(mag_obs, perm_res) %>%
  mutate(
    q_high = p.adjust(p_high, method = "BH"),
    q_low  = p.adjust(p_low,  method = "BH")
  )

# global 95% CI for dashed lines
global_se <- sd(all_ratios, na.rm = TRUE) / sqrt(sum(!is.na(all_ratios)))
global_ci <- c(
  global_mean - 1.96 * global_se,
  global_mean + 1.96 * global_se
)


mag_perm_filt <- mag_perm %>%
  group_by(Genus) %>%
  filter(dplyr::n_distinct(MAG_id) > 5) %>%   # keep genera with > 5 MAGs
  ungroup()



ggplot(
  mag_perm_filt,
  aes(
    x = reorder(Genus, mean_pi_ratio, FUN = median),
    y = mean_pi_ratio
  )
) +
  geom_hline(
    yintercept = global_ci,
    linetype   = "dashed",
    linewidth  = 0.5
  ) +
  geom_hline(
    yintercept = 1,
    linetype   = "dotted",
    linewidth  = 0.5
  ) +
  geom_jitter(
    aes(
      color = case_when(
        q_high < 0.05 ~ "higher-than-null (weaker purifying)",
        q_low  < 0.05 ~ "lower-than-null (stronger purifying)",
        TRUE         ~ "not-significant"
      )
    ),
    width  = 0.2,
    height = 0,
    size   = 2,
    alpha  = 0.9
  ) +
  scale_y_log10(
    breaks = c(1e-2, 1e-1, 1),
    labels = label_math(10^.x)
    # optionally, if your data never goes below 1e-2:
    # limits = c(1e-2, max(mag_perm_genus$mean_pi_ratio, na.rm = TRUE))
  ) +
  scale_color_manual(
    values = c(
      "higher-than-null (weaker purifying)" = "#e41a1c",
      "lower-than-null (stronger purifying)" = "#377eb8",
      "not-significant"                      = "grey50"
    ),
    name = NULL
  ) +
  labs(
    x = "Genus (≥ 6 MAGs)",
    y = expression(paste("Mean MAG-level ", pi[N]/pi[S])),
    title = expression(paste("Per-MAG nonsynonymous diversity ratio by Genus")),
    subtitle = "Points = MAGs; dashed = global 95% null CI; dotted = π[N]/π[S] = 1"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  coord_flip()





# ----- prepare data -----
# ensure we have the MAG-level table we expect
plot_dt <- mag_perm_filt %>%
  select(MAG_id, Genus, mean_pi_ratio, n_genes) %>%
  filter(!is.na(MAG_id), !is.na(Genus), is.finite(mean_pi_ratio)) %>%
  distinct()     # defend against accidental duplicates

# produce a plotting Genus column: keep color for target genera, collapse all others to "other"
plot_dt <- plot_dt %>%
  mutate(
    Genus_plot = ifelse(tolower(Genus) %in% tolower(target_genera), Genus, "other")
  )

# if user didn't supply colors for targets, make a default palette
if (!exists("cols_targets") || is.null(cols_targets)) {
  # create named colors for the target genera
  pal <- RColorBrewer::brewer.pal(min(length(target_genera), 8), "Set2")
  # if more than 8 targets, extend with viridis
  if (length(target_genera) > length(pal)) {
    extra <- viridis::viridis(length(target_genera) - length(pal))
    pal <- c(pal, extra)
  }
  cols_targets <- setNames(pal[1:length(target_genera)], target_genera)
}

# ensure there's an entry for "other"
if (!"other" %in% names(cols_targets)) {
  cols_targets <- c(cols_targets, other = "grey80")
}

# order MAGs by mean_pi_ratio (descending), create factor for y-axis
mag_order <- plot_dt %>%
  group_by(MAG_id) %>%
  summarize(mean_pi_ratio = median(mean_pi_ratio, na.rm = TRUE)) %>%
  arrange(desc(mean_pi_ratio)) %>%
  pull(MAG_id)

plot_dt <- plot_dt %>%
  mutate(MAG_id_f = factor(MAG_id, levels = mag_order))

# ----- the plot -----
p_mag_by_genus <- ggplot(plot_dt, aes(x = mean_pi_ratio, y = MAG_id_f, color = Genus_plot)) +

  # points: one per MAG (you already aggregated to MAG-level)
  geom_point(size = 2.2, alpha = 0.85) +

  # x-axis: log10 with explicit breaks like pathway analysis
  scale_x_log10(
    breaks = c(1e-3, 1e-2, 1e-1, 1),
    labels = c("10^-3", "10^-2", "10^-1", "10^0"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  # color only highlighted genera (others grey)
  scale_color_manual(
    values = cols_targets,
    breaks = c(target_genera, "other") %>% unique(),
    labels = c(target_genera, "other"),
    name = "Genus (highlighted)"
  ) +

  labs(
    x = expression(pi[N]/pi[S]),
    y = NULL,    # hide label for y axis (you requested "no label")
    title = "Per-MAG π[N]/π[S] (MAG rows ordered by mean π[N]/π[S])",
    subtitle = "Points colored for highlighted genera; other genera shown in grey"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    axis.text.y     = element_blank(),    # hide MAG labels
    axis.ticks.y    = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(8, "pt")
  ) +

  # flip so MAGs read top->bottom (optional — remove coord_flip() if you prefer vertical)
  coord_flip()

# print it
p_mag_by_genus














# standardize Month variable names for joining
gene_core <- gene_dt %>%
  select(colony_id, MAG_id, Month, corresponding_gene_call,
         pi_nonsyn, pi_syn, pi_all_sites = pi_all_sites)  # adjust name if needed

sweep_div <- sweep_0D %>%
  # join diversity at t1 (before)
  left_join(
    gene_core,
    by = c("colony_id", "MAG_id", "corresponding_gene_call" = "corresponding_gene_call", "Month_t1" = "Month")
  ) %>%
  rename(
    piN_before  = pi_nonsyn,
    piS_before  = pi_syn,
    piAll_before = pi_all_sites
  ) %>%
  # join diversity at t2 (after)
  left_join(
    gene_core,
    by = c("colony_id", "MAG_id", "corresponding_gene_call" = "corresponding_gene_call", "Month_t2" = "Month")
  ) %>%
  rename(
    piN_after  = pi_nonsyn,
    piS_after  = pi_syn,
    piAll_after = pi_all_sites
  )


sweep_div <- sweep_div %>%
  mutate(
    d_piN   = piN_after  - piN_before,
    d_piS   = piS_after  - piS_before,
    d_piAll = piAll_after - piAll_before
  )

# Restrict to pairs where both before and after are observed
sweep_div_clean <- sweep_div %>%
  filter(!is.na(d_piN), !is.na(d_piS), !is.na(d_piAll))

# One-sided Wilcoxon signed-rank tests (expect negative Δπ if sweeps reduce diversity)
test_piN <- wilcox.test(sweep_div_clean$d_piN, alternative = "less", exact = FALSE)
test_piS <- wilcox.test(sweep_div_clean$d_piS, alternative = "less", exact = FALSE)
test_piAll <- wilcox.test(sweep_div_clean$d_piAll, alternative = "less", exact = FALSE)






# Long format for plotting
d_long <- sweep_div_clean %>%
  select(d_piN, d_piS, d_piAll) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "metric",
    values_to = "delta_pi"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("d_piN", "d_piS", "d_piAll"),
      labels = c("Δπ_N (0D)", "Δπ_S (4D)", "Δπ_all")
    )
  )

# Quick range check (optional)
summary(d_long$delta_pi)




# Histogram + density for each Δπ
ggplot(d_long, aes(x = delta_pi)) +
  geom_histogram(
    bins   = 80,
    alpha  = 0.5,
    colour = "grey20",
    fill   = "grey70"
  ) +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ metric, scales = "free_x") +
  labs(
    x = expression(Delta * pi),
    y = "Count / density",
    title = "Change in nucleotide diversity (π) in genes with 0D sweeps",
    subtitle = "Δπ = π_after − π_before; dashed line at 0"
  ) +
  theme_minimal(base_size = 12)



pi_long <- sweep_div_clean %>%
  transmute(
    colony_id, MAG_id, corresponding_gene_call,
    Month_t1, Month_t2,
    piN_before, piN_after,
    piS_before, piS_after,
    piAll_before, piAll_after
  ) %>%
  pivot_longer(
    cols = c(piN_before, piN_after, piS_before, piS_after, piAll_before, piAll_after),
    names_to = c("metric", "timepoint"),
    names_sep = "_",
    values_to = "pi"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("piN", "piS", "piAll"),
      labels = c("π_N (0D)", "π_S (4D)", "π_all")
    ),
    timepoint = factor(timepoint, levels = c("before", "after"))
  )

ggplot(pi_long, aes(x = timepoint, y = pi, group = interaction(MAG_id, corresponding_gene_call))) +
  geom_line(alpha = 0.02) +
  geom_boxplot(outlier.alpha = 0.1) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = NULL,
    y = "Nucleotide diversity π",
    title = "Nucleotide diversity before vs after 0D sweeps",
    subtitle = "Thin lines = individual genes; boxplots = marginal distributions"
  ) +
  theme_minimal(base_size = 12)




sig_events <- sweep_div_clean %>%
  filter(
    !is.na(d_piN), !is.na(d_piS), !is.na(d_piAll),
    d_piN   < 0,
    d_piS   < 0,
    d_piAll < 0
    
    
    
    
    
  )



sig_table <- sig_events %>%
  arrange(d_piAll) %>%  # most negative at top
  select(
    colony_id,
    MAG_id,
    corresponding_gene_call,
    Month_t1, Month_t2,
    piN_before, piN_after,  d_piN,
    piS_before, piS_after,  d_piS,
    piAll_before, piAll_after, d_piAll
  )

# maybe just inspect the first 50 rows
head(sig_table, 50)

# write to disk for inspection
data.table::fwrite(
  sig_table,
  file = "/path/to/0D_sweep_pi_drop_events.tsv",
  sep  = "\t"
)


sig_events_strong <- sweep_div_clean %>%
  filter(
    !is.na(d_piN), !is.na(d_piS), !is.na(d_piAll),
    d_piN   < 0, d_piS < 0, d_piAll < 0,
    d_piAll < quantile(d_piAll, 0.05, na.rm = TRUE)  # strongest 5% drops
  )







#PERMUTATIONNNNN

gene_dt2 <- gene_dt %>%
  rename(
    piN   = pi_nonsyn,
    piS   = pi_syn,
    piAll = pi_all_sites   # or whatever your "all-sites" π column is
  )

# enforce month order
month_levels <- c("May", "June", "July", "August", "September",
                  "October", "November", "January", "February")

gene_dt2 <- gene_dt2 %>%
  mutate(
    Month = factor(Month, levels = month_levels, ordered = TRUE)
  )

# Build all consecutive intervals and Δπ
all_intervals <- gene_dt2 %>%
  arrange(colony_id, MAG_id, corresponding_gene_call, Month) %>%
  group_by(colony_id, MAG_id, corresponding_gene_call) %>%
  mutate(
    Month_t1  = Month,
    Month_t2  = lead(Month),
    piN_t1    = piN,
    piN_t2    = lead(piN),
    piS_t1    = piS,
    piS_t2    = lead(piS),
    piAll_t1  = piAll,
    piAll_t2  = lead(piAll),
    d_piN     = piN_t2   - piN_t1,
    d_piS     = piS_t2   - piS_t1,
    d_piAll   = piAll_t2 - piAll_t1
  ) %>%
  # keep only valid intervals
  filter(!is.na(Month_t2)) %>%
  ungroup()

# Drop NAs in deltas
all_intervals <- all_intervals %>%
  filter(
    !is.na(d_piN),
    !is.na(d_piS),
    !is.na(d_piAll)
  )

dim(all_intervals)
head(all_intervals)



sweep_div_clean2 <- sweep_div_clean %>%
  filter(
    !is.na(d_piN),
    !is.na(d_piS),
    !is.na(d_piAll)
  )

n_sweeps <- nrow(sweep_div_clean2)
n_sweeps




set.seed(123)

perm_test_delta <- function(
    sweep_delta,      # numeric: Δπ for sweep intervals
    background_delta, # numeric: Δπ for all intervals (null pool)
    B = 1000L         # number of permutations
) {
  sweep_delta  <- sweep_delta[is.finite(sweep_delta)]
  background_delta <- background_delta[is.finite(background_delta)]
  
  n_sweeps <- length(sweep_delta)
  if (n_sweeps == 0L) {
    stop("No sweep delta values provided.")
  }
  if (sum(!is.na(background_delta)) < n_sweeps) {
    stop("Not enough background intervals to sample from.")
  }
  
  obs_mean <- mean(sweep_delta, na.rm = TRUE)
  
  perm_means <- replicate(B, {
    sample_vals <- sample(background_delta, size = n_sweeps, replace = FALSE)
    mean(sample_vals, na.rm = TRUE)
  })
  
  # One-tailed: we expect sweeps to have *more negative* Δπ
  p_val <- (sum(perm_means <= obs_mean) + 1) / (B + 1)
  
  list(
    obs_mean   = obs_mean,
    perm_mean  = mean(perm_means),
    perm_sd    = sd(perm_means),
    p_value    = p_val,
    perm_means = perm_means
  )
}





# Pull the background deltas
bg_d_piN   <- all_intervals$d_piN
bg_d_piS   <- all_intervals$d_piS
bg_d_piAll <- all_intervals$d_piAll

# Pull the sweep deltas
sw_d_piN   <- sweep_div_clean2$d_piN
sw_d_piS   <- sweep_div_clean2$d_piS
sw_d_piAll <- sweep_div_clean2$d_piAll

# Run permutation tests
res_piN   <- perm_test_delta(sw_d_piN,   bg_d_piN,   B = 1000)
res_piS   <- perm_test_delta(sw_d_piS,   bg_d_piS,   B = 1000)
res_piAll <- perm_test_delta(sw_d_piAll, bg_d_piAll, B = 1000)


res_piN
res_piS
res_piAll



plot_perm <- function(res, metric_label) {
  df <- data.frame(perm_mean = res$perm_means)
  
  ggplot(df, aes(x = perm_mean)) +
    geom_histogram(
      bins   = 50,
      fill   = "grey70",
      colour = "grey20"
    ) +
    geom_vline(
      xintercept = res$obs_mean,
      colour     = "red",
      linewidth  = 1
    ) +
    labs(
      x = paste0("Mean Δπ under null (", metric_label, ")"),
      y = "Count of permutations",
      title = paste0("Permutation test for Δπ in sweep genes (", metric_label, ")"),
      subtitle = paste0(
        "Observed mean = ", signif(res$obs_mean, 3),
        "; null mean = ", signif(res$perm_mean, 3),
        "; p(one-tailed) = ", signif(res$p_value, 3)
      )
    ) +
    theme_minimal(base_size = 12)
}
plot_perm(res_piN,   "π_N")
plot_perm(res_piS,   "π_S")
plot_perm(res_piAll, "π_all")





## starting from your per-gene table
gene_dt <- as.data.table(gene_dt)

month_levels <- c("May","June","July","August","September","October","November","January","February")
month_map <- data.table(
  Month     = month_levels,
  month_idx = seq_along(month_levels)
)

gene_dt <- merge(gene_dt, month_map, by = "Month", all.x = TRUE)

## sanity check:
stopifnot("month_idx" %in% names(gene_dt))

gene_small <- gene_dt[, .(
  colony_id,
  MAG_id,
  corresponding_gene_call,
  Month,
  month_idx,
  piAll = pi_all_sites,  # adjust to your actual column names
  piS   = pi_syn,
  piN   = pi_nonsyn
)]

## 1B. Build t1–t2 pairs (consecutive months) per gene
setkey(gene_small, colony_id, MAG_id, corresponding_gene_call, month_idx)

setorder(gene_small, colony_id, MAG_id, corresponding_gene_call, month_idx)

gene_pairs <- gene_small[, .(
  Month_t1       = Month,
  Month_t2       = shift(Month,       type = "lead"),
  month_idx_t1   = month_idx,
  month_idx_t2   = shift(month_idx,   type = "lead"),
  piAll_t1       = piAll,
  piAll_t2       = shift(piAll,       type = "lead"),
  piS_t1         = piS,
  piS_t2         = shift(piS,         type = "lead"),
  piN_t1         = piN,
  piN_t2         = shift(piN,         type = "lead")
), by = .(colony_id, MAG_id, corresponding_gene_call)]


gene_pairs <- gene_pairs[
  !is.na(Month_t2) & (month_idx_t2 == month_idx_t1 + 1L)
]

gene_pairs[, `:=`(
  d_piAll = piAll_t2 - piAll_t1,
  d_piS   = piS_t2   - piS_t1,
  d_piN   = piN_t2   - piN_t1
)]

sweep_dt <- as.data.table(sweep_dt)

gene_sweep_flags <- sweep_dt[
  ,
  .(
    has_0D = any(degeneracy == "0D"),
    has_4D = any(degeneracy == "4D")
  ),
  by = .(colony_id, MAG_id, corresponding_gene_call, Month_t1, Month_t2)
]

gene_pairs <- merge(
  gene_pairs,
  gene_sweep_flags,
  by = c("colony_id","MAG_id","corresponding_gene_call","Month_t1","Month_t2"),
  all.x = TRUE
)

gene_pairs[is.na(has_0D), has_0D := FALSE]
gene_pairs[is.na(has_4D), has_4D := FALSE]





df_dpi_all <- rbind(
  pairs_0D[,   .(status = "0D_sweep", d_piAll, d_piS, d_piN)],
  pairs_null[, .(status = "null",     d_piAll, d_piS, d_piN)]
)

## example: Δπ_all
wilcox.test(
  d_piAll ~ status,
  data = df_dpi_all,
  alternative = "less"  # H1: 0D_sweep has *more negative* Δπ than null
)



ggplot(df_dpi_all, aes(x = status, y = d_piAll)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "",
    y = expression(Delta * pi[all]),
    title = expression("Change in " * pi[all] * " in 0D sweep intervals vs null")
  ) +
  theme_minimal(base_size = 12)






perm_test_dpi <- function(sweep_df, bg_df, value_col = "d_piAll",
                          n_perm = 1000,
                          one_tailed = c("less","greater")) {
  one_tailed <- match.arg(one_tailed)
  
  # observed mean
  obs <- mean(sweep_df[[value_col]], na.rm = TRUE)
  
  # background values
  bg_vals <- bg_df[[value_col]]
  bg_vals <- bg_vals[!is.na(bg_vals)]
  
  n_sweep <- sum(!is.na(sweep_df[[value_col]]))
  
  perm_means <- replicate(
    n_perm,
    mean(sample(bg_vals, n_sweep, replace = TRUE))
  )
  
  if (one_tailed == "less") {
    pval <- mean(perm_means <= obs)
  } else {
    pval <- mean(perm_means >= obs)
  }
  
  list(
    obs_mean   = obs,
    perm_means = perm_means,
    p_value    = pval
  )
}

res_unmatched_piAll <- perm_test_dpi(
  sweep_df  = pairs_0D,
  bg_df     = pairs_null,
  value_col = "d_piAll",
  n_perm    = 2000,
  one_tailed = "less"  # expecting more negative Δπ in sweeps
)

res_unmatched_piS <- perm_test_dpi(
  sweep_df  = pairs_0D,
  bg_df     = pairs_null,
  value_col = "d_piS",
  n_perm    = 2000,
  one_tailed = "less"  # expecting more negative Δπ in sweeps
)

res_unmatched_piAll
res_unmatched_piS


plot_perm <- function(res, title = expression(Delta * pi[all] ~ " in 0D sweep intervals (unmatched null)")) {
  df <- data.frame(perm_means = res$perm_means)
  
  ggplot(df, aes(x = perm_means)) +
    geom_histogram(bins = 40, fill = "grey70", color = "grey30") +
    geom_vline(xintercept = res$obs_mean, color = "red", size = 1.1) +
    labs(
      x = expression("Mean " * Delta * pi[all] * " under null"),
      y = "Count of permutations",
      title = title,
      subtitle = sprintf(
        "Observed mean Δπ = %.3g; null mean = %.3g; p(one-tailed) = %.4g",
        res$obs_mean, mean(res$perm_means), res$p_value
      )
    ) +
    theme_minimal(base_size = 12)
}

plot_perm(res_unmatched_piAll)

plot_perm(res_unmatched_piS)









perm_test_dpi_matched <- function(sweep_df, bg_df,
                                  value_col = "d_piAll",
                                  match_col = "piAll_t1",
                                  n_bins = 10,
                                  n_perm = 1000,
                                  one_tailed = c("less","greater")) {
  one_tailed <- match.arg(one_tailed)
  
  # keep only rows with finite match & value
  sweep_df <- copy(sweep_df)
  bg_df    <- copy(bg_df)
  
  sweep_df <- sweep_df[is.finite(get(match_col)) & is.finite(get(value_col))]
  bg_df    <- bg_df[is.finite(get(match_col)) & is.finite(get(value_col))]
  
  if (nrow(sweep_df) == 0L || nrow(bg_df) == 0L) {
    stop("No valid sweep or background intervals after filtering.")
  }
  
  # define common bins on match_col
  all_match <- c(sweep_df[[match_col]], bg_df[[match_col]])
  breaks <- quantile(all_match, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  # widen edges so everything lands inside
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  
  sweep_df[, bin := cut(get(match_col), breaks = breaks, include.lowest = TRUE)]
  bg_df[,    bin := cut(get(match_col), breaks = breaks, include.lowest = TRUE)]
  
  # drop bins that have no background
  valid_bins <- intersect(
    unique(sweep_df$bin),
    unique(bg_df$bin)
  )
  sweep_df <- sweep_df[bin %in% valid_bins]
  bg_df    <- bg_df[bin %in% valid_bins]
  
  if (nrow(sweep_df) == 0L) {
    stop("After bin-matching, no sweep intervals remain with a valid background.")
  }
  
  # observed mean Δπ in sweeps
  obs <- sweep_df[, mean(get(value_col), na.rm = TRUE)]
  
  # pre-split BG by bin for speed
  bg_split <- split(bg_df, bg_df$bin)
  
  # how many sweeps per bin?
  bin_counts <- sweep_df[, .N, by = bin]
  
  perm_means <- replicate(n_perm, {
    sampled_vals <- lapply(seq_len(nrow(bin_counts)), function(i) {
      b  <- bin_counts$bin[i]
      nb <- bin_counts$N[i]
      pool <- bg_split[[as.character(b)]]
      if (is.null(pool) || nrow(pool) == 0L) return(numeric(0))
      sample(pool[[value_col]], nb, replace = TRUE)
    })
    sampled_vals <- unlist(sampled_vals)
    mean(sampled_vals)
  })
  
  if (one_tailed == "less") {
    pval <- mean(perm_means <= obs)
  } else {
    pval <- mean(perm_means >= obs)
  }
  
  list(
    obs_mean   = obs,
    perm_means = perm_means,
    p_value    = pval
  )
}



res_matched_piAll <- perm_test_dpi_matched(
  sweep_df  = pairs_0D,
  bg_df     = pairs_null,
  value_col = "d_piAll",
  match_col = "pi_all_t1",   # match on starting π_all
  n_bins    = 10,
  n_perm    = 2000,
  one_tailed = "less"
)

plot_perm(res_matched_piAll,
          title = expression(Delta * pi[all] ~ " in 0D sweep intervals (π[all](t1)-matched null)"))


res_matched_piAll


res_matched_piS <- perm_test_dpi_matched(
  sweep_df  = pairs_0D,
  bg_df     = pairs_null,
  value_col = "d_piS",
  match_col = "piS_t1",     # match on starting π_S
  n_bins    = 10,
  n_perm    = 2000,
  one_tailed = "less"
)
res_matched_piS








perm_means_two_groups <- function(x_sweep,
                                  x_null,
                                  n_sample = 100,
                                  n_perm   = 5000,
                                  seed     = 1L) {
  set.seed(seed)
  
  # Clean NAs
  x_sweep <- x_sweep[!is.na(x_sweep)]
  x_null  <- x_null[!is.na(x_null)]
  
  if (length(x_sweep) == 0L || length(x_null) == 0L) {
    stop("One of the groups has zero non-NA values.")
  }
  
  # If group has < n_sample, sample with replacement (still fine)
  replace_sweep <- length(x_sweep) < n_sample
  replace_null  <- length(x_null)  < n_sample
  
  sweep_means <- numeric(n_perm)
  null_means  <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    sweep_means[i] <- mean(sample(x_sweep, n_sample, replace = replace_sweep))
    null_means[i]  <- mean(sample(x_null,  n_sample, replace = replace_null))
  }
  
  diff_means <- sweep_means - null_means
  
  list(
    sweep_means = sweep_means,
    null_means  = null_means,
    diff_means  = diff_means,
    # one-sided p-values you might care about
    p_sweep_gt_null = mean(diff_means >  0),
    p_sweep_lt_null = mean(diff_means <  0)
  )
}


# 0D-sweep intervals
pairs_0D <- gene_pairs[has_0D == TRUE]

# Null intervals: no 0D and no 4D
pairs_null <- gene_pairs[has_0D == FALSE & has_4D == FALSE]

df_unmatched_piS <- rbind(
  pairs_0D[ , .(status = "0D_sweep", piS_t2)],
  pairs_null[, .(status = "null",      piS_t2)]
)

df_unmatched_piS[, .N, by = status]




# 0D sweep intervals
pairs_0D <- gene_pairs[has_0D == TRUE]

# Month-pair distribution for sweeps
sweep_month_pairs <- unique(pairs_0D[, .(colony_id, Month_t1, Month_t2)])

# Null intervals *restricted* to month-pairs used by sweeps
pairs_null_matched <- gene_pairs[
  has_0D == FALSE & has_4D == FALSE
][sweep_month_pairs, on = .(colony_id, Month_t1, Month_t2), nomatch = 0]

df_matched_piS <- rbind(
  pairs_0D[,           .(status = "0D_sweep", piS_t2)],
  pairs_null_matched[, .(status = "null",      piS_t2)]
)

res_unmatched_piS <- perm_means_two_groups(
  x_sweep = df_unmatched_piS[status == "0D_sweep", piS_t2],
  x_null  = df_unmatched_piS[status == "null",     piS_t2],
  n_sample = 500,
  n_perm   = 5000,
  seed     = 1
)

res_unmatched_piS$p_sweep_gt_null  # P(sweep mean > null mean)
res_unmatched_piS$p_sweep_lt_null  # P(sweep mean < null mean)



res_matched_piS <- perm_means_two_groups(
  x_sweep = df_matched_piS[status == "0D_sweep", piS_t2],
  x_null  = df_matched_piS[status == "null",     piS_t2],
  n_sample = 500,
  n_perm   = 5000,
  seed     = 1
)

res_matched_piS$p_sweep_gt_null
res_matched_piS$p_sweep_lt_null



plot_perm_means <- function(res, title = "Permutation distribution of mean πS") {
  df_plot <- rbind(
    data.frame(group = "0D_sweep", mean_piS = res$sweep_means),
    data.frame(group = "null",     mean_piS = res$null_means)
  )
  
  ggplot(df_plot, aes(x = mean_piS, fill = group)) +
    geom_density(alpha = 0.4) +
    labs(
      x = expression("Mean " * pi[S] * " from subsamples (n == 100, 5000 cycles)"),
      y = "Density",
      title = title
    ) +
    scale_fill_manual(values = c("0D_sweep" = "#d95f02", "null" = "#1b9e77")) +
    theme_minimal(base_size = 12)
}

plot_perm_means(res_unmatched_piS,
                title = expression("Mean " * pi[S] * " in 0D-sweep genes vs null (unmatched)"))

plot_perm_means(res_matched_piS,
                title = expression("Mean " * pi[S] * " in 0D-sweep genes vs null (month-pair matched)"))

















# 0D sweep intervals
pairs_0D <- gene_pairs[has_0D == TRUE]

# Null: no 0D, no 4D sweep
pairs_null <- gene_pairs[has_0D == FALSE & has_4D == FALSE]

# πS at t2
df_unmatched_piS_t2 <- rbind(
  pairs_0D[,   .(status = "0D_sweep", piS_t2)],
  pairs_null[, .(status = "null",      piS_t2)]
)

df_unmatched_piS_t1 <- rbind(
  pairs_0D[,   .(status = "0D_sweep", piS_t1)],
  pairs_null[, .(status = "null",      piS_t1)]
)
df_unmatched_piS_t2[, .N, by = status]
df_unmatched_piS_t1[, .N, by = status]



w_t2 <- wilcox.test(piS_t2 ~ status,
                    data = df_unmatched_piS_t2,
                    alternative = "two.sided")

w_t2
df_unmatched_piS_t2[, .(
  median_piS = median(piS_t2, na.rm = TRUE),
  mean_piS   = mean(piS_t2,   na.rm = TRUE),
  n          = .N
), by = status]



perm_means_two_groups <- function(x_sweep,
                                  x_null,
                                  n_sample = 500,
                                  n_perm   = 5000,
                                  seed     = 1L) {
  set.seed(seed)
  
  x_sweep <- x_sweep[!is.na(x_sweep)]
  x_null  <- x_null[!is.na(x_null)]
  
  if (length(x_sweep) == 0L || length(x_null) == 0L) {
    stop("One of the groups has zero non-NA values.")
  }
  
  replace_sweep <- length(x_sweep) < n_sample
  replace_null  <- length(x_null)  < n_sample
  
  sweep_means <- numeric(n_perm)
  null_means  <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    sweep_means[i] <- mean(sample(x_sweep, n_sample, replace = replace_sweep))
    null_means[i]  <- mean(sample(x_null,  n_sample, replace = replace_null))
  }
  
  diff_means <- sweep_means - null_means
  
  list(
    sweep_means     = sweep_means,
    null_means      = null_means,
    diff_means      = diff_means,
    p_sweep_gt_null = mean(diff_means > 0),
    p_sweep_lt_null = mean(diff_means < 0)
  )
}

res_piS_t2 <- perm_means_two_groups(
  x_sweep = df_unmatched_piS_t2[status == "0D_sweep", piS_t2],
  x_null  = df_unmatched_piS_t2[status == "null",     piS_t2],
  n_sample = 500,
  n_perm   = 5000,
  seed     = 1
)

res_piS_t2$p_sweep_gt_null  # P(mean piS_t2(sweep) > mean piS_t2(null))
res_piS_t2$p_sweep_lt_null  # P(mean piS_t2(sweep) < mean piS_t2(null))



res_piS_t1 <- perm_means_two_groups(
  x_sweep = df_unmatched_piS_t1[status == "0D_sweep", piS_t1],
  x_null  = df_unmatched_piS_t1[status == "null",     piS_t1],
  n_sample = 500,
  n_perm   = 5000,
  seed     = 1
)

res_piS_t1$p_sweep_gt_null
res_piS_t1$p_sweep_lt_null













gene_pairs[, d_piS := piS_t2 - piS_t1]

df_delta <- gene_pairs[
  ,
  .(status = ifelse(has_0D, "0D_sweep",
                    ifelse(!has_0D & !has_4D, "null", NA_character_)),
    d_piS, piS_t1)
][!is.na(status)]

w_delta <- wilcox.test(d_piS ~ status,
                       data = df_delta[status %in% c("0D_sweep", "null")],
                       alternative = "less")  # if you expect d_piS_sweep < d_piS_null

df_delta[status %in% c("0D_sweep","null"),
         .(median_d = median(d_piS, na.rm=TRUE),
           mean_d   = mean(d_piS, na.rm=TRUE)),
         by = status]




ggplot(df_delta[status %in% c("0D_sweep","null")],
       aes(x = status, y = d_piS, fill = status)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    x = "",
    y = expression(Delta*pi[S] ~ "(t[2] - t[1])"),
    title = expression("Change in " * pi[S] * " in 0D-sweep vs non-swept genes")
  ) +
  scale_fill_manual(values = c("0D_sweep" = "#d95f02", "null" = "#1b9e77")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")+
  scale_y_continuous(trans="log10")+
  coord_cartesian(ylim=c(0.00001,0.1))



fit_delta <- lm(d_piS ~ status + piS_t1,
                data = df_delta[status %in% c("0D_sweep","null")])

summary(fit_delta)




perm_delta_piS <- function(df, n_perm = 5000, n_sample = 500, seed = 1L) {
  set.seed(seed)
  df <- df[!is.na(d_piS) & status %in% c("0D_sweep","null")]
  
  x_sweep <- df[status == "0D_sweep", d_piS]
  x_null  <- df[status == "null",     d_piS]
  
  # same perm_means_two_groups idea but for d_piS
  perm_means_two_groups(x_sweep, x_null,
                        n_sample = n_sample,
                        n_perm   = n_perm,
                        seed     = seed)
}
perm_delta_piS(df_delta)
with(df_delta[status %in% c("0D_sweep","null")],
     tapply(d_piS, status, mean, na.rm = TRUE))










fit_good_model <- function(dat) {
  # Optionally drop rows with missing or zero pi_syn_2
  dat <- subset(dat, is.finite(pi_ratio_decoupled) & is.finite(pi_syn_2) & pi_syn_2 > 0)
  
  if (nrow(dat) < 5) return(NULL)
  
  tryCatch({
    nls(
      pi_ratio_decoupled ~ (1 - fd) +
        fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2)),
      data      = dat,
      start     = list(fd = 0.5, k = 1000),
      algorithm = "port",               # bounds
      lower     = c(fd = 0,   k = 1),
      upper     = c(fd = 1,   k = 1e6)
    )
  }, error = function(e) NULL)
}

# Check it still works on your original data:
fit_model <- fit_good_model(magcoltime_split)
if (!is.null(fit_model)) coef(fit_model)


set.seed(1)

B <- 500   # number of bootstrap replicates; bump to 1000 if fast

boot_estimates <- replicate(B, {
  idx <- sample(seq_len(nrow(magcoltime_split)), replace = TRUE)
  dat_b <- magcoltime_split[idx, , drop = FALSE]
  
  fit_b <- fit_good_model(dat_b)
  if (is.null(fit_b)) return(c(fd = NA, k = NA))
  
  coef(fit_b)[c("fd", "k")]
})

boot_estimates <- t(boot_estimates)
boot_estimates <- as.data.frame(boot_estimates)

# Drop failed fits
boot_estimates <- subset(boot_estimates, is.finite(fd) & is.finite(k))

# Bootstrap CIs
fd_ci <- quantile(boot_estimates$fd, c(0.025, 0.5, 0.975), na.rm = TRUE)
k_ci  <- quantile(boot_estimates$k,  c(0.025, 0.5, 0.975), na.rm = TRUE)

fd_ci
k_ci

# Optional: quick diagnostic plots
hist(boot_estimates$fd, breaks = 30, main = "Bootstrap fd", xlab = "fd")
hist(boot_estimates$k,  breaks = 30, main = "Bootstrap k",  xlab = "k")



make_split_and_fit <- function(dat) {
  # assumes dat has syn_diff_total and syn_sites and also whatever you need for piN
  
  # 1. Random split of synonymous differences into 2 independent Poisson thinned counts
  n_total <- dat$syn_diff_total
  
  # Poisson thinning ≡ Binomial split conditional on total
  n1 <- rbinom(length(n_total), size = n_total, prob = 0.5)
  n2 <- n_total - n1
  
  dat$pi_syn_1 <- n1 / dat$syn_sites
  dat$pi_syn_2 <- n2 / dat$syn_sites
  
  # 2. Recompute pi_ratio_decoupled using pi_syn_1 in the denominator.
  #    Change this to match however you defined pi_ratio_decoupled originally.
  #    Example if you have piN already:
  # dat$pi_ratio_decoupled <- dat$piN / dat$pi_syn_1
  
  # 3. Fit the Good-style model on this split
  fit_good_model(dat)
}





# choose target genera
target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium",
                   "commensalibacter", "frischella", "gilliamella",
                   "lactobacillus", "snodgrassella")

# make a plotting column: keep targets as-is, collapse others to "other"
summ_dt <- summ_dt %>%
  mutate(
    Genus_plot = ifelse(Genus %in% target_genera, Genus, "other"),
    # optional alpha column so 'other' is more transparent
    alpha_plot = ifelse(Genus_plot == "other", 0.35, 0.85)
  )

# color palette for targets + grey for 'other'
pal_targets <- brewer.pal(n = max(3, length(target_genera)), name = "Set2")
pal_targets <- pal_targets[1:length(target_genera)]         # trim if brewer returns more colors
cols <- c(setNames(pal_targets, target_genera), other = "grey70")

# control legend order: targets then "other"
legend_order <- c(target_genera, "other")

# binned line as before
binned <- summ_dt %>%
  mutate(bin = cut(log10(pi_syn), breaks = pretty(log10(pi_syn), n = 15))) %>%
  group_by(bin) %>%
  summarize(
    pi_syn_mid = 10^mean(log10(pi_syn), na.rm = TRUE),
    pi_ratio_mean = mean(pi_ratio, na.rm = TRUE),
    .groups = "drop"
  )

# plot
ggplot(summ_dt, aes(x = pi_syn, y = pi_ratio)) +
  # MAG points: colored only for target genera, others grey + lower alpha
  geom_point(
    aes(color = Genus_plot, alpha = alpha_plot),
    size = 1.6,
    position = position_jitter(width = 0, height = 0)
  ) +
  # trend line from binned summary
  geom_line(
    data = binned,
    aes(x = pi_syn_mid, y = pi_ratio_mean),
    color = "orange", size = 1
  ) +
  # linear fit in log-log space (same as before)
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "red"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = cols,
    breaks = legend_order,
    labels = legend_order,
    name = "Genus (highlighted)"
  ) +
  scale_alpha_identity() +       # use alpha values as provided (not shown in legend)
  guides(alpha = "none") +       # hide alpha legend
  labs(
    x = expression(pi[S] ~ "(synonymous nucleotide diversity)"),
    y = expression(pi[N]/pi[S] ~ "(nonsyn/syn diversity ratio)"),
    title = "Selection efficacy increases with neutral diversity",
    
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    axis.text = element_text(color = "black")
  )










# 1) Aggregate to your "summ_dt" (exactly your definition)
make_summ_dt <- function(gene_dt) {
  gene_dt %>%
    as_tibble() %>%
    group_by(MAG_id, colony_id, Month, Genus) %>%
    summarise(
      pi_syn    = sum(pi_syn    * effective_length_4D_sites, na.rm = TRUE) / sum(effective_length_4D_sites, na.rm = TRUE),
      pi_nonsyn = sum(pi_nonsyn * effective_length_0D_sites, na.rm = TRUE) / sum(effective_length_0D_sites, na.rm = TRUE),
      # weights: total "effective" sites contributing (these will be useful later)
      syn_sites_total    = sum(effective_length_4D_sites, na.rm = TRUE),
      nonsyn_sites_total = sum(effective_length_0D_sites, na.rm = TRUE),
      n_genes = dplyr::n(),
      .groups = "drop"
    ) %>%
    filter(is.finite(pi_syn), is.finite(pi_nonsyn), pi_syn > 0, pi_nonsyn > 0) %>%
    mutate(pi_ratio = pi_nonsyn / pi_syn)
}


# 2) Make pi_syn_1 and pi_syn_2 by random gene-level thinning inside each MAG/colony/Month/Genus
#    Then compute pi_ratio_decoupled = pi_nonsyn / pi_syn_1 and set w = total synonymous effective length.
make_decoupled_dt <- function(gene_dt,
                              summ_dt = NULL,
                              min_genes_per_split = 10,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # If summ_dt not provided, compute it.
  if (is.null(summ_dt)) summ_dt <- make_summ_dt(gene_dt)
  
  # Keep only rows needed for building pi_syn splits
  syn_gene <- gene_dt %>%
    as_tibble() %>%
    filter(
      is.finite(pi_syn),
      is.finite(effective_length_4D_sites),
      effective_length_4D_sites > 0
    ) %>%
    select(MAG_id, colony_id, Month, Genus, pi_syn, effective_length_4D_sites)
  
  # Randomly split genes into two sets within each group
  syn_split <- syn_gene %>%
    group_by(MAG_id, colony_id, Month, Genus) %>%
    mutate(split_group = sample.int(2L, size = dplyr::n(), replace = TRUE)) %>%
    group_by(MAG_id, colony_id, Month, Genus, split_group) %>%
    summarise(
      L4D = sum(effective_length_4D_sites, na.rm = TRUE),
      pi_syn_split = sum(pi_syn * effective_length_4D_sites, na.rm = TRUE) / L4D,
      n_genes_split = dplyr::n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = split_group,
      values_from = c(L4D, pi_syn_split, n_genes_split),
      names_glue  = "{.value}_{split_group}"
    ) %>%
    # enforce enough genes in both splits
    filter(
      n_genes_split_1 >= min_genes_per_split,
      n_genes_split_2 >= min_genes_per_split,
      is.finite(pi_syn_split_1), is.finite(pi_syn_split_2),
      pi_syn_split_1 > 0, pi_syn_split_2 > 0
    ) %>%
    transmute(
      MAG_id, colony_id, Month, Genus,
      pi_syn_1 = pi_syn_split_1,
      pi_syn_2 = pi_syn_split_2,
      n_genes_1 = n_genes_split_1,
      n_genes_2 = n_genes_split_2,
      w = L4D_1 + L4D_2   # total synonymous effective length used across both splits
    )
  
  # Join to get pi_nonsyn (from summ_dt) and compute decoupled ratio
  decoupled <- summ_dt %>%
    left_join(syn_split, by = c("MAG_id", "colony_id", "Month", "Genus")) %>%
    filter(
      is.finite(pi_nonsyn),
      is.finite(pi_syn_1), is.finite(pi_syn_2),
      pi_syn_1 > 0, pi_syn_2 > 0,
      is.finite(w), w > 0
    ) %>%
    mutate(
      pi_ratio_decoupled = pi_nonsyn / pi_syn_1
    ) %>%
    filter(is.finite(pi_ratio_decoupled))
  
  decoupled
}


fit_good_model <- function(dat) {
  dat <- subset(dat, is.finite(pi_ratio_decoupled) & is.finite(pi_syn_2) & pi_syn_2 > 0)
  if (nrow(dat) < 20) return(NULL)
  
  tryCatch({
    nls(
      pi_ratio_decoupled ~ (1 - fd) + fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2)),
      data      = dat,
      start     = list(fd = 0.5, k = 1000),
      algorithm = "port",
      lower     = c(fd = 0,   k = 1),
      upper     = c(fd = 1,   k = 1e6)
    )
  }, error = function(e) NULL)
}

fit_good_weighted <- function(dat) {
  dat <- subset(dat,
                is.finite(pi_ratio_decoupled) & is.finite(pi_syn_2) & pi_syn_2 > 0 &
                  is.finite(w) & w > 0)
  if (nrow(dat) < 20) return(NULL)
  
  tryCatch({
    nls(
      sqrt(w) * pi_ratio_decoupled ~ sqrt(w) * (
        (1 - fd) + fd * ((1 - exp(-k * pi_syn_2)) / (k * pi_syn_2))
      ),
      data      = dat,
      start     = list(fd = 0.5, k = 1000),
      algorithm = "port",
      lower     = c(fd = 0,   k = 1),
      upper     = c(fd = 1,   k = 1e6)
    )
  }, error = function(e) NULL)
}

fit_hill <- function(dat) {
  dat <- subset(dat, is.finite(pi_ratio_decoupled) & is.finite(pi_syn_2) & pi_syn_2 > 0)
  if (nrow(dat) < 20) return(NULL)
  
  tryCatch({
    nls(
      pi_ratio_decoupled ~ c + (1 - c) / (1 + (pi_syn_2 / x0)^h),
      data      = dat,
      start     = list(c = 0.2, x0 = median(dat$pi_syn_2, na.rm = TRUE), h = 1),
      algorithm = "port",
      lower     = c(c = 0,  x0 = min(dat$pi_syn_2), h = 0.1),
      upper     = c(c = 1,  x0 = max(dat$pi_syn_2), h = 10)
    )
  }, error = function(e) NULL)
}

fit_gam <- function(dat) {
  dat <- subset(dat, is.finite(pi_ratio_decoupled) & is.finite(pi_syn_2) & pi_syn_2 > 0)
  gam(pi_ratio_decoupled ~ s(pi_syn_2, k = 5), data = dat)
}






summ_dt <- make_summ_dt(gene_dt)

# If you have Genus_plot / alpha_plot and want to carry them through,
# join them into summ_dt here (or compute them before this step).
# (Leaving as-is since you didn't show how Genus_plot/alpha_plot are created.)

set.seed(1)
decoupled_dt <- make_decoupled_dt(
  gene_dt  = gene_dt,
  summ_dt  = summ_dt,
  min_genes_per_split = 1,
  seed = 1
)

dim(summ_dt)
dim(decoupled_dt)





fit_good   <- fit_good_model(decoupled_dt)
fit_good_w <- fit_good_weighted(decoupled_dt)
fit_h      <- fit_hill(decoupled_dt)
gam_fit    <- fit_gam(decoupled_dt)

coef(fit_good)
coef(fit_good_w)
coef(fit_h)
summary(gam_fit)

AIC(fit_good, fit_good_w, fit_h)


set.seed(2)

B <- 500
boot_est <- replicate(B, {
  idx <- sample.int(nrow(decoupled_dt), replace = TRUE)
  dat_b <- decoupled_dt[idx, , drop = FALSE]
  
  fb <- fit_good_model(dat_b)
  if (is.null(fb)) return(c(fd = NA, k = NA))
  coef(fb)[c("fd", "k")]
})

boot_est <- as.data.frame(t(boot_est))
boot_est <- subset(boot_est, is.finite(fd) & is.finite(k))

quantile(boot_est$fd, c(0.025, 0.5, 0.975))
quantile(boot_est$k,  c(0.025, 0.5, 0.975))




set.seed(3)

R_splits <- 100

split_est <- replicate(R_splits, {
  dat_r <- make_decoupled_dt(
    gene_dt = gene_dt,
    summ_dt = summ_dt,
    min_genes_per_split = 10
    # no seed: allow randomness each replicate
  )
  fr <- fit_good_model(dat_r)
  if (is.null(fr)) return(c(fd = NA, k = NA))
  coef(fr)[c("fd", "k")]
})

split_est <- as.data.frame(t(split_est))
split_est <- subset(split_est, is.finite(fd) & is.finite(k))

fd_split_ci <- quantile(split_est$fd, c(0.025, 0.5, 0.975))
k_split_ci  <- quantile(split_est$k,  c(0.025, 0.5, 0.975))

fd_split_ci
k_split_ci





x_grid <- seq(min(decoupled_dt$pi_syn_2, na.rm = TRUE),
              max(decoupled_dt$pi_syn_2, na.rm = TRUE),
              length.out = 300)

pred_good <- predict(fit_good, newdata = data.frame(pi_syn_2 = x_grid))
pred_hill <- if (!is.null(fit_h)) predict(fit_h, newdata = data.frame(pi_syn_2 = x_grid)) else NA
pred_gam  <- predict(gam_fit, newdata = data.frame(pi_syn_2 = x_grid))

plot(pi_ratio_decoupled ~ pi_syn_2, data = decoupled_dt, pch = 16, cex = 0.5,
     xlab = expression(pi[S]~"(split 2)"),
     ylab = expression(pi[N]/pi[S]~"(split 1 denom)"))

lines(x_grid, pred_gam,  lwd = 2)
lines(x_grid, pred_good, lwd = 2, lty = 2)
if (all(is.finite(pred_hill))) lines(x_grid, pred_hill, lwd = 2, lty = 3)

legend("topright",
       legend = c("GAM", "Good-style", "Hill/logistic"),
       lty    = c(1, 2, 3),
       lwd    = 2,
       bty    = "n")
decoupled_dt
ggplot(data = decoupled_dt,aes(x=pi_syn_2,y=pi_ratio_decoupled))+
  geom_point()+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(trans="log10")

coef(fit_good)
coef(fit_good_w)

min(decoupled_dt$pi_syn_2)








make_decoupled_dt_splitboth <- function(gene_dt,
                                        min_genes_per_split = 1,
                                        min_4D_sites_per_split = 0,
                                        min_0D_sites_per_split = 0,
                                        seed = NULL) {
  if (!is.null(seed) ) set.seed(seed)
  
  # Keep genes with needed columns
  g <- gene_dt %>%
    as_tibble() %>%
    filter(
      is.finite(pi_syn), is.finite(pi_nonsyn),
      is.finite(effective_length_4D_sites), effective_length_4D_sites > 0,
      is.finite(effective_length_0D_sites), effective_length_0D_sites > 0
    ) %>%
    select(MAG_id, colony_id, Month, Genus,
           pi_syn, pi_nonsyn,
           effective_length_4D_sites, effective_length_0D_sites)
  
  # Randomly split genes within each MAG/colony/Month/Genus
  split_sum <- g %>%
    group_by(MAG_id, colony_id, Month, Genus) %>%
    mutate(split_group = sample.int(2L, size = n(), replace = TRUE)) %>%
    group_by(MAG_id, colony_id, Month, Genus, split_group) %>%
    summarise(
      L4D = sum(effective_length_4D_sites, na.rm = TRUE),
      L0D = sum(effective_length_0D_sites, na.rm = TRUE),
      pi_syn_split = sum(pi_syn * effective_length_4D_sites, na.rm = TRUE) / L4D,
      pi_nonsyn_split = sum(pi_nonsyn * effective_length_0D_sites, na.rm = TRUE) / L0D,
      n_genes_split = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = split_group,
      values_from = c(L4D, L0D, pi_syn_split, pi_nonsyn_split, n_genes_split),
      names_glue  = "{.value}_{split_group}"
    ) %>%
    # Filters: ensure both splits exist and have enough genes/sites
    filter(
      n_genes_split_1 >= min_genes_per_split,
      n_genes_split_2 >= min_genes_per_split,
      L4D_1 >= min_4D_sites_per_split,
      L4D_2 >= min_4D_sites_per_split,
      L0D_1 >= min_0D_sites_per_split,
      L0D_2 >= min_0D_sites_per_split
    ) %>%
    transmute(
      MAG_id, colony_id, Month, Genus,
      n_genes_1 = n_genes_split_1,
      n_genes_2 = n_genes_split_2,
      L4D_1, L4D_2, L0D_1, L0D_2,
      pi_syn_1 = pi_syn_split_1,
      pi_syn_2 = pi_syn_split_2,
      pi_nonsyn_1 = pi_nonsyn_split_1,
      pi_nonsyn_2 = pi_nonsyn_split_2,
      # decoupled ratio uses split 1 numerator + denominator, x uses split 2
      pi_ratio_decoupled = pi_nonsyn_1 / pi_syn_1,
      # weights (choose what you like; I include a few)
      w4D = L4D_1 + L4D_2,
      w0D = L0D_1 + L0D_2,
      w_min4D = pmin(L4D_1, L4D_2),
      w_min0D = pmin(L0D_1, L0D_2)
    ) %>%
    filter(
      is.finite(pi_ratio_decoupled),
      is.finite(pi_syn_2),
      pi_syn_1 > 0, pi_syn_2 > 0,
      pi_ratio_decoupled > 0
    )
  
  split_sum
}

# Build once with minimal filtering (so we can inspect distributions)
dec_min <- make_decoupled_dt_splitboth(gene_dt,
                                       min_genes_per_split = 1,
                                       min_4D_sites_per_split = 0,
                                       min_0D_sites_per_split = 0,
                                       seed = 1)

# Key derived quantities
dec_min <- dec_min %>%
  mutate(
    n_genes_total = n_genes_1 + n_genes_2,
    n_genes_min = pmin(n_genes_1, n_genes_2)
  )

# Gene-count distributions
p1 <- ggplot(dec_min, aes(x = n_genes_total)) +
  geom_histogram(bins = 60) +
  theme_classic() +
  labs(x = "Total genes (split1 + split2)", y = "Count")

p2 <- ggplot(dec_min, aes(x = n_genes_min)) +
  geom_histogram(bins = 60) +
  theme_classic() +
  labs(x = "Min genes per split (min(n_genes_1, n_genes_2))", y = "Count")

# Effective-site distributions (often better than gene count)
p3 <- ggplot(dec_min, aes(x = w_min4D)) +
  geom_histogram(bins = 60) +
  scale_x_log10() +
  theme_classic() +
  labs(x = "Min 4D effective sites per split (log10)", y = "Count")

p4 <- ggplot(dec_min, aes(x = w_min0D)) +
  geom_histogram(bins = 60) +
  scale_x_log10() +
  theme_classic() +
  labs(x = "Min 0D effective sites per split (log10)", y = "Count")

p1; p2; p3; p4




quantile(dec_min$n_genes_min, probs = c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE)
quantile(dec_min$w_min4D,     probs = c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE)
quantile(dec_min$w_min0D,     probs = c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE)




threshold_grid <- 1:30

retain_by_genes <- lapply(threshold_grid, function(mg) {
  d <- make_decoupled_dt_splitboth(gene_dt, min_genes_per_split = mg, seed = 1)
  tibble(
    min_genes_per_split = mg,
    n_groups = nrow(d),
    frac_retained = nrow(d) / nrow(summ_dt)  # if you have summ_dt; else divide by max
  )
}) %>% bind_rows()

retain_by_genes

ggplot(retain_by_genes, aes(x = min_genes_per_split, y = frac_retained)) +
  geom_line() + geom_point() +
  theme_classic() +
  labs(x = "min_genes_per_split", y = "Fraction of (MAG, colony, Month, Genus) retained")




sites_grid <- c(0, 50, 100, 200, 500, 1000, 2000, 5000)

retain_by_sites <- lapply(sites_grid, function(ms) {
  d <- make_decoupled_dt_splitboth(gene_dt,
                                   min_genes_per_split = 1,
                                   min_4D_sites_per_split = ms,
                                   min_0D_sites_per_split = 0,
                                   seed = 1)
  tibble(
    min_4D_sites_per_split = ms,
    n_groups = nrow(d),
    frac_retained = nrow(d) / nrow(summ_dt)
  )
}) %>% bind_rows()

retain_by_sites

ggplot(retain_by_sites, aes(x = min_4D_sites_per_split, y = frac_retained)) +
  geom_line() + geom_point() +
  scale_x_log10() +
  theme_classic() +
  labs(x = "min_4D_sites_per_split (log10)", y = "Fraction retained")






thresholds <- tidyr::crossing(
  min_genes_per_split = c(1, 2, 5),
  min_4D_sites_per_split = c(0, 200, 500, 1000, 2000)
)


R <- 50  # random split replicates per threshold combo

stab_grid <- thresholds %>%
  rowwise() %>%
  do({
    mg <- .$min_genes_per_split
    ms <- .$min_4D_sites_per_split
    
    est <- replicate(R, {
      d <- make_decoupled_dt_splitboth(
        gene_dt = gene_dt,
        min_genes_per_split = mg,
        min_4D_sites_per_split = ms,
        seed = NULL
      )
      f <- fit_good_model(d)
      if (is.null(f)) return(c(fd = NA, k = NA, n = nrow(d)))
      c(coef(f)[c("fd","k")], n = nrow(d))
    })
    
    est <- as.data.frame(t(est))
    est$min_genes_per_split <- mg
    est$min_4D_sites_per_split <- ms
    est
  }) %>%
  ungroup() %>%
  filter(is.finite(fd), is.finite(k))

stab_summary <- stab_grid %>%
  group_by(min_genes_per_split, min_4D_sites_per_split) %>%
  summarise(
    n_med = median(n),
    fd_med = median(fd),
    fd_q025 = quantile(fd, 0.025),
    fd_q975 = quantile(fd, 0.975),
    k_med = median(k),
    k_q025 = quantile(k, 0.025),
    k_q975 = quantile(k, 0.975),
    .groups = "drop"
  )

stab_summary


















# Choose sweet-spot thresholds
MG <- 2
MS4D <- 500

set.seed(1)
dec_final <- make_decoupled_dt_splitboth(
  gene_dt = gene_dt,
  min_genes_per_split = MG,
  min_4D_sites_per_split = MS4D,
  seed = 1
)

# Fit models on this final dataset
fit_good   <- fit_good_model(dec_final)
fit_h      <- fit_hill(dec_final)
gam_fit    <- fit_gam(dec_final)

coef(fit_good)
coef(fit_h)
summary(gam_fit)
AIC(fit_good, fit_h)

# log-spaced x grid for smooth curves on log axis
x_min <- min(dec_final$pi_syn_2, na.rm = TRUE)
x_max <- max(dec_final$pi_syn_2, na.rm = TRUE)
x_grid <- exp(seq(log(x_min), log(x_max), length.out = 400))
newdat <- data.frame(pi_syn_2 = x_grid)

pred_df <- tibble(
  pi_syn_2 = x_grid,
  `Good-style`     = predict(fit_good, newdata = newdat),
  `Hill/logistic`  = predict(fit_h,    newdata = newdat),
  `GAM`            = predict(gam_fit,  newdata = newdat, type = "response")
) %>%
  pivot_longer(-pi_syn_2, names_to = "model", values_to = "pi_ratio_pred") %>%
  filter(is.finite(pi_ratio_pred), pi_ratio_pred > 0)

ggplot(dec_final, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_line(data = pred_df, aes(y = pi_ratio_pred, linetype = model), linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = expression(pi[S]~"(split 2)"),
    y = expression(pi[N]/pi[S]~"(split 1)"),
    linetype = "Fit"
  )

dec_final

ggplot(dec_final, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_line(data = pred_df,
            aes(y = pi_ratio_pred, linetype = model, linewidth = model)) +
  scale_x_log10() +
  scale_y_log10() +
  scale_linewidth_manual(values = c(`Good-style` = 1.3, `Hill/logistic` = 0.9, `GAM` = 1.0)) +
  theme_classic() +
  labs(
    x = expression(pi[S]~"(split 2)"),
    y = expression(pi[N]/pi[S]~"(split 1)"),
    linetype = "Fit",
    linewidth = "Fit"
  )


ggplot(dec_final, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  geom_point(data = dec_final, aes(color = Genus), alpha = 0.7, size = 1) +
  geom_line(data = pred_df,
            aes(y = pi_ratio_pred, linetype = model, linewidth = model)) +
  scale_x_log10() +
  scale_y_log10() +
  scale_linewidth_manual(values = c(`Good-style` = 1.3, `Hill/logistic` = 0.9, `GAM` = 1.0)) +
  theme_classic() +
  labs(
    x = expression(pi[S]~"(split 2)"),
    y = expression(pi[N]/pi[S]~"(split 1)"),
    linetype = "Fit",
    linewidth = "Fit"
  )



x_grid <- exp(seq(log(1e-4), log(1e-1), length.out = 400))
newdat <- data.frame(pi_syn_2 = x_grid)

pred_df <- tibble(
  pi_syn_2 = x_grid,
  `Good-style`     = predict(fit_good, newdata = newdat),
  `Hill/logistic`  = predict(fit_h,    newdata = newdat),
  `GAM`            = predict(gam_fit,  newdata = newdat, type = "response")
) |>
  tidyr::pivot_longer(-pi_syn_2, names_to = "model", values_to = "pi_ratio_pred") |>
  dplyr::filter(is.finite(pi_ratio_pred), pi_ratio_pred > 0)



ggplot(dec_final, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_line(data = pred_df, aes(y = pi_ratio_pred, linetype = model), linewidth = 1) +
  scale_x_log10(
    breaks = 10^(-4:-1),
    labels = scales::label_math(10^.x)
  ) +
  scale_y_log10() +
  coord_cartesian(xlim = c(1e-4, 1e-1)) +
  theme_classic()+
  geom_hline(yintercept=1)



ggplot(dec_final, aes(x = pi_syn_2, y = pi_ratio_decoupled)) +
  geom_point(data = dec_final, aes(color = Genus), alpha = 0.6, size = 1.4) +
  geom_line(data = pred_df, aes(y = pi_ratio_pred, linetype = model), linewidth = 1) +
  scale_x_log10(
    breaks = 10^(-4:-1),
    labels = scales::label_math(10^.x)
  ) +
  scale_y_log10() +
  coord_cartesian(xlim = c(1e-4, 1e-1)) +
  theme_classic()+
  geom_hline(yintercept=1)










# quick checks
fit0 <- gam(log10(pi_ratio) ~ s(log10(pi_syn), k = 10), data = summ_dt, method="REML")
summary(fit0)
gam.check(fit0)
mgcv::k.check(fit0)   # warns if k too low


# Find extreme leverage / extreme values
summ_dt %>% 
  mutate(log_pi_syn = log10(pi_syn), log_pi_ratio = log10(pi_ratio)) %>%
  arrange(desc(abs(log_pi_ratio))) %>% head(20)  # extreme response












# target genera list you provided
target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium",
                   "commensalibacter", "frischella", "gilliamella",
                   "lactobacillus", "snodgrassella", "bombella", "bombilactobacillus")



# 1) Prepare data explicitly with log10 columns so there is NO ambiguity
df <- summ_dt %>%
  filter(is.finite(pi_syn), pi_syn > 0,
         is.finite(pi_nonsyn), pi_nonsyn > 0) %>%
  mutate(
    pi_ratio = pi_nonsyn / pi_syn,
    log_pi_syn   = log10(pi_syn),
    log_pi_ratio = log10(pi_ratio),
    Genus_plot = ifelse(Genus %in% target_genera, Genus, "other")
  )

# 2) Optional: inspect extremes (you already did)
top_extremes <- df %>% arrange(log_pi_ratio) %>% head(30)
print(top_extremes[, c("MAG_id","colony_id","Genus","pi_syn","pi_ratio","log_pi_ratio")])

# 3) Optional: create a trimmed / winsorized dataset to test sensitivity
#    (either winsorize top/bottom 0.5% or remove them)
low_thr  <- quantile(df$log_pi_ratio, 0.005, na.rm = TRUE)
high_thr <- quantile(df$log_pi_ratio, 0.995, na.rm = TRUE)
df_wins <- df %>% mutate(log_pi_ratio_w = pmin(pmax(log_pi_ratio, low_thr), high_thr))

# 4) Fit GAM on log–log scale (explicit column names: log_pi_ratio ~ s(log_pi_syn))
#    Use REML, select=TRUE, gamma>1 to reduce oversmoothing/overfitting
set.seed(1)
gam_fit <- mgcv::gam(
  log_pi_ratio ~ s(log_pi_syn, k = 25, bs = "tp"),
  data = df,
  method = "REML",
  select = TRUE,
  gamma = 1.4
)

summary(gam_fit)
mgcv::k.check(gam_fit)
gam.check(gam_fit)   # scan residuals and qqplot

# 5) Build a prediction grid INSIDE observed range (no extrapolation)
logx_min <- min(df$log_pi_syn, na.rm = TRUE)
logx_max <- max(df$log_pi_syn, na.rm = TRUE)
grid_logx <- seq(logx_min, logx_max, length.out = 400)

# Important: create newdata column with the same name used in model.frame()
pred_var_name <- names(model.frame(gam_fit))[2]   # should be "log_pi_syn"
newdata <- setNames(data.frame(grid_logx), pred_var_name)

pred <- predict(gam_fit, newdata = newdata, se.fit = TRUE)

pred_df <- data.frame(
  log_pi_syn = grid_logx,
  fit_log10  = pred$fit,
  se_log10   = pred$se.fit
) %>%
  mutate(
    pi_syn = 10^log_pi_syn,
    fit   = 10^fit_log10,
    upper = 10^(fit_log10 + 1.96 * se_log10),
    lower = 10^(fit_log10 - 1.96 * se_log10)
  )

# 6) Plot: points (colored only for target genera) + smooth + ribbon
cols <- setNames(RColorBrewer::brewer.pal(10,"Set3"), target_genera)
cols["other"] <- "grey80"

p <- ggplot() +
  geom_point(data = df, aes(x = pi_syn, y = pi_ratio, color = Genus_plot),
             size = 1.6, alpha = 0.5) +
  geom_ribbon(data = pred_df, aes(x = pi_syn, ymin = lower, ymax = upper),
              fill = "black", alpha = 0.12) +
  geom_line(data = pred_df, aes(x = pi_syn, y = fit),
            color = "black", size = 1) +
  scale_x_log10(labels = scales::label_math(10^.x)) +
  scale_y_log10(labels = scales::label_math(10^.x)) +
  scale_color_manual(values = cols) +
  labs(x = expression(pi[S]), y = expression(pi[N]/pi[S]),
       title = "GAM (log–log) fit: piN/piS ~ s(piS)",
       subtitle = "Model fit on log10 scale, predictions back-transformed; ribbon = 95% CI") +
  theme_minimal(base_size = 13)

print(p)
















# Colors
cols_targets <- RColorBrewer::brewer.pal(10, "Set3")
names(cols_targets) <- tolower(target_genera)
cols_targets <- c(cols_targets, other = "grey70")

# plotting params
grid_n <- 400          # prediction grid density (x)
pt_alpha <- 0.45
pt_size  <- 1.5
curve_lwd <- 1.1

# PREP DATA -----------------------------------------------------------------
# summ_dt must be your aggregated MAG×colony×Month table (pi_syn, pi_nonsyn, pi_ratio)
# decoupled_dt must be the result of make_decoupled_dt() with pi_syn_1, pi_syn_2, pi_ratio_decoupled, w
stopifnot(exists("summ_dt"), exists("decoupled_dt"))

# basic clean:
summ_df <- summ_dt %>%
  filter(is.finite(pi_syn), is.finite(pi_nonsyn), pi_syn > 0, pi_nonsyn > 0) %>%
  mutate(pi_ratio = pi_nonsyn / pi_syn,
         log_pi_syn = log10(pi_syn),
         log_pi_ratio = log10(pi_ratio),
         Genus_lc = tolower(Genus),
         Genus_plot = ifelse(Genus_lc %in% tolower(target_genera), Genus_lc, "other"))

dec_df <- decoupled_dt %>%
  filter(is.finite(pi_ratio_decoupled), is.finite(pi_syn_2), pi_syn_2 > 0) %>%
  mutate(log_pi_syn2 = log10(pi_syn_2),
         log_pi_ratio_dec = log10(pi_ratio_decoupled))

# compute genus medians (from summ_df)
genus_medians <- summ_df %>%
  group_by(Genus_plot) %>%
  summarize(
    median_pi_syn = 10^(median(log_pi_syn, na.rm = TRUE)),
    median_pi_ratio = 10^(median(log_pi_ratio, na.rm = TRUE)),
    n = n()
  ) %>% ungroup()

# prediction grid (on same x axis, use range from summ_df and dec_df)
x_min <- min(c(summ_df$pi_syn, dec_df$pi_syn_2), na.rm = TRUE)
x_max <- max(c(summ_df$pi_syn, dec_df$pi_syn_2), na.rm = TRUE)
grid_x <- 10^seq(log10(x_min) - 0.1, log10(x_max) + 0.1, length.out = grid_n) # slight pad
grid_logx <- log10(grid_x)

# GAM (fitted on summ_df) predictions (if not already fit, fit now)
if (!exists("gam_fit_summ") || is.null(gam_fit_summ)) {
  message("Fitting GAM on summ_df (log10 scale)...")
  gam_fit_summ <- tryCatch(
    mgcv::gam(log_pi_ratio ~ s(log_pi_syn, k = 25, bs = "tp"),
              data = summ_df, method = "REML", select = TRUE, gamma = 1.2),
    error = function(e) { warning("GAM fit failed: ", e$message); NULL }
  )
}
pred_gam <- NULL
if (!is.null(gam_fit_summ)) {
  pred_gam_raw <- tryCatch(predict(gam_fit_summ, newdata = data.frame(log_pi_syn = grid_logx), se.fit = TRUE),
                           error = function(e) NULL)
  if (!is.null(pred_gam_raw)) {
    pred_gam <- tibble::tibble(
      pi_syn = grid_x,
      fit_log10 = pred_gam_raw$fit,
      se = pred_gam_raw$se.fit
    ) %>% mutate(
      fit = 10^fit_log10,
      fit_u = 10^(fit_log10 + 1.96*se),
      fit_l = 10^(fit_log10 - 1.96*se)
    )
  }
}

# Good model predictions (unweighted) => using params fd, k if fit exists
pred_good <- NULL
if (exists("fit_good") && inherits(fit_good, "nls")) {
  co <- tryCatch(coef(fit_good), error = function(e) NULL)
  if (!is.null(co) && all(c("fd", "k") %in% names(co))) {
    fd <- as.numeric(co["fd"]); k <- as.numeric(co["k"])
    # formula: (1 - fd) + fd * ((1 - exp(-k * x)) / (k * x))
    pred_good <- tibble::tibble(
      pi_syn = grid_x,
      fit = (1 - fd) + fd * ((1 - exp(-k * grid_x)) / (k * grid_x))
    )
  }
}

# Weighted Good model predictions (fit_good_w)
pred_good_w <- NULL
if (exists("fit_good_w") && inherits(fit_good_w, "nls")) {
  co <- tryCatch(coef(fit_good_w), error = function(e) NULL)
  if (!is.null(co) && all(c("fd", "k") %in% names(co))) {
    fd <- as.numeric(co["fd"]); k <- as.numeric(co["k"])
    pred_good_w <- tibble::tibble(
      pi_syn = grid_x,
      fit = (1 - fd) + fd * ((1 - exp(-k * grid_x)) / (k * grid_x))
    )
  }
}

# Hill fit predictions (if fit_h exists; params c, x0, h)
pred_hill <- NULL
if (exists("fit_h") && inherits(fit_h, "nls")) {
  co <- tryCatch(coef(fit_h), error = function(e) NULL)
  if (!is.null(co) && all(c("c", "x0", "h") %in% names(co))) {
    c0 <- as.numeric(co["c"]); x0 <- as.numeric(co["x0"]); h <- as.numeric(co["h"])
    pred_hill <- tibble::tibble(
      pi_syn = grid_x,
      fit = c0 + (1 - c0) / (1 + (grid_x / x0)^h)
    )
  }
}

# GAM on decoupled data (pi_ratio_decoupled ~ s(pi_syn_2))
pred_gam_dec <- NULL
if (exists("gam_fit") && inherits(gam_fit, "gam")) {
  # NOTE: this gam_fit may refer to earlier GAMs — if you fitted a GAM on decoupled_dt call it gam_fit_dec
  # If you didn't, fit a quick one here:
  if (!inherits(gam_fit, "gam") || (exists("gam_fit") && !inherits(gam_fit, "gam"))) {
    gam_fit_dec <- tryCatch(mgcv::gam(log10(pi_ratio_decoupled) ~ s(log10(pi_syn_2), k = 15),
                                      data = dec_df, method = "REML", gamma = 1.2),
                            error = function(e) NULL)
  } else {
    gam_fit_dec <- gam_fit
  }
  if (!is.null(gam_fit_dec)) {
    pr <- tryCatch(predict(gam_fit_dec, newdata = data.frame('log10(pi_syn_2)' = grid_logx), se.fit = TRUE), error = function(e) NULL)
    if (!is.null(pr)) {
      pred_gam_dec <- tibble::tibble(
        pi_syn = grid_x,
        fit_log10 = pr$fit,
        se = pr$se.fit
      ) %>% mutate(
        fit = 10^fit_log10,
        fit_u = 10^(fit_log10 + 1.96*se),
        fit_l = 10^(fit_log10 - 1.96*se)
      )
    }
  }
}

# PLOTTING ------------------------------------------------------------------
p <- ggplot() +
  # background points from summ_df (MAG snapshots)
  geom_point(data = summ_df, aes(x = pi_syn, y = pi_ratio, color = Genus_plot),
             alpha = pt_alpha, size = pt_size) +
  scale_x_log10(labels = label_math(10^.x), limits = c(x_min*0.9, x_max*1.1)) +
  scale_y_log10(labels = label_math(10^.x)) +
  scale_color_manual(values = cols_targets, guide = guide_legend(order = 1)) +
  labs(
    x = expression(pi[S]),
    y = expression(pi[N]/pi[S]),
    title = "Combined model overlay: GAM (summ), Good, weighted-Good, Hill, GAM(decoupled)",
    subtitle = "GAM (black) = fitted on summ_dt; other model curves colored; highlighted genus medians = filled circles with dark border"
  ) +
  theme_minimal(base_size = 12)

# add GAM (summ) ribbon & line
if (!is.null(pred_gam)) {
  p <- p +
    geom_ribbon(data = pred_gam, aes(x = pi_syn, ymin = fit_l, ymax = fit_u),
                inherit.aes = FALSE, fill = "black", alpha = 0.12) +
    geom_line(data = pred_gam, aes(x = pi_syn, y = fit), color = "black", size = curve_lwd)
}

# add GAM(decoupled) if available (thin dashed)
if (!is.null(pred_gam_dec)) {
  p <- p +
    geom_ribbon(data = pred_gam_dec, aes(x = pi_syn, ymin = fit_l, ymax = fit_u),
                inherit.aes = FALSE, fill = "black", alpha = 0.06) +
    geom_line(data = pred_gam_dec, aes(x = pi_syn, y = fit), color = "black", size = 0.7, linetype = "dashed")
}

# add Good (unweighted) curve
if (!is.null(pred_good)) {
  p <- p + geom_line(data = pred_good, aes(x = pi_syn, y = fit), color = "#D73027", size = curve_lwd, linetype = "solid")
}

# add Good (weighted) curve
if (!is.null(pred_good_w)) {
  p <- p + geom_line(data = pred_good_w, aes(x = pi_syn, y = fit), color = "#FC8D59", size = curve_lwd, linetype = "dotdash")
}

# add Hill curve
if (!is.null(pred_hill)) {
  p <- p + geom_line(data = pred_hill, aes(x = pi_syn, y = fit), color = "#4575B4", size = curve_lwd, linetype = "longdash")
}

# add genus medians for target genera (filled circle, black stroke)
p <- p +
  geom_point(data = genus_medians %>% filter(Genus_plot != "other"),
             aes(x = median_pi_syn, y = median_pi_ratio, fill = Genus_plot),
             shape = 21, color = "black", size = 3.2, stroke = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = cols_targets, guide = FALSE)

# add a vertical rug so you can see density of x
p <- p + geom_rug(data = summ_df, aes(x = pi_syn), inherit.aes = FALSE, sides = "b", alpha = 0.15)

# final theme tweaks
p <- p +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p)

# Save figure
ggsave("combined_piNpiS_models_overlay.png", p, width = 10, height = 9, dpi = 300)
min(summ_df$pi_syn)
