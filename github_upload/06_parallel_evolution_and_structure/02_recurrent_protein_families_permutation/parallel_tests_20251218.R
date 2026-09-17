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
library(ggrepel)
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










#------------------ Season map -----------------------------
  
  season_map <- c("May"="Spring","June"="Spring",
                  "July"="Summer","August"="Summer","September"="Summer",
                  "October"="Fall","November"="Fall",
                  "January"="Winter","February"="Winter")

## ------------------ Helpers for annotation strings --------

safe_unlist <- function(x){
  if (is.null(x)) return(NA_character_)
  y <- unlist(x)
  if (length(y) == 0) return(NA_character_)
  paste(unique(as.character(y)), collapse = ";;")
}

normalize_anno_string_one <- function(x) {
  # x should be length-1, may be NA/empty
  if (length(x) == 0L || is.na(x)) return(NA_character_)
  s <- trimws(as.character(x))
  if (s == "" || s == "<NA>") return(NA_character_)
  
  # Python-style ['a', 'b']
  s <- gsub("^\\[|\\]$", "", s)
  s <- gsub("',\\s*'", ";;", s)
  s <- gsub("'", "", s)
  
  # R-style c("a", "b")
  s <- gsub("^c\\(", "", s)
  s <- gsub("\\)$", "", s)
  s <- gsub('"', "", s)
  s <- gsub(",\\s*", ";;", s)
  
  s
}

parse_anno_terms_one <- function(x,
                                 drop_wrappers = TRUE,
                                 drop_go_ids   = TRUE) {
  s <- normalize_anno_string_one(x)
  if (is.na(s)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  parts <- strsplit(s, ";;", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  if (!length(parts)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  if (drop_wrappers) {
    wrappers <- c("molecular_function",
                  "biological_process",
                  "cellular_component")
    parts <- parts[!parts %in% wrappers]
  }
  
  if (drop_go_ids) {
    parts <- parts[
      !grepl("^GO:\\d{7}$", parts) &
        !grepl("^\\d{7}$",   parts)
    ]
  }
  
  if (!length(parts)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  list(primary = parts[1L], all_terms = parts)
}

## ------------------ Load base tables -----------------------

message("Loading gene summary, taxonomy, and annotations...")

gene_dt  <- fread(GENE_SUM_FILE)
tax_dt   <- fread(TAX_FILE)
ann_dt   <- readRDS(ANN_RDS)

# ensure data.table
gene_dt <- as.data.table(gene_dt)
if (!is.data.table(ann_dt)) ann_dt <- as.data.table(ann_dt)

# canonicalize ID types
for (c in c("MAG_id","colony_id","corresponding_gene_call","Month")) {
  if (c %in% names(gene_dt)) gene_dt[[c]] <- as.character(gene_dt[[c]])
}
tax_dt[, MAG_id := as.character(MAG_id)]

## ------------------ Taxonomy: attach Genus -----------------

tax_dt_unique <- unique(tax_dt[, .(MAG_id, Genus)])
gene_dt <- merge(
  gene_dt,
  tax_dt_unique,
  by = "MAG_id",
  all.x = TRUE
)

## ------------------ Collapse annotation strings ------------

if ("corresponding_gene_call" %in% names(ann_dt)) {
  ann_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]
}

# If GO_str / pfam_str not precomputed, build from list columns
if (!"GO_str" %in% names(ann_dt) && "GO_slim_name_list" %in% names(ann_dt)) {
  ann_dt[, GO_str := sapply(GO_slim_name_list, safe_unlist)]
}
if (!"pfam_str" %in% names(ann_dt) && "pfam_label_list" %in% names(ann_dt)) {
  ann_dt[, pfam_str := sapply(pfam_label_list, safe_unlist)]
}

ann_sub <- ann_dt[, .(corresponding_gene_call, GO_str, pfam_str)]

gene_dt <- merge(
  gene_dt,
  ann_sub,
  by = "corresponding_gene_call",
  all.x = TRUE
)

## ------------------ Parse GO / PFAM primary terms ----------

gene_dt[, c("GO_primary", "GO_terms") := {
  parsed <- lapply(GO_str, parse_anno_terms_one)
  list(
    vapply(parsed, `[[`, character(1), "primary"),
    lapply(parsed, `[[`, "all_terms")
  )
}]

gene_dt[, c("PFAM_primary", "PFAM_terms") := {
  parsed <- lapply(
    pfam_str,
    parse_anno_terms_one,
    drop_wrappers = FALSE,  # PFAM doesn't have wrapper terms
    drop_go_ids   = FALSE
  )
  list(
    vapply(parsed, `[[`, character(1), "primary"),
    lapply(parsed, `[[`, "all_terms")
  )
}]

## ------------------ Add Season, Colony, pi_ratio -----------

if ("Month" %in% names(gene_dt)) {
  gene_dt[, Season := season_map[Month]]
}

# convenience: Colony as char copy of colony_id
if ("colony_id" %in% names(gene_dt)) {
  gene_dt[, Colony := as.character(colony_id)]
}

# pi_ratio = pi_nonsyn / pi_syn, avoid weird infinities
if (all(c("pi_nonsyn","pi_syn") %in% names(gene_dt))) {
  gene_dt[, pi_ratio := pi_nonsyn / pi_syn]
  gene_dt[!is.finite(pi_ratio), pi_ratio := NA_real_]
}

## ------------------ Optional cleaning / filtering ----------

# Drop rows with missing MAG or colony if you want
gene_clean <- gene_dt[!is.na(MAG_id) & !is.na(Colony)]

# Set a key for fast joins later
setkey(gene_clean, corresponding_gene_call, MAG_id, Colony)

message("gene_clean rebuilt: ", nrow(gene_clean), " rows; columns: ",
        paste(names(gene_clean), collapse = ", "))



gene_clean



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

sweep_dt




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

par_zeroA <- parallel_zeroD_A
par_zeroA_expanded <- par_zeroA[, .(corresponding_gene_call = corresponding_gene_call,
                                    n_MAG, n_Genus, n_Col,
                                    MAGs = sapply(MAGs, function(x) paste(x, collapse=";")),
                                    Genus_list = sapply(Genus_list, function(x) paste(x, collapse=";")))]
# add annotation cols
par_zeroA_expanded <- merge(par_zeroA_expanded, unique(gene_dt[, .(corresponding_gene_call, GO_primary, PFAM_primary, Genus)]),
                            by = "corresponding_gene_call", all.x = TRUE)




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
# Instead we compute slice_ymin/ymax as: y - cell_height/2 + (y - y)?? simple# Filepaths (edit if needed)
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




normalize_anno_string_one <- function(x) {
  # x should be length-1, may be NA/empty
  if (length(x) == 0L || is.na(x)) return(NA_character_)
  s <- trimws(as.character(x))
  if (s == "" || s == "<NA>") return(NA_character_)
  
  # --- Handle Python-style list: ['a', 'b'] ---
  s <- gsub("^\\[|\\]$", "", s)
  s <- gsub("',\\s*'", ";;", s)
  s <- gsub("'", "", s)
  
  # --- Handle R c("a", "b") style ---
  s <- gsub("^c\\(", "", s)
  s <- gsub("\\)$", "", s)
  s <- gsub('"', "", s)
  s <- gsub(",\\s*", ";;", s)
  
  s
}

parse_anno_terms_one <- function(x,
                                 drop_wrappers = TRUE,
                                 drop_go_ids   = TRUE) {
  s <- normalize_anno_string_one(x)
  if (is.na(s)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  parts <- strsplit(s, ";;", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  if (!length(parts)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  if (drop_wrappers) {
    wrappers <- c("molecular_function",
                  "biological_process",
                  "cellular_component")
    parts <- parts[!parts %in% wrappers]
  }
  
  if (drop_go_ids) {
    parts <- parts[
      !grepl("^GO:\\d{7}$", parts) &
        !grepl("^\\d{7}$",   parts)
    ]
  }
  
  if (!length(parts)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  list(primary = parts[1L], all_terms = parts)
}



gene_dt <- as.data.table(gene_dt)

gene_dt[, c("GO_primary", "GO_terms") := {
  parsed <- lapply(GO_str, parse_anno_terms_one)
  list(
    vapply(parsed, `[[`, character(1), "primary"),
    lapply(parsed, `[[`, "all_terms")
  )
}]

gene_dt[, c("PFAM_primary", "PFAM_terms") := {
  parsed <- lapply(
    pfam_str,
    parse_anno_terms_one,
    drop_wrappers = FALSE,
    drop_go_ids   = FALSE
  )
  list(
    vapply(parsed, `[[`, character(1), "primary"),
    lapply(parsed, `[[`, "all_terms")
  )
}]


gene_dt[!is.na(GO_primary)][1:20,
                            .(corresponding_gene_call, GO_str, GO_primary)]

gene_dt[!is.na(PFAM_primary)][1:20,
                              .(corresponding_gene_call, pfam_str, PFAM_primary)]
r to map 0..1 to [y-cellh/2, y+cellh/2]
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




normalize_anno_string_one <- function(x) {
  # x should be length-1, may be NA/empty
  if (length(x) == 0L || is.na(x)) return(NA_character_)
  s <- trimws(as.character(x))
  if (s == "" || s == "<NA>") return(NA_character_)
  
  # --- Handle Python-style list: ['a', 'b'] ---
  s <- gsub("^\\[|\\]$", "", s)
  s <- gsub("',\\s*'", ";;", s)
  s <- gsub("'", "", s)
  
  # --- Handle R c("a", "b") style ---
  s <- gsub("^c\\(", "", s)
  s <- gsub("\\)$", "", s)
  s <- gsub('"', "", s)
  s <- gsub(",\\s*", ";;", s)
  
  s
}

parse_anno_terms_one <- function(x,
                                 drop_wrappers = TRUE,
                                 drop_go_ids   = TRUE) {
  s <- normalize_anno_string_one(x)
  if (is.na(s)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  parts <- strsplit(s, ";;", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  if (!length(parts)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  if (drop_wrappers) {
    wrappers <- c("molecular_function",
                  "biological_process",
                  "cellular_component")
    parts <- parts[!parts %in% wrappers]
  }
  
  if (drop_go_ids) {
    parts <- parts[
      !grepl("^GO:\\d{7}$", parts) &
        !grepl("^\\d{7}$",   parts)
    ]
  }
  
  if (!length(parts)) {
    return(list(primary = NA_character_, all_terms = character(0)))
  }
  
  list(primary = parts[1L], all_terms = parts)
}



gene_dt <- as.data.table(gene_dt)

gene_dt[, c("GO_primary", "GO_terms") := {
  parsed <- lapply(GO_str, parse_anno_terms_one)
  list(
    vapply(parsed, `[[`, character(1), "primary"),
    lapply(parsed, `[[`, "all_terms")
  )
}]

gene_dt[, c("PFAM_primary", "PFAM_terms") := {
  parsed <- lapply(
    pfam_str,
    parse_anno_terms_one,
    drop_wrappers = FALSE,
    drop_go_ids   = FALSE
  )
  list(
    vapply(parsed, `[[`, character(1), "primary"),
    lapply(parsed, `[[`, "all_terms")
  )
}]


gene_dt[!is.na(GO_primary)][1:20,
                            .(corresponding_gene_call, GO_str, GO_primary)]

gene_dt[!is.na(PFAM_primary)][1:20,
                              .(corresponding_gene_call, pfam_str, PFAM_primary)]





## --- Genus-level means per GO term ----------------------------------------

pathway_genus <- go_long %>%
  # keep only pathways that passed your n_obs / n_genes filters
  semi_join(pathway_perm_plot, by = "GO_term") %>%
  group_by(GO_term, Genus) %>%
  summarise(
    mean_pi_ratio = mean(pi_ratio, na.rm = TRUE),
    n_genes       = n_distinct(corresponding_gene_call),
    n_obs         = n(),
    .groups       = "drop"
  ) %>%
  # optional: drop tiny Genus–pathway combos
  filter(n_genes >= 5)

## --- Use a single GO-term ordering for both layers -----------------------


go_order <- pathway_perm_plot %>%
  arrange(mean_pi_ratio) %>%
  pull(GO_term)

pathway_perm_plot <- pathway_perm_plot %>%
  mutate(GO_term = factor(GO_term, levels = go_order))

pathway_genus <- pathway_genus %>%
  mutate(GO_term = factor(GO_term, levels = go_order))




## MAG-level table: create a color group
pathway_mag <- pathway_mag %>%
  mutate(
    Genus_plot = if_else(
      Genus %in% target_genera,
      Genus,
      "other"
    ),
    Genus_plot = factor(
      Genus_plot,
      levels = c(target_genera, "other")
    )
  )




## typical null CI as before
global_ci <- c(
  median(pathway_perm_plot$null_ci_low,  na.rm = TRUE),
  median(pathway_perm_plot$null_ci_high, na.rm = TRUE)
)

## make sure GO_term is ordered by ascending mean piN/piS
go_order <- pathway_perm_plot %>%
  arrange(mean_pi_ratio) %>%
  pull(GO_term) %>%
  unique()

pathway_perm_plot <- pathway_perm_plot %>%
  mutate(GO_term = factor(GO_term, levels = go_order))

pathway_mag <- pathway_mag %>%
  mutate(GO_term = factor(GO_term, levels = go_order))


# color palette for highlighted genera
cols_targets <- c(
  "apilactobacillus" = "#1b9e77",
  "bartonella"       = "#d95f02",
  "bifidobacterium"  = "#7570b3",
  "commensalibacter" = "#e7298a",
  "frischella"       = "#66a61e",
  "gilliamella"      = "#e6ab02",
  "lactobacillus"    = "#a6761d",
  "snodgrassella"    = "#1f78b4",
  "other"            = "grey80"
)







## --- 1. Add sig_group & compute "typical" null CI (Option 2) ----

pathway_perm_plot <- pathway_perm %>%
  mutate(
    sig_group = case_when(
      q_high < 0.05 ~ "higher-than-null (weaker purifying)",
      q_low  < 0.05 ~ "lower-than-null (stronger purifying)",
      TRUE         ~ "not-significant"
    )
  ) %>%
  filter(
    n_MAGs >= 5,
    !str_detect(
      GO_term,
      regex("ribosome|organelle|term tracker item", ignore_case = TRUE)
    )
  )
# restrict MAG×pathway table to the same filtered set
keep_go <- unique(pathway_perm_plot$GO_term)

pathway_mag <- pathway_mag %>%
  filter(GO_term %in% keep_go)

## --- Recompute null CI (if using the median-of-null-intervals approach) ---

global_ci <- c(
  median(pathway_perm_plot$null_ci_low,  na.rm = TRUE),
  median(pathway_perm_plot$null_ci_high, na.rm = TRUE)
)


go_order <- pathway_perm_plot %>%
  arrange(mean_pi_ratio) %>%
  pull(GO_term) %>%
  unique()

pathway_perm_plot <- pathway_perm_plot %>%
  mutate(GO_term = factor(GO_term, levels = go_order))

pathway_mag <- pathway_mag %>%
  mutate(GO_term = factor(GO_term, levels = go_order))


## --- 4. Y limits so everything fits ----

y_min <- min(
  pathway_mag$mag_mean_pi_ratio,
  pathway_perm_plot$mean_pi_ratio,
  global_ci,
  na.rm = TRUE
)
y_max <- max(
  pathway_mag$mag_mean_pi_ratio,
  pathway_perm_plot$mean_pi_ratio,
  global_ci,
  na.rm = TRUE
)


## --- 5. Final plot: fixed-size species-mean points ----

p_pathway <- ggplot() +
  # MAG-level points (small, colored by Genus)
  geom_point(
    data = pathway_mag,
    aes(
      x     = mag_mean_pi_ratio,
      y     = GO_term,
      color = Genus_plot
    ),
    size  = 1,
    alpha = 0.6
  ) +
  
  # vertical reference lines (null CI + piN/piS = 1)
  geom_vline(
    xintercept = global_ci,
    linetype   = "dashed",
    linewidth  = 0.4
  ) +
  geom_vline(
    xintercept = 1,
    linetype   = "dotted",
    linewidth  = 1
  ) +
  
  # pathway mean across all MAGs (large, black; shape = significance)
  geom_point(
    data = pathway_perm_plot,
    aes(
      x     = mean_pi_ratio,
      y     = GO_term,
      shape = sig_group
    ),
    size   = 2,
    color  = "black",
    fill = "black",
    stroke = 1,
    shape = 21
  ) +
  
  scale_x_log10(
   ) +
  
  scale_color_manual(
    values = cols_targets,
    breaks = target_genera,         # legend only for highlighted genera
    name   = "Genus (highlighted)"
  ) +

  
  labs(
    x = expression(pi[N]/pi[S]),
    y = "GO term (pathway)",
    title = "Pathway-level nonsynonymous diversity relative to purifying-selection null",
    subtitle = paste(
      "Small points = MAG × pathway (colored by highlighted genera; others grey);",
      "large black points = pathway means across all MAGs (shape = deviation from null);",
      "dashed = typical null CI; dotted = π[N]/π[S] = 1"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 8, angle = 30),
    axis.text.x        = element_text(angle = 45, hjust = 0.5)
  )
p_pathway






## Tunable parameters -------------------------------------------------
xmin      <- 1e-3   # lower end of data
xmax      <- 2      # upper end of data
mid_low   <- 0.1    # start of zoom region
mid_high  <- 1      # end of zoom region

left_frac  <- 0.10  # fraction of axis width for < mid_low
mid_frac   <- 0.80  # fraction for [mid_low, mid_high]
right_frac <- 0.10  # fraction for > mid_high

## Precompute logs
l_min      <- log10(xmin)
l_max      <- log10(xmax)
l_mid_low  <- log10(mid_low)
l_mid_high <- log10(mid_high)

## Forward transform: x -> squished coordinate ------------------------
squish_log_forward <- function(x) {
  out <- rep(NA_real_, length(x))        # preserve NA positions
  
  ok <- !is.na(x) & x > 0                # finite, positive values only
  if (!any(ok)) return(out)
  
  lx <- log10(x[ok])
  res <- numeric(length(lx))
  
  left_idx  <- lx <  l_mid_low
  mid_idx   <- lx >= l_mid_low  & lx <= l_mid_high
  right_idx <- lx >  l_mid_high
  
  # left part
  if (any(left_idx)) {
    res[left_idx] <- (lx[left_idx] - l_min) /
      (l_mid_low - l_min) * left_frac
  }
  
  # middle (zoomed) part
  if (any(mid_idx)) {
    res[mid_idx] <- left_frac +
      (lx[mid_idx] - l_mid_low) /
      (l_mid_high - l_mid_low) * mid_frac
  }
  
  # right part
  if (any(right_idx)) {
    res[right_idx] <- left_frac + mid_frac +
      (lx[right_idx] - l_mid_high) /
      (l_max - l_mid_high) * right_frac
  }
  
  out[ok] <- res
  out
}

## Inverse transform: squished coordinate -> x ------------------------
squish_log_inverse <- function(y) {
  out <- rep(NA_real_, length(y))        # preserve NA positions
  
  ok <- !is.na(y)
  if (!any(ok)) return(out)
  
  yy <- y[ok]
  res <- numeric(length(yy))
  
  left_idx  <- yy < left_frac
  mid_idx   <- yy >= left_frac & yy <= (left_frac + mid_frac)
  right_idx <- yy > (left_frac + mid_frac)
  
  # left
  if (any(left_idx)) {
    res[left_idx] <- 10^(
      l_min + (yy[left_idx] / left_frac) * (l_mid_low - l_min)
    )
  }
  
  # middle
  if (any(mid_idx)) {
    res[mid_idx] <- 10^(
      l_mid_low +
        ((yy[mid_idx] - left_frac) / mid_frac) * (l_mid_high - l_mid_low)
    )
  }
  
  # right
  if (any(right_idx)) {
    res[right_idx] <- 10^(
      l_mid_high +
        ((yy[right_idx] - left_frac - mid_frac) / right_frac) *
        (l_max - l_mid_high)
    )
  }
  
  out[ok] <- res
  out
}

squish_log_trans <- trans_new(
  name      = "squish_log",
  transform = squish_log_forward,
  inverse   = squish_log_inverse
)





p_pathway <- ggplot() +
  # MAG-level points (small, colored by Genus)
  geom_jitter(
    data = pathway_mag,
    aes(
      x     = mag_mean_pi_ratio,
      y     = GO_term,
      color = Genus_plot
    ),
    size  = 1.5,
    alpha = 0.4,
    height=0.1
  ) +
  
  # vertical reference lines (null CI + piN/piS = 1)
  geom_vline(
    xintercept = global_ci,
    linetype   = "dashed",
    linewidth  = 0.4
  ) +
  geom_vline(
    xintercept = 1,
    linetype   = "dotted",
    linewidth  = 1
  ) +
  
  # pathway mean across all MAGs (large, black; shape = significance)
  geom_point(
    data = pathway_perm_plot,
    aes(
      x     = mean_pi_ratio,
      y     = GO_term,
      shape = sig_group
    ),
    size   = 2,
    color  = "black",
    fill = "black",
    stroke = 1,
    shape = 21
  ) +
  
  scale_x_continuous(
    trans  = squish_log_trans,
    breaks = c(1e-1, 1, 10),
    limits = c(xmin, xmax)
  ) +
  
  scale_color_manual(
    values = cols_targets,
    breaks = target_genera,         # legend only for highlighted genera
    name   = "Genus (highlighted)"
  ) +
  
  
  labs(
    x = expression(pi[N]/pi[S]),
    y = "GO term (pathway)",
    title = "Pathway-level nonsynonymous diversity relative to purifying-selection null",
    subtitle = paste(
      "Small points = MAG × pathway (colored by highlighted genera; others grey);",
      "large black points = pathway means across all MAGs (shape = deviation from null);",
      "dashed = typical null CI; dotted = π[N]/π[S] = 1"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 8, angle = 30),
    axis.text.x        = element_text(angle = 45, hjust = 0.5)
  )
p_pathway












pathway_order <- pathway_perm_plot %>%
  arrange(mean_pi_ratio) %>%
  pull(GO_term)

pathway_perm_plot <- pathway_perm_plot %>%
  mutate(GO_term = factor(GO_term, levels = pathway_order))

pathway_mag <- pathway_mag %>%
  mutate(
    GO_term = factor(GO_term, levels = pathway_order),
    # highlight only target genera, collapse others to "other"
    Genus_plot = ifelse(Genus %in% target_genera, Genus, "other")
  )

# pathway_genus_means should contain: GO_term, Genus, mean_pi_ratio_genus
# ensure matching factor levels and collapsed colour category
pathway_genus_means <- pathway_genus_means %>%
  mutate(
    GO_term = factor(GO_term, levels = pathway_order),
    Genus_plot = ifelse(Genus %in% target_genera, Genus, "other")
  )

# ensure the named colours include "other"
if (!"other" %in% names(cols_targets)) {
  cols_targets <- c(cols_targets, other = "grey80")
}

# ---- plotting ----
p_pathway <- ggplot() +
  
  # MAG-level points (small, colored by Genus if in target_genera else grey)
  geom_jitter(
    data = pathway_mag,
    aes(
      x     = mag_mean_pi_ratio,
      y     = GO_term,
      color = Genus_plot
    ),
    size  = 1.2,
    alpha = 0.35,
    height = 0.10,
    width  = 0
  ) +
  
  # NEW: Genus-level means (triangles) — one point per Genus × pathway
  geom_point(
    data = pathway_genus_means,
    aes(
      x     = mean_pi_ratio_genus,
      y     = GO_term,
      fill = Genus_plot,
    ),
    shape = 21, 
    size  = 2,
    alpha = 0.95,
    stroke = 0.8,
  ) +
  
  # vertical reference lines (null CI + piN/piS = 1)
  geom_vline(
    xintercept = global_ci,
    linetype   = "dashed",
    linewidth  = 0.4
  ) +
  geom_vline(
    xintercept = 1,
    linetype   = "dotted",
    linewidth  = 0.8
  ) +
  
  # pathway mean across all MAGs (large, black; shape = significance)
  geom_point(
    data = pathway_perm_plot,
    aes(
      x     = mean_pi_ratio,
      y     = GO_term,
      shape = sig_group
    ),
    size   = 2.5,
    color  = "black",
    fill   = "black",
    stroke = 0.8
  ) +
  
  # X axis: squished log with explicit breaks at 1e-3,1e-2,1e-1,1
  # Use your squish transform; if not defined, replace with scale_x_log10()
  scale_x_continuous(
    trans  = squish_log_trans,                    # <- assumes squish_log_trans already defined
    breaks = c(1e-3, 1e-2, 1e-1, 1),
    labels = c("10^-3", "10^-2", "10^-1", "10^0"),
    limits = c(xmin, xmax)
  ) +
  
  scale_color_manual(
    values = cols_targets,
    breaks = c(target_genera, "other") %>% unique(),
    name   = "Genus (highlighted)"
  ) +
  
  scale_shape_manual(
    values = c(
      "higher-than-null (weaker purifying)" = 21, 
      "not-significant"                     = 21,
      "lower-than-null (stronger purifying)"= 21
    ),
    guide = guide_legend(override.aes = list(fill = "black")),
    name = NULL
  ) +
  
  labs(
    x = expression(pi[N]/pi[S]),
    y = "GO term (pathway)",
    title = "Pathway-level nonsynonymous diversity relative to purifying-selection null",
    subtitle = paste(
      "Small points = MAG × pathway (colored for highlighted genera; others grey);",
      "triangles = Genus means (≥5 MAGs);",
      "large black points = pathway means across all MAGs (shape = deviation from null);",
      "dashed = typical null CI; dotted = π[N]/π[S] = 1"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor.y = element_blank(),
    axis.text.y        = element_text(size = 9,angle=30, face = "bold", family = "sans"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    legend.position    = "right"
  )

p_pathway





























# ---------- user params (adjust) ----------
target_genera <- c("apilactobacillus","bartonella","bifidobacterium",
                   "commensalibacter","frischella","gilliamella",
                   "lactobacillus","snodgrassella")

# mid-range we want to emphasize
mid_low  <- 1e-1
mid_high <- 1
xmin_allowed <- 1e-3
xmax_allowed <- 10
mid_frac <- 0.80

# ensure global_ci exists (from your earlier code)
if (!exists("global_ci")) stop("please create 'global_ci' before running plotting block")

# ---------- corrected squish transform ----------
squish_log_trans <- function(lmin = xmin_allowed,
                             mid_low = mid_low,
                             mid_high = mid_high,
                             rmax = xmax_allowed,
                             mid_frac = mid_frac) {
  
  # checks
  if (!(lmin > 0 && mid_low > lmin && mid_high > mid_low && rmax > mid_high)) {
    stop("Invalid bounds for squish_log_trans.")
  }
  
  L0 <- log10(lmin); L1 <- log10(mid_low); L2 <- log10(mid_high); L3 <- log10(rmax)
  left_frac  <- (1 - mid_frac) / 2
  right_frac <- left_frac
  
  # forward transform: x -> y in [0,1]
  fwd <- function(x) {
    out <- rep(NA_real_, length(x))
    ok <- !is.na(x) & is.finite(x) & x > 0
    if (!any(ok)) return(out)
    xv <- x[ok]
    lx <- log10(xv)
    
    yv <- numeric(length(xv))
    
    # masks & indices
    left_idx  <- which(lx <= L1)
    mid_idx   <- which(lx > L1  & lx <= L2)
    right_idx <- which(lx > L2)
    
    if (length(left_idx) > 0) {
      if (L1 > L0) {
        prop <- (lx[left_idx] - L0) / (L1 - L0)
        yv[left_idx] <- 0 + prop * left_frac
      } else {
        yv[left_idx] <- left_frac / 2
      }
    }
    
    if (length(mid_idx) > 0) {
      if (L2 > L1) {
        prop <- (lx[mid_idx] - L1) / (L2 - L1)
        yv[mid_idx] <- left_frac + prop * mid_frac
      } else {
        yv[mid_idx] <- left_frac + mid_frac/2
      }
    }
    
    if (length(right_idx) > 0) {
      if (L3 > L2) {
        prop <- (lx[right_idx] - L2) / (L3 - L2)
        yv[right_idx] <- left_frac + mid_frac + prop * right_frac
      } else {
        yv[right_idx] <- left_frac + mid_frac + right_frac/2
      }
    }
    
    out[ok] <- yv
    out
  }
  
  # inverse: y -> x
  inv <- function(y) {
    out <- rep(NA_real_, length(y))
    ok <- !is.na(y) & is.finite(y)
    if (!any(ok)) return(out)
    yv <- y[ok]
    
    left_hi <- left_frac
    mid_hi  <- left_frac + mid_frac
    
    xv <- numeric(length(yv))
    
    left_idx  <- which(yv <= left_hi)
    mid_idx   <- which(yv > left_hi & yv <= mid_hi)
    right_idx <- which(yv > mid_hi)
    
    if (length(left_idx) > 0) {
      prop <- (yv[left_idx] - 0) / left_frac
      xv[left_idx] <- 10^(L0 + prop * (L1 - L0))
    }
    if (length(mid_idx) > 0) {
      prop <- (yv[mid_idx] - left_frac) / mid_frac
      xv[mid_idx] <- 10^(L1 + prop * (L2 - L1))
    }
    if (length(right_idx) > 0) {
      prop <- (yv[right_idx] - (left_frac + mid_frac)) / right_frac
      xv[right_idx] <- 10^(L2 + prop * (L3 - L2))
    }
    
    out[ok] <- xv
    out
  }
  
  trans_new(
    name = paste0("squish_log_", formatC(mid_frac, digits=2)),
    transform = fwd,
    inverse   = inv,
    domain = c(lmin, rmax)
  )
}


# ---------- prepare plot data ----------
# mag_perm_filt should contain MAG_id, Genus, mean_pi_ratio
plot_dt <- mag_perm_filt %>%
  filter(!is.na(mean_pi_ratio)) %>%
  mutate(Genus_plot = ifelse(tolower(Genus) %in% tolower(target_genera), Genus, "other"))

# color vector (targets colored, others grey)
target_palette <- RColorBrewer::brewer.pal(max(3, min(length(target_genera),8)), "Set1")
cols_targets <- setNames(target_palette[1:length(target_genera)], target_genera)
cols_targets <- c(cols_targets, other = "grey70")

# order MAGs by mean_pi_ratio (descending; change as you like)
mag_order <- plot_dt %>% group_by(MAG_id) %>% summarize(m = median(mean_pi_ratio, na.rm=TRUE)) %>% arrange(desc(m)) %>% pull(MAG_id)
plot_dt$MAG_id_f <- factor(plot_dt$MAG_id, levels = mag_order)

# create transform
squish_trans <- squish_log_trans(lmin = xmin_allowed, mid_low = mid_low, mid_high = mid_high, rmax = xmax_allowed, mid_frac = mid_frac)

# ---------- PLOT A: original orientation (MAGs as rows) using coord_flip() -------------
p1 <- ggplot(plot_dt, aes(x = mean_pi_ratio, y = MAG_id_f, color = Genus_plot)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_vline(xintercept = global_ci[1], linetype="dashed") +
  geom_vline(xintercept = global_ci[2], linetype="dashed") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  scale_x_continuous(trans = squish_trans, breaks = c(1e-3,1e-2,1e-1,1), labels = c("10^-3","10^-2","10^-1","10^0"), limits = c(xmin_allowed, xmax_allowed)) +
  scale_color_manual(values = cols_targets) +
  labs(x = expression(pi[N]/pi[S]), y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank()) 

# ---------- PLOT B: without coord_flip() — swap axes to get same visual result -------------
# you asked for "invert it now" — this is how to get same visual without coord_flip:
p2 <- ggplot(plot_dt, aes(x = MAG_id_f, y = mean_pi_ratio, color = Genus_plot)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept = global_ci[1], linetype = "dashed") +
  geom_hline(yintercept = global_ci[2], linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_y_continuous(trans = squish_trans, breaks = c(1e-3,1e-2,1e-1,1), labels = c("10^-3","10^-2","10^-1","10^0"), limits = c(xmin_allowed, xmax_allowed)) +
  scale_color_manual(values = cols_targets) +
  labs(y = expression(pi[N]/pi[S]), x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())

# pick whichever orientation you want:
print(p1)  # if you prefer the "MAG rows" display, use this (coord_flip)
# OR
print(p2) # if you prefer no coord_flip, use this (x axis = MAGs, y = pi ratio)

# ---------- end ----------







# assume summ_dt_clean has columns corresponding_gene_call, MAG_id, pi_ratio
set.seed(42)
B <- 1000

# prepare: per-MAG aggregated vector of gene ratios or per-MAG mean
mag_gene_list <- summ_dt_clean %>%
  group_by(MAG_id) %>%
  summarize(pi_list = list(pi_ratio), n_genes = n(), .groups = "drop")

mag_ids <- mag_gene_list$MAG_id
n_mag <- length(mag_ids)

boot_mag_means <- replicate(B, {
  # sample MAGs with replacement
  sampled_mags <- sample(mag_ids, size = n_mag, replace = TRUE)
  # for each sampled MAG, compute its mean (use precomputed list)
  mag_means <- sapply(sampled_mags, function(m) mean(unlist(mag_gene_list$pi_list[mag_gene_list$MAG_id==m]), na.rm=TRUE))
  # return the overall mean across the sampled MAG means
  mean(mag_means, na.rm = TRUE)
})

global_ci_mag_boot <- quantile(boot_mag_means, c(0.025, 0.975))
global_mean_mag_boot <- mean(boot_mag_means)

global_mean_mag_boot
global_ci_mag_boot








############### PARALLELISM ############




# --- 1. Gene × MAG presence and 0D sweeps -----------------------------------

# 0D sweeps per (MAG, gene) – any colony/month
zeroD_gene_mag <- unique(zeroD_dt[, .(MAG_id, corresponding_gene_call)])
zeroD_gene_mag[, swept0D := 1L]

# All gene × MAG presences (from gene_dt)
# (if gene_dt has multiple rows per MAG/gene across colonies/months, unique() fixes it)
gene_mag_presence <- unique(gene_dt[, .(MAG_id, corresponding_gene_call)])

# Gene-level: in how many MAGs is this gene present? swept?
gene_presence_counts <- gene_mag_presence[
  , .(MAG_present = uniqueN(MAG_id)),
  by = corresponding_gene_call
]

gene_sweep_counts <- zeroD_gene_mag[
  , .(MAG_swept = uniqueN(MAG_id)),
  by = corresponding_gene_call
]

gene_summary <- merge(
  gene_presence_counts,
  gene_sweep_counts,
  by = "corresponding_gene_call",
  all.x = TRUE
)

gene_summary[is.na(MAG_swept), MAG_swept := 0L]

gene_summary[, swept_any   := MAG_swept >= 1L]
gene_summary[, swept_multi := MAG_swept >= 2L]

# Empirical probabilities
prob_swept_any       <- gene_summary[, mean(swept_any)]          # P(gene swept at least once)
prob_swept_multi     <- gene_summary[, mean(swept_multi)]        # P(gene swept in ≥2 MAGs)
prob_multi_given_2p  <- gene_summary[MAG_present >= 2L,
                                     mean(swept_multi)]          # P(parallel | present in ≥2 MAGs)

message("P(swept at least once)               = ", signif(prob_swept_any, 3))
message("P(swept in ≥2 MAGs)                  = ", signif(prob_swept_multi, 3))
message("P(swept in ≥2 MAGs | present ≥2 MAG) = ", signif(prob_multi_given_2p, 3))

# Optional: save the gene-level summary
fwrite(gene_summary,
       file.path(OUTDIR, "gene_level_0D_sweep_summary_across_MAGs.csv"))







#!/usr/bin/env Rscript

## ================================================================
##  0D sweep parallelism + clustering analysis
##  - Uses gene-level opportunities (effective_length_0D_sites)
##  - Tests whether 0D sweeps are randomly distributed across genes
##    and whether the same genes sweep in multiple MAGs more than expected
## ================================================================

suppressPackageStartupMessages({
  library(data.table)
})

## ------------------ User paths ----------------------------------

SWEEP_FILE   <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/sweep_gene_merged.csv"
GENE_SUM_FILE <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/final_analysis_all_colonies_aggregated_polymorphic_trajectory_SNVs/master_per_gene_summary.csv"
OUTDIR       <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ------------------ Load data -----------------------------------

message("Loading data...")
sweep_dt <- fread(SWEEP_FILE)
gene_dt  <- fread(GENE_SUM_FILE)

# canonicalize IDs as character
for (c in c("MAG_id", "corresponding_gene_call")) {
  if (c %in% names(gene_dt))  gene_dt[[c]]  <- as.character(gene_dt[[c]])
  if (c %in% names(sweep_dt)) sweep_dt[[c]] <- as.character(sweep_dt[[c]])
}

## ------------------ 0D sweeps table -----------------------------

if (!"degeneracy" %in% names(sweep_dt)) {
  stop("sweep_dt does not have a 'degeneracy' column; cannot select 0D sweeps.")
}

zeroD_dt <- sweep_dt[degeneracy == "0D"]

message("Total rows in sweep_dt:  ", nrow(sweep_dt))
message("Total 0D sweep rows:     ", nrow(zeroD_dt))

if (nrow(zeroD_dt) == 0L) {
  stop("No 0D sweeps found (degeneracy == '0D').")
}

## ------------------ Gene × MAG presence & sweeps ---------------

# genes present in each MAG (from gene_dt)
gene_mag_presence <- unique(gene_dt[, .(MAG_id, corresponding_gene_call)])

# genes with 0D sweeps in each MAG
zeroD_gene_mag <- unique(zeroD_dt[, .(MAG_id, corresponding_gene_call)])

message("Distinct MAGs in gene_dt:     ", uniqueN(gene_mag_presence$MAG_id))
message("Distinct genes in gene_dt:    ", uniqueN(gene_mag_presence$corresponding_gene_call))
message("Distinct MAGs with 0D sweeps: ", uniqueN(zeroD_gene_mag$MAG_id))
message("Distinct genes with 0D sweeps:", uniqueN(zeroD_gene_mag$corresponding_gene_call))

# gene-level counts: in how many MAGs is each gene present / swept?
gene_presence_counts <- gene_mag_presence[
  , .(MAG_present = uniqueN(MAG_id)),
  by = corresponding_gene_call
]

gene_sweep_counts <- zeroD_gene_mag[
  , .(MAG_swept = uniqueN(MAG_id)),
  by = corresponding_gene_call
]

gene_summary <- merge(
  gene_presence_counts,
  gene_sweep_counts,
  by = "corresponding_gene_call",
  all.x = TRUE
)

gene_summary[is.na(MAG_swept), MAG_swept := 0L]

# indicators
gene_summary[, swept_any   := MAG_swept >= 1L]
gene_summary[, swept_multi := MAG_swept >= 2L]

# observed probabilities
prob_swept_any      <- gene_summary[, mean(swept_any)]
prob_swept_multi    <- gene_summary[, mean(swept_multi)]
prob_multi_given_2p <- gene_summary[MAG_present >= 2L,
                                    mean(swept_multi)]

message("----------------------------------------------------")
message("Observed probabilities (0D sweeps, gene-level):")
message("  P(swept at least once)               = ", signif(prob_swept_any, 3))
message("  P(swept in ≥2 MAGs)                  = ", signif(prob_swept_multi, 3))
message("  P(swept in ≥2 MAGs | present ≥2 MAG) = ", signif(prob_multi_given_2p, 3))
message("----------------------------------------------------")

# Save gene-level summary
fwrite(gene_summary,
       file.path(OUTDIR, "gene_level_0D_sweep_summary_across_MAGs.csv"))

## ------------------ Opportunity (0D sites) ----------------------

if (!"effective_length_0D_sites" %in% names(gene_dt)) {
  stop("gene_dt does not contain 'effective_length_0D_sites'; cannot weight by 0D opportunity.")
}

# Take one row per MAG×gene for the 0D opportunity
op_dt <- unique(
  gene_dt[, .(MAG_id, corresponding_gene_call, opportunity_0D = effective_length_0D_sites)]
)

# If duplicates remain with conflicting values, collapse explicitly
op_dt <- op_dt[
  , .(opportunity_0D = unique(opportunity_0D)[1]),  # take first unique
  by = .(MAG_id, corresponding_gene_call)
]

# Merge opportunity into presence table
gene_mag_presence <- merge(
  gene_mag_presence,
  op_dt,
  by = c("MAG_id", "corresponding_gene_call"),
  all.x = TRUE
)

# Fill missing opportunity with median across genes
med_op <- median(gene_mag_presence$opportunity_0D, na.rm = TRUE)
if (!is.finite(med_op)) {
  stop("Median opportunity_0D is not finite; check effective_length_0D_sites.")
}
gene_mag_presence[is.na(opportunity_0D), opportunity_0D := med_op]

## ------------------ Per-MAG swept counts ------------------------

k_dt <- zeroD_gene_mag[, .N, by = MAG_id]
setnames(k_dt, "N", "n_swept")
k_vec <- setNames(k_dt$n_swept, k_dt$MAG_id)

## ------------------ Build lists for permutation -----------------

mag_genes_list   <- split(gene_mag_presence$corresponding_gene_call,
                          gene_mag_presence$MAG_id)
mag_weights_list <- split(gene_mag_presence$opportunity_0D,
                          gene_mag_presence$MAG_id)

mag_ids <- names(mag_genes_list)

## ------------------ Weighted permutation function ---------------

# Returns:
#  - parallel: P(gene swept in ≥2 MAGs | present ≥2 MAGs)
#  - any:      P(gene swept at least once) across all genes
perm_once_weighted <- function() {
  perm_list <- vector("list", length(mag_ids))
  names(perm_list) <- mag_ids
  
  for (i in seq_along(mag_ids)) {
    m <- mag_ids[i]
    genes_m   <- mag_genes_list[[m]]
    weights_m <- mag_weights_list[[m]]
    k         <- k_vec[m]  # number of swept genes in MAG m (observed)
    
    if (is.na(k) || k == 0L || length(genes_m) == 0L) {
      perm_list[[i]] <- NULL
      next
    }
    
    # safety
    if (k > length(genes_m)) k <- length(genes_m)
    
    perm_list[[i]] <- data.table(
      MAG_id = m,
      corresponding_gene_call = sample(
        genes_m,
        size    = k,
        prob    = weights_m,
        replace = FALSE
      )
    )
  }
  
  perm_zeroD <- rbindlist(perm_list, use.names = TRUE, fill = TRUE)
  if (nrow(perm_zeroD) == 0L) {
    return(c(parallel = NA_real_, any = NA_real_))
  }
  
  perm_counts <- perm_zeroD[
    , .(MAG_swept = uniqueN(MAG_id)),
    by = corresponding_gene_call
  ]
  
  perm_gene <- merge(
    gene_summary[, .(corresponding_gene_call, MAG_present)],
    perm_counts,
    by = "corresponding_gene_call",
    all.x = TRUE
  )
  
  perm_gene[is.na(MAG_swept), MAG_swept := 0L]
  
  stat_parallel <- perm_gene[MAG_present >= 2L, mean(MAG_swept >= 2L)]
  stat_any      <- perm_gene[,                       mean(MAG_swept >= 1L)]
  
  c(parallel = stat_parallel, any = stat_any)
}

## ------------------ Run permutations ----------------------------

set.seed(1)
n_perm <- 1000L

message("Running ", n_perm, " weighted permutations...")
perm_mat <- replicate(n_perm, perm_once_weighted())
perm_mat <- as.matrix(perm_mat)

# Drop columns where everything is NA (shouldn't happen, but check)
keep_cols <- colSums(is.na(perm_mat)) < nrow(perm_mat)
perm_mat  <- perm_mat[, keep_cols, drop = FALSE]

perm_parallel <- perm_mat["parallel", ]
perm_any      <- perm_mat["any", ]

obs_prob_parallel <- prob_multi_given_2p
obs_prob_any      <- prob_swept_any

# Summaries under null
perm_mean_parallel <- mean(perm_parallel, na.rm = TRUE)
perm_sd_parallel   <- sd(perm_parallel, na.rm = TRUE)

perm_mean_any <- mean(perm_any, na.rm = TRUE)
perm_sd_any   <- sd(perm_any, na.rm = TRUE)

# Empirical p-values: P(null >= observed)
p_emp_parallel <- (sum(perm_parallel >= obs_prob_parallel, na.rm = TRUE) + 1) /
  (sum(!is.na(perm_parallel)) + 1)

p_emp_any <- (sum(perm_any >= obs_prob_any, na.rm = TRUE) + 1) /
  (sum(!is.na(perm_any)) + 1)

message("----------------------------------------------------")
message("Weighted null (opportunity_0D-based) results:")
message(" Parallelism (genes swept in ≥2 MAGs | present ≥2 MAGs)")
message("   Observed       = ", signif(obs_prob_parallel, 3))
message("   Null mean      = ", signif(perm_mean_parallel, 3),
        " ± ", signif(perm_sd_parallel, 3))
message("   Empirical p    = ", signif(p_emp_parallel, 3))
message("   Fold-enrichment= ",
        signif(obs_prob_parallel / perm_mean_parallel, 3))

message("")
message(" Clustering across genes (P(gene swept at least once))")
message("   Observed       = ", signif(obs_prob_any, 3))
message("   Null mean      = ", signif(perm_mean_any, 3),
        " ± ", signif(perm_sd_any, 3))
message("   Empirical p    = ", signif(p_emp_any, 3))
message("   Fold-enrichment= ",
        signif(obs_prob_any / perm_mean_any, 3))
message("----------------------------------------------------")

## Interpretation:
##  - Parallelism test:
##      If obs_prob_parallel > perm_mean_parallel with small p_emp_parallel,
##      then the same gene is swept in multiple MAGs more often than expected
##      under random distribution of sweeps across 0D sites.
##
##  - Clustering test:
##      If obs_prob_any < perm_mean_any, sweeps are concentrated into fewer
##      genes than expected by chance (more clustering into a subset of genes).
##      If obs_prob_any > perm_mean_any, sweeps are spread over more genes.

## ------------------ Save permutation distributions --------------

perm_dt <- data.table(
  perm_id               = seq_len(length(perm_parallel)),
  P_parallel_given_2MAG = perm_parallel,
  P_swept_any           = perm_any
)

fwrite(perm_dt,
       file.path(OUTDIR, "perm_weighted_parallel_and_any_0D.csv"))

message("Analysis complete. Results written to: ", OUTDIR)



gene_clean












## ------------------ Weighted permutation (prob + counts) -------------------

perm_once_weighted <- function() {
  perm_list <- vector("list", length(mag_ids))
  names(perm_list) <- mag_ids
  
  for (i in seq_along(mag_ids)) {
    m <- mag_ids[i]
    genes_m   <- mag_genes_list[[m]]
    weights_m <- mag_weights_list[[m]]
    k         <- k_vec[m]  # number of swept genes in MAG m (observed)
    
    if (is.na(k) || k == 0L || length(genes_m) == 0L) {
      perm_list[[i]] <- NULL
      next
    }
    
    if (k > length(genes_m)) k <- length(genes_m)
    
    perm_list[[i]] <- data.table(
      MAG_id = m,
      corresponding_gene_call = sample(
        genes_m,
        size    = k,
        prob    = weights_m,
        replace = FALSE
      )
    )
  }
  
  perm_zeroD <- rbindlist(perm_list, use.names = TRUE, fill = TRUE)
  if (nrow(perm_zeroD) == 0L) {
    return(c(prob_parallel = NA_real_,
             n_parallel    = NA_integer_,
             prob_any      = NA_real_))
  }
  
  perm_counts <- perm_zeroD[
    , .(MAG_swept = uniqueN(MAG_id)),
    by = corresponding_gene_call
  ]
  
  perm_gene <- merge(
    gene_summary[, .(corresponding_gene_call, MAG_present)],
    perm_counts,
    by = "corresponding_gene_call",
    all.x = TRUE
  )
  
  perm_gene[is.na(MAG_swept), MAG_swept := 0L]
  
  # stats:
  # P(swept in ≥2 MAGs | present ≥2 MAGs)
  prob_parallel <- perm_gene[MAG_present >= 2L, mean(MAG_swept >= 2L)]
  # count of genes with ≥2 sweeps (among genes present ≥2 MAGs)
  n_parallel    <- perm_gene[MAG_present >= 2L, sum(MAG_swept >= 2L)]
  # P(gene swept at least once) overall
  prob_any      <- perm_gene[, mean(MAG_swept >= 1L)]
  
  c(prob_parallel = prob_parallel,
    n_parallel    = n_parallel,
    prob_any      = prob_any)
}



set.seed(1)
n_perm <- 1000L

message("Running ", n_perm, " weighted permutations (prob + counts)...")
perm_mat <- replicate(n_perm, perm_once_weighted())
perm_mat <- t(perm_mat)  # rows = perms, cols = stats
perm_dt  <- as.data.table(perm_mat)
setnames(perm_dt, c("prob_parallel", "n_parallel", "prob_any"))

# drop NA perms if any
perm_dt  <- perm_dt[!is.na(prob_parallel)]

# observed values
obs_prob_parallel <- prob_multi_given_2p
obs_n_parallel    <- gene_summary[MAG_present >= 2L & MAG_swept >= 2L, .N]
obs_prob_any      <- prob_swept_any

# recompute p-values with counts if you want
p_emp_parallel <- (sum(perm_dt$prob_parallel >= obs_prob_parallel) + 1) /
  (nrow(perm_dt) + 1)
p_emp_any <- (sum(perm_dt$prob_any >= obs_prob_any) + 1) /
  (nrow(perm_dt) + 1)

message("Observed n_parallel genes (MAG_swept ≥2 | present≥2) = ", obs_n_parallel)
message("Empirical p for parallel prob    = ", signif(p_emp_parallel, 3))
message("Empirical p for 'any sweep' prob = ", signif(p_emp_any, 3))

fwrite(perm_dt,
       file.path(OUTDIR, "perm_weighted_parallel_prob_and_counts_0D.csv"))


p_counts <- ggplot(perm_dt, aes(x = n_parallel)) +
  geom_histogram(binwidth = 1, color = "black") +
  geom_vline(xintercept = obs_n_parallel, linetype = "dashed", size = 0.8) +
  labs(
    x = "Number of genes with 0D sweeps in ≥2 MAGs",
    y = "Permutations",
    title = "Null distribution of parallel 0D-swept genes (opportunity-weighted)",
    subtitle = paste0("Observed = ", obs_n_parallel,
                      " (p = ", signif(p_emp_parallel, 3), ")")
  ) +
  theme_classic()
p_counts
ggsave(file.path(OUTDIR, "parallel_0D_sweeps_count_null_vs_observed.png"),
       p_counts, width = 6, height = 4, dpi = 300)












# Use gene_clean as the annotated gene table
gene_ann <- as.data.table(gene_clean)

# make sure IDs are character
gene_ann[, corresponding_gene_call := as.character(corresponding_gene_call)]

# background: genes present in ≥2 MAGs (same as for prob_multi_given_2p)
bg <- merge(
  gene_summary[MAG_present >= 2L,
               .(corresponding_gene_call, MAG_present, MAG_swept)],
  unique(gene_ann[, .(corresponding_gene_call, GO_primary, PFAM_primary)]),
  by = "corresponding_gene_call",
  all.x = TRUE
)

bg[, is_parallel := MAG_swept >= 2L]  # TRUE/FALSE




# drop unannotated if you want GO-based enrichment
bg_go <- bg[!is.na(GO_primary)]

# N total parallel/non-parallel in this background (for sanity)
total_parallel    <- sum(bg_go$is_parallel)
total_nonparallel <- sum(!bg_go$is_parallel)

go_stats <- bg_go[
  ,
  {
    # a: parallel AND has this GO
    a <- sum(is_parallel)
    # b: not parallel AND has this GO
    b <- sum(!is_parallel)
    # c: parallel AND NOT this GO
    c <- total_parallel    - a
    # d: not parallel AND NOT this GO
    d <- total_nonparallel - b
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    ft  <- fisher.test(mat)
    
    list(
      a = a, b = b, c = c, d = d,
      odds_ratio = unname(ft$estimate),
      p_value    = ft$p.value
    )
  },
  by = GO_primary
]

go_stats[, p_adj := p.adjust(p_value, method = "BH")]
setorder(go_stats, p_value)

fwrite(go_stats, file.path(OUTDIR, "GO_enrichment_parallel_0D_swept_genes.csv"))




bg_pf <- bg[!is.na(PFAM_primary)]

total_parallel_pf    <- sum(bg_pf$is_parallel)
total_nonparallel_pf <- sum(!bg_pf$is_parallel)

bg_pf <- bg[!is.na(PFAM_primary)]

pfam_counts <- bg_pf[
  ,
  .(
    n_bg        = .N,
    n_parallel  = sum(is_parallel)
  ),
  by = PFAM_primary
]



pfam_stats[, p_adj := p.adjust(p_value, method = "BH")]
setorder(pfam_stats, p_value)



fwrite(pfam_stats, file.path(OUTDIR, "PFAM_enrichment_parallel_0D_swept_genes.csv"))









zeroD_dt <- as.data.table(zeroD_dt)
gene_clean <- as.data.table(gene_clean)
for (c in c("MAG_id","colony_id","corresponding_gene_call","Month","Month_t1","Month_t2")) {
  if (c %in% names(zeroD_dt)) zeroD_dt[[c]] <- as.character(zeroD_dt[[c]])
  if (c %in% names(gene_clean)) gene_clean[[c]] <- as.character(gene_clean[[c]])
}

# 0D sweeps summarized to (gene, MAG, colony, month_of_sweep)
zeroD_by_month <- unique(
  zeroD_dt[!is.na(Month_t2),
           .(corresponding_gene_call, MAG_id, colony_id, Month = Month_t2)]
)

# same-month sweeps in ≥2 MAGs
same_month_sweeps <- zeroD_by_month[
  ,
  .(
    n_MAG    = uniqueN(MAG_id),
    MAGs     = list(sort(unique(MAG_id))),
    Colonies = list(sort(unique(colony_id)))
  ),
  by = .(corresponding_gene_call, Month)
][n_MAG >= 2L]




# same-month pN/pS > 1 across MAGs
pnps_by_month <- gene_clean[
  !is.na(pN_pS_ratio) & pN_pS_ratio > 1,
  .(
    n_MAG_pnps    = uniqueN(MAG_id),
    MAGs_pnps     = list(sort(unique(MAG_id))),
    Colonies_pnps = list(sort(unique(colony_id)))
  ),
  by = .(corresponding_gene_call, Month)
]

pnps_same_month <- pnps_by_month[n_MAG_pnps >= 2L]





# intersect the same (gene, Month) that satisfy both conditions
narrative_candidates <- merge(
  zeroD_events,
  pnps_events,
  by = c("corresponding_gene_call"),
  all = FALSE,
  suffixes = c("_sweep", "_pnps")
)

message("Number of (gene, month) narrative candidates = ",
        nrow(narrative_candidates))







#### FST ###


OUTDIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel"

# Fst files from Python (gene-level)
FST_AVG_FILE <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/temporal_fst_gene_level/gene_average_temporal_fst.csv"
FST_VS0_FILE <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/temporal_fst_gene_level/gene_vs_initial_month_fst.csv"

# You should already have gene_clean and zeroD_dt in memory.
# If not, reload them here.
# source("script_that_builds_gene_clean.R")  # or paste from earlier
# zeroD_dt <- fread("...")  # if needed

## ------------------------------------------------------------
## 0. Load Fst and FILTER OUT negative values
## ------------------------------------------------------------

fst_avg_raw <- fread(FST_AVG_FILE)
fst_vs0_raw <- fread(FST_VS0_FILE)

# Ensure IDs are character
for (c in c("Colony","MAG_id","corresponding_gene_call")) {
  if (c %in% names(fst_avg_raw)) fst_avg_raw[[c]] <- as.character(fst_avg_raw[[c]])
  if (c %in% names(fst_vs0_raw)) fst_vs0_raw[[c]] <- as.character(fst_vs0_raw[[c]])
}

# Choose an Fst metric (vs initial month, Method1)
fst_avg_raw[, Fst := Average_Fst_vs_Initial_Month_Method1]

# Filter out negative Fst and NA
fst_avg <- fst_avg_raw[!is.na(Fst) & Fst >= 0]

# For per-comparison Fst (used for seasons), also filter negative
fst_vs0_raw[, Fst := Fst_Method1]
fst_vs0 <- fst_vs0_raw[!is.na(Fst) & Fst >= 0]

message("Rows in fst_avg_raw: ", nrow(fst_avg_raw),
        " -> after filtering (Fst>=0): ", nrow(fst_avg))
message("Rows in fst_vs0_raw: ", nrow(fst_vs0_raw),
        " -> after filtering (Fst>=0): ", nrow(fst_vs0))

## ------------------------------------------------------------
## 1. Rebuild fst_annot (Fst + gene annotations)
## ------------------------------------------------------------

gene_clean <- as.data.table(gene_clean)
gene_clean[, corresponding_gene_call := as.character(corresponding_gene_call)]

if ("colony_id" %in% names(gene_clean)) {
  gene_clean[, Colony := as.character(colony_id)]
}
gene_clean
fst_annot <- merge(
  fst_avg,
  unique(gene_clean[, .(Colony, MAG_id, corresponding_gene_call,
                        Genus, GO_primary, PFAM_primary,
                        pN_pS_ratio, GO_str, pfam_str)]),
  by = c("Colony","MAG_id","corresponding_gene_call"),
  all.x = TRUE
)

# Keep only rows with non-NA Fst
fst_annot <- fst_annot[!is.na(Fst)]
fst_annot
## ------------------------------------------------------------
## 2. Volcano plot: Z-scored Fst vs FDR
## ------------------------------------------------------------

mu  <- mean(fst_annot$Fst, na.rm = TRUE)
sig <- sd(fst_annot$Fst, na.rm = TRUE)

fst_annot[, Z_Fst := (Fst - mu) / sig]
fst_annot[, p_two_sided := 2 * pnorm(abs(Z_Fst), lower.tail = FALSE)]
fst_annot[, q_two_sided := p.adjust(p_two_sided, method = "BH")]
fst_annot[, neg_log10_q := -log10(pmax(q_two_sided, 1e-300))]

# label a few top outliers
fst_annot[, label_volcano := NA_character_]
fst_annot[rank(p_two_sided, ties.method = "first") <= 20,
          label_volcano := corresponding_gene_call]

p_volcano <- ggplot(fst_annot, aes(x = Z_Fst, y = neg_log10_q)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_text_repel(aes(label = label_volcano),
                  size = 2.2, max.overlaps = 20, na.rm = TRUE) +
  labs(
    x = "Z-scored temporal Fst (vs initial month)",
    y = "-log10(FDR)",
    title = "Gene-level temporal Fst outliers (Fst ≥ 0)"
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "volcano_temporal_Fst_genes_Fst_ge0.png"),
       p_volcano, width = 6, height = 5, dpi = 300)

## ------------------------------------------------------------
## 3. PFAM / GO group Fst outliers
## ------------------------------------------------------------

global_mu <- mu
global_sd <- sig

### 3.1 PFAM

min_genes_pf <- 20

pfam_stats <- fst_annot[!is.na(PFAM_primary),
                        .(
                          n_genes = .N,
                          mean_Fst = mean(Fst, na.rm = TRUE)
                        ),
                        by = PFAM_primary
][n_genes >= min_genes_pf]

pfam_stats[, Z_group := (mean_Fst - global_mu) / (global_sd / sqrt(n_genes))]
pfam_stats[, p_two_sided := 2 * pnorm(abs(Z_group), lower.tail = FALSE)]
pfam_stats[, q_BH := p.adjust(p_two_sided, method = "BH")]
pfam_stats <- pfam_stats[order(p_two_sided)]

fwrite(pfam_stats,
       file.path(OUTDIR, "PFAM_temporal_Fst_group_outliers_Fst_ge0.csv"))

### 3.2 GO

min_genes_go <- 20

go_stats <- fst_annot[!is.na(GO_primary),
                      .(
                        n_genes = .N,
                        mean_Fst = mean(Fst, na.rm = TRUE)
                      ),
                      by = GO_primary
][n_genes >= min_genes_go]

go_stats[, Z_group := (mean_Fst - global_mu) / (global_sd / sqrt(n_genes))]
go_stats[, p_two_sided := 2 * pnorm(abs(Z_group), lower.tail = FALSE)]
go_stats[, q_BH := p.adjust(p_two_sided, method = "BH")]
go_stats <- go_stats[order(p_two_sided)]

fwrite(go_stats,
       file.path(OUTDIR, "GO_temporal_Fst_group_outliers_Fst_ge0.csv"))

## ------------------------------------------------------------
## 4. MAG-level and Genus-level mean Fst
## ------------------------------------------------------------

# MAG-level
fst_mag <- fst_annot[
  ,
  .(
    n_genes = .N,
    mean_Fst = mean(Fst, na.rm = TRUE)
  ),
  by = .(Colony, MAG_id)
]

mu_mag <- mean(fst_mag$mean_Fst, na.rm = TRUE)
sd_mag <- sd(fst_mag$mean_Fst, na.rm = TRUE)
fst_mag[, Z_mag := (mean_Fst - mu_mag) / sd_mag]

fwrite(fst_mag,
       file.path(OUTDIR, "MAG_level_mean_temporal_Fst_Fst_ge0.csv"))

p_mag <- ggplot(fst_mag, aes(x = mean_Fst)) +
  geom_histogram(bins = 30, color = "black", fill = "grey70") +
  labs(
    x = "Mean temporal Fst (per MAG_id × Colony)",
    y = "Count",
    title = "Distribution of MAG-level temporal Fst (Fst ≥ 0)"
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "hist_MAG_temporal_Fst_Fst_ge0.png"),
       p_mag, width = 6, height = 4, dpi = 300)

# Genus-level
fst_genus <- fst_annot[
  !is.na(Genus),
  .(
    n_genes = .N,
    mean_Fst = mean(Fst, na.rm = TRUE)
  ),
  by = Genus
]

fst_genus <- fst_genus[order(-mean_Fst)]
fwrite(fst_genus,
       file.path(OUTDIR, "Genus_level_mean_temporal_Fst_Fst_ge0.csv"))

p_genus <- ggplot(fst_genus,
                  aes(x = reorder(Genus, mean_Fst), y = mean_Fst)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Genus",
    y = "Mean temporal Fst",
    title = "Genus-level temporal Fst (Fst ≥ 0)"
  ) +
  theme_classic()
p_genus
ggsave(file.path(OUTDIR, "Genus_temporal_Fst_barplot_Fst_ge0.png"),
       p_genus, width = 6, height = 5, dpi = 300)

## ------------------------------------------------------------
## 5. Fst vs 0D sweeps and pN/pS > 1
## ------------------------------------------------------------

zeroD_dt <- as.data.table(zeroD_dt)
zeroD_gene_mag <- unique(zeroD_dt[, .(MAG_id, corresponding_gene_call)])
zeroD_gene_mag[, corresponding_gene_call := as.character(corresponding_gene_call)]
zeroD_gene_mag[, has_0D_sweep := TRUE]

fst_sweep <- merge(
  fst_annot,
  zeroD_gene_mag,
  by = c("MAG_id", "corresponding_gene_call"),
  all.x = TRUE
)
fst_sweep[is.na(has_0D_sweep), has_0D_sweep := FALSE]

pnps_flag <- unique(
  gene_clean[!is.na(pN_pS_ratio) & pN_pS_ratio > 1,
             .(Colony, MAG_id, corresponding_gene_call)]
)
pnps_flag[, corresponding_gene_call := as.character(corresponding_gene_call)]
pnps_flag[, pnps_gt1 := TRUE]

fst_sweep <- merge(
  fst_sweep,
  pnps_flag,
  by = c("Colony","MAG_id","corresponding_gene_call"),
  all.x = TRUE
)
fst_sweep[is.na(pnps_gt1), pnps_gt1 := FALSE]

# 0D sweeps vs not
fst_sweep[, group_sweep := ifelse(has_0D_sweep, "0D_sweep", "no_0D_sweep")]

t_sweep   <- t.test(Fst ~ group_sweep, data = fst_sweep)
wil_sweep <- wilcox.test(Fst ~ group_sweep, data = fst_sweep)
print(t_sweep)
print(wil_sweep)

p_sweep <- ggplot(fst_sweep, aes(x = group_sweep, y = Fst)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  labs(
    x = "",
    y = "Temporal Fst",
    title = "Temporal Fst for genes with vs without 0D sweeps (Fst ≥ 0)"
  ) +
  theme_classic()+
  scale_y_continuous(trans="log10")

p_sweep
ggsave(file.path(OUTDIR, "Fst_0Dsweep_vs_not_Fst_ge0.png"),
       p_sweep, width = 5, height = 4, dpi = 300)

# pN/pS > 1 vs ≤ 1
fst_sweep[, group_pnps := ifelse(pnps_gt1, "pN/pS > 1", "pN/pS ≤ 1")]

t_pnps   <- t.test(Fst ~ group_pnps, data = fst_sweep)
wil_pnps <- wilcox.test(Fst ~ group_pnps, data = fst_sweep)
print(t_pnps)
print(wil_pnps)

p_pnps <- ggplot(fst_sweep, aes(x = group_pnps, y = Fst)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  labs(
    x = "",
    y = "Temporal Fst",
    title = "Temporal Fst for genes with vs without pN/pS > 1 (Fst ≥ 0)"
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "Fst_pnps_gt1_vs_not_Fst_ge0.png"),
       p_pnps, width = 5, height = 4, dpi = 300)

## ------------------------------------------------------------
## 6. Fst by Season (using per-comparison Fst, Fst>=0 only)
## ------------------------------------------------------------

season_map <- c("May"="Spring","June"="Spring",
                "July"="Summer","August"="Summer","September"="Summer",
                "October"="Fall","November"="Fall",
                "January"="Winter","February"="Winter")

fst_vs0[, Season := season_map[Month]]

fst_vs0_season <- merge(
  fst_vs0,
  unique(gene_clean[, .(Colony, MAG_id, corresponding_gene_call, Genus, GO_primary, PFAM_primary)]),
  by = c("Colony","MAG_id","corresponding_gene_call"),
  all.x = TRUE
)

fst_vs0_season <- fst_vs0_season[!is.na(Season) & !is.na(Fst)]

p_season <- ggplot(fst_vs0_season, aes(x = Season, y = Fst)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.15) +
  labs(
    x = "Season of later timepoint",
    y = "Temporal Fst (vs initial month)",
    title = "Temporal Fst by season (Fst ≥ 0)"
  ) +
  theme_classic()+
  scale_y_continuous(trans="log10")

ggsave(file.path(OUTDIR, "Fst_by_season_Fst_ge0.png"),
       p_season, width = 5.5, height = 4, dpi = 300)
p_season
anova_season <- aov(Fst ~ Season, data = fst_vs0_season)
print(summary(anova_season))

print(kruskal.test(Fst ~ Season, data = fst_vs0_season))

message("All analyses recomputed using Fst >= 0 only.")






season_levels <- c("Spring","Summer","Fall","Winter")

fst_vs0_season[, Season := factor(Season, levels = season_levels, ordered = TRUE)]

p_season_ordered <- ggplot(fst_vs0_season, aes(x = Season, y = Fst)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.15) +
  labs(
    x = "Season of later timepoint",
    y = "Temporal Fst (vs initial month)",
    title = "Temporal Fst by season (Fst ≥ 0)"
  ) +
  theme_classic()
p_season_ordered
ggsave(file.path(OUTDIR, "Fst_by_season_ordered_Fst_ge0.png"),
       p_season_ordered, width = 5.5, height = 4, dpi = 300)



target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium",
                   "commensalibacter", "frischella", "gilliamella",
                   "lactobacillus", "snodgrassella")

fst_target <- fst_annot[!is.na(Genus) & Genus %chin% target_genera]
fst_vs0_target <- fst_vs0_season[!is.na(Genus) & Genus %chin% target_genera]



var_partition <- function(dt, value_col, group_col) {
  dt <- dt[!is.na(get(value_col)) & !is.na(get(group_col))]
  if (nrow(dt) < 2) return(list(between = NA_real_, within = NA_real_, ratio = NA_real_))
  
  # group means
  group_means <- dt[, .(mean_val = mean(get(value_col), na.rm = TRUE)), by = group_col]
  between_var <- var(group_means$mean_val, na.rm = TRUE)
  
  # within variance per group, then average
  within_vals <- dt[, .(var_val = var(get(value_col), na.rm = TRUE)), by = group_col]
  within_var <- mean(within_vals$var_val, na.rm = TRUE)
  
  list(between = between_var, within = within_var,
       ratio = between_var / within_var)
}






# MAG-level (across genes and comparisons within MAG)
mag_part <- var_partition(fst_target, value_col = "Fst", group_col = "MAG_id")

# Genus-level
genus_part <- var_partition(fst_target, value_col = "Fst", group_col = "Genus")

# Colony-level
colony_part <- var_partition(fst_target, value_col = "Fst", group_col = "Colony")

# Season-level (using per-comparison Fst)
season_part <- var_partition(fst_vs0_target, value_col = "Fst", group_col = "Season")

mag_part
genus_part
colony_part
season_part






# per-MAG stats
fst_mag_target <- fst_target[
  ,
  .(
    n_genes = .N,
    mean_Fst = mean(Fst, na.rm = TRUE),
    sd_Fst   = sd(Fst, na.rm = TRUE)
  ),
  by = .(Genus, Colony, MAG_id)
]

# overall variance of MAG means
var_mag_means <- var(fst_mag_target$mean_Fst, na.rm = TRUE)
# average within-MAG variance
mean_within_mag_var <- mean(fst_mag_target$sd_Fst^2, na.rm = TRUE)

var_mag_means
mean_within_mag_var
var_mag_means / mean_within_mag_var





# Season factor ordered
fst_vs0_target[, Season := factor(Season, levels = c("Spring","Summer","Fall","Winter"), ordered = TRUE)]

fst_genus_season <- fst_vs0_target[
  ,
  .(
    n_obs   = .N,
    mean_Fst = mean(Fst, na.rm = TRUE),
    sd_Fst   = sd(Fst, na.rm = TRUE)
  ),
  by = .(Genus, Season)
][n_obs >= 10]

fwrite(fst_genus_season,
       file.path(OUTDIR, "Genus_Season_temporal_Fst_target_genera.csv"))
fst_genus_season



anova_genus <- aov(Fst ~ Genus + Season, data = fst_vs0_target)
summary(anova_genus)





fst_mag_colony <- fst_target[
  ,
  .(
    n_genes = .N,
    mean_Fst = mean(Fst, na.rm = TRUE)
  ),
  by = .(MAG_id, Colony, Genus)
]



mag_multi_colony <- fst_mag_colony[
  ,
  .N,
  by = MAG_id
][N > 1, MAG_id]

fst_mag_multi <- fst_mag_colony[MAG_id %chin% mag_multi_colony]

# For each such MAG, variance of mean_Fst across colonies
mag_between_col_var <- fst_mag_multi[
  ,
  .(var_between_colonies = var(mean_Fst, na.rm = TRUE),
    n_col = .N),
  by = MAG_id
]

# Typical variance across colonies for the same MAG:
mean(mag_between_col_var$var_between_colonies, na.rm = TRUE)





# reuse fst_mag_target or recompute for same subset
fst_mag_gene_sd <- fst_target[
  ,
  .(
    n_genes = .N,
    sd_Fst = sd(Fst, na.rm = TRUE)
  ),
  by = .(MAG_id, Colony)
]

# average within-MAG variance across (MAG, Colony)
mean_within_mag_var <- mean((fst_mag_gene_sd$sd_Fst)^2, na.rm = TRUE)

mean_between_col_var <- mean(mag_between_col_var$var_between_colonies, na.rm = TRUE)

mean_between_col_var
mean_within_mag_var
mean_between_col_var / mean_within_mag_var




fst_vs0_target <- fst_vs0_season[!is.na(Genus) & Genus %chin% target_genera]

# Your month-time function, adapted to data.table and usable for any month column
add_cumulative_time_from_month_nucdiv <- function(dt, month_col = "Month", new_col = "Time") {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34, 
                        "September" = 28, "October" = 22, "November" = 42, 
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  
  dt[, (new_col) := month_to_cumulative_time_map[as.character(get(month_col))]]
  dt
}

# Work on a copy
fst_time <- copy(fst_vs0_target)

# Add absolute times for the later and initial month of each comparison
fst_time <- add_cumulative_time_from_month_nucdiv(fst_time, month_col = "Month",         new_col = "Time_month")
fst_time <- add_cumulative_time_from_month_nucdiv(fst_time, month_col = "Initial_Month", new_col = "Time_initial")

# Temporal distance between samples in each Fst comparison
fst_time[, DeltaTime := Time_month - Time_initial]

# Sanity filter
fst_time <- fst_time[!is.na(DeltaTime) & DeltaTime >= 0 & !is.na(Fst)]
summary(fst_time$DeltaTime)





season_levels <- c("Spring","Summer","Fall","Winter")

fst_season <- copy(fst_vs0_target)
fst_season <- fst_season[!is.na(Fst) & !is.na(Season)]
fst_season[, Season := factor(Season, levels = season_levels, ordered = TRUE)]
fst_season[, MAG_id := factor(MAG_id)]





m_season_mag <- lmer(
  Fst ~ Season + (1 | MAG_id),
  data = fst_season
)

summary(m_season_mag)
anova(m_season_mag)      # season effect
VarCorr(m_season_mag)


vc1 <- as.data.frame(VarCorr(m_season_mag))
vc1
vc1$prop <- vc1$vcov / sum(vc1$vcov)
vc1


emm_season <- emmeans(m_season_mag, ~ Season)
emm_season
pairs(emm_season)   # seasonal contrasts










fst_season_pfam <- fst_season[!is.na(PFAM_primary)]
fst_season_pfam[, PFAM_primary := factor(PFAM_primary)]

m_season_mag_pfam <- lmer(
  Fst ~ Season + (1 | MAG_id) + (1 | PFAM_primary),
  data = fst_season_pfam
)

summary(m_season_mag_pfam)
VarCorr(m_season_mag_pfam)



vc2 <- as.data.frame(VarCorr(m_season_mag_pfam))
vc2
vc2$prop <- vc2$vcov / sum(vc2$vcov)
vc2



fst_season_go <- fst_season[!is.na(GO_primary)]
fst_season_go[, GO_primary := factor(GO_primary)]

m_season_mag_go <- lmer(
  Fst ~ Season + (1 | MAG_id) + (1 | GO_primary),
  data = fst_season_go
)

vc_go <- as.data.frame(VarCorr(m_season_mag_go))
vc_go$prop <- vc_go$vcov / sum(vc_go$vcov)
vc_go




m_no_season <- lmer(
  Fst ~ 1 + (1 | MAG_id) + (1 | PFAM_primary),
  data = fst_season_pfam,
  REML = FALSE
)

m_with_season <- update(m_no_season, . ~ Season + (1 | MAG_id) + (1 | PFAM_primary))

anova(m_no_season, m_with_season)













## ------------------------------------------------------------
## 7. Parallel 0D genes vs MAG-specific temporal Fst distribution
## ------------------------------------------------------------


# Make sure IDs are character
fst_annot[, MAG_id := as.character(MAG_id)]
fst_annot[, corresponding_gene_call := as.character(corresponding_gene_call)]
gene_clean[, MAG_id := as.character(MAG_id)]
gene_clean[, corresponding_gene_call := as.character(corresponding_gene_call)]

# 7.1 Collapse to MAG × gene mean temporal Fst
gene_mag_fst <- fst_annot[
  !is.na(corresponding_gene_call),
  .(
    mean_temporal_Fst = mean(Fst, na.rm = TRUE)
  ),
  by = .(MAG_id, corresponding_gene_call, Colony, Genus)
]

# 7.2 Mean temporal Fst per MAG
mag_mean <- gene_mag_fst[
  ,
  .(mean_MAG_temporal_Fst = mean(mean_temporal_Fst, na.rm = TRUE)),
  by = MAG_id
]

gene_mag_fst <- merge(
  gene_mag_fst,
  mag_mean,
  by = "MAG_id",
  all.x = TRUE
)

## Make sure both sides are character
gene_mag_fst[, corresponding_gene_call := as.character(corresponding_gene_call)]

parallel_zeroD_A <- fread(file.path(OUTDIR, "parallel_zeroD_A_genes.csv"))
parallel_zeroD_B <- fread(file.path(OUTDIR, "parallel_zeroD_B_genes.csv"))
parallel_zeroD_C <- fread(file.path(OUTDIR, "parallel_zeroD_C_genes.csv"))

parallel_zeroD_A[, corresponding_gene_call := as.character(corresponding_gene_call)]
parallel_zeroD_B[, corresponding_gene_call := as.character(corresponding_gene_call)]
parallel_zeroD_C[, corresponding_gene_call := as.character(corresponding_gene_call)]

parallel_zeroD_genes <- unique(c(
  parallel_zeroD_A$corresponding_gene_call,
  parallel_zeroD_B$corresponding_gene_call,
  parallel_zeroD_C$corresponding_gene_call
))

# Now this works:
gene_mag_fst[, parallel_zeroD := corresponding_gene_call %chin% parallel_zeroD_genes]

# 7.4 Restrict to genes where we know parallel-zeroD status
gene_mag_fst_use <- gene_mag_fst[!is.na(parallel_zeroD)]

# Quick summaries
gene_mag_fst_use[
  , .(
    n = .N,
    mean_centered = mean(Fst_centered, na.rm = TRUE),
    median_centered = median(Fst_centered, na.rm = TRUE),
    q25 = quantile(Fst_centered, 0.25, na.rm = TRUE),
    q75 = quantile(Fst_centered, 0.75, na.rm = TRUE)
  ),
  by = parallel_zeroD
]

# Wilcoxon: parallel vs non-parallel (centered Fst)
wilcox.test(
  Fst_centered ~ parallel_zeroD,
  data = gene_mag_fst_use
)
# Mixed model: MAG_id random effect (on already-centered Fst)
m_par0D <- lmer(
  Fst_centered ~ parallel_zeroD + (1 | MAG_id),
  data = gene_mag_fst_use
)
summary(m_par0D)

# Plot
p_parallel_centered <- ggplot(gene_mag_fst_use,
                              aes(x = parallel_zeroD, y = Fst_centered, fill = parallel_zeroD)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_classic() +
  labs(
    x = "Parallel 0D-sweep gene?",
    y = "Centered temporal Fst (gene mean – MAG mean)",
    title = "Parallel 0D genes relative to MAG-specific temporal Fst"
  ) +
  coord_cartesian(
    ylim = c(
      quantile(gene_mag_fst_use$Fst_centered, 0.01, na.rm = TRUE),
      quantile(gene_mag_fst_use$Fst_centered, 0.99, na.rm = TRUE)
    )
  )



ggsave(file.path(OUTDIR, "Fst_centered_parallel0D_vs_nonparallel.png"),
       p_parallel_centered, width = 6, height = 4, dpi = 300)






## ------------------------------------------------------------
## 8. Top 20 genes with largest mean temporal Fst relative to MAG average
## ------------------------------------------------------------

# Add some annotations from gene_clean (GO, PFAM, pN/pS, etc.)
gene_ann_sub <- unique(
  gene_clean[
    ,
    .(MAG_id, corresponding_gene_call,
      Genus, GO_primary, PFAM_primary,
      pN_pS_ratio, GO_str, pfam_str,
      unique_gene_callers_id)
  ]
)

gene_mag_fst_annot <- merge(
  gene_mag_fst,
  gene_ann_sub,
  by = c("MAG_id", "corresponding_gene_call"),
  all.x = TRUE
)

# Order by Fst_centered (largest positive = strongest temporal divergence)
top20_hot <- gene_mag_fst_annot[order(-Fst_centered)][1:20]

top20_hot[
  ,
  .(MAG_id, Colony,
    corresponding_gene_call, unique_gene_callers_id,
    Genus.y,
    mean_temporal_Fst,
    mean_MAG_temporal_Fst,
    Fst_centered,
    GO_primary, PFAM_primary,
    pN_pS_ratio)
]

fwrite(
  top20_hot,
  file.path(OUTDIR, "top20_genes_high_centered_temporal_Fst_Fst_ge0.csv")
)

# If you also want the most "conserved" over time:
top20_cold <- gene_mag_fst_annot[order(Fst_centered)][1:20]

fwrite(
  top20_cold,
  file.path(OUTDIR, "top20_genes_low_centered_temporal_Fst_Fst_ge0.csv")
)
































# Ensure types
gene_clean <- as.data.table(gene_clean)
gene_clean[, `:=`(
  MAG_id = as.character(MAG_id),
  corresponding_gene_call = as.character(corresponding_gene_call)
)]

# Helper to go from (MAGs=list) to long format
expand_parallel <- function(dt, type_label, criterion_label) {
  if (!"MAGs" %in% names(dt)) stop("dt must have a MAGs list column")
  dt[, .(
    MAG_id = unlist(MAGs),
    type    = type_label,
    criterion = criterion_label
  ), by = corresponding_gene_call]
}

# pN/pS parallel sets
pnps_A_long <- expand_parallel(parallel_pnps_A, "pnps_parallel", "A")
pnps_B_long <- expand_parallel(parallel_pnps_B, "pnps_parallel", "B")
pnps_C_long <- expand_parallel(parallel_pnps_C, "pnps_parallel", "C")

# 0D sweep parallel sets
zeroD_A_long <- expand_parallel(parallel_zeroD_A, "zeroD_parallel", "A")
zeroD_B_long <- expand_parallel(parallel_zeroD_B, "zeroD_parallel", "B")
zeroD_C_long <- expand_parallel(parallel_zeroD_C, "zeroD_parallel", "C")

parallel_long <- rbind(
  pnps_A_long, pnps_B_long, pnps_C_long,
  zeroD_A_long, zeroD_B_long, zeroD_C_long,
  use.names = TRUE, fill = TRUE
)

parallel_long[, `:=`(
  MAG_id = as.character(MAG_id),
  corresponding_gene_call = as.character(corresponding_gene_call)
)]

parallel_long <- unique(parallel_long)

# Attach gene-level annotations (Genus, GO, PFAM, etc.)
parallel_annot <- merge(
  parallel_long,
  unique(gene_clean[, .(MAG_id, corresponding_gene_call,
                        Genus, GO_primary, PFAM_primary,
                        GO_str, pfam_str, pN_pS_ratio)]),
  by = c("MAG_id","corresponding_gene_call"),
  all.x = TRUE
)

parallel_annot



















pfam_genus_counts <- parallel_annot[
  !is.na(PFAM_primary) & !is.na(Genus),
  .N,
  by = .(Genus, PFAM_primary)
]

# Keep PFAMs with at least a few parallel genes to avoid crazy sparsity
pfam_genus_counts <- pfam_genus_counts[N >= 2]

p_pfam_genus <- ggplot(pfam_genus_counts,
                       aes(x = Genus, y = PFAM_primary, fill = N)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  coord_flip() +
  labs(
    x = "Genus",
    y = "PFAM (primary)",
    fill = "# parallel genes",
    title = "Functional distribution of parallel genes (PFAM × Genus)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTDIR, "parallel_genes_PFAM_by_Genus_heatmap.png"),
       p_pfam_genus, width = 7, height = 6, dpi = 300)

p_pfam_genus









go_genus_counts <- parallel_annot[
  !is.na(GO_primary) & !is.na(Genus),
  .N,
  by = .(Genus, GO_primary)
][N >= 2]

p_go_genus <- ggplot(go_genus_counts,
                     aes(x = Genus, y = GO_primary, fill = N)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  coord_flip() +
  labs(
    x = "Genus",
    y = "GO primary term",
    fill = "# parallel genes",
    title = "Functional distribution of parallel genes (GO × Genus)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_go_genus
ggsave(file.path(OUTDIR, "parallel_genes_GO_by_Genus_heatmap.png"),
       p_go_genus, width = 7, height = 6, dpi = 300)

















# Wide-type indicators per MAG×gene
parallel_wide <- dcast(
  parallel_annot,
  MAG_id + corresponding_gene_call + Genus + GO_primary + PFAM_primary ~ type,
  fun.aggregate = length,
  value.var = "type"
)

# Turn counts into logical
for (col in c("pnps_parallel","zeroD_parallel")) {
  if (!col %in% names(parallel_wide)) parallel_wide[[col]] <- 0L
  parallel_wide[[col]] <- parallel_wide[[col]] > 0
}

parallel_wide[
  ,
  parallel_class := fifelse(pnps_parallel & zeroD_parallel, "both",
                            fifelse(pnps_parallel, "pnps_only",
                                    fifelse(zeroD_parallel, "zeroD_only","none")))
]

parallel_wide <- parallel_wide[parallel_class != "none"]

p_class <- ggplot(parallel_wide,
                  aes(x = Genus, fill = parallel_class)) +
  geom_bar(position = "fill") +
  labs(
    x = "Genus",
    y = "Proportion of parallel genes",
    fill = "Parallel type",
    title = "Relative contribution of pN/pS vs 0D sweeps to gene-level parallelism"
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "parallel_type_proportions_by_Genus.png"),
       p_class, width = 6, height = 4, dpi = 300)












# For each gene, how many distinct MAGs where it is parallel?
gene_parallel_spread <- parallel_annot[
  ,
  .(n_MAGs = uniqueN(MAG_id)),
  by = corresponding_gene_call
]

top_genes <- gene_parallel_spread[order(-n_MAGs)][1:50]$corresponding_gene_call

incidence <- unique(
  parallel_annot[corresponding_gene_call %chin% top_genes,
                 .(MAG_id, corresponding_gene_call, type)]
)

p_incidence <- ggplot(incidence,
                      aes(x = MAG_id, y = factor(corresponding_gene_call),
                          fill = type)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_manual(values = c(pnps_parallel = "#1b9e77",
                               zeroD_parallel = "#d95f02")) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "MAGs",
    y = "Parallel genes (top 50 by breadth across MAGs)",
    fill = "Parallel signal type",
    title = "Which MAGs share the same parallel genes?"
  )
p_incidence
ggsave(file.path(OUTDIR, "parallel_genes_MAG_incidence_top50.png"),
       p_incidence, width = 7, height = 6, dpi = 300)










# zeroD_dt assumed as in your sweep table, with Season_t1 / Season_t2
zeroD_dt <- as.data.table(zeroD_dt)
zeroD_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]
zeroD_dt[, MAG_id := as.character(MAG_id)]

season_map <- c("May"="Spring","June"="Spring",
                "July"="Summer","August"="Summer","September"="Summer",
                "October"="Fall","November"="Fall",
                "January"="Winter","February"="Winter")
if (!"Season_t1" %in% names(zeroD_dt) & "Month_t1" %in% names(zeroD_dt)) {
  zeroD_dt[, Season_t1 := season_map[Month_t1]]
}

# Restrict to *parallel* 0D genes and unpack to MAG×gene
parallel_zeroD_all <- rbind(parallel_zeroD_A, parallel_zeroD_B, parallel_zeroD_C, fill = TRUE)
parallel_zeroD_all[, corresponding_gene_call := as.character(corresponding_gene_call)]

zeroD_parallel_hits <- merge(
  zeroD_dt[Season_t1 %in% c("Spring","Summer","Fall","Winter")],
  unique(parallel_zeroD_all[, .(corresponding_gene_call)]),
  by = "corresponding_gene_call",
  all.x = FALSE
)

# Count parallel 0D sweeps per Season × Genus
zeroD_season_genus <- zeroD_parallel_hits[
  !is.na(Genus),
  .N,
  by = .(Genus, Season_t1)
]

zeroD_season_genus[, Season_t1 := factor(Season_t1,
                                         levels = c("Spring","Summer","Fall","Winter"),
                                         ordered = TRUE)]

p_zeroD_season <- ggplot(zeroD_season_genus,
                         aes(x = Season_t1, y = N, fill = Season_t1)) +
  geom_col() +
  facet_wrap(~ Genus, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Season of sweep (t1)",
    y = "# parallel 0D sweep events",
    title = "Seasonal distribution of parallel 0D sweeps by genus"
  )
p_zeroD_season
ggsave(file.path(OUTDIR, "parallel_zeroD_sweeps_by_season_genus.png"),
       p_zeroD_season, width = 8, height = 5, dpi = 300)


### Try to make matrix ###

target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium",
                   "commensalibacter", "frischella", "gilliamella",
                   "lactobacillus", "snodgrassella")


## ------------------------------------------------------------
## 1. Build core 0D-sweep matrix (MAG × gene × Month)
## ------------------------------------------------------------

zeroD_dt <- as.data.table(zeroD_dt)
gene_clean <- as.data.table(gene_clean)

# Make sure key IDs are character
zeroD_dt[, `:=`(
  MAG_id = as.character(MAG_id),
  corresponding_gene_call = as.character(corresponding_gene_call),
  colony_id = as.character(colony_id)
)]
gene_clean[, `:=`(
  MAG_id = as.character(MAG_id),
  corresponding_gene_call = as.character(corresponding_gene_call),
  colony_id = as.character(colony_id)
)]

# Pick a "sweep time" – here we use the *later* month (Month_t2) if present,
# otherwise fall back to whatever Month column you had
if ("Month_t2" %in% names(zeroD_dt)) {
  zeroD_dt[, Month_sweep := Month_t2]
} else if ("Month" %in% names(zeroD_dt)) {
  zeroD_dt[, Month_sweep := Month]
} else {
  stop("No Month or Month_t2 column in zeroD_dt")
}

# Order months
month_levels <- c("May","June","July","August","September","October",
                  "November","January","February")
zeroD_dt[, Month_sweep := factor(Month_sweep, levels = month_levels, ordered = TRUE)]

# Attach Genus from zeroD_dt (you already merged tax info earlier)
if (!"Genus" %in% names(zeroD_dt)) {
  stop("Genus column not found in zeroD_dt; make sure tax_dt was merged earlier")
}

# Restrict to target genera if you want core guys only
zeroD_dt_core <- zeroD_dt[Genus %chin% target_genera]

# Collapse to one record per MAG × Colony × gene × Month_sweep
zeroD_gene_month <- unique(
  zeroD_dt_core[
    ,
    .(MAG_id,
      Colony = colony_id,
      corresponding_gene_call,
      Genus,
      Month_sweep)
  ]
)

## --- Attach PFAM / GO / pN/pS info from gene_clean -------------------------

gene_ann <- unique(
  gene_clean[
    ,
    .(MAG_id, Colony = colony_id, corresponding_gene_call,
      Genus, GO_primary, PFAM_primary,
      GO_str, pfam_str, pN_pS_ratio)
  ]
)

zeroD_gene_month <- merge(
  zeroD_gene_month,
  gene_ann,
  by = c("MAG_id", "Colony", "corresponding_gene_call"),
  all.x = TRUE
)

## --- Flag pN/pS > 1 genes (within MAG×Colony×gene) ------------------------

pnps_flag <- unique(
  gene_clean[!is.na(pN_pS_ratio) & pN_pS_ratio > 1,
             .(MAG_id, Colony = colony_id, corresponding_gene_call)]
)
pnps_flag[, pnps_gt1 := TRUE]

zeroD_gene_month <- merge(
  zeroD_gene_month,
  pnps_flag,
  by = c("MAG_id","Colony","corresponding_gene_call"),
  all.x = TRUE
)
zeroD_gene_month[is.na(pnps_gt1), pnps_gt1 := FALSE]

## --- Flag "parallel zeroD" genes -------------------------------------------

parallel_zeroD_A <- fread(file.path(OUTDIR, "parallel_zeroD_A_genes.csv"))
parallel_zeroD_B <- fread(file.path(OUTDIR, "parallel_zeroD_B_genes.csv"))
parallel_zeroD_C <- fread(file.path(OUTDIR, "parallel_zeroD_C_genes.csv"))

parallel_zeroD_A[, corresponding_gene_call := as.character(corresponding_gene_call)]
parallel_zeroD_B[, corresponding_gene_call := as.character(corresponding_gene_call)]
parallel_zeroD_C[, corresponding_gene_call := as.character(corresponding_gene_call)]

parallel_zeroD_genes <- unique(c(
  parallel_zeroD_A$corresponding_gene_call,
  parallel_zeroD_B$corresponding_gene_call,
  parallel_zeroD_C$corresponding_gene_call
))

zeroD_gene_month[, parallel_zeroD := corresponding_gene_call %chin% parallel_zeroD_genes]

zeroD_gene_month[]
# columns: MAG_id, Colony, corresponding_gene_call, Genus, Month_sweep,
#          GO_primary, PFAM_primary, GO_str, pfam_str, pN_pS_ratio,
#          pnps_gt1, parallel_zeroD




## ------------------------------------------------------------
## PFAM × Month × Genus matrix
## ------------------------------------------------------------

pfam_mat <- zeroD_gene_month[!is.na(PFAM_primary)]

pfam_summary <- pfam_mat[
  ,
  .(
    n_hits = .N,                                  # number of 0D-sweep genes in this PFAM × Month × Genus
    n_MAGs = uniqueN(MAG_id),                     # how many MAGs contribute
    any_parallel_zeroD = any(parallel_zeroD),
    any_pnps_gt1      = any(pnps_gt1)
  ),
  by = .(Genus.x, PFAM_primary, Month_sweep)
]

pfam_summary[, Month_sweep := droplevels(Month_sweep)]
pfam_summary


p_pfam_mat <- ggplot(pfam_summary,
                     aes(x = Month_sweep, y = PFAM_primary, fill = n_MAGs)) +
  geom_tile(color = "grey80") +
  # overlay markers: triangle for pN/pS>1, circle border for parallel_zeroD
  geom_point(
    data = pfam_summary[any_pnps_gt1 == TRUE],
    aes(x = Month_sweep, y = PFAM_primary),
    shape = 24, size = 2, fill = "black", color = "black"
  ) +
  geom_point(
    data = pfam_summary[any_parallel_zeroD == TRUE],
    aes(x = Month_sweep, y = PFAM_primary),
    shape = 21, size = 2, fill = NA, color = "red", stroke = 0.6
  ) +
  scale_fill_viridis_c(option = "plasma", direction = 1) +
  facet_wrap(~ Genus.x, scales = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    x = "Month of 0D sweep",
    y = "PFAM (primary)",
    fill = "# MAGs with 0D sweep",
    title = "0D sweep matrix: PFAM × Month (faceted by Genus)",
    subtitle = "Red outline = parallel 0D gene(s); black triangle = any pN/pS > 1"
  )
p_pfam_mat
ggsave(file.path(OUTDIR, "matrix_0D_PFAM_by_Month_faceted_Genus.png"),
       p_pfam_mat, width = 8, height = 6, dpi = 300)












## ------------------------------------------------------------
## Gene × Month × Genus matrix (parallel 0D only)
## ------------------------------------------------------------

gene_mat <- zeroD_gene_month[parallel_zeroD == TRUE]

# optional: keep only genes with at least 2 MAGs or 3 months etc
gene_filter <- gene_mat[
  ,
  .(
    n_MAGs  = uniqueN(MAG_id),
    n_month = uniqueN(Month_sweep)
  ),
  by = corresponding_gene_call
]

keep_genes <- gene_filter[n_MAGs >= 2 | n_month >= 3]$corresponding_gene_call

gene_mat_filt <- gene_mat[corresponding_gene_call %chin% keep_genes]

gene_summary <- gene_mat_filt[
  ,
  .(
    n_hits = .N,
    n_MAGs = uniqueN(MAG_id),
    any_pnps_gt1 = any(pnps_gt1)
  ),
  by = .(Genus.x, corresponding_gene_call, Month_sweep)
]

gene_summary[, Month_sweep := droplevels(Month_sweep)]

p_gene_mat <- ggplot(gene_summary,
                     aes(x = Month_sweep,
                         y = factor(corresponding_gene_call),
                         fill = n_MAGs)) +
  geom_tile(color = "grey85") +
  geom_point(
    data = gene_summary[any_pnps_gt1 == TRUE],
    aes(x = Month_sweep,
        y = factor(corresponding_gene_call)),
    shape = 24, size = 2, fill = "black", color = "black"
  ) +
  scale_fill_viridis_c(option = "magma", direction = 1) +
  facet_wrap(~ Genus.x, scales = "free_y") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5)
  ) +
  labs(
    x = "Month of 0D sweep",
    y = "corresponding_gene_call",
    fill = "# MAGs",
    title = "Parallel 0D sweeps across genes and months",
    subtitle = "Triangles mark genes with pN/pS > 1 in that MAG/Colony"
  )
p_gene_mat
ggsave(file.path(OUTDIR, "matrix_0D_gene_by_Month_parallel_only.png"),
       p_gene_mat, width = 8, height = 6, dpi = 300)














# target genera as before
target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium",
                   "commensalibacter", "frischella", "gilliamella",
                   "lactobacillus", "snodgrassella")

fst_vs0_season <- as.data.table(fst_vs0_season)

# Keep only target genera and PFAM with non-NA
fst_pfam <- fst_vs0_season[
  !is.na(PFAM_primary) & !is.na(Fst) & !is.na(Season) &
    !is.na(MAG_id) & Genus %chin% target_genera
]

# Make sure Season is ordered Spring → Winter
season_levels <- c("Spring","Summer","Fall","Winter")
fst_pfam[, Season := factor(Season, levels = season_levels, ordered = TRUE)]
fst_pfam[, MAG_id := as.character(MAG_id)]


## ------------------------------------------------------------
## PFAM-level concordance: coherent seasonal response vs MAG heterogeneity
## ------------------------------------------------------------

pfam_concordance <- fst_pfam[
  ,
  {
    # cell means: one Fst per MAG × Season (within PFAM)
    cell <- .SD[, .(Fst_cell = mean(Fst, na.rm = TRUE)),
                by = .(MAG_id, Season)]
    
    if (nrow(cell) < 4) {
      # too little data to estimate anything meaningful
      list(
        n_obs   = .N,
        n_MAG   = uniqueN(MAG_id),
        n_season = uniqueN(Season),
        between_season_var = NA_real_,
        within_season_var  = NA_real_,
        concordance_index  = NA_real_
      )
    } else {
      # between-season variance of season means
      season_means <- cell[, .(mean_Fst = mean(Fst_cell, na.rm = TRUE)), by = Season]
      between_season_var <- var(season_means$mean_Fst, na.rm = TRUE)
      
      # within-season variance across MAGs, averaged over seasons
      within_season_var <- cell[
        ,
        .(var_within = var(Fst_cell, na.rm = TRUE)),
        by = Season
      ][, mean(var_within, na.rm = TRUE)]
      
      # concordance index
      denom <- between_season_var + within_season_var
      concordance_index <- if (is.finite(denom) && denom > 0) {
        between_season_var / denom
      } else NA_real_
      
      list(
        n_obs   = .N,
        n_MAG   = uniqueN(MAG_id),
        n_season = uniqueN(Season),
        between_season_var = between_season_var,
        within_season_var  = within_season_var,
        concordance_index  = concordance_index
      )
    }
  },
  by = PFAM_primary
]

# filter to reasonably well-sampled pathways
pfam_concordance_filt <- pfam_concordance[
  n_obs >= 50 & n_MAG >= 3 & n_season >= 2 & !is.na(concordance_index)
]

pfam_concordance_filt[order(-concordance_index)][1:20]


p_pfam_conc <- ggplot(pfam_concordance_filt,
                      aes(x = between_season_var,
                          y = within_season_var,
                          color = concordance_index)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    x = "Between-season variance in mean Fst",
    y = "Within-season variance across MAGs",
    color = "Concordance\nindex",
    title = "Pathway-level coherence of seasonal Fst response (PFAM)",
    subtitle = "High index: MAGs behave similarly across seasons in that pathway"
  ) +
  theme_classic(base_size = 10)
p_pfam_conc
ggsave(file.path(OUTDIR, "PFAM_concordance_season_vs_MAG_Fst.png"),
       p_pfam_conc, width = 6, height = 5, dpi = 300)



top_pfam <- pfam_concordance_filt[order(-concordance_index)][1:10]
bottom_pfam <- pfam_concordance_filt[order(concordance_index)][1:10]

top_pfam[, .(PFAM_primary, n_MAG, n_season, concordance_index)]
bottom_pfam[, .(PFAM_primary, n_MAG, n_season, concordance_index)]


top15 <- pfam_concordance_filt[order(-concordance_index)][1:15]

p_top15 <- ggplot(top15,
                  aes(x = reorder(PFAM_primary, concordance_index),
                      y = concordance_index)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "PFAM (primary)",
    y = "Seasonal concordance index",
    title = "Top PFAM pathways with coherent seasonal Fst responses"
  ) +
  theme_classic(base_size = 9)
p_top15
ggsave(file.path(OUTDIR, "PFAM_top15_seasonal_concordance.png"),
       p_top15, width = 6, height = 5, dpi = 300)





fst_go <- fst_vs0_season[
  !is.na(GO_primary) & !is.na(Fst) & !is.na(Season) &
    !is.na(MAG_id) & Genus %chin% target_genera
]

fst_go[, Season := factor(Season, levels = season_levels, ordered = TRUE)]
fst_go[, MAG_id := as.character(MAG_id)]

go_concordance <- fst_go[
  ,
  {
    cell <- .SD[, .(Fst_cell = mean(Fst, na.rm = TRUE)),
                by = .(MAG_id, Season)]
    if (nrow(cell) < 4) {
      list(n_obs=.N, n_MAG=uniqueN(MAG_id), n_season=uniqueN(Season),
           between_season_var=NA_real_, within_season_var=NA_real_,
           concordance_index=NA_real_)
    } else {
      season_means <- cell[, .(mean_Fst = mean(Fst_cell, na.rm = TRUE)),
                           by = Season]
      between_season_var <- var(season_means$mean_Fst, na.rm = TRUE)
      within_season_var <- cell[
        ,
        .(var_within = var(Fst_cell, na.rm = TRUE)),
        by = Season
      ][, mean(var_within, na.rm = TRUE)]
      denom <- between_season_var + within_season_var
      concordance_index <- if (is.finite(denom) && denom > 0) {
        between_season_var / denom
      } else NA_real_
      list(n_obs=.N, n_MAG=uniqueN(MAG_id), n_season=uniqueN(Season),
           between_season_var=between_season_var,
           within_season_var=within_season_var,
           concordance_index=concordance_index)
    }
  },
  by = GO_primary
]

go_concordance_filt <- go_concordance[
  n_obs >= 50 & n_MAG >= 3 & n_season >= 2 & !is.na(concordance_index)
]
go_concordance_filt[order(-concordance_index)][1:20]

# Pathways with high concordance
high_pfam <- pfam_concordance_filt[concordance_index > 0.7]$PFAM_primary

# Among 0D sweep genes, how many in high-concordance vs others?
zeroD_gene_month[, high_concordance_PFAM := PFAM_primary %chin% high_pfam]

table(zeroD_gene_month$high_concordance_PFAM, zeroD_gene_month$parallel_zeroD)










season_map <- c("May"="Spring","June"="Spring",
                "July"="Summer","August"="Summer","September"="Summer",
                "October"="Fall","November"="Fall",
                "January"="Winter","February"="Winter")
season_levels <- c("Spring","Summer","Fall","Winter")

target_genera <- c("apilactobacillus", "bartonella", "bifidobacterium",
                   "commensalibacter", "frischella", "gilliamella",
                   "lactobacillus", "snodgrassella")



fst_vs0_season <- as.data.table(fst_vs0_season)
fst_vs0_season[, Season := factor(Season, levels = season_levels, ordered = TRUE)]
fst_vs0_season[, MAG_id := as.character(MAG_id)]

fst_pfam <- fst_vs0_season[
  !is.na(PFAM_primary) & !is.na(Fst) & !is.na(Season) &
    !is.na(MAG_id) & !is.na(Genus) & Genus %chin% target_genera
]

pfam_Fst_cells <- fst_pfam[
  ,
  .(Fst_cell = mean(Fst, na.rm = TRUE)),
  by = .(PFAM_primary, MAG_id, Season)
]



gene_clean <- as.data.table(gene_clean)
gene_clean[, `:=`(
  MAG_id = as.character(MAG_id),
  Colony = as.character(colony_id),
  Season = season_map[Month]
)]
gene_clean[, Season := factor(Season, levels = season_levels, ordered = TRUE)]

gene_pfam <- gene_clean[
  !is.na(PFAM_primary) & !is.na(Season) & !is.na(MAG_id) &
    !is.na(Genus) & Genus %chin% target_genera
]

# (optional) tame insane pN/pS outliers:
# gene_pfam[pN_pS_ratio <= 0 | pN_pS_ratio > 50, pN_pS_ratio := NA]

pfam_gene_cells <- gene_pfam[
  ,
  .(
    n_gene_obs = .N,
    mean_pNpS  = mean(pN_pS_ratio, na.rm = TRUE),
    mean_piN   = mean(pi_nonsyn,    na.rm = TRUE)
  ),
  by = .(PFAM_primary, MAG_id, Season)
]



zeroD_dt <- as.data.table(zeroD_dt)
zeroD_dt[, `:=`(
  MAG_id = as.character(MAG_id),
  Colony = as.character(colony_id),
  corresponding_gene_call = as.character(corresponding_gene_call)
)]

# decide which month is "sweep month"
if ("Month_t2" %in% names(zeroD_dt)) {
  zeroD_dt[, Month_sweep := Month_t2]
} else if ("Month" %in% names(zeroD_dt)) {
  zeroD_dt[, Month_sweep := Month]
} else {
  stop("Need Month or Month_t2 in zeroD_dt for sweep timing.")
}

month_levels <- c("May","June","July","August","September","October",
                  "November","January","February")
zeroD_dt[, Month_sweep := factor(Month_sweep, levels = month_levels, ordered = TRUE)]

zeroD_gene_month <- unique(
  zeroD_dt[
    Genus %chin% target_genera &
      !is.na(pfam_str),
    .(MAG_id, Colony, corresponding_gene_call, Genus, Month_sweep, pfam_str)
  ]
)
zeroD_dt
zeroD_gene_month[, Season := season_map[as.character(Month_sweep)]]
zeroD_gene_month[, Season := factor(Season, levels = season_levels, ordered = TRUE)]

pfam_sweep_cells <- zeroD_gene_month[
  !is.na(Season),
  .(
    n_sweep_genes = uniqueN(corresponding_gene_call)
  ),
  by = .(pfam_str, MAG_id, Season)
]

pfam_sweep_cells
pfam_sweep_cells <- rename(pfam_sweep_cells, PFAM_primary = PFAM_PRIMARY)

pfam_cells <- Reduce(
  function(x, y) merge(x, y,
                       by = c("PFAM_primary", "MAG_id", "Season"),
                       all = TRUE),
  list(pfam_Fst_cells, pfam_gene_cells, pfam_sweep_cells)
)

# PFAM_primary, MAG_id, Season,
# Fst_cell, n_gene_obs, mean_pNpS, mean_piN, n_sweep_genes







compute_concordance <- function(dt, value_col) {
  v <- dt[[value_col]]
  if (all(is.na(v)) || length(na.omit(v)) < 4) {
    return(list(
      between_var = NA_real_,
      within_var  = NA_real_,
      concordance = NA_real_,
      mean_val    = NA_real_
    ))
  }
  
  # one metric value per MAG × Season
  cell <- dt[!is.na(get(value_col)) &
               !is.na(Season) &
               !is.na(MAG_id),
             .(value = mean(get(value_col), na.rm = TRUE)),
             by = .(MAG_id, Season)]
  
  if (nrow(cell) < 4 || uniqueN(cell$Season) < 2 || uniqueN(cell$MAG_id) < 2) {
    return(list(
      between_var = NA_real_,
      within_var  = NA_real_,
      concordance = NA_real_,
      mean_val    = mean(v, na.rm = TRUE)
    ))
  }
  
  # between-season variance of season means
  season_means <- cell[, .(mean_val = mean(value, na.rm = TRUE)), by = Season]
  between_var <- var(season_means$mean_val, na.rm = TRUE)
  
  # within-season variance across MAGs, averaged
  within_var <- cell[
    ,
    .(var_within = var(value, na.rm = TRUE)),
    by = Season
  ][, mean(var_within, na.rm = TRUE)]
  
  denom <- between_var + within_var
  concordance <- if (is.finite(denom) && denom > 0) between_var / denom else NA_real_
  
  list(
    between_var = between_var,
    within_var  = within_var,
    concordance = concordance,
    mean_val    = mean(cell$value, na.rm = TRUE)
  )
}






pfam_conc_all <- pfam_cells[
  ,
  {
    base_info <- list(
      n_rows   = .N,
      n_MAG    = uniqueN(MAG_id),
      n_season = uniqueN(Season)
    )
    
    c_Fst   <- compute_concordance(.SD, "Fst_cell")
    c_pNpS  <- compute_concordance(.SD, "mean_pNpS")
    c_piN   <- compute_concordance(.SD, "mean_piN")
    c_sweep <- compute_concordance(.SD, "n_sweep_genes")
    
    c(
      base_info,
      setNames(c_Fst,   paste0("Fst_",   names(c_Fst))),
      setNames(c_pNpS,  paste0("pNpS_",  names(c_pNpS))),
      setNames(c_piN,   paste0("piN_",   names(c_piN))),
      setNames(c_sweep, paste0("sweep_", names(c_sweep)))
    )
  },
  by = PFAM_primary
]

# filter to sensibly sampled PFAMs
pfam_conc_filt <- pfam_conc_all[
  n_rows >= 50 & n_MAG >= 3 & n_season >= 2
]

pfam_conc_filt[1:10]


## --- build GO cells -------------------------------------------

fst_go <- fst_vs0_season[
  !is.na(GO_primary) & !is.na(Fst) & !is.na(Season) &
    !is.na(MAG_id) & !is.na(Genus) & Genus %chin% target_genera
]

fst_go[, Season := factor(Season, levels = season_levels, ordered = TRUE)]
fst_go[, MAG_id := as.character(MAG_id)]

go_Fst_cells <- fst_go[
  ,
  .(Fst_cell = mean(Fst, na.rm = TRUE)),
  by = .(GO_primary, MAG_id, Season)
]

gene_go <- gene_clean[
  !is.na(GO_primary) & !is.na(Season) & !is.na(MAG_id) &
    !is.na(Genus) & Genus %chin% target_genera
]

go_gene_cells <- gene_go[
  ,
  .(
    n_gene_obs = .N,
    mean_pNpS  = mean(pN_pS_ratio, na.rm = TRUE),
    mean_piN   = mean(pi_nonsyn,    na.rm = TRUE)
  ),
  by = .(GO_primary, MAG_id, Season)
]

zeroD_gene_month_go <- unique(
  zeroD_dt[
    Genus %chin% target_genera &
      !is.na(GO_str),
    .(MAG_id, Colony, corresponding_gene_call, Genus, Month_sweep, GO_str)
  ]
)
zeroD_gene_month_go[, Season := season_map[as.character(Month_sweep)]]
zeroD_gene_month_go[, Season := factor(Season, levels = season_levels, ordered = TRUE)]

zeroD_gene_month_go <- rename(zeroD_gene_month_go, GO_primary = GO_str)

go_sweep_cells <- zeroD_gene_month_go[
  !is.na(Season),
  .(
    n_sweep_genes = uniqueN(corresponding_gene_call)
  ),
  by = .(GO_primary, MAG_id, Season)
]

go_cells <- Reduce(
  function(x, y) merge(x, y,
                       by = c("GO_primary", "MAG_id", "Season"),
                       all = TRUE),
  list(go_Fst_cells, go_gene_cells, go_sweep_cells)
)

## --- GO concordance per metric --------------------------------

go_conc_all <- go_cells[
  ,
  {
    base_info <- list(
      n_rows   = .N,
      n_MAG    = uniqueN(MAG_id),
      n_season = uniqueN(Season)
    )
    
    c_Fst   <- compute_concordance(.SD, "Fst_cell")
    c_pNpS  <- compute_concordance(.SD, "mean_pNpS")
    c_piN   <- compute_concordance(.SD, "mean_piN")
    c_sweep <- compute_concordance(.SD, "n_sweep_genes")
    
    c(
      base_info,
      setNames(c_Fst,   paste0("Fst_",   names(c_Fst))),
      setNames(c_pNpS,  paste0("pNpS_",  names(c_pNpS))),
      setNames(c_piN,   paste0("piN_",   names(c_piN))),
      setNames(c_sweep, paste0("sweep_", names(c_sweep)))
    )
  },
  by = GO_primary
]

go_conc_filt <- go_conc_all[
  n_rows >= 50 & n_MAG >= 3 & n_season >= 2
]

go_conc_filt[1:10]







top50_pf_Fst <- pfam_conc_filt[
  !is.na(Fst_concordance)
][order(-Fst_concordance)][1:50]

ggplot(top50_pf_Fst,
       aes(x = reorder(PFAM_primary, Fst_concordance),
           y = Fst_concordance)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "PFAM",
    y = "Seasonal concordance (Fst)",
    title = "Top 50 PFAM pathways by Fst concordance"
  ) +
  theme_classic(base_size = 9)










