library(data.table)


# ============================================================
# Paths
# ============================================================

ROOT <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel/parallel_with_CLUSTERED_PROTEINS/"

CLUSTER_FILE <- file.path(
  ROOT,
  "selection_gene_clusters",
  "eligible_gene_mmseqs_clusters.csv.gz"
)

GENE_SUM_FILE <- file.path(
  ROOT,
  "master_per_gene_summary.csv"
)

SWEEP_FILE <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/sweep_gene_merged.csv"

TAX_FILE <- "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv"

OUTDIR <- file.path(
  ROOT,
  "selection_gene_clusters",
  "selection_flags"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Helpers
# ============================================================

clean_colony <- function(x) {
  x <- as.character(x)
  x <- sub("^colony_", "", x)
  x <- sub("\\.0$", "", x)
  x
}


# ============================================================
# 1. Load MMseq physical-gene table
# ============================================================

message("Loading MMseq physical-gene table...")

cluster_dt <- fread(CLUSTER_FILE)

cluster_dt[, MAG_id := as.character(MAG_id)]
cluster_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]

# physical-gene identity
cluster_dt[, gene_key := paste(
  MAG_id,
  corresponding_gene_call,
  sep = "::"
)]

stopifnot(
  uniqueN(cluster_dt$gene_key) == nrow(cluster_dt)
)

message(
  "Eligible physical genes: ",
  format(nrow(cluster_dt), big.mark = ",")
)

message(
  "MMseq clusters: ",
  format(uniqueN(cluster_dt$cluster_id), big.mark = ",")
)


# ============================================================
# 2. Attach canonical Genus
# ============================================================

message("Attaching canonical taxonomy...")

tax_dt <- fread(
  TAX_FILE,
  select = c("MAG_id", "Genus")
)

tax_dt[, MAG_id := as.character(MAG_id)]
tax_dt[, Genus := as.character(Genus)]

tax_unique <- unique(
  tax_dt[, .(MAG_id, Genus)]
)

# sanity check: one genus per MAG
tax_conflicts <- tax_unique[
  !is.na(Genus),
  .(n_Genus = uniqueN(Genus)),
  by = MAG_id
][n_Genus > 1]

if (nrow(tax_conflicts) > 0) {
  fwrite(
    tax_conflicts,
    file.path(
      OUTDIR,
      "taxonomy_conflicts.csv"
    )
  )
  
  stop(
    "Some MAGs map to more than one Genus."
  )
}

tax_unique <- tax_unique[
  ,
  .(Genus = Genus[1]),
  by = MAG_id
]

cluster_dt <- merge(
  cluster_dt,
  tax_unique,
  by = "MAG_id",
  all.x = TRUE
)

message(
  "Genes with genus assignment: ",
  sum(!is.na(cluster_dt$Genus)),
  " / ",
  nrow(cluster_dt)
)




# ============================================================
# 3. Load gene summary
# ============================================================

message("Loading gene summary...")

gene_dt <- fread(
  GENE_SUM_FILE,
  select = c(
    "MAG_id",
    "colony_id",
    "Month",
    "corresponding_gene_call",
    "pN_pS_ratio"
  )
)

gene_dt[, MAG_id := as.character(MAG_id)]
gene_dt[, colony_id := clean_colony(colony_id)]
gene_dt[, Month := as.character(Month)]
gene_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]

gene_dt[, gene_key := paste(
  MAG_id,
  corresponding_gene_call,
  sep = "::"
)]

# restrict to final eligible MMseq physical genes
eligible_keys <- cluster_dt$gene_key

gene_dt <- gene_dt[
  gene_key %chin% eligible_keys
]


# ============================================================
# 4. pN/pS > 1 at MAG x colony x gene level
# ============================================================

message("Building pN/pS population flags...")

gene_dt[, pnps_gt1_month :=
          !is.na(pN_pS_ratio) &
          pN_pS_ratio > 1
]

pnps_pop <- gene_dt[
  ,
  .(
    n_months_observed = uniqueN(Month),
    
    n_months_with_pnps = sum(
      !is.na(pN_pS_ratio)
    ),
    
    n_months_pnps_gt1 = sum(
      pnps_gt1_month
    ),
    
    any_pnps_gt1 = any(
      pnps_gt1_month
    ),
    
    max_pnps_ratio = {
      x <- pN_pS_ratio[
        !is.na(pN_pS_ratio)
      ]
      
      if (length(x)) {
        max(x)
      } else {
        NA_real_
      }
    },
    
    median_pnps_ratio = {
      x <- pN_pS_ratio[
        !is.na(pN_pS_ratio)
      ]
      
      if (length(x)) {
        median(x)
      } else {
        NA_real_
      }
    }
  ),
  by = .(
    gene_key,
    MAG_id,
    colony_id,
    corresponding_gene_call
  )
]

message(
  "Gene-populations: ",
  format(nrow(pnps_pop), big.mark = ",")
)

message(
  "Gene-populations with pN/pS > 1: ",
  format(
    sum(pnps_pop$any_pnps_gt1),
    big.mark = ","
  )
)


# ============================================================
# 5. Load sweep table and isolate 0D sweeps
# ============================================================

message("Loading sweep table...")

sweep_dt <- fread(SWEEP_FILE)

sweep_dt[, MAG_id := as.character(MAG_id)]
sweep_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]

if ("colony_id" %in% names(sweep_dt)) {
  
  sweep_dt[, colony_id := clean_colony(
    colony_id
  )]
  
} else if ("Colony" %in% names(sweep_dt)) {
  
  sweep_dt[, colony_id := clean_colony(
    Colony
  )]
  
} else {
  
  stop(
    "Sweep table contains neither colony_id nor Colony."
  )
}

sweep_dt[, gene_key := paste(
  MAG_id,
  corresponding_gene_call,
  sep = "::"
)]

zeroD_dt <- sweep_dt[
  degeneracy == "0D"
]

n_zeroD_all <- nrow(zeroD_dt)

zeroD_dt <- zeroD_dt[
  gene_key %chin% eligible_keys
]

message(
  "0D sweep rows total: ",
  format(n_zeroD_all, big.mark = ",")
)

message(
  "0D sweep rows in eligible MMseq universe: ",
  format(nrow(zeroD_dt), big.mark = ",")
)

message(
  "Match fraction: ",
  round(
    nrow(zeroD_dt) / n_zeroD_all,
    4
  )
)


# ============================================================
# 6. Build population-level 0D sweep support
# ============================================================

message("Building 0D population flags...")

# Preserve timing information where available.
if (
  all(
    c(
      "Month_t1",
      "Month_t2"
    ) %in% names(zeroD_dt)
  )
) {
  
  zeroD_dt[, sweep_timepair :=
             paste(
               Month_t1,
               Month_t2,
               sep = "->"
             )
  ]
  
} else if (
  all(
    c(
      "t1",
      "t2"
    ) %in% names(zeroD_dt)
  )
) {
  
  zeroD_dt[, sweep_timepair :=
             paste(
               t1,
               t2,
               sep = "->"
             )
  ]
  
} else {
  
  zeroD_dt[, sweep_timepair := NA_character_]
  
}


zeroD_pop <- zeroD_dt[
  ,
  .(
    has_0D_sweep = TRUE,
    
    n_0D_sweep_rows = .N,
    
    n_0D_sweep_timepairs = uniqueN(
      sweep_timepair,
      na.rm = TRUE
    )
  ),
  by = .(
    gene_key,
    MAG_id,
    colony_id,
    corresponding_gene_call
  )
]

message(
  "Gene-populations with >=1 0D sweep: ",
  format(nrow(zeroD_pop), big.mark = ",")
)


# ============================================================
# 7. Combine pN/pS and sweep support
# ============================================================

message("Combining population-level selection signals...")

gene_pop <- merge(
  pnps_pop,
  zeroD_pop,
  by = c(
    "gene_key",
    "MAG_id",
    "colony_id",
    "corresponding_gene_call"
  ),
  all.x = TRUE
)

gene_pop[
  is.na(has_0D_sweep),
  has_0D_sweep := FALSE
]

gene_pop[
  is.na(n_0D_sweep_rows),
  n_0D_sweep_rows := 0L
]

gene_pop[
  is.na(n_0D_sweep_timepairs),
  n_0D_sweep_timepairs := 0L
]


# attach cluster + genus
cluster_sub <- unique(
  cluster_dt[
    ,
    .(
      gene_key,
      cluster_id,
      Genus
    )
  ]
)

gene_pop <- merge(
  gene_pop,
  cluster_sub,
  by = "gene_key",
  all.x = TRUE
)

stopifnot(
  !any(is.na(gene_pop$cluster_id))
)


# ============================================================
# 8. Save population-level table
# ============================================================

fwrite(
  gene_pop,
  file.path(
    OUTDIR,
    "gene_population_selection_flags.csv.gz"
  )
)


# ============================================================
# 9. Collapse back to one row per physical gene
# ============================================================

message("Collapsing to physical-gene support...")

gene_support <- gene_pop[
  ,
  .(
    n_populations_observed =
      uniqueN(colony_id),
    
    n_populations_with_pnps_data =
      sum(
        n_months_with_pnps > 0
      ),
    
    n_populations_pnps_gt1 =
      sum(any_pnps_gt1),
    
    any_pnps_gt1 =
      any(any_pnps_gt1),
    
    max_pnps_ratio = {
      x <- max_pnps_ratio[
        !is.na(max_pnps_ratio)
      ]
      
      if (length(x)) {
        max(x)
      } else {
        NA_real_
      }
    },
    
    n_populations_with_0D_sweep =
      sum(has_0D_sweep),
    
    any_0D_sweep =
      any(has_0D_sweep),
    
    total_0D_sweep_rows =
      sum(n_0D_sweep_rows),
    
    total_0D_sweep_timepairs =
      sum(n_0D_sweep_timepairs)
  ),
  by = gene_key
]


gene_flags <- merge(
  cluster_dt,
  gene_support,
  by = "gene_key",
  all.x = TRUE
)


# explicit zeros / FALSE
count_cols <- c(
  "n_populations_observed",
  "n_populations_with_pnps_data",
  "n_populations_pnps_gt1",
  "n_populations_with_0D_sweep",
  "total_0D_sweep_rows",
  "total_0D_sweep_timepairs"
)

for (cc in count_cols) {
  
  set(
    gene_flags,
    which(is.na(gene_flags[[cc]])),
    cc,
    0L
  )
}


gene_flags[
  is.na(any_pnps_gt1),
  any_pnps_gt1 := FALSE
]

gene_flags[
  is.na(any_0D_sweep),
  any_0D_sweep := FALSE
]


# ============================================================
# 10. Save physical-gene table
# ============================================================

fwrite(
  gene_flags,
  file.path(
    OUTDIR,
    "gene_cluster_selection_flags.csv.gz"
  )
)


# ============================================================
# 11. QC summary
# ============================================================

message("")
message("===== QC SUMMARY =====")

message(
  "Eligible physical genes: ",
  format(nrow(gene_flags), big.mark = ",")
)

message(
  "Physical genes with pN/pS > 1: ",
  format(
    sum(gene_flags$any_pnps_gt1),
    big.mark = ","
  )
)

message(
  "Physical genes with >=1 0D sweep: ",
  format(
    sum(gene_flags$any_0D_sweep),
    big.mark = ","
  )
)

message(
  "Physical genes with both signals: ",
  format(
    sum(
      gene_flags$any_pnps_gt1 &
        gene_flags$any_0D_sweep
    ),
    big.mark = ","
  )
)

message(
  "Gene-populations with pN/pS > 1: ",
  format(
    sum(gene_pop$any_pnps_gt1),
    big.mark = ","
  )
)

message(
  "Gene-populations with >=1 0D sweep: ",
  format(
    sum(gene_pop$has_0D_sweep),
    big.mark = ","
  )
)

message("Done.")





# How many populations support each signal per physical gene?

table(gene_flags$n_populations_pnps_gt1)

table(gene_flags$n_populations_with_0D_sweep)

# Only among flagged genes

table(
  gene_flags[
    any_pnps_gt1 == TRUE,
    n_populations_pnps_gt1
  ]
)

table(
  gene_flags[
    any_0D_sweep == TRUE,
    n_populations_with_0D_sweep
  ]
)



cluster_selection_summary <- gene_flags[
  ,
  .(
    n_genes = .N,
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(Genus, na.rm = TRUE),
    
    n_pnps_gt1_genes = sum(any_pnps_gt1),
    n_0D_sweep_genes = sum(any_0D_sweep),
    n_both_genes = sum(
      any_pnps_gt1 & any_0D_sweep
    ),
    
    prop_pnps_gt1 = mean(any_pnps_gt1),
    prop_0D_sweep = mean(any_0D_sweep)
  ),
  by = cluster_id
]

cluster_selection_summary[
  order(-n_0D_sweep_genes)
][1:30]

cluster_selection_summary[
  order(-n_pnps_gt1_genes)
][1:30]

















# ============================================================
# NEXT STEP:
# Independent support for selection signals by MMseq cluster
# ============================================================

# Start from the population-level table because it preserves colony.
gene_pop

# ------------------------------------------------------------
# 1. 0D sweep support
# ------------------------------------------------------------

cluster_0D_support <- gene_pop[
  has_0D_sweep == TRUE,
  .(
    n_0D_genes = uniqueN(gene_key),
    
    n_0D_MAGs = uniqueN(MAG_id),
    
    n_0D_colonies = uniqueN(colony_id),
    
    n_0D_populations = uniqueN(
      paste(MAG_id, colony_id, sep = "::")
    ),
    
    n_0D_genera = uniqueN(
      Genus[!is.na(Genus)]
    )
  ),
  by = cluster_id
]


# ------------------------------------------------------------
# 2. pN/pS > 1 support
# ------------------------------------------------------------

cluster_pnps_support <- gene_pop[
  any_pnps_gt1 == TRUE,
  .(
    n_pnps_genes = uniqueN(gene_key),
    
    n_pnps_MAGs = uniqueN(MAG_id),
    
    n_pnps_colonies = uniqueN(colony_id),
    
    n_pnps_populations = uniqueN(
      paste(MAG_id, colony_id, sep = "::")
    ),
    
    n_pnps_genera = uniqueN(
      Genus[!is.na(Genus)]
    )
  ),
  by = cluster_id
]


# ------------------------------------------------------------
# 3. Eligible background for each cluster
# ------------------------------------------------------------

cluster_background <- gene_pop[
  ,
  .(
    n_eligible_genes = uniqueN(gene_key),
    
    n_eligible_MAGs = uniqueN(MAG_id),
    
    n_eligible_colonies = uniqueN(colony_id),
    
    n_eligible_populations = uniqueN(
      paste(MAG_id, colony_id, sep = "::")
    ),
    
    n_eligible_genera = uniqueN(
      Genus[!is.na(Genus)]
    )
  ),
  by = cluster_id
]


# ------------------------------------------------------------
# 4. Merge
# ------------------------------------------------------------

cluster_support <- merge(
  cluster_background,
  cluster_0D_support,
  by = "cluster_id",
  all.x = TRUE
)

cluster_support <- merge(
  cluster_support,
  cluster_pnps_support,
  by = "cluster_id",
  all.x = TRUE
)


# Replace missing support counts with zero
count_cols <- setdiff(
  names(cluster_support),
  "cluster_id"
)

for (cc in count_cols) {
  
  set(
    cluster_support,
    which(is.na(cluster_support[[cc]])),
    cc,
    0L
  )
}


# ------------------------------------------------------------
# 5. Simple descriptive proportions
# ------------------------------------------------------------

cluster_support[
  ,
  prop_genes_0D :=
    n_0D_genes / n_eligible_genes
]

cluster_support[
  ,
  prop_genes_pnps :=
    n_pnps_genes / n_eligible_genes
]


# ============================================================
# Inspect strongest recurrent 0D families
# ============================================================

cluster_support[
  n_0D_genes > 0,
  .(
    cluster_id,
    n_eligible_genes,
    n_eligible_MAGs,
    n_eligible_genera,
    
    n_0D_genes,
    n_0D_MAGs,
    n_0D_colonies,
    n_0D_populations,
    n_0D_genera,
    
    prop_genes_0D,
    
    n_pnps_genes,
    n_pnps_MAGs,
    n_pnps_genera
  )
][
  order(
    -n_0D_genera,
    -n_0D_MAGs,
    -n_0D_genes
  )
][1:30]







# ============================================================
# Protein-family recurrence of selection signals across taxa
# ============================================================

# gene_pop already contains:
# gene_key
# MAG_id
# colony_id
# corresponding_gene_call
# cluster_id
# Genus
# any_pnps_gt1
# has_0D_sweep


# ------------------------------------------------------------
# 1. One row per physical gene / taxon
# ------------------------------------------------------------

cluster_gene_selection <- unique(
  gene_pop[
    ,
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      any_pnps_gt1,
      has_0D_sweep
    )
  ]
)

cluster_gene_selection[
  ,
  any_selection_signal :=
    any_pnps_gt1 | has_0D_sweep
]


# ------------------------------------------------------------
# 2. Summarize recurrence within each protein cluster
# ------------------------------------------------------------

cluster_recurrence <- cluster_gene_selection[
  ,
  .(
    # eligible background
    n_genes = uniqueN(gene_key),
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(Genus, na.rm = TRUE),
    n_colonies = uniqueN(colony_id),
    
    # 0D sweeps
    n_0D_genes = uniqueN(
      gene_key[has_0D_sweep]
    ),
    
    n_0D_MAGs = uniqueN(
      MAG_id[has_0D_sweep]
    ),
    
    n_0D_genera = uniqueN(
      Genus[
        has_0D_sweep &
          !is.na(Genus)
      ]
    ),
    
    n_0D_colonies = uniqueN(
      colony_id[has_0D_sweep]
    ),
    
    # pN/pS > 1
    n_pnps_genes = uniqueN(
      gene_key[any_pnps_gt1]
    ),
    
    n_pnps_MAGs = uniqueN(
      MAG_id[any_pnps_gt1]
    ),
    
    n_pnps_genera = uniqueN(
      Genus[
        any_pnps_gt1 &
          !is.na(Genus)
      ]
    ),
    
    n_pnps_colonies = uniqueN(
      colony_id[any_pnps_gt1]
    ),
    
    # either form of evidence
    n_selected_genes = uniqueN(
      gene_key[any_selection_signal]
    ),
    
    n_selected_MAGs = uniqueN(
      MAG_id[any_selection_signal]
    ),
    
    n_selected_genera = uniqueN(
      Genus[
        any_selection_signal &
          !is.na(Genus)
      ]
    ),
    
    n_selected_colonies = uniqueN(
      colony_id[any_selection_signal]
    ),
    
    # does the family contain both types of evidence?
    has_0D_signal = any(has_0D_sweep),
    
    has_pnps_signal = any(any_pnps_gt1),
    
    has_both_signal_types =
      any(has_0D_sweep) &
      any(any_pnps_gt1)
  ),
  by = cluster_id
]


cluster_recurrence[
  ,
  tested_primary :=
    n_genes >= 5 &
    n_MAGs >= 2
]

cluster_recurrence[
  ,
  .(
    n_total_clusters = .N,
    n_tested_primary = sum(tested_primary)
  )
]


recurrent_0D <- cluster_recurrence[
  n_0D_MAGs >= 2
]

recurrent_0D[
  order(
    -n_0D_genera,
    -n_0D_MAGs,
    -n_0D_genes
  )
][1:30]


recurrent_pnps <- cluster_recurrence[
  n_pnps_MAGs >= 2
]

recurrent_pnps[
  order(
    -n_pnps_genera,
    -n_pnps_MAGs,
    -n_pnps_genes
  )
][1:30]




both_signal_clusters <- cluster_recurrence[
  has_both_signal_types == TRUE
]

both_signal_clusters[
  order(
    -n_selected_genera,
    -n_selected_MAGs,
    -n_selected_genes
  )
][1:30]




candidate_parallel_clusters <- cluster_recurrence[
  n_0D_MAGs >= 2 &
    n_0D_genera >= 2
]

candidate_parallel_clusters[
  order(
    -n_0D_genera,
    -n_0D_MAGs,
    -n_0D_genes
  )
]






# How many protein families show recurrent 0D sweeps?

cluster_recurrence[
  ,
  .(
    any_0D = sum(n_0D_MAGs >= 1),
    recurrent_0D_2MAG = sum(n_0D_MAGs >= 2),
    recurrent_0D_2genus = sum(
      n_0D_MAGs >= 2 &
        n_0D_genera >= 2
    ),
    recurrent_0D_3genus = sum(
      n_0D_MAGs >= 3 &
        n_0D_genera >= 3
    )
  )
]


# How many show recurrent pN/pS > 1?

cluster_recurrence[
  ,
  .(
    any_pnps = sum(n_pnps_MAGs >= 1),
    recurrent_pnps_2MAG = sum(n_pnps_MAGs >= 2),
    recurrent_pnps_2genus = sum(
      n_pnps_MAGs >= 2 &
        n_pnps_genera >= 2
    ),
    recurrent_pnps_3genus = sum(
      n_pnps_MAGs >= 3 &
        n_pnps_genera >= 3
    )
  )
]


# Families with BOTH forms of cross-taxon evidence

cluster_recurrence[
  ,
  .(
    both_any = sum(
      has_0D_signal &
        has_pnps_signal
    ),
    
    both_recurrent_2genus = sum(
      n_0D_genera >= 2 &
        n_pnps_genera >= 2
    ),
    
    both_recurrent_3genus = sum(
      n_0D_genera >= 3 &
        n_pnps_genera >= 3
    )
  )
]





strong_candidates <- cluster_recurrence[
  n_0D_MAGs >= 2 &
    n_0D_genera >= 2
][
  order(
    -n_0D_genera,
    -n_0D_MAGs,
    -n_pnps_genera,
    -n_pnps_MAGs
  )
]

strong_candidates[
  1:30,
  .(
    cluster_id,
    
    n_genes,
    n_MAGs,
    n_genera,
    
    n_0D_genes,
    n_0D_MAGs,
    n_0D_genera,
    
    n_pnps_genes,
    n_pnps_MAGs,
    n_pnps_genera,
    
    has_both_signal_types
  )
]









# ============================================================
# 0D sweep permutation universe
# ============================================================



# gene_pop:
# one row per MAG x colony x eligible physical gene
#
# gene_flags:
# one row per eligible physical gene, contains
# effective_length_0D_sites

opp_dt <- unique(
  gene_flags[
    ,
    .(
      gene_key,
      effective_length_0D_sites
    )
  ]
)

perm_dt <- merge(
  gene_pop,
  opp_dt,
  by = "gene_key",
  all.x = TRUE
)

perm_dt[
  ,
  effective_length_0D_sites :=
    as.numeric(effective_length_0D_sites)
]

# eligible for a 0D event
perm_dt[
  ,
  callable_0D :=
    !is.na(effective_length_0D_sites) &
    is.finite(effective_length_0D_sites) &
    effective_length_0D_sites > 0
]

# sanity check
perm_dt[
  has_0D_sweep == TRUE,
  .(
    n_swept = .N,
    n_swept_not_callable =
      sum(!callable_0D)
  )
]




# ============================================================
# Observed cluster statistics
# ============================================================


obs_0D <- perm_dt[
  has_0D_sweep == TRUE,
  .(
    obs_0D_MAGs = uniqueN(MAG_id),
    obs_0D_genera = uniqueN(
      Genus[!is.na(Genus)]
    )
  ),
  by = cluster_id
]

obs <- cluster_recurrence[
  ,
  .(
    cluster_id,
    tested_primary
  )
]

obs <- merge(
  obs,
  obs_0D,
  by = "cluster_id",
  all.x = TRUE
)

obs[
  is.na(obs_0D_MAGs),
  obs_0D_MAGs := 0L
]

obs[
  is.na(obs_0D_genera),
  obs_0D_genera := 0L
]

# ============================================================
# Precompute populations containing observed sweeps
# ============================================================

# Only callable genes can receive a randomized 0D sweep
callable_dt <- perm_dt[
  callable_0D == TRUE
]

# number of swept genes to preserve per MAG x colony
sweep_burden <- perm_dt[
  ,
  .(
    k = sum(has_0D_sweep)
  ),
  by = .(
    MAG_id,
    colony_id
  )
][
  k > 0
]

# candidate row indices for each swept population
strata <- vector(
  "list",
  nrow(sweep_burden)
)

for (i in seq_len(nrow(sweep_burden))) {
  
  mag_i <- sweep_burden$MAG_id[i]
  col_i <- sweep_burden$colony_id[i]
  k_i   <- sweep_burden$k[i]
  
  idx <- which(
    callable_dt$MAG_id == mag_i &
      callable_dt$colony_id == col_i
  )
  
  if (length(idx) < k_i) {
    stop(
      paste(
        "Not enough callable genes for",
        mag_i,
        col_i
      )
    )
  }
  
  strata[[i]] <- list(
    idx = idx,
    k = k_i
  )
}





# ============================================================
# Pilot permutation
# ============================================================

set.seed(12345)

N_PERM <- 1000L

n_clusters <- length(all_clusters)

exceed_MAG <- integer(n_clusters)
exceed_genus <- integer(n_clusters)

sum_null_MAG <- numeric(n_clusters)
sum_null_genus <- numeric(n_clusters)

obs_MAG <- obs$obs_0D_MAGs
obs_genus <- obs$obs_0D_genera


for (b in seq_len(N_PERM)) {
  
  selected_idx <- unlist(
    lapply(
      strata,
      function(s) {
        
        w <- callable_dt[
          s$idx,
          effective_length_0D_sites
        ]
        
        sample(
          s$idx,
          size = s$k,
          replace = FALSE,
          prob = w
        )
      }
    ),
    use.names = FALSE
  )
  
  selected <- callable_dt[
    selected_idx,
    .(
      cluster_id,
      MAG_id,
      Genus
    )
  ]
  
  # unique MAG support per cluster
  mag_counts <- selected[
    ,
    .(
      null_MAGs = uniqueN(MAG_id)
    ),
    by = cluster_id
  ]
  
  # unique genus support per cluster
  genus_counts <- selected[
    ,
    .(
      null_genera = uniqueN(
        Genus[!is.na(Genus)]
      )
    ),
    by = cluster_id
  ]
  
  null_MAG <- integer(n_clusters)
  
  null_MAG[
    match(
      mag_counts$cluster_id,
      all_clusters
    )
  ] <- mag_counts$null_MAGs
  
  
  null_genus <- integer(n_clusters)
  
  null_genus[
    match(
      genus_counts$cluster_id,
      all_clusters
    )
  ] <- genus_counts$null_genera
  
  
  exceed_MAG <- exceed_MAG +
    as.integer(
      null_MAG >= obs_MAG
    )
  
  exceed_genus <- exceed_genus +
    as.integer(
      null_genus >= obs_genus
    )
  
  sum_null_MAG <- sum_null_MAG +
    null_MAG
  
  sum_null_genus <- sum_null_genus +
    null_genus
  
  
  if (b %% 100 == 0) {
    message(
      "Permutation ",
      b,
      " / ",
      N_PERM
    )
  }
}







# ============================================================
# Empirical statistics
# ============================================================

perm_results <- copy(obs)

perm_results[
  ,
  expected_0D_MAGs :=
    sum_null_MAG / N_PERM
]

perm_results[
  ,
  expected_0D_genera :=
    sum_null_genus / N_PERM
]

perm_results[
  ,
  p_empirical_MAG :=
    (exceed_MAG + 1) /
    (N_PERM + 1)
]

perm_results[
  ,
  p_empirical_genus :=
    (exceed_genus + 1) /
    (N_PERM + 1)
]

perm_results[
  ,
  FDR_MAG :=
    p.adjust(
      p_empirical_MAG,
      method = "BH"
    )
]

perm_results[
  ,
  FDR_genus :=
    p.adjust(
      p_empirical_genus,
      method = "BH"
    )
]




perm_results <- merge(
  perm_results,
  cluster_recurrence[
    ,
    .(
      cluster_id,
      n_genes,
      n_MAGs,
      n_genera,
      n_0D_genes,
      n_0D_MAGs,
      n_0D_genera,
      n_pnps_genes,
      n_pnps_MAGs,
      n_pnps_genera,
      has_both_signal_types
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)


perm_results[
  n_0D_MAGs >= 2 &
    n_0D_genera >= 2
][
  order(
    FDR_MAG,
    p_empirical_MAG,
    -n_0D_genera,
    -n_0D_MAGs
  )
][1:30]







#### QC


perm_results[
  ,
  p_empirical_MAG :=
    (exceed_MAG + 1) /
    (N_PERM + 1)
]

perm_results[
  ,
  p_empirical_genus :=
    (exceed_genus + 1) /
    (N_PERM + 1)
]



# ============================================================
# Effect-size / descriptive quantities
# ============================================================

perm_results[
  ,
  excess_0D_MAGs :=
    obs_0D_MAGs -
    expected_0D_MAGs
]

perm_results[
  ,
  fold_0D_MAGs :=
    (obs_0D_MAGs + 0.5) /
    (expected_0D_MAGs + 0.5)
]

perm_results[
  ,
  excess_0D_genera :=
    obs_0D_genera -
    expected_0D_genera
]



# ============================================================
# Multiple-testing correction
# ============================================================

perm_results[
  ,
  FDR_MAG_primary := NA_real_
]

perm_results[
  tested_primary == TRUE,
  FDR_MAG_primary :=
    p.adjust(
      p_empirical_MAG,
      method = "BH"
    )
]

perm_results[
  ,
  FDR_genus_primary := NA_real_
]

perm_results[
  tested_primary == TRUE,
  FDR_genus_primary :=
    p.adjust(
      p_empirical_genus,
      method = "BH"
    )
]




perm_results <- merge(
  perm_results,
  cluster_recurrence[
    ,
    .(
      cluster_id,
      n_genes,
      n_MAGs,
      n_genera,
      
      n_0D_genes,
      n_0D_MAGs,
      n_0D_genera,
      
      n_pnps_genes,
      n_pnps_MAGs,
      n_pnps_genera,
      
      has_both_signal_types,
      
      tested_primary
    )
  ],
  by = "cluster_id",
  all.x = TRUE,
  suffixes = c("", "_meta")
)





if ("tested_primary_meta" %in% names(perm_results)) {
  
  perm_results[
    ,
    tested_primary_meta := NULL
  ]
}


# ============================================================
# Final cross-genus recurrent candidates
# ============================================================

final_candidates <- perm_results[
  tested_primary == TRUE &
    n_0D_MAGs >= 2 &
    n_0D_genera >= 2
][
  order(
    FDR_MAG_primary,
    p_empirical_MAG,
    -n_0D_MAGs,
    -n_0D_genera
  )
]


final_candidates[
  1:30,
  .(
    cluster_id,
    
    n_genes,
    n_MAGs,
    n_genera,
    
    n_0D_genes,
    n_0D_MAGs,
    n_0D_genera,
    
    expected_0D_MAGs,
    excess_0D_MAGs,
    fold_0D_MAGs,
    
    p_empirical_MAG,
    FDR_MAG_primary,
    
    n_pnps_genes,
    n_pnps_MAGs,
    n_pnps_genera,
    
    has_both_signal_types
  )
]


final_candidates[
  ,
  .(
    n_cross_genus_candidates = .N,
    
    n_nominal_P05 =
      sum(
        p_empirical_MAG < 0.05
      ),
    
    n_FDR_05 =
      sum(
        FDR_MAG_primary < 0.05,
        na.rm = TRUE
      ),
    
    n_FDR_10 =
      sum(
        FDR_MAG_primary < 0.10,
        na.rm = TRUE
      )
  )
]





# ============================================================
# Final bookkeeping sanity check
# ============================================================

perm_results[
  ,
  .(
    n_all_clusters = .N,
    n_tested_primary = sum(tested_primary),
    
    n_tested_with_any_0D =
      sum(
        tested_primary &
          obs_0D_MAGs >= 1
      ),
    
    n_tested_recurrent_2MAG =
      sum(
        tested_primary &
          obs_0D_MAGs >= 2
      ),
    
    n_tested_cross_genus =
      sum(
        tested_primary &
          obs_0D_MAGs >= 2 &
          obs_0D_genera >= 2
      ),
    
    n_FDR05_all_tested =
      sum(
        tested_primary &
          FDR_MAG_primary < 0.05,
        na.rm = TRUE
      )
  )
]



significant_recurrent_families <- perm_results[
  tested_primary == TRUE &
    obs_0D_MAGs >= 2 &
    obs_0D_genera >= 2 &
    FDR_MAG_primary < 0.05
]

nrow(significant_recurrent_families)
# expected from what you just showed: 100


fwrite(
  perm_results,
  file.path(
    OUTDIR,
    "MMseqs_0D_sweep_permutation_all_clusters_10000.csv.gz"
  )
)

fwrite(
  significant_recurrent_families,
  file.path(
    OUTDIR,
    "MMseqs_significant_recurrent_0D_families_FDR05.csv"
  )
)






ANN_RDS <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/gene_level_annotations_all_colonies.rds"




# ============================================================
# 1. Load annotations
# ============================================================

ann_dt <- readRDS(ANN_RDS)

if (!is.data.table(ann_dt)) {
  ann_dt <- as.data.table(ann_dt)
}

ann_dt[
  ,
  corresponding_gene_call :=
    as.character(corresponding_gene_call)
]

names(ann_dt)




ann_dt[
  ,
  .(
    n_rows = .N,
    n_gene_calls = uniqueN(corresponding_gene_call)
  )
]
# What does the annotation table look like?
names(ann_dt)

ann_dt[
  ,
  .(
    n_rows = .N,
    n_gene_calls = uniqueN(corresponding_gene_call)
  )
]

# How many annotation records per corresponding_gene_call?
ann_dt[
  ,
  .N,
  by = corresponding_gene_call
][
  order(-N)
][1:20]




# Standardize type
ann_dt[
  ,
  corresponding_gene_call :=
    as.character(corresponding_gene_call)
]

gene_flags[
  ,
  corresponding_gene_call :=
    as.character(corresponding_gene_call)
]

# Annotation-key coverage
gene_flags[
  ,
  .(
    n_genes = .N,
    
    n_matching_annotation_key = sum(
      corresponding_gene_call %chin%
        ann_dt$corresponding_gene_call
    )
  )
]



safe_unlist <- function(x) {
  
  if (is.null(x)) {
    return(NA_character_)
  }
  
  y <- unlist(x)
  
  if (length(y) == 0) {
    return(NA_character_)
  }
  
  paste(
    unique(as.character(y)),
    collapse = ";;"
  )
}


if (
  !"GO_str" %in% names(ann_dt) &&
  "GO_slim_name_list" %in% names(ann_dt)
) {
  
  ann_dt[
    ,
    GO_str :=
      sapply(
        GO_slim_name_list,
        safe_unlist
      )
  ]
}


if (
  !"pfam_str" %in% names(ann_dt) &&
  "pfam_label_list" %in% names(ann_dt)
) {
  
  ann_dt[
    ,
    pfam_str :=
      sapply(
        pfam_label_list,
        safe_unlist
      )
  ]
}





ann_gene <- ann_dt[
  ,
  .(
    GO_annotation = paste(
      unique(
        GO_str[
          !is.na(GO_str) &
            GO_str != ""
        ]
      ),
      collapse = ";;"
    ),
    
    PFAM_annotation = paste(
      unique(
        pfam_str[
          !is.na(pfam_str) &
            pfam_str != ""
        ]
      ),
      collapse = ";;"
    )
  ),
  by = corresponding_gene_call
]

ann_gene[
  GO_annotation == "",
  GO_annotation := NA_character_
]

ann_gene[
  PFAM_annotation == "",
  PFAM_annotation := NA_character_
]



gene_annotated <- merge(
  gene_flags,
  ann_gene,
  by = "corresponding_gene_call",
  all.x = TRUE
)


nrow(gene_flags)
nrow(gene_annotated)

stopifnot(
  nrow(gene_annotated) ==
    nrow(gene_flags)
)

stopifnot(
  uniqueN(gene_annotated$gene_key) ==
    nrow(gene_annotated)
)



gene_annotated[
  ,
  .(
    n_genes = .N,
    
    n_PFAM =
      sum(!is.na(PFAM_annotation)),
    
    n_GO =
      sum(!is.na(GO_annotation)),
    
    pct_PFAM =
      mean(!is.na(PFAM_annotation)) * 100,
    
    pct_GO =
      mean(!is.na(GO_annotation)) * 100
  )
]






sig_cluster_ids <- significant_recurrent_families$cluster_id

sig_genes <- gene_annotated[
  cluster_id %chin% sig_cluster_ids
]
cluster_interp <- sig_genes[
  ,
  .(
    n_PFAM_annotated = sum(!is.na(PFAM_annotation)),
    n_GO_annotated   = sum(!is.na(GO_annotation)),
    
    n_unique_PFAM = uniqueN(
      PFAM_annotation[!is.na(PFAM_annotation)]
    ),
    
    n_unique_GO = uniqueN(
      GO_annotation[!is.na(GO_annotation)]
    ),
    
    PFAM_labels = paste(
      unique(
        PFAM_annotation[!is.na(PFAM_annotation)]
      ),
      collapse = " || "
    ),
    
    GO_labels = paste(
      unique(
        GO_annotation[!is.na(GO_annotation)]
      ),
      collapse = " || "
    )
  ),
  by = cluster_id
]

cluster_interp <- merge(
  significant_recurrent_families,
  cluster_interp,
  by = "cluster_id",
  all.x = TRUE
)
cluster_interp[
  order(
    -n_0D_genera,
    -n_0D_MAGs,
    -n_pnps_genera
  )
][1:30]


cluster_interp[
  ,
  support_class := fifelse(
    n_pnps_genera >= 2,
    "Recurrent 0D + recurrent pN/pS",
    fifelse(
      n_pnps_genes >= 1,
      "Recurrent 0D + pN/pS support",
      "Recurrent 0D only"
    )
  )
]

table(cluster_interp$support_class)






# FIGURE

library(ggplot2)

plot_dt <- perm_results[
  tested_primary == TRUE
]

plot_dt[
  ,
  significant_recurrent :=
    FDR_MAG_primary < 0.05 &
    obs_0D_MAGs >= 2 &
    obs_0D_genera >= 2
]

ggplot(
  plot_dt,
  aes(
    x = expected_0D_MAGs,
    y = obs_0D_MAGs
  )
) +
  geom_point(
    aes(
      shape = significant_recurrent
    ),
    alpha = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  labs(
    x = "Expected MAGs with 0D sweep",
    y = "Observed MAGs with 0D sweep"
  ) +
  theme_classic()













# ============================================================
# Build long-form 0D sweep table for significant families
# ============================================================

season_map <- c(
  "May"       = "Spring",
  "June"      = "Spring",
  "July"      = "Summer",
  "August"    = "Summer",
  "September" = "Summer",
  "October"   = "Fall",
  "November"  = "Fall",
  "January"   = "Winter",
  "February"  = "Winter"
)

sig_ids <- significant_recurrent_families$cluster_id

# Start from actual 0D sweep rows
sig_sweeps <- zeroD_dt[
  gene_key %chin% gene_flags[
    cluster_id %chin% sig_ids,
    gene_key
  ]
]

# Attach cluster/genus/pN-pS/annotation metadata
sig_meta <- unique(
  gene_annotated[
    cluster_id %chin% sig_ids,
    .(
      gene_key,
      cluster_id,
      MAG_id,
      Genus,
      any_pnps_gt1,
      PFAM_annotation,
      GO_annotation
    )
  ]
)

sig_sweeps <- merge(
  sig_sweeps,
  sig_meta,
  by = c("gene_key", "MAG_id"),
  all.x = TRUE
)

# Add seasons
sig_sweeps[, Season_t1 := season_map[as.character(Month_t1)]]
sig_sweeps[, Season_t2 := season_map[as.character(Month_t2)]]

sig_sweeps[, season_transition :=
             paste(Season_t1, Season_t2, sep = " -> ")
]






sig_sweeps[
  ,
  .N,
  by = season_transition
][
  order(-N)
]



sig_sweep_events <- unique(
  sig_sweeps[
    ,
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2,
      Season_t1,
      Season_t2,
      season_transition
    )
  ]
)

sig_sweep_events[
  ,
  .(
    n_events = .N,
    n_clusters = uniqueN(cluster_id),
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(Genus, na.rm = TRUE)
  ),
  by = season_transition
][
  order(-n_events)
]





family_season_concordance <- sig_sweep_events[
  ,
  {
    tab <- table(season_transition)
    
    top_transition <- names(tab)[which.max(tab)]
    top_n <- max(tab)
    
    .(
      n_events = .N,
      n_MAGs = uniqueN(MAG_id),
      n_genera = uniqueN(Genus, na.rm = TRUE),
      n_transitions = uniqueN(season_transition),
      dominant_transition = top_transition,
      dominant_transition_n = top_n,
      seasonal_concordance = top_n / .N
    )
  },
  by = cluster_id
]

family_season_concordance[
  order(
    -seasonal_concordance,
    -n_MAGs
  )
][1:30]







all_zeroD_events <- unique(
  zeroD_dt[
    ,
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

all_zeroD_events[, Season_t1 := season_map[as.character(Month_t1)]]
all_zeroD_events[, Season_t2 := season_map[as.character(Month_t2)]]

all_zeroD_events[, season_transition :=
                   paste(Season_t1, Season_t2, sep = " -> ")
]

all_zeroD_events[, recurrent_family :=
                   gene_key %chin% sig_meta$gene_key
]

season_test <- all_zeroD_events[
  !is.na(season_transition),
  .N,
  by = .(
    recurrent_family,
    season_transition
  )
]

season_test





null_recurrent_family_count <- integer(N_PERM)

for (b in seq_len(N_PERM)) {
  

  selected <- callable_dt[
    selected_idx,
    .(
      cluster_id,
      MAG_id,
      Genus
    )
  ]
  
  null_cluster_stats <- selected[
    ,
    .(
      n_MAGs = uniqueN(MAG_id),
      n_genera = uniqueN(
        Genus[!is.na(Genus)]
      )
    ),
    by = cluster_id
  ]
  
  null_recurrent_family_count[b] <- null_cluster_stats[
    n_MAGs >= 2 &
      n_genera >= 2,
    .N
  ]
}


observed_recurrent_family_count <- cluster_recurrence[
  tested_primary == TRUE &
    n_0D_MAGs >= 2 &
    n_0D_genera >= 2,
  .N
]

observed_recurrent_family_count



null_dt <- data.table(
  n_recurrent_families =
    null_recurrent_family_count
)

ggplot(
  null_dt,
  aes(x = n_recurrent_families)
) +
  geom_histogram(
    bins = 30,
    color = "black",
    fill = "grey80"
  ) +
  geom_vline(
    xintercept = observed_recurrent_family_count,
    linewidth = 1.2
  ) +
  labs(
    x = "Protein families with recurrent 0D sweeps\nacross >=2 MAGs and >=2 genera",
    y = "Permutations"
  ) +
  theme_classic()











# ============================================================
# Genome-wide null:
# number of recurrent cross-genus protein families
# ============================================================

set.seed(12345)

N_PERM <- 1000

null_recurrent_family_count <- rep(NA_integer_, N_PERM)

testable_ids <- cluster_recurrence[
  tested_primary == TRUE,
  cluster_id
]

for (b in seq_len(N_PERM)) {
  
  selected_idx <- unlist(
    lapply(
      strata,
      function(s) {
        
        w <- callable_dt[
          s$idx,
          effective_length_0D_sites
        ]
        
        sample(
          s$idx,
          size = s$k,
          replace = FALSE,
          prob = w
        )
      }
    ),
    use.names = FALSE
  )
  
  selected <- callable_dt[
    selected_idx,
    .(
      cluster_id,
      MAG_id,
      Genus
    )
  ]
  
  null_cluster_stats <- selected[
    cluster_id %chin% testable_ids,
    .(
      n_MAGs = uniqueN(MAG_id),
      n_genera = uniqueN(
        Genus[!is.na(Genus)]
      )
    ),
    by = cluster_id
  ]
  
  null_recurrent_family_count[b] <- null_cluster_stats[
    n_MAGs >= 2 &
      n_genera >= 2,
    .N
  ]
  
  if (b %% 1000 == 0) {
    message(b, " / ", N_PERM)
  }
}





summary(null_recurrent_family_count)

sum(is.na(null_recurrent_family_count))

table(null_recurrent_family_count)[1:20]



observed_recurrent_family_count <- cluster_recurrence[
  tested_primary == TRUE &
    n_0D_MAGs >= 2 &
    n_0D_genera >= 2,
  .N
]

observed_recurrent_family_count

p_global <- (
  sum(
    null_recurrent_family_count >=
      observed_recurrent_family_count
  ) + 1
) / (N_PERM + 1)

p_global


null_dt <- data.table(
  n_recurrent_families =
    null_recurrent_family_count
)

ggplot(
  null_dt,
  aes(x = n_recurrent_families)
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    xintercept =
      observed_recurrent_family_count,
    linewidth = 1.1
  ) +
  annotate(
    "text",
    x = observed_recurrent_family_count,
    y = Inf,
    label = paste0(
      "Observed = ",
      observed_recurrent_family_count
    ),
    angle = 90,
    vjust = 1.5
  ) +
  labs(
    x = paste0(
      "Protein families with recurrent 0D sweeps\n",
      "across >=2 MAGs and >=2 genera"
    ),
    y = "Permutations"
  ) +
  theme_classic()






season_mat <- dcast(
  season_test,
  recurrent_family ~ season_transition,
  value.var = "N",
  fill = 0
)

season_mat



season_counts <- as.matrix(
  season_mat[
    ,
    -"recurrent_family"
  ]
)

rownames(season_counts) <-
  season_mat$recurrent_family

chisq.test(season_counts)


season_chi <- chisq.test(season_counts)

round(
  season_chi$stdres,
  2
)








all_zeroD_events <- unique(
  zeroD_dt[
    ,
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

all_zeroD_events[
  ,
  Season_t1 :=
    season_map[as.character(Month_t1)]
]

all_zeroD_events[
  ,
  Season_t2 :=
    season_map[as.character(Month_t2)]
]

all_zeroD_events[
  ,
  season_transition :=
    paste(
      Season_t1,
      Season_t2,
      sep = " -> "
    )
]

all_zeroD_events <- all_zeroD_events[
  !is.na(Season_t1) &
    !is.na(Season_t2)
]

all_zeroD_events[
  ,
  recurrent_family :=
    gene_key %chin% sig_meta$gene_key
]





transition_levels <- sort(
  unique(
    all_zeroD_events$season_transition
  )
)

obs_season <- all_zeroD_events[
  recurrent_family == TRUE,
  .N,
  by = season_transition
]

obs_vec <- setNames(
  integer(length(transition_levels)),
  transition_levels
)

obs_vec[
  obs_season$season_transition
] <- obs_season$N





recurrent_burden <- all_zeroD_events[
  ,
  .(
    k = sum(recurrent_family)
  ),
  by = .(
    MAG_id,
    colony_id
  )
][
  k > 0
]




set.seed(12345)

N_SEASON_PERM <- 1000L

null_season <- matrix(
  0L,
  nrow = N_SEASON_PERM,
  ncol = length(transition_levels)
)

colnames(null_season) <- transition_levels


for (b in seq_len(N_SEASON_PERM)) {
  
  sampled_transitions <- character()
  
  for (i in seq_len(nrow(recurrent_burden))) {
    
    mag_i <- recurrent_burden$MAG_id[i]
    col_i <- recurrent_burden$colony_id[i]
    k_i   <- recurrent_burden$k[i]
    
    pool <- callable_pools[[timepair_burden$stratum_id[i]]]
    
    sampled <- pool[
      sample(
        .N,
        size = k_i,
        replace = FALSE
      )
    ]
    
    sampled_transitions <- c(
      sampled_transitions,
      sampled$season_transition
    )
  }
  
  tab <- table(
    factor(
      sampled_transitions,
      levels = transition_levels
    )
  )
  
  null_season[b, ] <- as.integer(tab)
}



season_results <- data.table(
  season_transition = transition_levels,
  observed = as.integer(
    obs_vec[transition_levels]
  ),
  expected = colMeans(
    null_season
  )
)

season_results[
  ,
  enrichment :=
    observed / expected
]

season_results[
  ,
  p_upper :=
    vapply(
      seq_along(transition_levels),
      function(j) {
        
        (
          sum(
            null_season[, j] >=
              observed[j]
          ) + 1
        ) /
          (N_SEASON_PERM + 1)
        
      },
      numeric(1)
    )
]

season_results[
  ,
  FDR :=
    p.adjust(
      p_upper,
      method = "BH"
    )
]

season_results[
  order(-enrichment)
]






family_season_concordance[
  n_MAGs >= 3
][
  order(
    -seasonal_concordance,
    -n_genera,
    -n_MAGs
  )
]



family_season_concordance[
  n_MAGs >= 3 &
    seasonal_concordance >= 0.66
]










# ============================================================
# Cross-genus recurrence score
# ============================================================

cross_genus_pair_score <- function(dt) {
  
  x <- unique(
    dt[
      !is.na(Genus),
      .(
        cluster_id,
        MAG_id,
        Genus
      )
    ]
  )
  
  # total MAG pairs in each protein family
  family_pairs <- x[
    ,
    .(
      n_MAG = uniqueN(MAG_id),
      total_pairs =
        choose(uniqueN(MAG_id), 2)
    ),
    by = cluster_id
  ]
  
  # subtract pairs in which both MAGs belong to same genus
  within_genus <- x[
    ,
    .(
      n_MAG_genus = uniqueN(MAG_id)
    ),
    by = .(
      cluster_id,
      Genus
    )
  ][
    ,
    .(
      same_genus_pairs =
        sum(
          choose(n_MAG_genus, 2)
        )
    ),
    by = cluster_id
  ]
  
  z <- merge(
    family_pairs,
    within_genus,
    by = "cluster_id",
    all.x = TRUE
  )
  
  z[
    is.na(same_genus_pairs),
    same_genus_pairs := 0
  ]
  
  z[
    ,
    cross_genus_pairs :=
      total_pairs -
      same_genus_pairs
  ]
  
  sum(z$cross_genus_pairs)
}




obs_pair_score <- cross_genus_pair_score(
  perm_dt[
    has_0D_sweep == TRUE &
      cluster_id %chin% testable_ids
  ]
)

obs_pair_score




set.seed(12345)

N_PERM <- 1000

null_pair_score <- numeric(N_PERM)

for (b in seq_len(N_PERM)) {
  
  selected_idx <- unlist(
    lapply(
      strata,
      function(s) {
        
        w <- callable_dt[
          s$idx,
          effective_length_0D_sites
        ]
        
        sample(
          s$idx,
          size = s$k,
          replace = FALSE,
          prob = w
        )
      }
    ),
    use.names = FALSE
  )
  
  selected <- callable_dt[
    selected_idx,
    .(
      cluster_id,
      MAG_id,
      Genus
    )
  ][
    cluster_id %chin% testable_ids
  ]
  
  null_pair_score[b] <-
    cross_genus_pair_score(selected)
  
  if (b %% 1000 == 0) {
    message(b, " / ", N_PERM)
  }
}

summary(null_pair_score)

p_pair <- (
  sum(
    null_pair_score >=
      obs_pair_score
  ) + 1
) / (N_PERM + 1)

p_pair




pair_null_dt <- data.table(
  cross_genus_pairs =
    null_pair_score
)

ggplot(
  pair_null_dt,
  aes(x = cross_genus_pairs)
) +
  geom_histogram(
    bins = 40,
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    xintercept = obs_pair_score,
    linewidth = 1.2
  ) +
  labs(
    x = "Cross-genus pairs of MAGs with 0D sweeps\nin the same protein family",
    y = "Permutations"
  ) +
  theme_classic()



temporal_test_families <- family_season_concordance[
  n_MAGs >= 3,
  cluster_id
]


obs_temporal_score <- family_season_concordance[
  cluster_id %chin% temporal_test_families,
  sum(
    dominant_transition_n - 1
  )
]




















# Start from significant recurrent families only
candidate_families <- copy(
  significant_recurrent_families
)

# Build sweep-event table for these families
strict_events <- unique(
  sig_sweep_events[
    ,
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2,
      season_transition
    )
  ]
)



strict_family_summary <- strict_events[
  ,
  {
    # distinct swept lineages
    n_mag <- uniqueN(MAG_id)
    n_genus <- uniqueN(Genus[!is.na(Genus)])
    n_colony <- uniqueN(colony_id)
    
    # support by exact month/timepoint pair
    tp <- paste(
      Month_t1,
      Month_t2,
      sep = "->"
    )
    
    tp_tab <- table(tp)
    
    max_same_transition <- if (length(tp_tab)) {
      max(tp_tab)
    } else {
      0L
    }
    
    # For the dominant transition, how broad is taxonomic/colony support?
    dom_transition <- if (length(tp_tab)) {
      names(tp_tab)[which.max(tp_tab)]
    } else {
      NA_character_
    }
    
    dom_rows <- tp == dom_transition
    
    .(
      n_0D_MAGs = n_mag,
      n_0D_genera = n_genus,
      n_0D_colonies = n_colony,
      
      n_time_transitions =
        uniqueN(tp),
      
      dominant_transition =
        dom_transition,
      
      n_sweeps_same_transition =
        max_same_transition,
      
      same_transition_MAGs =
        uniqueN(MAG_id[dom_rows]),
      
      same_transition_genera =
        uniqueN(
          Genus[
            dom_rows &
              !is.na(Genus)
          ]
        ),
      
      same_transition_colonies =
        uniqueN(colony_id[dom_rows])
    )
  },
  by = cluster_id
]



strict_family_summary <- merge(
  strict_family_summary,
  cluster_recurrence[
    ,
    .(
      cluster_id,
      n_pnps_genes,
      n_pnps_MAGs,
      n_pnps_genera
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)



strict_family_summary[
  ,
  evidence_tier := fifelse(
    same_transition_genera >= 2 &
      same_transition_colonies >= 2 &
      n_pnps_genera >= 2,
    
    "Tier 5: concurrent cross-genus sweeps + recurrent pN/pS",
    
    fifelse(
      same_transition_genera >= 2 &
        same_transition_colonies >= 2,
      
      "Tier 4: concurrent cross-genus sweeps",
      
      fifelse(
        n_0D_genera >= 2 &
          n_0D_colonies >= 2,
        
        "Tier 3: cross-genus, cross-colony sweeps",
        
        "Lower support"
      )
    )
  )
]



table(
  strict_family_summary$evidence_tier
)

strict_family_summary[
  order(
    -same_transition_genera,
    -same_transition_colonies,
    -n_pnps_genera
  )
]







# ============================================================
# Rebuild zeroD_clustered cleanly
# ============================================================

# Check current columns
names(zeroD_dt)
names(gene_flags)

# Metadata from our authoritative physical-gene table
gene_meta <- unique(
  gene_flags[
    ,
    .(
      gene_key,
      cluster_id,
      Genus
    )
  ]
)

stopifnot(
  uniqueN(gene_meta$gene_key) ==
    nrow(gene_meta)
)

# Remove any old cluster/genus columns before merging
zeroD_clustered <- copy(zeroD_dt)

drop_cols <- intersect(
  names(zeroD_clustered),
  c(
    "cluster_id",
    "Genus",
    "Genus.x",
    "Genus.y"
  )
)

if (length(drop_cols) > 0) {
  zeroD_clustered[
    ,
    (drop_cols) := NULL
  ]
}

# Attach cluster + genus by physical gene
zeroD_clustered <- merge(
  zeroD_clustered,
  gene_meta,
  by = "gene_key",
  all.x = TRUE
)

# Sanity checks
names(zeroD_clustered)

zeroD_clustered[
  ,
  .(
    n_rows = .N,
    n_with_cluster = sum(!is.na(cluster_id)),
    n_with_genus = sum(!is.na(Genus))
  )
]




obs_concurrent_pair_score <-
  concurrent_cross_genus_pair_score(
    zeroD_clustered[
      cluster_id %chin% testable_ids
    ]
  )

obs_concurrent_pair_score

zeroD_clustered[
  1:5,
  .(
    gene_key,
    cluster_id,
    MAG_id,
    Genus,
    Month_t1,
    Month_t2
  )
]




# ============================================================
# Concurrent cross-genus pair score
# ============================================================

concurrent_cross_genus_pair_score <- function(dt) {
  
  # dt needs:
  # cluster_id
  # MAG_id
  # Genus
  # Month_t1
  # Month_t2
  
  x <- unique(
    dt[
      !is.na(Genus) &
        !is.na(Month_t1) &
        !is.na(Month_t2),
      .(
        cluster_id,
        MAG_id,
        Genus,
        Month_t1,
        Month_t2
      )
    ]
  )
  
  x[
    ,
    timepair := paste(
      Month_t1,
      Month_t2,
      sep = "->"
    )
  ]
  
  # total MAG pairs within each family x timepair
  total_pairs <- x[
    ,
    .(
      total_pairs =
        choose(uniqueN(MAG_id), 2)
    ),
    by = .(
      cluster_id,
      timepair
    )
  ]
  
  # subtract same-genus MAG pairs
  within_genus <- x[
    ,
    .(
      n_MAG_genus =
        uniqueN(MAG_id)
    ),
    by = .(
      cluster_id,
      timepair,
      Genus
    )
  ][
    ,
    .(
      same_genus_pairs =
        sum(
          choose(n_MAG_genus, 2)
        )
    ),
    by = .(
      cluster_id,
      timepair
    )
  ]
  
  z <- merge(
    total_pairs,
    within_genus,
    by = c(
      "cluster_id",
      "timepair"
    ),
    all.x = TRUE
  )
  
  z[
    is.na(same_genus_pairs),
    same_genus_pairs := 0
  ]
  
  z[
    ,
    cross_genus_pairs :=
      total_pairs -
      same_genus_pairs
  ]
  
  sum(z$cross_genus_pairs)
}



obs_concurrent_pair_score <-
  concurrent_cross_genus_pair_score(
    zeroD_clustered[
      gene_key %chin% gene_flags[
        cluster_id %chin% testable_ids,
        gene_key
      ],
      .(
        cluster_id = gene_flags[
          match(
            gene_key,
            gene_flags$gene_key
          ),
          cluster_id
        ],
        MAG_id,
        Genus,
        Month_t1,
        Month_t2
      )
    ]
  )




# ============================================================
# Observed 0D sweep timepairs by population
# ============================================================

observed_sweep_events <- unique(
  zeroD_clustered[
    ,
    .(
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      gene_key
    )
  ]
)

time_strata <- observed_sweep_events[
  ,
  .(
    Month_t1 = list(Month_t1),
    Month_t2 = list(Month_t2)
  ),
  by = .(
    MAG_id,
    colony_id
  )
]



set.seed(12345)

N_PERM <- 1000L

null_concurrent_pair_score <- numeric(N_PERM)

for (b in seq_len(N_PERM)) {
  
  perm_events <- vector(
    "list",
    nrow(time_strata)
  )
  
  for (i in seq_len(nrow(time_strata))) {
    
    mag_i <- time_strata$MAG_id[i]
    col_i <- time_strata$colony_id[i]
    
    t1_i <- time_strata$Month_t1[[i]]
    t2_i <- time_strata$Month_t2[[i]]
    
    k_i <- length(t1_i)
    
    pool <- callable_dt[
      MAG_id == mag_i &
        colony_id == col_i
    ]
    
    if (nrow(pool) < k_i) {
      stop(
        paste(
          "Not enough callable genes:",
          mag_i,
          col_i
        )
      )
    }
    
    sampled_idx <- sample(
      seq_len(nrow(pool)),
      size = k_i,
      replace = FALSE,
      prob = pool$effective_length_0D_sites
    )
    
    sampled <- pool[
      sampled_idx,
      .(
        cluster_id,
        MAG_id,
        Genus
      )
    ]
    
    # preserve actual observed timepairs
    sampled[
      ,
      Month_t1 := t1_i
    ]
    
    sampled[
      ,
      Month_t2 := t2_i
    ]
    
    perm_events[[i]] <- sampled
  }
  
  perm_dt_b <- rbindlist(
    perm_events,
    use.names = TRUE
  )
  
  perm_dt_b <- perm_dt_b[
    cluster_id %chin% testable_ids
  ]
  
  null_concurrent_pair_score[b] <-
    concurrent_cross_genus_pair_score(
      perm_dt_b
    )
  
  if (b %% 1000 == 0) {
    message(
      b,
      " / ",
      N_PERM
    )
  }
}




summary(
  null_concurrent_pair_score
)

p_concurrent <- (
  sum(
    null_concurrent_pair_score >=
      obs_concurrent_pair_score
  ) + 1
) /
  (N_PERM + 1)

p_concurrent





null_concurrent_dt <- data.table(
  concurrent_cross_genus_pairs =
    null_concurrent_pair_score
)

ggplot(
  null_concurrent_dt,
  aes(
    x = concurrent_cross_genus_pairs
  )
) +
  geom_histogram(
    bins = 40,
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    xintercept =
      obs_concurrent_pair_score,
    linewidth = 1.2
  ) +
  labs(
    x = paste0(
      "Cross-genus MAG pairs with 0D sweeps\\n",
      "in the same protein family and time interval"
    ),
    y = "Permutations"
  ) +
  theme_classic()+
  scale_x_continuous(trans="log10")





# How many unique observed sweep units are entering the score?
obs_score_dt <- unique(
  zeroD_clustered[
    cluster_id %chin% testable_ids &
      !is.na(Genus) &
      !is.na(Month_t1) &
      !is.na(Month_t2),
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2
    )
  ]
)

obs_score_dt[
  ,
  .(
    n_rows = .N,
    n_genes = uniqueN(gene_key),
    n_MAGs = uniqueN(MAG_id),
    n_clusters = uniqueN(cluster_id),
    n_timepairs = uniqueN(
      paste(Month_t1, Month_t2)
    )
  )
]



obs_score_dt[
  ,
  .N,
  by = .(
    gene_key,
    MAG_id,
    colony_id
  )
][
  order(-N)
][1:30]



table(
  obs_score_dt[
    ,
    .N,
    by = .(
      gene_key,
      MAG_id,
      colony_id
    )
  ]$N
)



sum(
  lengths(time_strata$Month_t1)
)

nrow(obs_score_dt)



obs_concurrent_clusters <- obs_score_dt[
  ,
  .(
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(
      Genus[!is.na(Genus)]
    )
  ),
  by = .(
    cluster_id,
    Month_t1,
    Month_t2
  )
][
  n_MAGs >= 2 &
    n_genera >= 2
]

obs_n_concurrent_families <-
  uniqueN(obs_concurrent_clusters$cluster_id)

obs_n_concurrent_families









#### MAYB EFINAL NULLLL #######


# ============================================================
# Observed concurrent cross-taxon families
# SAME protein family
# SAME exact month pair
# >=2 MAGs
# >=2 genera
# >=2 colonies
# ============================================================

obs_events <- unique(
  zeroD_clustered[
    cluster_id %chin% testable_ids &
      !is.na(Genus) &
      !is.na(Month_t1) &
      !is.na(Month_t2),
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2
    )
  ]
)

obs_concurrent_by_time <- obs_events[
  ,
  .(
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(Genus),
    n_colonies = uniqueN(colony_id)
  ),
  by = .(
    cluster_id,
    Month_t1,
    Month_t2
  )
][
  n_MAGs >= 2 &
    n_genera >= 2 &
    n_colonies >= 2
]

# Number of unique protein families meeting the stringent criterion
obs_n_concurrent_families <-
  uniqueN(obs_concurrent_by_time$cluster_id)

obs_n_concurrent_families

# Number of qualifying family x timepair combinations
nrow(obs_concurrent_by_time)




# ============================================================
# Observed sweep burden by population x exact timepair
# ============================================================

all_obs_events <- unique(
  zeroD_clustered[
    !is.na(Month_t1) &
      !is.na(Month_t2),
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

timepair_burden <- all_obs_events[
  ,
  .(
    k = uniqueN(gene_key)
  ),
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2
  )
]

summary(timepair_burden$k)
sum(timepair_burden$k)

# This should equal:
nrow(all_obs_events)




set.seed(12345)

N_PERM <- 10000L

null_concurrent_family_count <-
  integer(N_PERM)

null_concurrent_timepair_count <-
  integer(N_PERM)


for (b in seq_len(N_PERM)) {
  
  perm_list <- vector(
    "list",
    nrow(timepair_burden)
  )
  
  for (i in seq_len(nrow(timepair_burden))) {
    
    mag_i <- timepair_burden$MAG_id[i]
    col_i <- timepair_burden$colony_id[i]
    
    t1_i <- timepair_burden$Month_t1[i]
    t2_i <- timepair_burden$Month_t2[i]
    
    k_i <- timepair_burden$k[i]
    
    # All genes that could have received a 0D sweep
    # in this MAG x colony population
    pool <- callable_dt[
      MAG_id == mag_i &
        colony_id == col_i
    ]
    
    if (nrow(pool) < k_i) {
      stop(
        paste(
          "Insufficient callable genes:",
          mag_i,
          col_i,
          t1_i,
          t2_i
        )
      )
    }
    
    picked <- pool[
      sample(
        seq_len(.N),
        size = k_i,
        replace = FALSE,
        prob = effective_length_0D_sites
      ),
      .(
        cluster_id,
        MAG_id,
        colony_id,
        Genus
      )
    ]
    
    picked[, Month_t1 := t1_i]
    picked[, Month_t2 := t2_i]
    
    perm_list[[i]] <- picked
  }
  
  perm_events <- rbindlist(
    perm_list,
    use.names = TRUE
  )
  
  # Same primary family universe
  perm_events <- perm_events[
    cluster_id %chin% testable_ids &
      !is.na(Genus)
  ]
  
  # Score each protein family within each exact timepair
  perm_concurrent <- perm_events[
    ,
    .(
      n_MAGs = uniqueN(MAG_id),
      n_genera = uniqueN(Genus),
      n_colonies = uniqueN(colony_id)
    ),
    by = .(
      cluster_id,
      Month_t1,
      Month_t2
    )
  ][
    n_MAGs >= 2 &
      n_genera >= 2 &
      n_colonies >= 2
  ]
  
  null_concurrent_family_count[b] <-
    uniqueN(
      perm_concurrent$cluster_id
    )
  
  null_concurrent_timepair_count[b] <-
    nrow(perm_concurrent)
  
  
  if (b %% 1000 == 0) {
    message(b, " / ", N_PERM)
  }
}





summary(
  null_concurrent_family_count
)

obs_n_concurrent_families

p_concurrent_family <- (
  sum(
    null_concurrent_family_count >=
      obs_n_concurrent_families
  ) + 1
) /
  (N_PERM + 1)

p_concurrent_family



quantile(
  null_concurrent_family_count,
  c(
    0.025,
    0.5,
    0.975,
    0.99,
    0.999
  )
)







# ============================================================
# Quick concurrent-family permutation figure
# ============================================================

null_plot_dt <- data.table(
  n_concurrent_families = null_concurrent_family_count
)

# Empirical P from current run
p_plot <- (
  sum(
    null_concurrent_family_count >=
      obs_n_concurrent_families
  ) + 1
) /
  (length(null_concurrent_family_count) + 1)

# Null summaries
null_median <- median(
  null_concurrent_family_count
)

null_ci <- quantile(
  null_concurrent_family_count,
  c(0.025, 0.975)
)

p_concurrent_null <- ggplot(
  null_plot_dt,
  aes(x = n_concurrent_families)
) +
  
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "grey80",
    color = "black",
    linewidth = 0.3
  ) +
  
  # 95% null interval
  geom_vline(
    xintercept = null_ci,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # Observed value
  geom_vline(
    xintercept = obs_n_concurrent_families,
    linewidth = 1.1
  ) +
  
  annotate(
    "text",
    x = obs_n_concurrent_families + 0.5,
    y = Inf,
    label = paste0(
      "Observed = ",
      obs_n_concurrent_families
    ),
    angle = 90,
    hjust = 1.1,
    vjust = -0.3,
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = min(null_plot_dt$n_concurrent_families),
    y = Inf,
    label = paste0(
      "Null median = ",
      null_median,
      "\n95% interval = ",
      null_ci[1],
      "\u2013",
      null_ci[2],
      "\nEmpirical P = ",
      signif(p_plot, 2)
    ),
    hjust = 0,
    vjust = 1.3,
    size = 3.4
  ) +
  
  labs(
    x = paste0(
      "Protein families with concurrent nonsynonymous sweeps\n",
      "across \u22652 genera and \u22652 colonies"
    ),
    y = "Permutations"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title.x = element_text(
      margin = margin(t = 8)
    )
  )

p_concurrent_null
















matrix_candidates <- strict_family_summary[
  same_transition_genera >= 2 &
    same_transition_colonies >= 2
]

matrix_candidates <- merge(
  matrix_candidates,
  cluster_interp[
    ,
    .(
      cluster_id,
      PFAM_labels,
      GO_labels,
      n_genes,
      n_MAGs,
      n_genera
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_candidates[
  ,
  matrix_score :=
    3 * same_transition_genera +
    2 * same_transition_MAGs +
    same_transition_colonies +
    2 * n_pnps_genera
]

matrix_candidates[
  order(
    -matrix_score,
    -same_transition_genera,
    -same_transition_MAGs
  ),
  .(
    cluster_id,
    n_genes,
    n_genera,
    same_transition_MAGs,
    same_transition_genera,
    same_transition_colonies,
    dominant_transition,
    n_pnps_MAGs,
    n_pnps_genera,
    PFAM_labels,
    GO_labels
  )
][1:20]



top20_ids <- matrix_candidates[
  order(-matrix_score),
  head(cluster_id, 20)
]

sig_sweep_events[
  cluster_id %chin% top20_ids,
  .(
    n_families = uniqueN(cluster_id),
    n_swept_MAGs = uniqueN(MAG_id)
  ),
  by = Genus
][order(-n_families)]





matrix_long <- sig_sweep_events[
  cluster_id %chin% matrix_candidates$cluster_id,
  .(
    n_swept_MAGs = uniqueN(MAG_id),
    n_swept_colonies = uniqueN(colony_id),
    
    dominant_timepair = {
      tt <- paste(Month_t1, Month_t2, sep = " -> ")
      names(sort(table(tt), decreasing = TRUE))[1]
    }
  ),
  by = .(
    cluster_id,
    Genus
  )
]





family_presence <- gene_annotated[
  cluster_id %chin% matrix_candidates$cluster_id,
  .(
    n_homologs = uniqueN(gene_key),
    any_pnps = any(any_pnps_gt1)
  ),
  by = .(
    cluster_id,
    Genus
  )
]

matrix_long <- merge(
  family_presence,
  matrix_long,
  by = c("cluster_id", "Genus"),
  all.x = TRUE
)

matrix_long[
  is.na(n_swept_MAGs),
  n_swept_MAGs := 0L
]

matrix_long[
  is.na(n_swept_colonies),
  n_swept_colonies := 0L
]

matrix_long








# ============================================================
# 1. Identify qualifying concurrent intervals for each family
# ============================================================

family_time_support <- unique(
  sig_sweep_events[
    cluster_id %chin% matrix_rows &
      !is.na(Genus) &
      !is.na(Month_t1) &
      !is.na(Month_t2),
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2
    )
  ]
)[
  ,
  .(
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(Genus),
    n_colonies = uniqueN(colony_id)
  ),
  by = .(
    cluster_id,
    Month_t1,
    Month_t2
  )
][
  n_MAGs >= 2 &
    n_genera >= 2 &
    n_colonies >= 2
]

family_time_support[
  ,
  timepair := paste(
    Month_t1,
    Month_t2,
    sep = " -> "
  )
]

# Pick the strongest qualifying concurrent interval per family.
# Tie-break: genera, then MAGs, then colonies.
focal_timepair <- family_time_support[
  order(
    cluster_id,
    -n_genera,
    -n_MAGs,
    -n_colonies,
    Month_t1,
    Month_t2
  ),
  .SD[1],
  by = cluster_id
][
  ,
  .(
    cluster_id,
    focal_timepair = timepair,
    focal_n_MAGs = n_MAGs,
    focal_n_genera = n_genera,
    focal_n_colonies = n_colonies
  )
]

focal_timepair




# ============================================================
# 2. Sweep support in each family x genus
# ============================================================

sweep_cell_dt <- unique(
  sig_sweep_events[
    cluster_id %chin% matrix_rows,
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2
    )
  ]
)

sweep_cell_dt[
  ,
  timepair := paste(
    Month_t1,
    Month_t2,
    sep = " -> "
  )
]

# Attach each family's focal concurrent interval
sweep_cell_dt <- merge(
  sweep_cell_dt,
  focal_timepair[
    ,
    .(
      cluster_id,
      focal_timepair
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

# Cell summaries:
# - total number swept MAGs
# - number swept MAGs specifically in focal interval
cell_sweep_summary <- sweep_cell_dt[
  ,
  .(
    n_swept_MAGs = uniqueN(MAG_id),
    n_swept_colonies = uniqueN(colony_id),
    
    n_focal_MAGs = uniqueN(
      MAG_id[
        timepair == focal_timepair
      ]
    ),
    
    n_focal_colonies = uniqueN(
      colony_id[
        timepair == focal_timepair
      ]
    )
  ),
  by = .(
    cluster_id,
    Genus,
    focal_timepair
  )
]

cell_sweep_summary[
  ,
  is_concurrent_sweep :=
    n_focal_MAGs > 0
]




# ============================================================
# 3. Homolog presence + pN/pS + sweeps
# ============================================================

matrix_genera <- c(
  "apilactobacillus",
  "bartonella",
  "bifidobacterium",
  "commensalibacter",
  "frischella",
  "gilliamella",
  "lactobacillus",
  "snodgrassella"
)

matrix_grid <- CJ(
  cluster_id = matrix_rows,
  Genus = matrix_genera,
  unique = TRUE
)

matrix_plot_dt <- merge(
  matrix_grid,
  family_presence[
    cluster_id %chin% matrix_rows &
      Genus %chin% matrix_genera
  ],
  by = c(
    "cluster_id",
    "Genus"
  ),
  all.x = TRUE
)

matrix_plot_dt <- merge(
  matrix_plot_dt,
  cell_sweep_summary,
  by = c(
    "cluster_id",
    "Genus"
  ),
  all.x = TRUE
)

matrix_plot_dt <- merge(
  matrix_plot_dt,
  focal_timepair[
    ,
    .(
      cluster_id,
      focal_n_MAGs,
      focal_n_genera,
      focal_n_colonies
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

# Fill missing values
matrix_plot_dt[
  is.na(n_homologs),
  n_homologs := 0L
]

matrix_plot_dt[
  is.na(any_pnps),
  any_pnps := FALSE
]

for (cc in c(
  "n_swept_MAGs",
  "n_swept_colonies",
  "n_focal_MAGs",
  "n_focal_colonies"
)) {
  matrix_plot_dt[
    is.na(get(cc)),
    (cc) := 0L
  ]
}

matrix_plot_dt[
  ,
  homolog_present := n_homologs > 0
]

matrix_plot_dt[
  ,
  has_sweep := n_swept_MAGs > 0
]

matrix_plot_dt[
  ,
  is_concurrent_sweep :=
    n_focal_MAGs > 0
]





# ============================================================
# 3. Homolog presence + pN/pS + sweeps
# ============================================================

matrix_genera <- c(
  "apilactobacillus",
  "bartonella",
  "bifidobacterium",
  "commensalibacter",
  "frischella",
  "gilliamella",
  "lactobacillus",
  "snodgrassella"
)

matrix_grid <- CJ(
  cluster_id = matrix_rows,
  Genus = matrix_genera,
  unique = TRUE
)

matrix_plot_dt <- merge(
  matrix_grid,
  family_presence[
    cluster_id %chin% matrix_rows &
      Genus %chin% matrix_genera
  ],
  by = c(
    "cluster_id",
    "Genus"
  ),
  all.x = TRUE
)

matrix_plot_dt <- merge(
  matrix_plot_dt,
  cell_sweep_summary,
  by = c(
    "cluster_id",
    "Genus"
  ),
  all.x = TRUE
)

matrix_plot_dt <- merge(
  matrix_plot_dt,
  focal_timepair[
    ,
    .(
      cluster_id,
      focal_n_MAGs,
      focal_n_genera,
      focal_n_colonies
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

# Fill missing values
matrix_plot_dt[
  is.na(n_homologs),
  n_homologs := 0L
]

matrix_plot_dt[
  is.na(any_pnps),
  any_pnps := FALSE
]

for (cc in c(
  "n_swept_MAGs",
  "n_swept_colonies",
  "n_focal_MAGs",
  "n_focal_colonies"
)) {
  matrix_plot_dt[
    is.na(get(cc)),
    (cc) := 0L
  ]
}

matrix_plot_dt[
  ,
  homolog_present := n_homologs > 0
]

matrix_plot_dt[
  ,
  has_sweep := n_swept_MAGs > 0
]

matrix_plot_dt[
  ,
  is_concurrent_sweep :=
    n_focal_MAGs > 0
]






# ============================================================
# 4. Provisional row ordering by strength
# ============================================================

row_order <- unique(
  matrix_plot_dt[
    ,
    .(
      focal_n_genera = max(
        focal_n_genera,
        na.rm = TRUE
      ),
      focal_n_MAGs = max(
        focal_n_MAGs,
        na.rm = TRUE
      ),
      focal_n_colonies = max(
        focal_n_colonies,
        na.rm = TRUE
      ),
      pnps_genera = uniqueN(
        Genus[
          any_pnps == TRUE
        ]
      )
    ),
    by = cluster_id
  ][
    order(
      -focal_n_genera,
      -focal_n_MAGs,
      -focal_n_colonies,
      -pnps_genera
    ),
    cluster_id
  ]
)

matrix_plot_dt[
  ,
  cluster_short :=
    substr(
      cluster_id,
      5,
      14
    )
]

short_order <- substr(
  row_order,
  5,
  14
)

# Reverse so strongest family appears at top
matrix_plot_dt[
  has_sweep == TRUE,
  .(
    cluster_short,
    Genus,
    n_swept_MAGs,
    n_focal_MAGs,
    is_concurrent_sweep,
    focal_timepair,
    any_pnps,
    plot_size
  )
][order(cluster_short, Genus)]

matrix_plot_dt[
  ,
  plot_size := fifelse(
    is_concurrent_sweep,
    n_focal_MAGs,
    n_swept_MAGs
  )
]

# sanity check
matrix_plot_dt[
  ,
  .(
    min_plot_size = min(plot_size),
    max_plot_size = max(plot_size),
    n_nonzero = sum(plot_size > 0)
  )
]



# ============================================================
# 5. Matrix
# ============================================================

p_matrix <- ggplot(
  matrix_plot_dt,
  aes(
    x = Genus,
    y = display_label
  )
) +
  
  # ----------------------------------------------------------
# Homolog presence
# ----------------------------------------------------------
geom_tile(
  data = matrix_plot_dt[
    homolog_present == TRUE
  ],
  fill = "grey85",
  color = "white",
  linewidth = 0.5
) +
  
  # ----------------------------------------------------------
# Non-focal 0D sweeps:
# visible, but deliberately muted
# ----------------------------------------------------------
geom_point(
  data = matrix_plot_dt[
    has_sweep == TRUE &
      is_concurrent_sweep == FALSE
  ],
  aes(
    size = n_swept_MAGs
  ),
  shape = 21,
  fill = "grey65",
  color = "grey30",
  stroke = 0.4
) +
  
  # ----------------------------------------------------------
# Focal concurrent sweeps:
# color = exact concurrent time interval
# size = number of MAGs sweeping in that interval
# ----------------------------------------------------------
geom_point(
  data = matrix_plot_dt[
    is_concurrent_sweep == TRUE
  ],
  aes(
    size = n_focal_MAGs,
    fill = focal_timepair
  ),
  shape = 21,
  color = "black",
  stroke = 0.4
) +
  
  # ----------------------------------------------------------
# pN/pS support WITH a sweep:
# thicker outline/ring
# ----------------------------------------------------------
geom_point(
  data = matrix_plot_dt[
    any_pnps == TRUE &
      has_sweep == TRUE
  ],
  aes(
    size = plot_size
  ),
  shape = 21,
  fill = NA,
  color = "black",
  stroke = 1.4
) +
  
  # ----------------------------------------------------------
# pN/pS support but NO observed sweep:
# small open circle
# ----------------------------------------------------------
geom_point(
  data = matrix_plot_dt[
    any_pnps == TRUE &
      has_sweep == FALSE
  ],
  shape = 21,
  fill = NA,
  color = "black",
  stroke = 1,
  size = 2.4
) +
  
  scale_size_continuous(
    range = c(3, 7),
    breaks = sort(
      unique(
        matrix_plot_dt[
          plot_size > 0,
          plot_size
        ]
      )
    ),
    name = "Swept MAGs"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = "Concurrent sweep interval"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.key.height = unit(
      0.45,
      "cm"
    )
  )

p_matrix












# ============================================================
# Sweep events by degeneracy
# ============================================================

zeroD_events <- unique(
  sweep_dt[
    degeneracy == "0D",
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

fourD_events <- unique(
  sweep_dt[
    degeneracy == "4D",
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

# Attach MMseq family + taxonomy
gene_meta <- unique(
  gene_flags[
    ,
    .(
      gene_key,
      cluster_id,
      MAG_id,
      Genus
    )
  ]
)

zeroD_events <- merge(
  zeroD_events,
  gene_meta,
  by = c("gene_key", "MAG_id"),
  all.x = TRUE
)

fourD_events <- merge(
  fourD_events,
  gene_meta,
  by = c("gene_key", "MAG_id"),
  all.x = TRUE
)

zeroD_events[, sweep_class := "0D"]
fourD_events[, sweep_class := "4D"]



rbind(
  zeroD_events[
    ,
    .(
      sweep_class = "0D",
      n_events = .N,
      n_genes = uniqueN(gene_key),
      n_MAGs = uniqueN(MAG_id),
      n_clusters = uniqueN(cluster_id)
    )
  ],
  fourD_events[
    ,
    .(
      sweep_class = "4D",
      n_events = .N,
      n_genes = uniqueN(gene_key),
      n_MAGs = uniqueN(MAG_id),
      n_clusters = uniqueN(cluster_id)
    )
  ]
)



gene_sweep_overlap <- gene_flags[
  ,
  .(
    gene_key,
    MAG_id,
    cluster_id,
    Genus
  )
]

gene_sweep_overlap[
  ,
  has_0D := gene_key %chin% zeroD_events$gene_key
]

gene_sweep_overlap[
  ,
  has_4D := gene_key %chin% fourD_events$gene_key
]

table(
  gene_sweep_overlap$has_0D,
  gene_sweep_overlap$has_4D
)

fisher.test(
  table(
    gene_sweep_overlap$has_0D,
    gene_sweep_overlap$has_4D
  )
)



gene_sweep_overlap[
  ,
  .(
    P_4D_given_0D =
      mean(has_4D[has_0D]),
    
    P_4D_given_no0D =
      mean(has_4D[!has_0D]),
    
    P_0D_given_4D =
      mean(has_0D[has_4D])
  )
]





zeroD_tp <- unique(
  zeroD_events[
    ,
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

fourD_tp <- unique(
  fourD_events[
    ,
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

same_gene_same_interval <- merge(
  zeroD_tp,
  fourD_tp,
  by = c(
    "gene_key",
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  )
)

nrow(same_gene_same_interval)

uniqueN(
  same_gene_same_interval$gene_key
)

nrow(same_gene_same_interval) /
  nrow(zeroD_tp)











# cluster_interp should already contain annotation summaries per MMseq family
matrix_labels <- cluster_interp[
  cluster_id %chin% matrix_rows,
  .(
    cluster_id,
    PFAM_labels,
    GO_labels
  )
]

# crude first-pass label:
matrix_labels[
  ,
  display_label := fifelse(
    !is.na(PFAM_labels) &
      PFAM_labels != "",
    PFAM_labels,
    "Uncharacterized protein"
  )
]

# if PFAM_labels contains multiple entries, keep first for now
matrix_labels[
  ,
  display_label := sub(
    ";.*$",
    "",
    display_label
  )
]

matrix_plot_dt <- merge(
  matrix_plot_dt,
  matrix_labels[
    ,
    .(
      cluster_id,
      display_label
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_sweeps <- merge(
  matrix_sweeps,
  matrix_labels[
    ,
    .(
      cluster_id,
      display_label
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)



label_order <- matrix_plot_dt[
  ,
  .SD[1],
  by = cluster_id
][
  match(row_order, cluster_id),
  display_label
]

matrix_plot_dt[
  ,
  display_label := factor(
    display_label,
    levels = rev(label_order)
  )
]

matrix_sweeps[
  ,
  display_label := factor(
    display_label,
    levels = levels(matrix_plot_dt$display_label)
  )
]




















# ============================================================
# Every 0D sweep represented in the matrix
# ============================================================

matrix_sweep_inventory <- unique(
  matrix_sweeps[
    ,
    .(
      cluster_id,
      Genus,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      timepair,
      concurrent,
      focal_timepair
    )
  ]
)

matrix_sweep_inventory[
  order(
    cluster_id,
    Genus,
    MAG_id,
    colony_id,
    Month_t1
  )
]






matrix_sweep_genes <- unique(
  sig_sweep_events[
    cluster_id %chin% matrix_rows,
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2
    )
  ]
)

matrix_sweep_genes[
  ,
  timepair := paste(
    Month_t1,
    Month_t2,
    sep = " -> "
  )
]

matrix_sweep_genes <- merge(
  matrix_sweep_genes,
  focal_timepair[
    ,
    .(
      cluster_id,
      focal_timepair
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_sweep_genes[
  ,
  concurrent := timepair == focal_timepair
]




matrix_sweep_genes[
  concurrent == TRUE,
  .(
    cluster_id,
    Genus,
    MAG_id,
    colony_id,
    gene_key,
    Month_t1,
    Month_t2
  )
][
  order(
    cluster_id,
    Genus,
    MAG_id
  )
]

matrix_allele_rows <- merge(
  sweep_dt[
    degeneracy == "0D"
  ],
  matrix_sweep_genes[
    ,
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      Genus,
      concurrent,
      focal_timepair
    )
  ],
  by = c(
    "gene_key",
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  ),
  all = FALSE
)


nrow(matrix_allele_rows)

names(matrix_allele_rows)


concurrent_alleles <- matrix_allele_rows[
  concurrent == TRUE
]

concurrent_alleles[
  order(
    cluster_id,
    Genus,
    MAG_id,
    colony_id,
    Month_t1
  )
]





trajectory_candidates <- unique(
  concurrent_alleles[
    ,
    .(
      cluster_id,
      Genus,
      MAG_id,
      colony_id,
      gene_key,
      corresponding_gene_call,
      Month_t1,
      Month_t2
    )
  ]
)

trajectory_candidates[
  order(
    cluster_id,
    Genus,
    MAG_id
  )
]





# ============================================================
# MAG-level 0D sweep burden
# ============================================================

mag_sweep_burden <- unique(
  zeroD_events[
    ,
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)[
  ,
  .(
    MAG_n_0D_sweep_events = .N,
    MAG_n_0D_swept_genes = uniqueN(gene_key),
    MAG_n_colonies_with_0D_sweep = uniqueN(colony_id),
    MAG_n_timepairs_with_0D_sweep = uniqueN(
      paste(Month_t1, Month_t2, sep = "->")
    )
  ),
  by = MAG_id
]


trajectory_candidates <- merge(
  trajectory_candidates,
  mag_sweep_burden,
  by = "MAG_id",
  all.x = TRUE
)




pop_sweep_burden <- unique(
  zeroD_events[
    ,
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)[
  ,
  .(
    population_n_0D_sweep_events = .N,
    population_n_0D_swept_genes = uniqueN(gene_key)
  ),
  by = .(
    MAG_id,
    colony_id
  )
]

trajectory_candidates <- merge(
  trajectory_candidates,
  pop_sweep_burden,
  by = c(
    "MAG_id",
    "colony_id"
  ),
  all.x = TRUE
)


trajectory_candidates








core_genera <- c(
  "lactobacillus",
  "gilliamella",
  "snodgrassella",
  "bartonella",
  "commensalibacter",
  "bifidobacterium",
  "frischella"
)

# ============================================================
# UpSet 1: all protein families containing >=1 0D sweep
# by genus
# ============================================================

all_0D_family_genus <- unique(
  zeroD_clustered[
    cluster_id %chin% testable_ids &
      Genus %chin% core_genera,
    .(
      cluster_id,
      Genus
    )
  ]
)

# Convert to family x genus presence/absence
all_0D_upset <- dcast(
  all_0D_family_genus[
    ,
    present := TRUE
  ],
  cluster_id ~ Genus,
  value.var = "present",
  fill = FALSE
)

# Make sure every genus column exists
for (g in core_genera) {
  if (!g %in% names(all_0D_upset)) {
    all_0D_upset[, (g) := FALSE]
  }
}

# sanity checks
nrow(all_0D_upset)

all_0D_family_genus[
  ,
  .(
    n_0D_families = uniqueN(cluster_id)
  ),
  by = Genus
][order(-n_0D_families)]





p_upset_all_0D <- upset(
  as.data.frame(all_0D_upset),
  intersect = core_genera,
  min_size = 1,
  width_ratio = 0.2,
  
  base_annotations = list(
    "Protein families" =
      intersection_size(
        counts = TRUE
      )
  ),
  
  set_sizes = upset_set_size()
) +
  labs(
    title = "Protein families containing 0D sweeps across core genera"
  )

p_upset_all_0D





# ============================================================
# UpSet 2: concurrent recurrent 0D-sweep families
# by genus
# ============================================================

concurrent_ids <- unique(
  obs_concurrent_by_time$cluster_id
)

length(concurrent_ids)
# 31




concurrent_0D_family_genus <- unique(
  zeroD_clustered[
    cluster_id %chin% concurrent_ids &
      Genus %chin% core_genera,
    .(
      cluster_id,
      Genus
    )
  ]
)


concurrent_0D_upset <- dcast(
  concurrent_0D_family_genus[
    ,
    present := TRUE
  ],
  cluster_id ~ Genus,
  value.var = "present",
  fill = FALSE
)

for (g in core_genera) {
  if (!g %in% names(concurrent_0D_upset)) {
    concurrent_0D_upset[, (g) := FALSE]
  }
}

nrow(concurrent_0D_upset)

concurrent_0D_family_genus[
  ,
  .(
    n_concurrent_families = uniqueN(cluster_id)
  ),
  by = Genus
][order(-n_concurrent_families)]






p_upset_concurrent <- upset(
  as.data.frame(concurrent_0D_upset),
  intersect = core_genera,
  min_size = 1,
  width_ratio = 0.2,
  
  base_annotations = list(
    "Protein families" =
      intersection_size(
        counts = TRUE
      )
  ),
  
  set_sizes = upset_set_size()
) +
  labs(
    title = "Concurrent recurrent 0D-sweep families across core genera"
  )

p_upset_concurrent


ANN_RDS <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/gene_level_annotations_all_colonies.rds"
ann_dt   <- readRDS(ANN_RDS)                # annotations (list-columns)



# ============================================================
# Broad annotation labels from ANN_RDS
# ============================================================

if (!is.data.table(ann_dt)) {
  ann_dt <- as.data.table(ann_dt)
}

ann_dt[
  ,
  corresponding_gene_call :=
    as.character(corresponding_gene_call)
]

# Helper: collapse list-column values safely
collapse_unique <- function(x) {
  
  x <- unlist(
    x,
    recursive = TRUE,
    use.names = FALSE
  )
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      x != "" &
      x != "NA"
  ]
  
  x <- unique(x)
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(
    x,
    collapse = "; "
  )
}


# Full annotation strings
ann_dt[
  ,
  pfam_str := vapply(
    pfam_label_list,
    collapse_unique,
    character(1)
  )
]

ann_dt[
  ,
  GO_str := vapply(
    GO_slim_name_list,
    collapse_unique,
    character(1)
  )
]



# ============================================================
# One compact display label per gene
# ============================================================

first_annotation <- function(x) {
  
  if (is.na(x) || x == "") {
    return(NA_character_)
  }
  
  trimws(
    strsplit(
      x,
      ";",
      fixed = TRUE
    )[[1]][1]
  )
}


ann_dt[
  ,
  pfam_primary := vapply(
    pfam_str,
    first_annotation,
    character(1)
  )
]

ann_dt[
  ,
  GO_primary := vapply(
    GO_str,
    first_annotation,
    character(1)
  )
]

ann_dt[
  ,
  broad_function := fifelse(
    !is.na(pfam_primary),
    pfam_primary,
    fifelse(
      !is.na(GO_primary),
      GO_primary,
      "Uncharacterized protein"
    )
  )
]



# ============================================================
# Collapse annotation table to one record per gene call
# ============================================================

ann_gene <- ann_dt[
  ,
  .(
    pfam_str = collapse_unique(pfam_str),
    GO_str = collapse_unique(GO_str),
    
    pfam_primary = {
      z <- pfam_primary[!is.na(pfam_primary)]
      if (length(z)) z[1] else NA_character_
    },
    
    GO_primary = {
      z <- GO_primary[!is.na(GO_primary)]
      if (length(z)) z[1] else NA_character_
    },
    
    broad_function = {
      z <- broad_function[
        !is.na(broad_function) &
          broad_function != "Uncharacterized protein"
      ]
      
      if (length(z)) {
        z[1]
      } else {
        "Uncharacterized protein"
      }
    }
  ),
  by = corresponding_gene_call
]

stopifnot(
  uniqueN(ann_gene$corresponding_gene_call) ==
    nrow(ann_gene)
)




ann_sub <- ann_gene[
  ,
  .(
    corresponding_gene_call,
    pfam_str,
    GO_str,
    pfam_primary,
    GO_primary,
    broad_function
  )
]

gene_dt <- merge(
  gene_dt,
  ann_sub,
  by = "corresponding_gene_call",
  all.x = TRUE
)

sweep_dt <- merge(
  sweep_dt,
  ann_sub,
  by = "corresponding_gene_call",
  all.x = TRUE
)




zeroD_annotated <- sweep_dt[
  degeneracy == "0D"
]

zeroD_annotated[
  ,
  .(
    n_sweep_rows = .N,
    n_genes = uniqueN(corresponding_gene_call),
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(Genus.y)
  ),
  by = broad_function.y
][
  order(-n_genes)
][1:30]


cluster_function_labels <- gene_flags[
  ,
  .(
    corresponding_gene_call,
    cluster_id
  )
][
  ann_gene,
  on = "corresponding_gene_call"
][
  !is.na(cluster_id),
  .(
    n_annotated_genes = .N,
    
    consensus_function = {
      z <- broad_function[
        !is.na(broad_function) &
          broad_function != "Uncharacterized protein"
      ]
      
      if (length(z) == 0) {
        "Uncharacterized protein"
      } else {
        names(
          sort(
            table(z),
            decreasing = TRUE
          )
        )[1]
      }
    },
    
    consensus_fraction = {
      z <- broad_function[
        !is.na(broad_function) &
          broad_function != "Uncharacterized protein"
      ]
      
      if (length(z) == 0) {
        NA_real_
      } else {
        max(table(z)) / length(z)
      }
    }
  ),
  by = cluster_id
]

















#### BUILDING

###

core_genera <- c(
  "lactobacillus",
  "gilliamella",
  "snodgrassella",
  "bartonella",
  "commensalibacter",
  "bifidobacterium",
  "frischella"
)

# Attach consensus annotation to cluster IDs
cluster_label_lookup <- unique(
  cluster_function_labels[
    ,
    .(
      cluster_id,
      consensus_function,
      consensus_fraction
    )
  ]
)

# ============================================================
# 1. ALL 0D-SWEEP FAMILIES BY GENUS
# ============================================================

all_0D_family_genus <- unique(
  zeroD_clustered[
    cluster_id %chin% testable_ids &
      Genus %chin% core_genera,
    .(
      cluster_id,
      Genus
    )
  ]
)

all_0D_upset <- dcast(
  all_0D_family_genus[
    ,
    present := TRUE
  ],
  cluster_id ~ Genus,
  value.var = "present",
  fill = FALSE
)

for (g in core_genera) {
  if (!g %in% names(all_0D_upset)) {
    all_0D_upset[, (g) := FALSE]
  }
}




all_0D_upset <- merge(
  all_0D_upset,
  cluster_label_lookup,
  by = "cluster_id",
  all.x = TRUE
)

all_0D_upset[
  is.na(consensus_function),
  consensus_function := "Uncharacterized protein"
]





p_upset_all_0D <- upset(
  as.data.frame(all_0D_upset),
  intersect = core_genera,
  min_size = 1,
  width_ratio = 0.20,
  
  base_annotations = list(
    "Protein families" =
      intersection_size(
        counts = TRUE
      )
  ),
  
  set_sizes =
    upset_set_size()
) +
  labs(
    title = "Protein families with nonsynonymous sweeps across core genera"
  ) +
  theme(
    text = element_text(size = 10)
  )

p_upset_all_0D


all_0D_upset[
  lactobacillus == TRUE &
    gilliamella == TRUE &
    snodgrassella == TRUE,
  .(
    cluster_id,
    consensus_function,
    consensus_fraction
  )
]




# ============================================================
# 2. CONCURRENT 0D-SWEEP FAMILIES BY GENUS
# ============================================================

concurrent_ids <- unique(
  obs_concurrent_by_time$cluster_id
)

concurrent_0D_family_genus <- unique(
  zeroD_clustered[
    cluster_id %chin% concurrent_ids &
      Genus %chin% core_genera,
    .(
      cluster_id,
      Genus
    )
  ]
)

concurrent_0D_upset <- dcast(
  concurrent_0D_family_genus[
    ,
    present := TRUE
  ],
  cluster_id ~ Genus,
  value.var = "present",
  fill = FALSE
)

for (g in core_genera) {
  if (!g %in% names(concurrent_0D_upset)) {
    concurrent_0D_upset[, (g) := FALSE]
  }
}

concurrent_0D_upset <- merge(
  concurrent_0D_upset,
  cluster_label_lookup,
  by = "cluster_id",
  all.x = TRUE
)

concurrent_0D_upset[
  is.na(consensus_function),
  consensus_function :=
    "Uncharacterized protein"
]


p_upset_concurrent <- upset(
  as.data.frame(concurrent_0D_upset),
  intersect = core_genera,
  min_size = 1,
  width_ratio = 0.20,
  
  base_annotations = list(
    "Protein families" =
      intersection_size(
        counts = TRUE
      )
  ),
  
  set_sizes =
    upset_set_size()
) +
  labs(
    title = "Concurrent nonsynonymous-sweep families across core genera"
  ) +
  theme(
    text = element_text(size = 10)
  )

p_upset_concurrent



matrix_labels <- cluster_label_lookup[
  cluster_id %chin% matrix_rows
]

matrix_labels[
  is.na(consensus_function) |
    consensus_function == "",
  consensus_function :=
    "Uncharacterized protein"
]




matrix_labels[
  ,
  n_same_label := .N,
  by = consensus_function
]

matrix_labels[
  ,
  display_label := fifelse(
    n_same_label > 1,
    paste0(
      consensus_function,
      " [",
      substr(cluster_id, 5, 10),
      "]"
    ),
    consensus_function
  )
]

matrix_plot_dt[
  ,
  display_label := NULL
]

matrix_sweeps[
  ,
  display_label := NULL
]

matrix_plot_dt <- merge(
  matrix_plot_dt,
  matrix_labels[
    ,
    .(
      cluster_id,
      display_label
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_sweeps <- merge(
  matrix_sweeps,
  matrix_labels[
    ,
    .(
      cluster_id,
      display_label
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)



label_order <- matrix_labels[
  match(row_order, cluster_id),
  display_label
]

matrix_plot_dt[
  ,
  display_label := factor(
    display_label,
    levels = rev(label_order)
  )
]

matrix_sweeps[
  ,
  display_label := factor(
    display_label,
    levels = levels(
      matrix_plot_dt$display_label
    )
  )
]

matrix_plot_dt[
  ,
  Genus := factor(
    Genus,
    levels = matrix_genera
  )
]

genus_lookup <- data.table(
  Genus = matrix_genera,
  genus_x = seq_along(matrix_genera)
)

# matrix_sweeps may already contain these;
# rebuild only if necessary
if (!"genus_x" %in% names(matrix_sweeps)) {
  
  matrix_sweeps <- merge(
    matrix_sweeps,
    genus_lookup,
    by = "Genus",
    all.x = TRUE
  )
}

if (!"plot_x" %in% names(matrix_sweeps)) {
  
  matrix_sweeps[
    ,
    sweep_index := seq_len(.N),
    by = .(
      cluster_id,
      Genus
    )
  ]
  
  matrix_sweeps[
    ,
    n_cell_sweeps := .N,
    by = .(
      cluster_id,
      Genus
    )
  ]
  
  matrix_sweeps[
    ,
    x_offset :=
      (
        sweep_index -
          (n_cell_sweeps + 1) / 2
      ) * 0.16
  ]
  
  matrix_sweeps[
    ,
    plot_x := genus_x + x_offset
  ]
}



p_matrix <- ggplot() +
  
  # homolog presence
  geom_tile(
    data = matrix_plot_dt[
      homolog_present == TRUE
    ],
    aes(
      x = as.numeric(Genus),
      y = display_label
    ),
    fill = "grey85",
    color = "white",
    linewidth = 0.5
  ) +
  
  # pN/pS support at family x genus level
  geom_tile(
    data = matrix_plot_dt[
      any_pnps == TRUE
    ],
    aes(
      x = as.numeric(Genus),
      y = display_label
    ),
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # every observed 0D sweep
  geom_point(
    data = matrix_sweeps,
    aes(
      x = plot_x,
      y = display_label,
      fill = timepair
    ),
    shape = 21,
    size = 2.9,
    color = "black",
    stroke = 0.35
  ) +
  
  # concurrent sweep marker
  geom_point(
    data = matrix_sweeps[
      concurrent == TRUE
    ],
    aes(
      x = plot_x,
      y = display_label
    ),
    shape = 16,
    size = 1.3,
    color = "black"
  ) +
  
  scale_x_continuous(
    breaks = seq_along(matrix_genera),
    labels = matrix_genera
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = "Sweep interval"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 9
    ),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.key.height =
      grid::unit(0.42, "cm")
  )

p_matrix








# ============================================================
# Concurrent 0D UpSet:
# stacked intersection bars by protein-function label
# ============================================================

# clean annotation labels
concurrent_0D_upset[
  is.na(consensus_function) |
    consensus_function == "",
  consensus_function := "Uncharacterized protein"
]

p_upset_concurrent <- upset(
  as.data.frame(concurrent_0D_upset),
  
  intersect = core_genera,
  min_size = 1,
  width_ratio = 0.20,
  
  base_annotations = list(
    
    "Protein families" = (
      ggplot(
        mapping = aes(
          fill = consensus_function
        )
      ) +
        geom_bar() +
        labs(
          y = "Protein families",
          fill = "Protein function"
        ) +
        theme(
          axis.title.x = element_blank()
        )
    )
  ),
  
  set_sizes =
    upset_set_size()
) +
  
  labs(
    title =
      "Concurrent nonsynonymous-sweep families across core genera"
  )

p_upset_concurrent



# Count functions among the 31 concurrent families
function_counts <- concurrent_0D_upset[
  ,
  .N,
  by = consensus_function
][
  order(-N)
]

function_counts





keep_functions <- function_counts[
  N >= 2 &
    consensus_function != "Uncharacterized protein",
  consensus_function
]

concurrent_0D_upset[
  ,
  function_group := fifelse(
    consensus_function == "Uncharacterized protein",
    "Uncharacterized protein",
    fifelse(
      consensus_function %chin% keep_functions,
      consensus_function,
      "Other"
    )
  )
]




p_upset_concurrent <- upset(
  as.data.frame(concurrent_0D_upset),
  
  intersect = core_genera,
  min_size = 1,
  width_ratio = 0.20,
  
  base_annotations = list(
    
    "Protein families" = (
      ggplot(
        mapping = aes(
          fill = function_group
        )
      ) +
        geom_bar() +
        labs(
          y = "Protein families",
          fill = "Protein function"
        ) +
        theme(
          axis.title.x = element_blank()
        )
    )
  ),
  
  set_sizes =
    upset_set_size()
)

p_upset_concurrent






# ============================================================
# Observed concurrent ABC-transporter families
# ============================================================

obs_concurrent_annot <- merge(
  unique(
    obs_concurrent_by_time[
      ,
      .(cluster_id)
    ]
  ),
  cluster_label_lookup[
    ,
    .(
      cluster_id,
      consensus_function
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

obs_n_ABC_concurrent <- obs_concurrent_annot[
  consensus_function == "ABC_tran",
  uniqueN(cluster_id)
]

obs_n_concurrent <- uniqueN(
  obs_concurrent_annot$cluster_id
)

obs_n_ABC_concurrent
obs_n_concurrent

# should be approximately:
# 14
# 31

obs_fraction_ABC <-
  obs_n_ABC_concurrent /
  obs_n_concurrent

obs_fraction_ABC


# ============================================================
# Precompute lookup objects ONCE
# ============================================================

cluster_to_function <- setNames(
  cluster_label_lookup$consensus_function,
  cluster_label_lookup$cluster_id
)

callable_dt[
  ,
  stratum_id := paste(
    MAG_id,
    colony_id,
    sep = "|"
  )
]

callable_pools <- split(
  callable_dt,
  callable_dt$stratum_id
)

timepair_burden[
  ,
  stratum_id := paste(
    MAG_id,
    colony_id,
    sep = "|"
  )
]
N_PERM_ABC <- 5000L

null_n_ABC_concurrent <- integer(N_PERM_ABC)
null_fraction_ABC_concurrent <- numeric(N_PERM_ABC)

set.seed(12345)

for (b in seq_len(N_PERM_ABC)) {
  
  perm_list <- vector(
    "list",
    nrow(timepair_burden)
  )
  
  for (i in seq_len(nrow(timepair_burden))) {
    
    t1_i <- timepair_burden$Month_t1[i]
    t2_i <- timepair_burden$Month_t2[i]
    k_i  <- timepair_burden$k[i]
    
    # FAST lookup instead of repeatedly filtering callable_dt
    pool <- callable_pools[[timepair_burden$stratum_id[i]]]
    
    if (is.null(pool)) {
      stop(
        paste(
          "No callable pool for",
          timepair_burden$stratum_id[i]
        )
      )
    }
    
    if (nrow(pool) < k_i) {
      stop(
        paste(
          "Insufficient callable genes for",
          timepair_burden$stratum_id[i]
        )
      )
    }
    
    idx <- sample(
      seq_len(nrow(pool)),
      size = k_i,
      replace = FALSE,
      prob = pool$effective_length_0D_sites
    )
    
    picked <- pool[
      idx,
      .(
        cluster_id,
        MAG_id,
        colony_id,
        Genus
      )
    ]
    
    picked[
      ,
      `:=`(
        Month_t1 = t1_i,
        Month_t2 = t2_i
      )
    ]
    
    perm_list[[i]] <- picked
  }
  
  perm_events <- rbindlist(
    perm_list,
    use.names = TRUE
  )
  
  # Same locked protein-family universe
  perm_events <- perm_events[
    cluster_id %chin% testable_ids &
      !is.na(Genus)
  ]
  
  # Same concurrency definition
  perm_concurrent <- perm_events[
    ,
    .(
      n_MAGs = uniqueN(MAG_id),
      n_genera = uniqueN(Genus),
      n_colonies = uniqueN(colony_id)
    ),
    by = .(
      cluster_id,
      Month_t1,
      Month_t2
    )
  ][
    n_MAGs >= 2 &
      n_genera >= 2 &
      n_colonies >= 2
  ]
  
  # Unique concurrent protein families
  perm_concurrent_ids <- unique(
    perm_concurrent$cluster_id
  )
  
  n_total <- length(
    perm_concurrent_ids
  )
  
  # Fast annotation lookup
  functions <- cluster_to_function[
    perm_concurrent_ids
  ]
  
  n_ABC <- sum(
    functions == "ABC_tran",
    na.rm = TRUE
  )
  
  null_n_ABC_concurrent[b] <- n_ABC
  
  null_fraction_ABC_concurrent[b] <-
    if (n_total > 0) {
      n_ABC / n_total
    } else {
      NA_real_
    }
  
  if (b %% 100 == 0) {
    message(
      b,
      " / ",
      N_PERM_ABC
    )
  }
}



summary(null_n_ABC_concurrent)

quantile(
  null_n_ABC_concurrent,
  c(0.025, 0.5, 0.975, 0.99)
)

p_ABC_count <- (
  sum(
    null_n_ABC_concurrent >=
      obs_n_ABC_concurrent
  ) + 1
) /
  (N_PERM_ABC + 1)

p_ABC_count



summary(
  null_fraction_ABC_concurrent
)

p_ABC_fraction <- (
  sum(
    null_fraction_ABC_concurrent >=
      obs_fraction_ABC,
    na.rm = TRUE
  ) + 1
) /
  (
    sum(!is.na(
      null_fraction_ABC_concurrent
    )) + 1
  )

p_ABC_fraction


















# ============================================================
# Curate matrix rows:
# strongest families, max 2 per PFAM label
# ============================================================

matrix_candidate_rank <- matrix_candidates[
  cluster_id %chin% concurrent_ids
]

matrix_candidate_rank <- merge(
  matrix_candidate_rank,
  cluster_label_lookup[
    ,
    .(
      cluster_id,
      consensus_function,
      consensus_fraction
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_candidate_rank[
  is.na(consensus_function) |
    consensus_function == "",
  consensus_function := "Uncharacterized protein"
]

# Rank by strength of recurrent/concurrent evidence
setorder(
  matrix_candidate_rank,
  -same_transition_genera,
  -same_transition_MAGs,
  -same_transition_colonies,
  -n_pnps_genera
)

# Within each PFAM, retain at most the two strongest families
matrix_candidate_capped <- matrix_candidate_rank[
  ,
  head(.SD, 2),
  by = consensus_function
]

expected_ABC <- mean(null_n_ABC_concurrent)

fold_ABC <- obs_n_ABC_concurrent / expected_ABC

expected_ABC
fold_ABC




setorder(
  matrix_candidate_capped,
  -same_transition_genera,
  -same_transition_MAGs,
  -same_transition_colonies,
  -n_pnps_genera
)

matrix_rows <- head(
  matrix_candidate_capped$cluster_id,
  16
)

matrix_candidate_capped[
  cluster_id %chin% matrix_rows,
  .N,
  by = consensus_function
][order(-N)]
















# ============================================================
# 1. CURATE FINAL MATRIX CANDIDATES
#    - concurrent candidate families only
#    - max 2 per consensus PFAM label
#    - 16 rows total
# ============================================================

matrix_candidate_rank <- matrix_candidates[
  cluster_id %chin% concurrent_ids
]

matrix_candidate_rank <- merge(
  matrix_candidate_rank,
  cluster_label_lookup[
    ,
    .(
      cluster_id,
      consensus_function,
      consensus_fraction
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_candidate_rank[
  is.na(consensus_function) |
    consensus_function == "",
  consensus_function := "Uncharacterized protein"
]

# strongest concurrent evidence first
setorder(
  matrix_candidate_rank,
  -same_transition_genera,
  -same_transition_MAGs,
  -same_transition_colonies,
  -n_pnps_genera
)

# max 2 representatives per PFAM
# allow up to 3 uncharacterized families if desired
matrix_candidate_capped <- matrix_candidate_rank[
  ,
  if (consensus_function == "Uncharacterized protein") {
    head(.SD, 3)
  } else {
    head(.SD, 2)
  },
  by = consensus_function
]

setorder(
  matrix_candidate_capped,
  -same_transition_genera,
  -same_transition_MAGs,
  -same_transition_colonies,
  -n_pnps_genera
)

matrix_rows <- head(
  matrix_candidate_capped$cluster_id,
  16
)

# audit
matrix_candidate_capped[
  cluster_id %chin% matrix_rows,
  .N,
  by = consensus_function
][order(-N)]



# ============================================================
# 2. GENUS ORDER
# ============================================================

matrix_genera <- c(
  "apilactobacillus",
  "bartonella",
  "bifidobacterium",
  "commensalibacter",
  "frischella",
  "gilliamella",
  "lactobacillus",
  "snodgrassella"
)



# ============================================================
# 3. FOCAL CONCURRENT INTERVAL PER FAMILY
# ============================================================

family_time_support <- unique(
  sig_sweep_events[
    cluster_id %chin% matrix_rows &
      !is.na(Genus) &
      !is.na(Month_t1) &
      !is.na(Month_t2),
    .(
      cluster_id,
      gene_key,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2
    )
  ]
)[
  ,
  .(
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(Genus),
    n_colonies = uniqueN(colony_id)
  ),
  by = .(
    cluster_id,
    Month_t1,
    Month_t2
  )
][
  n_MAGs >= 2 &
    n_genera >= 2 &
    n_colonies >= 2
]

family_time_support[
  ,
  timepair := paste(
    Month_t1,
    Month_t2,
    sep = " -> "
  )
]

# strongest qualifying interval per family
focal_timepair <- family_time_support[
  order(
    cluster_id,
    -n_genera,
    -n_MAGs,
    -n_colonies,
    Month_t1,
    Month_t2
  ),
  .SD[1],
  by = cluster_id
][
  ,
  .(
    cluster_id,
    focal_timepair = timepair,
    focal_n_MAGs = n_MAGs,
    focal_n_genera = n_genera,
    focal_n_colonies = n_colonies
  )
]






# ============================================================
# 4. FAMILY x GENUS BACKGROUND GRID
# ============================================================

family_presence_matrix <- gene_flags[
  cluster_id %chin% matrix_rows &
    Genus %chin% matrix_genera,
  .(
    n_homologs = uniqueN(gene_key),
    any_pnps = any(any_pnps_gt1, na.rm = TRUE)
  ),
  by = .(
    cluster_id,
    Genus
  )
]

matrix_grid <- CJ(
  cluster_id = matrix_rows,
  Genus = matrix_genera,
  unique = TRUE
)

matrix_plot_dt <- merge(
  matrix_grid,
  family_presence_matrix,
  by = c("cluster_id", "Genus"),
  all.x = TRUE
)

matrix_plot_dt[
  is.na(n_homologs),
  n_homologs := 0L
]

matrix_plot_dt[
  is.na(any_pnps),
  any_pnps := FALSE
]

matrix_plot_dt[
  ,
  homolog_present := n_homologs > 0
]







# ============================================================
# 5. INDIVIDUAL 0D SWEEP SYMBOLS
# ============================================================

matrix_sweeps <- unique(
  sig_sweep_events[
    cluster_id %chin% matrix_rows &
      Genus %chin% matrix_genera,
    .(
      cluster_id,
      Genus,
      MAG_id,
      colony_id,
      gene_key,
      Month_t1,
      Month_t2
    )
  ]
)

matrix_sweeps[
  ,
  timepair := paste(
    Month_t1,
    Month_t2,
    sep = " -> "
  )
]

matrix_sweeps <- merge(
  matrix_sweeps,
  focal_timepair[
    ,
    .(
      cluster_id,
      focal_timepair
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_sweeps[
  ,
  concurrent :=
    !is.na(focal_timepair) &
    timepair == focal_timepair
]





# ============================================================
# 6. OFFSET MULTIPLE SWEEPS WITHIN EACH CELL
# ============================================================

genus_lookup <- data.table(
  Genus = matrix_genera,
  genus_x = seq_along(matrix_genera)
)

matrix_sweeps <- merge(
  matrix_sweeps,
  genus_lookup,
  by = "Genus",
  all.x = TRUE
)

matrix_sweeps[
  ,
  sweep_index := seq_len(.N),
  by = .(
    cluster_id,
    Genus
  )
]

matrix_sweeps[
  ,
  n_cell_sweeps := .N,
  by = .(
    cluster_id,
    Genus
  )
]

matrix_sweeps[
  ,
  x_offset :=
    (
      sweep_index -
        (n_cell_sweeps + 1) / 2
    ) * 0.16
]

matrix_sweeps[
  ,
  plot_x := genus_x + x_offset
]




# ============================================================
# 7. FINAL ROW LABELS
# ============================================================

matrix_labels <- unique(
  matrix_candidate_capped[
    cluster_id %chin% matrix_rows,
    .(
      cluster_id,
      consensus_function,
      consensus_fraction
    )
  ]
)

matrix_labels[
  ,
  duplicate_n := .N,
  by = consensus_function
]

matrix_labels[
  ,
  display_label := fifelse(
    duplicate_n > 1,
    paste0(
      consensus_function,
      " [",
      substr(cluster_id, 5, 10),
      "]"
    ),
    consensus_function
  )
]

matrix_plot_dt <- merge(
  matrix_plot_dt,
  matrix_labels[
    ,
    .(
      cluster_id,
      display_label
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

matrix_sweeps <- merge(
  matrix_sweeps,
  matrix_labels[
    ,
    .(
      cluster_id,
      display_label
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)






# ============================================================
# 8. ROW ORDER
# ============================================================

row_order <- matrix_candidate_capped[
  cluster_id %chin% matrix_rows
][
  order(
    -same_transition_genera,
    -same_transition_MAGs,
    -same_transition_colonies,
    -n_pnps_genera
  ),
  cluster_id
]

label_order <- matrix_labels[
  match(row_order, cluster_id),
  display_label
]

matrix_plot_dt[
  ,
  display_label := factor(
    display_label,
    levels = rev(label_order)
  )
]

matrix_sweeps[
  ,
  display_label := factor(
    display_label,
    levels = levels(matrix_plot_dt$display_label)
  )
]

matrix_plot_dt[
  ,
  Genus := factor(
    Genus,
    levels = matrix_genera
  )
]





# ============================================================
# 9. FINAL MATRIX FIGURE
# ============================================================

p_matrix_final <- ggplot() +
  
  # homolog present
  geom_tile(
    data = matrix_plot_dt[
      homolog_present == TRUE
    ],
    aes(
      x = as.numeric(Genus),
      y = display_label
    ),
    fill = "grey85",
    color = "white",
    linewidth = 0.5
  ) +
  
  # pN/pS > 1 support at family x genus level
  geom_tile(
    data = matrix_plot_dt[
      any_pnps == TRUE
    ],
    aes(
      x = as.numeric(Genus),
      y = display_label
    ),
    fill = NA,
    color = "black",
    linewidth = 0.85
  ) +
  
  # every observed 0D sweep
  geom_point(
    data = matrix_sweeps,
    aes(
      x = plot_x,
      y = display_label,
      fill = timepair
    ),
    shape = 21,
    size = 2.9,
    color = "black",
    stroke = 0.35
  ) +
  
  # central dot = member of focal concurrent event
  geom_point(
    data = matrix_sweeps[
      concurrent == TRUE
    ],
    aes(
      x = plot_x,
      y = display_label
    ),
    shape = 16,
    size = 1.25,
    color = "black"
  ) +
  
  scale_x_continuous(
    breaks = seq_along(matrix_genera),
    labels = tools::toTitleCase(matrix_genera),
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = "0D sweep interval"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 9
    ),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.height = unit(
      0.38,
      "cm"
    )
  )

p_matrix_final















# Replace this with the exact table used for Fig. 5B
fig5_dt <- sweep_transition_summary %>%
  filter(
    S > 0,
    is.finite(abs_delta_B),
    is.finite(B_t2)
  ) %>%
  mutate(
    resident_preserving =
      abs_delta_B <= 0.10 &
      B_t2 <= 0.20
  )

# ------------------------------------------------------------
# Observed effects
# ------------------------------------------------------------

transition_obs <- fig5_dt %>%
  summarise(
    n = n(),
    n_resident = sum(resident_preserving),
    estimate = mean(resident_preserving)
  )

population_obs <- fig5_dt %>%
  group_by(population_id) %>%
  summarise(
    resident_preserving =
      any(resident_preserving),
    .groups = "drop"
  ) %>%
  summarise(
    n = n(),
    n_resident = sum(resident_preserving),
    estimate = mean(resident_preserving)
  )





# Replace this with the exact table used for Fig. 5B
fig5_dt <- sweep_transition_summary %>%
  filter(
    S > 0,
    is.finite(abs_delta_B),
    is.finite(B_t2)
  ) %>%
  mutate(
    resident_preserving =
      abs_delta_B <= 0.10 &
      B_t2 <= 0.20
  )

# ------------------------------------------------------------
# Observed effects
# ------------------------------------------------------------

transition_obs <- fig5_dt %>%
  summarise(
    n = n(),
    n_resident = sum(resident_preserving),
    estimate = mean(resident_preserving)
  )

population_obs <- fig5_dt %>%
  group_by(population_id) %>%
  summarise(
    resident_preserving =
      any(resident_preserving),
    .groups = "drop"
  ) %>%
  summarise(
    n = n(),
    n_resident = sum(resident_preserving),
    estimate = mean(resident_preserving)
  )






gene_pop
gene_flags
cluster_recurrence
zeroD_clustered
callable_dt
testable_ids
cluster_label_lookup









































# ============================================================
# FINAL FIGURE 6 ANALYSIS
# ============================================================

# ------------------------------------------------------------
# 1. Observed 0D events
# ------------------------------------------------------------

obs_events_final <- unique(
  zeroD_clustered[
    !is.na(Month_t1) &
      !is.na(Month_t2),
    .(
      gene_key,
      cluster_id,
      MAG_id,
      colony_id,
      Genus,
      Month_t1,
      Month_t2
    )
  ]
)

# Number of unique swept genes to preserve within each
# MAG x colony x exact timepair.
timepair_burden_final <- obs_events_final[
  ,
  .(
    k = uniqueN(gene_key)
  ),
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2
  )
]


# ------------------------------------------------------------
# 2. Callable pools
# ------------------------------------------------------------

callable_final <- copy(
  callable_dt[
    callable_0D == TRUE
  ]
)

callable_final[
  ,
  stratum_id := paste(
    MAG_id,
    colony_id,
    sep = "|"
  )
]

timepair_burden_final[
  ,
  stratum_id := paste(
    MAG_id,
    colony_id,
    sep = "|"
  )
]

callable_pools_final <- split(
  callable_final,
  callable_final$stratum_id
)

get_concurrent_families <- function(events) {
  
  x <- unique(
    events[
      cluster_id %chin% testable_ids &
        !is.na(Genus) &
        !is.na(Month_t1) &
        !is.na(Month_t2),
      .(
        cluster_id,
        MAG_id,
        colony_id,
        Genus,
        Month_t1,
        Month_t2
      )
    ]
  )
  
  support <- x[
    ,
    .(
      n_MAGs = uniqueN(MAG_id),
      n_genera = uniqueN(Genus),
      n_colonies = uniqueN(colony_id)
    ),
    by = .(
      cluster_id,
      Month_t1,
      Month_t2
    )
  ]
  
  qualifying <- support[
    n_MAGs >= 2 &
      n_genera >= 2 &
      n_colonies >= 2
  ]
  
  list(
    support = support,
    qualifying = qualifying,
    ids = unique(qualifying$cluster_id)
  )
}


obs_concurrent <- get_concurrent_families(
  obs_events_final
)

obs_concurrent_ids <- obs_concurrent$ids

length(obs_concurrent_ids)





# ============================================================
# ORTHOGONAL pN/pS SUPPORT
# ============================================================

pnps_independent <- gene_pop[
  any_pnps_gt1 == TRUE &
    has_0D_sweep == FALSE,
  .(
    n_independent_pnps_genes =
      uniqueN(gene_key),
    
    n_independent_pnps_MAGs =
      uniqueN(MAG_id),
    
    n_independent_pnps_colonies =
      uniqueN(colony_id),
    
    n_independent_pnps_genera =
      uniqueN(
        Genus[!is.na(Genus)]
      )
  ),
  by = cluster_id
]

pnps_independent <- merge(
  data.table(
    cluster_id = testable_ids
  ),
  pnps_independent,
  by = "cluster_id",
  all.x = TRUE
)

for (
  cc in setdiff(
    names(pnps_independent),
    "cluster_id"
  )
) {
  
  set(
    pnps_independent,
    which(is.na(
      pnps_independent[[cc]]
    )),
    cc,
    0L
  )
}

# Strong independent support:
# elevated pN/pS in non-sweeping populations
# from >=2 genera
pnps_independent[
  ,
  independent_pnps_support :=
    n_independent_pnps_genera >= 2
]

independent_pnps_ids <-
  pnps_independent[
    independent_pnps_support == TRUE,
    cluster_id
  ]

obs_n_concurrent <-
  length(obs_concurrent_ids)

obs_n_pnps_supported <-
  sum(
    obs_concurrent_ids %chin%
      independent_pnps_ids
  )

obs_pnps_fraction <-
  obs_n_pnps_supported /
  obs_n_concurrent

c(
  n_concurrent = obs_n_concurrent,
  n_independent_pnps =
    obs_n_pnps_supported,
  fraction =
    obs_pnps_fraction
)














# ============================================================
# INDEPENDENT CONVERGENCE PAIR SCORE
# ============================================================

independent_pair_score <- function(events) {
  
  x <- unique(
    events[
      cluster_id %chin% testable_ids &
        !is.na(Genus) &
        !is.na(Month_t1) &
        !is.na(Month_t2),
      .(
        cluster_id,
        MAG_id,
        colony_id,
        Genus,
        Month_t1,
        Month_t2
      )
    ]
  )
  
  x[
    ,
    population_id := paste(
      MAG_id,
      colony_id,
      sep = "::"
    )
  ]
  
  x <- unique(
    x[
      ,
      .(
        cluster_id,
        population_id,
        colony_id,
        Genus,
        Month_t1,
        Month_t2
      )
    ]
  )
  
  
  # all population pairs
  total <- x[
    ,
    .(
      all_pairs =
        choose(.N, 2)
    ),
    by = .(
      cluster_id,
      Month_t1,
      Month_t2
    )
  ]
  
  
  # pairs from same genus
  same_genus <- x[
    ,
    .N,
    by = .(
      cluster_id,
      Month_t1,
      Month_t2,
      Genus
    )
  ][
    ,
    .(
      same_genus_pairs =
        sum(choose(N, 2))
    ),
    by = .(
      cluster_id,
      Month_t1,
      Month_t2
    )
  ]
  
  
  # pairs from same colony
  same_colony <- x[
    ,
    .N,
    by = .(
      cluster_id,
      Month_t1,
      Month_t2,
      colony_id
    )
  ][
    ,
    .(
      same_colony_pairs =
        sum(choose(N, 2))
    ),
    by = .(
      cluster_id,
      Month_t1,
      Month_t2
    )
  ]
  
  
  # added back because subtracted twice
  same_both <- x[
    ,
    .N,
    by = .(
      cluster_id,
      Month_t1,
      Month_t2,
      Genus,
      colony_id
    )
  ][
    ,
    .(
      same_both_pairs =
        sum(choose(N, 2))
    ),
    by = .(
      cluster_id,
      Month_t1,
      Month_t2
    )
  ]
  
  
  z <- Reduce(
    function(a, b) {
      merge(
        a,
        b,
        by = c(
          "cluster_id",
          "Month_t1",
          "Month_t2"
        ),
        all = TRUE
      )
    },
    list(
      total,
      same_genus,
      same_colony,
      same_both
    )
  )
  
  for (
    cc in c(
      "same_genus_pairs",
      "same_colony_pairs",
      "same_both_pairs"
    )
  ) {
    
    z[
      is.na(get(cc)),
      (cc) := 0
    ]
  }
  
  z[
    ,
    independent_pairs :=
      all_pairs -
      same_genus_pairs -
      same_colony_pairs +
      same_both_pairs
  ]
  
  sum(z$independent_pairs)
}


obs_pair_score <- independent_pair_score(
  obs_events_final
)

obs_pair_score






obs_family_depth <- obs_concurrent$qualifying[
  ,
  .(
    max_concurrent_MAGs =
      max(n_MAGs),
    
    max_concurrent_genera =
      max(n_genera),
    
    n_concurrent_intervals =
      .N
  ),
  by = cluster_id
]

max_depth <- max(
  obs_family_depth$max_concurrent_MAGs
)

depth_thresholds <- seq(
  2L,
  max_depth
)

obs_depth_counts <- sapply(
  depth_thresholds,
  function(k) {
    
    sum(
      obs_family_depth$
        max_concurrent_MAGs >= k
    )
  }
)

obs_depth_counts


MIN_FUNCTION_BACKGROUND <- 5L

function_levels <- function_background_summary[
  n_testable_families >=
    MIN_FUNCTION_BACKGROUND &
    consensus_function !=
    "Uncharacterized protein",
  consensus_function
]

length(function_levels)



# ============================================================
# FUNCTION LOOKUP
# ============================================================

function_meta <- unique(
  cluster_label_lookup[
    ,
    .(
      cluster_id,
      consensus_function
    )
  ]
)

function_meta[
  is.na(consensus_function) |
    consensus_function == "",
  consensus_function :=
    "Uncharacterized protein"
]


function_background <- merge(
  data.table(
    cluster_id = testable_ids
  ),
  function_meta,
  by = "cluster_id",
  all.x = TRUE
)

function_background[
  is.na(consensus_function),
  consensus_function :=
    "Uncharacterized protein"
]

function_background_summary <-
  function_background[
    ,
    .(
      n_testable_families = .N
    ),
    by = consensus_function
  ][
    order(-n_testable_families)
  ]

function_background_summary





# ============================================================
# FINAL FIGURE 6 PERMUTATION ENGINE
# ============================================================

N_PERM_FIG6 <- 1000L
# FINAL:
# N_PERM_FIG6 <- 10000L

set.seed(12345)

null_n_families <-
  integer(N_PERM_FIG6)

null_pair_score <-
  numeric(N_PERM_FIG6)

null_pnps_count <-
  integer(N_PERM_FIG6)

null_pnps_fraction <-
  rep(NA_real_, N_PERM_FIG6)

null_depth_counts <- matrix(
  0L,
  nrow = N_PERM_FIG6,
  ncol = length(depth_thresholds)
)

colnames(null_depth_counts) <-
  paste0("MAG", depth_thresholds)


null_function_count <- matrix(
  0L,
  nrow = N_PERM_FIG6,
  ncol = length(function_levels)
)

colnames(null_function_count) <-
  function_levels


for (
  b in seq_len(N_PERM_FIG6)
) {
  
  perm_list <- vector(
    "list",
    nrow(timepair_burden_final)
  )
  
  for (
    i in seq_len(
      nrow(timepair_burden_final)
    )
  ) {
    
    pool <- callable_pools_final[[
      timepair_burden_final$
        stratum_id[i]
    ]]
    
    k_i <-
      timepair_burden_final$k[i]
    
    if (
      is.null(pool) ||
      nrow(pool) < k_i
    ) {
      
      stop(
        paste(
          "Insufficient callable genes:",
          timepair_burden_final$
            stratum_id[i]
        )
      )
    }
    
    picked <- pool[
      sample(
        seq_len(.N),
        size = k_i,
        replace = FALSE,
        prob =
          effective_length_0D_sites
      ),
      .(
        cluster_id,
        MAG_id,
        colony_id,
        Genus
      )
    ]
    
    picked[
      ,
      Month_t1 :=
        timepair_burden_final$
        Month_t1[i]
    ]
    
    picked[
      ,
      Month_t2 :=
        timepair_burden_final$
        Month_t2[i]
    ]
    
    perm_list[[i]] <- picked
  }
  
  
  perm_events <- rbindlist(
    perm_list,
    use.names = TRUE
  )
  
  
  # ----------------------------------------
  # Concurrent families
  # ----------------------------------------
  
  conc <- get_concurrent_families(
    perm_events
  )
  
  ids <- conc$ids
  
  null_n_families[b] <-
    length(ids)
  
  
  # ----------------------------------------
  # Independent convergence pairs
  # ----------------------------------------
  
  null_pair_score[b] <-
    independent_pair_score(
      perm_events
    )
  
  
  # ----------------------------------------
  # independent pN/pS support
  # ----------------------------------------
  
  null_pnps_count[b] <-
    sum(
      ids %chin%
        independent_pnps_ids
    )
  
  if (length(ids) > 0) {
    
    null_pnps_fraction[b] <-
      null_pnps_count[b] /
      length(ids)
  }
  
  
  # ----------------------------------------
  # recurrence depth
  # ----------------------------------------
  
  if (
    nrow(conc$qualifying) > 0
  ) {
    
    fam_depth <- conc$qualifying[
      ,
      .(
        max_MAGs = max(n_MAGs)
      ),
      by = cluster_id
    ]
    
    null_depth_counts[b, ] <-
      sapply(
        depth_thresholds,
        function(k) {
          sum(
            fam_depth$max_MAGs >= k
          )
        }
      )
  }
  
  
  # ----------------------------------------
  # function composition
  # ----------------------------------------
  
  if (length(ids) > 0) {
    
    this_functions <- function_meta[
      cluster_id %chin% ids
    ]
    
    this_functions[
      is.na(consensus_function),
      consensus_function :=
        "Uncharacterized protein"
    ]
    
    ff <- this_functions[
      ,
      .N,
      by = consensus_function
    ]
    
    jj <- match(
      ff$consensus_function,
      function_levels
    )
    
    keep <- !is.na(jj)
    
    null_function_count[
      b,
      jj[keep]
    ] <- ff$N[keep]
  }
  
  
  if (b %% 100 == 0) {
    
    message(
      b,
      " / ",
      N_PERM_FIG6
    )
  }
}




# ============================================================
# FIGURE 6B
# Independent pN/pS concordance
# ============================================================

null_pnps_valid <-
  null_pnps_fraction[
    is.finite(
      null_pnps_fraction
    )
  ]

pnps_null_median <-
  median(null_pnps_valid)

pnps_null_ci <-
  quantile(
    null_pnps_valid,
    c(0.025, 0.975)
  )

p_pnps <- (
  sum(
    null_pnps_valid >=
      obs_pnps_fraction
  ) + 1
) /
  (
    length(null_pnps_valid) + 1
  )

c(
  observed = obs_pnps_fraction,
  null_median = pnps_null_median,
  null_lo = pnps_null_ci[1],
  null_hi = pnps_null_ci[2],
  P = p_pnps
)




fig6B_dt <- data.table(
  y = "Independent elevated pN/pS support",
  
  observed =
    obs_pnps_fraction,
  
  expected =
    pnps_null_median,
  
  lo =
    pnps_null_ci[1],
  
  hi =
    pnps_null_ci[2]
)


p6B <- ggplot(
  fig6B_dt,
  aes(y = y)
) +
  
  geom_segment(
    aes(
      x = lo,
      xend = hi,
      yend = y
    ),
    linewidth = 2,
    color = "grey75"
  ) +
  
  geom_point(
    aes(x = expected),
    shape = 21,
    fill = "white",
    size = 3
  ) +
  
  geom_point(
    aes(x = observed),
    shape = 18,
    size = 4
  ) +
  
  annotate(
    "text",
    x = obs_pnps_fraction,
    y = 1.18,
    label = paste0(
      "Observed ",
      percent(
        obs_pnps_fraction,
        accuracy = 0.1
      ),
      "\nPermutation P = ",
      signif(p_pnps, 2)
    ),
    size = 3.3
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    labels = percent
  ) +
  
  labs(
    x = "Concurrent protein families\nwith independent pN/pS > 1 support",
    y = NULL
  ) +
  
  theme_classic(base_size = 12)

p6B










depth_null_summary <- data.table(
  depth = depth_thresholds,
  
  expected = apply(
    null_depth_counts,
    2,
    median
  ),
  
  lo = apply(
    null_depth_counts,
    2,
    quantile,
    probs = 0.025
  ),
  
  hi = apply(
    null_depth_counts,
    2,
    quantile,
    probs = 0.975
  ),
  
  observed =
    obs_depth_counts
)


p6C_depth <- ggplot(
  depth_null_summary,
  aes(x = depth)
) +
  
  geom_ribbon(
    aes(
      ymin = lo,
      ymax = hi
    ),
    alpha = 0.18
  ) +
  
  geom_line(
    aes(y = expected),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(y = expected),
    shape = 21,
    fill = "white",
    size = 2.5
  ) +
  
  geom_line(
    aes(y = observed),
    linewidth = 1
  ) +
  
  geom_point(
    aes(y = observed),
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = depth_thresholds
  ) +
  
  labs(
    x = "Minimum number of independently\nswept MAG populations",
    y = "Concurrent protein families"
  ) +
  
  theme_classic(base_size = 12)

p6C_depth








pair_null_ci <- quantile(
  null_pair_score,
  c(0.025, 0.975)
)

pair_null_median <-
  median(null_pair_score)

p_pair <- (
  sum(
    null_pair_score >=
      obs_pair_score
  ) + 1
) /
  (length(null_pair_score) + 1)

c(
  observed = obs_pair_score,
  expected = pair_null_median,
  lo = pair_null_ci[1],
  hi = pair_null_ci[2],
  P = p_pair
)





p6C_pairs <- ggplot(
  data.table(
    pair_score =
      null_pair_score
  ),
  aes(x = pair_score)
) +
  
  geom_histogram(
    bins = 40,
    fill = "grey85",
    color = "white"
  ) +
  
  geom_vline(
    xintercept =
      obs_pair_score,
    linewidth = 1.1
  ) +
  
  labs(
    x = paste0(
      "Independent cross-genus, cross-colony population pairs\n",
      "with concurrent 0D sweeps in the same protein family"
    ),
    y = "Permutations"
  ) +
  
  annotate(
    "text",
    x = obs_pair_score,
    y = Inf,
    label = paste0(
      "Observed = ",
      obs_pair_score,
      "\nP = ",
      signif(p_pair, 2)
    ),
    hjust = 1.05,
    vjust = 1.4,
    size = 3.3
  ) +
  
  theme_classic(base_size = 12)

p6C_pairs











# ============================================================
# OBSERVED FUNCTION COMPOSITION
# ============================================================

obs_function <- function_meta[
  cluster_id %chin%
    obs_concurrent_ids
][
  ,
  .N,
  by = consensus_function
]


function_results <- data.table(
  consensus_function =
    function_levels
)

function_results[
  ,
  observed := 0L
]

function_results[
  obs_function,
  on = "consensus_function",
  observed := i.N
]

function_results[
  ,
  observed_fraction :=
    observed /
    obs_n_concurrent
]




null_function_fraction <- matrix(
  NA_real_,
  nrow = N_PERM_FIG6,
  ncol = length(function_levels)
)

colnames(null_function_fraction) <-
  function_levels

for (
  b in seq_len(N_PERM_FIG6)
) {
  
  if (null_n_families[b] > 0) {
    
    null_function_fraction[b, ] <-
      null_function_count[b, ] /
      null_n_families[b]
  }
}





function_results[
  ,
  expected_fraction :=
    colMeans(
      null_function_fraction,
      na.rm = TRUE
    )
]

function_results[
  ,
  null_lo :=
    apply(
      null_function_fraction,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    )
]

function_results[
  ,
  null_hi :=
    apply(
      null_function_fraction,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    )
]


function_results[
  ,
  p_upper :=
    vapply(
      seq_len(.N),
      function(j) {
        
        z <-
          null_function_fraction[, j]
        
        z <- z[
          is.finite(z)
        ]
        
        (
          sum(
            z >=
              function_results$
              observed_fraction[j]
          ) + 1
        ) /
          (length(z) + 1)
      },
      numeric(1)
    )
]


function_results[
  ,
  FDR :=
    p.adjust(
      p_upper,
      method = "BH"
    )
]


# Smoothed fold enrichment for display
function_results[
  ,
  log2_enrichment :=
    log2(
      (
        observed + 0.5
      ) /
        (
          expected_fraction *
            obs_n_concurrent +
            0.5
        )
    )
]


function_results <- merge(
  function_results,
  function_background_summary,
  by = "consensus_function",
  all.x = TRUE
)

function_results[
  order(
    FDR,
    -log2_enrichment
  )
]




function_plot_dt <- function_results[
  observed >= 2
][
  order(log2_enrichment)
]





p6D <- ggplot(
  function_plot_dt,
  aes(
    x = log2_enrichment,
    y = reorder(
      consensus_function,
      log2_enrichment
    )
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_point(
    aes(
      size = observed,
      fill = -log10(FDR)
    ),
    shape = 21,
    stroke = 0.5
  ) +
  
  scale_size_continuous(
    name = "Concurrent\nfamilies"
  ) +
  
  labs(
    x = expression(
      log[2]~
        "enrichment among recurrent families"
    ),
    y = NULL,
    fill = expression(
      -log[10]~FDR
    )
  ) +
  
  theme_classic(base_size = 12)

p6D







# ============================================================
# INDEPENDENT pN/pS FAMILY-ENRICHMENT ANALYSIS
# ============================================================

pnps_universe <- gene_pop[
  n_months_with_pnps > 0 &
    has_0D_sweep == FALSE
]

pnps_universe[
  ,
  pnps_positive :=
    any_pnps_gt1 == TRUE
]

pnps_universe[
  ,
  population_id :=
    paste(
      MAG_id,
      colony_id,
      sep = "::"
    )
]









# ============================================================
# SWEEP-INDEPENDENT pN/pS UNIVERSE
# ============================================================

pnps_universe <- copy(
  gene_pop[
    n_months_with_pnps > 0 &
      has_0D_sweep == FALSE
  ]
)

pnps_universe[
  ,
  pnps_positive :=
    any_pnps_gt1 == TRUE
]

pnps_universe[
  ,
  population_id := paste(
    MAG_id,
    colony_id,
    sep = "::"
  )
]

# Preserve MAG x colony background AND
# opportunity to observe pN/pS > 1
pnps_universe[
  ,
  perm_stratum := paste(
    MAG_id,
    colony_id,
    n_months_with_pnps,
    sep = "|"
  )
]

stopifnot(
  nrow(pnps_universe) > 0,
  length(pnps_universe$perm_stratum) ==
    nrow(pnps_universe),
  !anyNA(pnps_universe$perm_stratum)
)

pnps_strata <- split(
  seq_len(nrow(pnps_universe)),
  pnps_universe$perm_stratum
)

pnps_k <- pnps_universe[
  ,
  .(
    k = sum(pnps_positive),
    n = .N
  ),
  by = perm_stratum
]

k_lookup <- setNames(
  pnps_k$k,
  pnps_k$perm_stratum
)


family_ids <- pnps_testable_ids

obs_vec <- pnps_family_obs[
  match(
    family_ids,
    cluster_id
  ),
  n_positive
]

obs_vec[is.na(obs_vec)] <- 0L





set.seed(12345)

N_PERM_PNPS <- 1000L
# final:
# N_PERM_PNPS <- 10000L

null_sum <- numeric(
  length(family_ids)
)

null_exceed <- integer(
  length(family_ids)
)

# useful for global recurrence test
null_recurrent_2genus <-
  integer(N_PERM_PNPS)

null_recurrent_3pop <-
  integer(N_PERM_PNPS)







for (b in seq_len(N_PERM_PNPS)) {
  
  positive_idx <- integer()
  
  for (s in names(pnps_strata)) {
    
    idx <- pnps_strata[[s]]
    
    k <- k_lookup[[s]]
    
    if (k == 0)
      next
    
    positive_idx <- c(
      positive_idx,
      sample(
        idx,
        size = k,
        replace = FALSE
      )
    )
  }
  
  
  perm_positive <-
    pnps_universe[
      positive_idx,
      .(
        cluster_id,
        MAG_id,
        colony_id,
        Genus
      )
    ]
  
  
  # -----------------------------------------
  # family positive counts
  # -----------------------------------------
  
  cc <- perm_positive[
    cluster_id %chin%
      family_ids,
    .(
      null_positive = .N
    ),
    by = cluster_id
  ]
  
  null_vec <- integer(
    length(family_ids)
  )
  
  jj <- match(
    cc$cluster_id,
    family_ids
  )
  
  null_vec[jj] <-
    cc$null_positive
  
  null_sum <-
    null_sum + null_vec
  
  null_exceed <-
    null_exceed +
    as.integer(
      null_vec >= obs_vec
    )
  
  
  # -----------------------------------------
  # recurrence breadth
  # -----------------------------------------
  
  breadth <- perm_positive[
    cluster_id %chin%
      family_ids,
    .(
      n_positive_populations = .N,
      
      n_positive_genera =
        uniqueN(
          Genus[
            !is.na(Genus)
          ]
        )
    ),
    by = cluster_id
  ]
  
  null_recurrent_2genus[b] <-
    breadth[
      n_positive_genera >= 2,
      .N
    ]
  
  null_recurrent_3pop[b] <-
    breadth[
      n_positive_populations >= 3,
      .N
    ]
  
  
  if (b %% 100 == 0) {
    message(
      b,
      " / ",
      N_PERM_PNPS
    )
  }
}





pnps_family_test <- pnps_family_obs[
  cluster_id %chin%
    family_ids
]

pnps_family_test[
  ,
  expected_positive :=
    null_sum /
    N_PERM_PNPS
]

pnps_family_test[
  ,
  p_empirical :=
    (
      null_exceed + 1
    ) /
    (
      N_PERM_PNPS + 1
    )
]

pnps_family_test[
  ,
  FDR :=
    p.adjust(
      p_empirical,
      method = "BH"
    )
]

pnps_family_test[
  ,
  log2_enrichment :=
    log2(
      (n_positive + 0.5) /
        (expected_positive + 0.5)
    )
]
pnps_family_test <- merge(
  pnps_family_test,
  cluster_label_lookup[
    ,
    .(
      cluster_id,
      consensus_function,
      consensus_fraction
    )
  ],
  by = "cluster_id",
  all.x = TRUE
)

pnps_family_test[
  order(
    FDR,
    -log2_enrichment,
    -n_positive
  )
][1:50]




obs_recurrent_2genus <-
  pnps_family_test[
    n_positive_genera >= 2,
    .N
  ]

obs_recurrent_3pop <-
  pnps_family_test[
    n_positive >= 3,
    .N
  ]

p_pnps_2genus <- (
  sum(
    null_recurrent_2genus >=
      obs_recurrent_2genus
  ) + 1
) /
  (N_PERM_PNPS + 1)

p_pnps_3pop <- (
  sum(
    null_recurrent_3pop >=
      obs_recurrent_3pop
  ) + 1
) /
  (N_PERM_PNPS + 1)

c(
  obs_2genus =
    obs_recurrent_2genus,
  
  null_median_2genus =
    median(
      null_recurrent_2genus
    ),
  
  P_2genus =
    p_pnps_2genus,
  
  obs_3pop =
    obs_recurrent_3pop,
  
  null_median_3pop =
    median(
      null_recurrent_3pop
    ),
  
  P_3pop =
    p_pnps_3pop
)




pnps_enriched_ids <- pnps_family_test[
  FDR < 0.05 &
    log2_enrichment > 0,
  cluster_id
]

sweep_recurrent_ids <-
  obs_concurrent_ids

observed_overlap <- length(
  intersect(
    pnps_enriched_ids,
    sweep_recurrent_ids
  )
)

observed_overlap

















library(data.table)

# ============================================================
# UNIQUE PHYSICAL SWEEP EVENTS
# ============================================================

sweep_sites <- unique(
  sweep_dt[
    !is.na(contig_name) &
      !is.na(pos_in_contig) &
      !is.na(Month_t1) &
      !is.na(Month_t2) &
      degeneracy %chin% c("0D", "4D"),
    .(
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      contig_name,
      pos_in_contig,
      unique_SNV_identifier,
      degeneracy,
      freq_t1,
      freq_t2,
      gene_key,
      corresponding_gene_call
    )
  ]
)

sweep_sites[
  ,
  delta_AF := freq_t2 - freq_t1
]

sweep_sites[
  ,
  abs_delta_AF := abs(delta_AF)
]






c(
  original_rows = nrow(
    sweep_dt[
      degeneracy %chin% c("0D", "4D")
    ]
  ),
  unique_events = nrow(sweep_sites)
)







zeroD_sites <- sweep_sites[
  degeneracy == "0D"
]

fourD_sites <- sweep_sites[
  degeneracy == "4D"
]

c(
  n_0D = nrow(zeroD_sites),
  n_4D = nrow(fourD_sites)
)






sweep_transition_counts <- sweep_sites[
  ,
  .(
    n_0D = sum(degeneracy == "0D"),
    n_4D = sum(degeneracy == "4D"),
    n_total = .N
  ),
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2
  )
]

summary(sweep_transition_counts$n_0D)
summary(sweep_transition_counts$n_4D)

sweep_transition_counts[
  ,
  .(
    transitions = .N,
    with_0D = sum(n_0D > 0),
    with_4D = sum(n_4D > 0),
    with_both = sum(
      n_0D > 0 &
        n_4D > 0
    )
  )
]












# Give every focal nonsynonymous event an ID
zeroD_sites[
  ,
  focal_0D_id := .I
]

zeroD_for_pair <- zeroD_sites[
  ,
  .(
    focal_0D_id,
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_0D = pos_in_contig,
    freq_0D_t1 = freq_t1,
    freq_0D_t2 = freq_t2,
    delta_0D = delta_AF,
    abs_delta_0D = abs_delta_AF,
    gene_0D = gene_key
  )
]

fourD_for_pair <- fourD_sites[
  ,
  .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_4D = pos_in_contig,
    freq_4D_t1 = freq_t1,
    freq_4D_t2 = freq_t2,
    delta_4D = delta_AF,
    abs_delta_4D = abs_delta_AF,
    gene_4D = gene_key,
    SNV_4D = unique_SNV_identifier
  )
]








zero_four_pairs <- merge(
  zeroD_for_pair,
  fourD_for_pair,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2",
    "contig_name"
  ),
  allow.cartesian = TRUE
)

zero_four_pairs[
  ,
  distance_bp :=
    abs(pos_4D - pos_0D)
]





nearest_4D <- zero_four_pairs[
  order(distance_bp),
  .SD[1],
  by = focal_0D_id
]






zeroD_nearest <- merge(
  zeroD_for_pair,
  nearest_4D[
    ,
    .(
      focal_0D_id,
      nearest_4D_bp = distance_bp,
      nearest_4D_pos = pos_4D,
      nearest_4D_delta = delta_4D,
      nearest_4D_abs_delta = abs_delta_4D
    )
  ],
  by = "focal_0D_id",
  all.x = TRUE
)




zeroD_nearest[
  ,
  .(
    n_0D = .N,
    
    prop_same_contig =
      mean(!is.na(nearest_4D_bp)),
    
    prop_within_1kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 1000
      ),
    
    prop_within_5kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 5000
      ),
    
    prop_within_10kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 10000
      ),
    
    prop_within_25kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 25000
      )
  )
]







nearest_plot_dt <- zeroD_nearest[
  !is.na(nearest_4D_bp)
]

ggplot(
  nearest_plot_dt,
  aes(
    x = nearest_4D_bp / 1000
  )
) +
  geom_histogram(
    bins = 40,
    color = "white"
  ) +
  scale_x_log10(
    labels = scales::label_number()
  ) +
  labs(
    x = "Distance to nearest contemporaneous\n4D sweep (kb; log scale)",
    y = "0D sweep events"
  ) +
  theme_classic(base_size = 12)




distance_thresholds <- c(
  100,
  250,
  500,
  1000,
  2500,
  5000,
  10000,
  25000,
  50000,
  100000
)

distance_curve <- data.table(
  distance_bp = distance_thresholds
)

distance_curve[
  ,
  fraction_0D :=
    vapply(
      distance_bp,
      function(d) {
        
        mean(
          !is.na(
            zeroD_nearest$
              nearest_4D_bp
          ) &
            zeroD_nearest$
            nearest_4D_bp <= d
        )
        
      },
      numeric(1)
    )
]




ggplot(
  distance_curve,
  aes(
    x = distance_bp / 1000,
    y = fraction_0D
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_x_log10() +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = "Maximum distance to 4D sweep (kb)",
    y = "0D sweeps with a contemporaneous\n4D sweep within distance"
  ) +
  theme_classic(base_size = 12)











zeroD_nearest[
  ,
  same_direction :=
    sign(delta_0D) ==
    sign(nearest_4D_delta)
]
zeroD_nearest[
  !is.na(nearest_4D_bp),
  .(
    n = .N,
    same_direction_fraction =
      mean(same_direction)
  )
]





ggplot(
  zero_four_pairs,
  aes(
    x = distance_bp / 1000,
    y = abs_delta_4D
  )
) +
  geom_point(
    alpha = 0.08,
    size = 1
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  scale_x_log10() +
  labs(
    x = "Distance from focal 0D sweep (kb)",
    y = expression(
      "|AF| of contemporaneous 4D sweep"
    )
  ) +
  theme_classic(base_size = 12)









REGION_GAP <- 5000




sweep_regions <- sweep_sites[
  order(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_in_contig
  ),
  {
    
    p <- sort(
      unique(pos_in_contig)
    )
    
    if (length(p) == 0) {
      
      .(
        n_sites = 0L,
        n_regions = 0L
      )
      
    } else {
      
      region_id <- cumsum(
        c(
          TRUE,
          diff(p) >
            REGION_GAP
        )
      )
      
      .(
        n_sites = length(p),
        n_regions =
          uniqueN(region_id)
      )
    }
    
  },
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    contig_name
  )
]




transition_regions <- sweep_regions[
  ,
  .(
    n_sweep_sites =
      sum(n_sites),
    
    n_spatial_regions =
      sum(n_regions),
    
    n_affected_contigs =
      sum(n_sites > 0)
  ),
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2
  )
]

summary(
  transition_regions$
    n_spatial_regions
)




# 1
c(
  original_rows = nrow(
    sweep_dt[
      degeneracy %chin%
        c("0D", "4D")
    ]
  ),
  unique_events =
    nrow(sweep_sites),
  n_0D =
    nrow(zeroD_sites),
  n_4D =
    nrow(fourD_sites)
)

table(
  transition_regions$
    n_spatial_regions
)

# 3
zeroD_nearest[
  ,
  .(
    n_0D = .N,
    prop_same_contig =
      mean(!is.na(nearest_4D_bp)),
    prop_within_1kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 1000
      ),
    prop_within_5kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 5000
      ),
    prop_within_10kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 10000
      ),
    prop_within_25kb =
      mean(
        !is.na(nearest_4D_bp) &
          nearest_4D_bp <= 25000
      )
  )
]








