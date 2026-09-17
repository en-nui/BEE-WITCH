#!/usr/bin/env Rscript

# ============================================================
# Protein-family representative screen for Figure 6E / future
# AlphaFold visualization
#
# Purpose
# -------
# Starting from the finalized MMseqs protein-family universe and
# observed 0D Delta RAF events, identify protein families that are:
#   1) strongly and independently recurrent,
#   2) temporally concurrent across genera/colonies,
#   3) supported by many observed 0D events, and
#   4) sufficiently similar in protein sequence to support a common
#      representative structure.
#
# Workflow A-F
# ------------
# A. Load and QC the authoritative physical-gene / cluster universe
# B. Reconstruct strict recurrent/concurrent family candidates
# C. Extract candidate member sequences and event support
# D. Align top candidate families with MAFFT and calculate within-family
#    identity / coverage from the MSA
# E. Select a representative (MSA medoid) and summarize conservation
# F. Export ranked candidates, representative FASTAs, MSAs, and event
#    inventories for downstream AlphaFold / structural mapping
#
# IMPORTANT
# ---------
# This is an exploratory candidate-selection workflow. Sequence similarity
# is NOT used to define the recurrence test or alter Figure 6 inferential
# statistics. It is used only after the recurrence analysis to identify a
# biologically interpretable representative family for visualization.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(ape)
  library(ips)
})

# ------------------------------------------------------------
# 0. Paths and controls
# ------------------------------------------------------------

ROOT <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel/parallel_with_CLUSTERED_PROTEINS/"

CLUSTER_FILE <- file.path(
  ROOT,
  "selection_gene_clusters",
  "eligible_gene_mmseqs_clusters.csv.gz"
)

SWEEP_FILE <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/sweep_gene_merged.csv"

ECOLOGY_FILE <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"

ALL_FASTA <- file.path(
  ROOT,
  "selection_gene_clusters",
  "eligible_bacterial_mmseqs_all_seqs.fasta"
)

REP_FASTA <- file.path(
  ROOT,
  "selection_gene_clusters",
  "eligible_bacterial_mmseqs_rep_seq.fasta"
)

PERM_FILE <- file.path(
  ROOT,
  "selection_gene_clusters",
  "selection_flags",
  "MMseqs_0D_sweep_permutation_all_clusters_10000.csv.gz"
)

OUTDIR <- file.path(
  ROOT,
  "selection_gene_clusters",
  "protein_representative_screen"
)
SNV_DIR <- paste0(
  "/home/robinch/projects/BEE-WITCH/",
  "40_internal_pnps_calculations_from_anvio/",
  "colony_outputs"
)

snv_files <- list.files(
  SNV_DIR,
  pattern = "^colony_[0-9]+_snps_with_degeneracy\\.csv$",
  full.names = TRUE
)

length(snv_files)
snv_files[1:5]


message("Loading codon-position annotations...")

snv_codon_map <- rbindlist(
  lapply(
    snv_files,
    function(f) {
      
      x <- fread(
        f,
        select = c(
          "Colony",
          "MAG_id",
          "contig_name",
          "pos_in_contig",
          "corresponding_gene_call",
          "codon_order_in_gene",
          "base_pos_in_codon",
          "unique_pos_identifier"
        ),
        showProgress = FALSE
      )
      
      x[
        ,
        .(
          colony_id =
            as.character(Colony),
          
          MAG_id =
            as.character(MAG_id),
          
          contig_name =
            as.character(contig_name),
          
          pos_in_contig =
            as.integer(pos_in_contig),
          
          corresponding_gene_call =
            as.character(corresponding_gene_call),
          
          codon_order_in_gene =
            as.integer(codon_order_in_gene),
          
          base_pos_in_codon =
            as.integer(base_pos_in_codon),
          
          unique_pos_identifier =
            as.character(unique_pos_identifier)
        )
      ]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

message(
  "Loaded ",
  format(nrow(snv_codon_map), big.mark = ","),
  " SNV annotation rows"
)

snv_codon_map[
  ,
  .(
    n_rows = .N,
    n_gene_calls =
      uniqueN(corresponding_gene_call),
    n_codon_orders =
      uniqueN(codon_order_in_gene),
    n_base_positions =
      uniqueN(base_pos_in_codon)
  ),
  by = .(
    colony_id,
    MAG_id,
    contig_name,
    pos_in_contig
  )
][
  n_gene_calls > 1 |
    n_codon_orders > 1 |
    n_base_positions > 1
][
  1:20
]


snv_codon_map_qc <- snv_codon_map[
  !is.na(colony_id) &
    !is.na(MAG_id) &
    !is.na(contig_name) &
    !is.na(pos_in_contig)
]


dup_codon_sites <- snv_codon_map_qc[
  ,
  .(
    n_rows = .N,
    n_gene_calls = uniqueN(corresponding_gene_call),
    n_codon_orders = uniqueN(codon_order_in_gene),
    n_base_positions = uniqueN(base_pos_in_codon)
  ),
  by = .(
    colony_id,
    MAG_id,
    contig_name,
    pos_in_contig
  )
][
  n_gene_calls > 1 |
    n_codon_orders > 1 |
    n_base_positions > 1
]

nrow(dup_codon_sites)


site_multiplicity <- snv_codon_map_qc[
  ,
  .(
    n_rows = .N
  ),
  by = .(
    colony_id,
    MAG_id,
    contig_name,
    pos_in_contig
  )
]

site_multiplicity[
  ,
  .(
    n_sites = .N,
    singletons = sum(n_rows == 1),
    duplicated = sum(n_rows > 1),
    max_rows = max(n_rows)
  )
]


snv_codon_map_unique <- unique(
  snv_codon_map_qc[
    ,
    .(
      colony_id,
      MAG_id,
      contig_name,
      pos_in_contig,
      corresponding_gene_call,
      codon_order_in_gene,
      base_pos_in_codon
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_in_contig"
  )
)






snv_codon_map_unique <- unique(
  snv_codon_map,
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_in_contig"
  )
)


snv_codon_map_unique[
  ,
  .N,
  by = .(
    colony_id,
    MAG_id,
    contig_name,
    pos_in_contig
  )
][
  N > 1
]
sweep_0D_unique <- unique(
  sweep_dt[
    degeneracy == "0D",
    .(
      colony_id,
      MAG_id,
      contig_name,
      pos_in_contig,
      unique_SNV_identifier,
      gene_key,
      corresponding_gene_call
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_in_contig"
  )
)

nrow(sweep_0D_unique)


sweep_0D_unique[
  ,
  .(
    n_sites = .N,
    n_with_codon =
      sum(
        snv_codon_map_unique[
          ,
          .(key = paste(
            colony_id,
            MAG_id,
            contig_name,
            pos_in_contig,
            sep = "|"
          ))
        ]$key %chin%
          paste(
            colony_id,
            MAG_id,
            contig_name,
            pos_in_contig,
            sep = "|"
          )
      )
  )
]


sweep_0D_unique[
  ,
  site_key := paste(
    colony_id,
    MAG_id,
    contig_name,
    pos_in_contig,
    sep = "|"
  )
]

snv_codon_map_unique[
  ,
  site_key := paste(
    colony_id,
    MAG_id,
    contig_name,
    pos_in_contig,
    sep = "|"
  )
]

sweep_0D_unique[
  ,
  has_codon_annotation :=
    site_key %chin% snv_codon_map_unique$site_key
]

sweep_0D_unique[
  ,
  .(
    n_sites = .N,
    n_matched = sum(has_codon_annotation),
    prop_matched = mean(has_codon_annotation)
  )
]





DIR_ALIGNMENT <- file.path(OUTDIR, "alignments")
DIR_MEMBER_FASTA <- file.path(OUTDIR, "member_fastas")
DIR_REP_FASTA <- file.path(OUTDIR, "representatives")
DIR_EVENTS <- file.path(OUTDIR, "event_inventories")
DIR_MSA_TABLES <- file.path(OUTDIR, "msa_tables")

for (d in c(
  OUTDIR,
  DIR_ALIGNMENT,
  DIR_MEMBER_FASTA,
  DIR_REP_FASTA,
  DIR_EVENTS,
  DIR_MSA_TABLES
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Number of families to take forward for MSA.
# Start small; increase after inspecting the first run.
N_ALIGN_FAMILIES <- 30L

# Minimum support for the structural-screen candidate pool.
MIN_CONCURRENT_GENERA <- 2L
MIN_CONCURRENT_COLONIES <- 2L

# Sequence-quality criteria used ONLY for ranking / interpretation.
# These do not redefine recurrence.
MIN_MEDIAN_IDENTITY_FOR_PREFERRED <- 0.70
MIN_MEDIAN_COVERAGE_FOR_PREFERRED <- 0.80

# Optional use of final permutation significance as a flag.
USE_PERMUTATION_FLAG <- TRUE

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

clean_colony <- function(x) {
  x <- as.character(x)
  x <- sub("^colony_", "", x)
  x <- sub("\\.0$", "", x)
  x
}

safe_num <- function(x) {
  as.numeric(as.character(x))
}

write_fasta <- function(seqs, filepath) {
  if (length(seqs) == 0L) {
    stop("No sequences supplied for FASTA: ", filepath)
  }
  Biostrings::writeXStringSet(seqs, filepath)
}

# Return a symmetric pairwise identity matrix and a symmetric
# pairwise alignment-coverage matrix from an MSA.
# Identity is calculated only over alignment columns in which both
# sequences contain residues; coverage is the fraction of the shorter
# ungapped sequence represented in the pairwise aligned overlap.
msa_pairwise_metrics <- function(aln) {
  mat <- as.matrix(aln)
  ids <- rownames(mat)
  if (is.null(ids)) {
    ids <- names(aln)
  }
  
  n <- nrow(mat)
  ungapped_len <- rowSums(mat != "-")
  
  identity <- matrix(
    NA_real_,
    nrow = n,
    ncol = n,
    dimnames = list(ids, ids)
  )
  
  coverage <- matrix(
    NA_real_,
    nrow = n,
    ncol = n,
    dimnames = list(ids, ids)
  )
  
  diag(identity) <- 1
  diag(coverage) <- 1
  
  if (n <= 1L) {
    return(list(
      identity = identity,
      coverage = coverage
    ))
  }
  
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      both <- mat[i, ] != "-" & mat[j, ] != "-"
      overlap <- sum(both)
      
      if (overlap == 0L) {
        next
      }
      
      pid <- sum(mat[i, both] == mat[j, both]) / overlap
      cov <- overlap / min(ungapped_len[i], ungapped_len[j])
      
      identity[i, j] <- pid
      identity[j, i] <- pid
      coverage[i, j] <- cov
      coverage[j, i] <- cov
    }
  }
  
  list(
    identity = identity,
    coverage = coverage
  )
}

summarize_msa <- function(aln, cluster_id) {
  ids <- names(aln)
  m <- msa_pairwise_metrics(aln)
  
  sim_long <- data.table()
  if (length(ids) > 1L) {
    pairs <- which(upper.tri(m$identity), arr.ind = TRUE)
    if (nrow(pairs) > 0L) {
      sim_long <- data.table(
        cluster_id = cluster_id,
        seq_i = ids[pairs[, 1]],
        seq_j = ids[pairs[, 2]],
        pair_identity = m$identity[pairs],
        pair_coverage = m$coverage[pairs]
      )
    }
  }
  
  # Per-sequence medoid statistics. We use pairs with reasonably strong
  # overlap for representative selection, but retain unrestricted summaries.
  per_seq <- rbindlist(
    lapply(seq_along(ids), function(i) {
      pid <- m$identity[i, ]
      cov <- m$coverage[i, ]
      keep <- seq_along(ids) != i & is.finite(pid) & is.finite(cov)
      strong <- keep & cov >= MIN_MEDIAN_COVERAGE_FOR_PREFERRED
      
      data.table(
        cluster_id = cluster_id,
        sequence_id = ids[i],
        median_identity = if (any(keep)) median(pid[keep], na.rm = TRUE) else NA_real_,
        mean_identity = if (any(keep)) mean(pid[keep], na.rm = TRUE) else NA_real_,
        min_identity = if (any(keep)) min(pid[keep], na.rm = TRUE) else NA_real_,
        median_coverage = if (any(keep)) median(cov[keep], na.rm = TRUE) else NA_real_,
        min_coverage = if (any(keep)) min(cov[keep], na.rm = TRUE) else NA_real_,
        n_members_compared = sum(keep),
        n_strong_pairs = sum(strong)
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  list(
    pairwise = sim_long,
    per_sequence = per_seq
  )
}

# Given a candidate MSA and sequence metadata, choose a medoid based on
# median identity among pairs with sufficient overlap. Tie-break on median
# coverage, then MMseqs representative status, then total evidence later.
choose_medoid <- function(seq_summary, member_meta) {
  x <- merge(
    seq_summary,
    member_meta[
      ,
      .(
        sequence_id,
        is_mmseq_rep,
        n_physical_genes,
        n_MAGs
      )
    ],
    by = "sequence_id",
    all.x = TRUE
  )
  
  x[
    ,
    preferred :=
      is.finite(median_identity) &
      is.finite(median_coverage) &
      median_identity >= MIN_MEDIAN_IDENTITY_FOR_PREFERRED &
      median_coverage >= MIN_MEDIAN_COVERAGE_FOR_PREFERRED
  ]
  
  setorderv(
    x,
    c(
      "preferred",
      "median_identity",
      "median_coverage",
      "is_mmseq_rep",
      "n_physical_genes",
      "n_MAGs"
    ),
    c(-1, -1, -1, -1, -1, -1),
    na.last = TRUE
  )
  
  x[1L]
}

# Parse a minimal FASTA string robustly when headers contain spaces.
read_named_aa_fasta <- function(path) {
  x <- Biostrings::readAAStringSet(path)
  nm <- names(x)
  nm <- sub("\\s+.*$", "", nm)
  names(x) <- nm
  x
}

message("============================================================")
message("Protein-family representative screen")
message("============================================================")

# ------------------------------------------------------------
# A. Load authoritative physical-gene universe and event data
# ------------------------------------------------------------

stopifnot(
  file.exists(CLUSTER_FILE),
  file.exists(SWEEP_FILE),
  file.exists(ALL_FASTA),
  file.exists(REP_FASTA)
)
ecology_dt <- fread(ECOLOGY_FILE)
cluster_dt <- fread(CLUSTER_FILE)
sweep_dt <- fread(SWEEP_FILE)

cluster_dt[, MAG_id := as.character(MAG_id)]
cluster_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]
cluster_dt[, cluster_id := as.character(cluster_id)]
cluster_dt[, sequence_id := as.character(sequence_id)]
cluster_dt[, cluster_rep := as.character(cluster_rep)]
cluster_dt[, gene_key := paste(MAG_id, corresponding_gene_call, sep = "::")]

mag_genus_map <- unique(ecology_dt[!is.na(Genus), .(MAG_id, Genus)])
cluster_dt[mag_genus_map, genus := i.Genus, on = "MAG_id"]
stopifnot(uniqueN(cluster_dt$gene_key) == nrow(cluster_dt))
stopifnot(!anyDuplicated(cluster_dt$gene_key))

cluster_dt
if (!"genus" %in% names(cluster_dt)) {
  stop("Cluster table lacks Genus. Attach canonical taxonomy before running.")
}

sweep_dt[, MAG_id := as.character(MAG_id)]
sweep_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]

if ("colony_id" %in% names(sweep_dt)) {
  sweep_dt[, colony_id := clean_colony(colony_id)]
} else if ("Colony" %in% names(sweep_dt)) {
  sweep_dt[, colony_id := clean_colony(Colony)]
} else {
  stop("Sweep table contains neither colony_id nor Colony.")
}

sweep_dt[, gene_key := paste(MAG_id, corresponding_gene_call, sep = "::")]

zeroD <- sweep_dt[
  degeneracy == "0D" & gene_key %chin% cluster_dt$gene_key
]

# Collapse row-level duplicates to one physical-gene x population x interval event.
if (!all(c("Month_t1", "Month_t2") %in% names(zeroD))) {
  stop("Expected Month_t1 and Month_t2 in sweep table for temporal concurrence.")
}

zeroD_events <- unique(
  zeroD[
    !is.na(Month_t1) & !is.na(Month_t2),
    .(
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2
    )
  ]
)

zeroD_events <- merge(
  zeroD_events,
  unique(
    cluster_dt[
      ,
      .(
        gene_key,
        cluster_id,
        sequence_id,
        cluster_rep,
        genus,
        is_cluster_rep
      )
    ]
  ),
  by = "gene_key",
  all.x = TRUE
)

stopifnot(!any(is.na(zeroD_events$cluster_id)))

zeroD[mag_genus_map, genus := i.Genus, on = "MAG_id"]


zeroD
# ------------------------------------------------------------
# B. Reconstruct strict recurrent / temporally concurrent families
# ------------------------------------------------------------

# Family eligibility for the primary recurrence universe.
cluster_background <- cluster_dt[
  ,
  .(
    n_genes = uniqueN(gene_key),
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(genus, na.rm = TRUE),
    n_colonies = uniqueN(MAG_id)
  ),
  by = cluster_id
]
cluster_dt
zeroD_events
# Exact family x interval support.
family_time_support <- zeroD_events[
  ,
  .(
    n_0D_genes = uniqueN(gene_key),
    n_MAGs = uniqueN(MAG_id),
    n_genera = uniqueN(genus, na.rm = TRUE),
    n_colonies = uniqueN(colony_id),
    n_event_records = .N
  ),
  by = .(
    cluster_id,
    Month_t1,
    Month_t2
  )
]

qualifying_intervals <- family_time_support[
  n_MAGs >= 2 &
    n_genera >= MIN_CONCURRENT_GENERA &
    n_colonies >= MIN_CONCURRENT_COLONIES
]

strict_ids <- unique(qualifying_intervals$cluster_id)

message("Strictly concurrent protein families: ", length(strict_ids))

# Whole-family recurrence summary.
cluster_recurrence_screen <- zeroD_events[
  ,
  .(
    total_0D_event_records = .N,
    n_0D_genes = uniqueN(gene_key),
    n_0D_MAGs = uniqueN(MAG_id),
    n_0D_genera = uniqueN(genus, na.rm = TRUE),
    n_0D_colonies = uniqueN(colony_id),
    n_timepairs_with_0D = uniqueN(
      paste(Month_t1, Month_t2, sep = "->")
    )
  ),
  by = cluster_id
]

# Strongest exact-time support per family.
concurrent_family_summary <- qualifying_intervals[
  ,
  .(
    n_concurrent_intervals = .N,
    max_concurrent_MAGs = max(n_MAGs),
    max_concurrent_genera = max(n_genera),
    max_concurrent_colonies = max(n_colonies),
    total_concurrent_gene_records = sum(n_event_records),
    dominant_concurrent_interval = {
      ord <- order(
        -n_genera,
        -n_MAGs,
        -n_colonies,
        Month_t1,
        Month_t2
      )
      paste(
        Month_t1[ord[1L]],
        Month_t2[ord[1L]],
        sep = " -> "
      )
    }
  ),
  by = cluster_id
]

candidate_family_summary <- merge(
  cluster_background,
  cluster_recurrence_screen,
  by = "cluster_id",
  all.x = TRUE
)

candidate_family_summary <- merge(
  candidate_family_summary,
  concurrent_family_summary,
  by = "cluster_id",
  all.x = TRUE
)

candidate_family_summary[
  is.na(n_0D_event_records),
  `:=`(
    n_0D_event_records = 0L,
    n_0D_genes = 0L,
    n_0D_MAGs = 0L,
    n_0D_genera = 0L,
    n_0D_colonies = 0L,
    n_timepairs_with_0D = 0L
  )
]
candidate_family_summary

candidate_family_summary[
  is.na(n_concurrent_intervals),
  `:=`(
    n_concurrent_intervals = 0L,
    max_concurrent_MAGs = 0L,
    max_concurrent_genera = 0L,
    max_concurrent_colonies = 0L,
    total_concurrent_gene_records = 0L,
    dominant_concurrent_interval = NA_character_
  )
]

candidate_family_summary[
  ,
  strict_concurrent := cluster_id %chin% strict_ids
]

# Keep only the families that satisfy the strict concurrence definition.
candidate_family_summary <- candidate_family_summary[
  strict_concurrent == TRUE
]

# Add optional final permutation significance flag if the finalized table exists.
if (USE_PERMUTATION_FLAG && file.exists(PERM_FILE)) {
  perm_results <- fread(PERM_FILE)
  perm_results[, cluster_id := as.character(cluster_id)]
  
  perm_meta <- unique(
    perm_results[
      ,
      .(
        cluster_id,
        p_empirical = if ("p_empirical_MAG" %in% names(perm_results)) p_empirical_MAG else NA_real_,
        FDR = if ("FDR_MAG_primary" %in% names(perm_results)) FDR_MAG_primary else NA_real_
      )
    ]
  )
  
  candidate_family_summary <- merge(
    candidate_family_summary,
    perm_meta,
    by = "cluster_id",
    all.x = TRUE
  )
  
  candidate_family_summary[
    ,
    FDR05 := is.finite(FDR) & FDR < 0.05
  ]
} else {
  candidate_family_summary[
    ,
    `:=`(
      p_empirical = NA_real_,
      FDR = NA_real_,
      FDR05 = FALSE
    )
  ]
}

# Independent pN/pS support, intentionally NOT including genes that already
# contain a 0D event. This makes it genuinely orthogonal to sweep recurrence.
gene_level_flags <- unique(
  cluster_dt[
    ,
    .(
      gene_key,
      cluster_id,
      MAG_id,
      genus
    )
  ]
)

# Load only the pN/pS information needed from the gene summary.
GENE_SUM_FILE <- file.path(
  ROOT,
  "master_per_gene_summary.csv"
)

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
gene_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]
gene_dt[, gene_key := paste(MAG_id, corresponding_gene_call, sep = "::")]

gene_dt <- gene_dt[gene_key %chin% cluster_dt$gene_key]

gene_dt[
  ,
  pnps_gt1 := !is.na(pN_pS_ratio) & pN_pS_ratio > 1
]

# Population-level genes with any pN/pS > 1.
pnps_pop <- gene_dt[
  ,
  .(
    any_pnps_gt1 = any(pnps_gt1)
  ),
  by = .(
    gene_key,
    MAG_id,
    colony_id,
    corresponding_gene_call
  )
]

# Mark whether a gene-population carries any 0D event.
swept_gene_pops <- unique(
  zeroD_events[
    ,
    .(
      gene_key,
      MAG_id,
      colony_id
    )
  ]
)

pnps_pop <- merge(
  pnps_pop,
  swept_gene_pops,
  by = c("gene_key", "MAG_id", "colony_id"),
  all.x = TRUE,
  suffixes = c("", "_sweep")
)
pnps_independent_family <- merge(
  unique(
    pnps_pop[
      any_pnps_gt1 == TRUE,  # Removed missing column check
      .(gene_key, MAG_id, colony_id)
    ]
  ),
  unique(
    cluster_dt[
      ,
      .(gene_key, cluster_id, genus)  # Changed Genus to genus
    ]
  ),
  by = "gene_key",
  all.x = TRUE
)[
  ,
  .(
    n_independent_pnps_genes = uniqueN(gene_key),
    n_independent_pnps_MAGs = uniqueN(MAG_id),
    n_independent_pnps_colonies = uniqueN(colony_id),
    n_independent_pnps_genera = uniqueN(genus, na.rm = TRUE)
  ),
  by = cluster_id
]



candidate_family_summary <- merge(
  candidate_family_summary,
  pnps_independent_family,
  by = "cluster_id",
  all.x = TRUE
)

for (cc in c(
  "n_independent_pnps_genes",
  "n_independent_pnps_MAGs",
  "n_independent_pnps_colonies",
  "n_independent_pnps_genera"
)) {
  if (!cc %in% names(candidate_family_summary)) {
    candidate_family_summary[, (cc) := 0L]
  }
  candidate_family_summary[
    is.na(get(cc)),
    (cc) := 0L
  ]
}

candidate_family_summary[
  ,
  independent_pnps_support :=
    n_independent_pnps_genera >= 2
]

# ------------------------------------------------------------
# C. Load protein FASTAs and extract candidate family members
# ------------------------------------------------------------

message("Loading eligible amino-acid sequences...")
aa_all <- read_named_aa_fasta(ALL_FASTA)
aa_rep <- read_named_aa_fasta(REP_FASTA)

message("Unique all-sequence proteins: ", length(aa_all))
message("Representative proteins: ", length(aa_rep))

# Sequence mapping QC.
seq_map <- cluster_dt[
  ,
  .(
    cluster_id = cluster_id[1L],
    cluster_rep = cluster_rep[1L],
    is_mmseq_rep = any(is_cluster_rep %in% TRUE)
  ),
  by = sequence_id
]

seq_map[
  ,
  sequence_in_fasta := sequence_id %chin% names(aa_all)
]

seq_map[
  ,
  rep_in_rep_fasta := cluster_rep %chin% names(aa_rep)
]

message(
  "Cluster-linked unique sequences present in all-seq FASTA: ",
  sum(seq_map$sequence_in_fasta),
  " / ", nrow(seq_map),
  " (", round(mean(seq_map$sequence_in_fasta) * 100, 2), "%)"
)

message(
  "Cluster representatives present in representative FASTA: ",
  sum(seq_map$rep_in_rep_fasta),
  " / ", nrow(seq_map)
)

if (mean(seq_map$sequence_in_fasta) < 0.99) {
  warning("Unexpectedly many cluster-linked sequences missing from all-seq FASTA.")
}

# Only families that can be aligned are taken forward.
candidate_family_summary <- candidate_family_summary[
  cluster_id %chin% unique(seq_map$cluster_id)
]

# ------------------------------------------------------------
# Rank before sequence alignment
# ------------------------------------------------------------

# Evidence-oriented components. These are descriptive ranking quantities only.
candidate_family_summary[
  ,
  recurrence_rank_component :=
    frank(
      -(
        4 * max_concurrent_genera +
          2 * max_concurrent_MAGs +
          1 * max_concurrent_colonies +
          log1p(n_concurrent_intervals)
      ),
      ties.method = "average"
    )
]

candidate_family_summary[
  ,
  sweep_burden_rank := frank(
    -log1p(total_0D_event_records),
    ties.method = "average"
  )
]

# Add function labels later if available; keep the scientific ranking independent
# of annotation to avoid circularity.

setorder(
  candidate_family_summary,
  -FDR05,
  -max_concurrent_genera,
  -max_concurrent_MAGs,
  -max_concurrent_colonies,
  -n_concurrent_intervals,
  -total_0D_event_records
)

prealignment_candidates <- head(
  candidate_family_summary,
  max(N_ALIGN_FAMILIES * 3L, 60L)
)


fwrite(
  candidate_family_summary,
  file.path(OUTDIR, "candidate_family_summary_pre_alignment.csv.gz")
)

# Member-level metadata for these families.
member_meta <- seq_map[
  cluster_id %chin% prealignment_candidates$cluster_id &
    sequence_in_fasta == TRUE,
  .(
    cluster_id,
    sequence_id,
    cluster_rep,
    is_mmseq_rep,
    sequence_length = width(aa_all[sequence_id])
  )
]
member_meta
# Number of physical genes / independent populations represented by each unique
# protein sequence.

# 1. Extract distinct gene_key -> colony_id mapping from sweep_dt
gene_colony_map <- unique(sweep_dt[!is.na(colony_id), .(gene_key, colony_id)])

# 2. Merge mapping into cluster_dt and aggregate
physical_support <- merge(
  cluster_dt[sequence_id %chin% member_meta$sequence_id],
  gene_colony_map,
  by = "gene_key",
  all.x = TRUE
)[
  , colony_key := fifelse(is.na(colony_id), NA_character_, paste(MAG_id, colony_id, sep = "::"))
][
  , .(
    n_physical_genes = uniqueN(gene_key),
    n_MAGs           = uniqueN(MAG_id),
    n_colonies       = uniqueN(colony_key, na.rm = TRUE)
  ),
  by = sequence_id
]




member_meta <- merge(
  member_meta,
  physical_support,
  by = "sequence_id",
  all.x = TRUE
)

# ------------------------------------------------------------
# C continued: write family-member FASTAs and event inventories
# ------------------------------------------------------------

for (cid in prealignment_candidates$cluster_id) {
  mem <- member_meta[cluster_id == cid]
  if (nrow(mem) == 0L) next
  
  seqs <- aa_all[mem$sequence_id]
  names(seqs) <- mem$sequence_id
  
  fasta_out <- file.path(
    DIR_MEMBER_FASTA,
    paste0(cid, ".faa")
  )
  write_fasta(seqs, fasta_out)
  
  ev <- zeroD_events[
    cluster_id == cid,
    .(
      cluster_id,
      gene_key,
      sequence_id,
      MAG_id,
      colony_id,
      genus,
      Month_t1,
      Month_t2,
      is_cluster_rep
    )
  ]
  
  if (nrow(ev) > 0L) {
    fwrite(
      ev,
      file.path(
        DIR_EVENTS,
        paste0(cid, "_0D_events.csv")
      )
    )
  }
}


# ============================================================
# Attach protein-cluster information to zeroD
# ============================================================

cluster_map <- unique(
  cluster_dt[
    ,
    .(
      gene_key,
      cluster_id,
      sequence_id,
      contig,
      gene_start = start,
      gene_stop  = stop,
      direction,
      sequence_length
    )
  ]
)



sweep_dt
grep(
  "codon_order|base_pos",
  names(sweep_dt),
  value = TRUE
)





target_clusters <- c(
  "seq_fa06c9c3d284d238e540bc215d89e7f0a05ec050",
  "seq_68ab79ba456ef8093145192b017fd1656be5672b",
  "seq_480042084cedb0f09b60f3020bacd9ef0f42fa8f",
  "seq_4fadd4da5fa72844bb0ef8c53a7565b9ec5708fd"
)

target_genes <- cluster_map[
  cluster_id %chin% target_clusters,
  unique(gene_key)
]

focal_target <- unique(
  sweep_dt[
    degeneracy == "0D" &
      gene_key %chin% target_genes &
      !is.na(corresponding_gene_call) &
      !is.na(ancestral_codon) &
      !is.na(derived_codon),
    .(
      colony_id = as.character(colony_id),
      MAG_id = as.character(MAG_id),
      gene_key,
      corresponding_gene_call =
        as.character(corresponding_gene_call),
      contig_name,
      pos_in_contig,
      unique_SNV_identifier,
      Month_t1,
      Month_t2,
      ancestral_codon,
      derived_codon,
      ancestral_aa,
      derived_aa,
      freq_t1,
      freq_t2
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_in_contig"
  )
)

nrow(focal_target)
focal_target




CODON_AA_DIR <- paste0(
  "/home/robinch/projects/BEE-WITCH/",
  "40_internal_pnps_calculations_from_anvio/",
  "colony_outputs"
)

codon_files <- list.files(
  CODON_AA_DIR,
  pattern = "^colony_[0-9]+_snp_codon_aa_changes\\.csv$",
  full.names = TRUE
)

target_colonies <- unique(focal_target$colony_id)

target_codon_files <- codon_files[
  sub(
    "^colony_([0-9]+)_.*$",
    "\\1",
    basename(codon_files)
  ) %chin%
    target_colonies
]

target_codon_files

codon_aa_target <- rbindlist(
  lapply(
    target_codon_files,
    function(f) {
      
      x <- fread(
        f,
        showProgress = FALSE
      )
      
      colony_id <- sub(
        "^colony_([0-9]+)_.*$",
        "\\1",
        basename(f)
      )
      
      x[
        ,
        colony_id := colony_id
      ]
      
      x
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

names(codon_aa_target)








target_gene_calls <- unique(
  focal_target[
    ,
    .(
      colony_id,
      MAG_id,
      corresponding_gene_call
    )
  ]
)

target_snv_rows <- rbindlist(
  lapply(
    snv_files,
    function(f) {
      
      colony_from_file <- sub(
        "^colony_([0-9]+)_.*$",
        "\\1",
        basename(f)
      )
      
      if (!colony_from_file %chin%
          target_gene_calls$colony_id) {
        return(NULL)
      }
      
      x <- fread(
        f,
        select = c(
          "Colony",
          "Month",
          "MAG_id",
          "contig_name",
          "pos_in_contig",
          "unique_pos_identifier",
          "corresponding_gene_call",
          "codon_order_in_gene",
          "base_pos_in_codon",
          "polarized_major_allele",
          "polarized_minor_allele",
          "degeneracy"
        ),
        showProgress = FALSE
      )
      
      x[
        ,
        `:=`(
          Colony = as.character(Colony),
          MAG_id = as.character(MAG_id),
          corresponding_gene_call =
            as.character(corresponding_gene_call)
        )
      ]
      
      wanted_calls <- target_gene_calls[
        colony_id == colony_from_file,
        corresponding_gene_call
      ]
      
      x[
        corresponding_gene_call %chin%
          wanted_calls
      ]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

nrow(target_snv_rows)




target_gene_presence <- target_gene_calls[
  ,
  {
    
    x <- target_snv_rows[
      Colony == .BY$colony_id &
        MAG_id == .BY$MAG_id &
        corresponding_gene_call ==
        .BY$corresponding_gene_call
    ]
    
    .(
      n_rows = nrow(x),
      n_sites =
        uniqueN(x$unique_pos_identifier),
      n_codon_positions =
        uniqueN(x$codon_order_in_gene)
    )
  },
  by = .(
    colony_id,
    MAG_id,
    corresponding_gene_call
  )
]

target_gene_presence





FINAL_INPUT_DIR <- paste0(
  "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/",
  "ALLSPECIES_intrapop_diversity/",
  "colony_polarized_outputs/",
  "pooled_by_month_global_consensus"
)

target_colonies <- unique(
  focal_target$colony_id
)

final_input_files <- file.path(
  FINAL_INPUT_DIR,
  paste0(
    "colony_",
    target_colonies,
    "_pooled_final_input.csv"
  )
)

final_input_files <- final_input_files[
  file.exists(final_input_files)
]






map_cols <- c(
  "MAG_id",
  "contig_name",
  "pos_in_contig",
  "corresponding_gene_call",
  "codon_order_in_gene",
  "base_pos_in_codon",
  "polarized_major_allele",
  "polarized_minor_allele"
)




snv_position_map <- rbindlist(
  lapply(
    final_input_files,
    function(f) {
      
      colony_id <- sub(
        "^colony_([0-9]+)_.*$",
        "\\1",
        basename(f)
      )
      
      x <- fread(
        f,
        select = map_cols,
        showProgress = FALSE
      )
      
      x[
        ,
        `:=`(
          colony_id = colony_id,
          MAG_id =
            as.character(MAG_id),
          contig_name =
            as.character(contig_name),
          pos_in_contig =
            as.integer(pos_in_contig),
          corresponding_gene_call =
            as.character(
              corresponding_gene_call
            ),
          codon_order_in_gene =
            as.integer(
              codon_order_in_gene
            ),
          base_pos_in_codon =
            as.integer(
              base_pos_in_codon
            )
        )
      ]
      
      x
    }
  ),
  use.names = TRUE,
  fill = TRUE
)







snv_position_map <- unique(
  snv_position_map[
    ,
    .(
      colony_id,
      MAG_id,
      contig_name,
      pos_in_contig,
      corresponding_gene_call,
      codon_order_in_gene,
      base_pos_in_codon,
      polarized_major_allele,
      polarized_minor_allele
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_in_contig"
  )
)




snv_position_map[
  ,
  .(
    n_gene_calls =
      uniqueN(corresponding_gene_call),
    n_codon_orders =
      uniqueN(codon_order_in_gene),
    n_base_positions =
      uniqueN(base_pos_in_codon),
    n_major =
      uniqueN(polarized_major_allele),
    n_minor =
      uniqueN(polarized_minor_allele)
  ),
  by = .(
    colony_id,
    MAG_id,
    contig_name,
    pos_in_contig
  )
][
  n_gene_calls > 1 |
    n_codon_orders > 1 |
    n_base_positions > 1 |
    n_major > 1 |
    n_minor > 1
]



candidate_events <- merge(
  focal_target,
  snv_position_map[
    ,
    .(
      colony_id,
      MAG_id,
      contig_name,
      pos_in_contig,
      codon_order_in_gene,
      base_pos_in_codon
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_in_contig"
  ),
  all.x = TRUE
)



candidate_events[
  ,
  aa_position :=
    as.integer(codon_order_in_gene)
]

candidate_events[
  ,
  .(
    n_events = .N,
    n_missing_position =
      sum(is.na(aa_position))
  )
]

codon_validation <- codon_aa_target[
  degeneracy == "0D",
  .(
    colony_id,
    corresponding_gene_call,
    codon_order,
    ancestral_codon,
    derived_codon,
    ancestral_aa,
    derived_aa
  )
]



candidate_events[
  ,
  key_no_position := paste(
    colony_id,
    corresponding_gene_call,
    ancestral_codon,
    derived_codon,
    sep = "|"
  )
]

codon_validation[
  ,
  key_no_position := paste(
    colony_id,
    corresponding_gene_call,
    ancestral_codon,
    derived_codon,
    sep = "|"
  )
]

candidate_events[
  ,
  .(
    n = .N,
    matched =
      sum(
        key_no_position %chin%
          codon_validation$key_no_position
      )
  )
]



candidate_events <- merge(
  candidate_events,
  cluster_map[
    ,
    .(
      gene_key,
      cluster_id,
      sequence_id,
      sequence_length
    )
  ],
  by = "gene_key",
  all.x = TRUE
)



candidate_events[
  ,
  .(
    n_events = .N,
    n_missing_aa =
      sum(is.na(aa_position)),
    n_out_of_bounds =
      sum(
        aa_position < 1 |
          aa_position > sequence_length,
        na.rm = TRUE
      ),
    n_missing_sequence =
      sum(is.na(sequence_id)),
    n_missing_cluster =
      sum(is.na(cluster_id))
  )
]

candidate_events[
  order(cluster_id, aa_position),
  .(
    cluster_id,
    gene_key,
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    unique_SNV_identifier,
    sequence_id,
    sequence_length,
    aa_position,
    base_pos_in_codon,
    ancestral_codon,
    derived_codon,
    ancestral_aa,
    derived_aa
  )
]



candidate_events
fwrite(
  candidate_events[
    ,
    .(
      cluster_id,
      corresponding_gene_call,
      gene_key,
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      unique_SNV_identifier,
      sequence_id,
      sequence_length,
      aa_position,
      base_pos_in_codon,
      ancestral_codon,
      derived_codon,
      ancestral_aa,
      derived_aa
    )
  ],
  "candidate_protein_0D_deltaRAF_events.csv"
)

