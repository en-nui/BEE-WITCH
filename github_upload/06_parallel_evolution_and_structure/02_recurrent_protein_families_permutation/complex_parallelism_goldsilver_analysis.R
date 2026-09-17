library(data.table)
library(ggplot2)

# ============================================================
# 1. UNIQUE FOCAL 0D SWEEPS
# ============================================================

focal_0D <- unique(
  sweep_dt[
    degeneracy == "0D" &
      !is.na(contig_name) &
      !is.na(pos_in_contig) &
      !is.na(Month_t1) &
      !is.na(Month_t2),
    .(
      colony_id = as.character(colony_id),
      MAG_id = as.character(MAG_id),
      
      Month_t1 = as.character(Month_t1),
      Month_t2 = as.character(Month_t2),
      
      contig_name = as.character(contig_name),
      pos_0D = as.integer(pos_in_contig),
      
      unique_SNV_identifier,
      
      freq_0D_t1 = freq_t1,
      freq_0D_t2 = freq_t2,
      
      gene_key,
      corresponding_gene_call
    )
  ]
)

focal_0D[
  ,
  delta_0D := freq_0D_t2 - freq_0D_t1
]

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

snv_files
length(snv_files)

focal_0D[
  ,
  abs_delta_0D := abs(delta_0D)
]

focal_0D[
  ,
  focal_id := .I
]

nrow(focal_0D)

snv_cols <- c(
  "Colony",
  "Month",
  "contig_name",
  "pos_in_contig",
  "coverage",
  "unique_pos_identifier",
  "corresponding_gene_call",
  "MAG_id",
  "polarized_major_allele",
  "polarized_minor_allele",
  "departure_from_polarized_consensus",
  "degeneracy"
)


message("Loading 4D SNVs...")

fourD_raw <- rbindlist(
  lapply(
    snv_files,
    function(f) {
      
      x <- fread(
        f,
        select = snv_cols,
        showProgress = FALSE
      )
      
      x[
        degeneracy == "4D"
      ]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

message(
  "Loaded ",
  format(nrow(fourD_raw), big.mark = ","),
  " 4D observations"
)



fourD_raw[
  ,
  Colony := as.character(Colony)
]

fourD_raw[
  ,
  MAG_id := as.character(MAG_id)
]

fourD_raw[
  ,
  Month := as.character(Month)
]

fourD_raw[
  ,
  contig_name := as.character(contig_name)
]

fourD_raw[
  ,
  pos_in_contig :=
    as.integer(pos_in_contig)
]

fourD_raw[
  ,
  AF :=
    as.numeric(
      departure_from_polarized_consensus
    )
]

summary(fourD_raw$AF)

range(
  fourD_raw$AF,
  na.rm = TRUE
)




dup_check <- fourD_raw[
  ,
  .N,
  by = .(
    Colony,
    MAG_id,
    Month,
    contig_name,
    pos_in_contig
  )
][
  N > 1
]

nrow(dup_check)

if (nrow(dup_check) > 0) {
  print(
    dup_check[1:20]
  )
}




polarization_check <- fourD_raw[
  ,
  .(
    n_major =
      uniqueN(
        polarized_major_allele
      ),
    
    n_minor =
      uniqueN(
        polarized_minor_allele
      )
  ),
  by = .(
    Colony,
    MAG_id,
    contig_name,
    pos_in_contig
  )
][
  n_major > 1 |
    n_minor > 1
]

nrow(polarization_check)




#dupCHECK

dup_detail <- fourD_raw[
  ,
  .(
    N = .N,
    
    n_AF = uniqueN(AF),
    n_cov = uniqueN(coverage),
    
    AF_min = min(AF, na.rm = TRUE),
    AF_max = max(AF, na.rm = TRUE),
    
    cov_min = min(coverage, na.rm = TRUE),
    cov_max = max(coverage, na.rm = TRUE),
    
    n_gene_calls =
      uniqueN(corresponding_gene_call),
    
    n_unique_pos_ids =
      uniqueN(unique_pos_identifier)
  ),
  by = .(
    Colony,
    MAG_id,
    Month,
    contig_name,
    pos_in_contig
  )
]

dup_detail[N > 1][1:30]



dup_detail[
  N > 1,
  .(
    duplicated_sites = .N,
    
    exact_same_AF =
      mean(n_AF == 1),
    
    multiple_AF_values =
      mean(n_AF > 1),
    
    exact_same_coverage =
      mean(n_cov == 1),
    
    multiple_coverages =
      mean(n_cov > 1),
    
    multiple_gene_calls =
      mean(n_gene_calls > 1),
    
    multiple_unique_pos_ids =
      mean(n_unique_pos_ids > 1)
  )
]




fourD_raw[
  Colony == "100" &
    MAG_id == "mmag_302" &
    Month == "July" &
    contig_name == "c_107328" &
    pos_in_contig == 2016
]
fourD_raw[
  Colony == "100" &
    MAG_id == "mmag_302" &
    Month == "May" &
    contig_name == "c_107332" &
    pos_in_contig == 38660
]

fourD_raw[
  Colony == "100" &
    MAG_id == "mmag_302" &
    Month == "May" &
    contig_name == "c_107332" &
    pos_in_contig == 38660,
  .(
    coverage,
    AF,
    unique_pos_identifier,
    corresponding_gene_call,
    polarized_major_allele,
    polarized_minor_allele
  )
]



fourD_clean <- unique(
  fourD_raw,
  by = c(
    "Colony",
    "MAG_id",
    "Month",
    "contig_name",
    "pos_in_contig"
  )
)
stopifnot(
  fourD_clean[
    ,
    .N,
    by = .(
      Colony,
      MAG_id,
      Month,
      contig_name,
      pos_in_contig
    )
  ][
    N > 1,
    .N
  ] == 0
)
fourD_clean[
  ,
  AF := as.numeric(
    departure_from_polarized_consensus
  )
]

summary(fourD_clean$AF)
range(fourD_clean$AF, na.rm = TRUE)


needed_endpoints <- unique(
  rbind(
    focal_0D[
      ,
      .(
        Colony = colony_id,
        MAG_id,
        Month = Month_t1
      )
    ],
    
    focal_0D[
      ,
      .(
        Colony = colony_id,
        MAG_id,
        Month = Month_t2
      )
    ]
  )
)

setkey(
  fourD_clean,
  Colony,
  MAG_id,
  Month
)

setkey(
  needed_endpoints,
  Colony,
  MAG_id,
  Month
)

fourD_sub <- fourD_clean[
  needed_endpoints,
  nomatch = 0
]



c(
  all_4D_rows = nrow(fourD_clean),
  focal_endpoint_rows = nrow(fourD_sub)
)




focal_transitions <- unique(
  focal_0D[
    ,
    .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2
    )
  ]
)

focal_transitions[
  ,
  transition_id := .I
]



fourD_t1 <- merge(
  focal_transitions,
  
  fourD_sub[
    ,
    .(
      colony_id = Colony,
      MAG_id,
      Month_t1 = Month,
      contig_name,
      pos_in_contig,
      
      AF_t1 = AF,
      coverage_t1 = coverage
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1"
  ),
  
  allow.cartesian = TRUE
)



fourD_t2 <- fourD_sub[
  ,
  .(
    colony_id = Colony,
    MAG_id,
    Month_t2 = Month,
    contig_name,
    pos_in_contig,
    
    AF_t2 = AF,
    coverage_t2 = coverage
  )
]




fourD_pairs <- merge(
  fourD_t1,
  fourD_t2,
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t2",
    "contig_name",
    "pos_in_contig"
  ),
  
  all = FALSE
)

fourD_pairs[
  ,
  delta_4D := AF_t2 - AF_t1
]

fourD_pairs[
  ,
  abs_delta_4D := abs(delta_4D)
]


summary(fourD_pairs$abs_delta_4D)

range(
  fourD_pairs$abs_delta_4D,
  na.rm = TRUE
)




fourD_transition_summary <- fourD_pairs[
  ,
  .(
    n_trackable_4D = .N,
    
    median_abs_delta_4D =
      median(
        abs_delta_4D,
        na.rm = TRUE
      )
  ),
  by = .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2
  )
]

summary(
  fourD_transition_summary$
    n_trackable_4D
)




quantile(
  fourD_transition_summary$
    n_trackable_4D,
  c(
    0,
    .1,
    .25,
    .5,
    .75,
    .9,
    1
  )
)




focal_for_join <- focal_0D[
  ,
  .(
    focal_id,
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_0D,
    delta_0D,
    abs_delta_0D
  )
]






local_4D <- merge(
  focal_for_join,
  
  fourD_pairs[
    ,
    .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      
      pos_4D = pos_in_contig,
      
      AF_4D_t1 = AF_t1,
      AF_4D_t2 = AF_t2,
      
      delta_4D,
      abs_delta_4D
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name"
  ),
  
  allow.cartesian = TRUE
)

local_4D[
  ,
  distance_bp :=
    abs(pos_4D - pos_0D)
]


c(
  total_focal_0D =
    nrow(focal_0D),
  
  focal_0D_with_trackable_4D =
    uniqueN(local_4D$focal_id),
  
  local_4D_comparisons =
    nrow(local_4D)
)



MIN_TRACKABLE_4D <- 50

eligible_transitions <- fourD_transition_summary[
  n_trackable_4D >= MIN_TRACKABLE_4D,
  .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2
  )
]

fourD_pairs_filt <- merge(
  fourD_pairs,
  eligible_transitions,
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2"
  ),
  all = FALSE
)



c(
  total_transitions =
    nrow(fourD_transition_summary),
  
  eligible_transitions =
    nrow(eligible_transitions),
  
  fraction_retained =
    nrow(eligible_transitions) /
    nrow(fourD_transition_summary)
)



transition_4D_background <- fourD_pairs_filt[
  ,
  .(
    genome_median_4D =
      median(
        abs_delta_4D,
        na.rm = TRUE
      ),
    
    genome_mean_4D =
      mean(
        abs_delta_4D,
        na.rm = TRUE
      ),
    
    n_genome_4D = .N
  ),
  by = .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2
  )
]




local_4D <- merge(
  focal_for_join,
  fourD_pairs_filt[
    ,
    .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      pos_4D = pos_in_contig,
      AF_4D_t1 = AF_t1,
      AF_4D_t2 = AF_t2,
      delta_4D,
      abs_delta_4D
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name"
  ),
  allow.cartesian = TRUE
)

local_4D[
  ,
  distance_bp :=
    abs(pos_4D - pos_0D)
]

local_4D <- merge(
  local_4D,
  transition_4D_background,
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)

local_4D[
  ,
  excess_abs_delta_4D :=
    abs_delta_4D -
    genome_median_4D
]






distance_breaks <- c(
  0,
  100,
  250,
  500,
  1000,
  2500,
  5000,
  10000,
  25000,
  50000,
  100000,
  Inf
)

distance_labels <- c(
  "<100 bp",
  "0.1-0.25 kb",
  "0.25-0.5 kb",
  "0.5-1 kb",
  "1-2.5 kb",
  "2.5-5 kb",
  "5-10 kb",
  "10-25 kb",
  "25-50 kb",
  "50-100 kb",
  ">100 kb"
)

local_4D[
  ,
  distance_bin := cut(
    distance_bp,
    breaks = distance_breaks,
    labels = distance_labels,
    right = FALSE
  )
]

focal_distance <- local_4D[
  ,
  .(
    local_median_4D =
      median(
        abs_delta_4D,
        na.rm = TRUE
      ),
    
    local_excess_4D =
      median(
        excess_abs_delta_4D,
        na.rm = TRUE
      ),
    
    n_4D = .N
  ),
  by = .(
    focal_id,
    distance_bin
  )
]

distance_summary <- focal_distance[
  ,
  .(
    median_abs_delta =
      median(
        local_median_4D,
        na.rm = TRUE
      ),
    
    median_excess =
      median(
        local_excess_4D,
        na.rm = TRUE
      ),
    
    q25_excess =
      quantile(
        local_excess_4D,
        0.25,
        na.rm = TRUE
      ),
    
    q75_excess =
      quantile(
        local_excess_4D,
        0.75,
        na.rm = TRUE
      ),
    
    n_focal =
      uniqueN(focal_id)
  ),
  by = distance_bin
]

distance_summary




p_local_4D <- ggplot(
  distance_summary,
  aes(
    x = distance_bin,
    y = median_excess,
    group = 1
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_ribbon(
    aes(
      ymin = q25_excess,
      ymax = q75_excess
    ),
    alpha = 0.18
  ) +
  
  geom_line(
    linewidth = 1.1
  ) +
  
  geom_point(
    size = 2.7
  ) +
  
  labs(
    x = "Distance from focal 0D sweep",
    y = expression(
      "Excess synonymous " *
        "|" * Delta * "AF|"
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

p_local_4D






local_4D[
  ,
  same_direction :=
    sign(delta_0D) ==
    sign(delta_4D)
]
direction_by_focal <- local_4D[
  abs(delta_4D) > 0,
  .(
    same_direction_fraction =
      mean(same_direction)
  ),
  by = .(
    focal_id,
    distance_bin
  )
]
direction_summary <- direction_by_focal[
  ,
  .(
    median_same_direction =
      median(
        same_direction_fraction,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        same_direction_fraction,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        same_direction_fraction,
        0.75,
        na.rm = TRUE
      ),
    
    n_focal =
      uniqueN(focal_id)
  ),
  by = distance_bin
]

direction_summary
distance_summary[
  ,
  .(
    distance_bin,
    median_excess,
    q25_excess,
    q75_excess,
    n_focal
  )
]





#REBIO:DOMG TP CPMTOMIE


# ============================================================
# REBUILD 4D PAIRS WITH GENE IDENTITY
# ============================================================

fourD_t1 <- merge(
  focal_transitions,
  
  fourD_sub[
    ,
    .(
      colony_id = Colony,
      MAG_id,
      Month_t1 = Month,
      contig_name,
      pos_in_contig,
      
      gene_4D =
        as.character(corresponding_gene_call),
      
      AF_t1 = AF,
      coverage_t1 = coverage
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1"
  ),
  
  allow.cartesian = TRUE
)


fourD_t2 <- fourD_sub[
  ,
  .(
    colony_id = Colony,
    MAG_id,
    Month_t2 = Month,
    contig_name,
    pos_in_contig,
    
    gene_4D =
      as.character(corresponding_gene_call),
    
    AF_t2 = AF,
    coverage_t2 = coverage
  )
]

fourD_pairs <- merge(
  fourD_t1,
  fourD_t2,
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t2",
    "contig_name",
    "pos_in_contig",
    "gene_4D"
  ),
  
  all = FALSE
)

fourD_pairs[
  ,
  delta_4D := AF_t2 - AF_t1
]

fourD_pairs[
  ,
  abs_delta_4D := abs(delta_4D)
]





focal_for_join <- focal_0D[
  ,
  .(
    focal_id,
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_0D,
    
    gene_0D =
      as.character(corresponding_gene_call),
    
    delta_0D,
    abs_delta_0D
  )
]



fourD_transition_summary <- fourD_pairs[
  ,
  .(
    n_trackable_4D = .N,
    
    median_abs_delta_4D =
      median(
        abs_delta_4D,
        na.rm = TRUE
      )
  ),
  by = .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2
  )
]


MIN_TRACKABLE_4D <- 50

eligible_transitions <- fourD_transition_summary[
  n_trackable_4D >= MIN_TRACKABLE_4D,
  .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2
  )
]


fourD_pairs_filt <- merge(
  fourD_pairs,
  eligible_transitions,
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2"
  ),
  
  all = FALSE
)

local_4D <- merge(
  focal_for_join,
  
  fourD_pairs_filt[
    ,
    .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      
      pos_4D = pos_in_contig,
      gene_4D,
      
      AF_4D_t1 = AF_t1,
      AF_4D_t2 = AF_t2,
      
      delta_4D,
      abs_delta_4D
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name"
  ),
  
  allow.cartesian = TRUE
)


local_4D[
  ,
  distance_bp :=
    abs(pos_4D - pos_0D)
]


local_4D[
  ,
  same_gene :=
    !is.na(gene_0D) &
    !is.na(gene_4D) &
    gene_0D == gene_4D
]




local_4D[
  ,
  .(
    n_pairs = .N,
    n_same_gene = sum(same_gene),
    n_different_gene = sum(!same_gene),
    prop_same_gene = mean(same_gene)
  )
]



local_4D[
  ,
  distance_bin := cut(
    distance_bp,
    breaks = distance_breaks,
    labels = distance_labels,
    right = FALSE
  )
]


gene_context_summary <- local_4D[
  ,
  .(
    n_pairs = .N,
    n_focal = uniqueN(focal_id),
    
    prop_same_gene =
      mean(same_gene),
    
    n_same_gene =
      sum(same_gene),
    
    n_different_gene =
      sum(!same_gene)
  ),
  by = distance_bin
]

gene_context_summary






transition_4D_background <- fourD_pairs_filt[
  ,
  .(
    genome_median_4D =
      median(
        abs_delta_4D,
        na.rm = TRUE
      ),
    
    genome_mean_4D =
      mean(
        abs_delta_4D,
        na.rm = TRUE
      ),
    
    n_genome_4D = .N
  ),
  by = .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2
  )
]


local_4D <- merge(
  local_4D,
  transition_4D_background,
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2"
  ),
  
  all.x = TRUE
)


local_4D[
  ,
  excess_abs_delta_4D :=
    abs_delta_4D -
    genome_median_4D
]






summarize_spatial_response <- function(dt) {
  
  by_focal <- dt[
    ,
    .(
      local_excess_4D =
        median(
          excess_abs_delta_4D,
          na.rm = TRUE
        ),
      
      local_abs_delta_4D =
        median(
          abs_delta_4D,
          na.rm = TRUE
        ),
      
      n_4D = .N
    ),
    by = .(
      focal_id,
      distance_bin
    )
  ]
  
  by_focal[
    ,
    .(
      median_excess =
        median(
          local_excess_4D,
          na.rm = TRUE
        ),
      
      q25_excess =
        quantile(
          local_excess_4D,
          0.25,
          na.rm = TRUE
        ),
      
      q75_excess =
        quantile(
          local_excess_4D,
          0.75,
          na.rm = TRUE
        ),
      
      n_focal =
        uniqueN(focal_id)
    ),
    by = distance_bin
  ]
}




spatial_all <- summarize_spatial_response(
  local_4D
)

spatial_same_gene <- summarize_spatial_response(
  local_4D[
    same_gene == TRUE
  ]
)

spatial_diff_gene <- summarize_spatial_response(
  local_4D[
    same_gene == FALSE
  ]
)



spatial_same_gene
spatial_diff_gene


spatial_all[
  ,
  context := "All same-contig"
]

spatial_same_gene[
  ,
  context := "Same gene"
]

spatial_diff_gene[
  ,
  context := "Different gene, same contig"
]


spatial_context <- rbindlist(
  list(
    spatial_all,
    spatial_same_gene,
    spatial_diff_gene
  ),
  use.names = TRUE
)




plot_bins <- distance_labels[
  1:8
]

ggplot(
  spatial_context[
    distance_bin %in% plot_bins
  ],
  aes(
    x = distance_bin,
    y = median_excess,
    group = context,
    linetype = context,
    shape = context
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2.5
  ) +
  
  labs(
    x = "Distance from focal nonsynonymous sweep",
    y = expression(
      "Excess synonymous " *
        "|" * Delta * "AF|"
    ),
    linetype = NULL,
    shape = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )






message("Loading 0D SNVs...")

zeroD_raw_all <- rbindlist(
  lapply(
    snv_files,
    function(f) {
      
      x <- fread(
        f,
        select = snv_cols,
        showProgress = FALSE
      )
      
      x[
        degeneracy == "0D"
      ]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

message(
  "Loaded ",
  format(
    nrow(zeroD_raw_all),
    big.mark = ","
  ),
  " 0D observations"
)



zeroD_raw_all[
  ,
  `:=`(
    Colony =
      as.character(Colony),
    
    MAG_id =
      as.character(MAG_id),
    
    Month =
      as.character(Month),
    
    contig_name =
      as.character(contig_name),
    
    pos_in_contig =
      as.integer(pos_in_contig),
    
    corresponding_gene_call =
      as.character(corresponding_gene_call),
    
    AF =
      as.numeric(
        departure_from_polarized_consensus
      )
  )
]

zeroD_clean_all <- unique(
  zeroD_raw_all,
  by = c(
    "Colony",
    "MAG_id",
    "Month",
    "contig_name",
    "pos_in_contig"
  )
)

zeroD_clean_all[
  ,
  .N,
  by = .(
    Colony,
    MAG_id,
    Month,
    contig_name,
    pos_in_contig
  )
][
  N > 1
]


setkey(
  zeroD_clean_all,
  Colony,
  MAG_id,
  Month
)

zeroD_sub_all <- zeroD_clean_all[
  needed_endpoints,
  nomatch = 0
]



zeroD_t1_all <- merge(
  focal_transitions,
  
  zeroD_sub_all[
    ,
    .(
      colony_id = Colony,
      MAG_id,
      Month_t1 = Month,
      contig_name,
      pos_in_contig,
      
      gene_0D_control =
        corresponding_gene_call,
      
      AF_t1 = AF
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1"
  ),
  
  allow.cartesian = TRUE
)


zeroD_t2_all <- zeroD_sub_all[
  ,
  .(
    colony_id = Colony,
    MAG_id,
    Month_t2 = Month,
    contig_name,
    pos_in_contig,
    
    gene_0D_control =
      corresponding_gene_call,
    
    AF_t2 = AF
  )
]




zeroD_pairs_all <- merge(
  zeroD_t1_all,
  zeroD_t2_all,
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t2",
    "contig_name",
    "pos_in_contig",
    "gene_0D_control"
  ),
  
  all = FALSE
)


zeroD_pairs_all[
  ,
  `:=`(
    delta_0D_control =
      AF_t2 - AF_t1,
    
    abs_delta_0D_control =
      abs(AF_t2 - AF_t1)
  )
]



focal_0D[
  ,
  sweep_site_key := paste(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_0D,
    sep = "|"
  )
]


zeroD_pairs_all[
  ,
  site_key := paste(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_in_contig,
    sep = "|"
  )
]


zeroD_pairs_all[
  ,
  is_focal_sweep :=
    site_key %chin%
    focal_0D$sweep_site_key
]





zeroD_pairs_all[
  ,
  .(
    n_trackable_0D = .N,
    n_focal_sweeps =
      sum(is_focal_sweep),
    n_nonsweep_controls =
      sum(!is_focal_sweep)
  )
]



control_0D <- zeroD_pairs_all[
  is_focal_sweep == FALSE
]





control_availability <- focal_0D[
  ,
  .(
    focal_id,
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_0D,
    freq_0D_t1
  )
][
  
  control_0D,
  on = .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2,
    contig_name
  ),
  allow.cartesian = TRUE,
  nomatch = NA
  
][
  ,
  .(
    n_same_contig_controls =
      sum(!is.na(pos_in_contig)),
    
    nearest_control_AF_diff = {
      x <- abs(
        AF_t1 - freq_0D_t1
      )
      
      if (all(is.na(x))) {
        NA_real_
      } else {
        min(
          x,
          na.rm = TRUE
        )
      }
    }
  ),
  by = focal_id
]




summary(
  control_availability$
    n_same_contig_controls
)

quantile(
  control_availability$
    n_same_contig_controls,
  c(
    0,
    .1,
    .25,
    .5,
    .75,
    .9,
    1
  ),
  na.rm = TRUE
)

mean(
  control_availability$
    n_same_contig_controls > 0
)




focal_match_dt <- focal_0D[
  ,
  .(
    focal_id,
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2,
    contig_name,
    pos_0D,
    freq_0D_t1
  )
]

control_match_dt <- unique(
  control_0D[
    ,
    .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      
      pos_control =
        pos_in_contig,
      
      AF_control_t1 =
        AF_t1
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name",
    "pos_control"
  )
)
focal_control_pairs <- merge(
  focal_match_dt,
  control_match_dt,
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name"
  ),
  
  all.x = TRUE,
  allow.cartesian = TRUE
)
nrow(focal_control_pairs)

sum(
  is.na(
    focal_control_pairs$focal_id
  )
)


focal_control_pairs[
  ,
  control_AF_diff :=
    abs(
      AF_control_t1 -
        freq_0D_t1
    )
]


control_availability <- focal_control_pairs[
  ,
  .(
    n_same_contig_controls =
      sum(
        !is.na(pos_control)
      ),
    
    nearest_control_AF_diff = {
      x <- control_AF_diff[
        !is.na(pos_control)
      ]
      
      if (length(x) == 0) {
        NA_real_
      } else {
        min(
          x,
          na.rm = TRUE
        )
      }
    }
  ),
  by = focal_id
]


summary(
  control_availability$
    n_same_contig_controls
)



quantile(
  control_availability$
    n_same_contig_controls,
  c(
    0,
    .1,
    .25,
    .5,
    .75,
    .9,
    1
  ),
  na.rm = TRUE
)
mean(
  control_availability$
    n_same_contig_controls > 0
)


summary(
  control_availability$
    nearest_control_AF_diff
)

quantile(
  control_availability$
    nearest_control_AF_diff,
  c(
    0,
    .1,
    .25,
    .5,
    .75,
    .9,
    .95,
    1
  ),
  na.rm = TRUE
)


focal_control_pairs[
  !is.na(pos_control),
  .(
    n_pairs = .N,
    
    same_position =
      sum(
        pos_control == pos_0D
      )
  )
]






matchable_focal_ids <- control_availability[
  n_same_contig_controls > 0,
  focal_id
]

length(matchable_focal_ids)

length(matchable_focal_ids) /
  nrow(focal_0D)





spatial_matchable_observed <- summarize_spatial_response(
  local_4D[
    focal_id %in% matchable_focal_ids
  ]
)

spatial_matchable_diff_gene <- summarize_spatial_response(
  local_4D[
    focal_id %in% matchable_focal_ids &
      same_gene == FALSE
  ]
)

merge(
  spatial_diff_gene[
    ,
    .(
      distance_bin,
      all_focal = median_excess,
      n_all = n_focal
    )
  ],
  spatial_matchable_diff_gene[
    ,
    .(
      distance_bin,
      matchable_focal = median_excess,
      n_matchable = n_focal
    )
  ],
  by = "distance_bin",
  all = TRUE
)
















control_sites <- unique(
  control_0D[
    ,
    .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      
      pos_control =
        pos_in_contig,
      
      gene_control =
        as.character(
          gene_0D_control
        ),
      
      AF_control_t1 =
        AF_t1,
      
      AF_control_t2 =
        AF_t2,
      
      delta_control =
        delta_0D_control
    )
  ],
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name",
    "pos_control"
  )
)

control_sites[
  ,
  control_id := .I
]

nrow(control_sites)




stopifnot(
  control_sites[
    ,
    .N,
    by = .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      pos_control
    )
  ][N > 1, .N] == 0
)




control_local_4D <- merge(
  control_sites,
  
  fourD_pairs_filt[
    ,
    .(
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      
      pos_4D =
        pos_in_contig,
      
      gene_4D,
      
      abs_delta_4D
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name"
  ),
  
  allow.cartesian = TRUE
)


control_local_4D[
  ,
  distance_bp :=
    abs(
      pos_4D -
        pos_control
    )
]


control_local_4D <- merge(
  control_local_4D,
  transition_4D_background,
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2"
  ),
  
  all.x = TRUE
)


control_local_4D[
  ,
  excess_abs_delta_4D :=
    abs_delta_4D -
    genome_median_4D
]

control_local_4D[
  ,
  same_gene :=
    !is.na(gene_control) &
    !is.na(gene_4D) &
    gene_control == gene_4D
]

control_local_4D[
  ,
  distance_bin :=
    cut(
      distance_bp,
      breaks = distance_breaks,
      labels = distance_labels,
      right = FALSE
    )
]

control_response_diff_gene <- control_local_4D[
  same_gene == FALSE,
  .(
    local_excess_4D =
      median(
        excess_abs_delta_4D,
        na.rm = TRUE
      ),
    
    n_4D =
      .N
  ),
  by = .(
    control_id,
    distance_bin
  )
]


focal_candidates <- merge(
  focal_0D[
    focal_id %in%
      matchable_focal_ids,
    .(
      focal_id,
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      pos_0D,
      freq_0D_t1
    )
  ],
  
  control_sites[
    ,
    .(
      control_id,
      colony_id,
      MAG_id,
      Month_t1,
      Month_t2,
      contig_name,
      pos_control,
      AF_control_t1
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "Month_t1",
    "Month_t2",
    "contig_name"
  ),
  
  allow.cartesian = TRUE
)

focal_candidates[
  ,
  AF_diff :=
    abs(
      AF_control_t1 -
        freq_0D_t1
    )
]

setorder(
  focal_candidates,
  focal_id,
  AF_diff
)

focal_candidates[
  ,
  AF_rank :=
    seq_len(.N),
  by = focal_id
]

setorder(
  focal_candidates,
  focal_id,
  AF_diff
)

focal_candidates[
  ,
  AF_rank :=
    seq_len(.N),
  by = focal_id
]


focal_candidates[
  AF_rank <= 5,
  .(
    n_candidates =
      .N,
    
    best_AF_diff =
      min(AF_diff),
    
    worst_top5_AF_diff =
      max(AF_diff)
  ),
  by = focal_id
][
  ,
  .(
    median_candidates =
      median(n_candidates),
    
    median_best =
      median(best_AF_diff),
    
    median_worst_top5 =
      median(worst_top5_AF_diff)
  )
]



TOP_K <- 5L

candidate_pool <- focal_candidates[
  AF_rank <= TOP_K
]


controls_with_response <- unique(
  control_response_diff_gene$
    control_id
)

candidate_pool <- candidate_pool[
  control_id %in%
    controls_with_response
]

candidate_pool[
  ,
  uniqueN(focal_id)
]

candidate_pool[
  ,
  uniqueN(focal_id)
] /
  length(matchable_focal_ids)



inference_focal_ids <- unique(
  candidate_pool$focal_id
)



obs_inference <- summarize_spatial_response(
  local_4D[
    focal_id %in%
      inference_focal_ids &
      same_gene == FALSE
  ]
)

obs_inference



set.seed(12345)

N_PERM <- 1000L



analysis_bins <- c(
  "<100 bp",
  "0.1-0.25 kb",
  "0.25-0.5 kb",
  "0.5-1 kb",
  "1-2.5 kb",
  "2.5-5 kb",
  "5-10 kb"
)


null_mat <- matrix(
  NA_real_,
  nrow = N_PERM,
  ncol = length(analysis_bins)
)

colnames(null_mat) <-
  analysis_bins

sample_one_control <- function(dt) {
  
  dt[
    ,
    .SD[
      sample.int(
        .N,
        1L
      )
    ],
    by = focal_id
  ]
}


for (b in seq_len(N_PERM)) {
  
  picked <- sample_one_control(
    candidate_pool
  )
  
  perm_response <- merge(
    picked[
      ,
      .(
        focal_id,
        control_id
      )
    ],
    
    control_response_diff_gene,
    
    by = "control_id",
    all = FALSE,
    allow.cartesian = TRUE
  )
  
  perm_summary <- perm_response[
    distance_bin %in%
      analysis_bins,
    .(
      median_excess =
        median(
          local_excess_4D,
          na.rm = TRUE
        )
    ),
    by = distance_bin
  ]
  
  idx <- match(
    analysis_bins,
    as.character(
      perm_summary$distance_bin
    )
  )
  
  null_mat[b, ] <-
    perm_summary$
    median_excess[idx]
  
  if (b %% 100L == 0L) {
    
    message(
      b,
      " / ",
      N_PERM
    )
  }
}


obs_for_test <- obs_inference[
  match(
    analysis_bins,
    as.character(
      distance_bin
    )
  )
]



matched_results <- data.table(
  distance_bin =
    factor(
      analysis_bins,
      levels = analysis_bins
    ),
  
  observed =
    obs_for_test$
    median_excess,
  
  null_median =
    apply(
      null_mat,
      2,
      median,
      na.rm = TRUE
    ),
  
  null_lo =
    apply(
      null_mat,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
  
  null_hi =
    apply(
      null_mat,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    )
)

matched_results


matched_results[
  ,
  p_upper :=
    vapply(
      seq_len(nrow(.SD)),
      function(j) {
        
        (
          sum(
            null_mat[, j] >=
              observed[j],
            na.rm = TRUE
          ) + 1
        ) /
          (
            sum(
              is.finite(
                null_mat[, j]
              )
            ) + 1
          )
        
      },
      numeric(1)
    )
]




obs_score <- mean(
  matched_results$observed,
  na.rm = TRUE
)

null_score <- rowMeans(
  null_mat,
  na.rm = TRUE
)

p_global <- (
  sum(
    null_score >= obs_score
  ) + 1
) /
  (
    length(null_score) + 1
  )


c(
  observed =
    obs_score,
  
  null_median =
    median(
      null_score,
      na.rm = TRUE
    ),
  
  null_lo =
    quantile(
      null_score,
      0.025,
      na.rm = TRUE
    ),
  
  null_hi =
    quantile(
      null_score,
      0.975,
      na.rm = TRUE
    ),
  
  P =
    p_global
)










focal_4D_percentile <- fourD_pairs_filt[
  focal_for_join,
  on = .(
    colony_id,
    MAG_id,
    Month_t1,
    Month_t2
  ),
  allow.cartesian = TRUE
][
  ,
  .(
    focal_percentile =
      mean(
        abs_delta_4D <
          abs_delta_0D
      ),
    
    n_4D = .N
  ),
  by = focal_id
]



# ============================================================
# FOCAL SWEEP VS LINKED SYNONYMOUS BACKGROUND
# ============================================================

local_4D[
  ,
  neighborhood := fifelse(
    same_gene,
    "Same gene",
    fifelse(
      distance_bp <= 500,
      "Different gene <=0.5 kb",
      fifelse(
        distance_bp <= 1000,
        "Different gene <=1 kb",
        fifelse(
          distance_bp <= 5000,
          "Different gene <=5 kb",
          fifelse(
            distance_bp <= 10000,
            "Different gene <=10 kb",
            NA_character_
          )
        )
      )
    )
  )
]

neighborhood_defs <- list(
  same_gene = function(x) {
    x$same_gene == TRUE
  },
  
  diff_gene_500bp = function(x) {
    x$same_gene == FALSE &
      x$distance_bp <= 500
  },
  
  diff_gene_1kb = function(x) {
    x$same_gene == FALSE &
      x$distance_bp <= 1000
  },
  
  diff_gene_5kb = function(x) {
    x$same_gene == FALSE &
      x$distance_bp <= 5000
  },
  
  diff_gene_10kb = function(x) {
    x$same_gene == FALSE &
      x$distance_bp <= 10000
  }
)




calc_focal_percentile <- function(
    dt,
    keep
) {
  
  x <- dt[
    keep(dt)
  ]
  
  x[
    ,
    .(
      focal_percentile =
        mean(
          abs_delta_4D <
            abs_delta_0D,
          na.rm = TRUE
        ),
      
      focal_abs_delta =
        first(
          abs_delta_0D
        ),
      
      median_4D =
        median(
          abs_delta_4D,
          na.rm = TRUE
        ),
      
      n_4D =
        .N
    ),
    by = focal_id
  ]
}




percentile_same_gene <-
  calc_focal_percentile(
    local_4D,
    neighborhood_defs$same_gene
  )

percentile_500bp <-
  calc_focal_percentile(
    local_4D,
    neighborhood_defs$diff_gene_500bp
  )

percentile_1kb <-
  calc_focal_percentile(
    local_4D,
    neighborhood_defs$diff_gene_1kb
  )

percentile_5kb <-
  calc_focal_percentile(
    local_4D,
    neighborhood_defs$diff_gene_5kb
  )

percentile_10kb <-
  calc_focal_percentile(
    local_4D,
    neighborhood_defs$diff_gene_10kb
  )


percentile_same_gene[
  ,
  context := "Same gene"
]

percentile_500bp[
  ,
  context := "Different gene <=0.5 kb"
]

percentile_1kb[
  ,
  context := "Different gene <=1 kb"
]

percentile_5kb[
  ,
  context := "Different gene <=5 kb"
]

percentile_10kb[
  ,
  context := "Different gene <=10 kb"
]

focal_percentile_spatial <- rbindlist(
  list(
    percentile_same_gene,
    percentile_500bp,
    percentile_1kb,
    percentile_5kb,
    percentile_10kb
  ),
  use.names = TRUE,
  fill = TRUE
)


genome_percentile <- copy(
  focal_4D_percentile
)

genome_percentile[
  ,
  context := "Whole MAG"
]

percentile_all_scales <- rbindlist(
  list(
    focal_percentile_spatial,
    genome_percentile
  ),
  use.names = TRUE,
  fill = TRUE
)


percentile_all_scales[
  ,
  context := factor(
    context,
    levels = c(
      "Same gene",
      "Different gene <=0.5 kb",
      "Different gene <=1 kb",
      "Different gene <=5 kb",
      "Different gene <=10 kb",
      "Whole MAG"
    )
  )
]



percentile_summary <- percentile_all_scales[
  ,
  .(
    n_focal =
      uniqueN(focal_id),
    
    median_percentile =
      median(
        focal_percentile,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        focal_percentile,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        focal_percentile,
        0.75,
        na.rm = TRUE
      ),
    
    prop_gt_90 =
      mean(
        focal_percentile >= 0.90,
        na.rm = TRUE
      ),
    
    prop_gt_95 =
      mean(
        focal_percentile >= 0.95,
        na.rm = TRUE
      ),
    
    prop_gt_99 =
      mean(
        focal_percentile >= 0.99,
        na.rm = TRUE
      )
  ),
  by = context
]

percentile_summary

calc_focal_contrast <- function(
    dt,
    keep
) {
  
  x <- dt[
    keep(dt)
  ]
  
  x[
    ,
    .(
      focal_abs_delta =
        first(
          abs_delta_0D
        ),
      
      median_4D =
        median(
          abs_delta_4D,
          na.rm = TRUE
        ),
      
      contrast =
        first(
          abs_delta_0D
        ) -
        median(
          abs_delta_4D,
          na.rm = TRUE
        ),
      
      n_4D = .N
    ),
    by = focal_id
  ]
}

contrast_same_gene <-
  calc_focal_contrast(
    local_4D,
    neighborhood_defs$same_gene
  )

contrast_500bp <-
  calc_focal_contrast(
    local_4D,
    neighborhood_defs$diff_gene_500bp
  )

contrast_1kb <-
  calc_focal_contrast(
    local_4D,
    neighborhood_defs$diff_gene_1kb
  )

contrast_5kb <-
  calc_focal_contrast(
    local_4D,
    neighborhood_defs$diff_gene_5kb
  )

contrast_10kb <-
  calc_focal_contrast(
    local_4D,
    neighborhood_defs$diff_gene_10kb
  )

contrast_same_gene[, context := "Same gene"]
contrast_500bp[, context := "Different gene <=0.5 kb"]
contrast_1kb[, context := "Different gene <=1 kb"]
contrast_5kb[, context := "Different gene <=5 kb"]
contrast_10kb[, context := "Different gene <=10 kb"]

focal_contrast_spatial <- rbindlist(
  list(
    contrast_same_gene,
    contrast_500bp,
    contrast_1kb,
    contrast_5kb,
    contrast_10kb
  )
)


contrast_summary <- focal_contrast_spatial[
  ,
  .(
    n_focal =
      uniqueN(focal_id),
    
    median_contrast =
      median(
        contrast,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        contrast,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        contrast,
        0.75,
        na.rm = TRUE
      ),
    
    prop_focal_larger =
      mean(
        contrast > 0,
        na.rm = TRUE
      )
  ),
  by = context
]

contrast_summary



# ============================================================
# FIG 6B
# Focal nonsynonymous sweep percentile relative to 4D variation
# ============================================================

fig6B_dt <- percentile_summary[
  ,
  .(
    context,
    n_focal,
    prop_gt_90,
    prop_gt_95
  )
]

fig6B_dt <- melt(
  fig6B_dt,
  id.vars = c(
    "context",
    "n_focal"
  ),
  measure.vars = c(
    "prop_gt_90",
    "prop_gt_95"
  ),
  variable.name = "threshold",
  value.name = "fraction"
)

fig6B_dt[
  ,
  threshold := factor(
    threshold,
    levels = c(
      "prop_gt_90",
      "prop_gt_95"
    ),
    labels = c(
      ">90th percentile",
      ">95th percentile"
    )
  )
]

fig6B_dt[
  ,
  context := factor(
    context,
    levels = c(
      "Same gene",
      "Different gene <=0.5 kb",
      "Different gene <=1 kb",
      "Different gene <=5 kb",
      "Different gene <=10 kb",
      "Whole MAG"
    ),
    labels = c(
      "Same gene",
      "Different gene\n≤0.5 kb",
      "Different gene\n≤1 kb",
      "Different gene\n≤5 kb",
      "Different gene\n≤10 kb",
      "Whole MAG"
    )
  )
]


fig6B_dt[
  ,
  `:=`(
    successes = round(
      fraction * n_focal
    ),
    failures = n_focal -
      round(
        fraction * n_focal
      )
  )
]

fig6B_dt[
  ,
  c("ci_lo", "ci_hi") := {
    
    bt <- binom.test(
      successes,
      n_focal
    )
    
    list(
      bt$conf.int[1],
      bt$conf.int[2]
    )
    
  },
  by = .(
    context,
    threshold
  )
]








contrast_plot_dt <- copy(
  focal_contrast_spatial
)

contrast_plot_dt[
  ,
  context := factor(
    context,
    levels = c(
      "Same gene",
      "Different gene <=0.5 kb",
      "Different gene <=1 kb",
      "Different gene <=5 kb",
      "Different gene <=10 kb"
    ),
    labels = c(
      "Same gene",
      "Different gene\n≤0.5 kb",
      "Different gene\n≤1 kb",
      "Different gene\n≤5 kb",
      "Different gene\n≤10 kb"
    )
  )
]




p6B_contrast <- ggplot(
  contrast_plot_dt,
  aes(
    x = context,
    y = contrast
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  

  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

p6B_contrast







spatial_plot_dt <- rbindlist(
  list(
    copy(
      spatial_same_gene
    ),
    copy(
      spatial_diff_gene
    )
  ),
  use.names = TRUE
)

spatial_plot_dt <- spatial_plot_dt[
  distance_bin %in%
    c(
      "<100 bp",
      "0.1-0.25 kb",
      "0.25-0.5 kb",
      "0.5-1 kb",
      "1-2.5 kb",
      "2.5-5 kb",
      "5-10 kb",
      "10-25 kb"
    )
]

spatial_plot_dt[
  ,
  distance_bin :=
    factor(
      distance_bin,
      levels = c(
        "<100 bp",
        "0.1-0.25 kb",
        "0.25-0.5 kb",
        "0.5-1 kb",
        "1-2.5 kb",
        "2.5-5 kb",
        "5-10 kb",
        "10-25 kb"
      )
    )
]






p6C_spatial_diffgene <- ggplot(
  spatial_diff_gene[
    distance_bin %in%
      c(
        "0.1-0.25 kb",
        "0.25-0.5 kb",
        "0.5-1 kb",
        "1-2.5 kb",
        "2.5-5 kb",
        "5-10 kb",
        "10-25 kb"
      )
  ],
  aes(
    x = distance_bin,
    y = median_excess,
    group = 1
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2.7
  ) +
  
  labs(
    x = paste0(
      "Distance from focal nonsynonymous sweep\n",
      "(different genes on same contig)"
    ),
    y = expression(
      "Excess synonymous " *
        "|" * Delta * "AF|"
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 40,
        hjust = 1
      )
  )

p6C_spatial_diffgene

p6C_spatial




# ============================================================
# FIG 6B CANDIDATE:
# percentile exceedance curves
# ============================================================

pct_plot <- copy(percentile_all_scales)

pct_plot[
  ,
  context := as.character(context)
]

# Main-figure scales only
pct_plot <- pct_plot[
  context %in% c(
    "Same gene",
    "Different gene <=1 kb",
    "Whole MAG"
  )
]

# Require a minimally informative local neighborhood.
# Whole-MAG comparisons already contain many sites.
pct_plot <- pct_plot[
  context == "Whole MAG" |
    n_4D >= 3
]






pct_thresholds <- seq(
  0.50,
  0.99,
  by = 0.01
)

exceedance_curve <- rbindlist(
  lapply(
    pct_thresholds,
    function(q) {
      
      pct_plot[
        ,
        .(
          fraction_exceeding =
            mean(
              focal_percentile >= q,
              na.rm = TRUE
            ),
          
          n_focal =
            uniqueN(focal_id)
        ),
        by = context
      ][
        ,
        threshold := q
      ]
    }
  )
)



null_curve <- data.table(
  threshold = pct_thresholds,
  fraction_exceeding = 1 - pct_thresholds
)


exceedance_curve[
  ,
  context := factor(
    context,
    levels = c(
      "Same gene",
      "Different gene <=1 kb",
      "Whole MAG"
    ),
    labels = c(
      "Same gene",
      "Different gene ≤1 kb",
      "Whole MAG"
    )
  )
]




p6B_exceedance <- ggplot() +
  
  # Random expectation
  geom_line(
    data = null_curve,
    aes(
      x = threshold,
      y = fraction_exceeding
    ),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # Observed focal sweeps
  geom_line(
    data = exceedance_curve,
    aes(
      x = threshold,
      y = fraction_exceeding,
      linetype = context
    ),
    linewidth = 1.15
  ) +
  
  # Highlight familiar thresholds
  geom_vline(
    xintercept = c(0.90, 0.95),
    linetype = "dotted",
    linewidth = 0.45
  ) +
  
  scale_x_continuous(
    limits = c(0.50, 0.99),
    breaks = c(
      0.50,
      0.70,
      0.80,
      0.90,
      0.95,
      0.99
    ),
    labels = scales::percent
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  
  labs(
    x = "Synonymous |ΔAF| percentile threshold",
    y = paste0(
      "Focal nonsynonymous sweeps\n",
      "exceeding threshold"
    ),
    linetype = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "top"
  )

p6B_exceedance




enrichment_dt <- copy(fig6B_dt)

enrichment_dt[
  ,
  expected := fifelse(
    threshold == ">90th percentile",
    0.10,
    0.05
  )
]

enrichment_dt[
  ,
  fold_enrichment :=
    fraction / expected
]

ggplot(
  enrichment_dt,
  aes(
    x = context,
    y = fold_enrichment,
    shape = threshold
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3,
    position = position_dodge(
      width = 0.35
    )
  ) +
  labs(
    x = NULL,
    y = "Enrichment of extreme focal changes",
    shape = NULL
  ) +
  theme_classic()




contrast_forest <- contrast_summary[
  context %in% c(
    "Same gene",
    "Different gene <=0.5 kb",
    "Different gene <=1 kb",
    "Different gene <=5 kb",
    "Different gene <=10 kb"
  )
]

contrast_forest[
  ,
  context := factor(
    context,
    levels = rev(
      c(
        "Same gene",
        "Different gene <=0.5 kb",
        "Different gene <=1 kb",
        "Different gene <=5 kb",
        "Different gene <=10 kb"
      )
    ),
    labels = rev(
      c(
        "Same gene",
        "Different gene ≤0.5 kb",
        "Different gene ≤1 kb",
        "Different gene ≤5 kb",
        "Different gene ≤10 kb"
      )
    )
  )
]



spatial_diff_plot <- copy(
  spatial_diff_gene
)

bin_midpoints_kb <- c(
  "<100 bp" = 0.05,
  "0.1-0.25 kb" = 0.175,
  "0.25-0.5 kb" = 0.375,
  "0.5-1 kb" = 0.75,
  "1-2.5 kb" = 1.75,
  "2.5-5 kb" = 3.75,
  "5-10 kb" = 7.5,
  "10-25 kb" = 17.5
)

spatial_diff_plot[
  ,
  distance_kb :=
    bin_midpoints_kb[
      as.character(
        distance_bin
      )
    ]
]

spatial_diff_plot <- spatial_diff_plot[
  !is.na(distance_kb)
]




p6C_spatial_log <- ggplot(
  spatial_diff_plot,
  aes(
    x = distance_kb,
    y = median_excess
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_line(
    linewidth = 1.1
  ) +
  
  geom_point(
    aes(
      size = n_focal
    )
  ) +
  
  scale_x_log10(
    breaks = c(
      0.1,
      0.25,
      0.5,
      1,
      2.5,
      5,
      10,
      25
    )
  ) +
  
  scale_size_continuous(
    range = c(2, 5)
  ) +
  
  labs(
    x = paste0(
      "Distance from focal nonsynonymous sweep (kb)\n",
      "different genes on same contig"
    ),
    y = expression(
      "Excess synonymous " *
        "|" * Delta * "AF|"
    ),
    size = "Focal sweeps"
  ) +
  
  theme_classic(
    base_size = 12
  )

p6C_spatial_log











pct_thresholds <- seq(
  0.50,
  0.95,
  by = 0.01
)

exceedance_curve <- rbindlist(
  lapply(
    pct_thresholds,
    function(q) {
      
      pct_plot[
        ,
        .(
          fraction_exceeding =
            mean(
              focal_percentile >= q,
              na.rm = TRUE
            ),
          
          n_focal =
            uniqueN(focal_id)
        ),
        by = context
      ][
        ,
        threshold := q
      ]
    }
  )
)

null_curve <- data.table(
  threshold = pct_thresholds,
  fraction_exceeding = 1 - pct_thresholds
)




p6B_exceedance <- ggplot() +
  
  geom_line(
    data = null_curve,
    aes(
      x = threshold,
      y = fraction_exceeding
    ),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_line(
    data = exceedance_curve,
    aes(
      x = threshold,
      y = fraction_exceeding,
      linetype = context
    ),
    linewidth = 1.15
  ) +
  
  geom_vline(
    xintercept = c(
      0.90,
      0.95
    ),
    linetype = "dotted",
    linewidth = 0.45
  ) +
  
  scale_x_continuous(
    limits = c(
      0.50,
      0.95
    ),
    breaks = c(
      0.50,
      0.60,
      0.70,
      0.80,
      0.90,
      0.95
    ),
    labels = scales::percent
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    labels = scales::percent
  ) +
  
  labs(
    x = "Synonymous |ΔAF| percentile threshold",
    y = paste0(
      "Focal nonsynonymous sweeps\n",
      "exceeding threshold"
    ),
    linetype = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "top"
  )

p6B_exceedance

















month_order <- c(
  "May",
  "June",
  "July",
  "August",
  "September",
  "October",
  "November",
  "January",
  "February"
)

transition_lookup <- data.table(
  Month_t1 = head(
    month_order,
    -1
  ),
  Month_t2 = tail(
    month_order,
    -1
  )
)

transition_lookup

LOCAL_RADIUS <- 5000



focal_neighborhood <- merge(
  focal_0D[
    ,
    .(
      focal_id,
      colony_id,
      MAG_id,
      contig_name,
      pos_0D,
      
      gene_0D =
        as.character(
          corresponding_gene_call
        ),
      
      focal_Month_t1 =
        Month_t1,
      
      focal_Month_t2 =
        Month_t2
    )
  ],
  
  unique(
    fourD_clean[
      ,
      .(
        colony_id = Colony,
        MAG_id,
        contig_name,
        pos_4D =
          pos_in_contig,
        
        gene_4D =
          as.character(
            corresponding_gene_call
          )
      )
    ]
  ),
  
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name"
  ),
  
  allow.cartesian = TRUE
)


focal_neighborhood[
  ,
  distance_bp :=
    abs(
      pos_4D - pos_0D
    )
]

focal_neighborhood <- focal_neighborhood[
  distance_bp <= LOCAL_RADIUS &
    !is.na(gene_0D) &
    !is.na(gene_4D) &
    gene_0D != gene_4D
]

neighborhood_support <- focal_neighborhood[
  ,
  .(
    n_local_4D =
      uniqueN(pos_4D)
  ),
  by = focal_id
]

summary(
  neighborhood_support$
    n_local_4D
)

eligible_gold_ids <- neighborhood_support[
  n_local_4D >= 3,
  focal_id
]




fourD_long <- fourD_clean[
  ,
  .(
    colony_id =
      as.character(Colony),
    
    MAG_id =
      as.character(MAG_id),
    
    Month =
      as.character(Month),
    
    contig_name =
      as.character(contig_name),
    
    pos_4D =
      as.integer(pos_in_contig),
    
    AF =
      as.numeric(AF)
  )
]

fourD_all_transitions <- rbindlist(
  lapply(
    seq_len(
      nrow(transition_lookup)
    ),
    function(i) {
      
      m1 <-
        transition_lookup$
        Month_t1[i]
      
      m2 <-
        transition_lookup$
        Month_t2[i]
      
      x1 <- fourD_long[
        Month == m1,
        .(
          colony_id,
          MAG_id,
          contig_name,
          pos_4D,
          AF_t1 = AF
        )
      ]
      
      x2 <- fourD_long[
        Month == m2,
        .(
          colony_id,
          MAG_id,
          contig_name,
          pos_4D,
          AF_t2 = AF
        )
      ]
      
      out <- merge(
        x1,
        x2,
        by = c(
          "colony_id",
          "MAG_id",
          "contig_name",
          "pos_4D"
        ),
        all = FALSE
      )
      
      out[
        ,
        `:=`(
          Month_t1 = m1,
          Month_t2 = m2,
          delta_4D =
            AF_t2 - AF_t1,
          abs_delta_4D =
            abs(
              AF_t2 - AF_t1
            )
        )
      ]
      
      out
    }
  ),
  use.names = TRUE
)




gold_long <- merge(
  focal_neighborhood[
    focal_id %in%
      eligible_gold_ids,
    .(
      focal_id,
      colony_id,
      MAG_id,
      contig_name,
      pos_4D,
      focal_Month_t1,
      focal_Month_t2
    )
  ],
  
  fourD_all_transitions,
  
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_4D"
  ),
  
  all = FALSE,
  allow.cartesian = TRUE
)

gold_long[
  ,
  is_focal_interval :=
    Month_t1 ==
    focal_Month_t1 &
    Month_t2 ==
    focal_Month_t2
]


gold_region_transition <- gold_long[
  ,
  .(
    local_4D_movement =
      median(
        abs_delta_4D,
        na.rm = TRUE
      ),
    
    n_trackable_local_4D =
      .N,
    
    is_focal_interval =
      first(
        is_focal_interval
      )
  ),
  by = .(
    focal_id,
    Month_t1,
    Month_t2
  )
]

MIN_LOCAL_4D <- 3

gold_region_transition <- gold_region_transition[
  n_trackable_local_4D >=
    MIN_LOCAL_4D
]




gold_specificity <- gold_region_transition[
  ,
  {
    
    focal_val <-
      local_4D_movement[
        is_focal_interval
      ]
    
    other_vals <-
      local_4D_movement[
        !is_focal_interval
      ]
    
    if (
      length(focal_val) == 1 &
      length(other_vals) >= 2
    ) {
      
      .(
        focal_local_movement =
          focal_val,
        
        other_median =
          median(
            other_vals,
            na.rm = TRUE
          ),
        
        temporal_excess =
          focal_val -
          median(
            other_vals,
            na.rm = TRUE
          ),
        
        focal_percentile_in_time =
          mean(
            other_vals <
              focal_val,
            na.rm = TRUE
          ),
        
        n_other_intervals =
          length(
            other_vals
          )
      )
      
    } else {
      
      NULL
    }
  },
  by = focal_id
]





summary(
  gold_specificity$
    temporal_excess
)

quantile(
  gold_specificity$
    focal_percentile_in_time,
  c(
    0,
    .25,
    .5,
    .75,
    .9,
    1
  ),
  na.rm = TRUE
)

mean(
  gold_specificity$
    temporal_excess > 0
)

mean(
  gold_specificity$
    focal_percentile_in_time >=
    0.9
)




# ============================================================
# GOLD NULL:
# Is local 4D movement specifically aligned with the
# actual focal 0D sweep interval?
# ============================================================

gold_null_input <- gold_region_transition[
  ,
  .(
    focal_id,
    Month_t1,
    Month_t2,
    local_4D_movement,
    is_focal_interval
  )
]

# Only focal regions with:
#   - one actual focal interval
#   - at least 2 comparison intervals
gold_valid_ids <- gold_null_input[
  ,
  .(
    n_intervals = .N,
    n_focal = sum(is_focal_interval)
  ),
  by = focal_id
][
  n_intervals >= 3 &
    n_focal == 1,
  focal_id
]

gold_null_input <- gold_null_input[
  focal_id %in% gold_valid_ids
]






gold_obs <- gold_null_input[
  ,
  {
    focal_val <-
      local_4D_movement[
        is_focal_interval
      ]
    
    other_vals <-
      local_4D_movement[
        !is_focal_interval
      ]
    
    .(
      temporal_excess =
        focal_val -
        median(
          other_vals,
          na.rm = TRUE
        ),
      
      focal_is_max =
        focal_val >=
        max(
          local_4D_movement,
          na.rm = TRUE
        )
    )
  },
  by = focal_id
]

obs_median_excess <-
  median(
    gold_obs$temporal_excess
  )

obs_prop_positive <-
  mean(
    gold_obs$temporal_excess > 0
  )

obs_prop_max <-
  mean(
    gold_obs$focal_is_max
  )

c(
  median_excess =
    obs_median_excess,
  prop_positive =
    obs_prop_positive,
  prop_max =
    obs_prop_max
)





set.seed(12345)

N_PERM_GOLD <- 1000L

null_median_excess <-
  numeric(N_PERM_GOLD)

null_prop_positive <-
  numeric(N_PERM_GOLD)

null_prop_max <-
  numeric(N_PERM_GOLD)

gold_split <- split(
  gold_null_input,
  by = "focal_id",
  keep.by = TRUE
)

for (b in seq_len(N_PERM_GOLD)) {
  
  perm_stats <- rbindlist(
    lapply(
      gold_split,
      function(x) {
        
        pseudo_idx <-
          sample.int(
            nrow(x),
            1L
          )
        
        pseudo_val <-
          x$local_4D_movement[
            pseudo_idx
          ]
        
        other_vals <-
          x$local_4D_movement[
            -pseudo_idx
          ]
        
        data.table(
          temporal_excess =
            pseudo_val -
            median(
              other_vals,
              na.rm = TRUE
            ),
          
          pseudo_is_max =
            pseudo_val >=
            max(
              x$local_4D_movement,
              na.rm = TRUE
            )
        )
      }
    )
  )
  
  null_median_excess[b] <-
    median(
      perm_stats$temporal_excess
    )
  
  null_prop_positive[b] <-
    mean(
      perm_stats$temporal_excess > 0
    )
  
  null_prop_max[b] <-
    mean(
      perm_stats$pseudo_is_max
    )
}




p_excess <- (
  sum(
    null_median_excess >=
      obs_median_excess
  ) + 1
) /
  (
    N_PERM_GOLD + 1
  )

p_positive <- (
  sum(
    null_prop_positive >=
      obs_prop_positive
  ) + 1
) /
  (
    N_PERM_GOLD + 1
  )

p_max <- (
  sum(
    null_prop_max >=
      obs_prop_max
  ) + 1
) /
  (
    N_PERM_GOLD + 1
  )

c(
  observed_median_excess =
    obs_median_excess,
  
  null_median_excess =
    median(
      null_median_excess
    ),
  
  null_excess_lo =
    quantile(
      null_median_excess,
      0.025
    ),
  
  null_excess_hi =
    quantile(
      null_median_excess,
      0.975
    ),
  
  P_excess =
    p_excess,
  
  observed_prop_positive =
    obs_prop_positive,
  
  null_prop_positive =
    median(
      null_prop_positive
    ),
  
  P_positive =
    p_positive,
  
  observed_prop_max =
    obs_prop_max,
  
  null_prop_max =
    median(
      null_prop_max
    ),
  
  P_max =
    p_max
)




# ============================================================
# FIGURE 6C
# Event-centered temporal localization
# ============================================================

transition_lookup <- data.table(
  Month_t1 = c(
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "January"
  ),
  Month_t2 = c(
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "January",
    "February"
  ),
  transition_index = 1:8
)

# attach chronological index to every measured regional transition
gold_event_time <- merge(
  gold_region_transition,
  transition_lookup,
  by = c(
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)

# index of the actual focal sweep transition
focal_index <- merge(
  focal_0D[
    ,
    .(
      focal_id,
      Month_t1,
      Month_t2
    )
  ],
  transition_lookup,
  by = c(
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)[
  ,
  .(
    focal_id,
    focal_transition_index =
      transition_index
  )
]

gold_event_time <- merge(
  gold_event_time,
  focal_index,
  by = "focal_id",
  all.x = TRUE
)

gold_event_time[
  ,
  relative_transition :=
    transition_index -
    focal_transition_index
]


gold_event_time[
  ,
  region_temporal_baseline :=
    median(
      local_4D_movement[
        relative_transition != 0
      ],
      na.rm = TRUE
    ),
  by = focal_id
]

gold_event_time[
  ,
  temporal_excess :=
    local_4D_movement -
    region_temporal_baseline
]


gold_event_summary <- gold_event_time[
  relative_transition >= -3 &
    relative_transition <= 3,
  .(
    median_excess =
      median(
        temporal_excess,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        temporal_excess,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        temporal_excess,
        0.75,
        na.rm = TRUE
      ),
    
    n_focal =
      uniqueN(focal_id)
  ),
  by = relative_transition
][
  order(relative_transition)
]

gold_event_summary





p6C_temporal <- ggplot(
  gold_event_summary,
  aes(
    x = relative_transition,
    y = median_excess
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75
    ),
    alpha = 0.18
  ) +
  
  geom_line(
    linewidth = 1.1
  ) +
  
  geom_point(
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = -3:3,
    labels = c(
      "-3",
      "-2",
      "-1",
      "Sweep",
      "+1",
      "+2",
      "+3"
    )
  ) +
  
  labs(
    x = "Transition relative to focal nonsynonymous sweep",
    y = expression(
      "Local synonymous " *
        "|" * Delta * "AF|" *
        " relative to temporal baseline"
    )
  ) +
  
  theme_classic(
    base_size = 12
  )

p6C_temporal





p6C_temporal <- p6C_temporal +
  
  annotate(
    "text",
    x = 2.9,
    y = Inf,
    hjust = 1,
    vjust = 1.4,
    label = paste0(
      "Median excess = ",
      sprintf(
        "%.3f",
        obs_median_excess
      ),
      "\n",
      "Sweep interval most dynamic: ",
      round(
        100 * obs_prop_max
      ),
      "% vs ",
      round(
        100 *
          median(
            null_prop_max
          )
      ),
      "% expected",
      "\nP < 0.001"
    ),
    size = 3.2
  )

p6C_temporal






gold_null_plot <- data.table(
  null_median_excess =
    null_median_excess
)

ggplot(
  gold_null_plot,
  aes(
    x = null_median_excess
  )
) +
  
  geom_histogram(
    bins = 35
  ) +
  
  geom_vline(
    xintercept =
      obs_median_excess,
    linewidth = 1.1
  ) +
  
  annotate(
    "text",
    x = obs_median_excess,
    y = Inf,
    label = paste0(
      "Observed = ",
      sprintf(
        "%.3f",
        obs_median_excess
      )
    ),
    hjust = 1.1,
    vjust = 1.5
  ) +
  
  labs(
    x = paste0(
      "Median temporal excess in local\n",
      "synonymous allele-frequency change"
    ),
    y = "Permutations"
  ) +
  
  theme_classic(
    base_size = 12
  )



gold_event_summary_plot <- gold_event_summary[
  relative_transition >= -2 &
    relative_transition <= 2
]

p6C_temporal_clean <- ggplot(
  gold_event_summary_plot,
  aes(
    x = relative_transition,
    y = median_excess
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = q25,
      ymax = q75
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  
  geom_line(
    linewidth = 1.05
  ) +
  
  geom_point(
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = -2:2,
    labels = c(
      "-2",
      "-1",
      "Sweep",
      "+1",
      "+2"
    )
  ) +
  
  labs(
    x = "Transition relative to focal nonsynonymous sweep",
    y = expression(
      "Local synonymous " *
        "|" * Delta * "AF|" *
        " relative to temporal baseline"
    )
  ) +
  
  theme_classic(
    base_size = 12
  )

p6C_temporal_clean













#SILVER BABY


# ============================================================
# SILVER
# Full-timecourse trajectory coherence
# ============================================================

focal_0D_long <- merge(
  focal_0D[
    ,
    .(
      focal_id,
      colony_id,
      MAG_id,
      contig_name,
      pos_0D
    )
  ],
  
  zeroD_clean_all[
    ,
    .(
      colony_id = Colony,
      MAG_id,
      Month,
      contig_name,
      pos_0D =
        pos_in_contig,
      
      focal_AF =
        AF
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_0D"
  ),
  
  all = FALSE
)




silver_pairs <- merge(
  focal_neighborhood[
    ,
    .(
      focal_id,
      colony_id,
      MAG_id,
      contig_name,
      pos_4D,
      distance_bp
    )
  ],
  
  fourD_clean[
    ,
    .(
      colony_id =
        as.character(Colony),
      
      MAG_id =
        as.character(MAG_id),
      
      Month =
        as.character(Month),
      
      contig_name =
        as.character(contig_name),
      
      pos_4D =
        as.integer(pos_in_contig),
      
      local_4D_AF =
        AF
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_4D"
  ),
  
  allow.cartesian = TRUE
)

silver_pairs <- merge(
  silver_pairs,
  
  focal_0D_long[
    ,
    .(
      focal_id,
      Month,
      focal_AF
    )
  ],
  
  by = c(
    "focal_id",
    "Month"
  ),
  
  all = FALSE
)


silver_cor <- silver_pairs[
  ,
  {
    
    ok <-
      is.finite(focal_AF) &
      is.finite(local_4D_AF)
    
    if (sum(ok) >= 4) {
      
      r <- cor(
        focal_AF[ok],
        local_4D_AF[ok],
        method = "spearman"
      )
      
      .(
        n_months =
          sum(ok),
        
        rho = r,
        
        abs_rho =
          abs(r),
        
        distance_bp =
          first(
            distance_bp
          )
      )
      
    } else {
      
      NULL
    }
  },
  by = .(
    focal_id,
    pos_4D
  )
]



silver_cor[
  ,
  distance_bin := cut(
    distance_bp,
    breaks = c(
      0,
      500,
      1000,
      2500,
      5000
    ),
    labels = c(
      "<0.5 kb",
      "0.5-1 kb",
      "1-2.5 kb",
      "2.5-5 kb"
    ),
    right = FALSE
  )
]



silver_by_focal <- silver_cor[
  !is.na(distance_bin),
  .(
    median_abs_rho =
      median(
        abs_rho,
        na.rm = TRUE
      ),
    
    median_signed_rho =
      median(
        rho,
        na.rm = TRUE
      ),
    
    n_local_4D =
      .N
  ),
  by = .(
    focal_id,
    distance_bin
  )
]



silver_summary <- silver_by_focal[
  ,
  .(
    median_abs_rho =
      median(
        median_abs_rho,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        median_abs_rho,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        median_abs_rho,
        0.75,
        na.rm = TRUE
      ),
    
    median_signed_rho =
      median(
        median_signed_rho,
        na.rm = TRUE
      ),
    
    n_focal =
      uniqueN(focal_id)
  ),
  by = distance_bin
]

silver_summary





# local result already exists:
# silver_cor: focal_id, pos_4D, rho, abs_rho, distance_bp

silver_local_by_focal <- silver_cor[
  distance_bp <= 5000,
  .(
    local_abs_rho = median(abs_rho, na.rm = TRUE),
    local_signed_rho = median(rho, na.rm = TRUE)
  ),
  by = focal_id
]
N_BG_PER_FOCAL <- 200L



silver_compare <- merge(
  silver_local_by_focal,
  silver_bg_by_focal,
  by = "focal_id"
)

silver_compare[
  ,
  delta_abs_rho :=
    local_abs_rho - bg_abs_rho
]

summary(silver_compare$delta_abs_rho)
mean(silver_compare$delta_abs_rho > 0)





N_BG_PER_FOCAL <- 200L
set.seed(12345)

# focal metadata
silver_focal_meta <- unique(
  focal_0D[
    ,
    .(
      focal_id,
      colony_id,
      MAG_id,
      contig_name,
      pos_0D
    )
  ]
)

# candidate background 4D physical sites
fourD_sites_bg <- unique(
  fourD_clean[
    ,
    .(
      colony_id = as.character(Colony),
      MAG_id = as.character(MAG_id),
      contig_name = as.character(contig_name),
      pos_4D = as.integer(pos_in_contig)
    )
  ]
)

# join within same colony x MAG
silver_bg_candidates <- merge(
  silver_focal_meta,
  fourD_sites_bg,
  by = c(
    "colony_id",
    "MAG_id"
  ),
  allow.cartesian = TRUE
)

# exclude local <=5 kb sites on same contig
silver_bg_candidates[
  ,
  distance_bp := fifelse(
    contig_name.x == contig_name.y,
    abs(pos_4D - pos_0D),
    Inf
  )
]

silver_bg_candidates <- silver_bg_candidates[
  contig_name.x != contig_name.y |
    distance_bp > 5000
]

setnames(
  silver_bg_candidates,
  c("contig_name.x", "contig_name.y"),
  c("focal_contig", "contig_name")
)

silver_bg_sample <- silver_bg_candidates[
  ,
  .SD[
    sample(
      .N,
      min(.N, N_BG_PER_FOCAL)
    )
  ],
  by = focal_id
]




silver_bg_long <- merge(
  silver_bg_sample[
    ,
    .(
      focal_id,
      colony_id,
      MAG_id,
      contig_name,
      pos_4D
    )
  ],
  
  fourD_clean[
    ,
    .(
      colony_id = as.character(Colony),
      MAG_id = as.character(MAG_id),
      Month = as.character(Month),
      contig_name = as.character(contig_name),
      pos_4D = as.integer(pos_in_contig),
      bg_4D_AF = AF
    )
  ],
  
  by = c(
    "colony_id",
    "MAG_id",
    "contig_name",
    "pos_4D"
  ),
  allow.cartesian = TRUE
)



silver_bg_long <- merge(
  silver_bg_long,
  
  focal_0D_long[
    ,
    .(
      focal_id,
      Month,
      focal_AF
    )
  ],
  
  by = c(
    "focal_id",
    "Month"
  ),
  
  all = FALSE
)



silver_bg_cor <- silver_bg_long[
  ,
  {
    
    ok <-
      is.finite(focal_AF) &
      is.finite(bg_4D_AF)
    
    if (sum(ok) >= 4) {
      
      r <- cor(
        focal_AF[ok],
        bg_4D_AF[ok],
        method = "spearman"
      )
      
      .(
        n_months = sum(ok),
        rho = r,
        abs_rho = abs(r)
      )
      
    } else {
      
      NULL
    }
  },
  by = .(
    focal_id,
    contig_name,
    pos_4D
  )
]



silver_local_by_focal <- silver_cor[
  distance_bp <= 5000,
  .(
    local_abs_rho =
      median(abs_rho, na.rm = TRUE),
    
    local_signed_rho =
      median(rho, na.rm = TRUE)
  ),
  by = focal_id
]

silver_bg_by_focal <- silver_bg_cor[
  ,
  .(
    bg_abs_rho =
      median(abs_rho, na.rm = TRUE),
    
    bg_signed_rho =
      median(rho, na.rm = TRUE)
  ),
  by = focal_id
]



silver_compare <- merge(
  silver_local_by_focal,
  silver_bg_by_focal,
  by = "focal_id"
)

silver_compare[
  ,
  `:=`(
    delta_abs_rho =
      local_abs_rho - bg_abs_rho,
    
    delta_signed_rho =
      local_signed_rho - bg_signed_rho
  )
]

summary(
  silver_compare$delta_abs_rho
)

mean(
  silver_compare$delta_abs_rho > 0
)


wilcox.test(
  silver_compare$local_abs_rho,
  silver_compare$bg_abs_rho,
  paired = TRUE,
  alternative = "greater"
)
median(
  silver_compare$local_abs_rho
)





c(
  n_focal = nrow(silver_compare),
  
  local_median =
    median(
      silver_compare$local_abs_rho
    ),
  
  bg_median =
    median(
      silver_compare$bg_abs_rho
    ),
  
  median_delta =
    median(
      silver_compare$delta_abs_rho
    ),
  
  prop_local_greater =
    mean(
      silver_compare$delta_abs_rho > 0
    )
)
median(
  silver_compare$bg_abs_rho
)












# ============================================================
# BRONZE
# Breadth of linked synonymous participation
# ============================================================

calc_bronze <- function(dt, radius_bp, min_sites = 3) {
  
  out <- dt[
    same_gene == FALSE &
      distance_bp <= radius_bp,
    .(
      focal_abs_delta =
        first(abs_delta_0D),
      
      focal_delta =
        first(delta_0D),
      
      local_median_abs_delta =
        median(
          abs_delta_4D,
          na.rm = TRUE
        ),
      
      linked_ratio =
        median(
          abs_delta_4D,
          na.rm = TRUE
        ) /
        first(abs_delta_0D),
      
      # Any local 4D site moving at least half
      # as strongly as the focal sweep
      fraction_half_focal =
        mean(
          abs_delta_4D >=
            0.5 * first(abs_delta_0D),
          na.rm = TRUE
        ),
      
      # Same, but also moving in the same direction
      fraction_half_focal_concordant =
        mean(
          abs_delta_4D >=
            0.5 * first(abs_delta_0D) &
            sign(delta_4D) ==
            sign(first(delta_0D)),
          na.rm = TRUE
        ),
      
      # Stricter: nearly as large as focal
      fraction_80pct_focal_concordant =
        mean(
          abs_delta_4D >=
            0.8 * first(abs_delta_0D) &
            sign(delta_4D) ==
            sign(first(delta_0D)),
          na.rm = TRUE
        ),
      
      n_local_4D = .N
    ),
    by = focal_id
  ]
  
  out <- out[
    n_local_4D >= min_sites
  ]
  
  out[
    ,
    radius_kb := radius_bp / 1000
  ]
  
  out
}


bronze_05kb <- calc_bronze(
  local_4D,
  radius_bp = 500
)

bronze_1kb <- calc_bronze(
  local_4D,
  radius_bp = 1000
)

bronze_5kb <- calc_bronze(
  local_4D,
  radius_bp = 5000
)

bronze_10kb <- calc_bronze(
  local_4D,
  radius_bp = 10000
)

bronze_all <- rbindlist(
  list(
    bronze_05kb,
    bronze_1kb,
    bronze_5kb,
    bronze_10kb
  ),
  use.names = TRUE
)




bronze_summary <- bronze_all[
  ,
  .(
    n_focal =
      uniqueN(focal_id),
    
    median_linked_ratio =
      median(
        linked_ratio,
        na.rm = TRUE
      ),
    
    q25_linked_ratio =
      quantile(
        linked_ratio,
        0.25,
        na.rm = TRUE
      ),
    
    q75_linked_ratio =
      quantile(
        linked_ratio,
        0.75,
        na.rm = TRUE
      ),
    
    median_fraction_half =
      median(
        fraction_half_focal,
        na.rm = TRUE
      ),
    
    median_fraction_concordant =
      median(
        fraction_half_focal_concordant,
        na.rm = TRUE
      ),
    
    median_fraction_80pct_concordant =
      median(
        fraction_80pct_focal_concordant,
        na.rm = TRUE
      ),
    
    prop_majority_concordant =
      mean(
        fraction_half_focal_concordant >= 0.5,
        na.rm = TRUE
      )
  ),
  by = radius_kb
][
  order(radius_kb)
]

bronze_summary





summary(
  bronze_1kb$linked_ratio
)

summary(
  bronze_1kb$fraction_half_focal_concordant
)

summary(
  bronze_5kb$linked_ratio
)

summary(
  bronze_5kb$fraction_half_focal_concordant
)



data.table(
  metric = c(
    "1 kb: median 4D >= 50% focal",
    "1 kb: majority of 4D strongly concordant",
    "5 kb: median 4D >= 50% focal",
    "5 kb: majority of 4D strongly concordant"
  ),
  
  fraction = c(
    mean(
      bronze_1kb$linked_ratio >= 0.5,
      na.rm = TRUE
    ),
    
    mean(
      bronze_1kb$
        fraction_half_focal_concordant >= 0.5,
      na.rm = TRUE
    ),
    
    mean(
      bronze_5kb$linked_ratio >= 0.5,
      na.rm = TRUE
    ),
    
    mean(
      bronze_5kb$
        fraction_half_focal_concordant >= 0.5,
      na.rm = TRUE
    )
  )
)



ggplot(
  bronze_1kb,
  aes(
    x = linked_ratio,
    y = fraction_half_focal_concordant
  )
) +
  geom_point(
    alpha = 0.45
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed"
  ) +
  labs(
    x = paste0(
      "Median local synonymous |ΔAF| / ",
      "focal nonsynonymous |ΔAF|"
    ),
    y = paste0(
      "Fraction of local synonymous sites\n",
      "with ≥50% focal change, same direction"
    )
  ) +
  theme_classic()






bronze_plot_dt <- rbindlist(
  list(
    copy(bronze_1kb)[, radius := "≤1 kb"],
    copy(bronze_5kb)[, radius := "≤5 kb"]
  )
)

bronze_plot_dt[
  ,
  radius := factor(
    radius,
    levels = c(
      "≤1 kb",
      "≤5 kb"
    )
  )
]

p_bronze <- ggplot(
  bronze_plot_dt,
  aes(
    x = linked_ratio,
    linetype = radius
  )
) +
  
  geom_density(
    linewidth = 1.1,
    adjust = 1.1
  ) +
  
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  
  coord_cartesian(
    xlim = c(0, 1.5)
  ) +
  
  labs(
    x = paste0(
      "Median local synonymous |ΔAF| / ",
      "focal nonsynonymous |ΔAF|"
    ),
    y = "Density",
    linetype = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "top"
  )

p_bronze



bronze_1kb[
  ,
  architecture := fifelse(
    linked_ratio >= 0.8,
    "Broad linked",
    fifelse(
      linked_ratio >= 0.3,
      "Intermediate",
      "Localized"
    )
  )
][
  ,
  .N,
  by = architecture
][
  ,
  fraction := N / sum(N)
]




resident_key <- unique(
  evo_sweep[
    ,
    .(
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      resident_preserving
    )
  ]
)

local_4D_resident <- merge(
  local_4D,
  resident_key,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)

spatial_by_background <- local_4D_resident[
  same_gene == FALSE &
    !is.na(distance_bin) &
    !is.na(resident_preserving),
  .(
    focal_bin_excess =
      median(
        excess_abs_delta_4D,
        na.rm = TRUE
      )
  ),
  by = .(
    focal_id,
    resident_preserving,
    distance_bin
  )
][
  ,
  .(
    median_excess =
      median(
        focal_bin_excess,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        focal_bin_excess,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        focal_bin_excess,
        0.75,
        na.rm = TRUE
      ),
    
    n_focal =
      uniqueN(focal_id)
  ),
  by = .(
    resident_preserving,
    distance_bin
  )
]


spatial_by_background










run_parallel_null <- function(
    sweep_events,
    degeneracy_keep,
    B = 1000,
    seed = 1
) {
  
  x <- sweep_events[
    degeneracy == degeneracy_keep
  ]
  
  # -------------------------------------------------
  # Everything below should be IDENTICAL to your
  # existing pair-score/null machinery:
  #
  # 1. map genes -> MMseqs homologous family
  # 2. define independent cross-colony/cross-genus
  #    concurrent sweep pairs
  # 3. calculate observed score
  # 4. randomize according to existing null
  # 5. return null score vector
  # -------------------------------------------------
  
  list(
    observed = observed_score,
    null = null_scores
  )
}


