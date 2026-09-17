library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)

pm_dir <- paste0(
  "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/",
  "ALLSPECIES_intrapop_diversity/",
  "private_marker_results_resident_oriented"
)

pm_sum <- fread(
  file.path(
    pm_dir,
    "private_marker_summary_all_colonies.csv"
  )
)

pm_month <- fread(
  file.path(
    pm_dir,
    "private_marker_monthly_trends_all_colonies.csv"
  )
)

pm_trans <- fread(
  file.path(
    pm_dir,
    "private_marker_transition_summary_all_colonies.csv"
  )
)

pm_snv <- fread(
  file.path(
    pm_dir,
    "private_marker_SNV_trajectories_all_colonies.csv"
  )
)


for (x in list(pm_sum, pm_trans, pm_snv)) {
  x[, Colony := as.character(Colony)]
  x[, MAG_id := as.character(MAG_id)]
}


#How many private markers support each SNV trajectory?

summary(pm_sum$Total_private)
quantile(
  pm_sum$Total_private,
  probs = c(
    0, 0.05, 0.10, 0.25,
    0.50, 0.75, 0.90, 0.95, 1
  ),
  na.rm = TRUE
)

p_private_n <- ggplot(
  pm_sum,
  aes(x = Total_private)
) +
  geom_histogram(bins = 50)+
  scale_x_log10() +
  theme_classic()


#next, we'll check to for sources of turnover in snv trajectories

summary(pm_month$callable_fraction)
quantile(
  pm_month$callable_fraction,
  probs = c(
    0, 0.05, 0.10, 0.25,
    0.50, 0.75, 0.90, 0.95, 1
  ),
  na.rm = TRUE
)


p_callable <- ggplot(
  pm_month,
  aes(x = callable_fraction)
) +
  geom_histogram(bins = 50)+
  scale_x_log10() +
  theme_classic()
p_callable

#p_callable looks good, median near 0.986 

#now we'll look at the resident-background turnover distribution to garner whether:
#we see a mass approx 0 and a mass approx 1
#a nice continuous distribution
#a bunch of intermediates

summary(pm_sum$max_background_turnover)
quantile(
  pm_sum$max_background_turnover,
  probs = c(
    0, 0.05, 0.10, 0.25,
    0.50, 0.75, 0.90, 0.95, 0.99, 1
  ),
  na.rm = TRUE
)


p_turnover <- ggplot(
  pm_sum,
  aes(x = max_background_turnover)
)+
  geom_histogram(
    bins=50
  )+
  theme_classic()
p_turnover

#looks to be mostly intermediate values, median near ~0.3

p_compare <- ggplot(
  pm_sum,
  aes(
    x=Frac_preserved,
    y=max_background_turnover
  ))+
    geom_point(
      aes(size= Total_private),
      alpha = 0.5
    )+
    scale_size_continuous(
      trans = "log10"
    )+
    theme_classic()
p_compare

#interesting distribution, negative correlation... helps us differentiate between edge cases
# where a MAG with different private alleles independently fluctuate throughout the year
#These alleles could be labeled as disrupted but if they never move together, then we never see lineage replacement


summary(pm_trans$delta_background_turnover)
quantile(
  pm_trans$delta_background_turnover,
  probs = c(
    0, 0.01, 0.05, 0.10,
    0.25, 0.50, 0.75, 0.90,
    0.95, 0.99, 1
  ),
  na.rm = TRUE
)


pm_trans[
  ,
  abs_delta_background :=
    abs(delta_background_turnover)
]

p_delta <- ggplot(
  pm_trans,
  aes(x=delta_background_turnover)
)+
  geom_histogram(bins = 50)+
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  )+
  theme_classic()
p_delta  

#this figure (above ) shows the frequency of transitions rather than entire snv trajectories


pm_trans_good <- pm_trans[
  n_private_baseline >= 20 &
    n_private_callable_t1 >= 20 &
    n_private_callable_t2 >= 20 &
    min_callable_fraction_pair >= 0.8 &
    is.finite(background_turnover_t1) &
    is.finite(background_turnover_t2)
]

nrow(pm_trans)
nrow(pm_trans_good)

uniqueN(
  pm_trans_good[
    ,
    paste(
      MAG_id,
      Colony,
      sep = "::"
    )
  ]
)


#sensitivity analysis to ask "how many snvs do we need to detect these transitions?

marker_thresholds <- c(
  10, 20, 50, 100
)

threshold_summary <- rbindlist(
  lapply(
    marker_thresholds,
    function(k) {
      
      x <- pm_trans[
        n_private_baseline >= k &
          min_callable_fraction_pair >= 0.8
      ]
      
      data.table(
        min_private = k,
        n_transitions = nrow(x),
        n_MAG_colony = uniqueN(
          paste(
            x$MAG_id,
            x$Colony
          )
        ),
        median_turnover_t2 =
          median(
            x$background_turnover_t2,
            na.rm = TRUE
          ),
        median_abs_delta =
          median(
            abs(
              x$delta_background_turnover
            ),
            na.rm = TRUE
          )
      )
    }
  )
)

threshold_summary


#revisiting known cases
focus_mags <- c(
  "mmag_254",
  "mmag_191",
  "mmag_411",
  "mmag_648",
  "mmag_283",
  "mmag_139",
  "mmag_385",
  "mmag_688"
)

pm_sum[
  MAG_id %chin% focus_mags,
  .(
    Colony,
    MAG_id,
    Total_private,
    Frac_preserved,
    Frac_disrupted,
    max_background_turnover,
    max_frac_lost,
    median_monthly_callability
  )
][
  order(MAG_id, Colony)
]



pm_month_focus <- pm_month[
  MAG_id %chin% focus_mags
]

p_focus <- ggplot(
  pm_month_focus,
  aes(
    x = month_index,
    y = background_turnover_median,
    group = interaction(
      MAG_id,
      Colony
    )
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(
    ~ MAG_id + Colony,
    scales = "free_x"
  ) +
  scale_x_continuous(
    breaks = 0:8,
    labels = c(
      "May",
      "Jun",
      "Jul",
      "Aug",
      "Sep",
      "Oct",
      "Nov",
      "Jan",
      "Feb"
    )
  ) +
  labs(
    x = NULL,
    y = "Resident-background turnover"
  ) +
  theme_classic() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

print(p_focus)
#very promising! This new metric very well recapitulates what we observed using our existing dataset


#large takeaways:
#callability is very good, ensuring that resident-background turnover is not because 
#large fractions of the original marker became unobservable

#turnover distribution is very continuous, suggesting lots of mixed populations -- as expected
#correlation between frac_preserved and max_background_turnover are due in part
#that they derive from the same measure, but they capture different aspects of the data

"""
Now we build our second axis: genomic breadth of change 
This basically asks: how much of the genome exhibits turnover?

For every MAG x Colony x pairwise timepoint (t_1 -> t_2)
we ask:
  B_ict = resident-background turnover
  
and independently check:
  G_ict = number of genes containing significant AF changes/number opf genes eligable/callable
"""


SWEEP_FILE <- paste0(
  "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/",
  "ALLSPECIES_intrapop_diversity/",
  "sweep_analysis/sweep_gene_merged.csv"
)






SWEEP_FILE <- paste0(
  "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/",
  "ALLSPECIES_intrapop_diversity/",
  "sweep_analysis/sweep_gene_merged.csv"
)

sweep_dt <- fread(SWEEP_FILE)

sweep_dt[, MAG_id := as.character(MAG_id)]
sweep_dt[, colony_id := sub("\\.0$", "", as.character(colony_id))]
sweep_dt[, Month_t1 := tools::toTitleCase(as.character(Month_t1))]
sweep_dt[, Month_t2 := tools::toTitleCase(as.character(Month_t2))]
sweep_dt[, corresponding_gene_call := as.character(corresponding_gene_call)]

sweep_dt[, gene_key := paste(
  MAG_id,
  corresponding_gene_call,
  sep = "::"
)]

sweep_breadth_num <- sweep_dt[
  !is.na(Month_t1) &
    !is.na(Month_t2),
  .(
    n_sweep_SNVS = .N,
    
    n_swept_genes =
      uniqueN(corresponding_gene_call),
    
    n_0D_swept_genes =
      uniqueN(
        corresponding_gene_call[
          degeneracy == "0D"
        ]
      ),
    
    n_4D_swept_genes =
      uniqueN(
        corresponding_gene_call[
          degeneracy == "4D"
        ]
      )
  ),
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2
  )
]

sweep_breadth_num


focus_mags <- c(
  "mmag_254",
  "mmag_191",
  "mmag_411",
  "mmag_648",
  "mmag_283",
  "mmag_139"
)

sweep_breadth_num[
  MAG_id %chin% focus_mags
][
  order(MAG_id, colony_id, Month_t1)
]









# ============================================================
# Paths
# ============================================================

TRAJ_DIR <- paste0(
  "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/",
  "ALLSPECIES_intrapop_diversity/",
  "snv_frequency_trajectories"
)

# pm_trans_good should already exist from the private-marker analysis
# Required columns:
# MAG_id, Colony, Month_t1, Month_t2


# ============================================================
# Standardize transition table
# ============================================================

pairs <- unique(
  pm_trans_good[
    ,
    .(
      MAG_id = as.character(MAG_id),
      colony_id = as.character(Colony),
      Month_t1 = tools::toTitleCase(as.character(Month_t1)),
      Month_t2 = tools::toTitleCase(as.character(Month_t2))
    )
  ]
)

pairs[
  ,
  pair_id := .I
]

cat(
  "Transitions to analyze:",
  nrow(pairs),
  "\n"
)


# ============================================================
# Process one colony at a time
# ============================================================

denominator_list <- vector(
  "list",
  length = uniqueN(pairs$colony_id)
)

colonies <- sort(
  unique(pairs$colony_id)
)


for (ii in seq_along(colonies)) {
  
  colony_i <- colonies[ii]
  
  cat(
    "\nProcessing colony",
    colony_i,
    "...\n"
  )
  
  
  # ----------------------------------------------------------
  # Find relevant transition pairs
  # ----------------------------------------------------------
  
  pairs_i <- pairs[
    colony_id == colony_i
  ]
  
  mags_i <- unique(
    pairs_i$MAG_id
  )
  
  
  # ----------------------------------------------------------
  # Load colony trajectory file
  # ----------------------------------------------------------
  
  f <- file.path(
    TRAJ_DIR,
    paste0(
      "snv_frequency_trajectory_",
      colony_i,
      ".csv"
    )
  )
  
  if (!file.exists(f)) {
    
    warning(
      "Missing trajectory file: ",
      f
    )
    
    next
  }
  
  
  traj <- fread(
    f,
    select = c(
      "MAG_id",
      "Month",
      "contig_name",
      "pos_in_contig",
      "departure_from_polarized_consensus"
    )
  )
  
  
  traj[
    ,
    MAG_id :=
      as.character(MAG_id)
  ]
  
  traj[
    ,
    Month :=
      tools::toTitleCase(
        as.character(Month)
      )
  ]
  
  
  # Only retain MAGs involved in transitions we care about
  traj <- traj[
    MAG_id %chin% mags_i
  ]
  
  
  # ----------------------------------------------------------
  # Define physical SNV identity
  # ----------------------------------------------------------
  
  traj[
    ,
    SNV_ID := paste(
      MAG_id,
      contig_name,
      pos_in_contig,
      sep = "::"
    )
  ]
  
  
  # An SNV counts as observed at a month only if an AF estimate exists.
  traj <- traj[
    !is.na(
      departure_from_polarized_consensus
    )
  ]
  
  
  # We need only presence/absence of each SNV at each month.
  snv_month <- unique(
    traj[
      ,
      .(
        MAG_id,
        Month,
        SNV_ID
      )
    ]
  )
  
  
  # ----------------------------------------------------------
  # Count overlap for each exact transition
  # ----------------------------------------------------------
  
  out_i <- vector(
    "list",
    nrow(pairs_i)
  )
  
  
  for (jj in seq_len(nrow(pairs_i))) {
    
    p <- pairs_i[jj]
    
    
    snv_t1 <- snv_month[
      MAG_id == p$MAG_id &
        Month == p$Month_t1,
      SNV_ID
    ]
    
    
    snv_t2 <- snv_month[
      MAG_id == p$MAG_id &
        Month == p$Month_t2,
      SNV_ID
    ]
    
    
    n_t1 <- length(
      unique(snv_t1)
    )
    
    n_t2 <- length(
      unique(snv_t2)
    )
    
    n_both <- length(
      intersect(
        snv_t1,
        snv_t2
      )
    )
    
    
    # Useful QC only:
    #
    # 1 means essentially every SNV from the smaller endpoint
    # was also observable at the other endpoint.
    overlap_fraction_min <- if (
      min(n_t1, n_t2) > 0
    ) {
      
      n_both /
        min(n_t1, n_t2)
      
    } else {
      
      NA_real_
    }
    
    
    out_i[[jj]] <- data.table(
      
      pair_id =
        p$pair_id,
      
      MAG_id =
        p$MAG_id,
      
      colony_id =
        p$colony_id,
      
      Month_t1 =
        p$Month_t1,
      
      Month_t2 =
        p$Month_t2,
      
      n_SNVS_t1 =
        n_t1,
      
      n_SNVS_t2 =
        n_t2,
      
      # THIS is our proposed denominator
      n_trackable_SNVS =
        n_both,
      
      overlap_fraction_min =
        overlap_fraction_min
    )
  }
  
  
  denominator_list[[ii]] <- rbindlist(
    out_i
  )
  
  
  rm(
    traj,
    snv_month,
    out_i
  )
  
  gc()
}


# ============================================================
# Combine colonies
# ============================================================

trackable_snv <- rbindlist(
  denominator_list,
  use.names = TRUE,
  fill = TRUE
)


setorder(
  trackable_snv,
  colony_id,
  MAG_id,
  Month_t1,
  Month_t2
)


# ============================================================
# Basic QC
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "TRACKABLE-SNV QC\n"
)

cat(
  "============================================================\n"
)

cat(
  "Transitions:",
  nrow(trackable_snv),
  "\n\n"
)


cat(
  "Trackable SNVs per transition:\n"
)

print(
  summary(
    trackable_snv$n_trackable_SNVS
  )
)


cat(
  "\nQuantiles:\n"
)

print(
  quantile(
    trackable_snv$n_trackable_SNVS,
    probs = c(
      0,
      0.01,
      0.05,
      0.10,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)


cat(
  "\nEndpoint overlap fraction:\n"
)

print(
  summary(
    trackable_snv$overlap_fraction_min
  )
)


cat(
  "\nTransitions with zero trackable SNVs:",
  sum(
    trackable_snv$n_trackable_SNVS == 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "============================================================\n"
)





check_dt <- merge(
  trackable_snv,
  sweep_breadth_num[
    ,
    .(
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      n_sweep_SNVS
    )
  ],
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)

check_dt[
  is.na(n_sweep_SNVS),
  n_sweep_SNVS := 0L
]


summary(check_dt$n_trackable_SNVS)

summary(check_dt$n_sweep_SNVS)

check_dt[
  n_sweep_SNVS > n_trackable_SNVS
]



names(sweep_dt)

sweep_dt[
  ,
  sweep_SNV_ID := paste(
    MAG_id,
    contig_name,
    pos_in_contig,
    sep = "::"
  )
]

sweep_check <- sweep_dt[
  ,
  .(
    n_rows = .N,
    n_unique_SNVS = uniqueN(sweep_SNV_ID)
  ),
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2
  )
]

summary(
  sweep_check$n_rows /
    sweep_check$n_unique_SNVS
)

sweep_check[
  n_rows != n_unique_SNVS
][1:20]



bad_pairs <- check_dt[
  n_sweep_SNVS > n_trackable_SNVS
]

merge(
  bad_pairs[
    ,
    .(
      MAG_id,
      colony_id,
      Month_t1,
      Month_t2,
      n_trackable_SNVS,
      n_sweep_SNVS
    )
  ],
  sweep_check,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)




#rebuilding numerator
# Physical SNV identity
sweep_dt[
  ,
  sweep_SNV_ID := paste(
    MAG_id,
    contig_name,
    pos_in_contig,
    sep = "::"
  )
]

sweep_breadth_num <- sweep_dt[
  !is.na(Month_t1) &
    !is.na(Month_t2),
  .(
    n_sweep_SNVS =
      uniqueN(sweep_SNV_ID),
    
    n_0D_sweep_SNVS =
      uniqueN(
        sweep_SNV_ID[
          degeneracy == "0D"
        ]
      ),
    
    n_4D_sweep_SNVS =
      uniqueN(
        sweep_SNV_ID[
          degeneracy == "4D"
        ]
      )
  ),
  by = .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2
  )
]



check_dt <- merge(
  trackable_snv,
  sweep_breadth_num,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)

for (v in c(
  "n_sweep_SNVS",
  "n_0D_sweep_SNVS",
  "n_4D_sweep_SNVS"
)) {
  check_dt[
    is.na(get(v)),
    (v) := 0L
  ]
}

check_dt[
  n_sweep_SNVS > n_trackable_SNVS
]




summary(check_dt$sweep_breadth)

quantile(
  check_dt$sweep_breadth,
  probs = c(
    0,
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  ),
  na.rm = TRUE
)

# Remove old object if it exists
if (exists("sweep_breadth")) {
  rm(sweep_breadth)
}

# Create the actual column explicitly
check_dt[
  ,
  sweep_breadth_frac :=
    n_sweep_SNVS / n_trackable_SNVS
]

summary(check_dt$sweep_breadth_frac)

quantile(
  check_dt$sweep_breadth_frac,
  probs = c(
    0,
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  ),
  na.rm = TRUE
)
focus_mags <- c(
  "mmag_254",
  "mmag_191",
  "mmag_411",
  "mmag_648",
  "mmag_283",
  "mmag_139"
)

check_dt[
  MAG_id %chin% focus_mags,
  .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    n_trackable_SNVS,
    n_sweep_SNVS,
    sweep_breadth_frac
  )
][
  order(
    MAG_id,
    colony_id,
    Month_t1
  )
]





#making private-marker turnover!

# Standardize IDs just to be safe
pm_trans_good[, MAG_id := as.character(MAG_id)]
pm_trans_good[, Colony := as.character(Colony)]

check_dt[, MAG_id := as.character(MAG_id)]
check_dt[, colony_id := as.character(colony_id)]

# Only take the sweep-related columns we need from check_dt
sweep_metric <- check_dt[
  ,
  .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    n_trackable_SNVS,
    n_sweep_SNVS,
    sweep_breadth_frac
  )
]

# Add colony_id to private-marker table
pm_trans_good[
  ,
  colony_id := as.character(Colony)
]

# Merge by exact population transition
evo_dt <- merge(
  pm_trans_good,
  sweep_metric,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE
)



nrow(pm_trans_good)
nrow(evo_dt)
sum(is.na(evo_dt$n_trackable_SNVS))
sum(is.na(evo_dt$sweep_breadth_frac))


evo_sweep <- evo_dt[
  n_sweep_SNVS > 0 &
    is.finite(sweep_breadth_frac) &
    is.finite(delta_background_turnover) &
    is.finite(background_turnover_t2)
]

cat(
  "ALL QC-PASSING TRANSITIONS:",
  nrow(evo_dt),
  "\n"
)

cat(
  "SWEEP-POSITIVE TRANSITIONS:",
  nrow(evo_sweep),
  "\n"
)

cat(
  "MAG X COLONY POPULATIONS REPRESENTED:",
  uniqueN(
    paste(
      evo_sweep$MAG_id,
      evo_sweep$colony_id
      
    )
  ),
  "\n"
)


#figure a, landscape of events
# x = deltaB
# y = S = N_sweep/N_trackable

p_event <- ggplot(
  evo_sweep,
  aes(
    x = delta_background_turnover,
    y = sweep_breadth_frac
  )
)+
  geom_point(
    aes(size = n_sweep_SNVS),
    alpha = 0.5
  )+
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  )+
  scale_y_log10(
    labels = label_percent(
      accuracy = 0.001
    )
  )+
  scale_size_continuous(
    range = c(1.5, 7)
  ) +
  theme_classic()

p_event


#figure B - endpoint, how does it look at end?

p_state <- ggplot(
  evo_sweep,
  aes(
    x = background_turnover_t2,
    y = sweep_breadth_frac
  )
) +
  geom_point(
    aes(size = n_sweep_SNVS),
    alpha = 0.45
  ) +
  scale_y_log10(
    labels = label_percent(
      accuracy = 0.001
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1)
  ) +
  scale_size_continuous(
    range = c(1.5, 7)
  ) +
  labs(
    x = "Resident-background turnover at t2",
    y = "Fraction of trackable SNVs with significant AF change",
    size = "Sweep SNVs"
  ) +
  theme_classic(
    base_size = 13
  )

print(p_state)



focus_mags <- c(
  "mmag_254",
  "mmag_191",
  "mmag_411",
  "mmag_648",
  "mmag_283",
  "mmag_139"
)

focus_labels <- evo_sweep[
  MAG_id %chin% focus_mags,
  .SD[
    which.max(sweep_breadth_frac)
  ],
  by = MAG_id
]

focus_labels[
  ,
  plot_label := paste0(
    MAG_id,
    "\nC",
    colony_id,
    " ",
    Month_t1,
    "\u2192",
    Month_t2
  )
]


p_event_labeled <- p_event +
  geom_point(
    data = focus_labels,
    aes(
      x = delta_background_turnover,
      y = sweep_breadth_frac
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_text_repel(
    data = focus_labels,
    aes(
      x = delta_background_turnover,
      y = sweep_breadth_frac,
      label = plot_label
    ),
    inherit.aes = FALSE,
    size = 3.3,
    max.overlaps = Inf
  )

print(p_event_labeled)



p_state_labeled <- p_state +
  geom_point(
    data = focus_labels,
    aes(
      x = background_turnover_t2,
      y = sweep_breadth_frac
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_text_repel(
    data = focus_labels,
    aes(
      x = background_turnover_t2,
      y = sweep_breadth_frac,
      label = plot_label
    ),
    inherit.aes = FALSE,
    size = 3.3,
    max.overlaps = Inf
  )

print(p_state_labeled)





ggplot(
  evo_sweep,
  aes(
    x = delta_background_turnover,
    y = sweep_breadth_frac
  )
) +
  geom_point(
    aes(
      color = background_turnover_t2
    ),
    alpha = 0.55,
    size = 2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_y_log10(
    labels = scales::label_percent(
      accuracy = 0.001
    )
  ) +
  scale_color_viridis_c(
    limits = c(0, 1),
    name = "Resident-background\nturnover at t2"
  ) +
  labs(
    x = expression(
      Delta*" resident-background turnover"
    ),
    y = "Fraction of trackable SNVs with significant AF change"
  ) +
  theme_classic()












# ------------------------------------------------------------
# Prepare analysis table
# ------------------------------------------------------------

evo_sweep[
  ,
  population_id := paste(
    MAG_id,
    colony_id,
    sep = "::"
  )
]

evo_sweep[
  ,
  abs_delta_B :=
    abs(delta_background_turnover)
]


# ------------------------------------------------------------
# Cluster bootstrap
#
# Resample entire MAG x colony trajectories with replacement.
# ------------------------------------------------------------

cluster_bootstrap <- function(
    dt,
    stat_fun,
    B = 1000,
    seed = 123
) {
  
  clusters <- unique(dt$population_id)
  
  set.seed(seed)
  
  boot_values <- vapply(
    seq_len(B),
    function(b) {
      
      sampled_clusters <- sample(
        clusters,
        size = length(clusters),
        replace = TRUE
      )
      
      boot_dt <- rbindlist(
        lapply(
          sampled_clusters,
          function(z) {
            dt[population_id == z]
          }
        ),
        use.names = TRUE
      )
      
      stat_fun(boot_dt)
    },
    numeric(1)
  )
  
  boot_values <- boot_values[
    is.finite(boot_values)
  ]
  
  estimate <- stat_fun(dt)
  
  ci <- quantile(
    boot_values,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  list(
    estimate = estimate,
    lower_95 = unname(ci[1]),
    upper_95 = unname(ci[2]),
    boot_values = boot_values
  )
}




stat_abs_delta <- function(d) {
  
  cor(
    d$abs_delta_B,
    d$sweep_breadth_frac,
    method = "spearman",
    use = "complete.obs"
  )
}


stat_state <- function(d) {
  
  cor(
    d$background_turnover_t2,
    d$sweep_breadth_frac,
    method = "spearman",
    use = "complete.obs"
  )
}



rho_abs_delta <- cluster_bootstrap(
  evo_sweep,
  stat_abs_delta,
  B = 1000,
  seed = 123
)

rho_state <- cluster_bootstrap(
  evo_sweep,
  stat_state,
  B = 1000,
  seed = 456
)

rho_abs_delta[
  c(
    "estimate",
    "lower_95",
    "upper_95"
  )
]

rho_state[
  c(
    "estimate",
    "lower_95",
    "upper_95"
  )
]




m_delta <- lmer(
  log10(sweep_breadth_frac) ~
    abs_delta_B +
    (1 | population_id),
  data = evo_sweep
)

summary(m_delta)
confint(
  m_delta,
  parm = "abs_delta_B",
  method = "Wald"
)



PRIMARY_DELTA <- 0.10
PRIMARY_STATE <- 0.20

evo_sweep[
  ,
  resident_preserving :=
    abs_delta_B <= PRIMARY_DELTA &
    background_turnover_t2 <= PRIMARY_STATE
]




evo_sweep[
  ,
  .(
    n_sweep_positive_transitions = .N,
    
    n_resident_preserving =
      sum(resident_preserving),
    
    prop_resident_preserving =
      mean(resident_preserving),
    
    median_sweep_breadth_retained =
      median(
        sweep_breadth_frac[
          resident_preserving
        ],
        na.rm = TRUE
      ),
    
    q90_sweep_breadth_retained =
      quantile(
        sweep_breadth_frac[
          resident_preserving
        ],
        0.90,
        na.rm = TRUE
      )
  )
]



stat_prevalence <- function(d) {
  
  mean(
    d$abs_delta_B <= 0.10 &
      d$background_turnover_t2 <= 0.20,
    na.rm = TRUE
  )
}

prev_primary <- cluster_bootstrap(
  evo_sweep,
  stat_prevalence,
  B = 1000,
  seed = 789
)

prev_primary[
  c(
    "estimate",
    "lower_95",
    "upper_95"
  )
]


#sensitivity
cut_grid <- CJ(
  delta_cut = c(
    0.05,
    0.10,
    0.20
  ),
  state_cut = c(
    0.10,
    0.20
  )
)

prevalence_results <- rbindlist(
  lapply(
    seq_len(nrow(cut_grid)),
    function(i) {
      
      dc <- cut_grid$delta_cut[i]
      sc <- cut_grid$state_cut[i]
      
      stat_i <- function(d) {
        
        mean(
          d$abs_delta_B <= dc &
            d$background_turnover_t2 <= sc,
          na.rm = TRUE
        )
      }
      
      result <- cluster_bootstrap(
        evo_sweep,
        stat_i,
        B = 2000,
        seed = 1000 + i
      )
      
      data.table(
        delta_cut = dc,
        state_cut = sc,
        
        prevalence =
          result$estimate,
        
        lower_95 =
          result$lower_95,
        
        upper_95 =
          result$upper_95
      )
    }
  )
)

prevalence_results[
  ,
  `:=`(
    prevalence_pct =
      100 * prevalence,
    
    lower_95_pct =
      100 * lower_95,
    
    upper_95_pct =
      100 * upper_95
  )
]

prevalence_results






population_summary <- evo_sweep[
  ,
  .(
    n_sweep_transitions = .N,
    
    any_resident_preserving =
      any(resident_preserving),
    
    n_resident_preserving =
      sum(resident_preserving),
    
    max_sweep_breadth =
      max(
        sweep_breadth_frac,
        na.rm = TRUE
      )
  ),
  by = .(
    MAG_id,
    colony_id,
    population_id
  )
]

population_summary[
  ,
  .(
    n_sweep_populations = .N,
    
    n_with_resident_preserving_sweep =
      sum(any_resident_preserving),
    
    prop_with_resident_preserving_sweep =
      mean(any_resident_preserving)
  )
]




mag_summary <- population_summary[
  ,
  .(
    any_resident_preserving =
      any(any_resident_preserving)
  ),
  by = MAG_id
]

mag_summary[
  ,
  .(
    n_MAGs_with_sweeps = .N,
    
    n_MAGs_with_resident_preserving_sweep =
      sum(any_resident_preserving),
    
    proportion =
      mean(any_resident_preserving)
  )
]





#panel C: do substantial genetic changes require substantial ecological abundance changes?
  #for each MAG x colony transition found in evo_sweep, calculate the replicate-averaged relative abundance at t1 and t2 then:
    #A = |log_2*(mean(p_t2)/mean(p_t1))|



ecology <- read.csv('/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv')

abund <- as.data.table(ecology)

abund[, Colony := as.character(Colony)]
abund[, MAG_id := as.character(MAG_id)]
abund[, Month := tools::toTitleCase(as.character(Month))]

# Relative abundance within each biological replicate/sample
abund[
  ,
  rel_abund :=
    mean / sum(mean, na.rm = TRUE),
  by = sample_name
]


abund_month <- abund[
  ,
  .(
    mean_rel_abund =
      mean(rel_abund, na.rm = TRUE),
    
    sd_rel_abund =
      sd(rel_abund, na.rm = TRUE),
    
    n_rep =
      uniqueN(Rep)
  ),
  by = .(
    MAG_id,
    Colony,
    Month
  )
]



abund_t1 <- copy(abund_month)

setnames(
  abund_t1,
  c(
    "Colony",
    "Month",
    "mean_rel_abund",
    "sd_rel_abund",
    "n_rep"
  ),
  c(
    "colony_id",
    "Month_t1",
    "abundance_t1",
    "abundance_sd_t1",
    "n_rep_t1"
  )
)


abund_t2 <- copy(abund_month)

setnames(
  abund_t2,
  c(
    "Colony",
    "Month",
    "mean_rel_abund",
    "sd_rel_abund",
    "n_rep"
  ),
  c(
    "colony_id",
    "Month_t2",
    "abundance_t2",
    "abundance_sd_t2",
    "n_rep_t2"
  )
)



evo_ecol <- merge(
  evo_sweep,
  abund_t1,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1"
  ),
  all.x = TRUE
)

evo_ecol <- merge(
  evo_ecol,
  abund_t2,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t2"
  ),
  all.x = TRUE
)



c(
  n_total = nrow(evo_ecol),
  missing_t1 = sum(is.na(evo_ecol$abundance_t1)),
  missing_t2 = sum(is.na(evo_ecol$abundance_t2)),
  zero_t1 = sum(evo_ecol$abundance_t1 == 0, na.rm = TRUE),
  zero_t2 = sum(evo_ecol$abundance_t2 == 0, na.rm = TRUE)
)



evo_ecol[
  abundance_t1 > 0 &
    abundance_t2 > 0,
  log2_abundance_change :=
    log2(abundance_t2 / abundance_t1)
]

evo_ecol[
  ,
  abs_log2_abundance_change :=
    abs(log2_abundance_change)
]




p_C <- ggplot(
  evo_ecol[
    is.finite(abs_log2_abundance_change)
  ],
  aes(
    x = abs_log2_abundance_change,
    y = sweep_breadth_frac
  )
) +
  geom_point(
    aes(
      color = background_turnover_t2
    ),
    alpha = 0.55,
    size = 2
  ) +
  scale_y_log10(
    labels = label_percent(
      accuracy = 0.001
    )
  ) +
  scale_color_viridis_c(
    limits = c(0, 1),
    name = "Resident-background\nturnover at t2",
    option = "cividis"
  ) +
  labs(
    x = expression(
      "|"*log[2]*" abundance fold-change|"
    ),
    y = "Fraction of trackable SNVs\nwith significant AF change"
  ) +
  theme_classic(
    base_size = 13
  )

print(p_C)

focus_C <- evo_ecol[
  MAG_id %chin% c(
    "mmag_191",
    "mmag_254",
    "mmag_668"
  ) &
    is.finite(abs_log2_abundance_change)
]

p_C_focus <- p_C +
  geom_point(
    data = focus_C,
    aes(
      x = abs_log2_abundance_change,
      y = sweep_breadth_frac
    ),
    inherit.aes = FALSE,
    size = 3
  )

print(p_C_focus)



evo_ecol[
  ,
  population_id := paste(
    MAG_id,
    colony_id,
    sep = "::"
  )
]

evo_ecol_good <- evo_ecol[
  is.finite(abs_log2_abundance_change) &
    is.finite(sweep_breadth_frac)
]

stat_abundance_sweep <- function(d) {
  
  cor(
    d$abs_log2_abundance_change,
    d$sweep_breadth_frac,
    method = "spearman",
    use = "complete.obs"
  )
}

rho_abundance_sweep <- cluster_bootstrap(
  evo_ecol[
    is.finite(abs_log2_abundance_change)
  ],
  stat_abundance_sweep,
  B = 1000,
  seed = 321
)



ABUND_STABLE <- 1.0
DELTA_STABLE <- 0.10
STATE_RETAINED <- 0.20

evo_ecol_good[
  ,
  abundance_stable :=
    abs_log2_abundance_change <= ABUND_STABLE
]

evo_ecol_good[
  ,
  background_retained :=
    abs_delta_B <= DELTA_STABLE &
    background_turnover_t2 <= STATE_RETAINED
]

evo_ecol_good[
  ,
  dynamic_class := fcase(
    
    abundance_stable & background_retained,
    "Stable abundance + retained background",
    
    abundance_stable & !background_retained,
    "Stable abundance + background turnover",
    
    !abundance_stable & background_retained,
    "Abundance change + retained background",
    
    !abundance_stable & !background_retained,
    "Abundance change + background turnover"
  )
]


class_summary <- evo_ecol_good[
  ,
  .(
    n = .N,
    proportion = .N / nrow(evo_ecol_good),
    
    median_sweep_breadth =
      median(sweep_breadth_frac, na.rm = TRUE),
    
    median_abs_delta_B =
      median(abs_delta_B, na.rm = TRUE),
    
    median_abundance_change =
      median(abs_log2_abundance_change, na.rm = TRUE)
  ),
  by = dynamic_class
][
  order(-proportion)
]

class_summary[
  ,
  proportion_pct := 100 * proportion
]

class_summary





evo_ecol_good[
  ,
  cryptic_positive_turnover :=
    abs_log2_abundance_change <= 1 &
    delta_background_turnover > 0.10 &
    background_turnover_t2 > 0.20
]

evo_ecol_good[
  ,
  cryptic_recovery :=
    abs_log2_abundance_change <= 1 &
    delta_background_turnover < -0.10
]

evo_ecol_good[
  ,
  .(
    n_total = .N,
    
    n_cryptic_positive_turnover =
      sum(cryptic_positive_turnover),
    
    prop_cryptic_positive_turnover =
      mean(cryptic_positive_turnover),
    
    n_cryptic_recovery =
      sum(cryptic_recovery),
    
    prop_cryptic_recovery =
      mean(cryptic_recovery)
  )
]



# Make sure population_id exists
evo_ecol_good[
  ,
  population_id := paste(
    MAG_id,
    colony_id,
    sep = "::"
  )
]

# Statistic:
# association between magnitude of abundance change
# and sweep breadth
stat_abundance_sweep <- function(d) {
  
  cor(
    d$abs_log2_abundance_change,
    d$sweep_breadth_frac,
    method = "spearman",
    use = "complete.obs"
  )
}

rho_abundance_sweep <- cluster_bootstrap(
  evo_ecol_good,
  stat_abundance_sweep,
  B = 1000,
  seed = 321
)

rho_abundance_sweep[
  c(
    "estimate",
    "lower_95",
    "upper_95"
  )
]



stat_abundance_background <- function(d) {
  
  cor(
    d$abs_log2_abundance_change,
    d$abs_delta_B,
    method = "spearman",
    use = "complete.obs"
  )
}

rho_abundance_background <- cluster_bootstrap(
  evo_ecol_good,
  stat_abundance_background,
  B = 1000,
  seed = 322
)

rho_abundance_background[
  c(
    "estimate",
    "lower_95",
    "upper_95"
  )
]




abundance_cuts <- c(0.5, 1, 2)

abundance_sensitivity <- rbindlist(
  lapply(
    seq_along(abundance_cuts),
    function(i) {
      
      ac <- abundance_cuts[i]
      
      stat_stable <- function(d) {
        mean(
          d$abs_log2_abundance_change <= ac,
          na.rm = TRUE
        )
      }
      
      stat_stable_retained <- function(d) {
        mean(
          d$abs_log2_abundance_change <= ac &
            d$abs_delta_B <= 0.10 &
            d$background_turnover_t2 <= 0.20,
          na.rm = TRUE
        )
      }
      
      stat_stable_not_retained <- function(d) {
        mean(
          d$abs_log2_abundance_change <= ac &
            !(
              d$abs_delta_B <= 0.10 &
                d$background_turnover_t2 <= 0.20
            ),
          na.rm = TRUE
        )
      }
      
      b1 <- cluster_bootstrap(
        evo_ecol_good,
        stat_stable,
        B = 2000,
        seed = 2000 + i
      )
      
      b2 <- cluster_bootstrap(
        evo_ecol_good,
        stat_stable_retained,
        B = 2000,
        seed = 3000 + i
      )
      
      b3 <- cluster_bootstrap(
        evo_ecol_good,
        stat_stable_not_retained,
        B = 2000,
        seed = 4000 + i
      )
      
      data.table(
        abundance_cut = ac,
        approximate_fold_change = 2^ac,
        
        prop_abundance_stable = b1$estimate,
        stable_lower = b1$lower_95,
        stable_upper = b1$upper_95,
        
        prop_stable_retained = b2$estimate,
        retained_lower = b2$lower_95,
        retained_upper = b2$upper_95,
        
        prop_stable_not_retained = b3$estimate,
        not_retained_lower = b3$lower_95,
        not_retained_upper = b3$upper_95
      )
    }
  )
)

abundance_sensitivity[
  ,
  c(
    "prop_abundance_stable",
    "stable_lower",
    "stable_upper",
    "prop_stable_retained",
    "retained_lower",
    "retained_upper",
    "prop_stable_not_retained",
    "not_retained_lower",
    "not_retained_upper"
  ) :=
    lapply(
      .SD,
      function(x) 100 * x
    ),
  .SDcols = c(
    "prop_abundance_stable",
    "stable_lower",
    "stable_upper",
    "prop_stable_retained",
    "retained_lower",
    "retained_upper",
    "prop_stable_not_retained",
    "not_retained_lower",
    "not_retained_upper"
  )
]

abundance_sensitivity





deg_metric <- check_dt[
  ,
  .(
    MAG_id,
    colony_id,
    Month_t1,
    Month_t2,
    n_sweep_SNVS,
    n_0D_sweep_SNVS,
    n_4D_sweep_SNVS
  )
]

evo_deg <- merge(
  evo_sweep,
  deg_metric,
  by = c(
    "MAG_id",
    "colony_id",
    "Month_t1",
    "Month_t2"
  ),
  all.x = TRUE,
  suffixes = c("", "_deg")
)




evo_deg[
  ,
  n_0D4D :=
    n_0D_sweep_SNVS +
    n_4D_sweep_SNVS
]

evo_deg[
  n_0D4D > 0,
  prop_4D_among_0D4D :=
    n_4D_sweep_SNVS / n_0D4D
]




cor.test(
  evo_deg$abs_delta_B,
  evo_deg$prop_4D_among_0D4D,
  method = "spearman",
  exact = FALSE
)




stat_4D_turnover <- function(d) {
  
  x <- d[
    is.finite(prop_4D_among_0D4D)
  ]
  
  cor(
    x$abs_delta_B,
    x$prop_4D_among_0D4D,
    method = "spearman",
    use = "complete.obs"
  )
}

rho_4D_turnover <- cluster_bootstrap(
  evo_deg,
  stat_4D_turnover,
  B = 1000,
  seed = 5001
)

rho_4D_turnover[
  c("estimate", "lower_95", "upper_95")
]









focus_mags <- c(
  "mmag_254",
  "mmag_385",
  "mmag_688"
)

pB <- ggplot(
  evo_sweep,
  aes(
    x = delta_background_turnover,
    y = sweep_breadth_frac,
    color = background_turnover_t2
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.5
  ) +
  geom_point(
    data = evo_sweep[
      MAG_id %chin% focus_mags
    ],
    size = 3.2,
    shape = 21,
    stroke = 0.8
  ) +
  geom_text_repel(
    data = evo_sweep[
      MAG_id %chin% focus_mags
    ],
    aes(label = MAG_id),
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_y_log10(
    labels = label_percent(
      accuracy = 0.001
    )
  ) +
  scale_color_viridis_c(
    limits = c(0, 1),
    name = "Resident-background\nturnover at t2",
    option="cividis"
  ) +
  labs(
    x = expression(
      Delta*" resident-background turnover"
    ),
    y = "Trackable SNVs with\nsignificant AF change"
  ) +
  theme_classic(base_size = 12)

pB
