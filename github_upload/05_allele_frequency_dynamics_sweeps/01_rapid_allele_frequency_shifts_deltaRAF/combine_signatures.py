#!/usr/bin/env python3
"""
integration_master_script.py
----------------------------

Integrates:
    Δp sweeps
    π statistics
    SNV-density states
    Poisson SNV-density sweeps
    Seasonal classification
    Evolutionary pattern classification
    Sweep rates (within vs between)
    Poisson rate test

Output:
    combined_gene_selection_signatures.parquet
    evolutionary_patterns.csv
"""

import os
import json
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import norm
from tqdm import tqdm

# --------------------------------------------------------
# CONFIG
# --------------------------------------------------------
OUTDIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs"

PI = f"{OUTDIR}/pi_statistics.parquet"
DP = f"{OUTDIR}/gene_month_deltas_with_sweepflags.parquet"
SNVD = f"{OUTDIR}/gene_snv_density_states.parquet"
POISS = f"{OUTDIR}/poisson_timeseries_sweep_results.csv"
ALLELIC = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/significant_sweep_snvs_ALL.csv"

MONTH_ORDER = ["May","June","July","August","September","October","November","January","February"]
SEASON_MAP = {
    "May":"Spring","June":"Spring",
    "July":"Summer","August":"Summer",
    "September":"Summer","October":"Autumn","November":"Autumn",
    "January":"Winter","February":"Winter"
}

# --------------------------------------------------------
# LOAD ALL FRAMEWORKS
# --------------------------------------------------------
pi = pd.read_parquet(PI)
delta = pd.read_parquet(DP)
snvd = pd.read_parquet(SNVD)
poiss = pd.read_csv(POISS)
allelic = pd.read_csv(ALLELIC)

# --------------------------------------------------------
# CLEAN TYPES
# --------------------------------------------------------
for df in (pi, delta, snvd, poiss):
    df['colony_id'] = df['colony_id'].astype(str)
    df['MAG_id'] = df['MAG_id'].astype(str)
    df['unique_gene_callers_id'] = df['unique_gene_callers_id'].astype(str)
    df['Month'] = df['Month'].astype(str)

delta['Month'] = pd.Categorical(delta['Month'], categories=MONTH_ORDER, ordered=True)
snvd['Month'] = pd.Categorical(snvd['Month'], categories=MONTH_ORDER, ordered=True)
poiss['Month'] = pd.Categorical(poiss['Month'], categories=MONTH_ORDER, ordered=True)

# --------------------------------------------------------
# MERGE EVERYTHING
# --------------------------------------------------------
combined = delta.merge(
    snvd[['colony_id','MAG_id','Month','unique_gene_callers_id',
          'snv_density','delta_snv_density',
          'z_delta','snv_state','snv_state_confidence']],
    on=['colony_id','MAG_id','Month','unique_gene_callers_id'], how='left'
)

combined = combined.merge(
    pi[['colony_id','MAG_id','Month','unique_gene_callers_id','pi_1D','pi_4D','pN_pS_ratio']],
    on=['colony_id','MAG_id','Month','unique_gene_callers_id'], how='left'
)

combined = combined.merge(
    poiss[['MAG_id','Month','unique_gene_callers_id','observed_snps','expected_snps']],
    on=['MAG_id','Month','unique_gene_callers_id'], how='left'
)

# --------------------------------------------------------
# ADD SEASONS AND TRANSITION TYPES
# --------------------------------------------------------
combined = combined.sort_values(['colony_id','MAG_id','unique_gene_callers_id','Month'])
combined['prev_month'] = combined.groupby(
    ['colony_id','MAG_id','unique_gene_callers_id']
)['Month'].shift(1)

combined['season'] = combined['Month'].map(SEASON_MAP)
combined['prev_season'] = combined['prev_month'].map(SEASON_MAP)

combined['transition_type'] = np.where(
    combined['prev_season'].isna(),
    None,
    np.where(combined['season']==combined['prev_season'], 'within','between')
)

# --------------------------------------------------------
# EVOLUTIONARY PATTERN CLASSIFIER
# --------------------------------------------------------
def classify_pattern(row, af_thresh=0.1):
    sweep = row['sweep_delta_flag']
    snv  = row['snv_state']
    ad   = row['abs_delta']

    if sweep and snv == 2:
        return "Sweep_with_Diversity_Collapse"
    if sweep and snv == 0:
        return "Sweep_Only"
    if sweep and snv == 1:
        return "Sweep_with_Diversity_Increase"
    
    if not sweep and snv == 1 and ad < af_thresh:
        return "Diversity_Increase_No_AF_Change"
    if not sweep and snv == 2 and ad < af_thresh:
        return "Diversity_Decrease_No_AF_Change"
    
    if not sweep and snv == 0 and ad >= af_thresh:
        return "AF_Change_No_Diversity_Change"

    if not sweep and snv == 1 and ad >= af_thresh:
        return "Complex_Increase"
    if not sweep and snv == 2 and ad >= af_thresh:
        return "Complex_Decrease"

    return "Neutral"

combined['evolutionary_pattern'] = combined.apply(classify_pattern, axis=1)

# --------------------------------------------------------
# SAVE
# --------------------------------------------------------
combined.to_parquet(f"{OUTDIR}/combined_gene_selection_signatures.parquet", index=False)
combined.to_csv(f"{OUTDIR}/evolutionary_patterns.csv", index=False)

print("Done.")
