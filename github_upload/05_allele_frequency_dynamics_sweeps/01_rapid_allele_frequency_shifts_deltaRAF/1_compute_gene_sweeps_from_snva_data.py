#!/usr/bin/env python3
# compute_gene_sweeps_from_snvs.py
"""
Compute per-gene per-month:
    - mean departure_from_polarized_consensus
    - Δp (delta)
    - abs_delta

Define sweep-like events using:
    abs_delta >= max(0.5, q95_of_MAG_specific_abs_delta)

Output:
    - per-colony gene_month_means_<colony>.csv
    - combined gene_month_deltas_with_sweepflags.parquet
"""

import os, glob
import pandas as pd
import numpy as np
from tqdm import tqdm

SNV_DIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/snv_frequency_trajectories"
OUTDIR = "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs"
os.makedirs(OUTDIR, exist_ok=True)

MONTH_ORDER = ["May","June","July","August","September","October","November","January","February"]
ABS_MIN_DELTA = 0.5        # absolute minimum Δp
REL_QUANTILE = 0.95        # relative Δp threshold

# ------------------------------------------------------------
# STEP 1 — Compute per-gene per-month mean_p and deltas
# ------------------------------------------------------------
files = sorted(glob.glob(os.path.join(SNV_DIR,"snv_frequency_trajectory_*.csv")))
all_means = []

for fp in tqdm(files, desc="Per-colony SNV -> per-gene summarize"):
    colony = os.path.basename(fp).split('_')[-1].replace('.csv','')

    usecols = ['MAG_id','Month','unique_gene_callers_id','departure_from_polarized_consensus']
    df = pd.read_csv(fp, usecols=lambda c: c in usecols, low_memory=False)
    df = df.dropna(subset=['unique_gene_callers_id','departure_from_polarized_consensus'])

    gm = df.groupby(
        ['MAG_id','unique_gene_callers_id','Month'], as_index=False
    )['departure_from_polarized_consensus'].agg(['mean','count']).reset_index()

    gm = gm.rename(columns={'mean':'mean_departure','count':'n_snvs'})
    gm['colony_id'] = colony
    gm.to_csv(os.path.join(OUTDIR, f"gene_month_means_{colony}.csv"), index=False)
    all_means.append(gm)

combined = pd.concat(all_means, ignore_index=True)
combined['Month'] = pd.Categorical(combined['Month'], categories=MONTH_ORDER, ordered=True)
combined = combined.sort_values(['colony_id','MAG_id','unique_gene_callers_id','Month'])

# Compute Δp for each gene time series
rows = []
for (colony, mag, gid), sub in combined.groupby(['colony_id','MAG_id','unique_gene_callers_id']):
    sub = sub.sort_values('Month').reset_index(drop=True)
    sub['prev_departure'] = sub['mean_departure'].shift(1)
    sub['delta'] = sub['mean_departure'] - sub['prev_departure']
    sub['abs_delta'] = sub['delta'].abs()
    rows.append(sub)

combined_deltas = pd.concat(rows, ignore_index=True)

# ------------------------------------------------------------
# STEP 2 — Add MAG-specific relative Δp threshold (q95)
# ------------------------------------------------------------
mag_q95 = combined_deltas.groupby('MAG_id')['abs_delta'].quantile(REL_QUANTILE)
mag_q95 = mag_q95.rename('mag_q95')

combined_deltas = combined_deltas.merge(mag_q95, on='MAG_id', how='left')

# Sweep definition:
# abs_delta >= max(ABS_MIN_DELTA, MAG_specific_q95)
combined_deltas['delta_thresh'] = np.maximum(ABS_MIN_DELTA, combined_deltas['mag_q95'])
combined_deltas['sweep_delta_flag'] = combined_deltas['abs_delta'] >= combined_deltas['delta_thresh']

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------
outp = os.path.join(OUTDIR, "gene_month_deltas_with_sweepflags.parquet")
combined_deltas.to_parquet(outp, index=False)
print("Wrote:", outp)

