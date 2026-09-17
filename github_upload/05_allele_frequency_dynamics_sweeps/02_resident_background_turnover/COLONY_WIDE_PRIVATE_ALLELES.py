#!/usr/bin/env python3
"""
colony_private_marker_all_colonies_with_monthly_change.py
---------------------------------------------------------

Extends the private marker analysis to compute month-by-month
changes in preserved and disrupted alleles per MAG and colony.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import glob

# ---------------- CONFIG ----------------
TRAJ_DIR = Path("/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/snv_frequency_trajectories")
OUT_DIR = TRAJ_DIR.parent / "private_marker_results_all_colonies"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MONTH_ORDER = ["May","June","July","August","September","October","November",'January',"February"]
HIGH_FREQ = 0.8
PROP_HIGH_REQ = 0.8
MIN_FREQ_REQ = 0.5
MIN_SNVS_PER_COLONY = 500
# ----------------------------------------

def classify_snv(sub):
    freqs = sub["major_freq"].dropna()
    n_total = len(freqs)
    n_high = (freqs >= HIGH_FREQ).sum()
    prop_high = n_high / n_total if n_total else np.nan
    min_freq = freqs.min() if n_total else np.nan

    if n_total == 0:
        status = "Missing"
    elif (prop_high >= PROP_HIGH_REQ) and (min_freq >= MIN_FREQ_REQ):
        status = "Preserved_private"
    else:
        status = "Disrupted_private"

    return pd.Series({
        "n_months": n_total,
        "prop_high": prop_high,
        "min_freq": min_freq,
        "status": status
    })

# ---------- Main loop ----------
trajectory_files = sorted(TRAJ_DIR.glob("snv_frequency_trajectory_*.csv"))
all_summaries = []
monthly_trends = []  # <-- new container for per-month changes

for traj_file in trajectory_files:
    colony_id = traj_file.stem.split("_")[-1]
    print(f"\n▶ Processing colony {colony_id}")

    df = pd.read_csv(traj_file)
    df["Month"] = pd.Categorical(df["Month"].str.title(), categories=MONTH_ORDER, ordered=True)
    df["major_freq"] = 1 - df["departure_from_polarized_consensus"]
    df["SNV_ID"] = (
        df["MAG_id"].astype(str) + "_" +
        df["contig_name"].astype(str) + "_" +
        df["pos_in_contig"].astype(str)
    )

    # identify May private SNVs (departure 0 or 1)
    may_df = df[df["Month"] == "May"]
    private_ids = may_df.loc[
        may_df["departure_from_polarized_consensus"].isin([0.0, 1.0]),
        "SNV_ID"
    ].unique()

    n_private = len(private_ids)
    print(f"  Found {n_private:,} May-fixed private SNVs.")
    if n_private == 0:
        continue

    df_private = df[df["SNV_ID"].isin(private_ids)].copy()

    # ---- Classify preserved/disrupted (overall) ----
    classified = (
        df_private.groupby(["MAG_id","contig_name","pos_in_contig"], group_keys=False)
        .apply(classify_snv)
        .reset_index()
    )

    summary = (
        classified.groupby(["MAG_id","status"], as_index=False)
        .size()
        .pivot(index="MAG_id", columns="status", values="size")
        .fillna(0)
    )
    summary["Colony"] = colony_id
    num_cols = summary.select_dtypes(include=[np.number]).columns
    summary["Total_private"] = summary[num_cols].sum(axis=1)
    summary["Frac_preserved"] = summary.get("Preserved_private",0) / summary["Total_private"]
    summary["Frac_disrupted"] = summary.get("Disrupted_private",0) / summary["Total_private"]
    all_summaries.append(summary.reset_index())

    # ---- Month-by-month trend ----
    # For each MAG, count how many private SNVs are high vs disrupted each month
    for mag, sub in df_private.groupby("MAG_id"):
        baseline_ids = sub["SNV_ID"].unique()
        for month in MONTH_ORDER:
            sub_month = sub[sub["Month"] == month]
            if sub_month.empty:
                continue
            freqs = sub_month.groupby("SNV_ID")["major_freq"].mean()  # average freq per SNV that month
            n_total = len(baseline_ids)
            n_high = (freqs >= HIGH_FREQ).sum()
            n_disrupted = (freqs < MIN_FREQ_REQ).sum()
            monthly_trends.append({
                "Colony": colony_id,
                "MAG_id": mag,
                "Month": month,
                "n_private_baseline": n_total,
                "n_high_freq": n_high,
                "n_disrupted": n_disrupted,
                "prop_high_freq": n_high / n_total if n_total else np.nan,
                "prop_disrupted": n_disrupted / n_total if n_total else np.nan
            })

# ---------- Save colony summaries ----------
all_summary_df = pd.concat(all_summaries, ignore_index=True)
all_summary_out = OUT_DIR / "private_marker_summary_all_colonies.csv"
all_summary_df.to_csv(all_summary_out, index=False)
print(f"\n✅ Wrote all-colony summary: {all_summary_out}")

# ---------- Save monthly trends ----------
monthly_df = pd.DataFrame(monthly_trends)
monthly_out = OUT_DIR / "private_marker_monthly_trends_all_colonies.csv"
monthly_df.to_csv(monthly_out, index=False)
print(f"✅ Wrote monthly trend data: {monthly_out}")

# ---------- Cross-colony MAGs ----------
print("\nScanning for MAGs shared across colonies (≥500 SNVs per colony)...")
mag_colony_counts = (
    all_summary_df[["MAG_id","Colony","Total_private"]]
    .groupby(["MAG_id","Colony"], as_index=False)
    .sum()
)
mag_colony_counts = mag_colony_counts[mag_colony_counts["Total_private"] >= MIN_SNVS_PER_COLONY]
multi_colony = (
    mag_colony_counts.groupby("MAG_id")
    .filter(lambda x: len(x["Colony"].unique()) > 1)
)
if len(multi_colony) > 0:
    out_path = OUT_DIR / "MAGs_in_multiple_colonies_over500.csv"
    multi_colony.to_csv(out_path, index=False)
    print(f"⚠️  {len(multi_colony['MAG_id'].unique())} MAGs appear in >1 colony (≥500 SNVs) → {out_path}")
else:
    print("No MAGs appear in multiple colonies with ≥500 SNVs.")
