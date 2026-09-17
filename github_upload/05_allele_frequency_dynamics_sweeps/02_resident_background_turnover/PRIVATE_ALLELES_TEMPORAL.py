#!/usr/bin/env python3
"""
colony_private_marker_all_colonies_with_monthly_change_and_SNV_export.py
-------------------------------------------------------------------------

Runs private marker analysis across ALL colonies:

1. Identify private SNVs (fixed in May: departure ∈ {0,1})
2. Temporal classification:
       Preserved_private   (>=0.8 in >=80% months AND min>=0.5)
       Disrupted_private   (everything else)
3. Monthly trends for private SNVs
4. Cross-colony MAGs with >=500 private SNVs
5. NEW: Export ALL SNV-level trajectories in one file for plotting:
       private_marker_SNV_trajectories_all_colonies.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import glob

# ---------------- CONFIG ----------------
TRAJ_DIR = Path(
    "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/"
    "ALLSPECIES_intrapop_diversity/snv_frequency_trajectories"
)
OUT_DIR = TRAJ_DIR.parent / "private_marker_results_all_colonies"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MONTH_ORDER = ["May", "June", "July", "August", "September", "October","November","January","February"]
HIGH_FREQ = 0.8
PROP_HIGH_REQ = 0.8
MIN_FREQ_REQ = 0.5
MIN_SNVS_PER_COLONY = 500
# ----------------------------------------

def classify_snv(sub):
    """Temporal fate classification per-SNV."""
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


trajectory_files = sorted(TRAJ_DIR.glob("snv_frequency_trajectory_*.csv"))
all_summaries = []
monthly_trends = []
snv_records = []   # NEW: SNV-level export

print(f"Found {len(trajectory_files)} trajectory files.")

# ---------------- MAIN LOOP ----------------
for traj_file in trajectory_files:
    colony_id = traj_file.stem.split("_")[-1]
    print(f"\n▶ Processing colony {colony_id}")

    df = pd.read_csv(traj_file)
    df["Month"] = df["Month"].astype(str).str.title()
    df["Month"] = pd.Categorical(df["Month"], categories=MONTH_ORDER, ordered=True)
    df["major_freq"] = 1 - df["departure_from_polarized_consensus"]

    df["SNV_ID"] = (
        df["MAG_id"].astype(str) + "_" +
        df["contig_name"].astype(str) + "_" +
        df["pos_in_contig"].astype(str)
    )

    # ---------------- Identify private SNVs in May ----------------
    may_df = df[df["Month"] == "May"]
    private_ids = may_df.loc[
        may_df["departure_from_polarized_consensus"].isin([0.0, 1.0]),
        "SNV_ID"
    ].unique()

    n_private = len(private_ids)
    print(f"  Found {n_private:,} private SNVs (fixed in May).")
    if n_private == 0:
        continue

    df_private = df[df["SNV_ID"].isin(private_ids)].copy()

    # ---------------- Temporal classification ----------------
    classified = (
        df_private.groupby(["MAG_id", "contig_name", "pos_in_contig"], group_keys=False)
        .apply(classify_snv)
        .reset_index()
    )

    # attach MAG_id/contig/pos-level classification back onto each row
    df_private_merged = df_private.merge(
        classified[["MAG_id", "contig_name", "pos_in_contig", "status"]],
        on=["MAG_id", "contig_name", "pos_in_contig"],
        how="left"
    )

    # ---------------- SNV-LEVEL EXPORT (NEW) ----------------
    for _, row in df_private_merged.iterrows():
        snv_records.append({
            "Colony": colony_id,
            "MAG_id": row["MAG_id"],
            "SNV_ID": row["SNV_ID"],
            "contig_name": row["contig_name"],
            "pos_in_contig": row["pos_in_contig"],
            "Month": row["Month"],
            "major_freq": row["major_freq"],
            "temporal_status": row["status"]
        })

    # ---------------- MAG × Colony summary ----------------
    summary = (
        classified.groupby(["MAG_id","status"], as_index=False)
        .size()
        .pivot(index="MAG_id", columns="status", values="size")
        .fillna(0)
    )

    summary["Colony"] = colony_id
    num_cols = summary.select_dtypes(include=[np.number]).columns
    summary["Total_private"] = summary[num_cols].sum(axis=1)
    summary["Frac_preserved"] = summary.get("Preserved_private", 0) / summary["Total_private"]
    summary["Frac_disrupted"] = summary.get("Disrupted_private", 0) / summary["Total_private"]

    all_summaries.append(summary.reset_index())

    # ---------------- Month-by-month trend ----------------
    for mag, sub in df_private.groupby("MAG_id"):
        baseline_ids = sub["SNV_ID"].unique()
        for month in MONTH_ORDER:
            sub_month = sub[sub["Month"] == month]
            if sub_month.empty:
                continue

            freqs = sub_month.groupby("SNV_ID")["major_freq"].mean()
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


# ---------------- SAVE OUTPUTS ----------------

print("\nSaving outputs...")

# SNV-level trajectories
snv_df = pd.DataFrame(snv_records)
snv_out = OUT_DIR / "private_marker_SNV_trajectories_all_colonies.csv"
snv_df.to_csv(snv_out, index=False)
print(f"✔️ SNV trajectories saved: {snv_out}")

# Colony-level MAG summaries
all_summary_df = pd.concat(all_summaries, ignore_index=True)
summary_out = OUT_DIR / "private_marker_summary_all_colonies.csv"
all_summary_df.to_csv(summary_out, index=False)
print(f"✔️ Summary saved: {summary_out}")

# Monthly trends
monthly_df = pd.DataFrame(monthly_trends)
monthly_out = OUT_DIR / "private_marker_monthly_trends_all_colonies.csv"
monthly_df.to_csv(monthly_out, index=False)
print(f"✔️ Monthly trends saved: {monthly_out}")

# Identify MAGs shared across colonies with >=500 private SNVs
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
    cross_out = OUT_DIR / "MAGs_in_multiple_colonies_over500.csv"
    multi_colony.to_csv(cross_out, index=False)
    print(f"⚠️ MAGs shared across colonies (>=500 private SNVs): {cross_out}")
else:
    print("No MAGs found in multiple colonies with >=500 private SNVs.")

print("\n🎉 DONE")
