#!/usr/bin/env python3

"""
PRIVATE_ALLELES_TEMPORAL_RESIDENT_ORIENTED.py
=============================================

Track the fate of the resident genomic background through time.

For each MAG x colony:

1. Identify the first observed timepoint for that MAG.
2. Identify SNVs fixed at that initial timepoint:
       departure_from_polarized_consensus ~= 0 or ~= 1
3. Orient every marker toward the allele that was fixed in the
   resident population at the initial timepoint.

   Therefore:

       resident_freq(t0) = 1

   for every private marker, regardless of whether the original
   polarized departure frequency was 0 or 1.

4. Track each resident marker through time.

5. Classify each private marker as:

   Preserved_private:
       resident allele frequency >= 0.8 in >= 80% of callable months
       AND
       resident allele frequency never falls below 0.5

   Disrupted_private:
       does not satisfy the above criteria

6. Generate MAG x colony x month measures of resident-background
   retention and turnover.

7. Generate MAG x colony x transition summaries that can later be
   merged directly with allele-sweep timepoint pairs.

IMPORTANT:
Missing/un-callable private markers are NOT automatically treated as
disrupted. Callability is reported separately from allele-frequency
change.
"""

import numpy as np
import pandas as pd
from pathlib import Path


# ============================================================
# CONFIG
# ============================================================

TRAJ_DIR = Path(
    "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/"
    "ALLSPECIES_intrapop_diversity/snv_frequency_trajectories"
)

# Use a new directory so we do not overwrite the old private-marker results.
OUT_DIR = (
    TRAJ_DIR.parent /
    "private_marker_results_resident_oriented"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


MONTH_ORDER = [
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "January",
    "February"
]

MONTH_INDEX = {
    month: i
    for i, month in enumerate(MONTH_ORDER)
}


# ------------------------------------------------------------
# Private-marker thresholds
# ------------------------------------------------------------

# Preserved marker:
# frequency >= 0.8 for >=80% of callable months
PRESERVED_FREQ = 0.80
PRESERVED_PROP = 0.80

# If resident allele ever falls below this,
# it cannot be classified as preserved.
MIN_PRESERVED_FREQ = 0.50

# Strong loss of original resident allele.
# This is NOT used in the binary Preserved/Disrupted classification.
# It is useful for the evolution/replacement analysis.
LOST_FREQ = 0.20

# Existing QC threshold for cross-colony MAG summaries.
MIN_SNVS_PER_COLONY = 500

# Floating-point tolerance when deciding whether t0 frequency
# is fixed at exactly 0 or 1.
FIXED_ATOL = 1e-8


# ------------------------------------------------------------
# Optional SNV metadata to retain if present
# ------------------------------------------------------------

# Keeping degeneracy information here will let us perform
# the 4D-private-marker sensitivity analysis later if the
# trajectory files already contain this information.
OPTIONAL_METADATA = [
    "degeneracy",
    "site_degeneracy",
    "codon_degeneracy",
    "corresponding_gene_call",
    "gene_callers_id"
]


# ============================================================
# HELPERS
# ============================================================

def first_non_null(x):
    """
    Return first non-null value in a pandas Series.
    Used for metadata when collapsing duplicate SNV-month rows.
    """

    x = x.dropna()

    if len(x) == 0:
        return np.nan

    return x.iloc[0]


def normalize_months(df):
    """
    Standardize month names and confirm all months have an
    explicitly defined chronological position.
    """

    df = df.copy()

    df["Month"] = (
        df["Month"]
        .astype(str)
        .str.strip()
        .str.title()
    )

    observed = set(
        df.loc[
            df["Month"].notna(),
            "Month"
        ]
    )

    unknown = sorted(
        observed - set(MONTH_ORDER)
    )

    if unknown:
        raise ValueError(
            "Found months not present in MONTH_ORDER: "
            + ", ".join(unknown)
            + "\nAdd them to MONTH_ORDER before running."
        )

    df["month_index"] = df["Month"].map(MONTH_INDEX)

    return df


def build_snv_id(df):
    """
    Construct stable SNV identifier.
    """

    return (
        df["MAG_id"].astype(str)
        + "_"
        + df["contig_name"].astype(str)
        + "_"
        + df["pos_in_contig"].astype(str)
    )


def collapse_snv_month_rows(df):
    """
    Collapse any duplicate observations of the same SNV within
    a MAG x month.

    departure_from_polarized_consensus is averaged.

    Optional metadata fields are retained using their first
    non-null value.
    """

    metadata_cols = [
        c
        for c in OPTIONAL_METADATA
        if c in df.columns
    ]

    group_cols = [
        "MAG_id",
        "SNV_ID",
        "contig_name",
        "pos_in_contig",
        "Month",
        "month_index"
    ]

    agg_dict = {
        "departure_from_polarized_consensus": "mean"
    }

    for col in metadata_cols:
        agg_dict[col] = first_non_null

    out = (
        df
        .groupby(
            group_cols,
            as_index=False,
            dropna=False
        )
        .agg(agg_dict)
    )

    return out


def get_sampled_months(mag_df):
    """
    Return observed months for a MAG in chronological order.
    """

    observed = set(
        mag_df["Month"]
        .dropna()
        .astype(str)
    )

    return [
        month
        for month in MONTH_ORDER
        if month in observed
    ]


def classify_private_markers(
    private_df,
    sampled_months
):
    """
    One row per private SNV.

    Classification is based only on months in which that marker
    was callable.

    Missing markers are therefore not automatically interpreted
    as disrupted.
    """

    records = []

    n_expected_months = len(sampled_months)

    for snv_id, sub in private_df.groupby(
        "SNV_ID",
        sort=False
    ):

        freqs = (
            sub["resident_freq"]
            .dropna()
            .astype(float)
        )

        n_callable = len(freqs)

        if n_callable == 0:

            n_high = 0
            prop_high = np.nan
            min_freq = np.nan
            mean_freq = np.nan
            median_freq = np.nan
            status = "Missing"

        else:

            n_high = int(
                (freqs >= PRESERVED_FREQ).sum()
            )

            prop_high = (
                n_high /
                n_callable
            )

            min_freq = freqs.min()
            mean_freq = freqs.mean()
            median_freq = freqs.median()

            if (
                prop_high >= PRESERVED_PROP
                and
                min_freq >= MIN_PRESERVED_FREQ
            ):
                status = "Preserved_private"

            else:
                status = "Disrupted_private"

        prop_months_callable = (
            n_callable / n_expected_months
            if n_expected_months > 0
            else np.nan
        )

        records.append({
            "SNV_ID": snv_id,
            "n_months_expected": n_expected_months,
            "n_months_callable": n_callable,
            "prop_months_callable": prop_months_callable,
            "n_months_high": n_high,
            "prop_high": prop_high,
            "min_resident_freq": min_freq,
            "mean_resident_freq": mean_freq,
            "median_resident_freq": median_freq,
            "status": status
        })

    return pd.DataFrame(records)


def summarize_private_month(
    sub,
    n_private_baseline
):
    """
    Summarize the resident genomic background at one month.

    All proportions describing allele state use only callable
    private markers as the denominator.

    Callability itself is reported separately.
    """

    if sub.empty:

        return {
            "n_private_callable": 0,
            "callable_fraction": 0.0,

            "n_retained_ge_0p8": 0,
            "frac_retained_callable": np.nan,

            "n_below_0p5": 0,
            "frac_below_0p5_callable": np.nan,

            "n_lost_le_0p2": 0,
            "frac_lost_callable": np.nan,

            "mean_resident_freq": np.nan,
            "median_resident_freq": np.nan,

            "background_turnover_mean": np.nan,
            "background_turnover_median": np.nan
        }

    # There should normally be one value per SNV/month after collapse,
    # but group again defensively.
    freqs = (
        sub
        .groupby("SNV_ID")["resident_freq"]
        .mean()
        .dropna()
    )

    n_callable = len(freqs)

    if n_callable == 0:

        return {
            "n_private_callable": 0,
            "callable_fraction": 0.0,

            "n_retained_ge_0p8": 0,
            "frac_retained_callable": np.nan,

            "n_below_0p5": 0,
            "frac_below_0p5_callable": np.nan,

            "n_lost_le_0p2": 0,
            "frac_lost_callable": np.nan,

            "mean_resident_freq": np.nan,
            "median_resident_freq": np.nan,

            "background_turnover_mean": np.nan,
            "background_turnover_median": np.nan
        }

    n_retained = int(
        (freqs >= PRESERVED_FREQ).sum()
    )

    n_below_half = int(
        (freqs < MIN_PRESERVED_FREQ).sum()
    )

    n_lost = int(
        (freqs <= LOST_FREQ).sum()
    )

    mean_freq = freqs.mean()
    median_freq = freqs.median()

    return {
        "n_private_callable": n_callable,

        "callable_fraction":
            n_callable / n_private_baseline,

        "n_retained_ge_0p8":
            n_retained,

        "frac_retained_callable":
            n_retained / n_callable,

        "n_below_0p5":
            n_below_half,

        "frac_below_0p5_callable":
            n_below_half / n_callable,

        "n_lost_le_0p2":
            n_lost,

        "frac_lost_callable":
            n_lost / n_callable,

        "mean_resident_freq":
            mean_freq,

        "median_resident_freq":
            median_freq,

        # 0 = original resident background fully retained
        # 1 = original resident background entirely lost
        "background_turnover_mean":
            1.0 - mean_freq,

        "background_turnover_median":
            1.0 - median_freq
    }


# ============================================================
# FIND TRAJECTORY FILES
# ============================================================

trajectory_files = sorted(
    TRAJ_DIR.glob(
        "snv_frequency_trajectory_*.csv"
    )
)

if len(trajectory_files) == 0:
    raise FileNotFoundError(
        f"No trajectory files found in:\n{TRAJ_DIR}"
    )


print(
    f"Found {len(trajectory_files)} trajectory files."
)


# ============================================================
# OUTPUT CONTAINERS
# ============================================================

all_snv_trajectories = []
all_marker_classifications = []
all_mag_summaries = []
all_monthly_summaries = []

qc_baseline_zero = 0
qc_baseline_one = 0


# ============================================================
# MAIN LOOP: ONE FILE / COLONY AT A TIME
# ============================================================

for traj_file in trajectory_files:

    colony_id = (
        traj_file
        .stem
        .split("_")[-1]
    )

    print(
        f"\nProcessing colony {colony_id}"
    )

    df = pd.read_csv(
        traj_file,
        low_memory=False
    )


    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

    required_cols = [
        "MAG_id",
        "Month",
        "contig_name",
        "pos_in_contig",
        "departure_from_polarized_consensus"
    ]

    missing_cols = [
        c
        for c in required_cols
        if c not in df.columns
    ]

    if missing_cols:

        raise ValueError(
            f"{traj_file.name} is missing required columns: "
            + ", ".join(missing_cols)
        )


    # --------------------------------------------------------
    # Standardize data
    # --------------------------------------------------------

    df["MAG_id"] = (
        df["MAG_id"]
        .astype(str)
    )

    df[
        "departure_from_polarized_consensus"
    ] = pd.to_numeric(
        df[
            "departure_from_polarized_consensus"
        ],
        errors="coerce"
    )

    df = normalize_months(df)

    df["SNV_ID"] = build_snv_id(df)


    # Check allele-frequency bounds.
    bad_freq = (
        df[
            "departure_from_polarized_consensus"
        ].notna()
        &
        (
            (
                df[
                    "departure_from_polarized_consensus"
                ] < 0
            )
            |
            (
                df[
                    "departure_from_polarized_consensus"
                ] > 1
            )
        )
    )

    if bad_freq.any():

        raise ValueError(
            f"{traj_file.name} contains "
            f"{bad_freq.sum()} allele frequencies outside [0,1]."
        )


    # Collapse any duplicate SNV x month entries.
    df_collapsed = collapse_snv_month_rows(df)


    # --------------------------------------------------------
    # Analyze each MAG separately
    # --------------------------------------------------------

    for mag_id, mag_df in df_collapsed.groupby(
        "MAG_id",
        sort=False
    ):

        sampled_months = get_sampled_months(
            mag_df
        )

        if len(sampled_months) == 0:
            continue

        baseline_month = sampled_months[0]


        # ----------------------------------------------------
        # Identify private markers at first observed month
        # ----------------------------------------------------

        baseline_df = mag_df[
            mag_df["Month"] == baseline_month
        ].copy()

        baseline_freq = (
            baseline_df[
                "departure_from_polarized_consensus"
            ]
        )

        fixed_zero = np.isclose(
            baseline_freq,
            0.0,
            atol=FIXED_ATOL,
            rtol=0
        )

        fixed_one = np.isclose(
            baseline_freq,
            1.0,
            atol=FIXED_ATOL,
            rtol=0
        )

        baseline_private = (
            baseline_df[
                fixed_zero |
                fixed_one
            ]
            .copy()
        )

        if baseline_private.empty:
            continue


        # One baseline state per SNV.
        baseline_private = (
            baseline_private
            .drop_duplicates(
                subset="SNV_ID"
            )
        )

        baseline_private[
            "baseline_departure"
        ] = baseline_private[
            "departure_from_polarized_consensus"
        ]


        n_private = len(
            baseline_private
        )

        n_fixed_zero = int(
            np.isclose(
                baseline_private[
                    "baseline_departure"
                ],
                0.0,
                atol=FIXED_ATOL,
                rtol=0
            ).sum()
        )

        n_fixed_one = int(
            np.isclose(
                baseline_private[
                    "baseline_departure"
                ],
                1.0,
                atol=FIXED_ATOL,
                rtol=0
            ).sum()
        )

        qc_baseline_zero += n_fixed_zero
        qc_baseline_one += n_fixed_one


        print(
            f"  {mag_id}: "
            f"{n_private:,} private markers "
            f"at {baseline_month} "
            f"(departure=0: {n_fixed_zero:,}; "
            f"departure=1: {n_fixed_one:,})"
        )


        # ----------------------------------------------------
        # Attach baseline state to all subsequent observations
        # ----------------------------------------------------

        baseline_orientation = (
            baseline_private[
                [
                    "SNV_ID",
                    "baseline_departure"
                ]
            ]
        )

        private_df = mag_df.merge(
            baseline_orientation,
            on="SNV_ID",
            how="inner"
        )


        # ----------------------------------------------------
        # Orient toward the RESIDENT allele
        # ----------------------------------------------------

        # If departure(t0) = 0:
        #     resident allele = consensus allele
        #     resident_freq = 1 - departure
        #
        # If departure(t0) = 1:
        #     resident allele = departure allele
        #     resident_freq = departure

        private_df[
            "resident_freq"
        ] = np.where(
            private_df[
                "baseline_departure"
            ] < 0.5,

            1.0 -
            private_df[
                "departure_from_polarized_consensus"
            ],

            private_df[
                "departure_from_polarized_consensus"
            ]
        )


        # Numerical protection only.
        private_df[
            "resident_freq"
        ] = private_df[
            "resident_freq"
        ].clip(
            lower=0,
            upper=1
        )


        private_df[
            "resident_allele_orientation"
        ] = np.where(
            private_df[
                "baseline_departure"
            ] < 0.5,
            "consensus_at_t0",
            "departure_at_t0"
        )


        private_df[
            "baseline_month"
        ] = baseline_month


        # ----------------------------------------------------
        # QC:
        # every callable private marker at t0 MUST now equal 1
        # ----------------------------------------------------

        baseline_oriented = private_df[
            private_df["Month"] ==
            baseline_month
        ]["resident_freq"].dropna()

        if not np.allclose(
            baseline_oriented,
            1.0,
            atol=FIXED_ATOL,
            rtol=0
        ):

            raise RuntimeError(
                f"Resident allele orientation failed for "
                f"{mag_id} in colony {colony_id}."
            )


        # ----------------------------------------------------
        # Marker-level temporal classification
        # ----------------------------------------------------

        classified = classify_private_markers(
            private_df,
            sampled_months
        )


        # Add marker genomic metadata.
        metadata_cols = [
            c
            for c in OPTIONAL_METADATA
            if c in baseline_private.columns
        ]

        marker_meta_cols = [
            "SNV_ID",
            "contig_name",
            "pos_in_contig",
            "baseline_departure"
        ] + metadata_cols

        marker_meta = (
            baseline_private[
                marker_meta_cols
            ]
            .drop_duplicates(
                subset="SNV_ID"
            )
        )

        classified = classified.merge(
            marker_meta,
            on="SNV_ID",
            how="left"
        )

        classified["Colony"] = colony_id
        classified["MAG_id"] = mag_id
        classified[
            "baseline_month"
        ] = baseline_month


        all_marker_classifications.append(
            classified
        )


        # ----------------------------------------------------
        # Attach classification to SNV-level trajectories
        # ----------------------------------------------------

        private_export = private_df.merge(
            classified[
                [
                    "SNV_ID",
                    "status",
                    "n_months_expected",
                    "n_months_callable",
                    "prop_months_callable",
                    "prop_high",
                    "min_resident_freq"
                ]
            ],
            on="SNV_ID",
            how="left"
        )

        private_export[
            "Colony"
        ] = colony_id

        private_export[
            "Trajectory_file"
        ] = traj_file.name

        all_snv_trajectories.append(
            private_export
        )


        # ----------------------------------------------------
        # MAG x colony summary
        # ----------------------------------------------------

        status_counts = (
            classified[
                "status"
            ]
            .value_counts()
        )

        n_preserved = int(
            status_counts.get(
                "Preserved_private",
                0
            )
        )

        n_disrupted = int(
            status_counts.get(
                "Disrupted_private",
                0
            )
        )

        n_missing = int(
            status_counts.get(
                "Missing",
                0
            )
        )


        all_mag_summaries.append({
            "Colony": colony_id,
            "MAG_id": mag_id,

            "baseline_month":
                baseline_month,

            "n_months_sampled":
                len(sampled_months),

            "sampled_months":
                ";".join(sampled_months),

            "Total_private":
                n_private,

            "Private_fixed_departure0":
                n_fixed_zero,

            "Private_fixed_departure1":
                n_fixed_one,

            "Preserved_private":
                n_preserved,

            "Disrupted_private":
                n_disrupted,

            "Missing_private":
                n_missing,

            "Frac_preserved":
                (
                    n_preserved /
                    n_private
                    if n_private > 0
                    else np.nan
                ),

            "Frac_disrupted":
                (
                    n_disrupted /
                    n_private
                    if n_private > 0
                    else np.nan
                ),

            "median_marker_callability":
                classified[
                    "prop_months_callable"
                ].median()
        })


        # ----------------------------------------------------
        # MAG x colony x month summaries
        # ----------------------------------------------------

        for month in sampled_months:

            month_sub = private_df[
                private_df["Month"] ==
                month
            ]

            metrics = summarize_private_month(
                month_sub,
                n_private
            )

            row = {
                "Colony": colony_id,
                "MAG_id": mag_id,
                "baseline_month":
                    baseline_month,
                "Month": month,
                "month_index":
                    MONTH_INDEX[month],
                "n_private_baseline":
                    n_private
            }

            row.update(metrics)

            all_monthly_summaries.append(
                row
            )


# ============================================================
# CONCATENATE OUTPUT TABLES
# ============================================================

print(
    "\nCombining results..."
)


if not all_mag_summaries:

    raise RuntimeError(
        "No private-marker trajectories were generated."
    )


summary_df = pd.DataFrame(
    all_mag_summaries
)

monthly_df = pd.DataFrame(
    all_monthly_summaries
)


if all_snv_trajectories:

    snv_df = pd.concat(
        all_snv_trajectories,
        ignore_index=True,
        sort=False
    )

else:

    snv_df = pd.DataFrame()


if all_marker_classifications:

    marker_df = pd.concat(
        all_marker_classifications,
        ignore_index=True,
        sort=False
    )

else:

    marker_df = pd.DataFrame()


# ============================================================
# SORT TABLES
# ============================================================

summary_df = summary_df.sort_values(
    [
        "Colony",
        "MAG_id"
    ]
)


monthly_df = monthly_df.sort_values(
    [
        "Colony",
        "MAG_id",
        "month_index"
    ]
)


if not snv_df.empty:

    snv_df = snv_df.sort_values(
        [
            "Colony",
            "MAG_id",
            "month_index",
            "contig_name",
            "pos_in_contig"
        ]
    )


if not marker_df.empty:

    marker_df = marker_df.sort_values(
        [
            "Colony",
            "MAG_id",
            "contig_name",
            "pos_in_contig"
        ]
    )


# ============================================================
# BUILD TRANSITION TABLE
# ============================================================

transition_records = []


for (
    colony,
    mag_id
), sub in monthly_df.groupby(
    [
        "Colony",
        "MAG_id"
    ],
    sort=False
):

    sub = (
        sub
        .sort_values(
            "month_index"
        )
        .reset_index(
            drop=True
        )
    )


    if len(sub) < 2:
        continue


    for i in range(
        1,
        len(sub)
    ):

        t1 = sub.iloc[i - 1]
        t2 = sub.iloc[i]


        turnover_t1 = t1[
            "background_turnover_median"
        ]

        turnover_t2 = t2[
            "background_turnover_median"
        ]


        retained_t1 = t1[
            "frac_retained_callable"
        ]

        retained_t2 = t2[
            "frac_retained_callable"
        ]


        lost_t1 = t1[
            "frac_lost_callable"
        ]

        lost_t2 = t2[
            "frac_lost_callable"
        ]


        transition_records.append({
            "Colony":
                colony,

            "MAG_id":
                mag_id,

            "baseline_month":
                t1["baseline_month"],

            "Month_t1":
                t1["Month"],

            "Month_t2":
                t2["Month"],

            "n_private_baseline":
                t1["n_private_baseline"],

            "n_private_callable_t1":
                t1["n_private_callable"],

            "n_private_callable_t2":
                t2["n_private_callable"],

            "callable_fraction_t1":
                t1["callable_fraction"],

            "callable_fraction_t2":
                t2["callable_fraction"],

            "min_callable_fraction_pair":
                np.nanmin([
                    t1["callable_fraction"],
                    t2["callable_fraction"]
                ]),

            "frac_retained_t1":
                retained_t1,

            "frac_retained_t2":
                retained_t2,

            "delta_frac_retained":
                (
                    retained_t2 -
                    retained_t1
                    if (
                        pd.notna(retained_t1)
                        and
                        pd.notna(retained_t2)
                    )
                    else np.nan
                ),

            "frac_lost_t1":
                lost_t1,

            "frac_lost_t2":
                lost_t2,

            "delta_frac_lost":
                (
                    lost_t2 -
                    lost_t1
                    if (
                        pd.notna(lost_t1)
                        and
                        pd.notna(lost_t2)
                    )
                    else np.nan
                ),

            "background_turnover_t1":
                turnover_t1,

            "background_turnover_t2":
                turnover_t2,

            # Positive =
            # loss of original resident background
            # between the two timepoints.
            "delta_background_turnover":
                (
                    turnover_t2 -
                    turnover_t1
                    if (
                        pd.notna(turnover_t1)
                        and
                        pd.notna(turnover_t2)
                    )
                    else np.nan
                )
        })


transition_df = pd.DataFrame(
    transition_records
)


# ============================================================
# MAG-LEVEL MAXIMUM TURNOVER SUMMARY
# ============================================================

turnover_summary = (
    monthly_df
    .groupby(
        [
            "Colony",
            "MAG_id"
        ],
        as_index=False
    )
    .agg(
        max_background_turnover=(
            "background_turnover_median",
            "max"
        ),

        max_frac_lost=(
            "frac_lost_callable",
            "max"
        ),

        min_frac_retained=(
            "frac_retained_callable",
            "min"
        ),

        median_monthly_callability=(
            "callable_fraction",
            "median"
        )
    )
)


summary_df = summary_df.merge(
    turnover_summary,
    on=[
        "Colony",
        "MAG_id"
    ],
    how="left"
)


# ============================================================
# CROSS-COLONY MAG QC
# ============================================================

mag_colony_counts = (
    summary_df[
        [
            "MAG_id",
            "Colony",
            "Total_private"
        ]
    ]
    .copy()
)


mag_colony_counts = mag_colony_counts[
    mag_colony_counts[
        "Total_private"
    ] >= MIN_SNVS_PER_COLONY
]


multi_colony = (
    mag_colony_counts
    .groupby(
        "MAG_id",
        group_keys=False
    )
    .filter(
        lambda x:
        x["Colony"].nunique() > 1
    )
)


# ============================================================
# SAVE
# ============================================================

print(
    "\nSaving outputs..."
)


# ------------------------------------------------------------
# SNV-level trajectories
# ------------------------------------------------------------

snv_out = (
    OUT_DIR /
    "private_marker_SNV_trajectories_all_colonies.csv"
)

snv_df.to_csv(
    snv_out,
    index=False
)

print(
    f"SNV trajectories:\n  {snv_out}"
)


# ------------------------------------------------------------
# One row per private marker
# ------------------------------------------------------------

marker_out = (
    OUT_DIR /
    "private_marker_SNV_classification_all_colonies.csv"
)

marker_df.to_csv(
    marker_out,
    index=False
)

print(
    f"SNV classifications:\n  {marker_out}"
)


# ------------------------------------------------------------
# One row per MAG x colony
# ------------------------------------------------------------

summary_out = (
    OUT_DIR /
    "private_marker_summary_all_colonies.csv"
)

summary_df.to_csv(
    summary_out,
    index=False
)

print(
    f"MAG summaries:\n  {summary_out}"
)


# ------------------------------------------------------------
# One row per MAG x colony x month
# ------------------------------------------------------------

monthly_out = (
    OUT_DIR /
    "private_marker_monthly_trends_all_colonies.csv"
)

monthly_df.to_csv(
    monthly_out,
    index=False
)

print(
    f"Monthly summaries:\n  {monthly_out}"
)


# ------------------------------------------------------------
# One row per MAG x colony x consecutive timepoint pair
# ------------------------------------------------------------

transition_out = (
    OUT_DIR /
    "private_marker_transition_summary_all_colonies.csv"
)

transition_df.to_csv(
    transition_out,
    index=False
)

print(
    f"Transition summaries:\n  {transition_out}"
)


# ------------------------------------------------------------
# Cross-colony MAGs with >=500 private markers
# ------------------------------------------------------------

cross_out = (
    OUT_DIR /
    "MAGs_in_multiple_colonies_over500.csv"
)

multi_colony.to_csv(
    cross_out,
    index=False
)

print(
    f"Cross-colony MAG QC:\n  {cross_out}"
)


# ============================================================
# FINAL QC
# ============================================================

print(
    "\n"
    "============================================================"
)

print(
    "PRIVATE-MARKER QC"
)

print(
    "============================================================"
)


print(
    f"MAG x colony trajectories: "
    f"{len(summary_df):,}"
)


print(
    f"Total baseline private markers: "
    f"{summary_df['Total_private'].sum():,}"
)


print(
    f"Median private markers per trajectory: "
    f"{summary_df['Total_private'].median():.1f}"
)


print(
    f"Baseline markers originally at departure = 0: "
    f"{qc_baseline_zero:,}"
)


print(
    f"Baseline markers originally at departure = 1: "
    f"{qc_baseline_one:,}"
)


if (
    qc_baseline_zero +
    qc_baseline_one
) > 0:

    pct_one = (
        qc_baseline_one /
        (
            qc_baseline_zero +
            qc_baseline_one
        )
        * 100
    )

    print(
        f"Percent requiring departure=1 orientation: "
        f"{pct_one:.2f}%"
    )


print(
    f"Median monthly private-marker callability: "
    f"{monthly_df['callable_fraction'].median():.3f}"
)


print(
    f"Median fraction preserved per MAG trajectory: "
    f"{summary_df['Frac_preserved'].median():.3f}"
)


print(
    f"Median maximum background turnover: "
    f"{summary_df['max_background_turnover'].median():.3f}"
)


print(
    f"MAGs represented in >1 colony with >= "
    f"{MIN_SNVS_PER_COLONY} private markers: "
    f"{multi_colony['MAG_id'].nunique():,}"
)


print(
    "============================================================"
)

print(
    "\nDONE"
)
