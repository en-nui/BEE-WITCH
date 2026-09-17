#!/usr/bin/env python3
"""
Align candidate MMseqs protein families and rank structural representatives.

Inputs
------
1) Candidate family summary CSV produced by the R recurrence pipeline.
2) Directory containing one FASTA per candidate family. Expected naming:
       <cluster_id>.fasta   or   <cluster_id>.fa   or   <cluster_id>.faa
3) MAFFT available on PATH (or pass --mafft-bin).

Outputs
-------
<outdir>/alignments/<cluster_id>.aligned.fasta
<outdir>/family_sequence_metrics.csv
<outdir>/structural_candidate_rank.csv
<outdir>/representative_sequences.fasta
<outdir>/run_summary.txt

The script does NOT redefine the statistical recurrence universe. It starts
from the candidate families supplied by the R pipeline and uses sequence-level
alignment only to assess within-family coherence and select an observed
representative sequence (medoid) for downstream structural work.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import median, mean
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--candidate-csv",
        default="/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel/parallel_with_CLUSTERED_PROTEINS/selection_gene_clusters/protein_representative_screen/candidate_family_summary_pre_alignment.csv",
        help="Candidate family summary CSV from the R pipeline.",
    )
    p.add_argument(
        "--member-fastas",
        default=(
            "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/"
            "ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel/"
            "parallel_with_CLUSTERED_PROTEINS/selection_gene_clusters/"
            "protein_representative_screen/member_fastas"
        ),
        help="Directory containing candidate-family FASTA files.",
    )
    p.add_argument(
        "--outdir",
        default=(
            "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/"
            "ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel/"
            "parallel_with_CLUSTERED_PROTEINS/selection_gene_clusters/"
            "protein_representative_screen/alignment_results"
        ),
        help="Output directory.",
    )
    p.add_argument(
        "--mafft-bin",
        default="mafft",
        help="MAFFT executable or absolute path.",
    )
    p.add_argument(
        "--max-family-size",
        type=int,
        default=0,
        help="Optional cap on sequences per family for alignment; 0 = no cap.",
    )
    return p.parse_args()


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name = None
    seq_parts: List[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).replace(" ", "")))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).replace(" ", "")))
    return records


def write_fasta(records: List[Tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def parse_mafft_stdout(text: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name = None
    seq_parts: List[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq_parts)))
            name = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts)))
    return records


def pairwise_metrics(seq_a: str, seq_b: str) -> Tuple[float, float]:
    """Return identity and aligned coverage.

    Identity is calculated over columns where both sequences have residues.
    Coverage is the number of jointly occupied aligned columns divided by the
    length of the shorter ungapped sequence.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("Aligned sequences have unequal lengths")

    compared = 0
    identical = 0
    nongap_a = 0
    nongap_b = 0

    for a, b in zip(seq_a, seq_b):
        if a != "-":
            nongap_a += 1
        if b != "-":
            nongap_b += 1
        if a != "-" and b != "-":
            compared += 1
            if a == b:
                identical += 1

    denom = min(nongap_a, nongap_b)
    coverage = compared / denom if denom else float("nan")
    identity = identical / compared if compared else float("nan")
    return identity, coverage


def summarize_alignment(records: List[Tuple[str, str]]) -> Dict[str, object]:
    if len(records) < 2:
        return {
            "n_sequences": len(records),
            "median_pairwise_identity": float("nan"),
            "mean_pairwise_identity": float("nan"),
            "min_pairwise_identity": float("nan"),
            "q25_pairwise_identity": float("nan"),
            "q75_pairwise_identity": float("nan"),
            "median_pairwise_coverage": float("nan"),
            "min_pairwise_coverage": float("nan"),
            "medoid": records[0][0] if records else "",
            "medoid_median_identity": float("nan"),
            "medoid_mean_identity": float("nan"),
        }

    identities: List[float] = []
    coverages: List[float] = []
    per_seq: Dict[str, List[float]] = {name: [] for name, _ in records}

    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            iden, cov = pairwise_metrics(records[i][1], records[j][1])
            if not math.isnan(iden):
                identities.append(iden)
                per_seq[records[i][0]].append(iden)
                per_seq[records[j][0]].append(iden)
            if not math.isnan(cov):
                coverages.append(cov)

    ranked = sorted(
        (
            median(vals) if vals else float("nan"),
            mean(vals) if vals else float("nan"),
            name,
        )
        for name, vals in per_seq.items()
    )
    medoid_median_identity, medoid_mean_identity, medoid = ranked[-1]

    def pct(vals: List[float], q: float) -> float:
        if not vals:
            return float("nan")
        vals = sorted(vals)
        idx = (len(vals) - 1) * q
        lo = math.floor(idx)
        hi = math.ceil(idx)
        if lo == hi:
            return vals[lo]
        return vals[lo] + (vals[hi] - vals[lo]) * (idx - lo)

    return {
        "n_sequences": len(records),
        "median_pairwise_identity": median(identities) if identities else float("nan"),
        "mean_pairwise_identity": mean(identities) if identities else float("nan"),
        "min_pairwise_identity": min(identities) if identities else float("nan"),
        "q25_pairwise_identity": pct(identities, 0.25),
        "q75_pairwise_identity": pct(identities, 0.75),
        "median_pairwise_coverage": median(coverages) if coverages else float("nan"),
        "min_pairwise_coverage": min(coverages) if coverages else float("nan"),
        "medoid": medoid,
        "medoid_median_identity": medoid_median_identity,
        "medoid_mean_identity": medoid_mean_identity,
    }


def quantile_rank(values: List[float], x: float) -> float:
    vals = [v for v in values if not math.isnan(v)]
    if not vals or math.isnan(x):
        return float("nan")
    vals.sort()
    # Empirical percentile rank, conservative and easy to interpret.
    return sum(v <= x for v in vals) / len(vals)


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def main() -> int:
    args = parse_args()
    candidate_csv = Path(args.candidate_csv)
    member_dir = Path(args.member_fastas)
    outdir = Path(args.outdir)
    align_dir = outdir / "alignments"
    input_dir = outdir / "candidate_inputs"

    if not candidate_csv.exists():
        raise FileNotFoundError(f"Candidate CSV not found: {candidate_csv}")
    if not member_dir.exists():
        raise FileNotFoundError(
            f"Member FASTA directory not found: {member_dir}\n"
            "This path must exist on the HPC filesystem where the script is run."
        )

    mafft = shutil.which(args.mafft_bin) if not os.path.isabs(args.mafft_bin) else args.mafft_bin
    if not mafft or not Path(mafft).exists():
        raise RuntimeError(
            f"MAFFT was not found: {args.mafft_bin}\n"
            "Load MAFFT or provide --mafft-bin /absolute/path/to/mafft."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    align_dir.mkdir(exist_ok=True)
    input_dir.mkdir(exist_ok=True)

    with candidate_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError("Candidate CSV contains no rows.")
    required = {"cluster_id", "n_genes", "n_MAGs", "n_genera", "n_colonies", "n_concurrent_intervals", "max_concurrent_genera", "max_concurrent_colonies"}
    missing = required - set(rows[0])
    if missing:
        raise RuntimeError(f"Candidate CSV is missing columns: {sorted(missing)}")

    summary_rows: List[Dict[str, object]] = []
    rep_records: List[Tuple[str, str]] = []
    missing_fastas: List[str] = []
    failed_alignments: List[Tuple[str, str]] = []

    fasta_extensions = [".fasta", ".fa", ".faa"]

    for row in rows:
        cluster_id = row["cluster_id"].strip()
        fasta = None
        for ext in fasta_extensions:
            candidate = member_dir / f"{cluster_id}{ext}"
            if candidate.exists():
                fasta = candidate
                break
        if fasta is None:
            missing_fastas.append(cluster_id)
            continue

        records = read_fasta(fasta)
        if not records:
            failed_alignments.append((cluster_id, "empty FASTA"))
            continue

        # Deduplicate sequence IDs while preserving first occurrence.
        seen = set()
        deduped: List[Tuple[str, str]] = []
        for name, seq in records:
            if name in seen:
                continue
            seen.add(name)
            deduped.append((name, seq))
        records = deduped

        if args.max_family_size and len(records) > args.max_family_size:
            # Deterministic truncation for diagnostic runs only.
            records = records[:args.max_family_size]

        input_fasta = input_dir / f"{cluster_id}.fasta"
        aligned_fasta = align_dir / f"{cluster_id}.aligned.fasta"
        write_fasta(records, input_fasta)

        if len(records) == 1:
            aligned = records
        else:
            try:
                proc = subprocess.run(
                    [mafft, "--auto", str(input_fasta)],
                    check=True,
                    text=True,
                    capture_output=True,
                )
            except subprocess.CalledProcessError as exc:
                failed_alignments.append((cluster_id, exc.stderr.strip() or "MAFFT failed"))
                continue
            aligned = parse_mafft_stdout(proc.stdout)

        if not aligned:
            failed_alignments.append((cluster_id, "MAFFT returned no alignment"))
            continue

        write_fasta(aligned, aligned_fasta)
        metrics = summarize_alignment(aligned)

        out: Dict[str, object] = dict(row)
        out.update(metrics)
        out["alignment_file"] = str(aligned_fasta)
        out["input_fasta"] = str(input_fasta)
        summary_rows.append(out)
        medoid_name = str(metrics["medoid"])
        medoid_seq = next(seq for name, seq in aligned if name == medoid_name)
        rep_records.append((f"{cluster_id}|medoid|{medoid_name}", medoid_seq.replace("-", "")))

    if not summary_rows:
        raise RuntimeError("No candidate families were successfully aligned.")

    # Rank structural candidates. Sequence coherence is NOT used to redefine
    # recurrence; it is a secondary visualization/representative-selection axis.
    for r in summary_rows:
        r["recurrence_strength"] = (
            3.0 * safe_float(str(r.get("max_concurrent_genera", "nan")))
            + 2.0 * safe_float(str(r.get("max_concurrent_MAGs", "nan")))
            + 1.0 * safe_float(str(r.get("max_concurrent_colonies", "nan")))
            + 2.0 * safe_float(str(r.get("n_concurrent_intervals", "nan")))
            + 0.5 * safe_float(str(r.get("total_0D_event_records", "nan")))
        )

    id_values = [safe_float(str(r["median_pairwise_identity"])) for r in summary_rows]
    cov_values = [safe_float(str(r["median_pairwise_coverage"])) for r in summary_rows]
    for r in summary_rows:
        r["family_identity_percentile"] = quantile_rank(id_values, safe_float(str(r["median_pairwise_identity"])))
        r["family_coverage_percentile"] = quantile_rank(cov_values, safe_float(str(r["median_pairwise_coverage"])))
        # Balanced secondary rank: recurrence remains dominant.
        r["structural_rank_score"] = (
            float(r["recurrence_strength"])
            + 5.0 * r["family_identity_percentile"]
            + 2.0 * r["family_coverage_percentile"]
        )

    summary_rows.sort(
        key=lambda r: (
            -safe_float(str(r["structural_rank_score"])),
            -safe_float(str(r["max_concurrent_genera"])),
            -safe_float(str(r["max_concurrent_MAGs"])),
            -safe_float(str(r["max_concurrent_colonies"])),
        )
    )

    fieldnames = list(summary_rows[0].keys())
    with (outdir / "family_sequence_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    with (outdir / "structural_candidate_rank.csv").open("w", newline="") as handle:
        rank_fields = [
            "cluster_id",
            "n_genes",
            "n_MAGs",
            "max_concurrent_MAGs",
            "n_genera",
            "max_concurrent_genera",
            "n_colonies",
            "max_concurrent_colonies",
            "n_concurrent_intervals",
            "dominant_concurrent_interval",
            "total_0D_event_records",
            "n_0D_genes",
            "n_0D_MAGs",
            "n_0D_genera",
            "n_0D_colonies",
            "independent_pnps_support",
            "median_pairwise_identity",
            "min_pairwise_identity",
            "median_pairwise_coverage",
            "min_pairwise_coverage",
            "medoid",
            "recurrence_strength",
            "family_identity_percentile",
            "family_coverage_percentile",
            "structural_rank_score",
            "input_fasta",
            "alignment_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=rank_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary_rows)

    write_fasta(rep_records, outdir / "representative_sequences.fasta")

    with (outdir / "run_summary.txt").open("w") as handle:
        handle.write(f"Candidate families in CSV: {len(rows)}\n")
        handle.write(f"Successfully aligned: {len(summary_rows)}\n")
        handle.write(f"Missing family FASTAs: {len(missing_fastas)}\n")
        handle.write(f"Failed alignments: {len(failed_alignments)}\n")
        handle.write(f"MAFFT: {mafft}\n")
        if missing_fastas:
            handle.write("\nMissing FASTAs:\n")
            handle.write("\n".join(missing_fastas) + "\n")
        if failed_alignments:
            handle.write("\nFailed alignments:\n")
            for cid, reason in failed_alignments:
                handle.write(f"{cid}\t{reason}\n")
        handle.write("\nTop candidate families:\n")
        for r in summary_rows[:10]:
            handle.write(
                f"{r['cluster_id']}\t"
                f"concurrent_genera={r['max_concurrent_genera']}\t"
                f"concurrent_MAGs={r['max_concurrent_MAGs']}\t"
                f"concurrent_colonies={r['max_concurrent_colonies']}\t"
                f"identity={r['median_pairwise_identity']:.4f}\t"
                f"coverage={r['median_pairwise_coverage']:.4f}\t"
                f"medoid={r['medoid']}\n"
            )

    print(f"Candidate families in CSV: {len(rows)}")
    print(f"Successfully aligned: {len(summary_rows)}")
    print(f"Missing family FASTAs: {len(missing_fastas)}")
    print(f"Failed alignments: {len(failed_alignments)}")
    print(f"Outputs: {outdir}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
