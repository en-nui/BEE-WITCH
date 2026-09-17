#!/usr/bin/env python3
"""
Map 0D Delta RAF events from the original anvi'o gene-coordinate system onto
MMseqs protein-family MAFFT alignments.

Key idea:
  aa_position is defined on the original gene sequence used by the pN/pS
  annotation pipeline. That sequence is translated here, then aligned to the
  corresponding MMseqs member protein before mapping into the family MAFFT
  alignment. We therefore do NOT assume that raw amino-acid numbering is
  identical between the original gene protein and MMseqs protein.

Inputs:
  --events CSV with at least:
      cluster_id, gene_key, MAG_id, colony_id, Month_t1, Month_t2,
      sequence_id, aa_position, ancestral_aa, derived_aa,
      unique_SNV_identifier, corresponding_gene_call
  --alignment-dir directory containing <cluster_id>.aligned.fasta
  --member-fasta eligible_bacterial_mmseqs_all_seqs.fasta
  --nt-gene-fasta original nt_gene_seq.fa

Outputs:
  candidate_mutation_alignment_mapping.csv
  candidate_mutation_alignment_summary.csv
  candidate_structure_positions.csv

The final mapping is:
  original gene AA position
      -> translated original gene protein
      -> MMseqs member residue position
      -> family alignment column
      -> medoid position
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import SeqIO

TARGETS = {
    "seq_fa06c9c3d284d238e540bc215d89e7f0a05ec050": "seq_fa06c9c3d284d238e540bc215d89e7f0a05ec050",
    "seq_68ab79ba456ef8093145192b017fd1656be5672b": "seq_aaafb6913940cadc446568042466eafdddd77e53",
    "seq_480042084cedb0f09b60f3020bacd9ef0f42fa8f": "seq_6fb9a660011902e06c96e91585d951e7a03dafc7",
    "seq_4fadd4da5fa72844bb0ef8c53a7565b9ec5708fd": "seq_81f39b249899cddab87cf42f4b8530150db9e9b9",
}

REQUIRED = {
    "cluster_id", "sequence_id", "aa_position", "MAG_id", "colony_id",
    "Month_t1", "Month_t2", "ancestral_aa", "derived_aa",
    "unique_SNV_identifier", "corresponding_gene_call",
}


def read_fasta(path: Path) -> Dict[str, str]:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}


def load_events(path: Path) -> List[dict]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("Event CSV has no header")
        missing = REQUIRED - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing required event columns: {sorted(missing)}")
        return list(reader)


def build_aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -5.0
    a.extend_gap_score = -0.5
    return a


def alignment_query_to_target_map(query: str, target: str, aligner: PairwiseAligner) -> Dict[int, int]:
    """Map 1-based positions in query to 1-based positions in target."""
    if not query or not target:
        return {}
    aln = aligner.align(query, target)[0]
    q_blocks, t_blocks = aln.aligned
    mapping: Dict[int, int] = {}
    for (q0, q1), (t0, t1) in zip(q_blocks, t_blocks):
        n = min(q1 - q0, t1 - t0)
        for i in range(n):
            mapping[q0 + i + 1] = t0 + i + 1
    return mapping


def member_position_to_alignment_column(aligned_seq: str, member_pos: int) -> int:
    count = 0
    for col, aa in enumerate(aligned_seq, start=1):
        if aa != "-":
            count += 1
            if count == member_pos:
                return col
    raise ValueError(
        f"Member residue {member_pos} exceeds alignment ungapped length {count}"
    )


def medoid_position(aligned_rep: str, alignment_col: int) -> Optional[int]:
    if alignment_col < 1 or alignment_col > len(aligned_rep):
        return None
    if aligned_rep[alignment_col - 1] == "-":
        return None
    return sum(aa != "-" for aa in aligned_rep[:alignment_col])


def translate_gene(nt_seq: str) -> str:
    nt_seq = nt_seq.upper().replace("-", "")
    if len(nt_seq) < 3:
        return ""
    # Original pipeline translates codons directly; retain complete codons.
    usable = len(nt_seq) - (len(nt_seq) % 3)
    return str(Seq(nt_seq[:usable]).translate(to_stop=False)).rstrip("*")


def normalize_id(x: str) -> str:
    return x.strip().split()[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True, type=Path)
    ap.add_argument("--alignment-dir", required=True, type=Path)
    ap.add_argument(
        "--member-fasta",
        type=Path,
        default=Path(
            "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/"
            "selection_pipeline_outputs_parallel/parallel_with_CLUSTERED_PROTEINS/"
            "selection_gene_clusters/eligible_bacterial_mmseqs_all_seqs.fasta"
        ),
    )
    ap.add_argument(
        "--nt-gene-fasta",
        type=Path,
        default=Path(
            "/home/robinch/projects/BEE-WITCH/40_internal_pnps_calculations_from_anvio/nt_gene_seq.fa"
        ),
    )
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    events = load_events(args.events)
    member_seqs = read_fasta(args.member_fasta)
    nt_genes = read_fasta(args.nt_gene_fasta)
    aligner = build_aligner()

    family_rows: Dict[str, List[dict]] = defaultdict(list)
    for row in events:
        if row["cluster_id"] in TARGETS:
            family_rows[row["cluster_id"]].append(row)

    all_mapped: List[dict] = []
    summary: List[dict] = []

    for cluster_id, medoid_id in TARGETS.items():
        aln_path = args.alignment_dir / f"{cluster_id}.aligned.fasta"
        rows = family_rows.get(cluster_id, [])
        if not aln_path.exists():
            summary.append({
                "cluster_id": cluster_id,
                "status": "missing_alignment",
                "n_alignment_sequences": 0,
                "alignment_length": 0,
                "medoid_id": medoid_id,
                "n_events": len(rows),
                "n_events_mapped": 0,
                "n_unique_alignment_positions": 0,
                "alignment_positions": "",
            })
            continue

        aln = read_fasta(aln_path)
        if medoid_id not in aln:
            raise ValueError(f"Medoid {medoid_id} not found in {aln_path}")
        rep_aln = aln[medoid_id]
        lengths = {len(s) for s in aln.values()}
        if len(lengths) != 1:
            raise ValueError(f"Unequal alignment lengths in {aln_path}: {sorted(lengths)}")

        events_mapped = 0
        event_columns: List[int] = []
        status_counts = defaultdict(int)

        # Cache translated original genes and original->member mappings.
        gene_proteins: Dict[str, str] = {}
        gene_to_member_map: Dict[Tuple[str, str], Dict[int, int]] = {}

        for row in rows:
            try:
                gene_call = normalize_id(row["corresponding_gene_call"])
                member_id = normalize_id(row["sequence_id"])
                orig_pos = int(row["aa_position"])

                if gene_call not in nt_genes:
                    raise ValueError(f"Gene call {gene_call} not found in nt gene FASTA")
                original_protein = gene_proteins.setdefault(
                    gene_call, translate_gene(nt_genes[gene_call])
                )
                if orig_pos < 1 or orig_pos > len(original_protein):
                    raise ValueError(
                        f"Original AA position {orig_pos} outside translated gene length {len(original_protein)}"
                    )

                if member_id not in member_seqs:
                    raise ValueError(f"MMseqs member {member_id} not found in member FASTA")
                member_protein = member_seqs[member_id].replace("-", "")
                if not member_protein:
                    raise ValueError(f"MMseqs member {member_id} has empty sequence")

                key = (gene_call, member_id)
                if key not in gene_to_member_map:
                    gene_to_member_map[key] = alignment_query_to_target_map(
                        original_protein, member_protein, aligner
                    )
                q2m = gene_to_member_map[key]
                if orig_pos not in q2m:
                    raise ValueError(
                        f"Could not map original protein position {orig_pos} from gene {gene_call} to member {member_id}"
                    )
                member_pos = q2m[orig_pos]

                # Validate the translated ancestral-state position when possible.
                original_residue = original_protein[orig_pos - 1]
                anc = (row.get("ancestral_aa") or "").upper()
                der = (row.get("derived_aa") or "").upper()
                state_check = "OK" if original_residue in {anc, der} else "WARNING_ORIGINAL_FASTA_RESIDUE"

                # Standard path: map the raw member sequence to its existing family alignment row.
                if member_id not in aln:
                    raise ValueError(f"sequence_id {member_id} not present in family alignment")
                member_aln = aln[member_id]
                nongap = member_aln.replace("-", "")
                if nongap != member_protein:
                    status_counts["alignment_sequence_differs_from_member_fasta"] += 1
                    # Still proceed by using the alignment sequence's residue numbering.
                    if member_pos > len(nongap):
                        raise ValueError(
                            f"Member position {member_pos} exceeds aligned ungapped length {len(nongap)}"
                        )
                if len(nongap) == 0:
                    raise ValueError("Family alignment row is all gaps")

                aln_col = member_position_to_alignment_column(member_aln, member_pos)
                rep_pos = medoid_position(rep_aln, aln_col)
                if rep_pos is None:
                    raise ValueError(
                        f"Alignment column {aln_col} corresponds to a gap in medoid {medoid_id}"
                    )

                member_residue = member_aln[aln_col - 1]
                medoid_residue = rep_aln[aln_col - 1]

                out = dict(row)
                out.update({
                    "original_protein_position": orig_pos,
                    "original_protein_residue": original_residue,
                    "member_protein_position": member_pos,
                    "alignment_column": aln_col,
                    "medoid_position": rep_pos,
                    "medoid_residue": medoid_residue,
                    "member_residue": member_residue,
                    "original_state_check": state_check,
                    "mapping_check": "OK",
                })
                all_mapped.append(out)
                events_mapped += 1
                event_columns.append(aln_col)

            except Exception as exc:
                out = dict(row)
                out.update({
                    "original_protein_position": row.get("aa_position", "NA"),
                    "original_protein_residue": "NA",
                    "member_protein_position": "NA",
                    "alignment_column": "NA",
                    "medoid_position": "NA",
                    "medoid_residue": "NA",
                    "member_residue": "NA",
                    "original_state_check": "NA",
                    "mapping_check": f"ERROR: {exc}",
                })
                all_mapped.append(out)
                status_counts["mapping_error"] += 1

        unique_cols = sorted(set(event_columns))
        summary.append({
            "cluster_id": cluster_id,
            "status": "mapped" if events_mapped == len(rows) else "partial",
            "n_alignment_sequences": len(aln),
            "alignment_length": len(rep_aln),
            "medoid_id": medoid_id,
            "n_events": len(rows),
            "n_events_mapped": events_mapped,
            "n_unique_alignment_positions": len(unique_cols),
            "alignment_positions": ";".join(map(str, unique_cols)),
            "mapping_warnings": ";".join(f"{k}={v}" for k, v in sorted(status_counts.items())),
        })

    event_out = args.outdir / "candidate_mutation_alignment_mapping.csv"
    if all_mapped:
        fields = list(all_mapped[0].keys())
        with event_out.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(all_mapped)
    else:
        event_out.write_text("")

    summary_out = args.outdir / "candidate_mutation_alignment_summary.csv"
    if summary:
        fields = sorted({k for r in summary for k in r})
        with summary_out.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(summary)

    by_family_pos: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for row in all_mapped:
        if row.get("medoid_position", "NA") != "NA":
            by_family_pos[(row["cluster_id"], str(row["medoid_position"]))].append(row)

    structure_rows = []
    for (cluster_id, rep_pos), rs in sorted(by_family_pos.items()):
        genera = sorted({r.get("Genus", "") for r in rs if r.get("Genus")})
        mags = sorted({r.get("MAG_id", "") for r in rs if r.get("MAG_id")})
        colonies = sorted({r.get("colony_id", "") for r in rs if r.get("colony_id")})
        intervals = sorted({f"{r['Month_t1']}->{r['Month_t2']}" for r in rs})
        substitutions = sorted({f"{r.get('ancestral_aa','?')}>{r.get('derived_aa','?')}" for r in rs})
        members = sorted({r.get("sequence_id", "") for r in rs if r.get("sequence_id")})
        structure_rows.append({
            "cluster_id": cluster_id,
            "medoid_position": rep_pos,
            "n_events": len(rs),
            "n_genera": len(genera),
            "genera": ";".join(genera),
            "n_MAGs": len(mags),
            "MAGs": ";".join(mags),
            "n_colonies": len(colonies),
            "colonies": ";".join(colonies),
            "intervals": ";".join(intervals),
            "substitutions": ";".join(substitutions),
            "members": ";".join(members),
        })

    structure_out = args.outdir / "candidate_structure_positions.csv"
    if structure_rows:
        fields = list(structure_rows[0].keys())
        with structure_out.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(structure_rows)
    else:
        structure_out.write_text("")

    print(f"Families targeted: {len(TARGETS)}")
    print(f"Event rows loaded: {len(events)}")
    print(f"Event rows in target families: {sum(len(v) for v in family_rows.values())}")
    print(f"Wrote: {event_out}")
    print(f"Wrote: {summary_out}")
    print(f"Wrote: {structure_out}")


if __name__ == "__main__":
    main()
