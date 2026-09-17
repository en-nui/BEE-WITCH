#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:46:30 2026

@author: robinch
"""

from pathlib import Path
from itertools import combinations

import numpy as np
from Bio.PDB import PDBParser


BASE = Path("/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel/parallel_with_CLUSTERED_PROTEINS/selection_gene_clusters/protein_representative_screen/a59cc_medoid_prediction/")

PATTERN = (
    "seq_b95da49716b9118ad697139dc9799dd86adc97e9_"
    "unrelaxed_rank_{rank:03d}_alphafold2_ptm_model_"
)

MODEL_SUFFIXES = {
    1: "2_seed_000.pdb",
    2: "3_seed_000.pdb",
    3: "5_seed_000.pdb",
    4: "4_seed_000.pdb",
    5: "1_seed_000.pdb",
}

FOCAL = [33, 39, 52, 77]
CLUSTER = [33, 39, 52]
WALKER_A = range(43, 51)


def heavy_atoms(residue):
    return [
        atom
        for atom in residue.get_atoms()
        if atom.element != "H"
    ]


def min_heavy_distance(r1, r2):
    return min(
        a - b
        for a in heavy_atoms(r1)
        for b in heavy_atoms(r2)
    )


print("rank\t33-39\t33-52\t39-52\t33-77\t39-77\t52-77\tcluster_diam\tL52-WA")


for rank in range(1, 6):

    pdb = BASE / (
        PATTERN.format(rank=rank)
        + MODEL_SUFFIXES[rank]
    )

    structure = PDBParser(QUIET=True).get_structure(
        f"rank{rank}",
        pdb
    )

    model = next(structure.get_models())
    chain = model["A"]

    residues = {
        p: chain[(" ", p, " ")]
        for p in FOCAL
    }

    pairwise = {}

    for a, b in combinations(FOCAL, 2):
        pairwise[(a, b)] = (
            residues[a]["CA"] -
            residues[b]["CA"]
        )

    cluster_diameter = max(
        pairwise[(a, b)]
        for a, b in combinations(CLUSTER, 2)
    )

    walker_residues = [
        chain[(" ", p, " ")]
        for p in WALKER_A
    ]

    l52_walker = min(
        min_heavy_distance(
            residues[52],
            walker
        )
        for walker in walker_residues
    )

    print(
        f"{rank}\t"
        f"{pairwise[(33,39)]:.2f}\t"
        f"{pairwise[(33,52)]:.2f}\t"
        f"{pairwise[(39,52)]:.2f}\t"
        f"{pairwise[(33,77)]:.2f}\t"
        f"{pairwise[(39,77)]:.2f}\t"
        f"{pairwise[(52,77)]:.2f}\t"
        f"{cluster_diameter:.2f}\t"
        f"{l52_walker:.2f}"
    )