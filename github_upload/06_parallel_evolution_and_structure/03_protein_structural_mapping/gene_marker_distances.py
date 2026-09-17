#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
import Bio.PDB

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

CIF_FILE = Path("/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/selection_pipeline_outputs_parallel/parallel_with_CLUSTERED_PROTEINS/selection_gene_clusters/protein_representative_screen/alphafold_mdb_files/AF-A0A316J9J2-F1-model_v6.cif")
CHAIN_ID = "A"

# PT2 medoid -> AlphaFold homolog
MAPPING = {
    33: 37,
    39: 43,
    52: 56,
    77: 81,
}

AF_POSITIONS = list(MAPPING.values())
# Define the 20 standard 3-letter amino acid codes
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", 
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}

# Your loop
if residue.get_resname() not in STANDARD_AA:
    print(f"Non-standard residue or ligand found: {residue.get_resname()}")


# ------------------------------------------------------------
# Load structure
# ------------------------------------------------------------

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("a59cc_AF", CIF_FILE)

model = next(structure.get_models())

print("Available chains:", [c.id for c in model])

chain = model[CHAIN_ID]


# ------------------------------------------------------------
# Validate focal residues
# ------------------------------------------------------------

selected = {}

print("\nAlphaFold structural coordinates")
print("================================")

for medoid_pos, af_pos in MAPPING.items():

    key = (" ", af_pos, " ")

    if key not in chain:
        raise ValueError(
            f"Residue {CHAIN_ID}:{af_pos} not found."
        )

    residue = chain[key]
    selected[af_pos] = residue

    print(
        f"Medoid {medoid_pos:>3} "
        f"-> AF {af_pos:>3} "
        f"{residue.get_resname()}"
    )


# ------------------------------------------------------------
# pLDDT
# ------------------------------------------------------------

print("\nAlphaFold confidence")
print("====================")

for medoid_pos, af_pos in MAPPING.items():

    residue = selected[af_pos]

    values = [
        atom.get_bfactor()
        for atom in residue.get_atoms()
        if atom.element != "H"
    ]

    print(
        f"Medoid {medoid_pos:>3} "
        f"AF {af_pos:>3} "
        f"{residue.get_resname():>3} "
        f"mean pLDDT = {np.mean(values):6.2f}"
    )


# ------------------------------------------------------------
# SASA
# ------------------------------------------------------------

sr = ShrakeRupley()
sr.compute(model, level="R")

print("\nSolvent accessibility")
print("======================")

for medoid_pos, af_pos in MAPPING.items():

    residue = selected[af_pos]

    print(
        f"Medoid {medoid_pos:>3} "
        f"AF {af_pos:>3} "
        f"{residue.get_resname():>3} "
        f"SASA = {residue.sasa:8.2f} A^2"
    )


# ------------------------------------------------------------
# Distance helpers
# ------------------------------------------------------------

def heavy_atoms(residue):
    return [
        atom
        for atom in residue.get_atoms()
        if atom.element != "H"
    ]


def minimum_heavy_atom_distance(res1, res2):

    atoms1 = heavy_atoms(res1)
    atoms2 = heavy_atoms(res2)

    minimum = np.inf
    closest_pair = None

    for atom1 in atoms1:
        for atom2 in atoms2:

            d = atom1 - atom2

            if d < minimum:
                minimum = d
                closest_pair = (
                    atom1.get_name(),
                    atom2.get_name()
                )

    return minimum, closest_pair


# ------------------------------------------------------------
# Pairwise distances
# ------------------------------------------------------------

print("\nPairwise structural distances")
print("=============================")

for medoid1, medoid2 in combinations(MAPPING, 2):

    af1 = MAPPING[medoid1]
    af2 = MAPPING[medoid2]

    residue1 = selected[af1]
    residue2 = selected[af2]

    ca_distance = residue1["CA"] - residue2["CA"]

    heavy_distance, pair = (
        minimum_heavy_atom_distance(
            residue1,
            residue2
        )
    )

    print(
        f"Medoid {medoid1:>3} "
        f"(AF {af1:>3}) - "
        f"Medoid {medoid2:>3} "
        f"(AF {af2:>3}) : "
        f"C-alpha = {ca_distance:6.2f} A, "
        f"heavy = {heavy_distance:6.2f} A "
        f"({pair[0]} - {pair[1]})"
    )
    
    
    
    
    
WALKER_A = range(47, 55)

walker_residues = {
    pos: chain[(" ", pos, " ")]
    for pos in WALKER_A
}

print("\nDistance to Walker A motif")
print("==========================")

for medoid_pos, af_pos in MAPPING.items():

    focal = selected[af_pos]

    best_distance = np.inf
    best_residue = None
    best_atoms = None

    for motif_pos, motif_residue in walker_residues.items():

        distance, atom_pair = minimum_heavy_atom_distance(
            focal,
            motif_residue
        )

        if distance < best_distance:
            best_distance = distance
            best_residue = motif_pos
            best_atoms = atom_pair

    print(
        f"Medoid {medoid_pos:>3} "
        f"(AF {af_pos:>3}) -> "
        f"Walker A {best_residue}: "
        f"{best_distance:6.2f} A "
        f"({best_atoms[0]} - {best_atoms[1]})"
    )
    
    
    
    
    
coords = np.array([
    selected[37]["CA"].coord,
    selected[43]["CA"].coord,
    selected[56]["CA"].coord,
])

centroid = coords.mean(axis=0)

distances_to_centroid = np.linalg.norm(
    coords - centroid,
    axis=1
)

print("\nThree-site structural cluster")
print("=============================")

print(
    "Centroid:",
    " ".join(f"{x:.3f}" for x in centroid)
)

print(
    "Distances to centroid:",
    ", ".join(
        f"{d:.2f} A"
        for d in distances_to_centroid
    )
)

print(
    "Cluster diameter:",
    max(
        selected[a]["CA"] - selected[b]["CA"]
        for a, b in combinations([37, 43, 56], 2)
    )
)













walker_positions = list(range(47, 55))

walker_residues = [
    chain[(" ", pos, " ")]
    for pos in walker_positions
]

walker_distances = []

for residue in chain:

    if residue.get_resname() not in STANDARD_AA:
        continue

    position = residue.id[1]

    # Don't include Walker A itself
    if position in walker_positions:
        continue

    minimum = np.inf

    for walker_residue in walker_residues:

        distance, _ = minimum_heavy_atom_distance(
            residue,
            walker_residue
        )

        minimum = min(minimum, distance)

    walker_distances.append(
        (position, residue.get_resname(), minimum)
    )


walker_distances.sort(key=lambda x: x[2])

print("\nResidues closest to Walker A")
print("============================")

for position, aa, distance in walker_distances[:30]:

    marker = ""

    if position in MAPPING.values():
        marker = "  <-- focal"

    print(
        f"{position:>3} "
        f"{aa:>3} "
        f"{distance:6.2f} A"
        f"{marker}"
    )
    
    
    
all_distances = np.array(
    [x[2] for x in walker_distances]
)

print("\nFocal-site Walker-A proximity")
print("=============================")

for medoid_pos, af_pos in MAPPING.items():

    focal_distance = next(
        d
        for pos, aa, d in walker_distances
        if pos == af_pos
    )

    percentile = (
        np.mean(all_distances <= focal_distance) * 100
    )

    print(
        f"Medoid {medoid_pos:>3} "
        f"(AF {af_pos:>3}): "
        f"{focal_distance:6.2f} A, "
        f"{percentile:5.1f}th distance percentile"
    )