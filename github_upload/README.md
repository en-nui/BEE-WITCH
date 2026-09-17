# Ecological stability masks rapid, convergent evolution within the honey bee worker microbiome

**Authors:** Chris R. P. Robinson, Adam G. Dolezal, Irene L. G. Newton  
**Affiliations:** Department of Biology, Indiana University Bloomington & University of Illinois Urbana-Champaign  

---

NOTE: MUCH OF THIS IS HARDCODED AND MADE FOR CRPR. FEEL FREE TO REACH OUT WITH ANY QUESTIONS FOR YOUR PARTICULAR USE CASE. THIS CODE WILL LIKELY NOT WORK IF  YOU DOWNLOAD IT AND TRY TO RUN IT. 


## Overview

This repository contains the bioinformatic, population-genetic, and ecological analysis pipelines supporting the manuscript: **"Ecological stability masks rapid, convergent evolution within the honey bee worker microbiome"**.

The analysis couples deeply sequenced longitudinal metagenomics (298 metagenomic libraries across 12 colonies and 3 apiaries over 12 months) with population-genetic modeling to track bacterial lineage turnover, standing variation, purifying selection, rapid allele-frequency shifts ($\Delta\text{RAF}$), and convergent protein evolution within the honey bee gut microbiome.

---

## Directory Organization

```text
github_upload/
├── 01_upstream_pipeline/                  # Metagenome assembly, binning, dereplication, annotation & profiling
│   ├── 00_preprocessing/                 # Read QC, adapter trimming (bbduk), sample organization
│   ├── 01_assembly/                      # MEGAHIT co-assembly across biological replicates
│   ├── 02_mapping/                       # Bowtie2 read mapping and coverage estimation
│   ├── 05_binning/                       # MetaBAT2, VAMB, and Binette bin integration
│   ├── 07_taxonomy_gtdbtk/               # GTDB-Tk v2 taxonomic classification
│   ├── 08_dereplication_drep/            # dRep (95% ANI, 85% breadth) dereplication into rMAGs
│   ├── 10_annotation_dram/               # DRAM functional annotation
│   └── 15_anvio_profiling/               # Anvi'o contigs/profile database generation & single-base variability
│
├── 02_ecological_analyses/                # Ecological dynamics, alpha/beta diversity, turnover & synchrony (Figs 1-2)
│   ├── 01_rmag_abundance_and_persistence/# Coverage filtering, occupancy ($O_{ic}$), and persistence classification
│   ├── 02_alpha_beta_diversity/          # Shannon alpha diversity, Bray-Curtis ordination, dbRDA, PERMANOVA
│   ├── 03_turnover_synchrony_and_lorenz/ # Consecutive turnover ($D_{c,t}$), apiary synchrony, Lorenz curves ($N_{50}$)
│   └── 04_larval_analyses/               # Larval metagenome microbiome quantification and comparison
│
├── 03_population_genetics_and_diversity/  # SNV calling, nucleotide diversity, temporal differentiation (Fig 3)
│   ├── 01_snv_profiling_and_filtering/   # Core genome identification, copy-number filtering, SNV QC
│   ├── 02_nucleotide_diversity_pnps/     # Nucleotide diversity ($\pi, \pi_S, \pi_N$), Watterson's $\theta_W$, Tajima's $D$
│   └── 03_temporal_fst/                  # Genome-wide and genic temporal $F_{ST}$
│
├── 04_purifying_selection_models/         # Purifying selection null models and functional annotations (Fig 4)
│   ├── 01_piN_piS_decoupling_null_models/# Gene-splitting decoupling models ($\pi_N/\pi_S$ vs $\pi_S$)
│   └── 02_functional_pathway_annotations/# Pfam-to-GO mapping, KEGG pathway mapping, and genic selection
│
├── 05_allele_frequency_dynamics_sweeps/   # Rapid allele-frequency shifts and background tracking (Fig 5)
│   ├── 01_rapid_allele_frequency_shifts_deltaRAF/ # Binomial error models for $\Delta\text{RAF}$ detection & trajectory plots
│   ├── 02_resident_background_turnover/  # High-frequency marker alleles, resident background displacement ($B_t$)
│   └── 03_neutral_drift_wright_fisher/   # Wright-Fisher neutral diffusion drift model ($N_{e,\min}, Z_{\max}$)
│
└── 06_parallel_evolution_and_structure/   # Spatial linkage, recurrence across populations & structural mapping (Fig 6)
    ├── 01_spatial_temporal_covariance/   # Synonymous 4D SNV spatial (10–25 kb) & temporal covariance around 0D sweeps
    ├── 02_recurrent_protein_families_permutation/ # MMseqs2 protein clustering, concurrent $\Delta\text{RAF}$, permutation nulls
    └── 03_protein_structural_mapping/    # MAFFT alignment, AlphaFold2-ptm / ColabFold, heavy-atom distance & ATP ligand modeling
```

---

## Code Mapping by Manuscript Methods & Figures

### 1. Metagenomic Assembly, Binning, and Dereplication (Methods 5.2 – 5.3)
- **Quality Control & Read Cleaning:** `01_upstream_pipeline/00_preprocessing/` (`processing_data.py`, `rename.py`, `makedirectory_and_move.py`) using `bbduk`.
- **Co-Assembly:** `01_upstream_pipeline/01_assembly/` (`coassembly.py`, `assembly_report.py`, `larval_coassembly.py`) using `MEGAHIT` per colony timepoint.
- **Read Mapping:** `01_upstream_pipeline/02_mapping/` (`mapping.py`, `rMAG_mapping.py`, `cov_file_gen.py`) using `bowtie2` and `samtools`.
- **Binning & Bin Refinement:** `01_upstream_pipeline/05_binning/` (`metabat2_v2.py`, `vamb.py`, `vamb_write_fasta.py`, `binette.py`) integrating `MetaBAT2` and `VAMB` with `Binette`.
- **Taxonomic Assignment:** `01_upstream_pipeline/07_taxonomy_gtdbtk/` (`gtdbtk_script_v2.py`, `gtdbk_count_and_concat.py`) using `GTDB-Tk v2`.
- **Dereplication & rMAG Databases:** `01_upstream_pipeline/08_dereplication_drep/` (`mMAG_dRep_processing.py`, `mMAG_dRep_metadata_enrichment.py`, `MAG_trimming.py`, `ani_plots.R`, `fig_gen_drep_counts.R`) using `dRep` at 95% ANI and 85% breadth.
- **Functional Annotation & Profiling:** `01_upstream_pipeline/10_annotation_dram/` (`mMAG_DRAM_v2.py`) and `01_upstream_pipeline/15_anvio_profiling/` (`anvio_processing_parallel.py`, `anvio_gen_var_script.py`, `anvio_processing_v2.py`).

---

### 2. Ecological Analyses & Diversity Dynamics (Methods 5.4.1 – 5.4.3; Figures 1 & 2)
- **rMAG Quality Control & Persistence (Figure 3A):**  
  `02_ecological_analyses/01_rmag_abundance_and_persistence/`  
  - `mmag_relabund_script.py`, `mmag_presabs_script.py`, `mosdepth_deconvolution.py`: Depth filtering ($\text{median} \ge 5\times$), contig breadth ($\ge 25\%$), relative abundances.
  - `mMAG_relabund.R`, `mAMG_prevalence_plot.R`: Longitudinal occupancy ($O_{ic}$), persistent ($\ge 75\%$ eligible months) vs intermittent classification, and replicate-corrected temporal variance ($\sigma^2_{bio}$).
- **Alpha & Beta Diversity (Figures 1C–E, 2A–C):**  
  `02_ecological_analyses/02_alpha_beta_diversity/`  
  - `alpha_diversity.R`: Subsampled ($B=5000$) Shannon diversity trajectories across months and seasons.
  - `NMDS_Plots.R`, `NMDS_Plots_Genus.R`, `comeback_20250925_BETADIVERSITY_AND_GLMGAM.R`: Genus-level Bray-Curtis dissimilarity, condition-conditioned ordination (`Condition(colony)` dbRDA), PERMANOVA, time lag slopes, and circular shift seasonal permutations.
  - `We_actually_did_VEGAN_scripts.R`: Multivariate ecological modeling via `vegan`.
- **Community Turnover & Apiary Synchrony (Figures 2D–F):**  
  `02_ecological_analyses/03_turnover_synchrony_and_lorenz/`  
  - `synchrony_scripts.R`: Consecutive turnover ($D_{c,t}$), pairwise Spearman trajectories ($r_{ij}$), intra- vs inter-apiary synchrony bootstrap ($B=5000$).
  - `comeback_20250917.R`: Decomposition into individual MAG fractional contributions ($C_i$), Lorenz curves ($N_{50}, N_{80}$ turnover), and genus permutation tests with BH FDR correction.
- **Larval Microbiome Comparisons (Figure 1E):**  
  `02_ecological_analyses/04_larval_analyses/` (`larval_comparison_code.R`, `larval_contig_vs_coverage.R`).

---

### 3. Population Genetics & Standing Diversity (Methods 5.4.4; Figure 3)
- **Core Genome & SNV QC Filters:**  
  `03_population_genetics_and_diversity/01_snv_profiling_and_filtering/`  
  - `CORE_SNP_COVERAGE_FILTERS.py`, `PANGENOME_CORE_SNP_FILTER.py`, `PANGENOME_CNV_CORE_CALLER.py`: Copy-number ratio bounds ($0.3 \le C_{i,t} \le 3$, $<10\%$ failure rate) and coverage bounds ($0.3\bar{D} < D < 3\bar{D}$).
  - `POLARIZED_SNP_FILTERING.py`, `CLEAN_ANVIO_GENE_LIST.py`, `GENE_BLACKLIST_GENERATOR.py`, `GENE_MAP_GENERATOR.py`, `anvio_variability_parser.py`: Polarization, blacklist filtering ($\ge 95\%$ ANI across distinct rMAGs), and SNV extraction.
  - `high_snvs_low_snvs_PARSED_R_NO_VARIABLES.R`: Monthly SNV accumulation by population (Figure 3B).
- **Nucleotide Diversity & $p_N/p_S$ (Figure 3C):**  
  `03_population_genetics_and_diversity/02_nucleotide_diversity_pnps/`  
  - `LIVE_calculate_pnps_from_polarized_anvio_snps.py`: Computes polarized 0-fold (0D) and 4-fold (4D) degenerate nucleotide diversity ($\pi, \pi_S, \pi_N$), Watterson's $\theta_W$, and Tajima's $D$.
  - `nucDiv_by_gene.py`, `ADDING_METADATA_TO_NUCDIV_CSVs.py`, `nucdiv_plots.R`, `SummaryStatsAnalysis.R`, `core_heterozygosity_call.py`.
- **Temporal $F_{ST}$ (Figure 3D):**  
  `03_population_genetics_and_diversity/03_temporal_fst/`  
  - `Fst_single_colony_multiple_month.py`, `Fst_multiple_colony_multiple_month.py`, `GENIC-FST_single_colony_mutliple_month.py`, `temporal_fst.py`: Genome-wide and genic temporal Wright's $F_{ST}$.

---

### 4. Purifying Selection Models & Pathway Enrichment (Methods 5.4.6; Figure 4)
- **Gene-Splitting Decoupling Null Model (Figure 4A):**  
  `04_purifying_selection_models/01_piN_piS_decoupling_null_models/`  
  - `purifying_selection_and_model_code_R.R`, `modelCode.R`, `model_code_20251019.R`: Decoupled gene-splitting null model to assess $\pi_N/\pi_S$ as a function of genealogical depth ($\pi_S$) without shared-denominator artifacts.
- **Functional Annotations & Pathway Selection (Figures 4B & 4C):**  
  `04_purifying_selection_models/02_functional_pathway_annotations/`  
  - `pfam_to_go_mapping_genes.py`, `pfam_to_go_map.py`, `UpSET_PLOT_R.R`: Pfam-to-GO / KEGG metabolic pathway mappings, pathway-wide $\pi_N/\pi_S$ confidence intervals, and genic $F_{ST}$ vs $\pi$ deviations.

---

### 5. Allele-Frequency Dynamics & Selective Sweeps (Methods 5.4.5 – 5.4.6; Figure 5)
- **Rapid Allele Frequency Shifts ($\Delta\text{RAF}$) (Figures 5A, 5E):**  
  `05_allele_frequency_dynamics_sweeps/01_rapid_allele_frequency_shifts_deltaRAF/`  
  - `allele_frequency_R.R`, `allele_trajectories_and_enrichments.R`: Identifies transitions ($<20\%$ to $>70\%$) using binomial sampling error models ($\text{FDR} < 0.1, q < 0.05$) and validates against biological replicates.
  - `1_compute_gene_sweeps_from_snva_data.py`, `combine_signatures.py`, `EXTREME_SNP_TRAJECTORY.py`, `FREQUENCY_TRAJECTORY_SNV_SUBSET.py`.
- **Genomic Background Displacement ($B_t$) (Figures 5B–D):**  
  `05_allele_frequency_dynamics_sweeps/02_resident_background_turnover/`  
  - `private_marker_v2.R`, `COLONY_WIDE_PRIVATE_ALLELES.py`, `PRIVATE_ALLELES_TEMPORAL.py`, `PRIVATE_ALLELES_TEMPORAL_V2.py`: Tracks baseline marker alleles ($f_{\text{resident}} > 0.9$ at $t_0$), calculates displacement score $B_t = 1 - \text{median}(f_{\text{resident},t})$, and classifies resident-preserving sweeps ($B_t \le 0.20, |\Delta B| \le 0.10$).
- **Neutral Drift Wright-Fisher Model (Figure 5C):**  
  `05_allele_frequency_dynamics_sweeps/03_neutral_drift_wright_fisher/`  
  - `sweep_drift_null_model.R`: Estimates minimum effective population size ($N_{e,\min}$) from 4D variance under Wright-Fisher diffusion and generates standardized trajectory scores $Z = \Delta f / \sqrt{V_{tot}}$ and $Z_{\max}$.

---

### 6. Convergent Evolution, Linkage, and Structural Mapping (Methods 5.4.6; Figure 6)
- **Spatial & Temporal Covariance of Synonymous Sites (Figures 6A & 6B):**  
  `06_parallel_evolution_and_structure/01_spatial_temporal_covariance/`  
  - `sweep_piS_delta.R`: Measures excess change in 4D synonymous SNVs as a function of physical distance ($\le 5\text{ kb}$ to $25\text{ kb}$) from focal 0D $\Delta\text{RAF}$ loci, and compares focal ($t=0$) intervals against region-specific permutation nulls.
- **Recurrent & Concurrent Protein Families (Figures 6C & 6D):**  
  `06_parallel_evolution_and_structure/02_recurrent_protein_families_permutation/`  
  - `alphafold_sweep_pipeline.R`, `parallel_tests_20251218.R`, `complex_parallelism_goldsilver_analysis.R`, `matrix_and_upset.R`, `parallelism_sweep_pnps_code.R`: MMseqs2 protein clustering ($\ge 50\%$ identity, $\ge 85\%$ coverage), identification of temporally concurrent 0D $\Delta\text{RAF}$ events across $\ge 2$ MAGs, genera, and colonies, tested against gene-length-weighted permutation nulls ($B = 5000$).
- **Structural Mapping & AlphaFold Modeling (Figure 6E):**  
  `06_parallel_evolution_and_structure/03_protein_structural_mapping/`  
  - `sweep_mafft.py`: MAFFT alignment of candidate protein homologs and mapping of population-specific variants onto representative medoid sequences.
  - `protein_mutation_to_alignment.py`, `gene_marker_distances.py`, `protein_prediction_resilience.py`: Structure coordinate extraction from AlphaFold2-ptm / ColabFold v1.6.2 models, calculation of pairwise $\text{C}\alpha$ and heavy-atom distances (e.g. 13.24 Å clustering around conserved ATP-binding loop), and modeling distances to the Walker A motif and ATP ligand (from aligned structure PDB 1L2T via UCSF ChimeraX).

---

## Dependencies & Environment

### Bioinformatic Tools
- **Assembly & Mapping:** `bbmap/bbduk` (v38+), `MEGAHIT` (v1.2.9), `Bowtie2` (v2.5+), `samtools` (v1.17+)
- **Binning & Quality:** `MetaBAT2` (v2.15), `VAMB` (v3.0+), `Binette` (v1.0+), `CheckM2` (v1.0.2)
- **Taxonomy & Dereplication:** `GTDB-Tk` (v2.3+ / Release 214+), `dRep` (v3.4.0+)
- **Annotation & Pangenomics:** `DRAM` (v1.4+), `MMseqs2` (v14+), `anvi'o` (v7.1+)
- **Phylogenetics & Structure:** `MAFFT` (v7.5+), `ColabFold` / `AlphaFold2-ptm` (v1.6.2), `UCSF ChimeraX` (v1.6+)

### R Environment & Key Packages
- R $\ge$ 4.2.0
- `vegan`, `tidyverse` (`ggplot2`, `dplyr`, `tidyr`, `readr`, `purrr`), `data.table`, `lme4`, `UpSetR`, `ape`, `Biostrings`, `cowplot`, `ggrepel`, `pheatmap`, `FSA`

### Python Environment & Key Packages
- Python $\ge$ 3.9
- `numpy`, `scipy`, `pandas`, `scikit-learn`, `biopython`, `tqdm`, `pyarrow`, `fastparquet`

---

## Citation & Contact

If you use code or workflows from this repository, please cite:
> Robinson, C. R. P., Dolezal, A. G., & Newton, I. L. G. (2026). *Ecological stability masks rapid, convergent evolution within the honey bee worker microbiome.*

For questions regarding data or implementation, please contact Chris R. P. Robinson (Indiana University Bloomington).
