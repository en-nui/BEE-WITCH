library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(purrr)
library(broom)
library(Hmisc)
library(rcorr)
r#!/usr/bin/env Rscript

# ===============================================================
#  Sweep localization analysis
#  MAG_96, Colony_104
# ===============================================================

# ----------------- Paths & Parameters -----------------
MAG_ID <- "mmag_191"
COLONY_ID <- "106"

base_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity"
sweep_path <- file.path(base_dir, "sweep_analysis/significant_sweep_snvs_ALL.csv")
traj_path  <- file.path(base_dir, "snv_frequency_trajectories_polymorphic_only",
                        paste0("snv_frequency_trajectory_", COLONY_ID, ".csv"))
out_dir    <- file.path(base_dir, "sweep_analysis/figures_sweep_localization", paste0(MAG_ID, "_", COLONY_ID))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# -------------------------------------------------------

# ----------------- Load & Filter Data -----------------
# 1. Sweeping SNVs
sweeps <- read_csv(sweep_path, show_col_types = FALSE) %>%
  filter(MAG_id == MAG_ID, colony_id == as.numeric(COLONY_ID)) %>%
  mutate(unique_SNV_identifier = paste(MAG_id, contig, position, sep = "_"))

cat(nrow(sweeps), "sweeping SNVs loaded for", MAG_ID, "in colony", COLONY_ID, "\n")

# 2. Trajectory file (to get gene info)
traj <- read_csv(traj_path, show_col_types = FALSE) %>%
  filter(MAG_id == MAG_ID) %>%
  select(MAG_id, contig_name, pos_in_contig,
         unique_gene_callers_id, corresponding_gene_call) %>%
  distinct(MAG_id, contig_name, pos_in_contig, .keep_all = TRUE)

cat(nrow(traj), "unique SNV positions in trajectory file for", MAG_ID, "\n")




# Deduplicate sweeps so each unique SNV is represented once
sweeps_unique <- sweeps %>%
  distinct(MAG_id, contig, position, .keep_all = TRUE)

cat(nrow(sweeps_unique), "unique sweeping SNV positions in sweep file for", MAG_ID, "\n")

# Clean, one-to-one join between unique sweeps and unique trajectory positions
annotated <- sweeps_unique %>%
  left_join(traj,
            by = c("MAG_id" = "MAG_id",
                   "contig" = "contig_name",
                   "position" = "pos_in_contig"),
            relationship = "one-to-one")



# Save annotated table for record
write_csv(annotated, file.path(out_dir, "sweep_snvs_with_gene_info.csv"))









gene_summary <- annotated %>%
  drop_na(unique_gene_callers_id) %>%
  group_by(unique_gene_callers_id, corresponding_gene_call) %>%
  summarise(n_sweeps = n(), .groups = "drop")

ggplot(gene_summary, aes(x = n_sweeps)) +
  geom_histogram(binwidth = 1, fill = "#1f78b4", color = "white", boundary = 0.5) +
  scale_x_continuous(breaks = seq(1, max(gene_summary$n_sweeps), 1)) +
  labs(x = "Sweep SNVs per gene", y = "Number of genes",
       title = paste(MAG_ID, "Colony", COLONY_ID)) +
  theme_minimal(base_size = 12)


total_genes <- traj %>%
  distinct(unique_gene_callers_id) %>%
  nrow()

genes_with_sweeps <- nrow(gene_summary)
prop_swept <- genes_with_sweeps / total_genes

prop_df <- tibble(
  category = c("Swept genes", "Other genes"),
  count = c(genes_with_sweeps, total_genes - genes_with_sweeps)
)

ggplot(prop_df, aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = c("#e31a1c", "gray80")) +
  labs(title = sprintf("%s (Colony %s): %.1f%% genes swept",
                       MAG_ID, COLONY_ID, 100 * prop_swept)) +
  theme_void(base_size = 12) +
  theme(legend.position = "bottom")


# approximate Poisson enrichment
lambda <- nrow(annotated) / total_genes
gene_summary <- gene_summary %>%
  mutate(p_poisson = ppois(n_sweeps - 1, lambda = lambda, lower.tail = FALSE),
         p_adj = p.adjust(p_poisson, method = "BH")) %>%
  arrange(p_adj)
head(gene_summary, 10)







# ===============================================================
# 1️⃣ Sweep SNVs per gene per time window
# ===============================================================
gene_time <- annotated %>%
  drop_na(unique_gene_callers_id) %>%
  mutate(window = paste0(t1, "_to_", t2)) %>%
  group_by(window, unique_gene_callers_id, corresponding_gene_call) %>%
  summarise(n_sweeps = n(), .groups = "drop")

write_csv(gene_time, file.path(out_dir, "gene_sweep_counts_by_window.csv"))
cat("Saved: gene_sweep_counts_by_window.csv\n")

# ===============================================================
# 2️⃣ Top genes per time window
# ===============================================================
top_genes <- gene_time %>%
  group_by(window) %>%
  arrange(desc(n_sweeps), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write_csv(top_genes, file.path(out_dir, "top_genes_per_window.csv"))
cat("Saved: top_genes_per_window.csv\n")

# ===============================================================
# 3️⃣ Proportion of genes affected per time window
# ===============================================================
total_genes <- traj %>%
  distinct(unique_gene_callers_id) %>%
  nrow()

gene_prop <- gene_time %>%
  group_by(window) %>%
  summarise(
    genes_with_sweeps = n(),
    total_genes = total_genes,
    prop_swept = genes_with_sweeps / total_genes
  )

write_csv(gene_prop, file.path(out_dir, "prop_genes_with_sweeps_by_window.csv"))
cat("Saved: prop_genes_with_sweeps_by_window.csv\n")

# ===============================================================
# 4️⃣ Quick visualization (optional)
# ===============================================================

# Histogram: sweeps per gene (across all time windows)
p1 <- ggplot(gene_time, aes(x = n_sweeps)) +
  geom_histogram(binwidth = 1, fill = "#1f78b4", color = "white", boundary = 0.5) +
  facet_wrap(~ window, ncol = 2) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Sweep SNVs per gene",
    y = "Number of genes",
    title = paste(MAG_ID, "Colony", COLONY_ID, "\nSweep SNVs per gene across time windows")
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "fig1_histogram_sweeps_per_gene_by_window.png"),
       p1, width = 6, height = 4, dpi = 300)

# Bar: proportion of genes affected per window
p2 <- ggplot(gene_prop, aes(x = window, y = prop_swept)) +
  geom_col(fill = "#e31a1c") +
  geom_text(aes(label = scales::percent(prop_swept, accuracy = 0.1)),
            vjust = -0.5, size = 3) +
  labs(
    x = "Time window (t1→t2)",
    y = "Proportion of genes with ≥1 sweep SNV",
    title = paste(MAG_ID, "Colony", COLONY_ID)
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "fig2_prop_genes_with_sweeps_by_window.png"),
       p2, width = 5, height = 3.5, dpi = 300)

cat("✅ Figures and summary tables saved in:\n", out_dir, "\n")
p1
p2



### BIG LOOP ####
base_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity"
sweep_path <- file.path(base_dir, "sweep_analysis/significant_sweep_snvs_ALL.csv")
traj_dir <- file.path(base_dir, "snv_frequency_trajectories_polymorphic_only")
out_root <- file.path(base_dir, "sweep_analysis/all_MAGs_summaries")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# month ordering - edit if you have other names
month_order <- c("May","June","July","August","September","October","November","January","February")

sweeps_all <- read_csv(sweep_path, show_col_types = FALSE)

MAG_list <- unique(sweeps_all$MAG_id)

for (MAG_ID in MAG_list) {
  # identify colonies for which sweeps exist for this MAG
  sweeps_mag <- sweeps_all %>% filter(MAG_id == MAG_ID)
  colonies <- unique(sweeps_mag$colony_id)
  
  # we will process each colony separately (trajectory files are per colony)
  for (COLONY_ID in colonies) {
    out_dir <- file.path(out_root, paste0(MAG_ID, "_", COLONY_ID))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # load sweeps and trajectory; skip if traj missing
    sweeps <- sweeps_mag %>% filter(colony_id == COLONY_ID) %>%
      mutate(window = paste0(t1, "_to_", t2)) %>%
      distinct(MAG_id, contig, position, t1, t2, .keep_all = TRUE)
    traj_file <- file.path(traj_dir, paste0("snv_frequency_trajectory_", COLONY_ID, ".csv"))
    if (!file.exists(traj_file)) {
      message("Missing traj for colony ", COLONY_ID, " skipping ", MAG_ID)
      next
    }
    traj <- read_csv(traj_file, show_col_types = FALSE) %>%
      filter(MAG_id == MAG_ID) %>%
      distinct(MAG_id, contig_name, pos_in_contig, unique_gene_callers_id, corresponding_gene_call)
    
    # annotate sweeps
    annotated <- sweeps %>%
      left_join(traj, by = c("MAG_id" = "MAG_id",
                             "contig" = "contig_name",
                             "position" = "pos_in_contig"))
    
    # gene counts by window
    gene_time <- annotated %>%
      drop_na(unique_gene_callers_id) %>%
      mutate(window = paste0(t1, "_to_", t2)) %>%
      group_by(window, unique_gene_callers_id, corresponding_gene_call) %>%
      summarise(n_sweeps = n(), .groups = "drop")
    
    write_csv(gene_time, file.path(out_dir, paste0("gene_sweep_counts_by_window_", MAG_ID, "_", COLONY_ID, ".csv")))
    
    # FIGURE 1: histogram faceted by window, order months
    p1 <- ggplot(gene_time, aes(x = n_sweeps)) +
      geom_histogram(binwidth = 1, fill = "#1f78b4", color = "white", boundary = 0.5) +
      facet_wrap(~ window, scales = "free_y", ncol = 2) +
      labs(title = paste(MAG_ID, "colony", COLONY_ID),
           x = "Sweep SNVs per gene", y = "Number of genes") +
      theme_minimal(base_size = 10)
    
    ggsave(file.path(out_dir, paste0("hist_sweeps_per_gene_", MAG_ID, "_", COLONY_ID, ".png")),
           p1, width = 6, height = 4, dpi = 300)
    
    # FIGURE 2: proportion of genes affected per window
    total_genes <- traj %>% distinct(unique_gene_callers_id) %>% nrow()
    gene_prop <- gene_time %>%
      group_by(window) %>%
      summarise(genes_with_sweeps = n(), .groups = "drop") %>%
      mutate(total_genes = total_genes, prop_swept = genes_with_sweeps / total_genes) %>%
      mutate(window = factor(window, levels = unique(window)))  # preserve observed order
    
    p2 <- ggplot(gene_prop, aes(x = window, y = prop_swept)) +
      geom_col(fill = "#e31a1c") +
      geom_text(aes(label = scales::percent(prop_swept, accuracy = 0.1)), vjust = -0.5, size = 3) +
      labs(x = "Time window (t1→t2)", y = "Proportion of genes with ≥1 sweep SNV") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(out_dir, paste0("prop_genes_with_sweeps_by_window_", MAG_ID, "_", COLONY_ID, ".png")),
           p2, width = 5, height = 3.5, dpi = 300)
    
    message("Saved plots for ", MAG_ID, " colony ", COLONY_ID, " to ", out_dir)
  }
} # end loop MAGs






# -------------------- Parameters ----------------
nperm <- 2000     # permutations per window
min_snvs <- 100   # minimum sweeps per MAG×colony to include
alpha_sig <- 0.05 # significance threshold
# ------------------------------------------------

# -------------------- Load sweep data ----------------
sweeps_all <- read_csv(sweep_path, show_col_types = FALSE) %>%
  mutate(window = paste0(t1, "_to_", t2))

cat("Total sweep SNVs:", nrow(sweeps_all), "\n")

mag_summary <- sweeps_all %>%
  group_by(MAG_id, colony_id) %>%
  summarise(n_snvs = n(), .groups = "drop") %>%
  filter(n_snvs >= min_snvs)

cat("Eligible MAG×colony pairs:", nrow(mag_summary), "\n")

# ===============================================================
#  Main loop
# ===============================================================

all_hits <- list()  # collect all significant genes across MAGs

for (i in seq_len(nrow(mag_summary))) {
  
  MAG_ID    <- mag_summary$MAG_id[i]
  COLONY_ID <- mag_summary$colony_id[i]
  out_dir   <- file.path(out_root, paste0(MAG_ID, "_", COLONY_ID))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\n▶ Running weighted permutation for ", MAG_ID, " (Colony ", COLONY_ID, ")")
  
  # ---------------- Load per-MAG data ----------------
  sweeps <- sweeps_all %>%
    filter(MAG_id == MAG_ID, colony_id == COLONY_ID) %>%
    distinct(MAG_id, contig, position, t1, t2, .keep_all = TRUE)
  
  traj_file <- file.path(traj_dir, paste0("snv_frequency_trajectory_", COLONY_ID, ".csv"))
  if (!file.exists(traj_file)) {
    message("⚠️ Missing trajectory file for colony ", COLONY_ID, " — skipping.")
    next
  }
  
  traj <- read_csv(traj_file, show_col_types = FALSE) %>%
    filter(MAG_id == MAG_ID) %>%
    select(MAG_id, contig_name, pos_in_contig,
           unique_gene_callers_id, corresponding_gene_call, gene_length) %>%
    mutate(unique_gene_callers_id = as.character(unique_gene_callers_id)) %>%
    distinct(MAG_id, contig_name, pos_in_contig, .keep_all = TRUE)
  
  # Weight genes by gene_length
  gene_table <- traj %>%
    distinct(unique_gene_callers_id, corresponding_gene_call, gene_length) %>%
    mutate(
      gene_length = ifelse(is.na(gene_length) | gene_length <= 0, median(gene_length, na.rm = TRUE), gene_length),
      weight = gene_length / sum(gene_length, na.rm = TRUE)
    )
  
  genes <- gene_table$unique_gene_callers_id
  weights <- gene_table$weight
  G <- length(genes)
  if (G == 0) {
    message("⚠️ No genes found for ", MAG_ID, " ", COLONY_ID)
    next
  }
  
  # ------------------------------------------------------
  #  Permutation per time window
  # ------------------------------------------------------
  windows <- unique(sweeps$t1 %>% paste0("_to_", sweeps$t2))
  perm_results <- list()
  
  for (w in windows) {
    sw <- sweeps %>% mutate(window = paste0(t1, "_to_", t2)) %>% filter(window == w)
    m <- nrow(sw)
    if (m == 0) next
    
    message("  • Window ", w, " (", m, " sweeps)")
    
    # Observed counts per gene
    obs <- sw %>%
      left_join(traj, by = c("MAG_id" = "MAG_id",
                             "contig" = "contig_name",
                             "position" = "pos_in_contig")) %>%
      mutate(unique_gene_callers_id = as.character(unique_gene_callers_id)) %>%
      group_by(unique_gene_callers_id) %>%
      summarise(obs_n = n(), .groups = "drop") %>%
      right_join(tibble(unique_gene_callers_id = genes), by = "unique_gene_callers_id") %>%
      mutate(obs_n = replace_na(obs_n, 0)) %>%
      arrange(unique_gene_callers_id)
    
    # Permutation: weighted sampling by gene length
    perm_counts <- matrix(0L, nrow = G, ncol = nperm)
    for (p in seq_len(nperm)) {
      draws <- sample(genes, size = m, replace = TRUE, prob = weights)
      tcount <- table(draws)
      idx <- match(names(tcount), genes)
      perm_counts[idx, p] <- as.integer(tcount)
    }
    
    p_emp <- sapply(seq_len(G), function(j) mean(perm_counts[j, ] >= obs$obs_n[j]))
    perm_df <- tibble(
      MAG_id = MAG_ID,
      colony_id = COLONY_ID,
      window = w,
      unique_gene_callers_id = genes,
      obs_n = obs$obs_n,
      p_emp = p_emp
    ) %>%
      mutate(
        p_adj = p.adjust(p_emp, method = "BH"),
        is_sig = p_adj < alpha_sig
      )
    
    perm_results[[w]] <- perm_df
  }
  
  perm_all <- bind_rows(perm_results)
  if (nrow(perm_all) == 0) {
    message("⚠️ No results for ", MAG_ID, " ", COLONY_ID)
    next
  }
  
  # Save per-MAG results
  write_csv(perm_all, file.path(out_dir, paste0("perm_enrichment_weighted_", MAG_ID, "_", COLONY_ID, ".csv")))
  
  # Append significant hits to master list
  hits <- perm_all %>%
    filter(is_sig) %>%
    left_join(gene_table, by = "unique_gene_callers_id") %>%
    select(MAG_id, colony_id, window, unique_gene_callers_id, corresponding_gene_call,
           obs_n, p_emp, p_adj, gene_length)
  if (nrow(hits) > 0) all_hits[[paste(MAG_ID, COLONY_ID)]] <- hits
  
  # Temporal recurrence summary (how many windows a gene is enriched)
  temp_summary <- perm_all %>%
    filter(is_sig) %>%
    group_by(unique_gene_callers_id) %>%
    summarise(n_windows_sig = n_distinct(window),
              windows = paste(window, collapse = "; "),
              .groups = "drop") %>%
    arrange(desc(n_windows_sig))
  if (nrow(temp_summary) > 0) {
    write_csv(temp_summary, file.path(out_dir, "temporal_recurrence_summary.csv"))
  }
  
  message("✅ Completed ", MAG_ID, " (Colony ", COLONY_ID, ") — ",
          sum(perm_all$is_sig), " significant genes")
}

# ===============================================================
# Combine all hits into one master table
# ===============================================================
if (length(all_hits) > 0) {
  all_hits_df <- bind_rows(all_hits, .id = "MAG_Colony")
  write_csv(all_hits_df, file.path(out_root, "significant_genes_allMAGs.csv"))
  
  cat("\n✅ Combined significant hits written to:\n",
      file.path(out_root, "significant_genes_allMAGs.csv"), "\n")
  
  cat("Genes appearing in >1 time window (temporal recurrence):\n")
  print(all_hits_df %>%
          group_by(MAG_id, unique_gene_callers_id) %>%
          summarise(n_windows = n_distinct(window), .groups = "drop") %>%
          filter(n_windows > 1) %>%
          arrange(desc(n_windows)))
} else {
  cat("\n⚠️ No significant hits detected.\n")
}





