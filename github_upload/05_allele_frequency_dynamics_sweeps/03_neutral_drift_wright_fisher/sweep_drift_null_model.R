#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)
## ---------- GLOBAL PARAMETERS ----------

gens_per_day   <- 12        # ~2h doubling
coverage_min   <- 10L
min_intervals  <- 3L       # min neutral intervals per MAG for Ne

month_base_times <- c(
  "May"       = 0,
  "June"      = 36,
  "July"      = 31,
  "August"    = 34,
  "September" = 28,
  "October"   = 22,
  "November"  = 42,
  "January"   = 47,
  "February"  = 31
)
cumulative_times <- cumsum(month_base_times)
month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))

traj_dir   <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/snv_frequency_trajectories_polymorphic_only"
deg_dir    <- "/home/robinch/projects/BEE-WITCH/40_internal_pnps_calculations_from_anvio/colony_outputs"
sweep_file <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/significant_sweep_snvs_ALL.csv"

## ---------- HELPERS ----------

add_time_from_month <- function(dt) {
  dt[, Time := month_to_cumulative_time_map[Month]]
  dt
}

# numerically robust drift test for a single SNV
test_snv_drift <- function(dt_snv, Ne_hat, gens_per_day, coverage_min = 20L) {
  dt_snv <- copy(dt_snv)
  dt_snv <- add_time_from_month(dt_snv)
  dt_snv[, GenTime := Time * gens_per_day]
  setorder(dt_snv, GenTime)
  
  dt_snv[, f := departure_from_polarized_consensus]
  dt_snv[, `:=`(
    f_next       = shift(f,        type = "lead"),
    cov_next     = shift(coverage, type = "lead"),
    GenTime_next = shift(GenTime,  type = "lead")
  )]
  
  dt_int <- dt_snv[
    !is.na(f_next) &
      coverage >= coverage_min &
      cov_next >= coverage_min
  ]
  
  if (nrow(dt_int) < 1L) {
    return(list(
      T_stat      = NA_real_,
      df          = 0L,
      p           = NA_real_,
      Z_max       = NA_real_,
      n_intervals = 0L
    ))
  }
  
  # clamp away from 0 and 1 to avoid 0 variance
  eps <- 1e-6
  dt_int[, f_clamp      := pmin(pmax(f,      eps), 1 - eps)]
  dt_int[, f_next_clamp := pmin(pmax(f_next, eps), 1 - eps)]
  
  dt_int[, tau := GenTime_next - GenTime]
  dt_int[, df  := f_next - f]
  
  dt_int[, V_sample := f_clamp * (1 - f_clamp) * (1/coverage + 1/cov_next)]
  dt_int[, V_drift  := f_clamp * (1 - f_clamp) * (1 - exp(-tau / (2 * Ne_hat)))]
  dt_int[, V_tot    := V_sample + V_drift]
  
  dt_int <- dt_int[is.finite(V_tot) & V_tot > 0]
  if (nrow(dt_int) < 1L) {
    return(list(
      T_stat      = NA_real_,
      df          = 0L,
      p           = NA_real_,
      Z_max       = NA_real_,
      n_intervals = 0L
    ))
  }
  
  dt_int[, Z := df / sqrt(V_tot)]
  dt_int <- dt_int[is.finite(Z)]
  if (nrow(dt_int) < 1L) {
    return(list(
      T_stat      = NA_real_,
      df          = 0L,
      p           = NA_real_,
      Z_max       = NA_real_,
      n_intervals = 0L
    ))
  }
  
  T_stat <- sum(dt_int$Z^2)
  df_chi <- nrow(dt_int)
  p_val  <- 1 - pchisq(T_stat, df = df_chi)
  
  list(
    T_stat      = T_stat,
    df          = df_chi,
    p           = p_val,
    Z_max       = max(abs(dt_int$Z)),
    n_intervals = df_chi
  )
}

## ---------- PER-COLONY FUNCTION ----------

run_colony_all_modes <- function(colony_id) {
  cat("\n======= Colony", colony_id, "=======\n")
  
  traj_file <- file.path(traj_dir, paste0("snv_frequency_trajectory_", colony_id, ".csv"))
  deg_file  <- file.path(deg_dir,  paste0("colony_", colony_id, "_snps_with_degeneracy.csv"))
  
  if (!file.exists(traj_file)) {
    cat("  Trajectory file missing, skipping.\n")
    return(NULL)
  }
  if (!file.exists(deg_file)) {
    cat("  Degeneracy file missing, skipping.\n")
    return(NULL)
  }
  
  ## --- Load trajectories ---
  traj <- fread(traj_file, showProgress = TRUE)
  traj <- traj[, .(MAG_id, Month, contig_name, pos_in_contig,
                   coverage, departure_from_polarized_consensus)]
  traj <- add_time_from_month(traj)
  traj[, key := paste(MAG_id, contig_name, pos_in_contig)]
  
  ## --- Load degeneracy ---
  deg <- fread(deg_file, showProgress = TRUE)
  deg <- deg[Colony == colony_id]
  deg_pos <- unique(
    deg[, .(MAG_id, contig_name, pos_in_contig, degeneracy)],
    by = c("MAG_id", "contig_name", "pos_in_contig")
  )
  
  ## --- Join degeneracy onto traj ---
  traj_deg <- merge(
    traj,
    deg_pos,
    by = c("MAG_id", "contig_name", "pos_in_contig"),
    all.x = TRUE
  )
  traj_deg[, key := paste(MAG_id, contig_name, pos_in_contig)]
  
  ## --- Load sweeps for this colony ---
  sweeps_all <- fread(sweep_file, showProgress = FALSE)
  sweeps_c   <- sweeps_all[colony_id == colony_id]
  setnames(sweeps_c,
           old = c("contig", "position"),
           new = c("contig_name", "pos_in_contig"))
  sweeps_c[, key := paste(MAG_id, contig_name, pos_in_contig)]
  
  ## --- Local helper: compute Ne for a neutral_mode ---
  compute_Ne_for_mode <- function(neutral_mode) {
    cat("  [", neutral_mode, "] computing Ne...\n", sep = "")
    
    neutral_traj <- traj_deg[!key %in% sweeps_c$key]
    
    if (neutral_mode == "4D") {
      neutral_traj_use <- neutral_traj[degeneracy == "4D"]
    } else if (neutral_mode == "non0D") {
      neutral_traj_use <- neutral_traj[is.na(degeneracy) | degeneracy != "0D"]
    } else if (neutral_mode == "all") {
      neutral_traj_use <- neutral_traj
    } else {
      stop("Unknown neutral_mode: ", neutral_mode)
    }
    
    cat("    Neutral rows: ", nrow(neutral_traj_use), "\n", sep = "")
    
    if (nrow(neutral_traj_use) == 0L) {
      return(data.table(MAG_id = character(),
                        Ne = numeric(),
                        n_intervals = integer(),
                        S1 = numeric(), S2 = numeric(), S3 = numeric()))
    }
    
    neutral_traj_use[, GenTime := Time * gens_per_day]
    neutral_traj_use[, f := departure_from_polarized_consensus]
    setorder(neutral_traj_use, MAG_id, contig_name, pos_in_contig, GenTime)
    
    neutral_traj_use[, `:=`(
      GenTime_next = shift(GenTime, type = "lead"),
      f_next       = shift(f,       type = "lead"),
      cov_next     = shift(coverage, type = "lead")
    ), by = .(MAG_id, contig_name, pos_in_contig)]
    
    intervals <- neutral_traj_use[
      !is.na(f_next) &
        coverage >= coverage_min &
        cov_next >= coverage_min
    ]
    cat("    Neutral intervals: ", nrow(intervals), "\n", sep = "")
    
    if (nrow(intervals) == 0L) {
      return(data.table(MAG_id = character(),
                        Ne = numeric(),
                        n_intervals = integer(),
                        S1 = numeric(), S2 = numeric(), S3 = numeric()))
    }
    
    intervals[, tau := GenTime_next - GenTime]
    intervals[, df  := f_next - f]
    intervals[, V_sample := f * (1 - f) * (1/coverage + 1/cov_next)]
    
    Ne_by_MAG <- intervals[, {
      S1 <- sum(f * (1 - f) * tau)
      S2 <- sum(V_sample)
      S3 <- sum(df^2)
      n_int <- .N
      
      Ne_hat <- if (S3 > S2 && n_int >= min_intervals) {
        S1 / (2 * (S3 - S2))
      } else {
        NA_real_
      }
      
      .(Ne = Ne_hat,
        n_intervals = n_int,
        S1 = S1, S2 = S2, S3 = S3)
    }, by = MAG_id]
    
    cat("    MAGs with finite Ne: ", sum(!is.na(Ne_by_MAG$Ne)), "\n", sep = "")
    Ne_by_MAG
  }
  
  ## --- Compute Ne for all 3 neutral modes ---
  Ne_4D    <- compute_Ne_for_mode("4D")
  Ne_non0D <- compute_Ne_for_mode("non0D")
  Ne_all   <- compute_Ne_for_mode("all")
  
  # Merge Ne tables
  Ne_4D_m <- Ne_4D[,    .(MAG_id, Ne_4D    = Ne, n_int_4D    = n_intervals)]
  Ne_n0_m <- Ne_non0D[, .(MAG_id, Ne_non0D = Ne, n_int_non0D = n_intervals)]
  Ne_all_m<- Ne_all[,   .(MAG_id, Ne_all   = Ne, n_int_all   = n_intervals)]
  
  Ne_all_modes <- Reduce(
    function(x, y) merge(x, y, by = "MAG_id", all = TRUE),
    list(Ne_4D_m, Ne_n0_m, Ne_all_m)
  )
  
  # Ne_min per MAG (minimum of the three, ignoring NA)
  if (nrow(Ne_all_modes) > 0L) {
    Ne_all_modes[, Ne_min := {
      xs <- unlist(.SD)
      xs <- xs[!is.na(xs)]
      if (length(xs) == 0) NA_real_ else min(xs)
    }, by = MAG_id,
    .SDcols = c("Ne_4D", "Ne_non0D", "Ne_all")]
  }
  
  cat("  MAGs with finite Ne_min: ",
      sum(!is.na(Ne_all_modes$Ne_min)), "\n", sep = "")
  
  ## --- Sweep tests using Ne_min ---
  dt_c <- traj_deg
  
  sweeps_c_unique <- unique(
    sweeps_c,
    by = c("MAG_id", "contig_name", "pos_in_contig")
  )
  sweeps_c_unique <- sweeps_c_unique[key %in% dt_c$key]
  cat("  Sweeps with trajectories: ", nrow(sweeps_c_unique), "\n", sep = "")
  
  # Attach Ne_min
  sweeps_c_unique <- merge(
    sweeps_c_unique,
    Ne_all_modes[, .(MAG_id, Ne_min)],
    by = "MAG_id",
    all.x = TRUE
  )
  
  sweep_tests <- sweeps_c_unique[, {
    snv_traj <- dt_c[key == .BY$key]
    
    res <- if (!is.na(Ne_min) && Ne_min > 0) {
      test_snv_drift(
        dt_snv       = snv_traj,
        Ne_hat       = Ne_min,
        gens_per_day = gens_per_day,
        coverage_min = coverage_min
      )
    } else {
      list(
        T_stat      = NA_real_,
        df          = 0L,
        p           = NA_real_,
        Z_max       = NA_real_,
        n_intervals = 0L
      )
    }
    
    .(T_stat      = res$T_stat,
      df          = res$df,
      p           = res$p,
      Z_max       = res$Z_max,
      n_intervals = res$n_intervals)
  }, by = .(MAG_id, contig_name, pos_in_contig, key)]
  
  sweeps_c_tests <- merge(
    sweeps_c_unique,
    sweep_tests,
    by = c("MAG_id", "contig_name", "pos_in_contig", "key"),
    all.x = TRUE
  )
  
  ## --- Write outputs for this colony ---
  Ne_out_file     <- paste0("colony_", colony_id, "_Ne_by_MAG_allmodes.tsv")
  sweeps_out_file <- paste0("colony_", colony_id, "_sweep_drift_tests_Nemin.tsv")
  
  fwrite(Ne_all_modes,   Ne_out_file,     sep = "\t")
  fwrite(sweeps_c_tests, sweeps_out_file, sep = "\t")
  
  invisible(list(
    Ne_all_modes = Ne_all_modes,
    sweeps       = sweeps_c_tests
  ))
}

#driver

traj_files <- list.files(
  traj_dir,
  pattern = "^snv_frequency_trajectory_\\d+\\.csv$",
  full.names = TRUE
)
colony_ids <- as.integer(sub("^.*snv_frequency_trajectory_(\\d+)\\.csv$", "\\1", traj_files))

all_results <- lapply(colony_ids, run_colony_all_modes)




#### summary report gen ####



## 1. Discover colony output files
files <- list.files(out_dir, full.names = TRUE)

ne_files    <- files[grepl(ne_pattern, basename(files))]
sweep_files <- files[grepl(sweep_pattern, basename(files))]

get_colony_id <- function(x, pattern) {
  as.integer(sub(pattern, "\\1", basename(x)))
}

ne_colonies    <- get_colony_id(ne_files,    ne_pattern)
sweep_colonies <- get_colony_id(sweep_files, sweep_pattern)

common_colonies <- intersect(ne_colonies, sweep_colonies)
cat("Colonies with both Ne and sweep outputs:", paste(common_colonies, collapse = ", "), "\n")

## 2. Load and combine Ne tables and sweep tables

all_Ne <- rbindlist(lapply(common_colonies, function(c_id) {
  f <- ne_files[ne_colonies == c_id][1]
  dt <- fread(f)
  dt[, colony_id := c_id]
  dt
}), fill = TRUE)

all_sweeps <- rbindlist(lapply(common_colonies, function(c_id) {
  f <- sweep_files[sweep_colonies == c_id][1]
  dt <- fread(f)
  dt[, colony_id := c_id]
  dt
}), fill = TRUE)

cat("Total MAG × colony Ne rows:", nrow(all_Ne), "\n")
cat("Total sweep SNV rows:", nrow(all_sweeps), "\n")

## 3. Add FDR-corrected p-values to sweep tests

all_sweeps_fdr <- copy(all_sweeps)

# FDR within each MAG × colony (this is the main family we care about)
all_sweeps_fdr[!is.na(p),
               p_adj_magcol := p.adjust(p, method = "BH"),
               by = .(colony_id, MAG_id)]

# Optional: FDR across all sweeps in each colony
all_sweeps_fdr[!is.na(p),
               p_adj_colony := p.adjust(p, method = "BH"),
               by = colony_id]

## 4. MAG × colony summary table

sweep_summary_mag_colony <- all_sweeps_fdr[!is.na(p), .(
  n_sweeps_tested   = .N,
  n_sweeps_p_0.05   = sum(p < 0.05, na.rm = TRUE),
  n_sweeps_p_0.01   = sum(p < 0.01, na.rm = TRUE),
  n_sweeps_p_0.001  = sum(p < 0.001, na.rm = TRUE),
  # FDR-based counts using BH within MAG × colony
  n_sweeps_q_0.10   = sum(p_adj_magcol < 0.10, na.rm = TRUE),
  n_sweeps_q_0.05   = sum(p_adj_magcol < 0.05, na.rm = TRUE),
  n_sweeps_q_0.01   = sum(p_adj_magcol < 0.01, na.rm = TRUE),
  median_p          = median(p, na.rm = TRUE),
  median_q_magcol   = median(p_adj_magcol, na.rm = TRUE),
  median_Zmax       = median(Z_max, na.rm = TRUE)
), by = .(colony_id, MAG_id)]

# Merge with Ne_all_modes (which should have Ne_4D, Ne_non0D, Ne_all, Ne_min, etc.)
mag_colony_report <- merge(
  all_Ne,
  sweep_summary_mag_colony,
  by = c("colony_id", "MAG_id"),
  all.x = TRUE
)

# Replace NA counts with 0 where appropriate
for (col in c("n_sweeps_tested", "n_sweeps_p_0.05",
              "n_sweeps_p_0.01", "n_sweeps_p_0.001",
              "n_sweeps_q_0.10", "n_sweeps_q_0.05", "n_sweeps_q_0.01")) {
  mag_colony_report[is.na(get(col)), (col) := 0L]
}

fwrite(
  mag_colony_report,
  file = file.path(out_dir, "MAG_colony_drift_summary.tsv"),
  sep  = "\t"
)

## 5. Overall per-MAG summary across colonies

mag_overall_summary <- mag_colony_report[, .(
  n_colonies_with_Ne = sum(!is.na(Ne_min)),
  median_Ne_min      = median(Ne_min, na.rm = TRUE),
  q25_Ne_min         = quantile(Ne_min, 0.25, na.rm = TRUE),
  q75_Ne_min         = quantile(Ne_min, 0.75, na.rm = TRUE),
  
  total_sweeps_tested = sum(n_sweeps_tested,   na.rm = TRUE),
  total_sweeps_p_0.05 = sum(n_sweeps_p_0.05,   na.rm = TRUE),
  total_sweeps_p_0.01 = sum(n_sweeps_p_0.01,   na.rm = TRUE),
  total_sweeps_p_0.001= sum(n_sweeps_p_0.001,  na.rm = TRUE),
  total_sweeps_q_0.10 = sum(n_sweeps_q_0.10,   na.rm = TRUE),
  total_sweeps_q_0.05 = sum(n_sweeps_q_0.05,   na.rm = TRUE),
  total_sweeps_q_0.01 = sum(n_sweeps_q_0.01,   na.rm = TRUE)
), by = MAG_id]

mag_overall_summary[, frac_sweeps_p_0.05  := ifelse(total_sweeps_tested > 0,
                                                    total_sweeps_p_0.05 / total_sweeps_tested, NA_real_)]
mag_overall_summary[, frac_sweeps_q_0.05  := ifelse(total_sweeps_tested > 0,
                                                    total_sweeps_q_0.05 / total_sweeps_tested, NA_real_)]

fwrite(
  mag_overall_summary,
  file = file.path(out_dir, "MAG_overall_drift_summary.tsv"),
  sep  = "\t"
)

## 6. Per-MAG sweep reports (with raw & FDR p-values)

unique_MAGS <- sort(unique(all_sweeps_fdr$MAG_id))

dir.create(file.path(out_dir, "MAG_sweep_reports"), showWarnings = FALSE)

for (mag in unique_MAGS) {
  dt_mag <- all_sweeps_fdr[MAG_id == mag]
  dt_mag <- merge(
    dt_mag,
    all_Ne[, .(colony_id, MAG_id, Ne_min)],
    by = c("colony_id", "MAG_id"),
    all.x = TRUE
  )
  out_path <- file.path(out_dir, "null_model_for_drift/MAG_sweep_reports",
                        paste0("MAG_", mag, "_sweep_drift_report.tsv"))
  fwrite(dt_mag, out_path, sep = "\t")
}







mc <- fread("MAG_colony_drift_summary.tsv")

mc[, frac_sweeps_q_0.05 := ifelse(
  n_sweeps_tested > 0,
  n_sweeps_q_0.05 / n_sweeps_tested,
  NA_real_
)]

# Optional: order MAGs by overall selection intensity
mc[, MAG_id := factor(MAG_id)]
mag_order <- mc[, .(
  overall_frac_q0.05 = sum(n_sweeps_q_0.05, na.rm = TRUE) /
    sum(n_sweeps_tested,   na.rm = TRUE)
), by = MAG_id][order(-overall_frac_q0.05)]$MAG_id

mc[, MAG_id := factor(MAG_id, levels = mag_order)]
mc[, colony_id := factor(colony_id)]  # treat colonies as categorical


ggplot(mc, aes(x = colony_id, y = MAG_id, fill = frac_sweeps_q_0.05)) +
  geom_tile(color = "grey80") +
  scale_fill_viridis_c(
    option = "C",
    na.value = "white",
    name = "Frac sweeps\nq<0.05"
  ) +
  labs(
    x = "Colony (host)",
    y = "MAG (species)",
    title = "Fraction of sweep SNVs per MAG × colony that reject drift (FDR q<0.05)",
    subtitle = "Ne_min per MAG used as maximally generous drift baseline"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )





mo <- fread("MAG_overall_drift_summary.tsv")

ggplot(mo, aes(x = median_Ne_min, y = frac_sweeps_q_0.05)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  labs(
    x = "Median Ne_min across colonies (log scale)",
    y = "Fraction of sweep SNVs with q<0.05 (MAG-level)",
    title = "Species-level summary of selection vs drift"
  ) +
  theme_minimal(base_size = 12)+
  coord_flip()





report_dir <- "MAG_sweep_reports"
report_files <- list.files(report_dir, pattern = "^MAG_.*_sweep_drift_report.tsv$", full.names = TRUE)
all_reports <- rbindlist(lapply(report_files, fread), fill = TRUE)

# Filter to sweeps with valid p
all_reports_valid <- all_reports[!is.na(p)]
dt <- as.data.table(all_reports_valid)
# Clamp zeros to a small value so log10 works, but label them as "0"
eps <- 1e-6
dt[, p_adj_clamped := ifelse(p_adj_magcol == 0, eps, p_adj_magcol)]



ggplot(dt, aes(x = p_adj_clamped)) +
  geom_histogram(bins = 40, boundary = 0, closed = "left") +
  scale_x_log10(
    breaks = c(eps, 1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = c("0", "0.0001", "0.001", "0.01", "0.1", "1"),
    limits = c(eps, 1)
  ) +
  labs(
    x = "FDR-adjusted p (BH within MAG × colony)",
    y = "Count of sweep SNVs",
    title = "Distribution of q-values for drift test across all sweep SNVs"
  ) +
  theme_minimal(base_size = 12)+
  geom_vline(xintercept = 0.05)

holder <- dt[dt$p_adj_magcol<0.05,]
print(19522/25888)

top_mags <- mo[order(-total_sweeps_tested)][1:6, MAG_id]

plot_data <- all_reports_valid[MAG_id %in% top_mags]

ggplot(plot_data, aes(x = p, colour = MAG_id)) +
  stat_ecdf() +
  labs(
    x = "Raw drift test p-value",
    y = "ECDF",
    title = "Per-MAG ECDF of drift test p-values for sweep SNVs"
  ) +
  theme_minimal(base_size = 12)



volc_df <- read.csv('/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/private_marker_results_all_colonies/private_marker_summary_all_colonies.csv')

volc_df <- as.data.table(volc_df)
# Compute log2 ratio (effect size), guarding against zeros
eps <- 1e-6
volc_df[, log2_ratio := log2((Frac_preserved + eps) / (Frac_disrupted + eps))]


ggplot(volc_df, aes(x = Frac_preserved, y = Total_private)) +
  geom_vline(xintercept = 0.5, colour = "grey60", linetype = "dashed") +
  geom_point(size = 3) +
  scale_colour_viridis_c(option = "C", name = "Frac preserved") +
  labs(
    x = "Fraction of private sweeps preserved",
    y = "Total number of private sweeps",
    title = "Balance of preserved vs disrupted private sweeps per MAG",
    subtitle = "Points to the right of 0.5: more preserved than disrupted"
  ) +
  theme_minimal(base_size = 12)+
  scale_y_continuous(trans="log10")

















add_time_from_month <- function(dt) {
  dt[, Time := month_to_cumulative_time_map[Month]]
  dt
}

## --- Read all sweeps once ---
sweeps_all <- fread(sweep_file)

## Helper: build per-interval table for *sweep* SNVs in one colony
build_sweep_intervals_for_colony <- function(col_id, out_dir) {
  cat("\n==== Colony", col_id, "====\n")
  
  ## ---------- File paths ----------
  traj_file  <- file.path(
    traj_dir,
    sprintf("snv_frequency_trajectory_%d.csv", col_id)
  )
  ne_file    <- file.path(
    out_dir,
    sprintf("colony_%d_Ne_by_MAG_allmodes.tsv", col_id)
  )
  sweep_file <- file.path(
    out_dir,
    sprintf("colony_%d_sweep_drift_tests_Nemin.tsv", col_id)
  )
  
  if (!file.exists(traj_file)) {
    cat("  Trajectory file missing, skipping colony", col_id, "\n")
    return(NULL)
  }
  if (!file.exists(ne_file)) {
    cat("  Ne file missing, skipping colony", col_id, "\n")
    return(NULL)
  }
  if (!file.exists(sweep_file)) {
    cat("  Sweep drift file missing, skipping colony", col_id, "\n")
    return(NULL)
  }
  
  ## ---------- Load trajectories ----------
  traj <- fread(traj_file)
  # keep only what we need
  traj <- traj[
    ,
    .(MAG_id,
      Month,
      contig_name,
      pos_in_contig,
      coverage,
      departure_from_polarized_consensus)
  ]
  traj <- add_time_from_month(traj)
  traj[, GenTime := Time * gens_per_day]
  setorder(traj, MAG_id, contig_name, pos_in_contig, GenTime)
  
  ## ---------- Load Ne (per MAG) ----------
  Ne_tab <- fread(ne_file)
  if (!"MAG_id" %in% names(Ne_tab)) {
    stop("Ne file for colony ", col_id, " has no MAG_id column. Names: ",
         paste(names(Ne_tab), collapse = ", "))
  }
  if (!"Ne_min" %in% names(Ne_tab)) {
    stop("Ne file for colony ", col_id, " has no Ne_min column. Names: ",
         paste(names(Ne_tab), collapse = ", "))
  }
  Ne_tab <- Ne_tab[, .(MAG_id, Ne_min)]
  
  ## ---------- Load sweeps ----------
  sweeps_all <- fread(sweep_file)
  
  # We do NOT need colony_id from the file; this file is already for a single colony.
  # Just make sure the locus columns exist:
  locus_cols <- c("MAG_id", "contig_name", "pos_in_contig")
  # In your per-colony sweep files, these may still be called 'contig' and 'position':
  if ("contig" %in% names(sweeps_all) && "position" %in% names(sweeps_all)) {
    setnames(sweeps_all,
             old = c("contig", "position"),
             new = c("contig_name", "pos_in_contig"))
  }
  
  if (!all(locus_cols %in% names(sweeps_all))) {
    stop("Sweep file for colony ", col_id,
         " is missing locus columns. Names: ",
         paste(names(sweeps_all), collapse = ", "))
  }
  
  # Unique sweep SNVs
  sweeps_c <- unique(
    sweeps_all[, .(MAG_id, contig_name, pos_in_contig)],
    by = c("MAG_id", "contig_name", "pos_in_contig")
  )
  
  if (!nrow(sweeps_c)) {
    cat("  No sweeps in sweep file for colony", col_id, "\n")
    return(NULL)
  }
  
  ## ---------- Attach Ne_min per MAG ----------
  sweeps_c <- merge(
    sweeps_c,
    Ne_tab,
    by = "MAG_id",
    all.x = TRUE
  )
  
  sweeps_c <- sweeps_c[!is.na(Ne_min) & Ne_min > 0]
  if (!nrow(sweeps_c)) {
    cat("  No sweeps with finite Ne_min for colony", col_id, "\n")
    return(NULL)
  }
  
  ## ---------- Build per-interval df, Z, etc. for each sweep SNV ----------
  res <- sweeps_c[
    ,
    {
      # subset the trajectory for this SNV
      snv_traj <- traj[
        MAG_id == .BY$MAG_id &
          contig_name == .BY$contig_name &
          pos_in_contig == .BY$pos_in_contig
      ]
      
      if (nrow(snv_traj) < 2L) {
        return(NULL)
      }
      
      snv_traj[, f := departure_from_polarized_consensus]
      snv_traj[
        ,
        `:=`(
          GenTime_next = shift(GenTime, type = "lead"),
          f_next       = shift(f,       type = "lead"),
          cov_next     = shift(coverage, type = "lead")
        )
      ]
      
      dt_int <- snv_traj[
        !is.na(GenTime_next) &
          coverage >= coverage_min &
          cov_next >= coverage_min
      ]
      if (!nrow(dt_int)) {
        return(NULL)
      }
      
      eps <- 1e-6
      dt_int[, f_clamp := pmin(pmax(f, eps), 1 - eps)]
      dt_int[, tau := GenTime_next - GenTime]
      dt_int[, df  := f_next - f]
      dt_int[, V_sample := f_clamp * (1 - f_clamp) * (1/coverage + 1/cov_next)]
      dt_int[, V_drift  := f_clamp * (1 - f_clamp) *
               (1 - exp(-tau / (2 * Ne_min)))]
      dt_int[, V_tot := V_sample + V_drift]
      
      dt_int <- dt_int[is.finite(V_tot) & V_tot > 0]
      if (!nrow(dt_int)) {
        return(NULL)
      }
      
      dt_int[, Z := df / sqrt(V_tot)]
      
      dt_int[
        ,
        .(
          colony_id     = col_id,
          MAG_id        = .BY$MAG_id,
          contig_name   = .BY$contig_name,
          pos_in_contig = .BY$pos_in_contig,
          Month,
          GenTime,
          GenTime_next,
          df,
          Z,
          tau,
          V_sample,
          V_drift,
          V_tot,
          Ne_min       = Ne_min
        )
      ]
    },
    by = .(MAG_id, contig_name, pos_in_contig)
  ]
  
  res
}


## --- Build one big table across all colonies ---

# Colonies where you actually have Ne & sweeps already summarized
# (re-use whatever 'common_colonies' and 'out_dir' you set above)
# If you still have 'common_colonies' in memory from your summary code, use that; 
# otherwise reconstruct it:
out_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/null_model_for_drift"

ne_pattern    <- "^colony_(\\d+)_Ne_by_MAG_allmodes.tsv$"
sweep_pattern <- "^colony_(\\d+)_sweep_drift_tests_Nemin.tsv$"

files <- list.files(out_dir, full.names = TRUE)
tax <- fread("/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv")[
  , .(MAG_id, Genus)
]
tax <- unique(tax, by = "MAG_id")

ne_files    <- files[grepl(ne_pattern,    basename(files))]
sweep_files <- files[grepl(sweep_pattern, basename(files))]

get_colony_id <- function(x, pattern) {
  as.integer(sub(pattern, "\\1", basename(x)))
}

ne_colonies    <- get_colony_id(ne_files,    ne_pattern)
sweep_colonies <- get_colony_id(sweep_files, sweep_pattern)

common_colonies <- intersect(ne_colonies, sweep_colonies)
cat("Colonies with both Ne and sweep files:",
    paste(common_colonies, collapse = ", "), "\n")

intervals_sweeps <- rbindlist(
  lapply(common_colonies, build_sweep_intervals_for_colony, out_dir = out_dir),
  fill = TRUE
)

cat("Total sweep intervals:", nrow(intervals_sweeps), "\n")
# 1. Identify which column names are duplicates
dup_colnames <- duplicated(colnames(intervals_sweeps))

# 2. Identify the indices (numbers) of the columns to KEEP (those that are NOT duplicated)
keep_indices <- which(!dup_colnames)

# 3. Subsetting the data table using the indices
intervals_sweeps <- intervals_sweeps[, keep_indices, with = FALSE] 
# 'with = FALSE' tells data.table to treat the numbers as column indices, not names.

intervals_sweeps[
  # i: The right-hand side table ('tax') is joined to the left-hand side (intervals_sweeps)
  tax, 
  
  # j: We create a new column called 'Genus' (:=) in intervals_sweeps, 
  # and assign it the value of Genus from the 'tax' table (i.Genus).
  Genus := i.Genus, 
  
  # on: The column to join on. This is where the lookup happens.
  on = "MAG_id"
]

head(intervals_sweeps)



intervals_sweeps[, abs_df := abs(df)]
intervals_sweeps[, exp_sd := sqrt(V_tot)]
intervals_sweeps <- intervals_sweeps[is.finite(exp_sd) & is.finite(abs_df)]



ggplot(intervals_sweeps, aes(x = Z)) +
  geom_histogram(aes(y = ..density..),
                 bins = 60,
                 boundary = 0,
                 closed = "left") +
  stat_function(
    fun  = dnorm,
    args = list(mean = 0, sd = 1),
    colour = "red",
    linewidth = 1
  ) +
  labs(
    x = expression(paste(Delta, "f / ", sqrt(V[tot]))),
    y = "Density",
    title = "Standardized allele-frequency changes for sweep SNVs",
    subtitle = "Black = observed; red = N(0,1) drift-only expectation"
  ) +
  theme_minimal(base_size = 12)



ggplot(intervals_sweeps, aes(x = exp_sd, y = abs_df)) +
  geom_hex() +  # or geom_point(alpha = 0.1) if you prefer
  # Under drift, E|Z| = sqrt(2/pi), so E|Δf| ~ sqrt(2/pi) * exp_sd
  geom_abline(
    slope     = sqrt(2 / pi),
    intercept = 0,
    linetype  = "dashed",
    colour    = "red"
  ) +
  scale_x_continuous(name = "Expected SD of Δf under drift + sampling (sqrt(V_tot))") +
  scale_y_continuous(name = "Observed |Δf| per interval") +
  labs(
    title    = "Observed allele-frequency changes vs drift-only expectation",
    subtitle = "Dashed line = expected mean |Δf| under drift (sqrt(2/pi) * sqrt(V_tot))"
  ) +
  theme_minimal(base_size = 12)




genus_summary <- intervals_sweeps[
  ,
  .(
    n_intervals   = .N,
    med_abs_df    = median(abs_df, na.rm = TRUE),
    med_exp_sd    = median(exp_sd, na.rm = TRUE),
    med_ratio     = median(abs_df / exp_sd, na.rm = TRUE),
    frac_Z_gt_2   = mean(abs(Z) > 2, na.rm = TRUE)  # optional
  ),
  by = Genus
][n_intervals >= 50]  # require some minimum data per Genus

ggplot(genus_summary,
       aes(x = reorder(Genus, med_ratio), y = med_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Genus",
    y = "Median |Δf| / expected SD under drift",
    title = "Excess allele-frequency change relative to drift-only null, by Genus",
    subtitle = "Values > 1 indicate more change than expected from drift + sampling"
  ) +
  theme_minimal(base_size = 12)









# 1. Calculate Expected Delta f (simulated)
# Repeat each expected standard deviation (exp_sd) for simulation
intervals_sweeps[, abs_df := abs(df)] # Ensure abs_df is still present

# For a drift-only null model, the change in allele frequency (df) is normally distributed
# with mean 0 and variance V_tot. We can simulate the magnitude of this change.
# Let's simulate one expected delta_f (magnitide) for every observed interval:
intervals_sweeps[, 
                 expected_abs_df := abs(rnorm(
                   .N, 
                   mean = 0, 
                   sd = exp_sd
                 ))
]

# 2. Reshape the data for ggplot (long format)
plot_data <- data.table(
  Observed = intervals_sweeps$abs_df,
  Drift_Expected = intervals_sweeps$expected_abs_df
)
plot_data_long <- melt(plot_data, 
                       measure.vars = c("Observed", "Drift_Expected"),
                       variable.name = "Model",
                       value.name = "Absolute_Delta_f")

# 3. Create the Density Plot
ggplot(plot_data_long, aes(x = Absolute_Delta_f, fill = Model, colour = Model)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values = c("Observed" = "black", "Drift_Expected" = "red")) +
  scale_colour_manual(values = c("Observed" = "black", "Drift_Expected" = "red")) +
  xlim(0, quantile(plot_data_long$Absolute_Delta_f, 0.99)) + # Focus on the main distribution
  labs(
    x = expression(paste("Absolute Allele Frequency Change ", "|", Delta, "f|")),
    y = "Density",
    title = "Observed vs. Drift-Only Expected Rates of Allele Frequency Change"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


plot_data_long$Absolute_Delta_f <- plot_data_long$Absolute_Delta_f + 0.0000001
plot_data_long$ <- plot_data_long$Absolute_Delta_f + 0.0000001

# Reshape the data for ggplot (long format) and include Genus
plot_data_long <- data.table(
  Genus = intervals_sweeps$Genus,
  Observed = intervals_sweeps$abs_df,
  Drift_Expected = intervals_sweeps$expected_abs_df
)
plot_data_long <- melt(plot_data_long, 
                       id.vars = "Genus",
                       measure.vars = c("Observed", "Drift_Expected"),
                       variable.name = "Model",
                       value.name = "Absolute_Delta_f")

# Filter out NAs if any Genus data is missing (optional cleanup)
plot_data_long <- plot_data_long[!is.na(Genus)]

# 2. Create the Faceted Density Plot with Log-Scale
ggplot(plot_data_long, aes(x = Absolute_Delta_f, fill = Model, colour = Model)) +
  
  # Use geom_density for smooth distribution curves
  geom_density(alpha = 0.5, linewidth = 0.8) +
  
  # Apply Log10 Transformation to the x-axis
  scale_x_log10(
    labels = scales::label_number(), # Use scientific notation or standard numbers
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1) 
  ) +
  
  # Separate plots for each Genus
  facet_wrap(~ Genus, scales = "free_y") + 
  
  scale_fill_manual(values = c("Observed" = "black", "Drift_Expected" = "red")) +
  scale_colour_manual(values = c("Observed" = "black", "Drift_Expected" = "red")) +
  
  labs(
    x = expression(paste("Absolute Allele Frequency Change ", "|", Delta, "f| (log scale)")),
    # Update y-axis label to reflect the log transformation
    y = "Density (log scale)", 
    title = "Observed vs. Drift-Only Expected Rates by Genus",
    subtitle = "Deviation (black) from neutral drift (red) indicates selective sweeps"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        # Shrink facet titles slightly if many genera are present
        strip.text = element_text(size = 9, face = "bold"))




ks_result <- ks.test(
  x = intervals_sweeps$abs_df, 
  y = intervals_sweeps$expected_abs_df
)

print(ks_result)



genus_ks_results <- intervals_sweeps %>%
  group_by(Genus) %>%
  # Use summarise to apply the KS test to each group
  dplyr::summarise(
    N_intervals = n(),
    KS_D_statistic = ks.test(abs_df, expected_abs_df)$statistic,
    KS_P_value = ks.test(abs_df, expected_abs_df)$p.value
  ) %>%
  # Convert to data.table for final viewing if you prefer
  data.table()

print(genus_ks_results)
