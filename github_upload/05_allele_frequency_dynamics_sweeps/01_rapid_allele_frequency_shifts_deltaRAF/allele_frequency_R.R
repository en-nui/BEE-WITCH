# --- Libraries ---
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(patchwork)
library(tidyr)
library(grid)

# --- 0. Configuration ---
SNV_DIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/snv_frequency_trajectories/"
OUTPUT_PLOT_DIR <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/plots/"
if (!dir.exists(OUTPUT_PLOT_DIR)) dir.create(OUTPUT_PLOT_DIR, recursive = TRUE)

# --- 1. Color palette ---
gruvbox_palette <- c(
  dark_bg = "#282828",
  blue    = "#83a598",
  orange  = "#fe8019",
  gray_bg = "#928374"
)

month_order <- c("May","June","July","August","September","October","November","January","February")

# --- 2. Load Ecology Data ---
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"

ecology_df <- read_csv(metadata_file) %>%
  filter(type == "mmag", Caste == "Worker")

mag_abundances <- ecology_df %>%
  group_by(Colony, Month, sample_name, MAG_id, Genus) %>%
  summarise(mean_coverage = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(Colony, Month, sample_name) %>%
  mutate(total_coverage = sum(mean_coverage),
         rep_rel_abund = mean_coverage / total_coverage) %>%
  ungroup()

message("✅ Ecology data processed with replicate information.")

# --- 3. Load Significant SNVs ---
significant_snvs_file <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/sweep_analysis/significant_sweep_snvs_ALL.csv"
significant_snvs_df <- read_csv(significant_snvs_file) %>%
  mutate(t1 = factor(t1, levels = month_order),
         t2 = factor(t2, levels = month_order))
message(paste("✅ Loaded", nrow(significant_snvs_df), "significant SNVs."))

# --- 4. Unique MAG × Colony combinations ---
all_candidates <- significant_snvs_df %>%
  distinct(colony_id, MAG_id) %>%
  arrange(colony_id, MAG_id)
message(paste0("\n📊 Generating plots for ", nrow(all_candidates), " MAG × colony combinations."))

global_max_rel_abund <- max(mag_abundances$rep_rel_abund, na.rm = TRUE)

# --- 5. Continuous time helpers ---
add_cumulative_time_from_month <- function(data) {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34,
                        "September" = 28, "October" = 22, "November" = 42,
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  data$Time <- month_to_cumulative_time_map[data$Month]
  return(data)
}

highlight_colors <- c(June="#fabd2f", October="#fe8019", January="#83a598")
get_month_annotations <- function(month_order, highlight_colors, min_time, max_time) {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34,
                        "September" = 28, "October" = 22, "November" = 42,
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  annos <- list()
  for (m in names(highlight_colors)) {
    if (m %in% names(month_to_cumulative_time_map)) {
      xmid <- month_to_cumulative_time_map[m]
      if (xmid >= min_time & xmid <= max_time) {
        annos[[length(annos) + 1]] <-
          annotate("rect", xmin = xmid - 15, xmax = xmid + 15,
                   ymin = -Inf, ymax = Inf,
                   fill = highlight_colors[m], alpha = 0.3)
      }
    }
  }
  annos
}

# --- 6. Main loop ---
plot_list <- list()

for (i in 1:nrow(all_candidates)) {
  target_colony_id <- all_candidates$colony_id[i]
  target_mag_id <- all_candidates$MAG_id[i]
  message(paste0("\n[", i, "/", nrow(all_candidates), "] Colony ", target_colony_id, " — MAG ", target_mag_id))
  
  # --- A. Abundance ---
  mag_rel_abund_data <- mag_abundances %>%
    filter(Colony == target_colony_id, MAG_id == target_mag_id) %>%
    mutate(Month = factor(Month, levels = month_order, ordered = TRUE),
           rel_abund_mag = rep_rel_abund + 1e-6)
  if (nrow(mag_rel_abund_data) == 0) next
  
  mag_rel_abund_data <- add_cumulative_time_from_month(mag_rel_abund_data)
  mag_rel_abund_summary <- mag_rel_abund_data %>%
    group_by(Month, Time) %>%
    summarise(mean_rel_abund = mean(rel_abund_mag, na.rm = TRUE),
              sd_rel_abund   = sd(rel_abund_mag, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(sd_rel_abund = ifelse(is.na(sd_rel_abund), 0, sd_rel_abund),
           ymin = pmax(mean_rel_abund - sd_rel_abund, 1e-6),
           ymax = pmax(mean_rel_abund + sd_rel_abund, 1e-6))
  
  min_time <- min(mag_rel_abund_summary$Time, na.rm = TRUE)
  max_time <- max(mag_rel_abund_summary$Time, na.rm = TRUE)
  month_annos <- get_month_annotations(month_order, highlight_colors, min_time, max_time)
  
  p1 <- ggplot(mag_rel_abund_summary,
               aes(x = Time, y = mean_rel_abund, group = 1)) +
    month_annos +
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
                fill = gruvbox_palette["orange"], alpha = 0.3) +
    geom_line(color = "black", linewidth = 0.8) +
    geom_point(color = "black", size = 2) +
    scale_x_continuous(
      breaks = seq(min_time, max_time, by = 20),
      labels = seq(min_time, max_time, by = 20),
      expand = expansion(mult = c(0.02, 0.02))   # <-- small buffer (~2%)
    )+
    scale_y_continuous(
      trans = "log10",
      breaks = c(1, 1e-1, 1e-2, 1e-3),
      labels = c(expression(10^0), expression(10^-1),
                 expression(10^-2), expression(10^-3))
    ) +
    coord_cartesian(ylim = c(1e-3, 1)) +
    labs(x = NULL, y = "Relative Abundance (log10)") +
    theme_minimal() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      plot.margin  = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  # --- B. SNVs ---
  snv_file_path <- file.path(SNV_DIR, paste0("snv_frequency_trajectory_", target_colony_id, ".csv"))
  if (!file.exists(snv_file_path)) next
  
  full_colony_snv_df <- read_csv(snv_file_path, show_col_types = FALSE) %>%
    filter(MAG_id == target_mag_id) %>%
    mutate(Month = factor(Month, levels = month_order, ordered = TRUE),
           unique_SNV_identifier = paste(MAG_id, contig_name, pos_in_contig, sep = "_")) %>%
    filter(!is.na(departure_from_polarized_consensus), coverage > 0)
  if (nrow(full_colony_snv_df) == 0) next
  
  full_colony_snv_df <- add_cumulative_time_from_month(full_colony_snv_df)
  
  mag_sig_snvs <- significant_snvs_df %>%
    filter(colony_id == target_colony_id, MAG_id == target_mag_id) %>%
    pull(unique_SNV_identifier) %>% unique()
  
  significant_snvs_for_plot <- full_colony_snv_df %>%
    filter(unique_SNV_identifier %in% mag_sig_snvs)
  
  non_significant_snvs <- full_colony_snv_df %>%
    filter(!(unique_SNV_identifier %in% mag_sig_snvs),
           departure_from_polarized_consensus < 0.5)
  if (n_distinct(non_significant_snvs$unique_SNV_identifier) > 400) {
    sampled_non_sig <- non_significant_snvs %>%
      distinct(unique_SNV_identifier) %>%
      sample_n(400) %>%
      pull(unique_SNV_identifier)
    non_significant_snvs_for_plot <- non_significant_snvs %>%
      filter(unique_SNV_identifier %in% sampled_non_sig)
  } else non_significant_snvs_for_plot <- non_significant_snvs
  
  p2 <- ggplot() +
    month_annos +
    geom_line(data = non_significant_snvs_for_plot,
              aes(x = Time, y = departure_from_polarized_consensus,
                  group = unique_SNV_identifier),
              color = gruvbox_palette["gray_bg"], alpha = 0.25, linewidth = 0.25) +
    geom_line(data = significant_snvs_for_plot,
              aes(x = Time, y = departure_from_polarized_consensus,
                  group = unique_SNV_identifier),
              color = gruvbox_palette["blue"], alpha = 0.35, linewidth = 0.4) +
    scale_x_continuous(
      breaks = seq(min_time, max_time, by = 20),
      labels = seq(min_time, max_time, by = 20),
      expand = expansion(mult = c(0.02, 0.02))   # <-- small buffer (~2%)
    ) +
    ylim(0, 1) +
    labs(x = "Days since start of season", y = "Allele Frequency") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      plot.margin  = margin(t = 0, r = 5, b = 5, l = 5)
    )
  
  # --- C. Spacer and combine ---
  spacer <- ggplot() + theme_void() + theme(plot.margin = margin(0, 2, 0, 2))
  combined_plot <- patchwork::wrap_plots(p1, spacer, p2, ncol = 1,
                                         heights = c(1.2, 0.15, 2.3)) +
    patchwork::plot_annotation(
      title = paste0("Colony ", target_colony_id, " - MAG ", target_mag_id, ": Abundance & SNV Sweeps"),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    ) &
    theme(plot.spacing = unit(0, "cm"), plot.margin = margin(0, 0, 0, 0))
  
  # --- D. Save ---
  svg_file <- file.path(
    OUTPUT_PLOT_DIR,
    paste0("colony_", target_colony_id, "_mag_", target_mag_id, "_sweep_visualization.svg")
  )
  grDevices::svg(svg_file, width = 8, height = 10)
  print(combined_plot)
  dev.off()
  message(paste("  ✅ Saved:", svg_file))
  
  plot_list[[i]] <- combined_plot
}

message("\n🎉 All MAG × colony plots generated with dynamic continuous time axes, adaptive shading, and fixed tick spacing (20 days).")

