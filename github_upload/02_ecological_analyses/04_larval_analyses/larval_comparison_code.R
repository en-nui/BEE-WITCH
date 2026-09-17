library(ggplot2)
library(dplyr)
library(tidyr)


cat("--- Script Started ---\n")

# --- 2. LOAD AND INITIAL DATA PROCESSING (GLOBAL) ---
# Load the concatenated dataframe
data_path <- '/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output/all_filtered_mmag_contigs_concatenated.csv'
df <- read_csv(data_path)

# Clean column names to lowercase for consistency
df <- df %>%
  rename_with(tolower)
df <- df[df$caste=="Larvae",]
df<-df[df$contig_presence_ratio<0.5,]
# --- Define month order for consistent plotting ---
full_month_order <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")

# Initial filtering and addition of the 'apiary' column
df_processed <- df %>%
  filter(type == "mmag") %>%
  mutate(apiary = case_when(
    colony > 500 ~ "HB",
    colony <= 103 ~ "PT",
    TRUE ~ "OB"
  ))

# IMPORTANT: Ensure each MAG is counted only once per sample by summarizing coverage.
# This step prevents double-counting and creates a clean base for abundance calculations.
df_processed <- df_processed %>%
  group_by(mag_id, month, colony, rep, caste, genus, apiary, contig_presence_ratio) %>%
  summarise(
    avg_cov = first(average_genome_wide_cov),
    .groups = 'drop'
  )


# --- 3. DEFINE TOP GENERA FOR CONSISTENT COLORING (GLOBAL) ---
# We will identify the most abundant genera to give them unique colors in the plots.
# For this example, we'll define top genera based on overall abundance across the dataset.
top_genera <- df_processed %>%
  filter(!is.na(genus) & genus != "no match") %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(avg_cov, na.rm = TRUE)) %>%
  slice_max(order_by = total_abundance, n = 10) %>% # Select top 10 most abundant
  pull(genus)

cat(paste("\nTop 10 genera selected for unique coloring:", paste(top_genera, collapse = ", "), "\n"))

# --- 4. CREATE FINAL TAXONOMIC CATEGORIES (GLOBAL) ---
# Group less abundant or low-quality MAGs into categories for cleaner plots.
df_processed <- df_processed %>%
  mutate(
    final_genus_category = case_when(
      contig_presence_ratio < 0.5 ~ "MAG Below Threshold",
      genus %in% top_genera ~ genus,
      is.na(genus) | genus == "" | genus == "no match" ~ "Unknown/Other",
      TRUE ~ "Unknown/Other" # All other genera are grouped
    )
  )

# --- 5. SETUP PLOTTING ENVIRONMENT ---
# Define the base directory for saving figures
figure_output_base_dir <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/figures/larvae_by_colony_rep_test"
dir.create(figure_output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Define a color palette
# Create a named list of colors for consistent plotting
# Manually assign colors to your top genera and fixed categories
# We'll generate a palette for the top genera and add fixed colors for the rest
color_palette_dynamic <- scales::hue_pal()(length(top_genera))
names(color_palette_dynamic) <- top_genera

# Combine with fixed colors for special categories
plot_colors <- c(
  color_palette_dynamic,
  "MAG Below Threshold" = "grey50",
  "Unknown/Other" = "grey80"
)



# --- 6. FILTER FOR LARVAE AND GENERATE COMPARATIVE PLOTS ---

# Filter the dataset for Larvae only for this analysis
df_larvae <- df_processed %>%
  filter(caste == "Larvae")

# Get unique Colonies to loop through
unique_colonies <- unique(df_larvae$colony) %>% sort()

cat("\n--- Starting Comparative Plot Generation for Larvae ---\n")

for (current_colony in unique_colonies) {
  
  cat(paste("Processing: Colony", current_colony, "\n"))

  # Filter data for the current colony
  df_subset <- df_larvae %>%
    filter(colony == current_colony)

  if (nrow(df_subset) == 0) {
    cat("  No data. Skipping.\n")
    next
  }

  # --- Calculate Relative Abundance (now grouped by month AND rep) ---
  df_rel_abundance <- df_subset %>%
    # Group by both month and rep to get the correct total for each individual bar
    group_by(month, rep) %>%
    mutate(total_sample_abundance = sum(avg_cov, na.rm = TRUE)) %>%
    ungroup() %>%
    # Calculate relative abundance for each genus category within each sample
    group_by(month, rep, final_genus_category, total_sample_abundance) %>%
    summarise(genus_sample_abundance = sum(avg_cov, na.rm = TRUE), .groups = 'drop') %>%
    # Avoid division by zero if a sample has zero coverage
    filter(total_sample_abundance > 0) %>% 
    mutate(relative_abundance = genus_sample_abundance / total_sample_abundance)

  # Ensure month and rep are ordered factors for plotting
  df_rel_abundance$month <- factor(df_rel_abundance$month, levels = full_month_order)
  df_rel_abundance$rep <- as.factor(df_rel_abundance$rep)

  # --- Generate Faceted Bar Plot ---
  # `facet_wrap` will create a panel for each month
  comparative_barplot <- ggplot(df_rel_abundance, aes(x = rep, y = relative_abundance, fill = final_genus_category)) +
    # --- MODIFICATION IS HERE ---
    geom_bar(stat = "identity", position = "stack", width = 0.7) + # Adjust this value as needed
    # ----------------------------
  facet_wrap(~ month, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = plot_colors) +
    labs(
      title = paste("Larval Microbiome Comparison by Replicate - Colony:", current_colony),
      subtitle = "Each panel displays a different month",
      x = "Replicate ID",
      y = "Relative Abundance",
      fill = "Genus / Category"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 12),
      panel.spacing = unit(1.5, "lines"),
      legend.position = "right"
    )
  
  # --- Save the Plot ---
  # The output directory is the same, but the filename is simpler
  plot_filename <- paste0("larvae_colony_", current_colony, "_comparative_abundance.png")
  plot_filepath <- file.path(figure_output_base_dir, plot_filename)
  ggsave(plot_filepath, plot = comparative_barplot, width = 14, height = 7, units = "in", bg = "white")
}

cat("--- Finding Overlapping MAGs with High-Quality Larval Observations ---\n")

overlapping_mags_filtered <- df_processed %>%
  # 1. Filter for the specific colonies and castes
  filter(
    colony %in% c(100, 999),
    caste %in% c("Larvae", "Worker"),
    avg_cov > 5
  ) %>%
  # 2. Group by the main identifiers
  group_by(colony, month, mag_id) %>%
  # 3. Ensure the MAG is in both castes AND has at least one high-quality larval hit
  filter(
    n_distinct(caste) >= 2 &
      any(caste == "Larvae" & contig_presence_ratio > 0.5)
  ) %>%
  # 4. Summarise the results, now only including high-quality larval data
  summarise(
    genus = first(genus),
    larvae_details = {
      # This logic now only operates on the pre-filtered high-quality data
      larvae_df <- distinct(tibble(
        rep = rep[caste == "Larvae" & contig_presence_ratio > 0.5],
        ratio = contig_presence_ratio[caste == "Larvae" & contig_presence_ratio > 0.5]
      ))
      toString(paste(larvae_df$rep, round(larvae_df$ratio, 2), sep = ": "))
    },
    worker_ratios = toString(unique(round(contig_presence_ratio[caste == "Worker"], 2))),
    .groups = 'drop'
  ) %>%
  # 5. Arrange for readability
  arrange(colony, month, mag_id)

# --- Display the Results ---
print(overlapping_mags_filtered, n = 50)


# 1. Add a source label to each dataframe
shared_df <- overlapping_mags_filtered %>% mutate(group = "Shared (Larvae & Worker)")
larvae_only_df <- larvae_only_mags %>% mutate(group = "Larvae-Only")

# 2. Combine them into a single dataframe for plotting
combined_plot_data <- bind_rows(shared_df, larvae_only_df)

# 3. Create a summary of counts for the plot
genus_counts <- combined_plot_data %>%
  filter(!is.na(genus) & genus != "no match") %>%
  group_by(colony, genus, group) %>%
  summarise(mag_count = n_distinct(mag_id), .groups = 'drop')

# 4. Generate the plot
library(forcats) # For reordering factors

genus_histogram <- ggplot(genus_counts, aes(x = mag_count, y = fct_reorder(genus, mag_count, .fun = sum), fill = group)) +
  geom_col(position = "stack") +
  facet_wrap(~ colony, scales = "free") +
  labs(
    title = "Count of Unique MAGs by Genus and Colony",
    subtitle = "Comparing MAGs shared between castes vs. those found only in larvae",
    x = "Number of Unique MAGs",
    y = "Genus",
    fill = "MAG Origin"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

# --- Display the Plot ---
print(genus_histogram)




# --- 1. Prepare Data ---
test_data <- df_processed %>%
  filter(colony %in% c(100, 999), contig_presence_ratio > 0.5, avg_cov > 5)

# --- CORRECTED CODE FOR GETTING TOTAL SAMPLE COUNTS ---
total_larva_samples <- test_data %>%
  filter(caste == "Larvae") %>%
  summarise(n = n_distinct(rep, colony, month)) %>%
  pull(n)

total_worker_samples <- test_data %>%
  filter(caste == "Worker") %>%
  summarise(n = n_distinct(rep, colony, month)) %>%
  pull(n)

# Get a list of all unique MAGs to test
all_mags <- unique(test_data$mag_id)
enrichment_results <- list()

# --- 2. Loop and Test Each MAG ---
for (mag in all_mags) {
  # --- CORRECTED CODE FOR COUNTING PRESENCE ---
  present_in_larvae <- test_data %>%
    filter(mag_id == mag, caste == "Larvae") %>%
    summarise(n = n_distinct(rep, colony, month)) %>%
    pull(n)
  
  present_in_workers <- test_data %>%
    filter(mag_id == mag, caste == "Worker") %>%
    summarise(n = n_distinct(rep, colony, month)) %>%
    pull(n)
  
  # Build the 2x2 contingency table
  contingency_table <- matrix(c(
    present_in_larvae,                                  # a
    present_in_workers,                                 # b
    total_larva_samples - present_in_larvae,            # c
    total_worker_samples - present_in_workers           # d
  ), nrow = 2, byrow = TRUE)
  
  # Perform the test
  test_result <- fisher.test(contingency_table, alternative = "greater")
  
  # Store the results
  enrichment_results[[mag]] <- tibble(
    mag_id = mag,
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate
  )
}

# --- 3. View Results ---
final_results_df <- bind_rows(enrichment_results) %>%
  arrange(p_value) %>%
  left_join(df_processed %>% distinct(mag_id, genus), by = "mag_id")

# Adjust p-values for multiple comparisons
final_results_df$p_adj <- p.adjust(final_results_df$p_value, method = "BH")

print(head(final_results_df))


mags_to_investigate <- c(
  "mmag_984", "mmag_968", "mmag_934", "mmag_849", "mmag_84",
  "mmag_805", "mmag_721", "mmag_608", "mmag_605", "mmag_411"
)

# --- 2. Filter the 'df' dataframe ---
# The '%in%' operator checks if the value in the 'mag_id' column
# is present in our 'mags_to_investigate' list.
investigation_df <- df %>%
  filter(mag_id %in% mags_to_investigate)
investigation_df <- investigation_df[investigation_df$genome_wide_median_cov>0,]
investigation_df <- investigation_df[investigation_df$caste == "Larvae",]
investigation_df <- investigation_df[investigation_df$average_genome_wide_cov >= 5,]
