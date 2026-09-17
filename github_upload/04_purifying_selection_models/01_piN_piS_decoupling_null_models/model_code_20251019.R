
# --- 1. Define Paths and Constants ---
# MODIFIED: Path to the new master summary file
master_summary_file <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/final_analysis_all_colonies_aggregated/master_per_gene_summary.csv"

metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"
base_output_dir_png <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/figures/abundance_pi_relationships_png_corrected"
base_output_dir_svg <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/figures/abundance_pi_relationships_svg_corrected"

month_order <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")

# Create output directories if they don't exist
dir.create(base_output_dir_png, recursive = TRUE, showWarnings = FALSE)
dir.create(base_output_dir_svg, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load and Prepare All Necessary Data ---
cat("--- Loading and Preparing Data ---\n")

# a. Load the NEW master gene-level population genetics summary
gene_level_data_raw <- read_csv(master_summary_file)

# Check if essential columns exist
required_cols <- c("colony_id", "MAG_id", "Month", "corresponding_gene_call",
                   "pi_1D", "pi_4D", "Watterson_theta_per_bp_1D_sites",
                   "Watterson_theta_per_bp_4D_sites", "TajimaD_1D_sites",
                   "TajimaD_4D_sites", "unique_gene_callers_id") # Check unique_gene_callers_id too
missing_cols <- setdiff(required_cols, colnames(gene_level_data_raw))
if (length(missing_cols) > 0) {
  stop(paste("Error: The following required columns are missing from master_per_gene_summary.csv:", paste(missing_cols, collapse=", ")))
}

# Rename colony_id for consistency
gene_level_data <- gene_level_data_raw %>% rename(Colony = colony_id)

# b. Load ecology data and create Genus map
ecology_df <- read_csv(metadata_file) %>%
  filter(type == "mmag", Caste == "Worker") # Filter earlier
genus_map <- ecology_df %>% distinct(MAG_id, Genus)

# c. Add Genus to gene_level_data and filter out "no match"
# MODIFIED: Join Genus info here and apply filter
gene_level_data <- left_join(gene_level_data, genus_map, by = "MAG_id") %>%
  filter(Genus != "no match") # Filter out 'no match' Genus entries

# d. Aggregate gene-level data to the MAG-level
# MODIFIED: Add aggregation for pi_all_sites
mag_level_popgen <- gene_level_data %>%
  group_by(MAG_id, Colony, Month, Genus) %>% # Keep Genus here
  summarise(
    # Check for the correct pivoted column name for all sites pi
    # Common names from pivot: pi_all_sites, pi_all_sites (if values_from was just 'pi')
    # If the column name is different, adjust here.
    mean_pi_all_sites = mean(pi_all_sites, na.rm = TRUE), # ADDED THIS LINE
    mean_pi_1D = mean(pi_1D, na.rm = TRUE),
    mean_pi_4D = mean(pi_4D, na.rm = TRUE),
    mean_Watterson_theta_1D = mean(Watterson_theta_per_bp_1D_sites, na.rm = TRUE),
    mean_Watterson_theta_4D = mean(Watterson_theta_per_bp_4D_sites, na.rm = TRUE),
    mean_TajimaD_1D = mean(TajimaD_1D_sites, na.rm = TRUE),
    mean_TajimaD_4D = mean(TajimaD_4D_sites, na.rm = TRUE),
    .groups = 'drop'
  )

# e. Add numeric Time variable
mag_level_popgen$TimePoint <- match(mag_level_popgen$Month, month_order) # Renamed to avoid conflict

# f. Load and calculate relative abundance data
abundance_data <- ecology_df %>%
  group_by(Colony, Month, MAG_id, Genus, Season) %>%
  summarise(mean_coverage = mean(mean, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Colony, Month) %>%
  mutate(rel_abund_mag = mean_coverage / sum(mean_coverage)) %>%
  ungroup() %>%
  # Filter out "no match" here too for consistency
  filter(Genus != "no match")

# g. Merge datasets
full_data <- inner_join(mag_level_popgen,
                        abundance_data,
                        by = c("Colony", "Month", "MAG_id", "Genus")) %>%
  filter(!is.na(Season)) # Ensure Season is not NA

# h. Add cumulative time
add_cumulative_time_from_month_nucdiv <- function(data) {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34,
                        "September" = 28, "October" = 22, "November" = 42,
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  data <- data %>% mutate(CumulativeTime = month_to_cumulative_time_map[Month])
  return(data)
}
full_data <- add_cumulative_time_from_month_nucdiv(full_data)

# i. Filter groups with too few timepoints
full_data <- full_data %>%
  group_by(Colony, MAG_id) %>%
  filter(n() >= 3) %>%
  ungroup()

# j. Calculate Alpha Diversity (Shannon and Richness)
# Use abundance_data as it's already grouped correctly before merging
abundance_wide <- abundance_data %>%
  select(Colony, Month, MAG_id, rel_abund_mag) %>%
  tidyr::pivot_wider(
    names_from = MAG_id,
    values_from = rel_abund_mag,
    values_fill = 0 # Fill missing MAGs in a sample with 0 abundance
  )

shannon_values <- diversity(abundance_wide %>% select(-Colony, -Month), "shannon")
richness_values <- specnumber(abundance_wide %>% select(-Colony, -Month)) # Use vegan::specnumber for richness

alpha_div_df <- abundance_wide %>%
  select(Colony, Month) %>%
  mutate(shannon_diversity = shannon_values,
         richness = richness_values)

# Join alpha diversity back to the main 'full_data' table
full_data <- left_join(full_data, alpha_div_df, by = c("Colony", "Month"))

# k. Convert Season to factor and add genome size
full_data$Season <- as.factor(full_data$Season)
genome_size_map <- ecology_df %>% distinct(MAG_id, genome_size)
full_data <- full_data %>% left_join(genome_size_map, by = "MAG_id")

# l. Final filtering for complete cases in key variables
full_data <- full_data %>%
  filter(!is.na(mean_pi_1D), !is.na(mean_pi_4D), !is.na(shannon_diversity),
         !is.na(richness), !is.na(mean_coverage), !is.na(Season), !is.na(genome_size))


cat("--- Data Loading and Preparation Complete ---\n")
cat("Final data dimensions:", dim(full_data), "\n")
cat("First few rows of the final data for analysis:\n\n")
print(head(full_data))

# --- 3. Analysis 1: Linear Regression (Stability vs. Shift) ---
cat("\n--- 1. Linear Regression (Mean Diversity vs. Time) ---\n")

regression_results_1D <- full_data %>%
  group_by(Colony, MAG_id, Genus) %>%
  # Use tryCatch in case lm fails (e.g., only one unique Time value)
  do(tryCatch({ broom::tidy(lm(mean_pi_1D ~ CumulativeTime, data = .)) }, error = function(e) NULL)) %>%
  filter(!is.null(term), term == "CumulativeTime") %>% # Filter for successful models and the Time term
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p_adj) %>%
  mutate(site_type = "1D")

regression_results_4D <- full_data %>%
  group_by(Colony, MAG_id, Genus) %>%
  do(tryCatch({ broom::tidy(lm(mean_pi_4D ~ CumulativeTime, data = .)) }, error = function(e) NULL)) %>%
  filter(!is.null(term), term == "CumulativeTime") %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p_adj) %>%
  mutate(site_type = "4D")

regression_results <- bind_rows(regression_results_1D, regression_results_4D)
print(head(regression_results))

# --- 4. Analysis 2: Change Point Analysis ---
cat("\n--- 2. Change Point Analysis (Abrupt Shifts in Mean Diversity) ---\n")
changepoint_results_1D <- full_data %>%
  arrange(Colony, MAG_id, CumulativeTime) %>%
  group_by(Colony, MAG_id, Genus) %>%
  summarise(
    cpt_object = list(tryCatch({
      cpt.mean(mean_pi_1D, method = "PELT", minseglen = 2)
    }, error = function(e) { NA })),
    .groups = 'drop'
  ) %>% filter(!is.na(cpt_object)) %>%
  mutate(num_changepoints = purrr::map_int(cpt_object, ~(nseg(.x) - 1))) %>%
  select(-cpt_object) %>% filter(num_changepoints > 0) %>%
  arrange(desc(num_changepoints)) %>% mutate(site_type = "1D")

changepoint_results_4D <- full_data %>%
  arrange(Colony, MAG_id, CumulativeTime) %>%
  group_by(Colony, MAG_id, Genus) %>%
  summarise(
    cpt_object = list(tryCatch({
      cpt.mean(mean_pi_4D, method = "PELT", minseglen = 2)
    }, error = function(e) { NA })),
    .groups = 'drop'
  ) %>% filter(!is.na(cpt_object)) %>%
  mutate(num_changepoints = purrr::map_int(cpt_object, ~(nseg(.x) - 1))) %>%
  select(-cpt_object) %>% filter(num_changepoints > 0) %>%
  arrange(desc(num_changepoints)) %>% mutate(site_type = "4D")

changepoint_results <- bind_rows(changepoint_results_1D, changepoint_results_4D)
print(head(changepoint_results))

# --- 5. Analysis 3: ANCOVA ---
cat("\n--- 3. ANCOVA (Comparing Mean Diversity Trends within Genera) ---\n")
ancova_results_1D <- full_data %>%
  group_by(Colony, Genus) %>%
  filter(n_distinct(MAG_id) > 1) %>%
  do(tryCatch({ broom::tidy(aov(mean_pi_1D ~ CumulativeTime * MAG_id, data = .)) }, error = function(e) NULL)) %>%
  filter(!is.null(term), grepl("CumulativeTime:MAG_id", term)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% arrange(p_adj) %>% mutate(site_type = "1D")

ancova_results_4D <- full_data %>%
  group_by(Colony, Genus) %>%
  filter(n_distinct(MAG_id) > 1) %>%
  do(tryCatch({ broom::tidy(aov(mean_pi_4D ~ CumulativeTime * MAG_id, data = .)) }, error = function(e) NULL)) %>%
  filter(!is.null(term), grepl("CumulativeTime:MAG_id", term)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>% arrange(p_adj) %>% mutate(site_type = "4D")

ancova_results <- bind_rows(ancova_results_1D, ancova_results_4D)
print(head(ancova_results))

# --- 6. Analysis 4: Effect of Season (LMM) ---
cat("\n--- 4. Testing for Significant Seasonal Effects (LMM) ---\n")
lmm_pi_1D <- lmer(log(mean_pi_1D + 1e-9) ~ Season + (1 | Colony) + (1 | MAG_id), data = full_data)
lmm_pi_4D <- lmer(log(mean_pi_4D + 1e-9) ~ Season + (1 | Colony) + (1 | MAG_id), data = full_data)
lmm_abund <- lmer(log(mean_coverage + 1) ~ Season + (1 | Colony) + (1 | MAG_id), data = full_data)
cat("\n--- LMM Summary: Effect of Season on pi_1D ---\n"); print(summary(lmm_pi_1D))
cat("\n--- LMM Summary: Effect of Season on pi_4D ---\n"); print(summary(lmm_pi_4D))
cat("\n--- LMM Summary: Effect of Season on Absolute Abundance ---\n"); print(summary(lmm_abund))

# --- 7. Analysis 5: Correlation of Diversity and Relative Abundance (Raw) ---
cat("\n--- 5. Pearson Correlation (Raw Diversity vs. Relative Abundance) ---\n")
correlation_results_raw <- full_data %>%
  group_by(Colony, MAG_id, Genus) %>%
  filter(n() > 4) %>%
  summarise(
    pearson_cor_1D = cor(mean_pi_1D, rel_abund_mag, method = "pearson", use = "pairwise.complete.obs"),
    pearson_cor_4D = cor(mean_pi_4D, rel_abund_mag, method = "pearson", use = "pairwise.complete.obs"),
    n_months = n(), .groups = 'drop' )
print(head(correlation_results_raw))

# --- 8. Analysis 6: Detrended Correlation (Focus on 4D sites, Dynamic Filter) ---
cat("\n--- 6. Detrended Correlation (4D Diversity vs. Relative Abundance) ---\n")

# --- NEW: Calculate threshold for each colony ---
# Count unique months sampled per colony in the full dataset
colony_month_counts <- full_data %>%
  distinct(Colony, Month) %>%
  group_by(Colony) %>%
  summarise(total_months_in_colony = n(), .groups = 'drop')
full_data
# Join this count back to the main data
full_data_with_counts <- left_join(full_data, colony_month_counts, by = "Colony")

full_data
# --- Step A: Detrending (Only 4D sites, Dynamic Filter) ---
detrended_data <- full_data_with_counts %>%
  group_by(Colony, MAG_id) %>%
  # --- MODIFIED FILTER ---
  # Keep groups present in at least 25% of the colony's sampled months
  # Ensure total_months_in_colony is available and > 0
  filter(!is.na(total_months_in_colony) & total_months_in_colony > 2) %>%
  do({
    group_data <- .
    # Models for detrending (only 4D and abundance)
    pi_4D_model <- tryCatch(lm(mean_pi_4D ~ CumulativeTime, data = group_data, na.action = na.exclude), error = function(e) NULL)
    abund_model <- tryCatch(lm(rel_abund_mag ~ CumulativeTime, data = group_data, na.action = na.exclude), error = function(e) NULL)

    # Add residuals
    group_data$pi_4D_residuals <- if (!is.null(pi_4D_model)) resid(pi_4D_model) else NA_real_
    group_data$abund_residuals <- if (!is.null(abund_model)) resid(abund_model) else NA_real_

    group_data
  }) %>%
  ungroup()

# --- Step B & C: Calculate Correlation & p-values (Only 4D sites) ---
correlation_p_values_4D <- detrended_data %>%
  # Filter out groups where residuals couldn't be calculated
  filter(!is.na(pi_4D_residuals) & !is.na(abund_residuals)) %>%
  group_by(Colony, MAG_id, Genus) %>% # Group by Genus too if available
  # Apply cor.test independently to each group using do()
  do({
    group_data <- .
    # Run cor.test for 4D, wrapped in tryCatch and tidy
    test_4D <- tryCatch({
      broom::tidy(cor.test(~ pi_4D_residuals + abund_residuals, data = group_data, use = "pairwise.complete.obs"))
    }, error = function(e) {
      tibble() # Return empty tibble if test fails
    })
    test_4D # Return the result (or empty tibble)
  }) %>%
  ungroup() %>% # Ungroup after do()
  # Rename columns and calculate adjusted p-values
  rename(pearson_cor_detrended = estimate, p.value = p.value) %>%
  filter(!is.na(pearson_cor_detrended)) %>% # Ensure test ran successfully
  # P-value adjustment (no need to group by site_type anymore)
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  # Select and arrange final columns
  select(Colony, MAG_id, Genus, pearson_cor_detrended, p.value, p_adj, statistic, parameter, conf.low, conf.high, method, alternative) %>%
  arrange(p_adj)

cat("\n--- Detrended Correlation (4D only) with p-values ---\n")
print(correlation_p_values_4D, n = 20) # Print more rows

significant_correlations_4D <- correlation_p_values_4D %>%
  filter(p_adj < 0.05)

cat("\n--- Only Significant Detrended Correlations (4D only, p_adj < 0.05) ---\n")
print(significant_correlations_4D)







# --- 8b. Analysis 6b: Detrended Correlation (ALL SITES) ---
cat("\n--- 6b. Detrended Correlation (ALL SITES Diversity vs. Relative Abundance) ---\n")

# --- Step A: Detrending (All sites, Dynamic Filter) ---
# Check if the 'pi_all_sites' column exists from the pivot
# The pivot creates columns like 'pi_all_sites', 'pi_1D_sites', 'pi_4D_sites'
# Let's verify and potentially rename for clarity if needed.
# If your pivot resulted in just 'pi' for all sites, adjust accordingly.
pi_all_sites_col <- "mean_pi_all_sites" # Default assumption from pivot
if (!pi_all_sites_col %in% colnames(full_data_with_counts)) {
  # If the pivot used a different naming, try common alternatives
  if ("pi" %in% colnames(full_data_with_counts)) {
    pi_all_sites_col <- "pi"
    warning("Using column 'pi' as pi_all_sites. Verify this is correct.")
  } else if ("mean_pi" %in% colnames(full_data_with_counts)) {
    # Check if the old mean_pi aggregation still exists inadvertently
    pi_all_sites_col <- "mean_pi"
    warning("Using column 'mean_pi' as pi_all_sites. Verify this is correct.")
  } else {
    stop("Could not find a column representing pi calculated from all sites (expected 'pi_all_sites' or similar). Please check the column names in 'full_data_with_counts'.")
  }
}

detrended_data_allsites <- full_data_with_counts %>%
  group_by(Colony, MAG_id) %>%
  # Use the same dynamic filter
  filter(!is.na(total_months_in_colony) & total_months_in_colony > 2) %>%
  do({
    group_data <- .
    # Models for detrending (pi_all_sites and abundance)
    # Use .data[[pi_all_sites_col]] to dynamically use the correct column name
    pi_all_model <- tryCatch(lm(as.formula(paste(pi_all_sites_col, "~ CumulativeTime")), data = group_data, na.action = na.exclude), error = function(e) NULL)
    abund_model <- tryCatch(lm(rel_abund_mag ~ CumulativeTime, data = group_data, na.action = na.exclude), error = function(e) NULL)
    
    # Add residuals
    group_data$pi_all_sites_residuals <- if (!is.null(pi_all_model)) resid(pi_all_model) else NA_real_
    group_data$abund_residuals <- if (!is.null(abund_model)) resid(abund_model) else NA_real_
    
    group_data
  }) %>%
  ungroup()

# --- Step B & C: Calculate Correlation & p-values (All sites) ---
correlation_p_values_allsites <- detrended_data_allsites %>%
  # Filter out groups where residuals couldn't be calculated
  filter(!is.na(pi_all_sites_residuals) & !is.na(abund_residuals)) %>%
  group_by(Colony, MAG_id, Genus) %>% # Group by Genus too if available
  # Apply cor.test independently to each group using do()
  do({
    group_data <- .
    # Run cor.test for all sites, wrapped in tryCatch and tidy
    test_all <- tryCatch({
      broom::tidy(cor.test(~ pi_all_sites_residuals + abund_residuals, data = group_data, use = "pairwise.complete.obs"))
    }, error = function(e) {
      tibble() # Return empty tibble if test fails
    })
    test_all # Return the result (or empty tibble)
  }) %>%
  ungroup() %>% # Ungroup after do()
  # Rename columns and calculate adjusted p-values
  rename(pearson_cor_detrended = estimate, p.value = p.value) %>%
  filter(!is.na(pearson_cor_detrended)) %>% # Ensure test ran successfully
  # P-value adjustment
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  # Select and arrange final columns
  select(Colony, MAG_id, Genus, pearson_cor_detrended, p.value, p_adj, statistic, parameter, conf.low, conf.high, method, alternative) %>%
  arrange(p_adj)

cat("\n--- Detrended Correlation (ALL SITES) with p-values ---\n")
print(correlation_p_values_allsites, n = 20) # Print more rows

significant_correlations_allsites <- correlation_p_values_allsites %>%
  filter(p_adj < 0.1)

cat("\n--- Only Significant Detrended Correlations (ALL SITES, p_adj < 0.05) ---\n")
print(significant_correlations_allsites)

# You can now compare 'correlation_p_values_allsites' with 'correlation_p_values_4D'

hist(correlation_p_values_allsites$pearson_cor_detrended)






# --- 1. Define Thresholds ---
bh_fdr_threshold <- 0.05 # Adjusted p-value threshold for significance (e.g., 0.05 or 0.1)
pos_threshold <- 0.3
neg_threshold <- -0.3

# --- 2. Create Categories using p_adj_BH ---
plot_summary_allsites_bh <- correlation_p_values_allsites %>%
  # Handle potential NA p_adj values (replace with 1 for categorization)
  mutate(p_adj_clean = ifelse(is.na(p_adj), 1, p_adj)) %>%
  mutate(
    correlation_type = case_when(
      # Significant Positive (Stronger than threshold)
      p_adj_clean < bh_fdr_threshold & pearson_cor_detrended >= pos_threshold  ~ "Significant Positive",
      # Significant Negative (Stronger than threshold)
      p_adj_clean < bh_fdr_threshold & pearson_cor_detrended <= neg_threshold  ~ "Significant Negative",
      # Significant Weak/None (Between thresholds, but significant)
      p_adj_clean < bh_fdr_threshold & pearson_cor_detrended > neg_threshold & pearson_cor_detrended < pos_threshold ~ "Significant Weak/None",
      # Not Significant Positive (Stronger than threshold)
      p_adj_clean >= bh_fdr_threshold & pearson_cor_detrended >= pos_threshold ~ "Not Significant Positive",
      # Not Significant Negative (Stronger than threshold)
      p_adj_clean >= bh_fdr_threshold & pearson_cor_detrended <= neg_threshold ~ "Not Significant Negative",
      # Not Significant Weak/None (Between thresholds)
      p_adj_clean >= bh_fdr_threshold & pearson_cor_detrended > neg_threshold & pearson_cor_detrended < pos_threshold ~ "Not Significant Weak/None",
      # Catch unexpected NAs in pearson_cor_detrended if any
      TRUE                                                  ~ "Check Data"
    )
  ) %>%
  # Count the number of MAGs in each category for each colony
  group_by(Colony, correlation_type) %>%
  summarise(n = n(), .groups = 'drop') %>%
  # Calculate the proportion of each category within each colony
  group_by(Colony) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() # Ungroup for plotting

# --- 3. Define Plotting Ele# --- 3. Define Plotting Ele# --- 3. Define Plotting Elements ---
# Define the color palette (same as before)
color_palette <- c(
  "Significant Positive"      = "#cc241d", # Strong Red
  "Significant Negative"      = "#458588", # Strong Blue
  "Significant Weak/None"     = "#fabd2f", # Yellow/Orange
  "Not Significant Positive"  = "#ebdbb2", # Lighter Orange/Red
  "Not Significant Negative"  = "#665c54", # Lighter Blue/Grey
  "Not Significant Weak/None" = "#32302f", # Dark Grey
  "Check Data"                = "purple"
)

# Define the logical order for the legend/stacking
category_levels <- c(
  "Significant Positive", "Not Significant Positive",
  "Significant Weak/None", "Not Significant Weak/None",
  "Not Significant Negative", "Significant Negative"
)

# Ensure the categories are ordered logically, handling missing categories
plot_summary_allsites_bh$correlation_type <- factor(
  plot_summary_allsites_bh$correlation_type,
  levels = intersect(category_levels, unique(plot_summary_allsites_bh$correlation_type))
)

# --- 4. Build the ggplot ---
plot_allsites_summary_bh <- ggplot(plot_summary_allsites_bh, aes(x = as.factor(Colony), y = proportion, fill = correlation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette, drop = FALSE) + # drop=FALSE keeps all colors
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportions of Detrended Correlations (All Sites Pi vs Abundance)",
    # MODIFIED: Updated subtitle to mention p_adj
    subtitle = paste("Significance: BH p_adj <", bh_fdr_threshold, "; Strength: |r| >=", pos_threshold),
    x = "Colony",
    y = "Proportion of MAGs",
    fill = "Correlation Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 5. Display and Save the Plot ---
print(plot_allsites_summary_bh)

# Save the plot (using BH in filename)
ggsave(file.path(base_output_dir_png, "correlation_summary_allsites_by_colony_BH_thresh.png"), plot = plot_allsites_summary_bh, width = 10, height = 6, dpi = 150, bg = "white")
ggsave(file.path(base_output_dir_svg, "correlation_summary_allsites_by_colony_BH_thresh.svg"), plot = plot_allsites_summary_bh, width = 10, height = 6)

cat("\n--- All sites correlation summary plot (using BH p_adj) generated. ---\n")






# --- XI. Plotting Correlation Summaries (All Sites, Magnitude Cutoff Only) ---
cat("\n--- XI. Generating Correlation Summary Plots (Magnitude Cutoff Only) ---\n")

# Ensure 'correlation_results_allsites' dataframe exists

# --- 1. Define Thresholds ---
# Using the correlation magnitude cutoffs only, ignoring p_adj/q_value for categorization
pos_threshold <- 0.3
neg_threshold <- -0.3

# --- 2. Create Simplified Categories for All Sites Data ---
correlation_categories_allsites_magnitude <- correlation_p_values_allsites %>%
  mutate(
    correlation_type_magnitude = case_when(
      pearson_cor_detrended >= pos_threshold  ~ "Positive Correlation",
      pearson_cor_detrended <= neg_threshold  ~ "Negative Correlation",
      TRUE                                    ~ "No/Weak Correlation" # Catches everything between -0.3 and 0.3
    )
  )

df = correlation_categories_allsites_magnitude[correlation_categories_allsites_magnitude$correlation_type_magnitude=="Positive Correlation",]
# --- 3. Define Plotting Elements (Colors & Levels) ---
# Define a custom color palette for the 3 categories
color_palette_magnitude <- c(
  "Positive Correlation"  = "#d79921", # Yellow for positive
  "Negative Correlation"  = "#458588", # Blue for negative
  "No/Weak Correlation"   = "#32302f"  # Dark Grey for no/weak
)

# Define the logical order for the legend/stacking
category_levels_magnitude <- c(
  "Positive Correlation",
  "No/Weak Correlation",
  "Negative Correlation"
)

# --- Plot 1: Grouped by Colony ---

# Prepare data for Colony plot
plot_summary_colony <- correlation_categories_allsites_magnitude %>%
  group_by(Colony, correlation_type_magnitude) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Colony) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Ensure categories are ordered
plot_summary_colony$correlation_type_magnitude <- factor(
  plot_summary_colony$correlation_type_magnitude,
  levels = intersect(category_levels_magnitude, unique(plot_summary_colony$correlation_type_magnitude))
)

# Build the ggplot for Colony
plot_allsites_colony_magnitude <- ggplot(plot_summary_colony, aes(x = as.factor(Colony), y = proportion, fill = correlation_type_magnitude)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_magnitude, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Correlations (All Sites Pi vs Abundance) by Colony",
    subtitle = paste("Correlation Strength Threshold: |r| >=", pos_threshold, " (Ignoring p-value)"),
    x = "Colony",
    y = "Proportion of MAGs",
    fill = "Correlation Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_allsites_colony_magnitude)
ggsave(file.path(base_output_dir_png, "correlation_summary_allsites_by_colony_magnitude.png"), plot = plot_allsites_colony_magnitude, width = 10, height = 6, dpi = 150, bg = "white")
ggsave(file.path(base_output_dir_svg, "correlation_summary_allsites_by_colony_magnitude.svg"), plot = plot_allsites_colony_magnitude, width = 10, height = 6)


# --- Plot 2: Grouped by Genus ---

# Prepare data for Genus plot
plot_summary_genus <- correlation_categories_allsites_magnitude %>%
  # Filter out Genus == NA if you don't want them in the plot
  filter(!is.na(Genus) & Genus != "") %>% # Also remove empty strings
  group_by(Genus, correlation_type_magnitude) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Genus) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Filter for genera with at least N MAGs (e.g., 5 MAGs) to avoid too many small bars
# This is optional, but often makes genus plots more readable
min_mags_per_genus <- 5
genera_to_plot <- plot_summary_genus %>%
  group_by(Genus) %>%
  summarise(total_mags = sum(n), .groups = 'drop') %>%
  filter(total_mags >= min_mags_per_genus) %>%
  pull(Genus)

plot_summary_genus_filtered <- plot_summary_genus %>%
  filter(Genus %in% genera_to_plot)

# Ensure categories are ordered
plot_summary_genus_filtered$correlation_type_magnitude <- factor(
  plot_summary_genus_filtered$correlation_type_magnitude,
  levels = intersect(category_levels_magnitude, unique(plot_summary_genus_filtered$correlation_type_magnitude))
)

# Build the ggplot for Genus
plot_allsites_genus_magnitude <- ggplot(plot_summary_genus_filtered, aes(x = as.factor(Genus), y = proportion, fill = correlation_type_magnitude)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_magnitude, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Correlations (All Sites Pi vs Abundance) by Genus",
    subtitle = paste("Correlation Strength Threshold: |r| >=", pos_threshold, " (Ignoring p-value). Only genera with >=", min_mags_per_genus, "MAGs shown."),
    x = "Genus",
    y = "Proportion of MAGs",
    fill = "Correlation Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_allsites_genus_magnitude)
ggsave(file.path(base_output_dir_png, "correlation_summary_allsites_by_genus_magnitude.png"), plot = plot_allsites_genus_magnitude, width = 12, height = 7, dpi = 150, bg = "white") # Wider plot for genera
ggsave(file.path(base_output_dir_svg, "correlation_summary_allsites_by_genus_magnitude.svg"), plot = plot_allsites_genus_magnitude, width = 12, height = 7)

cat("\n--- All sites correlation summary plots (magnitude only) generated. ---\n")






# --- Save Plot 1 (Colony) ---

# Save PNG (no change needed)
ggsave(file.path(base_output_dir_png, "correlation_summary_allsites_by_colony_magnitude.png"),
       plot = plot_allsites_colony_magnitude, width = 10, height = 6, dpi = 150, bg = "white")

# Save SVG using base R's SVG device via ggsave's device argument
ggsave(file.path(base_output_dir_svg, "correlation_summary_allsites_by_colony_magnitude.svg"),
       plot = plot_allsites_colony_magnitude, width = 10, height = 6,
       device = grDevices::svg) # Specify the device explicitly

cat("\n--- Saved Colony summary plots (PNG & SVG without svglite). ---\n")


# --- Save Plot 2 (Genus) ---

# Save PNG (no change needed)
ggsave(file.path(base_output_dir_png, "correlation_summary_allsites_by_genus_magnitude.png"),
       plot = plot_allsites_genus_magnitude, width = 12, height = 7, dpi = 150, bg = "white")

# Save SVG using base R's SVG device via ggsave's device argument
ggsave(file.path(base_output_dir_svg, "correlation_summary_allsites_by_genus_magnitude.svg"),
       plot = plot_allsites_genus_magnitude, width = 12, height = 7,
       device = grDevices::svg) # Specify the device explicitly

cat("\n--- Saved Genus summary plots (PNG & SVG without svglite). ---\n")





























# --- 9. Analysis 7: Test for Effect of Season (per MAG) ---
cat("\n--- 7. Testing for Seasonal Effects per MAG ---\n")
test_season_effect <- function(data, metric_name) {
  if (nrow(data) < 4 || n_distinct(data$Season) < 2) { return(NA_real_) }
  formula_str <- paste(metric_name, "~ Season") # Ensure metric_name is correct
  model <- tryCatch(lm(as.formula(formula_str), data = data), error = function(e) NULL)
  if (is.null(model)) return(NA_real_)
  summary_stats <- broom::glance(model)
  return(summary_stats$p.value)
}

mag_seasonal_summary <- full_data %>%
  group_by(MAG_id) %>% nest() %>%
  mutate(
    p_value_pi_1D = map_dbl(data, ~test_season_effect(.x, "mean_pi_1D")), # Use correct column name
    p_value_pi_4D = map_dbl(data, ~test_season_effect(.x, "mean_pi_4D")), # Use correct column name
    p_value_abund = map_dbl(data, ~test_season_effect(.x, "rel_abund_mag"))
  ) %>% select(-data) %>% ungroup() %>%
  mutate(
    p_adj_pi_1D = p.adjust(p_value_pi_1D, method = "BH"),
    p_adj_pi_4D = p.adjust(p_value_pi_4D, method = "BH"),
    p_adj_abund = p.adjust(p_value_abund, method = "BH") )

seasonal_pi_mags_1D <- mag_seasonal_summary %>% filter(p_adj_pi_1D < 0.1) %>% arrange(p_adj_pi_1D)
seasonal_pi_mags_4D <- mag_seasonal_summary %>% filter(p_adj_pi_4D < 0.1) %>% arrange(p_adj_pi_4D)
seasonal_abund_mags <- mag_seasonal_summary %>% filter(p_adj_abund < 0.1) %>% arrange(p_adj_abund)
cat("--- MAGs with Sig Seasonal Pi (1D) ---\n"); print(seasonal_pi_mags_1D)
cat("\n--- MAGs with Sig Seasonal Pi (4D) ---\n"); print(seasonal_pi_mags_4D)
cat("\n--- MAGs with Sig Seasonal Abundance ---\n"); print(seasonal_abund_mags)


# --- 10. Analysis 8: GLMM Temporal (Alpha Div, Season, Abundance) ---
cat("\n--- 8. GLMM: Temporal Fluctuation (Diversity vs. Alpha Div, Season, Abundance) ---\n")

# Scale predictors for GLMM
full_data_scaled <- full_data %>%
  mutate(across(c(shannon_diversity, richness, mean_coverage, genome_size), scale))

# Run separate models for pi_1D and pi_4D using Gamma family
glmm_temporal_1D <- glmer(
  mean_pi_1D ~ shannon_diversity + richness + Season + mean_coverage + genome_size + (1 | Colony) + (1 | MAG_id),
  data = full_data_scaled, family = Gamma(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)) # Increase iterations
)
glmm_temporal_4D <- glmer(
  mean_pi_4D ~ shannon_diversity + richness + Season + mean_coverage + genome_size + (1 | Colony) + (1 | MAG_id),
  data = full_data_scaled, family = Gamma(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
cat("\n--- GLMM Summary: Temporal Fluctuation (pi_1D) ---\n"); print(summary(glmm_temporal_1D))
cat("\n--- GLMM Summary: Temporal Fluctuation (pi_4D) ---\n"); print(summary(glmm_temporal_4D))

# --- 11. Analysis 9: Beta GLMM ---
cat("\n--- 9. Beta GLMM ---\n")

# a. Transform pi for Beta distribution (requires values > 0 and < 1)
n_obs <- nrow(full_data_scaled)
full_data_beta <- full_data_scaled %>%
  # Apply transformation separately for 1D and 4D
  mutate(
    pi_beta_1D = (mean_pi_1D * (n() - 1) + 0.5) / n(), # Use n() within group if needed, or n_obs if global
    pi_beta_4D = (mean_pi_4D * (n() - 1) + 0.5) / n()
  ) %>%
  # Filter out any rows where transformation might fail (if pi was exactly 0 or 1, though unlikely)
  filter(pi_beta_1D > 0, pi_beta_1D < 1, pi_beta_4D > 0, pi_beta_4D < 1)

# b. Fit Beta GLMM models
glmm_beta_1D <- glmmTMB(
  pi_beta_1D ~ shannon_diversity + richness + Season + mean_coverage + genome_size + (1 | Colony) + (1 | MAG_id),
  data = full_data_beta, family = beta_family(link = "logit")
)
glmm_beta_4D <- glmmTMB(
  pi_beta_4D ~ shannon_diversity + richness + Season + mean_coverage + genome_size + (1 | Colony) + (1 | MAG_id),
  data = full_data_beta, family = beta_family(link = "logit")
)
cat("\n--- Beta GLMM Summary (pi_1D) ---\n"); print(summary(glmm_beta_1D))
cat("\n--- Beta GLMM Summary (pi_4D) ---\n"); print(summary(glmm_beta_4D))

# --- 12. Model Diagnostics (Example for Beta GLMM 1D) ---
cat("\n--- 10. Diagnostics for Beta GLMM (pi_1D) ---\n")

# a. R-squared
cat("\n--- R-squared ---\n"); print(MuMIn::r.squaredGLMM(glmm_beta_1D))

# b. DHARMa Diagnostics
cat("\n--- DHARMa Residual Plots ---\n")
simulationOutput_beta_1D <- simulateResiduals(fittedModel = glmm_beta_1D, plot = TRUE)
cat("\n--- DHARMa Dispersion Test ---\n"); testDispersion(simulationOutput_beta_1D)

# c. Null Model Comparison
cat("\n--- Full vs. Null Model Comparison ---\n")
glmm_beta_null_1D <- glmmTMB(
  pi_beta_1D ~ 1 + (1 | Colony) + (1 | MAG_id),
  data = full_data_beta, family = beta_family(link = "logit")
)
print(anova(glmm_beta_1D, glmm_beta_null_1D))



# b. Fit Beta GLMM models

)
glmm_beta_4D <- glmmTMB(
  pi_beta_4D ~ shannon_diversity + richness + Season + mean_coverage + genome_size + (1 | Colony) + (1 | MAG_id),
  data = full_data_beta, family = beta_family(link = "logit")
)
cat("\n--- Beta GLMM Summary (pi_4D) ---\n"); print(summary(glmm_beta_4D))

# --- 12. Model Diagnostics (Example for Beta GLMM 1D) ---
cat("\n--- 10. Diagnostics for Beta GLMM (pi_4D) ---\n")

# a. R-squared
cat("\n--- R-squared ---\n"); print(MuMIn::r.squaredGLMM(glmm_beta_4D))

# b. DHARMa Diagnostics
cat("\n--- DHARMa Residual Plots ---\n")
simulationOutput_beta_4D <- simulateResiduals(fittedModel = glmm_beta_4D, plot = TRUE)
cat("\n--- DHARMa Dispersion Test ---\n"); testDispersion(simulationOutput_beta_4D)

# c. Null Model Comparison
cat("\n--- Full vs. Null Model Comparison ---\n")
glmm_beta_null_4D <- glmmTMB(
  pi_beta_4D ~ 1 + (1 | Colony) + (1 | MAG_id),
  data = full_data_beta, family = beta_family(link = "logit")
)
print(anova(glmm_beta_4D, glmm_beta_null_4D))

# --- 13. Plotting Predicted Effects (Example for Beta GLMM 1D Shannon Effect) ---
cat("\n--- 11. Plotting Predicted Effects (Beta GLMM pi_1D vs Shannon) ---\n")

predicted_effects_beta_1D <- ggpredict(glmm_beta_1D, terms = "shannon_diversity [all]")
predicted_effects_beta_4D <- ggpredict(glmm_beta_4D, terms = "shannon_diversity [all]")

# Get original mean/sd for back-transforming axis labels
shannon_mean_orig <- mean(full_data$shannon_diversity, na.rm = TRUE)
shannon_sd_orig <- sd(full_data$shannon_diversity, na.rm = TRUE)

plot_beta_shannon <- ggplot(predicted_effects_beta_1D, aes(x = x, y = predicted)) +
  geom_point(data = full_data_beta, aes(x = shannon_diversity, y = pi_beta_1D), color = "#a89984", alpha = 0.3, shape = 16) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#83a598", alpha = 0.3) +
  geom_line(color = "#282828", linewidth = 1.2) +
  labs(title = "Effect of Community Diversity on Genetic Diversity (pi_1D)",
       subtitle = "Prediction from Beta GLMM with 95% Confidence Interval",
       x = "Shannon Diversity of Bee Gut Community",
       y = expression(paste("Predicted Nucleotide Diversity (", pi ["1D"], ")"))) +
  scale_x_continuous(labels = ~round((.x * shannon_sd_orig) + shannon_mean_orig, 1)) + # Back-transform labels
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), axis.title = element_text())

print(plot_beta_shannon)
ggsave(file.path(base_output_dir_png, "beta_glmm_shannon_effect_1D.png"), plot = plot_beta_shannon, width = 7, height = 5)
ggsave(file.path(base_output_dir_svg, "beta_glmm_shannon_effect_1D.svg"), plot = plot_beta_shannon, width = 7, height = 5)

# --- 14. Plotting Individual Gene/MAG Dynamics (Example for one MAG) ---
# Select a MAG_id to plot
mag_to_plot <- "mmag_691" # Or choose another MAG

cat(paste("\n--- 12. Generating Dynamics Plot for", mag_to_plot, "---\n"))

# Prepare data for the specific MAG
plot_mag_data <- gene_level_data %>%
  filter(MAG_id == mag_to_plot) %>%
  mutate(Month = factor(Month, levels = month_order)) %>%
  filter(!is.na(Month)) %>%
  mutate(CumulativeTime = match(Month, month_order)) # Simple time index for plotting

plot_avg_data <- full_data %>%
  filter(MAG_id == mag_to_plot) %>%
  mutate(Month = factor(Month, levels = month_order)) %>%
  filter(!is.na(Month))

# Plotting function (adapted)
plot_mag_dynamics <- function(mag_data, avg_data) {
  p_pi_1D <- ggplot() +
    geom_line(data = mag_data, aes(x = CumulativeTime, y = pi_1D, group = corresponding_gene_call), color = "#83a598", alpha = 0.15) +
    geom_line(data = avg_data, aes(x = CumulativeTime, y = mean_pi_1D), color = "#282828", linewidth = 1.2) +
    facet_wrap(~ Colony, nrow = 1) +
    labs(x = "Month", y = expression(pi["1D"])) +
    scale_x_continuous(breaks = 1:length(month_order), labels = month_order) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_pi_4D <- ggplot() +
    geom_line(data = mag_data, aes(x = CumulativeTime, y = pi_4D, group = corresponding_gene_call), color = "#fb4934", alpha = 0.15) +
    geom_line(data = avg_data, aes(x = CumulativeTime, y = mean_pi_4D), color = "#282828", linewidth = 1.2) +
    facet_wrap(~ Colony, nrow = 1) +
    labs(x = "Month", y = expression(pi["4D"])) +
    scale_x_continuous(breaks = 1:length(month_order), labels = month_order) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots
  (p_pi_1D / p_pi_4D) + plot_annotation(title = paste("Temporal Dynamics for", unique(avg_data$MAG_id)))
}

if(nrow(plot_mag_data) > 0 && nrow(plot_avg_data) > 0) {
  dynamics_plot <- plot_mag_dynamics(plot_mag_data, plot_avg_data)
  print(dynamics_plot)
  ggsave(file.path(base_output_dir_png, paste0(mag_to_plot, "_dynamics.png")), plot = dynamics_plot, width = 12, height = 8)
  ggsave(file.path(base_output_dir_svg, paste0(mag_to_plot, "_dynamics.svg")), plot = dynamics_plot, width = 12, height = 8)
} else {
  cat(paste("Skipping dynamics plot for", mag_to_plot, "- insufficient data.\n"))
}


cat("\n--- Analysis and Plotting Complete ---\n")








# --- X. Analysis X: Correlation Plots by Genus and MAG ---
cat("\n--- X. Generating Correlation Plots (Diversity vs. Alpha Diversity) ---\n")

# a. Define the Genus-level plotting function (MODIFIED)
plot_genus_correlation <- function(target_genus, data, pi_metric_col, metric_label, output_subdir) {
  
  # Ensure the metric column exists
  if (!pi_metric_col %in% colnames(data)) {
    cat(paste("Skipping Genus", target_genus, "- Metric column", pi_metric_col, "not found.\n"))
    return()
  }
  
  genus_data <- data %>% filter(Genus == target_genus)
  
  # Check if there's enough non-NA data for correlation
  valid_pairs <- sum(!is.na(genus_data[[pi_metric_col]]) & !is.na(genus_data$shannon_diversity))
  if (valid_pairs < 3) {
    cat(paste("Skipping Genus", target_genus, "- fewer than 3 valid data points for", metric_label, ".\n"))
    return()
  }
  
  corr_test <- tryCatch(
    cor.test(genus_data$shannon_diversity, genus_data[[pi_metric_col]]),
    error = function(e) NULL # Return NULL if cor.test fails
  )
  
  if (is.null(corr_test)) {
    cat(paste("Skipping Genus", target_genus, "- correlation test failed for", metric_label, ".\n"))
    return()
  }
  
  corr_label <- paste0("r = ", round(corr_test$estimate, 2), ", p = ", format.pval(corr_test$p.value, digits = 2))
  
  line_color <- "gray50" # Default grey
  if (!is.na(corr_test$p.value) && corr_test$p.value < 0.05) {
    if (!is.na(corr_test$estimate) && corr_test$estimate > 0) {
      line_color <- "#fb4934" # Red
    } else if (!is.na(corr_test$estimate) && corr_test$estimate < 0) {
      line_color <- "#83a598" # Blue/Green
    }
  }
  
  
  p <- ggplot(genus_data, aes(x = shannon_diversity, y = .data[[pi_metric_col]])) +
    geom_point(color = "black", alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = line_color, linewidth = 1.5) +
    annotate("text", x = -Inf, y = Inf, label = corr_label, hjust = -0.1, vjust = 1.5, size = 4) +
    labs(title = paste("Genus:", target_genus, "-", metric_label),
         subtitle = "Nucleotide Diversity vs. Community (Shannon) Diversity",
         x = "Shannon Diversity", y = paste("Nucleotide Diversity (", metric_label, ")")) +
    theme_bw()
  
  output_dir_full <- file.path(base_output_dir_svg, output_subdir) # Use SVG base dir
  dir.create(output_dir_full, showWarnings = FALSE, recursive = TRUE)
  svg_path <- file.path(output_dir_full, paste0(target_genus, "_", metric_label, "_correlation.svg"))
  
  svg(filename = svg_path, width = 6, height = 5)
  print(p)
  dev.off()
  
  cat(paste("Saved SVG plot for Genus:", target_genus, "-", metric_label, "\n"))
}

# b. Run the loop for Genus plots (for 1D and 4D)
unique_genera <- unique(full_data$Genus)
purrr::walk(unique_genera, ~plot_genus_correlation(.x, full_data, "mean_pi_1D", "pi_1D", "by_genus_svg_1D"))
purrr::walk(unique_genera, ~plot_genus_correlation(.x, full_data, "mean_pi_4D", "pi_4D", "by_genus_svg_4D"))

# c. Define the MAG-level plotting function (MODIFIED)
plot_mag_correlation <- function(target_mag, data, pi_metric_col, metric_label, output_subdir) {
  
  if (!pi_metric_col %in% colnames(data)) {
    cat(paste("Skipping MAG", target_mag, "- Metric column", pi_metric_col, "not found.\n"))
    return()
  }
  
  mag_data <- data %>% filter(MAG_id == target_mag)
  
  valid_pairs <- sum(!is.na(mag_data[[pi_metric_col]]) & !is.na(mag_data$shannon_diversity))
  if (valid_pairs < 3) {
    cat(paste("Skipping MAG", target_mag, "- fewer than 3 valid data points for", metric_label, ".\n"))
    return()
  }
  
  corr_test <- tryCatch(
    cor.test(mag_data$shannon_diversity, mag_data[[pi_metric_col]]),
    error = function(e) NULL
  )
  
  if (is.null(corr_test)) {
    cat(paste("Skipping MAG", target_mag, "- correlation test failed for", metric_label, ".\n"))
    return()
  }
  
  corr_label <- paste0("r = ", round(corr_test$estimate, 2), ", p = ", format.pval(corr_test$p.value, digits = 2))
  
  line_color <- "#665c54" # Default gruvbox grey
  if (!is.na(corr_test$p.value) && corr_test$p.value < 0.05) {
    if (!is.na(corr_test$estimate) && corr_test$estimate > 0) {
      line_color <- "#fb4934" # Red
    } else if (!is.na(corr_test$estimate) && corr_test$estimate < 0) {
      line_color <- "#458588" # Blue
    }
  }
  
  
  p <- ggplot(mag_data, aes(x = shannon_diversity, y = .data[[pi_metric_col]])) +
    geom_point(aes(color = as.factor(Colony)), alpha = 0.8) + # Color points by Colony
    geom_smooth(method = "lm", se = FALSE, color = line_color, linewidth = 1.2) +
    annotate("text", x = -Inf, y = Inf, label = corr_label, hjust = -0.1, vjust = 1.5, size = 4) +
    labs(title = paste("MAG:", target_mag, "-", metric_label),
         subtitle = "Nucleotide Diversity vs. Community (Shannon) Diversity",
         x = "Shannon Diversity", y = paste("Nucleotide Diversity (", metric_label, ")"),
         color = "Colony") + # Add legend title
    theme_bw()
  
  output_dir_full <- file.path(base_output_dir_png, output_subdir) # Use PNG base dir
  dir.create(output_dir_full, showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = file.path(output_dir_full, paste0(target_mag, "_", metric_label, "_correlation.png")),
         plot = p, width = 7, height = 5, dpi = 150) # Increased width slightly for legend
  
  cat(paste("Saved PNG plot for MAG:", target_mag, "-", metric_label, "\n"))
}

# d. Run the loop for MAG plots (for 1D and 4D)
unique_mags <- unique(full_data$MAG_id)
purrr::walk(unique_mags, ~plot_mag_correlation(.x, full_data, "mean_pi_1D", "pi_1D", "by_mag_id_png_1D"))
purrr::walk(unique_mags, ~plot_mag_correlation(.x, full_data, "mean_pi_4D", "pi_4D", "by_mag_id_png_4D"))

cat("\n--- All correlation plots generated successfully! ---\n")





# --- XIII. Plotting Faceted Temporal Dynamics by Genus within Colony (Detailed) ---
# Gruvbox palette definition
# Gruvbox palette definition
gruvbox_palette <- c(
  dark_bg = "#282828", blue = "#83a598", red = "#fb4934", orange = "#fe8019", grey = "#a89984"
)

# a. Define the faceted plotting function (Revised for Faceting, Points, SD)
plot_genus_dynamics_faceted <- function(target_colony, target_genus,
                                        gene_data, mag_data, ecology_data, month_order,
                                        output_dir_png, output_dir_svg) {
  
  # --- Filter Data for the specific Colony and Genus ---
  genus_mag_ids <- mag_data %>%
    filter(Colony == target_colony, Genus == target_genus) %>%
    distinct(MAG_id) %>% pull(MAG_id)
  
  if (length(genus_mag_ids) == 0) {
    cat(paste("Skipping:", target_colony, "-", target_genus, "- No MAGs found.\n"))
    return()
  }
  
  # Filter primary data frames
  genus_gene_data <- gene_data %>%
    filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>%
    filter(!is.na(Month)) %>%
    mutate(TimePoint = as.integer(Month))
  
  genus_mag_data <- mag_data %>%
    filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>%
    filter(!is.na(Month)) %>%
    mutate(TimePoint = as.integer(Month))
  
  ecology_genus_data <- ecology_data %>%
    filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids)
  
  # Check if enough data points exist *after* filtering
  if(nrow(genus_mag_data) == 0 || nrow(genus_gene_data) == 0 || nrow(ecology_genus_data) == 0) {
    cat(paste("Skipping:", target_colony, "-", target_genus, "- Insufficient data after initial filtering.\n"))
    return()
  }
  # Further check: Ensure at least one MAG has >= 3 timepoints in the averaged data
  min_timepoints_check <- genus_mag_data %>% group_by(MAG_id) %>% filter(n() >= 3)
  if (nrow(min_timepoints_check) == 0) {
    cat(paste("Skipping:", target_colony, "-", target_genus, "- No MAG has >= 3 time points.\n"))
    return()
  }
  
  # --- Prepare Abundance Data ---
  replicate_totals <- ecology_data %>% # Calculate totals across the whole colony/month
    filter(Colony == target_colony, Caste == "Worker") %>%
    group_by(Month, sample_name) %>%
    summarise(total_rep_coverage = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    filter(total_rep_coverage > 0)
  
  abundance_reps <- ecology_genus_data %>%
    inner_join(replicate_totals, by = c("Month", "sample_name"), relationship = "many-to-many") %>%
    mutate(rep_rel_abund = mean / total_rep_coverage) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>%
    mutate(TimePoint = as.integer(Month))
  
  abundance_summary <- abundance_reps %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>% # Ensure MAG_id is included
    summarise(
      mean_abund = mean(rep_rel_abund, na.rm = TRUE),
      sd_abund = sd(rep_rel_abund, na.rm = TRUE),
      .groups = 'drop' ) %>%
    mutate(sd_abund = ifelse(is.na(sd_abund), 0, sd_abund)) # Handle single replicate cases
  
  # --- Prepare Diversity Data (Mean +/- SD from Genes) ---
  diversity_summary <- genus_gene_data %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>% # Group by MAG_id
    summarise(
      mean_pi_1D_genes = mean(pi_1D, na.rm = TRUE),
      sd_pi_1D_genes = sd(pi_1D, na.rm = TRUE),
      mean_pi_4D_genes = mean(pi_4D, na.rm = TRUE),
      sd_pi_4D_genes = sd(pi_4D, na.rm = TRUE),
      .groups = 'drop' ) %>%
    # Handle single gene cases for SD
    mutate(sd_pi_1D_genes = ifelse(is.na(sd_pi_1D_genes), 0, sd_pi_1D_genes),
           sd_pi_4D_genes = ifelse(is.na(sd_pi_4D_genes), 0, sd_pi_4D_genes))
  
  # --- Determine Shared Y-Axis Limits ---
  # Abundance limits (based on mean +/- sd for visual consistency with ribbon)
  abund_range_data <- abundance_summary %>%
    mutate(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund)
  y_limits_abund <- c(max(0, min(abund_range_data$ymin, na.rm=T)), max(abund_range_data$ymax, na.rm=T))
  # Add buffer if desired
  buffer_abund <- (y_limits_abund[2] - y_limits_abund[1]) * 0.05
  y_limits_abund <- c(max(0, y_limits_abund[1] - buffer_abund), y_limits_abund[2] + buffer_abund)
  
  
  # Combined Pi limits (considering individual points and mean +/- sd range)
  pi_range_data_summary <- diversity_summary %>%
    mutate(ymin_1D = mean_pi_1D_genes - sd_pi_1D_genes, ymax_1D = mean_pi_1D_genes + sd_pi_1D_genes,
           ymin_4D = mean_pi_4D_genes - sd_pi_4D_genes, ymax_4D = mean_pi_4D_genes + sd_pi_4D_genes)
  
  all_pi_values <- c(genus_gene_data$pi_1D, genus_gene_data$pi_4D,
                     pi_range_data_summary$ymin_1D, pi_range_data_summary$ymax_1D,
                     pi_range_data_summary$ymin_4D, pi_range_data_summary$ymax_4D)
  valid_pi_values <- all_pi_values[!is.na(all_pi_values) & is.finite(all_pi_values)]
  
  if (length(valid_pi_values) > 0) {
    min_pi <- min(valid_pi_values); max_pi <- max(valid_pi_values)
    buffer_pi <- (max_pi - min_pi) * 0.05
    y_limits_pi <- c(max(0, min_pi - buffer_pi), max_pi + buffer_pi)
  } else { y_limits_pi <- c(0, 0.01) } # Default
  
  # --- Filter Data for Plotting (Remove NAs needed for plotting) ---
  plot_abundance_summary <- abundance_summary %>% filter(!is.na(TimePoint), !is.na(mean_abund))
  plot_diversity_summary <- diversity_summary %>% filter(!is.na(TimePoint))
  plot_gene_pi <- genus_gene_data %>% filter(!is.na(TimePoint)) # Keep points even if pi is NA
  
  # --- Build Plots ---
  x_axis_scale <- scale_x_continuous(breaks = 1:length(month_order), labels = month_order)
  point_alpha <- 0.3
  point_size <- 0.8
  line_width <- 0.8
  ribbon_alpha <- 0.3
  
  # PLOT 1: Abundance
  p_abund <- ggplot(plot_abundance_summary, aes(x = TimePoint, y = mean_abund)) +
    geom_ribbon(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund),
                fill = gruvbox_palette["orange"], alpha = ribbon_alpha) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = line_width) +
    # geom_point(color = gruvbox_palette["dark_bg"], size = point_size + 0.2) + # Optional points for mean
    facet_wrap(~ MAG_id, nrow = 1) + # Facet by MAG
    coord_cartesian(ylim = y_limits_abund) + # Shared Y axis for abundance
    labs(y = "Relative Abundance") + x_axis_scale + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey90"),
          panel.spacing.x = unit(0.1, "lines")) # Reduce space between facets

  
  # --- Build Plots ---
  x_axis_scale <- scale_x_continuous(breaks = 1:length(month_order), labels = month_order)
  
  # --- ADD THESE DEFINITIONS ---
  gene_line_alpha <- 0.1   # Faint lines for individual genes
  gene_line_width <- 0.3   # Thin lines for individual genes
  mean_line_width <- 0.8   # Thicker line for the mean
  ribbon_alpha    <- 0.3   # Transparency for the SD ribbon
  # --- END OF ADDITIONS ---
  
  
  # PLOT 1: Abundance (Use mean_line_width, ribbon_alpha)
  p_abund <- ggplot(plot_abundance_summary, aes(x = TimePoint, y = mean_abund)) +
    geom_ribbon(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund),
                fill = gruvbox_palette["orange"], alpha = ribbon_alpha) + # Use variable
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) + # Use variable
    geom_point(color = gruvbox_palette["dark_bg"], size = 1.5) +
    facet_wrap(~ MAG_id, nrow = 1) +
    coord_cartesian(ylim = y_limits_abund) +
    labs(y = "Relative Abundance") + x_axis_scale + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey90"),
          panel.spacing.x = unit(0.1, "lines"))
  
  # PLOT 2: Diversity (1D) (Use gene_line_alpha, gene_line_width, mean_line_width, ribbon_alpha)
  p_pi_1D <- ggplot(plot_diversity_summary, aes(x = TimePoint, y = mean_pi_1D_genes)) +
    geom_line(data = plot_gene_pi, aes(y = pi_1D, group = corresponding_gene_call),
              color = gruvbox_palette["blue"], alpha = gene_line_alpha, linewidth = gene_line_width, na.rm = TRUE) + # Use variables
    geom_ribbon(aes(ymin = mean_pi_1D_genes - sd_pi_1D_genes, ymax = mean_pi_1D_genes + sd_pi_1D_genes),
                fill = gruvbox_palette["blue"], alpha = ribbon_alpha) + # Use variable
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) + # Use variable
    coord_cartesian(ylim = y_limits_pi) + facet_wrap(~ MAG_id, nrow = 1) +
    labs(x = "Month", y = expression(paste(pi ["1D"]))) +
    x_axis_scale + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x = element_blank(),
          panel.spacing.x = unit(0.1, "lines"))
  
  # PLOT 3: Diversity (4D) (Use gene_line_alpha, gene_line_width, mean_line_width, ribbon_alpha)
  p_pi_4D <- ggplot(plot_diversity_summary, aes(x = TimePoint, y = mean_pi_4D_genes)) +
    geom_line(data = plot_gene_pi, aes(y = pi_4D, group = corresponding_gene_call),
              color = gruvbox_palette["red"], alpha = gene_line_alpha, linewidth = gene_line_width, na.rm = TRUE) + # Use variables
    geom_ribbon(aes(ymin = mean_pi_4D_genes - sd_pi_4D_genes, ymax = mean_pi_4D_genes + sd_pi_4D_genes),
                fill = gruvbox_palette["red"], alpha = ribbon_alpha) + # Use variable
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) + # Use variable
    coord_cartesian(ylim = y_limits_pi) + facet_wrap(~ MAG_id, nrow = 1) +
    labs(x = "Month", y = expression(paste(pi ["4D"]))) +
    x_axis_scale + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x = element_blank(),
          panel.spacing.x = unit(0.1, "lines"))
  
  
   #Also update the subtitle in the plot_annotation call:
   final_plot <- (p_abund / p_pi_1D / p_pi_4D) +
    plot_annotation(
      title = paste("Dynamics in Colony", target_colony, "- Genus:", target_genus),
      subtitle = "Facets show individual MAGs. Faint lines = gene pi; Line/Ribbon = MAG mean +/- SD." # Updated subtitle
    ) &
    theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))
  
  # Use target_colony and target_genus for path
  colony_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", target_colony)
  genus_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", target_genus)
  
  # Adjust width based on number of MAGs (facets)
  n_facets <- length(genus_mag_ids)
  plot_width <- max(6, 2 * n_facets) # Adjust base width and multiplier as needed
  
  # Save PNG
  output_dir_p_nested <- file.path(output_dir_png, "genus_dynamics_faceted", colony_name_safe)
  dir.create(output_dir_p_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_p <- paste0(genus_name_safe, "_dynamics_faceted.png")
  ggsave(file.path(output_dir_p_nested, file_name_p), plot = final_plot, width = plot_width, height = 8, dpi = 150, bg = "white") # Slightly less height
  
  # Save SVG
  output_dir_s_nested <- file.path(output_dir_svg, "genus_dynamics_faceted", colony_name_safe)
  dir.create(output_dir_s_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_s <- paste0(genus_name_safe, "_dynamics_faceted.svg")
  svg(filename = file.path(output_dir_s_nested, file_name_s), width = plot_width, height = 8)
  print(final_plot); dev.off()
  
  cat(paste("Saved faceted dynamics plots for Genus:", target_genus, "in Colony:", target_colony, "\n"))
}

# b. Run the loop for faceted dynamics plots
cat("\n--- Running loop to generate faceted MAG dynamics plots ---\n")
# Get unique Colony-Genus pairs present in the data
colony_genus_pairs <- full_data %>% # Use full_data as it's aggregated and filtered
  distinct(Colony, Genus) %>%
  filter(!is.na(Genus), Genus != "") # Ensure Genus is valid

# Ensure gene_level_data, full_data, ecology_df are prepared correctly
purrr::pwalk(colony_genus_pairs, ~plot_genus_dynamics_faceted(target_colony = ..1, target_genus = ..2,
                                                              gene_data = gene_level_data,
                                                              mag_data = full_data,
                                                              ecology_data = ecology_df,
                                                              month_order = month_order,
                                                              output_dir_png = base_output_dir_png,
                                                              output_dir_svg = base_output_dir_svg))

cat("\n--- All faceted MAG dynamics plots generation attempted! ---\n")




















########### SAME 4D SCALE, ALL MAGS
# --- XIV. Calculate Global Y-Limits for Faceted Plots ---
cat("\n--- XIV. Calculating Global Y-Limits ---\n")

# Ensure gene_level_data, full_data, ecology_df are prepared and filtered

# a. Global Abundance Limits (based on Mean +/- SD)
# Recalculate abundance summary across ALL relevant data
global_abundance_reps <- ecology_df %>% # Use the filtered ecology_df
  filter(MAG_id %in% unique(full_data$MAG_id)) %>% # Only MAGs passing initial filters
  inner_join(
    ecology_df %>% filter(Caste == "Worker") %>%
      group_by(Colony, Month, sample_name) %>%
      summarise(total_rep_coverage = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
      filter(total_rep_coverage > 0),
    by = c("Colony", "Month", "sample_name"), relationship = "many-to-many"
  ) %>%
  mutate(rep_rel_abund = mean / total_rep_coverage)

global_abundance_summary <- global_abundance_reps %>%
  group_by(MAG_id, Colony, Month) %>% # Group to get mean/sd per timepoint
  summarise(
    mean_abund = mean(rep_rel_abund, na.rm = TRUE),
    sd_abund = sd(rep_rel_abund, na.rm = TRUE),
    .groups = 'drop' ) %>%
  mutate(sd_abund = ifelse(is.na(sd_abund), 0, sd_abund)) %>%
  mutate(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund)

# Determine global range
global_y_limits_abund <- c(0, max(global_abundance_summary$ymax, na.rm = TRUE))
# Add buffer
buffer_abund_global <- (global_y_limits_abund[2] - global_y_limits_abund[1]) * 0.05
global_y_limits_abund <- c(0, global_y_limits_abund[2] + buffer_abund_global)
cat("Global Abundance Y-Limits:", global_y_limits_abund, "\n")


# b. Global 4D Pi Limits (based on Genes + Mean +/- SD)
# Calculate Mean +/- SD for ALL relevant gene data
global_diversity_summary <- gene_level_data %>% # Use the filtered gene_level_data
  filter(MAG_id %in% unique(full_data$MAG_id)) %>%
  group_by(MAG_id, Colony, Month) %>%
  summarise(
    mean_pi_4D_genes = mean(pi_4D, na.rm = TRUE),
    sd_pi_4D_genes = sd(pi_4D, na.rm = TRUE),
    .groups = 'drop' ) %>%
  mutate(sd_pi_4D_genes = ifelse(is.na(sd_pi_4D_genes), 0, sd_pi_4D_genes)) %>%
  mutate(ymin_4D = mean_pi_4D_genes - sd_pi_4D_genes,
         ymax_4D = mean_pi_4D_genes + sd_pi_4D_genes)

# Combine gene points and summary ranges
global_all_pi4D_values <- c(
  gene_level_data %>% filter(MAG_id %in% unique(full_data$MAG_id)) %>% pull(pi_4D),
  global_diversity_summary$ymin_4D,
  global_diversity_summary$ymax_4D
)
global_valid_pi4D_values <- global_all_pi4D_values[!is.na(global_all_pi4D_values) & is.finite(global_all_pi4D_values)]

if (length(global_valid_pi4D_values) > 0) {
  min_pi4D_global <- min(global_valid_pi4D_values); max_pi4D_global <- max(global_valid_pi4D_values)
  buffer_pi4D_global <- (max_pi4D_global - min_pi4D_global) * 0.05
  global_y_limits_pi4D <- c(max(0, min_pi4D_global - buffer_pi4D_global), max_pi4D_global + buffer_pi4D_global)
} else { global_y_limits_pi4D <- c(0, 0.01) } # Default
cat("Global 4D Pi Y-Limits:", global_y_limits_pi4D, "\n")



# c. Define the faceted plotting function (4D only, Global Y)
plot_genus_dynamics_faceted_4D_global <- function(target_colony, target_genus,
                                                  gene_data, mag_data, ecology_data, month_order,
                                                  global_ylim_abund, global_ylim_pi4D, # NEW arguments
                                                  output_dir_png, output_dir_svg) {
  
  # --- Filter Data (same as before) ---
  genus_mag_ids <- mag_data %>% filter(Colony == target_colony, Genus == target_genus) %>% distinct(MAG_id) %>% pull(MAG_id)
  if (length(genus_mag_ids) == 0) { return() }
  genus_gene_data <- gene_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  genus_mag_data <- mag_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  ecology_genus_data <- ecology_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids)
  if(nrow(genus_mag_data) == 0 || nrow(genus_gene_data) == 0 || nrow(ecology_genus_data) == 0) { return() }
  min_timepoints_check <- genus_mag_data %>% group_by(MAG_id) %>% filter(n() >= 3)
  if (nrow(min_timepoints_check) == 0) { return() }
  
  # --- Prepare Abundance Data (same as before) ---
  replicate_totals <- ecology_data %>% filter(Colony == target_colony, Caste == "Worker") %>%
    group_by(Month, sample_name) %>% summarise(total_rep_coverage = sum(mean, na.rm = TRUE), .groups = 'drop') %>% filter(total_rep_coverage > 0)
  abundance_reps <- ecology_genus_data %>%
    inner_join(replicate_totals, by = c("Month", "sample_name"), relationship = "many-to-many") %>%
    mutate(rep_rel_abund = mean / total_rep_coverage) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  abundance_summary <- abundance_reps %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>%
    summarise(mean_abund = mean(rep_rel_abund, na.rm = TRUE), sd_abund = sd(rep_rel_abund, na.rm = TRUE), .groups = 'drop' ) %>%
    mutate(sd_abund = ifelse(is.na(sd_abund), 0, sd_abund))
  
  # --- Prepare Diversity Data (Mean +/- SD from Genes - 4D only needed now) ---
  diversity_summary <- genus_gene_data %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>%
    summarise(mean_pi_4D_genes = mean(pi_4D, na.rm = TRUE), sd_pi_4D_genes = sd(pi_4D, na.rm = TRUE), .groups = 'drop' ) %>%
    mutate(sd_pi_4D_genes = ifelse(is.na(sd_pi_4D_genes), 0, sd_pi_4D_genes))
  
  # --- Filter Data for Plotting ---
  plot_abundance_summary <- abundance_summary %>% filter(!is.na(TimePoint), !is.na(mean_abund))
  plot_diversity_summary <- diversity_summary %>% filter(!is.na(TimePoint))
  plot_gene_pi <- genus_gene_data %>% filter(!is.na(TimePoint), !is.na(pi_4D), is.finite(pi_4D)) # Filter for 4D lines
  
  # --- Build Plots (Abundance and 4D Pi only) ---
  x_axis_scale <- scale_x_continuous(breaks = 1:length(month_order), labels = month_order)
  gene_line_alpha <- 0.1; gene_line_width <- 0.3; mean_line_width <- 0.8; ribbon_alpha <- 0.3
  
  # PLOT 1: Abundance (Apply Global Limit)
  p_abund <- ggplot(plot_abundance_summary, aes(x = TimePoint, y = mean_abund)) +
    geom_ribbon(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), fill = gruvbox_palette["orange"], alpha = ribbon_alpha) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) +
    geom_point(color = gruvbox_palette["dark_bg"], size = 1.5) +
    facet_wrap(~ MAG_id, nrow = 1) +
    coord_cartesian(ylim = global_ylim_abund) + # USE GLOBAL LIMIT
    labs(y = "Relative Abundance") + x_axis_scale + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey90"), panel.spacing.x = unit(0.1, "lines"))
  
  # PLOT 2: Diversity (4D) (Apply Global Limit)
  p_pi_4D <- ggplot(plot_diversity_summary, aes(x = TimePoint, y = mean_pi_4D_genes)) +
    geom_line(data = plot_gene_pi, aes(y = pi_4D, group = corresponding_gene_call), color = gruvbox_palette["red"], alpha = gene_line_alpha, linewidth = gene_line_width, na.rm = TRUE) +
    geom_ribbon(aes(ymin = mean_pi_4D_genes - sd_pi_4D_genes, ymax = mean_pi_4D_genes + sd_pi_4D_genes), fill = gruvbox_palette["red"], alpha = ribbon_alpha) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) +
    facet_wrap(~ MAG_id, nrow = 1) +
    coord_cartesian(ylim = global_ylim_pi4D) + # USE GLOBAL LIMIT
    labs(x = "Month", y = expression(paste(pi ["4D"]))) +
    x_axis_scale + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x = element_blank(),
          panel.spacing.x = unit(0.1, "lines"))
  
  
  # --- Combine and Save (Abundance / Pi_4D) ---
  final_plot <- (p_abund / p_pi_4D) + # Combine only two plots
    plot_annotation(
      title = paste("Dynamics in Colony", target_colony, "- Genus:", target_genus, "(Global Y-Scale)"), # Updated title
      subtitle = "Facets show individual MAGs. Top: Rel Abund (Avg +/- SD). Bottom: Pi_4D (Avg +/- SD, Faint lines = gene pi)." # Updated subtitle
    ) &
    theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))
  
  # Get Colony/Genus info and sanitize names (same as before)
  colony_name <- target_colony; genus_name <- target_genus
  colony_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", colony_name)
  genus_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", genus_name)
  
  # Determine plot width (same as before)
  n_facets <- length(genus_mag_ids)
  plot_width <- max(6, 2 * n_facets)
  plot_height <- 6 # Reduced height for 2 panels
  
  # Save PNG (update directory/filename)
  output_dir_p_nested <- file.path(output_dir_png, "genus_dynamics_faceted_4D_global", colony_name_safe)
  dir.create(output_dir_p_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_p <- paste0(genus_name_safe, "_dynamics_4D_global.png")
  ggsave(file.path(output_dir_p_nested, file_name_p), plot = final_plot, width = plot_width, height = plot_height, dpi = 150, bg = "white")
  
  # Save SVG (update directory/filename)
  output_dir_s_nested <- file.path(output_dir_svg, "genus_dynamics_faceted_4D_global", colony_name_safe)
  dir.create(output_dir_s_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_s <- paste0(genus_name_safe, "_dynamics_4D_global.svg")
  svg(filename = file.path(output_dir_s_nested, file_name_s), width = plot_width, height = plot_height)
  print(final_plot); dev.off()
  
  cat(paste("Saved 4D/Global faceted plots for Genus:", target_genus, "in Colony:", target_colony, "\n"))
}



# c. Define the faceted plotting function (4D only, Global Y)
plot_genus_dynamics_faceted_4D_global <- function(target_colony, target_genus,
                                                  gene_data, mag_data, ecology_data, month_order,
                                                  global_ylim_abund, global_ylim_pi4D, # NEW arguments
                                                  output_dir_png, output_dir_svg) {
  
  # --- Filter Data (same as before) ---
  genus_mag_ids <- mag_data %>% filter(Colony == target_colony, Genus == target_genus) %>% distinct(MAG_id) %>% pull(MAG_id)
  if (length(genus_mag_ids) == 0) { return() }
  genus_gene_data <- gene_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  genus_mag_data <- mag_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  ecology_genus_data <- ecology_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids)
  if(nrow(genus_mag_data) == 0 || nrow(genus_gene_data) == 0 || nrow(ecology_genus_data) == 0) { return() }
  min_timepoints_check <- genus_mag_data %>% group_by(MAG_id) %>% filter(n() >= 3)
  if (nrow(min_timepoints_check) == 0) { return() }
  
  # --- Prepare Abundance Data (same as before) ---
  replicate_totals <- ecology_data %>% filter(Colony == target_colony, Caste == "Worker") %>%
    group_by(Month, sample_name) %>% summarise(total_rep_coverage = sum(mean, na.rm = TRUE), .groups = 'drop') %>% filter(total_rep_coverage > 0)
  abundance_reps <- ecology_genus_data %>%
    inner_join(replicate_totals, by = c("Month", "sample_name"), relationship = "many-to-many") %>%
    mutate(rep_rel_abund = mean / total_rep_coverage) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  abundance_summary <- abundance_reps %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>%
    summarise(mean_abund = mean(rep_rel_abund, na.rm = TRUE), sd_abund = sd(rep_rel_abund, na.rm = TRUE), .groups = 'drop' ) %>%
    mutate(sd_abund = ifelse(is.na(sd_abund), 0, sd_abund))
  
  # --- Prepare Diversity Data (Mean +/- SD from Genes - 4D only needed now) ---
  diversity_summary <- genus_gene_data %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>%
    summarise(mean_pi_4D_genes = mean(pi_4D, na.rm = TRUE), sd_pi_4D_genes = sd(pi_4D, na.rm = TRUE), .groups = 'drop' ) %>%
    mutate(sd_pi_4D_genes = ifelse(is.na(sd_pi_4D_genes), 0, sd_pi_4D_genes))
  
  # --- Filter Data for Plotting ---
  plot_abundance_summary <- abundance_summary %>% filter(!is.na(TimePoint), !is.na(mean_abund))
  plot_diversity_summary <- diversity_summary %>% filter(!is.na(TimePoint))
  plot_gene_pi <- genus_gene_data %>% filter(!is.na(TimePoint), !is.na(pi_4D), is.finite(pi_4D)) # Filter for 4D lines
  
  # --- Build Plots (Abundance and 4D Pi only) ---
  x_axis_scale <- scale_x_continuous(breaks = 1:length(month_order), labels = month_order)
  gene_line_alpha <- 0.1; gene_line_width <- 0.3; mean_line_width <- 0.8; ribbon_alpha <- 0.3
  
  # PLOT 1: Abundance (Apply Global Limit)
  p_abund <- ggplot(plot_abundance_summary, aes(x = TimePoint, y = mean_abund)) +
    geom_ribbon(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), fill = gruvbox_palette["orange"], alpha = ribbon_alpha) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) +
    geom_point(color = gruvbox_palette["dark_bg"], size = 1.5) +
    facet_wrap(~ MAG_id, nrow = 1) +
    coord_cartesian(ylim = global_ylim_abund) + # USE GLOBAL LIMIT
    labs(y = "Relative Abundance") + x_axis_scale + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey90"), panel.spacing.x = unit(0.1, "lines"))
  
  # PLOT 2: Diversity (4D) (Apply Global Limit)
  p_pi_4D <- ggplot(plot_diversity_summary, aes(x = TimePoint, y = mean_pi_4D_genes)) +
    geom_line(data = plot_gene_pi, aes(y = pi_4D, group = corresponding_gene_call), color = gruvbox_palette["red"], alpha = gene_line_alpha, linewidth = gene_line_width, na.rm = TRUE) +
    geom_ribbon(aes(ymin = mean_pi_4D_genes - sd_pi_4D_genes, ymax = mean_pi_4D_genes + sd_pi_4D_genes), fill = gruvbox_palette["red"], alpha = ribbon_alpha) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) +
    facet_wrap(~ MAG_id, nrow = 1) +
    coord_cartesian(ylim = global_ylim_pi4D) + # USE GLOBAL LIMIT
    labs(x = "Month", y = expression(paste(pi ["4D"]))) +
    x_axis_scale + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x = element_blank(),
          panel.spacing.x = unit(0.1, "lines"))
  
  
  # --- Combine and Save (Abundance / Pi_4D) ---
  final_plot <- (p_abund / p_pi_4D) + # Combine only two plots
    plot_annotation(
      title = paste("Dynamics in Colony", target_colony, "- Genus:", target_genus, "(Global Y-Scale)"), # Updated title
      subtitle = "Facets show individual MAGs. Top: Rel Abund (Avg +/- SD). Bottom: Pi_4D (Avg +/- SD, Faint lines = gene pi)." # Updated subtitle
    ) &
    theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))
  
  # Get Colony/Genus info and sanitize names (same as before)
  colony_name <- target_colony; genus_name <- target_genus
  colony_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", colony_name)
  genus_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", genus_name)
  
  # Determine plot width (same as before)
  n_facets <- length(genus_mag_ids)
  plot_width <- max(6, 2 * n_facets)
  plot_height <- 6 # Reduced height for 2 panels
  
  # Save PNG (update directory/filename)
  output_dir_p_nested <- file.path(output_dir_png, "genus_dynamics_faceted_4D_global", colony_name_safe)
  dir.create(output_dir_p_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_p <- paste0(genus_name_safe, "_dynamics_4D_global.png")
  ggsave(file.path(output_dir_p_nested, file_name_p), plot = final_plot, width = plot_width, height = plot_height, dpi = 150, bg = "white")
  
  # Save SVG (update directory/filename)
  output_dir_s_nested <- file.path(output_dir_svg, "genus_dynamics_faceted_4D_global", colony_name_safe)
  dir.create(output_dir_s_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_s <- paste0(genus_name_safe, "_dynamics_4D_global.svg")
  svg(filename = file.path(output_dir_s_nested, file_name_s), width = plot_width, height = plot_height)
  print(final_plot); dev.off()
  
  cat(paste("Saved 4D/Global faceted plots for Genus:", target_genus, "in Colony:", target_colony, "\n"))
}


# d. Run the loop for 4D/Global faceted dynamics plots
cat("\n--- Running loop to generate 4D/Global faceted MAG dynamics plots ---\n")
# Get unique Colony-Genus pairs (same as before)
colony_genus_pairs <- full_data %>% distinct(Colony, Genus) %>% filter(!is.na(Genus), Genus != "")

# Rename columns to match function arguments (same as before)
colony_genus_pairs_renamed <- colony_genus_pairs %>%
  rename(target_colony = Colony, target_genus = Genus)

# Ensure gene_level_data, full_data, ecology_df are prepared correctly
# Use the NEW function and pass the global limits
purrr::pwalk(colony_genus_pairs_renamed,
             plot_genus_dynamics_faceted_4D_global, # Use the new function
             # Pass other arguments
             gene_data = gene_level_data,
             mag_data = full_data,
             ecology_data = ecology_df,
             month_order = month_order,
             global_ylim_abund = global_y_limits_abund, # Pass global abund limits
             global_ylim_pi4D = global_y_limits_pi4D,    # Pass global pi 4D limits
             output_dir_png = base_output_dir_png,
             output_dir_svg = base_output_dir_svg)

cat("\n--- All 4D/Global faceted MAG dynamics plots generation attempted! ---\n")













#new 4D mean pi and FIXED

# --- XV. Plotting Faceted Dynamics (4D Global Y, Abundance Variable Y) ---
cat("\n--- XV. Generating Faceted Plots (4D Global Y, Abundance Variable Y) ---\n")

# a. Define the faceted plotting function (Variable Abund Y)
plot_genus_dynamics_faceted_4D_variable_abund <- function(target_colony, target_genus,
                                                          gene_data, mag_data, ecology_data, month_order,
                                                          global_ylim_pi4D, # Keep global Pi limit
                                                          output_dir_png, output_dir_svg) {
  
  # --- Filter Data (same as before) ---
  genus_mag_ids <- mag_data %>% filter(Colony == target_colony, Genus == target_genus) %>% distinct(MAG_id) %>% pull(MAG_id)
  if (length(genus_mag_ids) == 0) { return() }
  genus_gene_data <- gene_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  genus_mag_data <- mag_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  ecology_genus_data <- ecology_data %>% filter(Colony == target_colony, Genus == target_genus, MAG_id %in% genus_mag_ids)
  if(nrow(genus_mag_data) == 0 || nrow(genus_gene_data) == 0 || nrow(ecology_genus_data) == 0) { return() }
  min_timepoints_check <- genus_mag_data %>% group_by(MAG_id) %>% filter(n() >= 3)
  if (nrow(min_timepoints_check) == 0) { return() }
  
  # --- Prepare Abundance Data (same as before) ---
  replicate_totals <- ecology_data %>% filter(Colony == target_colony, Caste == "Worker") %>%
    group_by(Month, sample_name) %>% summarise(total_rep_coverage = sum(mean, na.rm = TRUE), .groups = 'drop') %>% filter(total_rep_coverage > 0)
  abundance_reps <- ecology_genus_data %>%
    inner_join(replicate_totals, by = c("Month", "sample_name"), relationship = "many-to-many") %>%
    mutate(rep_rel_abund = mean / total_rep_coverage) %>%
    mutate(Month = factor(Month, levels = month_order)) %>% filter(!is.na(Month)) %>% mutate(TimePoint = as.integer(Month))
  abundance_summary <- abundance_reps %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>%
    summarise(mean_abund = mean(rep_rel_abund, na.rm = TRUE), sd_abund = sd(rep_rel_abund, na.rm = TRUE), .groups = 'drop' ) %>%
    mutate(sd_abund = ifelse(is.na(sd_abund), 0, sd_abund))
  
  # --- Prepare Diversity Data (4D only needed - same as before) ---
  diversity_summary <- genus_gene_data %>%
    group_by(Colony, Month, TimePoint, MAG_id) %>%
    summarise(mean_pi_4D_genes = mean(pi_4D, na.rm = TRUE), sd_pi_4D_genes = sd(pi_4D, na.rm = TRUE), .groups = 'drop' ) %>%
    mutate(sd_pi_4D_genes = ifelse(is.na(sd_pi_4D_genes), 0, sd_pi_4D_genes))
  
  # --- Determine LOCAL Abundance Y-Axis Limits ---
  # Calculate limits based only on the current group's abundance data
  abund_range_data_local <- abundance_summary %>%
    mutate(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund)
  y_limits_abund_local <- c(0, max(abund_range_data_local$ymax, na.rm=T))
  # Add buffer if desired
  buffer_abund_local <- (y_limits_abund_local[2] - y_limits_abund_local[1]) * 0.05
  y_limits_abund_local <- c(0, y_limits_abund_local[2] + buffer_abund_local)
  
  # --- Filter Data for Plotting ---
  plot_abundance_summary <- abundance_summary %>% filter(!is.na(TimePoint), !is.na(mean_abund))
  plot_diversity_summary <- diversity_summary %>% filter(!is.na(TimePoint))
  plot_gene_pi <- genus_gene_data %>% filter(!is.na(TimePoint), !is.na(pi_4D), is.finite(pi_4D))
  
  # --- Build Plots (Abundance uses LOCAL Limit, Pi uses GLOBAL Limit) ---
  x_axis_scale <- scale_x_continuous(breaks = 1:length(month_order), labels = month_order)
  gene_line_alpha <- 0.1; gene_line_width <- 0.3; mean_line_width <- 0.8; ribbon_alpha <- 0.3
  
  # PLOT 1: Abundance (Apply LOCAL Limit)
  p_abund <- ggplot(plot_abundance_summary, aes(x = TimePoint, y = mean_abund)) +
    geom_ribbon(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), fill = gruvbox_palette["orange"], alpha = ribbon_alpha) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) +
    geom_point(color = gruvbox_palette["dark_bg"], size = 1.5) +
    facet_wrap(~ MAG_id, nrow = 1) +
    coord_cartesian(ylim = y_limits_abund_local) + # USE LOCAL LIMIT
    labs(y = "Relative Abundance") + x_axis_scale + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey90"), panel.spacing.x = unit(0.1, "lines"))
  
  # PLOT 2: Diversity (4D) (Apply GLOBAL Limit)
  p_pi_4D <- ggplot(plot_diversity_summary, aes(x = TimePoint, y = mean_pi_4D_genes)) +
    geom_line(data = plot_gene_pi, aes(y = pi_4D, group = corresponding_gene_call), color = gruvbox_palette["red"], alpha = gene_line_alpha, linewidth = gene_line_width, na.rm = TRUE) +
    geom_ribbon(aes(ymin = mean_pi_4D_genes - sd_pi_4D_genes, ymax = mean_pi_4D_genes + sd_pi_4D_genes), fill = gruvbox_palette["red"], alpha = ribbon_alpha) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = mean_line_width) +
    facet_wrap(~ MAG_id, nrow = 1) +
    coord_cartesian(ylim = global_ylim_pi4D) + # USE GLOBAL LIMIT
    labs(x = "Month", y = expression(paste(pi ["4D"]))) +
    x_axis_scale + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x = element_blank(),
          panel.spacing.x = unit(0.1, "lines"))
  
  # --- Combine and Save ---
  final_plot <- (p_abund / p_pi_4D) +
    plot_annotation(
      title = paste("Dynamics in Colony", target_colony, "- Genus:", target_genus, "(Pi 4D Global Y)"), # Updated title
      subtitle = "Facets show individual MAGs. Top: Rel Abund (Avg +/- SD). Bottom: Pi_4D (Avg +/- SD, Faint lines = gene pi)."
    ) &
    theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))
  
  # Get Colony/Genus info and sanitize names (same as before)
  colony_name <- target_colony; genus_name <- target_genus
  colony_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", colony_name)
  genus_name_safe <- gsub("[^A-Za-z0-9_.-]", "_", genus_name)
  
  # Determine plot width (same as before)
  n_facets <- length(genus_mag_ids)
  plot_width <- max(6, 2 * n_facets)
  plot_height <- 6 # Reduced height for 2 panels
  
  # Save PNG (update directory/filename)
  output_dir_p_nested <- file.path(output_dir_png, "genus_dynamics_faceted_4Dglobal_AbundVar", colony_name_safe) # New subdir
  dir.create(output_dir_p_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_p <- paste0(genus_name_safe, "_dynamics_4Dglobal_AbundVar.png") # New filename part
  ggsave(file.path(output_dir_p_nested, file_name_p), plot = final_plot, width = plot_width, height = plot_height, dpi = 150, bg = "white")
  
  # Save SVG (update directory/filename)
  output_dir_s_nested <- file.path(output_dir_svg, "genus_dynamics_faceted_4Dglobal_AbundVar", colony_name_safe) # New subdir
  dir.create(output_dir_s_nested, showWarnings = FALSE, recursive = TRUE)
  file_name_s <- paste0(genus_name_safe, "_dynamics_4Dglobal_AbundVar.svg") # New filename part
  svg(filename = file.path(output_dir_s_nested, file_name_s), width = plot_width, height = plot_height)
  print(final_plot); dev.off()
  
  cat(paste("Saved 4D-Global/Abund-Var faceted plots for Genus:", target_genus, "in Colony:", target_colony, "\n"))
}

# b. Run the loop for 4D-Global/Abund-Variable faceted dynamics plots
cat("\n--- Running loop to generate 4D-Global/Abund-Variable faceted MAG dynamics plots ---\n")
# Get unique Colony-Genus pairs (same as before)
colony_genus_pairs <- full_data %>% distinct(Colony, Genus) %>% filter(!is.na(Genus), Genus != "")

# Rename columns to match function arguments (same as before)
colony_genus_pairs_renamed <- colony_genus_pairs %>%
  rename(target_colony = Colony, target_genus = Genus)

# Ensure gene_level_data, full_data, ecology_df, and global_y_limits_pi4D are ready
# Use the NEW function and pass only the global Pi limit
purrr::pwalk(colony_genus_pairs_renamed,
             plot_genus_dynamics_faceted_4D_variable_abund, # Use the new function
             # Pass other arguments
             gene_data = gene_level_data,
             mag_data = full_data,
             ecology_data = ecology_df,
             month_order = month_order,
             # global_ylim_abund = global_y_limits_abund, # DO NOT PASS global abund limits
             global_ylim_pi4D = global_y_limits_pi4D,    # Pass global pi 4D limits
             output_dir_png = base_output_dir_png,
             output_dir_svg = base_output_dir_svg)

cat("\n--- All 4D-Global/Abund-Variable faceted MAG dynamics plots generation attempted! ---\n")
