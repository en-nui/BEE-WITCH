library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
file_path <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output/all_filtered_mmag_contigs_concatenated.csv"
tryCatch({
  df <- read.csv(file_path)
  print("File loaded successfully.")
}, error = function(e) {
  stop(paste("Error loading file:", e$message, "\nPlease check the file path and permissions."))
})

# 3. Identify all possible Colony,Month combinations
# We use a distinct() call to get all unique pairs of Colony and Month
unique_combinations <- df %>%
  distinct(Colony, Month)

total_combinations <- nrow(unique_combinations)
print(paste("Total number of unique Colony,Month combinations:", total_combinations))

# 4. For each MAG_id, identify the number of Colony,Month combinations it's associated with
# We group by MAG_id and then count the number of distinct Colony,Month pairs
mag_counts <- df %>%
  group_by(MAG_id) %>%
  summarise(
    Genus = first(Genus), # Get the Genus value for the MAG_id
    count_associations = n_distinct(Colony, Month)
  )

# 5. Calculate the ratio
results_df <- mag_counts %>%
  mutate(ratio = count_associations / total_combinations) %>%
  ungroup()

# 6. Reorder and select the final columns for the output dataframe
results_df <- results_df %>%
  select(MAG_id, Genus, count_associations, ratio)

# Display the first few rows of the final dataframe
print("Final results dataframe (first 10 rows):")
print(head(results_df, 10))


# --- Second Analysis: Colony,Month,Rep combinations ---
# 3b. Identify all possible Colony,Month,Rep combinations
unique_combinations_cmr <- df %>%
  distinct(Colony, Month, Rep)

total_combinations_cmr <- nrow(unique_combinations_cmr)
print(paste("Total number of unique Colony,Month,Rep combinations:", total_combinations_cmr))

# 4b. For each MAG_id, identify the number of Colony,Month,Rep combinations it's associated with
mag_counts_cmr <- df %>%
  group_by(MAG_id) %>%
  summarise(
    Genus = first(Genus),
    count_associations_with_reps = n_distinct(Colony, Month, Rep)
  )

# 5b. Calculate the ratio
results_df_cmr <- mag_counts_cmr %>%
  mutate(ratio_with_reps = count_associations_with_reps / total_combinations_cmr) %>%
  ungroup()

# 6b. Select the final columns
results_df_cmr <- results_df_cmr %>%
  select(MAG_id, Genus, count_associations_with_reps, ratio_with_reps)

results_df$MAG_id <- factor(results_df$MAG_id, levels = results_df$MAG_id[order(results_df$ratio, decreasing = TRUE)])

# Create the plot
p <- ggplot(data = results_df, aes(x = MAG_id, y = ratio,color=Genus)) +
  geom_point(shape = 20, size = 3) + # Use geom_point for a scatter/dot plot with small dots
  theme_minimal() + # A clean theme for a better look
  labs(
    title = "MAG Association Ratio by MAG ID",
    y = "Ratio of Unique Colony,Month Combinations",
    x = NULL # Remove the x-axis label
  ) +
  theme(
    axis.text.y = element_blank(), # Remove the MAG_id axis labels
    axis.ticks.y = element_blank() # Remove the x-axis ticks
  )+
  coord_flip()
p


# 1. Calculate the weighted average ratio for each Genus
genus_weighted_means <- results_df %>%
  group_by(Genus) %>%
  summarise(
    weighted_mean_ratio = sum(ratio * count_associations, na.rm = TRUE) / sum(count_associations, na.rm = TRUE)
  ) %>%
  ungroup()

# 2. Order the Genus factor levels by the new weighted mean ratio in descending order
# This ensures ggplot plots the genera in the desired order
results_df$Genus <- factor(results_df$Genus, levels = genus_weighted_means$Genus[order(genus_weighted_means$weighted_mean_ratio, decreasing = FALSE)])

# 3. Create the plot with the newly ordered Genus axis
p <- ggplot(data = results_df, aes(x = Genus, y = ratio)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  theme_minimal() + # A clean theme
  labs(
    title = "MAG Association Ratio by Genus (Ordered by Weighted Average Ratio)",
    subtitle = "Each point represents a unique MAG_id",
    x = "Genus", # x and y labels are swapped for the flipped plot
    y = "Frequency of observations relative to total number of possible observations"
  ) +
  # Use coord_flip() to swap the axes
  coord_flip()

# Print the plot
print(p)
# --- 1. Calculate the Unique High-Coverage Sample Counts for Each MAG ---
# (Assuming 'df' is your main dataframe)
high_coverage_unique_counts <- df %>%
  # First, keep only the observations with high coverage
  filter(average_genome_wide_cov > 5) %>%
  # Group by MAG_id
  group_by(MAG_id) %>%
  # THE KEY CHANGE: Count the number of *distinct* sample_names for each MAG
  summarise(
    count_unique_high_cov_samples = n_distinct(sample_name)
  )

# --- 2. Join the New Count Column to the Results Dataframe ---
# (Assuming 'results_df_cmr' is your results dataframe)
results_df_cmr <- results_df_cmr %>%
  # Use a left_join to add the new counts
  left_join(high_coverage_unique_counts, by = "MAG_id") %>%
  # If a MAG had zero high-coverage samples, the join creates an NA. Replace NAs with 0.
  mutate(
    count_unique_high_cov_samples = replace_na(count_unique_high_cov_samples, 0)
  )

# --- 3. View the Updated Dataframe ---
print("Updated results_df_cmr with new unique count column:")
print(head(results_df_cmr))

