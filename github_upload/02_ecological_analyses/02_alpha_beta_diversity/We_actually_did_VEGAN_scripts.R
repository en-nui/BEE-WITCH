# Ensure you have the necessary packages installed:
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("readr")
# install.packages("vegan")
# install.packages("data.table")
# install.packages("ggplot2")

library(dplyr)
library(tidyr)
library(readr)
library(vegan)
library(data.table)
library(ggplot2) # Load ggplot2 for plotting

# --- 1. Define File Paths ---
abundance_metadata_file <- "/home/robinch/projects/BEE-WITCH/16_rMAG_mapping/colony_rMAG_mapping_with_metadata/concat_input/all_colonies_RETAINED_r_vMAGs.csv"
taxonomy_file <- "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_HOME_NO_TOUCH_20250327.csv"

message("Starting data preparation for vegan community analysis...")

# --- 2. Process Abundance Data ---
message("Processing abundance data...")
df_abundance_raw <- read_csv(abundance_metadata_file, show_col_types = FALSE)

# Select relevant columns and rename 'mean' to 'Abundance'
# Ensure sample_name and r_MAG_id are unambiguous character types
otu_df_long <- df_abundance_raw %>%
  select(sample_name, r_MAG_id, mean) %>%
  rename(Abundance = mean) %>%
  mutate(sample_name = as.character(sample_name),
         r_MAG_id = as.character(r_MAG_id))

# Explicitly convert 'Abundance' to integer and ensure it's a simple vector
otu_df_long$Abundance <- as.vector(as.integer(floor(otu_df_long$Abundance)))

# --- DIAGNOSTIC & HANDLING: Check for duplicate (sample_name, r_MAG_id) combinations ---
message("\n--- DIAGNOSTIC & HANDLING: Checking for duplicate (sample_name, r_MAG_id) combinations ---")
duplicates_check <- otu_df_long %>%
  dplyr::count(sample_name, r_MAG_id) %>%
  dplyr::filter(n > 1)

if(nrow(duplicates_check) > 0) {
  message("WARNING: Duplicate (sample_name, r_MAG_id) combinations found! Summing abundances.")
  message("Example duplicates (first 5):")
  print(head(duplicates_check, 5))
  message("Proceeding to SUM 'Abundance' for these duplicate combinations. ")
  message("If this is not the desired behavior, stop the script and review your data.")
  
  # Sum the Abundance for duplicate combinations
  otu_df_long <- otu_df_long %>%
    group_by(sample_name, r_MAG_id) %>%
    summarise(Abundance = sum(Abundance), .groups = 'drop') # Ensure groups are dropped after summarizing
  
  message("Duplicates summed. Data now unique for reshaping.")
} else {
  message("No duplicate (sample_name, r_MAG_id) combinations found. Data is unique for reshaping.")
}
message("-----------------------------------------------------------------------------------------")

# --- GENERATE SHORT IDs FOR r_MAG_id TO AVOID LONG COLUMN NAMES ---
message("\nGenerating short IDs for r_MAG_id to avoid 'variable names limited' error...")
unique_r_MAG_ids <- sort(unique(otu_df_long$r_MAG_id))
mag_id_mapping <- data.frame(
  r_MAG_id = unique_r_MAG_ids,
  short_MAG_id = paste0("MAG_", sprintf(paste0("%0", nchar(length(unique_r_MAG_ids)), "d"), 1:length(unique_r_MAG_ids))),
  stringsAsFactors = FALSE
)

# Apply short IDs to the abundance data for pivoting
otu_df_long_short_id <- otu_df_long %>%
  left_join(mag_id_mapping, by = "r_MAG_id") %>%
  select(sample_name, short_MAG_id, Abundance) # Select the new short ID for pivoting

message(paste("Generated", nrow(mag_id_mapping), "short MAG IDs. Example:", head(mag_id_mapping, 3)))
message("------------------------------------------------------------------")

# --- DATA.TABLE DCAST ALTERNATIVE FOR PIVOTING USING SHORT IDs ---
message("\nAttempting to pivot using data.table::dcast() with short IDs...")

# Convert to data.table for efficient pivoting
setDT(otu_df_long_short_id)

# Perform dcast (similar to pivot_wider/reshape)
abundance_dt <- dcast(otu_df_long_short_id, sample_name ~ short_MAG_id,
                      value.var = "Abundance",
                      fill = 0)

# Convert the data.table to a matrix for vegan.
abundance_matrix <- as.matrix(abundance_dt, rownames = "sample_name")

message("Pivoting with data.table::dcast() completed successfully using short IDs.")

message(paste("Abundance matrix dimensions:", dim(abundance_matrix)[1], "samples,", dim(abundance_matrix)[2], "features"))
message(paste("Example Abundance matrix (first 5 rows, first 5 cols):"))
print(abundance_matrix[1:min(5, nrow(abundance_matrix)), 1:min(5, ncol(abundance_matrix))])

# --- 3. Process Taxonomy Data ---
message("\nProcessing Taxonomy data...")
tax_df_raw <- read_csv(taxonomy_file, show_col_types = FALSE)

# Select desired taxonomic levels and ensure r_MAG_id is unique for row names
tax_df <- tax_df_raw %>%
  select(r_MAG_id, Family, Genus, Species) %>%
  # Handle potential duplicate r_MAG_ids by taking the first entry (or refine as needed)
  distinct(r_MAG_id, .keep_all = TRUE) %>%
  # Join with the short ID mapping to get the short_MAG_id
  left_join(mag_id_mapping, by = "r_MAG_id") %>%
  # Ensure we have a short_MAG_id for all taxa that might appear in the OTU table
  # (though filtering later will handle missing ones)
  filter(!is.na(short_MAG_id)) %>%
  # Set short_MAG_id as row names for the taxonomy matrix
  tibble::column_to_rownames(var = "short_MAG_id") %>%
  select(-r_MAG_id) # Remove the original long r_MAG_id column

# Convert to a data.frame (or keep as a matrix if preferred for taxonomy in vegan)
taxonomy_table <- as.data.frame(tax_df) # Storing as data.frame is often more flexible for taxonomy

message(paste("Taxonomy table dimensions:", dim(taxonomy_table)[1], "features,", dim(taxonomy_table)[2], "ranks"))
message(paste("Example Taxonomy table (first 5 rows):"))
print(taxonomy_table[1:min(5, nrow(taxonomy_table)), ])

# --- 4. Process Metadata ---
message("\nProcessing Metadata...")
# Start with the raw abundance/metadata file, select distinct sample metadata
metadata_df <- df_abundance_raw %>%
  select(sample_name, State, Month, Caste, Colony) %>%
  distinct(sample_name, .keep_all = TRUE) # Ensure one row per unique sample_name

# Derive Apiary
metadata_df <- metadata_df %>%
  mutate(Apiary = case_when(
    Colony %in% c(100, 101, 102, 103) ~ "PT",
    Colony %in% c(104, 105, 106, 107) ~ "OB",
    Colony %in% c(555, 777, 888, 999) ~ "HB",
    TRUE ~ NA_character_ # Assign NA if Colony does not match any rule
  ))

# Derive Season
metadata_df <- metadata_df %>%
  mutate(Season = case_when(
    Month %in% c("May", "June") ~ "Spring",
    Month %in% c("July", "August", "September") ~ "Summer",
    Month %in% c("October", "November") ~ "Fall",
    Month %in% c("January", "February") ~ "Winter",
    TRUE ~ NA_character_ # Assign NA if Month does not match any rule
  )) %>%
  tibble::column_to_rownames(var = "sample_name")

# Convert relevant columns to factors if they are to be used as grouping variables in vegan
metadata_df <- metadata_df %>%
  mutate(across(c(State, Month, Caste, Apiary, Season), as.factor))

# Define the order of months explicitly for consistent plotting on X-axis
ordered_months_for_factor <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")

# Convert Month to an ordered factor in metadata_df before adding Time (crucial for plotting order)
metadata_df <- metadata_df %>%
  mutate(Month = factor(Month, levels = ordered_months_for_factor, ordered = TRUE))

sample_metadata <- metadata_df

message(paste("Sample metadata table dimensions:", dim(sample_metadata)[1], "samples,", dim(sample_metadata)[2], "variables"))
message(paste("Example Sample metadata (first 5 rows):"))
print(head(sample_metadata, 5))


# --- 5. Align and Filter Tables ---
message("\nAligning and filtering tables...")

# Identify features (short_MAG_ids) present in the abundance matrix
features_in_abundance <- colnames(abundance_matrix)
samples_in_abundance <- rownames(abundance_matrix)

# Filter taxonomy table to include only those features present in abundance matrix
taxonomy_table_filtered <- taxonomy_table[rownames(taxonomy_table) %in% features_in_abundance, , drop = FALSE]

# Order taxonomy table rows to match abundance matrix columns (important for consistency)
taxonomy_table_filtered <- taxonomy_table_filtered[match(features_in_abundance, rownames(taxonomy_table_filtered)), , drop = FALSE]

# Filter sample metadata to include only samples present in abundance matrix
sample_metadata_filtered <- sample_metadata[rownames(sample_metadata) %in% samples_in_abundance, , drop = FALSE]

# Order sample metadata rows to match abundance matrix rows (important for consistency)
sample_metadata_filtered <- sample_metadata_filtered[match(samples_in_abundance, rownames(sample_metadata_filtered)), , drop = FALSE]

# Final check for matching dimensions
if (!all(rownames(abundance_matrix) == rownames(sample_metadata_filtered))) {
  stop("Sample names in abundance matrix and sample metadata do not match after filtering/ordering!")
}
if (!all(colnames(abundance_matrix) == rownames(taxonomy_table_filtered))) {
  stop("Feature names in abundance matrix and taxonomy table do not match after filtering/ordering!")
}

# Assign filtered objects back to original names for clarity
abundance_matrix <- abundance_matrix
taxonomy_table <- taxonomy_table_filtered
sample_metadata <- sample_metadata_filtered


message("Alignment complete. Data is consistent across all tables.")

# --- 6. Add Cumulative Time to Metadata ---
message("\nAdding cumulative time from month to metadata...")

# Define the function to add cumulative time
add_cumulative_time_from_month <- function(data_df) {
  # Define the base time values for each month in order
  # These are interpreted as increments from the previous month's start for cumulative sum
  month_base_times <- c(
    "May" = 0, # Starting point
    "June" = 36,
    "July" = 31,
    "August" = 34,
    "September" = 28,
    "October" = 22,
    "November" = 42,
    "January" = 47,
    "February" = 31
  )
  
  # Calculate the cumulative time values based on the ordered months
  ordered_months <- names(month_base_times)
  # Ensure the internal cumulative_times calculation uses the correct order
  cumulative_times <- cumsum(month_base_times)
  
  # Create a mapping between Month and Cumulative Time
  month_to_cumulative_time_map <- setNames(cumulative_times, ordered_months)
  
  # Add the "Time" column based on the "Month" column in the metadata
  # Ensure Month is treated as a character for mapping even if it's a factor
  data_df$Time <- month_to_cumulative_time_map[as.character(data_df$Month)]
  
  # Check if any months in the data were not found in the mapping
  if (any(is.na(data_df$Time))) {
    warning("Warning: Some months in the 'Month' column were not found in the mapping and have been assigned NA in the 'Time' column.")
  }
  
  return(data_df)
}

# Apply the function to your sample_metadata (which now has Month as ordered factor)
sample_metadata <- add_cumulative_time_from_month(sample_metadata)

message("Cumulative 'Time' column added to sample_metadata.")
message("Example sample_metadata with Time (first 5 rows):")
print(head(sample_metadata, 5))


# --- 7. Calculate Richness and Diversity Indices ---
message("\nCalculating vMAG richness, Shannon, and Inverse Simpson diversity indices...")

# Calculate richness (number of non-zero MAGs) for each sample
# Create a data frame where 'sample_name' is a regular column, not row names yet
sample_richness_df <- data.frame(
  sample_name = rownames(abundance_matrix),
  richness = specnumber(abundance_matrix),
  stringsAsFactors = FALSE
)

# Calculate Shannon and Inverse Simpson diversity
diversity_indices_df <- data.frame(
  sample_name = rownames(abundance_matrix),
  shannon = diversity(abundance_matrix, "shannon"),
  simpson_inverse = diversity(abundance_matrix, "invsimpson"), # Using inverse Simpson for intuitive interpretation (higher value = more diverse)
  stringsAsFactors = FALSE
)

# Combine richness, diversity, and sample_metadata into one comprehensive data frame
combined_diversity_data <- sample_richness_df %>%
  left_join(diversity_indices_df, by = "sample_name") %>%
  left_join(sample_metadata %>% tibble::rownames_to_column(var = "sample_name"), by = "sample_name") %>%
  # Filter out samples where Time is NA, as they cannot be plotted on the time axis
  filter(!is.na(Time))

message(paste("Combined diversity data prepared. Total samples with valid time:", nrow(combined_diversity_data)))
message("Example combined diversity data (first 5 rows):")
print(head(combined_diversity_data, 5))


# --- 8. Generate Diversity Plots ---
message("\nGenerating diversity plots (Richness, Shannon, Inverse Simpson)...")

# Define a common plotting function to reduce redundancy
# 'y_limits' and 'y_breaks' are optional for custom axis scales
plot_diversity <- function(data, metric_col, metric_name, y_limits = NULL, y_breaks = NULL) {
  
  plot_list <- list()
  
  # --- Plot by Apiary (all samples) ---
  p_apiary <- ggplot(data, aes_string(x = "Time", y = metric_col, color = "Apiary")) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = paste0(metric_name, " Over Time by Apiary"),
      x = "Cumulative Time (Days)",
      y = metric_name,
      color = "Apiary"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(y_limits)) {
    p_apiary <- p_apiary + scale_y_continuous(limits = y_limits, breaks = y_breaks)
  }
  plot_list[[paste0("apiary_", tolower(gsub(" ", "_", metric_name)))]] <- p_apiary
  
  # --- Plot by Colony (Worker Caste Only) - Scatterplot with LOESS Line ---
  current_data_colony_worker <- data %>% filter(Caste == "Worker")
  
  if (nrow(current_data_colony_worker) == 0) {
    message(paste0("No 'Worker' caste samples found for plotting ", metric_name, " by Colony. Skipping plot."))
  } else {
    p_colony_scatterplot <- ggplot(current_data_colony_worker, aes_string(x = "Time", y = metric_col, color = "Colony")) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = FALSE, color = "darkblue") + # LOESS smooth line for trend
      facet_wrap(~ Colony, scales = "fixed", ncol = 4) + # Fixed scales for consistent Y-axis
      labs(
        title = paste0(metric_name, " Over Time by Colony (Worker Caste)"),
        x = "Cumulative Time (Days)",
        y = metric_name,
        color = "Colony"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.text = element_text(face = "bold"),
            legend.position = "none")
    
    if (!is.null(y_limits)) {
      p_colony_scatterplot <- p_colony_scatterplot + scale_y_continuous(limits = y_limits, breaks = y_breaks)
    }
    plot_list[[paste0("colony_worker_scatterplot_", tolower(gsub(" ", "_", metric_name)))]] <- p_colony_scatterplot
  }
  
  return(plot_list)
}


# --- Plot Richness ---
message("\n--- Generating Richness Plots ---")
richness_plots <- plot_diversity(
  data = combined_diversity_data,
  metric_col = "richness",
  metric_name = "vMAG Richness (Number of Unique vMAGs)",
  y_limits = c(0, 100),
  y_breaks = seq(0, 100, by = 20)
)
print(richness_plots$`apiary_vmag_richness_(number_of_unique_vmags)`)
# Only the scatterplot version for colony now
print(richness_plots$`colony_worker_scatterplot_vmag_richness_(number_of_unique_vmags)`)


# --- Plot Richness ---
message("\n--- Generating Richness Plots ---")

# --- Step 1: Define the correct column name for time points ---
time_column_name <- "Time"

# --- Step 2: Basic validation of the 'Time' column ---
if (!time_column_name %in% colnames(combined_diversity_data)) {
  stop(paste("Error: Column '", time_column_name, "' not found in combined_diversity_data. Please check column names.", sep = ""))
}
if (all(is.na(combined_diversity_data[[time_column_name]]))) {
  stop(paste("Error: Column '", time_column_name, "' contains only missing values (NA). Cannot plot time points.", sep = ""))
}

# --- Step 3: Generate the base richness plot ---
# IMPORTANT: Ensure your plot_diversity function uses the 'Time' column for the x-axis.
# If your plot_diversity function has an 'x_col' or similar argument, you'd pass it like this:
richness_base_plot <- plot_diversity(
  data = combined_diversity_data,
  metric_col = "richness",
  # x_col = time_column_name, # Uncomment and use if your plot_diversity function supports an x_col argument
  metric_name = "vMAG Richness (Number of Unique vMAGs)",
  y_limits = c(0, 100),
  y_breaks = seq(0, 100, by = 20)
)

# --- Step 4: Prepare data for vertical lines ---
x_intercepts_all_timepoints <- unique(combined_diversity_data[[time_column_name]])
x_intercepts_all_timepoints <- sort(x_intercepts_all_timepoints)

# Convert to numeric if 'Time' is a Date/POSIXct type, otherwise use as-is (assuming it's numeric)
if (inherits(combined_diversity_data[[time_column_name]], c("Date", "POSIXct"))) {
  message("Converting 'Time' column to numeric for xintercepts.")
  x_intercepts_all_timepoints_for_plot <- as.numeric(x_intercepts_all_timepoints)
} else if (is.numeric(combined_diversity_data[[time_column_name]])) {
  message("'Time' column is already numeric. Using as-is for xintercepts.")
  x_intercepts_all_timepoints_for_plot <- x_intercepts_all_timepoints
} else {
  # This block will handle cases where 'Time' might be a factor or character that *should* be numeric.
  # It attempts a conversion but will warn if it encounters issues.
  message(paste("Warning: '", time_column_name, "' is neither Date/POSIXct nor numeric. Attempting conversion to numeric.", sep = ""))
  x_intercepts_all_timepoints_for_plot <- suppressWarnings(as.numeric(as.character(x_intercepts_all_timepoints)))
  x_intercepts_all_timepoints_for_plot <- na.omit(x_intercepts_all_timepoints_for_plot) # Remove NAs from conversion
}

message(paste("Number of unique x-intercepts for '", time_column_name, "':", sep = ""))
print(length(x_intercepts_all_timepoints_for_plot))
if (length(x_intercepts_all_timepoints_for_plot) > 0) {
  message("First few x-intercepts:")
  print(head(x_intercepts_all_timepoints_for_plot))
}

# --- Step 5: Add faceting and vertical lines to the plot ---
richness_plots_faceted_and_vlines <- richness_base_plot$`apiary_vmag_richness_(number_of_unique_vmags)` +
  facet_wrap(~ Apiary, scales = "free_x") + # 'free_x' allows each facet to have its own x-axis range if needed
  geom_vline(xintercept = x_intercepts_all_timepoints_for_plot,
             linetype = "dotted",
             color = "black",
             linewidth = 0.5)

print(richness_plots_faceted_and_vlines)




# --- Plot Shannon Diversity ---
message("\n--- Generating Shannon Diversity Plots ---")
max_shannon <- max(combined_diversity_data$shannon, na.rm = TRUE)
shannon_plots <- plot_diversity(
  data = combined_diversity_data,
  metric_col = "shannon",
  metric_name = "Shannon Diversity Index",
  y_limits = c(0, ceiling(max_shannon * 1.1)), # Round up for a clean upper limit
  y_breaks = seq(0, ceiling(max_shannon * 1.1), by = 1)
)
print(shannon_plots$apiary_shannon_diversity_index)
print(shannon_plots$colony_worker_scatterplot_shannon_diversity_index)


# --- Plot Inverse Simpson Diversity ---
message("\n--- Generating Inverse Simpson Diversity Plots ---")
max_simpson_inverse <- max(combined_diversity_data$simpson_inverse, na.rm = TRUE)
simpson_plots <- plot_diversity(
  data = combined_diversity_data,
  metric_col = "simpson_inverse",
  metric_name = "Inverse Simpson Diversity Index",
  y_limits = c(0, ceiling(max_simpson_inverse * 1.1)), # Round up for a clean upper limit
  y_breaks = seq(0, ceiling(max_simpson_inverse * 1.1), by = max(1, floor(ceiling(max_simpson_inverse * 1.1) / 5))) # Aim for ~5 breaks
)
print(simpson_plots$apiary_inverse_simpson_diversity_index)
print(simpson_plots$colony_worker_scatterplot_inverse_simpson_diversity_index)


message("\nAll diversity plots generated!")


# ==============================================================================
# --- Bray-Curtis Dissimilarity Analysis ---
# ==============================================================================
message("\n--- Starting Bray-Curtis Dissimilarity Analysis ---")

# Ensure abundance_matrix and sample_metadata are aligned before proceeding
# (This was handled in step 5, but good to be mindful)
abundance_matrix_aligned <- abundance_matrix[rownames(sample_metadata), ]

# Calculate Bray-Curtis dissimilarity matrix for all samples
bray_dist_all <- vegdist(abundance_matrix_aligned, method = "bray")


# --- Part 1: Compare Bray-Curtis Dissimilarity Between Colony and Apiary ---
message("\n--- Starting Bray-Curtis Dissimilarity & Similarity Analysis ---")

# Ensure abundance_matrix and sample_metadata are aligned before proceeding
# (This was handled in step 5, but good to be mindful)
abundance_matrix_aligned <- abundance_matrix[rownames(sample_metadata), ]

# Calculate Bray-Curtis dissimilarity matrix for all samples
bray_dist_all <- vegdist(abundance_matrix_aligned, method = "bray")


# --- Part 1: Compare Bray-Curtis Dissimilarity Between Colony and Apiary (PERMANOVA & NMDS) ---
message("\n--- Part 1: Community Structure Comparison (PERMANOVA & NMDS) ---")
message("Performing PERMANOVA (adonis2) to compare community structure by Apiary, Colony, and Caste...")

# Ensure metadata factors are correctly set for adonis2
sample_metadata_for_adonis <- sample_metadata %>%
  mutate(
    Apiary = factor(Apiary),
    Colony = factor(Colony),
    Caste = factor(Caste)
  )

# Perform PERMANOVA using adonis2
# Using sequential sums of squares (Type I) which means order matters.
# Apiary first, then Colony (nested within Apiary implicitly if Colony IDs are unique across apiaries, otherwise needs specific nesting syntax if Colony IDs repeat).
# Caste is also a major grouping. Time as a continuous covariate.
set.seed(123) # For reproducibility of permutations
permanova_result <- adonis2(bray_dist_all ~ Apiary + Colony + Caste + Time, 
                            data = sample_metadata_for_adonis,
                            permutations = 999,
                            by = "term") # "by = term" gives results for each term adjusted for terms before it.

message("\n--- PERMANOVA Results (Bray-Curtis Dissimilarity) ---")
print(permanova_result)


# --- Part 2: Bray-Curtis Similarity over Time (Pairwise Comparisons & Mantel Test) ---
message("\n--- Part 2: Bray-Curtis Similarity over Time (Pairwise & Mantel Test) ---")

# Filter to Worker caste for this analysis (as per previous colony plots and user request)
worker_data_for_bray_time <- combined_diversity_data %>% filter(Caste == "Worker")
worker_abundance_matrix_time <- abundance_matrix[rownames(abundance_matrix) %in% worker_data_for_bray_time$sample_name, ]
# Ensure samples are in the same order
worker_abundance_matrix_time <- worker_abundance_matrix_time[worker_data_for_bray_time$sample_name, ]

# Calculate Bray-Curtis distance for all worker samples once
full_worker_bray_dist_time <- vegdist(worker_abundance_matrix_time, method = "bray")
# Convert to similarity for plotting and Mantel test
full_worker_bray_sim_time <- 1 - full_worker_bray_dist_time

# --- ADD THIS NEW LINE ---
# Convert the similarity 'dist' object to a full matrix for direct indexing
full_worker_bray_sim_matrix <- as.matrix(full_worker_bray_sim_time)
# --- END NEW LINE ---


# Create a Time Difference Matrix for all worker samples
time_values_worker <- worker_data_for_bray_time$Time[match(rownames(worker_abundance_matrix_time), worker_data_for_bray_time$sample_name)]
time_diff_matrix <- as.matrix(dist(time_values_worker, method = "euclidean"))
colnames(time_diff_matrix) <- rownames(worker_abundance_matrix_time)
rownames(time_diff_matrix) <- rownames(worker_abundance_matrix_time)


# --- Perform Mantel Tests ---
message("\n--- Performing Mantel Tests for Bray-Curtis Similarity vs. Time Difference ---")

# Overall Mantel test for all worker samples
if (length(time_values_worker) > 1) {
  message("\nOverall Mantel test (All Worker Samples):")
  mantel_overall <- mantel(full_worker_bray_sim_time, time_diff_matrix, permutations = 999, method = "spearman")
  print(mantel_overall)
  message(paste("Interpretation: A positive r means higher similarity with longer time difference, a negative r means lower similarity (higher dissimilarity) with longer time difference."))
} else {
  message("Not enough worker samples for overall Mantel test.")
}


# Mantel tests per Apiary (Worker Caste only)
message("\nMantel tests per Apiary (Worker Caste):")
unique_apiaries <- unique(worker_data_for_bray_time$Apiary)
mantel_apiary_results <- list()
for (apiary in unique_apiaries) {
  apiary_samples <- worker_data_for_bray_time %>% filter(Apiary == apiary) %>% pull(sample_name)
  if (length(apiary_samples) > 1) {
    bray_sim_sub <- as.dist(as.matrix(full_worker_bray_sim_time)[apiary_samples, apiary_samples])
    time_diff_sub <- as.dist(time_diff_matrix[apiary_samples, apiary_samples])
    
    if(length(bray_sim_sub) > 0 && length(time_diff_sub) > 0) { # Ensure there are actual distances to test
      mantel_res <- mantel(bray_sim_sub, time_diff_sub, permutations = 999, method = "pearson")
      mantel_apiary_results[[apiary]] <- mantel_res
      message(paste0("Apiary ", apiary, ": r = ", round(mantel_res$statistic, 3), ", p = ", round(mantel_res$signif, 3)))
    } else {
      message(paste0("Apiary ", apiary, ": Not enough unique pairwise comparisons for Mantel test."))
    }
  } else {
    message(paste0("Apiary ", apiary, ": Not enough worker samples for Mantel test."))
  }
}

# Mantel tests per Colony (Worker Caste only)
message("\nMantel tests per Colony (Worker Caste):")
unique_colonies <- unique(worker_data_for_bray_time$Colony)
mantel_colony_results <- list()
for (colony in unique_colonies) {
  colony_samples <- worker_data_for_bray_time %>% filter(Colony == colony) %>% pull(sample_name)
  if (length(colony_samples) > 1) {
    bray_sim_sub <- as.dist(as.matrix(full_worker_bray_sim_time)[colony_samples, colony_samples])
    time_diff_sub <- as.dist(time_diff_matrix[colony_samples, colony_samples])
    
    if(length(bray_sim_sub) > 0 && length(time_diff_sub) > 0) {
      mantel_res <- mantel(bray_sim_sub, time_diff_sub, permutations = 999, method = "pearson")
      mantel_colony_results[[as.character(colony)]] <- mantel_res
      message(paste0("Colony ", colony, ": r = ", round(mantel_res$statistic, 3), ", p = ", round(mantel_res$signif, 3)))
    } else {
      message(paste0("Colony ", colony, ": Not enough unique pairwise comparisons for Mantel test."))
    }
  } else {
    message(paste0("Colony ", colony, ": Not enough worker samples for Mantel test."))
  }
}


# Prepare dataframe for pairwise results (Now with Bray-Curtis Similarity)
pairwise_dissimilarity_results <- data.frame(
  sample1 = character(),
  sample2 = character(),
  time_diff = numeric(),
  bray_curtis_dissimilarity = numeric(), # Keep this for reference
  bray_curtis_similarity = numeric(), # New column for similarity
  group_type = character(), # "Apiary" or "Colony"
  group_name = character(),
  stringsAsFactors = FALSE
)

# Function to process pairs and extract data (modified to directly use similarity matrix)
process_pairs_for_time_and_similarity <- function(data_subset, group_name, group_type, bray_similarity_matrix, time_df) { # Changed argument name
  pairs_data <- data.frame(
    sample1 = character(),
    sample2 = character(),
    time_diff = numeric(),
    bray_curtis_dissimilarity = numeric(),
    bray_curtis_similarity = numeric(),
    group_type = character(),
    group_name = character(),
    stringsAsFactors = FALSE
  )
  
  samples_in_group <- data_subset$sample_name
  
  if (length(samples_in_group) < 2) {
    return(pairs_data) # Cannot form pairs with less than 2 samples
  }
  
  for (i in 1:(length(samples_in_group) - 1)) {
    for (j in (i + 1):length(samples_in_group)) {
      s1 <- samples_in_group[i]
      s2 <- samples_in_group[j]
      
      bc_sim <- bray_similarity_matrix[s1, s2] # Directly use similarity from the passed matrix
      bc_diss <- 1 - bc_sim # Calculate dissimilarity for logging if needed
      
      time1 <- time_df$Time[time_df$sample_name == s1]
      time2 <- time_df$Time[time_df$sample_name == s2]
      time_diff <- abs(time1 - time2)
      
      pairs_data <- rbind(pairs_data, data.frame(
        sample1 = s1,
        sample2 = s2,
        time_diff = time_diff,
        bray_curtis_dissimilarity = bc_diss,
        bray_curtis_similarity = bc_sim,
        group_type = group_type,
        group_name = group_name,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(pairs_data)
}


# --- Process by Apiary (Worker Caste only) ---
message("Recalculating pairwise Bray-Curtis results with similarity...")
unique_apiaries <- unique(worker_data_for_bray_time$Apiary)
for (apiary in unique_apiaries) {
  apiary_data <- worker_data_for_bray_time %>% filter(Apiary == apiary)
  if (nrow(apiary_data) > 1) {
    pairwise_dissimilarity_results <- rbind(pairwise_dissimilarity_results, 
                                            process_pairs_for_time_and_similarity(apiary_data, apiary, "Apiary", full_worker_bray_sim_matrix, worker_data_for_bray_time)) # Updated here
  }
}

# --- Process by Colony (Worker Caste only) ---
unique_colonies <- unique(worker_data_for_bray_time$Colony)
for (colony in unique_colonies) {
  colony_data <- worker_data_for_bray_time %>% filter(Colony == colony)
  if (nrow(colony_data) > 1) {
    pairwise_dissimilarity_results <- rbind(pairwise_dissimilarity_results,
                                            process_pairs_for_time_and_similarity(colony_data, as.character(colony), "Colony", full_worker_bray_sim_matrix, worker_data_for_bray_time)) # Updated here
  }
}

message("Pairwise Bray-Curtis similarity results prepared.")
message("Example pairwise results (first 5 rows):")
print(head(pairwise_dissimilarity_results, 5))

# --- Plotting Part 2 (Now with Bray-Curtis Similarity on Y-axis) ---
message("\nPlotting pairwise Bray-Curtis Similarity over time...")

# Plot for Apiary
plot_bray_sim_apiary_time <- ggplot(
  pairwise_dissimilarity_results %>% filter(group_type == "Apiary"),
  aes(x = time_diff, y = bray_curtis_similarity, color = group_name)
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE) + # Show confidence interval
  facet_wrap(~ group_name, scales = "free_x") + # Free x-scales might be useful for different time ranges across apiaries
  labs(
    title = "Bray-Curtis Similarity vs. Time Difference within Apiaries (Worker Caste)",
    x = "Absolute Time Difference (Days)",
    y = "Bray-Curtis Similarity",
    color = "Apiary"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold"))

print(plot_bray_sim_apiary_time)

# Plot for Colony
plot_bray_sim_colony_time <- ggplot(
  pairwise_dissimilarity_results %>% filter(group_type == "Colony"),
  aes(x = time_diff, y = bray_curtis_similarity, color = group_name)
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) + # Show confidence interval
  facet_wrap(~ group_name, scales = "free_x", ncol = 4) + # Free x-scales, fixed cols for facet layout
  labs(
    title = "Bray-Curtis Similarity vs. Time Difference within Colonies (Worker Caste)",
    x = "Absolute Time Difference (Days)",
    y = "Bray-Curtis Similarity",
    color = "Colony"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none") # Hide legend as color is also facet label

print(plot_bray_sim_colony_time)


message("\nBray-Curtis dissimilarity and similarity analysis complete. All results and plots generated!")


# (Your existing code up to the point of creating worker_data_for_bray_time)

message("\n--- Diagnosing Mantel Test p=1 issue ---")
message("Number of worker samples (n_distinct) per Apiary in the filtered dataset:")
worker_data_for_bray_time %>%
  group_by(Apiary) %>%
  summarise(n_samples = n_distinct(sample_name)) %>%
  print()

message("\nNumber of worker samples (n_distinct) per Colony in the filtered dataset:")
worker_data_for_bray_time %>%
  group_by(Colony) %>%
  summarise(n_samples = n_distinct(sample_name)) %>%
  print()

# Identify groups that might be problematic (i.e., too few samples for robust permutation)
problem_apiaries <- worker_data_for_bray_time %>%
  group_by(Apiary) %>%
  summarise(n_samples = n_distinct(sample_name)) %>%
  filter(n_samples < 4) %>% pull(Apiary)

problem_colonies <- worker_data_for_bray_time %>%
  group_by(Colony) %>%
  summarise(n_samples = n_distinct(sample_name)) %>%
  filter(n_samples < 4) %>% pull(Colony)

if (length(problem_apiaries) > 0) {
  message(paste("Warning: Apiaries with fewer than 4 worker samples (p=1 or NaN expected for Mantel):", paste(problem_apiaries, collapse = ", ")))
}
if (length(problem_colonies) > 0) {
  message(paste("Warning: Colonies with fewer than 4 worker samples (p=1 or NaN expected for Mantel):", paste(problem_colonies, collapse = ", ")))
}














