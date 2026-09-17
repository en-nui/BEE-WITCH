# Load necessary libraries
library(dplyr)
library(ggplot2)

# --- Load the concatenated dataframe ---
df <- read.csv('/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output/all_filtered_mmag_contigs_concatenated.csv')

# First, convert all column names to lowercase
df <- df %>%
  rename_with(tolower)

# Filter data for df_one (using lowercase column names and values)
df <- df %>%
  filter(caste %in% c("Worker"))

# Base directory for figures (now organized by apiary)
figure_output_base_dir <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/figures/figures_by_apiary"

# Define month order (using lowercase column names)
full_month_order <- c("May", "June", "July", "August", "September", "October", "November","January","February")
df$month <- factor(df$month, levels = full_month_order)

castes_to_plot <- c("Worker", "Larvae") # Specify castes for plotting (using lowercase values)

# --- Filter for type == "mmag" ---
df <- df %>%
  filter(type == "mmag")

# Add apiary column (using lowercase colony)
df <- df %>%
  mutate(apiary = case_when(
    colony > 500 ~ "HB",
    colony <= 103 ~ "PT",
    TRUE ~ "OB"
  ))

# --- IMPORTANT: Handle unique average_genome_wide_cov per mag_id ---
# To ensure each MAG_id's coverage is counted only once, we group by unique identifiers
# and take the first 'average_genome_wide_cov' value.
# DO NOT create final_genus_category here yet, as it depends on global top_genera.
df_processed_initial <- df %>%
  group_by(mag_id, month, colony, rep, caste, genus, type, apiary, contig_presence_ratio) %>%
  summarise(
    avg_cov = first(average_genome_wide_cov), # Use average_genome_wide_cov as the input for relative abundance
    .groups = 'drop'
  )

cat("\n--- Script Complete (Initial Data Processing) ---\n")

# --- Define a consistent color palette for genera (Gruvbox-inspired) ---
gruvbox_colors <- c(
  "#a89984", "#fabd2f", "#fe8019", "#d7ba7d", "#cc241d", "#b16286",
  "#98971a", "#689d6a", "#458588", "#83a598", "#d3869b", "#8ec07c"
)

# --- GLOBAL: Identify Genera present in EVERY Apiary-Month combination (for consistent coloring) ---
# This will be the set of genera that get unique colors across ALL plots.
total_unique_apiary_month_combinations <- df_processed_initial %>%
  distinct(apiary, month) %>%
  nrow()

genera_in_all_apiary_months_globally_raw <- df_processed_initial %>%
  filter(genus != "unknown_mmag", genus != "no match") %>% # Exclude both "unknown_mmag" AND "no match"
  group_by(genus) %>%
  summarise(apiary_month_count = n_distinct(paste(apiary, month)), .groups = 'drop') %>%
  filter(apiary_month_count == total_unique_apiary_month_combinations) %>%
  pull(genus)

cat(paste("\n--- Genera present in EVERY Apiary-Month combination (Globally for coloring):",
          ifelse(length(genera_in_all_apiary_months_globally_raw) > 0,
                 paste(genera_in_all_apiary_months_globally_raw, collapse = ", "),
                 "None"), "---\n"))

# Calculate overall abundance for these globally consistent genera to pick the top 12
overall_abundance_for_globally_consistent <- df_processed_initial %>%
  filter(genus %in% genera_in_all_apiary_months_globally_raw) %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(avg_cov, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(total_abundance))

# Define 'top_genera' for coloring: top 12 most abundant among those consistently present globally
top_genera <- head(overall_abundance_for_globally_consistent$genus, 12)
top_n <- length(top_genera) # Actual number of top genera selected for specific colors (will be <= 12)

# --- GLOBAL: Define the color palette for all plots ---
top_genus_colors <- list()
if (top_n > 0) {
  colors_to_use <- head(gruvbox_colors, top_n)
  names(colors_to_use) <- top_genera
  top_genus_colors <- as.list(colors_to_use)
}
# Add color for "MAG below length threshold (0.5)"
top_genus_colors[["MAG below length threshold (0.5)"]] <- "#888888" # Dark grey
# Add color for "unknown_mmag"
top_genus_colors[["unknown_mmag"]] <- "#888888" # A slightly lighter grey for unknown
# Add color for "Other Genera" (must be last to catch everything else)
top_genus_colors[["Other Genera"]] <- "#888888"


# --- NOW: Add the final_genus_category to df_processed based on global top_genera and new rules ---
df_processed <- df_processed_initial %>%
  mutate(
    final_genus_category = case_when(
      contig_presence_ratio < 0.5 ~ "MAG below length threshold (0.5)",
      genus == "no match" ~ "Other Genera", # Explicitly handle "no match"
      genus %in% top_genera ~ genus, # Assign specific genus name if it's one of the global top
      is.na(genus) | genus == "" ~ "unknown_mmag", # Handle explicit unknown
      TRUE ~ "Other Genera" # Catch all remaining genera
    )
  )

cat("\n--- Script Complete (Global Setup) ---\n")

# --- Get unique Apiary values for the outer loop ---
unique_apiaries <- unique(df_processed$apiary)

# --- Start Apiary Loop ---
for (current_apiary in unique_apiaries) {
  cat(paste("\n--- Processing Apiary Directory:", current_apiary, "---\n"))
  
  # --- Create apiary-specific directory within the base figure directory ---
  apiary_figure_output_dir <- file.path(figure_output_base_dir, paste0("apiary_", current_apiary))
  dir.create(apiary_figure_output_dir, recursive = TRUE, showWarnings = FALSE) # Create apiary dir
  
  # --- Loop through each Caste for the current Apiary ---
  for (current_caste in castes_to_plot) {
    cat(paste("  Processing Caste:", current_caste, "---\n"))
    
    # --- Filter data for the current Apiary and Caste, aggregating across all Reps ---
    df_apiary_caste <- df_processed %>%
      filter(apiary == current_apiary, caste == current_caste, type == "mmag")
    
    # --- Skip to next Caste if no data for this Apiary-Caste combination ---
    if(nrow(df_apiary_caste) == 0) {
      cat(paste("  No data for Apiary", current_apiary, "Caste", current_caste, ". Skipping to next caste.\n"))
      next
    }
    
    # --- Calculate total abundance per month for this Apiary-Caste ---
    month_totals_apiary_caste <- df_apiary_caste %>%
      group_by(month) %>%
      summarise(total_month_abundance = sum(avg_cov)) # Using avg_cov
    
    # --- Calculate relative abundance per month and final_genus_category for this Apiary-Caste ---
    df_rel_abundance_apiary_caste <- df_apiary_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(genus_month_abundance = sum(avg_cov)) %>% # Using avg_cov
      left_join(month_totals_apiary_caste, by = "month") %>%
      mutate(relative_abundance = genus_month_abundance / total_month_abundance)
    
    # --- Re-aggregate data by Month and final_genus_category ---
    # This step is still needed to sum up relative abundances for categories like "Other Genera"
    df_rel_abundance_apiary_caste_aggregated <- df_rel_abundance_apiary_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(relative_abundance = sum(relative_abundance))
    
    # --- Convert Month to ordered factor for plotting ---
    df_rel_abundance_apiary_caste_aggregated$month <- factor(df_rel_abundance_apiary_caste_aggregated$month, levels = full_month_order)
    
    # --- Dynamic Month Order for the current Apiary-Caste ---
    current_apiary_caste_months_in_data <- unique(as.character(df_rel_abundance_apiary_caste_aggregated$month))
    current_apiary_caste_month_order <- full_month_order[full_month_order %in% current_apiary_caste_months_in_data]
    
    # --- Create Area Plot for the current Apiary-Caste ---
    apiary_caste_area_plot <- ggplot(df_rel_abundance_apiary_caste_aggregated,
                                     aes(x = month, y = relative_abundance, fill = final_genus_category, group = final_genus_category, color = final_genus_category)) +
      geom_area(position = "stack") +
      scale_x_discrete(limits = current_apiary_caste_month_order) +
      scale_fill_manual(values = top_genus_colors) + # Use the GLOBAL color palette
      scale_color_manual(values = top_genus_colors) + # Ensure consistent colors for outlines too
      labs(title = paste("Microbial Genus Relative Abundance - Apiary:", current_apiary, "- Caste:", current_caste),
           x = "Month",
           y = "Relative Abundance",
           fill = "Bacterial Genus",
           color = "Bacterial Genus") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5), legend.position="right")
    
    print(apiary_caste_area_plot)
    
    # --- Save the Apiary-Caste Area plot to a file (optional) ---
    if (!dir.exists(apiary_figure_output_dir)) {
      dir.create(apiary_figure_output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    apiary_caste_area_plot_file_name <- paste0("apiary_", current_apiary, "_caste_", current_caste, "_genus_area_plot_global_consistent_colors_aggregated_reps.png") # Updated filename
    apiary_caste_area_plot_file_path <- file.path(apiary_figure_output_dir, apiary_caste_area_plot_file_name)
    ggsave(apiary_caste_area_plot_file_path, plot = apiary_caste_area_plot, width = 10, height = 6, units = "in")
    cat(paste("  Apiary & Caste mMAG Genus Area plot (Global Consistent Colors, Aggregated Reps) saved for Apiary", current_apiary, "Caste", current_caste, "to:", apiary_caste_area_plot_file_path, "\n"))
  } # --- END Caste Loop ---
} # --- END Apiary Loop ---

cat("\n--- Script Complete ---\n")

###################### COLONIES ####################3

df <- df %>%
  rename_with(tolower)

# Filter data for df_one (using lowercase column names and values)
df <- df %>%
  filter(caste %in% c("Worker", "Larvae"))

# Base directory for figures (now organized by colony)
figure_output_base_dir <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/figures/figures_by_colony" # Changed to figures_by_colony

# Define month order (using lowercase column names)
full_month_order <- c("May", "June", "July", "August", "September", "October", "November","January","February")
df$month <- factor(df$month, levels = full_month_order)

castes_to_plot <- c("Worker", "Larvae") # Specify castes for plotting (using lowercase values)

# --- Filter for type == "mmag" ---
df <- df %>%
  filter(type == "mmag")

# Add apiary column (using lowercase colony)
df <- df %>%
  mutate(apiary = case_when(
    colony > 500 ~ "HB",
    colony <= 103 ~ "PT",
    TRUE ~ "OB"
  ))

# --- IMPORTANT: Handle unique average_genome_wide_cov per mag_id ---
# To ensure each MAG_id's coverage is counted only once, we group by unique identifiers
# and take the first 'average_genome_wide_cov' value.
# DO NOT create final_genus_category here yet, as it depends on global top_genera.
df_processed_initial <- df %>%
  group_by(mag_id, month, colony, rep, caste, genus, type, apiary, contig_presence_ratio) %>%
  summarise(
    avg_cov = first(average_genome_wide_cov), # Use average_genome_wide_cov as the input for relative abundance
    .groups = 'drop'
  )

cat("\n--- Script Complete (Initial Data Processing) ---\n")

# --- Define a consistent color palette for genera (Gruvbox-inspired) ---
gruvbox_colors <- c(
  "#a89984", "#fabd2f", "#fe8019", "#d7ba7d", "#cc241d", "#b16286",
  "#98971a", "#689d6a", "#458588", "#83a598", "#d3869b", "#8ec07c"
)

# --- GLOBAL: Identify Genera present in EVERY Apiary-Month combination (for consistent coloring) ---
# This will be the set of genera that get unique colors across ALL plots.
total_unique_apiary_month_combinations <- df_processed_initial %>%
  distinct(apiary, month) %>%
  nrow()

genera_in_all_apiary_months_globally_raw <- df_processed_initial %>%
  filter(genus != "unknown_mmag", genus != "no match") %>% # Exclude both "unknown_mmag" AND "no match"
  group_by(genus) %>%
  summarise(apiary_month_count = n_distinct(paste(apiary, month)), .groups = 'drop') %>%
  filter(apiary_month_count == total_unique_apiary_month_combinations) %>%
  pull(genus)

cat(paste("\n--- Genera present in EVERY Apiary-Month combination (Globally for coloring):",
          ifelse(length(genera_in_all_apiary_months_globally_raw) > 0,
                 paste(genera_in_all_apiary_months_globally_raw, collapse = ", "),
                 "None"), "---\n"))

# Calculate overall abundance for these globally consistent genera to pick the top 12
overall_abundance_for_globally_consistent <- df_processed_initial %>%
  filter(genus %in% genera_in_all_apiary_months_globally_raw) %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(avg_cov, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(total_abundance))

# Define 'top_genera' for coloring: top 12 most abundant among those consistently present globally
top_genera <- head(overall_abundance_for_globally_consistent$genus, 12)
top_n <- length(top_genera) # Actual number of top genera selected for specific colors (will be <= 12)

# --- GLOBAL: Define the color palette for all plots ---
top_genus_colors <- list()
if (top_n > 0) {
  colors_to_use <- head(gruvbox_colors, top_n)
  names(colors_to_use) <- top_genera
  top_genus_colors <- as.list(colors_to_use)
}
# Add color for "MAG below length threshold (0.5)"
top_genus_colors[["MAG below length threshold (0.5)"]] <- "#888888" # Dark grey
# Add color for "unknown_mmag"
top_genus_colors[["unknown_mmag"]] <- "#888888" # A slightly lighter grey for unknown
# Add color for "Other Genera" (must be last to catch everything else)
top_genus_colors[["Other Genera"]] <- "#888888"


# --- NOW: Add the final_genus_category to df_processed based on global top_genera and new rules ---
df_processed <- df_processed_initial %>%
  mutate(
    final_genus_category = case_when(
      contig_presence_ratio < 0.5 ~ "MAG below length threshold (0.5)",
      genus == "no match" ~ "Other Genera", # Explicitly handle "no match"
      genus %in% top_genera ~ genus, # Assign specific genus name if it's one of the global top
      is.na(genus) | genus == "" ~ "unknown_mmag", # Handle explicit unknown
      TRUE ~ "Other Genera" # Catch all remaining genera
    )
  )

cat("\n--- Script Complete (Global Setup) ---\n")

# --- Get unique Colony values for the outer loop ---
unique_colonies <- unique(df_processed$colony)

# --- Start Colony Loop ---
for (current_colony in unique_colonies) {
  cat(paste("\n--- Processing Colony Directory:", current_colony, "---\n"))
  
  # --- Create colony-specific directory within the base figure directory ---
  colony_figure_output_dir <- file.path(figure_output_base_dir, paste0("colony_", current_colony))
  dir.create(colony_figure_output_dir, recursive = TRUE, showWarnings = FALSE) # Create colony dir
  
  # --- Loop through each Caste for the current Colony ---
  for (current_caste in castes_to_plot) {
    cat(paste("  Processing Caste:", current_caste, "---\n"))
    
    # --- Filter data for the current Colony and Caste, aggregating across all Reps ---
    df_colony_caste <- df_processed %>%
      filter(colony == current_colony, caste == current_caste, type == "mmag")
    
    # --- Skip to next Caste if no data for this Colony-Caste combination ---
    if(nrow(df_colony_caste) == 0) {
      cat(paste("  No data for Colony", current_colony, "Caste", current_caste, ". Skipping to next caste.\n"))
      next
    }
    
    # --- Calculate total abundance per month for this Colony-Caste ---
    month_totals_colony_caste <- df_colony_caste %>%
      group_by(month) %>%
      summarise(total_month_abundance = sum(avg_cov)) # Using avg_cov
    
    # --- Calculate relative abundance per month and final_genus_category for this Colony-Caste ---
    df_rel_abundance_colony_caste <- df_colony_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(genus_month_abundance = sum(avg_cov)) %>% # Using avg_cov
      left_join(month_totals_colony_caste, by = "month") %>%
      mutate(relative_abundance = genus_month_abundance / total_month_abundance)
    
    # --- Re-aggregate data by Month and final_genus_category ---
    # This step is still needed to sum up relative abundances for categories like "Other Genera"
    df_rel_abundance_colony_caste_aggregated <- df_rel_abundance_colony_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(relative_abundance = sum(relative_abundance))
    
    # --- Convert Month to ordered factor for plotting ---
    df_rel_abundance_colony_caste_aggregated$month <- factor(df_rel_abundance_colony_caste_aggregated$month, levels = full_month_order)
    
    # --- Dynamic Month Order for the current Colony-Caste ---
    current_colony_caste_months_in_data <- unique(as.character(df_rel_abundance_colony_caste_aggregated$month))
    current_colony_caste_month_order <- full_month_order[full_month_order %in% current_colony_caste_months_in_data]
    
    # --- Create Area Plot for the current Colony-Caste ---
    colony_caste_area_plot <- ggplot(df_rel_abundance_colony_caste_aggregated,
                                     aes(x = month, y = relative_abundance, fill = final_genus_category, group = final_genus_category, color = final_genus_category)) +
      geom_area(position = "stack") +
      scale_x_discrete(limits = current_colony_caste_month_order) +
      scale_fill_manual(values = top_genus_colors) + # Use the GLOBAL color palette
      scale_color_manual(values = top_genus_colors) + # Ensure consistent colors for outlines too
      labs(title = paste("Microbial Genus Relative Abundance - Colony:", current_colony, "- Caste:", current_caste),
           x = "Month",
           y = "Relative Abundance",
           fill = "Bacterial Genus",
           color = "Bacterial Genus") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5), legend.position="right")
    
    print(colony_caste_area_plot)
    
    # --- Save the Colony-Caste Area plot to a file (optional) ---
    if (!dir.exists(colony_figure_output_dir)) {
      dir.create(colony_figure_output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    colony_caste_area_plot_file_name <- paste0("colony_", current_colony, "_caste_", current_caste, "_genus_area_plot_global_consistent_colors_aggregated_reps.png") # Updated filename
    colony_caste_area_plot_file_path <- file.path(colony_figure_output_dir, colony_caste_area_plot_file_name)
    ggsave(colony_caste_area_plot_file_path, plot = colony_caste_area_plot, width = 10, height = 6, units = "in")
    cat(paste("  Colony & Caste mMAG Genus Area plot (Global Consistent Colors, Aggregated Reps) saved for Colony", current_colony, "Caste", current_caste, "to:", colony_caste_area_plot_file_path, "\n"))
  } # --- END Caste Loop ---
} # --- END Colony Loop ---

cat("\n--- Script Complete ---\n")



##### COLONY - REP COMBINATION ########



castes_to_plot <- c("Worker") # Specify castes for plotting (using lowercase values)

# --- Filter for type == "mmag" ---
df <- df %>%
  filter(type == "mmag")

# Add apiary column (using lowercase colony)
df <- df %>%
  mutate(apiary = case_when(
    colony > 500 ~ "HB",
    colony <= 103 ~ "PT",
    TRUE ~ "OB"
  ))

# --- IMPORTANT: Handle unique average_genome_wide_cov per mag_id ---
# To ensure each MAG_id's coverage is counted only once, we group by unique identifiers
# and take the first 'average_genome_wide_cov' value.
# DO NOT create final_genus_category here yet, as it depends on global top_genera.
df_processed_initial <- df %>%
  group_by(mag_id, month, colony, rep, caste, genus, type, apiary, contig_presence_ratio) %>%
  summarise(
    avg_cov = first(average_genome_wide_cov), # Use average_genome_wide_cov as the input for relative abundance
    .groups = 'drop' # Ensure ungrouping after this summarise
  )

cat("\n--- Script Complete (Initial Data Processing) ---\n")

# --- Define a consistent color palette for genera (Gruvbox-inspired) ---
gruvbox_colors <- c(
  "#a89984", "#fabd2f", "#fe8019", "#d7ba7d", "#cc241d", "#b16286",
  "#98971a", "#689d6a", "#458588", "#83a598", "#d3869b", "#8ec07c"
)

# --- GLOBAL: Identify Genera present in EVERY Apiary-Month combination (for consistent coloring) ---
# This will be the set of genera that get unique colors across ALL plots.
total_unique_apiary_month_combinations <- df_processed_initial %>%
  distinct(apiary, month) %>%
  nrow()

genera_in_all_apiary_months_globally_raw <- df_processed_initial %>%
  filter(genus != "unknown_mmag", genus != "no match") %>% # Exclude both "unknown_mmag" AND "no match"
  group_by(genus) %>%
  summarise(apiary_month_count = n_distinct(paste(apiary, month)), .groups = 'drop') %>% # Ensure ungrouping
  filter(apiary_month_count == total_unique_apiary_month_combinations) %>%
  pull(genus)

cat(paste("\n--- Genera present in EVERY Apiary-Month combination (Globally for coloring):",
          ifelse(length(genera_in_all_apiary_months_globally_raw) > 0,
                 paste(genera_in_all_apiary_months_globally_raw, collapse = ", "),
                 "None"), "---\n"))

# Calculate overall abundance for these globally consistent genera to pick the top 12
overall_abundance_for_globally_consistent <- df_processed_initial %>%
  filter(genus %in% genera_in_all_apiary_months_globally_raw) %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(avg_cov, na.rm = TRUE), .groups = 'drop') %>% # Ensure ungrouping
  arrange(desc(total_abundance))

# Define 'top_genera' for coloring: top 12 most abundant among those consistently present globally
top_genera <- head(overall_abundance_for_globally_consistent$genus, 12)
top_n <- length(top_genera) # Actual number of top genera selected for specific colors (will be <= 12)

# --- GLOBAL: Define the color palette for all plots ---
top_genus_colors <- list()
if (top_n > 0) {
  colors_to_use <- head(gruvbox_colors, top_n)
  names(colors_to_use) <- top_genera
  top_genus_colors <- as.list(colors_to_use)
}
# Add color for "MAG below length threshold (0.5)"
top_genus_colors[["MAG below length threshold (0.5)"]] <- "#333333" # Dark grey
# Add color for "unknown_mmag"
top_genus_colors[["unknown_mmag"]] <- "#888888" # A slightly lighter grey for unknown
# Add color for "Other Genera" (must be last to catch everything else)
top_genus_colors[["Other Genera"]] <- "lightgrey"


# --- NOW: Add the final_genus_category to df_processed based on global top_genera and new rules ---
df_processed <- df_processed_initial %>%
  mutate(
    final_genus_category = case_when(
      contig_presence_ratio < 0.5 ~ "MAG below length threshold (0.5)",
      genus == "no match" ~ "Other Genera", # Explicitly handle "no match"
      genus %in% top_genera ~ genus, # Assign specific genus name if it's one of the global top
      is.na(genus) | genus == "" ~ "unknown_mmag", # Handle explicit unknown
      TRUE ~ "Other Genera" # Catch all remaining genera
    )
  )

cat("\n--- Script Complete (Global Setup) ---\n")

# --- Get unique Colony values and Rep combinations for looping ---
unique_colony_reps <- df_processed %>%
  distinct(colony, rep) %>%
  arrange(colony, rep)

# --- Get unique Colony values for the outer loop ---
unique_colonies <- unique(unique_colony_reps$colony)

# --- Start Colony Loop --- # Colony loop for directory structure
for (current_colony in unique_colonies) {
  cat(paste("\n--- Processing Colony Directory:", current_colony, "---\n"))
  
  # --- Create colony-specific directory within the base figure directory ---
  colony_figure_output_dir <- file.path(figure_output_base_dir, paste0("colony_", current_colony))
  dir.create(colony_figure_output_dir, recursive = TRUE, showWarnings = FALSE) # Create colony dir
  
  # --- Loop through each Rep within the current Colony --- # NEW: Rep loop within Colony loop
  for (i in 1:nrow(unique_colony_reps)) {
    current_colony_rep_row <- unique_colony_reps[i,]
    current_colony_for_rep <- current_colony_rep_row$colony
    current_rep <- current_colony_rep_row$rep
    
    # Only process if the current Colony in the Rep loop matches the current Colony directory loop
    if(current_colony_for_rep == current_colony) {
      cat(paste("\n--- Processing Colony-Rep Combination: Colony", current_colony, "- Rep:", current_rep, "---\n"))
      
      # --- Loop through each Caste for the current Colony-Rep combination ---
      for (current_caste in castes_to_plot) {
        cat(paste("  Processing Caste:", current_caste, "---\n"))
        
        # --- Filter data for the current Colony, Rep, and Caste ---
        df_colony_rep_caste <- df_processed %>%
          filter(colony == current_colony, rep == current_rep, caste == current_caste, type == "mmag")
        
        # --- Skip to next Caste if no data for this Colony-Rep-Caste combination ---
        if(nrow(df_colony_rep_caste) == 0) {
          cat(paste("  No data for Colony", current_colony, "Rep", current_rep, "Caste", current_caste, ". Skipping to next caste.\n"))
          next
        }
        
        # --- Calculate total abundance per month for this Colony-Rep-Caste ---
        month_totals_colony_rep_caste <- df_colony_rep_caste %>%
          group_by(month) %>%
          summarise(total_month_abundance = sum(avg_cov), .groups = 'drop') # Ensure ungrouping
        
        # --- Calculate relative abundance per month and final_genus_category for this Colony-Rep-Caste ---
        df_rel_abundance_colony_rep_caste <- df_colony_rep_caste %>%
          group_by(month, final_genus_category) %>%
          summarise(genus_month_abundance = sum(avg_cov), .groups = 'drop') %>% # Ensure ungrouping
          left_join(month_totals_colony_rep_caste, by = "month") %>%
          mutate(relative_abundance = genus_month_abundance / total_month_abundance)
        
        # --- Re-aggregate data by Month and final_genus_category ---
        df_rel_abundance_colony_rep_caste_aggregated <- df_rel_abundance_colony_rep_caste %>%
          group_by(month, final_genus_category) %>%
          summarise(relative_abundance = sum(relative_abundance), .groups = 'drop') # Ensure ungrouping
        
        # --- Convert Month to ordered factor for plotting ---
        df_rel_abundance_colony_rep_caste_aggregated$month <- factor(df_rel_abundance_colony_rep_caste_aggregated$month, levels = full_month_order)
        
        # --- Dynamic Month Order for the current Colony-Rep-Caste ---
        current_colony_rep_caste_months_in_data <- unique(as.character(df_rel_abundance_colony_rep_caste_aggregated$month))
        current_colony_rep_caste_month_order <- full_month_order[full_month_order %in% current_colony_rep_caste_months_in_data]
        
        # --- DEBUGGING: Check categories and colors before plotting ---
        unique_categories_in_data <- unique(df_rel_abundance_colony_rep_caste_aggregated$final_genus_category)
        missing_colors_for_categories <- setdiff(unique_categories_in_data, names(top_genus_colors))
        if(length(missing_colors_for_categories) > 0) {
          cat(paste("  WARNING: Missing colors for categories in plot data:", paste(missing_colors_for_categories, collapse = ", "), "\n"))
          cat(paste("  Available colors for:", paste(names(top_genus_colors), collapse = ", "), "\n"))
        }
        # --- END DEBUGGING ---
        
        # --- Create Area Plot for the current Colony-Rep-Caste ---
        colony_rep_caste_area_plot <- ggplot(df_rel_abundance_colony_rep_caste_aggregated,
                                             aes(x = month, y = relative_abundance, fill = final_genus_category, group = final_genus_category, color = final_genus_category)) +
          geom_area(position = "stack") +
          scale_x_discrete(limits = current_colony_rep_caste_month_order) +
          scale_fill_manual(values = top_genus_colors) + # Use the GLOBAL color palette
          scale_color_manual(values = top_genus_colors) + # Ensure consistent colors for outlines too
          labs(title = paste("Microbial Genus Relative Abundance - Colony:", current_colony, "- Rep:", current_rep, "- Caste:", current_caste),
               x = "Month",
               y = "Relative Abundance",
               fill = "Bacterial Genus",
               color = "Bacterial Genus") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5), legend.position="right")
        
        print(colony_rep_caste_area_plot)
        
        # --- Save the Colony-Rep-Caste Area plot to a file (optional) ---
        if (!dir.exists(colony_figure_output_dir)) {
          dir.create(colony_figure_output_dir, recursive = TRUE, showWarnings = FALSE)
        }
        colony_rep_caste_area_plot_file_name <- paste0("colony_", current_colony, "_rep_", current_rep, "_caste_", current_caste, "_genus_area_plot_global_consistent_colors.png") # Updated filename
        colony_rep_caste_area_plot_file_path <- file.path(colony_figure_output_dir, colony_rep_caste_area_plot_file_name)
        ggsave(colony_rep_caste_area_plot_file_path, plot = colony_rep_caste_area_plot, width = 10, height = 6, units = "in")
        cat(paste("  Colony-Rep & Caste mMAG Genus Area plot (Global Consistent Colors) saved for Colony", current_colony, "Rep", current_rep, "Caste", current_caste, "to:", colony_rep_caste_area_plot_file_path, "\n"))
      } # --- END Caste Loop ---
    } # --- END if current_colony_for_rep == current_colony ---
  } # --- END Rep Loop ---
} # --- END Colony Directory Loop ---

cat("\n--- Script Complete ---\n")



#### LARVAE #####
library(scales) # For hue_pal() to generate distinct colors


# First, convert all column names to lowercase
df <- df %>%
  rename_with(tolower)
df <- df %>%
  filter(caste %in% c("Larvae"))

# Filter data for df_one (using lowercase column names and values)
# Keep only "larvae" as per new requirement
castes_to_plot <- c("Larvae")

# Define month order (using lowercase column names)
full_month_order <- c("May", "June", "July", "August")
df$month <- factor(df$month, levels = full_month_order)

# --- Filter for type == "mmag" ---
df <- df %>%
  filter(type == "mmag")

# Add apiary column (using lowercase colony)
df <- df %>%
  mutate(apiary = case_when(
    colony > 500 ~ "HB",
    colony <= 103 ~ "PT",
    TRUE ~ "OB"
  ))

# --- IMPORTANT: Handle unique average_genome_wide_cov per mag_id ---
# To ensure each MAG_id's coverage is counted only once, we group by unique identifiers
# and take the first 'average_genome_wide_cov' value.
df_processed_initial <- df %>%
  group_by(mag_id, month, colony, rep, caste, genus, type, apiary, contig_presence_ratio) %>%
  summarise(
    avg_cov = first(average_genome_wide_cov), # Use average_genome_wide_cov as the input for relative abundance
    .groups = 'drop' # Ensure ungrouping after this summarise
  )

cat("\n--- Script Complete (Initial Data Processing) ---\n")

# --- Filter df_processed_initial to only include relevant colonies for Larvae ---
# As stated: "Only Colony == 100 and Colony == 999 contain Caste == Larvae"
df_processed_larvae_only <- df_processed_initial %>%
  filter(colony %in% c(100, 999), caste == "Larvae")

# --- Define final_genus_category for Larvae-specific coloring ---
df_processed_larvae <- df_processed_larvae_only %>%
  mutate(
    final_genus_category = case_when(
      contig_presence_ratio < 0.5 ~ "MAG below length threshold (0.5)",
      TRUE ~ as.character(genus) 
    )
  )

cat("\n--- Script Complete (Larvae-Specific Data Preparation) ---\n")

# --- Define fixed colors for special categories ---
fixed_special_colors <- c(
  "MAG below length threshold (0.5)" = "#333333", # Dark grey for below threshold
  "unknown_mmag" = "#888888" # A slightly lighter grey for unknown
)

# --- SET 1: Generate plots for each Colony (Larvae only, aggregating Reps) ---
figure_output_base_dir_colony <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/figures/figures_by_colony_larvae_aggregated"
unique_colonies <- unique(df_processed_larvae$colony)

for (current_colony in unique_colonies) {
  cat(paste("\n--- Processing Colony Directory (Aggregated Reps):", current_colony, "---\n"))
  
  colony_figure_output_dir <- file.path(figure_output_base_dir_colony, paste0("colony_", current_colony))
  dir.create(colony_figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (current_caste in castes_to_plot) { # This loop will only run for "larvae"
    cat(paste("  Processing Caste:", current_caste, "---\n"))
    
    df_colony_caste <- df_processed_larvae %>%
      filter(colony == current_colony, caste == current_caste)
    
    if(nrow(df_colony_caste) == 0) {
      cat(paste("  No data for Colony", current_colony, "Caste", current_caste, ". Skipping to next caste.\n"))
      next
    }
    
    month_totals_colony_caste <- df_colony_caste %>%
      group_by(month) %>%
      summarise(total_month_abundance = sum(avg_cov), .groups = 'drop')
    
    df_rel_abundance_colony_caste <- df_colony_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(genus_month_abundance = sum(avg_cov), .groups = 'drop') %>%
      left_join(month_totals_colony_caste, by = "month") %>%
      mutate(relative_abundance = genus_month_abundance / total_month_abundance)
    
    df_rel_abundance_colony_caste_aggregated <- df_rel_abundance_colony_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(relative_abundance = sum(relative_abundance), .groups = 'drop')
    
    df_rel_abundance_colony_caste_aggregated$month <- factor(df_rel_abundance_colony_caste_aggregated$month, levels = full_month_order)
    
    current_colony_caste_months_in_data <- unique(as.character(df_rel_abundance_colony_caste_aggregated$month))
    current_colony_caste_month_order <- full_month_order[full_month_order %in% current_colony_caste_months_in_data]
    
    # --- Dynamic Color Palette for current plot ---
    unique_categories_in_plot <- unique(df_rel_abundance_colony_caste_aggregated$final_genus_category)
    
    # Separate fixed categories from dynamic Genus names
    dynamic_genera <- setdiff(unique_categories_in_plot, names(fixed_special_colors))
    
    plot_colors <- fixed_special_colors
    if (length(dynamic_genera) > 0) {
      # Generate a distinct palette for Genus names
      genus_colors <- scales::hue_pal()(length(dynamic_genera))
      names(genus_colors) <- dynamic_genera
      plot_colors <- c(plot_colors, genus_colors)
    }
    
    colony_caste_area_plot <- ggplot(df_rel_abundance_colony_caste_aggregated,
                                     aes(x = month, y = relative_abundance, fill = final_genus_category, group = final_genus_category, color = final_genus_category)) +
      geom_bar(stat = "identity", position = "stack") + # Changed to geom_bar
      scale_x_discrete(limits = current_colony_caste_month_order) +
      scale_fill_manual(values = plot_colors) +
      scale_color_manual(values = plot_colors) +
      labs(title = paste("Larval Microbiome Relative Abundance - Colony:", current_colony, "- Caste:", current_caste, "(Aggregated Reps)"),
           x = "Month",
           y = "Relative Abundance",
           fill = "Genus / Category",
           color = "Genus / Category") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5), legend.position="right")
    
    print(colony_caste_area_plot)
    
    colony_caste_area_plot_file_name <- paste0("colony_", current_colony, "_caste_", current_caste, "_larvae_aggregated_reps.png")
    colony_caste_area_plot_file_path <- file.path(colony_figure_output_dir, colony_caste_area_plot_file_name)
    ggsave(colony_caste_area_plot_file_path, plot = colony_caste_area_plot, width = 10, height = 6, units = "in")
    cat(paste("  Colony Larval Microbiome plot (Aggregated Reps) saved for Colony", current_colony, "Caste", current_caste, "to:", colony_caste_area_plot_file_path, "\n"))
  }
}

# --- SET 2: Generate plots for each Colony-Rep combination (Larvae only) ---
figure_output_base_dir_colony_rep <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/figures/figures_by_colony_rep_larvae"
unique_colony_reps <- df_processed_larvae %>%
  distinct(colony, rep) %>%
  arrange(colony, rep)

for (i in 1:nrow(unique_colony_reps)) {
  current_colony_rep_row <- unique_colony_reps[i,]
  current_colony_for_rep <- current_colony_rep_row$colony
  current_rep <- current_colony_rep_row$rep
  
  colony_rep_figure_output_dir <- file.path(figure_output_base_dir_colony_rep, paste0("colony_", current_colony_for_rep, "_rep_", current_rep))
  dir.create(colony_rep_figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(paste("\n--- Processing Colony-Rep Combination (Larvae): Colony", current_colony_for_rep, "- Rep:", current_rep, "---\n"))
  
  for (current_caste in castes_to_plot) { # This loop will only run for "larvae"
    cat(paste("  Processing Caste:", current_caste, "---\n"))
    
    df_colony_rep_caste <- df_processed_larvae %>%
      filter(colony == current_colony_for_rep, rep == current_rep, caste == current_caste)
    
    if(nrow(df_colony_rep_caste) == 0) {
      cat(paste("  No data for Colony", current_colony_for_rep, "Rep", current_rep, "Caste", current_caste, ". Skipping to next caste.\n"))
      next
    }
    
    month_totals_colony_rep_caste <- df_colony_rep_caste %>%
      group_by(month) %>%
      summarise(total_month_abundance = sum(avg_cov), .groups = 'drop')
    
    df_rel_abundance_colony_rep_caste <- df_colony_rep_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(genus_month_abundance = sum(avg_cov), .groups = 'drop') %>%
      left_join(month_totals_colony_rep_caste, by = "month") %>%
      mutate(relative_abundance = genus_month_abundance / total_month_abundance)
    
    df_rel_abundance_colony_rep_caste_aggregated <- df_rel_abundance_colony_rep_caste %>%
      group_by(month, final_genus_category) %>%
      summarise(relative_abundance = sum(relative_abundance), .groups = 'drop')
    
    df_rel_abundance_colony_rep_caste_aggregated$month <- factor(df_rel_abundance_colony_rep_caste_aggregated$month, levels = full_month_order)
    
    current_colony_rep_caste_months_in_data <- unique(as.character(df_rel_abundance_colony_rep_caste_aggregated$month))
    current_colony_rep_caste_month_order <- full_month_order[full_month_order %in% current_colony_rep_caste_months_in_data]
    
    # --- Dynamic Color Palette for current plot ---
    unique_categories_in_plot <- unique(df_rel_abundance_colony_rep_caste_aggregated$final_genus_category)
    
    # Separate fixed categories from dynamic Genus names
    dynamic_genera <- setdiff(unique_categories_in_plot, names(fixed_special_colors))
    
    plot_colors <- fixed_special_colors
    if (length(dynamic_genera) > 0) {
      genus_colors <- scales::hue_pal()(length(dynamic_genera))
      names(genus_colors) <- dynamic_genera
      plot_colors <- c(plot_colors, genus_colors)
    }
    
    colony_rep_caste_area_plot <- ggplot(df_rel_abundance_colony_rep_caste_aggregated,
                                         aes(x = month, y = relative_abundance, fill = final_genus_category, group = final_genus_category, color = final_genus_category)) +
      geom_bar(stat = "identity", position = "stack") + # Changed to geom_bar
      scale_x_discrete(limits = current_colony_rep_caste_month_order) +
      scale_fill_manual(values = plot_colors) +
      scale_color_manual(values = plot_colors) +
      labs(title = paste("Larval Microbiome Relative Abundance - Colony:", current_colony_for_rep, "- Rep:", current_rep, "- Caste:", current_caste),
           x = "Month",
           y = "Relative Abundance",
           fill = "Genus / Category",
           color = "Genus / Category") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5), legend.position="right")
    
    print(colony_rep_caste_area_plot)
    
    colony_rep_caste_area_plot_file_name <- paste0("colony_", current_colony_for_rep, "_rep_", current_rep, "_caste_", current_caste, "_larvae.png")
    colony_rep_caste_area_plot_file_path <- file.path(colony_rep_figure_output_dir, colony_rep_caste_area_plot_file_name)
    ggsave(colony_rep_caste_area_plot_file_path, plot = colony_rep_caste_area_plot, width = 10, height = 6, units = "in")
    cat(paste("  Colony-Rep Larval Microbiome plot saved for Colony", current_colony_for_rep, "Rep", current_rep, "Caste", current_caste, "to:", colony_rep_caste_area_plot_file_path, "\n"))
  }
}

cat("\n--- Script Complete ---\n")

