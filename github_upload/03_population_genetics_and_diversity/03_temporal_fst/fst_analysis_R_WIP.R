library(dplyr)
library(ggplot2)
library(rstatix)
library(dunn.test)
library(tidyr)
library(ggpubr)
library(purrr)

df = read.csv('normalized_fst_multiple_colony.csv')
df$Fst_Method1[df$Fst_Method1 < 0] <- 0


### Helper function to add DaysSinceSampling to csvs

add_cumulative_time_from_month <- function(data) {
  # Define the base time values for each month in order
  month_base_times <- c(
    "May" = 0,
    "June" = 36,
    "July" = 31,
    "August" = 34,
    "September" = 28,
    "October" = 22,
    "November" = 42,
    "January" = 47,
    "February" = 31
  )
  
  # Get the order of the months as provided
  ordered_months <- names(month_base_times)
  
  # Calculate the cumulative time values
  cumulative_times <- cumsum(month_base_times)
  
  # Create a mapping between Month and Cumulative Time
  month_to_cumulative_time_map <- setNames(cumulative_times, ordered_months)
  
  # Add the "Time" column based on the "Month" column
  data$Time <- month_to_cumulative_time_map[data$Month2]
  
  # Check if any months in the data were not found in the mapping
  if (any(is.na(data$Time))) {
    warning("Warning: Some months in the 'Month' column were not found in the mapping and have been assigned NA in the 'Time' column.")
  }
  
  return(data)
}

#####
# Create the base boxplot
boxplot <- ggplot(data = df, aes(x = lifestyle, y = Fst_Method1)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(title = "Gene-Level Normalized Fst by Lifestyle",
       y = "Fst_Method1",
       x = "Lifestyle")+
  theme_minimal()

# Perform pairwise Dunn's test
pairwise_comparisons <- df %>%
  dunn_test(Fst_Method1 ~ lifestyle, p.adjust.method = "bonferroni") %>% # You can change the p.adjust.method
  adjust_pvalue() %>%
  add_xy_position(x = "lifestyle")
# Manually adjust y.position for staggering
# You might need to adjust these values based on your specific lifestyles and data range
pairwise_comparisons <- pairwise_comparisons %>%
  mutate(
    y.position = case_when(
      group1 == levels(factor(df$lifestyle))[1] & group2 == levels(factor(df$lifestyle))[2] ~ 0.45,
      group1 == levels(factor(df$lifestyle))[1] & group2 == levels(factor(df$lifestyle))[3] ~ 0.5,
      group1 == levels(factor(df$lifestyle))[2] & group2 == levels(factor(df$lifestyle))[3] ~ 0.55,
      TRUE ~ y.position # Keep default if not specified
    )
  )

# Add p-values to the plot
boxplot_with_pvalues <- boxplot +
  stat_pvalue_manual(pairwise_comparisons,
                     label = "p.adj")

# Print the plot
print(boxplot_with_pvalues)


###SHARED COMPARISONS, BETWEEN COLONY AND WITHIN COLONY

# Read the first CSV file
between_colony_df <- read.csv('gene_level_between_colony_fst_multiple_colony.csv')
between_colony_df$Source <- "Between Colony"

# Read the second CSV file (assuming this is the intended one)
temporal_df <- read.csv('gene_level_temporal_fst_multiple_colony.csv')
temporal_df$Source <- "Temporal"

# Combine the two dataframes
combined_df <- bind_rows(between_colony_df, temporal_df)

# Create the base boxplot
boxplot2 <- ggplot(combined_df, aes(x = Source, y = Fst_Method1, fill = Source)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.75)) +
  labs(title = "Comparison of Gene-Level Fst Distributions",
       y = "Fst_Method1",
       x = "Data Source",
       fill = "Data Source") +
  theme_minimal()

# Perform pairwise Dunn's test
pairwise_comparisons_shared_within_between <- combined_df %>%
  dunn_test(Fst_Method1 ~ Source, p.adjust.method = "bonferroni") %>%
  adjust_pvalue()

boxplot2
# Prepare the data frame for stat_pvalue_manual
pvalue_df <- pairwise_comparisons_shared_within_between %>%
  mutate(group1 = Source[1], group2 = Source[2]) %>% # Assuming only one comparison
  mutate(y.position = 0.65) # Adjust this value to position the p-value
boxplot2


########################SAME AS ABOVE BUT NOT SHARED #####################3

vMAG_within = read.csv('gene_level_month_to_month_fst_single_colony.csv')
vMAG_between = read.csv('gene_level_pairwise_fst_results.csv')

vMAG_within$Fst_Method1[vMAG_within$Fst_Method1 < 0] <- 0
vMAG_between$Fst_Method1[vMAG_between$Fst_Method1 < 0] <- 0

# Add a Source column to each dataframe
vMAG_between <- vMAG_between %>% mutate(Source = "vMAG Between")
vMAG_within <- vMAG_within %>% mutate(Source = "vMAG Within")

# Combine the two dataframes
combined_vMAG_df <- bind_rows(vMAG_between, vMAG_within)

# Create the base boxplot
boxplot_vMAG <- ggplot(combined_vMAG_df, aes(x = Source, y = Fst_Method1, fill = Source)) +
  geom_boxplot() +
  labs(title = "Comparison of Fst_Method1 Distributions",
       y = "Fst_Method1",
       x = "Data Source",
       fill = "Data Source") +
  theme_minimal()
boxplot_vMAG
# Perform pairwise Dunn's test
pairwise_comparisons_vMAG <- combined_vMAG_df %>%
  dunn_test(Fst_Method1 ~ Source, p.adjust.method = "bonferroni") %>%
  adjust_pvalue()

# Prepare the data frame for stat_pvalue_manual
pvalue_vMAG_df <- pairwise_comparisons_vMAG %>%
  mutate(group1 = Source[1], group2 = Source[2]) %>% # Assuming "vMAG Between" is the first level
  mutate(y.position = 0.65) # Adjust this value to position the p-value

# Add p-values to the plot
boxplot_with_pvalues_vMAG <- boxplot_vMAG +
  stat_pvalue_manual(pvalue_vMAG_df,
                     label = "p.adj")

# Print the plot
print(boxplot_with_pvalues_vMAG)



######################### SHARED, SAME MONTH #############3

dfsm_df = read.csv('gene_level_pairwise_fst_results_with_annotations.csv')
dfsm_df$Fst_Method1[dfsm_df$Fst_Method1 < 0] <- 0


ggplot(data=dfsm_df,aes(x=lifestyle,y=Fst_Method1))+
  geom_boxplot()#+
#scale_y_continuous(limits = c(0, 0.075))


ggplot(data=dfsm_df,aes(x=functional_category,y=Fst_Method1))+
  geom_boxplot()#+
#scale_y_continuous(limits = c(0, 0.075))

pairwise_dfsm_testtest <- dfsm_df %>%
  dunn_test(Fst_Method1 ~ lifestyle, p.adjust.method = "BH") %>%
  adjust_pvalue()

pairwise_dfsm_testtest









################ WITHIN COLONY, SINGLE COLONY ACROSS TIME ###########


genic_im = read.csv('gene_level_vs_initial_month_fst_single_colony_with_annotations.csv')
genic_mm = read.csv('gene_level_month_to_month_fst_single_colony_with_annotations.csv')
genic_im$Fst_Method1[genic_im$Fst_Method1 < 0] <- 0
genic_mm$Fst_Method1[genic_mm$Fst_Method1 < 0] <- 0

ggplot(data=genic_im,aes(x=lifestyle,y=Fst_Method1))+
  geom_boxplot()#+
 #scale_y_continuous(limits = c(0, 0.075))

ggplot(data=genic_im, aes(x=Fst_Method1,color=lifestyle))+
  geom_density()

pairwise_genic_im_fst_dunn_test <- genic_im %>%
  dunn_test(Fst_Method1 ~ lifestyle, p.adjust.method = "BH") %>%
  adjust_pvalue()
pairwise_genic_im_fst_dunn_test

pairwise_genic_mm_fst_dunn_test <- genic_mm %>%
  dunn_test(Fst_Method1 ~ lifestyle, p.adjust.method = "BH") %>%
  adjust_pvalue()
lytic_median = genic_im[genic_im$lifestyle == "Lytic",]
lysogenic_median = genic_im[genic_im$lifestyle == "Temperate",]
prophage_median = genic_im[genic_im$lifestyle == "Prophage",]

median(lytic_median$Fst_Method1, na.rm = TRUE)
median(lysogenic_median$Fst_Method1, na.rm = TRUE)
median(prophage_median$Fst_Method1, na.rm = TRUE)


########### COMPARISONS OF WITHIN AND BETWEEN DISTRIBUTIONS
genic_btwn = read.csv('gene_level_pairwise_fst_results_with_annotations.csv')
genic_btwn$Fst_Method1[genic_btwn$Fst_Method1 < 0] <- 0

# --- Prepare Data ---
# Add a 'Source' column to distinguish the dataframes
genic_im <- genic_im %>% mutate(Source = "Temporal")
genic_btwn <- genic_btwn %>% mutate(Source = "Distance")

# Combine the dataframes
combined_genic_df <- bind_rows(genic_im, genic_btwn)

# Find unique lifestyles
unique_lifestyles <- unique(combined_genic_df$lifestyle)
unique_lifestyles <- unique_lifestyles[!is.na(unique_lifestyles)] # Remove potential NA lifestyles

# --- Loop Through Lifestyles and Generate Plots/Tests ---

for (current_lifestyle in unique_lifestyles) {
  # Filter data for the current lifestyle
  subset_data <- combined_genic_df %>%
    filter(lifestyle == current_lifestyle)
  
  # --- Perform Dunn's Test ---
  # Dunn's test for pairwise comparison (Temporal vs Distance) within the lifestyle
  dunn_results <- subset_data %>%
    dunn_test(Fst_Method1 ~ Source, p.adjust.method = "BH") %>% # Using BH for p-value adjustment
    adjust_pvalue() # Ensure adjusted p-values are calculated
  
  # --- Create Plot ---
  # Prepare the data frame for stat_pvalue_manual (since it's a simple pairwise comparison)
  # We need to ensure the group names match the 'Source' levels in the plot
  pvalue_df <- dunn_results %>%
    mutate(group1 = as.character(group1), group2 = as.character(group2)) %>% # Ensure group names are characters
    mutate(
      # Assign y.position - you might need to adjust this value based on your data range
      y.position = max(subset_data$Fst_Method1, na.rm = TRUE) * 1.1 # Place above max Fst value
    )
  
  # Create the boxplot
  plot <- ggplot(subset_data, aes(x = Source, y = Fst_Method1, fill = Source)) +
    geom_boxplot() +
    # You might need to adjust y-axis limits based on your data
    # scale_y_continuous(limits = c(0, max(subset_data$Fst_Method1, na.rm = TRUE) * 1.2)) +
    labs(
      title = paste("Fst Distribution for Lifestyle:", current_lifestyle),
      subtitle = paste("Dunn Test p.adj:", format.pval(pvalue_df$p.adj[1], digits = 3)), # Display p.adj in subtitle
      y = "Fst_Method1",
      x = "Data Source",
      fill = "Data Source"
    ) +
    theme_minimal() # Use a minimal theme
  
  # Add p-values to the plot using stat_pvalue_manual
  # We add it directly if the dunn_results dataframe is structured correctly
  # If only one comparison per plot, can add manually or via subtitle
  # For this case (only 2 groups), adding to subtitle is clean
  # If you prefer the line annotation, uncomment the stat_pvalue_manual part and adjust y.position
  
  # plot <- plot +
  # stat_pvalue_manual(pvalue_df,
  #                      label = "p.adj",
  #                      y.position = pvalue_df$y.position)
  
  
  # Print the plot
  print(plot)
  
  # Optional: Print the Dunn test results for the current lifestyle
  cat("\nDunn test results for Lifestyle:", current_lifestyle, "\n")
  print(dunn_results)
  cat("---\n")
  
}







# --- Loop Through Functional Categories and Generate Plots/Tests ---

for (current_functional_category in unique_functional_categories) {
  # Filter data for the current functional category
  subset_data_func <- combined_genic_df %>%
    filter(functional_category == current_functional_category)

  # Check if there's data for both sources in this category
  if (length(unique(subset_data_func$Source)) < 2) {
    cat("Skipping functional category:", current_functional_category, "- Data for both sources not available.\n")
    next # Skip to the next category if only one source is present
  }

  # --- Perform Dunn's Test ---
  # Dunn's test for pairwise comparison (Temporal vs Distance) within the functional category
  dunn_results_func <- subset_data_func %>%
    dunn_test(Fst_Method1 ~ Source, p.adjust.method = "BH") %>% # Using BH for p-value adjustment
    adjust_pvalue() # Ensure adjusted p-values are calculated

  # --- Create Plot ---
  # Format the p-value for display in the subtitle
  p_value_text <- "Dunn Test p.adj: NA" # Default in case test fails or no groups
  if(nrow(dunn_results_func) > 0) {
      p_value_text <- paste("Dunn Test p.adj:", format.pval(dunn_results_func$p.adj[1], digits = 3))
  }

  # Create the boxplot
  plot_func <- ggplot(subset_data_func, aes(x = Fst_Method1, color = Source)) +
    geom_density() +
    # You might need to adjust y-axis limits based on your data
    # scale_y_continuous(limits = c(0, max(subset_data_func$Fst_Method1, na.rm = TRUE) * 1.2)) +
    labs(
      title = paste("Fst Distribution for Functional Category:", current_functional_category),
      subtitle = p_value_text, # Display p.adj in subtitle
      y = "Fst_Method1",
      x = "Data Source",
      fill = "Data Source"
    ) +
    theme_minimal() # Use a minimal theme

  # Print the plot
  print(plot_func)

  # Optional: Print the Dunn test results for the current functional category
  cat("\nDunn test results for Functional Category:", current_functional_category, "\n")
  print(dunn_results_func)
  cat("---\n")
}




# Find unique functional categories
unique_functional_categories <- unique(combined_genic_df$functional_category)
unique_functional_categories <- unique_functional_categories[!is.na(unique_functional_categories)] # Remove potential NA functional categories

# --- Loop Through Functional Categories and Generate Plots/Tests ---

for (current_functional_category in unique_functional_categories) {
  # Filter data for the current functional category
  subset_data_func <- combined_genic_df %>%
    filter(functional_category == current_functional_category)
  
  # Check if there's data for both sources in this category
  if (length(unique(subset_data_func$Source)) < 2) {
    cat("Skipping functional category:", current_functional_category, "- Data for both sources not available.\n")
    next # Skip to the next category if only one source is present
  }
  
  # --- Perform Dunn's Test ---
  # Dunn's test for pairwise comparison (Temporal vs Distance) within the functional category
  dunn_results_func <- subset_data_func %>%
    dunn_test(Fst_Method1 ~ Source, p.adjust.method = "BH") %>% # Using BH for p-value adjustment
    adjust_pvalue() # Ensure adjusted p-values are calculated
  
  # --- Create Plot ---
  # Format the p-value for display in the subtitle
  p_value_text <- "Dunn Test p.adj: NA" # Default in case test fails or no groups
  if(nrow(dunn_results_func) > 0) {
    p_value_text <- paste("Dunn Test p.adj:", format.pval(dunn_results_func$p.adj[1], digits = 3))
  }
  
  # Create the boxplot
  plot_func <- ggplot(subset_data_func, aes(x = Source, y = Fst_Method1, fill = Source)) +
    geom_boxplot() +
    # You might need to adjust y-axis limits based on your data
    # scale_y_continuous(limits = c(0, max(subset_data_func$Fst_Method1, na.rm = TRUE) * 1.2)) +
    labs(
      title = paste("Fst Distribution for Functional Category:", current_functional_category),
      subtitle = p_value_text, # Display p.adj in subtitle
      y = "Fst_Method1",
      x = "Data Source",
      fill = "Data Source"
    ) +
    theme_minimal() # Use a minimal theme
  
  # Print the plot
  print(plot_func)
  
  # Optional: Print the Dunn test results for the current functional category
  cat("\nDunn test results for Functional Category:", current_functional_category, "\n")
  print(dunn_results_func)
  cat("---\n")
}






















############# WITHIN COLONY AND BETWEEN, SAME VMAG ###################3
test_df = read.csv('gene_level_normalized_fst_multiple_colony_with_annotations.csv')
test_df2 = read.csv('gene_level_between_colony_fst_multiple_colony_with_annotations.csv')
test_df$Fst_Method1[test_df$Fst_Method1 < 0] <- 0
test_df2$Fst_Method1[test_df2$Fst_Method1 < 0] <- 0


ggplot(data=test_df,aes(x=lifestyle, y = Fst_Method1))+
  geom_boxplot()
dunn_results_test <- test_df %>%
  dunn_test(Fst_Method1 ~ lifestyle, p.adjust.method = "BH") %>% # Using BH for p-value adjustment
  adjust_pvalue() # Ensure adjusted p-values are calculated

dunn_results_test
# --- Prepare Data ---
# Add a 'Source' column to distinguish the dataframes
test_df <- test_df %>% mutate(Source = "Temporal")
test_df2 <- test_df2 %>% mutate(Source = "Distance")

# Combine the dataframes
combined_test_df <- bind_rows(test_df, test_df2)

# Find unique lifestyles
unique_lifestyles_test <- unique(combined_test_df$lifestyle)
unique_lifestyles_test <- unique_lifestyles_test[!is.na(unique_lifestyles_test)] # Remove potential NA lifestyles

# --- Loop Through Lifestyles and Generate Plots/Tests ---

for (current_lifestyle in unique_lifestyles_test) {
  # Filter data for the current lifestyle
  subset_data_test <- combined_test_df %>%
    filter(lifestyle == current_lifestyle)
  
  # --- Perform Dunn's Test ---
  # Dunn's test for pairwise comparison (Temporal vs Distance) within the lifestyle
  dunn_results_test <- subset_data_test %>%
    dunn_test(Fst_Method1 ~ Source, p.adjust.method = "BH") %>% # Using BH for p-value adjustment
    adjust_pvalue() # Ensure adjusted p-values are calculated
  
  # --- Create Plot ---
  # Format the p-value for display in the subtitle
  p_value_text <- "Dunn Test p.adj: NA" # Default in case test fails or no groups
  if(nrow(dunn_results_test) > 0) {
    p_value_text <- paste("Dunn Test p.adj:", format.pval(dunn_results_test$p.adj[1], digits = 3))
  }
  
  
  # Create the boxplot
  plot_test <- ggplot(subset_data_test, aes(x = Source, y = Fst_Method1, fill = Source)) +
    geom_boxplot() +
    # You might need to adjust y-axis limits based on your data
    # scale_y_continuous(limits = c(0, max(subset_data_test$Fst_Method1, na.rm = TRUE) * 1.2)) +
    labs(
      title = paste("Fst Distribution for Lifestyle:", current_lifestyle),
      subtitle = p_value_text, # Display p.adj in subtitle
      y = "Fst_Method1",
      x = "Data Source",
      fill = "Data Source"
    ) +
    theme_minimal() # Use a minimal theme
  
  # Print the plot
  print(plot_test)
  
  # Optional: Print the Dunn test results for the current lifestyle
  cat("\nDunn test results for Lifestyle:", current_lifestyle, "\n")
  print(dunn_results_test)
  cat("---\n")
}



# Find unique functional categories
unique_functional_categories_test <- unique(combined_test_df$functional_category)
unique_functional_categories_test <- unique_functional_categories_test[!is.na(unique_functional_categories_test)] # Remove potential NA functional categories

# --- Loop Through Functional Categories and Generate Plots/Tests ---

for (current_functional_category in unique_functional_categories_test) {
  # Filter data for the current functional category
  subset_data_func_test <- combined_test_df %>%
    filter(functional_category == current_functional_category)
  
  # Check if there's data for both sources in this category
  if (length(unique(subset_data_func_test$Source)) < 2) {
    cat("Skipping functional category:", current_functional_category, "- Data for both sources not available in test_df/test_df2 combined data.\n")
    next # Skip to the next category if only one source is present
  }
  
  # --- Perform Dunn's Test ---
  # Dunn's test for pairwise comparison (Temporal vs Distance) within the functional category
  dunn_results_func_test <- subset_data_func_test %>%
    dunn_test(Fst_Method1 ~ Source, p.adjust.method = "BH") %>% # Using BH for p-value adjustment
    adjust_pvalue() # Ensure adjusted p-values are calculated
  
  # --- Create Plot ---
  # Format the p-value for display in the subtitle
  p_value_text <- "Dunn Test p.adj: NA" # Default in case test fails or no groups
  if(nrow(dunn_results_func_test) > 0) {
    p_value_text <- paste("Dunn Test p.adj:", format.pval(dunn_results_func_test$p.adj[1], digits = 3))
  }
  
  # Create the boxplot
  plot_func_test <- ggplot(subset_data_func_test, aes(x = Source, y = Fst_Method1, fill = Source)) +
    geom_boxplot() +
    # You might need to adjust y-axis limits based on your data
    # scale_y_continuous(limits = c(0, max(subset_data_func_test$Fst_Method1, na.rm = TRUE) * 1.2)) +
    labs(
      title = paste("Fst Distribution for Functional Category:", current_functional_category),
      subtitle = p_value_text, # Display p.adj in subtitle
      y = "Fst_Method1",
      x = "Data Source",
      fill = "Data Source"
    ) +
    theme_minimal() # Use a minimal theme
  
  # Print the plot
  print(plot_func_test)
  
  # Optional: Print the Dunn test results for the current functional category
  cat("\nDunn test results for Functional Category:", current_functional_category, " (test_df/test_df2 combined data)\n")
  print(dunn_results_func_test)
  cat("---\n")
}

genome_wide_normalized = read.csv('normalized_fst_multiple_colony.csv')


TIME_MCMM_normalized_fst_multiple_colony <- add_cumulative_time_from_month(genome_wide_normalized)

TIME_MCMM_normalized_fst_multiple_colony <- TIME_MCMM_normalized_fst_multiple_colony %>%
  mutate(
    normalizedTransform = case_when(
      Normalized_Fst_Method1 < 0 ~ 0, # New condition: If Normalized_Fst_Method1 < 0, set to 0
      normalized_flag == FALSE & Fst_Method1 > 0 ~ 1,
      TRUE ~ Normalized_Fst_Method1
    ),
    normalizedTransform = ifelse(normalizedTransform > 3, 3, normalizedTransform)
  )


# Select the columns from test_df that define the combinations to keep
combinations_to_keep <- TIME_MCMM_normalized_fst_multiple_colony %>%
  select(MAG_id, Colony, Comparison) %>%
  distinct() # Get unique combinations

# Filter TIME_MCMM_normalized_fst_multiple_colony based on these combinations
filtered_TIME_MCMM_normalized_fst <- TIME_MCMM_normalized_fst_multiple_colony %>%
  semi_join(combinations_to_keep, by = c("MAG_id", "Colony", "Comparison"))

# --- Original Plotting Script (using the filtered data) ---

# Define the directory to save the plots
# Consider changing the output directory name to reflect the filtering
output_directory <- "/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/multiple_colony_multiple_month/r_figures/line_plots_filtered"


# Find unique combinations of MAG_id and Colony in the FILTERED data
unique_combinations_filtered <- unique(filtered_TIME_MCMM_normalized_fst[, c("MAG_id", "Colony")])

# Define the Time values for the vertical lines
seasonal_time_points <- c(July = 67, October = 151, January = 240)

# Define the fixed x-axis limits
fixed_x_limits <- c(0, 271)

# Loop through each unique combination in the FILTERED data and generate and save a plot
for (i in 1:nrow(unique_combinations_filtered)) {
  current_mag_id <- unique_combinations_filtered$MAG_id[i]
  current_colony <- unique_combinations_filtered$Colony[i]
  
  # Filter the FILTERED data for the current combination
  # Note: This second filter is technically redundant if unique_combinations_filtered comes from the filtered data,
  # but it's kept for clarity and consistency with the original loop structure.
  subset_data <- filtered_TIME_MCMM_normalized_fst[
    filtered_TIME_MCMM_normalized_fst$MAG_id == current_mag_id &
      filtered_TIME_MCMM_normalized_fst$Colony == current_colony,
  ]
  
  # Check if subset_data is empty
  if (nrow(subset_data) == 0) {
    cat(paste("Skipping plot for MAG_id:", current_mag_id, ", Colony:", current_colony, "- No data after filtering.\n"))
    next
  }
  
  # Get the lifestyle for the subtitle (assuming it's consistent for each MAG_id-Colony combination)
  current_lifestyle <- unique(subset_data$lifestyle)
  if (length(current_lifestyle) > 1) {
    warning(paste("Multiple lifestyles found for MAG_id:", current_mag_id, ", Colony:", current_colony))
    current_lifestyle <- current_lifestyle[1] # Take the first one if there are multiple
  } else if (length(current_lifestyle) == 0) {
    current_lifestyle <- "Unknown" # Handle case where no lifestyle is found
  }
  
  
  # Generate the line plot
  plot <- ggplot(data = subset_data, aes(x = Time, y = normalizedTransform)) + # Using normalizedTransform on y-axis
    geom_line() +
    scale_y_continuous(limits = c(0, 3)) +
    scale_x_continuous(limits = fixed_x_limits) + # Setting fixed x-axis limits
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") + # Adding the horizontal line
    geom_vline(xintercept = seasonal_time_points["July"], linetype = "dotted", color = "blue") + # V-line for July
    geom_vline(xintercept = seasonal_time_points["October"], linetype = "dotted", color = "blue") + # V-line for October
    geom_vline(xintercept = seasonal_time_points["January"], linetype = "dotted", color = "blue") + # V-line for January
    labs(
      title = paste("MAG_id:", current_mag_id, ", Colony:", current_colony),
      subtitle = paste("Lifestyle:", current_lifestyle), # Adding the subtitle
      x = "Time",
      y = "Normalized Fst (Method 1)" # Keeping the y-axis label consistent
    )
  
  # Create a filename for the plot
  filename <- paste0("MAG_", gsub("[^[:alnum:]]", "_", current_mag_id),
                     "_Colony_", gsub("[^[:alnum:]]", "_", current_colony),
                     "_filtered.png") # Added _filtered to filename
  
  # Save the plot
  ggsave(filename = filename, plot = plot, path = output_directory, width = 8, height = 6, units = "in")
  
  # Optional: Print a message to the console to track progress
  cat(paste("Saved plot:", filename, "\n"))
}

cat("All filtered plots have been saved to:", output_directory, "\n")


############################################PERMUTATIONS FOR WITHIN AND BETWEEN


# Find unique lifestyles
unique_lifestyles <- unique(combined_genic_df$lifestyle)
unique_lifestyles <- unique_lifestyles[!is.na(unique_lifestyles)] # Remove potential NA lifestyles

# --- Parameters for Permutation Test ---
n_permutations <- 1000 # Number of permutations

# --- Perform Permutation Test and Plot for Each Lifestyle ---

cat("--- Permutation Test Results (Mean Comparison) and Plots ---\n")

results_list <- map_dfr(unique_lifestyles, function(current_lifestyle) {
  cat("\nAnalyzing Lifestyle:", current_lifestyle, "\n")
  
  # Filter data for the current lifestyle
  subset_data <- combined_genic_df %>%
    filter(lifestyle == current_lifestyle)
  
  # Separate data by source
  distance_data <- subset_data %>% filter(Source == "Distance") %>% pull(Fst_Method1) %>% na.omit()
  temporal_data <- subset_data %>% filter(Source == "Temporal") %>% pull(Fst_Method1) %>% na.omit()
  
  n_distance <- length(distance_data)
  n_temporal <- length(temporal_data)
  
  cat("Sample size (Distance):", n_distance, "\n")
  cat("Sample size (Temporal):", n_temporal, "\n")
  
  # Determine the smaller and larger group for the test comparison
  if (n_distance <= n_temporal) {
    smaller_group_name <- "Distance"
    larger_group_name <- "Temporal"
    smaller_data_for_test <- distance_data
    larger_data_for_test <- temporal_data
    sample_size_for_resampling <- n_distance
  } else {
    smaller_group_name <- "Temporal"
    larger_group_name <- "Distance"
    smaller_data_for_test <- temporal_data
    larger_data_for_test <- distance_data
    sample_size_for_resampling <- n_temporal
  }
  
  # Skip if sample size for resampling is zero or either group has no data
  if (sample_size_for_resampling == 0 || length(smaller_data_for_test) == 0 || length(larger_data_for_test) == 0) {
    cat("Skipping analysis and plot for lifestyle", current_lifestyle, "- Insufficient data for test.\n")
    return(tibble(lifestyle = current_lifestyle,
                  smaller_group = smaller_group_name,
                  larger_group = larger_group_name,
                  n_smaller = length(smaller_data_for_test),
                  n_larger = length(larger_data_for_test),
                  observed_mean_smaller = NA,
                  permutation_p_value = NA))
  }
  
  
  # Calculate the observed mean of the smaller group for the test
  observed_mean_smaller <- mean(smaller_data_for_test)
  
  # Perform permutations
  permutation_means <- replicate(n_permutations, {
    resample <- sample(larger_data_for_test, size = sample_size_for_resampling, replace = TRUE)
    mean(resample)
  })
  
  # Calculate the observed mean of the larger group's original data (for plotting line)
  observed_mean_larger <- mean(larger_data_for_test)
  
  # Calculate p-value (proportion of permutation means as extreme as or more extreme than observed mean)
  # Two-tailed test
  p_value <- mean(abs(permutation_means - observed_mean_larger) >= abs(observed_mean_smaller - observed_mean_larger))
  
  # If the p-value is 0 (meaning no permutation mean was more extreme), use a small value for reporting
  if (p_value == 0) {
    p_value <- 1 / (n_permutations + 1)
  }
  
  cat("Observed Mean (", smaller_group_name, "):", observed_mean_smaller, "\n")
  cat("Observed Mean (", larger_group_name, "):", observed_mean_larger, "\n")
  cat("Permutation P-value:", format.pval(p_value, digits = 3), "\n")
  
  # --- Create Plot ---
  # Prepare data for plotting (need both original Distance and Temporal data)
  plot_data <- subset_data %>%
    select(Fst_Method1, Source) %>%
    filter(!is.na(Fst_Method1)) # Remove NAs for plotting
  
  plot <- ggplot(plot_data, aes(x = Fst_Method1, fill = Source)) +
    geom_density(alpha = 0.5) + # Use density plots
    geom_vline(xintercept = median(distance_data), color = "steelblue", linetype = "dashed", linewidth = 1) + # Median of Distance
    geom_vline(xintercept = median(temporal_data), color = "pink3", linetype = "dashed", linewidth = 1) + # Median of Temporal
    geom_vline(xintercept = mean(distance_data), color = "steelblue", linetype = "solid", linewidth = 1) + # Mean of Distance
    geom_vline(xintercept = mean(temporal_data), color = "pink3", linetype = "solid", linewidth = 1) + # Mean of Temporal
    labs(
      title = paste("Fst Distribution for Lifestyle:", current_lifestyle),
      subtitle = paste("Permutation Test p =", format.pval(p_value, digits = 3)), # Display p-value
      x = "Fst_Method1",
      y = "Density",
      fill = "Data Source"
    ) +
    theme_minimal()
  
  # Print the plot
  print(plot)
  
  # Return results for this lifestyle
  tibble(lifestyle = current_lifestyle,
         smaller_group = smaller_group_name,
         larger_group = larger_group_name,
         n_smaller = length(smaller_data_for_test),
         n_larger = length(larger_data_for_test),
         observed_mean_smaller = observed_mean_smaller,
         permutation_p_value = p_value)
})

cat("\n--- Summary of Permutation Test Results ---\n")
print(results_list)




### CATEGORY ###

unique_functional_categories <- unique(combined_genic_df$functional_category)
unique_functional_categories <- unique_functional_categories[!is.na(unique_functional_categories)] # Remove potential NA functional categories

# --- Parameters for Permutation Test ---
n_permutations <- 1000 # Number of permutations

# --- Perform Permutation Test and Plot for Each Functional Category ---

cat("--- Permutation Test Results (Mean Comparison) by Functional Category and Plots ---\n")

results_list_func_cat <- map_dfr(unique_functional_categories, function(current_functional_category) {
  cat("\nAnalyzing Functional Category:", current_functional_category, "\n")
  
  # Filter data for the current functional category
  subset_data <- combined_genic_df %>%
    filter(functional_category == current_functional_category)
  
  # Separate data by source
  distance_data <- subset_data %>% filter(Source == "Distance") %>% pull(Fst_Method1) %>% na.omit()
  temporal_data <- subset_data %>% filter(Source == "Temporal") %>% pull(Fst_Method1) %>% na.omit()
  
  n_distance <- length(distance_data)
  n_temporal <- length(temporal_data)
  
  cat("Sample size (Distance):", n_distance, "\n")
  cat("Sample size (Temporal):", n_temporal, "\n")
  
  # Determine the smaller and larger group for the test comparison
  if (n_distance <= n_temporal) {
    smaller_group_name <- "Distance"
    larger_group_name <- "Temporal"
    smaller_data_for_test <- distance_data
    larger_data_for_test <- temporal_data
    sample_size_for_resampling <- n_distance
  } else {
    smaller_group_name <- "Temporal"
    larger_group_name <- "Distance"
    smaller_data_for_test <- temporal_data
    larger_data_for_test <- distance_data
    sample_size_for_resampling <- n_temporal
  }
  
  # Skip if sample size for resampling is zero or either group has no data
  if (sample_size_for_resampling == 0 || length(smaller_data_for_test) == 0 || length(larger_data_for_test) == 0) {
    cat("Skipping analysis and plot for functional category", current_functional_category, "- Insufficient data for test.\n")
    return(tibble(functional_category = current_functional_category,
                  smaller_group = smaller_group_name,
                  larger_group = larger_group_name,
                  n_smaller = length(smaller_data_for_test),
                  n_larger = length(larger_data_for_test),
                  observed_mean_smaller = NA,
                  permutation_p_value = NA))
  }
  
  
  # Calculate the observed mean of the smaller group for the test
  observed_mean_smaller <- mean(smaller_data_for_test)
  
  # Perform permutations
  permutation_means <- replicate(n_permutations, {
    resample <- sample(larger_data_for_test, size = sample_size_for_resampling, replace = TRUE)
    mean(resample)
  })
  
  # Calculate the observed mean of the larger group's original data (for plotting line)
  observed_mean_larger <- mean(larger_data_for_test)
  
  # Calculate p-value (proportion of permutation means as extreme as or more extreme than observed mean)
  # Two-tailed test
  p_value <- mean(abs(permutation_means - observed_mean_larger) >= abs(observed_mean_smaller - observed_mean_larger))
  
  # If the p-value is 0 (meaning no permutation mean was more extreme), use a small value for reporting
  if (p_value == 0) {
    p_value <- 1 / (n_permutations + 1)
  }
  
  cat("Observed Mean (", smaller_group_name, "):", observed_mean_smaller, "\n")
  cat("Observed Mean (", larger_group_name, "):", observed_mean_larger, "\n")
  cat("Permutation P-value:", format.pval(p_value, digits = 3), "\n")
  
  # --- Create Plot ---
  # Prepare data for plotting (need both original Distance and Temporal data)
  plot_data <- subset_data %>%
    select(Fst_Method1, Source) %>%
    filter(!is.na(Fst_Method1)) # Remove NAs for plotting
  
  plot <- ggplot(plot_data, aes(x = Fst_Method1, fill = Source)) +
    geom_density(alpha = 0.5) + # Use density plots
    geom_vline(xintercept = median(distance_data), color = "steelblue", linetype = "dashed", linewidth = 1) + # Median of Distance
    geom_vline(xintercept = median(temporal_data), color = "pink3", linetype = "dashed", linewidth = 1) + # Median of Temporal
    geom_vline(xintercept = mean(distance_data), color = "steelblue", linetype = "solid", linewidth = 1) + # Mean of Distance
    geom_vline(xintercept = mean(temporal_data), color = "pink3", linetype = "solid", linewidth = 1) + # Mean of Temporal
    labs(
      title = paste("Fst Distribution for Functional Category:", current_functional_category),
      subtitle = paste("Permutation Test p =", format.pval(p_value, digits = 3)), # Display p-value
      x = "Fst_Method1",
      y = "Density",
      fill = "Data Source"
    ) +
    theme_minimal()
  
  # Print the plot
  print(plot)
  
  # Return results for this functional category
  tibble(functional_category = current_functional_category,
         smaller_group = smaller_group_name,
         larger_group = larger_group_name,
         n_smaller = length(smaller_data_for_test),
         n_larger = length(larger_data_for_test),
         observed_mean_smaller = observed_mean_smaller,
         permutation_p_value = p_value)
})

cat("\n--- Summary of Permutation Test Results by Functional Category ---\n")
print(results_list_func_cat)

ggplot(genic_im, aes(x = Fst_Method1, fill = lifestyle)) +
  geom_density(alpha = 0.5) + # Use density plots

  labs(
    title = paste("Fst Distribution for Functional Category:", current_functional_category),
    x = "Fst_Method1",
    y = "Density",
    fill = "Data Source"
  ) +
  theme_minimal()

##################### POPULATION GENETIC STATISTICS OVER TIME ###################

# --- Important: Make sure 'TIME_MCMM_normalized_fst_multiple_colony' dataframe is loaded in your R environment. ---
# For demonstration purposes, if you don't have this dataframe loaded,
# you can uncomment and run the following lines to create a dummy dataset:
# set.seed(123) # For reproducibility
# TIME_MCMM_normalized_fst_multiple_colony <- data.frame(
#   Time = rep(seq(1, 10, by = 1), each = 100),
#   lifestyle = sample(c("Lytic", "Prophage", "Temperate", "Lysogenic"), 1000, replace = TRUE),
#   Normalized_Fst_Method1 = runif(1000, 0.05, 0.8) + rep(c(0, 0.1, 0.2, -0.05), each = 250) * rep(seq(1, 10, by = 1), each = 100)[1:1000]
# )
# # Ensure no Fst values are negative if runif produces slightly below 0 due to addition/subtraction
# TIME_MCMM_normalized_fst_multiple_colony$Normalized_Fst_Method1[
#   TIME_MCMM_normalized_fst_multiple_colony$Normalized_Fst_Method1 < 0
#   ] <- 0.01


# Get all unique lifestyle values from your dataframe
unique_lifestyles <- unique(TIME_MCMM_normalized_fst_multiple_colony$lifestyle)

TIME_MCMM_normalized_fst_multiple_colony$Normalized_Fst_Method1 <-
  pmin(TIME_MCMM_normalized_fst_multiple_colony$Normalized_Fst_Method1, 3)
# Loop through each unique lifestyle to generate and display its plot
for (ls in unique_lifestyles) {
  
  # Filter the data for the current lifestyle
  current_lifestyle_data <- TIME_MCMM_normalized_fst_multiple_colony %>%
    filter(lifestyle == ls)
  
  # Calculate the mean and median of Normalized_Fst_Method1 for each 'Time' point
  summary_data <- current_lifestyle_data %>%
    group_by(Time) %>%
    summarise(
      mean_fst = mean(Normalized_Fst_Method1, na.rm = TRUE),
      median_fst = median(Normalized_Fst_Method1, na.rm = TRUE),
      .groups = 'drop' # Ensures the output is a simple dataframe
    )
  
  # Create the ggplot for the current lifestyle
  p <- ggplot(current_lifestyle_data, aes(x = Time, y = Normalized_Fst_Method1)) +
    # Add raw data points with high transparency to show density without overplotting
    geom_point(alpha = 0.05, size = 0.5, color = "grey50") + # Faint grey points
    
    # Add the mean trend line
    geom_line(data = summary_data, aes(y = mean_fst, color = "Mean"), linewidth = 1) +
    
    # Add the median trend line (dashed for distinction)
    geom_line(data = summary_data, aes(y = median_fst, color = "Median"), linewidth = 1, linetype = "dashed") +
    
    # Add the horizontal reference line at Normalized_Fst_Method1 == 1
    geom_hline(yintercept = 1, linetype = "dotted", color = "red", linewidth = 1) +
    
    # Define labels and titles for the plot
    labs(
      title = paste("Normalized Fst Method1 Over Time for Lifestyle:", ls),
      subtitle = "Mean and Median Trends with Reference Line at 1",
      x = "Time",
      y = "Normalized Fst Method1",
      color = "Statistic" # Legend title for mean/median lines
    ) +
    
    # Manually assign colors for the mean and median lines for clarity
    scale_color_manual(values = c("Mean" = "steelblue", "Median" = "darkgreen")) +
    
    # Set the y-axis scale. Fst values are typically between 0 and 1,
    # so a linear scale is usually appropriate unless extreme skew is present.
    # The limits are set to ensure the reference line at 1 is always visible,
    # and to start from 0 or the minimum observed value.

    
    # Apply a clean, minimal theme and further aesthetic adjustments
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # Center and bold main title
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"), # Center and lighten subtitle
      axis.title = element_text(size = 10), # Axis title font size
      axis.text = element_text(size = 9), # Axis text font size
      legend.title = element_text(size = 10, face = "bold"), # Legend title styling
      legend.text = element_text(size = 9), # Legend item text size
      legend.position = "right", # Position the legend to the right of the plot
      panel.grid.major = element_line(color = "grey90", linewidth = 0.15), # Lighter major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines for a cleaner look
      panel.background = element_rect(fill = "white", color = NA), # White plot panel background
      plot.background = element_rect(fill = "white", color = NA) # White overall plot background
  )+
    scale_y_continuous(limits = c(0, 3)) 
  # Print the plot to display it. If running in an interactive R environment
  # (like RStudio), this will display each plot as it's generated.
  print(p)
}



