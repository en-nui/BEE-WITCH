library(ggplot2)
library(dplyr)
library(readr)
library(broom)
library(purrr)
library(changepoint)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyr)
library(patchwork)
library(vegan)
library(lme4)
library(lmerTest)
library(ggeffects)
library(glmmTMB)
library(mclust)
library(gt)
library(broom.mixed)
library(MuMIn)
library(DHARMa)
#--- 2. Define Gruvbox Color Palette ---
  # We'll define a few key colors to use in our plot
  gruvbox_palette <- c(
    dark_bg    = "#282828", # A dark gray that will look black
    light_fg   = "#ebdbb2",
    faded_red  = "#cc241d",
    bright_red = "#fb4934",
    blue       = "#83a598",
    aqua       = "#8ec07c",
    orange     = "#fe8019"
  )
ecology_df = read.csv('/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv')
nucdiv_df = read.csv('/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/popgen_summaries_final_filtered/colony_103_filtered_snps_cleaned_pooled_by_gene.csv')
# b. Filter ecology data for the target colony
ecology_col101 <- ecology_df %>%
  filter(Colony == 106)

# c. Filter nucdiv data for the target MAG
nucdiv_mag173 <- nucdiv_df %>%
  filter(MAG_id == "mmag_254")

# d. Define the correct month order and apply it
month_order <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")
nucdiv_mag173$Month <- factor(nucdiv_mag173$Month, levels = month_order)

# --- 3. Calculate the Monthly Average 'pi' ---
avg_pi_mag173 <- nucdiv_mag173 %>%
  group_by(Month) %>%
  summarise(average_pi = mean(pi, na.rm = TRUE))

# --- 4. Create the Plot ---
ggplot() +
  # Plot the thin blue lines for each unique gene
  geom_line(
    data = nucdiv_mag173, 
    aes(x = Month, y = pi, group = unique_gene_id), 
    color = "blue", 
    alpha = 0.3
  ) +
  # Plot the thick black line for the monthly average
  geom_line(
    data = avg_pi_mag173, 
    aes(x = Month, y = average_pi, group = 1), # group = 1 is essential for a single line
    color = "black", 
    linewidth = 1.2
  ) +
  labs(
    title = "Nucleotide Diversity (pi) for MAG mmag_173 in Colony 101",
    subtitle = "Blue lines are individual genes; Black line is the monthly average",
    x = "Month",
    y = "Nucleotide Diversity (pi)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(unique(nucdiv_df$MAG_id))
#mmag_72 is interesting


head(ecology_df)












################ SCRIPT GEN ############
# --- 2. Load Data ---
# Define the path to your metadata file
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"

# Load the data
ecology_df <- read_csv(metadata_file)
ecology_df <- ecology_df[ecology_df$type == "mmag",]
ecology_df <- ecology_df[ecology_df$Caste == "Worker",]

# Step A: Aggregate technical replicates (sample_name) to the biological sample level (Colony, Month)
mag_abundances <- ecology_df %>%
  group_by(Colony, Month, MAG_id, Genus) %>%
  summarise(
    # Calculate the average coverage for each MAG in each biological sample
    mean_coverage = mean(mean, na.rm = TRUE),
    .groups = 'drop'
  )

# Step B: Calculate MAG-level relative abundance
# This adds the 'rel_abund_mag' column
mag_abundances <- mag_abundances %>%
  group_by(Colony, Month) %>%
  mutate(
    # For each sample, calculate the total coverage of all MAGs
    total_coverage_in_sample = sum(mean_coverage),
    # Divide the MAG's coverage by the total to get its relative abundance
    rel_abund_mag = mean_coverage / total_coverage_in_sample
  ) %>%
  ungroup() %>%
  # Remove the temporary total coverage column
  select(-total_coverage_in_sample)

# Step C: Calculate Genus-level relative abundance
# This creates a new summary table for genera
genus_abundances <- mag_abundances %>%
  group_by(Colony, Month, Genus) %>%
  # For each Genus, sum the relative abundances of all its member MAGs
  summarise(rel_abund_genus = sum(rel_abund_mag), .groups = 'drop')

# Step D: Join the genus-level data back to the main table
final_abundances <- left_join(mag_abundances, genus_abundances, by = c("Colony", "Month", "Genus"))





# --- 2. Load and Prepare All Necessary Data ---

# a. Load and combine all GENE-level population genetics summaries
popgen_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/popgen_summaries_final_filtered"
popgen_files <- list.files(path = popgen_dir, pattern = "_by_gene.csv", full.names = TRUE, recursive = TRUE)

gene_level_data <- popgen_files %>%
  map_dfr(read_csv)

# b. Aggregate gene-level data to the MAG-level
mag_level_popgen <- gene_level_data %>%
  # MODIFIED: Removed Genus from this initial grouping
  group_by(MAG_id, Colony, Month) %>% 
  summarise(
    mean_pi = mean(pi, na.rm = TRUE),
    mean_Watterson_theta = mean(Watterson_theta_per_bp, na.rm = TRUE),
    mean_TajimaD = mean(TajimaD, na.rm = TRUE),
    .groups = 'drop'
  )

# c. Add a numeric Time variable for regression
month_order <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")
mag_level_popgen$Time <- match(mag_level_popgen$Month, month_order)

# d. Load ecology data and create a map from MAG_id to Genus
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"
ecology_df <- read_csv(metadata_file)

# --- NEW: Create a simple lookup table for Genus ---
genus_map <- ecology_df %>%
  distinct(MAG_id, Genus)

# --- NEW: Join the Genus information to the popgen data ---
mag_level_popgen <- left_join(mag_level_popgen, genus_map, by = "MAG_id")

# e. Load and calculate relative abundance data
abundance_data <- ecology_df %>%
  group_by(Colony, Month, MAG_id) %>%
  summarise(mean_coverage = mean(mean, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Colony, Month) %>%
  mutate(rel_abund_mag = mean_coverage / sum(mean_coverage)) %>%
  ungroup() %>%
  select(Colony, Month, MAG_id, rel_abund_mag)

# f. Merge the datasets to create the final analysis table
full_data <- inner_join(mag_level_popgen, abundance_data, by = c("Colony", "Month", "MAG_id"))
add_cumulative_time_from_month_nucdiv <- function(data) {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34, 
                        "September" = 28, "October" = 22, "November" = 42, 
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  data$Time <- month_to_cumulative_time_map[data$Month]
  return(data)
}
full_data <- add_cumulative_time_from_month_nucdiv(full_data)
# Filter out any groups with too few timepoints for analysis
full_data <- full_data %>%
  group_by(Colony, MAG_id) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("--- Data Loading and Aggregation Complete ---\n")
cat("First few rows of the final data for analysis:\n\n")
print(head(full_data))
print(unique(full_data$MAG_id))

# --- 3. Analysis 1: Linear Regression (Stability vs. Shift) ---
cat("\n--- 1. Linear Regression (Mean Diversity vs. Time) ---\n")

regression_results <- full_data %>%
  group_by(Colony, MAG_id, Genus) %>%
  do(broom::tidy(lm(mean_pi ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p_adj)

print(head(regression_results))

# --- 4. Analysis 2: Change Point Analysis (MODIFIED) ---
cat("\n--- 2. Change Point Analysis (Abrupt Shifts in Mean Diversity) ---\n")

changepoint_results <- full_data %>%
  arrange(Colony, MAG_id, Time) %>%
  group_by(Colony, MAG_id, Genus) %>%
  summarise(
    cpt_object = list(tryCatch({
      cpt.mean(mean_pi, method = "PELT", minseglen = 2)
    }, error = function(e) {
      NA
    })),
    .groups = 'drop'
  ) %>%
  filter(!is.na(cpt_object)) %>%
  # --- MODIFIED: Replaced npts() with the correct nseg() - 1 logic ---
  mutate(num_changepoints = purrr::map_int(cpt_object, ~(nseg(.x) - 1))) %>%
  select(-cpt_object) %>%
  filter(num_changepoints > 0) %>%
  arrange(desc(num_changepoints))


# --- 5. Analysis 3: ANCOVA (Comparing Species within a Genus) ---
cat("\n--- 3. ANCOVA (Comparing Mean Diversity Trends within Genera) ---\n")

ancova_results <- full_data %>%
  group_by(Colony, Genus) %>%
  filter(n_distinct(MAG_id) > 1) %>%
  do(broom::tidy(aov(mean_pi ~ Time * MAG_id, data = .))) %>%
  ungroup() %>%
  filter(grepl("Time:MAG_id", term)) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p_adj)

print(head(ancova_results))






# --- 2. Consolidated Data Preparation ---

# a. Load and filter ecology data
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"
ecology_df <- read_csv(metadata_file) %>%
  filter(type == "mmag", Caste == "Worker")

# b. Calculate mean coverage and relative abundance
abundance_data <- ecology_df %>%
  group_by(Colony, Month, MAG_id, Genus, Season) %>%
  summarise(mean_coverage = mean(mean, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Colony, Month) %>%
  mutate(rel_abund_mag = mean_coverage / sum(mean_coverage)) %>%
  ungroup()

# c. Load and aggregate population genetics data
popgen_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/popgen_summaries_final_filtered"
popgen_files <- list.files(path = popgen_dir, pattern = "_by_gene.csv", full.names = TRUE, recursive = TRUE)

gene_level_data <- popgen_files %>%
  map_dfr(read_csv)

mag_level_popgen <- gene_level_data %>%
  group_by(MAG_id, Colony, Month) %>% 
  summarise(
    mean_pi = mean(pi, na.rm = TRUE),
    mean_Watterson_theta = mean(Watterson_theta_per_bp, na.rm = TRUE),
    mean_TajimaD = mean(TajimaD, na.rm = TRUE),
    .groups = 'drop'
  )

# d. Merge pop-gen and abundance data into one final table
full_data <- inner_join(
  mag_level_popgen, 
  abundance_data, 
  by = c("Colony", "Month", "MAG_id")
) %>%
  filter(!is.na(Season)) # Ensure Season is not NA

# Convert Season to a factor for modeling
full_data$Season <- as.factor(full_data$Season)
add_cumulative_time_from_month_nucdiv <- function(data) {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34, 
                        "September" = 28, "October" = 22, "November" = 42, 
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  data$Time <- month_to_cumulative_time_map[data$Month]
  return(data)
}
full_data <- add_cumulative_time_from_month_nucdiv(full_data)
cat("--- Data Preparation Complete ---\n")
cat("First few rows of the final analysis table:\n\n")
print(head(full_data))


# --- 3. Analysis 1: Test for Effect of Season (LMM) ---
cat("\n--- 1. Testing for Significant Seasonal Effects ---\n")

# Model 1: Test if Season affects nucleotide diversity (mean_pi)
# We use a log transform because pi is always positive and often skewed
# Random effects for Colony and MAG_id account for repeated measures
lmm_pi <- lmer(log(mean_pi + 1e-9) ~ Season + (1 | Colony) + (1 | MAG_id), data = full_data)

cat("\n--- LMM Summary: Effect of Season on Nucleotide Diversity (pi) ---\n")
print(summary(lmm_pi))


# Model 2: Test if Season affects absolute abundance (mean_coverage)
lmm_abund <- lmer(log(mean_coverage + 1) ~ Season + (1 | Colony) + (1 | MAG_id), data = full_data)

cat("\n--- LMM Summary: Effect of Season on Absolute Abundance (Mean Coverage) ---\n")
print(summary(lmm_abund))


# --- 4. Analysis 2: Correlation of Diversity and Relative Abundance ---
cat("\n--- 2. Pearson Correlation (Diversity vs. Relative Abundance) ---\n")
# Note: This is performed per-MAG, across months

correlation_results <- full_data %>%
  group_by(Colony, MAG_id, Genus) %>%
  filter(n() > 2) %>% 
  summarise(
    pearson_cor = cor(mean_pi, rel_abund_mag, method = "pearson"),
    n_months = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(abs(pearson_cor)))

print(head(correlation_results))
hist(correlation_results$pearson_cor)



# Extract the random effects for MAG_id from the nucleotide diversity model
mag_effects <- ranef(lmm_pi)$MAG_id %>%
  as.data.frame() %>%
  tibble::rownames_to_column("MAG_id") %>%
  rename(intercept_deviation = `(Intercept)`) %>%
  arrange(desc(abs(intercept_deviation)))

cat("--- Top MAGs with the Largest Deviations from Mean Diversity ---\n")
print(head(mag_effects))



# --- 4. Correlation Follow-up: Detrending the Data (MODIFIED) ---

# Step A: For each MAG, calculate the residuals after removing the effect of Time
detrended_data <- full_data %>%
  group_by(Colony, MAG_id) %>%
  filter(n() > 2) %>% # Ensure there are enough points for regression
  # --- MODIFIED: Replaced mutate() with do() for this operation ---
  do({
    # '.' refers to the data for the current group
    group_data <- .
    
    # Fit models for the current group
    pi_model <- lm(mean_pi ~ Time, data = group_data)
    abund_model <- lm(rel_abund_mag ~ Time, data = group_data)
    
    # Add residuals back to the group's data
    group_data$pi_residuals <- resid(pi_model)
    group_data$abund_residuals <- resid(abund_model)
    
    # Return the modified data for the group
    group_data
  }) %>%
  ungroup()

# Step B: Calculate the correlation on the DETRENDED data
detrended_correlation_results <- detrended_data %>%
  group_by(Colony, MAG_id, Genus) %>%
  summarise(
    pearson_cor_detrended = cor(pi_residuals, abund_residuals, method = "pearson"),
    n_months = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(abs(pearson_cor_detrended)))

cat("\n--- Detrended Correlation (Diversity vs. Relative Abundance) ---\n")
print(detrended_correlation_results,n=150)
hist(detrended_correlation_results$pearson_cor_detrended)





# --- 2. Load and Prepare Data ---
# (This section is unchanged)
# Assuming 'full_data' is in your environment

# --- 3. Identify MAGs with a Significant Seasonal Pattern (MODIFIED) ---

# a. Define a function to run the model and extract the p-value for Season
test_season_effect <- function(data, metric_name) {
  if (nrow(data) < 4 || n_distinct(data$Season) < 2) {
    return(NA_real_)
  }
  formula <- as.formula(paste(metric_name, "~ Season"))
  model <- lm(formula, data = data)
  summary_stats <- broom::glance(model)
  return(summary_stats$p.value)
}

# b. Run the analysis for each MAG
mag_seasonal_summary <- full_data %>%
  group_by(MAG_id) %>%
  nest() %>%
  mutate(
    p_value_pi = map_dbl(data, ~test_season_effect(.x, "mean_pi")),
    p_value_abund = map_dbl(data, ~test_season_effect(.x, "rel_abund_mag"))
  ) %>%
  select(-data) %>%
  ungroup()

# --- NEW: Apply FDR (Benjamini-Hochberg) correction to the p-values ---
mag_seasonal_summary <- mag_seasonal_summary %>%
  mutate(
    p_adj_pi = p.adjust(p_value_pi, method = "BH"),
    p_adj_abund = p.adjust(p_value_abund, method = "BH")
  )

# c. Filter for MAGs with significant seasonal nucleotide diversity (using adjusted p-value)
seasonal_pi_mags <- mag_seasonal_summary %>%
  filter(p_adj_pi < 0.1) %>%
  arrange(p_adj_pi)

cat("--- MAGs with Significant Seasonal Differences in Nucleotide Diversity (pi) ---\n")
print(seasonal_pi_mags)


# d. Filter for MAGs with significant seasonal relative abundance (using adjusted p-value)
seasonal_abund_mags <- mag_seasonal_summary %>%
  filter(p_adj_abund < 0.1) %>%
  arrange(p_adj_abund)

cat("\n--- MAGs with Significant Seasonal Differences in Relative Abundance ---\n")
print(seasonal_abund_mags)



correlation_p_values <- detrended_data %>%
  group_by(Colony, MAG_id, Genus) %>%
  # Use do() with broom::tidy to run cor.test() on each group
  do(broom::tidy(cor.test(~ pi_residuals + abund_residuals, data = .))) %>%
  ungroup() %>%
  # Rename columns for clarity and add an adjusted p-value
  rename(pearson_cor_detrended = estimate) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p_adj)

cat("\n--- Detrended Correlation with p-values ---\n")
print(head(correlation_p_values))

# Now you can filter based on significance
significant_correlations <- correlation_p_values %>%
  filter(p_adj < 0.05)

cat("\n--- Only Significant Correlations (p_adj < 0.05) ---\n")
print(significant_correlations)



# --- 2. Prepare Data for Plotting ---
# Assuming 'correlation_p_values' dataframe from our previous step is in your environment.
# If not, you would first need to run the code to generate it.

# Create the categories for plotting using the thresholds
plot_summary <- correlation_p_values %>%
  mutate(
    # Use case_when for clear, conditional categorization
    correlation_type = case_when(
      p_adj < 0.05 & pearson_cor_detrended >= 0.1  ~ "Significant Positive",
      p_adj < 0.05 & pearson_cor_detrended <= -0.1 ~ "Significant Negative",
      p_adj >= 0.05 & pearson_cor_detrended >= 0.1  ~ "Not Significant Positive",
      p_adj >= 0.05 & pearson_cor_detrended <= -0.1 ~ "Not Significant Negative",
      TRUE                                       ~ "No Correlation" # Catches everything between -0.1 and 0.1
    )
  ) %>%
  # Count the number of MAGs in each category for each colony
  group_by(Colony, correlation_type) %>%
  summarise(n = n(), .groups = 'drop') %>%
  # Calculate the proportion of each category within each colony
  group_by(Colony) %>%
  mutate(proportion = n / sum(n))

# --- 3. Create the Figure ---
gruvbox_palette <- c(
  dark_bg    = "#282828", # A dark gray that will look black
  light_fg   = "#ebdbb2",
  faded_red  = "#cc241d",
  bright_red = "#fb4934",
  blue       = "#83a598",
  aqua       = "#8ec07c",
  orange     = "#fe8019"
)
# Define a custom color palette to make significant results stand out
color_palette <- c(
  "Significant Positive" = "#cc241d",      # Strong Red
  "Not Significant Positive" = "#ebdbb2",  # Lighter Orange/Red
  "Significant Negative" = "#689d6a",      # Strong Blue
  "Not Significant Negative" = "#665c54",  # Lighter Blue
  "No Correlation" = "#32302f"                # Light Grey
)

# Ensure the categories are ordered logically in the plot legend
plot_summary$correlation_type <- factor(plot_summary$correlation_type, levels = c(
  "Significant Positive", "Not Significant Positive", 
  "No Correlation", 
  "Not Significant Negative", "Significant Negative"
))

# Build the ggplot
ggplot(plot_summary, aes(x = as.factor(Colony), y = proportion, fill = correlation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  # Use the custom color palette
  scale_fill_manual(values = color_palette) +
  # Format the y-axis as percentages
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Correlations between Diversity and Abundance",
    subtitle = "Based on detrended time-series data for each MAG",
    x = "Colony",
    y = "Proportion of MAGs",
    fill = "Correlation Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels if needed

























# --- 2. Prepare Data for Plotting ---
# Assuming 'correlation_p_values' dataframe from our previous step is in your environment.

# Define the list of core genera you want to include in the plot
core_genera <- c(
  "bartonella", "apilactobacillus", "bifidobacterium", "bombilactobacillus",
  "lactobacillus", "frischella", "gilliamella", "snodgrassella"
)

# Create the categories for plotting, this time grouping by Genus
plot_summary_genus <- correlation_p_values %>%
  # --- MODIFIED: Filter to include only the specified core genera ---
  filter(Genus %in% core_genera) %>%
  mutate(
    correlation_type = case_when(
      p_adj < 0.05 & pearson_cor_detrended >= 0.1  ~ "Significant Positive",
      p_adj < 0.05 & pearson_cor_detrended <= -0.1 ~ "Significant Negative",
      p_adj >= 0.05 & pearson_cor_detrended >= 0.1  ~ "Not Significant Positive",
      p_adj >= 0.05 & pearson_cor_detrended <= -0.1 ~ "Not Significant Negative",
      TRUE                                       ~ "No Correlation"
    )
  ) %>%
  # --- MODIFIED: Group by Genus instead of Colony ---
  group_by(Genus, correlation_type) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Genus) %>%
  mutate(proportion = n / sum(n))

# --- 3. Create the Figure ---

# Define the same custom color palette
color_palette <- c(
  "Significant Positive" = "#cc241d",      # Strong Red
  "Not Significant Positive" = "#ebdbb2",  # Lighter Orange/Red
  "Significant Negative" = "#689d6a",      # Strong Blue
  "Not Significant Negative" = "#665c54",  # Lighter Blue
  "No Correlation" = "#32302f"                # Light Grey
)

# Ensure the categories are ordered logically in the plot legend
plot_summary_genus$correlation_type <- factor(plot_summary_genus$correlation_type, levels = c(
  "Significant Positive", "Not Significant Positive", 
  "No Correlation", 
  "Not Significant Negative", "Significant Negative"
))

# Build the ggplot, now with Genus on the x-axis
ggplot(plot_summary_genus, aes(x = Genus, y = proportion, fill = correlation_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Correlations between Diversity and Abundance by Genus",
    subtitle = "Based on detrended time-series data for each MAG",
    x = "Genus", # Updated label
    y = "Proportion of MAGs",
    fill = "Correlation Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Keep rotation for long names
 

# --- 2. Define Target and Load Data ---
mag_to_plot <- "mmag_691"
# Assuming 'gene_level_data', 'ecology_df', and 'full_data' are in your environment

# --- 3. Prepare Combined Data for Plotting ---

# a. Prepare Abundance Data (Replicate and Average)
# This part is unchanged
replicate_abund <- ecology_df %>%
  filter(Caste == "Worker") %>%
  group_by(Colony, Month, sample_name) %>%
  summarise(total_rep_coverage = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
  right_join(ecology_df, by = c("Colony", "Month", "sample_name")) %>%
  filter(MAG_id == mag_to_plot) %>%
  mutate(Value = mean / total_rep_coverage) %>%
  select(Colony, Month, Value) %>%
  mutate(Metric = "Relative Abundance", Type = "Replicate")

avg_abund <- full_data %>%
  filter(MAG_id == mag_to_plot) %>%
  select(Colony, Month, rel_abund_mag) %>%
  mutate(Metric = "Relative Abundance", Type = "Average") %>%
  rename(Value = rel_abund_mag)

# b. Prepare Diversity Data (Gene and Average)
# --- MODIFIED: We now include 'unique_gene_id' to group the lines ---
gene_pi <- gene_level_data %>%
  filter(MAG_id == mag_to_plot) %>%
  select(Colony, Month, pi, unique_gene_id) %>% # Added unique_gene_id
  mutate(Metric = "Nucleotide Diversity (π)", Type = "Gene") %>%
  rename(Value = pi)

avg_pi <- full_data %>%
  filter(MAG_id == mag_to_plot) %>%
  select(Colony, Month, mean_pi) %>%
  mutate(Metric = "Nucleotide Diversity (π)", Type = "Average") %>%
  rename(Value = mean_pi)

# c. Bind all data into one long-format dataframe
final_plot_data <- bind_rows(replicate_abund, avg_abund, gene_pi, avg_pi) %>%
  mutate(Time = match(Month, month_order)) %>%
  filter(!is.na(Time))

# --- 4. Build the Final Plot ---
ggplot() +
  
  # --- Top Panel: Abundance ---
  # Plot points for each replicate
  geom_point(data = final_plot_data %>% filter(Metric == "Relative Abundance", Type == "Replicate"),
             aes(x = Time, y = Value), color = "skyblue", alpha = 0.5, size = 2) +
  # Plot the average line for abundance
  geom_line(data = final_plot_data %>% filter(Metric == "Relative Abundance", Type == "Average"),
            aes(x = Time, y = Value), color = "navy", linewidth = 1.2) +
  
  # --- Bottom Panel: Diversity ---
  # --- MODIFIED: Plot lines for each individual gene ---
  geom_line(data = final_plot_data %>% filter(Metric == "Nucleotide Diversity (π)", Type == "Gene"),
            aes(x = Time, y = Value, group = unique_gene_id), # Group by gene ID to connect correctly
            color = "lightcoral", alpha = 0.2, linewidth = 0.5) +
  # Plot the average line for diversity
  geom_line(data = final_plot_data %>% filter(Metric == "Nucleotide Diversity (π)", Type == "Average"),
            aes(x = Time, y = Value), color = "darkred", linewidth = 1.2) +
  
  # Create the grid of plots, faceted by Metric (row) and Colony (column)
  facet_grid(Metric ~ Colony, scales = "free_y") +
  
  # Customize scales and labels
  scale_x_continuous(
    breaks = unique(final_plot_data$Time),
    labels = unique(final_plot_data$Month)
  ) +
  labs(
    title = paste("Temporal Dynamics of", mag_to_plot),
    subtitle = "Average (dark line) vs. Raw Data (faint points/lines)",
    x = "Month",
    y = "Value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"), # Facet label background
    legend.position = "none"
  )


# --- 3. Define Target and File Paths ---
# (This section is unchanged)
mag_to_plot <- "mmag_691"
month_order <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")

popgen_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/popgen_summaries_final_filtered"
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"


# --- 4. Load and Prepare Data ---
# (This section is unchanged)
popgen_files <- list.files(path = popgen_dir, pattern = "_by_gene.csv", full.names = TRUE, recursive = TRUE)
gene_level_data <- popgen_files %>% map_dfr(read_csv)

ecology_df <- read_csv(metadata_file)

mag_level_popgen <- gene_level_data %>%
  group_by(MAG_id, Colony, Month) %>% 
  summarise(mean_pi = mean(pi, na.rm = TRUE), .groups = 'drop')

full_data <- mag_level_popgen %>%
  mutate(Time = match(Month, month_order))

replicate_totals <- ecology_df %>%
  filter(Caste == "Worker") %>%
  group_by(sample_name) %>%
  summarise(total_coverage = sum(mean, na.rm = TRUE), .groups = 'drop')

target_mag_replicates <- ecology_df %>%
  filter(MAG_id == mag_to_plot, Caste == "Worker") %>%
  left_join(replicate_totals, by = "sample_name") %>%
  mutate(rep_rel_abund = mean / total_coverage)

abundance_summary <- target_mag_replicates %>%
  group_by(Colony, Month) %>%
  summarise(
    mean_abund = mean(rep_rel_abund, na.rm = TRUE),
    sd_abund = sd(rep_rel_abund, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(Time = match(Month, month_order)) %>%
  filter(!is.na(Time))

gene_pi_data <- gene_level_data %>%
  filter(MAG_id == mag_to_plot) %>%
  mutate(Time = match(Month, month_order)) %>%
  filter(!is.na(Time))

avg_pi_data <- full_data %>%
  filter(MAG_id == mag_to_plot, !is.na(Time))

# --- 5. Build the Plots with Gruvbox Colors ---

# PLOT 1: Abundance
p_abund <- ggplot(abundance_summary, aes(x = Time, y = mean_abund)) +
  geom_ribbon(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund),
              fill = gruvbox_palette["orange"], alpha = 0.5) +
  geom_line(color = gruvbox_palette["dark_bg"], linewidth = 1) +
  geom_point(color = gruvbox_palette["dark_bg"], size = 2) +
  facet_wrap(~ Colony, nrow = 1) +
  labs(y = "Relative Abundance") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# PLOT 2: Diversity
p_pi <- ggplot() +
  geom_line(data = gene_pi_data,
            aes(x = Time, y = pi, group = unique_gene_id),
            color = gruvbox_palette["blue"], alpha = 0.15) +
  geom_line(data = avg_pi_data,
            aes(x = Time, y = mean_pi),
            color = gruvbox_palette["dark_bg"], linewidth = 2) +
  facet_wrap(~ Colony, nrow = 1) +
  labs(x = "Month", y = "Nucleotide Diversity (π)") +
  scale_x_continuous(
    breaks = avg_pi_data$Time,
    labels = avg_pi_data$Month
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 6. Combine Plots ---
final_plot <- p_abund / p_pi +
  plot_annotation(
    title = paste("Temporal Dynamics of", mag_to_plot),
    subtitle = "Top: Mean Abundance +/- 1 SD of Replicates. Bottom: Mean Diversity vs. Individual Gene Diversity."
  )

# Display the final combined plot
final_plot







#### BEGIN ITERATING!!! ####
# --- 2. Define Palettes, Paths, and Constants ---
# (This section is unchanged)
gruvbox_palette <- c(
  dark_bg = "#282828", blue = "#83a598", red = "#fb4934", orange = "#fe8019"
)
month_order <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")
popgen_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/popgen_summaries_final_filtered"
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"
base_output_dir_png <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/figures/abundance_pi_relationships_png"
base_output_dir_svg <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/figures/abundance_pi_relationships_svg"

# --- 3. Load and Prepare All Data ---
# (This section is unchanged)
popgen_files <- list.files(path = popgen_dir, pattern = "_by_gene.csv", full.names = TRUE, recursive = TRUE)
gene_level_data <- popgen_files %>% map_dfr(read_csv)
ecology_df <- read_csv(metadata_file)
genus_map <- ecology_df %>% distinct(MAG_id, Genus)
mag_level_popgen <- gene_level_data %>%
  group_by(MAG_id, Colony, Month) %>% 
  summarise(mean_pi = mean(pi, na.rm = TRUE), .groups = 'drop')
full_data <- mag_level_popgen %>%
  left_join(genus_map, by = "MAG_id") %>%
  mutate(Month = factor(Month, levels = month_order)) %>%
  filter(!is.na(Month))

# --- 4. Define the Genus-Level Plotting Function ---
create_genus_plot <- function(Colony, Genus) {
  
  # (Data prep inside the function is unchanged)
  mag_ids_in_group <- full_data %>%
    filter(Colony == !!Colony, Genus == !!Genus) %>%
    distinct(MAG_id) %>% pull(MAG_id)
  if (length(mag_ids_in_group) == 0) return()
  replicate_totals <- ecology_df %>%
    filter(Caste == "Worker", Colony == !!Colony) %>%
    group_by(sample_name) %>%
    summarise(total_coverage = sum(mean, na.rm = TRUE), .groups = 'drop')
  abundance_summary <- ecology_df %>%
    filter(MAG_id %in% mag_ids_in_group, Colony == !!Colony, Caste == "Worker") %>%
    left_join(replicate_totals, by = "sample_name") %>%
    mutate(rep_rel_abund = mean / total_coverage) %>%
    group_by(MAG_id, Month) %>%
    summarise(mean_abund = mean(rep_rel_abund, na.rm = TRUE), sd_abund = sd(rep_rel_abund, na.rm = TRUE), .groups = 'drop') %>%
    mutate(Month = factor(Month, levels = month_order)) %>%
    complete(MAG_id, Month)
  gene_pi_data <- gene_level_data %>%
    filter(MAG_id %in% mag_ids_in_group, Colony == !!Colony) %>%
    mutate(Month = factor(Month, levels = month_order)) %>%
    rename(Value = pi) %>%
    complete(MAG_id, Month, unique_gene_id)
  avg_pi_data <- full_data %>%
    filter(MAG_id %in% mag_ids_in_group, Colony == !!Colony) %>%
    rename(Value = mean_pi) %>%
    complete(MAG_id, Month)
  if (nrow(abundance_summary) == 0 || nrow(avg_pi_data) == 0) {
    cat(paste("Skipping", Genus, "in Colony", Colony, "- no plot data after filtering.\n"))
    return()
  }
  
  # (Plot building is unchanged)
  p_abund <- ggplot(abundance_summary, aes(x = Month, y = mean_abund, group = 1)) +
    geom_ribbon(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), fill = gruvbox_palette["orange"], alpha = 0.5) +
    geom_line(color = gruvbox_palette["dark_bg"], linewidth = 1) +
    geom_point(color = gruvbox_palette["dark_bg"], size = 2) +
    facet_wrap(~ MAG_id, nrow = 1) + labs(y = "Relative Abundance") + scale_x_discrete(drop = FALSE) +
    theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_pi <- ggplot() +
    geom_line(data = gene_pi_data, aes(x = Month, y = Value, group = unique_gene_id), color = gruvbox_palette["blue"], alpha = 0.15) +
    geom_line(data = avg_pi_data, aes(x = Month, y = Value, group = 1), color = gruvbox_palette["dark_bg"], linewidth = 1.2) +
    facet_wrap(~ MAG_id, nrow = 1) + labs(x = "Month", y = "Nucleotide Diversity (π)") + scale_x_discrete(drop = FALSE) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks.x = element_blank())
  final_plot <- p_abund / p_pi + plot_annotation(title = paste("Temporal Dynamics of Genus", Genus, "in Colony", Colony))
  
  # --- MODIFIED: Saving process ---
  plot_width <- 4 + 2 * length(mag_ids_in_group)
  
  # 1. Save PNG (unchanged)
  output_dir_png <- file.path(base_output_dir_png, Colony, Genus)
  dir.create(output_dir_png, showWarnings = FALSE, recursive = TRUE)
  file_name_png <- paste0(Genus, "_dynamics.png")
  ggsave(file.path(output_dir_png, file_name_png), plot = final_plot, width = plot_width, height = 6, dpi = 300)
  
  # 2. Save SVG using base R graphics device
  output_dir_svg <- file.path(base_output_dir_svg, Colony, Genus)
  dir.create(output_dir_svg, showWarnings = FALSE, recursive = TRUE)
  file_name_svg <- paste0(Genus, "_dynamics.svg")
  
  # Open the SVG device
  svg(filename = file.path(output_dir_svg, file_name_svg), width = plot_width, height = 6)
  # You MUST explicitly print the ggplot object
  print(final_plot)
  # Close the device to finalize the file
  dev.off()
  
  cat(paste("Saved PNG and SVG for Genus", Genus, "in Colony", Colony, "\n"))
}

# --- 5. Run the Iteration ---
# (This section is unchanged)
plot_list <- full_data %>%
  filter(!is.na(Genus), Genus != "unclassified") %>%
  distinct(Colony, Genus)
purrr::pwalk(plot_list, .f = create_genus_plot)
cat("\n--- All plots generated successfully! ---\n")











########## correlations with alpha div metrics

abundance_wide <- full_data %>%
  select(Colony, Month, MAG_id, rel_abund_mag) %>%
  tidyr::pivot_wider(
    names_from = MAG_id,
    values_from = rel_abund_mag,
    values_fill = 0
  )

# b. Calculate Shannon diversity for each row (sample)
shannon_values <- diversity(abundance_wide %>% select(-Colony, -Month), "shannon")

# c. Create a dataframe with the diversity scores
shannon_df <- abundance_wide %>%
  select(Colony, Month) %>%
  mutate(shannon_diversity = shannon_values)

# d. Join the Shannon diversity back to our main 'full_data' table
full_data <- left_join(full_data, shannon_df, by = c("Colony", "Month"))


# a. Aggregate data to the colony level
colony_avg_shannon <- full_data %>%
  group_by(Colony) %>%
  summarise(avg_shannon = mean(shannon_diversity, na.rm = TRUE))

mag_avg_pi <- full_data %>%
  group_by(Colony, MAG_id) %>%
  summarise(avg_pi = mean(mean_pi, na.rm = TRUE), .groups = 'drop')

# b. Join the aggregated data
q1_data <- left_join(mag_avg_pi, colony_avg_shannon, by = "Colony")

# c. Fit the model
lmm_colony_level <- lmer(avg_pi ~ avg_shannon + (1 | MAG_id), data = q1_data)

cat("--- Model 1: Colony-Level Diversity vs. MAG Diversity ---\n")
summary(lmm_colony_level)

# Fit the GLMM on the sample-level data
glmm_temporal <- glmer(
  mean_pi ~ shannon_diversity + Season + mean_coverage + (1 | Colony) + (1 | MAG_id),
  data = full_data,
  family = Gamma(link = "log")
)

cat("\n--- Model 2: Temporal Fluctuation in Diversity vs. MAG Diversity ---\n")
summary(glmm_temporal)

# --- 4. Calculate Predicted Effects from the Model ---
# 'ggpredict' calculates the predicted 'mean_pi' across a range of 'shannon_diversity'
predicted_effects <- ggpredict(glmm_temporal, terms = "shannon_diversity")

# --- 5. Create the Plot ---
ggplot(predicted_effects, aes(x = x, y = predicted)) +
  # Add the raw data points in the background for context
  geom_point(
    data = full_data, 
    aes(x = shannon_diversity, y = mean_pi), 
    alpha = 1, 
    color = "gray"
  ) +
  # Add the model's confidence interval as a shaded ribbon
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#458588") +
  # Add the model's predicted regression line
  geom_line(color = "#282828", linewidth = 1.2) +
  labs(
    title = "Effect of Community Diversity on MAG Nucleotide Diversity",
    subtitle = "Prediction from GLMM with 95% Confidence Interval",
    x = "Shannon Diversity of Bee Gut Community",
    y = "Predicted Nucleotide Diversity (π) of MAG"
  ) +
  theme_bw()+
  scale_y_continuous(trans=("log10"))


# --- 1. Get the ID of the single most influential point ---
# We'll create a unique ID for each row to filter it out
data_with_id <- full_data %>%
  mutate(row_id = row_number())

most_influential_id <- data_with_id %>%
  mutate(cooks_distance = cooks.distance(glmm_temporal)) %>%
  arrange(desc(cooks_distance)) %>%
  slice(1) %>% # Get the top row
  pull(row_id) # Get its unique ID

# --- 2. Create a new dataset excluding that one point ---
data_filtered <- data_with_id %>%
  filter(row_id != most_influential_id)

# --- 3. Re-fit the GLMM on the filtered data ---
glmm_filtered <- glmer(
  mean_pi ~ shannon_diversity + (1 | Colony) + (1 | MAG_id) +,
  data = data_filtered,
  family = Gamma(link = "log")
)

# --- 4. Compare the Results ---
cat("\n--- Original Model Summary ---\n")
summary(glmm_temporal)

cat("\n--- Model Summary with Most Influential Point Removed ---\n")
summary(glmm_filtered)



####### RICHNESS ######
# --- 3. Calculate Richness ---
# For each sample (Colony-Month), count the number of unique MAGs present
richness_df <- full_data %>%
  group_by(Colony, Month) %>%
  summarise(richness = n_distinct(MAG_id), .groups = 'drop')

# Join the richness data back to the main dataframe
full_data <- left_join(full_data, richness_df, by = c("Colony", "Month"))

# --- 3. Fit the GLMM with Richness and Abundance ---
glmm_richness_abundance <- glmer(
  mean_pi ~ richness + mean_coverage + Season + (1 | Colony) + (1 | MAG_id),
  data = full_data,
  family = Gamma(link = "log")
)

# --- 4. View the Results ---
cat("--- GLMM Results: Controlling for Abundance ---\n")
summary(glmm_richness_abundance)





# --- 2. Calculate R-squared ---
r2_glmm <- r.squaredGLMM(glmm_temporal)

cat("--- GLMM Coefficient of Determination ---\n")
print(r2_glmm)


simulationOutput <- simulateResiduals(fittedModel = glmm_temporal, plot = TRUE)

testDispersion(simulationOutput)


glmm_null <- glmer(
  mean_pi ~ 1 + (1 | Colony) + (1 | MAG_id),
  data = full_data,
  family = Gamma(link = "log")
)

# --- 2. Compare the Models ---
# A Likelihood Ratio Test
model_comparison <- anova(glmm_temporal, glmm_null)

cat("\n--- Full vs. Null Model Comparison ---\n")
print(model_comparison)


genome_size_map <- ecology_df %>%
  distinct(MAG_id, genome_size)

# b. Join genome_size to the main analysis dataframe
full_data <- full_data %>%
  left_join(genome_size_map, by = "MAG_id")

# c. Scale all continuous predictor variables for model stability
full_data_scaled <- full_data %>%
  mutate(
    shannon_scaled = scale(shannon_diversity),
    coverage_scaled = scale(mean_coverage),
    genome_size_scaled = scale(genome_size) # Scale the new variable
  )

# --- 3. Fit the GLMM with a Different Optimizer ---
glmm_full_stable <- glmer(
  mean_pi ~ shannon_scaled + coverage_scaled + Season + genome_size_scaled + 
    (1 | Colony) + (1 | MAG_id),
  data = full_data_scaled,
  family = Gamma(link = "log"),
  # --- ADD THIS LINE to switch the optimizer ---
  control = glmerControl(optimizer = "bobyqa")
)

# --- 4. View the Results ---
# Check the summary for any warnings
summary(glmm_full_stable)



# --- 2. R-squared: How much variance does the model explain? ---
cat("--- 1. R-squared Values ---\n")
print(r.squaredGLMM(glmm_full_stable))

# --- 3. DHARMa Diagnostics: Check model assumptions ---
cat("\n--- 2. DHARMa Residual Plots ---\n")
# This simulates residuals from the final model to check its fit
simulationOutput <- simulateResiduals(fittedModel = glmm_full_stable, plot = TRUE)

# --- 4. Dispersion Test: Is the variance modeled correctly? ---
cat("\n--- 3. DHARMa Dispersion Test ---\n")
testDispersion(simulationOutput)

# --- 5. Null Model Comparison: Are the fixed effects significant overall? ---
glmm_null <- glmer(
  mean_pi ~ 1 + (1 | Colony) + (1 | MAG_id),
  data = full_data_scaled,
  family = Gamma(link = "log"),
  control = glmerControl(optimizer = "bobyqa") # Use same optimizer for fair comparison
)
cat("\n--- 4. Full vs. Null Model Comparison ---\n")
print(anova(glmm_full_stable, glmm_null))



# Create the final analysis dataframe
full_data <- inner_join(mag_level_popgen, abundance_data, by = c("Colony", "Month", "MAG_id")) %>%
  left_join(ecology_df %>% distinct(Month, Season), by = "Month") %>%
  left_join(ecology_df %>% distinct(MAG_id, genome_size), by = "MAG_id") %>%
  left_join(shannon_df, by = c("Colony", "Month")) %>%
  filter(!is.na(Season.y) & !is.na(genome_size) & !is.na(shannon_diversity))

# --- 4. Prepare Data for Beta GLMM ---
n_obs <- nrow(full_data)
full_data_scaled <- full_data %>%
  # --- NEW: Transform pi for Beta distribution (cannot handle exact 0s or 1s) ---
  mutate(pi_beta = (mean_pi * (n_obs - 1) + 0.5) / n_obs) %>%
  # Scale continuous predictors
  mutate(
    shannon_scaled = scale(shannon_diversity),
    coverage_scaled = scale(mean_coverage),
    genome_size_scaled = scale(genome_size)
  )
# --- 5. Fit the Full Beta GLMM ---
glmm_beta <- glmmTMB(
  pi_beta ~ shannon_scaled + coverage_scaled + Season.y + genome_size_scaled + 
    (1 | Colony) + (1 | MAG_id),
  data = full_data_scaled,
  family = beta_family(link = "logit")
)
cat("--- FINAL MODEL: Beta GLMM Summary ---\n")
summary(glmm_beta)


#--- 6. Run Full Diagnostic Suite on the Beta Model ---

# a. R-squared
cat("\n\n--- 1. R-squared Values ---\n")
print(MuMIn::r.squaredGLMM(glmm_beta))

# b. DHARMa Diagnostics
cat("\n\n--- 2. DHARMa Residual Plots and Tests ---\n")
simulationOutput_beta <- simulateResiduals(fittedModel = glmm_beta, plot = TRUE)
testDispersion(simulationOutput_beta)

# c. Null Model Comparison
cat("\n\n--- 3. Full vs. Null Model Comparison ---\n")
glmm_beta_null <- glmmTMB(
  pi_beta ~ 1 + (1 | Colony) + (1 | MAG_id),
  data = full_data_scaled,
  family = beta_family(link = "logit")
)
print(anova(glmm_beta, glmm_beta_null))





# --- 4. Calculate Predicted Effects ---
# Use ggpredict to get the marginal effect of shannon_scaled
predicted_effects <- ggpredict(glmm_beta, terms = "shannon_scaled [all]")

# --- 5. Create Custom Labels for the X-Axis ---
# This will transform the scaled x-axis back to the original Shannon values for interpretation
shannon_mean <- mean(full_data$shannon_diversity, na.rm = TRUE)
shannon_sd <- sd(full_data$shannon_diversity, na.rm = TRUE)

# Calculate what those values are on the original scale
original_labels <- (scaled_breaks * shannon_sd) + shannon_mean
# --- 2. Define Gruvbox Palette ---
gruvbox_palette <- c(
  dark_bg = "#282828",
  blue    = "#83a598",
  grey    = "#a89984"
)

# --- 3. Ensure Final Model and Data are in Environment ---
# (Assuming 'glmm_beta' and 'full_data_scaled' from the previous step are ready)

# --- 4. Calculate Predicted Effects ---
predicted_effects <- ggpredict(glmm_beta, terms = "shannon_scaled [all]")

# --- 5. Build the Plot ---
ggplot(predicted_effects, aes(x = x, y = predicted)) +
  
  geom_point(
    data = full_data_scaled, 
    aes(x = shannon_scaled, y = pi_beta), 
    color = gruvbox_palette["grey"], 
    alpha = 0.3, 
    shape = 16
  ) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high), 
    fill = gruvbox_palette["blue"], 
    alpha = 0.3
  ) +
  
  geom_line(color = gruvbox_palette["dark_bg"], linewidth = 1.2) +
  
  labs(
    title = "Higher Community Diversity is Associated with Higher Genetic Diversity",
    subtitle = "Prediction from Beta GLMM with 95% Confidence Interval",
    x = "Shannon Diversity of Bee Gut Community",
    y = expression(paste("Predicted Nucleotide Diversity (", pi, ")"))
  ) +
  
  # --- MODIFIED: Set specific breaks on the x-axis ---
  scale_x_continuous(
    breaks = {
      # Define the breaks you want in the ORIGINAL data scale
      original_breaks <- c(1, 1.5, 2, 2.5, 3)
      # Calculate the mean and sd of the original data
      shannon_mean <- mean(full_data$shannon_diversity, na.rm = TRUE)
      shannon_sd <- sd(full_data$shannon_diversity, na.rm = TRUE)
      # Convert your desired breaks to the SCALED axis
      (original_breaks - shannon_mean) / shannon_sd
    },
    # Use the original values as the labels
    labels = c("1.0", "1.5", "2.0", "2.5", "3.0")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )+
  scale_y_continuous(trans="log10")
# --- 3. Extract All Model Components ---
# (This section is unchanged)
model_tidied <- broom.mixed::tidy(glmm_beta)
fixed_effects <- model_tidied %>%
  filter(effect == "fixed") %>%
  select(term, estimate, std.error, statistic, p.value)
random_effects <- model_tidied %>%
  filter(effect == "ran_pars" & term == "sd__(Intercept)") %>%
  select(group, estimate)
r2_values <- MuMIn::r.squaredGLMM(glmm_beta)
r2m <- r2_values[1, "R2m"]
r2c <- r2_values[1, "R2c"]

# --- 4. Build the Summary Dataframe ---
# (This section is unchanged)
summary_table <- bind_rows(
  fixed_effects %>%
    mutate(group = "**Fixed Effects**") %>%
    rename(Parameter = term),
  tibble(
    group = "**Model Summary**",
    Parameter = c(
      paste("SD (", random_effects$group, ")", sep = ""),
      "Marginal R² (Fixed Effects)",
      "Conditional R² (Full Model)"
    ),
    estimate = c(random_effects$estimate, r2m, r2c)
  )
)

# --- 5. Create and Print the 'gt' Table ---
summary_gt <- summary_table %>%
  gt(groupname_col = "group") %>%
  tab_header(
    title = "Beta GLMM Results",
    subtitle = "Predicting Nucleotide Diversity (π)"
  ) %>%
  cols_label(
    Parameter = "Parameter",
    estimate = "Estimate",
    std.error = "Std. Error",
    statistic = "z value",
    p.value = "p-value"
  ) %>%
  fmt_number(
    columns = c(estimate, std.error, statistic),
    decimals = 3
  ) %>%
  # This function should now exist after you update
  fmt_pvalue(columns = p.value, decimals = 3) %>%
  
  # --- MODIFIED: Use the new function name 'sub_missing()' ---
  sub_missing(
    columns = c(std.error, statistic, p.value),
    missing_text = ""
  )

# Print the final table
summary_gt



# --- 2. Fit Simpler Models for Each Predictor ---
# Model with only Shannon diversity
model_shannon <- glmmTMB(pi_beta ~ shannon_scaled + (1|Colony) + (1|MAG_id), data = full_data_scaled, family = beta_family())

# Model with only Season
model_season <- glmmTMB(pi_beta ~ Season + (1|Colony) + (1|MAG_id), data = full_data_scaled, family = beta_family())

# --- 3. Get Predictions for Each Model ---
preds_shannon <- ggpredict(model_shannon, terms = "shannon_scaled [all]")
preds_season <- ggpredict(model_season, terms = "Season")

# --- 4. Create a Plot for Each Effect ---

# Plot A: Shannon Effect
plot_a <- ggplot(preds_shannon, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#83a598", alpha = 0.5) +
  geom_point()+
  geom_line(color = "#282828") +
  labs(subtitle = "Effect of Community Diversity", x = "Shannon Diversity", y = expression(pi)) +
  scale_x_continuous(labels = ~round((.x * sd(full_data$shannon_diversity)) + mean(full_data$shannon_diversity), 1)) +
  theme_bw()

# Plot B: Season Effect
plot_b <- ggplot(preds_season, aes(x = x, y = predicted)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "#fb4934") +
  geom_point(size = 4, color = "#282828") +
  labs(subtitle = "Effect of Season", x = "Season", y = NULL) +
  theme_bw()

# --- 5. Combine Plots ---
(plot_a | plot_b) + 
  plot_annotation(title = "Comparing the Relative Strength of Predictors on Nucleotide Diversity")



# --- 2. Ensure Models and Data are in Environment ---
# (Assuming 'glmm_beta', 'glmm_beta_null', and 'full_data_scaled' are ready)

# --- 3. Calculate Predictions for the FULL Model ---
predicted_full <- ggpredict(glmm_beta, terms = "shannon_scaled [all]")

# --- 4. Get the Grand Mean from the NULL Model (MODIFIED) ---
# --- MODIFIED: Access the conditional model's intercept specifically ---
grand_mean <- plogis(fixef(glmm_beta_null)$cond[1])


# --- 5. Build the Plot ---
ggplot(predicted_full, aes(x = x, y = predicted)) +
  
  # a. Raw data points (as before)
  geom_point(
    data = full_data_scaled, 
    aes(x = shannon_scaled, y = pi_beta), 
    color = "#a89984", alpha = 0.3, shape = 16
  ) +
  
  # b. Full model's confidence interval and line (as before)
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#83a598", alpha = 0.3) +
  geom_line(color = "#282828", linewidth = 1.2) +
  
  # c. --- NEW: Add the null model's prediction as a horizontal line ---
  geom_hline(
    yintercept = grand_mean, 
    color = "#fb4934",       # Gruvbox red
    linewidth = 1.2, 
    linetype = "dashed"
  ) +
  
  # d. Labels and theme (as before)
  labs(
    title = "Effect of Community Diversity vs. Grand Mean",
    subtitle = "Solid Line: Full Model Prediction. Dashed Line: Null Model (Grand Mean).",
    x = "Shannon Diversity of Bee Gut Community",
    y = expression(paste("Predicted Nucleotide Diversity (", pi, ")"))
  ) +
  scale_x_continuous(labels = ~round((.x * sd(full_data$shannon_diversity)) + mean(full_data$shannon_diversity), 1)) +
  theme_classic(base_size = 14)

ggplot(predicted_effects, aes(x = x, y = predicted)) +
  
  geom_point(
    data = full_data_scaled, 
    aes(x = shannon_scaled, y = pi_beta), 
    color = gruvbox_palette["grey"], 
    alpha = 0.3, 
    shape = 16
  ) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high), 
    fill = gruvbox_palette["blue"], 
    alpha = 0.3
  ) +
  
  geom_line(color = gruvbox_palette["dark_bg"], linewidth = 1.2) +
  
  labs(
    title = "Higher Community Diversity is Associated with Higher Genetic Diversity",
    subtitle = "Prediction from Beta GLMM with 95% Confidence Interval",
    x = "Shannon Diversity of Bee Gut Community",
    y = expression(paste("Predicted Nucleotide Diversity (", pi, ")"))
  ) +
  
  scale_x_continuous(
    breaks = {
      original_breaks <- c(1, 1.5, 2, 2.5, 3)
      shannon_mean <- mean(full_data$shannon_diversity, na.rm = TRUE)
      shannon_sd <- sd(full_data$shannon_diversity, na.rm = TRUE)
      (original_breaks - shannon_mean) / shannon_sd
    },
    labels = c("1.0", "1.5", "2.0", "2.5", "3.0")
  ) +
  
  # --- MODIFIED: Use coord_cartesian() to set y-axis limits ---
  coord_cartesian(ylim = c(0.001, 0.01)) +
  # --- REMOVED: scale_y_continuous(trans="log10") ---
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, margin = margin(b = 15)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )



# --- 2. Load and Prepare Data ---
# This section reproduces the 'full_data' dataframe with all necessary variables
popgen_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/popgen_summaries_final_filtered"
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"

popgen_files <- list.files(path = popgen_dir, pattern = "_by_gene.csv", full.names = TRUE, recursive = TRUE)
gene_level_data <- popgen_files %>% map_dfr(read_csv)
ecology_df <- read_csv(metadata_file)

mag_level_popgen <- gene_level_data %>%
  group_by(MAG_id, Colony, Month) %>% 
  summarise(mean_pi = mean(pi, na.rm = TRUE), .groups = 'drop')

abundance_data <- ecology_df %>%
  filter(Caste == "Worker") %>%
  group_by(Colony, Month, MAG_id) %>%
  summarise(mean_coverage = mean(mean, na.rm = TRUE), .groups = 'drop')

shannon_df <- inner_join(mag_level_popgen, abundance_data) %>%
  group_by(Colony, Month) %>%
  mutate(total_abund = sum(mean_coverage)) %>%
  ungroup() %>%
  mutate(rel_abund = mean_coverage / total_abund) %>%
  group_by(Colony, Month) %>%
  summarise(shannon_diversity = vegan::diversity(rel_abund, "shannon"), .groups = 'drop')

full_data <- inner_join(mag_level_popgen, abundance_data, by = c("Colony", "Month", "MAG_id")) %>%
  left_join(ecology_df %>% distinct(MAG_id, Genus), by = "MAG_id") %>%
  left_join(shannon_df, by = c("Colony", "Month")) %>%
  filter(!is.na(Genus), !is.na(shannon_diversity), !is.na(mean_pi))

# --- 3. Define Output Directory ---
base_output_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/figures/correlation_plots"


# --- 4. Analysis 1: Plotting by Genus (MODIFIED) ---

# a. Define the plotting function
plot_genus_correlation <- function(target_genus, data) {
  
  genus_data <- data %>% filter(Genus == target_genus)
  
  if (nrow(genus_data) < 3) {
    cat(paste("Skipping Genus", target_genus, "- not enough data points.\n"))
    return()
  }
  
  corr_test <- cor.test(genus_data$shannon_diversity, genus_data$mean_pi)
  corr_label <- paste0("r = ", round(corr_test$estimate, 2), ", p = ", format.pval(corr_test$p.value, digits = 2))
  
  line_color <- "gray50"
  if (corr_test$p.value < 0.05) {
    if (corr_test$estimate > 0) { line_color <- "#fb4934" } else { line_color <- "#83a598" }
  }
  
  p <- ggplot(genus_data, aes(x = shannon_diversity, y = mean_pi)) +
    # --- MODIFIED: Set point color to black ---
    geom_point(color = "black", alpha = 0.6) +
    # --- MODIFIED: Make line thicker ---
    geom_smooth(method = "lm", se = FALSE, color = line_color, linewidth = 1.5) +
    annotate("text", x = -Inf, y = Inf, label = corr_label, hjust = -0.1, vjust = 1.5, size = 4) +
    labs(title = paste("Genus:", target_genus), subtitle = "Nucleotide Diversity vs. Community (Shannon) Diversity", x = "Shannon Diversity", y = "Nucleotide Diversity (π)") +
    theme_bw()
  
  # --- MODIFIED: Save as SVG using base R device ---
  output_dir <- file.path(base_output_dir, "by_genus_svg")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  svg_path <- file.path(output_dir, paste0(target_genus, "_correlation.svg"))
  
  svg(filename = svg_path, width = 6, height = 5)
  print(p)
  dev.off()
  
  cat(paste("Saved SVG plot for Genus:", target_genus, "\n"))
}

# b. Run the loop for Genus plots
unique_genera <- unique(full_data$Genus)
purrr::walk(unique_genera, ~plot_genus_correlation(.x, full_data))

  # --- 5. Analysis 2: Plotting by MAG_id ---

# a. Define the plotting function
plot_mag_correlation <- function(target_mag, data) {
  
  mag_data <- data %>% filter(MAG_id == target_mag)
  
  if (nrow(mag_data) < 3) {
    cat(paste("Skipping MAG", target_mag, "- not enough data points.\n"))
    return()
  }
  
  corr_test <- cor.test(mag_data$shannon_diversity, mag_data$mean_pi)
  corr_label <- paste0("r = ", round(corr_test$estimate, 2), ", p = ", format.pval(corr_test$p.value, digits = 2))
  
  # --- NEW: Choose color based on significance and direction ---
  line_color <- "#665c54" # Default for non-significant
  if (corr_test$p.value < 0.05) {
    if (corr_test$estimate > 0) {
      line_color <- "#fb4934" # Red for positive
    } else {
      line_color <- "#458588" # Blue for negative
    }
  }
  
  p <- ggplot(mag_data, aes(x = shannon_diversity, y = mean_pi)) +
    geom_point(aes(color = Colony), alpha = 0.8) +
    # --- MODIFIED: Use the new line_color variable ---
    geom_smooth(method = "lm", se = FALSE, color = line_color) +
    annotate("text", x = -Inf, y = Inf, label = corr_label, hjust = -0.1, vjust = 1.5, size = 4) +
    labs(title = paste("MAG:", target_mag), subtitle = "Nucleotide Diversity vs. Community (Shannon) Diversity", x = "Shannon Diversity", y = "Nucleotide Diversity (π)") +
    theme_bw()
  
  output_dir <- file.path(base_output_dir, "by_mag_id")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = file.path(output_dir, paste0(target_mag, "_correlation.png")), plot = p, width = 6, height = 5, dpi = 150)
  cat(paste("Saved plot for MAG:", target_mag, "\n"))
}

# b. Run the loop
unique_mags <- unique(full_data$MAG_id)
purrr::walk(unique_mags, ~plot_mag_correlation(.x, full_data))

cat("\n--- All correlation plots generated successfully! ---\n")

