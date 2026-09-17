library(dplyr)
library(tidyr)
library(stats)
library(ggplot2)
library(mgcv)
library(tweedie)
library(rstatix)
# --- Adapted Function to Add Cumulative Time ---
add_cumulative_time_from_month_nucdiv <- function(data) {
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
  # *** MODIFIED HERE to use data$Month instead of data$Month2 ***
  data$Time <- month_to_cumulative_time_map[data$Month]
  
  # Check if any months in the data were not found in the mapping
  if (any(is.na(data$Time))) {
    warning("Warning: Some months in the 'Month' column were not found in the mapping and have been assigned NA in the 'Time' column.")
  }
  
  return(data)
}


# --- Define File Paths ---
nucdiv_file <- "/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/genic/all_nucdiv_by_gene_RETAINED_filtered.csv"
pairwise_fst_file <- "/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/genic/multiple_colony_same_month/gene_level_pairwise_fst_results_with_annotations.csv"
vs_initial_fst_file <- "/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/genic/single_colony_multiple_month/gene_level_vs_initial_month_fst_single_colony_with_annotations.csv"
between_colony_fst_file <- "/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/genic/multiple_colony_multiple_month/gene_level_between_colony_fst_multiple_colony_with_annotations.csv"
temporal_fst_file <- "/home/robinch/projects/BEE-WITCH/22_anvio/fst_analyses/genic/multiple_colony_multiple_month/gene_level_temporal_fst_multiple_colony_with_annotations.csv"

# --- Load Data ---
nucdiv_df <- read.csv(nucdiv_file)
pairwise_fst_df <- read.csv(pairwise_fst_file)
vs_initial_fst_df <- read.csv(vs_initial_fst_file)
between_colony_fst_df <- read.csv(between_colony_fst_file)
temporal_fst_df <- read.csv(temporal_fst_file)
# Assuming nucdiv_df and vs_initial_fst_df are already loaded
nucdiv_df$nucdiv_by_gene_length <- nucdiv_df$nucDiv_by_gene/nucdiv_df$gene_length



# --- Identify unique MAG_id and Month combinations from vs_initial_fst_df ---

# Get MAG_id and Month1 combinations
mag_month1_combinations <- vs_initial_fst_df %>%
  select(MAG_id, Month1) %>%
  rename(Month = Month1)

# Get MAG_id and Month2 combinations
mag_month2_combinations <- vs_initial_fst_df %>%
  select(MAG_id, Month2) %>%
  rename(Month = Month2)

# Combine and get unique MAG_id and Month combinations
unique_mag_month_combinations_from_vs_initial <- bind_rows(mag_month1_combinations, mag_month2_combinations) %>%
  distinct() # Ensure uniqueness

# --- Filter nucdiv_df based on these combinations ---

# Filter nucdiv_df to keep rows where the MAG_id and Month combination is in the unique combinations
filtered_nucdiv_df <- nucdiv_df %>%
  semi_join(unique_mag_month_combinations_from_vs_initial, by = c("MAG_id", "Month"))

# Remove duplicate rows based on MAG_id, gene_pos, Month, and Rep
unique_filtered_nucdiv_df <- filtered_nucdiv_df %>%
  distinct(MAG_id, gene_pos, Month, Rep, .keep_all = TRUE)
unique_filtered_nucdiv_df_with_time <- add_cumulative_time_from_month_nucdiv(unique_filtered_nucdiv_df)

# --- Fit a Linear Model with Interaction Term ---
# Model: nucDiv_by_gene changes with Time, and the way it changes differs by lifestyle
lm_model <- lm(Wtheta_per_SNV ~ Time * lifestyle, data = unique_filtered_nucdiv_df_with_time)

# --- Summarize the Model Results ---
summary(lm_model)
plot(lm_model)

# --- Find the minimum number of unique Time points per lifestyle ---
unique_time_points_per_lifestyle <- unique_filtered_nucdiv_df_with_time %>%
  group_by(lifestyle) %>%
  summarize(n_unique_time = n_distinct(Time))

min_unique_time <- min(unique_time_points_per_lifestyle$n_unique_time, na.rm = TRUE) # Use na.rm in case of NA lifestyles

ggplot(nucdiv_df,aes(x=Wtheta_per_SNV))+
  geom_histogram()


ggplot(data = unique_filtered_nucdiv_df_with_time, aes(x = Time, y = nucdiv_by_gene_length, color = lifestyle)) +
  geom_point(alpha = 0.3, size = 1) + # Show individual data points with transparency
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) + # Add fitted linear regression lines for each lifestyle
  labs(
    title = "Linear Model Fit of Nucleotide Diversity by Gene Length Over Time by Lifestyle",
    x = "Time (Days Since May 1st)",
    y = "Nucleotide Diversity by Gene Length",
    color = "Lifestyle" # Legend title for color
  ) +
  theme_minimal()+
  scale_y_continuous(trans="log10")
# Optional: Adjust x-axis scale if needed
# scale_x_continuous(limits = c(0, max(unique_filtered_nucdiv_df_with_time$Time, na.rm = TRUE))) +
# Optional: Adjust y-axis scale if needed
# scale_y_continuous(limits = c(0, max(unique_filtered_nucdiv_df_with_time$nucDiv_by_gene, na.rm = TRUE)))



# --- Fit a Generalized Additive Model (GAM) with reduced 'k' ---
# We set k to be less than or equal to the minimum number of unique time points.
# A common choice is min_unique_time - 1, or a small value like 3 or 4
# that still allows for a curved relationship.
# Let's try min_unique_time - 1, but ensure it's at least 2 (for a line) or 3 (for a curve).
# If min_unique_time is very small (e.g., 1 or 2), a smooth might not be possible,
# and a simpler model (e.g., linear) might be necessary for that term.

# --- Ensure 'lifestyle' is a factor ---
# It's good practice to ensure grouping variables are factors
unique_filtered_nucdiv_df_with_time <- unique_filtered_nucdiv_df_with_time %>%
  mutate(lifestyle = as.factor(lifestyle))

# --- Check the dataframe structure and column presence just before gam call ---
cat("Checking data structure before gam call:\n")
print(str(unique_filtered_nucdiv_df_with_time))
cat("\nChecking for 'lifestyle' column presence:\n")
print("lifestyle" %in% colnames(unique_filtered_nucdiv_df_with_time))
cat("\nFirst few rows of relevant columns:\n")
print(head(unique_filtered_nucdiv_df_with_time %>% select(Time, lifestyle, nucDiv_by_gene)))
cat("\nUnique values in lifestyle column:\n")
print(unique(unique_filtered_nucdiv_df_with_time$lifestyle))


# --- Find the minimum number of unique Time points per lifestyle (re-calculating just in case) ---
unique_time_points_per_lifestyle <- unique_filtered_nucdiv_df_with_time %>%
  group_by(lifestyle) %>%
  summarize(n_unique_time = n_distinct(Time))

min_unique_time <- min(unique_time_points_per_lifestyle$n_unique_time, na.rm = TRUE)

cat("Minimum number of unique Time points within any lifestyle:", min_unique_time, "\n")

# Choose a value for k - ensuring k is at least 2 or 3 for a curve
chosen_k <- max(3, min_unique_time - 1)
cat("Using k =", chosen_k, "for smooth terms.\n")


# --- Fit the Generalized Additive Model (GAM) with reduced 'k' ---
# Specify the reduced k in the s() terms
cat("\nAttempting to fit GAM...\n")
gam_model <- gam(nucdiv_by_gene_length ~ s(Time, k = chosen_k) + s(Time, by = lifestyle, k = chosen_k),
                 data = unique_filtered_nucdiv_df_with_time)

# --- Summarize the Model Results ---
summary(gam_model)

# --- Interpret the Interaction Term ---
# As before, look at the p-values for the s(Time):lifestyle terms.

# --- Optional: Plot the GAM Smooths ---
 plot(gam_model, pages = 1, all.terms = TRUE) # plots all smooth terms

# --- Optional: Check Model Assumptions ---
# Check assumptions of the GAM
gam.check(gam_model) # Provides diagnostic plots and tests for checking assumptions


unique_filtered_nucdiv_df_with_time$nucdiv_by_gene_length = unique_filtered_nucdiv_df_with_time$nucdiv_by_gene_length + 0.00000000000000000001
glm_model <- glm(nucdiv_by_gene_length ~ Time * lifestyle, family = Gamma(link = "log"), data = unique_filtered_nucdiv_df_with_time)
summary(glm_model)
plot(glm_model, pages = 1, all.terms = TRUE)

# --- Create a dataframe for plotting the fitted lines ---
# Generate a sequence of Time values
time_range <- range(unique_filtered_nucdiv_df_with_time$Time, na.rm = TRUE)
new_time_data <- seq(time_range[1], time_range[2], length.out = 100) # Generate 100 points across the time range

# Get all unique lifestyle levels
lifestyle_levels <- levels(unique_filtered_nucdiv_df_with_time$lifestyle)

# Create a new dataframe with combinations of Time and lifestyle for prediction
predict_data <- expand.grid(Time = new_time_data, lifestyle = lifestyle_levels)

# Add predicted values from the GLM (on the response scale)
predict_data$predicted_nucdiv <- predict(glm_model, newdata = predict_data, type = "response")

# --- Plot original data and fitted GLM lines ---

ggplot(data = unique_filtered_nucdiv_df_with_time, aes(x = Time, y = nucdiv_by_gene_length, color = lifestyle)) +
  geom_point(alpha = 0.3, size = 1) + # Show original data points
  geom_line(data = predict_data, aes(x = Time, y = predicted_nucdiv, color = lifestyle), linewidth = 1.2) + # Add fitted lines
  labs(
    title = "GLM (Gamma) Fit of Nucleotide Diversity by Gene Length Over Time by Lifestyle",
    x = "Time (Days Since May 1st)",
    y = "Nucleotide Diversity by Gene Length",
    color = "Lifestyle" # Legend title for color
  ) +
  theme_minimal()+
  scale_y_continuous(trans="log10", limits = c(1e-4, 1e-2))
# Optional: Adjust x-axis scale if needed
# scale_x_continuous(limits = time_range) +
# Optional: Adjust y-axis scale if needed
# scale_y_continuous(limits = c(0, max(unique_filtered_nucdiv_df_with_time$nucdiv_by_gene_length, na.rm = TRUE)))