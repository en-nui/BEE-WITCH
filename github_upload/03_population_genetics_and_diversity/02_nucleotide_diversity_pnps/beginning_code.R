# --- 1. Load Libraries ---
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
# --- 2. Load and Prepare Initial Data ---

# a. Your function to create a numeric Time variable
add_cumulative_time_from_month_nucdiv <- function(data) {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34, 
                        "September" = 28, "October" = 22, "November" = 42, 
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  data$Time <- month_to_cumulative_time_map[data$Month]
  return(data)
}

# b. Load and aggregate metadata, then add the Time column
metadata_file <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/ecology_analyses/ecology_metadata_df.csv"
ecology_metadata <- read.csv(metadata_file)

unique_metadata <- ecology_metadata %>%
  group_by(sample_name, MAG_id) %>%
  summarise(
    mean_coverage = mean(mean, na.rm = TRUE),
    Colony = first(Colony),
    Month = first(Month),
    .groups = 'drop'
  ) %>%
  add_cumulative_time_from_month_nucdiv() # Apply your time function

# c. Load and prepare popgen data
popgen_dir <- "/home/robinch/projects/BEE-WITCH/34_mMAG_popgen/ALLSPECIES_intrapop_diversity/cnv_filtered_snps_by_colony/polarized_file_by_file/popgen_summaries"
file_list <- list.files(path = popgen_dir, pattern = "colony_.*_filtered_snps_by_gene\\.csv", full.names = TRUE)
all_popgen_data <- purrr::map_dfr(file_list, read.csv) %>%
  rename(sample_name = sample_id)

# d. Join the datasets
merged_data <- inner_join(all_popgen_data, unique_metadata, by = c("sample_name", "MAG_id"))
cat("Initial MAG count before any filtering:", n_distinct(merged_data$MAG_id), "\n")

# --- 3. Identify Core MAGs (The New Filter) ---

# a. For each Colony, find the total number of months it was sampled
colony_month_totals <- merged_data %>%
  group_by(Colony) %>%
  summarise(total_months = n_distinct(Month), .groups = 'drop')

# b. For each MAG in each Colony, count the number of months it was observed
mag_persistence_counts <- merged_data %>%
  group_by(Colony, MAG_id) %>%
  summarise(mag_months = n_distinct(Month), .groups = 'drop')

# c. Create a "key" of core MAGs by finding where the counts match
core_mag_key <- inner_join(mag_persistence_counts, colony_month_totals, by = "Colony") %>%
  filter(mag_months == total_months) %>%
  select(Colony, MAG_id)

# --- 4. Filter The Main Dataset ---
# Use a semi_join to keep only the data for the core MAGs we identified
core_filtered_data <- semi_join(merged_data, core_mag_key, by = c("Colony", "MAG_id"))

cat("MAG count after 100% persistence filter:", n_distinct(core_filtered_data$MAG_id), "\n")

# --- 5. Run Correlation on Filtered Data ---

# a. Filter by coverage and aggregate to get median diversity
mag_diversity_over_time <- core_filtered_data %>%
  filter(mean_coverage >= 5) %>%
  group_by(MAG_id, sample_name, Colony, Time) %>%
  summarise(median_nucDiv = median(nucDiv_per_bp, na.rm = TRUE), .groups = 'drop')

# b. Perform the correlation
correlation_table <- mag_diversity_over_time %>%
  group_by(MAG_id, Colony) %>%
  filter(n() >= 3) %>%
  do(broom::tidy(cor.test(~ median_nucDiv + Time, data = ., method = "pearson"))) %>%
  ungroup() %>%
  select(MAG_id, Colony, estimate, p.value) %>%
  rename(pearson_r = estimate) %>%
  arrange(p.value)

# --- 6. View the Final Results ---
cat("\n--- Correlation of Median Nucleotide Diversity vs. Time (Core MAGs) ---\n")
print(as.data.frame(correlation_table))

ggplot(correlation_table,aes(x=pearson_r))+
  geom_histogram()


