library(ggplot2)
library(tidyr)
library(dplyr)
library(fs)
library(purrr)
library(readr)


# 1. Define the directory path
data_dir <- "/home/robinch/projects/BEE-WITCH/16b_rMAG_mapping/analysis_output"
file_paths <- dir_ls(data_dir, glob = "*_Larvae*.filtered_mmag_contigs.csv")
combined_df <- map_dfr(file_paths, read_csv)
combined_df$genome_wide_median_cov <- combined_df$genome_wide_median_cov + 0.000001

# Define your actual median coverage cutoff value
# For this example, we'll just use > 0
median_cutoff <- 1

# Modify the dataframe to add a new 'presence_status' column
# and then create the plot.
larvae_mag_plot_v2 <- combined_df %>%
  mutate(presence_status = if_else(genome_wide_median_cov > median_cutoff, 
                                   "Present", 
                                   "Absent")) %>%
  ggplot(aes(x = contig_presence_ratio, 
             y = average_genome_wide_cov, 
             color = presence_status)) +
  geom_point(alpha = 0.8, size = 2) +
  # Use a log scale if average coverage also has a wide range (optional)
  # scale_y_log10() + 
  scale_color_manual(values = c("Present" = "dodgerblue", "Absent" = "grey70")) +
  labs(
    title = "MAG Average Coverage vs. Presence in Honey Bee Larvae",
    subtitle = "Color indicates presence/absence based on median coverage >= 1",
    x = "Contig Presence Ratio",
    y = "Genome-Wide Average Coverage",
    color = "Status (by Median Cov)" # This changes the legend title
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")+
  geom_vline(xintercept=0.5)

# Display the new plot
print(larvae_mag_plot_v2)
