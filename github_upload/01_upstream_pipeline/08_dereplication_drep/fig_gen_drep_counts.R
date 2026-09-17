library(ggplot2)
library(dplyr)
library(tidyr)
df = read.csv('RENAMED_CONTIG_TAXONOMY_20250528.csv')
df = df[df$type=="mmag",]
mag_counts_by_genus <- df %>%
  group_by(Genus) %>%
  summarise(
    unique_MAG_count = n_distinct(MAG_id),
    .groups = 'drop'
  )
print(mag_counts_by_genus)
r_mag_counts_by_genus <- df %>%
  group_by(Genus) %>%
  summarise(
    unique_rMAG_count = n_distinct(r_MAG_id),
    .groups = 'drop'
  )

cat("\n--- Unique r_MAG_id Counts per Genus ---\n")
print(r_mag_counts_by_genus)

# --- 3. Prepare the Data for Plotting ---
# Join the two datasets and calculate the "compressed" portion
plot_data_prep <- full_join(mag_counts_by_genus, r_mag_counts_by_genus, by = "Genus") %>%
  # For any genus with no compression, the rMAG count equals the MAG count
  mutate(unique_rMAG_count = ifelse(is.na(unique_rMAG_count), unique_MAG_count, unique_rMAG_count)) %>%
  # Calculate the number of MAGs that were compressed
  mutate(compressed_count = unique_MAG_count - unique_rMAG_count) %>%
  # Keep only the necessary columns for stacking
  select(Genus, unique_rMAG_count, compressed_count)

# Pivot the data into a long format suitable for ggplot
plot_data_long <- plot_data_prep %>%
  pivot_longer(
    cols = c(unique_rMAG_count, compressed_count),
    names_to = "Category",
    values_to = "Count"
  )

# --- 4. Generate the Stacked Bar Plot ---
compression_plot <- ggplot(
  # Use the original prep data to reorder the Genus factor by total MAG count
  data = plot_data_long,
  aes(x = reorder(Genus, -Count, sum), y = Count, fill = Category)
) +
  geom_bar(stat = "identity", position = "stack") +
  # Set custom colors and labels
  scale_fill_manual(
    name = "MAG Type",
    values = c("unique_rMAG_count" = "#2c3e50", "compressed_count" = "#bdc3c7"), # Dark blue-grey and light grey
    labels = c("unique_rMAG_count" = "Representative MAGs (r_MAGs)", "compressed_count" = "Additional (Compressed) MAGs")
  ) +
  labs(
    title = "MAG Compression into Representative MAGs by Genus",
    x = "Genus",
    y = "Number of Unique MAGs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


# --- 1. Create the new plot ordered by r_MAG_id count ---
# First, determine the descending order of Genera based on r_MAG_id count
genus_order_by_rMAG <- plot_data_prep %>%
  arrange(desc(unique_rMAG_count)) %>%
  pull(Genus)

# Generate the plot, applying the new order to the Genus factor
plot_ordered_by_rMAG <- ggplot(
  plot_data_long,
  # Apply the custom order to the x-axis
  aes(x = factor(Genus, levels = genus_order_by_rMAG), y = Count, fill = Category)
) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    name = "MAG Type",
    values = c("unique_rMAG_count" = "#2c3e50", "compressed_count" = "#bdc3c7"),
    labels = c("unique_rMAG_count" = "Representative MAGs (r_MAGs)", "compressed_count" = "Additional (Compressed) MAGs")
  ) +
  labs(
    title = "MAG Compression Ordered by r_MAG Count",
    x = "Genus",
    y = "Number of Unique MAGs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(plot_ordered_by_rMAG)


Yes, we can apply those modifications. The process involves a few data preparation steps to combine the specified genera and filter by the MAG count before feeding the data into the plotting script you've settled on.

The coord_flip() you've added is a great choice for making the genus labels more readable.

The R Script

This script takes the count dataframes, performs the requested cleaning and filtering, and then generates the final plot.
R

# --- 1. Load Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 2. Use the count dataframes from the previous step ---
# (Using a more complex sample dataset to demonstrate the changes)
mag_counts_by_genus <- data.frame(
  Genus = c("Lactobacillus", "Bifidobacterium", "Gilliamella", "unknown", "Snodgrassella", "no match", "Frischella"),
  unique_MAG_count = c(25, 18, 12, 8, 7, 5, 2)
)
r_mag_counts_by_genus <- data.frame(
  Genus = c("Lactobacillus", "Bifidobacterium", "Gilliamella", "unknown", "Snodgrassella", "no match", "Frischella"),
  unique_rMAG_count = c(8, 15, 12, 4, 6, 3, 2)
)

# --- 3. Data Preparation ---

# A. Combine 'unknown' and 'no match' into 'Unidentified' and re-aggregate counts
mag_counts_combined <- mag_counts_by_genus %>%
  mutate(Genus = if_else(Genus %in% c("unknown", "no match"), "Unidentified", Genus)) %>%
  group_by(Genus) %>%
  summarise(unique_MAG_count = sum(unique_MAG_count), .groups = 'drop')

r_mag_counts_combined <- r_mag_counts_by_genus %>%
  mutate(Genus = if_else(Genus %in% c("unknown", "no match"), "Unidentified", Genus)) %>%
  group_by(Genus) %>%
  summarise(unique_rMAG_count = sum(unique_rMAG_count), .groups = 'drop')

# B. Join the datasets and filter for genera with > 10 total MAGs
plot_data_prep <- full_join(mag_counts_combined, r_mag_counts_combined, by = "Genus") %>%
  filter(unique_rMAG_count > 2)

# C. Calculate the 'compressed' portion and pivot to a long format
plot_data_long <- plot_data_prep %>%
  mutate(compressed_count = unique_MAG_count - unique_rMAG_count) %>%
  select(Genus, unique_rMAG_count, compressed_count) %>%
  pivot_longer(
    cols = c(unique_rMAG_count, compressed_count),
    names_to = "Category",
    values_to = "Count"
  )

# --- 4. Generate the Final Plot ---
# This is the ggplot script you provided, now using our filtered and cleaned data
plot_ordered_by_MAG <- ggplot(
  plot_data_long,
  # reorder(Genus, Count, sum) sorts Genus by the sum of Count in ascending order
  aes(x = reorder(Genus, Count, sum), y = Count, fill = Category)
) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    name = "MAG Type",
    values = c("unique_rMAG_count" = "#2c3e50", "compressed_count" = "#bdc3c7"),
    labels = c("unique_rMAG_count" = "Representative MAGs (r_MAGs)", "compressed_count" = "Assembled MAGs")
  ) +
  labs(
    title = "MAG Dereplication for Genera with >1 Representative MAGs",
    x = "Genus",
    y = "Number of Unique MAGs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top"
  ) +
  coord_flip() # Flip coordinates to make bars horizontal

print(plot_ordered_by_MAG)

