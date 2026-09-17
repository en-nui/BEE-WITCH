# --- 1. Filter 'pointrange_data' for Core Genera ---
abundance_plot_data_core <- pointrange_data %>%
  filter(Color_Genus %in% target_genera)

# --- 2. Create the Filtered Abundance Plot ---
range_abundance_plot_core <- ggplot(
  abundance_plot_data_core, # Use the filtered data
  aes(
    x = reorder(mag_colony_id, -median_rel_abund),
    color = Color_Genus 
  )
) +
  geom_pointrange(
    aes(
      y = median_rel_abund, 
      ymin = min_rel_abund, 
      ymax = max_rel_abund
    ),
    fatten = 0.5 
  ) +
  scale_color_manual(values = bar_color_palette) +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "Median, Min, and Max Relative Abundance (Core Genera)",
    subtitle = "Point is median, line is range (min/max) of monthly abundances.",
    x = "MAG ID (Colony)",
    y = "Relative Abundance (log10 scale)",
    color = "Genus" 
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# --- 3. Print the plot ---
print(range_abundance_plot_core)


# --- 1. Filter 'plot_data_snv_range' for Core Genera ---
snv_plot_data_core <- plot_data_snv_range %>%
  filter(Color_Genus %in% target_genera)

# --- 2. Create the Filtered SNV Density Plot ---
snv_range_plot_core <- ggplot(
  snv_plot_data_core, # Use the filtered data
  aes(
    x = reorder(mag_colony_id, -median_rel_abund),
    color = Color_Genus
  )
) +
  geom_pointrange(
    aes(
      y = median_snv_density, 
      ymin = min_snv_density, 
      ymax = max_snv_density
    ),
    fatten = 0.5 
  ) +
  scale_color_manual(values = bar_color_palette) +
  # You might want the log scale here too, if the data is skewed
  # scale_y_continuous(trans = "log10") + 
  labs(
    title = "Median, Min, and Max SNV Density (Core Genera)",
    subtitle = "Point is median, line is range (min/max) of *monthly* density. Sorted by median abundance.",
    x = "MAG ID (Colony)",
    y = "SNVs per Mbp",
    color = "Genus"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )+
  scale_y_continuous(trans="log10")

# --- 3. Print the plot ---
print(snv_range_plot_core)





#### QUICK GROUPINGS OF REGRESSION 

# --- 1. Create ONE Master Data Frame for Analysis ---
# We join your two main data objects to get all columns in one place.
# 'pointrange_data' has the Abundance min/median/max
# 'plot_data_snv_range' has the SNV min/median/max and the abundance-based sort order

# We select the key SNV columns from plot_data_snv_range
snv_data_to_join <- plot_data_snv_range %>%
  select(
    mag_colony_id, 
    median_snv_density, min_snv_density, max_snv_density
  )

# Now, join this to your full 'pointrange_data'
master_plot_data <- pointrange_data %>%
  inner_join(snv_data_to_join, by = "mag_colony_id")

cat("--- Master data frame created. First few rows: ---\n")
print(head(master_plot_data))

# --- 2. Run Regression and Get Residuals ---
# Now we run the model on this complete data frame
snv_model <- lm(median_snv_density ~ median_rel_abund, data = master_plot_data)

# Use broom::augment() and join the residuals back
# This ensures we keep ALL original columns
data_with_residuals <- broom::augment(snv_model, data = master_plot_data)

# --- 3. Identify Group 1: "Relatively Normal" ---
cat("\n--- Group 1: Top 25 'Normal' Relationship MAGs ---\n")

group1_normal_relationship <- data_with_residuals %>%
  arrange(abs(.resid)) %>%
  slice_head(n = 25)
# No 'select()' needed, we just keep everything

print(head(group1_normal_relationship, 2))

# --- 4. Identify Group 2: "Normal Abundance, Low SNVs" ---
cat("\n--- Group 2: Top 25 'Normal Abundance, Low SNV' MAGs ---\n")

abundance_quantiles <- quantile(data_with_residuals$median_rel_abund, probs = c(0.25, 0.75), na.rm = TRUE)
q25_abund <- abundance_quantiles[1]
q75_abund <- abundance_quantiles[2]

group2_low_snvs <- data_with_residuals %>%
  filter(
    median_rel_abund >= q25_abund,
    median_rel_abund <= q75_abund
  ) %>%
  arrange(.resid) %>%
  slice_head(n = 25)
# No 'select()' needed

print(head(group2_low_snvs, 2))

# --- 5. Prepare Data for Faceting ---
# This is much simpler now
group1_data <- group1_normal_relationship %>%
  mutate(group = "Group 1: Normal Relationship")

group2_data <- group2_low_snvs %>%
  mutate(group = "Group 2: Low SNVs")

# This single data frame has ALL columns for BOTH plots
faceted_plot_data <- bind_rows(group1_data, group2_data)

# --- 6. Plot 1: Faceted Abundance (This will work) ---
faceted_abundance_plot <- ggplot(
  faceted_plot_data,
  aes(
    x = reorder(mag_colony_id, -median_rel_abund),
    color = Color_Genus
  )
) +
  geom_pointrange(
    aes(
      y = median_rel_abund,
      ymin = min_rel_abund,
      ymax = max_rel_abund
    ),
    fatten = 0.5
  ) +
  scale_color_manual(values = bar_color_palette) +
  scale_y_continuous(trans = "log10") 
  facet_wrap(~ group, scales = "free_x") +
  labs(
    title = "Relative Abundance for Selected MAGs",
    subtitle = "Point is median, line is range (min/max) of monthly abundances.",
    x = "MAG ID (Colony)",
    y = "Relative Abundance (log10 scale)",
    color = "Genus"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+

print(faceted_abundance_plot)

# --- 7. Plot 2: Faceted SNV Density (This will work) ---
faceted_snv_plot <- ggplot(
  faceted_plot_data, 
  aes(
    x = reorder(mag_colony_id, -median_rel_abund),
    color = Color_Genus
  )
) +
  geom_pointrange(
    aes(
      y = median_snv_density,
      ymin = min_snv_density, # This column now exists
      ymax = max_snv_density  # This column now exists
    ),
    fatten = 0.5
  ) +
  scale_color_manual(values = bar_color_palette) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~ group, scales = "free_x") +
  labs(
    title = "SNV Density for Selected MAGs",
    subtitle = "Point is median, line is range (min/max) of *monthly* density. Sorted by median abundance.",
    x = "MAG ID (Colony)",
    y = "SNVs per Mbp",
    color = "Genus"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
print(faceted_snv_plot)



######### SAME AS ABOVE, NOT PLOTTING MEDIAN ########



# --- 1. Define Core Genera (Unchanged) ---
target_genera <- c(
  "bartonella", "snodgrassella", "gilliamella", 
  "lactobacillus", "bombilactobacillus", "bifidobacterium"
)

# --- 2. Find "Complete" MAGs (Unchanged) ---
months_in_colony <- final_abundances %>%
  group_by(Colony) %>%
  summarise(total_months = n_distinct(Month), .groups = 'drop')

mag_month_counts <- final_abundances %>%
  group_by(Colony, MAG_id) %>%
  summarise(mag_months = n_distinct(Month), .groups = 'drop')

completeness_key <- mag_month_counts %>%
  left_join(months_in_colony, by = "Colony") %>%
  filter(mag_months == total_months) %>%
  select(Colony, MAG_id) %>%
  mutate(complete_key = paste0(Colony, "_", MAG_id))

cat(paste("--- Found", nrow(completeness_key), "complete MAG/Colony pairs ---\n"))

# --- 3. Create One Filtered *Monthly* Data Frame ---
# Combine abundance and SNV data, then filter
monthly_data_clean <- final_abundances %>%
  # Join with monthly SNV data
  left_join(all_monthly_snv_density, by = c("MAG_id", "Colony", "Month")) %>%
  # Handle 0-SNV months
  mutate(snv_per_mbp = replace_na(snv_per_mbp, 0)) %>%
  # Create Color_Genus and mag_colony_id
  mutate(
    Color_Genus = case_when(
      tolower(Genus) %in% target_genera ~ tolower(Genus),
      TRUE ~ "Other"
    ),
    Color_Genus = factor(Color_Genus, levels = c(target_genera, "Other")),
    mag_colony_id = paste0(MAG_id, " (Colony ", Colony, ")")
  ) %>%
  # a. Filter for Core Genera
  filter(Color_Genus %in% target_genera) %>%
  # b. Filter for Completeness
  mutate(complete_key = paste0(Colony, "_", MAG_id)) %>%
  filter(complete_key %in% completeness_key$complete_key)

cat(paste("--- Created clean monthly dataset with", n_distinct(monthly_data_clean$mag_colony_id), "MAGs ---\n"))



# --- 1. Create a Summary of the Clean Data ---
# (This is from your previous step, unchanged)
summary_data_clean <- monthly_data_clean %>%
  group_by(mag_colony_id, Color_Genus) %>%
  summarise(
    median_rel_abund = median(rel_abund_mag, na.rm = TRUE),
    min_rel_abund = min(rel_abund_mag, na.rm = TRUE),
    max_rel_abund = max(rel_abund_mag, na.rm = TRUE),
    median_snv_density = median(snv_per_mbp, na.rm = TRUE),
    min_snv_density = min(snv_per_mbp, na.rm = TRUE),
    max_snv_density = max(snv_per_mbp, na.rm = TRUE),
    .groups = 'drop'
  )

# --- NEW FILTERING STEP ---
# Filter the summary data before building the model
# We use max_snv_density > 0 to keep MAGs that had *at least one* SNV
# in *at least one* month.
summary_data_clean_filtered <- summary_data_clean %>%
  filter(max_snv_density > 0)

cat(paste("--- Kept", nrow(summary_data_clean_filtered), 
          "clean, complete, core MAGs that have > 0 intermediate SNVs ---\n"))

# --- 2. Re-build the Model (on filtered data) ---
snv_model_clean <- lm(median_snv_density ~ median_rel_abund, data = summary_data_clean_filtered)
data_with_residuals_clean <- broom::augment(snv_model_clean, data = summary_data_clean_filtered)

# --- 3. Re-Identify Group 1: "Relatively Normal" ---
cat("--- Group 1: Top 25 'Normal' (from SNV > 0 data) ---\n")
group1_new <- data_with_residuals_clean %>%
  arrange(abs(.resid)) %>%
  slice_head(n = 25) %>%
  mutate(group = "Group 1: Normal Relationship")

# --- 4. Re-Identify Group 2: "Normal Abundance, Low SNVs" ---
cat("--- Group 2: Top 25 'Low SNV' (from SNV > 0 data) ---\n")
# Define "normal abundance" from this new filtered dataset
abundance_quantiles_clean <- quantile(data_with_residuals_clean$median_rel_abund, probs = c(0.25, 0.75), na.rm = TRUE)
q25_abund_clean <- abundance_quantiles_clean[1]
q75_abund_clean <- abundance_quantiles_clean[2]

group2_new <- data_with_residuals_clean %>%
  filter(
    median_rel_abund >= q25_abund_clean,
    median_rel_abund <= q75_abund_clean
  ) %>%
  # We still arrange by .resid (most negative)
  arrange(.resid) %>% 
  slice_head(n = 25) %>%
  mutate(group = "Group 2: Low SNVs")

# --- 5. Create Final Datasets for Plotting ---
# a. The summary data for the linerange
linerange_data_final <- bind_rows(group1_new, group2_new)

# b. The monthly data for the points
plot_data_final_filtered <- monthly_data_clean %>%
  filter(mag_colony_id %in% linerange_data_final$mag_colony_id) %>%
  left_join(
    linerange_data_final %>% select(mag_colony_id, group, median_rel_abund), 
    by = "mag_colony_id"
  )


faceted_abundance_plot_final <- ggplot(
  linerange_data_final,
  aes(
    x = reorder(mag_colony_id, -median_rel_abund),
    color = Color_Genus
  )
) +
  geom_linerange(
    aes(ymin = min_rel_abund, ymax = max_rel_abund),
    alpha = 0.3, linewidth = 1 
  ) +
  geom_point(
    data = plot_data_final_filtered,
    aes(y = rel_abund_mag),
    alpha = 0.7
  ) +
  scale_color_manual(values = bar_color_palette) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~ group, scales = "free_x") +
  labs(
    title = "Monthly Relative Abundance (SNV > 0 Filtered)",
    subtitle = "Points are monthly measurements; vertical line shows full range (min/max).",
    x = "MAG ID (Colony)",
    y = "Relative Abundance (log10 scale)",
    color = "Genus"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(faceted_abundance_plot_final)



faceted_snv_plot_final <- ggplot(
  linerange_data_final,
  aes(
    x = reorder(mag_colony_id, -median_rel_abund),
    color = Color_Genus
  )
) +
  geom_linerange(
    aes(ymin = min_snv_density, ymax = max_snv_density),
    alpha = 0.3, linewidth = 1
  ) +
  geom_point(
    data = plot_data_final_filtered,
    aes(y = snv_per_mbp),
    alpha = 0.7
  ) +
  scale_color_manual(values = bar_color_palette) +
  facet_wrap(~ group, scales = "free_x") +
  labs(
    title = "Monthly SNV Density (SNV > 0 Filtered)")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_continuous(trans="log10")
  
faceted_snv_plot_final
