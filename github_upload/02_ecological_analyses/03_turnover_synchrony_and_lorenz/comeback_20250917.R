# --- 1. Load Libraries and Prepare Base Data ---
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstatix)
library(ggpubr)
library(scales)
library(purrr)
library(rstatix)
library(mclust)
library(corrplot)
ecology_data <- read.csv('ecology_metadata_df.csv')
ecology_data <- ecology_data[ecology_data$Caste=="Worker",]
# Apply the time function (assuming 'ecology_data' is loaded)
add_cumulative_time_from_month_nucdiv <- function(data) {
  month_base_times <- c("May" = 0, "June" = 36, "July" = 31, "August" = 34, 
                        "September" = 28, "October" = 22, "November" = 42, 
                        "January" = 47, "February" = 31)
  cumulative_times <- cumsum(month_base_times)
  month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))
  data$Time <- month_to_cumulative_time_map[data$Month]
  return(data)
}
ecology_data <- add_cumulative_time_from_month_nucdiv(ecology_data)


# --- 2. Create the list of data subsets ---
data_subsets <- list(
  Microbial = filter(ecology_data, type == "mmag"),
  Viral = filter(ecology_data, type == "vmag"),
  Combined = filter(ecology_data, type %in% c("mmag", "vmag"))
)



# --- 3. Loop Through Each Community Type to Run the Analysis ---
# We will store the final results for plots and tables in these lists
all_community_decay_data <- list()
all_community_mantel_results <- list()

for (community_name in names(data_subsets)) {
  
  message(paste("--- Analyzing:", community_name, "Community ---"))
  
  current_data <- data_subsets[[community_name]]
  
  # A. Create the community matrix
  community_matrix <- current_data %>%
    group_by(sample_name, MAG_id) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = MAG_id, values_from = total_abundance, values_fill = 0) %>%
    tibble::column_to_rownames("sample_name")
  
  # B. Align metadata
  metadata <- current_data %>%
    select(sample_name, Colony, Time) %>%
    distinct(sample_name, .keep_all = TRUE) %>%
    filter(sample_name %in% rownames(community_matrix))
  
  # C. Calculate decay data and run Mantel tests for each colony
  unique_colonies <- unique(metadata$Colony)
  colony_decay_data <- list()
  colony_mantel_results <- list()
  
  for (col in unique_colonies) {
    meta_subset <- metadata %>% filter(Colony == col)
    if (nrow(meta_subset) < 3) next
    
    comm_subset <- community_matrix[meta_subset$sample_name, ]
    bray_dist <- vegdist(comm_subset, method = "bray")
    time_lag_dist <- dist(meta_subset$Time, method = "euclidean")
    
    # Run Mantel test
    mantel_res <- mantel(bray_dist, time_lag_dist, permutations = 999)
    colony_mantel_results[[as.character(col)]] <- data.frame(Colony = col, Mantel_r = mantel_res$statistic, p_value = mantel_res$signif)
    
    # Tidy data for plotting
    bray_sim <- 1 - as.matrix(bray_dist)
    bray_sim[upper.tri(bray_sim, diag = TRUE)] <- NA
    time_lag_matrix <- as.matrix(time_lag_dist)
    time_lag_matrix[upper.tri(time_lag_matrix, diag = TRUE)] <- NA
    
    decay_df <- data.frame(Colony = col, Similarity = as.vector(bray_sim), TimeLag = as.vector(time_lag_matrix)) %>% na.omit()
    colony_decay_data[[as.character(col)]] <- decay_df
  }
  
  # D. Store the final results for this community type
  all_community_decay_data[[community_name]] <- bind_rows(colony_decay_data)
  all_community_mantel_results[[community_name]] <- bind_rows(colony_mantel_results)
}

# --- 4. View Results ---
# You can now access the results for each community type
print("--- Mantel Results for Microbial Community ---")
print(all_community_mantel_results$Microbial)

print("--- Mantel Results for Viral Community ---")
print(all_community_mantel_results$Viral)

# --- 5. Plotting Example (for the 'Combined' community) ---
ggplot(all_community_decay_data$Microbial, aes(x = TimeLag, y = Similarity)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ Colony, scales = "free") +
  labs(
    title = "Similarity Decay of the Combined Microbial",
    x = "Time Lag Between Samples (Days)",
    y = "Bray-Curtis Similarity"
  ) +
  theme_bw()

plot_data <- all_community_decay_data$Microbial %>%
  mutate(Dissimilarity = 1 - Similarity)

# --- 3. Prepare the Annotation Data ---
# We'll use the Mantel test results to create labels for each facet
label_data <- all_community_mantel_results$Microbial %>%
  mutate(
    # Create a clean text label for the plot
    label = paste0("Mantel R = ", round(Mantel_r, 2), "\np = ", round(p_value, 3)),
    # Add a column to check for significance at the p < 0.05 level
    is_significant = ifelse(p_value < 0.05, TRUE, FALSE)
  )

# --- 4. Create the Plot ---
similarity_decay_plot <- ggplot(plot_data, aes(x = TimeLag, y = Dissimilarity)) +
  # Add the points for each pairwise comparison
  geom_point(alpha = 1) +
  # Add a linear trend line. We'll color it based on the Mantel test significance.
  geom_smooth(method = "lm", se = FALSE, aes(color = "black"), linetype = "solid") +
  # Add the text annotations from our label_data
  geom_text(
    data = label_data,
    aes(x = Inf, y = Inf, label = label), # Place text in top-right corner
    hjust = 1.1, vjust = 1.2, size = 3
  ) +
  # Create a separate panel for each colony
  facet_wrap(~ Colony, scales = "free") +
  # Manually set the color for significant (red) vs. not significant (black)
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
  labs(
    title = "Temporal Decay of Microbial Communities",
    subtitle = "Bray-Curtis Dissimilarity vs. Time Lag for each Colony",
    x = "Time Lag Between Samples (Days)",
    y = "Bray-Curtis Dissimilarity (Higher = More Different)"
  ) +
  theme_bw()

print(similarity_decay_plot)

library(mgcv)
# --- 2. Prepare the Data ---
# (Assuming 'all_community_decay_data' is in your environment)
gam_decay_data <- all_community_decay_data$Microbial %>%
  mutate(Dissimilarity = 1 - Similarity)

# For Beta regression, the response must be > 0 and < 1.
# We "squeeze" the data slightly away from the boundaries.
n <- nrow(gam_decay_data)
gam_decay_data <- gam_decay_data %>%
  mutate(
    Dissimilarity_squeezed = (Dissimilarity * (n - 1) + 0.5) / n,
    Colony = as.factor(Colony)
  )

# --- 3. Fit the GAM ---
decay_gam <- gam(
  Dissimilarity_squeezed ~ s(TimeLag) + s(Colony, bs = "re"),
  data = gam_decay_data,
  family = betar(link = "logit"),
  method = "REML"
)

# --- 4. View Results ---
summary(decay_gam)
plot(decay_gam, pages = 1)



  
  # --- 2. Create the Plot ---
  # (Assuming 'mmag_alpha_df' is in your environment)
  alpha_season_plot <- ggplot(mmag_alpha_df, aes(x = Season, y = Shannon, fill = Season)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.4) +
    labs(
      title = "Alpha Diversity Across Seasons",
      x = "Season",
      y = "Shannon Diversity"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  print(alpha_season_plot)
  
  # --- 3. Run the Statistical Tests ---
  # Kruskal-Wallis test
  kruskal_result <- kruskal.test(Shannon ~ Season, data = mmag_alpha_df)
  cat("\n--- Kruskal-Wallis Test for Shannon Diversity by Season ---\n")
  print(kruskal_result)
  
  # Dunn's post-hoc test (if Kruskal-Wallis is significant)
  dunn_result <- mmag_alpha_df %>%
    dunn_test(Shannon ~ Season, p.adjust.method = "fdr")
  cat("\n--- Dunn's Post-Hoc Test ---\n")
  print(dunn_result)
  
  
  
  
  
  # --- 2. Prepare the Genus-Level Community Matrix ---
  # (Assuming 'ecology_data' is loaded)
  mmag_data <- filter(ecology_data, type == "mmag")
  
  # The key change is here: we group by 'Genus' instead of 'MAG_id'
  community_matrix_genus <- mmag_data %>%
    group_by(sample_name, Genus) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Genus, values_from = total_abundance, values_fill = 0) %>%
    tibble::column_to_rownames("sample_name")
  
  # Create the corresponding metadata file
  metadata <- mmag_data %>%
    select(sample_name, Season) %>%
    distinct() %>%
    filter(sample_name %in% rownames(community_matrix_genus))
  
  # --- 3. Alpha Diversity Analysis (by Genus) ---
  
  # a. Calculate Shannon diversity on the new genus-level matrix
  shannon_genus <- diversity(community_matrix_genus, index = "shannon")
  alpha_div_genus_df <- data.frame(
    Shannon = shannon_genus,
    sample_name = rownames(community_matrix_genus)
  ) %>%
    inner_join(metadata, by = "sample_name")
  
  # b. Create the violin plot
  alpha_season_plot_genus <- ggplot(alpha_div_genus_df, aes(x = Season, y = Shannon, fill = Season)) +
    geom_boxplot() +
    geom_jitter(width = 0.05, alpha = 0.9) +
    labs(title = "Alpha Diversity (Genus Level) Across Seasons", x = "Season", y = "Shannon Diversity") +
    theme_bw() +
    theme(legend.position = "none")
  
  print(alpha_season_plot_genus)
  
  # c. Run statistical tests
  kruskal_result_genus <- kruskal.test(Shannon ~ Season, data = alpha_div_genus_df)
  cat("\n--- Kruskal-Wallis Test (Genus Level) ---\n")
  print(kruskal_result_genus)
  
  dunn_result_genus <- alpha_div_genus_df %>%
    dunn_test(Shannon ~ Season, p.adjust.method = "fdr")
  cat("\n--- Dunn's Post-Hoc Test (Genus Level) ---\n")
  print(dunn_result_genus)
  
  # --- 4. Beta Diversity Analysis (by Genus) ---
  
  # a. Calculate Bray-Curtis distance and run PERMANOVA
  bray_dist_genus <- vegdist(community_matrix_genus, method = "bray")
  permanova_result_genus <- adonis2(bray_dist_genus ~ Season, data = metadata, permutations = 999)
  cat("\n--- PERMANOVA Test (Genus Level) ---\n")
  print(permanova_result_genus)
  
  # b. Create NMDS plot
  nmds_genus <- metaMDS(community_matrix_genus, distance = "bray")
  nmds_scores_genus <- as.data.frame(scores(nmds_genus))
  nmds_scores_genus$sample_name <- rownames(nmds_scores_genus)
  nmds_plot_data_genus <- inner_join(nmds_scores_genus, metadata, by = "sample_name")
  
  beta_season_plot_genus <- ggplot(nmds_plot_data_genus, aes(x = NMDS1, y = NMDS2, color = Season)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Season), type = "t") +
    labs(title = "NMDS of Communities (Genus Level) by Season") +
    theme_bw()
  
  print(beta_season_plot_genus)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  --- Prepare Data ---
    # Define the desired order for the seasons
    season_order <- c("Spring", "Summer", "Fall", "Winter")
  
  # Reorder the 'Season' factor in the dataframe
  plot_data_mag <- mmag_alpha_df %>%
    mutate(Season = factor(Season, levels = season_order))
  
  # --- Run Statistical Tests ---
  # We need the Dunn's test results to add to the plot
  dunn_result_mag <- plot_data_mag %>%
    dunn_test(Shannon ~ Season, p.adjust.method = "fdr")
  
  # --- Create the Plot with Significance Brackets ---
  alpha_plot_mag <- ggplot(plot_data_mag, aes(x = Season, y = Shannon, fill = Season)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    # Add significance brackets using the Dunn's test results
    stat_pvalue_manual(
      dunn_result_mag, 
      y.position = 4.0, # Adjust this value to move brackets up or down
      label = "p.adj.signif",
      tip.length = 0.01
    ) +
    labs(
      title = "Alpha Diversity (MAG Level) Across Seasons",
      x = "Season",
      y = "Shannon Diversity"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  print(alpha_plot_mag)# --- Prepare Data ---
# (Assuming 'alpha_div_genus_df' is in your environment)
plot_data_genus <- alpha_div_genus_df %>%
  mutate(Season = factor(Season, levels = season_order))

# --- Run Statistical Tests ---
dunn_result_genus <- plot_data_genus %>%
  dunn_test(Shannon ~ Season, p.adjust.method = "fdr")

# --- Create the Plot with Significance Brackets ---
alpha_plot_genus <- ggplot(plot_data_genus, aes(x = Season, y = Shannon, fill = Season)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  # Add significance brackets
  stat_pvalue_manual(
    dunn_result_genus, 
    y.position = 3.5, # Adjust this value as needed
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  labs(
    title = "Alpha Diversity (Genus Level) Across Seasons",
    x = "Season",
    y = "Shannon Diversity"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(alpha_plot_genus)
  
  
  
  
  
 








# --- Prepare Data ---
# Define the desired order for the seasons
season_order <- c("Spring", "Summer", "Fall", "Winter")

# Reorder the 'Season' factor in the dataframe
plot_data_mag <- mmag_alpha_df %>%
  mutate(Season = factor(Season, levels = season_order))

# --- Run Statistical Tests ---
# We need the Dunn's test results to add to the plot
dunn_result_mag <- plot_data_mag %>%
  dunn_test(Shannon ~ Season, p.adjust.method = "fdr")
# --- Create the Plot with Significance Brackets ---
alpha_plot_mag <- ggplot(plot_data_mag, aes(x = Season, y = Shannon, fill = Season)) +
  geom_boxplot(trim = FALSE) +
  geom_jitter(width = 0.05, alpha = 0.9) +
  # Add significance brackets using the Dunn's test results
  stat_pvalue_manual(
    dunn_result_mag, 
    y.position = 4.0, # Adjust this value to move brackets up or down
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  labs(
    title = "Alpha Diversity (MAG Level) Across Seasons",
    x = "Season",
    y = "Shannon Diversity"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(alpha_plot_mag)


# --- Prepare Data ---
# (Assuming 'alpha_div_genus_df' is in your environment)
plot_data_genus <- alpha_div_genus_df %>%
  mutate(Season = factor(Season, levels = season_order))

# --- Run Statistical Tests ---
dunn_result_genus <- plot_data_genus %>%
  dunn_test(Shannon ~ Season, p.adjust.method = "fdr")

# --- Create the Plot with Beeswarm Overlay ---
alpha_plot_genus_beeswarm <- ggplot(plot_data_genus, aes(x = Season, y = Shannon, fill = Season)) +
  geom_boxplot(alpha = 1) +
  # Replace geom_jitter() with geom_quasirandom()
  geom_quasirandom(width = 0.1, alpha = 0.9) + 
  stat_pvalue_manual(
    dunn_result_mag, 
    y.position = 4.0, 
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  labs(
    title = "Alpha Diversity (Genus Level) Across Seasons",
    x = "Season",
    y = "Shannon Diversity"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(alpha_plot_genus_beeswarm)

  
  
  
  
  
# --- Create the Plot with Beeswarm Overlay ---
alpha_plot_mag_beeswarm <- ggplot(plot_data_mag, aes(x = Season, y = Shannon, fill = Season)) +
  geom_boxplot(alpha = 1) +
  # Replace geom_jitter() with geom_quasirandom()
  geom_quasirandom(width = 0.1, alpha = 0.9) + 
  stat_pvalue_manual(
    dunn_result_mag, 
    y.position = 4.0, 
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  labs(
    title = "Alpha Diversity (MAG Level) Across Seasons",
    x = "Season",
    y = "Shannon Diversity"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(alpha_plot_mag_beeswarm)
  
  

# --- 1. Load Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan) # For Bray-Curtis distance

# --- 2. Calculate Genus Relative Abundance ---
genus_abund <- ecology_data %>%
  # <-- FIX: Remove rows with empty or NA Genus names
  filter(Genus != "" & !is.na(Genus)) %>% 
  group_by(sample_name, Genus) %>%
  summarise(absolute_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
  group_by(sample_name) %>%
  mutate(relative_abundance = absolute_abundance / sum(absolute_abundance)) %>%
  ungroup()

# --- 3. Identify Top 6 Genera and Prepare for Plotting ---
top_6_genera <- genus_abund %>%
  group_by(Genus) %>%
  summarise(mean_rel_abund = mean(relative_abundance)) %>%
  arrange(desc(mean_rel_abund)) %>%
  head(6) %>%
  pull(Genus)

plot_data <- genus_abund %>%
  mutate(Genus_colored = ifelse(Genus %in% top_6_genera, Genus, "Other"))

# --- 4. Hierarchical Clustering to Order Samples ---
community_matrix <- plot_data %>%
  pivot_wider(
    id_cols = sample_name, 
    names_from = Genus, 
    values_from = relative_abundance, 
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("sample_name")

bray_dist <- vegdist(community_matrix, method = "bray")
hclust_result <- hclust(bray_dist, method = "ward.D2")
sample_order <- hclust_result$label[hclust_result$order]

# --- 5. Create the Final Plot ---
plot_data$sample_name <- factor(plot_data$sample_name, levels = sample_order)

color_palette <- c(
  "commensalibacter" = "#a6cee3", "lactobacillus" = "#1f78b4", "snodgrassella" = "#b2df8a",
  "gilliamella" = "#33a02c", "bartonella" = "#fb9a99", "bombilactobacillus" = "#e31a1c",
  "Other" = "grey80"
)

abundance_plot <- ggplot(plot_data, aes(x = sample_name, y = relative_abundance, fill = Genus_colored)) +
  geom_col(width = 1) +
  scale_fill_manual(values = color_palette, name = "Genus") +
  labs(
    title = "Relative Abundance of Microbial Genera",
    subtitle = "Samples are ordered by hierarchical clustering (Ward.D2 method)",
    x = "Sample",
    y = "Relative Abundance"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(abundance_plot)



# --- 3. Create Comprehensive Annotation Data ---
# Select all metadata variables we need and apply the sample order
annotation_data <- ecology_data %>%
  select(sample_name, Apiary, Colony, Season, Caste) %>%
  distinct() %>%
  mutate(sample_name = factor(sample_name, levels = sample_order))

# --- 4. Create the Individual Plots ---
# a. Main abundance plot (modified to have no legends initially)
plot_data$sample_name <- factor(plot_data$sample_name, levels = sample_order)
color_palette <- c("commensalibacter"="#a6cee3", "lactobacillus"="#1f78b4", "snodgrassella"="#b2df8a", "gilliamella"="#33a02c", "bartonella"="#fb9a99", "bombilactobacillus"="#e31a1c", "Other"="grey80")

abundance_plot <- ggplot(plot_data, aes(x = sample_name, y = relative_abundance, fill = Genus_colored)) +
  geom_col(width = 1) +
  scale_fill_manual(values = color_palette, name = "Genus") +
  labs(x = "Sample", y = "Relative Abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

# b. Create a plot for each annotation bar
p_apiary <- ggplot(annotation_data, aes(x = sample_name, y = 1, fill = Apiary)) + geom_tile() + theme_void() + scale_fill_brewer(palette = "Set2", name = "Apiary")
p_season <- ggplot(annotation_data, aes(x = sample_name, y = 1, fill = Season)) + geom_tile() + theme_void() + scale_fill_brewer(palette = "Set3", name = "Season")
p_caste <- ggplot(annotation_data, aes(x = sample_name, y = 1, fill = Caste)) + geom_tile() + theme_void() + scale_fill_brewer(palette = "Dark2", name = "Caste")

# --- 5. Combine All Plots with Patchwork ---
final_plot_annotated <- p_caste + p_season  + p_apiary + abundance_plot +
  plot_layout(
    ncol = 1, 
    heights = c(0.5, 0.5, 0.5, 12), # Give each bar a height of 1 and the main plot a height of 12
    guides = 'collect' # Collect all legends into a single panel
  ) +
  plot_annotation(
    title = 'Relative Abundance of Microbial Genera',
    subtitle = 'Samples are ordered by hierarchical clustering'
  ) & # The '&' applies the theme to all subplots
  theme(legend.position = "bottom")

print(final_plot_annotated)
# --- 3. Loop Through Each Community Type to Run Analyses ---
anosim_results <- list()

for (community_name in names(data_subsets)) {
  
  message(paste("\n--- Analyzing Phenotype Dissimilarity for:", community_name, "Community ---"))
  
  current_data <- data_subsets[[community_name]]
  
  # A. Create community matrix
  community_matrix <- current_data %>%
    group_by(sample_name, MAG_id) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = MAG_id, values_from = total_abundance, values_fill = 0) %>%
    tibble::column_to_rownames("sample_name")
  
  # B. Align metadata and ensure Phenotype is a factor
  metadata <- current_data %>%
    select(sample_name, Phenotype) %>%
    distinct(sample_name, .keep_all = TRUE) %>%
    filter(sample_name %in% rownames(community_matrix)) %>%
    mutate(Phenotype = factor(Phenotype))
  
  # C. Calculate the full dissimilarity matrix
  bray_dist <- vegdist(community_matrix, method = "bray")
  
  # D. Run ANOSIM
  # Ensure there are at least two groups to compare
  if (nlevels(metadata$Phenotype) > 1) {
    anosim_res <- anosim(bray_dist, metadata$Phenotype, permutations = 999)
    
    anosim_results[[community_name]] <- data.frame(
      Community = community_name,
      ANOSIM_R = anosim_res$statistic,
      p_value = anosim_res$signif
    )
  } else {
    message("Skipping ANOSIM: Not enough phenotype groups in the data.")
  }
}

# --- Print the final ANOSIM results ---
final_anosim_results <- bind_rows(anosim_results)
print(final_anosim_results)


ancova_results <- list()

# Re-using the main loop structure and data from previous steps...
for (community_name in names(data_subsets)) {
  
  message(paste("\n--- Analyzing Differences in Slopes for:", community_name, "Community ---"))
  
  # Use the tidy pairwise data we generated before (re-creating it here for clarity)
  # This dataframe has one row per unique pair of samples
  tidy_pairwise_data <- all_community_decay_data[[community_name]] %>%
    # We want Dissimilarity on the y-axis, not similarity
    mutate(Dissimilarity = 1 - Similarity) %>%
    # Ensure Colony is a factor for the model
    mutate(Colony = factor(Colony))
  
  # Fit the ANCOVA model
  # We are testing if Dissimilarity can be explained by TimeLag, Colony, and their interaction
  if (nrow(tidy_pairwise_data) > 0 && nlevels(tidy_pairwise_data$Colony) > 1) {
    ancova_model <- lm(Dissimilarity ~ TimeLag * Colony, data = tidy_pairwise_data)
    
    # Get the ANOVA table for the model
    ancova_table <- anova(ancova_model)
    
    # Store the result for the interaction term
    ancova_results[[community_name]] <- list(
      Community = community_name,
      ANCOVA_Table = ancova_table
    )
  }
}

# --- Print and interpret the ANCOVA results ---
for (community_name in names(ancova_results)) {
  cat(paste("\n--- ANCOVA Results for", community_name, "Community ---\n"))
  print(ancova_results[[community_name]]$ANCOVA_Table)
  cat("Look for the 'TimeLag:Colony' interaction term. A p-value < 0.05 suggests that the rates of temporal change are significantly different among the colonies.\n")
}



slope_results <- list()

for (community_name in names(data_subsets)) {
  
  message(paste("\n--- Calculating Slopes by Apiary & Overall for:", community_name, "Community ---"))
  
  # We need metadata with Apiary information
  metadata_full <- data_subsets[[community_name]] %>%
    select(sample_name, Colony, Apiary, Time) %>%
    distinct(sample_name, .keep_all = TRUE)
  
  # Use the tidy pairwise data and add Apiary info
  tidy_pairwise_data <- all_community_decay_data[[community_name]] %>%
    mutate(Dissimilarity = 1 - Similarity)
  
  # A. Overall Analysis (Agnostic to Group)
  overall_model <- lm(Dissimilarity ~ TimeLag, data = tidy_pairwise_data)
  overall_slope <- coef(overall_model)[2]
  overall_p_value <- summary(overall_model)$coefficients[2, 4]
  
  # B. Apiary-level Analysis
  message("Correcting join and running Apiary-level analysis...")
  
  # Create a clean mapping of each Colony to its Apiary
  colony_to_apiary_map <- metadata_full %>%
    select(Colony, Apiary) %>%
    distinct()
  
  # Join the Apiary info to the pairwise data using the shared 'Colony' column
  # This adds the correct 'Apiary' label to each pairwise comparison
  pairs_with_apiary <- tidy_pairwise_data %>%
    left_join(colony_to_apiary_map, by = "Colony")
  
  # Now, proceed with the analysis using the corrected dataframe
  apiary_slopes <- pairs_with_apiary %>%
    filter(!is.na(Apiary)) %>% # Safeguard against any failed matches
    group_by(Apiary) %>%
    summarise(
      slope = coef(lm(Dissimilarity ~ TimeLag))[2],
      p_value = summary(lm(Dissimilarity ~ TimeLag))$coefficients[2, 4],
      .groups = 'drop' 
    )
  
  slope_results[[community_name]] <- list(
    Overall = data.frame(Group = "Overall", slope = overall_slope, p_value = overall_p_value),
    By_Apiary = apiary_slopes
  )
}

# --- Print the slope results ---
for(community_name in names(slope_results)) {
  cat(paste("\n--- Temporal Change Slopes for", community_name, "Community ---\n"))
  print(slope_results[[community_name]]$Overall)
  print(slope_results[[community_name]]$By_Apiary)
}


# --- 1. Prepare the Data for Plotting ---
# Combine the pairwise results for Microbial and Viral communities into one dataframe
# The .id argument creates a new column with the names from the list ('Microbial', 'Viral')
comparison_data <- bind_rows(
  all_community_decay_data[c("Microbial", "Viral")], 
  .id = "CommunityType"
) %>%
  # Calculate Dissimilarity from the existing Similarity column
  mutate(Dissimilarity = 1 - Similarity)

# --- 2. Create the Comparison Plot ---
dissimilarity_comparison_plot <- ggplot(
  comparison_data, 
  aes(x = TimeLag, y = Dissimilarity, color = CommunityType)
) +
  # Add the individual data points with some transparency
  geom_point(alpha = 0.3) +
  # Add a linear regression line for each community type without the confidence interval
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  # Use a color-blind friendly palette
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Rate of Community Change: Viral vs. Microbial",
    subtitle = "Comparing the increase in Bray-Curtis Dissimilarity over time",
    x = "Time Lag Between Samples (Days)",
    y = "Bray-Curtis Dissimilarity",
    color = "Community Type"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(dissimilarity_comparison_plot)




# --- 3. Initialize lists to store all results ---
dispersion_plots <- list()
dispersion_summary_stats <- list()
dispersion_test_results <- list()

# --- 4. Main Analysis Loop ---
for (community_name in names(data_subsets)) {
  
  message(paste("\n--- Analyzing Dispersion for:", community_name, "Community ---"))
  
  current_data <- data_subsets[[community_name]]
  
  # A. Create community matrix
  community_matrix <- current_data %>%
    group_by(sample_name, MAG_id) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = MAG_id, values_from = total_abundance, values_fill = 0) %>%
    tibble::column_to_rownames("sample_name")
  
  # B. Data Cleaning: Remove samples with zero total abundance
  community_matrix <- community_matrix[rowSums(community_matrix) > 0, ]
  
  # C. Align metadata
  metadata <- current_data %>%
    select(sample_name, Colony) %>%
    distinct(sample_name, .keep_all = TRUE) %>%
    filter(sample_name %in% rownames(community_matrix)) %>%
    mutate(Colony = factor(Colony))
  
  # D. Robustness Check
  if (nrow(community_matrix) < 3 || nlevels(metadata$Colony) < 2) {
    message(paste("Skipping", community_name, "- not enough samples or groups for analysis."))
    next
  }
  
  # E. Calculate Bray-Curtis dissimilarity
  bray_dist <- vegdist(community_matrix, method = "bray")
  
  # F. Perform multivariate dispersion analysis
  betadisper_result <- betadisper(bray_dist, group = metadata$Colony)
  
  # G. Test for significant differences in dispersion
  permutest_result <- permutest(betadisper_result, permutations = 999)
  dispersion_test_results[[community_name]] <- permutest_result
  
  # H. Calculate summary statistics (Mean distance & CV)
  distances_df <- data.frame(distance_to_centroid = betadisper_result$distances, Colony = betadisper_result$group)
  summary_stats <- distances_df %>%
    group_by(Colony) %>%
    summarise(
      mean_distance = mean(distance_to_centroid),
      cv_percent = (sd(distance_to_centroid) / mean(distance_to_centroid)) * 100,
      .groups = 'drop'
    ) %>%
    arrange(mean_distance)
  dispersion_summary_stats[[community_name]] <- summary_stats
  
  # --- I. Generate the Ordination "Spider Plot" (CORRECTED METHOD) ---
  
  # 1. Perform PCoA explicitly to get sample coordinates
  pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
  sample_scores <- as.data.frame(pcoa_result$points)
  colnames(sample_scores) <- c("PCoA1", "PCoA2")
  sample_scores$Colony <- metadata$Colony
  
  # 2. Extract centroid coordinates from the betadisper result
  centroid_scores <- as.data.frame(betadisper_result$centroids[, c("PCoA1", "PCoA2")])
  centroid_scores$Colony <- rownames(centroid_scores)
  
  # 3. Join the two dataframes for plotting
  segment_data <- left_join(sample_scores, centroid_scores, by = "Colony", suffix = c("_sample", "_centroid"))
  
  # 4. Create the ggplot object
  p <- ggplot() +
    geom_segment(data = segment_data, aes(x = PCoA1_sample, y = PCoA2_sample, xend = PCoA1_centroid, yend = PCoA2_centroid, color = Colony), alpha = 0.5) +
    geom_point(data = sample_scores, aes(x = PCoA1, y = PCoA2, color = Colony), size = 2) +
    geom_point(data = centroid_scores, aes(x = PCoA1, y = PCoA2, color = Colony), size = 5, shape = 18) +
    labs(
      title = paste("Community Dispersion by Colony:", community_name),
      subtitle = "Samples (circles) are connected to their group centroid (diamond)",
      x = "Principal Coordinate 1", y = "Principal Coordinate 2", color = "Colony"
    ) +
    theme_bw() +
    coord_equal()
  
  # 5. Store the completed plot in the list
  dispersion_plots[[community_name]] <- p
}

# --- 5. Now you can view your results and plots ---
message("\n--- Dispersion Analysis Complete ---")

# Example: View the summary statistics for the Microbial community
print(dispersion_summary_stats$Microbial)

# Example: View the ordination plot for the Microbial community
print(dispersion_plots$Microbial)




# --- Function to calculate alpha diversity on a data subset ---
calculate_alpha_diversity_subset <- function(input_df, mag_type_filter, analysis_name) {
  
  print(paste("--- Running Alpha Diversity for:", analysis_name, "Community ---"))
  
  # 1. Filter the data for the specific MAG type
  filtered_data <- input_df %>%
    filter(type == mag_type_filter)
  
  if (nrow(filtered_data) == 0) {
    warning(paste("No data found for type:", mag_type_filter, ". Skipping."), call. = FALSE)
    return(NULL)
  }
  
  # 2. Create the community matrix (Samples x Genera)
  community_matrix <- filtered_data %>%
    group_by(sample_name, MAG_id) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = MAG_id, values_from = total_abundance, values_fill = 0) %>%
    tibble::column_to_rownames("sample_name")
  
  # Check for empty matrix
  if (nrow(community_matrix) == 0 || ncol(community_matrix) == 0) {
    warning(paste("Community matrix for", mag_type_filter, "is empty. Skipping."), call. = FALSE)
    return(NULL)
  }
  
  # 3. Calculate vegan alpha diversity metrics
  richness <- specnumber(community_matrix)
  shannon <- diversity(community_matrix, index = "shannon")
  simpson <- diversity(community_matrix, index = "simpson")
  
  # 4. Combine and merge with metadata
  alpha_div_df <- data.frame(
    Observed = richness,
    Shannon = shannon,
    Simpson = simpson,
    sample_name = rownames(community_matrix)
  )
  
  metadata_df <- filtered_data %>%
    select(sample_name, State, Month, Caste, Colony, Rep, Season, Phenotype, Apiary) %>%
    distinct(sample_name, .keep_all = TRUE)
  
  final_df <- merge(alpha_div_df, metadata_df, by = "sample_name")
  
  print(paste("--- Finished analysis for:", analysis_name, "---"))
  return(final_df)
}

# --- Execute the function for vMAGs and mMAGs ---

# Analyze the viral community 🔬
vmag_alpha_df <- calculate_alpha_diversity_subset(
  input_df = ecology_data, 
  mag_type_filter = "vmag", 
  analysis_name = "Viral"
)

# Analyze the microbial community 🦠
mmag_alpha_df <- calculate_alpha_diversity_subset(
  input_df = ecology_data, 
  mag_type_filter = "mmag", 
  analysis_name = "Microbial"
)


# --- View the results ---
cat("\n\n### Viral Community Alpha Diversity ###\n")
head(vmag_alpha_df)

cat("\n### Microbial Community Alpha Diversity ###\n")
head(mmag_alpha_df)


# --- 1. Define metrics and the specific factor to test ---
metrics_to_test <- c("Observed", "Shannon", "Simpson")
factor_to_test <- "Apiary"

# A list to store the results
results_list <- list()

# --- 2. Loop through metrics and perform Kruskal-Wallis tests ---
for (metric in metrics_to_test) {
  test_formula <- as.formula(paste(metric, "~", factor_to_test))
  
  # Ensure the factor has more than one level
  if (length(unique(vmag_alpha_df[[factor_to_test]])) > 1) {
    kruskal_result <- kruskal.test(test_formula, data = mmag_alpha_df)
    
    results_list[[metric]] <- list(
      H_statistic = kruskal_result$statistic,
      p_value = kruskal_result$p.value
    )
  }
}

# --- 3. Format results and adjust p-values ---
# Convert the list to a dataframe
results_df <- do.call(rbind, lapply(results_list, as.data.frame))
results_df$Metric <- rownames(results_df)

# Adjust p-values using the Benjamini-Hochberg (FDR) method
results_df$Adjusted_P_Value <- p.adjust(results_df$p_value, method = "fdr")

# Print the final results table
print(results_df[, c("Metric", "H_statistic", "p_value", "Adjusted_P_Value")])


# --- 1. Perform Dunn's post-hoc test ---
# This test will compare every pair of colonies for Shannon diversity
dunn_results <- mmag_alpha_df %>%
  dunn_test(Observed ~ Apiary, p.adjust.method = "fdr")

# --- 2. View the results ---
# We are most interested in pairs with a low p.adj value
print(dunn_results)



combined_alpha_data <- bind_rows(
  list(Viral = vmag_alpha_df, Microbial = mmag_alpha_df),
  .id = "Community"
)
combined_alpha_data <- add_cumulative_time_from_month_nucdiv(combined_alpha_data)

# Pivot the data into a tidy, long format for plotting
tidy_alpha_over_time <- combined_alpha_data %>%
  select(Community, Time, Colony, Apiary, Observed, Shannon, Simpson) %>%
  pivot_longer(
    cols = c(Observed, Shannon, Simpson),
    names_to = "Metric",
    values_to = "Value"
  )


# --- Plot for the Microbial Community by Colony ---
plot_microbial_colony <- tidy_alpha_over_time %>%
  filter(Community == "Microbial") %>%
  ggplot(aes(x = Time, y = Value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  facet_grid(Metric ~ Colony, scales = "free_y") +
  labs(
    title = "Microbial Alpha Diversity Over Time by Colony",
    x = "Time (Cumulative Days)",
    y = "Alpha Diversity Value"
  ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 8))

# --- Plot for the Viral Community by Colony ---
plot_viral_colony <- tidy_alpha_over_time %>%
  filter(Community == "Viral") %>%
  ggplot(aes(x = Time, y = Value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  facet_grid(Metric ~ Colony, scales = "free_y") +
  labs(
    title = "Viral Alpha Diversity Over Time by Colony",
    x = "Time (Cumulative Days)",
    y = "Alpha Diversity Value"
  ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 8))

# Print the plots
print(plot_microbial_colony)
print(plot_viral_colony)


# --- Plot for the Microbial Community by Apiary ---
plot_microbial_apiary <- tidy_alpha_over_time %>%
  filter(Community == "Microbial") %>%
  ggplot(aes(x = Time, y = Value)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue") +
  facet_grid(Metric ~ Apiary, scales = "free_y") +
  labs(
    title = "Microbial Alpha Diversity Over Time by Apiary",
    x = "Time (Cumulative Days)",
    y = "Alpha Diversity Value"
  ) +
  theme_bw()

# --- Plot for the Viral Community by Apiary ---
plot_viral_apiary <- tidy_alpha_over_time %>%
  filter(Community == "Viral") %>%
  ggplot(aes(x = Time, y = Value)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "firebrick") +
  facet_grid(Metric ~ Apiary, scales = "free_y") +
  labs(
    title = "Viral Alpha Diversity Over Time by Apiary",
    x = "Time (Cumulative Days)",
    y = "Alpha Diversity Value"
  ) +
  theme_bw()

# Print the plots
print(plot_microbial_apiary)
print(plot_viral_apiary)



# --- 2. Prepare the Data ---
# Define our two community types and filter the main dataframe
prevalence_base_data <- ecology_data %>%
  filter(type %in% c("vmag", "mmag")) %>%
  mutate(CommunityType = if_else(type == "vmag", "Viral", "Microbial"))

# --- 3. Calculate Prevalence for Each MAG (This part remains the same) ---
mag_prevalence <- prevalence_base_data %>%
  group_by(CommunityType, MAG_id) %>%
  summarise(num_samples = n_distinct(sample_name), .groups = 'drop')

# --- 4. Calculate the Cumulative Distribution (CORRECTED, EXPLICIT METHOD) ---
prevalence_distribution <- mag_prevalence %>%
  # Group by the desired variable
  group_by(CommunityType) %>%
  # Physically split the dataframe into a list of dataframes (one for Microbial, one for Viral)
  group_split() %>%
  # Use map_dfr to apply a function to each dataframe in the list and row-bind the results
  map_dfr(~{
    # Inside this block, '.x' refers to the dataframe for the current group
    
    current_community <- .x$CommunityType[1]
    counts <- .x$num_samples
    total_mags_in_group <- length(counts)
    
    k_samples <- 1:max(counts)
    num_mags_at_or_above_k <- sapply(k_samples, function(k) sum(counts >= k))
    
    proportion_mags <- num_mags_at_or_above_k / total_mags_in_group
    
    # Return a tidy dataframe for this one group
    tibble(
      CommunityType = current_community,
      k_samples = k_samples,
      proportion_mags = proportion_mags
    )
  })

# --- 5. Generate the Plot with Log10 Y-Axis (This part remains the same) ---
prevalence_plot_lm <- ggplot(
  prevalence_distribution,
  aes(x = k_samples, y = proportion_mags, color = CommunityType)
) +
  # Add the points for the actual calculated data
  geom_point(alpha = 0.5) +
  # Add a linear model best-fit line instead of connecting the dots
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  
  # The rest of the plot code remains the same
  scale_y_log10(
    breaks = c(1, 0.1, 0.01, 0.001, 0.0001),
    labels = percent
  ) +
  scale_x_continuous(
    breaks = seq(0, max(prevalence_distribution$k_samples) + 25, by = 25),
    limits = c(0, max(prevalence_distribution$k_samples) + 5)
  ) +
  annotation_logticks(sides = "l") +
  labs(
    title = "Prevalence of Microbial and Viral MAGs (Log Scale)",
    subtitle = "Trendline shows linear model fit (lm)",
    x = "Number of Samples (k)",
    y = "Proportion of MAGs Found in at Least 'k' Samples",
    color = "Community Type"
  ) +
  theme_bw() +
  theme(legend.position = "top")

print(prevalence_plot_lm)


ecology_data_timed <- add_cumulative_time_from_month_nucdiv(ecology_data)

accumulation_base_data <- ecology_data_timed %>%
  filter(type %in% c("mmag", "vmag")) %>%
  mutate(CommunityType = if_else(type == "vmag", "Viral", "Microbial")) %>%
  select(CommunityType, Colony, Time, MAG_id, sample_name)

# --- 3. Calculate Cumulative MAG Count (CORRECTED METHOD) ---
accumulation_curves_data <- accumulation_base_data %>%
  # Group by the factors we want to analyze separately
  group_by(CommunityType, Colony) %>%
  # Physically split the data into a list of tables
  group_split() %>%
  # Use purrr::map_dfr to run a function on each table and row-bind the results
  map_dfr(~{
    # Here, '.x' is correctly defined as the dataframe for the current group
    
    # Arrange the group's samples by time
    group_data_ordered <- .x %>% arrange(Time)
    
    # Get the unique time points in their chronological order
    unique_times <- unique(group_data_ordered$Time)
    
    # Calculate the cumulative count of unique MAGs at each time point
    cumulative_mags <- map_int(unique_times, function(current_time) {
      group_data_ordered %>%
        filter(Time <= current_time) %>%
        summarise(n_distinct(MAG_id)) %>%
        pull()
    })
    
    # Return a tidy dataframe with the results for this group
    tibble(
      CommunityType = .x$CommunityType[1],
      Colony = .x$Colony[1],
      Time = unique_times,
      cumulative_mags = cumulative_mags
    )
  })

# --- 4. Generate the Plot (This part remains the same) ---
accumulation_plot <- ggplot(
  accumulation_curves_data,
  aes(x = Time, y = cumulative_mags, color = CommunityType)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Colony, scales = "free") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Cumulative MAG Discovery Over Time by Colony",
    x = "Time (Cumulative Days)",
    y = "Cumulative Number of Unique MAGs Observed",
    color = "Community Type"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(accumulation_plot)




library(ggrepel) # For non-overlapping text labels

# --- 2. Filter data and calculate metrics for each MAG ---
mag_stability_metrics <- ecology_data %>%
  # Filter for microbial MAGs only
  filter(type == "mmag") %>%
  # Group by Genus and MAG to calculate metrics per MAG
  group_by(Genus, MAG_id) %>%
  summarise(
    # Count how many unique samples the MAG is in
    prevalence = n_distinct(sample_name),
    # Calculate mean and standard deviation of abundance (only where present)
    mean_abundance = mean(mean, na.rm = TRUE),
    sd_abundance = sd(mean, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Calculate the Coefficient of Variation (Volatility)
  # Replace NA sd (for MAGs in only 1 sample) with 0, as there's no variation
  mutate(
    sd_abundance = if_else(is.na(sd_abundance), 0, sd_abundance),
    volatility_cv = sd_abundance / mean_abundance
  )

# View the results
# Each row is a MAG with its stability metrics
head(mag_stability_metrics)


# --- 3. Prepare for plotting: Find the most speciose genera ---
# This helps create a cleaner, more focused plot
top_genera <- mag_stability_metrics %>%
  count(Genus, sort = TRUE) %>%
  slice_head(n = 9) %>% # Let's plot the top 9 genera with the most MAGs
  pull(Genus)

# Filter the metrics data to only these top genera
plot_data <- mag_stability_metrics %>%
  filter(Genus %in% top_genera)

# --- 4. Create the quadrant plot ---
stability_plot <- ggplot(plot_data, aes(x = prevalence, y = volatility_cv)) +
  geom_point(aes(color = Genus), alpha = 0.7, show.legend = FALSE) +
  # Use ggrepel to label some of the most dynamic/stable MAGs without overlap
  geom_text_repel(
    data = . %>% group_by(Genus) %>% filter(
      prevalence == max(prevalence) | volatility_cv == max(volatility_cv)
    ),
    aes(label = MAG_id),
    size = 2.5,
    max.overlaps = 5
  ) +
  # Facet by genus to create a separate panel for each
  facet_wrap(~ Genus, scales = "free") +
  # Add lines for the median to create quadrants
  geom_hline(data = . %>% group_by(Genus) %>% summarise(med = median(volatility_cv)),
             aes(yintercept = med), linetype = "dashed", color = "grey50") +
  geom_vline(data = . %>% group_by(Genus) %>% summarise(med = median(prevalence)),
             aes(xintercept = med), linetype = "dashed", color = "grey50") +
  scale_y_log10() + # Log scale helps spread out the volatility values
  annotation_logticks(sides = "l") +
  labs(
    title = "MAG Stability within Top Genera",
    subtitle = "Comparing Prevalence and Abundance Volatility",
    x = "Prevalence (Number of Samples)",
    y = "Volatility (Coefficient of Variation)"
  ) +
  theme_bw()

print(stability_plot)


# --- 1. Prepare Data ---
# Filter for only microbial worker data
microbial_data <- ecology_data %>%
  filter(type %in% c("mmag"), Caste == "Worker")

# Create a list to store the tidy decay data for each level
decay_data_by_level <- list()

# --- 2. Loop through the two taxonomic levels ---
for (taxonomic_level in c("MAG_id", "Genus")) {
  
  message(paste("--- Analyzing temporal decay at the", taxonomic_level, "level ---"))
  
  # A. Create the community matrix at the appropriate level
  community_matrix <- microbial_data %>%
    group_by(sample_name, !!sym(taxonomic_level)) %>% # Use !!sym() to evaluate the variable
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = !!sym(taxonomic_level), values_from = total_abundance, values_fill = 0) %>%
    tibble::column_to_rownames("sample_name")
  
  # B. Get metadata for this subset
  metadata <- microbial_data %>%
    select(sample_name, Colony, Time) %>%
    distinct() %>%
    filter(sample_name %in% rownames(community_matrix))
  
  # C. Calculate pairwise dissimilarities within each colony
  unique_colonies <- unique(metadata$Colony)
  colony_decay_data <- list()
  for (col in unique_colonies) {
    meta_subset <- metadata %>% filter(Colony == col)
    if (nrow(meta_subset) < 2) next
    
    comm_subset <- community_matrix[meta_subset$sample_name, ]
    bray_dist <- vegdist(comm_subset, method = "bray", na.rm = TRUE)
    time_lag_dist <- dist(meta_subset$Time)
    
    # Tidy the matrices
    bray_matrix <- as.matrix(bray_dist)
    bray_matrix[upper.tri(bray_matrix, diag = TRUE)] <- NA
    time_lag_matrix <- as.matrix(time_lag_dist)
    time_lag_matrix[upper.tri(time_lag_matrix, diag = TRUE)] <- NA
    
    decay_df <- data.frame(Dissimilarity = as.vector(bray_matrix), TimeLag = as.vector(time_lag_matrix)) %>% na.omit()
    colony_decay_data[[as.character(col)]] <- decay_df
  }
  
  # D. Store the combined tidy data for this taxonomic level
  decay_data_by_level[[taxonomic_level]] <- bind_rows(colony_decay_data)
}

# --- 3. Combine and Plot ---
# Bind the results from the Genus and MAG_id levels into one dataframe
comparison_data <- bind_rows(decay_data_by_level, .id = "TaxonomicLevel")

# Create the comparison plot
dissimilarity_level_plot <- ggplot(comparison_data, aes(x = TimeLag, y = Dissimilarity, color = TaxonomicLevel)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Rate of Community Change: Genus vs. MAG Level",
    subtitle = "Comparing microbial community turnover at different taxonomic resolutions",
    x = "Time Lag Between Samples (Days)",
    y = "Bray-Curtis Dissimilarity",
    color = "Taxonomic Level"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(dissimilarity_level_plot)


# --- 2. Prepare Data ---
microbial_data_timed <- ecology_data %>%
  filter(type == "mmag") %>%
  add_cumulative_time_from_month_nucdiv()

# --- 3. Calculate Dynamism Score for each Genus in each Colony ---
genus_colony_dynamism <- microbial_data_timed %>%
  group_by(Genus) %>%
  filter(n_distinct(MAG_id) >= 2) %>%
  ungroup() %>%
  group_by(Colony, Genus) %>%
  group_split() %>%
  map_dfr(~{
    if (n_distinct(.x$Time) < 2) return(NULL)
    
    # Create a time-series matrix (Time x MAGs)
    time_matrix <- .x %>%
      select(Time, MAG_id, mean) %>%
      
      # --- FIX: Summarize duplicates by summing them ---
      group_by(Time, MAG_id) %>%
      summarise(total_mean = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
      
      # Now pivot the clean, summarized data
      pivot_wider(names_from = MAG_id, values_from = total_mean, values_fill = 0) %>%
      arrange(Time)
    
    abund_matrix <- time_matrix %>% select(-Time)
    
    consecutive_dissimilarities <- c()
    for (i in 1:(nrow(abund_matrix) - 1)) {
      dissim <- vegdist(abund_matrix[i:(i+1), ], method = "bray")
      consecutive_dissimilarities <- c(consecutive_dissimilarities, dissim)
    }
    
    dynamism_score <- mean(consecutive_dissimilarities, na.rm = TRUE)
    
    tibble(
      Colony = .x$Colony[1],
      Genus = .x$Genus[1],
      Dynamism_Score = dynamism_score
    )
  })

# --- 4. View and Visualize the Results ---
# View the most dynamic genera
print(head(arrange(genus_colony_dynamism, desc(Dynamism_Score))))

# Visualize the distribution of scores
ggplot(genus_colony_dynamism, aes(x = Dynamism_Score)) +
  geom_histogram(bins = 30, fill = "darkred", alpha = 0.7) +
  geom_vline(xintercept = median(genus_colony_dynamism$Dynamism_Score), 
             linetype = "dashed", color = "black") +
  labs(
    title = "Distribution of Genus Dynamism Scores",
    x = "Genus Dynamism Score (Low = Stable, High = Dynamic)",
    y = "Frequency"
  ) +
  theme_bw()

# Visualize the distribution of scores
ggplot(genus_colony_dynamism, aes(x = Dynamism_Score)) +
  geom_histogram(bins = 30, fill = "darkred", alpha = 0.7) +
  geom_vline(xintercept = median(genus_colony_dynamism$Dynamism_Score), 
             linetype = "dashed", color = "black") +
  labs(
    title = "Distribution of Genus Dynamism Scores",
    x = "Genus Dynamism Score (Low = Stable, High = Dynamic)",
    y = "Frequency"
  ) +
  theme_bw()

ggplot(data=classified_genera,aes(x=Dynamism_Score,fill=Stability_Class))+
  geom_histogram()

median_dynamism <- median(genus_colony_dynamism$Dynamism_Score, na.rm = TRUE)

# Add a classification column based on the median
classified_genera <- genus_colony_dynamism %>%
  mutate(Stability_Class = if_else(Dynamism_Score > median_dynamism, "Dynamic", "Stable"))

# View the most dynamic genus-colony pairs
classified_genera %>%
  filter(Stability_Class == "Dynamic") %>%
  arrange(desc(Dynamism_Score))


bic_values <- mclustBIC(genus_colony_dynamism$Dynamism_Score)
summary(bic_values)
plot(bic_values)

# --- 3. Fit the Gaussian Mixture Model ---
# We tell the model to look for G=2 groups (Stable and Dynamic)
gmm_model <- Mclust(genus_colony_dynamism$Dynamism_Score, G = 2)
# --- 1. Ensure the necessary objects are loaded ---
# (Assuming 'gmm_model' and 'genus_colony_dynamism' are in your environment)

# --- 2. Re-create the classified data ---
# This ensures we start from the correct base object
genus_colony_dynamism$classification <- gmm_model$classification
cluster_means <- gmm_model$parameters$mean
stable_cluster_number <- which.min(cluster_means)
classified_genera_gmm <- genus_colony_dynamism %>%
  mutate(Stability_Class = if_else(
    classification == stable_cluster_number, "Stable", "Dynamic"
  ))

# --- 3. Run the Summary Calculation with the FIX ---

# Add a check to see the object's class before the operation
cat("The class of 'classified_genera_gmm' is:", class(classified_genera_gmm), "\n")

# THE FIX: Use bind_rows() to ensure the input is a data frame
genus_summary_table <- bind_rows(classified_genera_gmm) %>%
  count(Genus, Stability_Class) %>%
  tidyr::pivot_wider(names_from = Stability_Class, values_from = n, values_fill = 0)

# --- 4. Place the final table into a list ---
summary_list <- list(
  GenusStabilitySummary = genus_summary_table
)

# View the results with the new classification column
head(classified_genera_gmm)



# --- Create a visualization of the model fit ---
ggplot(classified_genera_gmm, aes(x = Dynamism_Score)) +
  # Plot the histogram of the original data
  geom_histogram(aes(y = ..density..), bins = 30, fill = "grey", alpha = 0.6) +
  
  # Overlay the two fitted Gaussian curves from the model
  stat_function(
    fun = function(x, mean, sd, lambda) { lambda * dnorm(x, mean, sd) },
    args = list(mean = gmm_model$parameters$mean[1], sd = sqrt(gmm_model$parameters$variance$sigmasq[1]), lambda = gmm_model$parameters$pro[1]),
    geom = "area", fill = "#1f78b4", alpha = 0.5 # Blue for one cluster
  ) +
  stat_function(
    fun = function(x, mean, sd, lambda) { lambda * dnorm(x, mean, sd) },
    args = list(mean = gmm_model$parameters$mean[2], sd = sqrt(gmm_model$parameters$variance$sigmasq[2]), lambda = gmm_model$parameters$pro[2]),
    geom = "area", fill = "#e31a1c", alpha = 0.5 # Red for the other
  ) +
  
  labs(
    title = "GMM Classification of Genus Dynamism",
    subtitle = "Data distribution overlaid with Stable and Dynamic clusters",
    x = "Genus Dynamism Score",
    y = "Density"
  ) +
  theme_bw()


classified_genera_gmm %>%
  filter(Stability_Class == "Dynamic") %>%
  arrange(desc(Dynamism_Score))

# Count how many times each Genus was classified as Stable vs. Dynamic
counts <- classified_genera_gmm %>%
  dplyr::count(Genus, Stability_Class) %>%
  tidyr::pivot_wider(names_from = Stability_Class, values_from = n, values_fill = 0)





# --- 3. Perform the Shapiro-Wilk Test for each group ---
normality_test_results <- classified_genera_gmm %>%
  group_by(Stability_Class) %>%
  summarise(
    shapiro_statistic = shapiro.test(Dynamism_Score)$statistic,
    p_value = shapiro.test(Dynamism_Score)$p.value
  )

cat("--- Shapiro-Wilk Normality Test Results ---\n")
print(normality_test_results)

# --- 4. Create Q-Q Plots for each group ---
qq_plots <- ggplot(classified_genera_gmm, aes(sample = Dynamism_Score)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Stability_Class, scales = "free") +
  labs(
    title = "Q-Q Plots for Stability Classes",
    subtitle = "Points should fall on the line if data is normally distributed",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw()

print(qq_plots)

# install.packages("ggridges")
library(ggridges)
library(ggplot2)
library(dplyr)

# We use the 'classified_genera_gmm' dataframe
# Order the Genus factor by the median dynamism score for a clean look
dynamism_ridges <- ggplot(
  classified_genera_gmm,
  aes(x = Dynamism_Score, y = reorder(Genus, Dynamism_Score, median), fill = Stability_Class)
) +
  # Create the ridgeline densities
  geom_density_ridges(
    aes(point_color = Stability_Class, point_fill = Stability_Class),
    alpha = 0.7,
    jittered_points = TRUE, # Show the individual data points
    point_size = 1
  ) +
  scale_fill_manual(values = c("Stable" = "#1f78b4", "Dynamic" = "#e31a1c"), name = "GMM Class")+
  labs(
    title = "Distribution of Dynamism Scores by Genus",
    subtitle = "Each ridge shows the stability profile of a genus across all colonies",
    x = "Genus Dynamism Score (Higher = More Dynamic)",
    y = "Genus"
  ) +
  theme_ridges()

print(dynamism_ridges)


# We use the 'classified_genera_gmm' dataframe
# Order the Genus factor by the median dynamism score for a clean look
ridges_with_points <- ggplot(
  classified_genera_gmm,
  aes(x = Dynamism_Score, y = reorder(Genus, Dynamism_Score, median))
) +
  geom_density_ridges(
    alpha = 0.5, fill = "lightgrey",
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    # We only need to map point_color, as the shape "|" has no fill
    aes(point_color = Stability_Class),
    point_shape = "|", point_size = 3, point_alpha = 0.8
  ) +
  # THE FIX: Use the correct ggplot2 function name
  scale_color_manual(
    values = c("Stable" = "#1f78b4", "Dynamic" = "#e31a1c"), 
    name = "GMM Class"
  ) +
  labs(
    title = "Distribution of Dynamism Scores by Genus",
    subtitle = "Points under each ridge show the classification for each colony",
    x = "Genus Dynamism Score (Higher = More Dynamic)",
    y = "Genus"
  ) +
  theme_ridges()

print(ridges_with_points)















# --- 2. Prepare Data ---
microbial_data_timed <- ecology_data %>%
  filter(type == "mmag") %>%
  add_cumulative_time_from_month_nucdiv()

# --- 3. Flexible Persistence Filter ---

# a. Count total months sampled per colony
colony_month_counts <- microbial_data_timed %>%
  group_by(Colony) %>%
  summarise(total_months_sampled = n_distinct(Month), .groups = 'drop')

# b. Calculate persistence percentage for each genus in each colony
genus_persistence <- microbial_data_timed %>%
  group_by(Colony, Genus) %>%
  summarise(genus_months_observed = n_distinct(Month), .groups = 'drop') %>%
  inner_join(colony_month_counts, by = "Colony") %>%
  mutate(persistence_pct = genus_months_observed / total_months_sampled)

# c. Create the two filtered datasets based on the thresholds
persistent_data_25pct <- semi_join(
  microbial_data_timed, 
  filter(genus_persistence, persistence_pct >= 0.25), 
  by = c("Colony", "Genus")
)
persistent_data_50pct <- semi_join(
  microbial_data_timed, 
  filter(genus_persistence, persistence_pct >= 0.50), 
  by = c("Colony", "Genus")
)

cat("Original rows:", nrow(microbial_data_timed), "\n")
cat("Rows kept at >= 25% persistence:", nrow(persistent_data_25pct), "\n")
cat("Rows kept at >= 50% persistence:", nrow(persistent_data_50pct), "\n")


# --- 4. Create an Analysis Function ---
# This function contains the code for our dynamism analysis.
# By making it a function, we can easily run it on different datasets without copying code.
calculate_dynamism <- function(input_data) {
  
  dynamism_results <- input_data %>%
    group_by(Genus) %>%
    filter(n_distinct(MAG_id) >= 2) %>%
    ungroup() %>%
    group_by(Colony, Genus) %>%
    group_split() %>%
    map_dfr(~{
      if (n_distinct(.x$Time) < 2) return(NULL)
      time_matrix <- .x %>%
        select(Time, MAG_id, mean) %>%
        group_by(Time, MAG_id) %>%
        summarise(total_mean = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
        pivot_wider(names_from = MAG_id, values_from = total_mean, values_fill = 0) %>%
        arrange(Time)
      abund_matrix <- time_matrix %>% select(-Time)
      consecutive_dissimilarities <- c()
      for (i in 1:(nrow(abund_matrix) - 1)) {
        dissim <- vegdist(abund_matrix[i:(i+1), ], method = "bray")
        consecutive_dissimilarities <- c(consecutive_dissimilarities, dissim)
      }
      dynamism_score <- mean(consecutive_dissimilarities, na.rm = TRUE)
      tibble(
        Colony = .x$Colony[1],
        Genus = .x$Genus[1],
        Dynamism_Score = dynamism_score
      )
    })
  
  return(dynamism_results)
}

# --- 5. Run the Analysis on Both Datasets ---
# Run the analysis on the 25% persistence dataset
dynamism_results_25pct <- calculate_dynamism(persistent_data_25pct)

# Run the analysis on the 50% persistence dataset
dynamism_results_50pct <- calculate_dynamism(persistent_data_50pct)

cat("\n--- Dynamism Score Results (>= 25% Persistence) ---\n")
print(head(arrange(dynamism_results_25pct, desc(Dynamism_Score))))

cat("\n--- Dynamism Score Results (>= 50% Persistence) ---\n")
print(head(arrange(dynamism_results_50pct, desc(Dynamism_Score))))


# --- 3. Define the Complete Analysis Pipeline as a Function ---
run_gmm_and_normality_pipeline <- function(input_df, title_prefix = "") {
  
  cat(paste("\n\n--- Running Analysis for:", title_prefix, "---\n"))
  
  # --- a. Fit the Gaussian Mixture Model ---
  gmm_model <- Mclust(input_df$Dynamism_Score, G = 2)
  
  # --- b. Check BIC values for model validation ---
  cat("\n--- BIC Model Selection ---\n")
  print(summary(mclustBIC(input_df$Dynamism_Score)))
  
  # --- c. Classify the data based on the model ---
  input_df$classification <- gmm_model$classification
  cluster_means <- gmm_model$parameters$mean
  stable_cluster_number <- which.min(cluster_means)
  classified_data <- input_df %>%
    mutate(Stability_Class = if_else(
      classification == stable_cluster_number, "Stable", "Dynamic"
    ))
  
  # --- d. Perform Shapiro-Wilk Normality Tests on the results ---
  normality_test_results <- classified_data %>%
    group_by(Stability_Class) %>%
    summarise(
      shapiro_statistic = shapiro.test(Dynamism_Score)$statistic,
      p_value = shapiro.test(Dynamism_Score)$p.value,
      .groups = 'drop'
    )
  cat("\n--- Shapiro-Wilk Normality Test Results ---\n")
  print(normality_test_results)
  
  # --- e. Generate Q-Q Plots ---
  qq_plot <- ggplot(classified_data, aes(sample = Dynamism_Score)) +
    stat_qq() +
    stat_qq_line() +
    facet_wrap(~ Stability_Class, scales = "free") +
    labs(
      title = paste(title_prefix, "Q-Q Plots for Stability Classes")
    ) +
    theme_bw()
  
  # --- f. Generate GMM Visualization ---
  gmm_plot <- ggplot(classified_data, aes(x = Dynamism_Score)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "grey", alpha = 0.6) +
    stat_function(
      fun = function(x, m, s, l) { l * dnorm(x, m, s) },
      args = list(m = gmm_model$parameters$mean[1], s = sqrt(gmm_model$parameters$variance$sigmasq[1]), l = gmm_model$parameters$pro[1]),
      geom = "area", fill = "#1f78b4", alpha = 0.5
    ) +
    stat_function(
      fun = function(x, m, s, l) { l * dnorm(x, m, s) },
      args = list(m = gmm_model$parameters$mean[2], s = sqrt(gmm_model$parameters$variance$sigmasq[2]), l = gmm_model$parameters$pro[2]),
      geom = "area", fill = "#e31a1c", alpha = 0.5
    ) +
    labs(
      title = paste(title_prefix, "GMM Classification of Genus Dynamism"),
      x = "Genus Dynamism Score", y = "Density"
    ) +
    theme_bw()
  
  # --- g. Print plots and return a list of results ---
  print(gmm_plot)
  print(qq_plot)
  
  return(list(
    classified_data = classified_data,
    gmm_plot = gmm_plot,
    qq_plot = qq_plot,
    normality_test = normality_test_results
  ))
}

# --- 4. Run the Pipeline on Both Datasets ---
# Run for the 25% persistence dataset
results_25pct <- run_gmm_and_normality_pipeline(
  dynamism_results_25pct, 
  title_prefix = ">= 25% Persistence"
)

# Run for the 50% persistence dataset
results_50pct <- run_gmm_and_normality_pipeline(
  dynamism_results_50pct, 
  title_prefix = ">= 50% Persistence"
)


# --- 2. Use the Classified Data from the 25% Persistence Analysis ---
# The analysis function returned a list, so we extract the data frame
classified_data_25pct <- results_25pct$classified_data
classified_data_50pct <- results_50pct$classified_data

# --- 3. Create the Ridgeline Plot ---
# We order the Genus factor by the median dynamism score for a clean look
ridges_plot_25pct <- ggplot(
  classified_data_25pct,
  aes(x = Dynamism_Score, y = reorder(Genus, Dynamism_Score, median))
) +
  # Create the density ridges with a neutral fill
  geom_density_ridges(
    alpha = 0.5, fill = "lightgrey",
    # Add the individual data points (each is a colony)
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    # Color the points by their GMM classification
    aes(point_color = Stability_Class),
    point_shape = "|", point_size = 3, point_alpha = 0.8
  ) +
  # Define the custom colors for the points
  scale_color_manual(
    values = c("Stable" = "#1f78b4", "Dynamic" = "#e31a1c"), 
    name = "GMM Class"
  ) +
  labs(
    title = "Distribution of Dynamism Scores by Genus (>= 25% Persistence)",
    subtitle = "Points under each ridge show the classification for each colony",
    x = "Genus Dynamism Score (Higher = More Dynamic)",
    y = "Genus"
  ) +
  theme_ridges()

print(ridges_plot_25pct)


dynamism_ridges_25 <- ggplot(
  classified_data_25pct,
  aes(x = Dynamism_Score, y = reorder(Genus, Dynamism_Score, median), fill = Stability_Class)
) +
  # Create the ridgeline densities
  geom_density_ridges(
    aes(point_color = Stability_Class, point_fill = Stability_Class),
    alpha = 0.7,
    jittered_points = TRUE, # Show the individual data points
    point_size = 1
  ) +
  scale_fill_manual(values = c("Stable" = "#1f78b4", "Dynamic" = "#e31a1c"), name = "GMM Class")+
  labs(
    title = "Distribution of Dynamism Scores by Genus",
    subtitle = "Each ridge shows the stability profile of a genus across all colonies",
    x = "Genus Dynamism Score (Higher = More Dynamic)",
    y = "Genus"
  ) +
  theme_ridges()+
  geom_vline(xintercept=1,color="black")

dynamism_ridges_50 <- ggplot(
  classified_data_50pct,
  aes(x = Dynamism_Score, y = reorder(Genus, Dynamism_Score, median), fill = Stability_Class)
) +
  # Create the ridgeline densities
  geom_density_ridges(
    aes(point_color = Stability_Class, point_fill = Stability_Class),
    alpha = 0.7,
    jittered_points = TRUE, # Show the individual data points
    point_size = 1
  ) +
  scale_fill_manual(values = c("Stable" = "#1f78b4", "Dynamic" = "#e31a1c"), name = "GMM Class")+
  labs(
    title = "Distribution of Dynamism Scores by Genus",
    subtitle = "Each ridge shows the stability profile of a genus across all colonies",
    x = "Genus Dynamism Score (Higher = More Dynamic)",
    y = "Genus"
  ) +
  theme_ridges()
print(dynamism_ridges)
print(dynamism_ridges_25)
print(dynamism_ridges_50)


# --- 2. Prepare the Data ---
# Use the classified data from our 25% persistence analysis
classified_data_25pct <- results_25pct$classified_data

# Create a contingency table of observed counts: Genus vs. Stability Class
contingency_table_genus <- table(classified_data_25pct$Genus, classified_data_25pct$Stability_Class)
contingency_table_genus
 # Use the same contingency table from before
fisher_result_genus <- fisher.test(contingency_table_genus, simulate.p.value = TRUE)

cat("\n--- Fisher's Exact Test Results (with simulation) ---\n")
print(fisher_result_genus)

# --- 3. Extract and View the Standardized Residuals ---
# Standardized residuals tell us how many standard deviations each cell's
# observed count is from its expected count.
std_residuals <- round(chi_sq_result_genus$stdres, 2)

cat("--- Standardized Residuals ---\n")
print(std_residuals)

# --- 4. Visualize the Residuals ---
# This plot makes it easy to spot the significant cells.
corrplot(std_residuals, is.cor = FALSE, 
         method = "color", tl.col = "black", tl.srt = 45,
         title = "Standardized Residuals (Genus vs. Stability)", 
         mar = c(1,1,2,1))


# --- 1. Define the Genera to Exclude ---
genera_to_exclude <- c(
  "hafnia", "enterococcus", "bombella", "pelethocola", "escherichia", 
  "pantoea", "mobilisporobacter", "providencia", "fructobacillus", 
  "no match", "spiroplasma"
)

# --- 2. Filter Your Dataframe ---
# (Assuming 'classified_data_50pct' is your dataframe)
filtered_plot_data_25 <- classified_data_25pct %>%
  filter(!Genus %in% genera_to_exclude)

# You would do the same for your other dataframes, for example:
# filtered_plot_data_25 <- classified_data_25pct %>%
#   filter(!Genus %in% genera_to_exclude)


# --- 3. Re-create the Plot with Filtered Data ---
dynamism_ridges_25_filtered <- ggplot(
  filtered_plot_data_25,  # Use the new filtered dataframe
  aes(x = Dynamism_Score, y = reorder(Genus, Dynamism_Score, median), fill = Stability_Class)
) +
  geom_density_ridges(
    aes(point_color = Stability_Class, point_fill = Stability_Class),
    alpha = 0.7,
    jittered_points = TRUE,
    point_size = 1
  ) +
  scale_fill_manual(values = c("Stable" = "#1f78b4", "Dynamic" = "#e31a1c"), name = "GMM Class") +
  labs(
    title = "Distribution of Dynamism Scores by Genus (Filtered)",
    subtitle = "Each ridge shows the stability profile of a genus across all colonies",
    x = "Genus Dynamism Score (Higher = More Dynamic)",
    y = "Genus"
  ) +
  theme_ridges()

# Print the new, filtered plot
print(dynamism_ridges_25_filtered)





# Create the filtered dataframe
filtered_analysis_data <- classified_data_25pct %>%
  filter(!Genus %in% genera_to_exclude)

# --- 3. Create a New Contingency Table ---
contingency_table_filtered <- table(filtered_analysis_data$Genus, filtered_analysis_data$Stability_Class)

cat("--- Contingency Table (Filtered) ---\n")
print(contingency_table_filtered)

# --- 4. Run Fisher's Exact Test on Filtered Data ---
fisher_result_filtered <- fisher.test(contingency_table_filtered, simulate.p.value = TRUE)

cat("\n--- Fisher's Exact Test Results (Filtered) ---\n")
print(fisher_result_filtered)

# --- 5. Recreate Chi-Squared Results to Get New Residuals ---
# We run this not for the p-value, but to get the residuals.
chi_sq_result_filtered <- chisq.test(contingency_table_filtered)

# Extract and print the new standardized residuals
std_residuals_filtered <- round(chi_sq_result_filtered$stdres, 2)

cat("\n--- Standardized Residuals (Filtered) ---\n")
print(std_residuals_filtered)

# --- 6. Visualize the New Residuals ---
corrplot(std_residuals_filtered, is.cor = FALSE, 
         method = "color", tl.col = "black", tl.srt = 45,
         title = "Standardized Residuals (Filtered)", mar = c(1,1,2,1))



# --- 1. Sort the residuals matrix ---
# We'll sort the rows based on the value in the "Dynamic" column, in descending order.
sorted_residuals <- std_residuals_filtered[
  order(std_residuals_filtered[, "Dynamic"], decreasing = TRUE), 
]

# --- 2. Re-plot with the sorted data ---
# We add 'order = "original"' to tell corrplot to respect our new sorted order
# and 'cl.pos = "n"' to remove the color legend, which is redundant here.
corrplot(sorted_residuals, 
         is.cor = FALSE, 
         method = "color",
         order = "original",  # <-- The key change
         tl.col = "black", 
         tl.srt = 45,
         cl.pos = "n", # Hide the color legend for a cleaner look
         title = "Standardized Residuals (Sorted by Dynamic Predisposition)", 
         mar = c(1,1,2,1))


# a. Define your custom sorting order
colony_order <- c(100, 101, 102, 103, 104, 105, 106, 107, 555, 777, 888, 999)
month_order <- c("May", "June", "July", "August", "September", "October", "November", "December", "January","February")
my_colors <- colorRampPalette(c("#FF975C", "white", "#2A3D45"))(100)

# --- 2. Aggregate Replicates by Colony and Month ---
# (Assuming 'mmag_data' is loaded)
pooled_data <- mmag_data %>%
  group_by(Colony, Month, MAG_id) %>%
  # For each MAG, find its average abundance across replicates for that month
  summarise(mean_abundance = mean(mean, na.rm = TRUE), .groups = 'drop') %>%
  # Create a new, unique identifier for each Colony-Month combination
  mutate(colony_month = paste(Colony, Month, sep = "_"))

# --- 3. Create the New Community Matrix ---
pooled_community_matrix <- pooled_data %>%
  pivot_wider(
    id_cols = colony_month,
    names_from = MAG_id,
    values_from = mean_abundance,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("colony_month")

# --- 4. Calculate New Dissimilarity Matrix ---
pooled_bray_dist <- vegdist(pooled_community_matrix, method = "bray")
pooled_bray_matrix <- as.matrix(pooled_bray_dist)

# --- 5. Generate the New Heatmap ---
# This time, we LET pheatmap cluster the rows and columns to find patterns
pheatmap(
  pooled_bray_matrix,
  color = my_colors,
  cluster_rows = FALSE,      
  cluster_cols = FALSE,      
  show_rownames = FALSE,     
  show_colnames = FALSE, 
  annotation_col = annotation_df,
  main = "Bray-Curtis Dissimilarity Between Pooled Colony-Month Samples"
)

pheatmap(
  bray_matrix_ordered,
  cluster_rows = FALSE,      
  cluster_cols = FALSE,      
  show_rownames = FALSE,     
  show_colnames = FALSE,     
  annotation_col = annotation_df,
  color = my_colors, # <-- The only new argument
  main = "Bray-Curtis Dissimilarity (Custom Colors)"
)

