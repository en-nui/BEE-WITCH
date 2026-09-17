library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
ecology_data <- read.csv('ecology_metadata_df.csv')
ecology_data <- ecology_data[ecology_data$Caste == "Worker",]

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





















# --- 0. Setup: Define inputs ---
# (Assuming 'mmag_alpha_df' is already created from your previous step)
# Let's verify sample sizes first
sample_counts <- table(mmag_alpha_df$Season)
print("Original Sample Sizes per Season:")
print(sample_counts)

min_size <- min(sample_counts)
print(paste("Downsampling all groups to N =", min_size))

# --- 1. Define the Permutation Function ---
run_balanced_kruskal <- function(df, metric_col, group_col, n_perms = 1000) {
  
  # Get the minimum sample size across groups
  group_counts <- table(df[[group_col]])
  min_n <- min(group_counts)
  
  # Storage for results
  p_values <- numeric(n_perms)
  h_statistics <- numeric(n_perms)
  
  set.seed(123) # For reproducibility
  
  for (i in 1:n_perms) {
    # Stratified subsampling: 
    # Take 'min_n' random rows from EACH group
    balanced_subset <- df %>%
      group_by(.data[[group_col]]) %>%
      sample_n(min_n) %>%
      ungroup()
    
    # Run Kruskal-Wallis on this balanced subset
    # Formula: Metric ~ Group
    f <- as.formula(paste(metric_col, "~", group_col))
    kruskal <- kruskal.test(f, data = balanced_subset)
    
    p_values[i] <- kruskal$p.value
    h_statistics[i] <- kruskal$statistic
  }
  
  return(data.frame(P_Value = p_values, H_Stat = h_statistics))
}

# --- 2. Run the Analysis for Shannon and Observed ---

# Run for Shannon Diversity
shannon_perms <- run_balanced_kruskal(mmag_alpha_df, "Shannon", "Season", n_perms = 1000)

# Run for Observed Richness
observed_perms <- run_balanced_kruskal(mmag_alpha_df, "Observed", "Season", n_perms = 1000)

# --- 3. Summarize and Visualize Results ---

summarize_perms <- function(perm_df, metric_name) {
  median_p <- median(perm_df$P_Value)
  percent_sig <- mean(perm_df$P_Value < 0.05) * 100
  
  cat(paste0("\n--- Results for ", metric_name, " (1000 Permutations) ---\n"))
  cat(paste("Median P-value:", format(median_p, digits=4), "\n"))
  cat(paste("Significant Iterations:", percent_sig, "%\n"))
  
  return(c(Median_P = median_p, Sig_Pct = percent_sig))
}

shannon_summary <- summarize_perms(shannon_perms, "Shannon Diversity")
observed_summary <- summarize_perms(observed_perms, "Observed Richness")

# --- 4. Plot the P-value Distribution ---
# This visualizes how robust your result is. 
# If the histogram is smashed against 0 (left side), the result is robustly significant.
# If it's spread out, the sample imbalance was likely driving the original result.

p_val_plot_data <- rbind(
  data.frame(Metric = "Shannon", P_Value = shannon_perms$P_Value),
  data.frame(Metric = "Observed", P_Value = observed_perms$P_Value)
)

ggplot(p_val_plot_data, aes(x = P_Value, fill = Metric)) +
  geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~Metric) +
  labs(
    title = "Robustness of Seasonality Effect",
    subtitle = paste("Distribution of P-values across 1000 balanced subsamples (N =", min_size, "per group)"),
    x = "Kruskal-Wallis P-value",
    y = "Count of Iterations"
  ) +
  theme_bw()























# --- 1. Define the Bootstrap Function ---
run_bootstrapped_means <- function(df, metric_col, group_col, n_perms = 1000) {
  
  # Get minimum sample size
  min_n <- min(table(df[[group_col]]))
  
  # Create a list to store the means from each iteration
  results_list <- list()
  
  set.seed(123)
  
  for (i in 1:n_perms) {
    # Subsample to balance groups
    balanced_subset <- df %>%
      group_by(.data[[group_col]]) %>%
      sample_n(min_n) %>%
      ungroup()
    
    # Calculate the MEAN metric for each group in this subset
    # (You can switch to median() if your data is very skewed)
    iteration_means <- balanced_subset %>%
      group_by(.data[[group_col]]) %>%
      summarise(Mean_Value = mean(.data[[metric_col]]), .groups = 'drop') %>%
      mutate(Iteration = i)
    
    results_list[[i]] <- iteration_means
  }
  
  # Combine all iterations into one big dataframe
  final_df <- bind_rows(results_list)
  return(final_df)
}

# --- 2. Run the Bootstrap ---
# (Assuming 'mmag_alpha_df' is your data)
boot_data <- run_bootstrapped_means(mmag_alpha_df, "Shannon", "Season", n_perms = 1000)

# Preview the data (It will look like: Season | Mean_Value | Iteration)
head(boot_data)





ggplot(boot_data, aes(x = Mean_Value, fill = Season)) +
  # alpha = 0.6 makes it transparent so you can see the overlap
  geom_density(alpha = 0.6, color = "black", size = 0.3) + 
  scale_fill_brewer(palette = "Set1") + # Nice distinct colors
  labs(
    title = "Bootstrapped Seasonal Diversity",
    subtitle = paste("Distribution of Mean Shannon Diversity across 1000 balanced subsamples"),
    x = "Mean Shannon Index",
    y = "Density (Frequency of Outcome)"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    legend.position = "top"
  )


ggplot(boot_data, aes(x = Mean_Value, y = Season, fill = Season)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) + # scale > 1 makes them overlap slightly
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Robustness of Seasonal Differences",
    x = "Mean Shannon Index (Bootstrapped)",
    y = "Season"
  ) +
  theme_bw() +
  theme(legend.position = "none") # Legend is redundant here
















# --- 1. Define the Bootstrap Function (With Replacement) ---
run_bootstrapped_means <- function(df, metric_col, group_col, n_perms = 1000) {
  
  # Determine the size to subsample to.
  # We use the minimum group size to keep it fair.
  min_n <- min(table(df[[group_col]]))
  
  # Optional: If min_n is tiny (e.g., < 5), you might want to hardcode a size
  # based on your knowledge, but min_n is usually the safest default.
  
  print(paste("Bootstrapping with N =", min_n, "per group (with replacement)."))
  
  results_list <- list()
  set.seed(123)
  
  for (i in 1:n_perms) {
    # 1. Stratified Bootstrap:
    # Group by Season -> Sample N rows WITH replacement
    boot_subset <- df %>%
      group_by(.data[[group_col]]) %>%
      sample_n(size = min_n, replace = TRUE) %>% # <--- CRITICAL FIX
      ungroup()
    
    # 2. Calculate the MEAN metric for each group in this bootstrap
    iteration_means <- boot_subset %>%
      group_by(.data[[group_col]]) %>%
      summarise(Mean_Value = mean(.data[[metric_col]]), .groups = 'drop') %>%
      mutate(Iteration = i)
    
    results_list[[i]] <- iteration_means
  }
  
  return(bind_rows(results_list))
}

# --- 2. Run Analysis ---
# (Assuming 'mmag_alpha_df' is loaded)
boot_data <- run_bootstrapped_means(mmag_alpha_df, "Shannon", "Season", n_perms = 1000)

# --- 3. Plot with ggridges (The Clean Look) ---
ggplot(boot_data, aes(x = Mean_Value, y = Season, fill = Season)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = 0.7, scale = 1.2) + 
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Bootstrapped Seasonal Diversity",
    subtitle = "Distribution of Mean Shannon Diversity (1000 Bootstraps)",
    x = "Mean Shannon Index",
    y = NULL # Season labels are already on the axis
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold")
  )



























# --- 1. Setup Data & Order ---
# Assuming 'mmag_alpha_df' is your source data
# We MUST set the factor levels so the plot is chronological, not alphabetical
mmag_alpha_df$Season <- factor(mmag_alpha_df$Season, 
                               levels = c("Spring", "Summer", "Fall", "Winter"))

# --- 2. The Bootstrap Function (With Replacement) ---
run_bootstrapped_means <- function(df, metric_col, group_col, n_perms = 5000) {
  
  # We use the smallest sample size (Winter = 18) to standardize variance
  min_n <- min(table(df[[group_col]]))
  #min_n <- 5
  
  results_list <- list()
  set.seed(123) # Reproducibility
  
  for (i in 1:n_perms) {
    # Stratified Bootstrap WITH Replacement
    boot_subset <- df %>%
      group_by(.data[[group_col]]) %>%
      sample_n(size = min_n, replace = TRUE) %>% 
      ungroup()
    
    # Calculate Mean
    iteration_means <- boot_subset %>%
      group_by(.data[[group_col]]) %>%
      summarise(Mean_Value = mean(.data[[metric_col]]), .groups = 'drop') %>%
      mutate(Iteration = i)
    
    results_list[[i]] <- iteration_means
  }
  
  return(bind_rows(results_list))
}

# Run the bootstrap
boot_data <- run_bootstrapped_means(mmag_alpha_df, "Shannon", "Season", n_perms = 5000)

# --- 3. The Visualization (Ridgeline Plot) ---
# We use Spring at the bottom and Winter at the top (or vice versa) to show time.
# "Spring -> Winter" logic usually implies reading Top-to-Bottom or Bottom-to-Top.
# Let's put Spring at the TOP so it reads like a timeline going down? 
# Actually, standard plots usually put the first factor at the bottom. 
# Let's stick to the Factor Order defined above.

ridges_plot <- ggplot(boot_data, aes(x = Mean_Value, y = Season, fill = Season)) +
  # density_ridges creates the "Mountain" look
  geom_density_ridges(
    scale = 1.5,             # How much they overlap (1.5 is good for comparison)
    quantile_lines = TRUE,   # Adds the median line
    quantiles = 2,           # Draws a line at the 50th percentile (median)
    alpha = 0.8,             # Transparency
    color = "white",         # White border looks cleaner
    size = 0.5
  ) +
  scale_fill_manual(values = c(
    "Spring" = "#66C2A5", 
    "Summer" = "#FC8D62", 
    "Fall"   = "#8DA0CB", 
    "Winter" = "#E78AC3"
  )) +
  labs(
    title = "Seasonal Progression of Alpha Diversity",
    subtitle = paste("Distribution of Mean Shannon Index (Bootstrapped N =", 18, "samples)"),
    x = "Mean Shannon Diversity",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title.x = element_text(size = 12),
    panel.grid.major.y = element_blank() # Removes horizontal grid lines for cleaner look
  )

print(ridges_plot)
# --- 4. Interpretation Stats: "Probability of Direction" ---
# This calculates: "In what % of simulations was Fall > Spring?"

# Pivot to wide format to compare columns easily
wide_boot <- boot_data %>%
  pivot_wider(names_from = Season, values_from = Mean_Value)

# Calculate probabilities
prob_fall_gt_spring <- mean(wide_boot$Fall > wide_boot$Spring)
prob_fall_gt_summer <- mean(wide_boot$Fall > wide_boot$Summer)
prob_summer_gt_spring <- mean(wide_boot$Summer > wide_boot$Spring)
prob_fall_gt_winter <- mean(wide_boot$Fall > wide_boot$Winter)

cat("\n--- Probability of Direction (Support for Hypothesis) ---\n")
cat(paste("Probability Fall > Spring: ", prob_fall_gt_spring * 100, "%\n", sep=""))
cat(paste("Probability Summer > Spring: ", prob_summer_gt_spring * 100, "%\n", sep=""))
cat(paste("Probability Fall > Winter:   ", prob_fall_gt_winter * 100, "%\n", sep=""))




























# --- 1. Prepare the Data ---
# Ensure the Season factor is ordered chronologically for both datasets
order_seasons <- c("Spring", "Summer", "Fall", "Winter")
boot_data$Season <- factor(boot_data$Season, levels = order_seasons)
mmag_alpha_df$Season <- factor(mmag_alpha_df$Season, levels = order_seasons)

# --- 2. Calculate Raw Summary Statistics (The "Observed") ---
raw_stats <- mmag_alpha_df %>%
  group_by(Season) %>%
  summarise(
    Raw_Mean = mean(Shannon),
    Raw_SE = sd(Shannon) / sqrt(n()), # Standard Error
    .groups = 'drop'
  )

# --- 3. Create the Composite Plot ---
validation_plot <- ggplot() +
  
  # LAYER 1: The "Robust Prediction" (Bootstrapped Ridges)
  # This shows the range of likely means given equal sampling
  geom_density_ridges(
    data = boot_data, 
    aes(x = Mean_Value, y = Season, fill = Season),
    scale = 1.2, 
    alpha = 0.6, 
    color = NA # No border for the background ridges
  ) +
  
  # LAYER 2: The "Observed Reality" (Raw Mean + Error Bars)
  # We overlay the raw data summary on top
  geom_pointrange(
    data = raw_stats,
    aes(x = Raw_Mean, y = Season, xmin = Raw_Mean - Raw_SE, xmax = Raw_Mean + Raw_SE),
    color = "black",
    size = 0.8,    # Size of the line
    fatten = 4     # Size of the dot
  ) +
  
  # Formatting
  scale_fill_manual(values = c(
    "Spring" = "#66C2A5", 
    "Summer" = "#FC8D62", 
    "Fall"   = "#8DA0CB", 
    "Winter" = "#E78AC3"
  )) +
  labs(
    title = "Model Validation: Robust vs. Observed Diversity",
    subtitle = "Colored Ridges = Bootstrapped Means (N=18) | Black Dots = Raw Observed Means (+/- SE)",
    x = "Shannon Diversity Index",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    panel.grid.major.y = element_blank()
  )

print(validation_plot)





















# --- 1. The Contrast Function ---
# This calculates the difference between specific pairs of seasons
run_bootstrap_contrasts <- function(df, metric_col, group_col, n_perms = 5000) {
  
  min_n <- min(table(df[[group_col]]))
  contrast_list <- list()
  set.seed(123)
  
  for (i in 1:n_perms) {
    # Resample
    boot_subset <- df %>%
      group_by(.data[[group_col]]) %>%
      sample_n(size = min_n, replace = TRUE) %>%
      ungroup()
    
    # Calculate Means
    means <- boot_subset %>%
      group_by(.data[[group_col]]) %>%
      summarise(Mean = mean(.data[[metric_col]]), .groups = 'drop')
    
    # Pull values for calculation
    spring <- means$Mean[means[[group_col]] == "Spring"]
    summer <- means$Mean[means[[group_col]] == "Summer"]
    fall   <- means$Mean[means[[group_col]] == "Fall"]
    winter <- means$Mean[means[[group_col]] == "Winter"]
    
    # Calculate Deltas (Change from previous season)
    # You can customize these pairs based on what you want to show
    contrast_list[[i]] <- data.frame(
      Iteration = i,
      `Summer vs Spring` = summer - spring,
      `Fall vs Summer`   = fall - summer,
      `Winter vs Fall`   = winter - fall,
      `Winter vs Spring` = winter - spring # Optional: Full loop check
    )
  }
  
  return(bind_rows(contrast_list))
}

# --- 2. Run Analysis ---
# (Assuming 'mmag_alpha_df' is loaded)
contrast_data <- run_bootstrap_contrasts(mmag_alpha_df, "Shannon", "Season", n_perms = 5000)

# Pivot for plotting
plot_contrasts <- contrast_data %>%
  pivot_longer(cols = -Iteration, names_to = "Comparison", values_to = "Difference")

# Set order of comparisons to match the timeline
plot_contrasts$Comparison <- factor(plot_contrasts$Comparison, 
                                    levels = c("Summer vs Spring", "Fall vs Summer", "Winter vs Fall", "Winter vs Spring"))

# --- 3. The Plot ---
ggplot(plot_contrasts, aes(x = Difference, y = Comparison, fill = after_stat(x))) +
  # Add the vertical line at 0 (No Change)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  
  # The Ridges
  geom_density_ridges_gradient(
    scale = 2.0, 
    rel_min_height = 0.01,
    gradient_locus = 0.05 # Adjusts where the color shift happens
  ) +
  
  # Color the gradient: Blue = Decrease, Red = Increase
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  
  labs(
    title = "Inferred Seasonal Changes in Alpha Diversity",
    subtitle = "Distributions of differences (Delta) derived from 5000 bootstraps",
    x = "Magnitude of Change (Shannon Index)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )






















# --- 1. Calculate Prevalence Metrics ---
calculate_mag_fidelity <- function(df) {
  
  # A. How many months was each Colony sampled total?
  colony_effort <- df %>%
    group_by(Colony) %>%
    summarise(Total_Months_Sampled = n_distinct(Month), .groups = 'drop')
  
  # B. How many months was each MAG present in each Colony?
  # (Assuming abundance > 0 means "Present")
  mag_presence <- df %>%
    filter(mean > 0) %>% 
    group_by(MAG_id, Colony) %>%
    summarise(Months_Present = n_distinct(Month), .groups = 'drop')
  
  # C. Combine to get Fidelity Score
  fidelity_df <- mag_presence %>%
    left_join(colony_effort, by = "Colony") %>%
    mutate(
      Fidelity_Score = Months_Present / Total_Months_Sampled,
      Stability_Class = case_when(
        Fidelity_Score >= 0.7 ~ "Resident (High)",
        Fidelity_Score >= 0.3 ~ "Transient (Med)",
        TRUE ~ "Rare (Low)"
      )
    )
  
  return(fidelity_df)
}

# Run the calculation (Ensure you filter for just ONE type, e.g., mmag or vmag)
input_data <- ecology_data %>% filter(type == "mmag", Caste == "Worker")
mag_stats <- calculate_mag_fidelity(input_data)

# --- 2. Create the Network Object ---

# --- FIX IS HERE: Force Colony names to be characters ---
nodes <- bind_rows(
  data.frame(name = as.character(unique(mag_stats$MAG_id)), type = "MAG"),
  data.frame(name = as.character(unique(mag_stats$Colony)), type = "Colony")
)

# Create Edge List
# Ensure 'to' and 'from' are also characters to match the node list
edges <- mag_stats %>%
  mutate(
    from = as.character(MAG_id), 
    to = as.character(Colony)
  ) %>%
  select(from, to, weight = Fidelity_Score, class = Stability_Class)

# Build Graph
my_graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

# Calculate Degree (Breadth)
my_graph <- my_graph %>%
  activate(nodes) %>%
  mutate(
    degree = centrality_degree(),
    # Only label Colonies
    label = ifelse(type == "Colony", name, NA) 
  )

# --- 3. Visualize ---
set.seed(123) 

network_plot <- ggraph(my_graph, layout = "stress") + 
  
  # Edges
  geom_edge_link(aes(width = weight, alpha = weight, color = class)) +
  scale_edge_width(range = c(0.2, 2)) + 
  scale_edge_alpha(range = c(0.1, 0.8)) +
  scale_edge_color_manual(values = c(
    "Resident (High)" = "#d73027", 
    "Transient (Med)" = "#fdae61", 
    "Rare (Low)" = "grey80"
  )) +
  
  # Colony Nodes
  geom_node_point(aes(filter = type == "Colony", fill = name), 
                  shape = 21, size = 8, color = "black", stroke = 1.5) +
  
  # MAG Nodes
  geom_node_point(aes(filter = type == "MAG", size = degree), 
                  color = "grey20", alpha = 0.7) +
  
  # Labels
  geom_node_text(aes(label = label), repel = TRUE, fontface = "bold", size = 5) +
  
  theme_void() + 
  labs(
    title = "Colony-MAG Association Network",
    subtitle = "Red lines indicate high stability (>70% of sampled months)",
    edge_width = "Stability"
  ) +
  theme(legend.position = "right")

print(network_plot)








# Summarize per MAG
mag_summary <- mag_stats %>%
  group_by(MAG_id) %>%
  summarise(
    Mean_Fidelity = mean(Fidelity_Score), # On average, how stable is it where it occurs?
    Colony_Breadth = n_distinct(Colony)   # How many colonies does it inhabit?
  )

ggplot(mag_summary, aes(x = Colony_Breadth, y = Mean_Fidelity)) +
  geom_jitter(width = 0.2, height = 0.02, size = 3, alpha = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    title = "MAG Lifestyle Classification",
    x = "Breadth (Number of Colonies)",
    y = "Mean Fidelity (Stability within Colony)",
    caption = "Top Left: Stable Specialists | Top Right: Core Generalists\nBottom Right: Transient Generalists"
  ) +
  theme_bw()










# ---------------------------
# Parameters
# ---------------------------
alpha_metric <- "shannon"   # or "richness"
pseudocount  <- 1e-9

# ---------------------------
# Prep: filter to core genera
# ---------------------------
genus_df <- ecology_data %>%
  mutate(Genus = tolower(as.character(Genus))) %>%
  filter(Genus %in% target_genera) %>%
  mutate(time_cont = month_to_cumulative_time_map[as.character(Month)]) %>%
  filter(!is.na(time_cont))

# ---------------------------
# Build genus-specific alpha diversity over time
# ---------------------------
alpha_results <- list()
alpha_plots   <- list()



alpha_long_all <- list()

for (g in target_genera) {
  
  g_df <- genus_df %>% filter(Genus == g)
  
  g_comm <- g_df %>%
    select(sample_name, MAG_id, mean) %>%
    group_by(sample_name, MAG_id) %>%
    summarise(abund = sum(mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = MAG_id, values_from = abund, values_fill = 0)
  
  if (nrow(g_comm) < 3 || ncol(g_comm) < 3) next
  
  g_comm_mat <- g_comm %>%
    tibble::column_to_rownames("sample_name") %>%
    as.matrix()
  
  g_comm_mat[g_comm_mat == 0] <- pseudocount
  
  alpha_vals <- vegan::diversity(g_comm_mat, index = "shannon")
  
  alpha_df <- data.frame(
    sample_name = names(alpha_vals),
    alpha = as.numeric(alpha_vals)
  ) %>%
    left_join(
      g_df %>%
        select(sample_name, Month, time_cont, Colony, Apiary) %>%
        distinct(),
      by = "sample_name"
    ) %>%
    mutate(Genus = g)
  
  alpha_long_all[[g]] <- alpha_df
}

alpha_long <- bind_rows(alpha_long_all)



str(alpha_long)
# sample_name | alpha | Month | time_cont | Colony | Apiary | Genus



alpha_mean <- alpha_long %>%
  group_by(Genus, time_cont) %>%
  summarise(mean_alpha = mean(alpha, na.rm = TRUE),
            .groups = "drop")



p_colony <- ggplot(alpha_long, aes(x = time_cont, y = alpha)) +
  geom_line(
    aes(group = Colony),
    color = "grey60",
    alpha = 0.6,
    linewidth = 0.4
  ) +
  geom_line(
    data = alpha_mean,
    aes(x = time_cont, y = mean_alpha),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.1
  ) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 4) +
  labs(
    title = "Within-genus α-diversity over time (by Colony)",
    x = "Continuous time",
    y = "Shannon α-diversity (MAGs)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )



p_apiary <- ggplot(alpha_long, aes(x = time_cont, y = alpha)) +
  geom_line(
    aes(group = Apiary),
    color = "grey60",
    alpha = 0.6,
    linewidth = 0.4
  ) +
  geom_line(
    data = alpha_mean,
    aes(x = time_cont, y = mean_alpha),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.1
  ) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 4) +
  labs(
    title = "Within-genus α-diversity over time (by Apiary)",
    x = "Continuous time",
    y = "Shannon α-diversity (MAGs)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )



library(cowplot)

final_alpha_panel <- plot_grid(
  p_colony,
  p_apiary,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 14
)

print(final_alpha_panel)



p_colony_smooth <- ggplot(alpha_long, aes(x = time_cont, y = alpha)) +
  # thin colony-level smoothers
  geom_smooth(
    aes(group = Colony),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    color = "grey65",
    linewidth = 0.5
  ) +
  # bold global mean smoother
  geom_smooth(
    data = alpha_long,
    aes(x = time_cont, y = alpha),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k=4),
    se = FALSE,
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2
  ) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 3) +
  labs(
    title = "Within-genus α-diversity over time (by Colony)",
    x = "Continuous time",
    y = "Shannon α-diversity (MAGs)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )


p_apiary_smooth <- ggplot(alpha_long, aes(x = time_cont, y = alpha)) +
  # thin apiary-level smoothers
  geom_smooth(
    aes(group = Apiary),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k=3),
    se = FALSE,
    color = "grey65",
    linewidth = 0.5
  ) +
  # bold global mean smoother
  geom_smooth(
    data = alpha_long,
    aes(x = time_cont, y = alpha),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2
  ) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 4) +
  labs(
    title = "Within-genus α-diversity over time (by Apiary)",
    x = "Continuous time",
    y = "Shannon α-diversity (MAGs)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )





final_alpha_panel <- plot_grid(
  p_colony_smooth,
  p_apiary_smooth,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 14
)

print(final_alpha_panel)




#AESTHETIC TWEAKS

y_limits <- range(alpha_long$alpha, na.rm = TRUE)

# Optional: add a small buffer so lines don’t touch the panel edges
y_pad <- diff(y_limits) * 0.05
y_limits <- c(y_limits[1] - y_pad, y_limits[2] + y_pad)




p_colony_smooth <- ggplot(alpha_long, aes(x = time_cont, y = alpha)) +
  # colony-specific smoothers
  geom_line()+
  geom_smooth(
    data = alpha_long,
    aes(x = time_cont, y = alpha),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 3),
    se = FALSE,
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2
  ) +
  facet_wrap(~ Genus, ncol = 4) +   # NOTE: no free_y
  coord_cartesian(ylim = y_limits) +
  labs(
    title = "Within-genus α-diversity over time (by Colony)",
    x = "Continuous time",
    y = "Shannon α-diversity (MAGs)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )
p_colony_smooth







p_apiary_smooth <- ggplot(alpha_long, aes(x = time_cont, y = alpha)) +
  # apiary-specific smoothers
  geom_smooth(
    aes(group = Apiary),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 3),
    se = FALSE,
    color = "grey70",
    linewidth = 0.45
  ) +
  # global mean smoother (same as above)
  geom_smooth(
    data = alpha_long,
    aes(x = time_cont, y = alpha),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 3),
    se = FALSE,
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2
  ) +
  facet_wrap(~ Genus, ncol = 4) +
  coord_cartesian(ylim = y_limits) +
  labs(
    title = "Within-genus α-diversity over time (by Apiary)",
    x = "Continuous time",
    y = "Shannon α-diversity (MAGs)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )
p_apiary_smooth





month_levels <- c("May","June","July","August","September",
                  "October","November","January","February")

alpha_long$Month <- factor(alpha_long$Month, levels = month_levels)


p_core_genus_loess <- ggplot(alpha_long, aes(x = Month, y = alpha)) +
  
  # Thin, dotted, apiary-colored colony trajectories
  geom_smooth(
    aes(group = Colony, color = Apiary),
    method = "loess",
    span = 0.4,
    se = FALSE,
    linetype = "dotted",
    alpha = 0.6,
    linewidth = 0.6
  ) +
  
  # Bold, smooth, genus-wide mean
  geom_smooth(
    aes(group = Genus),
    method = "loess",
    span = 0.6,        # slightly smoother than colonies
    se = FALSE,
    color = "black",
    linewidth = 1.4
  ) +
  
  facet_wrap(~ Genus, ncol = 4, scales = "fixed") +
  
  labs(
    title = "Within-genus α-diversity over time (core genera)",
    subtitle = "Dotted lines: individual colonies (colored by apiary); solid black: genus-wide mean",
    x = "Month",
    y = "Shannon α-diversity (MAGs)",
    color = "Apiary"
  ) +
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )


print(p_core_genus_loess)
