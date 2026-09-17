# Load necessary libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# --- Load and Prepare Data ---
# Make sure your CSV file has a 'Genus' column.
ecology_data <- read.csv("ecology_metadata_df.csv", header = TRUE)
ecology_data$MAG_id <- as.factor(ecology_data$MAG_id)
ecology_data$Genus <- as.factor(ecology_data$Genus) # Ensure Genus is a factor

# --- Helper Functions for Data Derivation ---
get_season <- function(month) {
  if (month %in% c("May", "June")) {
    return('Spring')
  } else if (month %in% c("July", "August", "September")) {
    return('Summer')
  } else if (month %in% c("October", "November")) {
    return('Fall')
  } else if (month %in% c("January", "February")) {
    return('Winter')
  }
  return(NA)
}

get_apiary <- function(colony) {
  if (colony > 500) {
    return('HB')
  } else if (colony < 104) {
    return('PT')
  }
  return('OB')
}

# --- Derive Metadata Columns ---
ecology_data$Season <- sapply(ecology_data$Month, get_season)
ecology_data$Phenotype <- sapply(ecology_data$Season,
                                 function(x) ifelse(x == 'Winter', 'Winter Bee', 'Non-winter Bee'))
ecology_data$Apiary <- sapply(ecology_data$Colony, get_apiary)

# Ensure metadata columns are factors
ecology_data$Season <- as.factor(ecology_data$Season)
ecology_data$Apiary <- as.factor(ecology_data$Apiary)
ecology_data$Caste <- as.factor(ecology_data$Caste)
ecology_data$Month <- as.factor(ecology_data$Month)
ecology_data$State <- as.factor(ecology_data$State)
ecology_data$Colony <- as.factor(ecology_data$Colony)
ecology_data$Rep <- as.factor(ecology_data$Rep)
ecology_data$Phenotype <- as.factor(ecology_data$Phenotype)

# --- Define Custom Colors and Shapes (Global for all plots) ---
season_colors <- c("Spring" = "#6a994e", "Summer" = "#c08497", "Fall" = "#f7b801", "Winter" = "#118ab2")
apiary_shapes <- c("HB" = 22, "OB" = 21, "PT" = 24)
colony_colors <- c(
  "100" = "#a6cee3", "101" = "#1f78b4", "102" = "#cab2d6", "103" = "#6a3d9a",
  "104" = "#fb9a99", "105" = "#e31a1c", "106" = "#fdbf6f", "107" = "#ff7f00",
  "555" = "#f1c40f", "777" = "#b15928", "888" = "#b2df8a", "999" = "#33a02c"
)
phenotype_colors <- c("Winter Bee" = "darkred", "Non-winter Bee" = "darkblue")

# --- Global lists to store all generated plots and tables ---
all_nmds_plots <- list()
all_statistical_tables <- list()

# --- Helper function to run tests and collect results ---
run_tests_and_collect <- function(community_matrix, metadata_df, factors_to_test) {
  anosim_results_raw <- list()
  permanova_results_raw <- list()
  
  anosim_p_values <- numeric()
  anosim_p_names <- character()
  permanova_p_values <- numeric()
  permanova_p_names <- character()
  
  bray_curtis_dist <- vegdist(sqrt(community_matrix), method = "bray")
  
  for (factor_name in factors_to_test) {
    if (factor_name %in% colnames(metadata_df) && length(unique(metadata_df[[factor_name]])) > 1) {
      # ANOSIM
      anosim_res <- anosim(bray_curtis_dist, metadata_df[[factor_name]], permutations = 999)
      anosim_results_raw[[factor_name]] <- list(R = anosim_res$statistic, p = anosim_res$signif)
      anosim_p_values <- c(anosim_p_values, anosim_res$signif)
      anosim_p_names <- c(anosim_p_names, factor_name)
      
      # PERMANOVA
      formula_str <- as.formula(paste("bray_curtis_dist ~", factor_name))
      adonis_res <- adonis2(formula_str, data = metadata_df, permutations = 999)
      permanova_results_raw[[factor_name]] <- list(R2 = adonis_res$R2[1], p = adonis_res$`Pr(>F)`[1])
      permanova_p_values <- c(permanova_p_values, adonis_res$`Pr(>F)`[1])
      permanova_p_names <- c(permanova_p_names, factor_name)
      
    } else {
      # Add NA placeholders for skipped factors
      anosim_results_raw[[factor_name]] <- list(R = NA_real_, p = NA_real_)
      anosim_p_values <- c(anosim_p_values, NA_real_)
      anosim_p_names <- c(anosim_p_names, factor_name)
      
      permanova_results_raw[[factor_name]] <- list(R2 = NA_real_, p = NA_real_)
      permanova_p_values <- c(permanova_p_values, NA_real_)
      permanova_p_names <- c(permanova_p_names, factor_name)
    }
  }
  
  # Apply FDR correction
  all_raw_p_values <- c(anosim_p_values, permanova_p_values)
  
  valid_p_indices <- !is.na(all_raw_p_values)
  if (any(valid_p_indices)) {
    adjusted_p_values <- rep(NA_real_, length(all_raw_p_values))
    adjusted_p_values[valid_p_indices] <- p.adjust(all_raw_p_values[valid_p_indices], method = "fdr")
  } else {
    adjusted_p_values <- all_raw_p_values
  }
  
  adj_anosim_p_values <- adjusted_p_values[1:length(anosim_p_values)]
  adj_permanova_p_values <- adjusted_p_values[(length(anosim_p_values) + 1):length(adjusted_p_values)]
  
  # Build the result table
  results_table <- tibble(
    Factor = character(),
    Test = character(),
    Statistic = numeric(),
    Raw_P_value = numeric(),
    Adjusted_P_value_FDR = numeric()
  )
  
  for (i in seq_along(anosim_p_names)) {
    factor_name <- anosim_p_names[i]
    results_table <- bind_rows(results_table, tibble(
      Factor = factor_name,
      Test = "ANOSIM",
      Statistic = anosim_results_raw[[factor_name]]$R,
      Raw_P_value = anosim_results_raw[[factor_name]]$p,
      Adjusted_P_value_FDR = adj_anosim_p_values[i]
    ))
  }
  
  for (i in seq_along(permanova_p_names)) {
    factor_name <- permanova_p_names[i]
    results_table <- bind_rows(results_table, tibble(
      Factor = factor_name,
      Test = "PERMANOVA",
      Statistic = permanova_results_raw[[factor_name]]$R2,
      Raw_P_value = permanova_results_raw[[factor_name]]$p,
      Adjusted_P_value_FDR = adj_permanova_p_values[i]
    ))
  }
  
  return(list(anosim_results = anosim_results_raw, permanova_results = permanova_results_raw, results_table = results_table))
}


# --- MODIFIED Analysis Function: Workers Only (by Genus) ---
run_workers_only_analysis_genus <- function(mag_type_filter) {
  analysis_name <- "Workers_Only_by_Genus"
  print(paste("### Starting Analysis for:", analysis_name, " - MAG Type:", mag_type_filter, "###"))
  
  filtered_data <- ecology_data %>%
    filter(Caste == "Worker")
  
  if (mag_type_filter != "all") {
    filtered_data <- filtered_data %>% filter(type == mag_type_filter)
  }
  
  if (nrow(filtered_data) == 0) {
    print(paste("  -> No data for this filter. Skipping."))
    return(list(plots = list(), results_table = tibble()))
  }
  
  print("  -> Preparing community matrix at Genus level...")
  community_long <- filtered_data %>%
    select(sample_name, Genus, mean) %>%
    group_by(sample_name, Genus) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop') %>%
    ungroup()
  
  community_matrix_tibble <- community_long %>%
    pivot_wider(names_from = Genus, values_from = total_abundance, values_fill = 0)
  
  community_matrix <- as.data.frame(community_matrix_tibble)
  rownames(community_matrix) <- community_matrix$sample_name
  community_matrix <- community_matrix %>% select(-sample_name)
  community_matrix <- as.matrix(community_matrix)
  
  community_matrix <- community_matrix[rowSums(community_matrix) > 0, , drop = FALSE]
  
  if (nrow(community_matrix) < 3) {
    warning("  -> Not enough data points for NMDS. Skipping.")
    return(list(plots = list(), results_table = tibble()))
  }
  
  print("  -> Preparing metadata...")
  metadata_df_raw <- filtered_data %>%
    select(sample_name, State, Month, Caste, Colony, Rep, Season, Phenotype, Apiary) %>%
    distinct(sample_name, .keep_all = TRUE)
  metadata_df <- data.frame(sample_name = rownames(community_matrix)) %>%
    left_join(metadata_df_raw, by = "sample_name")
  rownames(metadata_df) <- metadata_df$sample_name
  
  print("  -> Running NMDS...")
  set.seed(123)
  nmds_result <- metaMDS(sqrt(community_matrix), distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
  stress_value <- nmds_result$stress
  
  plot_data <- as.data.frame(scores(nmds_result, display = "sites"))
  plot_data <- cbind(plot_data, metadata_df)
  
  factors_for_test <- c("Season", "Apiary", "Caste", "Colony", "State", "Phenotype")
  test_results <- run_tests_and_collect(community_matrix, metadata_df, factors_for_test)
  
  plots_for_this_analysis <- list()
  plot_title_suffix_base <- paste(analysis_name, " - MAG:", mag_type_filter)
  
  # Plot 1: Color by Colony, Shape by Apiary
  p1 <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = Colony, shape = Apiary)) +
    geom_point(size = 2) +
    scale_color_manual(values = colony_colors) +
    scale_shape_manual(values = apiary_shapes) +
    theme_bw() +
    labs(title = paste("NMDS -", plot_title_suffix_base), subtitle = paste("Stress:", formatC(stress_value, digits = 3)))
  plots_for_this_analysis[["by_Colony_Apiary"]] <- p1
  
  return(list(plots = plots_for_this_analysis, results_table = test_results$results_table))
}

# --- MODIFIED Analysis Function: Apiary OB & Workers Only (by Genus) ---
run_ob_workers_analysis_genus <- function(mag_type_filter) {
  analysis_name <- "Apiary_OB_Workers_Only_by_Genus"
  print(paste("### Starting Analysis for:", analysis_name, " - MAG Type:", mag_type_filter, "###"))
  
  filtered_data <- ecology_data %>%
    filter(Apiary == "OB" & Caste == "Worker")
  
  if (mag_type_filter != "all") {
    filtered_data <- filtered_data %>% filter(type == mag_type_filter)
  }
  
  if (nrow(filtered_data) == 0) {
    print(paste("  -> No data for this filter. Skipping."))
    return(list(plots = list(), results_table = tibble()))
  }
  
  print("  -> Preparing community matrix at Genus level...")
  community_long <- filtered_data %>%
    select(sample_name, Genus, mean) %>%
    group_by(sample_name, Genus) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop')
  
  community_matrix_tibble <- community_long %>%
    pivot_wider(names_from = Genus, values_from = total_abundance, values_fill = 0)
  
  community_matrix <- as.data.frame(community_matrix_tibble)
  rownames(community_matrix) <- community_matrix$sample_name
  community_matrix <- community_matrix %>% select(-sample_name)
  community_matrix <- as.matrix(community_matrix)
  
  community_matrix <- community_matrix[rowSums(community_matrix) > 0, , drop = FALSE]
  
  if (nrow(community_matrix) < 3) {
    warning("  -> Not enough data points for NMDS. Skipping.")
    return(list(plots = list(), results_table = tibble()))
  }
  
  print("  -> Preparing metadata...")
  metadata_df_raw <- filtered_data %>%
    select(sample_name, State, Month, Caste, Colony, Rep, Season, Phenotype, Apiary) %>%
    distinct(sample_name, .keep_all = TRUE)
  metadata_df <- data.frame(sample_name = rownames(community_matrix)) %>%
    left_join(metadata_df_raw, by = "sample_name")
  rownames(metadata_df) <- metadata_df$sample_name
  
  print("  -> Running NMDS...")
  set.seed(123)
  nmds_result <- metaMDS(sqrt(community_matrix), distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
  stress_value <- nmds_result$stress
  
  plot_data <- as.data.frame(scores(nmds_result, display = "sites"))
  plot_data <- cbind(plot_data, metadata_df)
  
  factors_for_test <- c("Season", "Phenotype", "Colony")
  test_results <- run_tests_and_collect(community_matrix, metadata_df, factors_for_test)
  
  plots_for_this_analysis <- list()
  plot_title_suffix_base <- paste(analysis_name, " - MAG:", mag_type_filter)
  
  # Plot 1: Color by Phenotype
  p1 <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = Phenotype)) +
    geom_point(size = 2) +
    scale_color_manual(values = phenotype_colors) +
    theme_bw() +
    labs(title = paste("NMDS -", plot_title_suffix_base), subtitle = paste("Stress:", formatC(stress_value, digits = 3)))
  plots_for_this_analysis[["by_Phenotype"]] <- p1
  
  return(list(plots = plots_for_this_analysis, results_table = test_results$results_table))
}

# --- MODIFIED Analysis Function: All Data (by Genus) ---
run_all_data_analysis_genus <- function(mag_type_filter) {
  analysis_name <- "All_Data_by_Genus"
  print(paste("### Starting Analysis for:", analysis_name, " - MAG Type:", mag_type_filter, "###"))
  
  filtered_data <- ecology_data
  
  if (mag_type_filter != "all") {
    filtered_data <- filtered_data %>% filter(type == mag_type_filter)
  }
  
  if (nrow(filtered_data) == 0) {
    print(paste("  -> No data for this filter. Skipping."))
    return(list(plots = list(), results_table = tibble()))
  }
  
  print("  -> Preparing community matrix at Genus level...")
  community_long <- filtered_data %>%
    select(sample_name, Genus, mean) %>%
    group_by(sample_name, Genus) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = 'drop')
  
  community_matrix_tibble <- community_long %>%
    pivot_wider(names_from = Genus, values_from = total_abundance, values_fill = 0)
  
  community_matrix <- as.data.frame(community_matrix_tibble)
  rownames(community_matrix) <- community_matrix$sample_name
  community_matrix <- community_matrix %>% select(-sample_name)
  community_matrix <- as.matrix(community_matrix)
  
  community_matrix <- community_matrix[rowSums(community_matrix) > 0, , drop = FALSE]
  
  if (nrow(community_matrix) < 3) {
    warning("  -> Not enough data points for NMDS. Skipping.")
    return(list(plots = list(), results_table = tibble()))
  }
  
  print("  -> Preparing metadata...")
  metadata_df_raw <- filtered_data %>%
    select(sample_name, State, Month, Caste, Colony, Rep, Season, Phenotype, Apiary) %>%
    distinct(sample_name, .keep_all = TRUE)
  metadata_df <- data.frame(sample_name = rownames(community_matrix)) %>%
    left_join(metadata_df_raw, by = "sample_name")
  rownames(metadata_df) <- metadata_df$sample_name
  
  print("  -> Running NMDS...")
  set.seed(123)
  nmds_result <- metaMDS(sqrt(community_matrix), distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
  stress_value <- nmds_result$stress
  
  plot_data <- as.data.frame(scores(nmds_result, display = "sites"))
  plot_data <- cbind(plot_data, metadata_df)
  
  factors_for_test <- c("Season", "Apiary", "Caste", "Colony", "State", "Phenotype")
  test_results <- run_tests_and_collect(community_matrix, metadata_df, factors_for_test)
  
  plots_for_this_analysis <- list()
  plot_title_suffix_base <- paste(analysis_name, " - MAG:", mag_type_filter)
  
  # Plot 1: Color by Season, Shape by Apiary
  p1 <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = Season, shape = Apiary)) +
    geom_point(size = 2) +
    scale_color_manual(values = season_colors) +
    scale_shape_manual(values = apiary_shapes) +
    theme_bw() +
    labs(title = paste("NMDS -", plot_title_suffix_base), subtitle = paste("Stress:", formatC(stress_value, digits = 3)))
  plots_for_this_analysis[["by_Season_Apiary"]] <- p1
  
  return(list(plots = plots_for_this_analysis, results_table = test_results$results_table))
}

# --- Main Execution Block ---

mag_types_to_run <- c("mmag") # Excludes "vMAG"

for (mag_type in mag_types_to_run) {
  # Run Workers Only Analysis
  workers_results <- run_workers_only_analysis_genus(mag_type)
  analysis_id <- paste("Workers_Only_Genus", mag_type, sep = "_")
  all_nmds_plots[[analysis_id]] <- workers_results$plots
  all_statistical_tables[[analysis_id]] <- workers_results$results_table
  
  # Run OB Workers Analysis
  ob_workers_results <- run_ob_workers_analysis_genus(mag_type)
  analysis_id <- paste("OB_Workers_Genus", mag_type, sep = "_")
  all_nmds_plots[[analysis_id]] <- ob_workers_results$plots
  all_statistical_tables[[analysis_id]] <- ob_workers_results$results_table
  
  # Run All Data Analysis
  all_data_results <- run_all_data_analysis_genus(mag_type)
  analysis_id <- paste("All_Data_Genus", mag_type, sep = "_")
  all_nmds_plots[[analysis_id]] <- all_data_results$plots
  all_statistical_tables[[analysis_id]] <- all_data_results$results_table
}

# --- To view your results ---
# Print a specific plot, for example:
print(all_nmds_plots$Workers_Only_Genus_mmag$by_Colony_Apiary)
print(all_nmds_plots$OB_Workers_Genus_mmag$by_Phenotype)
print(all_nmds_plots$All_Data_Genus_mmag$by_Season_Apiary)

# View a specific statistical table
print(all_statistical_tables$Workers_Only_Genus_mmag)
print(all_statistical_tables$OB_Workers_Genus_mmag)
print(all_statistical_tables$All_Data_Genus_mmag)

# You can also save plots and tables to files as needed.
