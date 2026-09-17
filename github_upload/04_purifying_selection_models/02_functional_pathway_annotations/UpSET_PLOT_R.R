library(ggplot2)
library(dplyr)
library(UpSetR)
library(cowplot)

df = read.csv('all_colonies_RETAINED_r_vMAGs.csv')

# Add the 'Apiary' column based on 'Colony' values
df <- df %>%
  mutate(
    Apiary = case_when(
      Colony < 104 ~ "PT",
      Colony > 500 ~ "HB",
      # The 'TRUE' catches all remaining conditions, which is Colony >= 104 and Colony <= 500
      TRUE ~ "OB" 
    )
  )

# Create the binary matrix
binary_matrix <- df %>%
  distinct(r_MAG_id, Colony) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Colony,
    values_from = present,
    values_fill = 0
  ) %>%
  select(-r_MAG_id) %>%
  as.data.frame()

# Define your desired custom order for the Colony sets
custom_colony_order <- c("100", "101", "102", "103", "104", "105", "106", "107", "555", "777", "888", "999")
custom_colony_order <- as.character(custom_colony_order)

# Filter the custom order to only include colonies that are actually present in your data
present_colony_order <- intersect(custom_colony_order, colnames(binary_matrix))

# Robustness Check
if (length(present_colony_order) == 0) {
  stop("No matching Colony IDs found between your custom order and the binary matrix columns. Upset plot cannot be generated.")
}

# Reorder the columns of the binary_matrix to match the desired display order.
binary_matrix <- binary_matrix[, present_colony_order, drop = FALSE]

# Generate the upset plot with custom set order and disabled automatic intersection ordering
# Removed 'keep.order = TRUE'
upset(binary_matrix, 
      sets = present_colony_order,       # Pass your custom ordered vector here
      nsets = min(ncol(binary_matrix), 12), # Display up to 12 sets (or fewer if less are present)
      mainbar.y.label = 'Intersection Size',
      sets.x.label = 'Set Size',
      text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),
      order.by = "freq",keep.order=TRUE) # Keep order.by = "null"






# --- Data Preparation for Apiary Intersection ---
# Create a binary matrix where columns are Apiary types and rows are r_MAG_ids
apiary_binary_matrix <- df %>%
  distinct(r_MAG_id, Apiary) %>% # Consider unique r_MAG_id within each Apiary
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Apiary,
    values_from = present,
    values_fill = 0
  ) %>%
  select(-r_MAG_id) %>%
  as.data.frame()

# Define a desired order for Apiary sets (optional, but good practice)
custom_apiary_order <- c("PT", "OB", "HB") # Example order

# Filter to only include Apiaries actually present in the data
present_apiary_order <- intersect(custom_apiary_order, colnames(apiary_binary_matrix))

# Reorder columns of the binary matrix
if (length(present_apiary_order) > 0) {
  apiary_binary_matrix <- apiary_binary_matrix[, present_apiary_order, drop = FALSE]
} else {
  message("No Apiary types found in the data to plot intersections.")
}

# --- Generate Upset Plot for Apiary Intersections ---
if (ncol(apiary_binary_matrix) > 0) {
  upset(apiary_binary_matrix,
        sets = present_apiary_order,
        nsets = min(ncol(apiary_binary_matrix), 3), # Limit to 3 sets (PT, OB, HB)
        mainbar.y.label = 'r_MAG_id Intersection Size (Apiary)',
        sets.x.label = 'Apiary Set Size',
        text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),
        order.by = "freq", # Order intersections by frequency
        keep.order = TRUE) # Keep the order of sets as specified
} else {
  message("Apiary binary matrix is empty or has no columns. Cannot generate plot.")
}



# --- Define custom colors for intersection degrees ---
degree_queries <- list(
  list(query = intersects, params = list(degree = 1), color = "gray", active = TRUE),
  list(query = intersects, params = list(degree = 2), color = "darkgray", active = TRUE),
  list(query = intersects, params = list(degree = 3), color = "black", active = TRUE)
  # You can add more queries here for higher degrees if needed,
  # or any intersections with degree > 3 will take the default main.bar.color.
)

# --- Generate Upset Plot for Colony Intersections, colored by Degree ---
# We will use 'group.by = "degree"' to color, and let UpSetR assign default colors
# since other direct coloring methods are causing errors.
if (ncol(binary_matrix) > 0) {
  upset(binary_matrix, 
        sets = present_colony_order,
        nsets = min(ncol(binary_matrix), 12),
        mainbar.y.label = 'r_MAG_id Intersection Size',
        sets.x.label = 'Colony Set Size',
        text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),
        order.by = "freq", # Order intersections by frequency
        group.by = "degree",
        color.pal = # Group and color intersections by their degree
        # Removed 'colors = degree_colors' as it caused an error
        # Removed 'queries' as it's causing an internal error
        keep.order = TRUE)
} else {
  message("Colony binary matrix is empty or has no columns. Cannot generate plot for Block 2.")
}


# Create a binary matrix where columns are Month values and rows are r_MAG_ids
month_binary_matrix <- df %>%
  distinct(r_MAG_id, Month) %>% # Consider unique r_MAG_id within each Month
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Month,
    values_from = present,
    values_fill = 0
  ) %>%
  select(-r_MAG_id) %>%
  as.data.frame()

# Define your desired custom order for the Month sets
# Ensure these match the month names as they appear in your data (case-sensitive)
custom_month_order <- c("May", "June", "July", "August", "September", "October", "November", "January", "February")

# Filter the custom order to only include months that are actually present in your data
present_month_order <- intersect(custom_month_order, colnames(month_binary_matrix))

# Robustness Check
if (length(present_month_order) == 0) {
  stop("No matching Month names found between your custom order and the binary matrix columns. Upset plot cannot be generated.")
}

# Reorder the columns of the month_binary_matrix to match the desired display order.
month_binary_matrix <- month_binary_matrix[, present_month_order, drop = FALSE]

# --- Generate Upset Plot for Month Intersections ---
if (ncol(month_binary_matrix) > 0) {
  upset(month_binary_matrix,
        sets = present_month_order,       # Pass your custom ordered Month vector here
        nsets = min(ncol(month_binary_matrix), 12), # Display up to 12 sets (or fewer if less are present)
        mainbar.y.label = 'r_MAG_id Intersection Size (Month)',
        sets.x.label = 'Month Set Size',
        text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),
        order.by = "freq", # Order intersections by frequency
        keep.order = TRUE) # Keep the order of sets as specified
} else {
  message("Month binary matrix is empty or has no columns. Cannot generate plot.")
}






# --- 1. Load Data (if not already in your environment) ---
# df = read.csv('all_colonies_RETAINED_r_vMAGs.csv')

# --- 2. Define Custom Colony Order (same as before) ---
custom_colony_order <- c("100", "101", "102", "103", "104", "105", "106", "107", "555", "777", "888", "999")
custom_colony_order <- as.character(custom_colony_order)

# --- 3. Create Separate Data Frames for Each Lifestyle ---

# Prophage Data Frame
df_prophage <- df %>%
  filter(lifestyle == "Prophage")

# Lytic Data Frame
df_lytic <- df %>%
  filter(lifestyle == "Lytic")

# Temperate Data Frame
df_temperate <- df %>%
  filter(lifestyle == "Temperate")

# --- 4. Function to Generate Upset Plot for a Given Lifestyle DF ---
# This function encapsulates the steps to avoid repetition
generate_lifestyle_upset <- function(data_frame, lifestyle_name, custom_colony_order_list) {
  message(paste("\n--- Generating Upset Plot for Lifestyle:", lifestyle_name, "---"))
  
  # Check if data frame is empty
  if (nrow(data_frame) == 0) {
    message(paste("  No data found for lifestyle:", lifestyle_name, ". Skipping plot generation."))
    return(invisible(NULL)) # Exit function if no data
  }
  
  # Create the binary matrix
  binary_matrix <- data_frame %>%
    distinct(r_MAG_id, Colony) %>%
    mutate(present = 1) %>%
    pivot_wider(
      names_from = Colony,
      values_from = present,
      values_fill = 0
    ) %>%
    select(-r_MAG_id) %>%
    as.data.frame()
  
  # Ensure columns are present and sum is not zero
  present_colony_order_subset <- intersect(custom_colony_order_list, colnames(binary_matrix))
  
  if (length(present_colony_order_subset) < 2 || sum(binary_matrix[, present_colony_order_subset, drop = FALSE]) == 0) {
    message(paste("  Insufficient data or no presences in binary matrix for lifestyle:", lifestyle_name, ". Skipping plot generation."))
    return(invisible(NULL))
  }
  
  # Reorder columns
  binary_matrix <- binary_matrix[, present_colony_order_subset, drop = FALSE]
  
  # Generate the plot
  upset(binary_matrix,
        sets = present_colony_order_subset,
        nsets = min(ncol(binary_matrix), 12),
        mainbar.y.label = paste0('r_MAG_id Intersection Size (', lifestyle_name, ')'),
        sets.x.label = 'Colony Set Size',
        text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),
        order.by = "freq",
        keep.order = TRUE)
}

# --- 5. Generate Plots for Each Lifestyle Data Frame ---

# Plot for Prophage
generate_lifestyle_upset(df_prophage, "Prophage", custom_colony_order)

# Plot for Lytic
generate_lifestyle_upset(df_lytic, "Lytic", custom_colony_order)

# Plot for Temperate
generate_lifestyle_upset(df_temperate, "Temperate", custom_colony_order)

# 1. Calculate the interaction degree for each unique r_MAG_id
r_mag_interaction_degree <- df %>%
  group_by(r_MAG_id) %>%
  summarise(
    num_colonies = n_distinct(Colony),
    .groups = 'drop'
  )

# 2. Get the lifestyle for each unique r_MAG_id
unique_r_mag_lifestyle_map <- df %>%
  select(r_MAG_id, lifestyle) %>%
  distinct(r_MAG_id, .keep_all = TRUE)

# Combine degree and lifestyle, ensuring unique r_MAG_id per row
combined_r_mag_data_for_plot <- r_mag_interaction_degree %>%
  left_join(unique_r_mag_lifestyle_map, by = "r_MAG_id")

# --- APPLY THE FILTER: num_colonies > 1 as requested ---
filtered_r_mags_data_colonies <- combined_r_mag_data_for_plot %>%
  filter(num_colonies > 1)

# 3. Summarize for plotting: Count r_MAG_ids by num_colonies and lifestyle (from filtered data)
plot_data_lifestyle_counts_filtered <- filtered_r_mags_data_colonies %>%
  group_by(num_colonies, lifestyle) %>%
  summarise(
    count_of_r_mags = n(),
    .groups = 'drop'
  )

# 4. Ensure 'num_colonies' is treated as a discrete factor for the x-axis categories
plot_data_lifestyle_counts_filtered$num_colonies <- as.factor(plot_data_lifestyle_counts_filtered$num_colonies)

# --- Calculate Total Bin Sizes for Labels ---
# Create a new dataframe to hold the total count for each 'num_colonies' category (after filtering).
total_bin_counts_colonies <- plot_data_lifestyle_counts_filtered %>%
  group_by(num_colonies) %>%
  summarise(
    total_count = sum(count_of_r_mags),
    .groups = 'drop'
  )

# --- Plotting ---
message("\nGenerating plot for r_MAG_id interaction with Colonies (filtered num_colonies > 1, with total bin labels)...")

if (nrow(plot_data_lifestyle_counts_filtered) == 0) {
  message("No r_MAG_ids with more than 1 colony interaction found for plotting. Skipping plot generation.")
} else {
  # Dynamically determine appropriate breaks for the Y-axis based on the filtered data's range.
  # The maximum count will be lower since num_colonies = 1 is excluded.
  max_filtered_colony_count_val <- max(plot_data_lifestyle_counts_filtered$count_of_r_mags)
  
  # Heuristic for breaks:
  if (max_filtered_colony_count_val > 500) { # Adjusted thresholds for potentially smaller numbers
    colony_y_breaks <- seq(0, ceiling(max_filtered_colony_count_val / 100) * 100, by = 100)
  } else if (max_filtered_colony_count_val > 50) {
    colony_y_breaks <- seq(0, ceiling(max_filtered_colony_count_val / 10) * 10, by = 10)
  } else {
    colony_y_breaks <- unique(c(0, pretty(c(0, max_filtered_colony_count_val), n = 5))) # Let pretty() determine for very small ranges
  }
  
  ggplot(plot_data_lifestyle_counts_filtered, aes(x = num_colonies, y = count_of_r_mags, fill = lifestyle)) +
    geom_col(position = "stack", color = "black", width = 0.7) +
    # Add geom_text layer to display the total count for each bar
    geom_text(
      data = total_bin_counts_colonies, # Use the new dataframe for labels
      aes(x = num_colonies, y = total_count, label = scales::comma(total_count)), # Y is the total height, label is formatted total
      inherit.aes = FALSE, # IMPORTANT: Prevents inheriting 'fill'
      vjust = -0.5, # Adjust vertical position (slightly above the bar)
      size = 3.5, # Text size
      color = "black" # Text color
    ) +
    labs(
      title = "Unique r_MAG_ids by Colony Interaction Count (num_colonies > 1) and Lifestyle",
      x = "Number of Colonies Interacted With",
      y = "Count of Unique r_MAG_ids",
      fill = "Lifestyle"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right"
    )
}



