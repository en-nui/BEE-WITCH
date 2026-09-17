library(dplyr)
library(tidyr)
library(vroom)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
# --- 1. CONFIGURATION ---
# Define the paths to your data
diamond_results_dir <- "/home/robinch/projects/BEE-WITCH/41_mmag_ani/diamond_results"
taxonomy_file <- "/home/robinch/projects/BEE-WITCH/15_MAG_dRep/RENAMED_CONTIG_TAXONOMY_20250528.csv"
output_plot_path <- "/home/robinch/projects/BEE-WITCH/41_mmag_ani/ani_by_genus_distribution.png"

# --- 2. LOAD AND PREPARE TAXONOMY DATA ---
# Load the taxonomy file and select only the columns we need
taxonomy_df <- read_csv(taxonomy_file) %>%
  select(MAG_id, Genus) %>%
  distinct() # Ensure one entry per MAG_id

taxonomy_df_family <- read_csv(taxonomy_file) %>%
  select(MAG_id, Family) %>%
  distinct() # Ensure one entry per MAG_id

# --- 3. PROCESS ALL DIAMOND HITS ---
# Load ALL hits, not just the single best ones per protein,
# as we need the full dataset for the new Step 4 logic.
diamond_files <- list.files(path = diamond_results_dir, pattern = "\\.tsv$", full.names = TRUE)
diamond_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# This data frame will be large
all_hits_df <- vroom(diamond_files, id = "filepath", col_names = diamond_cols, delim = "\t") %>%
  mutate(MAG_id = str_remove(basename(filepath), "\\.tsv$")) %>%
  select(MAG_id, qseqid, sseqid, pident)
# --- END MODIFIED LINE ---

# --- 4. CREATE THE GENOME-LEVEL 'ANI' PROXY (True Top Hit Method) ---

mag_level_ani_proxy <- all_hits_df %>%
  group_by(MAG_id) %>%
  nest() %>% # Nests all data for each MAG
  mutate(
    ANI_proxy = map_dbl(data, ~{
      # .x is the tibble of *all* hits for one MAG
      mag_data <- .x 
      
      # 1. Find the "global best hit" for each protein
      global_best_hits <- mag_data %>%
        group_by(qseqid) %>%
        slice_max(order_by = pident, n = 1, with_ties = FALSE) %>%
        ungroup()
      
      # 2. From this, find the single most frequent sseqid ("Genome A")
      top_sseqid <- global_best_hits %>%
        count(sseqid, sort = TRUE) %>%
        slice(1) %>%
        pull(sseqid)
      
      # Handle case where no hits are found
      if (length(top_sseqid) == 0) {
        return(NA_real_)
      }
      
      # 3. Go back to the *full* MAG data. Filter for *all* hits to "Genome A".
      genome_A_hits <- mag_data %>%
        filter(sseqid == top_sseqid)
      
      # 4. From this, find the *best hit* to "Genome A" for each protein
      #    (This fulfills your "keep only the highest ANI for that
      #     protein, genome A combination" requirement)
      best_hits_to_A <- genome_A_hits %>%
        group_by(qseqid) %>%
        slice_max(order_by = pident, n = 1, with_ties = FALSE) %>%
        ungroup()
      
      # 5. Calculate the mean pident of *all* these best hits.
      mean_pident <- mean(best_hits_to_A$pident, na.rm = TRUE)
      
      return(mean_pident)
    })
  ) %>%
  select(MAG_id, ANI_proxy) %>%
  ungroup()



# --- 5. CALCULATE RELATIVE DENSITY ---
relative_density_data <- best_hits_df %>%
  group_by(MAG_id) %>%
  nest() %>%
  mutate(
    density_calc = map(data, ~density(.x$pident, from = 0, to = 100, n = 500,adjust=2)),
    density_df = map(density_calc, ~tibble(pident = .x$x, density = .x$y)),
    relative_density_df = map(density_df, ~mutate(.x, relative_density = density / max(density, na.rm = TRUE)))
  ) %>%
  select(MAG_id, relative_density_df) %>%
  unnest(cols = relative_density_df)

# 1. Create the base plot_data
# This joins the densities with taxonomy and our new "Top Hit" ANI proxy
print("Creating base plot_data...")
plot_data <- relative_density_data %>%
  left_join(taxonomy_df, by = "MAG_id") %>%
  left_join(mag_level_ani_proxy, by = "MAG_id") %>%
  filter(!is.na(Genus) & !is.na(ANI_proxy))

# 2. Create the final binned and filtered data for plotting
# This takes plot_data and applies all our rules for the final plots
print("Creating binned and filtered plot_data_binned...")
plot_data_binned <- plot_data %>%
  # Round the ANI proxy to create integer bins
  mutate(ANI_bin = round(ANI_proxy)) %>%
  
  # Filter to only include bins in our 75-100 range
  filter(ANI_bin >= 75 & ANI_bin <= 100) %>%
  
  # Create the ordered factor for the discrete color scale
  mutate(ANI_factor = factor(ANI_bin, levels = factor_levels)) %>% 
  
  # Filter down to only the 10 genera we want to plot
  filter(Genus %in% genera_to_keep)

print("plot_data and plot_data_binned are ready.")

# --- NEW: CREATE THE GRUVBOX-INSPIRED PALETTE ---
# We need 26 colors to map to the 26 integer levels (75 to 100)
n_colors <- 26
factor_levels <- 100:75

# --- MODIFIED: New color ramp ---
# This palette interpolates along your new anchor points:
# Red -> Blue -> Deeper Blue -> Gray -> White
gruvbox_palette_ramp <- colorRampPalette(c(
  "#cc241d",  # Gruvbox Red (at 100%)
  "#83a598",  # Gruvbox Aqua (the "more blue" around 97%)
  "#458588",  # Gruvbox Darker Blue (the "very blue" at 95%)
  "#a89984",  # Gruvbox Neutral/Beige (the "gray")
  "#f0f0f0"   # Light Gray/White (at 75%)
))(n_colors)
# --- END MODIFIED PART ---
names(gruvbox_palette_ramp) <- factor_levels
# --- 6. JOIN ALL DATA and PLOT 🎨 ---
# Join the calculated density curves with taxonomy AND our new ANI proxy
plot_data <- relative_density_data %>%
  left_join(taxonomy_df, by = "MAG_id") %>%
  left_join(mag_level_ani_proxy, by = "MAG_id") %>%
  filter(!is.na(Genus) & !is.na(ANI_proxy))

# --- NEW: Define the list of genera you want to keep ---
genera_to_keep <- c(
  "apilactobacillus", "bartonella", "bifidobacterium", "bombella", 
  "bombilactobacillus", "commensalibacter", "frischella", "gilliamella", 
  "lactobacillus", "snodgrassella"
)
# --- END NEW PART ---

# Re-create plot_data_binned from plot_data (as in the previous step)
# AND add the new filter for specific genera
plot_data_binned <- plot_data %>%
  mutate(ANI_bin = round(ANI_proxy)) %>%
  filter(ANI_bin >= 25 & ANI_bin <= 100) %>%
  mutate(ANI_factor = factor(ANI_bin, levels = 100:25))# %>%
  
  # --- NEW: Filter for the specified genera ---
  #filter(Genus %in% genera_to_keep)
# --- END NEW PART ---


# --- 3. CREATE THE FINAL PLOT ---
# This ggplot code is the same as you provided,
# but it now uses the *filtered* 'plot_data_binned'
recreated_plot <- ggplot(
  plot_data_binned,
  aes(x = pident, y = relative_density, group = MAG_id, color = ANI_factor)
) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  
  # facet_wrap will now only create panels for the 10 genera
  facet_wrap(~ Genus) + 
  
  scale_color_manual(
    values = gruvbox_palette_ramp,
    name = "Avg. Protein ID (%)",
    drop = FALSE # Ensures all colors/levels appear in legend
  ) +
  
  scale_x_continuous(breaks = c(0, 100)) +
  xlim(75, 100) +
  labs(
    title = "Distribution of Best-Hit Protein Identity by Genus",
    x = "Protein Identity (%)",
    y = "Relative Density"
  ) +
  theme_bw()

# Display the plot
recreated_plot


ggsave(output_plot_path, plot = recreated_plot, width = 10, height = 8, dpi = 300)

print(paste("Plot saved successfully to:", output_plot_path))


# --- 7. SAVE INDIVIDUAL PLOTS (with alpha highlighting) ---

# Get the base directory from your original output path
output_dir <- dirname(output_plot_path)

# Loop through each genus in your list
print("Starting to save individual plots as SVG...")
for (current_genus in genera_to_keep) {
  
  # 1. Filter data AND create the new highlighting column
  genus_data <- plot_data_binned %>%
    filter(Genus == current_genus) %>%
    mutate(ANI_highlight = ifelse(ANI_proxy < 95, "Below 95%", "95% and Above"))
  
  # 2. Create the plot object
  genus_plot <- ggplot(
    genus_data,
    aes(
      x = pident, 
      y = relative_density, 
      group = MAG_id, 
      color = ANI_factor,
      alpha = ANI_highlight 
    )
  ) +
    geom_line(linewidth = 0.7) + 
    
    scale_color_manual(
      values = gruvbox_palette_ramp,
      name = "Avg. Protein ID (%)",
      drop = FALSE 
    ) +
    
    scale_alpha_manual(
      name = "ANI Highlight",
      values = c("Below 95%" = 1.0, "95% and Above" = 0.25),
      guide = "none" 
    ) +
    
    scale_x_continuous(breaks = c(75, 80, 85, 90, 95, 100)) +
    xlim(75, 100) +
    
    labs(
      title = paste("Distribution of Best-Hit Protein Identity:", current_genus),
      subtitle = "Curves with <95% Avg. Protein ID are highlighted (opaque)",
      x = "Protein Identity (%)",
      y = "Relative Density"
    ) +
    theme_bw() +
    theme(legend.position = "right") +
    
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  
  # --- MODIFIED PART ---
  # 3. Define the new output filename with .svg extension
  output_filename <- file.path(output_dir, paste0("ani_distribution_highlighted_", current_genus, "."))
  
  # 4. Save the individual plot as SVG
  # We just changed the extension to .svg and removed the 'dpi' argument.
  ggsave(output_filename, plot = genus_plot, width = 7, height = 5)
  # --- END MODIFIED PART ---
  
  print(paste("Saved plot:", output_filename))
}

print("All individual highlighted plots saved successfully as SVG.")



snodgrassella_mags <- plot_data_binned %>%
  filter(Genus == "snodgrassella") %>%
  select(MAG_id,ANI_bin) %>%
  distinct()

# Print the resulting tibble of unique MAG IDs
print(snodgrassella_mags)



