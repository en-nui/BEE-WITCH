# ---------------------------
# Produce requested plots IN MEMORY (no saving)
# ---------------------------
library(dplyr); library(tidyr); library(vegan); library(ggplot2); library(ggrepel); library(ape)

# ---------- Helpers ----------
compute_pcoa_and_percent <- function(community_matrix) {
  bray <- vegdist(sqrt(community_matrix), method = "bray")
  p <- ape::pcoa(as.matrix(bray))
  rel <- NULL
  if (!is.null(p$values$Relative_eig)) rel <- p$values$Relative_eig
  else if (!is.null(p$values$Eigenvalues)) {
    ev <- p$values$Eigenvalues; evpos <- ev[ev > 0]; rel <- evpos / sum(evpos)
  } else rel <- rep(NA, nrow(community_matrix))
  list(pcoa = p, rel_var = rel, scores = p$vectors)
}

compute_centroids <- function(plot_data, group_col, xcol = "NMDS1", ycol = "NMDS2") {
  plot_data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(across(all_of(group_col))) %>%
    summarise(centroid_x = mean(.data[[xcol]], na.rm = TRUE),
              centroid_y = mean(.data[[ycol]], na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    rename(group = 1)
}

plot_ord_with_centroids <- function(plot_data, xcol = "NMDS1", ycol = "NMDS2",
                                    color_by = NULL, shape_by = NULL,
                                    centroid_group = NULL, centroid_type = "centroid",
                                    connect = TRUE,
                                    point_size = 2.5, centroid_size = 4.5, centroid_label_size = 3,
                                    rel_var = NULL, zoom_limits = NULL) {
  # basic checks
  if (!(xcol %in% names(plot_data)) || !(ycol %in% names(plot_data))) stop("plot_data must contain x/y columns")
  if (!is.null(color_by) && !(color_by %in% names(plot_data))) color_by <- NULL
  if (!is.null(shape_by) && !(shape_by %in% names(plot_data))) shape_by <- NULL
  # build mapping
  aes_map <- aes_string(x = xcol, y = ycol)
  if (!is.null(color_by)) aes_map <- modifyList(aes_map, aes_string(color = color_by))
  if (!is.null(shape_by)) aes_map <- modifyList(aes_map, aes_string(shape = shape_by))
  p <- ggplot(plot_data, aes_map) +
    geom_point(size = point_size, stroke = 0.8, alpha = 0.95) +
    theme_bw() + theme(legend.position = "right") +
    coord_equal()
  # centroids
  if (!is.null(centroid_group) && centroid_group %in% names(plot_data)) {
    cent_df <- compute_centroids(plot_data, centroid_group, xcol = xcol, ycol = ycol)
    # optional connect segments (use full original plot_data so lines reflect full set)
    if (isTRUE(connect)) {
      seg_df <- plot_data %>% filter(!is.na(.data[[centroid_group]])) %>% left_join(cent_df, by = c(setNames("group", centroid_group)))
      p <- p + geom_segment(data = seg_df, aes_string(x = xcol, y = ycol, xend = "centroid_x", yend = "centroid_y"),
                            inherit.aes = FALSE, color = "gray40", alpha = 0.35, size = 0.35)
    }
    p <- p + geom_point(data = cent_df, aes(x = centroid_x, y = centroid_y),
                        inherit.aes = FALSE, shape = 21, fill = "white", size = centroid_size, stroke = 1)
    p <- p + ggrepel::geom_text_repel(data = cent_df, aes(x = centroid_x, y = centroid_y, label = group),
                                      inherit.aes = FALSE, size = centroid_label_size)
  }
  # axis labels (PCoA percent var)
  if (!is.null(rel_var) && length(rel_var) >= 2) {
    p <- p + labs(x = paste0(xcol, " (", round(rel_var[1]*100,2), "%)"),
                  y = paste0(ycol, " (", round(rel_var[2]*100,2), "%)"))
  } else {
    p <- p + labs(x = xcol, y = ycol)
  }
  # zoom (optional), zoom_limits should be list(xlim = c(min,max), ylim = c(min,max))
  if (!is.null(zoom_limits) && is.list(zoom_limits)) {
    p <- p + coord_cartesian(xlim = zoom_limits$xlim, ylim = zoom_limits$ylim)
  }
  p
}

# ---------- Utility: build community matrix & metadata, returns list ----------
build_comm_and_meta <- function(df) {
  # df must include sample_name, MAG_id, mean
  if (!all(c("sample_name","MAG_id","mean") %in% names(df))) stop("Input missing sample_name/MAG_id/mean")
  comm_long <- df %>% select(sample_name, MAG_id, mean) %>% group_by(sample_name, MAG_id) %>%
    summarise(total_abundance = sum(mean, na.rm = TRUE), .groups = "drop")
  comm_wide <- comm_long %>% pivot_wider(names_from = MAG_id, values_from = total_abundance, values_fill = 0)
  if (!("sample_name" %in% names(comm_wide))) stop("sample_name missing in community tibble")
  sample_names <- comm_wide$sample_name
  comm_mat <- comm_wide %>% select(-sample_name) %>% as.data.frame()
  rownames(comm_mat) <- sample_names
  comm_mat <- as.matrix(comm_mat); comm_mat[is.na(comm_mat)] <- 0
  comm_mat <- comm_mat[rowSums(comm_mat) > 0, , drop = FALSE]
  # metadata aligned to remaining samples
  meta_raw <- df %>% select(sample_name, State, Month, Caste, Colony, Rep, Season, Phenotype, Apiary) %>% distinct(sample_name, .keep_all = TRUE)
  meta_aligned <- data.frame(sample_name = rownames(comm_mat)) %>% left_join(meta_raw, by = "sample_name")
  list(community = comm_mat, metadata = meta_aligned)
}

# ---------- Main: compute NMDS/PCoA and create requested plots ----------
make_requested_plots <- function(ecology_data) {
  results <- list()
  # prepare derived columns (safe)
  if (!("Season" %in% names(ecology_data))) {
    get_season <- function(month) {
      if (month %in% c("May","June")) 'Spring' else if (month %in% c("July","August","September")) 'Summer' else if (month %in% c("October","November")) 'Fall' else if (month %in% c("January","February")) 'Winter' else NA
    }
    ecology_data$Season <- sapply(ecology_data$Month, get_season)
  }
  if (!("Phenotype" %in% names(ecology_data))) ecology_data$Phenotype <- sapply(ecology_data$Season, function(x) ifelse(x == 'Winter','Winter Bee','Non-winter Bee'))
  if (!("Apiary" %in% names(ecology_data))) {
    get_apiary <- function(colony) {
      colnum <- suppressWarnings(as.numeric(as.character(colony)))
      if (!is.na(colnum) && colnum > 500) 'HB' else if (!is.na(colnum) && colnum < 104) 'PT' else 'OB'
    }
    ecology_data$Apiary <- sapply(ecology_data$Colony, get_apiary)
  }
  
  # ensure factors
  for (c in c("Season","Apiary","Caste","Month","State","Colony","Rep","Phenotype")) if (c %in% names(ecology_data)) ecology_data[[c]] <- as.factor(ecology_data[[c]])
  
  ##### 1) Workers only, type == "mmag" (no larvae) - NMDS plots colored by Colony and by Season
  df_workers_mmag <- ecology_data %>% filter(Caste == "Worker", type == "mmag")
  bc1 <- build_comm_and_meta(df_workers_mmag)
  if (nrow(bc1$community) == 0) stop("No samples after building community matrix for Workers mmag")
  # NMDS
  set.seed(123); nmds1 <- metaMDS(sqrt(bc1$community), distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
  nm_plotdata <- as.data.frame(scores(nmds1, display = "sites")); nm_plotdata$sample_name <- rownames(nm_plotdata)
  nm_plotdata <- left_join(nm_plotdata, bc1$metadata, by = "sample_name")
  # also PCoA for %var (optional)
  pcoa1 <- compute_pcoa_and_percent(bc1$community)
  pcoa_scores1 <- as.data.frame(pcoa1$scores); pcoa_scores1$sample_name <- rownames(pcoa_scores1); pcoa_scores1 <- left_join(pcoa_scores1, bc1$metadata, by = "sample_name")
  # plots (NMDS)
  results$workers_mmag <- list()
  results$workers_mmag$nmds <- list()
  results$workers_mmag$nmds$by_colony <- plot_ord_with_centroids(nm_plotdata, xcol = "NMDS1", ycol = "NMDS2", color_by = "Colony", shape_by = NULL, centroid_group = NULL, connect = FALSE)
  results$workers_mmag$nmds$by_season <- plot_ord_with_centroids(nm_plotdata, xcol = "NMDS1", ycol = "NMDS2", color_by = "Season", shape_by = NULL, centroid_group = NULL, connect = FALSE)
  # PCoA versions (with percent var labels)
  results$workers_mmag$pcoa <- list()
  if (all(c("PCoA1","PCoA2") %in% names(pcoa_scores1))) {
    results$workers_mmag$pcoa$by_colony <- plot_ord_with_centroids(pcoa_scores1, xcol = "PCoA1", ycol = "PCoA2", color_by = "Colony", shape_by = NULL, centroid_group = NULL, connect = FALSE, rel_var = pcoa1$rel_var)
    results$workers_mmag$pcoa$by_season <- plot_ord_with_centroids(pcoa_scores1, xcol = "PCoA1", ycol = "PCoA2", color_by = "Season", shape_by = NULL, centroid_group = NULL, connect = FALSE, rel_var = pcoa1$rel_var)
  }
  
  ##### 2) type == "mmag" (includes Larvae; i.e. don't filter Caste) - plot with zoom/transformation
  df_all_mmag <- ecology_data %>% filter(type == "mmag")
  bc2 <- build_comm_and_meta(df_all_mmag)
  if (nrow(bc2$community) == 0) stop("No samples for type mmag (all castes)")
  set.seed(123); nmds2 <- metaMDS(sqrt(bc2$community), distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
  nm2_pd <- as.data.frame(scores(nmds2, display = "sites")); nm2_pd$sample_name <- rownames(nm2_pd); nm2_pd <- left_join(nm2_pd, bc2$metadata, by = "sample_name")
  pcoa2 <- compute_pcoa_and_percent(bc2$community)
  pcoa_scores2 <- as.data.frame(pcoa2$scores); pcoa_scores2$sample_name <- rownames(pcoa_scores2); pcoa_scores2 <- left_join(pcoa_scores2, bc2$metadata, by = "sample_name")
  # compute zoom limits that cut off extreme 1% tails on each axis (visual only)
  xlims <- quantile(nm2_pd$NMDS1, probs = c(0.01, 0.99), na.rm = TRUE)
  ylims <- quantile(nm2_pd$NMDS2, probs = c(0.01, 0.99), na.rm = TRUE)
  zoom_limits <- list(xlim = as.numeric(xlims), ylim = as.numeric(ylims))
  results$all_mmag <- list()
  results$all_mmag$nmds <- list()
  results$all_mmag$nmds$by_colony_zoom <- plot_ord_with_centroids(nm2_pd, xcol = "NMDS1", ycol = "NMDS2", color_by = "Colony", centroid_group = NULL, connect = FALSE, zoom_limits = zoom_limits)
  results$all_mmag$nmds$by_season_zoom <- plot_ord_with_centroids(nm2_pd, xcol = "NMDS1", ycol = "NMDS2", color_by = "Season", centroid_group = NULL, connect = FALSE, zoom_limits = zoom_limits)
  # also PCoA zoomed
  if (all(c("PCoA1","PCoA2") %in% names(pcoa_scores2))) {
    # build zoom limits on PCoA too
    xlims_p <- quantile(pcoa_scores2$PCoA1, probs = c(0.01, 0.99), na.rm = TRUE)
    ylims_p <- quantile(pcoa_scores2$PCoA2, probs = c(0.01, 0.99), na.rm = TRUE)
    zoom_limits_p <- list(xlim = as.numeric(xlims_p), ylim = as.numeric(ylims_p))
    results$all_mmag$pcoa <- list()
    results$all_mmag$pcoa$by_colony_zoom <- plot_ord_with_centroids(pcoa_scores2, xcol = "PCoA1", ycol = "PCoA2", color_by = "Colony", rel_var = pcoa2$rel_var, zoom_limits = zoom_limits_p)
    results$all_mmag$pcoa$by_season_zoom <- plot_ord_with_centroids(pcoa_scores2, xcol = "PCoA1", ycol = "PCoA2", color_by = "Season", rel_var = pcoa2$rel_var, zoom_limits = zoom_limits_p)
  }
  
  ##### 3) Centroids: for mmag (Workers-only OR full mmag choose which you want)
  # We'll produce centroids that connect samples -> Apiary and samples -> Season for Workers mmag (your original priority).
  # Use nm_plotdata (Workers mmag NMDS) as base; if you prefer 'nm2_pd' (all mmag including larvae) swap variable.
  base_for_centroids <- nm_plotdata  # change to nm2_pd if you want centroids for full mmag (includes larvae)
  results$centroids <- list()
  results$centroids$workers_mmag <- list()
  # connect to Apiary centroids
  if ("Apiary" %in% names(base_for_centroids)) {
    results$centroids$workers_mmag$connect_apiary <- plot_ord_with_centroids(base_for_centroids, xcol = "NMDS1", ycol = "NMDS2",
                                                                             color_by = "Colony", shape_by = "Season",
                                                                             centroid_group = "Apiary", connect = TRUE)
  }
  
  
  if ("Colony" %in% names(base_for_centroids)) {
    results$centroids$workers_mmag$connect_colony <- plot_ord_with_centroids(base_for_centroids, xcol = "NMDS1", ycol = "NMDS2",
                                                                             color_by = "Colony", shape_by = "Season",
                                                                             centroid_group = "Colony", connect = TRUE)
  }
  # connect to Season centroids
  if ("Season" %in% names(base_for_centroids)) {
    results$centroids$workers_mmag$connect_season <- plot_ord_with_centroids(base_for_centroids, xcol = "NMDS1", ycol = "NMDS2",
                                                                             color_by = "Colony", shape_by = "Apiary",
                                                                             centroid_group = "Season", connect = TRUE)
  }
  
  # store some metadata about ordinations
  results$info <- list(nmds_workers_mmag = nmds1$stress, nmds_all_mmag = nmds2$stress, pcoa_workers_rel = pcoa1$rel_var, pcoa_all_rel = pcoa2$rel_var)
  return(results)
}

# ---- Run (do not save) ----
plots_out <- make_requested_plots(ecology_data)

# ---- How to inspect / print ----
# Example: Workers-only NMDS colored by Colony
print(plots_out$workers_mmag$nmds$by_colony)

# Example: Workers-only NMDS colored by Season
print(plots_out$workers_mmag$nmds$by_season)

# Example: All mmag (includes larvae) NMDS zoomed (colored by Season)
print(plots_out$all_mmag$nmds)

# Example: Centroid connections (workers mmag) -> Apiary
if (!is.null(plots_out$centroids$workers_mmag$connect_apiary)) print(plots_out$centroids$workers_mmag$connect_apiary)
if (!is.null(plots_out$centroids$workers_mmag$connect_season)) print(plots_out$centroids$workers_mmag$connect_season)
if (!is.null(plots_out$centroids$workers_mmag$connect_season)) print(plots_out$centroids$workers_mmag$connect_colony)

# You can also toggle to use PCoA plots (with percent-variance labels), e.g.:
if (!is.null(plots_out$workers_mmag$pcoa$by_colony)) print(plots_out$workers_mmag$pcoa$by_colony)
pc <- compute_pcoa_and_percent(community_matrix)   # returns list(pcoa, rel_var, scores)
rel_var <- pc$rel_var
# rel_var[1], rel_var[2] are relative explained variance for first two axes (fraction)
x_pct <- round(rel_var[1] * 100, 2)
y_pct <- round(rel_var[2] * 100, 2)










# assume nmds_result from metaMDS and plot_data has NMDS1/NMDS2
stress_value <- nmds_result$stress
p_nmds <- plot_ord_with_centroids(plot_data, xcol="NMDS1", ycol="NMDS2", color_by="Season", shape_by="Apiary")
p_nmds <- p_nmds + labs(subtitle = paste0("NMDS stress: ", formatC(stress_value, digits=3)))
print(p_nmds)
