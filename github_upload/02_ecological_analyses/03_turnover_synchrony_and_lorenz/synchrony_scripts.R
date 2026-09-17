ecology_data <- read.csv('ecology_metadata_df.csv')
mmag_worker
#keep only worker mmags
mmag_worker <- ecology_data %>%
  filter(type == "mmag", Caste == "Worker")

# 0B) collapse contig-level rows -> sample × MAG
# choose "first(mean)" because it's genome-wide and duplicated; could also use mean(mean) for safety
mag_by_sample <- mmag_worker %>%
  group_by(sample_name, Colony, Month, Time, Rep, MAG_id, Genus, Season) %>%
  summarise(
    mag_mean = first(mean),
    .groups = "drop"
  )

# sanity check: should now be 1 row per sample_name × MAG_id
mag_by_sample %>%
  count(sample_name, MAG_id) %>%
  summarise(max_n = max(n))
head(mag_by_sample)
## relabund
mag_rel_by_sample <- mag_by_sample %>%
  group_by(sample_name) %>%
  mutate(
    total_microbial = sum(mag_mean, na.rm = TRUE),
    rel_abund = if_else(total_microbial > 0, mag_mean / total_microbial, NA_real_)
  ) %>%
  ungroup()


replicate_noise <- mag_rel_by_sample %>%
  group_by(Colony, Month, Time, MAG_id) %>%
  summarise(
    n_rep = n_distinct(sample_name),
    mean_rel = mean(rel_abund, na.rm = TRUE),
    var_rep  = var(rel_abund, na.rm = TRUE),
    sd_rep   = sd(rel_abund, na.rm = TRUE),
    cv_rep   = sd_rep / (mean_rel + 1e-12),
    .groups = "drop"
  ) %>%
  filter(n_rep >= 2)


ggplot(data=replicate_noise,aes(x=mean_rel, y=cv_rep))+
  geom_point()


# 3A) replicate means per month
mag_monthly <- mag_rel_by_sample %>%
  group_by(Colony, Month, Time, MAG_id) %>%
  summarise(
    n_rep = n_distinct(sample_name),
    x_bar = mean(rel_abund, na.rm = TRUE),
    s2_rep = var(rel_abund, na.rm = TRUE),
    .groups = "drop"
  )

# 3B) variance decomposition per colony×MAG
mag_noise_partition <- mag_monthly %>%
  group_by(Colony, MAG_id) %>%
  filter(n_distinct(Time) >= 3) %>%
  summarise(
    T = n_distinct(Time),
    var_between = var(x_bar, na.rm = TRUE),
    mean_noise_in_mean = mean(s2_rep / pmax(n_rep,1), na.rm = TRUE),
    proc_var = pmax(var_between - mean_noise_in_mean, 0),
    noise_fraction = if_else(var_between > 0, mean_noise_in_mean / var_between, NA_real_),
    .groups = "drop"
  )



eps <- 1e-8

mag_changes <- mag_monthly %>%
  arrange(Colony, MAG_id, Time) %>%
  group_by(Colony, MAG_id) %>%
  mutate(
    mu = log((x_bar + eps) / (lag(x_bar) + eps)),
    dt = Time - lag(Time)
  ) %>%
  ungroup() %>%
  filter(!is.na(mu), dt > 0)

# Quick Laplace MLE for b (assuming symmetric around 0): b_hat = mean(|mu|)
b_hat <- mean(abs(mag_changes$mu), na.rm = TRUE)

# Compare Laplace vs Normal log-likelihood (simple AIC-style check)
ll_laplace <- sum(-log(2*b_hat) - abs(mag_changes$mu)/b_hat, na.rm = TRUE)

sd_hat <- sd(mag_changes$mu, na.rm = TRUE)
ll_normal <- sum(dnorm(mag_changes$mu, mean=0, sd=sd_hat, log=TRUE), na.rm = TRUE)

c(ll_laplace=ll_laplace, ll_normal=ll_normal, b_hat=b_hat, sd_hat=sd_hat)





taylor_df <- mag_monthly %>%
  group_by(Colony, MAG_id) %>%
  filter(n_distinct(Time) >= 4) %>%
  summarise(
    mean_x = mean(x_bar, na.rm=TRUE),
    var_x  = var(x_bar, na.rm=TRUE),
    .groups="drop"
  ) %>%
  filter(mean_x > 0, var_x > 0)

taylor_fit <- lm(log(var_x) ~ log(mean_x), data=taylor_df)
summary(taylor_fit)
coef(taylor_fit)["log(mean_x)"]  # this is beta_hat



taylor_df2 <- mag_monthly %>%
  group_by(Colony, MAG_id) %>%
  filter(n_distinct(Time) >= 4) %>%
  summarise(
    mean_x = mean(x_bar, na.rm=TRUE),
    var_x_raw = var(x_bar, na.rm=TRUE),
    mean_noise_in_mean = mean(s2_rep / pmax(n_rep,1), na.rm=TRUE),
    var_x = pmax(var_x_raw - mean_noise_in_mean, 0),
    .groups="drop"
  ) %>%
  filter(mean_x > 0, var_x > 0)

taylor_fit2 <- lm(log(var_x) ~ log(mean_x), data=taylor_df2)
summary(taylor_fit2)




estimate_tau <- function(df_one, eps=1e-8){
  df_one <- df_one %>% arrange(Time)
  x <- df_one$x_bar
  t <- df_one$Time
  if (length(x) < 4) return(NULL)
  
  xbar <- mean(x, na.rm=TRUE)
  x0 <- x[-length(x)] - xbar
  x1 <- x[-1]         - xbar
  dt <- t[-1] - t[-length(t)]
  
  ok <- is.finite(x0) & is.finite(x1) & is.finite(dt) & dt > 0
  x0 <- x0[ok]; x1 <- x1[ok]; dt <- dt[ok]
  if (length(x0) < 3) return(NULL)
  
  # slope through origin: phi = sum(x0*x1)/sum(x0^2)
  phi <- sum(x0 * x1) / sum(x0^2)
  # only meaningful if 0<phi<1
  if (!is.finite(phi) || phi <= 0 || phi >= 1) return(tibble(phi=phi, tau=NA_real_))
  
  dt_med <- median(dt)
  tau <- -dt_med / log(phi)
  tibble(phi=phi, tau=tau)
}

tau_df <- mag_monthly %>%
  group_by(Colony, MAG_id) %>%
  group_modify(~{
    out <- estimate_tau(.x)
    if (is.null(out)) return(tibble())
    out
  }) %>%
  ungroup()
tau_df

median(mag_noise_partition$noise_fraction, na.rm = TRUE)
IQR(mag_noise_partition$noise_fraction, na.rm = TRUE)





sigma_df <- mag_monthly %>%
  group_by(Colony, MAG_id) %>%
  filter(n_distinct(Time) >= 4) %>%   # use >=4 for more stable moments
  summarise(
    mu = mean(x_bar, na.rm=TRUE),
    var_raw = var(x_bar, na.rm=TRUE),
    noise_mean = mean(s2_rep / pmax(n_rep,1), na.rm=TRUE),
    var_proc = pmax(var_raw - noise_mean, 0),
    cv2_proc = var_proc / (mu^2 + 1e-20),
    sigma_i  = (2*cv2_proc) / (1 + cv2_proc),
    .groups="drop"
  ) %>%
  filter(is.finite(sigma_i), mu > 0)

sigma_global

sigma_global <- median(sigma_df$sigma_i, na.rm=TRUE)
sigma_by_colony <- sigma_df %>%
  group_by(Colony) %>%
  summarise(sigma_med = median(sigma_i, na.rm=TRUE),
            sigma_iqr = IQR(sigma_i, na.rm=TRUE),
            n = n(), .groups="drop")




K_df <- sigma_df %>%
  mutate(
    K_i = mu / (1 - sigma_global/2)
  )


# A) noise dominance fractions
mean(mag_noise_partition$noise_fraction > 0.5, na.rm=TRUE)
mean(mag_noise_partition$noise_fraction < 0.2, na.rm=TRUE)

# B) sigma estimates (global median + IQR)
sigma_global <- median(sigma_df$sigma_i, na.rm=TRUE)
sigma_iqr <- IQR(sigma_df$sigma_i, na.rm=TRUE)
c(sigma_global=sigma_global, sigma_iqr=sigma_iqr)


sigma_df_filt <- sigma_df %>%
  filter(mu > 1e-4)  # tune as needed

sigma_global_filt <- median(sigma_df_filt$sigma_i, na.rm=TRUE)
IQR(sigma_df_filt$sigma_i, na.rm=TRUE)




summary_sigma <- sigma_df %>%
  mutate(rare = mu < 1e-4) %>%
  summarise(
    n_total = n(),
    n_rare  = sum(rare),
    sigma_med_all = median(sigma_i, na.rm=TRUE),
    sigma_iqr_all = IQR(sigma_i, na.rm=TRUE),
    sigma_med_nonrare = median(sigma_i[!rare], na.rm=TRUE),
    sigma_iqr_nonrare = IQR(sigma_i[!rare], na.rm=TRUE)
  )
summary_sigma


thr <- quantile(sigma_df_filt$sigma_i, 0.25, na.rm=TRUE)
stable_mag <- sigma_df_filt %>% mutate(stable = sigma_i <= thr)
mean(stable_mag$stable)  # should be 0.25 by construction globally













#main take away is that variation in relative abundance, only 16% of variance can be explained by sampling error 



# 1) genus-level abundances per sample (replicate-level)
genus_long <- mag_by_sample %>%
  group_by(sample_name, Colony, Time, Month, Season, Rep, Genus) %>%
  summarise(abund = sum(mag_mean, na.rm=TRUE), .groups="drop")

# 2) mag-level abundances per sample (replicate-level)
mag_long <- mag_by_sample %>%
  group_by(sample_name, Colony, Time, Month, Season, Rep, MAG_id) %>%
  summarise(abund = sum(mag_mean, na.rm=TRUE), .groups="drop")

# 3) collapse replicates within a colony-timepoint if you want a single trajectory
genus_long_ct <- genus_long %>%
  group_by(Colony, Time, Month, Season, Genus) %>%
  summarise(abund = mean(abund, na.rm=TRUE), .groups="drop")



make_mat <- function(dat, feature_col){
  dat %>%
    tidyr::pivot_wider(names_from = all_of(feature_col), values_from = abund, values_fill = 0) %>%
    tibble::column_to_rownames("sample_id")
}

# Create a sample_id column before pivot
genus_wide <- genus_long %>%
  mutate(sample_id = sample_name) %>%
  select(sample_id, Genus, abund) %>%
  pivot_wider(names_from=Genus, values_from=abund, values_fill=0)

# for colony-timepoint aggregated
genus_wide_ct <- genus_long_ct %>%
  mutate(sample_id = paste0("C", Colony, "_T", Time)) %>%
  select(sample_id, Genus, abund) %>%
  pivot_wider(names_from=Genus, values_from=abund, values_fill=0)




# metadata for the aggregated matrix
meta_ct <- genus_long_ct %>%
  distinct(Colony, Time, Month, Season) %>%
  mutate(sample_id = paste0("C", Colony, "_T", Time))

X <- genus_wide_ct %>% tibble::column_to_rownames("sample_id") %>% as.matrix()

bray <- vegdist(X, method="bray")
pcoa <- cmdscale(bray, k=2, eig=TRUE)

scores <- data.frame(sample_id = rownames(pcoa$points),
                     PC1 = pcoa$points[,1],
                     PC2 = pcoa$points[,2]) %>%
  left_join(meta_ct, by="sample_id") %>%
  arrange(Colony, Time)

# colony trajectories: correlation of PC1 among colonies (requires aligning timepoints)
pc1_wide <- scores %>%
  select(Colony, Time, PC1) %>%
  pivot_wider(names_from=Colony, values_from=PC1)

# pairwise colony synchrony (cor across aligned timepoints)
cors <- cor(pc1_wide %>% select(-Time), use="pairwise.complete.obs", method="spearman")

cors


print(genus_sync,n=41)
eps <- 1e-6

genus_anom <- genus_long_ct %>%
  mutate(y = log(abund + eps)) %>%
  group_by(Colony, Genus) %>%
  mutate(z = y - mean(y, na.rm=TRUE)) %>%
  ungroup()

# wide by colony for each genus; compute mean pairwise cor
genus_sync <- genus_anom %>%
  select(Genus, Colony, Time, z) %>%
  group_by(Genus) %>%
  group_modify(~{
    # 1. Pivot the data
    tmp <- .x %>% 
      select(Colony, Time, z) %>% 
      pivot_wider(names_from = Colony, values_from = z)
    
    # 2. Extract matrix (excluding Time)
    M <- as.matrix(tmp %>% select(-Time))
    
    # 3. SAFETY CHECK: Do we have at least 2 columns (Colonies) and some data?
    if(ncol(M) < 2) {
      return(data.frame(mean_pairwise_r = NA_real_))
    }
    
    # 4. Compute correlation
    C <- cor(M, use = "pairwise.complete.obs", method = "spearman")
    
    # 5. Handle cases where C might be all NA
    if(all(is.na(C))) {
      return(data.frame(mean_pairwise_r = NA_real_))
    }
    
    mean_cor <- mean(C[upper.tri(C)], na.rm = TRUE)
    data.frame(mean_pairwise_r = mean_cor)
  }) %>% 
  ungroup()





genus_sync2 <- genus_anom %>%
  select(Genus, Colony, Time, z) %>%
  group_by(Genus) %>%
  group_modify(~{
    tmp <- .x %>% pivot_wider(names_from=Colony, values_from=z)
    M <- as.matrix(tmp %>% select(-Time))
    
    # colonies that have enough non-NA points and nonzero variance
    ok <- apply(M, 2, function(v){
      sum(!is.na(v)) >= 6 && sd(v, na.rm=TRUE) > 0
    })
    M2 <- M[, ok, drop=FALSE]
    
    if (ncol(M2) < 3) {
      return(data.frame(mean_pairwise_r = NA_real_,
                        n_colonies_used = ncol(M2),
                        n_timepoints_overlap = sum(complete.cases(M2))))
    }
    
    C <- cor(M2, use="pairwise.complete.obs", method="spearman")
    mean_cor <- mean(C[upper.tri(C)], na.rm=TRUE)
    
    data.frame(mean_pairwise_r = mean_cor,
               n_colonies_used = ncol(M2),
               n_timepoints_overlap = sum(complete.cases(M2)))
  }) %>% ungroup() %>%
  arrange(desc(mean_pairwise_r))
genus_sync2

#alighned pc1 matrix (time by colony)
# scores must have: Colony, Time, PC1
pc1_wide <- scores %>%
  select(Colony, Time, PC1) %>%
  mutate(Colony = as.character(Colony)) %>%
  pivot_wider(names_from = Colony, values_from = PC1) %>%
  arrange(Time)

M <- as.matrix(pc1_wide %>% select(-Time))  # rows=timepoints, cols=colonies


#colony-colony correlation matrix + per-colony mean correlation

cors_obs <- cor(M, use="pairwise.complete.obs", method="spearman")

# per-colony mean correlation to all other colonies (exclude self)
mean_cor_per_colony <- sapply(seq_len(ncol(cors_obs)), function(i){
  mean(cors_obs[i, -i], na.rm=TRUE)
})
mean_cor_per_colony <- sort(mean_cor_per_colony, decreasing=TRUE)

# overall synchrony statistic (mean off-diagonal correlation)
overall_sync_obs <- mean(cors_obs[upper.tri(cors_obs)], na.rm=TRUE)

overall_sync_obs
mean_cor_per_colony



#permutation (shuffling time within colony)

set.seed(1)

perm_sync <- function(M){
  Mperm <- M
  for (j in seq_len(ncol(Mperm))){
    idx <- which(!is.na(Mperm[, j]))
    Mperm[idx, j] <- sample(Mperm[idx, j], size=length(idx), replace=FALSE)
  }
  C <- cor(Mperm, use="pairwise.complete.obs", method="spearman")
  mean(C[upper.tri(C)], na.rm=TRUE)
}

B <- 5000
sync_null <- replicate(B, perm_sync(M))

p_val <- (sum(sync_null >= overall_sync_obs) + 1) / (B + 1)

list(overall_sync_obs = overall_sync_obs,
     p_value = p_val,
     null_mean = mean(sync_null),
     null_sd = sd(sync_null))







col_meta <- tibble::tibble(
  Colony = c("555","777","888","999",
             "104","105","106","107",
             "100","101","102","103"),
  Apiary = c(rep("HB", 4),
             rep("OB", 4),
             rep("PT", 4))
)

# sanity check alignment
stopifnot(is.matrix(cors_obs), !is.null(colnames(cors_obs)))


cols <- colnames(cors_obs)

cors_df <- expand.grid(Col1 = cols, Col2 = cols, stringsAsFactors = FALSE) %>%
  mutate(r = cors_obs[cbind(match(Col1, cols), match(Col2, cols))]) %>%
  # keep each pair once (upper triangle excluding diagonal)
  filter(match(Col1, cols) < match(Col2, cols)) %>%
  left_join(col_meta, by = c("Col1" = "Colony")) %>%
  rename(Apiary1 = Apiary) %>%
  left_join(col_meta, by = c("Col2" = "Colony")) %>%
  rename(Apiary2 = Apiary) %>%
  mutate(within_apiary = (Apiary1 == Apiary2))

obs_within_mean  <- mean(cors_df$r[cors_df$within_apiary], na.rm=TRUE)
obs_between_mean <- mean(cors_df$r[!cors_df$within_apiary], na.rm=TRUE)
obs_delta <- obs_within_mean - obs_between_mean

c(within = obs_within_mean, between = obs_between_mean, delta = obs_delta)
table(cors_df$within_apiary, useNA="ifany")
summary(cors_df$r)




set.seed(2)
B <- 5000

apiary_perm_delta <- function(cors_df, col_meta){
  perm_meta <- col_meta
  perm_meta$Apiary <- sample(perm_meta$Apiary, replace=FALSE)
  
  tmp <- cors_df %>%
    select(Col1, Col2, r) %>%
    left_join(perm_meta, by = c("Col1" = "Colony")) %>% rename(Apiary1 = Apiary) %>%
    left_join(perm_meta, by = c("Col2" = "Colony")) %>% rename(Apiary2 = Apiary) %>%
    mutate(within_apiary = (Apiary1 == Apiary2))
  
  w <- mean(tmp$r[tmp$within_apiary], na.rm=TRUE)
  b <- mean(tmp$r[!tmp$within_apiary], na.rm=TRUE)
  w - b
}

delta_null <- replicate(B, apiary_perm_delta(cors_df, col_meta))
p_apiary <- (sum(delta_null >= obs_delta) + 1) / (B + 1)

list(obs_delta = obs_delta,
     p_value = p_apiary,
     null_mean = mean(delta_null),
     null_sd = sd(delta_null))


cors_df %>%
  mutate(group = ifelse(within_apiary, "Within apiary", "Between apiaries")) %>%
  ggplot(aes(x = r, fill = group)) +
  geom_histogram(bins = 15, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = obs_within_mean, linetype = 2) +
  geom_vline(xintercept = obs_between_mean, linetype = 2) +
  labs(x = "Colony–colony synchrony (Spearman corr of PC1)",
       y = "Number of colony pairs",
       title = "Are colonies within the same apiary more synchronized?")







min_tp <- 5      # minimum timepoints per colony-genus series
min_cols <- 6    # minimum colonies contributing per genus

# 1) keep only genus-colony series with enough data + nonzero variance
dat_ok <- genus_anom %>%
  group_by(Colony, Genus) %>%
  filter(sum(!is.na(z)) >= min_tp, sd(z, na.rm=TRUE) > 0) %>%
  ungroup()

# 2) keep genera observed in enough colonies
keep_genera <- dat_ok %>%
  distinct(Genus, Colony) %>%
  count(Genus, name="n_colonies") %>%
  filter(n_colonies >= min_cols) %>%
  pull(Genus)

dat_ok <- dat_ok %>% filter(Genus %in% keep_genera)

# 3) build a Time x Colony matrix per genus, compute leave-one-out consensus correlations
calc_r_for_genus <- function(g){
  tmp <- dat_ok %>% filter(Genus == g) %>%
    select(Time, Colony, z) %>%
    pivot_wider(names_from=Colony, values_from=z) %>%
    arrange(Time)
  
  M <- as.matrix(tmp %>% select(-Time))
  cols <- colnames(M)
  
  # compute leave-one-out consensus for each colony (rowwise mean excluding that column)
  r <- sapply(seq_along(cols), function(j){
    x <- M[, j]
    Y <- M[, -j, drop=FALSE]
    cons <- rowMeans(Y, na.rm=TRUE)
    if (sum(!is.na(x) & !is.na(cons)) < min_tp) return(NA_real_)
    suppressWarnings(cor(x, cons, use="pairwise.complete.obs", method="spearman"))
  })
  
  tibble(Genus = g, Colony = cols, r = as.numeric(r))
}

r_cg <- map_dfr(keep_genera, calc_r_for_genus)

# 4) make a genus x colony matrix for heatmap
Rmat <- r_cg %>%
  pivot_wider(names_from=Colony, values_from=r) %>%
  tibble::column_to_rownames("Genus") %>%
  as.matrix()

# 5) annotate columns by apiary
ann_col <- col_meta %>%
  mutate(Colony = as.character(Colony)) %>%
  filter(Colony %in% colnames(Rmat)) %>%
  tibble::column_to_rownames("Colony")

# reorder columns to match annotation rownames
Rmat <- Rmat[, rownames(ann_col), drop=FALSE]

pheatmap(Rmat,
         annotation_col = ann_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Genus seasonal synchrony by colony (leave-one-out consensus r)",
         na_col = "grey90")




outliers <- r_cg %>%
  group_by(Genus) %>%
  mutate(q1 = quantile(r, 0.25, na.rm=TRUE),
         q3 = quantile(r, 0.75, na.rm=TRUE),raki
         iqr = q3 - q1,
         is_outlier = r < (q1 - 1.5*iqr)) %>%
  filter(is_outlier) %>%
  arrange(Genus, r)

outliers

















library(mgcv)


# scores must include: Colony, Time, PC1
scores2 <- scores %>%
  mutate(Colony = as.character(Colony)) %>%
  left_join(col_meta, by=c("Colony"="Colony")) %>%
  mutate(Apiary = factor(Apiary),
         Colony = factor(Colony))

# hierarchical GAM: global smooth + random apiary + random colony
fit <- gam(PC1 ~ s(Time, k=6) + s(Apiary, bs="re") + s(Colony, bs="re"),
           data = scores2, method="REML")

scores2 <- scores2 %>%
  mutate(fitted = predict(fit, newdata = scores2),
         resid  = PC1 - fitted)

# helper: build correlation matrix from a simulated dataset
corr_from_scores <- function(df_sim){
  pc1_wide <- df_sim %>%
    select(Colony, Time, PC1_sim) %>%
    pivot_wider(names_from=Colony, values_from=PC1_sim) %>%
    arrange(Time)
  M <- as.matrix(pc1_wide %>% select(-Time))
  cor(M, use="pairwise.complete.obs", method="spearman")
}

# helper: apiary delta from correlation matrix
apiary_delta_from_C <- function(C){
  cols <- colnames(C)
  cors_df <- expand.grid(Col1 = cols, Col2 = cols, stringsAsFactors=FALSE) %>%
    mutate(r = C[cbind(match(Col1, cols), match(Col2, cols))]) %>%
    filter(match(Col1, cols) < match(Col2, cols)) %>%
    left_join(col_meta, by=c("Col1"="Colony")) %>% rename(Apiary1=Apiary) %>%
    left_join(col_meta, by=c("Col2"="Colony")) %>% rename(Apiary2=Apiary) %>%
    mutate(within = (Apiary1 == Apiary2))
  
  mean(cors_df$r[cors_df$within], na.rm=TRUE) - mean(cors_df$r[!cors_df$within], na.rm=TRUE)
}

set.seed(123)
B <- 2000

sim_stats <- replicate(B, {
  # resample residuals within colony (preserves colony variance)
  df_sim <- scores2 %>%
    group_by(Colony) %>%
    mutate(resid_sim = sample(resid, size=n(), replace=TRUE)) %>%
    ungroup() %>%
    mutate(PC1_sim = fitted + resid_sim)
  
  C <- corr_from_scores(df_sim)
  c(overall_sync = mean(C[upper.tri(C)], na.rm=TRUE),
    apiary_delta = apiary_delta_from_C(C))
})

sim_stats <- t(sim_stats)




obs_overall <- overall_sync_obs
obs_delta <- obs_delta  # from your apiary test

p_overall_sim <- (sum(sim_stats[, "overall_sync"] >= obs_overall) + 1) / (B + 1)
p_delta_sim   <- (sum(sim_stats[, "apiary_delta"] >= obs_delta) + 1) / (B + 1)

list(p_overall_sim = p_overall_sim,
     p_delta_sim = p_delta_sim,
     sim_overall_mean = mean(sim_stats[, "overall_sync"]),
     sim_delta_mean = mean(sim_stats[, "apiary_delta"]))




dfp <- as.data.frame(sim_stats)

ggplot(dfp, aes(x=overall_sync)) +
  geom_histogram(bins=40) +
  geom_vline(xintercept=obs_overall, linetype=2) +
  labs(title="Overall synchrony: observed vs simulated null",
       x="Mean off-diagonal Spearman r", y="Count")

ggplot(dfp, aes(x=apiary_delta)) +
  geom_histogram(bins=40) +
  geom_vline(xintercept=obs_delta, linetype=2) +
  labs(title="Apiary structure: observed vs simulated null",
       x="Δ (within-apiary mean r − between-apiary mean r)", y="Count")















outlier_by_colony <- outliers %>%
  count(Colony, name="n_genus_outliers") %>%
  left_join(col_meta, by="Colony") %>%
  arrange(desc(n_genus_outliers))

outlier_by_colony

genus_summary <- r_cg %>%
  group_by(Genus) %>%
  summarise(mean_r = mean(r, na.rm=TRUE),
            frac_negative = mean(r < 0, na.rm=TRUE),
            n = sum(!is.na(r)),
            .groups="drop") %>%
  arrange(desc(mean_r))

genus_summary






shift_non_na <- function(v){
  idx <- which(!is.na(v))
  x <- v[idx]
  n <- length(x)
  if (n <= 1) return(v)
  k <- sample(0:(n-1), 1)
  x2 <- if (k == 0) x else c(tail(x, k), head(x, n-k))
  v[idx] <- x2
  v
}

perm_sync_shift <- function(M){
  Mperm <- M
  for (j in seq_len(ncol(Mperm))){
    Mperm[, j] <- shift_non_na(Mperm[, j])
  }
  C <- cor(Mperm, use="pairwise.complete.obs", method="spearman")
  mean(C[upper.tri(C)], na.rm=TRUE)
}

set.seed(10)
B <- 5000
sync_null_shift <- replicate(B, perm_sync_shift(M))
p_shift <- (sum(sync_null_shift >= overall_sync_obs) + 1) / (B + 1)

list(overall_sync_obs = overall_sync_obs,
     p_shift = p_shift,
     null_mean = mean(sync_null_shift),
     null_sd = sd(sync_null_shift))






fit2 <- gam(PC1 ~ s(Time, k=6) + s(Time, Apiary, bs="fs", k=6) + s(Colony, bs="re"),
            data = scores2, method="REML")

scores2 <- scores2 %>%
  mutate(fitted2 = predict(fit2, newdata = scores2),
         resid2  = PC1 - fitted2)




set.seed(11)
B <- 2000

sim_stats2 <- replicate(B, {
  df_sim <- scores2 %>%
    group_by(Colony) %>%
    mutate(resid_sim = shift_non_na(resid2)) %>%
    ungroup() %>%
    mutate(PC1_sim = fitted2 + resid_sim)
  
  C <- corr_from_scores(df_sim)
  c(overall_sync = mean(C[upper.tri(C)], na.rm=TRUE),
    apiary_delta = apiary_delta_from_C(C))
})

sim_stats2 <- t(sim_stats2)

p_overall_sim2 <- (sum(sim_stats2[, "overall_sync"] >= overall_sync_obs) + 1) / (B + 1)
p_delta_sim2   <- (sum(sim_stats2[, "apiary_delta"] >= obs_delta) + 1) / (B + 1)

list(p_overall_sim2 = p_overall_sim2,
     p_delta_sim2 = p_delta_sim2,
     sim_overall_mean2 = mean(sim_stats2[, "overall_sync"]),
     sim_delta_mean2 = mean(sim_stats2[, "apiary_delta"]))




AIC(fit, fit2)
anova(fit, fit2, test="F")




ph <- pheatmap(cors_obs,
               clustering_distance_rows = as.dist(sqrt(2*(1 - cors_obs))),
               clustering_distance_cols = as.dist(sqrt(2*(1 - cors_obs))),
               clustering_method = "average",
               silent = TRUE)

col_order <- ph$tree_col$labels[ph$tree_col$order]
Rmat2 <- Rmat[, col_order, drop = FALSE]

ann_col <- col_meta %>%
  mutate(Colony = as.character(Colony)) %>%
  filter(Colony %in% colnames(Rmat2)) %>%
  tibble::column_to_rownames("Colony")

# outlier counts per colony (optional but nice)
outlier_counts <- outliers %>%
  count(Colony, name="n_outliers") %>%
  mutate(Colony = as.character(Colony))

ann_col <- ann_col %>%
  left_join(outlier_counts, by = c("row.names" = "Colony")) %>%
  tibble::column_to_rownames("row.names")

pheatmap(Rmat2,
         annotation_col = ann_col[, "Apiary", drop=FALSE],
         cluster_rows = TRUE,
         cluster_cols = FALSE,   # keep same colony order as colony heatmap
         main = "Genus seasonal synchrony by colony (leave-one-out Spearman r)",
         na_col = "grey90")





# scores2 has Colony, Time, PC1, Apiary
# global mean trajectory
mu_t <- scores2 %>%
  group_by(Time) %>%
  summarise(mu = mean(PC1), .groups="drop")

# residuals around global trajectory
scores3 <- scores2 %>%
  left_join(mu_t, by="Time") %>%
  mutate(r0 = PC1 - mu)

# apiary-shared residual at each time: mean residual within apiary
A_at <- scores3 %>%
  group_by(Apiary, Time) %>%
  summarise(A = mean(r0), .groups="drop")

# colony-specific residual = remainder
scores3 <- scores3 %>%
  left_join(A_at, by=c("Apiary","Time")) %>%
  mutate(eps = r0 - A)

var_A  <- var(A_at$A)
var_eps <- var(scores3$eps)
c(var_apiary_shared = var_A, var_colony_noise = var_eps, ratio = var_A/var_eps)


shift_non_na <- function(v){
  idx <- which(!is.na(v))
  x <- v[idx]; n <- length(x)
  if (n <= 1) return(v)
  k <- sample(0:(n-1), 1)
  x2 <- if (k == 0) x else c(tail(x, k), head(x, n-k))
  v[idx] <- x2
  v
}

simulate_PC1 <- function(include_apiary=FALSE){
  # simulate eps by circularly shifting each colony’s eps series (keeps structure)
  sim <- scores3 %>%
    group_by(Colony) %>%
    mutate(eps_sim = shift_non_na(eps)) %>%
    ungroup()
  
  if (include_apiary){
    # also circularly shift apiary-shared component as a block per apiary
    A_sim <- A_at %>%
      group_by(Apiary) %>%
      mutate(A_shift = shift_non_na(A)) %>%
      ungroup()
    
    sim <- sim %>%
      left_join(A_sim %>% select(Apiary, Time, A_shift), by=c("Apiary","Time")) %>%
      mutate(PC1_sim = mu + A_shift + eps_sim)
  } else {
    sim <- sim %>% mutate(PC1_sim = mu + eps_sim)
  }
  
  sim
}

corr_from_scores <- function(df_sim){
  pc1_wide <- df_sim %>%
    select(Colony, Time, PC1_sim) %>%
    pivot_wider(names_from=Colony, values_from=PC1_sim) %>%
    arrange(Time)
  M <- as.matrix(pc1_wide %>% select(-Time))
  cor(M, use="pairwise.complete.obs", method="spearman")
}

apiary_delta_from_C <- function(C){
  cols <- colnames(C)
  cors_df <- expand.grid(Col1 = cols, Col2 = cols, stringsAsFactors=FALSE) %>%
    mutate(r = C[cbind(match(Col1, cols), match(Col2, cols))]) %>%
    filter(match(Col1, cols) < match(Col2, cols)) %>%
    left_join(col_meta, by=c("Col1"="Colony")) %>% rename(Apiary1=Apiary) %>%
    left_join(col_meta, by=c("Col2"="Colony")) %>% rename(Apiary2=Apiary) %>%
    mutate(within = (Apiary1 == Apiary2))
  
  mean(cors_df$r[cors_df$within], na.rm=TRUE) - mean(cors_df$r[!cors_df$within], na.rm=TRUE)
}

set.seed(99)
B <- 5000

sim_null <- replicate(B, {
  df <- simulate_PC1(include_apiary=FALSE)
  C <- corr_from_scores(df)
  c(overall = mean(C[upper.tri(C)], na.rm=TRUE),
    delta = apiary_delta_from_C(C))
})
sim_null <- t(sim_null)

sim_fit <- replicate(B, {
  df <- simulate_PC1(include_apiary=TRUE)
  C <- corr_from_scores(df)
  c(overall = mean(C[upper.tri(C)], na.rm=TRUE),
    delta = apiary_delta_from_C(C))
})
sim_fit <- t(sim_fit)




cors_df2 <- cors_df %>%
  mutate(group = ifelse(within_apiary, "Within apiary", "Between apiaries"))

overall_mean <- mean(cors_df2$r, na.rm=TRUE)

ggplot(cors_df2, aes(x = group, y = r)) +
  geom_violin(trim=TRUE, alpha=0.4) +
  geom_jitter(width=0.12, height=0, alpha=0.7, size=2) +
  stat_summary(fun=mean, geom="point", size=3) +
  geom_hline(yintercept = overall_mean, linetype=2) +
  labs(x = NULL,
       y = "Pairwise synchrony (Spearman r of PC1)",
       title = "Colony synchrony: within-apiary pairs are more similar") +
  theme_bw()



scores_plot <- scores %>%
  mutate(Colony = as.character(Colony)) %>%
  left_join(col_meta, by="Colony")

ggplot(scores_plot, aes(Time, PC1, group=Colony, color=Apiary)) +
  geom_line(alpha=0.35) +
  stat_summary(aes(group=Apiary), fun=mean, geom="line", linewidth=1.2) +
  theme_bw() +
  labs(title="Community trajectory (PC1) through time",
       y="PC1 (Bray–Curtis PCoA)")







cors_df_glmm <- cors_df %>%
  mutate(
    z = atanh(pmin(pmax(r, -0.999999), 0.999999)),
    Col1 = factor(Col1),
    Col2 = factor(Col2),
    within = within_apiary
  )

m0 <- lmer(z ~ 1 + (1|Col1) + (1|Col2), data=cors_df_glmm, REML=FALSE)
m1 <- lmer(z ~ within + (1|Col1) + (1|Col2), data=cors_df_glmm, REML=FALSE)

anova(m0, m1)
summary(m1)










obs_delta <- obs_delta
obs_overall <- overall_sync_obs

df_sim <- rbind(
  data.frame(model="Null (no apiary shared)", overall=sim_null[, "overall"], delta=sim_null[, "delta"]),
  data.frame(model="Fitted (apiary shared)", overall=sim_fit[, "overall"],  delta=sim_fit[, "delta"])
)

ggplot(df_sim, aes(x=delta, fill=model)) +
  geom_histogram(bins=60, alpha=0.5, position="identity") +
  geom_vline(xintercept=obs_delta, linetype=2, linewidth=1) +
  theme_bw() +
  labs(x="Δ = mean(within-apiary r) − mean(between-apiary r)",
       y="Count", title="Apiary structure: observed vs simulated")



ggplot(df_sim, aes(x=overall, fill=model)) +
  geom_histogram(bins=60, alpha=0.5, position="identity") +
  geom_vline(xintercept=obs_overall, linetype=2, linewidth=1) +
  theme_bw() +
  labs(x="Mean off-diagonal r", y="Count",
       title="Overall synchrony: observed vs simulated")

emmeans(m1, ~ within, type="response")  # response here still on z-scale


m1b <- lmer(z ~ within + (1|Col1), data=cors_df_glmm, REML=FALSE)
anova(m0, m1b)
summary(m1b)




# Desired fixed order
col_order <- c("100","101","102","103","104","105","106","107","555","777","888","999")

# Reorder the correlation matrix
cors_ord <- cors_obs[col_order, col_order]

# If you have an apiary annotation dataframe (rownames = colony IDs)
ann_col <- col_meta %>%
  mutate(Colony = as.character(Colony)) %>%
  filter(Colony %in% col_order) %>%
  tibble::column_to_rownames("Colony")

pheatmap(cors_ord,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = ann_col,
         annotation_col = ann_col,
         main = "Colony synchrony (Spearman corr of PC1)")






genus_ab <- mag_by_sample %>%
  group_by(sample_name, Colony, Time, Season, Rep, Genus) %>%
  summarise(gen_ab = sum(mag_mean, na.rm=TRUE), .groups="drop")

genus_ab <- genus_ab %>%
  group_by(Colony, Time, Season, Genus) %>%
  summarise(gen_ab = mean(gen_ab, na.rm=TRUE), .groups="drop")


tau <- 0  # or a small threshold
genus_season <- genus_ab %>%
  mutate(present = gen_ab > tau,
         log_ab = log1p(gen_ab)) %>%
  group_by(Genus, Season) %>%
  summarise(prev = mean(present),
            med_log_ab = median(log_ab[present], na.rm=TRUE),
            n_samples = n(),
            .groups="drop")





genus_bias <- genus_season %>%
  group_by(Genus) %>%
  mutate(
    top_season = Season[which.max(prev)],
    top_prev = max(prev, na.rm=TRUE),
    other_prev = mean(prev[Season != top_season], na.rm=TRUE),
    delta_prev = top_prev - other_prev,
    
    top_med = med_log_ab[Season == top_season][1],
    other_med = median(med_log_ab[Season != top_season], na.rm=TRUE),
    delta_med = top_med - other_med
  ) %>%
  distinct(Genus, top_season, top_prev, other_prev, delta_prev, delta_med) %>%
  arrange(desc(delta_prev), desc(delta_med))




strong_genus <- genus_bias %>%
  filter(delta_prev >= 0.2, delta_med >= 1)  # example thresholds



mag_ab <- mag_by_sample %>%
  group_by(Colony, Time, Season, Genus, MAG_id) %>%
  summarise(ab = mean(mag_mean, na.rm=TRUE), .groups="drop")


mag_ab <- mag_ab %>%
  group_by(Colony, Time, Season, Genus) %>%
  mutate(gen_total = sum(ab, na.rm=TRUE),
         frac = ifelse(gen_total > 0, ab / gen_total, NA_real_)) %>%
  ungroup()

# dominance per timepoint (top MAG fraction)
dom <- mag_ab %>%
  group_by(Colony, Time, Season, Genus) %>%
  summarise(dom_top = max(frac, na.rm=TRUE), .groups="drop")




dom_biased <- dom %>%
  left_join(genus_bias %>% select(Genus, top_season), by="Genus") %>%
  filter(Season == top_season)

dom_summary <- dom_biased %>%
  group_by(Genus, Colony) %>%
  summarise(med_dom_top = median(dom_top, na.rm=TRUE),
            .groups="drop")

all_community_decay_data$Microbial







# 1) Clean your decay data
decay <- all_community_decay_data$Microbial %>%
  mutate(dissim = 1 - Similarity) %>%
  filter(TimeLag > 0)   # drop replicate/zero-lag pairs for the decay slope

# 2) Recreate all within-colony timepoint pairs from PC1
pairs_pc1 <- scores %>%
  transmute(Colony = as.character(Colony), Time, PC1) %>%
  group_by(Colony) %>%
  tidyr::crossing(
    Time2 = Time,
    PC1_2 = PC1
  ) %>%
  ungroup() %>%
  rename(Time1 = Time, PC1_1 = PC1) %>%
  filter(Time2 >= Time1) %>%
  mutate(TimeLag = Time2 - Time1,
         dPC1 = abs(PC1_2 - PC1_1)) %>%
  filter(TimeLag > 0)

# 3) Summarize mean dPC1 by lag within colony (so it matches your decay data structure)
pc1_by_lag <- pairs_pc1 %>%
  group_by(Colony, TimeLag) %>%
  summarise(mean_dPC1 = mean(dPC1, na.rm=TRUE), .groups="drop")

pc1_by_lag
decay <- all_community_decay_data$Microbial %>%
  mutate(
    Colony = as.character(Colony),
    bray = Similarity  # IMPORTANT: this is Bray–Curtis dissimilarity
  ) %>%
  filter(TimeLag > 0)

bridge <- decay %>%
  left_join(pc1_by_lag, by = c("Colony", "TimeLag"))


# quick check: does dissimilarity increase with dPC1?
summary(lm(bray ~ mean_dPC1, data = bridge))



ggplot(bridge, aes(x = mean_dPC1, y = bray)) +
  geom_point(alpha = 0.15, position = position_jitter(width = 0.002, height = 0)) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  labs(
    x = "Mean |ΔPC1| for that lag (within colony)",
    y = "Bray–Curtis dissimilarity",
    title = "Temporal decay scales with movement along the seasonal axis (PC1)"
  )


ggplot(bridge, aes(mean_dPC1, bray)) +
  geom_bin2d() +
  geom_smooth(method="lm", se=FALSE, color="black") +
  theme_bw() +
  labs(x="Mean |ΔPC1| (within colony, by lag)", y="Bray–Curtis dissimilarity")



bridge$Colony <- factor(bridge$Colony)

m <- gam(bray ~ s(mean_dPC1, k=6) + s(TimeLag, k=8) + s(Colony, bs="re"),
         data = bridge, method = "REML")
summary(m)

decay2 <- decay %>%
  mutate(Colony = factor(Colony),
         logLag = log1p(TimeLag))

m_decay <- lmer(dissim ~ logLag + (1 + logLag | Colony), data=decay2, REML=FALSE)
summary(m_decay)



genus_long <- mag_by_sample %>%
  group_by(sample_name, Colony, Time, Season, Genus) %>%
  summarise(ab = sum(mag_mean, na.rm=TRUE), .groups="drop")


genus_wide <- genus_long %>%
  pivot_wider(names_from = Genus, values_from = ab, values_fill = 0)

meta <- genus_wide %>% select(sample_name, Colony, Time, Season)
X <- genus_wide %>% select(-sample_name, -Colony, -Time, -Season)

X_rel <- sweep(as.matrix(X), 1, rowSums(X) + 1e-12, "/")


bc <- vegdist(X_rel, method="bray")

bd <- betadisper(bc, group = meta$Season)  # distances to season centroid

# Classic ANOVA (not permutation-based)
anova(bd)

# Permutation test, stratified within Colony (respects repeated measures)
set.seed(1)
perm <- permutest(bd, permutations = 9999, strata = meta$Colony)
perm




bd_df <- data.frame(sample = meta$sample_name,
                    Colony = meta$Colony,
                    Season = meta$Season,
                    dist_to_centroid = bd$distances)


ggplot(bd_df, aes(x=Season, y=dist_to_centroid)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Colony), width=0.15, alpha=0.55, size=1.6) +
  theme_bw() +
  labs(y="Distance to season centroid (Bray–Curtis)",
       title="Seasonal convergence test (betadisper)")


TukeyHSD(bd)



ord <- cmdscale(bc, k=2, eig=TRUE)
ord_df <- data.frame(PCo1 = ord$points[,1],
                     PCo2 = ord$points[,2],
                     Season = meta$Season,
                     Colony = meta$Colony)

ggplot(ord_df, aes(PCo1, PCo2, color=Season)) +
  geom_point(alpha=0.8) +
  theme_bw() +
  labs(title="PCoA of genus composition (Bray–Curtis)")


range(all_community_decay_data$Microbial$Similarity)






# population-level smooth (exclude colony random effects)
grid <- expand.grid(
  mean_dPC1 = seq(min(bridge$mean_dPC1, na.rm=TRUE),
                  max(bridge$mean_dPC1, na.rm=TRUE), length.out=200),
  TimeLag   = median(bridge$TimeLag, na.rm=TRUE),
  Colony    = levels(bridge$Colony)[1]   # placeholder; we'll exclude it
)

pred <- predict(m, newdata=grid, se.fit=TRUE, exclude="s(Colony)")
grid$fit <- pred$fit
grid$lo  <- pred$fit - 1.96*pred$se.fit
grid$hi  <- pred$fit + 1.96*pred$se.fit

ggplot(grid, aes(mean_dPC1, fit)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.2) +
  geom_line(linewidth=1) +
  theme_bw() +
  labs(x="Mean |ΔPC1| (within-colony, by lag)",
       y="Predicted Bray–Curtis dissimilarity",
       title="Seasonal-axis movement predicts turnover (controlling for time lag and colony)")





grid2 <- expand.grid(
  mean_dPC1 = seq(min(bridge$mean_dPC1, na.rm=TRUE),
                  max(bridge$mean_dPC1, na.rm=TRUE), length.out=80),
  TimeLag = seq(min(bridge$TimeLag, na.rm=TRUE),
                max(bridge$TimeLag, na.rm=TRUE), length.out=80),
  Colony = levels(bridge$Colony)[1]
)

grid2$fit <- predict(m, newdata=grid2, exclude="s(Colony)")




m0 <- gam(bray ~ s(TimeLag, k=8) + s(Colony, bs="re"),
          data=bridge, method="ML")

m1 <- gam(bray ~ s(mean_dPC1, k=6) + s(TimeLag, k=8) + s(Colony, bs="re"),
          data=bridge, method="ML")

anova(m0, m1, test="Chisq")
summary(m0)$dev.expl
summary(m1)$dev.expl







# 1) Partial residuals for s(mean_dPC1):
# subtract prediction from the model WITHOUT that term (so what's left is f(mean_dPC1)+noise)
bridge$partial_y <- bridge$bray - predict(m1, newdata=bridge, exclude="s(mean_dPC1)")

# 2) Smooth for s(mean_dPC1) alone (with SE)
grid <- data.frame(
  mean_dPC1 = seq(min(bridge$mean_dPC1, na.rm=TRUE),
                  max(bridge$mean_dPC1, na.rm=TRUE), length.out=200),
  TimeLag = median(bridge$TimeLag, na.rm=TRUE),
  Colony = levels(bridge$Colony)[1]
)

pr <- predict(m1, newdata=grid, type="terms", terms="s(mean_dPC1)", se.fit=TRUE)
grid$fit <- pr$fit[,1]
grid$se  <- pr$se.fit[,1]
grid$lo  <- grid$fit - 1.96*grid$se
grid$hi  <- grid$fit + 1.96*grid$se

# 3) Plot: partial residual points + partial smooth

ggplot(bridge, aes(mean_dPC1, partial_y)) +
  geom_point(alpha=0.10, size=0.7) +
  geom_ribbon(
    data = grid,
    aes(x = mean_dPC1, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  geom_line(
    data = grid,
    aes(x = mean_dPC1, y = fit),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    x="Mean |ΔPC1| (within-colony, by lag)",
    y="Partial Bray–Curtis dissimilarity\n(adjusted for time lag + colony)",
    title="Seasonal-axis displacement predicts turnover beyond time lag"
  )














eps <- 1e-8

# ------------------------------------------------------------
# 1. Zero-complete MAG table within each colony
# ------------------------------------------------------------

sample_meta <- mag_by_sample %>%
  distinct(sample_name, Colony, Month, Time, Rep, Season)

mag_colony_meta <- mag_by_sample %>%
  distinct(Colony, MAG_id, Genus)

mag_by_sample0 <- sample_meta %>%
  left_join(mag_colony_meta, by = "Colony", relationship = "many-to-many") %>%
  left_join(
    mag_by_sample %>% select(sample_name, MAG_id, mag_mean),
    by = c("sample_name", "MAG_id")
  ) %>%
  mutate(mag_mean = coalesce(mag_mean, 0))

# ------------------------------------------------------------
# 2. Relative abundance per sample
# ------------------------------------------------------------

mag_rel_by_sample <- mag_by_sample0 %>%
  group_by(sample_name) %>%
  mutate(
    total_microbial = sum(mag_mean, na.rm = TRUE),
    rel_abund = if_else(total_microbial > 0, mag_mean / total_microbial, 0)
  ) %>%
  ungroup()





# ------------------------------------------------------------
# Replicate-level variation
# ------------------------------------------------------------

replicate_noise <- mag_rel_by_sample %>%
  group_by(Colony, Month, Time, MAG_id, Genus) %>%
  summarise(
    n_rep = n_distinct(sample_name),
    mean_rel = mean(rel_abund, na.rm = TRUE),
    var_rep  = var(rel_abund, na.rm = TRUE),
    sd_rep   = sd(rel_abund, na.rm = TRUE),
    cv_rep   = sd_rep / (mean_rel + eps),
    .groups = "drop"
  ) %>%
  filter(n_rep >= 2)

# ------------------------------------------------------------
# Replicate-averaged abundance per colony-time-MAG
# ------------------------------------------------------------

mag_monthly <- mag_rel_by_sample %>%
  group_by(Colony, Month, Time, Season, MAG_id, Genus) %>%
  summarise(
    n_rep = n_distinct(sample_name),
    x_bar = mean(rel_abund, na.rm = TRUE),
    s2_rep = var(rel_abund, na.rm = TRUE),
    .groups = "drop"
  )








# ------------------------------------------------------------
# Variance partition per Colony x MAG
# ------------------------------------------------------------

mag_noise_partition <- mag_monthly %>%
  group_by(Colony, MAG_id, Genus) %>%
  filter(n_distinct(Time) >= 4) %>%
  summarise(
    T = n_distinct(Time),
    mu = mean(x_bar, na.rm = TRUE),
    var_obs = var(x_bar, na.rm = TRUE),
    mean_noise_in_mean = mean(s2_rep / pmax(n_rep, 1), na.rm = TRUE),
    var_bio = pmax(var_obs - mean_noise_in_mean, 0),
    noise_fraction = if_else(var_obs > 0, mean_noise_in_mean / var_obs, NA_real_),
    .groups = "drop"
  ) %>%
  filter(is.finite(mu), mu > 0, is.finite(var_obs))




noise_summary <- mag_noise_partition %>%
  summarise(
    n = n(),
    median_noise_fraction = median(noise_fraction, na.rm = TRUE),
    IQR_noise_fraction = IQR(noise_fraction, na.rm = TRUE),
    frac_noise_gt_0.5 = mean(noise_fraction > 0.5, na.rm = TRUE),
    frac_noise_lt_0.2 = mean(noise_fraction < 0.2, na.rm = TRUE)
  )

noise_summary




# ------------------------------------------------------------
# Taylor's Law: raw variance
# ------------------------------------------------------------

taylor_raw <- mag_noise_partition %>%
  filter(mu > 0, var_obs > 0)

taylor_fit_raw <- lm(log(var_obs) ~ log(mu), data = taylor_raw)

# ------------------------------------------------------------
# Taylor's Law: de-noised biological variance
# ------------------------------------------------------------

taylor_bio <- mag_noise_partition %>%
  filter(mu > 0, var_bio > 0)

taylor_fit_bio <- lm(log(var_bio) ~ log(mu), data = taylor_bio)

taylor_summary <- tibble::tibble(
  model = c("raw", "denoised"),
  beta = c(
    coef(taylor_fit_raw)[["log(mu)"]],
    coef(taylor_fit_bio)[["log(mu)"]]
  ),
  intercept = c(
    coef(taylor_fit_raw)[["(Intercept)"]],
    coef(taylor_fit_bio)[["(Intercept)"]]
  ),
  r2 = c(
    summary(taylor_fit_raw)$r.squared,
    summary(taylor_fit_bio)$r.squared
  ),
  n = c(nrow(taylor_raw), nrow(taylor_bio))
)

taylor_summary





# ------------------------------------------------------------
# SLM-inspired fluctuation parameter and effective K
# ------------------------------------------------------------

slm_df <- mag_noise_partition %>%
  mutate(
    cv2_bio = var_bio / (mu^2 + eps),
    sigma_ic = (2 * cv2_bio) / (1 + cv2_bio)
  ) %>%
  filter(is.finite(sigma_ic), sigma_ic >= 0)

sigma_global <- median(slm_df$sigma_ic, na.rm = TRUE)

slm_df <- slm_df %>%
  mutate(
    sigma_global = sigma_global,
    K_ic = mu / (1 - sigma_global / 2),
    stability_rank = percent_rank(sigma_ic),
    stability_class = case_when(
      sigma_ic <= quantile(sigma_ic, 0.25, na.rm = TRUE) ~ "low fluctuation",
      sigma_ic >= quantile(sigma_ic, 0.75, na.rm = TRUE) ~ "high fluctuation",
      TRUE ~ "intermediate"
    )
  )

slm_summary <- slm_df %>%
  summarise(
    n = n(),
    sigma_median = median(sigma_ic, na.rm = TRUE),
    sigma_IQR = IQR(sigma_ic, na.rm = TRUE),
    K_median = median(K_ic, na.rm = TRUE),
    K_IQR = IQR(K_ic, na.rm = TRUE)
  )

slm_summary






genus_slm_summary <- slm_df %>%
  group_by(Genus) %>%
  summarise(
    n_colony_mag = n(),
    n_colonies = n_distinct(Colony),
    median_mu = median(mu, na.rm = TRUE),
    median_sigma = median(sigma_ic, na.rm = TRUE),
    IQR_sigma = IQR(sigma_ic, na.rm = TRUE),
    median_K = median(K_ic, na.rm = TRUE),
    frac_low_fluctuation = mean(stability_class == "low fluctuation", na.rm = TRUE),
    frac_high_fluctuation = mean(stability_class == "high fluctuation", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_colony_mag))

genus_slm_summary





# ------------------------------------------------------------
# Month-to-month log abundance changes
# ------------------------------------------------------------

mag_changes <- mag_monthly %>%
  arrange(Colony, MAG_id, Time) %>%
  group_by(Colony, MAG_id, Genus) %>%
  mutate(
    dlog = log((x_bar + eps) / (lag(x_bar) + eps)),
    dt = Time - lag(Time)
  ) %>%
  ungroup() %>%
  filter(!is.na(dlog), is.finite(dlog), dt > 0)

b_hat <- mean(abs(mag_changes$dlog), na.rm = TRUE)
sd_hat <- sd(mag_changes$dlog, na.rm = TRUE)

ll_laplace <- sum(-log(2 * b_hat) - abs(mag_changes$dlog) / b_hat, na.rm = TRUE)
ll_normal <- sum(dnorm(mag_changes$dlog, mean = 0, sd = sd_hat, log = TRUE), na.rm = TRUE)

change_fit_summary <- tibble::tibble(
  model = c("Laplace", "Normal"),
  logLik = c(ll_laplace, ll_normal),
  parameter = c(b_hat, sd_hat)
) %>%
  mutate(
    delta_logLik_vs_best = max(logLik) - logLik
  )

change_fit_summary










# ------------------------------------------------------------
# Plot labels
# ------------------------------------------------------------

noise_med <- median(mag_noise_partition$noise_fraction, na.rm = TRUE)
noise_iqr <- IQR(mag_noise_partition$noise_fraction, na.rm = TRUE)

beta_bio <- coef(taylor_fit_bio)[["log(mu)"]]
r2_bio <- summary(taylor_fit_bio)$r.squared

sigma_med <- median(slm_df$sigma_ic, na.rm = TRUE)
sigma_iqr <- IQR(slm_df$sigma_ic, na.rm = TRUE)
sigma_q25 <- quantile(slm_df$sigma_ic, 0.25, na.rm = TRUE)
sigma_q75 <- quantile(slm_df$sigma_ic, 0.75, na.rm = TRUE)

# ------------------------------------------------------------
# A. Replicate CV vs mean relative abundance
# ------------------------------------------------------------

pA <- ggplot(replicate_noise, aes(x = mean_rel, y = cv_rep)) +
  geom_point(alpha = 0.35, size = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Mean relative abundance",
    y = "Replicate CV",
    title = "A. Replicate variability"
  ) +
  theme_bw()
pA
# ------------------------------------------------------------
# B. Noise fraction
# ------------------------------------------------------------

pB <- ggplot(mag_noise_partition, aes(x = noise_fraction)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = noise_med, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  coord_cartesian(xlim = c(0, 1.25)) +
  labs(
    x = "Replicate-noise fraction",
    y = "Colony × rMAG trajectories",
    title = "B. Replicate noise explains a minority of temporal variance",
    subtitle = paste0("Median = ", round(noise_med, 3), 
                      "; IQR = ", round(noise_iqr, 3))
  ) +
  theme_bw()
pB
# ------------------------------------------------------------
# C. Taylor's Law: de-noised biological variance
# ------------------------------------------------------------

pC <- ggplot(taylor_bio, aes(x = log(mu), y = log(var_bio))) +
  geom_point(alpha = 0.45, size = 0.9) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "log mean relative abundance",
    y = "log de-noised temporal variance",
    title = "C. Taylor's Law in rMAG abundance fluctuations",
    subtitle = paste0(
      "Slope = ", round(beta_bio, 2),
      "; R² = ", round(r2_bio, 2)
    )
  ) +
  theme_bw()
pC
# ------------------------------------------------------------
# D. SLM-derived sigma
# ------------------------------------------------------------

pD <- ggplot(slm_df, aes(x = sigma_ic)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = sigma_med, linetype = "dashed") +
  geom_vline(xintercept = sigma_q25, linetype = "dotted") +
  geom_vline(xintercept = sigma_q75, linetype = "dotted") +
  labs(
    x = expression("SLM-derived fluctuation parameter " * sigma[ic]),
    y = "Colony × rMAG trajectories",
    title = "D. Heterogeneity in rMAG fluctuation magnitude",
    subtitle = paste0(
      "Median = ", round(sigma_med, 3),
      "; IQR = ", round(sigma_iqr, 3)
    )
  ) +
  theme_bw()
pD
# ------------------------------------------------------------
# E. Effective K distribution
# ------------------------------------------------------------

pE <- ggplot(slm_df, aes(x = K_ic)) +
  geom_histogram(bins = 40) +
  scale_x_log10() +
  labs(
    x = expression("Effective relative carrying capacity " * K[ic]),
    y = "Colony × rMAG trajectories",
    title = "E. SLM-derived effective relative carrying capacity"
  ) +
  theme_bw()
pE
# ------------------------------------------------------------
# F. Genus-level fluctuation summary
# ------------------------------------------------------------

pF <- genus_slm_plot_df %>%
  mutate(Genus = fct_reorder(Genus, median_sigma)) %>%
  ggplot(aes(x = Genus, y = median_sigma)) +
  geom_point(aes(size = n_colony_mag), alpha = 0.8) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Median " * sigma[ic]),
    size = "Trajectories",
    title = "F. rMAG fluctuation magnitude differs among genera"
  ) +
  theme_bw()

supp_slm_fig <- (pA | pB) / (pC | pD) / (pE | pF)

supp_slm_fig

ggsave(
  "supplemental_rMAG_noise_Taylor_SLM.pdf",
  supp_slm_fig,
  width = 12,
  height = 13,
  units = "in"
)

















# -----------------------------
# Simple summary statistics
# -----------------------------

simple_summary <- mag_noise_partition %>%
  mutate(
    noise_dominated = mean_noise_in_mean >= var_obs | var_bio == 0,
    low_noise = noise_fraction < 0.2,
    high_noise = noise_fraction > 0.5
  ) %>%
  summarise(
    n = n(),
    median_noise_fraction = median(noise_fraction, na.rm = TRUE),
    IQR_noise_fraction = IQR(noise_fraction, na.rm = TRUE),
    frac_low_noise = mean(low_noise, na.rm = TRUE),
    frac_high_noise = mean(high_noise, na.rm = TRUE),
    frac_noise_dominated = mean(noise_dominated, na.rm = TRUE)
  )

simple_summary

# -----------------------------
# Taylor's Law, denoised only
# -----------------------------

taylor_bio <- mag_noise_partition %>%
  filter(mu > 0, var_bio > 0)

taylor_fit_bio <- lm(log(var_bio) ~ log(mu), data = taylor_bio)

taylor_beta <- coef(taylor_fit_bio)[["log(mu)"]]
taylor_r2 <- summary(taylor_fit_bio)$r.squared

# -----------------------------
# SLM-inspired sigma as fluctuation index
# -----------------------------

slm_df <- mag_noise_partition %>%
  filter(mu > 0, is.finite(mu), is.finite(var_bio)) %>%
  mutate(
    cv2_bio = if_else(var_bio > 0, var_bio / mu^2, 0),
    sigma_ic = (2 * cv2_bio) / (1 + cv2_bio),
    low_fluctuation = sigma_ic < 1
  )

sigma_summary <- slm_df %>%
  summarise(
    n = n(),
    median_sigma = median(sigma_ic, na.rm = TRUE),
    IQR_sigma = IQR(sigma_ic, na.rm = TRUE),
    frac_sigma_lt_1 = mean(sigma_ic < 1, na.rm = TRUE),
    frac_sigma_eq_0 = mean(sigma_ic == 0, na.rm = TRUE)
  )

sigma_summary


  noise_med <- median(mag_noise_partition$noise_fraction, na.rm = TRUE)
noise_iqr <- IQR(mag_noise_partition$noise_fraction, na.rm = TRUE)

sigma_med <- median(slm_df$sigma_ic, na.rm = TRUE)
sigma_lt1 <- mean(slm_df$sigma_ic < 1, na.rm = TRUE)

# A. Observed variance vs replicate-noise variance
pA <- ggplot(
  mag_noise_partition,
  aes(x = mean_noise_in_mean, y = var_obs)
) +
  geom_point(alpha = 0.35, size = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Expected replicate-noise variance",
    y = "Observed temporal variance",
    title = "A. Most rMAG trajectories exceed replicate noise"
  ) +
  theme_bw()

# B. Noise fraction
pB <- ggplot(mag_noise_partition, aes(x = noise_fraction)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = noise_med, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  coord_cartesian(xlim = c(0, 1.25)) +
  labs(
    x = "Replicate-noise fraction",
    y = "Colony × rMAG trajectories",
    title = "B. Replicate noise explains a minority of temporal variance",
    subtitle = paste0(
      "Median = ", round(noise_med, 3),
      "; IQR = ", round(noise_iqr, 3)
    )
  ) +
  theme_bw()

# C. Taylor's Law
pC <- ggplot(taylor_bio, aes(x = log(mu), y = log(var_bio))) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "log mean relative abundance",
    y = "log de-noised temporal variance",
    title = "C. rMAG abundance fluctuations follow Taylor's Law",
    subtitle = paste0(
      "Slope = ", round(taylor_beta, 2),
      "; R² = ", round(taylor_r2, 2)
    )
  ) +
  theme_bw()

# D. Sigma fluctuation index
pD <- ggplot(slm_df, aes(x = sigma_ic)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = sigma_med, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  labs(
    x = expression("SLM-derived fluctuation index " * sigma[ic]),
    y = "Colony × rMAG trajectories",
    title = "D. rMAGs vary in de-noised fluctuation magnitude",
    subtitle = paste0(
      "Median = ", round(sigma_med, 3),
      "; fraction ", expression(sigma[ic] < 1), " = ", round(sigma_lt1, 3)
    )
  ) +
  theme_bw()
pA
pB
pC
pD


ggsave(
  "supplemental_rMAG_stability_noise_Taylor.pdf",
  supp_simple,
  width = 11,
  height = 8.5,
  units = "in"
)




pA <- ggplot(mag_noise_partition, aes(x = mean_noise_in_mean, y = var_obs)) +
  geom_point(alpha = 0.35, size = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Expected replicate-level variance",
    y = "Observed temporal variance",
    title = "A. Most rMAG trajectories exceed replicate-level variance"
  ) +
  theme_bw()

pB <- ggplot(mag_noise_partition, aes(x = noise_fraction)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = median(mag_noise_partition$noise_fraction, na.rm = TRUE),
             linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  coord_cartesian(xlim = c(0, 1.25)) +
  labs(
    x = "Fraction of temporal variance attributable to replicate-level variation",
    y = "Colony × rMAG trajectories",
    title = "B. Replicate-level variation explains a minority of temporal variance"
  ) +
  theme_bw()

print(pA)
pB
plot(1:10)
R.version.string
RStudio.Version()$version
