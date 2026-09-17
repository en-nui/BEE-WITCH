library(data.table)
library(ggplot2)

# If not already loaded:
# gene_dt  <- fread("master_per_gene_summary.csv")
# sweep_dt <- fread("sweep_gene_merged.csv")

setDT(gene_dt)
setDT(sweep_dt)

# Keep only columns we need
gene_small <- gene_dt[
  , .(
    colony_id,
    MAG_id,
    corresponding_gene_call,
    Month,
    piAll = pi_all_sites,
    piS   = pi_syn,
    piN   = pi_nonsyn
  )
]

# Drop rows with all pi == NA
gene_small <- gene_small[
  !(is.na(piAll) & is.na(piS) & is.na(piN))
]

month_levels <- c(
  "May","June","July","August",
  "September","October","November",
  "January","February"
)

gene_small[, Month := factor(Month, levels = month_levels)]
gene_small <- gene_small[!is.na(Month)]

gene_small[, month_idx := as.integer(Month)]



# Make two copies: one for t1, one for t2
g1 <- copy(gene_small)
g2 <- copy(gene_small)

setnames(
  g1,
  old = c("Month",     "month_idx",      "piAll",      "piS",      "piN"),
  new = c("Month_t1",  "month_idx_t1",   "piAll_t1",   "piS_t1",   "piN_t1")
)

setnames(
  g2,
  old = c("Month",     "month_idx",      "piAll",      "piS",      "piN"),
  new = c("Month_t2",  "month_idx_t2",   "piAll_t2",   "piS_t2",   "piN_t2")
)

setkey(g1, colony_id, MAG_id, corresponding_gene_call, month_idx_t1)
setkey(g2, colony_id, MAG_id, corresponding_gene_call, month_idx_t2)

gene_pairs <- g1[g2, nomatch = 0][
  month_idx_t2 == month_idx_t1 + 1L
]

gene_pairs[, `:=`(
  d_piAll = piAll_t2 - piAll_t1,
  d_piS   = piS_t2   - piS_t1,
  d_piN   = piN_t2   - piN_t1
)]



setDT(sweep_dt)

gene_sweep_flags <- sweep_dt[
  ,
  .(
    has_0D = any(degeneracy == "0D"),
    has_4D = any(degeneracy == "4D")
  ),
  by = .(colony_id, MAG_id, corresponding_gene_call, Month_t1, Month_t2)
]

setkey(gene_pairs,
       colony_id, MAG_id, corresponding_gene_call, Month_t1, Month_t2)
setkey(gene_sweep_flags,
       colony_id, MAG_id, corresponding_gene_call, Month_t1, Month_t2)

gene_pairs <- gene_sweep_flags[gene_pairs]

gene_pairs[is.na(has_0D), has_0D := FALSE]
gene_pairs[is.na(has_4D), has_4D := FALSE]


gene_pairs[
  ,
  status := fifelse(has_0D,
                    "0D_sweep",
                    fifelse(!has_0D & !has_4D, "null", "other"))
]

gene_pairs2 <- gene_pairs[status %in% c("0D_sweep", "null")]


df_unmatched_piS_t1 <- gene_pairs2[
  ,
  .(status, piS_t1)
]

df_unmatched_piS_t2 <- gene_pairs2[
  ,
  .(status, piS_t2)
]

w_t1 <- wilcox.test(piS_t1 ~ status, data = df_unmatched_piS_t1)
w_t2 <- wilcox.test(piS_t2 ~ status, data = df_unmatched_piS_t2)

w_t1
w_t2

df_unmatched_piS_t1[
  ,
  .(
    median_piS_t1 = median(piS_t1, na.rm = TRUE),
    mean_piS_t1   = mean(piS_t1,   na.rm = TRUE),
    n             = .N
  ),
  by = status
]

df_unmatched_piS_t2[
  ,
  .(
    median_piS_t2 = median(piS_t2, na.rm = TRUE),
    mean_piS_t2   = mean(piS_t2,   na.rm = TRUE),
    n             = .N
  ),
  by = status
]

ggplot(df_unmatched_piS_t1, aes(x = status, y = piS_t1, fill = status)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = c("0D_sweep" = "#d95f02", "null" = "#1b9e77")) +
  labs(
    x = "",
    y = expression(pi[S](t[1])),
    title = expression(pi[S] ~ "before interval: 0D-sweep vs null")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggplot(df_unmatched_piS_t2, aes(x = status, y = piS_t2, fill = status)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = c("0D_sweep" = "#d95f02", "null" = "#1b9e77")) +
  labs(
    x = "",
    y = expression(pi[S](t[2])),
    title = expression(pi[S] ~ "after interval: 0D-sweep vs null")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")


df_delta <- gene_pairs2[
  ,
  .(status, d_piS, piS_t1, piS_t2)
]

df_delta[
  ,
  .(mean_d_piS = mean(d_piS, na.rm = TRUE),
    median_d_piS = median(d_piS, na.rm = TRUE),
    n = .N),
  by = status
]

ggplot(df_delta, aes(x = status, y = d_piS, fill = status)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_manual(values = c("0D_sweep" = "#d95f02", "null" = "#1b9e77")) +
  labs(
    x = "",
    y = expression(Delta * pi[S]),
    title = expression("Change in " * pi[S] ~ " between consecutive months")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")+
  scale_y_continuous(trans="log10")+
  coord_cartesian(ylim=c(0.1,0.0001))

## ----------------------------------------------------------
## 6. Linear model: ΔπS ~ status + πS(t1)
## ----------------------------------------------------------

fit_delta <- lm(
  d_piS ~ status + piS_t1,
  data = df_delta[status %in% c("0D_sweep", "null")]
)

summary(fit_delta)

## ----------------------------------------------------------
## 7. Residual ΔπS after conditioning on πS(t1)
## ----------------------------------------------------------

df_delta_res <- df_delta[
  status %in% c("0D_sweep","null") &
    !is.na(d_piS) & !is.na(piS_t1)
]

fit_res <- lm(d_piS ~ piS_t1, data = df_delta_res)
df_delta_res[, resid := resid(fit_res)]

ggplot(df_delta_res, aes(x = status, y = resid, fill = status)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_manual(values = c("0D_sweep" = "#d95f02", "null" = "#1b9e77")) +
  labs(
    x = "",
    y = expression("Residual " * Delta * pi[S] ~ " (after conditioning on " * pi[S](t[1]) * ")"),
    title = expression("Excess loss of " * pi[S] ~ " in 0D-sweep intervals")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")+
  scale_y_continuous(trans="log10")+
  coord_cartesian(ylim=c(0.1,0.0001))



perm_means_two_groups <- function(x_sweep,
                                  x_null,
                                  n_sample = 500,
                                  n_perm   = 5000,
                                  seed     = 1L) {
  set.seed(seed)
  
  x_sweep <- x_sweep[!is.na(x_sweep)]
  x_null  <- x_null[!is.na(x_null)]
  
  if (length(x_sweep) == 0L || length(x_null) == 0L) {
    stop("One of the groups has zero non-NA values.")
  }
  
  replace_sweep <- length(x_sweep) < n_sample
  replace_null  <- length(x_null)  < n_sample
  
  sweep_means <- numeric(n_perm)
  null_means  <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    sweep_means[i] <- mean(sample(x_sweep, n_sample, replace = replace_sweep))
    null_means[i]  <- mean(sample(x_null,  n_sample, replace = replace_null))
  }
  
  diff_means <- sweep_means - null_means
  
  list(
    sweep_means     = sweep_means,
    null_means      = null_means,
    diff_means      = diff_means,
    p_sweep_gt_null = mean(diff_means > 0),
    p_sweep_lt_null = mean(diff_means < 0)
  )
}

res_piS_t2 <- perm_means_two_groups(
  x_sweep = df_unmatched_piS_t2[status == "0D_sweep", piS_t2],
  x_null  = df_unmatched_piS_t2[status == "null",     piS_t2],
  n_sample = 500,
  n_perm   = 5000,
  seed     = 1
)

res_piS_t2$p_sweep_gt_null  # P(mean πS_t2(sweep) > mean πS_t2(null))
res_piS_t2$p_sweep_lt_null  # P(mean πS_t2(sweep) < mean πS_t2(null))

res_d_piS <- perm_means_two_groups(
  x_sweep = df_delta[status == "0D_sweep", d_piS],
  x_null  = df_delta[status == "null",     d_piS],
  n_sample = 500,
  n_perm   = 5000,
  seed     = 1
)

res_d_piS$p_sweep_gt_null   # P(mean ΔπS(sweep) > mean ΔπS(null))
res_d_piS$p_sweep_lt_null   # P(mean ΔπS(sweep) < mean ΔπS(null))

plot_perm_means <- function(res, title = "Permutation distribution of mean πS") {
  df_plot <- rbind(
    data.frame(group = "0D_sweep", mean_val = res$sweep_means),
    data.frame(group = "null",     mean_val = res$null_means)
  )
  
  ggplot(df_plot, aes(x = mean_val, fill = group)) +
    geom_density(alpha = 0.4) +
    labs(
      x = "Mean value from subsamples",
      y = "Density",
      title = title
    ) +
    scale_fill_manual(values = c("0D_sweep" = "#d95f02", "null" = "#1b9e77")) +
    theme_minimal(base_size = 12)
}

plot_perm_means(
  res_piS_t2,
  title = expression("Mean " * pi[S](t[2]) * " in 0D-sweep genes vs null")
)

plot_perm_means(
  res_d_piS,
  title = expression("Mean " * Delta * pi[S] * " (0D-sweep vs null)")
)


