library(zoo)
# Make sure both are data.table
sweep_dt <- as.data.table(sweep_dt)
gene_clean <- as.data.table(gene_clean)

# Harmonize types
sweep_dt[, `:=`(
  MAG_id = as.character(MAG_id),
  colony_id = as.character(colony_id),
  corresponding_gene_call = as.character(corresponding_gene_call)
)]

gene_clean[, `:=`(
  MAG_id = as.character(MAG_id),
  Colony = as.character(Colony),
  corresponding_gene_call = as.character(corresponding_gene_call)
)]

# Build a unique mapping table from gene_clean
pfam_map <- unique(
  gene_clean[!is.na(PFAM_primary) & PFAM_primary != "",
             .(MAG_id,
               colony_id = Colony,
               corresponding_gene_call,
               PFAM_primary)]
)

# Join onto sweep_dt
sweep_dt <- merge(
  sweep_dt,
  pfam_map,
  by = c("MAG_id", "colony_id", "corresponding_gene_call"),
  all.x = TRUE
)

# quick QC
cat("Rows in sweep_dt:", nrow(sweep_dt), "\n")
cat("Missing PFAM_primary after join:", sweep_dt[is.na(PFAM_primary) | PFAM_primary=="", .N], "\n")
sweep_dt



# ---- Choose signal: 0D sweeps as "changes" ----
changes <- as.data.table(sweep_dt)

# Require these columns; rename if yours differ
# MAG_id, colony_id, Month_t1, Month_t2, corresponding_gene_call, PFAM_primary, degeneracy
stopifnot(all(c("MAG_id","colony_id","corresponding_gene_call","PFAM_primary","degeneracy") %in% names(changes)))

# If you have Month_t1/Month_t2 use them; otherwise set a unit another way
stopifnot(all(c("Month_t1","Month_t2") %in% names(changes)))

# Filter to 0D changes (your definition)
changes <- changes[degeneracy == "0D"]
changes <- changes[!is.na(PFAM_primary) & PFAM_primary != ""]

# Define unit_id (QP analog)
changes[, unit_id := paste(MAG_id, colony_id, Month_t1, Month_t2, sep="|")]

# Observed number of changes per unit
n_by_unit <- changes[, .N, by = unit_id]
setnames(n_by_unit, "N", "n_changes")

# Observed counts per PFAM
obs_pfam <- changes[, .N, by = PFAM_primary]
setnames(obs_pfam, "N", "obs_changes")







g <- as.data.table(gene_clean)
g
# Ensure key columns exist
stopifnot(all(c("MAG_id","Colony","Month","corresponding_gene_call","PFAM_primary") %in% names(g)))

# Start from gene_clean but keep only scalar columns needed for lookup + sampling
WCOL <- if ("effective_length_0D_sites" %in% names(gene_clean)) "effective_length_0D_sites" else "gene_length"

g_key <- as.data.table(gene_clean)[
  , .(
    MAG_id = as.character(MAG_id),
    Colony = as.character(Colony),
    Month  = as.character(Month),
    corresponding_gene_call = as.character(corresponding_gene_call),
    PFAM_primary = as.character(PFAM_primary),
    w = as.numeric(get(WCOL))
  )
]

# Clean
g_key <- g_key[
  !is.na(PFAM_primary) & PFAM_primary != "" &
    is.finite(w) & w > 0
]

# Remove duplicate gene rows within a MAG×Colony×Month
g_key <- unique(g_key, by = c("MAG_id","Colony","Month","corresponding_gene_call","PFAM_primary"))

# Now key works (no list columns)
setkey(g_key, MAG_id, Colony, Month)
g_key

changes <- as.data.table(sweep_dt)

changes <- changes[
  degeneracy == "0D" &
    !is.na(PFAM_primary) & PFAM_primary != ""
]

changes[, `:=`(
  MAG_id    = as.character(MAG_id),
  colony_id = as.character(colony_id),
  Month_t1  = as.character(Month_t1),
  Month_t2  = as.character(Month_t2)
)]

changes[, unit_id := paste(MAG_id, colony_id, Month_t1, Month_t2, sep = "|")]

unit_meta <- unique(
  changes[, .(
    unit_id,
    MAG_id,
    Colony = colony_id,
    m1 = Month_t1,
    m2 = Month_t2
  )]
)
changes
unit_meta
# Build per-unit universes using g_key (NOT g)
unit_univ_list <- vector("list", nrow(unit_meta))
names(unit_univ_list) <- unit_meta$unit_id

for (i in seq_len(nrow(unit_meta))) {
  ui  <- unit_meta$unit_id[i]
  mag <- unit_meta$MAG_id[i]
  col <- unit_meta$Colony[i]
  m1  <- unit_meta$m1[i]
  m2  <- unit_meta$m2[i]
  
  dt1 <- g_key[J(mag, col, m1), nomatch = 0L]
  dt2 <- g_key[J(mag, col, m2), nomatch = 0L]
  u   <- unique(rbind(dt1, dt2))
  
  unit_univ_list[[ui]] <- u[, .(corresponding_gene_call, PFAM_primary, w)]
}






















stopifnot(all(c("PFAM_primary","MAG_id","colony_id","unit_id") %in% names(changes)))

# observed PFAM stats (analog of "gene class" in paper)
obs_pfam <- changes[!is.na(PFAM_primary) & PFAM_primary != "", .(
  obs_changes = .N,                      # total SNV changes
  n_units  = uniqueN(unit_id),
  n_cols   = uniqueN(colony_id),         # "hosts"
  n_MAGs   = uniqueN(MAG_id),
  n_genus  = uniqueN(Genus)              # if Genus present; else join it first
), by = PFAM_primary]












bootstrap_null1_pfam <- function(B = 2000, seed = 1L,
                                 unit_meta, n_by_unit, unit_univ_list,
                                 pfams, pfam_index, obs_vec) {
  set.seed(seed)
  
  # unit -> n_changes
  nmap <- n_by_unit$n_changes
  names(nmap) <- n_by_unit$unit_id
  
  null_mat <- matrix(0L, nrow = B, ncol = length(pfams))
  colnames(null_mat) <- pfams
  
  for (b in seq_len(B)) {
    counts <- integer(length(pfams))
    
    for (i in seq_len(nrow(unit_meta))) {
      ui <- unit_meta$unit_id[i]
      n_u <- nmap[[ui]]
      if (is.null(n_u) || n_u <= 0) next
      
      U <- unit_univ_list[[ui]]
      if (is.null(U) || nrow(U) == 0) next
      
      idx <- sample.int(nrow(U), size = n_u, replace = TRUE, prob = U$w)
      pf <- U$PFAM_primary[idx]
      tab <- table(pf)
      
      jj <- pfam_index[names(tab)]
      ok <- !is.na(jj)
      if (any(ok)) counts[jj[ok]] <- counts[jj[ok]] + as.integer(tab[ok])
    }
    
    null_mat[b, ] <- counts
  }
  
  # empirical p-values: P(null >= obs)
  ge <- colSums(null_mat >= matrix(obs_vec, nrow = B, ncol = length(obs_vec), byrow = TRUE))
  p_emp <- (ge + 1) / (B + 1)
  
  data.table(
    PFAM_primary = pfams,
    obs_changes  = obs_vec,
    null_mean    = colMeans(null_mat),
    null_sd      = apply(null_mat, 2, sd),
    p_null       = p_emp,
    q_BH         = p.adjust(p_emp, method = "BH")
  )
}
# PFAM set must cover observable PFAMs under the null
pfams <- sort(unique(c(
  obs_pfam$PFAM_primary,
  unique(unlist(lapply(unit_univ_list, function(u) u$PFAM_primary)))
)))
pfam_index <- setNames(seq_along(pfams), pfams)

# observed vector aligned
obs_vec <- integer(length(pfams))
obs_vec[pfam_index[obs_pfam$PFAM_primary]] <- obs_pfam$obs_changes

# run null
B <- 500
res_null <- bootstrap_null1_pfam(
  B = B, seed = 1L,
  unit_meta = unit_meta,
  n_by_unit = n_by_unit,
  unit_univ_list = unit_univ_list,
  pfams = pfams,
  pfam_index = pfam_index,
  obs_vec = obs_vec
)




res
res_null




res <- merge(res_null, obs_pfam, by = "PFAM_primary", all.x = TRUE)

# enrichment / z for reporting
res[, FE := obs_changes.x / pmax(null_mean, 1e-12)]
res[, Z  := (obs_changes.x - null_mean) / pmax(null_sd, 1e-12)]




recur_thr <- 2L
cand <- res[n_cols >= recur_thr & PFAM_primary != "hypothetical protein"]

# rank by q then FE
setorder(cand, q_BH, -FE)









plot_dt <- res[!is.na(null_mean)]
plot_dt[, recurrent := (n_cols >= recur_thr)]

pA <- ggplot(plot_dt, aes(x = null_mean, y = obs_changes.x)) +
  geom_point(aes(alpha = recurrent), size = 1.8) +
  scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.9), guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "Expected # 0D changes per PFAM under Null 1 (log10)",
    y = "Observed # 0D changes per PFAM (log10)",
    title = "Observed vs expected 0D changes per PFAM",
    subtitle = glue("Recurrent PFAMs: ≥{recur_thr} colonies highlighted")
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "FigA_obs_vs_null1_expected_perPFAM.png"),
       pA, width = 7, height = 6, dpi = 300)



pA





topN <- 20
top <- cand[1:min(topN, .N)]
top[, PFAM_primary := factor(PFAM_primary, levels = rev(PFAM_primary))]

pB <- ggplot(top, aes(x = PFAM_primary, y = obs_changes.x)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL, y = "# 0D SNV changes",
    title = glue("Top recurrent PFAMs (≥{recur_thr} colonies)"),
    subtitle = "Ranked by q_BH then enrichment"
  ) +
  theme_classic()
pB













copy(res)[!is.na(null_mean) & null_mean > 0]
plot_dt[, log2FE := log2( (obs_changes.x + 1) / (null_mean + 1) )]   # +1 stabilizes tiny means
plot_dt[, recurrent := (n_cols >= 4)]

# take top N by obs (or by q)
N <- 200
top <- plot_dt[order(-obs_changes.x)][1:min(N, .N)]
top[, PFAM_rank := .I]

p <- ggplot(top, aes(x = PFAM_rank, y = log2FE)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(size = obs_changes.x, color = n_cols, alpha = recurrent)) +
  scale_alpha_manual(values = c(`FALSE`=0.3, `TRUE`=0.95), guide="none") +
  scale_size_continuous(range=c(1.2, 5)) +
  labs(x = "PFAMs (ranked by observed # changes)",
       y = "log2 fold-enrichment vs Null 1  ( (obs+1)/(E+1) )",
       color = "# colonies",
       size = "# changes",
       title = "PFAM classes with excess 0D changes under Null 1") +
  theme_classic()

p









top <- res[!is.na(null_mean)][order(-(obs_changes.x - null_mean))][1:30]
top[, PFAM_primary := factor(PFAM_primary, levels = rev(PFAM_primary))]
top[, delta := obs_changes.x - null_mean]

p <- ggplot(top, aes(x = PFAM_primary, y = obs_changes.x)) +
  geom_col() +
  geom_point(aes(y = null_mean), shape = 21, fill = "white", size = 2) +
  geom_errorbar(aes(ymin = pmax(null_mean - 1.96*null_sd, 0),
                    ymax = null_mean + 1.96*null_sd),
                width = 0.2) +
  coord_flip() +
  labs(x=NULL,
       y="# 0D changes",
       title="Top PFAMs: observed vs Null 1 expectation",
       subtitle="Bars = observed; open dot = null mean; whiskers = ±1.96 SD") +
  theme_classic()

p












plot_dt <- res[!is.na(obs_changes.x)]
plot_dt[, recurrent := (n_cols >= 4)]

p <- ggplot(plot_dt, aes(x = n_cols, y = obs_changes.x)) +
  geom_point(aes(alpha = recurrent), size = 2) +
  scale_alpha_manual(values=c(`FALSE`=0.25, `TRUE`=0.9), guide="none") +
  scale_y_continuous(trans="log10") +
  labs(x = "# colonies with ≥1 change",
       y = "# 0D changes (log10)",
       title="Recurrent PFAMs accrue many 0D changes",
       subtitle="PFAMs with changes in ≥4 colonies highlighted") +
  theme_classic()

p





d1 <- res[!is.na(null_mean), .(value = null_mean, type = "Null expected")]
d2 <- res[!is.na(obs_changes.x), .(value = obs_changes.x, type = "Observed")]
dd <- rbind(d1, d2)

p <- ggplot(dd, aes(x = value)) +
  geom_histogram(bins=60) +
  facet_wrap(~type, scales="free_y") +
  scale_x_continuous(trans="log10") +
  theme_classic() +
  labs(x="# changes per PFAM (log10)", y="Count of PFAMs",
       title="Null expects <1 change per PFAM; observed has heavy tail")

p












# ----------------------------
# 1) Build tail over PFAMs
# ----------------------------
tail_pfam <- function(x) {
  # x is integer vector of obs_changes per PFAM (including zeros if you want)
  x <- as.integer(x)
  x <- x[is.finite(x)]
  kmax <- max(x)
  data.table(
    k = 1:kmax,
    S = vapply(1:kmax, function(kk) mean(x >= kk), numeric(1))
  )
}

# If you want ONLY PFAMs that were ever hit (exclude zeros), use res[obs_changes>0]
# If you want all PFAMs in universe (include zeros), use res (with zeros present).
# I recommend "hit PFAMs" for a cleaner tail, like your original.
tail_dt <- tail_pfam(res[obs_changes.x > 0]$obs_changes.x)

# ----------------------------
# 2) Pick top 20 PFAMs with >=4 colonies
# ----------------------------
topN <- 20
top_pfam <- res[n_cols >= 2 & obs_changes.x > 0][order(-obs_changes.x)][1:min(topN, .N)]
top_pfam[, PFAM_label := PFAM_primary]

# Map each PFAM’s y-position as S(k=obs_changes)
setkey(tail_dt, k)
top_pfam[, k_plot := pmin(obs_changes.x, max(tail_dt$k))]
top_pfam <- tail_dt[top_pfam, on = .(k = k_plot)]
# top_pfam now has: PFAM_primary, obs_changes, n_cols, S
top_pfam
# ----------------------------
# 3) Plot
# ----------------------------
p_tail <- ggplot(tail_dt, aes(x = k, y = S)) +
  geom_step(linewidth = 1, color = "black") +
  geom_point(
    data = top_pfam,
    aes(x = obs_changes.x, y = S),
    inherit.aes = FALSE,
    size = 2.6,
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = top_pfam,
    aes(x = obs_changes.x, y = S, label = PFAM_label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  scale_x_continuous(trans="log10")+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "k (0D changes per PFAM)",
    y = "Fraction of PFAMs with ≥ k changes",
    title = "Tail of 0D changes per PFAM",
    subtitle = "Red points = top 20 PFAMs with ≥4 colonies"
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "tail_pfam_top20_recurrent_ge4colonies.png"),
       p_tail, width = 8, height = 5.5, dpi = 300)

p_tail

test <- gene_clean[gene_clean$pN_pS_ratio>1,]



gene_clean






g_key <- as.data.table(gene_clean)[
  , .(
    MAG_id = as.character(MAG_id),
    Colony = as.character(Colony),
    Month  = as.character(Month),
    corresponding_gene_call = as.character(corresponding_gene_call),
    PFAM_primary = as.character(PFAM_primary),
    Genus = as.character(Genus),
    w = as.numeric(get(WCOL))
  )
]

g_key <- g_key[
  !is.na(PFAM_primary) & PFAM_primary != "" &
    is.finite(w) & w > 0
]

g_key <- unique(
  g_key,
  by = c("MAG_id","Colony","Month","corresponding_gene_call","PFAM_primary","Genus")
)

setkey(g_key, MAG_id, Colony, Month)



unit_univ_list <- vector("list", nrow(unit_meta))
names(unit_univ_list) <- unit_meta$unit_id

for (i in seq_len(nrow(unit_meta))) {
  ui  <- unit_meta$unit_id[i]
  mag <- unit_meta$MAG_id[i]
  col <- unit_meta$Colony[i]
  m1  <- unit_meta$m1[i]
  m2  <- unit_meta$m2[i]
  
  dt1 <- g_key[J(mag, col, m1), nomatch = 0L]
  dt2 <- g_key[J(mag, col, m2), nomatch = 0L]
  u   <- unique(rbind(dt1, dt2))
  
  unit_univ_list[[ui]] <- u[, .(
    corresponding_gene_call,
    PFAM_primary,
    MAG_id,
    Colony,
    Genus,
    w
  )]
}




obs_pfam_breadth <- changes[
  !is.na(PFAM_primary) & PFAM_primary != "",
  .(
    obs_changes = .N,
    obs_n_units = uniqueN(unit_id),
    obs_n_cols  = uniqueN(colony_id),
    obs_n_MAGs  = uniqueN(MAG_id),
    obs_n_genus = uniqueN(Genus)
  ),
  by = PFAM_primary
]


bootstrap_null_breadth <- function(
    B = 500,
    seed = 1L,
    unit_meta,
    n_by_unit,
    unit_univ_list,
    pfams
) {
  set.seed(seed)
  
  n_pf <- length(pfams)
  pfam_index <- setNames(seq_along(pfams), pfams)
  
  # unit -> n_changes
  nmap <- n_by_unit$n_changes
  names(nmap) <- n_by_unit$unit_id
  
  # store matrices
  null_changes <- matrix(0L, nrow = B, ncol = n_pf, dimnames = list(NULL, pfams))
  null_n_units <- matrix(0L, nrow = B, ncol = n_pf, dimnames = list(NULL, pfams))
  null_n_cols  <- matrix(0L, nrow = B, ncol = n_pf, dimnames = list(NULL, pfams))
  null_n_MAGs  <- matrix(0L, nrow = B, ncol = n_pf, dimnames = list(NULL, pfams))
  null_n_genus <- matrix(0L, nrow = B, ncol = n_pf, dimnames = list(NULL, pfams))
  
  for (b in seq_len(B)) {
    sim_events <- vector("list", nrow(unit_meta))
    
    for (i in seq_len(nrow(unit_meta))) {
      ui <- unit_meta$unit_id[i]
      n_u <- nmap[[ui]]
      if (is.null(n_u) || n_u <= 0) next
      
      U <- unit_univ_list[[ui]]
      if (is.null(U) || nrow(U) == 0) next
      
      idx <- sample.int(nrow(U), size = n_u, replace = TRUE, prob = U$w)
      samp <- copy(U[idx])
      samp[, unit_id := ui]
      sim_events[[i]] <- samp[, .(PFAM_primary, unit_id, MAG_id, Colony, Genus)]
    }
    
    sim_dt <- rbindlist(sim_events, use.names = TRUE, fill = TRUE)
    if (nrow(sim_dt) == 0) next
    
    sim_sum <- sim_dt[, .(
      obs_changes = .N,
      n_units = uniqueN(unit_id),
      n_cols  = uniqueN(Colony),
      n_MAGs  = uniqueN(MAG_id),
      n_genus = uniqueN(Genus)
    ), by = PFAM_primary]
    
    jj <- pfam_index[sim_sum$PFAM_primary]
    ok <- !is.na(jj)
    
    null_changes[b, jj[ok]] <- sim_sum$obs_changes[ok]
    null_n_units[b, jj[ok]] <- sim_sum$n_units[ok]
    null_n_cols[b, jj[ok]]  <- sim_sum$n_cols[ok]
    null_n_MAGs[b, jj[ok]]  <- sim_sum$n_MAGs[ok]
    null_n_genus[b, jj[ok]] <- sim_sum$n_genus[ok]
  }
  
  list(
    null_changes = null_changes,
    null_n_units = null_n_units,
    null_n_cols  = null_n_cols,
    null_n_MAGs  = null_n_MAGs,
    null_n_genus = null_n_genus
  )
}

pfams <- sort(unique(c(
  obs_pfam_breadth$PFAM_primary,
  unique(unlist(lapply(unit_univ_list, function(u) u$PFAM_primary)))
)))

null_breadth <- bootstrap_null_breadth(
  B = 500,
  seed = 1L,
  unit_meta = unit_meta,
  n_by_unit = n_by_unit,
  unit_univ_list = unit_univ_list,
  pfams = pfams
)

tail_from_counts <- function(x, K = NULL) {
  x <- as.integer(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) stop("No finite values in x")
  if (is.null(K)) K <- seq_len(max(x))
  
  data.table(
    k = K,
    S = vapply(K, function(kk) mean(x >= kk), numeric(1))
  )
}


obs_counts <- obs_pfam_breadth$obs_changes
K_changes <- seq_len(max(obs_counts))

obs_tail_changes <- tail_from_counts(obs_counts, K = K_changes)

null_tail_changes <- rbindlist(
  lapply(seq_len(nrow(null_breadth$null_changes)), function(b) {
    tt <- tail_from_counts(null_breadth$null_changes[b, ], K = K_changes)
    tt[, b := b]
    tt
  })
)

band_changes <- null_tail_changes[, .(
  lo  = quantile(S, 0.025, na.rm = TRUE),
  hi  = quantile(S, 0.975, na.rm = TRUE),
  mid = median(S, na.rm = TRUE)
), by = k]



topN <- 20
top_pfam <- obs_pfam_breadth[obs_n_cols >= 2 & obs_changes > 0][order(-obs_changes)][1:min(topN, .N)]
top_pfam[, PFAM_label := PFAM_primary]

setkey(obs_tail_changes, k)
top_pfam[, k_plot := pmin(obs_changes, max(obs_tail_changes$k))]
top_pfam <- obs_tail_changes[top_pfam, on = .(k = k_plot)]

p_tail_null <- ggplot() +
  geom_ribbon(data = band_changes, aes(x = k, ymin = lo, ymax = hi), alpha = 0.20) +
  geom_line(data = band_changes, aes(x = k, y = mid), linewidth = 0.8, linetype = 2) +
  geom_step(data = obs_tail_changes, aes(x = k, y = S), linewidth = 1, color = "black") +
  geom_point(
    data = top_pfam,
    aes(x = obs_changes, y = S),
    inherit.aes = FALSE,
    size = 2.6,
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = top_pfam,
    aes(x = obs_changes, y = S, label = PFAM_label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "k (0D changes per PFAM)",
    y = "Fraction of PFAMs with ≥ k changes",
    title = "Tail of 0D changes per PFAM",
    subtitle = "Black = observed; dashed line/ribbon = unit-aware null; red = top PFAMs"
  ) +
  theme_classic()

p_tail_null




obs_cols <- obs_pfam_breadth$obs_n_cols
K_cols <- seq_len(max(obs_cols))

obs_tail_cols <- tail_from_counts(obs_cols, K = K_cols)

null_tail_cols <- rbindlist(
  lapply(seq_len(nrow(null_breadth$null_n_cols)), function(b) {
    tt <- tail_from_counts(null_breadth$null_n_cols[b, ], K = K_cols)
    tt[, b := b]
    tt
  })
)

band_cols <- null_tail_cols[, .(
  lo  = quantile(S, 0.025, na.rm = TRUE),
  hi  = quantile(S, 0.975, na.rm = TRUE),
  mid = median(S, na.rm = TRUE)
), by = k]

p_cols_tail <- ggplot() +
  geom_ribbon(data = band_cols, aes(x = k, ymin = lo, ymax = hi), alpha = 0.20) +
  geom_line(data = band_cols, aes(x = k, y = mid), linewidth = 0.8, linetype = 2) +
  geom_step(data = obs_tail_cols, aes(x = k, y = S), linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(
    x = "k (colonies with ≥1 0D sweep hit in PFAM)",
    y = "Fraction of PFAMs with ≥ k colonies",
    title = "Breadth of PFAM sharing across colonies",
    subtitle = "Observed vs unit-aware null"
  )

p_cols_tail






obs_genus <- obs_pfam_breadth$obs_n_genus
K_genus <- seq_len(max(obs_genus))

obs_tail_genus <- tail_from_counts(obs_genus, K = K_genus)

null_tail_genus <- rbindlist(
  lapply(seq_len(nrow(null_breadth$null_n_genus)), function(b) {
    tt <- tail_from_counts(null_breadth$null_n_genus[b, ], K = K_genus)
    tt[, b := b]
    tt
  })
)

band_genus <- null_tail_genus[, .(
  lo  = quantile(S, 0.025, na.rm = TRUE),
  hi  = quantile(S, 0.975, na.rm = TRUE),
  mid = median(S, na.rm = TRUE)
), by = k]

p_genus_tail <- ggplot() +
  geom_ribbon(data = band_genus, aes(x = k, ymin = lo, ymax = hi), alpha = 0.20) +
  geom_line(data = band_genus, aes(x = k, y = mid), linewidth = 0.8, linetype = 2) +
  geom_step(data = obs_tail_genus, aes(x = k, y = S), linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  labs(
    x = "k (genera with ≥1 0D sweep hit in PFAM)",
    y = "Fraction of PFAMs with ≥ k genera",
    title = "Breadth of PFAM sharing across genera",
    subtitle = "Observed vs unit-aware null"
  )

p_genus_tail



pfam_index <- setNames(seq_along(pfams), pfams)

obs_aligned <- data.table(PFAM_primary = pfams)
obs_aligned <- merge(obs_aligned, obs_pfam_breadth, by = "PFAM_primary", all.x = TRUE)

for (v in c("obs_changes","obs_n_units","obs_n_cols","obs_n_MAGs","obs_n_genus")) {
  obs_aligned[is.na(get(v)), (v) := 0L]
}

p_cols <- (colSums(null_breadth$null_n_cols >= matrix(obs_aligned$obs_n_cols,
                                                      nrow = nrow(null_breadth$null_n_cols),
                                                      ncol = length(pfams),
                                                      byrow = TRUE)) + 1) /
  (nrow(null_breadth$null_n_cols) + 1)

p_genus <- (colSums(null_breadth$null_n_genus >= matrix(obs_aligned$obs_n_genus,
                                                        nrow = nrow(null_breadth$null_n_genus),
                                                        ncol = length(pfams),
                                                        byrow = TRUE)) + 1) /
  (nrow(null_breadth$null_n_genus) + 1)

breadth_res <- copy(obs_aligned)
breadth_res[, p_cols  := p_cols]
breadth_res[, q_cols  := p.adjust(p_cols, method = "BH")]
breadth_res[, p_genus := p_genus]
breadth_res[, q_genus := p.adjust(p_genus, method = "BH")]






changes[, month_pair := paste(Month_t1, Month_t2, sep = "->")]




obs_shared <- changes[
  !is.na(PFAM_primary) & PFAM_primary != "",
  .(
    n_events = .N,
    n_cols   = uniqueN(colony_id),
    n_MAGs   = uniqueN(MAG_id),
    n_genus  = uniqueN(Genus)
  ),
  by = .(PFAM_primary, month_pair)
]



obs_pfam_shared <- obs_shared[, .(
  T_genus2   = sum(n_genus >= 2),
  T_genus3   = sum(n_genus >= 3),
  T_colony2  = sum(n_cols  >= 2),
  T_colony3  = sum(n_cols  >= 3),
  max_genus  = max(n_genus),
  max_colony = max(n_cols)
), by = PFAM_primary]




bootstrap_null_shared <- function(
    B = 500,
    seed = 1L,
    unit_meta,
    n_by_unit,
    unit_univ_list,
    pfams
) {
  set.seed(seed)
  
  pfam_index <- setNames(seq_along(pfams), pfams)
  
  nmap <- n_by_unit$n_changes
  names(nmap) <- n_by_unit$unit_id
  
  null_T_genus2   <- matrix(0L, nrow = B, ncol = length(pfams), dimnames = list(NULL, pfams))
  null_T_colony2  <- matrix(0L, nrow = B, ncol = length(pfams), dimnames = list(NULL, pfams))
  null_max_genus  <- matrix(0L, nrow = B, ncol = length(pfams), dimnames = list(NULL, pfams))
  null_max_colony <- matrix(0L, nrow = B, ncol = length(pfams), dimnames = list(NULL, pfams))
  
  for (b in seq_len(B)) {
    sim_list <- vector("list", nrow(unit_meta))
    
    for (i in seq_len(nrow(unit_meta))) {
      ui <- unit_meta$unit_id[i]
      n_u <- nmap[[ui]]
      if (is.null(n_u) || n_u <= 0) next
      
      U <- unit_univ_list[[ui]]
      if (is.null(U) || nrow(U) == 0) next
      
      idx <- sample.int(nrow(U), size = n_u, replace = TRUE, prob = U$w)
      samp <- copy(U[idx])
      
      samp[, `:=`(
        unit_id    = ui,
        month_pair = paste(unit_meta$m1[i], unit_meta$m2[i], sep = "->")
      )]
      
      sim_list[[i]] <- samp[, .(PFAM_primary, month_pair, Colony, Genus)]
    }
    
    sim_dt <- rbindlist(sim_list, use.names = TRUE, fill = TRUE)
    if (nrow(sim_dt) == 0) next
    
    sim_shared <- sim_dt[, .(
      n_cols  = uniqueN(Colony),
      n_genus = uniqueN(Genus)
    ), by = .(PFAM_primary, month_pair)]
    
    sim_pfam <- sim_shared[, .(
      T_genus2   = sum(n_genus >= 2),
      T_colony2  = sum(n_cols  >= 2),
      max_genus  = max(n_genus),
      max_colony = max(n_cols)
    ), by = PFAM_primary]
    
    jj <- pfam_index[sim_pfam$PFAM_primary]
    ok <- !is.na(jj)
    
    null_T_genus2[b, jj[ok]]   <- sim_pfam$T_genus2[ok]
    null_T_colony2[b, jj[ok]]  <- sim_pfam$T_colony2[ok]
    null_max_genus[b, jj[ok]]  <- sim_pfam$max_genus[ok]
    null_max_colony[b, jj[ok]] <- sim_pfam$max_colony[ok]
  }
  
  list(
    null_T_genus2   = null_T_genus2,
    null_T_colony2  = null_T_colony2,
    null_max_genus  = null_max_genus,
    null_max_colony = null_max_colony
  )
}

pfams <- sort(unique(c(
  obs_pfam_shared$PFAM_primary,
  unique(unlist(lapply(unit_univ_list, function(u) u$PFAM_primary)))
)))

null_shared <- bootstrap_null_shared(
  B = 500,
  seed = 1L,
  unit_meta = unit_meta,
  n_by_unit = n_by_unit,
  unit_univ_list = unit_univ_list,
  pfams = pfams
)




obs_aligned <- data.table(PFAM_primary = pfams)
obs_aligned <- merge(obs_aligned, obs_pfam_shared, by = "PFAM_primary", all.x = TRUE)

for (v in c("T_genus2","T_colony2","max_genus","max_colony")) {
  obs_aligned[is.na(get(v)), (v) := 0L]
}

p_T_genus2 <- (colSums(null_shared$null_T_genus2 >= matrix(
  obs_aligned$T_genus2,
  nrow = nrow(null_shared$null_T_genus2),
  ncol = length(pfams),
  byrow = TRUE
)) + 1) / (nrow(null_shared$null_T_genus2) + 1)

p_T_colony2 <- (colSums(null_shared$null_T_colony2 >= matrix(
  obs_aligned$T_colony2,
  nrow = nrow(null_shared$null_T_colony2),
  ncol = length(pfams),
  byrow = TRUE
)) + 1) / (nrow(null_shared$null_T_colony2) + 1)

shared_res <- copy(obs_aligned)
shared_res[, p_T_genus2  := p_T_genus2]
shared_res[, q_T_genus2  := p.adjust(p_T_genus2, method = "BH")]
shared_res[, p_T_colony2 := p_T_colony2]
shared_res[, q_T_colony2 := p.adjust(p_T_colony2, method = "BH")]


tail_from_counts <- function(x, K = NULL) {
  x <- as.integer(x)
  x <- x[is.finite(x)]
  if (is.null(K)) K <- 0:max(x)
  data.table(
    k = K,
    S = vapply(K, function(kk) mean(x >= kk), numeric(1))
  )
}
obs_tail <- tail_from_counts(shared_res$T_genus2)

null_tail <- rbindlist(lapply(seq_len(nrow(null_shared$null_T_genus2)), function(b) {
  tt <- tail_from_counts(null_shared$null_T_genus2[b, ], K = obs_tail$k)
  tt[, b := b]
  tt
}))

band <- null_tail[, .(
  lo  = quantile(S, 0.025),
  hi  = quantile(S, 0.975),
  mid = median(S)
), by = k]

ggplot() +
  geom_ribbon(data = band, aes(x = k, ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(data = band, aes(x = k, y = mid), linetype = 2, linewidth = 0.8) +
  geom_step(data = obs_tail, aes(x = k, y = S), linewidth = 1) +
  theme_classic() +
  labs(
    x = "k (month pairs where PFAM is hit in >=2 genera)",
    y = "Fraction of PFAMs with >= k shared month pairs",
    title = "Same-month-pair shared 0D sweeps across genera",
    subtitle = "Observed vs unit-aware null"
  )















































































g <- copy(gene_clean)
g[, Colony := as.character(Colony)]   # you already have this
g[, Month  := as.character(Month)]

# Filters to stabilize pN/pS
g <- g[
  !is.na(PFAM_primary) & PFAM_primary != "" &
    is.finite(pN_pS_ratio) &
    is.finite(nN_gene_reference) & is.finite(nS_gene_reference) &
    nN_gene_reference > 0 & nS_gene_reference > 0
]

# Binary signal and a continuous transform
eps <- 1e-6
g[, pnps_gt1 := as.integer(pN_pS_ratio > 1)]
g[, log_pnps := log(pN_pS_ratio + eps)]

setkey(g, MAG_id, Colony, Month)





# unit_meta must have: unit_id, MAG_id, Colony, m1, m2
stopifnot(all(c("unit_id","MAG_id","Colony","m1","m2") %in% names(unit_meta)))

unit_univ_list <- vector("list", nrow(unit_meta))
names(unit_univ_list) <- unit_meta$unit_id

for (i in seq_len(nrow(unit_meta))) {
  ui  <- unit_meta$unit_id[i]
  mag <- unit_meta$MAG_id[i]
  col <- unit_meta$Colony[i]
  m1  <- unit_meta$m1[i]
  m2  <- unit_meta$m2[i]
  
  dt1 <- g[J(mag, col, m1), nomatch = 0L]
  dt2 <- g[J(mag, col, m2), nomatch = 0L]
  u   <- unique(rbind(dt1, dt2), by = c("corresponding_gene_call","MAG_id","Colony"))
  
  # Keep what we need for null draws
  unit_univ_list[[ui]] <- u[, .(
    corresponding_gene_call,
    PFAM_primary,
    pnps_gt1,
    log_pnps,
    w = pmax(nS_gene_reference, 1)   # weight option (see note below)
  )]
}

cat("Units:", length(unit_univ_list), "\n")
cat("Median genes per unit universe:", median(vapply(unit_univ_list, nrow, integer(1))), "\n")







obs_pfam_pnps <- g[, .(
  obs_n        = .N,
  obs_gt1      = sum(pnps_gt1, na.rm = TRUE),
  obs_mean_log = mean(log_pnps, na.rm = TRUE),
  n_cols       = uniqueN(Colony),
  n_MAGs       = uniqueN(MAG_id),
  n_genus      = uniqueN(Genus)
), by = PFAM_primary]






# observed K_u per unit: count of pnps_gt1 within the unit's observed universe
# (but we need unit-specific observed set; easiest: compute from your changes table if unit_id exists in gene_clean;
# if unit_id is not in gene_clean, compute K_u from dt1/dt2 union inside the loop below)

# Precompute K_u by reusing unit_univ_list:
n_by_unit <- data.table(
  unit_id = unit_meta$unit_id,
  K_u = vapply(unit_univ_list, function(U) sum(U$pnps_gt1, na.rm = TRUE), integer(1))
)

# PFAM universe
pfams <- sort(unique(c(obs_pfam_pnps$PFAM_primary,
                       unique(unlist(lapply(unit_univ_list, function(u) u$PFAM_primary))))))

pfam_index <- setNames(seq_along(pfams), pfams)

# Observed PFAM counts
obs_vec <- integer(length(pfams))
tmp <- obs_pfam_pnps[, .(PFAM_primary, obs_gt1)]
obs_vec[pfam_index[tmp$PFAM_primary]] <- tmp$obs_gt1

bootstrap_null_pnps_gt1 <- function(B = 2000, seed = 1L) {
  set.seed(seed)
  null_mat <- matrix(0L, nrow = B, ncol = length(pfams))
  colnames(null_mat) <- pfams
  
  Kmap <- n_by_unit$K_u
  names(Kmap) <- n_by_unit$unit_id
  
  for (b in seq_len(B)) {
    counts <- integer(length(pfams))
    
    for (i in seq_len(nrow(unit_meta))) {
      ui <- unit_meta$unit_id[i]
      K  <- Kmap[[ui]]
      if (is.null(K) || K <= 0) next
      
      U <- unit_univ_list[[ui]]
      if (is.null(U) || nrow(U) == 0) next
      
      # sample K genes from the unit universe (no replacement is fine here)
      idx <- sample.int(nrow(U), size = min(K, nrow(U)), replace = FALSE, prob = U$w)
      pf <- U$PFAM_primary[idx]
      tab <- table(pf)
      
      jj <- pfam_index[names(tab)]
      ok <- !is.na(jj)
      if (any(ok)) counts[jj[ok]] <- counts[jj[ok]] + as.integer(tab[ok])
    }
    
    null_mat[b, ] <- counts
  }
  
  ge <- colSums(null_mat >= matrix(obs_vec, nrow = B, ncol = length(obs_vec), byrow = TRUE))
  p_emp <- (ge + 1) / (B + 1)
  
  data.table(
    PFAM_primary = pfams,
    obs_gt1      = obs_vec,
    null_mean    = colMeans(null_mat),
    null_sd      = apply(null_mat, 2, sd),
    p_perm       = p_emp,
    q_BH         = p.adjust(p_emp, method = "BH")
  )
}

res_pnps_gt1 <- bootstrap_null_pnps_gt1(B = 2000, seed = 1)
res_pnps_gt1 <- merge(res_pnps_gt1, obs_pfam_pnps, by = "PFAM_primary", all.x = TRUE)
res_pnps_gt1[, FE := obs_gt1 / pmax(null_mean, 1e-12)]
res_pnps_gt1[, Z  := (obs_gt1 - null_mean) / pmax(null_sd, 1e-12)]



res_pnps_gt1







#CLEANING UP BEFORE MERGE

pnps_clean <- unique(res_pnps_gt1[, .(
  PFAM_primary,
  pnps_obs_gt1   = obs_gt1.x,   # use the enrichment obs count
  pnps_null_mean = null_mean,
  pnps_null_sd   = null_sd,
  pnps_p         = p_perm,
  pnps_q         = q_BH,
  pnps_FE        = FE,
  pnps_Z         = Z,
  pnps_obs_n     = obs_n.x,
  pnps_mean_log  = obs_mean_log.x,
  pnps_n_cols    = n_cols.x,
  pnps_n_MAGs    = n_MAGs.x,
  pnps_n_genus   = n_genus.x
)], by = "PFAM_primary")

sweep_clean <- unique(top_pfam[, .(
  PFAM_primary,
  sweep_obs_changes = obs_changes.x,  # or obs_changes.y; here they match in your output
  sweep_null_mean   = null_mean,
  sweep_null_sd     = null_sd,
  sweep_p           = p_null,
  sweep_q           = q_BH,
  sweep_FE          = FE,
  sweep_Z           = Z,
  sweep_n_units     = n_units,
  sweep_n_cols      = n_cols,
  sweep_n_MAGs      = n_MAGs,
  sweep_n_genus     = n_genus,
  tail_k            = k,
  tail_S            = S
)], by = "PFAM_primary")





combo <- merge(sweep_clean, pnps_clean, by = "PFAM_primary", all.x = TRUE)

# Rank by sweep significance and show pN/pS enrichment alongside
combo_ranked <- combo[order(sweep_q, pnps_q, -sweep_obs_changes)]

combo_ranked[, .(
  PFAM_primary,
  sweep_obs_changes, sweep_FE, sweep_q, sweep_n_cols, sweep_n_genus,
  pnps_obs_gt1 = pnps_obs_gt1, pnps_FE, pnps_q, pnps_obs_n
)][1:30]



plot_dt <- combo[!is.na(pnps_FE)]
plot_dt[, log2FE_sweep := log2(sweep_FE)]
plot_dt[, log2FE_pnps  := log2(pnps_FE)]

# label candidates: significant in either axis
plot_dt[, label_me := (sweep_q <= 0.05 | pnps_q <= 0.05) & sweep_n_cols >= 4]

p <- ggplot(plot_dt, aes(x = log2FE_sweep, y = log2FE_pnps)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(size = sweep_obs_changes, color = sweep_n_cols), alpha = 0.75) +
  ggrepel::geom_text_repel(
    data = plot_dt[label_me == TRUE],
    aes(label = PFAM_primary),
    size = 3,
    max.overlaps = 40
  ) +
  scale_color_viridis_c(name = "# colonies\n(sweep)") +
  scale_size_continuous(name = "# 0D changes\n(sweep)", range = c(2, 10)) +
  labs(
    x = expression(log[2]*" fold-enrichment of 0D changes (Null 1)"),
    y = expression(log[2]*" fold-enrichment of pN/pS > 1 genes (Null 1)"),
    title = "Do parallel-sweep PFAMs also show excess high pN/pS genes?"
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "sweep_vs_pnps_enrichment_scatter.png"),
       p, width = 7.5, height = 6, dpi = 300)



p



combo[, pnps_rate_gt1 := pnps_obs_gt1 / pmax(pnps_obs_n, 1L)]
combo[order(-pnps_rate_gt1)][1:20, .(PFAM_primary, pnps_obs_gt1, pnps_obs_n, pnps_rate_gt1, pnps_q)]






#weird assumptions in parallel sweep plot, so let's improve it in null 1. 

#need to add sifnificance scaling or indicator
#need to also add count-scale floor on the null to prevent odd outliers






#---------------------------
# 1) Standardize column names
#---------------------------
sweep_tab <- as.data.table(top_pfam)
pnps_tab  <- as.data.table(res_pnps_gt1)

# Ensure consistent names (edit these if your tables differ)
if ("obs_changes.x" %in% names(sweep_tab)) setnames(sweep_tab, "obs_changes.x", "obs_changes")
if ("null_mean" %in% names(sweep_tab) == FALSE && "null1_mean" %in% names(sweep_tab)) setnames(sweep_tab, "null1_mean", "null_mean")
if ("q_BH" %in% names(sweep_tab) == FALSE && "q_null1" %in% names(sweep_tab)) setnames(sweep_tab, "q_null1", "q_BH")

if ("obs_gt1.x" %in% names(pnps_tab)) setnames(pnps_tab, "obs_gt1.x", "obs_gt1")
if ("null_mean" %in% names(pnps_tab) == FALSE && "null_mean" %in% names(pnps_tab)) {} # ok
if ("q_BH" %in% names(pnps_tab) == FALSE && "q_null1" %in% names(pnps_tab)) setnames(pnps_tab, "q_null1", "q_BH")


# ---- build sweep_keep with required cols ----
sweep_keep <- sweep_tab[, .(
  PFAM_primary,
  sweep_obs  = obs_changes,
  sweep_null = null_mean,
  sweep_q    = q_BH
)]

# add optional cols if present
if ("n_cols" %in% names(sweep_tab))  sweep_keep[, sweep_n_cols  := sweep_tab$n_cols]
if ("n_units" %in% names(sweep_tab)) sweep_keep[, sweep_n_units := sweep_tab$n_units]
if (!("n_cols" %in% names(sweep_tab)))  sweep_keep[, sweep_n_cols := NA_integer_]
if (!("n_units" %in% names(sweep_tab))) sweep_keep[, sweep_n_units := NA_integer_]

# ---- build pnps_keep with required cols ----
pnps_keep <- pnps_tab[, .(
  PFAM_primary,
  pnps_obs_gt1 = obs_gt1,
  pnps_null    = null_mean,
  pnps_q       = q_BH
)]



# ---- merge ----
combo <- merge(sweep_keep, pnps_keep, by = "PFAM_primary", all = TRUE)

# sanity check
print(combo[1:5])

#---------------------------
# 2) FE floor for plotting (NOT for stats)
#---------------------------
# Count-scale floor: avoids absurd log2FE from null_mean ~ 0.
# For counts, floor at 1 is usually cleanest and easiest to explain in figure caption.
FE_FLOOR <- 1

combo[, sweep_FE_plot := sweep_obs / pmax(sweep_null, FE_FLOOR)]
combo[, pnps_FE_plot  := pnps_obs_gt1 / pmax(pnps_null, FE_FLOOR)]

combo[, log2_sweep_FE := log2(sweep_FE_plot)]
combo[, log2_pnps_FE  := log2(pnps_FE_plot)]

#---------------------------
# 3) Significance encoding
#---------------------------
ALPHA_Q <- 0.10
combo[, sig_sweep := !is.na(sweep_q) & sweep_q <= ALPHA_Q]
combo[, sig_pnps  := !is.na(pnps_q)  & pnps_q  <= ALPHA_Q]

# Plotting aesthetics
combo[, pt_alpha := ifelse(sig_sweep | sig_pnps, 0.85, 0.25)]
combo[, pt_shape := ifelse(sig_sweep | sig_pnps, 21, 21)]   # use shape 21 so fill works
combo[, pt_fill  := ifelse(sig_sweep & sig_pnps, "black",
                           ifelse(sig_sweep, "grey20",
                                  ifelse(sig_pnps, "grey50", "white")))]
combo[, pt_color := ifelse(sig_sweep | sig_pnps, "black", "grey70")]

# Label set: only double-sig, or customize
label_dt <- combo[sig_sweep & sig_pnps][
  order(-sweep_obs)][1:15]   # label top by sweep_obs, tweak

#---------------------------
# 4) Plot
#---------------------------
p <- ggplot(combo, aes(x = log2_sweep_FE, y = log2_pnps_FE)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(size = sweep_obs,
                 alpha = pt_alpha),
             shape = 21,
             stroke = 0.4,
             fill = combo$pt_fill,
             color = combo$pt_color) +
  ggrepel::geom_text_repel(
    data = label_dt,
    aes(label = PFAM_primary),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_alpha_identity() +
  scale_size_continuous(name = "# 0D changes") +
  labs(
    x = glue::glue("log2 fold-enrichment of 0D changes (Null 1; FE floor={FE_FLOOR})"),
    y = glue::glue("log2 fold-enrichment of pN/pS>1 genes (Null 1; FE floor={FE_FLOOR})"),
    title = "Do parallel-sweep PFAMs also show excess high pN/pS genes?",
    subtitle = glue::glue("Filled points indicate q≤{ALPHA_Q} for sweep and/or pN/pS enrichment; open points are not significant")
  ) +
  theme_classic()

print(p)








#xpand to all well-defined
#by this we mean, sweep_null is not tiny ( > 0.5)
#pnps_null is also not tiny
#both values are finite

combo_all <- merge(sweep_keep, pnps_keep, by="PFAM_primary", all=TRUE)

# fold enrichments (add tiny epsilon to avoid log2(0); better is to filter)
eps <- 1e-9
combo_all[, sweep_FE := sweep_obs / pmax(sweep_null, eps)]
combo_all[, pnps_FE  := pnps_obs_gt1 / pmax(pnps_null, eps)]

combo_all[, x := log2(sweep_FE)]
combo_all[, y := log2(pnps_FE)]

# define significance flags (you can tune threshold)
combo_all[, sig_any := (sweep_q <= 0.1) | (pnps_q <= 0.1)]
combo_all[, sig_both := (sweep_q <= 0.1) & (pnps_q <= 0.1)]

# filter to “stable” region (tune cutoffs)
plot_dt <- combo_all[
  is.finite(x) & is.finite(y) &
    !is.na(sweep_null) & !is.na(pnps_null) &
    sweep_null >= 1 & pnps_null >= 0.5
]

# label only the most extreme positive pnps_FE (or both)
lab_dt <- plot_dt[order(-y)][1:10]

ggplot(plot_dt, aes(x=x, y=y)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(aes(size = sweep_obs, alpha = sig_any, shape = sig_any),
             color="black") +
  scale_alpha_manual(values=c(`FALSE`=0.25, `TRUE`=0.8), guide="none") +
  scale_shape_manual(values=c(`FALSE`=1, `TRUE`=16), guide="none") +
  geom_text_repel(data=lab_dt, aes(label=PFAM_primary), size=3, max.overlaps=50) +
  labs(
    x = "log2 fold-enrichment of 0D changes (Null 1)",
    y = "log2 fold-enrichment of pN/pS > 1 genes (Null 1)",
    size = "# 0D changes"
  ) +
  theme_classic()

nrow(sweep_keep)
nrow(pnps_keep)








#--------- helper: rename first column that exists ----------
  rename_first_present <- function(dt, candidates, newname) {
    cand <- candidates[candidates %in% names(dt)]
    if (length(cand) == 0L) stop("None of these columns exist: ", paste(candidates, collapse=", "))
    setnames(dt, cand[1], newname)
    invisible(dt)
  }

# ---------- build sweep_tab_full ----------
sweep_tab_full <- merge(
  copy(res_null),      # PFAM_primary, obs_changes, null_mean, null_sd, p_null, q_BH  (but may be suffixed by merge)
  copy(obs_pfam),      # PFAM_primary, obs_changes, n_units, n_cols, n_MAGs, n_genus
  by = "PFAM_primary",
  all.x = TRUE
)

# --- rename obs changes (prefer the "null table" one) ---
# after merge it might be obs_changes.x / obs_changes.y, or obs_changes_null/obs_changes_obs, etc.
rename_first_present(
  sweep_tab_full,
  c("obs_changes", "obs_changes.x", "obs_changes_null", "obs_changes_res_null", "obs_changes_res"),
  "sweep_obs_changes"
)

# --- rename null summary columns ---
rename_first_present(sweep_tab_full, c("null_mean", "null_mean.x", "null_mean_null"), "sweep_null_mean")
rename_first_present(sweep_tab_full, c("null_sd",   "null_sd.x",   "null_sd_null"),   "sweep_null_sd")
rename_first_present(sweep_tab_full, c("p_null",    "p_null.x",    "p_null_null"),    "sweep_p")
rename_first_present(sweep_tab_full, c("q_BH",      "q_BH.x",      "q_BH_null"),      "sweep_q")

# --- rename metadata count columns if present ---
if ("n_units" %in% names(sweep_tab_full)) setnames(sweep_tab_full, "n_units", "sweep_n_units")
if ("n_cols"  %in% names(sweep_tab_full)) setnames(sweep_tab_full, "n_cols",  "sweep_n_cols")
if ("n_MAGs"  %in% names(sweep_tab_full)) setnames(sweep_tab_full, "n_MAGs",  "sweep_n_MAGs")
if ("n_genus" %in% names(sweep_tab_full)) setnames(sweep_tab_full, "n_genus", "sweep_n_genus")

# If merge created obs_changes.y etc, keep it around but don't use it.
# (Optional) sanity check:
# names(sweep_tab_full)

# --- derived stats ---
eps <- 1e-9
sweep_tab_full[, sweep_FE := sweep_obs_changes / pmax(sweep_null_mean, eps)]
sweep_tab_full[, sweep_Z  := (sweep_obs_changes - sweep_null_mean) / pmax(sweep_null_sd, eps)]

cat("PFAMs:", nrow(sweep_tab_full), "\n")
cat("Missing sweep_n_cols:", sum(is.na(sweep_tab_full$sweep_n_cols)), "\n")





#EXPANDING DATA


# ----- pick how many PFAMs to include in the enrichment scatter -----
TOPN <- 200  # try 200 first; 500 is fine too

setorder(sweep_tab_full, -sweep_obs_changes)
sweep_keep <- sweep_tab_full[1:TOPN]

# ----- tail curves (survival: fraction PFAMs with >= k changes) -----
# Observed tail over ALL PFAMs
obs_counts <- sweep_tab_full$sweep_obs_changes
k_grid_obs <- sort(unique(obs_counts[is.finite(obs_counts)]))

tail_df_obs <- data.table(
  k = k_grid_obs,
  S = sapply(k_grid_obs, function(kk) mean(obs_counts >= kk))
)

# Null "tail": using expected mean per PFAM (this is not a discrete count distribution,
# but it shows how concentrated the null expectation is vs observed)
null_mu <- sweep_tab_full$sweep_null_mean
k_grid_null <- sort(unique(round(null_mu, 6)))
tail_df_null <- data.table(
  k = k_grid_null,
  S = sapply(k_grid_null, function(kk) mean(null_mu >= kk))
)

# Labels: top 20 among those with >=4 colonies (if sweep_n_cols exists)
if ("sweep_n_cols" %in% names(sweep_tab_full)) {
  label_df <- sweep_tab_full[!is.na(sweep_n_cols) & sweep_n_cols >= 4]
  setorder(label_df, -sweep_obs_changes)
  label_df <- label_df[1:min(20, .N)]
} else {
  label_df <- sweep_tab_full[order(-sweep_obs_changes)][1:20]
}







p_tail <- ggplot() +
  geom_step(data = tail_df_obs,  aes(x = k, y = S), linewidth = 0.8) +
  geom_step(data = tail_df_null, aes(x = k, y = S), linewidth = 0.8, linetype = "dashed") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "k (# 0D changes per PFAM)",
    y = "Fraction of PFAMs with ≥ k changes",
    title = "Heavy tail of observed 0D changes vs Null 1 expectation",
    subtitle = "Solid = observed counts; dashed = PFAM-specific expected means under Null 1"
  ) +
  theme_classic()

print(p_tail)






sweep_keep[, rank_obs := frank(-sweep_obs_changes, ties.method = "first")]
sweep_keep[, log2FE := log2((sweep_obs_changes + 1) / (sweep_null_mean + 1))]

p_excess <- ggplot(sweep_keep, aes(x = rank_obs, y = log2FE)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(size = sweep_obs_changes,
                 alpha = ifelse(!is.na(sweep_n_cols) & sweep_n_cols >= 2, 0.9, 0.35)),
             show.legend = TRUE) +
  scale_alpha_identity() +
  labs(
    x = "PFAMs (ranked by observed # changes)",
    y = "log2 fold-enrichment vs Null 1  ( (obs+1)/(E+1) )",
    title = "PFAM classes with excess 0D changes under Null 1",
    subtitle = "More opaque points indicate PFAMs observed in ≥2 colonies"
  ) +
  theme_classic()

# label top 20 with >=4 colonies (if available)
if (nrow(label_df) > 0) {
  label_df2 <- merge(label_df[, .(PFAM_primary)], sweep_keep[, .(PFAM_primary, rank_obs, log2FE)], by="PFAM_primary", all.x=TRUE)
  p_excess <- p_excess +
    ggrepel::geom_text_repel(
      data = label_df2,
      aes(x = rank_obs, y = log2FE, label = PFAM_primary),
      size = 3,
      max.overlaps = Inf
    )
}

print(p_excess)




































# FULL sweep table from res (NOT just top 20)
sweep_tab_full <- copy(res)

obs_col <- if ("obs_changes" %in% names(sweep_tab_full)) {
  "obs_changes"
} else if ("obs_changes.x" %in% names(sweep_tab_full)) {
  "obs_changes.x"
} else {
  stop("Can't find obs_changes column in sweep_tab_full")
}



sweep_tab_full[, sweep_obs_changes := get(obs_col)]
sweep_tab_full[, sweep_null_mean  := null_mean]
sweep_tab_full[, sweep_null_sd    := null_sd]
sweep_tab_full[, sweep_p          := p_null]
sweep_tab_full[, sweep_q          := q_BH]


# carry optional metadata if present
for (v in c("n_units","n_cols","n_MAGs","n_genus")) {
  newv <- paste0("sweep_", v)
  if (v %in% names(sweep_tab_full)) {
    sweep_tab_full[, (newv) := get(v)]
  } else {
    sweep_tab_full[, (newv) := NA_integer_]
  }
}

# effect sizes
sweep_tab_full[, sweep_FE := sweep_obs_changes / pmax(sweep_null_mean, 1e-9)]
sweep_tab_full[, sweep_log2FE := log2((sweep_obs_changes + 1) / (sweep_null_mean + 1))]
sweep_tab_full[, sweep_Z := (sweep_obs_changes - sweep_null_mean) / pmax(sweep_null_sd, 1e-9)]

sweep_keep_full <- sweep_tab_full[, .(
  PFAM_primary,
  sweep_obs_changes, sweep_null_mean, sweep_null_sd,
  sweep_p, sweep_q,
  sweep_FE, sweep_log2FE, sweep_Z,
  sweep_n_units, sweep_n_cols, sweep_n_MAGs, sweep_n_genus
)]




pnps_tab_full <- copy(res_pnps_gt1)

obs_gt1_col <- if ("obs_gt1" %in% names(pnps_tab_full)) {
  "obs_gt1"
} else if ("obs_gt1.x" %in% names(pnps_tab_full)) {
  "obs_gt1.x"
} else if ("obs_gt1.y" %in% names(pnps_tab_full)) {
  "obs_gt1.y"
} else {
  stop("Can't find obs_gt1 column in pnps_tab_full")
}

obs_n_col <- if ("obs_n" %in% names(pnps_tab_full)) {
  "obs_n"
} else if ("obs_n.x" %in% names(pnps_tab_full)) {
  "obs_n.x"
} else if ("obs_n.y" %in% names(pnps_tab_full)) {
  "obs_n.y"
} else {
  NA_character_
}

pnps_tab_full[, pnps_obs_gt1 := get(obs_gt1_col)]
pnps_tab_full[, pnps_null_mean := null_mean]
pnps_tab_full[, pnps_null_sd := null_sd]
pnps_tab_full[, pnps_p := p_perm]
pnps_tab_full[, pnps_q := q_BH]

if (!is.na(obs_n_col)) {
  pnps_tab_full[, pnps_obs_n := get(obs_n_col)]
} else {
  pnps_tab_full[, pnps_obs_n := NA_integer_]
}

pnps_tab_full[, pnps_FE := pnps_obs_gt1 / pmax(pnps_null_mean, 1e-9)]
pnps_tab_full[, pnps_log2FE := log2((pnps_obs_gt1 + 1) / (pnps_null_mean + 1))]
pnps_tab_full[, pnps_Z := (pnps_obs_gt1 - pnps_null_mean) / pmax(pnps_null_sd, 1e-9)]

pnps_keep_full <- pnps_tab_full[, .(
  PFAM_primary,
  pnps_obs_gt1, pnps_obs_n,
  pnps_null_mean, pnps_null_sd,
  pnps_p, pnps_q,
  pnps_FE, pnps_log2FE, pnps_Z
)]


combo_full <- merge(sweep_keep_full, pnps_keep_full, by="PFAM_primary", all=TRUE)

combo_full[, sweep_sig := !is.na(sweep_q) & sweep_q < 0.1]
combo_full[, pnps_sig  := !is.na(pnps_q)  & pnps_q  < 0.1]
combo_full[, any_sig := sweep_sig | pnps_sig]

plot_dt <- combo_full[!is.na(sweep_log2FE) & !is.na(pnps_log2FE)]

cat("sweep_keep_full rows:", nrow(sweep_keep_full), "\n")
cat("pnps_keep_full rows:", nrow(pnps_keep_full), "\n")
cat("overlap PFAMs:", length(intersect(sweep_keep_full$PFAM_primary, pnps_keep_full$PFAM_primary)), "\n")

# confirm effect-size cols exist
stopifnot(all(c("sweep_log2FE","pnps_log2FE") %in% names(combo_full)))







ggplot(plot_dt, aes(x = sweep_log2FE, y = pnps_log2FE)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(alpha=0.15, size=1) +
  geom_point(
    data = plot_dt[any_sig == TRUE],
    aes(size = pmin(sweep_obs_changes, 2000)),
    alpha = 0.8
  ) +
  geom_point(
    data = plot_dt[!is.na(sweep_n_cols) & sweep_n_cols >= 4],
    alpha = 0.9, size = 2.2
  ) +
  scale_size_continuous(name = "# 0D changes (cap)") +
  labs(
    x = "log2 fold-enrichment of 0D changes (Null 1; +1 pseudocount)",
    y = "log2 fold-enrichment of pN/pS>1 genes (null; +1 pseudocount)"
  ) +
  theme_classic()


gene_clean








#distribution shift




dt <- copy(gene_clean)

# Keep genes with PFAM + usable pN/pS
dt <- dt[!is.na(PFAM_primary)]
dt <- dt[!is.na(pN_pS_ratio) & is.finite(pN_pS_ratio)]

# Optional: drop zeros if you consider pN/pS=0 as "no nonsyn observed"
# dt <- dt[pN_pS_ratio > 0]

# log transform with small pseudocount (robust, avoids -Inf)
eps <- 1e-6
dt[, log_pnps := log(pN_pS_ratio + eps)]

# define permutation strata (this is the key choice)
dt[, stratum := interaction(MAG_id, Colony, Month, drop=TRUE)]

# Precompute per-stratum row indices for fast shuffling
idx_by_stratum <- split(seq_len(nrow(dt)), dt$stratum)

# Tail thresholds
thr_gt1  <- 1
thr_q90  <- as.numeric(quantile(dt$log_pnps, probs = 0.90, na.rm = TRUE))










obs <- dt[, .(
  n_genes      = .N,
  mean_log     = mean(log_pnps),
  median_log   = median(log_pnps),
  frac_gt1     = mean(pN_pS_ratio > thr_gt1),
  frac_q90     = mean(log_pnps > thr_q90)
), by = PFAM_primary]

# Optional: require enough genes per PFAM for stability
min_genes <- 30
obs <- obs[n_genes >= min_genes]



perm_once <- function() {
  pf_perm <- dt$PFAM_primary
  
  # shuffle PFAM labels within each stratum
  for (ii in idx_by_stratum) {
    if (length(ii) > 1L) pf_perm[ii] <- sample(pf_perm[ii], length(ii), replace = FALSE)
  }
  
  tmp <- data.table(
    PFAM_primary = pf_perm,
    log_pnps     = dt$log_pnps,
    gt1          = dt$pN_pS_ratio > thr_gt1,
    q90          = dt$log_pnps > thr_q90
  )
  
  tmp[, .(
    mean_log   = mean(log_pnps),
    median_log = median(log_pnps),
    frac_gt1   = mean(gt1),
    frac_q90   = mean(q90),
    n_genes    = .N
  ), by = PFAM_primary][n_genes >= min_genes]
}



B <- 500  # start with 200–500; scale up later if needed

null_list <- vector("list", B)
for (b in seq_len(B)) {
  null_list[[b]] <- perm_once()
}

null <- rbindlist(null_list, idcol = "perm")







# Helper: compute permutation p-values PFAM-by-PFAM for a single metric
perm_pvals <- function(metric) {
  # observed vector
  obs_v <- obs[, .(PFAM_primary, obs = get(metric))]
  
  # null distribution: PFAM x perm -> value
  null_v <- null[, .(PFAM_primary, perm, val = get(metric))]
  
  # join and compute p = (1 + #{null >= obs})/(B+1)
  m <- merge(null_v, obs_v, by = "PFAM_primary", allow.cartesian = TRUE)
  p <- m[, .(p = (1 + sum(val >= obs, na.rm = TRUE)) / (B + 1)), by = PFAM_primary]
  p[, q := p.adjust(p, method = "BH")]
  setnames(p, c("p","q"), paste0(c("p_","q_"), metric))
  p
}

p_mean   <- perm_pvals("mean_log")
p_median <- perm_pvals("median_log")
p_gt1    <- perm_pvals("frac_gt1")
p_q90    <- perm_pvals("frac_q90")

pnps_shift <- Reduce(function(x,y) merge(x,y, by="PFAM_primary", all=TRUE),
                     list(obs, p_mean, p_median, p_gt1, p_q90))






sweep <- copy(res_null)  # or res
sweep[, sweep_FE := obs_changes / pmax(null_mean, 1e-9)]
sweep[, sweep_logFE := log2((obs_changes + 1) / (null_mean + 1))]
setnames(sweep, c("p_null","q_BH"), c("sweep_p","sweep_q"))

combo <- merge(sweep, pnps_shift, by="PFAM_primary", all.x=TRUE)

# Example overlap calls
combo[, sweep_hit := sweep_q <= 0.1]
combo[, pnps_hit_gt1 := q_frac_gt1 <= 0.1]
combo[, pnps_hit_mean := q_mean_log <= 0.1]

combo[, overlap_gt1  := sweep_hit & pnps_hit_gt1]
combo[, overlap_mean := sweep_hit & pnps_hit_mean]



# make hits non-NA booleans
combo[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.2]
combo[, pnps_hit_gt1 := !is.na(q_frac_gt1) & q_frac_gt1 <= 0.2]

# overlap is now TRUE/FALSE only
combo[, overlap_gt1 := sweep_hit & pnps_hit_gt1]

# NOW this will work (no NA “phantom rows”)
combo[overlap_gt1 == TRUE][order(sweep_q, q_frac_gt1)][1:30,
                                                       .(PFAM_primary, obs_changes, null_mean, sweep_FE, sweep_q,
                                                         n_genes, frac_gt1, pNPS_q = q_frac_gt1)]





# how many PFAMs have pnps stats?
cat("pnps available:", sum(!is.na(combo$q_frac_gt1)), "of", nrow(combo), "\n")

# how many are sweep hits?
cat("sweep hits:", sum(combo$sweep_hit), "\n")

# overlap count
cat("overlap:", sum(combo$overlap_gt1), "\n")













g <- copy(gene_clean)
g <- g[!is.na(PFAM_primary) & PFAM_primary != ""]

# sweep results table: use res_null (Key: PFAM_primary; has q_BH, FE, Z)
sweep_res <- copy(res_null)
setkey(sweep_res, PFAM_primary)

# choose threshold for "sweep-enriched"
alpha <- 0.10
sweep_hits <- sweep_res[q_BH <= alpha & is.finite(obs_changes) & obs_changes > 0,
                        unique(PFAM_primary)]

g[, sweep_hit := PFAM_primary %chin% sweep_hits]
table(g$sweep_hit)






# Keep finite values
g_pi <- g[is.finite(pi_nonsyn) & is.finite(pi_syn) & pi_syn > 0 & pi_nonsyn >= 0]

# (A) πN distribution
p_piN <- ggplot(g_pi, aes(x = pi_nonsyn, fill = sweep_hit)) +
  geom_density(alpha = 0.35) +
  scale_x_continuous(trans = "log10") +
  theme_classic() +
  labs(x = expression(pi[N]), y = "Density", fill = "Sweep-enriched PFAM",
       title = expression("Gene-level " * pi[N] * " shift"))

# (B) πS distribution
p_piS <- ggplot(g_pi, aes(x = pi_syn, fill = sweep_hit)) +
  geom_density(alpha = 0.35) +
  scale_x_continuous(trans = "log10") +
  theme_classic() +
  labs(x = expression(pi[S]), y = "Density", fill = "Sweep-enriched PFAM",
       title = expression("Gene-level " * pi[S] * " shift"))

p_piN; p_piS


wilcox.test(log10(pi_nonsyn) ~ sweep_hit, data = as.data.frame(g_pi[pi_nonsyn > 0]))
wilcox.test(log10(pi_syn)    ~ sweep_hit, data = as.data.frame(g_pi[pi_syn > 0]))







pfam_sum <- g_pi[, .(
  n_genes = .N,
  med_piN = median(pi_nonsyn, na.rm = TRUE),
  med_piS = median(pi_syn, na.rm = TRUE),
  med_ratio = median(pi_nonsyn / pmax(pi_syn, 1e-12), na.rm = TRUE)
), by = .(PFAM_primary)]

pfam_sum[, sweep_hit := PFAM_primary %chin% sweep_hits]

# Plot PFAM medians
ggplot(pfam_sum[n_genes >= 10], aes(x = med_piS, y = med_piN, color = sweep_hit)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  labs(x = "PFAM median πS", y = "PFAM median πN", color = "Sweep-enriched PFAM",
       title = "PFAM-level median diversity")

# Tests at PFAM-level (less PFAM-size bias)
wilcox.test(log10(med_piN) ~ sweep_hit, data = as.data.frame(pfam_sum[med_piN > 0 & n_genes >= 10]))
wilcox.test(log10(med_piS) ~ sweep_hit, data = as.data.frame(pfam_sum[med_piS > 0 & n_genes >= 10]))











g2 <- g[!is.na(PFAM_primary) & PFAM_primary != ""]
g2[, is_gt1 := is.finite(pN_pS_ratio) & (pN_pS_ratio > 1)]

pfam_evt <- g2[, .(n_genes = .N, n_gt1 = sum(is_gt1, na.rm = TRUE)), by = PFAM_primary]
pfam_evt[, sweep_hit := PFAM_primary %chin% sweep_hits]

# observed tail curve: fraction PFAMs with >= k
K <- 0:10
obs_frac <- data.table(k = K, frac = sapply(K, function(k) mean(pfam_evt$n_gt1 >= k)))

# Null: shuffle is_gt1 labels across genes (preserves PFAM sizes and overall event rate)
set.seed(1)
B <- 500
null_mat <- replicate(B, {
  shuffled <- sample(g2$is_gt1)
  pfam_evt_null <- g2[, .(n_gt1 = sum(shuffled, na.rm = TRUE)), by = PFAM_primary]
  sapply(K, function(k) mean(pfam_evt_null$n_gt1 >= k))
})

null_dt <- data.table(
  k = rep(K, times = B),
  frac = as.vector(null_mat),
  b = rep(seq_len(B), each = length(K))
)

band <- null_dt[, .(
  lo = quantile(frac, 0.025),
  hi = quantile(frac, 0.975),
  mid = median(frac)
), by = k]

ggplot() +
  geom_ribbon(data = band, aes(x = k, ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(data = band, aes(x = k, y = mid), linewidth = 0.8) +
  geom_point(data = obs_frac, aes(x = k, y = frac), size = 2) +
  theme_classic() +
  labs(x = "k (genes with pN/pS>1 in PFAM)", y = "Fraction of PFAMs with ≥ k",
       title = "Observed vs null recurrence of pN/pS>1 within PFAMs")
band





dt <- copy(combo)





dt_plot <- dt[!is.na(sweep_logFE) & !is.na(frac_gt1)]

# significance flags
dt_plot[, sig_sweep := !is.na(sweep_q) & sweep_q <= 0.1]
dt_plot[, sig_pnps  := !is.na(q_frac_gt1) & q_frac_gt1 <= 0.1]

dt_plot[, sig_class := fifelse(sig_sweep & sig_pnps, "Both",
                               fifelse(sig_sweep, "Sweep only",
                                       fifelse(sig_pnps, "pN/pS only", "Neither")))]
dt_plot[, sig_class := factor(sig_class, levels=c("Neither","Sweep only","pN/pS only","Both"))]

# y axis: recurrence signal (transform fraction; avoids hard 0/1 problems)
eps <- 1e-4
dt_plot[, y_recur := qlogis(pmin(pmax(frac_gt1, eps), 1-eps))]  # logit(frac_gt1)

# label top 10 in "Both" among top-right (x>0 and high recurrence)
lab_dt <- dt_plot[sig_class=="Both" & sweep_logFE>0][order(-y_recur)][1:10]

ggplot(dt_plot, aes(x = sweep_logFE, y = y_recur)) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(aes(fill = sig_class, size = pmin(obs_changes, 1600)),
             shape=21, color="black", alpha=0.7, stroke=0.2) +
  ggrepel::geom_text_repel(data=lab_dt, aes(label=PFAM_primary),
                           size=3, max.overlaps=Inf) +
  theme_classic() +
  labs(x="log2 fold-enrichment of 0D changes (Null 1)",
       y="logit( fraction genes with pN/pS>1 )",
       fill="Significance class", size="# 0D changes (cap)")










# pnps_tab_full should be your full PFAM-level pnps>1 null results (like res_pnps_gt1)
pnps_tab_full <- copy(res_pnps_gt1)

# standardize column names
setnames(pnps_tab_full,
         old = c("obs_gt1.x","null_mean","null_sd","q_BH"),
         new = c("pnps_obs_gt1","pnps_null_mean","pnps_null_sd","pnps_q"),
         skip_absent=TRUE)

# compute pnps FE and LFC (same spirit as sweep_logFE)
pnps_tab_full[, pnps_FE  := pnps_obs_gt1 / pmax(pnps_null_mean, 1e-9)]
pnps_tab_full[, pnps_LFC := log2(pnps_FE + 1e-9)]  # tiny floor

# merge into dt (dt already has sweep_logFE etc)
dt3 <- merge(dt, pnps_tab_full[, .(PFAM_primary, pnps_q, pnps_LFC, pnps_FE)],
             by="PFAM_primary", all.x=TRUE)






dt3[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_logFE > 0]
dt3[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_LFC   > 0]

dt3[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                           fifelse(sweep_hit, "Sweep only",
                                   fifelse(pnps_hit, "pN/pS only","Neither")))]
dt3[, sig_class := factor(sig_class, levels=c("Neither","Sweep only","pN/pS only","Both"))]









plot_dt <- dt3[!is.na(sweep_logFE) & !is.na(pnps_LFC)]

lab <- plot_dt[sig_class=="Both" & sweep_logFE>0 & pnps_LFC>0][
  order(-(sweep_logFE + pnps_LFC))
][1:10]

ggplot(plot_dt, aes(x=sweep_logFE, y=pnps_LFC)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(aes(fill=sig_class, size=pmin(obs_changes,1600)),
             shape=21, color="black", alpha=0.7, stroke=0.2) +
  ggrepel::geom_text_repel(data=lab, aes(label=PFAM_primary),
                           size=3, max.overlaps=Inf) +
  theme_classic() +
  labs(x="log2 fold-enrichment of 0D changes (Null 1)",
       y="log2 fold-enrichment of genes with pN/pS>1 (null-based)",
       fill="Significance class", size="# 0D changes (cap)")






tab <- table(Sweep = plot_dt$sweep_hit, pNPS = plot_dt$pnps_hit)

# guard against degenerate tables
if (nrow(tab) < 2 || ncol(tab) < 2) {
  print(tab)
  stop("Degenerate 2x2: pnps_hit or sweep_hit has no TRUEs. Check thresholds or inputs.")
}

ft <- fisher.test(tab)
ft$estimate
ft$p.value
tab
plot_dt[, .N, by=sig_class][order(sig_class)]
plot_dt[, .(sweep_true=sum(sweep_hit, na.rm=TRUE),
            pnps_true=sum(pnps_hit, na.rm=TRUE),
            both_true=sum(sweep_hit & pnps_hit, na.rm=TRUE))]
summary(plot_dt$pnps_LFC)
summary(plot_dt$pnps_q)




# pnps_LFC = log2 fold-enrichment of fraction(pN/pS>1) vs null
# edit `pnps_null_mean_col` to whatever your null mean column is called
pnps_null_mean_col <- NULL
for (cand in c("pnps_null_mean","null_mean_frac_gt1","frac_gt1_null_mean","null_mean")) {
  if (cand %in% names(dt)) { pnps_null_mean_col <- cand; break }
}

if (!"pnps_LFC" %in% names(dt)) {
  if (!is.null(pnps_null_mean_col)) {
    dt[, pnps_LFC := log2((frac_gt1 + 1e-6) / (get(pnps_null_mean_col) + 1e-6))]
  } else {
    # fallback: plot observed frac_gt1 directly (NOT enrichment)
    dt[, pnps_LFC := log2(frac_gt1 + 1e-6)]
    message("No pnps null mean found; pnps_LFC is log2(observed frac_gt1), not enrichment vs null.")
  }
}
plot_dt <- dt[
  !is.na(sweep_LFC) & !is.na(pnps_LFC),
  .(PFAM_primary,
    obs_changes, null_mean, null_sd,
    sweep_p, sweep_q, sweep_FE, sweep_LFC,
    frac_gt1, q_frac_gt1, pnps_LFC,
    sig_class)
]




dt[, .(
  has_sweep_LFC = "sweep_LFC" %in% names(dt),
  has_sweep_logFE = "sweep_logFE" %in% names(dt),
  has_pnps_LFC = "pnps_LFC" %in% names(dt),
  has_frac_gt1 = "frac_gt1" %in% names(dt)
)]
summary(dt$sweep_LFC)
summary(dt$pnps_LFC)

y_floor <- -5   # pick -4, -5, or -6 depending on how aggressive you want

ggplot(plot_dt, aes(x = sweep_LFC, y = pnps_LFC)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = sig_class, size = pmin(obs_changes, 1600)),
             shape = 21, color = "black", alpha = 0.7, stroke = 0.2) +
  coord_cartesian(ylim = c(y_floor, NA)) +   # <<< zoom without dropping data
  theme_classic() +
  labs(x = "log2 fold-enrichment of 0D changes (Null 1)",
       y = "log2 fold-enrichment of pN/pS>1 genes (Null 1; +1 pseudocount)",
       fill = "Significance class", size = "# 0D changes (cap)")



plot_dt

y_floor <- -5
plot_dt_zoom <- plot_dt[pnps_LFC >= y_floor | is.na(pnps_LFC)]

ggplot(plot_dt_zoom, aes(x = sweep_LFC, y = pnps_LFC)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = sig_class, size = pmin(obs_changes, 1600)),
             shape = 21, color = "black", alpha = 0.7, stroke = 0.2) +
  theme_classic()




y_floor <- -5

ggplot(plot_dt, aes(x = sweep_LFC, y = pnps_LFC)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = sig_class, size = pmin(obs_changes, 1600)),
             shape = 21, color = "black", alpha = 0.7, stroke = 0.2) +
  coord_cartesian(ylim = c(y_floor, NA)) +
  annotate("text", x = min(plot_dt$sweep_LFC, na.rm=TRUE), y = y_floor,
           label = paste0("Truncated below ", y_floor), hjust = 0, vjust = -0.5, size = 3) +
  theme_classic()








if (!"sweep_LFC" %in% names(dt)) {
  if ("sweep_logFE" %in% names(dt)) {
    dt[, sweep_LFC := sweep_logFE]
  } else if ("sweep_FE" %in% names(dt)) {
    dt[, sweep_LFC := log2(pmax(sweep_FE, 1e-12))]
  } else if (all(c("obs_changes","null_mean") %in% names(dt))) {
    dt[, sweep_LFC := log2((obs_changes + 1) / (null_mean + 1))]  # +1 pseudocount version
  } else {
    stop("Can't construct sweep_LFC: need sweep_logFE, sweep_FE, or (obs_changes & null_mean)")
  }
}
names(dt)






















dt <- copy(combo)

# sweep_LFC: log2 fold-enrichment of 0D changes vs null (with +1 pseudocount)
if (!"sweep_LFC" %in% names(dt)) {
  if ("sweep_logFE" %in% names(dt)) {
    dt[, sweep_LFC := sweep_logFE]
  } else if (all(c("obs_changes","null_mean") %in% names(dt))) {
    dt[, sweep_LFC := log2((obs_changes + 1) / (null_mean + 1))]
  } else if (all(c("sweep_obs_changes","sweep_null_mean") %in% names(dt))) {
    dt[, sweep_LFC := log2((sweep_obs_changes + 1) / (sweep_null_mean + 1))]
  } else {
    stop("Can't construct sweep_LFC: need sweep_logFE or (obs_changes,null_mean) or (sweep_obs_changes,sweep_null_mean).")
  }
}










pnps_tab_full <- copy(res_pnps_gt1)

# find obs_gt1 column robustly
obs_gt1_col <- if ("obs_gt1" %in% names(pnps_tab_full)) "obs_gt1"
else if ("obs_gt1.x" %in% names(pnps_tab_full)) "obs_gt1.x"
else stop("Can't find obs_gt1 column in res_pnps_gt1")

# standardize key cols
pnps_tab_full[, pnps_obs_gt1  := get(obs_gt1_col)]
pnps_tab_full[, pnps_null_mean := null_mean]
pnps_tab_full[, pnps_q := q_BH]

# null-based enrichment
pnps_tab_full[, pnps_FE  := pnps_obs_gt1 / pmax(pnps_null_mean, 1e-12)]
pnps_tab_full[, pnps_LFC := log2(pnps_FE + 1e-12)]

# merge into dt
dt <- merge(
  dt,
  pnps_tab_full[, .(PFAM_primary, pnps_obs_gt1, pnps_null_mean, pnps_LFC, pnps_q)],
  by = "PFAM_primary",
  all.x = TRUE
)





dt[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_LFC > 0]
dt[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_LFC  > 0]

dt[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                          fifelse(sweep_hit, "Sweep only",
                                  fifelse(pnps_hit, "pN/pS only", "Neither")))]
dt[, sig_class := factor(sig_class, levels=c("Neither","Sweep only","pN/pS only","Both"))]




plot_dt <- dt[!is.na(sweep_LFC) & !is.na(pnps_LFC)]

y_floor <- -5

lab <- plot_dt[sig_class=="Both" & sweep_LFC>0 & pnps_LFC>0][
  order(-(sweep_LFC + pnps_LFC))
][1:10]

ggplot(plot_dt, aes(x=sweep_LFC, y=pnps_LFC)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(aes(fill=sig_class, size=pmin(obs_changes,1600)),
             shape=21, color="black", alpha=0.7, stroke=0.2) +
  ggrepel::geom_text_repel(data=lab, aes(label=PFAM_primary),
                           size=3, max.overlaps=Inf) +
  coord_cartesian(ylim=c(y_floor, NA)) +
  theme_classic() +
  labs(x="log2 fold-enrichment of 0D changes (Null 1; +1 pseudocount)",
       y="log2 fold-enrichment of pN/pS>1 genes (null-based)",
       fill="Significance class", size="# 0D changes (cap)")


g2 <- g[!is.na(PFAM_primary) & PFAM_primary != ""]
g2[, is_gt1 := is.finite(pN_pS_ratio) & (pN_pS_ratio > 1)]

K <- 0:10
obs_frac <- data.table(k=K, frac=sapply(K, function(k) mean(g2[, sum(is_gt1), by=PFAM_primary]$V1 >= k)))

set.seed(1)
B <- 500

null_mat <- replicate(B, {
  g2[, is_gt1_shuf := sample(is_gt1)]
  tmp <- g2[, .(n_gt1 = sum(is_gt1_shuf, na.rm=TRUE)), by=PFAM_primary]
  sapply(K, function(k) mean(tmp$n_gt1 >= k))
})

null_dt <- data.table(k=rep(K, times=B), frac=as.vector(null_mat))
band <- null_dt[, .(
  lo  = quantile(frac, 0.025),
  hi  = quantile(frac, 0.975),
  mid = median(frac)
), by=k]

ggplot() +
  geom_ribbon(data=band, aes(x=k, ymin=lo, ymax=hi), alpha=0.2) +
  geom_line(data=band, aes(x=k, y=mid), linewidth=0.8) +
  geom_point(data=obs_frac, aes(x=k, y=frac), size=2) +
  theme_classic() +
  labs(x="k (genes with pN/pS>1 in PFAM)",
       y="Fraction of PFAMs with ≥ k",
       title="Observed vs null recurrence of pN/pS>1 within PFAMs")













# g2 = gene-level table
g2 <- as.data.table(gene_clean)
g2 <- g2[!is.na(PFAM_primary) & PFAM_primary != ""]

# Define gene-level "high pN/pS" event
g2[, is_gt1 := is.finite(pN_pS_ratio) & (pN_pS_ratio > 1)]

# PFAM-level observed recurrence
pfam_obs <- g2[, .(
  n_genes = .N,
  n_gt1   = sum(is_gt1, na.rm = TRUE),
  frac_gt1 = sum(is_gt1, na.rm = TRUE) / .N
), by = PFAM_primary]






Kmax <- 10
K <- 0:Kmax

# Observed curve: fraction of PFAMs with >=k events
obs_curve <- pfam_obs[, .(k = K, frac = sapply(K, function(k) mean(n_gt1 >= k)))]

# Null bootstraps
set.seed(1)
B <- 500

# pre-allocate: B x length(K)
null_mat <- matrix(NA_real_, nrow = B, ncol = length(K))

for (b in seq_len(B)) {
  shuffled <- sample(g2$is_gt1)  # shuffle event labels
  pfam_null <- g2[, .(n_gt1 = sum(shuffled, na.rm = TRUE)), by = PFAM_primary]
  null_mat[b, ] <- sapply(K, function(k) mean(pfam_null$n_gt1 >= k))
}

band <- data.table(
  k = K,
  lo = apply(null_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
  hi = apply(null_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
  mid = apply(null_mat, 2, median, na.rm = TRUE)
)

p_band <- ggplot() +
  geom_ribbon(data = band, aes(x = k, ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(data = band, aes(x = k, y = mid), linewidth = 0.8) +
  geom_point(data = obs_curve, aes(x = k, y = frac), size = 2) +
  theme_classic() +
  labs(
    x = "k (genes with pN/pS > 1 in PFAM)",
    y = "Fraction of PFAMs with ≥ k",
    title = "Observed vs null recurrence of pN/pS>1 within PFAMs"
  )

p_band






func_props <- plot_dt2[pfam_str %chin% top_cats,
                       .(N=.N),
                       by=.(pfam_str, sig_class)]
func_props[, total := sum(N), by=pfam_str]
func_props[, prop := N / total]

p_inset_props <- ggplot(func_props, aes(x=reorder(pfam_str, prop, sum), y=prop, fill=sig_class)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  theme_classic() +
  labs(x=NULL, y="Within-category fraction", fill=NULL,
       title="Composition of significance classes within functional categories")

p_inset_props






tab <- table(Sweep=plot_dt2$sweep_hit, pNPS=plot_dt2$pnps_hit)
ft <- fisher.test(tab)
tab
ft$estimate
ft$p.value

p_quadrant <- ggplot(plot_dt2[, .N, by=sig_class], aes(x=sig_class, y=N, fill=sig_class)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, hjust=1), legend.position="none") +
  labs(x=NULL, y="PFAMs", title="Counts per class")+
  scale_y_log10()

p_quadrant







ggplot(plot_dt2[!is.na(sweep_LFC)], aes(x=sweep_LFC, fill=pnps_hit)) +
  geom_density(alpha=0.4) +
  theme_classic() +
  labs(x="sweep_LFC", fill="pnps_hit")

ggplot(plot_dt2[!is.na(pnps_LFC)], aes(x=pnps_LFC, fill=sweep_hit)) +
  geom_density(alpha=0.4) +
  theme_classic() +
  labs(x="pnps_LFC", fill="sweep_hit")
































### -------------------------
### 0) Helpers: robust column pickers
pick_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) stop("Could not find any of: ", paste(candidates, collapse=", "))
  hit
}

### -------------------------
### 1) Sweep table (full)
sweep_tab_full <- copy(res_null)

sweep_tab_full[, sweep_obs_changes := obs_changes]
sweep_tab_full[, sweep_null_mean  := null_mean]
sweep_tab_full[, sweep_null_sd    := null_sd]
sweep_tab_full[, sweep_p          := p_null]
sweep_tab_full[, sweep_q          := q_BH]

# +1 pseudocount log2 fold enrichment
sweep_tab_full[, sweep_log2FE := log2((sweep_obs_changes + 1) / (sweep_null_mean + 1))]

# Optional: effect size without pseudo (can be unstable when null_mean~0)
sweep_tab_full[, sweep_FE := sweep_obs_changes / pmax(sweep_null_mean, 1e-9)]

# If you have these columns elsewhere (obs_pfam), merge them in; otherwise skip
# e.g. obs_pfam has n_units, n_cols, n_MAGs, n_genus
if (exists("obs_pfam")) {
  sweep_tab_full <- merge(
    sweep_tab_full,
    obs_pfam[, .(PFAM_primary, n_units, n_cols, n_MAGs, n_genus)],
    by="PFAM_primary", all.x=TRUE
  )
  setnames(sweep_tab_full,
           c("n_units","n_cols","n_MAGs","n_genus"),
           c("sweep_n_units","sweep_n_cols","sweep_n_MAGs","sweep_n_genus"))
}

sweep_keep_full <- sweep_tab_full[, .(
  PFAM_primary,
  obs_changes = sweep_obs_changes,
  null_mean   = sweep_null_mean,
  null_sd     = sweep_null_sd,
  sweep_p, sweep_q,
  sweep_FE, sweep_log2FE
)]

### -------------------------
### 2) pN/pS>1 table (full)
pnps_tab_full <- copy(res_pnps_gt1)

obs_gt1_col <- pick_col(pnps_tab_full, c("obs_gt1","obs_gt1.x","obs_gt1.y","obs_gt1.x"))
obs_n_col   <- (c("obs_n","obs_n.x","obs_n.y")[c("obs_n","obs_n.x","obs_n.y") %in% names(pnps_tab_full)][1])

pnps_tab_full[, pnps_obs_gt1 := get(obs_gt1_col)]
pnps_tab_full[, pnps_null_mean := null_mean]
pnps_tab_full[, pnps_null_sd   := null_sd]
pnps_tab_full[, pnps_p := if ("p_perm" %in% names(pnps_tab_full)) p_perm else if ("p_null" %in% names(pnps_tab_full)) p_null else NA_real_]
pnps_tab_full[, pnps_q := q_BH]

if (!is.na(obs_n_col)) pnps_tab_full[, pnps_obs_n := get(obs_n_col)] else pnps_tab_full[, pnps_obs_n := NA_integer_]

pnps_tab_full[, pnps_FE     := pnps_obs_gt1 / pmax(pnps_null_mean, 1e-9)]
pnps_tab_full[, pnps_log2FE := log2((pnps_obs_gt1 + 1) / (pnps_null_mean + 1))]
pnps_tab_full[, pnps_Z      := (pnps_obs_gt1 - pnps_null_mean) / pmax(pnps_null_sd, 1e-9)]

pnps_keep_full <- pnps_tab_full[, .(
  PFAM_primary,
  pnps_obs_gt1, pnps_obs_n,
  pnps_null_mean, pnps_null_sd,
  pnps_p, pnps_q,
  pnps_FE, pnps_log2FE, pnps_Z
)]

### -------------------------
### 3) Merge sweep + pnps
combo_full <- merge(sweep_keep_full, pnps_keep_full, by="PFAM_primary", all=TRUE)

# define hits (your requested definition)
combo_full[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_log2FE > 0]
combo_full[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_log2FE > 0]

combo_full[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                                  fifelse(sweep_hit, "Sweep only",
                                          fifelse(pnps_hit, "pN/pS only", "Neither")))]
combo_full[, sig_class := factor(sig_class, levels=c("Neither","Sweep only","pN/pS only","Both"))]



### -------------------------
### 4) Join PFAM functional label (pfam_str) for inset
pfam_func <- unique(gene_clean[!is.na(PFAM_primary) & PFAM_primary!="", .(PFAM_primary, pfam_str)])
combo_full <- merge(combo_full, pfam_func, by="PFAM_primary", all.x=TRUE)
combo_full[is.na(pfam_str) | pfam_str=="", pfam_str := "Unknown"]

### -------------------------
### 5) MAIN PLOT (Option A): sweep_log2FE vs pnps_log2FE with y-floor zoom
plot_dt <- combo_full[!is.na(sweep_log2FE) & !is.na(pnps_log2FE)]

# label top-right BOTH candidates
lab <- plot_dt[sig_class=="Both" & sweep_log2FE>0 & pnps_log2FE>0][
  order(-(sweep_log2FE + pnps_log2FE))
][1:10]

y_floor <- -5 


p_main <- ggplot(plot_dt, aes(x=sweep_log2FE, y=pnps_log2FE)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(aes(fill=sig_class, size=pmin(obs_changes,1600)),
             shape=21, color="black", alpha=0.7, stroke=0.2) +
  coord_cartesian(ylim=c(y_floor, NA)) +
  ggrepel::geom_text_repel(data=lab, aes(label=PFAM_primary),
                           size=3, max.overlaps=Inf) +
  theme_classic() +
  labs(
    x="log2 fold-enrichment of 0D changes (Null 1; +1 pseudocount)",
    y="log2 fold-enrichment of genes with pN/pS>1 (null; +1 pseudocount)",
    fill="Significance class", size="# 0D changes (cap)"
  )

p_main




### -------------------------
### 6) INSET: histogram-like bar counts of pfam_str by sig_class (top functions)
# count PFAMs per function × class
func_counts <- combo_full[, .N, by=.(pfam_str, sig_class)]

# choose top functions by total PFAMs
top_funcs <- combo_full[, .N, by=pfam_str][order(-N)][1:15, pfam_str]
func_counts <- func_counts[pfam_str %chin% top_funcs]

p_func <- ggplot(func_counts, aes(x=reorder(pfam_str, N, sum), y=N, fill=sig_class)) +
  geom_col(position="stack") +
  coord_flip() +
  theme_classic() +
  labs(x=NULL, y="# PFAMs (top functions)", fill="Significance class")

p_func












Yep — we can keep the “floor” and make it interpretable by switching the y-axis to a pseudocount log2 fold-enrichment that never goes to 
−∞
−∞, and by explicitly treating the “no events / no expectation” PFAMs as a separate regime.

Below is (i) what’s going on with the floor, (ii) what we’re estimating in each figure, (iii) clean bulk code to reproduce both figures (main scatter + a functional-class inset), and (iv) an optional inset that directly answers: among sweep-PFAMs, how many pN/pS>1 genes would we expect under a null?

1) Why the “floor” happens (and how the pseudocount fixes it)

If you use an enrichment like

log⁡2(obs_gt1null_mean)
log
2
	​

(
null_mean
obs_gt1
	​

)

then PFAMs with obs_gt1 = 0 go to 
−∞
−∞, and PFAMs with null_mean ≈ 0 can explode to 
+∞
+∞. That’s why your plot produces a huge negative shelf (“floor”).

Use the +1 pseudocount version instead:

pnps_log2FE=log⁡2(obs_gt1+1null_mean+1)
pnps_log2FE=log
2
	​

(
null_mean+1
obs_gt1+1
	​

)

This does two good things:

PFAMs with 0 events become 
log⁡2(1/(null+1))
log
2
	​

(1/(null+1)) — finite.

PFAMs with tiny expected counts don’t blow up infinitely.

So: keep the floor but make it meaningful. Values near the floor now mean “basically zero pN/pS>1 genes relative to expectation”, not “math artifact”.

2) What exactly each axis is estimating
X axis (sweep signal)

Observed enrichment of 0D changes per PFAM vs Null 1

sweep_log2FE=log⁡2(obs_changes+1null_mean+1)
sweep_log2FE=log
2
	​

(
null_mean+1
obs_changes+1
	​

)

Interpretation: PFAMs to the right have more 0D sweep changes than expected by your metagenome-randomization null.

Y axis (pN/pS>1 recurrence)

Observed enrichment of “genes with pN/pS>1” within PFAM vs a PFAM-aware null

pnps_log2FE=log⁡2(obs_gt1+1null_mean+1)
pnps_log2FE=log
2
	​

(
null_mean+1
obs_gt1+1
	​

)

Interpretation: PFAMs above 0 have more pN/pS>1 genes than expected given your null construction for pN/pS events.

Quadrants

Top-right: enriched for sweeps and enriched for pN/pS>1 recurrence → best “adaptive candidates”

Bottom-right: sweep-enriched but not pN/pS>1 enriched → still interesting (selection may not manifest as >1 events)

Top-left: pN/pS>1 enriched but no sweep enrichment → possibly selection without sweep-like 0D pattern (or different mechanism)

Bottom-left: neither

3) Bulk code: build the merged table + main scatter (Option A) + functional inset

This assumes you have:

res_null (full sweep PFAM null results: PFAM_primary, obs_changes, null_mean, null_sd, p_null, q_BH)

res_pnps_gt1 (PFAM pN/pS>1 null results: has PFAM_primary, an obs count column like obs_gt1 or obs_gt1.x, and null_mean, null_sd, p_perm, q_BH)

gene_clean containing PFAM_primary, pfam_str (or you can join from elsewhere)

library(data.table)
library(ggplot2)
library(ggrepel)

### -------------------------
### 0) Helpers: robust column pickers
pick_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) stop("Could not find any of: ", paste(candidates, collapse=", "))
  hit
}

### -------------------------
### 1) Sweep table (full)
sweep_tab_full <- copy(res_null)

sweep_tab_full[, sweep_obs_changes := obs_changes]
sweep_tab_full[, sweep_null_mean  := null_mean]
sweep_tab_full[, sweep_null_sd    := null_sd]
sweep_tab_full[, sweep_p          := p_null]
sweep_tab_full[, sweep_q          := q_BH]

# +1 pseudocount log2 fold enrichment
sweep_tab_full[, sweep_log2FE := log2((sweep_obs_changes + 1) / (sweep_null_mean + 1))]

# Optional: effect size without pseudo (can be unstable when null_mean~0)
sweep_tab_full[, sweep_FE := sweep_obs_changes / pmax(sweep_null_mean, 1e-9)]

# If you have these columns elsewhere (obs_pfam), merge them in; otherwise skip
# e.g. obs_pfam has n_units, n_cols, n_MAGs, n_genus
if (exists("obs_pfam")) {
  sweep_tab_full <- merge(
    sweep_tab_full,
    obs_pfam[, .(PFAM_primary, n_units, n_cols, n_MAGs, n_genus)],
    by="PFAM_primary", all.x=TRUE
  )
  setnames(sweep_tab_full,
           c("n_units","n_cols","n_MAGs","n_genus"),
           c("sweep_n_units","sweep_n_cols","sweep_n_MAGs","sweep_n_genus"))
}

sweep_keep_full <- sweep_tab_full[, .(
  PFAM_primary,
  obs_changes = sweep_obs_changes,
  null_mean   = sweep_null_mean,
  null_sd     = sweep_null_sd,
  sweep_p, sweep_q,
  sweep_FE, sweep_log2FE
)]

### -------------------------
### 2) pN/pS>1 table (full)
pnps_tab_full <- copy(res_pnps_gt1)

obs_gt1_col <- pick_col(pnps_tab_full, c("obs_gt1","obs_gt1.x","obs_gt1.y","obs_gt1.x"))
obs_n_col   <- (c("obs_n","obs_n.x","obs_n.y")[c("obs_n","obs_n.x","obs_n.y") %in% names(pnps_tab_full)][1])

pnps_tab_full[, pnps_obs_gt1 := get(obs_gt1_col)]
pnps_tab_full[, pnps_null_mean := null_mean]
pnps_tab_full[, pnps_null_sd   := null_sd]
pnps_tab_full[, pnps_p := if ("p_perm" %in% names(pnps_tab_full)) p_perm else if ("p_null" %in% names(pnps_tab_full)) p_null else NA_real_]
pnps_tab_full[, pnps_q := q_BH]

if (!is.na(obs_n_col)) pnps_tab_full[, pnps_obs_n := get(obs_n_col)] else pnps_tab_full[, pnps_obs_n := NA_integer_]

pnps_tab_full[, pnps_FE     := pnps_obs_gt1 / pmax(pnps_null_mean, 1e-9)]
pnps_tab_full[, pnps_log2FE := log2((pnps_obs_gt1 + 1) / (pnps_null_mean + 1))]
pnps_tab_full[, pnps_Z      := (pnps_obs_gt1 - pnps_null_mean) / pmax(pnps_null_sd, 1e-9)]

pnps_keep_full <- pnps_tab_full[, .(
  PFAM_primary,
  pnps_obs_gt1, pnps_obs_n,
  pnps_null_mean, pnps_null_sd,
  pnps_p, pnps_q,
  pnps_FE, pnps_log2FE, pnps_Z
)]

### -------------------------
### 3) Merge sweep + pnps
combo_full <- merge(sweep_keep_full, pnps_keep_full, by="PFAM_primary", all=TRUE)

# define hits (your requested definition)
combo_full[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_log2FE > 0]
combo_full[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_log2FE > 0]

combo_full[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                           fifelse(sweep_hit, "Sweep only",
                             fifelse(pnps_hit, "pN/pS only", "Neither")))]
combo_full[, sig_class := factor(sig_class, levels=c("Neither","Sweep only","pN/pS only","Both"))]

### -------------------------
### 4) Join PFAM functional label (pfam_str) for inset
pfam_func <- unique(gene_clean[!is.na(PFAM_primary) & PFAM_primary!="", .(PFAM_primary, pfam_str)])
combo_full <- merge(combo_full, pfam_func, by="PFAM_primary", all.x=TRUE)
combo_full[is.na(pfam_str) | pfam_str=="", pfam_str := "Unknown"]

### -------------------------
### 5) MAIN PLOT (Option A): sweep_log2FE vs pnps_log2FE with y-floor zoom
plot_dt <- combo_full[!is.na(sweep_log2FE) & !is.na(pnps_log2FE)]

# label top-right BOTH candidates
lab <- plot_dt[
  sig_class=="Both" &
    sweep_log2FE > 0 &
    pnps_log2FE  > 0 &
    !is.na(PFAM_primary) & PFAM_primary != ""
][
  order(-(sweep_log2FE + pnps_log2FE))
][
  , .SD[1], by=PFAM_primary   # <-- keep best row per PFAM
][
  order(-(sweep_log2FE + pnps_log2FE))
][1:15]

y_floor <- -5  # tweak this (e.g. -4, -5, -6)

p_main <- ggplot(plot_dt, aes(x=sweep_log2FE, y=pnps_log2FE)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(aes(fill=sig_class, size=pmin(obs_changes,1600)),
             shape=21, color="black", alpha=0.7, stroke=0.2) +
  coord_cartesian(ylim=c(y_floor, NA)) +
  ggrepel::geom_text_repel(
    data = lab,
    aes(label = PFAM_primary),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_classic() +
  labs(
    x="log2 fold-enrichment of 0D changes (Null 1; +1 pseudocount)",
    y="log2 fold-enrichment of genes with pN/pS>1 (null; +1 pseudocount)",
    fill="Significance class", size="# 0D changes (cap)"
  )

p_main

### -------------------------
### 6) INSET: histogram-like bar counts of pfam_str by sig_class (top functions)
# count PFAMs per function × class
func_counts <- combo_full[, .N, by=.(pfam_str, sig_class)]

# choose top functions by total PFAMs
top_funcs <- combo_full[, .N, by=pfam_str][order(-N)][1:15, pfam_str]
func_counts <- func_counts[pfam_str %chin% top_funcs]

p_func <- ggplot(func_counts, aes(x=reorder(pfam_str, N, sum), y=N, fill=sig_class)) +
  geom_col(position="stack") +
  coord_flip() +
  theme_classic() +
  labs(x=NULL, y="# PFAMs (top functions)", fill="Significance class")

p_func














g2 <- gene_clean[!is.na(PFAM_primary) & PFAM_primary != ""]
g2[, is_gt1 := is.finite(pN_pS_ratio) & (pN_pS_ratio > 1)]

pfam_evt <- g2[, .(n_genes=.N, n_gt1=sum(is_gt1, na.rm=TRUE)), by=PFAM_primary]

# choose sweep PFAM set from combo_full (or from sweep_q threshold)
sweep_hits <- combo_full[sweep_hit==TRUE, PFAM_primary]
pfam_evt[, sweep_hit := PFAM_primary %chin% sweep_hits]

K <- 0:10
obs_frac_sweep <- data.table(
  k = K,
  frac = sapply(K, function(k) mean(pfam_evt[sweep_hit==TRUE]$n_gt1 >= k))
)

set.seed(1)
B <- 500
null_mat <- replicate(B, {
  shuffled <- sample(g2$is_gt1)
  tmp <- g2[, .(n_gt1=sum(shuffled, na.rm=TRUE)), by=PFAM_primary]
  tmp[, sweep_hit := PFAM_primary %chin% sweep_hits]
  sapply(K, function(k) mean(tmp[sweep_hit==TRUE]$n_gt1 >= k))
})

null_dt <- data.table(k=rep(K, B), frac=as.vector(null_mat))
band <- null_dt[, .(
  lo=quantile(frac, 0.025),
  hi=quantile(frac, 0.975),
  mid=median(frac)
), by=k]

ggplot() +
  geom_ribbon(data=band, aes(x=k, ymin=lo, ymax=hi), alpha=0.2) +
  geom_line(data=band, aes(x=k, y=mid), linewidth=0.8) +
  geom_point(data=obs_frac_sweep, aes(x=k, y=frac), size=2) +
  theme_classic() +
  labs(x="k (# genes with pN/pS>1 in PFAM)",
       y="Fraction of sweep-PFAMs with ≥k",
       title="Sweep PFAMs: observed vs null recurrence of pN/pS>1")










combo_full
combo_full[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_log2FE > 0]
combo_full[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_log2FE  > 0]

combo_full[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                                  fifelse(sweep_hit, "Sweep only",
                                          fifelse(pnps_hit, "pN/pS only", "Neither")))]
combo_full[, sig_class := factor(sig_class,
                                 levels = c("Neither","Sweep only","pN/pS only","Both"))]


table(combo_full$sig_class, useNA = "ifany")


both_pfams <- combo_full[sig_class == "Both", unique(PFAM_primary)]
sweep_only_pfams <- combo_full[sig_class == "Sweep only", unique(PFAM_primary)]


g <- copy(gene_clean)
g <- g[!is.na(PFAM_primary) & PFAM_primary != ""]
g <- g[is.finite(pN_pS_ratio)]

eps <- 1e-6
g[, log_pnps := log(pN_pS_ratio + eps)]

g[, gene_class := fifelse(PFAM_primary %chin% both_pfams, "Both",
                          fifelse(PFAM_primary %chin% sweep_only_pfams, "Sweep only",
                                  "Background"))]

g[, gene_class := factor(gene_class, levels = c("Background","Sweep only","Both"))]

table(g$gene_class)



sum_dt <- g[, .(
  n = .N,
  mean_pnps = mean(pN_pS_ratio, na.rm = TRUE),
  median_pnps = median(pN_pS_ratio, na.rm = TRUE),
  mean_log_pnps = mean(log_pnps, na.rm = TRUE),
  median_log_pnps = median(log_pnps, na.rm = TRUE),
  frac_gt1 = mean(pN_pS_ratio > 1, na.rm = TRUE)
), by = gene_class]

sum_dt



combo_full[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_log2FE > 0]
combo_full[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_log2FE  > 0]

combo_full[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                                  fifelse(sweep_hit, "Sweep only",
                                          fifelse(pnps_hit, "pN/pS only", "Neither")))]
combo_full[, sig_class := factor(sig_class,
                                 levels = c("Neither","Sweep only","pN/pS only","Both"))]



tail_pfam <- function(x) {
  x <- as.integer(x)
  x <- x[is.finite(x)]
  kmax <- max(x)
  data.table(
    k = 1:kmax,
    S = vapply(1:kmax, function(kk) mean(x >= kk), numeric(1))
  )
}

tail_dt <- tail_pfam(res[obs_changes.x > 0]$obs_changes.x)




both_pfams_dt <- combo_full[sig_class == "Both",
                            .(PFAM_primary, sig_class)]

both_tail <- merge(
  both_pfams_dt,
  res[, .(PFAM_primary, obs_changes.x, n_cols, null_mean, null_sd)],
  by = "PFAM_primary",
  all.x = TRUE
)

both_tail <- both_tail[!is.na(obs_changes.x) & obs_changes.x > 0]
both_tail[, PFAM_label := PFAM_primary]


setkey(tail_dt, k)
both_tail[, k_plot := pmin(obs_changes.x, max(tail_dt$k))]
both_tail <- tail_dt[both_tail, on = .(k = k_plot)]





p_tail_both <- ggplot(tail_dt, aes(x = k, y = S)) +
  geom_step(linewidth = 1, color = "black") +
  geom_point(
    data = both_tail,
    aes(x = obs_changes.x, y = S),
    inherit.aes = FALSE,
    size = 2.8,
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = both_tail,
    aes(x = obs_changes.x, y = S, label = PFAM_label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "k (0D changes per PFAM)",
    y = "Fraction of PFAMs with ≥ k changes",
    title = "Tail of 0D changes per PFAM",
    subtitle = "Red points = PFAMs enriched for both 0D sweeps and pN/pS"
  ) +
  theme_classic()

p_tail_both












both_bar <- both_tail[order(obs_changes.x - null_mean)]
both_bar[, PFAM_primary := factor(PFAM_primary, levels = PFAM_primary)]

ggplot(both_bar, aes(x = PFAM_primary, y = obs_changes.x)) +
  geom_col() +
  geom_point(aes(y = null_mean), shape = 21, fill = "white", size = 2) +
  geom_errorbar(aes(ymin = pmax(null_mean - 1.96 * null_sd, 0),
                    ymax = null_mean + 1.96 * null_sd),
                width = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "# 0D sweep changes",
    title = "Observed vs expected 0D sweeps in PFAMs enriched for both signals",
    subtitle = "Bars = observed; open circles = null mean; whiskers = ±1.96 SD"
  )



















both_tail_all <- both_tail

lab_both <- both_tail[order(-obs_changes.x)][1:min(10, .N)]



lab_both <- merge(
  both_tail,
  combo_full[, .(PFAM_primary, sweep_q, pnps_q, sweep_log2FE, pnps_log2FE)],
  by = "PFAM_primary",
  all.x = TRUE
)

lab_both[, rank_score := frank(sweep_q, ties.method = "average") +
           frank(pnps_q, ties.method = "average")]

lab_both <- lab_both[order(rank_score, -obs_changes.x)][1:min(10, .N)]
p_tail_both <- ggplot(tail_dt, aes(x = k, y = S)) +
  geom_step(linewidth = 1, color = "black") +
  geom_point(
    data = both_tail,
    aes(x = obs_changes.x, y = S),
    inherit.aes = FALSE,
    size = 2.4,
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = lab_both,
    aes(x = obs_changes.x, y = S, label = PFAM_primary),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "k (0D changes per PFAM)",
    y = "Fraction of PFAMs with ≥ k changes",
    title = "Tail of 0D changes per PFAM",
    subtitle = "Red points = PFAMs enriched for both 0D sweeps and pN/pS"
  ) +
  theme_classic()

p_tail_both
both_tail












# "Both" PFAMs already merged with sweep stats in both_tail
# Keep the 10 most extreme by observed 0D sweep count
lab_both <- both_tail[order(-obs_changes.x)][1:min(10, .N)]

lab_both

lab_both <- both_tail[order(-(obs_changes.x - null_mean))][1:min(10, .N)]






# tail from observed PFAM sweep counts
tail_pfam <- function(x) {
  x <- as.integer(x)
  x <- x[is.finite(x)]
  kmax <- max(x)
  data.table(
    k = 1:kmax,
    S = vapply(1:kmax, function(kk) mean(x >= kk), numeric(1))
  )
}

tail_dt <- tail_pfam(res[obs_changes.x > 0]$obs_changes.x)

# all PFAMs enriched for both signals
both_tail <- combo_full[sig_class == "Both", .(PFAM_primary, sig_class)]
both_tail <- merge(
  both_tail,
  res[, .(PFAM_primary, obs_changes.x, n_cols, null_mean, null_sd)],
  by = "PFAM_primary",
  all.x = TRUE
)

both_tail <- both_tail[!is.na(obs_changes.x) & obs_changes.x > 0]
both_tail[, PFAM_label := PFAM_primary]

# map to tail curve
setkey(tail_dt, k)
both_tail[, k_plot := pmin(obs_changes.x, max(tail_dt$k))]
both_tail <- tail_dt[both_tail, on = .(k = k_plot)]

# label only the top 10 most extreme "Both" PFAMs
lab_both <- both_tail[order(-obs_changes.x)][1:min(10, .N)]

p_tail_both <- ggplot(tail_dt, aes(x = k, y = S)) +
  geom_step(linewidth = 1, color = "black") +
  geom_point(
    data = both_tail,
    aes(x = obs_changes.x, y = S),
    inherit.aes = FALSE,
    size = 2.4,
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = lab_both,
    aes(x = obs_changes.x, y = S, label = PFAM_label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "k (0D changes per PFAM)",
    y = "Fraction of PFAMs with ≥ k changes",
    title = "Tail of 0D changes per PFAM",
    subtitle = "Red points = PFAMs enriched for both 0D sweeps and pN/pS; labels show top 10 by observed sweep count"
  ) +
  theme_classic()

p_tail_both





# rank the same 10 PFAMs for inset/barplot
bar_both <- both_tail[order(-(obs_changes.x - null_mean))][1:min(10, .N)]

bar_both[, delta := obs_changes.x - null_mean]
bar_both[, PFAM_primary := factor(PFAM_primary, levels = rev(bar_both$PFAM_primary))]

p_both_obs_exp <- ggplot(bar_both, aes(x = PFAM_primary, y = obs_changes.x)) +
  geom_col() +
  geom_point(aes(y = null_mean), shape = 21, fill = "white", size = 2.3) +
  geom_errorbar(
    aes(
      ymin = pmax(null_mean - 1.96 * null_sd, 0),
      ymax = null_mean + 1.96 * null_sd
    ),
    width = 0.2
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "# 0D sweep changes",
    title = "Observed vs expected 0D sweeps in PFAMs enriched for both signals",
    subtitle = "Top 10 'Both' PFAMs ranked by observed-minus-expected sweeps; bars = observed, open circles = null mean, whiskers = ±1.96 SD"
  )

p_both_obs_exp













changes_shared <- copy(changes)

changes_shared[, month_pair := paste(Month_t1, Month_t2, sep = "->")]

shared_pfam_pairs <- changes_shared[
  !is.na(PFAM_primary) & PFAM_primary != "",
  .(n_cols = uniqueN(colony_id)),
  by = .(PFAM_primary, month_pair)
][n_cols >= 2]

changes_shared <- merge(
  changes_shared,
  shared_pfam_pairs[, .(PFAM_primary, month_pair)],
  by = c("PFAM_primary", "month_pair"),
  all = FALSE
)





changes_shared_genus <- copy(changes)

changes_shared_genus[, month_pair := paste(Month_t1, Month_t2, sep = "->")]

shared_pfam_pairs_genus <- changes_shared_genus[
  !is.na(PFAM_primary) & PFAM_primary != "" & !is.na(Genus),
  .(n_genus = uniqueN(Genus)),
  by = .(PFAM_primary, month_pair)
][n_genus >= 2]

changes_shared_genus <- merge(
  changes_shared_genus,
  shared_pfam_pairs_genus[, .(PFAM_primary, month_pair)],
  by = c("PFAM_primary", "month_pair"),
  all = FALSE
)





obs_pfam_shared0D <- changes_shared[, .(
  obs_changes = .N,
  n_units  = uniqueN(unit_id),
  n_cols   = uniqueN(colony_id),
  n_MAGs   = uniqueN(MAG_id),
  n_genus  = uniqueN(Genus)
), by = PFAM_primary]

n_by_unit_shared <- merge(
  unit_meta[, .(unit_id)],
  changes_shared[, .N, by = unit_id],
  by = "unit_id",
  all.x = TRUE
)

n_by_unit_shared[is.na(N), N := 0L]
setnames(n_by_unit_shared, "N", "n_changes")



pfams_shared <- sort(unique(c(
  obs_pfam_shared0D$PFAM_primary,
  unique(unlist(lapply(unit_univ_list, function(u) u$PFAM_primary)))
)))

pfam_index_shared <- setNames(seq_along(pfams_shared), pfams_shared)

obs_vec_shared <- integer(length(pfams_shared))
obs_vec_shared[pfam_index_shared[obs_pfam_shared0D$PFAM_primary]] <- obs_pfam_shared0D$obs_changes

res_null_shared <- bootstrap_null1_pfam(
  B = 500,
  seed = 1L,
  unit_meta = unit_meta,
  n_by_unit = n_by_unit_shared,
  unit_univ_list = unit_univ_list,
  pfams = pfams_shared,
  pfam_index = pfam_index_shared,
  obs_vec = obs_vec_shared
)

res_shared <- merge(res_null_shared, obs_pfam_shared0D, by = "PFAM_primary", all.x = TRUE)
res_shared[, FE := obs_changes.x / pmax(null_mean, 1e-12)]
res_shared[, sweep_log2FE := log2((obs_changes.x + 1) / (null_mean + 1))]



combo_shared <- merge(
  res_shared[, .(
    PFAM_primary,
    sweep_obs_changes = obs_changes.x,
    sweep_null_mean = null_mean,
    sweep_null_sd = null_sd,
    sweep_p = p_null,
    sweep_q = q_BH,
    sweep_log2FE
  )],
  pnps_keep_full,
  by = "PFAM_primary",
  all = TRUE
)

combo_shared[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_log2FE > 0]
combo_shared[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_log2FE > 0]

combo_shared[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                                    fifelse(sweep_hit, "Sweep only",
                                            fifelse(pnps_hit, "pN/pS only", "Neither")))]
combo_shared[, sig_class := factor(sig_class,
                                   levels = c("Neither","Sweep only","pN/pS only","Both"))]




tail_dt_shared <- tail_pfam(res_shared[obs_changes.x > 0]$obs_changes.x)

both_tail_shared <- combo_shared[sig_class == "Both", .(PFAM_primary, sig_class)]
both_tail_shared <- merge(
  both_tail_shared,
  res_shared[, .(PFAM_primary, obs_changes.x, n_cols, null_mean, null_sd)],
  by = "PFAM_primary",
  all.x = TRUE
)

both_tail_shared <- both_tail_shared[!is.na(obs_changes.x) & obs_changes.x > 0]
both_tail_shared[, PFAM_label := PFAM_primary]

setkey(tail_dt_shared, k)
both_tail_shared[, k_plot := pmin(obs_changes.x, max(tail_dt_shared$k))]
both_tail_shared <- tail_dt_shared[both_tail_shared, on = .(k = k_plot)]

lab_both_shared <- both_tail_shared[order(-obs_changes.x)][1:min(10, .N)]


bar_both_shared <- both_tail_shared[order(-(obs_changes.x - null_mean))][1:min(10, .N)]
bar_both_shared[, PFAM_primary := factor(PFAM_primary, levels = rev(bar_both_shared$PFAM_primary))]



both_pfams_shared <- combo_shared[sig_class == "Both", unique(PFAM_primary)]
sweep_only_pfams_shared <- combo_shared[sig_class == "Sweep only", unique(PFAM_primary)]



g <- copy(gene_clean)
g <- g[!is.na(PFAM_primary) & PFAM_primary != ""]
g <- g[is.finite(pN_pS_ratio)]

eps <- 1e-6
g[, log_pnps := log(pN_pS_ratio + eps)]

g[, gene_class := fifelse(PFAM_primary %chin% both_pfams_shared, "Both",
                          fifelse(PFAM_primary %chin% sweep_only_pfams_shared, "Sweep only",
                                  "Background"))]

g[, gene_class := factor(gene_class, levels = c("Background","Sweep only","Both"))]






combo_shared[, sweep_hit := !is.na(sweep_q) & sweep_q <= 0.1 & sweep_log2FE > 0]
combo_shared[, pnps_hit  := !is.na(pnps_q)  & pnps_q  <= 0.1 & pnps_log2FE  > 0]

combo_shared[, sig_class := fifelse(sweep_hit & pnps_hit, "Both",
                                    fifelse(sweep_hit, "Sweep only",
                                            fifelse(pnps_hit, "pN/pS only", "Neither")))]
combo_shared[, sig_class := factor(sig_class,
                                   levels = c("Neither","Sweep only","pN/pS only","Both"))]

table(combo_shared$sig_class, useNA = "ifany")



tail_pfam <- function(x) {
  x <- as.integer(x)
  x <- x[is.finite(x)]
  x <- x[x > 0]
  kmax <- max(x)
  data.table(
    k = 1:kmax,
    S = vapply(1:kmax, function(kk) mean(x >= kk), numeric(1))
  )
}

tail_dt_shared <- tail_pfam(res_shared[obs_changes.x > 0]$obs_changes.x)



both_tail_shared <- combo_shared[sig_class == "Both", .(PFAM_primary, sig_class)]

both_tail_shared <- merge(
  both_tail_shared,
  res_shared[, .(PFAM_primary, obs_changes.x, n_cols, null_mean, null_sd)],
  by = "PFAM_primary",
  all.x = TRUE
)

both_tail_shared <- both_tail_shared[!is.na(obs_changes.x) & obs_changes.x > 0]
both_tail_shared[, PFAM_label := PFAM_primary]

setkey(tail_dt_shared, k)
both_tail_shared[, k_plot := pmin(obs_changes.x, max(tail_dt_shared$k))]
both_tail_shared <- tail_dt_shared[both_tail_shared, on = .(k = k_plot)]

# label only top 10 by shared observed sweep count
lab_both_shared <- both_tail_shared[order(-obs_changes.x)][1:min(10, .N)]

both_tail_shared
lab_both_shared








p_tail_both_shared <- ggplot(tail_dt_shared, aes(x = k, y = S)) +
  geom_step(linewidth = 1, color = "black") +
  geom_point(
    data = both_tail_shared,
    aes(x = obs_changes.x, y = S),
    inherit.aes = FALSE,
    size = 2.4,
    color = "red"
  ) +
  ggrepel::geom_text_repel(
    data = lab_both_shared,
    aes(x = obs_changes.x, y = S, label = PFAM_label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0
  ) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "k (shared 0D changes per PFAM)",
    y = "Fraction of PFAMs with ≥ k shared changes",
    title = "Tail of shared 0D changes per PFAM",
    subtitle = "Red points = PFAMs enriched for both shared 0D sweeps and pN/pS"
  ) +
  theme_classic()

p_tail_both_shared

bar_both_shared <- both_tail_shared[order(-(obs_changes.x - null_mean))][1:min(10, .N)]
bar_both_shared[, delta := obs_changes.x - null_mean]
bar_both_shared[, PFAM_primary := factor(PFAM_primary, levels = rev(bar_both_shared$PFAM_primary))]

p_both_obs_exp_shared <- ggplot(bar_both_shared, aes(x = PFAM_primary, y = obs_changes.x)) +
  geom_col() +
  geom_point(aes(y = null_mean), shape = 21, fill = "white", size = 2.3) +
  geom_errorbar(
    aes(
      ymin = pmax(null_mean - 1.96 * null_sd, 0),
      ymax = null_mean + 1.96 * null_sd
    ),
    width = 0.2
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "# shared 0D sweep changes",
    title = "Observed vs expected shared 0D sweeps in PFAMs enriched for both signals",
    subtitle = "Top 10 'Both' PFAMs ranked by observed-minus-expected shared sweeps"
  )

p_both_obs_exp_shared











plot_dt_shared <- combo_shared[!is.na(sweep_log2FE) & !is.na(pnps_log2FE)]

lab_scatter_shared <- plot_dt_shared[
  sig_class == "Both" & sweep_log2FE > 0 & pnps_log2FE > 0
][order(-(sweep_log2FE + pnps_log2FE))][1:min(10, .N)]

p_scatter_shared <- ggplot(plot_dt_shared, aes(x = sweep_log2FE, y = pnps_log2FE)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(fill = sig_class, size = pmin(sweep_obs_changes, 1600)),
    shape = 21, color = "black", alpha = 0.7, stroke = 0.2
  ) +
  ggrepel::geom_text_repel(
    data = lab_scatter_shared,
    aes(label = PFAM_primary),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_classic() +
  labs(
    x = "log2 fold-enrichment of shared 0D changes",
    y = "log2 fold-enrichment of pN/pS>1 genes",
    fill = "Significance class",
    size = "# shared 0D changes"
  )

p_scatter_shared








both_pfams_shared <- combo_shared[sig_class == "Both", unique(PFAM_primary)]
sweep_only_pfams_shared <- combo_shared[sig_class == "Sweep only", unique(PFAM_primary)]

g_shared <- copy(gene_clean)
g_shared <- g_shared[!is.na(PFAM_primary) & PFAM_primary != ""]
g_shared <- g_shared[is.finite(pN_pS_ratio)]

eps <- 1e-6
g_shared[, log_pnps := log(pN_pS_ratio + eps)]

g_shared[, gene_class := fifelse(PFAM_primary %chin% both_pfams_shared, "Both",
                                 fifelse(PFAM_primary %chin% sweep_only_pfams_shared, "Sweep only",
                                         "Background"))]

g_shared[, gene_class := factor(gene_class, levels = c("Background","Sweep only","Both"))]




sum_dt_shared <- g_shared[, .(
  n = .N,
  mean_pnps = mean(pN_pS_ratio, na.rm = TRUE),
  median_pnps = median(pN_pS_ratio, na.rm = TRUE),
  mean_log_pnps = mean(log_pnps, na.rm = TRUE),
  median_log_pnps = median(log_pnps, na.rm = TRUE),
  frac_gt1 = mean(pN_pS_ratio > 1, na.rm = TRUE)
), by = gene_class]

sum_dt_shared






# use the same top 10 you already selected
dumbbell_dt <- copy(bar_both_shared)

# order by observed - expected
dumbbell_dt[, delta := obs_changes.x - null_mean]
dumbbell_dt <- dumbbell_dt[order(delta)]
dumbbell_dt[, PFAM_primary := factor(PFAM_primary, levels = dumbbell_dt$PFAM_primary)]

p_both_dumbbell_shared <- ggplot(dumbbell_dt, aes(y = PFAM_primary)) +
  geom_segment(
    aes(x = null_mean, xend = obs_changes.x, yend = PFAM_primary),
    linewidth = 0.8,
    color = "grey50"
  ) +
  geom_point(aes(x = null_mean), size = 2.5, shape = 21, fill = "white", stroke = 0.8) +
  geom_point(aes(x = obs_changes.x), size = 2.8, color = "red") +
  theme_classic() +
  labs(
    x = "# shared 0D sweep changes",
    y = NULL,
    title = "Observed vs expected shared 0D sweeps in PFAMs enriched for both signals",
    subtitle = "Top 10 'Both' PFAMs ranked by observed-minus-expected shared sweeps"
  )

p_both_dumbbell_shared













## ---------- GLOBAL PARAMETERS ----------
gens_per_day   <- 12
coverage_min   <- L
min_intervals  <- 1L

month_base_times <- c(
  "May"       = 0,
  "June"      = 36,
  "July"      = 31,
  "August"    = 34,
  "September" = 28,
  "October"   = 22,
  "November"  = 42,
  "January"   = 47,
  "February"  = 31
)

cumulative_times <- cumsum(month_base_times)
month_to_cumulative_time_map <- setNames(cumulative_times, names(month_base_times))




## =========================================================
## 0) START FROM changes
## =========================================================
dt <- as.data.table(changes)

## Build cumulative time from Month_t1 and Month_t2
dt[, Month_t1 := as.character(Month_t1)]
dt[, Month_t2 := as.character(Month_t2)]

dt[, t1_cum := unname(month_to_cumulative_time_map[Month_t1])]
dt[, t2_cum := unname(month_to_cumulative_time_map[Month_t2])]

## TimeLag = cumulative distance from baseline-derived month positions
dt[, TimeLag := abs(t2_cum - t1_cum)]

## Keep only usable 0D rows
dt <- dt[
  degeneracy == "0D" &
    !is.na(Genus) & Genus != "" &
    !is.na(unit_id) &
    is.finite(TimeLag)
]

## Define a unique 0D sweep event within each comparison unit
dt[, event_id := paste(unit_id, unique_SNV_identifier, sep = "|")]

## Keep one row per unique event
evt <- unique(
  dt[, .(
    event_id,
    unit_id,
    Genus,
    TimeLag
  )],
  by = "event_id"
)

## =========================================================
## 1) COUNT UNIQUE 0D SWEEPS PER COMPARISON UNIT
## =========================================================
unit_counts <- evt[, .(
  n_0D_sweeps = uniqueN(event_id),
  Genus = first(Genus),
  TimeLag = first(TimeLag)
), by = unit_id]

cat("Unique 0D sweep events:", nrow(evt), "\n")
cat("Comparison units:", nrow(unit_counts), "\n")
cat("Genera:", uniqueN(unit_counts$Genus), "\n")

## =========================================================
## 2) SUMMARIZE BY Genus × TimeLag
## =========================================================
genus_tlag <- unit_counts[, .(
  mean_0D_sweeps = mean(n_0D_sweeps, na.rm = TRUE),
  n_units = .N
), by = .(Genus, TimeLag)]

setorder(genus_tlag, Genus, TimeLag)

## =========================================================
## 3) GLOBAL SUMMARY ACROSS ALL GENERA
## =========================================================
global_tlag <- unit_counts[, .(
  mean_0D_sweeps = mean(n_0D_sweeps, na.rm = TRUE),
  n_units = .N
), by = TimeLag][order(TimeLag)]

## =========================================================
## 4) COMPLETE MISSING TimeLag BINS
## =========================================================
all_tlags <- sort(unique(unit_counts$TimeLag))
all_genera <- sort(unique(unit_counts$Genus))

genus_tlag_full <- CJ(Genus = all_genera, TimeLag = all_tlags)
genus_tlag_full <- merge(
  genus_tlag_full,
  genus_tlag,
  by = c("Genus", "TimeLag"),
  all.x = TRUE
)










genus_tlag_full[is.na(mean_0D_sweeps), mean_0D_sweeps := 0]
genus_tlag_full[is.na(n_units), n_units := 0L]

global_tlag_full <- data.table(TimeLag = all_tlags)
global_tlag_full <- merge(global_tlag_full, global_tlag, by = "TimeLag", all.x = TRUE)
global_tlag_full[is.na(mean_0D_sweeps), mean_0D_sweeps := 0]
global_tlag_full[is.na(n_units), n_units := 0L]

## =========================================================
## 5) COMPUTE RUNNING MEANS
## =========================================================
k <- 3

genus_tlag_full[, runmean_0D := zoo::rollmean(
  mean_0D_sweeps,
  k = k,
  fill = NA,
  align = "center"
), by = Genus]

global_tlag_full[, runmean_0D := zoo::rollmean(
  mean_0D_sweeps,
  k = k,
  fill = NA,
  align = "center"
)]

## =========================================================
## 6) OPTIONAL: RESTRICT TO TOP GENERA FOR READABILITY
## =========================================================
top_genera <- unit_counts[, .N, by = Genus][order(-N)][1:min(8, .N), Genus]
plot_genus <- genus_tlag_full[Genus %chin% top_genera]

cat("Top genera shown:\n")
print(top_genera)

## =========================================================
## 7) MAIN PLOT
## =========================================================
p_timelag_genus <- ggplot() +
  geom_line(
    data = plot_genus,
    aes(x = TimeLag, y = mean_0D_sweeps, color = Genus, group = Genus),
    linetype = "dashed",
    linewidth = 0.6,
    alpha = 0.45
  ) +
  geom_point(
    data = plot_genus[mean_0D_sweeps > 0],
    aes(x = TimeLag, y = mean_0D_sweeps, color = Genus),
    size = 1.6,
    alpha = 0.65
  ) +
  geom_line(
    data = plot_genus,
    aes(x = TimeLag, y = runmean_0D, color = Genus, group = Genus),
    linewidth = 0.95,
    alpha = 0.95
  ) +
  geom_line(
    data = global_tlag_full,
    aes(x = TimeLag, y = mean_0D_sweeps),
    color = "black",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_line(
    data = global_tlag_full,
    aes(x = TimeLag, y = runmean_0D),
    color = "black",
    linewidth = 1.3
  ) +
  theme_classic() +
  labs(
    x = "Cumulative TimeLag from baseline (days)",
    y = "Mean # unique 0D sweeps per comparison",
    color = "Genus",
    title = "Temporal dynamics of 0D sweeps by genus",
    subtitle = "Dashed lines = observed means at each TimeLag; solid lines = running means; black = global"
  )

p_timelag_genus




gene_sites <- unique(
  gene_dt[, .(
    MAG_id,
    colony_id,
    corresponding_gene_call,
    Season,
    Genus,
    nN_gene_reference,
    nS_gene_reference
  )],
  by = c("MAG_id", "colony_id", "corresponding_gene_call", "Season")
)


gene_sites <- gene_sites[
  !is.na(nN_gene_reference) &
    nN_gene_reference > 0
]



sweeps_gene_0D <- sweep_dt[
  degeneracy == "0D",
  .N,
  by = .(
    Season = Season_t2,
    MAG_id,
    colony_id,
    corresponding_gene_call
  )
]

setnames(sweeps_gene_0D, "N", "n_sweeps_0D")



gene_rates <- merge(
  gene_sites,
  sweeps_gene_0D,
  by = c("Season", "MAG_id", "colony_id", "corresponding_gene_call"),
  all.x = TRUE
)

gene_rates[is.na(n_sweeps_0D), n_sweeps_0D := 0L]


gene_rates[, rate_0D_per_1e6 :=
             n_sweeps_0D / nN_gene_reference * 1e6
]



fit_0D <- glm(
  n_sweeps_0D ~ Season,
  offset = log(nN_gene_reference),
  family = poisson,
  data = gene_rates
)

summary(fit_0D)
exp(coef(fit_0D))


dispersion <- sum(residuals(fit_0D, type="pearson")^2) / fit_0D$df.residual
dispersion



fit_0D <- glm(
  n_sweeps_0D ~ Season,
  offset = log(nN_gene_reference),
  family = quasipoisson,
  data = gene_rates
)


fit_0D_genus <- glm(
  n_sweeps_0D ~ Season + Genus,
  offset = log(nN_gene_reference),
  family = quasipoisson,
  data = gene_rates
)

summary(fit_0D_genus)
exp(coef(fit_0D_genus))

season_summary <- gene_rates[, .(
  sweeps = sum(n_sweeps_0D),
  sites  = sum(nN_gene_reference)
), by = Season]

season_summary[, rate_0D_per_1e6 := sweeps / sites * 1e6]

season_summary


ggplot(season_summary,
       aes(x = Season, y = rate_0D_per_1e6)) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  labs(
    y = "0D sweeps per 10^6 nonsynonymous sites",
    x = NULL,
    title = "Seasonal bias in 0D sweep rates"
  )





ggplot(gene_rates,
       aes(x = Season, y = rate_0D_per_1e6)) +
  geom_boxplot(outlier.alpha = 0.25) +
  geom_jitter(alpha = 0.1,width = 0.05)+
  scale_y_log10() +
  theme_classic() +
  labs(
    y = "0D sweeps per 10^6 nonsynonymous sites",
    title = "Distribution of gene-level 0D sweep rates by season"
  )


ggplot(gene_rates,
       aes(x = Season, y = rate_0D_per_1e6)) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_jitter(width = 0.05, alpha = 0.15, size = 0.5) +
  scale_y_log10() +
  theme_classic() +
  labs(
    y = "0D sweeps per 10^6 nonsynonymous sites",
    title = "Distribution of gene-level 0D sweep rates by season"
  )


ggplot(gene_rates,
       aes(x = Season, y = rate_0D_per_1e6)) +
  geom_violin(fill = "lightgray") +
  geom_boxplot(width = 0.15) +
  geom_jitter(width = 0.05, alpha = 0.15, size = 0.5) +
  scale_y_log10() +
  theme_classic()

fit_0D <- glm(
  n_sweeps_0D ~ Season + Genus,
  offset = log(nN_gene_reference),
  family = quasipoisson,
  data = gene_rates
)




pred_data <- data.table(
  Season = unique(gene_rates$Season),
  Genus = levels(factor(gene_rates$Genus))[1],   # hold genus constant
  nN_gene_reference = 1
)



pred <- predict(
  fit_0D,
  newdata = pred_data,
  type = "link",
  se.fit = TRUE
)

pred_data[, `:=`(
  log_rate = pred$fit,
  se = pred$se.fit
)]

pred_data[, `:=`(
  rate = exp(log_rate),
  lower = exp(log_rate - 1.96 * se),
  upper = exp(log_rate + 1.96 * se)
)]


pred_data[, `:=`(
  rate_per_1e6 = rate * 1e6,
  lower_per_1e6 = lower * 1e6,
  upper_per_1e6 = upper * 1e6
)]




ggplot(pred_data,
       aes(x = Season, y = rate_per_1e6)) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = lower_per_1e6, ymax = upper_per_1e6),
    width = 0.2
  ) +
  theme_classic() +
  labs(
    y = "0D sweeps per 10^6 nonsynonymous sites",
    x = NULL,
    title = "Seasonal bias in 0D sweep rates",
    subtitle = "Model-estimated rates with 95% confidence intervals"
  )+ geom_point(
    data = season_summary,
    aes(x = Season, y = rate_per_1e6),
    color = "red",
    size = 3
  )



season_summary <- gene_rates[, .(
  sweeps = sum(n_sweeps_0D),
  sites = sum(nN_gene_reference)
), by = Season]

season_summary[, rate_per_1e6 := sweeps / sites * 1e6]






snv_raw2 <- snv_raw[snv_raw$degeneracy == "0D",]

print(length(unique(snv_raw$corresponding_gene_call)))

print(30344/5)
