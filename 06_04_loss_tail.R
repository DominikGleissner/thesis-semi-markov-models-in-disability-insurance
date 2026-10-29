###############################################################################
##
##  06_04_loss_tail.R
##
##  Section 6.4. The distribution of the policy loss L = L^h_{(a,0)}(t_0) and
##  its tail through the Monte Carlo. Everything is reported under the misspecified 
##  pricing E2, so every scenario is evaluated with the baseline premium.
##  A pool analysis for independent and identically distributed losses is done as well.
##
##  Input
##    INPUT_DIR   dav2021i_invalid_male_2ndorder.csv
##    OUTPUT_DIR  baseline_calibration_grid.csv,
##                sojourn_data_fit.csv,
##                sojourn_data_fit_matched.csv
##
##  Output
##    FIG_DIR     fig_loss_tail.pdf,
##                fig_loss_pooling.pdf
##    TAB_DIR     tab_tail_summary.tex, 
##                tail_numbers.tex
##    OUTPUT_DIR  loss_tail_summary.csv
##
##  Needs ggplot2 for the figure.
##
###############################################################################

library(ggplot2)

source("06_00_common.R")
source("06_00_calibration.R")
source("06_00_scenarios.R")
source("06_00_engine.R")

## ---------------------------------------------------------------------------
## 0. Configuration
## ---------------------------------------------------------------------------

SRC        <- "06_04_loss_tail.R"
TEX_PREFIX <- "Tail"
TEX_SCI    <- FALSE

POLICY      <- list(t0 = 30, n = 35, beta = 1, f = 0, delta = 0.03, h = 1 / 24)
POLICY_LATE <- modifyList(POLICY, list(t0 = 45, n = 20))
POLICY_DEF  <- modifyList(POLICY, list(f = 1))

M         <- 1e6             # path number of Definition 6.1.1
ALPHA_VAR <- 0.995           # Solvency II level of the Value-at-Risk
ALPHA_ES  <- 0.99            # Swiss Solvency Test level of the shortfall
GAMMA     <- 0.05            # confidence level of (5.29)

SEED       <- 19690902
SEED_EXTRA <- 2              # further seeds of the noise diagnostic
SEED_KEYS  <- c("baseline", paste0("weibull", KAPPA),
                "lognormal", "empirical", "fittedmm", "fitted")
VAR_KEYS   <- SEED_KEYS
POOL_KEYS  <- c("weibull0.5", "weibull2", "lognormal", "empirical")

POOL_K   <- c(1, 10, 100, 1000)    # pool sizes of Proposition 6.4.2
CV_AGE   <- 45                     # age at which the baseline cv is reported
TAIL_LO  <- 0.90                   # left end of the tail figure
TAIL_PTS <- 500                    # resolution of the tail figure
TAIL_PAD <- 0.02                   # left pad of the tail figure, in range units


## ---------------------------------------------------------------------------
## 1. Estimators of Definition 5.4.3
## ---------------------------------------------------------------------------

## Mean and empirical variance of (5.26), skewness and excess kurtosis
## as of below the definition, empirical Value-at-Risk of (5.27) and empirical
## expected shortfall (5.28).
emp_stats <- function(x, a_var = ALPHA_VAR, a_es = ALPHA_ES) {
  n <- length(x)
  if (n <= 1 / (1 - a_es))
    stop("the sample is too small for the level ", a_es)
  m <- mean(x)
  s <- sd(x)
  o <- sort(x)
  k <- floor(n * (1 - a_es))
  c(mean = m, sd = s,
    skew = mean((x - m)^3) / s^3,
    kurt = mean((x - m)^4) / s^4 - 3,
    var  = o[ceiling(n * a_var)],
    es   = mean(o[(n - k + 1):n]))
}

## Half-width of the interval (5.29) for the mean.
half_width <- function(x, gamma = GAMMA)
  qnorm(1 - gamma / 2) * sd(x) / sqrt(length(x))

## Gap between the simulated mean bias and the prediction of assertion (3) of 
## Corollary 6.4.1. Returns the gap and the half-width, one column per scenario.
mean_check <- function(run, keys) {
  base  <- run[["baseline"]]$e2
  prem0 <- run[["baseline"]]$pr$premium
  vapply(keys, function(k) {
    d <- run[[k]]$e2 - base
    c(gap = mean(d) - run[[k]]$pr$a_bar * (run[[k]]$pr$premium - prem0),
      hw  = half_width(d))
  }, numeric(2))
}

## Per-policy shortfall of a pool of K independent copies, obtained by
## regrouping the sample into floor(length(x) / K) blocks.
pool_es <- function(x, K, a_es = ALPHA_ES) {
  n <- length(x) %/% K
  if (n <= 1 / (1 - a_es)) return(NA_real_)
  s <- colSums(matrix(x[seq_len(n * K)], nrow = K))
  o <- sort(s)
  k <- floor(n * (1 - a_es))
  mean(o[(n - k + 1):n]) / K
}

## Empirical survival function (5.25) on a grid.
tail_curve <- function(x, grid) 1 - findInterval(grid, sort(x)) / length(x)

## Coefficient of variation of the baseline sojourn law at the integer entry
## age s.  For the first two moments, we use E[T^n] = n * int_0^infty z^{n-1} P(T>z) dz
## as given in Durrett Rick, Probability: Theory and Examples, 2019 (Lemma 2.2.13),
## Folland Gerald, Real Analysis: Modern Techniques and Their Applications, 1999
## (Proposition 6.24) or Dickson David et al., Actuarial mathematics for life
## contingent risks, 2013 (Section 2.5, Equations (2.22) and (2.23)).
## Under the baseline the exit intensity is the constant q on every year of age.
## Recall that the baseline survival of (5.31) states S(z) = exp(- int_0^z q(s+v) dv).
## On the k-th year of the invalidity with z in [k, k + 1), the survival is therefore
## S(z) = S(k) * exp(- q(k) * (z - k)) with S(k) the survival at the start of that year.
## Similar to step 1 in the proof of Proposition 6.1.6, we deduce that
##   int_k^{k+1} S(z) dz = S_k (1 - e^{-q}) / q.
## A substitution w = z - k and integration by parts further yields
##   int_k^{k+1} z S(z) dz = S(k) * int_0^1 (k+w) exp(-qw) dw
##     = S(k) * [ k * int_0^1 exp(-qw) dw + int_0^1 w exp(-qw) dw ]
##     = ... = S_k (k (1 - e^{-q}) / q + (1 - e^{-q} (1 + q)) / q^2).
## Beyond the last calibrated age the intensity is continued constantly, so from
## the duration K on the two integrals are S_K / q and S_K (K / q + 1 / q^2).
## Both moments are thus exact, and neither a quadrature nor a truncation enters.
baseline_cv <- function(cal, s) {
  stopifnot(abs(s - round(s)) < 1e-9)
  x <- s:max(cal$age)
  q <- q_at_age(cal, x) # vector of constant rates for all ages
  e <- exp(-q)
  S <- c(1, exp(-cumsum(q)))[seq_along(q)] # survival at the start of each year
  k <- x - s
  m1 <- sum(S * (1 - e) / q)                                      # first moment
  m2 <- 2 * sum(S * (k * (1 - e) / q + (1 - e * (1 + q)) / q^2))  # second moment
  ## the tail from the duration K on, where the intensity stays at its last value
  ST <- exp(-sum(q))
  qT <- q[length(q)]
  K  <- length(q)
  m1 <- m1 + ST / qT
  m2 <- m2 + 2 * ST * (K / qT + 1 / qT^2)
  ## Consistency check: computed mean must reproduce the tabulated one.
  stopifnot(abs(m1 - cal$m_hat[match(s, cal$age)]) < 1e-6)
  ## Well-known formula Var(T) = E[T^2] - E[T]^2.
  sqrt(m2 - m1^2) / m1
}


## ---------------------------------------------------------------------------
## 2. The simulation
## ---------------------------------------------------------------------------

## Returns the losses of both experiments E1 an E2 together with the deterministic
## quantities the checks of the section compare them against.
## The baseline has to be part of keys, since it fixes the premium of E2.
## The loss is computed via the discounted occupancies as in (6.32) in Corollary 6.4.1.
run_set <- function(cal, pol, scen, keys, m_paths, seed) {
  if (!"baseline" %in% keys) stop("the baseline fixes the premium of E2")
  mods  <- lapply(scen[keys], function(g) build_general(cal, pol, g))
  ## Compute the deterministic premiums and extract the baseline premium.
  prs   <- lapply(mods, price_g, pol = pol)
  prem0 <- prs[["baseline"]]$premium

  out <- list()
  for (k in keys) {
    ## Reset the seed for each scenario for higher comparability.
    set.seed(seed)
    occ <- simulate_occupancy(mods[[k]], pol, m_paths)
    ## Extract the Monte-Carlo premium and the two losses for E1 and E2.
    out[[k]] <- list(pr = prs[[k]],
                     e2 = pol$beta * occ$i - prem0            * occ$a,
                     e1 = pol$beta * occ$i - prs[[k]]$premium * occ$a)
  }
  out
}

## The three bias functionals of (6.4) of one run.
bias_table <- function(run, scale, which = "e2") {
  ref <- emp_stats(run[["baseline"]][[which]])
  do.call(rbind, lapply(names(run), function(k) {
    st <- emp_stats(run[[k]][[which]])
    data.frame(scenario = k, mean = st[["mean"]], sd = st[["sd"]],
               skew = st[["skew"]], kurt = st[["kurt"]],
               var = st[["var"]], es = st[["es"]],
               d_sd  = (st[["sd"]]   - ref[["sd"]])   / scale,
               ## Scaling invariance of skewness and kurtosis, since
               ## mu_k / sigma^k cancel the scaling
               d_skew = st[["skew"]] - ref[["skew"]],
               d_kurt = st[["kurt"]] - ref[["kurt"]],
               d_e   = (st[["mean"]] - ref[["mean"]]) / scale,
               d_var = (st[["var"]]  - ref[["var"]])  / scale,
               d_es  = (st[["es"]]   - ref[["es"]])   / scale,
               stringsAsFactors = FALSE)
  }))
}


## ---------------------------------------------------------------------------
## 3. Figures 6.7 and 6.8
## ---------------------------------------------------------------------------

## Figure 6.7
plot_loss_tail <- function(curves, labels) {
  sc <- scenario_scales(labels)
  curves$label <- factor(unlist(labels[curves$key], use.names = FALSE),
                         levels = sc$levels)

  ggplot(curves, aes(x, p, colour = label, linetype = label)) +
    geom_hline(yintercept = c(1 - ALPHA_VAR, 1 - ALPHA_ES),
               colour = COL_LINE_D, linewidth = 0.55) +
    geom_step(linewidth = 0.4) +
    scale_y_log10(limits = c(5e-5, 1)) +
    scale_scenario(sc) +
    labs(x = "policy loss L", y = "empirical survival of L") +
    theme_thesis()
}

## Figure 6.8
plot_pooling <- function(pool, mean_lv, labels) {
  lev <- unlist(labels[POOL_KEYS], use.names = FALSE)
  sc  <- scenario_scales(labels)
  pool$label    <- factor(unlist(labels[pool$key],    use.names = FALSE), lev)
  mean_lv$label <- factor(unlist(labels[mean_lv$key], use.names = FALSE), lev)

  ggplot(pool, aes(K, 100 * bias, colour = label, linetype = label)) +
    geom_hline(data = mean_lv, aes(yintercept = 100 * d_e, colour = label),
               linetype = "dotted", linewidth = 0.55, show.legend = FALSE) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.3, show.legend = FALSE) +
    scale_x_log10(breaks = POOL_K) +
    scale_colour_manual(values = sc$colour[lev]) +
    scale_linetype_manual(values = sc$linetype[lev]) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    labs(x = "pool size K", y = "per-policy shortfall bias (%)") +
    theme_thesis()
}


## ---------------------------------------------------------------------------
## 4. Driver
## ---------------------------------------------------------------------------

main <- function() {
  ## Read the already computed scenarios and fits
  cal <- read.csv(file.path(OUTPUT_DIR, "baseline_calibration_grid.csv"),
                  stringsAsFactors = FALSE)
  cal <- cal[order(cal$age), ]
  need <- c("age", "sigma", "mu", "q_i", "p_ia_hat", "m_hat")
  if (!all(need %in% names(cal)))
    stop("the calibrated grid must have the columns ", paste(need, collapse = ", "))

  fit_file <- file.path(OUTPUT_DIR, "sojourn_data_fit.csv")
  mm_file  <- file.path(OUTPUT_DIR, "sojourn_data_fit_matched.csv")
  if (!file.exists(fit_file) || !file.exists(mm_file))
    stop("run 06_01_sojourn_comparison.R first, it writes ", fit_file, " and ", mm_file)
  fit    <- read.csv(fit_file, stringsAsFactors = FALSE)
  fit_mm <- read.csv(mm_file,  stringsAsFactors = FALSE)
  
  ## The empirical scenario reads the table itself and not a fitted parameter.
  tab <- read_invalid_table(file.path(INPUT_DIR, "dav2021i_invalid_male_2ndorder.csv"))
  
  # Compute and calibrate the scenarios.
  m_i <- mean_sojourn(cal)
  set <- build_scenarios(m_i, fit, fit_mm, tab)
  keys   <- names(set$scen)
  labels <- set$labels
  
  ## Additional consistency check.
  ## The two files have to belong to the current settings of 06_00_scenarios.R.
  ## Both fits are therefore computed once more and compared with what was read.
  par  <- c("age", "scale", "shape")
  same <- function(a, b)
    isTRUE(all.equal(unlist(a[par]), unlist(b[par]), tolerance = 1e-8))
  if (!same(fit, fit_weibull_grid(tab)) ||
      !same(fit_mm, fit_weibull_matched_grid(tab, m_i)))
    stop("the fits in ", OUTPUT_DIR, " do not belong to the current settings ",
         "of 06_00_scenarios.R, run 06_01_sojourn_comparison.R again")

  ## Main run for simulation under E2.
  run   <- run_set(cal, POLICY, set$scen, keys, M, SEED)
  
  ## Computation of all bias functionals under the two experiments E1 and E2.
  scale <- run[["baseline"]]$pr$epv_b
  tb    <- bias_table(run, scale, "e2")
  tb1   <- bias_table(run, scale, "e1")
  
  ## Assertion (3) of Corollary 6.4.1: the simulated mean must reproduce 
  ## a_bar(S) * (Pi(S) - Pi(S^(0))) without any further run or simulation.
  ## For the baseline, this is useless. Thus, drop it here.
  prem0 <- run[["baseline"]]$pr$premium
  chk   <- tb$scenario[tb$scenario != "baseline"]
  pair  <- mean_check(run, chk)
  dev   <- abs(pair["gap", ])
  miss  <- dev / pair["hw", ]
  issue <- abs(mean(run[["baseline"]]$e2)) / half_width(run[["baseline"]]$e2)
  message("largest mean gap : ", format(max(miss), digits = 3),
          " half-widths, issue condition : ", format(issue, digits = 3))

  ## The two endpoints or bounds of assertion (2) of Corollary 6.4.1.
  N         <- round(POLICY$n / POLICY$h)
  disc      <- exp(-POLICY$delta * POLICY$h * seq_len(N - 1))
  a_cert    <- POLICY$h * sum(disc)
  ceiling_l <- POLICY$beta * POLICY$h * sum(disc[-1])
  floor_l   <- -prem0 * a_cert
  atom      <- mean(run[["baseline"]]$e2 <= floor_l + 1e-12)
  
  ## Further analysis for the data-based scenarios. Beyond the disability year
  ## Y the tabulated law is continued by the ultimate cell, so its intensity
  ## stops reacting to the elapsed duration, while the hazard rate of a Weibull
  ## law keeps moving with it. The criterion (5.44) compares the
  ## two only on the years 1, ..., Y, so no assumption beyond Y are included.
  ## The present comparison is taken at the longest duration the contract can pay for,
  ## that is an invalidity entered at t0 and never left until maturity, which
  ## by assertion (2) of Corollary 6.4.1 is the path attaining the loss ceiling.
  Z_GAP     <- POLICY$n
  surv_tail <- vapply(c("empirical", "fittedmm", "fitted"),
                      function(k) set$scen[[k]](POLICY$t0, Z_GAP), numeric(1))
  gap_fit   <- surv_tail[["empirical"]] / surv_tail[["fitted"]]
  gap_mm    <- surv_tail[["empirical"]] / surv_tail[["fittedmm"]]
  message("survival ratio at duration ", format(Z_GAP),
          " : free ", format(gap_fit, digits = 3),
          ", mean-matched ", format(gap_mm, digits = 3))
  ## A ratio at or below one would mean that a fitted law outlives the tabulated
  ## one at the duration feeding the tail, which is impossible with the negative
  ## Delta_VaR of both fits. It is therefore some error that requires attention.
  if (min(gap_fit, gap_mm) <= 1)
    stop("a fitted sojourn law survives longer than the tabulated one at the ",
         "duration ", format(Z_GAP), "; the data-based item of Section 6.4 is wrong")
  
  ## How far the three models fall short of the tabulated 99.5% quantile. The
  ## bias functionals (6.4) are measured against the Markov baseline.
  ## Further, we want to make some statements about the comparison of the fits
  ## against the tabulated law to assess how well the fits match the real data.
  var_emp    <- tb$var[tb$scenario == "empirical"]
  var_base   <- tb$var[tb$scenario == "baseline"]
  var_fits   <- tb$var[match(c("fittedmm", "fitted"), tb$scenario)]
  short_base <- (var_emp - var_base) / var_emp
  short_fit  <- max((var_emp - var_fits) / var_emp)
  message("shortfall against the tabulated VaR : baseline ",
          format(100 * short_base, digits = 3), "%, worst fit ",
          format(100 * short_fit, digits = 3), "%")
  ## The prose asserts two directions: that the baseline already understates
  ## the tabulated quantile, and that both fits understate it by more. Neither
  ## is forced by the construction, so both are checked rather than assumed.
  if (short_base <= 0 || min((var_emp - var_fits) / var_emp) <= short_base)
    stop("the ordering of the three Value-at-Risk levels has changed; the ",
         "closing sentence of the data-based item has to be rewritten")

  ## Noise diagnostic for different seeds
  seed_es <- matrix(NA_real_, length(SEED_KEYS), 1 + SEED_EXTRA,
                    dimnames = list(SEED_KEYS, NULL))
  seed_es[, 1] <- tb$d_es[match(SEED_KEYS, tb$scenario)]
  ## Check whether the mean assertion still holds for the other seeds.
  ## Value measured in half-widths of its confidence interval from (5.29).
  seed_miss <- max(miss[setdiff(SEED_KEYS, "baseline")])
  for (j in seq_len(SEED_EXTRA)) {
    rj <- run_set(cal, POLICY, set$scen, SEED_KEYS, M, SEED + j)
    seed_es[, j + 1] <- bias_table(rj, scale, "e2")$d_es
    mj        <- mean_check(rj, setdiff(SEED_KEYS, "baseline"))
    seed_miss <- max(seed_miss, abs(mj["gap", ]) / mj["hw", ])
    ## Remove the computation to release storage.
    rm(rj)
  }
  ## Record how the tail results move between the seeds.
  seed_spread <- max(apply(seed_es[SEED_KEYS != "baseline", , drop = FALSE],
                           1, function(v) diff(range(v))))

  ## The loss and loss functionals for the policy variations for a late
  ## entry and a deferred period, and the deviation to the standard policy.
  run_def    <- run_set(cal, POLICY_DEF,  set$scen, VAR_KEYS, M, SEED)
  run_late   <- run_set(cal, POLICY_LATE, set$scen, VAR_KEYS, M, SEED)
  tb_def     <- bias_table(run_def,  run_def[["baseline"]]$pr$epv_b,  "e2")
  tb_late    <- bias_table(run_late, run_late[["baseline"]]$pr$epv_b, "e2")
  def_spread <- max(abs(tb_def$d_es - tb$d_es[match(tb_def$scenario, tb$scenario)]))
  ## Remove the computation to release storage.
  rm(run_def, run_late)

  ## Pool analysis of the portfolio with multiple i.i.d. contracts. 
  ## First, compute the baseline expected shortfall with K i.i.d. baseline contracts.
  pool_ref <- vapply(POOL_K, function(K)
    pool_es(run[["baseline"]]$e2, K), numeric(1))
  ## For all scenarios selected for the pool analysis, compute the deviation against
  ## the baseline 
  pool <- do.call(rbind, lapply(POOL_KEYS, function(k)
    data.frame(key = k, K = POOL_K,
               bias = (vapply(POOL_K, function(KK) pool_es(run[[k]]$e2, KK),
                              numeric(1)) - pool_ref) / scale,
               stringsAsFactors = FALSE)))
  ## Assertion (2) of Proposition 6.4.2 holds for every pool size here, since it is a 
  ## multiple of the previous one, so the estimated sequence should decrease.
  ## For each scenario, only store the difference for varying pool sizes.
  up <- max(vapply(POOL_KEYS, function(k)
    max(diff(vapply(POOL_K, function(KK) pool_es(run[[k]]$e2, KK), numeric(1)))),
    numeric(1)))
  message("largest increase of the pooled shortfall : ", format(up, digits = 3))
  ## Bias functionals for a single contract or pools of size K = 1.
  mean_lv     <- data.frame(key = POOL_KEYS, stringsAsFactors = FALSE)
  mean_lv$d_e <- tb$d_e[match(POOL_KEYS, tb$scenario)]
  ## Bias functionals for largest chosen pool size.
  pool_last   <- pool$bias[pool$K == max(POOL_K)]
  ## Maximum over all scenarios of the difference between the biases for a single
  ## policy and for a pool of largest size.
  pool_gap    <- max(abs(pool_last - mean_lv$d_e))
  
  ## How far the resolution of the simulation sits below the smallest bias of the table.
  ## First, write all biases in absolute values in one vector.
  bias_all     <- abs(unlist(tb[tb$scenario != "baseline", c("d_e", "d_var", "d_es")]))
  ## Second, the factor reads as the scale of the smallest reported bias against
  ## the maximal reported deviation measured in EPV.
  ## The factor 10 is more elegant for rounding.
  noise_factor <- floor(10 * min(bias_all) / (max(dev) / scale)) / 10

  ## -------------------------------------------------------------------------
  ## From here on nothing is computed, only written
  ## -------------------------------------------------------------------------

  write.csv(tb, file.path(OUTPUT_DIR, "loss_tail_summary.csv"),
            row.names = FALSE)

  ## Figure preparation.
  lo   <- min(vapply(run, function(r) min(r$e2), numeric(1)))
  hi   <- max(vapply(run, function(r) max(r$e2), numeric(1)))
  lo   <- lo - TAIL_PAD * (hi - lo)
  grid <- seq(lo, hi, length.out = TAIL_PTS)
  
  curves <- do.call(rbind, lapply(keys, function(k)
    data.frame(key = k, x = grid, p = tail_curve(run[[k]]$e2, grid),
               stringsAsFactors = FALSE)))
  curves <- curves[curves$p > 0, ]

  ## Figure 6.7.
  fig1 <- open_tex_figure("fig_loss_tail.pdf", height = 5.50)
  print(plot_loss_tail(curves, labels))
  dev.off()

  ## Figure 6.8.
  fig2 <- open_tex_figure("fig_loss_pooling.pdf", height = 2.60)
  print(plot_pooling(pool, mean_lv, labels))
  dev.off()

  ## Table 6.7.
  ## The row labels are grouped as in Section 6.2, which keeps the first column narrow.
  disp <- c(baseline  = "baseline",
            lognormal = "log-normal",
            empirical = "\\quad empirical",
            fittedmm  = "\\quad mean-matched",
            fitted    = "\\quad free")
  for (k in KAPPA)
    disp[paste0("weibull", k)] <- sprintf("\\quad kappa = %s", format(k))
  ord <- c("baseline", "@Weibull", paste0("weibull", KAPPA), "lognormal",
           "@data-based", "empirical", "fittedmm", "fitted")

  pct_cell <- function(v) paste0("$", tex_fixed(100 * v, 2), "\\%$")
  num_cell <- function(v, d) paste0("$", tex_fixed(v, d), "$")

  body <- NULL
  for (nm in ord) {
    if (substr(nm, 1, 1) == "@") {
      body <- rbind(body, c(substring(nm, 2), rep("", 8)))
    } else {
      r <- tb[match(nm, tb$scenario), ]
      body <- rbind(body, c(disp[nm],
                            num_cell(r$sd, 3), num_cell(r$skew, 2),
                            num_cell(r$kurt, 2), num_cell(r$var, 3),
                            num_cell(r$es, 3),
                            if (nm == "baseline") rep("-", 3) else
                              c(pct_cell(r$d_e), pct_cell(r$d_var),
                                pct_cell(r$d_es))))
    }
  }
  body <- as.data.frame(body, stringsAsFactors = FALSE)

  write_tex_tabular(
    body,
    headers = c("scenario", "$\\widehat{\\sigma}_M$", "skew.", "kurt.",
                "$\\widehat{\\mathrm{VaR}}_{0.995}$",
                "$\\widehat{\\mathrm{ES}}_{0.99}$",
                "$\\Delta_{E}$", "$\\Delta_{\\mathrm{VaR}}$",
                "$\\Delta_{\\mathrm{ES}}$"),
    digits  = rep(NA_integer_, ncol(body)),
    align   = "lrrrrrrrr",
    above   = c(paste("     & \\multicolumn{5}{c}{empirical law of $L$} &",
                      "\\multicolumn{3}{c}{bias} \\\\"),
                "    \\cmidrule(lr){2-6} \\cmidrule(lr){7-9}"),
    path    = file.path(TAB_DIR, "tab_tail_summary.tex"))

  ## Every scalar in Section 6.4.
  pct  <- function(v) paste0(tex_fixed(100 * v, 2), "\\%")
  pick <- function(nm, col) tb[[col]][match(nm, tb$scenario)]
  
  ## Record how often the observed ordering of the biases holds.
  es_below_var <- tb$scenario != "baseline" &
    abs(tb$d_es) <= abs(tb$d_var) + 1e-12 &
    abs(tb$d_e) <= abs(tb$d_var) + 1e-12

  write_tex_macros(c(
    tex_macro("Paths",        "10^{6}"),
    tex_macro("ConfLevel",    format(1-GAMMA)),
    tex_macro("Seed",         format(SEED)),
    tex_macro("SeedRuns",     format(1 + SEED_EXTRA)),
    tex_macro("SeedRunsAdded",format(SEED_EXTRA)),
    tex_macro("ScenarioCount",format(sum(tb$scenario != "baseline"))),
    tex_macro("Ceiling",      tex_fixed(ceiling_l, 3)),
    tex_macro("FloorValue",   tex_fixed(floor_l, 3)),
    tex_macro("Atom",         pct(atom)),
    tex_macro("MeanCheckAbs",   paste0(tex_ceil(100 * max(dev) / scale, 4),"\\%")),
    tex_macro("MeanCheckWorst", tex_ceil(max(miss), 2)),
    tex_macro("MeanCheckSeed",  tex_ceil(seed_miss, 2)),
    tex_macro("IssueCheck",     tex_ceil(issue, 2)),
    tex_macro("NoiseFactor",    tex_fixed(noise_factor, 1)),
    tex_macro("BaselineSd",   tex_fixed(pick("baseline", "sd"), 3)),
    tex_macro("BaselineSkew", tex_fixed(pick("baseline", "skew"), 2)),
    tex_macro("BaselineKurt", tex_fixed(pick("baseline", "kurt"), 2)),
    tex_macro("BaselineVaR",  tex_fixed(pick("baseline", "var"), 3)),
    tex_macro("BaselineES",   tex_fixed(pick("baseline", "es"), 3)),
    tex_macro("EsBelowVarHits",      format(sum(es_below_var))),
    tex_macro("DeltaEWeibullHalf",   pct(pick("weibull0.5", "d_e"))),
    tex_macro("DeltaVaRWeibullHalf", pct(pick("weibull0.5", "d_var"))),
    tex_macro("DeltaESWeibullHalf",  pct(pick("weibull0.5", "d_es"))),
    tex_macro("DeltaEWeibullTwo",    pct(pick("weibull2", "d_e"))),
    tex_macro("DeltaVaRWeibullTwo",  pct(pick("weibull2", "d_var"))),
    tex_macro("DeltaESWeibullTwo",   pct(pick("weibull2", "d_es"))),
    tex_macro("DeltaSdLogNormal",    pct(pick("lognormal", "d_sd"))),
    tex_macro("DeltaELogNormal",     pct(pick("lognormal", "d_e"))),
    tex_macro("DeltaVaRLogNormal",   pct(pick("lognormal", "d_var"))),
    tex_macro("DeltaESLogNormal",    pct(pick("lognormal", "d_es"))),
    tex_macro("DeltaEWeibullOne",    pct(pick("weibull1", "d_e"))),
    tex_macro("DeltaVaRWeibullOne",  pct(pick("weibull1", "d_var"))),
    tex_macro("DeltaESWeibullOne",   pct(pick("weibull1", "d_es"))),
    tex_macro("SdLogNormal",   pct(pick("lognormal", "sd") / pick("baseline", "sd") - 1)),
    tex_macro("EpvB",          tex_fixed(scale, 3)),
    tex_macro("DeltaSdMax",    pct(max(abs(tb$d_sd)))),
    tex_macro("SdRelMax",      pct(max(abs(tb$sd  / pick("baseline", "sd")  - 1)))),
    tex_macro("VaRRelMax",     pct(max(abs(tb$var / pick("baseline", "var") - 1)))),
    tex_macro("ESRelMax",      pct(max(abs(tb$es  / pick("baseline", "es")  - 1)))),
    tex_macro("SkewSpreadMax", tex_fixed(max(abs(tb$d_skew)), 2)),
    tex_macro("KurtSpreadMax", tex_fixed(max(abs(tb$d_kurt)), 2)),
    tex_macro("DeltaVaRMax",   pct(max(abs(tb$d_var)))),
    tex_macro("DeltaESMax",    pct(max(abs(tb$d_es)))),
    tex_macro("SkewLogNormal", tex_fixed(pick("lognormal", "skew"), 2)),
    tex_macro("KurtLogNormal", tex_fixed(pick("lognormal", "kurt"), 2)),
    tex_macro("CvAge",         format(CV_AGE)),
    tex_macro("CvBaseline",    tex_fixed(baseline_cv(cal, CV_AGE), 3)),
    tex_macro("DeltaEEmpirical",   pct(pick("empirical", "d_e"))),
    tex_macro("DeltaVaREmpirical", pct(pick("empirical", "d_var"))),
    tex_macro("DeltaESEmpirical",  pct(pick("empirical", "d_es"))),
    tex_macro("DeltaEFittedMatched",  pct(pick("fittedmm", "d_e"))),
    tex_macro("DeltaESFittedMatched", pct(pick("fittedmm", "d_es"))),
    tex_macro("DeltaEFitted",         pct(pick("fitted", "d_e"))),
    tex_macro("DeltaESFitted",        pct(pick("fitted", "d_es"))),
    tex_macro("DeltaVaRFittedMatched", pct(pick("fittedmm", "d_var"))),
    tex_macro("DeltaVaRFitted",        pct(pick("fitted",   "d_var"))),
    tex_macro("SurvGapFree",           tex_fixed(gap_fit, 2)),
    tex_macro("BaselineVaRShort",      pct(short_base)),
    tex_macro("FitVaRShort",           pct(short_fit)),
    tex_macro("EOneESMax",        pct(max(abs(tb1$d_es)))),
    tex_macro("DeferredESSpread", tex_ceil(100 * def_spread, 2)),
    tex_macro("LateEntryESMax",   pct(max(abs(tb_late$d_es)))),
    tex_macro("SeedSpreadES",     tex_ceil(100 * seed_spread, 2)),
    tex_macro("PoolLarge",        pct(pool_last[POOL_KEYS == "empirical"])),
    tex_macro("PoolGapMax",       tex_ceil(100 * pool_gap, 2))
  ), file.path(TAB_DIR, "tail_numbers.tex"))

  message("tail section finished!")
}

main()
