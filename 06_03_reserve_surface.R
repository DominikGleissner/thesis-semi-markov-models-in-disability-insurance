###############################################################################
##
##  06_03_reserve_surface.R
##
##  Section 6.3. The prospective reserve on the duration lattice D_h of (5.8).
##  Everything runs on the duration-aware deterministic engine of 06_00_engine.R
##  under the self-consistent pricing E1, so every scenario carries its own
##  equivalence premium and the Monte Carlo engine is not needed here.
##
##  Input
##    INPUT_DIR   dav2021i_invalid_male_2ndorder.csv
##    OUTPUT_DIR  baseline_calibration_grid.csv
##
##  Output
##    FIG_DIR     fig_reserve_duration.pdf,
##                fig_reserve_surface.pdf
##    TAB_DIR     tab_reserve_gap.tex, 
##                reserve_numbers.tex
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

SRC        <- "06_03_reserve_surface.R"
TEX_PREFIX <- "Res"
TEX_SCI    <- FALSE

POLICY      <- list(t0 = 30, n = 35, beta = 1, f = 0, delta = 0.03, h = 1 / 24)
POLICY_LATE <- modifyList(POLICY, list(t0 = 45, n = 20))
F_LONG      <- 1                     # deferred period of the comparison run

## For the figures
ENTRY_AGES <- c(30, 45, 60)          # the three panels of fig_reserve_duration
##  At h = 1/24 the lattice carries (N + 1)(N + 2)/2 points, far too many to
##  visit one by one. The surface is therefore read on the coarser grids
##  U_STEP, DRIFT_STEP and HEAT_STEP, while the engine itself always runs at h.
##  All of them are multiples of h, so every point read is a genuine lattice
##  point and nothing is interpolated.
U_STEP     <- 1 / 12                  # duration grid of the profiles
DRIFT_STEP <- 1 / 12                  # grid of the drift and the E2 diagnostic
HEAT_STEP  <- 1 / 12                  # grid of fig_reserve_surface
SURF_KEYS  <- c("weibull0.5", "weibull1", "weibull2", 
                "empirical", "fittedmm", "fitted")

FRESH_AGE  <- 45                     # the reserve level quoted in the text
SPELL_AGE  <- 30                     # entry age of the crossing statements
LATE_U     <- 20                     # the long spell quoted in the text

## The six lattice points of tab_reserve_gap, three per entry age.
GAP_POINTS <- data.frame(s = c(30, 30, 30, 45, 45, 45),
                         u = c( 0, 10, 20,  0,  5, 10))


## ---------------------------------------------------------------------------
## 1. Small helpers
## ---------------------------------------------------------------------------

## Grid index of an age on T_h.
idx <- function(pol, age) round((age - pol$t0) / pol$h)

## One scenario on the grid, priced by itself (experiment E1).
## surv = NULL yields the baseline S^(0).
## Uses the duration-aware versions build_general of build_model and price_g of price
build_case <- function(cal, pol, surv) {
  mod <- build_general(cal, pol, surv)
  list(mod = mod, pr = price_g(mod, pol))
}

## V^{+,h}_{(i,u)}(t) at the lattice points (t, u), given as two vectors.
reserve_at <- function(case, pol, t, u)
  mapply(function(tt, uu)
    reserve_g(case$mod, pol, case$pr, "i", idx(pol, tt), round(uu / pol$h),
              case$pr$premium),
    t, u)

## Surface value of the unit premium family, respectively the annuity factor of (6.29).
annuity_at <- function(case, pol, t, u)
  mapply(function(tt, uu) {
    lt <- idx(pol, tt)
    step_g(case$mod, pol, case$pr$ann$a, case$pr$ann$i, lt,
           lt - round(uu / pol$h), 1, 0)[["i"]]
  }, t, u)

## The lattice points of D_h on a grid of step size step.
lattice_points <- function(pol, step) {
  t <- seq(pol$t0 + step, pol$t0 + pol$n, by = step)
  do.call(rbind, lapply(t, function(tt)
    data.frame(t = tt, u = seq(0, tt - pol$t0, by = step))))
}

## Largest movement of the invalid reserve in u at a fixed attained age,
## relative to the scale EPV^h_B as recorded in (6.3).
drift <- function(case, pol, pts, scale) {
  v <- reserve_at(case, pol, pts$t, pts$u)
  d <- tapply(v, pts$t, function(x) diff(range(x))) / scale
  list(max = max(d), age = as.numeric(names(d))[which.max(d)])
}

## First sign change of a gap profile, linearly interpolated between the two
## grid points that bracket it. Returns NA if the profile keeps its sign.
first_crossing <- function(u, g) {
  k <- which(g[-1] * g[-length(g)] < 0)[1]
  if (is.na(k)) return(NA_real_)
  u[k] + (u[k + 1] - u[k]) * abs(g[k]) / (abs(g[k]) + abs(g[k + 1]))
}

## The zero contour of a gap surface. At a fixed attained age the gap is a
## profile in u, so first_crossing gives the duration at which it changes its sign.
## Attained ages without a sign change are dropped.
zero_contour <- function(pts, g, key) {
  u <- tapply(seq_along(g), pts$t, function(i) first_crossing(pts$u[i], g[i]))
  z <- data.frame(t = as.numeric(names(u)), u = as.numeric(u), key = key,
                  stringsAsFactors = FALSE)
  z[!is.na(z$u), ]
}

## u_0(kappa) of (6.30) in Proposition 6.3.2.
u_zero <- function(kappa, m)
  m * (kappa * gamma(1 + 1 / kappa)^kappa)^(1 / (1 - kappa))

## Percentage cell for a column that write_tex_tabular prints as it stands.
pct_cell <- function(v, digits = 2)
  paste0("$", tex_fixed(100 * v, digits), "\\%$")


## ---------------------------------------------------------------------------
## 2. Figures 6.5 and 6.6
## ---------------------------------------------------------------------------

## Figure 6.5
plot_reserve_duration <- function(prof, labels) {
  sc <- scenario_scales(labels)
  prof$label <- factor(unlist(labels[prof$key], use.names = FALSE),
                       levels = sc$levels)
  prof$panel <- factor(paste0("entry age x = ", prof$entry),
                       levels = paste0("entry age x = ", ENTRY_AGES))
  
  ggplot(prof, aes(u, v, colour = label, linetype = label)) +
    geom_line(linewidth = 0.55) +
    facet_wrap(~ panel, nrow = 1, scales = "free_x") +
    scale_scenario(sc) +
    labs(x = "elapsed disability duration u (years)",
         y = "invalid reserve") +
    theme_thesis()
}

## Figure 6.6
plot_reserve_surface <- function(hm, zero, labels) {
  lev <- unlist(labels[SURF_KEYS], use.names = FALSE)
  hm$label   <- factor(unlist(labels[hm$key],   use.names = FALSE),
                       levels = lev)
  zero$label <- factor(unlist(labels[zero$key], use.names = FALSE),
                       levels = lev)
  
  ggplot(hm, aes(t, u)) +
    geom_tile(aes(fill = gap), width = HEAT_STEP, height = HEAT_STEP) +
    geom_line(data = zero, linewidth = 0.4, colour = "grey25") +
    facet_wrap(~ label, nrow = 2) +
    scale_fill_gradient2(low = COL_TABLE, mid = "white", high = COL_ALT_DARK,
                         midpoint = 0,
                         labels = function(v) paste0(round(100 * v, 1), "%")) +
    labs(x = "attained age t", y = "elapsed disability duration u (years)") +
    theme_thesis() +
    theme(legend.key.height = unit(0.32, "cm"),
          legend.key.width  = unit(1.32, "cm"))
}


## ---------------------------------------------------------------------------
## 3. Driver
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
  
  ## Compute and calibrate the scenarios.
  m_i <- mean_sojourn(cal)
  set  <- build_scenarios(m_i, fit, fit_mm, tab)
  keys   <- names(set$scen)
  
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

  cases <- lapply(set$scen, function(g) build_case(cal, POLICY, g))
  ## Extract the baseline.
  base  <- cases$baseline
  scale <- base$pr$epv_b

  ## The general engine with surv = NULL must reproduce the plain one.
  ## Otherwise the sojourn matrix is not the calibrated kernel.
  plain <- build_model(cal, POLICY)
  stopifnot(abs(price(plain, POLICY)$premium - base$pr$premium) < 1e-10)

  ## Assertion (2) of Proposition 6.3.1:
  ## With the split read at the attained age the invalid reserve does not move in u.
  ref    <- build_model(cal, POLICY, split = "attained")
  pr_ref <- price(ref, POLICY)
  flat   <- vapply(0:10, function(k)
    reserve(ref, POLICY, pr_ref, "i", idx(POLICY, 50), k * round(1 / POLICY$h),
            pr_ref$premium), numeric(1))
  stopifnot(diff(range(flat)) < 1e-10)

  ## Profiles for Figure 6.5 fig_reserve_duration
  ## For each entry age and duration combination in each scenario, 
  ## compute the duration-aware reserve and store in one data.frame.
  prof <- do.call(rbind, lapply(ENTRY_AGES, function(s) {   # iteration entry ages
    du <- seq(0, POLICY$t0 + POLICY$n - s, by = U_STEP)
    do.call(rbind, lapply(keys, function(k)                 # iteration scenarios
      data.frame(entry = s, u = du, key = k,
                 v = reserve_at(cases[[k]], POLICY, s + du, du))))
  }))

  ## Preparation for Table 6.5 tab_reserve_gap
  ## Baseline reserve.
  v0  <- reserve_at(base, POLICY, GAP_POINTS$s + GAP_POINTS$u, GAP_POINTS$u)
  ## Reserve bias functional for all scenarios.
  gap <- lapply(cases, function(cs)
    (reserve_at(cs, POLICY, GAP_POINTS$s + GAP_POINTS$u, GAP_POINTS$u) - v0) /
      scale)

  gap_tab <- data.frame(scenario = unlist(set$labels, use.names = FALSE),
                        stringsAsFactors = FALSE)
  for (j in seq_len(nrow(GAP_POINTS)))
    gap_tab[[paste0("c", j)]] <- pct_cell(vapply(gap, function(g) g[j],
                                                 numeric(1)))

  ## The crossings along a spell entered at SPELL_AGE for scalar numbers
  ## Evaluated for the two extreme Weibull scenarios
  du   <- seq(0, POLICY$t0 + POLICY$n - SPELL_AGE, by = U_STEP)
  vb   <- reserve_at(base, POLICY, SPELL_AGE + du, du)
  gp   <- lapply(cases, function(cs)
    (reserve_at(cs, POLICY, SPELL_AGE + du, du) - vb) / scale)
  cross <- vapply(gp, function(g) first_crossing(du, g), numeric(1))
  u0    <- vapply(c(0.5, 2), function(k) u_zero(k, m_i(SPELL_AGE)), numeric(1))

  ## Preparation for Figure 6.6
  pts_d <- lattice_points(POLICY, DRIFT_STEP)
  dr    <- lapply(cases, function(cs) drift(cs, POLICY, pts_d, scale))

  ## Define the grid and compute the reserve bias for the baseline and all scenarios
  pts_h <- lattice_points(POLICY, HEAT_STEP)
  vh    <- reserve_at(base, POLICY, pts_h$t, pts_h$u)
  gh    <- lapply(cases[keys != "baseline"], function(cs)
    (reserve_at(cs, POLICY, pts_h$t, pts_h$u) - vh) / scale)

  ## Extract the scenario name and the positional index of the worst deviation
  worst_key <- names(gh)[which.max(vapply(gh, function(g) max(abs(g)), numeric(1)))]
  worst_pos <- which.max(abs(gh[[worst_key]]))

  ## Write the reserve biases in one data.frame.
  hm <- do.call(rbind, lapply(SURF_KEYS, function(k)
    data.frame(t = pts_h$t, u = pts_h$u, gap = gh[[k]], key = k)))
  ## Compute the zero contour for all scenarios.
  zero <- do.call(rbind, lapply(SURF_KEYS, function(k)
    zero_contour(pts_h, gh[[k]], k)))

  ## Premium biases for all scenarios and all policies.
  prem <- vapply(cases, function(cs)
    abs(cs$pr$premium / base$pr$premium - 1), numeric(1))
  
  ## Pricing difference between the two experiments E1 and E2.
  e2 <- vapply(cases, function(cs)
    max(abs((cs$pr$premium - base$pr$premium) *
              annuity_at(cs, POLICY, pts_d$t, pts_d$u))) / scale, numeric(1))

  ## Deferred period of one year, at a fresh disability entered at SPELL_AGE.
  pol_f  <- modifyList(POLICY, list(f = F_LONG))
  base_f <- build_case(cal, pol_f, set$scen$baseline)
  vf     <- reserve_at(base_f, pol_f, SPELL_AGE, 0)
  def <- vapply(keys, function(k) {
    cs <- build_case(cal, pol_f, set$scen[[k]])
    g  <- (reserve_at(cs, pol_f, SPELL_AGE, 0) - vf) / base_f$pr$epv_b
    abs(g - gap[[k]][1])
  }, numeric(1))

  ## The late-entry policy of Definition 6.1.1.
  pts_l   <- lattice_points(POLICY_LATE, HEAT_STEP)
  base_l  <- build_case(cal, POLICY_LATE, set$scen$baseline)
  vl      <- reserve_at(base_l, POLICY_LATE, pts_l$t, pts_l$u)
  late <- vapply(keys, function(k) {
    cs <- build_case(cal, POLICY_LATE, set$scen[[k]])
    max(abs(reserve_at(cs, POLICY_LATE, pts_l$t, pts_l$u) - vl)) /
      base_l$pr$epv_b
  }, numeric(1))

  ## -------------------------------------------------------------------------
  ## From here on nothing is computed, only written
  ## -------------------------------------------------------------------------
  
  ## Table 6.5.
  write_tex_tabular(
    gap_tab,
    headers = c("scenario", "$u = 0$", "$u = 10$", "$u = 20$",
                "$u = 0$", "$u = 5$", "$u = 10$"),
    digits  = rep(NA_integer_, ncol(gap_tab)),
    align   = "lrrrrrr",
    above   = c(paste("     & \\multicolumn{3}{c}{$x = 30$} &",
                      "\\multicolumn{3}{c}{$x = 45$} \\\\"),
                "    \\cmidrule(lr){2-4} \\cmidrule(lr){5-7}"),
    path    = file.path(TAB_DIR, "tab_reserve_gap.tex"))
  
  ## Figure 6.5.
  open_tex_figure("fig_reserve_duration.pdf", height = 3.4)
  print(plot_reserve_duration(prof, set$labels))
  dev.off()
  
  ## Figure 6.6.
  open_tex_figure("fig_reserve_surface.pdf", height = 5.40)
  print(plot_reserve_surface(hm, zero, set$labels))
  dev.off()
  
  ## Every scalar in Section 6.3.
  pct <- function(v, digits = 2) paste0(tex_fixed(100 * v, digits), "\\%")
  
  write_tex_macros(c(
    tex_macro("FreshAge",         FRESH_AGE),
    tex_macro("FreshBaseline",    tex_fixed(reserve_at(base, POLICY, FRESH_AGE, 0), 4)),
    tex_macro("SpellAge",         SPELL_AGE),
    tex_macro("LateDuration",     LATE_U),

    tex_macro("BaselineDriftMax", pct(dr$baseline$max, 3)),
    tex_macro("BaselineDriftAge", dr$baseline$age),
    tex_macro("ScenarioDriftMax", pct(max(vapply(dr, function(d) d$max, numeric(1))), 2)),

    tex_macro("GapFreshHalf",     pct(gap[["weibull0.5"]][1])),
    tex_macro("GapLateHalf",      pct(gap[["weibull0.5"]][3])),
    tex_macro("GapFreshTwo",      pct(gap[["weibull2"]][1])),
    tex_macro("GapLateTwo",       pct(gap[["weibull2"]][3])),
    tex_macro("GapEmpiricalFresh", pct(gap[["empirical"]][1])),
    tex_macro("GapEmpiricalLate",  pct(gap[["empirical"]][3])),

    tex_macro("CrossHalf",        tex_fixed(cross[["weibull0.5"]], 2)),
    tex_macro("CrossTwo",         tex_fixed(cross[["weibull2"]], 2)),
    tex_macro("UZeroHalf",        tex_fixed(u0[1], 2)),
    tex_macro("UZeroTwo",         tex_fixed(u0[2], 2)),

    tex_macro("GapMax",           pct(max(abs(gh[[worst_key]])))),
    tex_macro("GapMaxScenario",   set$labels[[worst_key]]),
    tex_macro("GapMaxAge",        tex_fixed(pts_h$t[worst_pos], 1)),
    tex_macro("GapMaxDuration",   tex_fixed(pts_h$u[worst_pos], 1)),

    tex_macro("PremiumMax",       pct(max(prem))),
    tex_macro("ExperimentGap",    pct(max(e2), 3)),
    tex_macro("DeferredSpread",   tex_ceil(100 * max(def), 2)),
    tex_macro("LateEntryGapMax",  pct(max(late)))),
    file.path(TAB_DIR, "reserve_numbers.tex"))

  ## A short summary on the console.
  cat("largest reserve bias  ", pct(max(abs(gh[[worst_key]]))),
      " (", worst_key, ")\n", sep = "")
  cat("baseline drift in u   ", pct(dr$baseline$max, 3), "\n", sep = "")
  cat("sign changes at age ", SPELL_AGE, "\n", sep = "")
  print(round(cross, 2))
  
  message("reserve section finished!")
}

main()
