###############################################################################
##
##  06_01_sojourn_comparison.R
##
##  Section 6.1.3. The scenarios of Chapter 6 and the way their sojourn
##  survival functions deviate from the calibrated baseline S^(0) including
##  the crossings of the sojourn survival functions reported in Proposition 6.1.9.
##
##  Input
#     INPUT_DIR   dav2021i_invalid_male_2ndorder.csv
##    OUTPUT_DIR  baseline_calibration_grid.csv
##
##  Output
##    FIG_DIR     fig_sojourn_window.pdf, 
##                fig_sojourn_annihilation.pdf
##    TAB_DIR     sojourn_numbers.tex
##    OUTPUT_DIR  sojourn_crossing.csv, 
##                sojourn_data_fit.csv,
##                sojourn_data_fit_matched.csv
##
##  Needs ggplot2 for the figure.
##
###############################################################################

library(ggplot2)

source("06_00_common.R")
source("06_00_calibration.R")
source("06_00_scenarios.R")


## ---------------------------------------------------------------------------
## 0. Configuration
## ---------------------------------------------------------------------------

SRC        <- "06_01_sojourn_comparison.R"
TEX_PREFIX <- "Soj"
TEX_SCI    <- FALSE

## Only the entry age and the term of Definition 6.1.1 are used here.
POLICIES <- list(early = list(t0 = 30, n = 35),
                 late  = list(t0 = 45, n = 20))

WINDOW_AGE <- 45                   # entry age shown in the survival figure
Z_MAX      <- 100                  # longest duration scanned for a crossing
DZ         <- 0.05                 # width of the duration grid

CROSS_AGES  <- c(30, 35, 36, 40)   # entry ages of the annihilation figure
CROSS_KAPPA <- 0.75                # the DFR shape whose first crossing jumps
CROSS_ZMAX  <- 40                  # durations shown there


## ---------------------------------------------------------------------------
## 1. The crossing study
## ---------------------------------------------------------------------------

## First sign change of the difference of the two survival functions after the
## entry age x, and whether it happens before the end of the payable window.
## This is the numerical form of Proposition 6.1.9 for the calibrated
## baseline, which is not exponential (text below Proposition 6.1.9).
crossing_study <- function(cal, scen, dz = DZ) {
  ## Define the duration grid.
  z    <- seq(dz, Z_MAX, by = dz)
  ## Build one data frame with all pairs of ages and window.
  pair <- NULL
  for (pn in names(POLICIES)) {
    pol  <- POLICIES[[pn]]
    x    <- pol$t0:(pol$t0 + pol$n - 1)
    pair <- rbind(pair, data.frame(age = x, window = pol$t0 + pol$n - x))
  }
  ## As some pairs occur multiple times, we select the unique pairs and order them.
  pair <- unique(pair)
  pair <- pair[order(pair$age, pair$window), ]

  ## Collect the values across each pair.
  out <- NULL
  for (i in seq_len(nrow(pair))) {
    x   <- pair$age[i]
    eta <- pair$window[i]
    b   <- baseline_survival(cal, x, z)
    for (nm in setdiff(names(scen), "baseline")) {
      ## Difference between the two sojourn survival functions.
      g  <- scen[[nm]](x, z) - b
      ## Determine all indices where sign of g changes happen,
      ## respectively where the two sojourn survival functions intersect.
      k  <- which(sign(g[-1]) * sign(g[-length(g)]) < 0)
      ## Take the first sign change if there is one.
      z1 <- if (length(k)) z[k[1]] else NA_real_
      out <- rbind(out, data.frame(
        age = x, scenario = nm, window = eta, first = z1,
        margin = z1 - eta, constant = is.na(z1) || z1 >= eta,
        stringsAsFactors = FALSE))
    }
  }
  out
}


## ---------------------------------------------------------------------------
## 2. Figures 6.2 and 6.3
## ---------------------------------------------------------------------------

plot_window <- function(cal, scen, labels) {
  sc  <- scenario_scales(labels)
  pol <- POLICIES[["early"]]
  z   <- seq(0, pol$t0 + pol$n - WINDOW_AGE, by = 1 / 12)
  df  <- do.call(rbind, lapply(names(scen), function(nm)
    data.frame(z = z, 
               surv = if (nm == "baseline") baseline_survival(cal, WINDOW_AGE, z)
                      else scen[[nm]](WINDOW_AGE, z),
               scenario = factor(labels[[nm]], levels = sc$levels))))
  
  ggplot(df, aes(x = z, y = surv, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.55) +
    scale_x_continuous(breaks = seq(0, max(z), by = 5)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    scale_scenario(sc) +
    labs(x = "elapsed duration z (years)", y = "sojourn survival") +
    theme_thesis()
}

## The difference of the two survival functions for the DFR shape CROSS_KAPPA at
## four entry ages. At this scale the two curves themselves are visually
## identical, so their difference is drawn instead, together with the end of the
## payable window. Between the entry ages 35 and 36 the positive lobe vanishes.
plot_annihilation <- function(cal, m_i) {
  z   <- seq(0, CROSS_ZMAX, by = 1 / 12)
  sur <- scenario_weibull(m_i, CROSS_KAPPA)
  mat <- POLICIES[["early"]]$t0 + POLICIES[["early"]]$n
  lab <- sprintf("entry age x = %d", CROSS_AGES)
  df  <- do.call(rbind, lapply(CROSS_AGES, function(x)
    data.frame(z = z, g = sur(x, z) - baseline_survival(cal, x, z),
               window = mat - x,
               panel = factor(sprintf("entry age x = %d", x), levels = lab))))

  ggplot(df, aes(x = z, y = g)) +
    geom_hline(yintercept = 0, colour = COL_LINE, linewidth = 0.3) +
    geom_vline(aes(xintercept = window), colour = COL_LINE,
               linetype = "22", linewidth = 0.3) +
    geom_line(colour = COL_TABLE, linewidth = 0.55) +
    facet_wrap(~ panel, nrow = 2) +
    scale_x_continuous(breaks = seq(0, CROSS_ZMAX, by = 10)) +
    labs(x = "elapsed duration z (years)",
         y = "scenario survival minus baseline survival") +
    theme_thesis()
}


## ---------------------------------------------------------------------------
## 3. Driver
## ---------------------------------------------------------------------------

main <- function() {
  cal <- read.csv(file.path(OUTPUT_DIR, "baseline_calibration_grid.csv"),
                  stringsAsFactors = FALSE)
  cal <- cal[order(cal$age), ]
  need <- c("age", "q_i", "m_hat")
  if (!all(need %in% names(cal)))
    stop("the calibrated grid must have the columns ", paste(need, collapse = ", "))

  tab <- read_invalid_table(
    file.path(INPUT_DIR, "dav2021i_invalid_male_2ndorder.csv"))
  if (min(FIT_AGES) < tab$x_min || max(FIT_AGES) > tab$omega)
    stop("FIT_AGES leaves the age range of the invalid table")

  m_i <- mean_sojourn(cal)

  ## the recomputed mean has to reproduce the tabulated one of Section 6.1.2
  dev_m <- max(abs(m_i(cal$age) - cal$m_hat))
  message("largest deviation of the recomputed mean sojourn : ",
          format(dev_m, digits = 3), " years")
  stopifnot(dev_m < 1e-8)

  fit    <- fit_weibull_grid(tab)
  fit_mm <- fit_weibull_matched_grid(tab, m_i)
  sc     <- build_scenarios(m_i, fit, fit_mm, tab)
  cr     <- crossing_study(cal, sc$scen)
  
  ## The empirical law carries the baseline mean on the integer ages by Corollary 6.1.8,
  ## and the subsequent remark records the deviation at real ages in-between.
  s_grid   <- seq(min(sapply(POLICIES, `[[`, "t0")),
                  max(sapply(POLICIES, function(p) p$t0 + p$n)), by = 1 / 24)
  m_emp    <- vapply(s_grid, function(s) empirical_mean(tab, s), numeric(1))
  mean_gap <- max(abs(m_emp - m_i(s_grid)))
  int      <- abs(s_grid - round(s_grid)) < 1e-9
  message("largest gap between the empirical and the baseline mean : ",
          format(mean_gap, digits = 3), " years")
  stopifnot(max(abs(m_emp[int] - m_i(s_grid[int]))) < 1e-6)
  tail_mass <- tabulated_survival(tab, WINDOW_AGE, Y_ULT)[Y_ULT + 1L]
  emp       <- cr[cr$scenario == "empirical", ]

  ## Preparation of the crossings analysis. The misses split into the three
  ## groups reported in Section 6.1.3, namely the remaining scenarios with a
  ## short reversed stretch, the tabulated law and the DFR shape CROSS_KAPPA.
  miss <- cr[!cr$constant, ]
  emp_miss <- miss[miss$scenario == "empirical", ]
  dfr  <- miss[miss$scenario == paste0("weibull", CROSS_KAPPA), ]
  near <- miss[!miss$scenario %in% c("empirical", paste0("weibull", CROSS_KAPPA)), ]
  stopifnot(nrow(near) + nrow(emp_miss) + nrow(dfr) == nrow(miss))
  
  ## Under the strict convention a crossing falling on the maturity keeps the
  ## sign constant, so the smallest margin among the constant combinations is
  ## reported over the strictly positive margins only.
  touch <- cr[cr$constant & !is.na(cr$margin) & cr$margin <= 0, ]
  free  <- cr$margin[cr$constant & !is.na(cr$margin) & cr$margin > 0]
  
  message("data-based shape between ", format(min(fit$shape), digits = 3),
          " and ", format(max(fit$shape), digits = 3),
          ", worst discrepancy ", format(max(fit$dstar), digits = 3))
  message("sign constant on the window in ", sum(cr$constant), " of ", nrow(cr),
          " combinations, of which ", nrow(touch), " cross on the maturity")
  message("of the remaining ", nrow(miss), " combinations, ", nrow(near),
          " have a short reversed stretch, ", nrow(emp_miss),
          " belong to the empirical law and ", nrow(dfr),
          " to the shape ", CROSS_KAPPA)

  ## The figure of the payable window is drawn at an age without a crossing
  stopifnot(all(cr$constant[cr$age == WINDOW_AGE]))

  ## -------------------------------------------------------------------------
  ## From here on nothing is computed, only written
  ## -------------------------------------------------------------------------

  write.csv(cr,  file.path(OUTPUT_DIR, "sojourn_crossing.csv"), row.names = FALSE)
  write.csv(fit, file.path(OUTPUT_DIR, "sojourn_data_fit.csv"), row.names = FALSE)
  write.csv(fit_mm, file.path(OUTPUT_DIR, "sojourn_data_fit_matched.csv"), row.names = FALSE)
  
  ## Table 6.2.
  ap <- intersect(seq(20, 65, by = 5), fit$age)
  ii <- match(ap, fit$age)
  jj <- match(ap, fit_mm$age)
  write_tex_tabular(
    data.frame(age = ap, alpha = fit$scale[ii], kappa = fit$shape[ii],
               dstar = fit$dstar[ii], mean = fit$mean[ii],
               kappa_m = fit_mm$shape[jj], dstar_m = fit_mm$dstar[jj]),
    headers = c("$x$", "$\\alpha$", "$\\kappa$", "$D^{\\ast}$",
                "mean", "$\\kappa_m$", "$D^{\\ast}_{m}$"),
    digits  = c(0, 2, 3, 4, 2, 3, 4),
    path    = file.path(TAB_DIR, "tab_data_fit.tex"))
  
  ## Figures 6.2 and 6.3.
  fig1 <- open_tex_figure("fig_sojourn_window.pdf", height = 3.60)
  print(plot_window(cal, sc$scen, sc$labels))
  dev.off()
  fig2 <- open_tex_figure("fig_sojourn_annihilation.pdf", height = 4.60)
  print(plot_annihilation(cal, m_i))
  dev.off()

  kk     <- KAPPA[KAPPA != 1]

  ## Every scalar in Section 6.1.3.
  write_tex_macros(c(
    tex_macro("Horizon",           format(Y_ULT)),
    tex_macro("ThetaScaleLo",      tex_fixed(THETA0$scale[1], 1)),
    tex_macro("ThetaScaleHi",      format(THETA0$scale[2])),
    tex_macro("ThetaShapeLo",      tex_fixed(THETA0$shape[1], 1)),
    tex_macro("ThetaShapeHi",      format(THETA0$shape[2])),
    tex_macro("WindowAge",         format(WINDOW_AGE)),
    tex_macro("MeanMatchGap",      tex_ceil(mean_gap, 3)),
    tex_macro("TailMass",          tex_fixed(tail_mass, 4)),
    tex_macro("EmpiricalRuns",     format(nrow(emp))),
    tex_macro("EmpiricalMisses",   format(nrow(emp_miss))),
    tex_macro("EmpiricalAgeLo",    format(min(emp_miss$age))),
    tex_macro("EmpiricalAgeHi",    format(max(emp_miss$age))),
    tex_macro("EmpiricalMissWorst",tex_fixed(max(abs(emp_miss$margin)), 2)),
    tex_macro("EmpiricalCrossLo",  tex_fixed(min(emp_miss$first), 2)),
    tex_macro("EmpiricalCrossHi",  tex_fixed(max(emp_miss$first), 2)),
    tex_macro("CrossRatioMin",     tex_fixed(min(gamma(1 + 1 / kk)^(kk / (1 - kk))), 4)),
    tex_macro("WindowRuns",        format(nrow(cr))),
    tex_macro("WindowHits",        format(sum(cr$constant))),
    tex_macro("WindowMisses",      format(nrow(miss))),
    tex_macro("WindowMargin",      tex_fixed(min(free), 2)),
    tex_macro("WindowTouchAge",    paste(touch$age, collapse = ", ")),
    tex_macro("NearMisses",        format(nrow(near))),
    tex_macro("WindowMissWorst",   tex_fixed(max(abs(near$margin)), 2)),
    tex_macro("DFRMisses",         format(nrow(dfr))),
    tex_macro("DFRShape",          format(CROSS_KAPPA)),
    tex_macro("DFRAgeLo",          format(min(dfr$age))),
    tex_macro("DFRAgeHi",          format(max(dfr$age))),
    tex_macro("DFRMissWorst",      tex_fixed(max(abs(dfr$margin)), 2)),
    tex_macro("DFRMissShare",      tex_fixed(100 * max(abs(dfr$margin / dfr$window)), 1))
  ), file.path(TAB_DIR, "sojourn_numbers.tex"))

  message("sojourn comparison finished!")
}

main()
