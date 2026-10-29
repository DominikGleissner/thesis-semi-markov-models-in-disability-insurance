###############################################################################
##
##  06_01_baseline_calibration.R
##
##  Section 6.1.2. Calibration of the Markov baseline S^(0) from the
##  DAV 2021 I second-order tables for male policyholders. The steps themselves
##  live in 06_00_calibration.R, this file runs them in the order of the thesis
##  and writes the figure, the tables and the calibrated grid.
##
##  Input
##    INPUT_DIR   dav2021i_active_male_2ndorder.csv
##                dav2021i_invalid_male_2ndorder.csv
##
##  Output
##    FIG_DIR     fig_baseline_sojourn_fit.pdf
##    TAB_DIR     tab_baseline_calibration.tex
##                calibration_numbers.tex
##    OUTPUT_DIR  baseline_calibration_grid.csv, 
##                baseline_appendix_table.csv,
##                diagnostic_dependent_inception.csv, 
##                diagnostic_splitting_ratio.csv
##
##  Needs ggplot2 for the figure.
##
###############################################################################

library(ggplot2)

source("06_00_common.R")
source("06_00_calibration.R")


## ---------------------------------------------------------------------------
## 0. Configuration
## ---------------------------------------------------------------------------

SRC        <- "06_01_baseline_calibration.R"
TEX_PREFIX <- "Calib"
TEX_SCI    <- TRUE

APPENDIX_AGES  <- seq(20, 65, by = 5)
FIG_ENTRY_AGES <- c(30, 45, 60)
FIG_Z_MAX      <- 35L        # horizontal range of the figure, in years


## ---------------------------------------------------------------------------
## 1. Figure 6.1
## ---------------------------------------------------------------------------

## The calibrated baseline survival against the tabulated one, one panel per
## entry age. The curve is drawn on a fine grid, the tabulated values at
## the integer durations where the table lives.
plot_sojourn_fit <- function(tab, cal, entry_ages = FIG_ENTRY_AGES,
                             z_max = FIG_Z_MAX) {
  panels <- paste("entry age x =", entry_ages)
  z_fine <- seq(0, z_max, by = 0.05)

  fitted <- do.call(rbind, lapply(seq_along(entry_ages), function(i)
    data.frame(panel = factor(panels[i], levels = panels),
               z     = z_fine,
               surv  = baseline_survival(cal, entry_ages[i], z_fine))))

  tabled <- do.call(rbind, lapply(seq_along(entry_ages), function(i)
    data.frame(panel = factor(panels[i], levels = panels),
               z     = 0:z_max,
               surv  = tabulated_survival(tab, entry_ages[i], z_max))))

  ggplot(fitted, aes(x = z, y = surv)) +
    geom_line(aes(colour = "Markov baseline"), linewidth = 0.7) +
    geom_point(data = tabled, aes(colour = "DAV 2021 I table"), size = 0.9) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = c("Markov baseline"   = COL_FIT,
                                   "DAV 2021 I table"  = COL_TABLE),
                        breaks = c("Markov baseline", "DAV 2021 I table")) +
    guides(colour = guide_legend(
      override.aes = list(linetype = c(1, 0), shape = c(NA, 16)))) +
    scale_x_continuous(breaks = seq(0, z_max, by = 5)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    labs(x = "elapsed duration z (years)", y = "sojourn survival") +
    theme_thesis()
}


## ---------------------------------------------------------------------------
## 2. Driver
## ---------------------------------------------------------------------------

main <- function() {
  stopifnot(MEAN_RULE %in% c("trapezoidal", "constant_force"))

  active_raw <- read_active_table(
    file.path(INPUT_DIR, "dav2021i_active_male_2ndorder.csv"))
  tab <- read_invalid_table(
    file.path(INPUT_DIR, "dav2021i_invalid_male_2ndorder.csv"))

  message("active table  : ages ", min(active_raw$age), "-", max(active_raw$age))
  message("invalid table : ages ", tab$x_min, "-", tab$omega,
          ", disability years 1-", Y_ULT, " (", Y_ULT, " = ultimate column 10+)")
  message("mean rule     : ", MEAN_RULE)

  ## Active state.
  act <- invert_active_intensities(active_raw)

  ## Invalid state.
  cal <- calibrate_invalid_intensity(tab)
  message("condition (6.15) holds on the whole age range")
  message("largest gap between the trapezoidal and the constant-force mean : ",
          format(max(abs(cal$m_hat_trapezoidal - cal$m_hat_constant_force)),
                 digits = 3), " years")

  ## Jump probabilities.
  pa <- active_jump_probabilities(act)
  pv <- invalid_jump_probabilities(tab, tab$ages)

  ## Checks reported in Section 6.1.2.
  dev_a <- max(abs(pa$p_ai + pa$p_ad - 1))
  dev_i <- max(abs(pv$p_ia_hat + pv$p_id_hat - 1))
  message("largest deviation of the active  row sum from one : ", format(dev_a, digits = 3))
  message("largest deviation of the invalid row sum from one : ", format(dev_i, digits = 3))
  stopifnot(dev_a < 1e-10, dev_i < 1e-10)
  stopifnot(all(cal$q_i > 0), all(act$sigma > 0), all(act$mu > 0))

  ## Mean preservation: run (6.18) forward again from the
  ## calibrated q_x and compare with the tabulated target.
  n <- nrow(cal)
  m_back <- numeric(n)
  m_back[n] <- 1 / cal$q_i[n]
  if (n >= 2)
    for (k in (n - 1L):1L)
      m_back[k] <- g_recursion(cal$q_i[k], m_back[k + 1L])
  dev_m <- max(abs(m_back - cal$m_hat))
  message("largest deviation between m_i and m_hat_i : ", format(dev_m, digits = 3))
  stopifnot(dev_m < 1e-7)

  dep <- dependent_inception_check(act)
  message("DAV correction factor, largest relative deviation over the ages : ",
          format(max(dep$rel_dev), digits = 3))
  spr <- splitting_ratio_spread(tab, tab$ages)
  message("splitting ratio, largest spread over the disability years : ",
          format(max(spr$w_range), digits = 3))
  message("mass carried by the constant continuation, active  state : ",
          format(max(pa$a_tail_mass), digits = 3))
  message("mass carried by the constant continuation, invalid state : ",
          format(max(pv$i_tail_mass), digits = 3))
  
  # Near-constancy of p_ai in Table 6.1. The smallest mass the constant
  # continuation carries over the tabulated ages, and the split it is frozen at.
  tail_a_min  <- min(pa$a_tail_mass)
  split_omega <- with(act[which.max(act$age), ], sigma / (sigma + mu))
  message("active state, smallest mass carried by the continuation : ",
          format(tail_a_min, digits = 3), ", split at the last age : ",
          format(split_omega, digits = 3))

  grid <- merge(merge(act, pa, by = "age"),
                merge(merge(cal, pv, by = "age"),
                      alternative_intensities(tab, tab$ages), by = "age"),
                by = "age")

  ## -------------------------------------------------------------------------
  ## From here on nothing is computed, only written
  ## -------------------------------------------------------------------------

  write.csv(grid, file.path(OUTPUT_DIR, "baseline_calibration_grid.csv"),
            row.names = FALSE)
  write.csv(dep, file.path(OUTPUT_DIR, "diagnostic_dependent_inception.csv"),
            row.names = FALSE)
  write.csv(spr, file.path(OUTPUT_DIR, "diagnostic_splitting_ratio.csv"),
            row.names = FALSE)

  ## Table 6.1.
  sel <- intersect(APPENDIX_AGES, grid$age)
  if (length(sel) < length(APPENDIX_AGES))
    warning("the appendix ages ",
            paste(setdiff(APPENDIX_AGES, grid$age), collapse = ", "),
            " lie outside the calibrated range and are dropped")
  appx <- grid[match(sel, grid$age),
               c("age", "sigma", "mu", "q_i", "m_hat", "p_ai", "p_ia_hat")]
  write.csv(appx, file.path(OUTPUT_DIR, "baseline_appendix_table.csv"),
            row.names = FALSE)

  write_tex_tabular(
    appx,
    headers = c("$x$", "$\\sigma_x$", "$\\mu_x$", "$q_x$",
                "$\\widehat{m}_{\\mathtt{i}}(x)$",
                "$\\widehat{p}_{\\mathtt{a}\\mathtt{i}}(x)$",
                "$\\widehat{p}_{\\mathtt{i}\\mathtt{a}}(x)$"),
    digits = c(0, 5, 5, 5, 4, 4, 4),
    path   = file.path(TAB_DIR, "tab_baseline_calibration.tex"))
  
  ## Figure 6.1.
  open_tex_figure("fig_baseline_sojourn_fit.pdf")
  print(plot_sojourn_fit(tab, cal))
  dev.off()

  ## Every scalar in Section 6.1.2.
  write_tex_macros(c(
    tex_macro("MeanRule",         MEAN_RULE),
    tex_macro("AgeMin",           format(min(grid$age))),
    tex_macro("AgeMax",           format(max(grid$age))),
    tex_macro("DisabilityYears",  format(Y_ULT)),
    tex_macro("MeanRuleGap",      tex_scalar(max(abs(cal$m_hat_trapezoidal -
                                                     cal$m_hat_constant_force)))),
    tex_macro("RowSumDevActive",  tex_scalar(dev_a)),
    tex_macro("RowSumDevInvalid", tex_scalar(dev_i)),
    tex_macro("MeanDev",          tex_scalar(dev_m)),
    tex_macro("DavRelDev",        tex_scalar(max(dep$rel_dev))),
    tex_macro("SplitSpread",      tex_scalar(max(spr$w_range))),
    tex_macro("TailMassActive",   tex_scalar(max(pa$a_tail_mass))),
    tex_macro("TailMassActiveMinPercent", paste0(tex_floor(100 * tail_a_min, 0), "\\%")),
    tex_macro("ActiveSplitOmega",         tex_fixed(split_omega, 2)),
    tex_macro("TailMassInvalid",  tex_scalar(max(pv$i_tail_mass))),
    tex_macro("RVersion",         R.version.string),
    tex_macro("RStudioVersion",   as.character(RStudio.Version()$version)),
    tex_macro("RunDate",          format(Sys.Date(), "%d %B %Y"))
  ), file.path(TAB_DIR, "calibration_numbers.tex"))

  message("calibration finished!")
}

main()
