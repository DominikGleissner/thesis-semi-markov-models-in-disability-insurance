###############################################################################
##
##  06_02_equivalence_premium.R
##
##  Section 6.2. The equivalence premium of every scenario under
##  the self-consistent pricing E1, together with the two sides of the bias
##  decomposition of (6.24) for the reference policies of Definition 6.1.1.
##
##  Input
##    INPUT_DIR   dav2021i_invalid_male_2ndorder.csv
##    OUTPUT_DIR  baseline_calibration_grid.csv,
##                sojourn_data_fit.csv,
##                sojourn_data_fit_matched.csv
##
##  Output
##    FIG_DIR     fig_premium_bias.pdf
##    TAB_DIR     tab_premium_scenarios.tex, 
##                tab_premium_deferred.tex,
##                premium_numbers.tex
##    OUTPUT_DIR  premium_experiments.csv
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

SRC        <- "06_02_equivalence_premium.R"
TEX_PREFIX <- "Prem"
TEX_SCI    <- FALSE

## Reference policies of Definition 6.1.1.
POLICIES <- list(
  early = list(t0 = 30, n = 35, beta = 1, f = 0, delta = 0.03, h = 1/24),
  late  = list(t0 = 45, n = 20, beta = 1, f = 0, delta = 0.03, h = 1/24))
F_GRID <- c(0, 1/4, 1/2, 1)

MEAN_AGE <- 45               # age at which the body table reports m_i


## ---------------------------------------------------------------------------
## 1. The experiments
## ---------------------------------------------------------------------------

## Premium, benefit value and premium annuity of every scenario, for every
## deferred period and both entry ages, each reported against the baseline of
## the same policy. The sojourn matrix does not depend on f, so it is built
## once per scenario and policy.
run_experiments <- function(cal, scen) {
  ## Collect all values in one data.frame.
  out <- NULL
  for (pn in names(POLICIES)) {
    pol0 <- POLICIES[[pn]]
    mods <- lapply(scen, function(s) build_general(cal, pol0, s))
    for (f in F_GRID) {
      pol  <- modifyList(pol0, list(f = f))
      base <- price_g(mods[["baseline"]], pol)
      for (nm in names(scen)) {
        pr <- price_g(mods[[nm]], pol)
        out <- rbind(out, data.frame(
          policy = pn, t0 = pol$t0, n = pol$n, f = f, scenario = nm,
          premium = pr$premium, epv_b = pr$epv_b, a_bar = pr$a_bar,
          d_pi  = pr$premium / base$premium - 1,
          d_epv = pr$epv_b   / base$epv_b   - 1,
          d_a   = pr$a_bar   / base$a_bar   - 1,
          stringsAsFactors = FALSE))
      }
    }
  }
  out
}


## ---------------------------------------------------------------------------
## 2. Figure 6.4
## ---------------------------------------------------------------------------

## The premium bias along the shape sweep, one curve per deferred period.
plot_bias <- function(ex) {
  df <- ex[ex$scenario %in% paste0("weibull", KAPPA), ]
  df$kappa  <- as.numeric(sub("weibull", "", df$scenario))
  panels    <- sprintf("entry age x = %d, term n = %d",
                       sapply(POLICIES, `[[`, "t0"), sapply(POLICIES, `[[`, "n"))
  df$panel  <- factor(panels[match(df$policy, names(POLICIES))], levels = panels)
  df$defer  <- factor(df$f, levels = F_GRID,
                      labels = c("f = 0", "f = 1/4", "f = 1/2", "f = 1"))

  ggplot(df, aes(x = kappa, y = 100 * d_pi, colour = defer, shape = defer)) +
    geom_hline(yintercept = 0, colour = COL_LINE, linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = colorRampPalette(c(COL_FIT, COL_TABLE))(length(F_GRID))) +
    scale_x_continuous(breaks = KAPPA) +
    labs(x = "Weibull shape kappa", y = "premium bias (%)") +
    theme_thesis()
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
  sc  <- build_scenarios(m_i, fit, fit_mm, tab)
  scen   <- sc$scen
  labels <- sc$labels
  
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

  ## The two engines must yield the same values for the premium.
  pol0  <- POLICIES[["early"]]
  dev_e <- abs(price_g(build_general(cal, pol0), pol0)$premium /
               price(build_model(cal, pol0), pol0)$premium - 1)
  message("relative gap between the general engine and Section 6.1.4 : ",
          format(dev_e, digits = 3))
  stopifnot(dev_e < 1e-10)
  
  ## The mean of the baseline sojourn must match the empirical one.
  dev_g <- abs(empirical_mean(tab, MEAN_AGE) - m_i(MEAN_AGE))
  message("mean gap of the empirical scenario at age ", MEAN_AGE, " : ",
          format(dev_g, digits = 3), " years")
  stopifnot(dev_g < 1e-6)

  ## Run the experiments
  ex <- run_experiments(cal, scen)

  ## The main baseline policy of Definition 6.1.1.
  b0 <- ex[ex$policy == "early" & ex$f == 0, ]
  b0 <- b0[match(names(scen), b0$scenario), ]
  message("baseline premium : ", format(b0$premium[1], digits = 6))
  for (i in seq_len(nrow(b0)))
    message("  ", formatC(b0$scenario[i], width = -12), " ",
            formatC(100 * b0$d_pi[i], format = "f", digits = 2), " %")

  ## -------------------------------------------------------------------------
  ## From here on nothing is computed, only written
  ## -------------------------------------------------------------------------

  write.csv(ex, file.path(OUTPUT_DIR, "premium_experiments.csv"), row.names = FALSE)

  ## Figure 6.4.
  fig1 <- open_tex_figure("fig_premium_bias.pdf", height = 2.80)
  print(plot_bias(ex))
  dev.off()

  ## Table 6.4.
  kk  <- KAPPA[KAPPA != 1]
  rat <- rep("-", length(scen))
  rat[match(paste0("weibull", kk), names(scen))] <-
    paste0("$", tex_fixed(gamma(1 + 1 / kk)^(kk / (1 - kk)), 4), "$")
  mean_col <- ifelse(names(scen) == "fitted",
                     fit$mean[match(MEAN_AGE, fit$age)], m_i(MEAN_AGE))
  
  ## The three deviations carry the per cent sign in the table itself.
  pct_cell <- function(v) paste0("$", tex_fixed(100 * v, 2), "\\%$")
  
  ## Write Table 6.4.
  write_tex_tabular(
    data.frame(scenario = unlist(labels, use.names = FALSE), mean = mean_col,
               ratio = rat, premium = b0$premium, d_pi = pct_cell(b0$d_pi),
               d_epv = pct_cell(b0$d_epv), d_a = pct_cell(b0$d_a),
               stringsAsFactors = FALSE),
    headers = c("scenario", sprintf("$m_{\\mathtt{i}}(%d)$", MEAN_AGE),
                "$z^{\\ast} / m$", "$\\Pi^{h}$", "$\\Delta_{\\Pi}$",
                "$\\Delta_{\\mathrm{EPV}}$", "$\\Delta_{\\bar{a}}$"),
    digits  = c(NA, 2, NA, 5, NA, NA, NA),
    path    = file.path(TAB_DIR, "tab_premium_scenarios.tex"),
    align   = "lrrrrrr")
  
  ## For Table 6.5, the row labels are grouped, which keeps the first column narrow.
  ## The baseline is left out, since its deviation is zero by definition.
  disp <- c(lognormal = "log-normal",
            empirical = "\\quad empirical",
            fittedmm  = "\\quad mean-matched",
            fitted    = "\\quad free")
  for (k in KAPPA)
    disp[paste0("weibull", k)] <- sprintf("\\quad kappa = %s", format(k))
  
  ## Define the order of the scenarios in the table.
  ord <- c("@Weibull", paste0("weibull", KAPPA), "lognormal",
           "@data-based", "empirical", "fittedmm", "fitted")
  
  ## Data.frame for the wider format empty rows for the category lines.
  wide <- NULL
  for (nm in ord) {
    if (substr(nm, 1, 1) == "@") {
      row <- c(substring(nm, 2), rep("", length(POLICIES) * length(F_GRID)))
    } else {
      row <- disp[nm]
      for (pn in names(POLICIES)) for (f in F_GRID) {
        s   <- ex[ex$policy == pn & ex$f == f, ]
        row <- c(row, pct_cell(s$d_pi[match(nm, s$scenario)]))
      }
    }
    wide <- rbind(wide, row)
  }
  wide <- as.data.frame(wide, stringsAsFactors = FALSE)

  ## Determine the top line of the table with the according spaces in-between.
  span <- length(F_GRID)
  grp  <- character(0)
  cmid <- character(0)
  for (k in seq_along(POLICIES)) {
    p    <- POLICIES[[k]]
    ## The actual content of the line including entry age and term.
    grp  <- c(grp, sprintf("\\multicolumn{%d}{c}{entry age $x = %d$, term $n = %d$}",
                           span, p$t0, p$n))
    ## The small lines below the above line and the header line.
    j    <- 1 + (k - 1) * span
    cmid <- c(cmid, sprintf("\\cmidrule(lr){%d-%d}", j + 1, j + span))
  }

  ## Write Table 6.5.
  write_tex_tabular(
    wide,
    headers = c("scenario", rep(c("$f = 0$", "$f = \\frac{1}{4}$",
                                  "$f = \\frac{1}{2}$", "$f = 1$"), 2)),
    digits  = rep(NA, ncol(wide)),
    path    = file.path(TAB_DIR, "tab_premium_deferred.tex"),
    align   = "lrrrrrrrr",
    above   = c(paste0("     & ", paste(grp, collapse = " & "), " \\\\"),
                paste0("    ", paste(cmid, collapse = " "))))
  
  ## Every scalar in Section 6.2.
  pct    <- function(v) paste0(tex_fixed(100 * v, 2), "\\%")
  pick   <- function(nm) b0$d_pi[match(nm, b0$scenario)]
  m_base <- cal$m_hat[cal$age %in% FIT_AGES]
  spread <- 0
  for (pn in names(POLICIES)) for (nm in names(scen))
    spread <- max(spread, diff(range(ex$d_pi[ex$policy == pn & ex$scenario == nm])))

  write_tex_macros(c(
    tex_macro("MeanAge",           format(MEAN_AGE)),
    tex_macro("BaselinePremium",   tex_fixed(b0$premium[1], 5)),
    tex_macro("DeltaWeibullHalf",  pct(pick("weibull0.5"))),
    tex_macro("DeltaWeibullOne",   pct(pick("weibull1"))),
    tex_macro("DeltaWeibullTwo",   pct(pick("weibull2"))),
    tex_macro("DeltaLogNormal",    pct(pick("lognormal"))),
    tex_macro("DeltaData",         pct(pick("fitted"))),
    tex_macro("DeltaEmpirical",        pct(pick("empirical"))),
    tex_macro("DeltaEmpiricalAbs",     pct(abs(pick("empirical")))),
    tex_macro("DeltaDataMatched",      pct(pick("fittedmm"))),
    tex_macro("DataKappaRangeMatched", sprintf("[%s,\\ %s]",
                                               tex_floor(min(fit_mm$shape), 2),
                                               tex_ceil(max(fit_mm$shape), 2))),
    tex_macro("DataFitMaxMatched",     tex_fixed(max(fit_mm$dstar), 4)),
    tex_macro("DeferralSpread",    tex_fixed(100 * spread, 2)),
    tex_macro("AnnuityMax",        pct(ex$d_a[which.max(abs(ex$d_a))])),
    tex_macro("DataKappaRange",    sprintf("[%s,\\ %s]", tex_floor(min(fit$shape), 2),
                                           tex_ceil(max(fit$shape), 2))),
    tex_macro("DataFitMax",        tex_fixed(max(fit$dstar), 4)),
    tex_macro("DataMeanRange",     sprintf("[%s,\\ %s]", tex_floor(min(fit$mean), 1),
                                           tex_ceil(max(fit$mean), 1))),
    tex_macro("BaselineMeanRange", sprintf("[%s,\\ %s]", tex_floor(min(m_base), 1),
                                           tex_ceil(max(m_base), 1)))
  ), file.path(TAB_DIR, "premium_numbers.tex"))

  message("premium section finished!")
}

main()
