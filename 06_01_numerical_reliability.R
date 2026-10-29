###############################################################################
##
##  06_01_numerical_reliability.R
##
##  Section 6.1.4. The three diagnostics, all run under the Markov baseline
##  S^(0) of Section 6.1.2 and the reference policy of Definition 6.1.1.
##
##    (1) cross-validation of the deterministic and the Monte Carlo engine at a
##        set of lattice points, as well as the premium gap between S^(0) 
##        and the Markov reference model M^(0)
##    (2) premium, benefit value and one invalid reserve along the step sizes
##    (3) standard error of the mean simulated loss along the path numbers
##
##  Input
##    OUTPUT_DIR  baseline_calibration_grid.csv
##
##  Output
##    TAB_DIR     tab_step_size.tex 
##                reliability_numbers.tex 
##    OUTPUT_DIR  reliability_cross_validation.csv, 
##                reliability_simulation_size.csv,
##                reliability_simulation_size.csv
##
###############################################################################

source("06_00_common.R")
source("06_00_engine.R")


## ---------------------------------------------------------------------------
## 0. Configuration
## ---------------------------------------------------------------------------

SRC        <- "06_01_numerical_reliability.R"
TEX_PREFIX <- "Reli"
TEX_SCI    <- FALSE

## Reference policies in Definition 6.1.1.
POLICY   <- list(t0 = 30, n = 35, beta = 1, f = 0, delta = 0.03, h = 1/24)
H_GRID   <- c(1/4, 1/12, 1/24, 1/48, 1/360)
M_GRID   <- c(1e4, 1e5, 1e6)
STEP_ROW <- length(H_GRID) - 1

## Cross-validation: the lattice points, the seeds and the level 1 - gamma.
VALID_AGES <- c(35, 45, 55)  # attained ages t
VALID_DUR  <- c(0, 2)        # elapsed durations u, in years
VALID_M    <- 1e6            # paths per run
SEED_BASE  <- 19690902
SEEDS      <- 1:5
GAMMA      <- 0.05

RESERVE_AGE <- 45            # age of the invalid reserve in the step-size table


## ---------------------------------------------------------------------------
## 1. The three diagnostics
## ---------------------------------------------------------------------------

## (1) Both engines at the same lattice points, with several seeds each. What
## is reported is the coverage of the confidence interval (5.29), since single misses
## are expected at the level 1 - gamma, and the largest miss measured in half-widths.
cross_validation <- function(mod, pol, pr) {
  z   <- qnorm(1 - GAMMA / 2)
  
  ## Create the set of lattice points with all combinations of the values.
  pts <- expand.grid(state = c("a", "i"), u = VALID_DUR, t = VALID_AGES,
                     stringsAsFactors = FALSE)
  ## Collect all values in one data.frame.
  out <- NULL
  for (i in seq_len(nrow(pts))) {
    l_t <- round((pts$t[i] - pol$t0) / pol$h)
    l_u <- round(pts$u[i] / pol$h)
    ## Computation of the deterministic reserve.
    det <- reserve(mod, pol, pr, pts$state[i], l_t, l_u, pr$premium)
    for (r in SEEDS) {
      ## Vary the seeds for each point and for the different seeds.
      ## This reduces the probability to pick a few "(un)lucky" seeds.
      set.seed(SEED_BASE + 10 * i + r)
      ## Computation of the Monte-Carlo reserve.
      L  <- simulate_loss(mod, pol, pts$state[i], l_t, l_u, VALID_M, pr$premium)
      ## Determine the half-width of the confidence interval (5.29)
      hw <- z * sd(L) / sqrt(VALID_M)
      out <- rbind(out, data.frame(
        state = pts$state[i], t = pts$t[i], u = pts$u[i], seed = r,
        deterministic = det, mc_mean = mean(L), half_width = hw,
        miss = abs(mean(L) - det) / hw, stringsAsFactors = FALSE))
    }
  }
  ## Define another boolean column to see whether a value lies inside the confidence
  ## interval, which is the case if and only if it differs by at most one half-width.
  out$covered <- out$miss <= 1
  
  out
}

## (2) The premium, the benefit value and one invalid reserve along the step
## sizes, each with the relative change against the next finer step.
step_size_study <- function(cal, pol) {
  ## Collect all values in one data.frame.
  out <- NULL
  for (h in H_GRID) {
    p  <- modifyList(pol, list(h = h))  # update the step size of the policy
    m  <- build_model(cal, p)           # build the model for the new step size 
    pr <- price(m, p)                   # compute the corresponding price structure
    out <- rbind(out, data.frame(
      h = h, premium = pr$premium, epv_b = pr$epv_b,
      reserve_i = reserve(m, p, pr, "i", round((RESERVE_AGE - p$t0) / h), 0,
                          pr$premium)))
  }
  ## Function to calculate the neighboring relative differences in a vector
  ## rel = v_j / v_{j+1} - 1 = (v_j - v_{j+1}) / v_{j+1}
  rel <- function(v) c(v[-length(v)] / v[-1] - 1, NA)
  ## Add the relative differences to the output as well for each step size
  cbind(out, rel_premium = rel(out$premium), rel_epv = rel(out$epv_b),
        rel_reserve = rel(out$reserve_i))
}

## (3) The standard error of the mean of the policy loss L^h_{(a, 0)}(t_0),
## priced with the baseline premium, along the simulation sizes.
simulation_size_study <- function(mod, pol, pr) {
  ## Collect all values in one data.frame.
  out <- NULL
  for (M in M_GRID) {
    ## Reset the (same!) seed between all simulations
    set.seed(SEED_BASE)
    L <- simulate_loss(mod, pol, "a", 0, 0, M, pr$premium)
    out <- rbind(out, data.frame(M = M, mean = mean(L), sd = sd(L),
                                 se = sd(L) / sqrt(M)))
  }
  out
}


## ---------------------------------------------------------------------------
## 2. Driver
## ---------------------------------------------------------------------------

main <- function() {
  cal <- read.csv(file.path(OUTPUT_DIR, "baseline_calibration_grid.csv"),
                  stringsAsFactors = FALSE)
  need <- c("age", "sigma", "mu", "q_i", "p_ia_hat")
  if (!all(need %in% names(cal)))
    stop("the calibrated grid must have the columns ", paste(need, collapse = ", "))
  cal <- cal[order(cal$age), ]
  if (min(cal$age) > POLICY$t0 || max(cal$age) < POLICY$t0 + POLICY$n)
    stop("the calibrated grid does not cover the contract period")

  ## Determine the models with standard reference policy and the prices.
  mod <- build_model(cal, POLICY)                      # baseline S^(0)
  ref <- build_model(cal, POLICY, split = "attained")  # Markov reference M^(0)
  pr  <- price(mod, POLICY)                            # baseline S^(0)
  prm <- price(ref, POLICY)                            # Markov reference M^(0)
  ## Premium difference between Markov baseline S^(0) and Markov reference M^(0).
  gap <- prm$premium / pr$premium - 1

  message("premium  : ", format(pr$premium, digits = 6))
  message("EPV_B    : ", format(pr$epv_b, digits = 6))
  message("premium gap S^(0) against M^(0) : ", format(gap, digits = 3))

  ## Check whether the equivalence principle holds at issue.
  stopifnot(abs(reserve(mod, POLICY, pr, "a", 0, 0, pr$premium)) < 1e-10)

  ## Run the three diagnostics.
  cv <- cross_validation(mod, POLICY, pr)
  ss <- step_size_study(cal, POLICY)
  ms <- simulation_size_study(mod, POLICY, pr)
  
  ## The working grid against the finest grid of H_GRID,
  ## in the sense (X(h_work) - X(h_min)) / X(h_min), maximised over the three
  ## quantities. In thesis "the quantities vary by less than ...".
  i_w <- which.min(abs(ss$h - POLICY$h))   # row of the working grid
  i_f <- which.min(ss$h)                   # row of the finest grid
  wf  <- max(abs(c(ss$premium[i_w]   / ss$premium[i_f]   - 1,
                   ss$epv_b[i_w]     / ss$epv_b[i_f]     - 1,
                   ss$reserve_i[i_w] / ss$reserve_i[i_f] - 1)))

  message("coverage : ", sum(cv$covered), " of ", nrow(cv),
          " runs, largest miss ", format(max(cv$miss), digits = 3), " half-widths")
  message("step size: largest change between the two finest steps ",
          format(max(abs(c(ss$rel_premium[STEP_ROW], ss$rel_epv[STEP_ROW],
                           ss$rel_reserve[STEP_ROW]))),
                 digits = 3))
  message("step size: largest change working against finest step ",
          format(wf, digits = 3))
  message("path number: standard error ", format(ms$se[1], digits = 3), " -> ",
          format(ms$se[nrow(ms)], digits = 3))

  ## -------------------------------------------------------------------------
  ## From here on nothing is computed, only written
  ## -------------------------------------------------------------------------

  write.csv(cv, file.path(OUTPUT_DIR, "reliability_cross_validation.csv"),
            row.names = FALSE)
  write.csv(ss, file.path(OUTPUT_DIR, "reliability_step_size.csv"),
            row.names = FALSE)
  write.csv(ms, file.path(OUTPUT_DIR, "reliability_simulation_size.csv"),
            row.names = FALSE)

  ## Table 6.3.
  rel_col <- function(v) ifelse(is.na(v), "",
                               paste0("$", sapply(v, tex_scalar), "$"))
  write_tex_tabular(
    data.frame(h       = sprintf("$1/%d$", round(1 / ss$h)),
               premium = ss$premium,   rel_p = rel_col(ss$rel_premium),
               epv     = ss$epv_b,     rel_e = rel_col(ss$rel_epv),
               res     = ss$reserve_i, rel_r = rel_col(ss$rel_reserve),
               stringsAsFactors = FALSE),
    headers = c("$h$", "$\\Pi^{h}$", "rel.", "$\\mathrm{EPV}^{h}_{B}$", "rel.",
                sprintf("$V^{+, h}_{(\\mathtt{i}, 0)}(%d)$", RESERVE_AGE), "rel."),
    digits  = c(NA, 6, NA, 5, NA, 5, NA),
    path    = file.path(TAB_DIR, "tab_step_size.tex"))

  ## Every scalar in Section 6.1.4
  write_tex_macros(c(
    tex_macro("Points",       format(nrow(cv) / length(SEEDS))),
    tex_macro("Seeds",        format(length(SEEDS))),
    tex_macro("Runs",         format(nrow(cv))),
    tex_macro("Level",        tex_fixed(1 - GAMMA, 2)),
    tex_macro("ValidPaths",   sprintf("10^{%d}", round(log10(VALID_M)))),
    tex_macro("Ages",         paste(VALID_AGES, collapse = ",\\ ")),
    tex_macro("Durations",    paste(VALID_DUR, collapse = ",\\ ")),
    tex_macro("CoverageHits", format(sum(cv$covered))),
    tex_macro("WorstMiss",    tex_fixed(max(cv$miss), 2)),
    tex_macro("Premium",      tex_fixed(pr$premium, 5)),
    tex_macro("EPVB",         tex_fixed(pr$epv_b, 5)),
    tex_macro("PremiumGap",   tex_scalar(gap)),
    tex_macro("PremiumGapPercent", paste0(tex_ceil(100 * gap, 2), "\\%")),
    tex_macro("ReserveAge",   format(RESERVE_AGE)),
    tex_macro("StepMax",      tex_scalar(max(abs(c(ss$rel_premium[STEP_ROW],
                                                   ss$rel_epv[STEP_ROW],
                                                   ss$rel_reserve[STEP_ROW]))))),
    tex_macro("StepMaxPercent", paste0(tex_ceil(100 * max(abs(c(ss$rel_premium[STEP_ROW],
                                                                ss$rel_epv[STEP_ROW],
                                                                ss$rel_reserve[STEP_ROW]))), 1), "\\%")),
    tex_macro("StepMaxWorkingToFinest",        tex_scalar(wf)),
    tex_macro("StepMaxWorkingToFinestPercent", paste0(tex_ceil(100 * wf, 1), "\\%")),
    tex_macro("SeCoarse",     tex_scalar(ms$se[1])),
    tex_macro("SeFine",       tex_scalar(ms$se[nrow(ms)])),
    tex_macro("SeFactor",     tex_fixed(ms$se[1] / ms$se[nrow(ms)], 2))
  ), file.path(TAB_DIR, "reliability_numbers.tex"))

  message("reliability finished!")
}

main()
