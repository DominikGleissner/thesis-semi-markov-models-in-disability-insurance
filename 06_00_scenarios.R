###############################################################################
##
##  06_00_scenarios.R
##
##  The scenarios of Chapter 6: the shape sweep across Weibull scenarios,
##  the log-normal, the Markov baseline, the empirical and the two fitted
##  Weibull scenarios. See Definitions 5.5.8, 5.5.9, 6.1.7.
##
###############################################################################


## ---------------------------------------------------------------------------
## 1. Basic parameters and constants
## ---------------------------------------------------------------------------

KAPPA    <- c(0.5, 0.75, 1, 1.5, 2)      # shape sweep of Definition 5.5.8
SIGMA_LN <- sqrt(log(2))                 # log-normal scenario with cv = 1
## Search region for Weibull parameters for fitting, see Definition 5.5.9
THETA0   <- list(scale = c(0.5, 250), shape = c(0.2, 5))   
FIT_AGES <- 20:70                        # entry ages of the data-based fit


## ---------------------------------------------------------------------------
## 2. The baseline expected sojourn
## ---------------------------------------------------------------------------

## Baseline expected sojourn m_i(s) of Definition 5.5.1. 
## On [x, x+1) the exit intensity is the constant q_x, so splitting
## the integral at the integer ages gives the recursion of Proposition 6.1.6.
mean_sojourn <- function(cal) {
  q <- cal$q_i
  K <- length(q)
  m <- numeric(K)
  ## constant continuation beyond omega
  m[K] <- 1 / q[K]
  for (k in (K - 1):1) m[k] <- (1 - exp(-q[k])) / q[k] + exp(-q[k]) * m[k+1]

  ## Closure: keeps the integer-age vector m from the recursion and 
  ## extends it to any real age via one partial-year step.
  ## The recursion thus runs once per calibration, 
  ## so later evaluations (e.g. inside numerical integrals) are cheap.
  function(age) {
    ## Either return age or one boundary of the ages
    x  <- pmin(pmax(age, min(cal$age)), max(cal$age))
    ## Match this age to the next-smaller integer-valued age
    i  <- match(as.integer(floor(x)), cal$age)
    ## Part of the year still to come to the next full age
    th <- cal$age[i] + 1 - x
    ## Returns m_i(s) even for real-valued s
    (1 - exp(-q[i] * th)) / q[i] + exp(-q[i] * th) * m[pmin(i + 1, K)]
  }
}


## ---------------------------------------------------------------------------
## 3. The single-channel scenarios
## ---------------------------------------------------------------------------

## Mean-matched Weibull
scenario_weibull <- function(m_i, kappa)
  function(x, z) pweibull(z, shape = kappa, scale = m_i(x) / gamma(1 + 1 / kappa),
                          lower.tail = FALSE)

## Mean-matched log-normal
scenario_lognormal <- function(m_i, sigma = SIGMA_LN)
  function(x, z) plnorm(z, meanlog = log(m_i(x)) - sigma^2 / 2, sdlog = sigma,
                        lower.tail = FALSE)

## The empirical scenario S_hat with the tabulated law itself
scenario_empirical <- function(tab)
  function(x, z) tabulated_survival_at(tab, x, z)

## Fitted free Weibull
scenario_fitted <- function(fit)
  function(x, z) pweibull(z, shape = fit$shape[which.min(abs(fit$age - floor(x)))],
                          scale = fit$scale[which.min(abs(fit$age - floor(x)))],
                          lower.tail = FALSE)

## Fitted mean-matched Weibull. The shape is the fitted one at the integer
## entry age, the scale is the mean-matched scale at the entry age.
scenario_fitted_matched <- function(fit, m_i)
  function(x, z) {
    k <- fit$shape[which.min(abs(fit$age - floor(x)))]
    pweibull(z, shape = k, scale = m_i(x) / gamma(1 + 1 / k), lower.tail = FALSE)
  }

## Minimum distance estimator in Definition 5.5.9 for the Weibull family,
## with the pair minimising the grid discrepancy (5.44) over the
## compact search region. Several starting values, because the criterion is a
## maximum and therefore not smooth.
## Optimise over p = (log scale, log shape): keeps both parameters positive 
## and gives Nelder-Mead comparable relative step sizes for scale (~10-100) 
## and shape (~1). The bijection exp() maps back inside d().
fit_weibull <- function(tab, x, theta0 = THETA0) {
  y  <- seq_len(Y_ULT)
  Hh <- 1 - tabulated_survival(tab, x, Y_ULT)[-1] # only look at y >= 1
  d  <- function(p) {
    th <- exp(p)
    if (th[1] < theta0$scale[1] || th[1] > theta0$scale[2] ||
        th[2] < theta0$shape[1] || th[2] > theta0$shape[2]) return(1)
    max(abs(Hh - pweibull(y, shape = th[2], scale = th[1])))
  }
  best <- NULL
  for (a0 in c(10, 30, 80)) for (k0 in c(0.6, 1, 1.6)) {
    o <- optim(log(c(a0, k0)), d, method = "Nelder-Mead",
               control = list(reltol = 1e-12, maxit = 2000))
    if (is.null(best) || o$value < best$value) best <- o
  }
  th <- exp(best$par)
  data.frame(age = x, scale = th[1], shape = th[2], dstar = best$value,
             mean = th[1] * gamma(1 + 1 / th[2]))
}

## Fits the Weibull model for all ages and returns a data frame with one row per age
## and columns age, scale, shape, dstar, mean
fit_weibull_grid <- function(tab, ages = FIT_AGES)
  do.call(rbind, lapply(ages, function(x) fit_weibull(tab, x)))

## Similar approach as for free Weibull fit, but instead of optimizing over 
## (shape, scale), do it only over shape. This already fixes the scale by (5.40).
## The criterion is a maximum and therefore not smooth,
## which is why optimize alone is not trusted with the whole region.
## This could yield only a local minimum of d instead of global.
## Therefore, we do the discretisation to find the rough region beforehand.
fit_weibull_matched <- function(tab, m_i, x, theta0 = THETA0) {
  y  <- seq_len(Y_ULT)
  Hh <- 1 - tabulated_survival(tab, x, Y_ULT)[-1]
  m  <- m_i(x)
  d  <- function(lk) {
    k <- exp(lk)
    a <- m / gamma(1 + 1 / k)
    if (a < theta0$scale[1] || a > theta0$scale[2]) return(1)
    max(abs(Hh - pweibull(y, shape = k, scale = a)))
  }
  ## Discretise the search interval to find the the smallest value roughly
  lk <- seq(log(theta0$shape[1]), log(theta0$shape[2]), length.out = 200)
  v  <- vapply(lk, d, numeric(1))
  j  <- which.min(v)
  ## Precise optimization around the rough minimal point
  o  <- optimize(d, c(lk[max(j - 1L, 1L)], lk[min(j + 1L, length(lk))]),
                 tol = 1e-10)
  k  <- exp(o$minimum)
  data.frame(age = x, scale = m / gamma(1 + 1 / k), shape = k,
             dstar = o$objective, mean = m)
}

## Fits the Weibull model for all ages and returns a data frame with one row per age
## and columns age, scale, shape, dstar, mean
fit_weibull_matched_grid <- function(tab, m_i, ages = FIT_AGES)
  do.call(rbind, lapply(ages, function(x) fit_weibull_matched(tab, m_i, x)))


## ---------------------------------------------------------------------------
## 4. The scenario set and its labels
##
##  The order fixed here is the order in which every table and figure reports them.
##  The baseline entry carries no survival function,
##  since the engine builds it from the calibrated kernel itself.
## ---------------------------------------------------------------------------

build_scenarios <- function(m_i, fit, fit_mm, tab) {
  scen <- c(list(baseline = NULL),
            setNames(lapply(KAPPA, function(k) scenario_weibull(m_i, k)),
                     paste0("weibull", KAPPA)),
            list(lognormal = scenario_lognormal(m_i),
                 empirical = scenario_empirical(tab),
                 fittedmm  = scenario_fitted_matched(fit_mm, m_i),
                 fitted    = scenario_fitted(fit)))
  labels <- c(list(baseline = "baseline"),
              setNames(as.list(sprintf("Weibull, kappa = %g", KAPPA)),
                       paste0("weibull", KAPPA)),
              list(lognormal = "log-normal",
                   empirical = "empirical",
                   fittedmm  = "fitted Weibull, mean-matched",
                   fitted    = "fitted Weibull, free"))
  list(scen = scen, labels = labels)
}


## ---------------------------------------------------------------------------
## 5. The appearance of a scenario
##
##  The colour and the line type belong to the scenario and not to the figure,
##  so that a scenario looks the same wherever Chapter 6 draws it.
##  Mostly, the colours are taken from the TUM Corporate Design.
## ---------------------------------------------------------------------------

SCEN_SWEEP <- paste0("weibull", KAPPA)

## The scenarios outside the shape sweep get their colour and their line type
## one by one. The names are the keys of build_scenarios.
SCEN_ALT_COL <- c(lognormal = "#F9BF4E",
                  empirical = "#9FBA36", 
                  fittedmm  = "#C680BB",
                  fitted    = "#9B468D",
                  baseline  = "#D95117")
## Keep the shape sweep solid and the others dashed for an easier distinguishing
SCEN_ALT_LTY <- c(lognormal = "42",        # 42 = dashed
                  empirical = "42",
                  fittedmm  = "42",
                  fitted    = "42",
                  baseline  = "solid")
SCEN_ALT     <- names(SCEN_ALT_COL)

## Assign a colour to each scenario
scenario_colours <- function(nm) {
  blue <- colorRampPalette(c(COL_FIT, COL_TABLE))(length(SCEN_SWEEP))
  col  <- setNames(rep("black", length(nm)), nm)
  col[nm %in% SCEN_SWEEP] <- blue[match(nm[nm %in% SCEN_SWEEP], SCEN_SWEEP)]
  col[nm %in% SCEN_ALT]   <- SCEN_ALT_COL[nm[nm %in% SCEN_ALT]]
  col
}

## Assign a line type to each scenario
scenario_types <- function(nm) {
  lty <- setNames(rep("solid", length(nm)), nm)
  lty[nm %in% SCEN_ALT] <- SCEN_ALT_LTY[nm[nm %in% SCEN_ALT]]
  lty
}

## Scenarios drawn last and hence lying on top of all others. The baseline is
## the reference the remaining curves are read against, so it must stay visible
## where the scenarios coincide with it.
SCEN_ON_TOP <- "baseline"

## The same two vectors keyed by the labels of build_scenarios, which is the
## form the manual scales of ggplot2 need. The level order fixes the drawing order.
scenario_scales <- function(labels) {
  nm   <- names(labels)
  lev  <- unlist(labels, use.names = FALSE)
  top  <- intersect(SCEN_ON_TOP, nm)
  list(levels   = unlist(labels[c(setdiff(nm, top), top)], use.names = FALSE),
       breaks   = lev,
       colour   = setNames(scenario_colours(nm), lev),
       linetype = setNames(scenario_types(nm), lev))
}

## The colour scale, the line type scale and the legend layout of a scenario
## figure in one component.
scale_scenario <- function(sc, nrow = 4) {
  list(scale_colour_manual(values = sc$colour, breaks = sc$breaks),
       scale_linetype_manual(values = sc$linetype, breaks = sc$breaks),
       guides(colour = guide_legend(nrow = nrow), linetype = guide_legend(nrow = nrow)))
}