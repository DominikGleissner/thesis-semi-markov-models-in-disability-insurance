###############################################################################
##
##  06_00_calibration.R
##
##  Everything that turns the DAV 2021 I second-order tables into the Markov
##  baseline S^(0) of Section 6.1.2. Precise steps:
##
##    (1) invert the tabulated active probabilities into sigma_x, mu_x
##                                     -> Proposition 6.1.4
##    (2) invert the tabulated invalid probabilities into rho_y(x), nu_y(x),
##        build h_hat_i and the tabulated sojourn survival 1 - H_hat_i
##                                     -> Approach as in Proposition 6.1.4
##                                     -> Definition 5.5.9
##    (3) evaluate the tabulated mean sojourn m_hat_i(x)
##                                     -> Definition 6.1.5
##    (4) check the mean-residual-life condition m_hat(x) < 1 + m_hat(x+1)
##                                     -> Equation (6.15) in Proposition 6.1.6
##    (5) solve the backward recursion for the duration-free exit intensity q_x
##                                     -> Proposition 6.1.6
##    (6) embedded jump probabilities for the active and the invalid state
##                                     -> Definition 5.5.1
##                                     -> Equations (6.11) and (6.19)
##
##  Input (in folder INPUT_DIR)
##     dav2021i_active_male_2ndorder.csv    columns: age, i_x, q_a_x
##     dav2021i_invalid_male_2ndorder.csv   columns: age, y, r, q
##  Both files must be complete rectangles over their age range; gaps and
##  duplicates are an error, not something to be filled silently.
##
###############################################################################


## ---------------------------------------------------------------------------
## 0. Settings of the calibration
## ---------------------------------------------------------------------------

Y_ULT     <- 10L               # index of the ultimate column "10+"
MEAN_RULE <- "constant_force"  # "trapezoidal" or "constant_force"
Q_BRACKET <- c(1e-10, 1e3)     # search bracket for (MeanPreservingRecursion)
ROOT_TOL  <- .Machine$double.eps^0.75


## ---------------------------------------------------------------------------
## 1. Reading and validating the tables
## ---------------------------------------------------------------------------

## Net inversion of a single decrement, p = 1 - exp(-alpha).
## log1p keeps the precision at the small probabilities of the young ages.
net_intensity <- function(p) -log1p(-p)

read_active_table <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("age", "i_x", "q_a_x") %in% names(d)))
    stop("the active table must have the columns age, i_x, q_a_x")

  d <- d[order(d$age), c("age", "i_x", "q_a_x")]
  d$age <- as.integer(d$age)
  if (anyDuplicated(d$age) > 0 || !all(diff(d$age) == 1L))
    stop("the active table must cover consecutive integer ages without repetition")

  ## Each value must lie in (0, 1) on its own
  bad <- which(!(d$i_x > 0 & d$i_x < 1 & d$q_a_x > 0 & d$q_a_x < 1))
  if (length(bad) > 0)
    stop("tabulated active probability outside (0, 1) at age(s) ",
         paste(d$age[bad], collapse = ", "))
  d
}

read_invalid_table <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  if ("n" %in% names(d) && !("y" %in% names(d))) names(d)[names(d) == "n"] <- "y"
  if (!all(c("age", "y", "r", "q") %in% names(d)))
    stop("the invalid table must have the columns age, y, r, q")
  d$age <- as.integer(d$age)
  d$y   <- as.integer(d$y)

  ages <- sort(unique(d$age))
  if (!all(diff(ages) == 1L))
    stop("the invalid table has gaps in the age column")
  if (length(unique(d$y)) != Y_ULT || !all(sort(unique(d$y)) == seq_len(Y_ULT)))
    stop("the invalid table must resolve exactly the disability years 1..", Y_ULT)
  if (anyDuplicated(d[, c("age", "y")]) > 0)
    stop("the invalid table contains duplicated age x disability-year cells")
  if (nrow(d) != length(ages) * Y_ULT)
    stop("the invalid table is not a complete age x disability-year rectangle")

  bad <- which(!(d$r >= 0 & d$r < 1 & d$q >= 0 & d$q < 1))
  if (length(bad) > 0)
    stop("tabulated invalid probability outside [0, 1) at ",
         paste0("(age ", d$age[bad], ", y ", d$y[bad], ")", collapse = ", "))

  ## Long format from csv -> age x disability-year matrix
  fill <- function(v) {
    M <- matrix(NA_real_, nrow = length(ages), ncol = Y_ULT,
                dimnames = list(ages, seq_len(Y_ULT)))
    M[cbind(match(d$age, ages), d$y)] <- v
    if (anyNA(M)) stop("the invalid table has missing cells")
    M
  }
  r <- fill(d$r)
  q <- fill(d$q)

  rho <- net_intensity(r)
  nu  <- net_intensity(q)
  cc  <- rho + nu                    # total exit intensity of the cell
  if (any(cc <= 0))
    stop("a cell of the invalid table has r = q = 0, so the tabulated sojourn ",
         "survival never decays there and m_hat_i would diverge")

  list(x_min = min(ages), omega = max(ages), ages = ages,
       r = r, q = q, rho = rho, nu = nu, cc = cc,
       h = -expm1(-cc),              # h_hat_i = 1 - (1 - r)(1 - q) = 1 - exp(-c)
       w = rho / cc)                 # recovery share rho / (rho + nu)
}

## Continuation rule of Definition TabulatedMeanSojourn: the rates are continued
## constantly in the attained age beyond omega, the disability year is capped.
cell <- function(tab, what, x, y) {
  ## Compute the indices for the matrix
  ## Ages are taken into [x_min,omega] and then suited to matrix indices starting at 1
  ix <- min(max(as.integer(x), tab$x_min), tab$omega) - tab$x_min + 1L
  ## Disability year rounded into [1,Y_ULT]
  jy <- min(as.integer(y), Y_ULT)
  tab[[what]][ix, jy]
}

## Vectorised form of cell(), with the same continuation rule.  The empirical
## scenario reads one cell per grid step, so the scalar version would be called
## way too many times.
cell_at <- function(tab, what, x, y) {
  ix <- pmin(pmax(as.integer(x), tab$x_min), tab$omega) - tab$x_min + 1L
  jy <- pmin(as.integer(y), Y_ULT)
  tab[[what]][cbind(ix, jy)]
}

## First disability year from which every rate along a invalidity entered at age x is
## constant, so that the rest of the series can be summed in closed form.
y_star <- function(tab, x) max(Y_ULT, tab$omega - as.integer(x) + 1L)


## ---------------------------------------------------------------------------
## 2. Step (1): inversion in the active state
## ---------------------------------------------------------------------------

invert_active_intensities <- function(active) {
  data.frame(age   = active$age,
             i_x   = active$i_x,
             q_a_x = active$q_a_x,
             sigma = net_intensity(active$i_x),
             mu    = net_intensity(active$q_a_x))
}


## ---------------------------------------------------------------------------
## 3. Steps (2) and (3): tabulated sojourn survival and mean sojourn
## ---------------------------------------------------------------------------

## 1 - H_hat_i(x, y) for y = 0, ..., y_max; entry [y + 1] holds y.
tabulated_survival <- function(tab, x, y_max) {
  out    <- numeric(y_max + 1L)
  out[1] <- 1
  s      <- 1
  for (y in seq_len(y_max)) {
    s <- s * exp(-cell(tab, "cc", x + y - 1L, y))
    out[y + 1L] <- s
  }
  out
}

## Both readings of the mean sojourn of a invalidity entered at the integer age x,
##    trapezoidal    : 1/2 + sum_{y >= 1} (1 - H_hat_i(x, y))
##    constant_force : sum_{y >= 1} (1 - H_hat_i(x, y-1)) (1 - e^{-c_y}) / c_y,
## From y_star on the tail is geometric in both cases; it is added in closed form.
tabulated_mean_sojourn <- function(tab, x) {
  ys <- y_star(tab, x)
  s <- 1; m_tr <- 0.5; m_cf <- 0
  for (y in seq_len(ys - 1L)) {
    cy   <- cell(tab, "cc", x + y - 1L, y)
    m_cf <- m_cf + s * (-expm1(-cy)) / cy
    s    <- s * exp(-cy)
    m_tr <- m_tr + s
  }
  c_inf <- cell(tab, "cc", tab$omega, Y_ULT)
  c(trapezoidal    = m_tr + s * exp(-c_inf) / (-expm1(-c_inf)),
    constant_force = m_cf + s / c_inf)
}

## 1 - H_hat_i(s, z) of (6.21) for one real entry age s and a non-decreasing vector of 
## durations z >= 0, returned in the order of z.
## The integrand v -> c_{ceil(v)}(floor(s + v)) of (6.21) is a step function.
## It can only jump where the disability year ceil(v) or the attained age floor(s + v) changes:
##   the integer durations      v = 1, 2, ..., floor(zmax)
##   the integer attained ages  v = k - s,  k = ceil(s), ..., floor(s + zmax).
## Together with 0 and the durations z themselves these points form the
## partition 0 = b_0 < b_1 < ... < b_K = zmax, held in brk. On every cell
## (b_{k-1}, b_k) the integrand is constant, so it is read once at the midpoint v_k.
tabulated_survival_at <- function(tab, s, z) {
  z <- as.numeric(z)
  if (is.unsorted(z) || any(z < 0))
    stop("the durations of the empirical scenario must be a non-negative increasing grid")
  zmax <- max(z)
  ages <- if (ceiling(s) <= s+zmax) seq(ceiling(s), s + zmax, by = 1) - s else numeric(0)
  brk <- sort(unique(c(0, z, seq_len(floor(zmax)), ages)))
  d   <- diff(brk)
  mid <- brk[-length(brk)] + d / 2
  Lam <- c(0, cumsum(d * cell_at(tab, "cc", floor(s + mid), pmax(ceiling(mid), 1))))
  exp(-Lam[match(z, brk)])
}

## Mean of the empirical sojourn law at a real entry age s.
empirical_mean <- function(tab, s, dz = 1 / 24) {
  z0 <- max(Y_ULT, tab$omega - as.integer(floor(s)) + 1L)
  z  <- seq(dz, z0, by = dz)
  cc <- cell_at(tab, "cc", floor(s + z - dz / 2), pmax(ceiling(z - dz / 2), 1))
  S  <- exp(-cumsum(dz * cc))
  Sm <- c(1, S[-length(S)])
  sum(Sm * (-expm1(-dz * cc)) / cc) +
    S[length(S)] / cell(tab, "cc", tab$omega, Y_ULT)
}


## ---------------------------------------------------------------------------
## 4. Steps (4) and (5): mean-preserving calibration of q_i
## ---------------------------------------------------------------------------

## g_m(q) = (1 - e^{-q}) / q + e^{-q} m, extended continuously by g_m(0) = 1 + m.
## Scalar q only.
g_recursion <- function(q, m) {
  if (q < 1e-12) 1 + m else -expm1(-q) / q + exp(-q) * m
}

## The unique positive root of g_m(q) = target, which exists precisely when
## 0 < target < 1 + m (Step 2 of the proof of Proposition 6.1.6).
solve_recursion <- function(m, target, age) {
  f  <- function(q) g_recursion(q, m) - target
  lo <- Q_BRACKET[1]
  hi <- Q_BRACKET[2]
  if (f(lo) <= 0 || f(hi) >= 0)
    stop("no sign change of (MeanPreservingRecursion) on the search bracket at ",
         "age ", age, "; widen Q_BRACKET or inspect the tabulated means")
  uniroot(f, interval = c(lo, hi), tol = ROOT_TOL)$root
}

calibrate_invalid_intensity <- function(tab) {
  ages <- tab$ages
  n    <- length(ages)

  m_both <- vapply(ages, function(x) tabulated_mean_sojourn(tab, x),
                   c(trapezoidal = 0, constant_force = 0))
  m_tr   <- m_both["trapezoidal", ]
  m_cf   <- m_both["constant_force", ]
  m_hat  <- if (MEAN_RULE == "trapezoidal") m_tr else m_cf

  ## Condition (6.16) in Proposition 6.1.6; beyond omega the mean is constant
  m_next <- c(m_hat[-1], m_hat[n])
  ok     <- m_hat > 0 & m_hat < 1 + m_next
  if (!all(ok))
    stop("condition (6.16) fails at age(s) ",
         paste(ages[!ok], collapse = ", "),
         "; Proposition 6.1.6 does not apply there")

  q_i    <- numeric(n)
  q_i[n] <- 1 / m_hat[n]                     # tail value, Step 3
  if (n >= 2)                                # backward substitution, Step 4
    for (k in (n - 1L):1L)
      q_i[k] <- solve_recursion(m_hat[k + 1L], m_hat[k], ages[k])

  data.frame(age = ages, m_hat = m_hat, m_hat_trapezoidal = m_tr,
             m_hat_constant_force = m_cf, q_i = q_i)
}

## The calibrated q_i at an attained age, continued constantly as in the table.
q_at_age <- function(cal, x) {
  x <- pmin(pmax(as.integer(floor(x)), min(cal$age)), max(cal$age))
  cal$q_i[match(x, cal$age)]
}

## Baseline sojourn survival 1 - H^(0)_i(s, z) for the calibrated piecewise-constant q_i,
## on a grid of integer entry ages s and real durations z >= 0.
baseline_survival <- function(cal, s, z) {
  stopifnot(abs(s - round(s)) < 1e-10)
  vapply(z, function(zz) {
    if (zz <= 0) return(1)
    k    <- floor(zz)
    ## full captures the full years such that integral becomes sum
    full <- if (k >= 1) sum(q_at_age(cal, s + seq_len(k) - 1L)) else 0
    ## (zz - k) * q_at_age(cal, s + k) is the share of the started year
    exp(-(full + (zz - k) * q_at_age(cal, s + k)))
  }, numeric(1))
}


## ---------------------------------------------------------------------------
## 5. Step (6): embedded jump probabilities
## ---------------------------------------------------------------------------

## Using (6.11), over the year of age x + k the contribution to p_ai is S(k) sigma_{x+k} / 
## c_{x+k} (1 - exp(-c_{x+k})), with S(k) the probability of still being active at age x + k.
## The two rows are accumulated separately, so that p_ai + p_ad - 1 stays a genuine check.
active_jump_probabilities <- function(act) {
  ages    <- act$age
  omega_a <- max(ages)
  ## Compute the age index for the matrix under continuation assumption
  idx     <- function(x) min(max(as.integer(x), min(ages)), omega_a) - min(ages) + 1L

  sig_inf <- act$sigma[idx(omega_a)]
  mu_inf  <- act$mu[idx(omega_a)]
  c_inf   <- sig_inf + mu_inf

  ## Initialize the vectors and them compute them step by step
  p_ai <- p_ad <- tail_mass <- numeric(length(ages))
  for (i in seq_along(ages)) {
    x <- ages[i]
    s <- 1; a <- 0; d <- 0
    for (k in seq_len(max(omega_a - x, 0L)) - 1L) {
      j  <- idx(x + k)
      cj <- act$sigma[j] + act$mu[j]
      a  <- a + s * act$sigma[j] / cj * (-expm1(-cj))
      d  <- d + s * act$mu[j]    / cj * (-expm1(-cj))
      s  <- s * exp(-cj)
    }
    p_ai[i] <- a + s * sig_inf / c_inf
    p_ad[i] <- d + s * mu_inf  / c_inf
    tail_mass[i] <- s              # mass the constant continuation carries
  }
  data.frame(age = ages, p_ai = p_ai, p_ad = p_ad, a_tail_mass = tail_mass)
}

## The exit composition of a whole invalidity entered at the integer age x,
## weighted with rho_y / (rho_y + nu_y) stems from (6.19). The tail
## from y_star on adds w_inf (1 - H_hat_i(x, y_star - 1)) to the recovery share.
invalid_jump_probabilities <- function(tab, ages) {
  w_inf <- cell(tab, "w", tab$omega, Y_ULT)

  ## Initialize the vectors and them compute them step by step
  p_ia <- p_id <- tail_mass <- numeric(length(ages))
  for (i in seq_along(ages)) {
    x  <- ages[i]
    ys <- y_star(tab, x)
    s <- 1; a <- 0; d <- 0
    for (y in seq_len(ys - 1L)) {
      xy <- x + y - 1L
      hy <- cell(tab, "h", xy, y)
      wy <- cell(tab, "w", xy, y)
      a  <- a + s * wy * hy
      d  <- d + s * (1 - wy) * hy
      s  <- s * (1 - hy)
    }
    p_ia[i] <- a + w_inf * s
    p_id[i] <- d + (1 - w_inf) * s
    tail_mass[i] <- s              # 1 - H_hat_i(x, y_star - 1)
  }
  data.frame(age = ages, p_ia_hat = p_ia, p_id_hat = p_id, i_tail_mass = tail_mass)
}


## ---------------------------------------------------------------------------
## 6. Alternative duration-free rules as recorded in remark after Proposition 6.1.6
## ---------------------------------------------------------------------------

## Instead of matching the mean sojourn one may simply declare one disability
## year to carry the duration-free rate.  Both rules change the mean sojourn of
## the baseline, which is why they are robustness runs and not the main model.
alternative_intensities <- function(tab, ages) {
  data.frame(age            = ages,
             q_i_first_year = sapply(ages, function(x) cell(tab, "cc", x, 1L)),
             q_i_ultimate   = sapply(ages, function(x) cell(tab, "cc", x, Y_ULT)))
}


## ---------------------------------------------------------------------------
## 7. Diagnostics
## ---------------------------------------------------------------------------

## Remark after Proposition 6.1.4: the exact probability of an
## a -> i transition within the year of age against the DAV correction factor
## i_hat_x (1 - q_hat^a_x / 2) as in (6.9) and (6.10).
dependent_inception_check <- function(act) {
  cx    <- act$sigma + act$mu
  exact <- act$sigma / cx * (-expm1(-cx))
  dav   <- act$i_x * (1 - act$q_a_x / 2)
  data.frame(age = act$age, exact = exact, dav = dav,
             abs_dev = abs(exact - dav), rel_dev = abs(exact - dav) / exact)
}

## Second remark below Proposition 6.1.6: how far the splitting ratio
## rho_y / (rho_y + nu_y) moves over the disability years, age by age.
splitting_ratio_spread <- function(tab, ages) {
  w_min <- w_max <- numeric(length(ages))
  for (i in seq_along(ages)) {
    w <- sapply(seq_len(Y_ULT), function(y) cell(tab, "w", ages[i], y))
    w_min[i] <- min(w)
    w_max[i] <- max(w)
  }
  data.frame(age = ages, w_min = w_min, w_max = w_max, w_range = w_max - w_min)
}
