###############################################################################
##
##  06_00_engine.R
##
##  The valuation engines of Chapter 6, all of them working on the grid D_h of
##  Definition 5.2.1. This file then computes:
##
##    (1) the model on the grid    -> build_model to prepare the setup
##    (2) the deterministic engine -> Definition 5.2.1 with (5.10) and (5.11),
##                                    Proposition 5.2.3
##    (3) the Monte Carlo engine   -> Definition 5.4.1
##    (4)-(5) the duration-aware versions of (2) and (3), which read the sojourn
##        survival of the invalid state off a matrix instead of a vector.
##    (6) Validation for the two implementation approaches
##
##  Although the duration-aware version also can compute the duration-independent
##  scenarios by setting the respective parameters to NULL, we keep both versions
##  to ensure the correctness of the entire implementation by comparing the results
##  for duration-independent scenarios computed on the two separate approaches.
##
##  Input. A calibrated grid cal with the columns age, sigma, mu, q_i and
##  p_ia_hat, as written by 06_01_baseline_calibration.R, and a policy list pol
##  with the entries t0, n, beta, f, delta and h.
##
###############################################################################


## ---------------------------------------------------------------------------
## 1. The baseline on the grid
## ---------------------------------------------------------------------------

## Everything the two engines need on the grid T_h: the survival factors
## of the two non-absorbing states and the probability of the destination given a jump.
build_model <- function(cal, pol, split = c("entry", "attained")) {
  split <- match.arg(split)
  N <- round(pol$n / pol$h)
  if (abs(N * pol$h - pol$n) > 1e-9)
    stop("the step size ", pol$h, " does not divide the term ", pol$n)

  grid <- pol$t0 + (0:N) * pol$h
  ## Left ends of all cells (t_{l-1}, t_l], could be attained via grid[1:N] as well
  cell <- grid[-1] - pol$h            

  ## Additional rounding before computing to avoid numerical issues
  at <- function(col, age) {
    k <- pmin(pmax(as.integer(floor(round(age, 9))), min(cal$age)), max(cal$age))
    cal[[col]][match(k, cal$age)]
  }

  q_a <- at("sigma", cell) + at("mu", cell)
  q_i <- at("q_i", cell)

  list(N = N, h = pol$h, grid = grid, split = split,
       ## Sojourn survivals such that Ea(l) = 1 - H_a(t_0, t_l-t_0)
       ## Cumulative hazard at t_0 is the exponent in discrete version:
       ## int_t0^tl q_a(w) dw = h * sum_{m=1}^l q_a((t_{m-1}, t_m]) =: A_a(t_l)
       ## Therefore, A_a(t_l) = -log(E_a(l))
       ## Cumulative hazard at any time is then difference A_a(s+z) - A_a(s)
       ## Then 1 - H_a(t_s, t_l-t_s) = exp( -(A_a(t_l) - A_a(t_s)) )
       ##                            = E_a(l) / E_a(s)
       ##                            = \bar{phi}_a
       Ea = c(1, exp(-cumsum(q_a * pol$h))),  # active sojourn survival
       Ei = c(1, exp(-cumsum(q_i * pol$h))),  # invalid sojourn survival
       w_a     = at("sigma", cell) / q_a,   # prob. a -> i, given a jump in the cell
       w_i     = at("p_ia_hat", cell),      # prob. i -> a, read at the attained age
       w_entry = at("p_ia_hat", grid))      # prob. i -> a, read at the entry age
       ## With this, the prob. of entering the dead state is always 1-w
}


## ---------------------------------------------------------------------------
## 2. Deterministic engine based on Definition 5.2.1 with (5.10) and (5.11),
##                                  Proposition 5.2.3
## ---------------------------------------------------------------------------

## Right-hand side of (5.11), which collapses to (5.10) for u=0 in this approach.
## Invalidity entered at grid point s <= a (current time grid point, t in formula),
## so that u=(a-s)*h. Further, mod and pol refer to the model and the policy,
## va and vi are the vectors of the reserves V. b_a and b_i are the payment families.
step <- function(mod, pol, va, vi, a, s, b_a, b_i) {
  if (a >= mod$N) return(c(a = 0, i = 0))
  l  <- (a + 1):mod$N                         # grid points of the future for sum
  d  <- exp(-pol$delta * mod$h * (l - a))     # discount factor for entire formula
  on <- l <= mod$N - 1                        # the rates vanish at t_0 + n
  df <- (l - s) >= round(pol$f / mod$h)       # deferred period
  wi <- if (mod$split == "entry") mod$w_entry[s + 1] else mod$w_i[l]

  ## Note that the sojourn \bar{phi}_j(t,tau,u) = Ej(l) / Ej(a) in this notation
  c(a = (mod$h * b_a * sum(d * on * mod$Ea[l + 1]) +
         sum(d * (mod$Ea[l] - mod$Ea[l + 1]) * mod$w_a[l] * vi[l + 1])) / mod$Ea[a + 1],
    i = (mod$h * b_i * sum(d * on * df * mod$Ei[l + 1]) +
         sum(d * (mod$Ei[l] - mod$Ei[l + 1]) * wi * va[l + 1])) / mod$Ei[a + 1])
}

## Fresh-entry curves of both states by backward substitution in t as in Proposition 5.2.2.
fresh_entry <- function(mod, pol, b_a, b_i) {
  va <- vi <- numeric(mod$N + 1)
  for (a in (mod$N - 1):0) {
    v <- step(mod, pol, va, vi, a, a, b_a, b_i)
    va[a + 1] <- v[["a"]]
    vi[a + 1] <- v[["i"]]
  }
  list(a = va, i = vi)
}

## The two runs of Proposition 5.2.4 with the benefits alone and a
## unit premium alone. This determines the equivalence premium.
price <- function(mod, pol) {
  ben <- fresh_entry(mod, pol, b_a = 0, b_i = pol$beta)
  ann <- fresh_entry(mod, pol, b_a = 1, b_i = 0)
  list(ben = ben, ann = ann, epv_b = ben$a[1], a_bar = ann$a[1],
       premium = ben$a[1] / ann$a[1])
}

## Reserve at the lattice point (j, u, t), received from two fresh-entry
## curves by the linearity shown in the proof of Proposition 5.2.4.
reserve <- function(mod, pol, pr, j, l_t, l_u, premium) {
  vb <- step(mod, pol, pr$ben$a, pr$ben$i, l_t, l_t - l_u, 0, pol$beta)
  vp <- step(mod, pol, pr$ann$a, pr$ann$i, l_t, l_t - l_u, 1, 0)
  unname(vb[[j]] - premium * vp[[j]])
}


## ---------------------------------------------------------------------------
## 3. Monte Carlo engine for the baseline scenario based on Definition 5.4.1
## ---------------------------------------------------------------------------

## M path losses L^h_{(j, u)}(t) started in state j at the grid point l_t with
## elapsed duration l_u * h. The pair (destination, jump time) of (5.20)
## is drawn in two steps, the jump time first and the destination afterwards,
## which is admissible by (5.33) in Proposition 5.5.2. For the jump time
## the inverse transform method is used with the cumulative hazard A_j.
## The cell of the jump is the first grid point at which A_j exceeds A_j(t_{l_t}) 
## plus a standard exponential. All M paths are advanced jump by jump at once.
simulate_loss <- function(mod, pol, j, l_t, l_u, M, premium) {
  N  <- mod$N
  lf <- round(pol$f / mod$h)
  ## Cumulative hazards of both states
  A  <- list(a = -log(mod$Ea), i = -log(mod$Ei))

  disc <- exp(-pol$delta * mod$h * (0:N)) # factor exp(-delta*(t_l-t_0)) for l=0:N
  cum  <- c(0, cumsum(disc * c(rep(1, N), 0)))   # the rates vanish at t_0 + n
  ## Sum of disc over lo, ..., hi, resp. sum_{l=lo}^li exp(-delta(t_l-t_0))
  sum_disc <- function(lo, hi) {
    lo <- pmin(pmax(lo, 0), N + 1)
    hi <- pmin(hi, N)
    ifelse(hi < lo, 0, cum[hi + 2] - cum[lo + 1])
  }

  st  <- rep(j, M)                 # current state
  pos <- rep(l_t, M)               # grid index of the last jump
  ent <- rep(l_t - l_u, M)         # grid index of the entry into that state
  acc <- numeric(M)                # accumulated payments, starting at value 0
  ## Boolean vector storing whether a simulation is done
  run <- rep(TRUE, M)

  ## Iteration per jump. 
  ## Run as long as there are runs that are not dead and before maturity
  while (any(run)) {
    for (s in c("a", "i")) {
      ## Collect indices of all runs still active and currently in state s.
      k <- which(run & st == s)
      ## If no such runs exist (no active runs in state s), skip the iteration of s.
      if (!length(k)) next
      
      ## Determine the jump time with approach as in inverse transform sampling
      ## stated and explained below Definition 5.4.1.
      ## Recall that for E~Exp(1), we have U := exp(-E)~Unif(0,1) and E = -log(U)
      ## The cumulative hazard is extracted via -log(E_s) as well. Thus, add E instead U.
      ## Let x := A_s(t_p) + rexp(1) for t_p reached as pos[k]+1 element. 
      ## Then we search the index b such that
      ## A_s[b] <= x < A_s[b+1], respectively A_s(t_{b-1}) <= x < A_s(t_b)
      ## This yields b = min { l in {p+1,..., N} such that A_s(t_l) > x }
      ## and further P(b > l) = P(Exp(1) >= A_s(t_l) - A_s(t_p)) 
      ##                                  = E_s(l) / E_s(p)
      ##                                  = never-left prob. 1 - H_a(t_p, t_l-t_p)
      b  <- findInterval(A[[s]][pos[k] + 1] + rexp(length(k)), A[[s]])

      ## Discounted to t_0 accumulated (sojourn) payments depending on the current state
      if (s == "a") {       # premium is paid in active state
        acc[k] <- acc[k] - premium * mod$h * sum_disc(pos[k] + 1, b - 1)
      } else {              # benefit is paid in invalid state
        acc[k] <- acc[k] + pol$beta * mod$h *
          sum_disc(pmax(pos[k] + 1, ent[k] + lf), b - 1) # incorporates deferred period
      }

      ## Destination of the jump, if there is one before the end of the term
      ## Set all currently active runs to inactive by standard
      run[k] <- FALSE
      ## b <= N means there is a jump, keep only those with jumps.
      ## b > N means that J_{m+1}=infty and no further jump occurs until maturity
      kk <- k[b <= N] 
      ## Only do something if there are real jumps left
      if (length(kk)) {
        bb <- b[b <= N]
        ## Read probability of destinations based on current state and type of split
        w  <- if (s == "a") mod$w_a[bb] else
              if (mod$split == "entry") mod$w_entry[ent[kk] + 1] else mod$w_i[bb]
        ## Prob. of entering dead state is always 1-w regardless of the kind of w
        ## Use Bernoulli trails with parameter w to determine this next state
        to <- ifelse(runif(length(kk)) < w, if (s == "a") "i" else "a", "d")
        ## Update the current state, the last jump index and the entrance time
        st[kk]  <- to
        pos[kk] <- bb
        ent[kk] <- bb
        ## Reactivate / keep those trails that are neither in dead state nor after maturity
        run[kk] <- to != "d" & bb < N
      }
    }
  }
  ## Returns the vector of losses disocunted to the valuation time t instead t_0
  ## Above we computed acc = sum_l exp(-delta * (t_l - t_0)) * c_l
  ## Using exp(-delta(t_l - t)) = exp(-delta(t_l - t_0)) / exp(-delta(t - t_0)),
  ## this explains the denominator disc[l_t + 1] which is just exp(-delta(t - t_0)).
  acc / disc[l_t + 1]
}


## ---------------------------------------------------------------------------
## 4. The duration-aware deterministic engine
##
##  The same scheme as in Section 2 of this file. 
##  The only change is that the survival and the kernel increment of the invalid state
##  are read off the row s of SUR instead of the vector Ei to include the duration.
##  Nothing in the active state changes.  A scenario is a function surv(x,z) 
##  returning 1 - G(x,z) for one entry age x and a vector of durations z.
##  Setting surv = NULL gives the baseline back, which is used to validate the implementation.
## ---------------------------------------------------------------------------

## Compute the duration-dependent sojourn survivals of the invalid state.
## Input surv(x,z) is the survival function of the invalid state.
## If surv == NULL, we use the baseline survival \bar{phi}_i from above.
sojourn_matrix <- function(mod, pol, surv = NULL) {
  S <- matrix(0, mod$N + 1, mod$N + 1)
  for (s in 0:mod$N) {
    k <- (s + 1):(mod$N + 1)
    S[s + 1, k] <- if (is.null(surv)) mod$Ei[k] / mod$Ei[s + 1]
                   else surv(pol$t0 + s * pol$h, (0:(mod$N - s)) * pol$h)
  }
  S
}

## Same function as above, but add the sojourn matrix of the invalid state as well.
build_general <- function(cal, pol, surv = NULL) {
  mod <- build_model(cal, pol)
  mod$Si <- sojourn_matrix(mod, pol, surv)
  mod
}

## Same function as above, but includes the sojourn survival Si instead Ei.
## Note that here, the invalid split wi always reads at the entry age of invalidity.
step_g <- function(mod, pol, va, vi, a, s, b_a, b_i) {
  if (a >= mod$N) return(c(a = 0, i = 0))
  l  <- (a + 1):mod$N
  d  <- exp(-pol$delta * mod$h * (l - a))
  on <- l <= mod$N - 1
  df <- (l - s) >= round(pol$f / mod$h)
  Si <- mod$Si[s + 1, ]

  c(a = (mod$h * b_a * sum(d * on * mod$Ea[l + 1]) +
         sum(d * (mod$Ea[l] - mod$Ea[l + 1]) * mod$w_a[l] * vi[l + 1])) / mod$Ea[a + 1],
    i = (mod$h * b_i * sum(d * on * df * Si[l + 1]) +
         sum(d * (Si[l] - Si[l + 1]) * mod$w_entry[s + 1] * va[l + 1])) / Si[a + 1])
}

## Same function as above, but calls the step_g function instead the step function.
fresh_entry_g <- function(mod, pol, b_a, b_i) {
  va <- vi <- numeric(mod$N + 1)
  for (a in (mod$N - 1):0) {
    v <- step_g(mod, pol, va, vi, a, a, b_a, b_i)
    va[a + 1] <- v[["a"]]
    vi[a + 1] <- v[["i"]]
  }
  list(a = va, i = vi)
}

## Same function as above, but calls fresh_entry_g instead of fresh_entry.
price_g <- function(mod, pol) {
  ben <- fresh_entry_g(mod, pol, b_a = 0, b_i = pol$beta)
  ann <- fresh_entry_g(mod, pol, b_a = 1, b_i = 0)
  list(ben = ben, ann = ann, epv_b = ben$a[1], a_bar = ann$a[1],
       premium = ben$a[1] / ann$a[1])
}

## Same function as above, but calls the step_g function instead the step function.
reserve_g <- function(mod, pol, pr, j, l_t, l_u, premium) {
  vb <- step_g(mod, pol, pr$ben$a, pr$ben$i, l_t, l_t - l_u, 0, pol$beta)
  vp <- step_g(mod, pol, pr$ann$a, pr$ann$i, l_t, l_t - l_u, 1, 0)
  unname(vb[[j]] - premium * vp[[j]])
}

## ---------------------------------------------------------------------------
## 5. The duration-aware Monte Carlo engine
##
##  M occupancies in a and i. As above, the pair (destination, jump time) of 
##  (5.20) is drawn in two steps. In the active state the inverse transform
##  method is applied to the cumulative hazard Aa, as in simulate_loss.  In the
##  invalid state it is applied to the survival row itself, so the jump cell is
##  the first grid point at which 1 - G(x, .) falls below a uniform draw. 
##  All changes to simulate_loss are documented in the function itself.
## ---------------------------------------------------------------------------

simulate_occupancy <- function(mod, pol, M) {
  N  <- mod$N
  lf <- round(pol$f / mod$h)
  Aa <- -log(mod$Ea)
  
  disc <- exp(-pol$delta * mod$h * (0:N))
  cum  <- c(0, cumsum(disc * c(rep(1, N), 0)))
  sum_disc <- function(lo, hi) {
    lo <- pmin(pmax(lo, 0), N + 1)
    hi <- pmin(hi, N)
    ifelse(hi < lo, 0, cum[hi + 2] - cum[lo + 1])
  }
  
  ## Always start in active state with a fresh entry.
  st  <- rep("a", M)
  pos <- rep(0L, M)
  ent <- rep(0L, M)
  ## Instead of tracking the accumulated payments, we record the discounted occupancies.
  ## These are just the discounted times in the respective state. It can be obtained
  ## by choosing premium and payment functions as 1 and tracking the two payment
  ## streams for the premiums and the benefits separately.
  ## Since the premium and benefit rate are just real constants in our analysis,
  ## one can just multiply the occupancies with the respective rate and add them
  ## to get the overall loss by linearity: loss = premium * occ$a - benefit * occ$i.
  occ <- list(a = numeric(M), i = numeric(M))
  run <- rep(TRUE, M)
  
  while (any(run)) {
    for (s in c("a", "i")) {
      k <- which(run & st == s)
      if (!length(k)) next
      
      ## Determine the jump cell, on the cumulative hazard in the active state,
      ## and on the survival row of the entry age in the invalid one.
      ## As for the active state, instead of doing it on the U-scale directly,
      ## transform with -log to the cumulative hazard scale and add an Exp(1).
      ## This makes all values positive and uses the same idea for active and invalid.
      if (s == "a") {
        b <- findInterval(Aa[pos[k] + 1] + rexp(length(k)), Aa)
      } else {
        b <- integer(length(k))
        z <- rexp(length(k))
        for (e in unique(ent[k])) {
          j <- which(ent[k] == e)
          b[j] <- e + findInterval(z[j], -log(mod$Si[e + 1, (e + 1):(N + 1)]))
        }
      }
      
      if (s == "a") {
        occ$a[k] <- occ$a[k] + sum_disc(pos[k] + 1, b - 1)
      } else {
        occ$i[k] <- occ$i[k] + sum_disc(pmax(pos[k] + 1, ent[k] + lf), b - 1)
      }
      
      run[k] <- FALSE
      kk <- k[b <= N]
      if (length(kk)) {
        bb <- b[b <= N]
        ## All transitions happen at fresh entries, so that the split is
        ## always read at the entry age, not the attained age.
        w  <- if (s == "a") mod$w_a[bb] else mod$w_entry[ent[kk] + 1]
        to <- ifelse(runif(length(kk)) < w, if (s == "a") "i" else "a", "d")
        st[kk]  <- to
        pos[kk] <- bb
        ent[kk] <- bb
        run[kk] <- to != "d" & bb < N
      }
    }
  }
  ## Scale all values to the underlying grid with step size h.
  list(a = pol$h * occ$a, i = pol$h * occ$i)
}


## ===========================================================================
## 6. Check to run once after the insertion, not part of any section
##    To run it, just set CHECK to TRUE, else to FALSE
##
##  With surv = NULL the sojourn matrix is the calibrated kernel itself, so
##  simulate_occupancy must reproduce simulate_loss of Section 3 path by path.
##  Both draw the active state from the same cumulative hazard, so the same seed
##  gives the same paths and the two losses have to agree to machine precision.
## ===========================================================================

CHECK = TRUE

if (CHECK) {
  source("06_00_common.R")
  pol <- list(t0 = 30, n = 35, beta = 1, f = 0, delta = 0.03, h = 1 / 24)
  cal <- read.csv(file.path(OUTPUT_DIR, "baseline_calibration_grid.csv"))
  mod <- build_general(cal, pol)
  pr  <- price_g(mod, pol)
  M <- 1e5
  set.seed(1)
  l1 <- simulate_loss(mod, pol, "a", 0, 0, M, pr$premium)
  set.seed(1)
  oc <- simulate_occupancy(mod, pol, M)
  deviance <- l1 - (pol$beta * oc$i - pr$premium * oc$a)
  stopifnot(max(abs(deviance)) < 1e-12)
  # This value should be very very small!
  cat("largest pathwise gap :", max(abs(deviance)), "\n")
  message("Validation of the two Monte-Carlo-engines finished.")
}
