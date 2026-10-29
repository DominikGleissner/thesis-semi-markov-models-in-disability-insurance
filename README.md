# Semi-Markov Models in Disability Insurance

R code for the numerical analysis (Chapter 6) of my master's thesis
*Semi-Markov Models in Disability Insurance: Sojourn Time Distributions, Actuarial Sensitivity and Tail Risk*
(Technical University of Munich, Department of Mathematics). To receive access to my master's thesis, you can contact me at the email address *dominik.gleissner@gmx.de*.

The code calibrates a Markov baseline from the DAV 2021 I tables and compares it with several semi-Markov scenarios for the sojourn time in the invalid state. The comparison covers the equivalence premium, the prospective reserve and the tail of the policy loss.

## Requirements

- R with RStudio (the versions I used are written in `tables/calibration_numbers.tex`)
- R package `ggplot2` solely for the figures

## Files

- `06_00_common.R`: directories, LaTeX output helpers and figure style
- `06_00_calibration.R`: reading the DAV tables and calibrating the baseline
- `06_00_scenarios.R`: the sojourn time scenarios
- `06_00_engine.R`: the deterministic and the Monte Carlo valuation engine
- `06_01_*.R` to `06_04_*.R`: one script per section of Chapter 6, see below

The `06_00_*.R` files only define functions and settings. They are sourced by the other scripts and are not run on their own.

## Input

The folder `data/` holds the DAV 2021 I second-order tables for male policyholders:

- `dav2021i_active_male_2ndorder.csv` with the columns `age, i_x, q_a_x`
- `dav2021i_invalid_male_2ndorder.csv` with the columns `age, y, r, q`, where `y = 1, ..., 10` is the disability year and `10` corresponds to the ultimate column "10+"

## How to reproduce the results

Set the working directory to the repository folder and make sure that the folders
`output/`, `figures/` and `tables/` exist. Then run the scripts in this order:

1. `06_01_baseline_calibration.R` (Section 6.1.2, calibration of the Markov baseline)
2. `06_01_sojourn_comparison.R` (Section 6.1.3, comparison of the sojourn survival functions)
3. `06_01_numerical_reliability.R` (Section 6.1.4, numerical reliability)
4. `06_02_equivalence_premium.R` (Section 6.2, equivalence premium)
5. `06_03_reserve_surface.R` (Section 6.3, prospective reserve)
6. `06_04_loss_tail.R` (Section 6.4, policy loss and its tail)

The order matters, because the later scripts read files that the earlier ones write to `output/`. All simulations start from fixed seeds, so a new run gives the same numbers.

## Output

- `figures/`: the figures of Chapter 6 as PDF
- `tables/`: the LaTeX tables and the `*_numbers.tex` files, which hold every number quoted in the text as a LaTeX macro
- `output/`: CSV working files that are not included in the thesis, but serve to understand in-between results and to store computations in order to reduce redundancy
