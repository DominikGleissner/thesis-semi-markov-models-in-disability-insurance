###############################################################################
##
##  05_05_WeibullFamilies.R
##  
##  Mean-matched Weibull shape sweep with the values for Table 5.1.
##
##  Recall the following formulas:
##  1)  Weibull(kappa, alpha):  E[X] = alpha * Gamma(1 + 1/kappa)
##  2)  Mean-matching condition:  
##      alpha(kappa) * Gamma(1 + 1/kappa) = m_i(s)
##      <=>  alpha(kappa) / m_i(s) = 1 / Gamma(1 + 1/kappa)
##  3)  Coefficient of variation:
##      cv(kappa) = sqrt( Gamma(1 + 2/kappa) / Gamma(1 + 1/kappa)^2 - 1 )
##
###############################################################################

main = function () {
  kappa <- c(0.50, 0.75, 1.00, 1.50, 2.00)
  
  scale_ratio <- 1 / gamma(1 + 1 / kappa)
  cv          <- sqrt(gamma(1 + 2 / kappa) / gamma(1 + 1 / kappa)^2 - 1)
  
  aging_class <- ifelse(kappa < 1, "DFR",
                        ifelse(kappa == 1, "IFR ∩ DFR", "IFR"))
  
  result <- data.frame(
    kappa       = kappa,
    scale_ratio = round(scale_ratio, 4),
    aging_class = aging_class,
    cv          = round(cv, 4)
  )
  
  print(result, row.names = FALSE)
}

main()
