###############################################################################
# Reverse-MR power calculation (CNSP -> GP-ANXD)
#
# Reproduces Table S4(b) of the revised manuscript and
# results/04_reverse_MR_power_analysis.csv.
#
# Method: Brion M-JA, Shakhbazov K, Visscher PM (2013). Calculating statistical
# power in Mendelian randomization studies. Int J Epidemiol 42(5):1497-1501.
# The formula and the bisection search are identical to those in the original
# analysis script; only the inputs are fixed to the archived values, so the
# calculation runs without a live OpenGWAS query or an API credential.
#
# Run from the repository root:  Rscript code/04_reverse_MR_power_calculation.R
#
# Correction (2026-11-08): an earlier deposit contained placeholder inputs
# (12 SNPs, R2 = 0.008) and a result table whose minimum detectable ORs
# did not follow from its own power curve. Both have been replaced with the
# recomputed values below. See CORRECTIONS.md.
###############################################################################

N_outcome <- 459560      # GP-ANXD (ukb-b-6991) total sample size
R2_total  <- 0.001152    # combined R2 of the 26 pre-harmonisation CNSP instruments
n_snps    <- 26          # selected at p < 1e-5, before harmonisation
alpha     <- 0.05        # two-sided

# Brion IVW power approximation, as used in the original analysis script.
mr_power_ivw <- function(n_outcome, beta, R2_total, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha / 2)
  n_eff   <- n_outcome * R2_total
  ncp     <- sqrt(n_eff * beta^2)
  pnorm(ncp - z_alpha) + pnorm(-ncp - z_alpha)
}

# Minimum detectable OR: bisection on log(OR) in [0, 2], width <= 0.001.
detectable_or <- function(n_outcome, R2_total, power_target, alpha = 0.05) {
  lo  <- 0
  hi  <- 2
  mid <- (lo + hi) / 2
  while ((hi - lo) > 0.001) {
    if (mr_power_ivw(n_outcome, mid, R2_total, alpha) < power_target) {
      lo <- mid
    } else {
      hi <- mid
    }
    mid <- (lo + hi) / 2
  }
  exp(mid)
}

min_or_80 <- detectable_or(N_outcome, R2_total, 0.8)
min_or_50 <- detectable_or(N_outcome, R2_total, 0.5)

cat("Reverse MR power calculation (Brion approximation)\n")
cat(sprintf("  N = %d, R2 = %.6f, SNPs = %d, two-sided alpha = %.2f\n",
            N_outcome, R2_total, n_snps, alpha))
cat(sprintf("  Minimum detectable OR (80%% power): %.3f\n", min_or_80))
cat(sprintf("  Minimum detectable OR (50%% power): %.3f\n", min_or_50))
for (orv in c(1.1, 1.2, 1.5)) {
  cat(sprintf("  Power at OR = %.1f: %.1f%%\n", orv,
              100 * mr_power_ivw(N_outcome, log(orv), R2_total)))
}

power_results <- data.frame(
  Analysis                     = "Reverse MR (CNSP -> GP-ANXD)",
  Threshold                    = "p<1e-5",
  N_SNPs                       = n_snps,
  Total_R2                     = R2_total,
  Mean_F                       = 21.4,
  Min_F                        = 19.5,
  Outcome_N                    = N_outcome,
  Actual_reverse_OR            = 1.025,
  Actual_reverse_P             = 0.238,
  Min_detectable_OR_80power    = round(min_or_80, 3),
  Min_detectable_OR_50power    = round(min_or_50, 3),
  Power_adequate               = sprintf("ADEQUATE_FOR_OR>=%.2f", round(min_or_80, 2))
)

or_range <- seq(1.0, 2.5, by = 0.01)
power_curve <- data.frame(
  OR    = or_range,
  Power = sapply(log(or_range), function(b) mr_power_ivw(N_outcome, b, R2_total))
)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(power_results, "results/04_reverse_MR_power_analysis.csv", row.names = FALSE)
write.csv(power_curve,   "results/04_reverse_MR_power_curve.csv",   row.names = FALSE)
cat("Wrote results/04_reverse_MR_power_analysis.csv and results/04_reverse_MR_power_curve.csv\n")
