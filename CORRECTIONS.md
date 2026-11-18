# Corrections

## 2026-11-08 — Reverse-MR minimum detectable OR

The values 3.10 / 2.97 in the manuscript and 3.095 / 2.971 in Table S4 and in `results/04_reverse_MR_power_analysis.csv` did not follow from the archived power calculation. The original run log prints 59.2% power at OR = 1.1 and 98.7% at OR = 1.2, which is incompatible with a minimum detectable OR of 3.095.

Recomputed with the same Brion approximation and the archived parameters (outcome n = 459,560; combined R2 = 0.001152; two-sided alpha = 0.05):

- 80% power: **1.129** (reported as 1.13)
- 50% power: **1.089** (reported as 1.09)

`code/04_reverse_MR_power_calculation.R` now uses those fixed parameters, reproduces the archived power curve exactly, and regenerates both result files without a live API call. The earlier placeholder version (12 SNPs, R2 = 0.008) has been replaced. The manuscript's interpretation of the null reverse result was revised accordingly. A separate correction to the primary-analysis figures was made on 2026-11-18 and is recorded below.

## 2026-11-18 — Primary-analysis figures

During the current revision the primary-analysis figures in the manuscript were checked against the archived analytical output `results/verify_estimates.log`, which reports 40 genome-wide significant (p < 5 × 10⁻⁸), non-ambiguous instruments. The corresponding estimates are:

- IVW OR = 1.249 (95% CI 1.126–1.385), p = 2.53 × 10⁻⁵
- Weighted median OR = 1.185
- MR-Egger OR = 1.249 (95% CI 0.356–4.384), p = 0.729; intercept p = 0.979
- Cochran's Q = 45.5 (df = 39), p = 0.220

The manuscript, Table 1, Table S1 and the legends of Figures 2 and S1 now use this 40-instrument set. Earlier figures based on a 41-instrument set (IVW OR = 1.232) are withdrawn; no archived analysis output supports that set. `results/verify_estimates.log` is included in this deposit as the source of the corrected values. The relaxed-threshold sensitivity analysis (182 instruments, IVW OR = 1.190) is unchanged. The reverse-MR minimum detectable odds ratios (1.129 at 80% power and 1.089 at 50% power) are a separate quantity and are not affected by this correction.
