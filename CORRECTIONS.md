# Corrections

## 2026-11-08 — Reverse-MR minimum detectable OR

The values 3.10 / 2.97 in the manuscript and 3.095 / 2.971 in Table S4 and in `results/04_reverse_MR_power_analysis.csv` did not follow from the archived power calculation. The original run log prints 59.2% power at OR = 1.1 and 98.7% at OR = 1.2, which is incompatible with a minimum detectable OR of 3.095.

Recomputed with the same Brion approximation and the archived parameters (outcome n = 459,560; combined R2 = 0.001152; two-sided alpha = 0.05):

- 80% power: **1.129** (reported as 1.13)
- 50% power: **1.089** (reported as 1.09)

`code/04_reverse_MR_power_calculation.R` now uses those fixed parameters, reproduces the archived power curve exactly, and regenerates both result files without a live API call. The earlier placeholder version (12 SNPs, R2 = 0.008) has been replaced. The manuscript's interpretation of the null reverse result was revised accordingly; the primary estimate (OR = 1.232) and all other results are unchanged.
