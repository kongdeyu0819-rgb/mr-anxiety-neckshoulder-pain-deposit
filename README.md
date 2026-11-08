# Analysis Code and Derived Results: Anxiety-Related Primary-Care Consultation and Chronic Neck-Shoulder Pain

Analysis code and derived results for:

> Is Genetic Liability to Anxiety-Related Primary-Care Consultation Associated with Chronic Neck-Shoulder Pain? A Two-Sample Mendelian Randomization Study with Implications for Physiotherapy Assessment

## Contents

```
code/     Analysis scripts (R and Python)
results/  Derived result tables (CSV) produced by these scripts
```

## Data sources

No individual-level data were used. All input data are publicly available GWAS summary statistics:

| Dataset | ID | Source |
|---|---|---|
| GP-ANXD (exposure, primary) | ukb-b-6991 | IEU OpenGWAS |
| Chronic neck-shoulder pain (outcome) | ukb-b-16118 | IEU OpenGWAS |
| PGC Major Depressive Disorder | ieu-a-1187 | IEU OpenGWAS |
| FinnGen anxiety disorders | finn-b-KRA_PSY_ANXIETY | IEU OpenGWAS |
| EBI anxiety/panic | ebi-a-GCST90038651 | IEU OpenGWAS |
| UKB anxiety/panic (self-reported) | ukb-b-17243 | IEU OpenGWAS |
| Insomnia | ukb-b-3957 | IEU OpenGWAS |
| C-reactive protein | ieu-b-35 | IEU OpenGWAS |
| Body mass index | ieu-a-2 | IEU OpenGWAS |
| FinnGen R12 cervicalgia (M54.2) | M13_CERVICALGIA | FinnGen public release |
| FinnGen R12 shoulder lesions (M75) | M13_SHOULDER | FinnGen public release |

The full FinnGen R12 summary-statistic files (approx. 800 MB each) are not redistributed here; download them from the FinnGen public release page (https://www.finngen.fi/en/researchers/clinical-endpoints) and place them in a local directory before running `code/07_cross_design_FinnGen_outcome.py`.

## Software

- R 4.6.0 with TwoSampleMR 0.7.8, MendelianRandomization, MRPRESSO, ieugwasr (see `code/install_packages_R460.R`)
- Python 3.13 (standard library only for the cross-design script)

## Analysis pipeline

| Script | Analysis |
|---|---|
| `code/01_strict_threshold_comparison.R` | Primary genome-wide (p < 5e-8) vs relaxed-threshold (p < 1e-5) comparison |
| `code/02_phenotype_comparison_PGC.R` | GP-ANXD vs PGC-MDD phenotype comparison |
| `code/03_MVMR_anxiety_depression.R` | Multivariable MR (GP-ANXD and PGC-MDD) |
| `code/04_reverse_MR_power_calculation.R` | Reverse MR and formal power calculation |
| `code/05_FinnGen_validation.R` | Exposure-side analyses with three further anxiety GWAS datasets |
| `code/06_mediation_MR_analysis.R` | Two-step mediation MR (insomnia, CRP, BMI) |
| `code/07_cross_design_FinnGen_outcome.py` | Exploratory outcome-side cross-design analysis into FinnGen R12 cervicalgia and shoulder lesions |

## Results files

`results/` contains the derived tables corresponding to the manuscript and its supplementary materials (Tables S1-S7 and the cross-design results). Column definitions follow the TwoSampleMR output conventions; `07_cross_design_FinnGen_outcome_results.csv` is produced by the Python script.

## Citation

If you use this code or the derived results, please cite the associated manuscript. The repository is archived on Zenodo under concept DOI 10.5281/zenodo.22764938, which resolves to the latest version; the manuscript cites the version DOI current at submission. See CORRECTIONS.md and RELEASE_NOTES.md for the changes included in this alignment.

## License

This repository is released under the MIT License (see `LICENSE`).


## Correspondence with the revised manuscript

This deposit matches the minor-revision manuscript (Manuscript ID 3082992). The
revised manuscript reports an association, not a causal relationship; the
FinnGen outcome-side analyses are labelled exploratory supportive evidence and
are not described as replication (see `RELEASE_NOTES.md`).

The reverse-MR power calculation was corrected on 2026-11-08. Running
`Rscript code/04_reverse_MR_power_calculation.R` reproduces
`results/04_reverse_MR_power_analysis.csv` (minimum detectable OR 1.129 at 80%
power and 1.089 at 50% power) and `results/04_reverse_MR_power_curve.csv`
without any network access or credentials.
