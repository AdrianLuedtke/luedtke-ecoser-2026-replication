# Replication material

Lüdtke, A., Burkart, S., & Lippert, C. Valuing Biodiversity-Linked Ecosystem Services in Agricultural Landscapes: A Systematic Review and Estimation of an Exploratory Meta-Analytic Value Function with Implications for Farm-Level Decisions and Agri-Environmental Policy. Submitted to *Ecosystem Services* (manuscript ECOSER-D-25-01027).

## Files

- `Revision_meta_regression_final_2026-10-23.R`
  Replication script for the meta-regression (Table 7 of the manuscript, Appendix A.4, A.7 and A.8). Reads the dataset, applies the Cook's distance filter (D > 4/n), estimates the levels model with HC1 standard errors, runs the robustness checks, and writes all tables to `Revision_tables_final.xlsx`.
- `Waste_Treatment_Data_set_final.xlsx`
  The 53 ESVD records of water purification values from European agricultural land (biome "Intensive land use", continent "Europe", service "Water purification"; extraction August 2025, ESVD version April 2024v.1.0) after removal of four extreme values, with the five regional moderators and their sources per record. The Bavaria 3 m record of Sieber et al. (2010) is corrected from the ESVD entry of 431 EUR to the 143 EUR reported in the study (Appendix A.15).
- `Supplement_S1_ESVD_observations_final.xlsx`
  All 57 extracted ESVD records with inclusion status, moderator values and statistical sources (Supplement S1 of the manuscript).
- `Supplement_S2_review_studies_final_2026-10-23.xlsx`
  The 71 studies included in the systematic review and the 51 studies excluded at full-text stage, with coding labels (service, method family, country, spatial scale, where discussed) and exclusion reasons (Supplement S2 of the manuscript).

## How to run

1. Place the R script and `Waste_Treatment_Data_set_final.xlsx` in the same folder, or set `base_dir` at the top of the script.
2. Run the script in R (tested with R 4.3.1). Missing packages (readxl, lmtest, sandwich, car, strucchange, broom, writexl) are installed automatically.
3. The output `Revision_tables_final.xlsx` contains the coefficient tables, diagnostics, robustness checks and the Baden-Württemberg and Bavaria predictions.

## Data source

ESVD records are from the Ecosystem Services Valuation Database (Brander et al., 2024, https://doi.org/10.1016/j.ecoser.2024.101606) and are redistributed here for replication purposes. Moderator sources are listed in Appendix A.14 of the manuscript.

## Contact

Adrian Lüdtke, University of Hohenheim, adrian.luedtke@uni-hohenheim.de
