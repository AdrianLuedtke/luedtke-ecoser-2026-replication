# ==================================================================
# Meta-regression of water purification values from European
# agricultural land (ESVD records, n = 48 after influence filtering)
#
# Replication script for:
# Lüdtke, Burkart & Lippert, "Valuing Biodiversity-Linked Ecosystem
# Services in Agricultural Landscapes: A Systematic Review and
# Estimation of an Exploratory Meta-Analytic Value Function with
# Implications for Farm-Level Decisions and Agri-Environmental
# Policy", Ecosystem Services (manuscript ECOSER-D-25-01027).
#
# Input : Waste_Treatment_Data_set_final.xlsx (53 observations after
#         removal of four extreme values; the Bavaria 3 m record of
#         Sieber et al. 2010 is already corrected from 431 to 143 EUR,
#         see Appendix A.15 and Supplement S1)
# Output: Revision_tables_final.xlsx (Table 7, Appendix A.4, A.7, A.8)
#
# Run from the folder that contains the input file, or set `base_dir`.
# Robust tables are built with lmtest::coeftest() (HC1), because
# broom::tidy() on an lm object ignores a vcov argument.
# ==================================================================

# --- Package installation (runs only for packages not yet installed)
required_packages <- c("readxl", "lmtest", "sandwich", "car",
                       "strucchange", "broom", "writexl")
missing_packages <- required_packages[!required_packages %in%
                                        rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

# --- Packages ------------------------------------------------------
library(readxl)      # read_excel
library(lmtest)      # coeftest, coefci, waldtest, bptest, dwtest
library(sandwich)    # vcovHC, NeweyWest
library(car)         # vif, ncvTest
library(strucchange) # sctest (Rainbow)
library(broom)       # tidy (applied to coeftest output only)
library(writexl)     # write_xlsx

# --- Paths ---------------------------------------------------------
base_dir <- "."                                   # folder with the data file
in_file  <- file.path(base_dir, "Waste_Treatment_Data_set_final.xlsx")
out_file <- file.path(base_dir, "Revision_tables_final.xlsx")

# --- Helper: robust coefficient table ------------------------------
robust_table <- function(model, vc, bonferroni = TRUE) {
  ct <- coeftest(model, vcov. = vc)
  ci <- coefci(model, vcov. = vc)
  out <- data.frame(
    term     = rownames(ct),
    estimate = ct[, 1],
    se       = ct[, 2],
    t        = ct[, 3],
    p        = ct[, 4],
    ci_low   = ci[, 1],
    ci_high  = ci[, 2],
    row.names = NULL
  )
  if (bonferroni) {
    k <- sum(out$term != "(Intercept)")
    out$p_bonferroni <- ifelse(out$term == "(Intercept)", NA,
                               pmin(out$p * k, 1))
  }
  out
}

fit_stats <- function(model, vc, label) {
  w <- waldtest(model, vcov = vc, test = "F")
  data.frame(
    specification = label,
    n             = nobs(model),
    adj_r2        = summary(model)$adj.r.squared,
    wald_F        = w$F[2],
    wald_df1      = abs(w$Df[2]),
    wald_df2      = w$Res.Df[1],
    wald_p        = w$`Pr(>F)`[2],
    max_vif       = if (length(coef(model)) > 2) max(vif(model)) else NA
  )
}

# --- Data import (53 observations after extreme-value removal) -----
# File: Waste_Treatment_Data_set_final.xlsx
#   - moderators rebuilt from the primary statistical sources listed in
#     Appendix A.14; sources per record in Farm_Size_Source and
#     Moderator_Notes
#   - Bavaria 3 m record (Sieber et al. 2010) corrected from the ESVD
#     entry 431 EUR to 143 EUR (Table 1 of the study), i.e. 15.1 instead
#     of 45.6 Int.$ 2020 ha-1 yr-1 (Appendix A.15)
d53 <- read_excel(in_file, sheet = 1)
stopifnot(nrow(d53) == 53)
stopifnot(any(d53$Location_Name == "Free State of Bavaria" &
              abs(d53$Value - 15.1) < 0.05))     # corrected record present

d53$lnOrg <- log(d53$OrganicLandShare + 1)

# --- Step 1: Cook's distance filter (D > 4/n) ----------------------
m0 <- lm(Value ~ GDP_per_Capita + lnOrg + Population_Density +
           Average_Farm_Size + Publication_Year, data = d53)
cooks  <- cooks.distance(m0)
drop   <- which(cooks > 4 / nrow(d53))
d      <- d53[-drop, ]
cat("Observations removed by Cook's D:", length(drop), "\n")
print(d53[drop, c("Location_Name", "Country", "Publication_Year", "Valuation_Method", "Value")])

yr_mean <- mean(d$Publication_Year)
d$Publication_Year_c <- d$Publication_Year - yr_mean
cat("N analytic sample:", nrow(d), " | mean publication year:", round(yr_mean, 2), "\n")

# Sanity check of the corrected farm-size variable
print(unique(d53[, c("Location_Name", "Average_Farm_Size", "Farm_Size_Source")]))
stopifnot(all(!is.na(d53$Average_Farm_Size)))

# Indicators used in the revision
d$method_other <- as.integer(d$Valuation_Method != "Opportunity Cost")
d$eco_other    <- as.integer(d$Ecosystem != "Crop Land")
d$scale_nat    <- as.integer(d$Study_Site_Scale == "National")
d$german_block <- as.integer(d$Country == "DEU" & d$Publication_Year == 2010)
cat("German block (Sieber et al., 2010) n =", sum(d$german_block), "\n")

# Distinct covariate profiles
n_prof_all <- nrow(unique(d[, c("GDP_per_Capita", "OrganicLandShare",
                                "Population_Density", "Average_Farm_Size")]))
n_prof_blk <- nrow(unique(d[d$german_block == 1,
                            c("GDP_per_Capita", "OrganicLandShare",
                              "Population_Density", "Average_Farm_Size")]))
cat("Distinct covariate profiles: full =", n_prof_all,
    "| German block =", n_prof_blk, "\n")

# --- Step 2: Core model (main-text Table 5) ------------------------
f_base <- Value ~ GDP_per_Capita + lnOrg + Population_Density +
  Average_Farm_Size + Publication_Year_c
m <- lm(f_base, data = d)

vc_hc1 <- vcovHC(m, type = "HC1")
vc_nw  <- NeweyWest(m, lag = 2, prewhite = FALSE)   # R default: no df adjustment

tab5_hc1 <- robust_table(m, vc_hc1)
tab5_nw  <- robust_table(m, vc_nw)
fit5     <- rbind(fit_stats(m, vc_hc1, "Base, HC1"),
                  fit_stats(m, vc_nw,  "Base, Newey-West lag 2"))
print(tab5_hc1, digits = 4)
print(fit5, digits = 4)
# Replication in Python with the corrected record: n = 48, adj. R2 = 0.397,
# GDP 0.00103 (HC1 p = 0.098), lnOrg -274.7 (p = 0.015), pub. year -6.10 (p < 0.001),
# population density and farm size insignificant; same five Cook's D exclusions.

# --- Step 3: Diagnostics -------------------------------------------
diag_out <- data.frame(
  test = c("Rainbow", "Breusch-Pagan", "NCV", "Shapiro-Wilk (residuals)"),
  p    = c(sctest(m, type = "Rainbow")$p.value,
           bptest(m)$p.value,
           ncvTest(m)$p,
           shapiro.test(residuals(m))$p.value)
)
vif_base <- data.frame(term = names(vif(m)), vif = vif(m), row.names = NULL)
print(diag_out); print(vif_base)

# --- Step 4: Additional study descriptors (Appendix A.4.5, Table Z) -
specs <- list(
  "Base"                 = f_base,
  "+ Method dummy"       = update(f_base, . ~ . + method_other),
  "+ Ecosystem dummy"    = update(f_base, . ~ . + eco_other),
  "+ Scale dummy"        = update(f_base, . ~ . + scale_nat),
  "Block replaces year"  = Value ~ GDP_per_Capita + lnOrg + Population_Density +
                             Average_Farm_Size + german_block
)
tabZ_coef <- do.call(rbind, lapply(names(specs), function(s) {
  mm <- lm(specs[[s]], data = d)
  cbind(specification = s, robust_table(mm, vcovHC(mm, type = "HC1"),
                                        bonferroni = FALSE))
}))
tabZ_vif <- do.call(rbind, lapply(names(specs), function(s) {
  mm <- lm(specs[[s]], data = d)
  v  <- vif(mm)
  data.frame(specification = s, term = names(v), vif = v, row.names = NULL)
}))
tabZ_fit <- do.call(rbind, lapply(names(specs), function(s) {
  mm <- lm(specs[[s]], data = d)
  fit_stats(mm, vcovHC(mm, type = "HC1"), s)
}))
print(tabZ_fit, digits = 4)
# Check the VIF column: publication year is collinear with the method,
# ecosystem and scale dummies, which is the reason for excluding them.

# --- Step 5: German study block alone (Appendix A.4.8) -------------
d_blk <- subset(d, german_block == 1)
m_blk <- lm(Value ~ GDP_per_Capita + lnOrg + Population_Density +
              Average_Farm_Size, data = d_blk)
tab_blk <- robust_table(m_blk, vcovHC(m_blk, type = "HC1"), bonferroni = FALSE)
fit_blk <- fit_stats(m_blk, vcovHC(m_blk, type = "HC1"), "German block, HC1")
print(tab_blk, digits = 4); print(fit_blk, digits = 4)
# German block: n = 37 in the corrected dataset.

# --- Step 6: Dependent variable in levels vs. log (Appendix Table 4) -
m_log <- lm(log(Value) ~ GDP_per_Capita + lnOrg + Population_Density +
              Average_Farm_Size + Publication_Year_c, data = d)
tab_log <- robust_table(m_log, vcovHC(m_log, type = "HC1"), bonferroni = FALSE)
fit_log <- rbind(fit_stats(m, vc_hc1, "Levels, HC1"),
                 fit_stats(m_log, vcovHC(m_log, type = "HC1"), "Log DV, HC1"))
fit_log$shapiro_p <- c(shapiro.test(residuals(m))$p.value,
                       shapiro.test(residuals(m_log))$p.value)

# --- Step 7: Organic share, ln(share + 1) vs. logit (Appendix Table 5)
eps <- 1e-6
d$logitOrg <- with(d, log(pmin(pmax(OrganicLandShare, eps), 1 - eps) /
                          (1 - pmin(pmax(OrganicLandShare, eps), 1 - eps))))
m_logit <- lm(Value ~ GDP_per_Capita + logitOrg + Population_Density +
                Average_Farm_Size + Publication_Year_c, data = d)
tab_logit <- robust_table(m_logit, vcovHC(m_logit, type = "HC1"), bonferroni = FALSE)
fit_logit <- fit_stats(m_logit, vcovHC(m_logit, type = "HC1"), "Logit organic share, HC1")

# --- Step 8: Log-transformed moderators (Appendix Table 6) ---------
m_lnmod <- lm(Value ~ log(GDP_per_Capita) + lnOrg + log(Population_Density) +
                log(Average_Farm_Size) + Publication_Year_c, data = d)
tab_lnmod <- robust_table(m_lnmod, vcovHC(m_lnmod, type = "HC1"), bonferroni = FALSE)
fit_lnmod <- fit_stats(m_lnmod, vcovHC(m_lnmod, type = "HC1"), "Log moderators, HC1")

# --- Step 9: Trimming 5th-95th percentile of the DV (Appendix Table 7)
q <- quantile(d$Value, c(0.05, 0.95))
d_trim <- subset(d, Value >= q[1] & Value <= q[2])
m_trim <- lm(f_base, data = d_trim)
tab_trim <- robust_table(m_trim, vcovHC(m_trim, type = "HC1"), bonferroni = FALSE)
fit_trim <- fit_stats(m_trim, vcovHC(m_trim, type = "HC1"), "Trimmed 5-95%, HC1")

# --- Step 10: Before vs. after influence filtering (Appendix Table 8)
d53$Publication_Year_c <- d53$Publication_Year - mean(d53$Publication_Year)
m_full53 <- lm(f_base, data = d53)
tab_full53 <- robust_table(m_full53, vcovHC(m_full53, type = "HC1"), bonferroni = FALSE)
fit_full53 <- fit_stats(m_full53, vcovHC(m_full53, type = "HC1"), "n = 53 before Cook's D, HC1")

# --- Step 11: Leave-one-out ranges (Appendix Table 9) --------------
loo <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  mi <- lm(f_base, data = d[-i, ])
  rt <- robust_table(mi, vcovHC(mi, type = "HC1"), bonferroni = FALSE)
  cbind(dropped = i, rt[, c("term", "estimate", "p")])
}))
tab_loo <- do.call(rbind, lapply(split(loo, loo$term), function(x) data.frame(
  term     = x$term[1],
  est_min  = min(x$estimate), est_max = max(x$estimate),
  p_min    = min(x$p),        p_max   = max(x$p),
  share_sig_5pct = mean(x$p < 0.05)
)))

# --- Step 12: Appendix A.7 descriptive statistics (n = 50) ---------
# 12a. Sample coverage by year (Table 13)
tab13 <- data.frame(
  n_obs    = nrow(d),
  year_min = min(d$Publication_Year),
  year_max = max(d$Publication_Year),
  year_mean = mean(d$Publication_Year)
)
year_hist <- as.data.frame(table(Publication_Year = d$Publication_Year))

# 12b. Outcome distribution (Table 14), incl. skewness and kurtosis
skewness <- function(x) {
  n <- length(x); m <- mean(x)
  (sum((x - m)^3) / n) / (sum((x - m)^2) / n)^(3/2)
}
kurtosis_excess <- function(x) {
  n <- length(x); m <- mean(x)
  (sum((x - m)^4) / n) / (sum((x - m)^2) / n)^2 - 3
}
q <- quantile(d$Value, c(0.10, 0.25, 0.50, 0.75, 0.90))
tab14 <- data.frame(
  n      = nrow(d),
  mean   = mean(d$Value),
  sd     = sd(d$Value),
  median = median(d$Value),
  min    = min(d$Value),
  max    = max(d$Value),
  IQR    = IQR(d$Value),
  skew   = skewness(d$Value),
  kurtosis_excess = kurtosis_excess(d$Value),
  p10 = q[1], p25 = q[2], p50 = q[3], p75 = q[4], p90 = q[5],
  row.names = NULL
)

# 12c. Covariate summary (Table 15)
cov_vars <- c("GDP_per_Capita", "OrganicLandShare",
              "Population_Density", "Average_Farm_Size")
tab15 <- do.call(rbind, lapply(cov_vars, function(v) data.frame(
  variable = v,
  n        = sum(!is.na(d[[v]])),
  missing  = sum(is.na(d[[v]])),
  mean     = mean(d[[v]], na.rm = TRUE),
  sd       = sd(d[[v]], na.rm = TRUE),
  min      = min(d[[v]], na.rm = TRUE),
  max      = max(d[[v]], na.rm = TRUE)
)))

# 12d. Correlations among covariates (Tables 16 and 17)
tab16 <- as.data.frame(round(cor(d[, cov_vars], method = "pearson"), 3))
tab16 <- cbind(variable = rownames(tab16), tab16); rownames(tab16) <- NULL
tab17 <- as.data.frame(round(cor(d[, cov_vars], method = "spearman"), 3))
tab17 <- cbind(variable = rownames(tab17), tab17); rownames(tab17) <- NULL

# 12e. Heterogeneity patterns (Tables 18 to 20)
tab18 <- as.data.frame.matrix(table(d$Valuation_Method, d$Country))
tab18 <- cbind(method = rownames(tab18), tab18); rownames(tab18) <- NULL
tab19 <- as.data.frame.matrix(table(d$Valuation_Method, d$Study_Site_Scale))
tab19 <- cbind(method = rownames(tab19), tab19); rownames(tab19) <- NULL
tab20 <- do.call(rbind, lapply(split(d$Value, d$Valuation_Method), function(v)
  data.frame(n = length(v), mean = mean(v), median = median(v), IQR = IQR(v))))
tab20 <- cbind(method = rownames(tab20), tab20); rownames(tab20) <- NULL

# 12f. Country medians (Section 4.2 text)
country_medians <- do.call(rbind, lapply(split(d$Value, d$Country), function(v)
  data.frame(n = length(v), median = median(v), mean = mean(v))))
country_medians <- cbind(country = rownames(country_medians), country_medians)
rownames(country_medians) <- NULL

print(tab14, digits = 4); print(tab15, digits = 5); print(tab20, digits = 4)

# --- Step 13: Composition of the analytic sample (Section 4.2) -----
composition <- as.data.frame(table(Location = d$Location_Name,
                                   Country = d$Country,
                                   Year    = d$Publication_Year,
                                   Method  = d$Valuation_Method,
                                   Scale   = d$Study_Site_Scale))
composition <- subset(composition, Freq > 0)
desc50 <- data.frame(
  n = nrow(d), mean = mean(d$Value), sd = sd(d$Value),
  median = median(d$Value), min = min(d$Value), max = max(d$Value)
)

# --- Step 14: Predictions for Baden-Württemberg and Bavaria (A.8) --
# --- Covariate bundles (Appendix Table 21) -------------------------
# GDP per capita: Eurostat nama_10r_3gdp, 2020
# Organic land share: Regionaldatenbank Deutschland, 2020
# Population density: Eurostat demo_r_d3dens, 2022
# Average farm size: Landwirtschaftszählung 2020 (BW 36.0 ha, BY 37.0 ha)
bundles <- data.frame(
  Region             = c("Baden-Württemberg", "Bavaria"),
  GDP_per_Capita     = c(46600, 48700),
  OrganicLandShare   = c(0.123, 0.121),
  Population_Density = c(313, 187),
  Average_Farm_Size  = c(36.0, 37.0)
)
bundles$lnOrg              <- log(bundles$OrganicLandShare + 1)
bundles$Publication_Year_c <- 0            # publication year fixed at the sample mean

# --- Point predictions ---------------------------------------------
X <- model.matrix(delete.response(terms(m)), data = bundles)
fit <- as.vector(X %*% coef(m))

# --- HC1-robust standard errors and intervals ----------------------
se_mean <- sqrt(rowSums((X %*% vc_hc1) * X))                 # SE of the mean prediction
s2      <- sum(residuals(m)^2) / df.residual(m)          # residual variance
se_obs  <- sqrt(se_mean^2 + s2)                          # SE of a single observation
tcrit   <- qt(0.975, df.residual(m))

pred <- data.frame(
  Region        = bundles$Region,
  Prediction    = fit,
  SE_mean       = se_mean,
  CI_mean_lower = fit - tcrit * se_mean,
  CI_mean_upper = fit + tcrit * se_mean,
  PI_obs_lower  = fit - tcrit * se_obs,
  PI_obs_upper  = fit + tcrit * se_obs
)
print(round(pred[, -1], 2))

# --- Decomposition of the difference between the two regions -------
b <- coef(m)
contrib <- c(
  GDP_per_capita     = b["GDP_per_Capita"]    * diff(bundles$GDP_per_Capita),
  Organic_land_share = b["lnOrg"]             * diff(bundles$lnOrg),
  Population_density = b["Population_Density"]* diff(bundles$Population_Density),
  Average_farm_size  = b["Average_Farm_Size"] * diff(bundles$Average_Farm_Size)
)
decomp <- data.frame(Moderator = names(contrib),
                     Contribution_BY_minus_BW = as.vector(contrib))
decomp <- rbind(decomp, data.frame(Moderator = "Total",
                                   Contribution_BY_minus_BW = sum(contrib)))
print(round(decomp$Contribution_BY_minus_BW, 2))


# --- Export all tables to one workbook -----------------------------
write_xlsx(list(
  Table5_HC1            = tab5_hc1,
  Table5_NeweyWest      = tab5_nw,
  Table5_fit            = fit5,
  Diagnostics           = diag_out,
  VIF_base              = vif_base,
  A4_5_TableZ_coef      = tabZ_coef,
  A4_5_TableZ_VIF       = tabZ_vif,
  A4_5_TableZ_fit       = tabZ_fit,
  A4_8_GermanBlock      = tab_blk,
  A4_8_GermanBlock_fit  = fit_blk,
  A4_2_LogDV            = tab_log,
  A4_2_fit              = fit_log,
  A4_3_Logit            = tab_logit,
  A4_3_Logit_fit        = fit_logit,
  A4_4_LogModerators    = tab_lnmod,
  A4_4_LogMod_fit       = fit_lnmod,
  A4_Trim               = tab_trim,
  A4_Trim_fit           = fit_trim,
  A4_Filter_n53         = tab_full53,
  A4_Filter_n53_fit     = fit_full53,
  A4_LeaveOneOut        = tab_loo,
  A7_T13_YearCoverage   = tab13,
  A7_T13_YearHistogram  = year_hist,
  A7_T14_Distribution   = tab14,
  A7_T15_Covariates     = tab15,
  A7_T16_Pearson        = tab16,
  A7_T17_Spearman       = tab17,
  A7_T18_MethodCountry  = tab18,
  A7_T19_MethodScale    = tab19,
  A7_T20_OutcomeMethod  = tab20,
  A7_CountryMedians     = country_medians,
  Composition           = composition,
  Descriptives          = desc50,
  A8_Predictions        = pred,
  A8_Decomposition      = decomp,
  A8_Bundles            = bundles
), out_file)

cat("Tables written to:", out_file, "\n")
