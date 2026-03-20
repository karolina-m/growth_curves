# =============================================================================
# Power Analysis for Bivariate Stationary LCSM
# Target coupling effect: gamma = 0.20 (both directions)
# Method: Monte Carlo simulation via lavaan
# =============================================================================
# Requirements: install.packages(c("lavaan", "tidyverse"))

library(lavaan)
library(tidyverse)

# =============================================================================
# SECTION 1: Population model
# =============================================================================
# Parameters taken from the N=100 stationary bivariate LCSM fit,
# with coupling set to the target effect size of 0.20

population_model <- '
  # ---- Measurement model ----
  vt1 =~ 1*log_vocab_T1_z;  vt2 =~ 1*log_vocab_T2_z;  vt3 =~ 1*log_vocab_T3_z
  gt1 =~ 1*gram_T1_z;       gt2 =~ 1*gram_T2_z;        gt3 =~ 1*gram_T3_z
  log_vocab_T1_z ~~ 0*log_vocab_T1_z
  log_vocab_T2_z ~~ 0*log_vocab_T2_z
  log_vocab_T3_z ~~ 0*log_vocab_T3_z
  gram_T1_z ~~ 0*gram_T1_z
  gram_T2_z ~~ 0*gram_T2_z
  gram_T3_z ~~ 0*gram_T3_z

  # ---- Latent change scores ----
  vt2 ~ 1*vt1;  vt3 ~ 1*vt2
  gt2 ~ 1*gt1;  gt3 ~ 1*gt2
  dv2 =~ 1*vt2; dv3 =~ 1*vt3
  dg2 =~ 1*gt2; dg3 =~ 1*gt3

  # ---- Self-feedback (from N=100 model) ----
  dv2 ~ -0.388*vt1;  dv3 ~ -0.388*vt2
  dg2 ~ -0.156*gt1;  dg3 ~ -0.156*gt2

  # ---- Coupling: set to target effect size ----
  dv2 ~ 0.20*gt1;  dv3 ~ 0.20*gt2   # gamma_gv = 0.20
  dg2 ~ 0.20*vt1;  dg3 ~ 0.20*vt2   # gamma_vg = 0.20

  # ---- Variances (from N=100 model) ----
  vt1 ~~ 0.990*vt1;  gt1 ~~ 0.990*gt1
  dv2 ~~ 0.161*dv2;  dv3 ~~ 0.559*dv3
  dg2 ~~ 1.485*dg2;  dg3 ~~ 1.192*dg3
  vt2 ~~ 0*vt2;  vt3 ~~ 0*vt3
  gt2 ~~ 0*gt2;  gt3 ~~ 0*gt3

  # ---- Covariances (from N=100 model) ----
  vt1 ~~ 0.072*gt1
  dv2 ~~ 0.071*dg2
  dv3 ~~ 0.272*dg3

  # ---- Means (from N=100 model) ----
  vt1 ~ 0*1;   gt1 ~ 0*1
  dv2 ~ 0.447*1;  dv3 ~ 0.750*1
  dg2 ~ 0.758*1;  dg3 ~ 0.839*1
  log_vocab_T1_z ~ 0*1;  log_vocab_T2_z ~ 0*1;  log_vocab_T3_z ~ 0*1
  gram_T1_z ~ 0*1;       gram_T2_z ~ 0*1;        gram_T3_z ~ 0*1
'

# =============================================================================
# SECTION 2: Analysis model
# =============================================================================
# This is the stationary bivariate LCSM we fit to each simulated dataset.
# Identical to the model in the main script.

analysis_model <- '
  vt1 =~ 1*log_vocab_T1_z;  vt2 =~ 1*log_vocab_T2_z;  vt3 =~ 1*log_vocab_T3_z
  gt1 =~ 1*gram_T1_z;       gt2 =~ 1*gram_T2_z;        gt3 =~ 1*gram_T3_z
  log_vocab_T1_z ~~ 0*log_vocab_T1_z
  log_vocab_T2_z ~~ 0*log_vocab_T2_z
  log_vocab_T3_z ~~ 0*log_vocab_T3_z
  gram_T1_z ~~ 0*gram_T1_z
  gram_T2_z ~~ 0*gram_T2_z
  gram_T3_z ~~ 0*gram_T3_z

  vt2 ~ 1*vt1;  vt3 ~ 1*vt2
  gt2 ~ 1*gt1;  gt3 ~ 1*gt2
  dv2 =~ 1*vt2; dv3 =~ 1*vt3
  dg2 =~ 1*gt2; dg3 =~ 1*gt3

  dv2 ~ beta_v*vt1;  dv3 ~ beta_v*vt2
  dg2 ~ beta_g*gt1;  dg3 ~ beta_g*gt2

  dv2 ~ gamma_gv*gt1;  dv3 ~ gamma_gv*gt2
  dg2 ~ gamma_vg*vt1;  dg3 ~ gamma_vg*vt2

  vt1 ~~ vt1;  gt1 ~~ gt1
  dv2 ~~ dv2;  dv3 ~~ dv3
  dg2 ~~ dg2;  dg3 ~~ dg3
  vt2 ~~ 0*vt2;  vt3 ~~ 0*vt3
  gt2 ~~ 0*gt2;  gt3 ~~ 0*gt3

  vt1 ~~ gt1
  dv2 ~~ dg2
  dv3 ~~ dg3

  vt1 ~ 1;  gt1 ~ 1
  dv2 ~ 1;  dv3 ~ 1
  dg2 ~ 1;  dg3 ~ 1
  log_vocab_T1_z ~ 0*1;  log_vocab_T2_z ~ 0*1;  log_vocab_T3_z ~ 0*1
  gram_T1_z ~ 0*1;       gram_T2_z ~ 0*1;        gram_T3_z ~ 0*1
'

# =============================================================================
# SECTION 3: Simulation function
# =============================================================================

run_simulation <- function(n, n_reps = 1000, seed = 42) {
  
  set.seed(seed)
  
  # Storage for results
  results <- data.frame(
    rep          = 1:n_reps,
    converged    = NA,
    p_gamma_gv   = NA,   # p-value: grammar -> vocab coupling
    p_gamma_vg   = NA    # p-value: vocab -> grammar coupling
  )
  
  for (i in 1:n_reps) {
    
    # Step 1: simulate dataset from population model
    sim_data <- tryCatch(
      simulateData(population_model, sample.nobs = n),
      error = function(e) NULL
    )
    if (is.null(sim_data)) next
    
    # Step 2: fit analysis model to simulated data
    fit <- tryCatch(
      lavaan(analysis_model, data = sim_data,
             estimator = "ML", missing = "FIML"),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        lavaan(analysis_model, data = sim_data,
               estimator = "ML", missing = "FIML")
      )
    )
    if (is.null(fit)) next
    
    # Step 3: check convergence
    if (!lavInspect(fit, "converged")) next
    results$converged[i] <- TRUE
    
    # Step 4: extract p-values for coupling parameters
    pe <- parameterEstimates(fit)
    results$p_gamma_gv[i] <- pe$pvalue[pe$label == "gamma_gv"]
    results$p_gamma_vg[i] <- pe$pvalue[pe$label == "gamma_vg"]
  }
  
  # Step 5: compute power
  converged    <- results %>% filter(converged == TRUE)
  n_converged  <- nrow(converged)
  
  power_gv     <- mean(converged$p_gamma_gv < .05, na.rm = TRUE)
  power_vg     <- mean(converged$p_gamma_vg < .05, na.rm = TRUE)
  power_both   <- mean(converged$p_gamma_gv < .05 &
                         converged$p_gamma_vg < .05, na.rm = TRUE)
  
  data.frame(
    n            = n,
    n_reps       = n_reps,
    n_converged  = n_converged,
    conv_rate    = n_converged / n_reps,
    power_gv     = power_gv,    # grammar -> vocab only
    power_vg     = power_vg,    # vocab -> grammar only
    power_either = 1 - (1 - power_gv) * (1 - power_vg),  # either direction
    power_both   = power_both   # PRIMARY: both directions simultaneously
  )
}

# =============================================================================
# SECTION 4: Run simulation across sample sizes
# =============================================================================
# This will take a few minutes — 1000 reps x 5 sample sizes = 5000 model fits
# Reduce n_reps to 500 for a quicker first run, then increase to 1000

sample_sizes <- c(100, 150, 200, 250, 300)

cat("Running power simulation... this may take a few minutes.\n")

power_results <- map_dfr(sample_sizes, function(n) {
  cat(sprintf("  N = %d...\n", n))
  run_simulation(n = n, n_reps = 1000)
})

# =============================================================================
# SECTION 5: Results table
# =============================================================================

cat("\n========================================\n")
cat("POWER ANALYSIS RESULTS\n")
cat("Coupling effect size: gamma = 0.20 (both directions)\n")
cat("Alpha = .05, 1000 replications per N\n")
cat("========================================\n\n")

power_results %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  print()

cat("\nPrimary criterion (power_both): both gamma_gv AND gamma_vg significant\n")
cat("Secondary criterion (power_either): at least one coupling significant\n")

# =============================================================================
# SECTION 6: Power curve plot
# =============================================================================

power_long <- power_results %>%
  select(n, power_gv, power_vg, power_either, power_both) %>%
  pivot_longer(
    cols      = starts_with("power"),
    names_to  = "criterion",
    values_to = "power"
  ) %>%
  mutate(criterion = recode(criterion,
                            power_gv     = "gram -> vocab only",
                            power_vg     = "vocab -> gram only",
                            power_either = "either direction",
                            power_both   = "both (PRIMARY)"
  ))

ggplot(power_long, aes(x = n, y = power, colour = criterion, linetype = criterion)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0.90, linetype = "dotted", colour = "grey40") +
  annotate("text", x = min(sample_sizes) + 2, y = 0.81,
           label = "80% power", hjust = 0, size = 3.5, colour = "grey40") +
  annotate("text", x = min(sample_sizes) + 2, y = 0.91,
           label = "90% power", hjust = 0, size = 3.5, colour = "grey40") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_manual(values = c(
    "gram -> vocab only" = "#4477AA",
    "vocab -> gram only" = "#EE6677",
    "either direction"   = "#CCBB44",
    "both (PRIMARY)"     = "#228833"
  )) +
  labs(
    title    = "Power analysis: Bivariate stationary LCSM",
    subtitle = "Coupling effect size gamma = 0.20, alpha = .05",
    x        = "Sample size (N)",
    y        = "Statistical power",
    colour   = "Criterion",
    linetype = "Criterion"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("lcsm_power_curve.png", width = 8, height = 5.5, dpi = 150)
cat("\nPower curve saved as lcsm_power_curve.png\n")

# =============================================================================
# SECTION 7: Find minimum N for 80% and 90% power (primary criterion)
# =============================================================================

cat("\n--- Minimum N for primary criterion (both couplings significant) ---\n")

n_80 <- power_results %>%
  filter(power_both >= 0.80) %>%
  slice_min(n) %>%
  pull(n)

n_90 <- power_results %>%
  filter(power_both >= 0.90) %>%
  slice_min(n) %>%
  pull(n)

if (length(n_80) > 0) {
  cat(sprintf("  80%% power achieved at N = %d\n", n_80))
} else {
  cat("  80% power not reached within tested range — extend to larger N\n")
}

if (length(n_90) > 0) {
  cat(sprintf("  90%% power achieved at N = %d\n", n_90))
} else {
  cat("  90% power not reached within tested range — extend to larger N\n")
}

# =============================================================================
# NOTES
# =============================================================================
# 1. Population parameters come from the N=100 stationary bivariate LCSM fit.
#    These are mock data estimates — replace with estimates from pilot or
#    literature when available for the real pre-registration.
#
# 2. Coupling is set symmetrically (gamma_gv = gamma_vg = 0.20).
#    If you expect asymmetric coupling (e.g. vocab -> gram stronger),
#    adjust the population model accordingly and re-run.
#
# 3. MLR estimator is used in the main analyses but ML is used here
#    for simulation speed (simulateData() generates normal data).
#    For the final pre-registration power analysis, consider switching
#    to MLR throughout and generating non-normal data to match real
#    vocab distributions.
#
# 4. The convergence rate column (conv_rate) tells you what proportion
#    of replications produced usable results. If this drops below ~0.90
#    at any N, that sample size is too small for reliable estimation.
