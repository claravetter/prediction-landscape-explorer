# Fully synthetic public demonstration; no real records or study outputs.
make_demo_data <- function() {
  old <- if(exists(".Random.seed",.GlobalEnv)) .Random.seed else NULL
  on.exit({if(is.null(old)) {if(exists(".Random.seed",.GlobalEnv))rm(".Random.seed",envir=.GlobalEnv)} else assign(".Random.seed",old,.GlobalEnv)})
  set.seed(20260500)

  N <- 1200
  cohorts <- c("Cohort_A", "Cohort_B", "Cohort_C", "Cohort_D", "Cohort_E")
  cohort_probs <- c(0.30, 0.22, 0.18, 0.18, 0.12)

  age <- runif(N, min = 16, max = 35)
  cohort <- sample(cohorts, size = N, replace = TRUE, prob = cohort_probs)

  cogdis <- pmin(pmax(round(rnorm(N, mean = 3.5, sd = 1.8)), 0L), 9L)
  sips_p <- pmin(pmax(round(rnorm(N, mean = 12,  sd = 4)),   0L), 30L)
  sips_n <- pmin(pmax(round(rnorm(N, mean = 10,  sd = 4)),   0L), 30L)
  psychosoz <- pmin(pmax(rnorm(N, mean = 55, sd = 12), 20), 100)
  brainage_sbc <- rnorm(N, mean = 0, sd = 3)

  latent <- rnorm(N, mean = 0, sd = 1)
  EXP_LABEL <- rbinom(N, size = 1, prob = plogis(-1.6 + latent))

  # Per-subject signal-to-noise: classifier signal is strong for older /
  # low-COGDIS subjects and weak for younger / high-COGDIS subjects.
  # snr ranges roughly 0.2 (weak) to 3.2 (strong); Mean_Score = snr*latent + N(0,1).
  # The seed above was chosen so the moderation pattern is visible in
  # windowed AUC trajectories even with finite-sample noise.
  age_norm <- (age - 16) / (35 - 16)
  cog_norm <- cogdis / 9
  snr_i <- 0.2 + 1.5 * age_norm + 1.5 * (1 - cog_norm)
  Mean_Score <- snr_i * latent + rnorm(N, mean = 0, sd = 1)

  platt_fit <- glm(EXP_LABEL ~ Mean_Score, family = binomial())
  Probs_platt <- as.numeric(predict(platt_fit, type = "response"))
  Probs_platt <- pmin(pmax(Probs_platt, 1e-4), 1 - 1e-4)

  df <- tibble(
    group              = cohort,
    event              = as.integer(EXP_LABEL),
    risk_score         = Mean_Score,
    calibrated_prob    = Probs_platt,
    age                = age,
    cognitive_score    = cogdis,
    positive_symptoms  = sips_p,
    negative_symptoms  = sips_n,
    psychosocial_function = psychosoz,
    brain_age_delta    = brainage_sbc
  )
  df
}
