#---- Fit turnover models ----
run_adj_glm_v2 <- 
  function(data, outcome_var, exposure_var, confounder_list, profession,
           cluster_var = "loc_2", family_type = "poisson", model_label,
           output_dir = "out/models/main") {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # --- Extended abbreviation maps ---
    
    out_short <- case_when(
      outcome_var == "phq_co_poi"      ~ "dep",
      outcome_var == "gad_co_poi"      ~ "anx",
      outcome_var == "suic_idea_poi"   ~ "suic",
      outcome_var == "cage_co_poi"     ~ "alc",
      outcome_var == "work_18_poi"     ~ "turn",
      TRUE ~ gsub("[^a-z0-9]", "", tolower(outcome_var))  # fallback
    )
    
    exp_short <- case_when(
      exposure_var == "work_8_dic"    ~ "hrs",
      exposure_var == "work_11_dic"   ~ "night",
      exposure_var == "work_12_dic"   ~ "shift",
      exposure_var == "work_19_dic"   ~ "influence",
      exposure_var == "work_20_dic"   ~ "breaks",
      exposure_var == "work_21_dic"   ~ "cols",
      exposure_var == "work_22_dic"   ~ "superiors",
      exposure_var == "work_23_dic"   ~ "deadline",
      exposure_var == "work_24_dic"   ~ "angry",
      exposure_var == "work_30_dic"   ~ "har",
      exposure_var == "work_31_dic"   ~ "threats",
      exposure_var == "work_32_dic"   ~ "viol",
      exposure_var == "work_33_dic"   ~ "bul",
      exposure_var == "phqads_z"      ~ "phqads",
      exposure_var == "phq_co"        ~ "dep",
      exposure_var == "gad_co"        ~ "anx",
      exposure_var == "wc_demands"    ~ "dem",
      exposure_var == "wc_hazards"    ~ "haz",
      exposure_var == "wc_resources"  ~ "res",
      exposure_var == "wc_demands_z"   ~ "dem_z",
      exposure_var == "wc_hazards_z"   ~ "haz_z",
      exposure_var == "wc_resources_z" ~ "res_z",
      TRUE ~ gsub("[^a-z0-9]", "", tolower(exposure_var))  # fallback
    )
    
    prof_short <- case_when(
      profession == "Doctor" ~ "doc",
      profession == "Nurse"  ~ "nur",
      profession == "All"    ~ "all",
      profession == "GP"     ~ "gp",
      profession == "Non-GP" ~ "nongp",
      TRUE ~ gsub("[^a-z0-9]", "", tolower(profession))
    )
    
    model_suffix <- case_when(
      model_label == "Partially adjusted" ~ "str",
      model_label == "Fully adjusted"     ~ "adj",
      TRUE ~ tolower(gsub(" ", "_", model_label))
    )
    
    message("out: ", out_short, " | exp: ", exp_short, " | prof: ", prof_short)
    
    # --- Detect exposure type ---
    is_binary_exposure <- is.factor(data[[exposure_var]]) && 
      all(levels(data[[exposure_var]]) %in% c("No", "Yes"))
    
    # --- Subset and clean ---
    # Data is pre-filtered externally; profession label is for file naming only
    ds_sub <- data
    
    active_confounders <- setdiff(confounder_list, "work_2")
    
    required_vars <- unique(c(outcome_var, exposure_var, cluster_var, 
                              active_confounders))
    
    ds_clean <- ds_sub |> select(all_of(required_vars)) |> drop_na()
    
    n_size <- nrow(ds_clean)          # <-- must come BEFORE the messages
    if (n_size < 10) return(NULL)
    
    # Now the messages are safe
    message("n_size: ", n_size)
    message("outcome range: ", paste(range(ds_clean[[outcome_var]], na.rm = TRUE), collapse = "-"))
    
    n_size <- nrow(ds_clean)
    if (n_size < 10) return(NULL)
    
    if (is_binary_exposure) {
      ds_clean[[exposure_var]] <- factor(ds_clean[[exposure_var]], 
                                         levels = c("No", "Yes"))
    }
    
    y_raw <- as.numeric(ds_clean[[outcome_var]])
    if (max(y_raw, na.rm = TRUE) > 1) ds_clean[[outcome_var]] <- y_raw - 1
    
    curr_family <- if (family_type == "gaussian") gaussian() else 
      poisson(link = "log")
    
    # --- Fit models (caching unchanged) ---
    file_name_ml <- paste0("glm_", out_short, "_", exp_short, "_", 
                           prof_short, "_", model_suffix, "_ml.rds")
    file_path_ml <- file.path(output_dir, file_name_ml)
    
    ml <- NULL
    if (file.exists(file_path_ml)) {
      ml <- tryCatch(readRDS(file_path_ml), error = function(e) {
        file.remove(file_path_ml); NULL })
    }
    if (is.null(ml)) {
      f_ml <- as.formula(paste0(
        outcome_var, " ~ ", exposure_var, " + ",
        paste(active_confounders, collapse = " + "),
        " + (1 | ", cluster_var, ")"))
      ml <- glmer(f_ml, data = ds_clean, family = curr_family,
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 2e5)))
      saveRDS(ml, file_path_ml)
    }
    
    file_name_rob <- paste0("glm_", out_short, "_", exp_short, "_", 
                            prof_short, "_", model_suffix, "_rob.rds")
    file_path_rob <- file.path(output_dir, file_name_rob)
    
    rob <- NULL
    if (file.exists(file_path_rob)) {
      rob <- tryCatch(readRDS(file_path_rob), error = function(e) {
        file.remove(file_path_rob); NULL })
    }
    if (is.null(rob)) {
      f_rob <- as.formula(paste0(
        outcome_var, " ~ ", exposure_var, " + ",
        paste(active_confounders, collapse = " + ")))
      rob <- glm(f_rob, data = ds_clean, family = curr_family)
      saveRDS(rob, file_path_rob)
    }
    
    # --- Extract coefficients (type-aware) ---
    coef_name_ml  <- if (is_binary_exposure) paste0(exposure_var, "Yes") else 
      exposure_var
    coef_name_rob <- coef_name_ml
    
    coef_ml <- summary(ml)$coefficients
    beta_ml  <- coef_ml[coef_name_ml, "Estimate"]
    se_ml    <- coef_ml[coef_name_ml, "Std. Error"]
    pr_ml    <- exp(beta_ml)
    ci_ml    <- exp(beta_ml + c(-1, 1) * 1.96 * se_ml)
    
    vcov_rob <- sandwich::vcovCL(rob, cluster = ~ loc_2, data = ds_clean)
    coef_rob  <- summary(rob)$coefficients
    beta_rb   <- coef_rob[coef_name_rob, "Estimate"]
    se_rb     <- sqrt(vcov_rob[coef_name_rob, coef_name_rob])
    pr_rb     <- exp(beta_rb)
    ci_rb     <- exp(beta_rb + c(-1, 1) * 1.96 * se_rb)
    
    # --- Marginal risks (binary exposure only) ---
    if (is_binary_exposure) {
      dat_unexp <- ds_clean
      dat_exp   <- ds_clean
      dat_unexp[[exposure_var]] <- factor("No",  levels = c("No", "Yes"))
      dat_exp[[exposure_var]]   <- factor("Yes", levels = c("No", "Yes"))
      risk_unexp_str <- sprintf("%.1f%%", mean(predict(ml, newdata = dat_unexp,
                                                       type = "response", re.form = NA)) * 100)
      risk_exp_str   <- sprintf("%.1f%%", mean(predict(ml, newdata = dat_exp,
                                                       type = "response", re.form = NA)) * 100)
    } else {
      risk_unexp_str <- NA_character_
      risk_exp_str   <- NA_character_
    }
    
    # --- Output tibbles (structure unchanged) ---
    ml_results <- tibble(
      model = "Multilevel Poisson", Outcome = outcome_var,
      Exposure = exposure_var, Profession = profession,
      Adjustment = model_label, beta = beta_ml, SE = se_ml,
      PR = pr_ml, CI_l = ci_ml[1], CI_u = ci_ml[2],
      N_Analysis = n_size, Risk_Unexp = risk_unexp_str, 
      Risk_Exp = risk_exp_str
    )
    
    rb_results <- tibble(
      model = "Poisson + RSVE", Outcome = outcome_var,
      Exposure = exposure_var, Profession = profession,
      Adjustment = model_label, beta = beta_rb, SE = se_rb,
      PR = pr_rb, CI_l = ci_rb[1], CI_u = ci_rb[2],
      N_Analysis = n_size, Risk_Unexp = risk_unexp_str,
      Risk_Exp = risk_exp_str
    )
    
    return(list(multilevel = ml_results, robust = rb_results))
  }

#---- Fit mediation model ----

run_mediation <- function(data, exposure_var, mediator_var, outcome_var,
                          confounder_list, profession = "GP",
                          cluster_var = "loc_2",
                          sims = 500,
                          output_dir = "out/models/mediation") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --- Abbreviation maps ---
  exp_short <- case_when(
    exposure_var == "wc_demands_z"   ~ "dem_z",
    exposure_var == "wc_hazards_z"   ~ "haz_z",
    exposure_var == "wc_resources_z" ~ "res_z",
    TRUE ~ gsub("[^a-z0-9]", "", tolower(exposure_var))
  )
  
  med_short <- case_when(
    mediator_var == "phqads_z" ~ "phqads",
    mediator_var == "phq_co"   ~ "dep",
    mediator_var == "gad_co"   ~ "anx",
    TRUE ~ gsub("[^a-z0-9]", "", tolower(mediator_var))
  )
  
  out_short <- case_when(
    outcome_var == "work_18_poi" ~ "turn",
    TRUE ~ gsub("[^a-z0-9]", "", tolower(outcome_var))
  )
  
  prof_short <- case_when(
    profession == "GP"     ~ "gp",
    profession == "Non-GP" ~ "nongp",
    profession == "All"    ~ "all",
    TRUE ~ gsub("[^a-z0-9]", "", tolower(profession))
  )
  
  file_name <- paste0("med_", exp_short, "_", med_short, "_", 
                      out_short, "_", prof_short, ".rds")
  file_path <- file.path(output_dir, file_name)
  
  message("Running: ", file_name)
  
  # --- Detect mediator type ---
  is_binary_mediator <- is.factor(data[[mediator_var]]) &&
    all(levels(data[[mediator_var]]) %in% c("No", "Yes"))
  
  is_binary_outcome <- is.factor(data[[outcome_var]]) ||
    all(na.omit(as.numeric(data[[outcome_var]])) %in% c(0, 1, 1, 2))
  
  # --- Clean data ---
  required_vars <- unique(c(outcome_var, exposure_var, mediator_var,
                            cluster_var, confounder_list))
  
  ds_clean <- data |>
    dplyr::select(dplyr::all_of(required_vars)) |>
    drop_na()
  
  n_size <- nrow(ds_clean)
  message("n: ", n_size)
  if (n_size < 10) return(NULL)
  
  # --- Recode outcome to 0/1 ---
  y_raw <- as.numeric(ds_clean[[outcome_var]])
  if (max(y_raw, na.rm = TRUE) > 1) ds_clean[[outcome_var]] <- y_raw - 1
  
  # --- Recode binary mediator to 0/1 numeric for mediation package ---
  if (is_binary_mediator) {
    ds_clean[[mediator_var]] <- as.integer(ds_clean[[mediator_var]] == "Yes")
  }
  
  # --- Build formulas ---
  confounders_str <- paste(confounder_list, collapse = " + ")
  
  f_med <- as.formula(paste0(
    mediator_var, " ~ ", exposure_var, " + ", confounders_str,
    " + (1 | ", cluster_var, ")"
  ))
  
  f_out <- as.formula(paste0(
    outcome_var, " ~ ", exposure_var, " + ", mediator_var, " + ",
    confounders_str, " + (1 | ", cluster_var, ")"
  ))
  
  # --- Fit mediator model ---
  fit_med <- if (is_binary_mediator) {
    glmer(
      f_med, data = ds_clean, family = binomial(link = "logit"),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  } else {
    lmer(
      f_med, data = ds_clean,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  }
  
  # --- Fit outcome model ---
  fit_out <- glmer(
    f_out, data = ds_clean, family = poisson(link = "log"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  
  # --- Run mediation ---
  med_result <- mediate(
    fit_med, fit_out,
    treat   = exposure_var,
    mediator = mediator_var,
    data    = ds_clean,
    sims    = sims,
    boot    = FALSE,
    boot.ci.type = "perc"
  )
  
  # --- Extract results ---
  out <- tibble(
    Exposure          = exposure_var,
    Mediator          = mediator_var,
    Outcome           = outcome_var,
    Profession        = profession,
    N                 = n_size,
    ACME              = med_result$d.avg,
    ACME_CI_l         = med_result$d.avg.ci[1],
    ACME_CI_u         = med_result$d.avg.ci[2],
    ACME_p            = med_result$d.avg.p,
    ADE               = med_result$z.avg,
    ADE_CI_l          = med_result$z.avg.ci[1],
    ADE_CI_u          = med_result$z.avg.ci[2],
    ADE_p             = med_result$z.avg.p,
    Total             = med_result$tau.coef,
    Total_CI_l        = med_result$tau.ci[1],
    Total_CI_u        = med_result$tau.ci[2],
    Total_p           = med_result$tau.p,
    Prop_mediated     = med_result$n.avg,
    Prop_med_CI_l     = med_result$n.avg.ci[1],
    Prop_med_CI_u     = med_result$n.avg.ci[2],
    Prop_med_p        = med_result$n.avg.p
  )
  
  saveRDS(list(result = out, med_object = med_result), file_path)
  
  return(out)
}