
# Pool association estimates from imputed dataset -------------------------

pool_rubins_rules <- 
  function(results_df) {
    
    # If the model failed and returned NULL, skip it
    if(nrow(results_df) == 0) return(NULL) 
    
    m <- nrow(results_df)                     # Number of imputations
    Q_bar <- mean(results_df$beta)            # Pooled Beta
    U_bar <- mean(results_df$SE^2)            # Within-imputation variance
    B <- var(results_df$beta)                 # Between-imputation variance
    T_var <- U_bar + (1 + 1/m) * B            # Total Variance
    SE_pool <- sqrt(T_var)                    # Pooled Standard Error
    
    # Grab metadata (names, exposures, risks) from the first imputation
    meta <- 
      
      results_df[1, ] |> 
      select(-beta, -SE, -PR, -CI_l, -CI_u)
    
    # Calculate the final pooled Prevalence Ratios and CIs
    pooled_stats <- 
      tibble(
        beta = Q_bar,
        SE = SE_pool,
        PR = exp(Q_bar),
        CI_l = exp(Q_bar - 1.96 * SE_pool),
        CI_u = exp(Q_bar + 1.96 * SE_pool)
      )
    
    return(bind_cols(meta, pooled_stats))
    
  }

# Run alternative models --------------------------------------------------

# Function for crude models

run_crude_mods <- 
  function(data, outcome_var, exposure_var, profession, 
           cluster_var, family_type = "poisson",
           output_dir = "out/models/main") {
    
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    out_short <- case_when(
      outcome_var == "phq_co_poi"      ~ "dep",
      outcome_var == "gad_co_poi"      ~ "anx",
      outcome_var == "suic_idea_poi"   ~ "suic",
      outcome_var == "cage_co_poi"     ~ "alc",
      outcome_var == "phq_sc"          ~ "dep_discr",
      outcome_var == "gad_sc"          ~ "anx_discr",
      outcome_var == "mh_phq_9"        ~ "suic_discr",
      outcome_var == "cage_sc"         ~ "alc_discr"
    )
    
    # exposure abbreviation
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
      exposure_var == "work_25_dic"   ~ "council",
      exposure_var == "work_26_dic"   ~ "feedb",
      exposure_var == "work_27_dic"   ~ "stressplan",
      exposure_var == "work_28_dic"   ~ "harprot",
      exposure_var == "work_29_dic"   ~ "vioprot",
      exposure_var == "work_30_dic"   ~ "har",
      exposure_var == "work_33_dic"   ~ "bul",
      exposure_var == "work_31_dic"   ~ "threats",
      exposure_var == "work_32_dic"   ~ "viol"
    )
    
    # profession abbreviation
    prof_short <- case_when(
      profession == "Doctor" ~ "doc",
      profession == "Nurse"  ~ "nur",
      profession == "All" ~ "all"
    )
    
    # stratify data (or not)
    ds_sub <- if (profession == "All") {
      data
    } else {
      data |> filter(work_2 == profession)
    }
    
    # Clean
    req_vars <- c(outcome_var, exposure_var, cluster_var)
    
    
    
    ds_clean <- ds_sub |> 
      select(all_of(req_vars)) |> 
      drop_na() 
    
    
    n_size <- nrow(ds_clean)
    
    # calculate absolute risks/means 
    y_vals <- as.numeric(ds_clean[[outcome_var]]) 
    x_vals <- as.character(ds_clean[[exposure_var]])
    
    if (family_type != "gaussian") {
      # ensure it is 0/1
      if(max(y_vals, na.rm = TRUE) > 1) y_vals <- y_vals - 1
      
      # calculations
      risk_unexp_val <- mean(y_vals[x_vals == "No"], na.rm = TRUE) * 100
      risk_exp_val   <- mean(y_vals[x_vals == "Yes"], na.rm = TRUE) * 100
      
      risk_unexp_str <- sprintf("%.1f%%", risk_unexp_val)
      risk_exp_str   <- sprintf("%.1f%%", risk_exp_val)
      
    } else {
      
      mean_unexp_val <- mean(y_vals[x_vals == "No"], na.rm = TRUE)
      mean_exp_val   <- mean(y_vals[x_vals == "Yes"], na.rm = TRUE)
      
      risk_unexp_str <- sprintf("%.2f", mean_unexp_val)
      risk_exp_str   <- sprintf("%.2f", mean_exp_val)
    }
    
    # Define family
    
    if (family_type == "gaussian") {
      curr_family <- gaussian()
    } else if (family_type == "log-binomial") {
      curr_family <- binomial(link = "log")
    } else if (family_type == "logistic") {
      curr_family <- binomial(link = "logit") # standard logistic
    } else {
      curr_family <- poisson(link = "log")
    }  
    
    # Define the estimate type to report
    estimate_type <- if(family_type == "logistic") "OR" else "PR"
    
    
    # Define names
    file_name_ml <- paste0("glm_", out_short, "_", exp_short, "_", prof_short, 
                           "_cr_", family_type, "_ml", ".rds")
    
    file_path_ml <- file.path(output_dir, file_name_ml)
    
    ml <- NULL
    
    # Run or load
    if(file.exists(file_path_ml)) {
      
      message(glue::glue("Loading ML model: {file_name_ml}"))
      ml <- readRDS(file_path_ml)
      
    } else {
      message(glue::glue("Fitting ML model: {file_name_ml}"))
      
      
      f_ml <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var, 
               " + (1 | ", cluster_var, ")"))
      
      ml <- tryCatch({
        glmer(
          f_ml, 
          data = ds_clean, 
          family = curr_family,
          control = glmerControl(
            optimizer = "bobyqa", 
            optCtrl = list(maxfun = 2e5)
          )
        )}, error = function(e) {
          
          message(glue::glue("Convergence failure in {file_name_ml}: {e$message}"))
          return(NULL)
        })
      
      
      if(!is.null(ml)) saveRDS(ml, file_path_ml)
    }
    
    
    # Names
    file_name_rob <-  paste0("glm_", out_short, "_", exp_short, "_", prof_short, 
                             "_cr_", family_type, "_rob", ".rds")
    
    file_path_rob <- file.path(output_dir, file_name_rob)
    
    rob <- NULL
    
    # Run or load
    if(file.exists(file_path_rob)) {
      
      message(glue::glue("Loading Robust model: {file_name_rob}"))
      rob <- readRDS(file_path_rob)
      
    } else {
      message(glue::glue("Fitting Robust model: {file_name_rob}"))
      
      f_rob <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var))
      
      rob <- tryCatch({
        glm(
          f_rob, 
          data = ds_clean, 
          family = curr_family)
      }, error = function(e) {
        
        message(glue::glue("Convergence failure in {file_name_rob}: {e$message}"))
        return(NULL)
      })
      
      if(!is.null(rob)) saveRDS(rob, file_path_rob)
      
    }
    
    ml_results <- NULL
    rob_results <- NULL
    
    # Conditional extraction if model converged
    
    if(!is.null(ml)) {
      coef_ml <- summary(ml)$coefficients
      
      beta_ml <- coef_ml[paste0(exposure_var, "Yes"), "Estimate"]
      
      se_ml   <- coef_ml[paste0(exposure_var, "Yes"), "Std. Error"]
      
      pr_ml <-  exp(beta_ml)
      
      ci_ml <- exp(beta_ml + c(-1, 1) * 1.96 * se_ml)
      
      ml_results <- tibble(
        
        model = "Multilevel",
        
        Family_Type = family_type,
        
        Measure = estimate_type,
        
        Outcome = outcome_var,
        
        Exposure = exposure_var,
        
        Profession = profession,
        
        Adjustment = "Crude",
        
        beta  = beta_ml,
        
        SE    = se_ml,
        
        PR    = pr_ml,
        
        CI_l  = ci_ml[1],
        
        CI_u  = ci_ml[2],
        
        N_Analysis = n_size,
        
        Risk_Unexp = risk_unexp_str,
        
        Risk_Exp = risk_exp_str
        
      )
    }
    
    if(!is.null(rob)) {
      
      vcov_rob <- sandwich::vcovCL(rob, cluster = ds_clean[[cluster_var]])
      
      coef_glm <- summary(rob)$coefficients
      
      beta_rb <- coef_glm[paste0(exposure_var, "Yes"), "Estimate"]
      
      se_rb <- sqrt(vcov_rob[paste0(exposure_var, "Yes"), 
                             paste0(exposure_var, "Yes")])
      
      pr_rb <- exp(beta_rb)
      
      ci_rb <- exp(beta_rb + c(-1, 1) * 1.96 * se_rb)
      
      rob_results <- tibble(
        
        model = "Robust SE",
        
        Family_Type = family_type,
        
        Measure = estimate_type,
        
        Outcome = outcome_var,
        
        Exposure = exposure_var,
        
        Profession = profession,
        
        Adjustment = "Crude",
        
        beta  = beta_rb,
        
        SE    = se_rb,
        
        PR    = pr_rb,
        
        CI_l  = ci_rb[1],
        
        CI_u  = ci_rb[2],
        
        N_Analysis = n_size 
      )
    }
    
    return(list(
      multilevel = ml_results,
      robust     = rob_results
    ))
    
  }



# Function for adjusted models

run_adj_mods <- 
  function(data, outcome_var, exposure_var, confounder_list, profession,
           cluster_var = "loc_2", family_type = "poisson", model_label,
           output_dir = "out/models/main") {
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    
    # stratify data (or not)
    ds_sub <- if (profession == "All") {
      data
    } else {
      data |> filter(work_2 == profession)
    }
    
    # remove work_2 cause it is a constant
    active_confounders <- if (profession == "All") {
      confounder_list
    } else {
      setdiff(confounder_list, "work_2")
    }
    
    # define clean set
    required_vars <- 
      unique(c(outcome_var, exposure_var, cluster_var, active_confounders))
    
    out_short <- case_when(
      outcome_var == "phq_co_poi"      ~ "dep",
      outcome_var == "gad_co_poi"      ~ "anx",
      outcome_var == "suic_idea_poi"   ~ "suic",
      outcome_var == "cage_co_poi"     ~ "alc"
    )
    
    # exposure abbreviation
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
      exposure_var == "work_25_dic"   ~ "council",
      exposure_var == "work_26_dic"   ~ "feedb",
      exposure_var == "work_27_dic"   ~ "stressplan",
      exposure_var == "work_28_dic"   ~ "harprot",
      exposure_var == "work_29_dic"   ~ "vioprot",
      exposure_var == "work_30_dic"   ~ "har",
      exposure_var == "work_33_dic"   ~ "bul",
      exposure_var == "work_31_dic"   ~ "threats",
      exposure_var == "work_32_dic"   ~ "viol"
    )
    
    
    
    # profession abbreviation
    prof_short <- case_when(
      profession == "Doctor" ~ "doc",
      profession == "Nurse"  ~ "nur",
      profession == "All" ~ "all"
    )
    
    
    
    ds_clean <- ds_sub |> 
      select(all_of(required_vars)) |> 
      drop_na()
    
    
    n_size <- nrow(ds_clean)
    
    if(n_size < 10) return(NULL)
    
    # Convert Exposure to Factor for prediction consistency
    ds_clean[[exposure_var]] <- 
      
      factor(ds_clean[[exposure_var]], 
             levels = c("No", "Yes"))
    
    # Convert Outcome to 0/1 for Poisson/Logit
    y_raw <- as.numeric(ds_clean[[outcome_var]])
    
    if(max(y_raw, na.rm = TRUE) > 1) {
      ds_clean[[outcome_var]] <- y_raw - 1
    }
    
    # Define family 
    if (family_type == "gaussian") {
      curr_family <- gaussian()
    } else if (family_type == "log-binomial") {
      curr_family <- binomial(link = "log")
    } else if (family_type == "logistic") {
      curr_family <- binomial(link = "logit") # standard logistic
    } else {
      curr_family <- poisson(link = "log")
    }  
    
    # Estimate type
    
    estimate_type <- if(family_type == "logistic") "OR" else "PR"
    
    
    # Naming logic for each adjustment
    
    model_suffix <- case_when(
      model_label == "Fully adjusted" ~ "_adj",
      TRUE ~ tolower(gsub(" ", "_", model_label))
    )
    
    
    # Define names
    file_name_ml <- 
      
      paste0(
        "glm_", 
        out_short, 
        "_", 
        exp_short , 
        "_", 
        prof_short, 
        "_", 
        model_suffix, 
        "_",
        family_type,
        "_ml.rds"
      )
    
    file_path_ml <- file.path(output_dir, file_name_ml)
    
    ml <- NULL 
    
    # run or load
    if(file.exists(file_path_ml)) {
      
      message(glue::glue("Loading ML model: {file_name_ml}"))
      ml <- tryCatch({
        readRDS(file_path_ml)
      }, error = function(e) {
        
        message(glue::glue("Corrupted file detected! Deleting {file_name_ml} and refitting..."))
        file.remove(file_path_ml)
        return(NULL)
      })
    } 
    
    if(is.null(ml)) {
      message(glue::glue("Fitting ML model: {file_name_ml}"))
      
      f_ml <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var, "+", 
               paste(active_confounders, collapse = "+"),
               " + (1 | ", cluster_var, ")"))
      
      ml <- tryCatch({
        glmer(
          f_ml, 
          data = ds_clean, 
          family = curr_family,
          control = glmerControl(
            optimizer = "bobyqa", 
            optCtrl = list(maxfun = 2e5)
            )
          )}, error = function(e) {
            
            message(glue::glue("Convergence failure in {file_name_ml}: {e$message}"))
            return(NULL)
            
            })
      
        if(!is.null(ml)) saveRDS(ml, file_path_ml)
    }
    
    
    # Define names
    file_name_rob <- paste0("glm_", out_short, "_", exp_short , "_", prof_short,
                            "_", model_suffix, "_", family_type, "_rob.rds")
    file_path_rob <- file.path(output_dir, file_name_rob)
    
    
    rob <- NULL 
    
    # run or load
    if(file.exists(file_path_rob)) {
      
      message(glue::glue("Loading Robust SE model: {file_name_rob}"))
      rob <- tryCatch({
        readRDS(file_path_rob)
      }, error = function(e) {
        
        message(glue::glue("Corrupted file detected! Deleting {file_name_rob} and refitting..."))
        file.remove(file_path_rob)
        return(NULL)
      })
    } 
    
    if(is.null(rob)) {
      message(glue::glue("Fitting Robust model: {file_name_rob}"))
      
      # 1. Removed the random effect "(1 | cluster_var)" from the formula
      f_rob <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var, "+", 
               paste(active_confounders, collapse = "+")))
      
      # 2. Changed glmer() back to base glm(), and removed the glmerControl arguments
      
      rob <- tryCatch({
        glm(
          f_rob, 
          data = ds_clean, 
          family = curr_family
          )
      }, error = function(e) {
        
        message(glue::glue("Convergence failure in {file_name_rob}: {e$message}"))
        return(NULL)
      
        })
      
      if(!is.null(rob)) saveRDS(rob, file_path_rob)
    }
      
    ml_results <- NULL
    rob_results <-  NULL
    
    # Calculate adjusted risks (Marginal standardization)
    # Use ML as it accounts for the cluster variance structure
    if(!is.null(ml)) {
      
      # synthetic ds
      dat_unexp <- ds_clean
      dat_unexp[[exposure_var]] <- factor("No", levels = c("No", "Yes"))
      
      dat_exp <- ds_clean
      dat_exp[[exposure_var]] <- factor("Yes", levels = c("No", "Yes"))
      
      # Predict: Average predicted probability if everyone UNEXPOSED
      pred_unexp <- 
        
        predict(
          ml, 
          newdata = dat_unexp, 
          type = "response", 
          re.form = NA
        )
      
      
      pred_exp   <- 
        
        predict(
          ml, 
          newdata = dat_exp,   
          type = "response", 
          re.form = NA
        )
      
      risk_unexp_str <- sprintf("%.1f%%", mean(pred_unexp, na.rm = TRUE) * 100)
      
      risk_exp_str   <- sprintf("%.1f%%", mean(pred_exp,   na.rm = TRUE) * 100)
      
      # Extract for ml
      coef_ml <- summary(ml)$coefficients
      
      beta_ml <- coef_ml[paste0(exposure_var, "Yes"), "Estimate"]
      
      se_ml   <- coef_ml[paste0(exposure_var, "Yes"), "Std. Error"]
      
      pr_ml <-  exp(beta_ml)
      
      ci_ml <- exp(beta_ml + c(-1, 1) * 1.96 * se_ml)
      
      # Generate data
      
      ml_results <- tibble(
        
        model = "Multilevel",
        
        Family_Type = family_type,
        
        Measure = estimate_type,
        
        Outcome = outcome_var,
        
        Exposure = exposure_var,
        
        Profession = profession,
        
        Adjustment = model_label,
        
        beta  = beta_ml,
        
        SE    = se_ml,
        
        PR    = pr_ml,
        
        CI_l  = ci_ml[1],
        
        CI_u  = ci_ml[2],
        
        N_Analysis = n_size, 
        
        Risk_Unexp = risk_unexp_str,
        
        Risk_Exp = risk_exp_str
        
      )
    }
    
    
    # Extract for robust SE
    if(!is.null(rob)){
      
      vcov_rob <- sandwich::vcovCL(rob, cluster = ds_clean[[cluster_var]])
    
      coef_glm <- summary(rob)$coefficients
      
      beta_rb <- coef_glm[paste0(exposure_var, "Yes"), "Estimate"]
      
      se_rb <- sqrt(vcov_rob[paste0(exposure_var, "Yes"), 
                             paste0(exposure_var, "Yes")])
      
      pr_rb <- exp(beta_rb)
      
      ci_rb <- exp(beta_rb + c(-1, 1) * 1.96 * se_rb)
      
      rob_results <- tibble(
        
        model = "Robust SE",
        
        Family_Type = family_type,
        
        Measure = estimate_type,
        
        Outcome = outcome_var,
        
        Exposure = exposure_var,
        
        Profession = profession,
        
        Adjustment = model_label,
        
        beta  = beta_rb,
        
        SE    = se_rb,
        
        PR    = pr_rb,
        
        CI_l  = ci_rb[1],
        
        CI_u  = ci_rb[2],
        
        N_Analysis = n_size # number of observations used
        )
      }
    
    return(list(
      multilevel = ml_results,
      robust     = rob_results
    ))
  }