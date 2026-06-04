# Unweighted associations -------------------------------------------------

# Function for GLM crude models

run_crude_glm <- 
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
    } else {
      curr_family <- poisson(link = "log")
    }
    
    # Define names
    file_name_ml <- paste0("glm_", out_short, "_", exp_short, "_", prof_short, 
                           "_cr_ml", ".rds")
    
    file_path_ml <- file.path(output_dir, file_name_ml)
    
    # Run or load
    if(file.exists(file_path_ml)) {
      
      message(glue::glue("Loading ML model: {file_name_ml}"))
      ml <- readRDS(file_path_ml)
      
    } else {
      message(glue::glue("Fitting ML model: {file_name_ml}"))
      
      
      f_ml <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var, 
               " + (1 | ", cluster_var, ")"))
      
      ml <- glmer(
        f_ml, 
        data = ds_clean, 
        family = curr_family,
        control = glmerControl(
          optimizer = "bobyqa", 
          optCtrl = list(maxfun = 2e5)
        )
      )
      
      
      saveRDS(ml, file_path_ml)
    }
    
    
    # Names
    file_name_rob <-  paste0("glm_", out_short, "_", exp_short, "_", prof_short, 
                             "_cr_rob", ".rds")
    
    file_path_rob <- file.path(output_dir, file_name_rob)
    
    # Run or load
    if(file.exists(file_path_rob)) {
      
      message(glue::glue("Loading Robust model: {file_name_rob}"))
      rob <- readRDS(file_path_rob)
      
    } else {
      message(glue::glue("Fitting Robust model: {file_name_rob}"))
      
      f_rob <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var))
      
      rob <- glm(
        f_rob, 
        data = ds_clean, 
        family = curr_family)
      
      
      saveRDS(rob, file_path_rob)
      
    }
    
    
    # Extract for ml
    coef_ml <- summary(ml)$coefficients
    
    beta_ml <- coef_ml[paste0(exposure_var, "Yes"), "Estimate"]
    
    se_ml   <- coef_ml[paste0(exposure_var, "Yes"), "Std. Error"]
    
    pr_ml <-  exp(beta_ml)
    
    ci_ml <- exp(beta_ml + c(-1, 1) * 1.96 * se_ml)
    
    # Generate data
    
    ml_results <- tibble(
      
      model = "Multilevel Poisson",
      
      Outcome = outcome_var,
      
      Exposure = exposure_var,
      
      Profession = profession,
      
      Adjustment = "Crude",
      
      beta  = beta_ml,
      
      SE    = se_ml,
      
      PR    = pr_ml,
      
      CI_l  = ci_ml[1],
      
      CI_u  = ci_ml[2],
      
      N_Analysis = n_size, # number of observations used
      
      Risk_Unexp = risk_unexp_str,
      
      Risk_Exp = risk_exp_str
      
    )
    
    
    # Extract for robust SE
    
    vcov_rob <- sandwich::vcovCL(rob, cluster = ds_clean[[cluster_var]])
    
    coef_glm <- summary(rob)$coefficients
    
    beta_rb <- coef_glm[paste0(exposure_var, "Yes"), "Estimate"]
    
    se_rb <- sqrt(vcov_rob[paste0(exposure_var, "Yes"), 
                           paste0(exposure_var, "Yes")])
    
    pr_rb <- exp(beta_rb)
    
    ci_rb <- exp(beta_rb + c(-1, 1) * 1.96 * se_rb)
    
    rb_results <- tibble(
      
      model = "Poisson + RSVE",
      
      Outcome = outcome_var,
      
      Exposure = exposure_var,
      
      Profession = profession,
      
      Adjustment = "Crude",
      
      beta  = beta_rb,
      
      SE    = se_rb,
      
      PR    = pr_rb,
      
      CI_l  = ci_rb[1],
      
      CI_u  = ci_rb[2],
      
      N_Analysis = n_size # number of observations used
    )
    
    return(list(
      multilevel = ml_results,
      robust     = rb_results
    ))
    
  }



# Function for adjusted models

run_adj_glm <- 
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
    } else {
      curr_family <- poisson(link = "log")
    }
    
    
    # Naming logic for each adjustment
    
    model_suffix <- case_when(
      model_label == "Partially adjusted" ~ "_str",
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
      
      ml <- glmer(
        f_ml, 
        data = ds_clean, 
        family = curr_family,
        control = glmerControl(
          optimizer = "bobyqa", 
          optCtrl = list(maxfun = 2e5)
        )
      )
      
      saveRDS(ml, file_path_ml)
    }
    
    
    # Define names
    file_name_rob <- paste0("glm_", out_short, "_", exp_short , "_", prof_short,
                            "_", model_suffix, "_rob.rds")
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
      rob <- glm(
        f_rob, 
        data = ds_clean, 
        family = curr_family
      )
      
      saveRDS(rob, file_path_rob)
    }
    
    # Calculate adjusted risks (Marginal standardization)
    # Use ML as it accounts for the cluster variance structure
    
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
      
      model = "Multilevel Poisson",
      
      Outcome = outcome_var,
      
      Exposure = exposure_var,
      
      Profession = profession,
      
      Adjustment = model_label,
      
      beta  = beta_ml,
      
      SE    = se_ml,
      
      PR    = pr_ml,
      
      CI_l  = ci_ml[1],
      
      CI_u  = ci_ml[2],
      
      N_Analysis = n_size, # number of observations used
      
      Risk_Unexp = risk_unexp_str,
      
      Risk_Exp = risk_exp_str
      
    )
    
    
    # Extract for robust SE
    
    vcov_rob <- sandwich::vcovCL(rob, cluster = ds_clean[[cluster_var]])
    
    coef_glm <- summary(rob)$coefficients
    
    beta_rb <- coef_glm[paste0(exposure_var, "Yes"), "Estimate"]
    
    se_rb <- sqrt(vcov_rob[paste0(exposure_var, "Yes"), 
                           paste0(exposure_var, "Yes")])
    
    pr_rb <- exp(beta_rb)
    
    ci_rb <- exp(beta_rb + c(-1, 1) * 1.96 * se_rb)
    
    rb_results <- tibble(
      
      model = "Poisson + RSVE",
      
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
    
    return(list(
      multilevel = ml_results,
      robust     = rb_results
    ))
  }


# Generate results table -------------------------------------------------------

make_stratified_table <- 
  function(data, title_text) {
    
    # Pivot to Wide Format
    # We ensure factors are respected so columns appear in the correct order (Crude -> Partial -> Full)
    wide_df <- 
      
      data |>
      # Use existing factors if available, or set them here to ensure order
      mutate(
        Adjustment = factor(Adjustment, levels = c("Crude", "Adjusted")),
        Profession = factor(Profession, levels = c("Doctor", "Nurse")),
        Exposure = factor(Exposure, levels = exposure_order_named),
        Outcome = factor(Outcome, levels = c("Depression", "Anxiety", "Alcohol dependence", "Suicidal thoughts"))
      ) |>
      arrange(Outcome, Exposure) |>
      select(Outcome, Exposure, Profession, Adjustment, Est_CI) |>
      pivot_wider(
        names_from = c(Profession, Adjustment),
        values_from = Est_CI,
        names_sep = "_"
      )
    
    # Create Grouping Rows
    # This inserts "Header Rows" for each Outcome (Depression, Anxiety, etc.)
    grouped_df <- as_grouped_data(wide_df, groups = "Outcome")
    
    # For the header rows (where Outcome is NOT NA), copy the text into the Exposure column.
    # We convert to character first to avoid Factor level errors.
    grouped_df$Exposure <- as.character(grouped_df$Exposure)
    grouped_df$Outcome  <- as.character(grouped_df$Outcome)
    
    grouped_df$Exposure <- ifelse(
      !is.na(grouped_df$Outcome), # If this is a header row...
      grouped_df$Outcome,         # ...use the Outcome name (e.g., "Depression")
      grouped_df$Exposure         # ...otherwise keep the Exposure name
    )
    
    # Build Flextable
    ft <- 
      
      flextable(grouped_df, 
                col_keys = c("Exposure", "Doctor_Crude", "Doctor_Adjusted",
                             "Nurse_Crude", "Nurse_Adjusted")) |>
      
      # Headers 
      # Rename the raw columns (e.g., "Doctor_Crude") to clean names ("Crude")
      set_header_labels(
        Exposure = "Exposure",
        `Doctor_Crude` = "Crude",
        `Doctor_Adjusted` = "Adjusted",
        `Nurse_Crude` = "Crude",
        `Nurse_Adjusted` = "Adjusted"
      ) |>
      # Add Spanner (Doctor / Nurse)
      add_header_row(
        values = c("", "Doctor", "Nurse"), 
        colwidths = c(1, 2, 2) # 1 col for Exposure, 2 for Doctor, 2 for Nurse
      ) |>
      
      # Formatting 
      theme_booktabs() |> # Clean scientific style
      autofit() |>
      align(align = "center", part = "all") |> # Center everything
      align(j = 1, align = "left", part = "all") |> # Left align the Exposure column
      
      # Group Row Styling (The "Outcome" headers) 
      # Make the grouping rows (Depression, Anxiety) bold and spanning all columns
      bold(j = 1, i = ~ !is.na(Outcome)) |>
      bold(part = "header") |>
      
      # Indentation for Exposures
      # Indent rows that are NOT group headers (where Outcome is NA)
      padding(j = 1, i = ~ is.na(Outcome), padding.left = 2) |>
      
      # Final Touches 
      add_header_lines(values = title_text) |>
      add_footer_lines("Prevalence Ratios (95% CI).")
    
    return(ft)
  }


# Dose response effects ---------------------------------------------------


run_glm_dose_emm <- 
  function(
    data, outcome_var, exposure_var, confounder_list,
    cluster_var = "loc_2", family_type = "gaussian",
    model_label, output_dir = "out/models/dose",
    force_fit = FALSE) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  required_vars <- unique(c(outcome_var, exposure_var, cluster_var, confounder_list))
  
  # Clean dataset
  ds_clean <- 
    
    data |> 
    select(all_of(required_vars)) |> 
    drop_na() 
  
  n_size <- nrow(ds_clean)
  if (n_size < 10) return(NULL)
  
  
  # Create categorical exposure variable
  exposure_cat <- paste0(exposure_var, "_cat")
  ds_clean[[exposure_cat]] <- as.factor(ds_clean[[exposure_var]])
  
  
  # Short labels
  
  out_short <- case_when(
    outcome_var == "phq_sc"   ~ "dep",
    outcome_var == "gad_sc"   ~ "anx",
    outcome_var == "mh_phq_9" ~ "suic",
    outcome_var == "cage_sc"  ~ "alc"
  )
  
  exp_short <- case_when(
    
    # hazards
    exposure_var == "work_30_num" ~ "har",
    exposure_var == "work_33_num" ~ "bul",
    exposure_var == "work_31_num" ~ "threats",
    exposure_var == "work_32_num" ~ "viol",
    
    # demands
    
    exposure_var == "work_8_num"  ~ "hrs",
    exposure_var == "work_11_num" ~ "night",
    exposure_var == "work_12_num" ~ "shift",
    exposure_var == "work_23_num" ~ "deadl",
    exposure_var == "work_24_num" ~ "angry",
    
    # resources
    
    exposure_var == "work_19_num" ~ "influ",
    exposure_var == "work_20_num" ~ "break",
    exposure_var == "work_21_num" ~ "colsupp",
    exposure_var == "work_22_num" ~ "supsupp",
    
    # protocols
    
    exposure_var == "work_25_num" ~ "union",
    exposure_var == "work_26_num" ~ "feedb",
    exposure_var == "work_27_num" ~ "stress",
    exposure_var == "work_28_num" ~ "bullprot",
    exposure_var == "work_29_num" ~ "violprot"
  )
  
  file_name <- paste0("dose_", out_short, "_", exp_short, "_",
                      tolower(gsub(" ", "_", model_label)), ".rds")
  file_path <- file.path(output_dir, file_name)
  
  # Fit or load model
  if(file.exists(file_path) & !force_fit){
    message(glue::glue("Loading saved bundle: {file_name}"))
    
    bundle <- readRDS(file_path)
    model_cont <- bundle$model_cont
    model_cat <- bundle$model_cat
  
    } else {
      message(glue::glue("Fitting Dose-response: {file_name}"))
      
      if (model_label == "Crude" || is.null(confounder_list)) {
        covariates_str <- ""
      } else {
        covariates_str <- paste0(" + ", paste(confounder_list, collapse = " + "))
      }
      
      form_cont <- as.formula(paste(outcome_var, "~", exposure_var, 
                                    covariates_str, "+ (1 |", cluster_var, ")"))
      form_cat  <- as.formula(paste(outcome_var, "~", exposure_cat, 
                                    covariates_str, "+ (1 |", cluster_var, ")"))
        
      if (family_type == "gaussian") {
        
        model_cont <- tryCatch({
          lmer(
            form_cont, 
            data = ds_clean, 
            control = lmerControl(optimizer = "bobyqa"))
        }, error = function(e) { message("Error cont: ", e$message); return(NULL) })
        
        model_cat <- tryCatch({
          lmer(
            form_cat, 
            data = ds_clean, 
            control = lmerControl(optimizer = "bobyqa"))
        }, error = function(e) { message("Error cat: ", e$message); return(NULL) })
        
      } else {
          
        model_cont <- tryCatch({
          glmer(
            form_cont, 
            data = ds_clean, 
            family = family_type, 
            control = glmerControl(optimizer = "bobyqa"))
        }, error = function(e) { message("Error cont: ", e$message); return(NULL) })
        
        model_cat <- tryCatch({
          glmer(
            form_cat, 
            data = ds_clean, 
            family = family_type, 
            control = glmerControl(optimizer = "bobyqa"))
        }, error = function(e) { message("Error cat: ", e$message); return(NULL) })
      }
      
      if(is.null(model_cont) || is.null(model_cat)) return(NULL)
      
      saveRDS(list(model_cont = model_cont, model_cat = model_cat), file_path)
      }
    
  # Extract continuous EMM
  unique_doses <- sort(unique(ds_clean[[exposure_var]]))
  at_list <- setNames(list(unique_doses), exposure_var)
  
  # Extract EMMs
  
  emm_cont_obj <- tryCatch({
    emmeans(
      model_cont,
      specs = exposure_var,
      type = "response",
      nuisance = confounder_list,
      at = at_list, 
     df = NULL)
     }, error = function(e) NULL)
  
  if (!is.null(emm_cont_obj)) {
    
    df_cont <- as.data.frame(emm_cont_obj) |> 
      rename(
        estimate = emmean,
        conf.low = any_of(c("asymp.LCL", "lower.CL")),
        conf.high = any_of(c("asymp.UCL", "upper.CL"))
      ) |> 
      mutate(Level = as.character(.data[[exposure_var]]), Fit_type = "Continuous")
    
    } else {
      
      df_cont <- tibble()
      
      }
  
  # Extract categorical EMM
  
  emm_cat_obj <- tryCatch({
    emmeans(
      model_cat,
      specs = exposure_cat,
      type = "response",
      nuisance = confounder_list,
      at = at_list, 
      df = NULL)
    }, error = function(e) NULL)
  
  if (!is.null(emm_cat_obj)) {
    
    df_cat <- as.data.frame(emm_cat_obj) |> 
      rename(
        estimate = emmean,
        conf.low = any_of(c("asymp.LCL", "lower.CL")),
        conf.high = any_of(c("asymp.UCL", "upper.CL"))
      ) |> 
      mutate(Level = as.character(.data[[exposure_cat]]), Fit_type = "Categorical")
  
    } else {
      
      df_cat <- tibble()
      
      }
  
  # Merge
  emm_df <- bind_rows(df_cont, df_cat)
  
  if(nrow(emm_df) == 0){
    message(glue::glue("Warning: No EMMs found for {exposure_var}"))
    return(NULL)
  }
  
  emm_df <- emm_df |> 
    mutate(
      Outcome = outcome_var,
      Exposure = exposure_var,
      Model = model_label,
      N_Analysis = n_size
    ) |> 
    select(Outcome, Exposure, Fit_type, Level, estimate, conf.low, conf.high, Model, N_Analysis)

  return(emm_df)
}

# Survey models------------------------------------------------------------
# Run and save every survey model

run_svy_models <- 
  function(outcome, exposure, profession, design_object, 
           output_dir = "out/svy_models") {
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # outcome abbreviation
    out_short <- case_when(
      outcome == "phq_co_poi"      ~ "dep",
      outcome == "gad_co_poi"      ~ "anx",
      outcome == "suic_idea_poi"   ~ "suic"
    )
    
    # exposure abbreviation
    exp_short <- case_when( 
      exposure == "work_30_dic"   ~ "har",
      exposure == "work_33_dic"   ~ "bul",
      exposure == "work_31_dic"   ~ "threats",
      exposure == "work_32_dic"   ~ "viol"
    )
    
    # profession abbreviation
    prof_short <- case_when(
      profession == "Doctor" ~ "doc",
      profession == "Nurse"  ~ "nur"
    )
    
    # determine model type
    if (outcome == "wb_who_sc") {
      current_family <- gaussian()
    } else {
      current_family <- quasipoisson(link = "log")
    }
    
    # This line centers or adjust in the grand mean of neighboring strata
    # when any of these has few observations
    
    options(survey.lonely.psu = "adjust")
    
    # subset design
    sub_design <- subset(design_object, work_2 == profession)
    
    
    # build formulas
    f_crude <- as.formula(paste(outcome, "~", exposure))
    
    f_adj <- as.formula(paste(outcome, "~", exposure, "+",
                              paste(covs, collapse = "+")))
    
    # run
    
    m_adj <- try(svyglm(
      f_adj,
      design = sub_design,
      family = current_family
    ))
    
    
    # construct systematic names
    base_name <- paste("svy", out_short, exp_short, prof_short, sep = "_")
    
    name_crude <- paste0(base_name, "_cr")
    file_crude <- file.path(output_dir, paste0(name_crude, ".rds"))
    
    if (file.exists(file_crude)) {
      
      message(glue::glue("Loading from disk: {name_crude}"))
      m_crude <- readRDS(file_crude)
      
    } else {
      
      message(glue::glue("Fitting new model: {name_crude}"))
      
      m_crude <- try(svyglm(
        f_crude,
        design = sub_design,
        family = current_family
      ))
      
      if (!inherits(m_crude, "try-error")) {
        saveRDS(m_crude, file_crude)
      }
    }
    
    assign(name_crude, m_crude, envir = .GlobalEnv)
    
    
    name_adj <- paste0(base_name, "_adj")
    file_adj <- file.path(output_dir, paste0(name_adj, ".rds"))
    
    if (file.exists(file_adj)) {
      
      message(glue::glue("Loading from disk: {name_adj}"))
      m_adj <- readRDS(file_adj)
      
    } else {
      
      message(glue::glue("Fitting new model: {name_adj}"))
      
      m_adj <- try(svyglm(
        f_adj, 
        design = sub_design, 
        family = current_family
      ))
      
      if (!inherits(m_adj, "try-error")) {
        saveRDS(m_adj, file_adj)
      }
    }
    
    assign(name_adj, m_adj, envir = .GlobalEnv)
    
    
    return(invisible(TRUE))
  }


# Moderation effects ---------------------------------------------------

run_glm_moder_emm <- function(
    data, outcome_var, exposure_var, confounder_list, moderator_var,
    cluster_var = "loc_2", family_type = "gaussian",
    model_label, output_dir = "out/models/inter",
    force_fit = FALSE
) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  required_vars <- 
    unique(
      c(outcome_var, exposure_var, cluster_var, confounder_list, moderator_var)
    )
  
  # Clean dataset
  ds_clean <- data |> 
    select(all_of(required_vars)) # |> 
    # drop_na() |> 
    # mutate(
    #   !!exposure_var := factor(
    #     .data[[exposure_var]],
    #     levels = c("No", "Yes, a few times", "Yes, monthly", "Yes, weekly", "Yes, daily"),
    #     ordered = TRUE
    #   )
    # )
  
  n_size <- nrow(ds_clean)
  if (n_size < 10) return(NULL)
  
  # Short labels
  out_short <- case_when(
    outcome_var == "phq_sc"   ~ "dep",
    outcome_var == "gad_sc"   ~ "anx",
    outcome_var == "mh_phq_9" ~ "suic",
    outcome_var == "cage_sc"  ~ "alc"
  )
  
  exp_short <- case_when(
    exposure_var == "work_30" ~ "har",
    exposure_var == "work_33" ~ "bul",
    exposure_var == "work_31" ~ "threats",
    exposure_var == "work_32" ~ "viol",
    exposure_var == "work_8"  ~ "hr",
    exposure_var == "work_11" ~ "night"
  )
  
  moder_short <- case_when(
    moderator_var == "work_28_dic" ~ "bull_harass",
    moderator_var == "work_29_dic" ~ "threats_abuse_assault",
    moderator_var == "work_16"     ~ "contract"
    
  )
  
  file_name <- 
    
    paste0(
      "moder_", 
      out_short, "_", 
      exp_short, "_", 
      moder_short, "_", 
      tolower(gsub(" ", "_", model_label)), 
      ".rds")
  
  file_path <- file.path(output_dir, file_name)
  
  # Fit or load model
  
  if(file.exists(file_path) & !force_fit){
    message(glue::glue("Loading saved bundle: {file_name}"))
    bundle <- readRDS(file_path)
    model <- bundle$model
    balance_report <- bundle$balance
  } else {
    message(glue::glue("Fitting moderation: {file_name}"))
    
    balance_report <- NULL
    
    model_formula <- as.formula(
      paste(
        outcome_var, "~",
        exposure_var, "*", moderator_var, "+",
        paste(confounder_list, collapse = " + "),
        "+ (1 |", cluster_var, ")"
      )
    )
    
    model <- lmer(
      model_formula,
      data = ds_clean,
      control = lmerControl(optimizer = "bobyqa")
    )
    
    saveRDS(list(model = model, balance = balance_report), file_path)
  }
  
  # Extract EMMs
  
  emm_obj <- emmeans(
    model,
    specs = c(exposure_var, moderator_var),
    type = "response",
    nuisance = confounder_list, # required for adj models with several covs
    df = NULL
  )
  
  emm_df <- as.data.frame(
    emm_obj
  ) %>%
    rename(estimate = emmean, std.error = SE) %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      Level = as.character(.data[[exposure_var]]),
      Moderator_Level = as.character(.data[[moderator_var]]),
      Outcome = outcome_var,
      Exposure = exposure_var,
      Outcome = outcome_var,
      Exposure = exposure_var,
      Moderator = moderator_var,
      Model = model_label,
      N_Analysis = n_size
    ) %>%
    select(Outcome, Exposure, Moderator, Moderator_Level, Level, estimate, conf.low, conf.high, Model, N_Analysis)
  
  if(is.null(emm_df) || nrow(emm_df) == 0){
    message(glue::glue("Warning: No EMMs found for {exposure_var}"))
    return(NULL)
  }
  return(emm_df)
}

run_glm_moder_lrt <- function(
    data,
    outcome_var,
    exposure_var,
    moderator_var,
    confounder_list,
    cluster_var = "loc_2",
    family_type = "gaussian",
    model_label
) {
  
  # Variables necesarias
  required_vars <- unique(
    c(outcome_var, exposure_var, moderator_var,
      confounder_list, cluster_var)
  )
  
  # Clean dataset (MISMA lógica que tu función EMM)
  ds_clean <- data |>
    select(all_of(required_vars)) |>
    drop_na() |>
    mutate(
      !!exposure_var := factor(
        .data[[exposure_var]],
        levels = c(
          "No",
          "Yes, a few times",
          "Yes, monthly",
          "Yes, weekly",
          "Yes, daily"
        ),
        ordered = TRUE
      )
    )
  
  n_size <- nrow(ds_clean)
  if (n_size < 10) return(NULL)
  
  # Fórmulas
  formula_no_int <- as.formula(
    paste(
      outcome_var, "~",
      exposure_var, "+", moderator_var, "+",
      paste(confounder_list, collapse = " + "),
      "+ (1 |", cluster_var, ")"
    )
  )
  
  formula_int <- as.formula(
    paste(
      outcome_var, "~",
      exposure_var, "*", moderator_var, "+",
      paste(confounder_list, collapse = " + "),
      "+ (1 |", cluster_var, ")"
    )
  )
  
  # Ajuste con ML
  model_no_int <- lmer(
    formula_no_int,
    data = ds_clean,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  model_int <- lmer(
    formula_int,
    data = ds_clean,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  # Likelihood Ratio Test
  lrt <- anova(model_no_int, model_int)
  
  # Salida limpia
  tibble(
    Outcome    = outcome_var,
    Exposure   = exposure_var,
    Moderator  = moderator_var,
    Model      = model_label,
    N_Analysis = n_size,
    chisq      = lrt$Chisq[2],
    df         = lrt$Df[2],
    p_value    = lrt$`Pr(>Chisq)`[2]
  )
}

function(data, outcome_var, exposure_var, confounder_list, profession,
         cluster_var = "loc_2", family_type = "poisson", model_label,
         output_dir = "out/models/inter") {
  
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
    outcome_var == "suic_idea_poi"   ~ "suic"
  )
  
  # exposure abbreviation
  exp_short <- case_when( 
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
  } else {
    curr_family <- poisson(link = "log")
  }
  
  
  # Naming logic for each adjustment
  
  model_suffix <- case_when(
    model_label == "Partially adjusted" ~ "_str",
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
      "_ml.rds"
    )
  
  file_path_ml <- file.path(output_dir, file_name_ml)
  
  
  # run or load
  if(file.exists(file_path_ml)) {
    
    message(glue::glue("Loading ML model: {file_name_ml}"))
    ml <- readRDS(file_path_ml)
    
  } else {
    message(glue::glue("Fitting ML model: {file_name_ml}"))
    
    
    f_ml <- as.formula(
      paste0(outcome_var, " ~ " , exposure_var, "+", 
             paste(active_confounders, collapse = "+"),
             " + (1 | ", cluster_var, ")"))
    
    ml <- glmer(
      f_ml, 
      data = ds_clean, 
      family = curr_family,
      control = glmerControl(
        optimizer = "bobyqa", 
        optCtrl = list(maxfun = 2e5)
      )
    )
    
    
    saveRDS(ml, file_path_ml)
  }
  
  
  # Define names
  file_name_rob <- paste0("glm_", out_short, "_", exp_short , "_", prof_short, 
                          model_suffix, "_rob.rds")
  file_path_rob <- file.path(output_dir, file_name_rob)
  
  
  # run or load
  if(file.exists(file_path_rob)) {
    
    message(glue::glue("Loading Robust SE model: {file_name_rob}"))
    rob <- readRDS(file_path_rob)
    
  } else {
    message(glue::glue("Fitting Robust SE model: {file_name_rob}"))
    
    
    f_rob <- as.formula(
      paste0(outcome_var, " ~ " , exposure_var, "+", 
             paste(active_confounders, collapse = "+"),
             " + factor(", cluster_var, ")"))
    
    rob <- glm(    # robust needs to be glm not glmer
      f_rob, 
      data = ds_clean, 
      family = curr_family
    )
    
    
    saveRDS(rob, file_path_rob)
  }
  
  # Calculate adjusted risks (Marginal standardization)
  # Use ML as it accounts for the cluster variance structure
  
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
    
    model = "Multilevel Poisson",
    
    Outcome = outcome_var,
    
    Exposure = exposure_var,
    
    Profession = profession,
    
    Adjustment = model_label,
    
    beta  = beta_ml,
    
    SE    = se_ml,
    
    PR    = pr_ml,
    
    CI_l  = ci_ml[1],
    
    CI_u  = ci_ml[2],
    
    N_Analysis = n_size, # number of observations used
    
    Risk_Unexp = risk_unexp_str,
    
    Risk_Exp = risk_exp_str
    
  )
  
  
  # Extract for robust SE
  
  vcov_rob <- sandwich::vcovCL(rob, cluster = ds_clean[[cluster_var]])
  
  coef_glm <- summary(rob)$coefficients
  
  beta_rb <- coef_glm[paste0(exposure_var, "Yes"), "Estimate"]
  
  se_rb <- sqrt(vcov_rob[paste0(exposure_var, "Yes"), 
                         paste0(exposure_var, "Yes")])
  
  pr_rb <- exp(beta_rb)
  
  ci_rb <- exp(beta_rb + c(-1, 1) * 1.96 * se_rb)
  
  rb_results <- tibble(
    
    model = "Poisson + RSVE",
    
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
  
  return(list(
    multilevel = ml_results,
    robust     = rb_results
  ))
}



run_adj_moder_glm <- 
  
  function(data, outcome_var, exposure_var, confounder_list, moderator, moder_val,
           cluster_var = "loc_2", family_type = "poisson", model_label,
           output_dir = "out/models/inter/") {
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    
    # stratify data (or not)
    ds_sub <- 
      
      data |> 
      filter(.data[[moderator]] == moder_val)
    
    
    # define clean set
    
    required_vars <- 
      unique(c(outcome_var, exposure_var, cluster_var, moderator, confounder_list))
    
    out_short <- case_when(
      outcome_var == "phq_co_poi"      ~ "dep",
      outcome_var == "gad_co_poi"      ~ "anx",
      outcome_var == "suic_idea_poi"   ~ "suic",
      outcome_var == "cage_co_poi"     ~ "alc"
    )
    
    # exposure abbreviation
    exp_short <- case_when( 
      exposure_var == "work_30_dic"   ~ "har",
      exposure_var == "work_33_dic"   ~ "bul",
      exposure_var == "work_31_dic"   ~ "threats",
      exposure_var == "work_32_dic"   ~ "viol",
      exposure_var == "work_8_dic"    ~ "hrs",
      exposure_var == "work_11_dic"   ~ "nights"
    )
    
    
    # moderator abbreviation
    moder_short <- case_when(
      moderator == "work_28_dic" ~ "prothar",
      moderator == "work_29_dic" ~ "protvio",
      moderator == "work_16"     ~ "contract"
    )
    
    # moderator values abbreviation
    moder_val_short <- case_when(
      moder_val == "Yes" ~ "yes",
      moder_val == "No"  ~ "no"
    )
    
    
    
    ds_clean <- 
      ds_sub |> 
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
    } else {
      curr_family <- poisson(link = "log")
    }
    
    
    # Naming logic for each adjustment
    
    model_suffix <- case_when(
      model_label == "Partially adjusted" ~ "_str",
      model_label == "Fully adjusted"     ~ "_adj",
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
        moder_short, 
        "_", 
        moder_val,
        "_",
        model_suffix, 
        "_ml.rds"
      )
    
    file_path_ml <- file.path(output_dir, file_name_ml)
    
    
    # run or load
    if(file.exists(file_path_ml)) {
      
      message(glue::glue("Loading ML model: {file_name_ml}"))
      ml <- readRDS(file_path_ml)
      
    } else {
      message(glue::glue("Fitting ML model: {file_name_ml}"))
      
      
      f_ml <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var, "+", 
               paste(confounder_list, collapse = "+"),
               " + (1 | ", cluster_var, ")"))
      
      ml <- glmer(
        f_ml, 
        data = ds_clean, 
        family = curr_family,
        control = glmerControl(
          optimizer = "bobyqa", 
          optCtrl = list(maxfun = 2e5)
        )
      )
      
      
      saveRDS(ml, file_path_ml)
    }
    
    
    # Define names
    file_name_rob <- 
      paste0("glm_", 
             out_short, 
             "_", 
             exp_short , 
             "_", 
             moder_short, 
             "_", 
             moder_val,
             "_",
             model_suffix, 
             "_rob.rds")
    
    file_path_rob <- file.path(output_dir, file_name_rob)
    
    
    # run or load
    if(file.exists(file_path_rob)) {
      
      message(glue::glue("Loading Robust SE model: {file_name_rob}"))
      rob <- readRDS(file_path_rob)
      
    } else {
      message(glue::glue("Fitting Robust SE model: {file_name_rob}"))
      
      
      f_rob <- as.formula(
        paste0(outcome_var, " ~ " , exposure_var, "+", 
               paste(confounder_list, collapse = "+"),
               " + factor(", cluster_var, ")"))
      
      rob <- glm(    # robust needs to be glm not glmer
        f_rob, 
        data = ds_clean, 
        family = curr_family
      )
      
      
      saveRDS(rob, file_path_rob)
    }
    
    # Calculate adjusted risks (Marginal standardization)
    # Use ML as it accounts for the cluster variance structure
    
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
      
      model = "Multilevel Poisson",
      
      Outcome = outcome_var,
      
      Exposure = exposure_var,
      
      Moderator = moderator,
      
      ModeratorLevel = moder_val,
      
      Adjustment = model_label,
      
      beta  = beta_ml,
      
      SE    = se_ml,
      
      PR    = pr_ml,
      
      CI_l  = ci_ml[1],
      
      CI_u  = ci_ml[2],
      
      N_Analysis = n_size, # number of observations used
      
      Risk_Unexp = risk_unexp_str,
      
      Risk_Exp = risk_exp_str
      
    )
    
    
    # Extract for robust SE
    
    vcov_rob <- sandwich::vcovCL(rob, cluster = ds_clean[[cluster_var]])
    
    coef_glm <- summary(rob)$coefficients
    
    beta_rb <- coef_glm[paste0(exposure_var, "Yes"), "Estimate"]
    
    se_rb <- sqrt(vcov_rob[paste0(exposure_var, "Yes"), 
                           paste0(exposure_var, "Yes")])
    
    pr_rb <- exp(beta_rb)
    
    ci_rb <- exp(beta_rb + c(-1, 1) * 1.96 * se_rb)
    
    rb_results <- tibble(
      
      model = "Poisson + RSVE",
      
      Outcome = outcome_var,
      
      Exposure = exposure_var,
      
      Moderator = moderator,
      
      ModeratorLevel = moder_val,
      
      Adjustment = model_label,
      
      beta  = beta_rb,
      
      SE    = se_rb,
      
      PR    = pr_rb,
      
      CI_l  = ci_rb[1],
      
      CI_u  = ci_rb[2],
      
      N_Analysis = n_size # number of observations used
    )
    
    return(list(
      multilevel = ml_results,
      robust     = rb_results
    ))
  }

run_adj_moder_prevs<-
  
  function(data, outcome_var, exposure_var, moderator_var,
           confounder_list,
           profession = "All",
           cluster_var = "loc_2") {
    
    # ---- Stratify ----
    ds_sub <- if (profession == "All") {
      data
    } else {
      data |> filter(work_2 == profession)
    }
    
    active_confounders <- if (profession == "All") {
      confounder_list
    } else {
      setdiff(confounder_list, "work_2")
    }
    
    required_vars <- unique(c(outcome_var,
                              exposure_var,
                              moderator_var,
                              cluster_var,
                              active_confounders))
    
    ds_clean <- ds_sub |> 
      select(all_of(required_vars)) |> 
      drop_na()
    
    if(nrow(ds_clean) < 50) return(NULL)
    
    # Binary coding
    ds_clean[[exposure_var]]  <- factor(ds_clean[[exposure_var]],
                                        levels = c("No","Yes"))
    
    ds_clean[[moderator_var]] <- factor(ds_clean[[moderator_var]],
                                        levels = c("No","Yes"))
    
    y_raw <- as.numeric(ds_clean[[outcome_var]])
    if(max(y_raw, na.rm = TRUE) > 1) {
      ds_clean[[outcome_var]] <- y_raw - 1
    }
    
    # ---- Model ----
    f_ml <- as.formula(
      paste0(
        outcome_var, " ~ ",
        exposure_var, " * ", moderator_var, " + ",
        paste(active_confounders, collapse = "+"),
        " + (1 | ", cluster_var, ")"
      )
    )
    
    ml <- glmer(
      f_ml,
      data = ds_clean,
      family = poisson(link = "log"),
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
    
    coefs <- summary(ml)$coefficients
    vc    <- vcov(ml)
    
    # Coefficient names
    b1_name <- paste0(exposure_var, "Yes")
    b2_name <- paste0(moderator_var, "Yes")
    b3_name <- paste0(exposure_var, "Yes:",
                      moderator_var, "Yes")
    
    beta1 <- coefs[b1_name, "Estimate"]
    beta2 <- coefs[b2_name, "Estimate"]
    beta3 <- coefs[b3_name, "Estimate"]
    
    # ---- Compute RERI ----
    RR_11 <- exp(beta1 + beta2 + beta3)
    RR_10 <- exp(beta1)
    RR_01 <- exp(beta2)
    
    RERI <- RR_11 - RR_10 - RR_01 + 1
    
    # ---- Delta method ----
    g1 <- RR_11 - RR_10
    g2 <- RR_11 - RR_01
    g3 <- RR_11
    
    grad <- c(g1, g2, g3)
    
    Sigma <- vc[c(b1_name, b2_name, b3_name),
                c(b1_name, b2_name, b3_name)]
    
    var_reri <- t(grad) %*% Sigma %*% grad
    se_reri  <- sqrt(var_reri)
    
    CI_l <- RERI - 1.96 * se_reri
    CI_u <- RERI + 1.96 * se_reri
    
    # ---- Addition: multiplicative interaction ----
    
    MI <- exp(beta3)
    se_beta3 <- coefs[b3_name, "Std. Error"]
    
    MI_CI_l <- exp(beta3 - 1.96 * se_beta3)
    MI_CI_u <- exp(beta3 + 1.96 * se_beta3)
    
    tibble(
      Outcome    = outcome_var,
      Exposure   = exposure_var,
      Moderator  = moderator_var,
      Profession = profession,
      RR11 = as.numeric(RR_11),
      RR10 = as.numeric(RR_10),
      RR01 = as.numeric(RR_01),
      RERI       = as.numeric(RERI),
      SE         = as.numeric(se_reri),
      CI_l       = as.numeric(CI_l),
      CI_u       = as.numeric(CI_u),
      MI       = as.numeric(MI),
      MI_CI_l  = as.numeric(MI_CI_l),
      MI_CI_u  = as.numeric(MI_CI_u),
      N          = nrow(ds_clean)
    )
  }
