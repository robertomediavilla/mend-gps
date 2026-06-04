# Weighted prevalences ----------------------------------------------------

# These functions estimate prevalence of mental health problems, 
# by job position, and across countries


calc_out_prev <- 
  
  function(var, outcome, design, level = "Yes") {
    
    
    form_outcome <- as.formula(paste0("~I(", outcome, " =='", level, "')"))
    form_by <- as.formula(paste0("~loc_2"))
    
    prev <- svyby(
      form_outcome,
      form_by,
      design = design,
      FUN = svyciprop,
      vartype = c("ci"),
      method = "beta",
      level = 0.95,
      na.rm = TRUE
    ) |> 
      as_tibble() 
    
    prop_col <- grep("^I\\(", names(prev), value = TRUE)
    
    prev <- 
      prev |> 
      mutate(
        Variable = loc_2,
        prev = 100 * .data[[prop_col]],
        ci_low = 100 * ci_l,
        ci_high = 100 * ci_u,
        Prevalence = sprintf("%.1f%% (%.1f-%.1f)", prev, ci_low, ci_high)
      ) |>  
      select(Variable, Prevalence)
    
    return(prev)
  }

get_column_stats <- 
  
  function(design_subset, outcome) {
    
    # calculate overall for this specific subgroup
    overall_stat <- svyciprop(
      as.formula(paste0("~I(", outcome, " == 'Yes')")),
      design = design_subset, 
      vartype = "ci",
      method = "beta",
      na.rm = TRUE
    )
    
    overall_row <- tibble(
      Variable = "Overall",
      Prevalence = 
        sprintf(
          "%.1f%% (%.1f-%.1f)", 100 * coef(overall_stat),
          100 * confint(overall_stat)[1], 100* confint(overall_stat)[2])
    )
    
    rows_stats <- map_dfr("loc_2", 
                          calc_out_prev, 
                          outcome = outcome,
                          design = design_subset)
    
    bind_rows(overall_row, rows_stats)
    
  }

generate_outcome_table <- 
  
  function(outcome_var, outcome_title) {
    
    # generate subsets
    ds_doc <-
      subset(svy_design,
             work_2 == "Doctor")
    
    ds_nur <-
      subset(svy_design,
             work_2 == "Nurse")
    
    
    # unweighted n for headers
    n_doc <- 
      ds |> 
      filter(
        work_2 == "Doctor" &
          !is.na(.data[[outcome_var]])) |> 
      count() |> 
      pull() |> 
      style_number()
    
    n_nur <- 
      ds |> 
      filter(
        work_2 == "Nurse" &
          !is.na(.data[[outcome_var]])) |> 
      count() |> 
      pull() |> 
      style_number()
    
    
    # calculations
    
    
    
    res_doc <- get_column_stats(ds_doc, outcome_var) |>
      rename(Doc = Prevalence)
    
    res_nur <- get_column_stats(ds_nur, outcome_var) |>
      rename(Nur = Prevalence)
    
    
    
    
    # merge
    final_df <- 
      
      res_doc |>
      left_join(res_nur, by = "Variable")
    
    final_df #remove
    
    # final_df |>
    #   gt() |>
    #   cols_label(
    #     Variable = md("**Country**"),
    #     Doc = md(glue("**Doctor**<br>(N = {n_doc})")),
    #     Nur = md(glue("**Nurse**<br>(N = {n_nur})")),
    #   ) |>
    #   tab_header(
    #     title = md(glue("**{outcome_title}**"))
    #   ) |>
    #   cols_align(align = "center", columns = contains("_")) |>
    #   tab_footnote(footnote = "Weighted Prevalence (95% CI)")
  }


## Exposures

# Calculate subgroup weighted prevalences

calc_w_prev_exp <- 
  function(design_subset, exposure, level = "Yes") {
  
  # Get label from attribute or fallback to name
  var_label <- tryCatch(attr(design_subset$variables[[exposure]], "label"),
                        error = function(e) exposure)
  if (is.null(var_label)) var_label <- exposure
  
  prev <- tryCatch({
    svyciprop(
      as.formula(paste0("~I(", exposure, " == '", level, "')")),
      design = design_subset,
      vartype = "ci",
      method = "beta",
      na.rm = TRUE
    )
  }, error = function(e) return(NULL)) 
  
  if(is.null(prev)) return(tibble(Variable = var_label, Prevalence = "-"))
  
  tibble(
    Variable = var_label,
    Prevalence = sprintf("%.1f%% (%.1f-%.1f)", 
                         100 * coef(prev),
                         100 * confint(prev)[1], 
                         100 * confint(prev)[2])
  )
}


# Get weighted prevalences

calc_prev_exposure_w <- 
  function(design_subset, exposures_list) {
  
    map_dfr(exposures_list, ~calc_w_prev_exp(design_subset, .x))

    }


# Build weighted exposure table 

gen_w_exp_tbl <- function(metadata_df) {
  
  # The list of raw column names from your metadata
  exposures_list <- metadata_df$Column
  
  # --- 1. Subsets (using your current objects) ---
  d_f_doc   <- subset(svy_design, scdm_2_rec == "Female" & work_2 == "Doctor")
  d_f_nur   <- subset(svy_design, scdm_2_rec == "Female" & work_2 == "Nurse")
  d_m_doc   <- subset(svy_design, scdm_2_rec == "Male" & work_2 == "Doctor")
  d_m_nur   <- subset(svy_design, scdm_2_rec == "Male" & work_2 == "Nurse")
  d_all_doc <- subset(svy_design, work_2 == "Doctor")
  d_all_nur <- subset(svy_design, work_2 == "Nurse")
  
  # --- 2. Helper for N counts ---
  get_n <- function(gen = NULL, work) {
    res <- ds_w %>% filter(work_2 == work)
    if (!is.null(gen)) res <- res %>% filter(scdm_2_rec == gen)
    # Using format for comma separation (e.g., 1,000)
    format(nrow(res), big.mark = ",")
  }
  
  # --- 3. Calculations ---
  # We use the list of column names for the calculation
  calc_cols <- list(
    M_Doc = d_m_doc, F_Doc = d_f_doc, All_Doc = d_all_doc,
    M_Nur = d_m_nur, F_Nur = d_f_nur, All_Nur = d_all_nur
  ) %>% 
    map(~ calc_prev_exposure_w(.x, exposures_list))
  
  # --- 4. Join and Format ---
  # We merge all columns back together
  final_df <- calc_cols$M_Doc %>% rename(M_Doc = Prevalence) %>%
    left_join(calc_cols$F_Doc %>% rename(F_Doc = Prevalence), by = "Variable") %>%
    left_join(calc_cols$All_Doc %>% rename(All_Doc = Prevalence), by = "Variable") %>%
    left_join(calc_cols$M_Nur %>% rename(M_Nur = Prevalence), by = "Variable") %>%
    left_join(calc_cols$F_Nur %>% rename(F_Nur = Prevalence), by = "Variable") %>%
    left_join(calc_cols$All_Nur %>% rename(All_Nur = Prevalence), by = "Variable") %>%
    # Map the internal column names to your "Label" names using metadata
    left_join(metadata_df, by = c("Variable" = "Column")) %>%
    select(Category, Label, everything(), -Variable) 
  
  # --- 5. Generate GT Table ---
  final_df %>%
    group_by(Category) %>%
    gt(rowname_col = "Label") %>%
    tab_header(title = "Weighted Prevalence of Exposures") %>%
    cols_label(
      M_Doc   = md(glue("**Male**<br>(N = {get_n('Male', 'Doctor')})")),
      F_Doc   = md(glue("**Female**<br>(N = {get_n('Female', 'Doctor')})")),
      All_Doc = md(glue("**Overall**<br>(N = {get_n(work='Doctor')})")),
      M_Nur   = md(glue("**Male**<br>(N = {get_n('Male', 'Nurse')})")),
      F_Nur   = md(glue("**Female**<br>(N = {get_n('Female', 'Nurse')})")),
      All_Nur = md(glue("**Overall**<br>(N = {get_n(work='Nurse')})"))
    ) %>%
    tab_spanner(label = md("**Doctor**"), columns = c(M_Doc, F_Doc, All_Doc)) %>%
    tab_spanner(label = md("**Nurse**"), columns = c(M_Nur, F_Nur, All_Nur)) %>%
    cols_align(align = "center", columns = contains("_")) %>%
    tab_options(row_group.font.weight = "bold") %>%
    tab_footnote(footnote = "Weighted Prevalence (95% CI) using Taylor Series Linearization.")
}

## Outcomes
calc_w_prev_out <- 
  function(design_subset, outcome, level = "Yes") {
    
    var_label <- tryCatch(attr(design_subset$variables[[outcome]], "label"),
                          error = function(e) outcome)
    
    if (is.null(var_label)) { 
      var_label <- outcome 
    }
    
    prev <-tryCatch({
      svyciprop(
        as.formula(paste0("~I(", outcome, " == '", level, "')")),
        design = design_subset,
        vartype = "ci",
        method = "beta",
        na.rm = TRUE
      )
    }, error = function(e) return(NULL)) 
    
    if(is.null(prev)) return(tibble(Variable = var_label, Prevalence = "-"))
    
    tibble(
      Variable = var_label,
      Prevalence = sprintf(
        "%.1f%% (%.1f-%.1f)", 
        100 * coef(prev),
        100 * confint(prev)[1], 
        100 * confint(prev)[2]
      ) 
    )
  }

# Get outcome weighted prevalences

calc_prev_outcome_w <- function(design_subset, out_dic_list) {
  map_dfr(out_dic_list, ~calc_w_prev_out(design_subset, .x))
}



gen_w_out_tbl <- function(out_dic_list, outcome_order = out_order_names) {
  
  # --- 1. Subsets (svy_design) ---
  d_f_doc   <- subset(svy_design, scdm_2_rec == "Female" & work_2 == "Doctor")
  d_f_nur   <- subset(svy_design, scdm_2_rec == "Female" & work_2 == "Nurse")
  d_m_doc   <- subset(svy_design, scdm_2_rec == "Male" & work_2 == "Doctor")
  d_m_nur   <- subset(svy_design, scdm_2_rec == "Male" & work_2 == "Nurse")
  d_all_doc <- subset(svy_design, work_2 == "Doctor")
  d_all_nur <- subset(svy_design, work_2 == "Nurse")
  
  # --- 2. Helper for N counts (ds_w) ---
  get_n <- function(gen = NULL, work) {
    res <- ds_w %>% filter(work_2 == work)
    if (!is.null(gen)) res <- res %>% filter(scdm_2_rec == gen)
    format(nrow(res), big.mark = ",")
  }
  
  # --- 3. Calculations ---
  calc_cols <- list(
    M_Doc = d_m_doc, F_Doc = d_f_doc, All_Doc = d_all_doc,
    M_Nur = d_m_nur, F_Nur = d_f_nur, All_Nur = d_all_nur
  ) %>% 
    map(~ calc_prev_outcome_w(.x, out_dic_list))
  
  # --- 4. Join, Re-label, and Order ---
  final_df <- calc_cols$M_Doc %>% rename(M_Doc = Prevalence) %>%
    left_join(calc_cols$F_Doc   %>% rename(F_Doc = Prevalence),   by = "Variable") %>%
    left_join(calc_cols$All_Doc %>% rename(All_Doc = Prevalence), by = "Variable") %>%
    left_join(calc_cols$M_Nur   %>% rename(M_Nur = Prevalence),   by = "Variable") %>%
    left_join(calc_cols$F_Nur   %>% rename(F_Nur = Prevalence),   by = "Variable") %>%
    left_join(calc_cols$All_Nur %>% rename(All_Nur = Prevalence), by = "Variable") %>%
    mutate(
      Variable = case_when(
        Variable == "phq_co"    ~ "Probable depression (PHQ-9)",
        Variable == "gad_co"    ~ "Probable anxiety (GAD-7)",
        Variable == "cage_co"   ~ "Probable alcohol dependence disorder (CAGE)",
        Variable == "suic_idea" ~ "Passive suicide thoughts (Item 9 of PHQ-9)",
        TRUE ~ Variable 
      ),
      # Apply the factor order based on your custom list
      Variable = factor(Variable, levels = outcome_order)
    ) %>% 
    arrange(Variable)
  
  # --- 5. Generate GT Table ---
  final_df %>%
    gt() %>%
    tab_header(title = "Weighted Prevalence of Mental Health Outcomes") %>%
    cols_label(
      Variable = md("**Outcome**"),
      M_Doc    = md(glue("**Male**<br>(N = {get_n('Male', 'Doctor')})")),
      F_Doc    = md(glue("**Female**<br>(N = {get_n('Female', 'Doctor')})")),
      All_Doc  = md(glue("**Overall**<br>(N = {get_n(work='Doctor')})")),
      M_Nur    = md(glue("**Male**<br>(N = {get_n('Male', 'Nurse')})")),
      F_Nur    = md(glue("**Female**<br>(N = {get_n('Female', 'Nurse')})")),
      All_Nur  = md(glue("**Overall**<br>(N = {get_n(work='Nurse')})"))
    ) %>%
    tab_spanner(label = md("**Doctor**"), columns = c(M_Doc, F_Doc, All_Doc)) %>%
    tab_spanner(label = md("**Nurse**"),  columns = c(M_Nur, F_Nur, All_Nur)) %>%
    cols_align(align = "center", columns = contains("_")) %>%
    tab_footnote(
      footnote = "Weighted Prevalence (95% CI). PHQ-9: 9-item Patient Health Questionnaire (cutoff ≥ 10); GAD-7: 7-item anxiety scale (cutoff ≥ 10); CAGE: 4-item scale for alcohol use. Passive suicide thoughts based on Item 9 of PHQ-9."
    )
}

# Unweighted exposures by subgroups (SCDM) --------------------------------
calc_unweighted_row_exp <-
  function(data, exposure_var) {
    
    # Clean vector
    vec <- data[[exposure_var]]
    vec <- vec[!is.na(vec)] # Remove NAs
    
    n_total <- length(vec)
    
    # Placeholders if data is empty (prevents errors)
    if (n_total == 0) return("NA (NA-NA)")
    
    n_yes <- sum(vec == "Yes")
    
    # CI using prop.test (Wilson score is standard for props)
    ptest <- prop.test(n_yes, n_total, conf.level = 0.95)
    
    est <- ptest$estimate * 100
    ci_low <- ptest$conf.int[1] * 100
    ci_high <- ptest$conf.int[2] * 100
    
    sprintf("%.1f%% (%.1f-%.1f)", est, ci_low, ci_high)
  }



# Loop through variables to calculate row for each subgroup

get_unweighted_stats_exp <-
  function(subset_df, exposures, scdm) {
    
    # Overall Row
    overall_val <- calc_unweighted_row_exp(subset_df, exposures)
    
    overall_row <- tibble(
      Variable = "Overall",
      Category = "Overall",
      Prevalence = overall_val
    )
    
    # Loop through scdm variables
    rows_stats <- map_dfr(scdm, function(var) {
      
      # Label
      var_lbl <- tryCatch(attr(subset_df[[var]], "label"),
                          error = function(e) var)
      
      if(is.null(var_lbl)) var_lbl <- var
      
      # group by the variable categories and calculate stats
      subset_df |>
        filter(!is.na(.data[[var]])) |> # remove missing categories
        group_by(Category = .data[[var]]) |>
        summarise(
          Prevalence = calc_unweighted_row_exp(pick(everything()), exposures)
        ) |>
        mutate(Variable = var_lbl) |>
        mutate(Category = as.character(Category)) |>
        select(Variable, Category, Prevalence)
    })
    
    # Combine
    bind_rows(overall_row, rows_stats)
  }


# Build unweighted table

get_exposure_data_long <-
  function(data, exposure_var, exposure_title, scdm) {
    
    # Subsets
    d_m_doc <-
      data |>
      filter(scdm_2_rec == "Male",
             work_2 == "Doctor")
    
    d_m_nur <-
      data |>
      filter(scdm_2_rec == "Male",
             work_2 == "Nurse")
    
    d_f_doc <-
      data |>
      filter(scdm_2_rec == "Female",
             work_2 == "Doctor")
    
    d_f_nur <-
      data |>
      filter(scdm_2_rec == "Female",
             work_2 == "Nurse")
    
    d_all_doc <-
      data |> filter(work_2 == "Doctor")
    
    d_all_nur <-
      data |> filter(work_2 == "Nurse")
    
    # Ns for each exposure
    n_total <-
      
      data |>
      filter(!is.na(.data[[exposure_var]])) |>
      nrow() |>
      style_number()
    
    # Create a label that includes N (e.g., "Depression (N = 1,203)")
    lbl_combined <- glue("{exposure_title} (N = {n_total})")
    
    
    
    # Calculate Stats
    res_m_doc <-
      get_unweighted_stats_exp(d_m_doc, exposure_var, scdm) |>
      rename(M_Doc = Prevalence)
    
    res_m_nur <-
      get_unweighted_stats_exp(d_m_nur, exposure_var, scdm) |>
      rename(M_Nur = Prevalence)
    
    res_f_doc <-
      get_unweighted_stats_exp(d_f_doc, exposure_var, scdm) |>
      rename(F_Doc = Prevalence)
    
    res_f_nur <-
      get_unweighted_stats_exp(d_f_nur, exposure_var, scdm) |>
      rename(F_Nur = Prevalence)
    
    res_all_doc <-
      get_unweighted_stats_exp(d_all_doc, exposure_var, scdm) |>
      rename(All_Doc = Prevalence)
    
    res_all_nur <-
      get_unweighted_stats_exp(d_all_nur, exposure_var, scdm) |>
      rename(All_Nur = Prevalence)
    
    # Merge
    res_m_doc |>
      left_join(res_f_doc, by = c("Variable", "Category")) |>
      left_join(res_all_doc, by = c("Variable", "Category")) |>
      left_join(res_m_nur, by = c("Variable", "Category")) |>
      left_join(res_f_nur, by = c("Variable", "Category")) |>
      left_join(res_all_nur, by = c("Variable", "Category")) |>
      mutate(Exposure_Label = lbl_combined) |> # Add the exposure identifier
      select(Exposure_Label, Variable, Category, everything())
    
    
  }

# table generator

generate_exposure_by_scdm_table <-
  function(data, exposure_list, scdm, desired_order) {
    
    # Loop through exposures and stack data
    long_df <- map_dfr(names(exposure_list), function(var_name) {
      title <- exposure_list[[var_name]]
      get_exposure_data_long(data, var_name, title, scdm)
    })
    
   
    # Force custom variable order
    long_df <- 
      
      long_df |> 
      mutate(
        Exposure_Label = factor(Exposure_Label, levels = unique(Exposure_Label)),
        Variable = factor(Variable, levels = desired_order)) |> 
      arrange(Exposure_Label, Variable)
    
    # Prepare for gt: Insert "Header Rows" for Variables
    # We want a row that says "Age Group" followed by rows "<30", "30-50"
    formatted_df <- 
      
      long_df |>
      group_by(Exposure_Label, Variable) |>
      group_split() |>
      map_dfr(function(chunk) {
        # Create a dummy header row
        header_row <- tibble(
          Exposure_Label = unique(chunk$Exposure_Label),
          Variable = unique(chunk$Variable),
          Category = as.character(unique(chunk$Variable)), # The Category column becomes the Header text
          M_Doc = NA_character_, F_Doc = NA_character_, All_Doc = NA_character_,
          M_Nur = NA_character_, F_Nur = NA_character_, All_Nur = NA_character_
        )
        
        # Determine if this is the "Overall" row (optional: skip header for Overall if desired)
        if(unique(chunk$Variable) == "Overall") {
          return(chunk)
        } else {
          return(bind_rows(header_row, chunk))
        }
      }) |> 
      # Add flags to tell flextable exactly what to bold and indent
      mutate(
        is_var_header = (Category == Variable & Variable != "Overall"),
        is_category   = (Category != Variable & Variable != "Overall")
      )
    
    # Convert to flextable grouped data
    grouped_df <- as_grouped_data(x = formatted_df, 
                                  groups = "Exposure_Label")
    
    # Copy exposure title into category column
    grouped_df$Category <- ifelse(
      !is.na(grouped_df$Exposure_Label),
      as.character(grouped_df$Exposure_Label),
      as.character(grouped_df$Category)
    )
    
    # Build ft
    ft <- 
      
      flextable(
        grouped_df,
        # Explicitly define which columns to show (this safely hides Variable and our flags)
        col_keys = c("Category", "M_Doc", "F_Doc", "All_Doc", "M_Nur", "F_Nur", "All_Nur")
      ) |>
      
      # Headers and spanners
      set_header_labels(
        Category = "Subgroup",
        M_Doc = "Male", F_Doc = "Female", All_Doc = "Overall",
        M_Nur = "Male", F_Nur = "Female", All_Nur = "Overall"
      ) |>
      add_header_row(
        values = c("", "Doctor", "Nurse"),
        colwidths = c(1, 3, 3)
      ) |> 
      
      # Alignment
      align(
        j = c("M_Doc", "F_Doc", "All_Doc", "M_Nur", "F_Nur", "All_Nur"),
        align = "center", 
        part = "all"
        ) |> 
      align( 
        j = "Category",
        align = "left",
        part = "all") |> 
      
      # Format and style
      bold(i = ~is.na(Exposure_Label), j = "Category") |> 
      bold(i = ~is_var_header == TRUE, j = "Category") |> 
      padding(i = ~is_category == TRUE, j = "Category", padding.left = 15) |> 
      
      # Display empty space instead of NA
      colformat_char(na_str = "") |> 
      bold(part = "header") |> 
      theme_booktabs() |> 
      add_header_lines(values = "Exposures by sociodemographic groups") |> 
      add_footer_lines(values = "Unweighted prevalence (95% CI). N represents total sample for that exposure.") |> 
      autofit()
    
    return(ft)
      
   
  }

# Unweighted outcomes by subgroups (SCDM) ------------------------------

# Calculation of unweighted rows

calc_unweighted_row_out <-
  function(data, outcome_var) {
    
    # Clean vector
    vec <- data[[outcome_var]]
    vec <- vec[!is.na(vec)] # Remove NAs
    
    n_total <- length(vec)
    
    # Placeholders if data is empty (prevents errors)
    if (n_total == 0) return("NA (NA-NA)")
    
    n_yes <- sum(vec == "Yes")
    
    # CI using prop.test (Wilson score is standard for props)
    ptest <- prop.test(n_yes, n_total, conf.level = 0.95)
    
    est <- ptest$estimate * 100
    ci_low <- ptest$conf.int[1] * 100
    ci_high <- ptest$conf.int[2] * 100
    
    sprintf("%.1f%% (%.1f-%.1f)", est, ci_low, ci_high)
  }



# Loop through variables to calculate row for each subgroup

get_unweighted_stats_out <-
  function(subset_df, outcome, scdm) {
    
    # Overall Row
    overall_val <- calc_unweighted_row_out(subset_df, outcome)
    
    overall_row <- tibble(
      Variable = "Overall",
      Category = "Overall",
      Prevalence = overall_val
    )
    
    # Loop through scdm variables
    rows_stats <- map_dfr(scdm, function(var) {
      
      # Label
      var_lbl <- tryCatch(attr(subset_df[[var]], "label"),
                          error = function(e) var)
      
      if(is.null(var_lbl)) var_lbl <- var
      
      # group by the variable categories and calculate stats
      subset_df |>
        filter(!is.na(.data[[var]])) |> # remove missing categories
        group_by(Category = .data[[var]]) |>
        summarise(
          Prevalence = calc_unweighted_row_out(pick(everything()), outcome)
        ) |>
        mutate(Variable = var_lbl) |>
        mutate(Category = as.character(Category)) |>
        select(Variable, Category, Prevalence)
    })
    
    # Combine
    bind_rows(overall_row, rows_stats)
  }


# Build unweighted table

get_outcome_data_long <-
  function(data, outcome_var, outcome_title, scdm) {
    
    # Subsets
    d_m_doc <-
      data |>
      filter(scdm_2_rec == "Male",
             work_2 == "Doctor")
    
    d_m_nur <-
      data |>
      filter(scdm_2_rec == "Male",
             work_2 == "Nurse")
    
    d_f_doc <-
      data |>
      filter(scdm_2_rec == "Female",
             work_2 == "Doctor")
    
    d_f_nur <-
      data |>
      filter(scdm_2_rec == "Female",
             work_2 == "Nurse")
    
    d_all_doc <-
      data |> filter(work_2 == "Doctor")
    
    d_all_nur <-
      data |> filter(work_2 == "Nurse")
    
    # Ns for each outcome
    n_total <-
      
      data |>
      filter(!is.na(.data[[outcome_var]])) |>
      nrow() |>
      style_number()
    
    # Create a label that includes N (e.g., "Depression (N = 1,203)")
    lbl_combined <- glue("{outcome_title} (N = {n_total})")
    
    
    
    # Calculate Stats
    res_m_doc <-
      get_unweighted_stats_out(d_m_doc, outcome_var, scdm) |>
      rename(M_Doc = Prevalence)
    
    res_m_nur <-
      get_unweighted_stats_out(d_m_nur, outcome_var, scdm) |>
      rename(M_Nur = Prevalence)
    
    res_f_doc <-
      get_unweighted_stats_out(d_f_doc, outcome_var, scdm) |>
      rename(F_Doc = Prevalence)
    
    res_f_nur <-
      get_unweighted_stats_out(d_f_nur, outcome_var, scdm) |>
      rename(F_Nur = Prevalence)
    
    res_all_doc <-
      get_unweighted_stats_out(d_all_doc, outcome_var, scdm) |>
      rename(All_Doc = Prevalence)
    
    res_all_nur <-
      get_unweighted_stats_out(d_all_nur, outcome_var, scdm) |>
      rename(All_Nur = Prevalence)
    
    # Merge
    res_m_doc |>
      left_join(res_f_doc, by = c("Variable", "Category")) |>
      left_join(res_all_doc, by = c("Variable", "Category")) |>
      left_join(res_m_nur, by = c("Variable", "Category")) |>
      left_join(res_f_nur, by = c("Variable", "Category")) |>
      left_join(res_all_nur, by = c("Variable", "Category")) |>
      mutate(Outcome_Label = lbl_combined) |> # Add the outcome identifier
      select(Outcome_Label, Variable, Category, everything())
    
    
  }

# table generator

generate_outcome_by_scdm_table <-
  function(data, outcome_list, scdm, desired_order) {
    
    # Loop through exposures and stack data
    long_df <- map_dfr(names(outcome_list), function(var_name) {
      title <- outcome_list[[var_name]]
      get_outcome_data_long(data, var_name, title, scdm)
    })
    
   
    # Define order of scdm vars
    long_df <- 
      
      long_df |> 
      mutate(
        Outcome_Label = factor(Outcome_Label, levels = unique(Outcome_Label)),
        Variable = factor(Variable, 
                          levels = unique(c(desired_order, unique(Variable))))
      ) |> 
      arrange(Outcome_Label, Variable)
    
    
    # Prepare for gt: Insert "Header Rows" for Variables
    # We want a row that says "Age Group" followed by rows "<30", "30-50"
    formatted_df <- 
      
      long_df |>
      group_by(Outcome_Label, Variable) |>
      group_split() |>
      map_dfr(function(chunk) {
        # Create a dummy header row
        header_row <- tibble(
          Outcome_Label = unique(chunk$Outcome_Label),
          Variable = unique(chunk$Variable),
          Category = unique(chunk$Variable), # The Category column becomes the Header text
          M_Doc = NA_character_, F_Doc = NA_character_, All_Doc = NA_character_,
          M_Nur = NA_character_, F_Nur = NA_character_, All_Nur = NA_character_
        )
        
        # Determine if this is the "Overall" row (optional: skip header for Overall if desired)
        if(unique(chunk$Variable) == "Overall") {
          return(chunk)
        } else {
          return(bind_rows(header_row, chunk))
        }
      }) |> 
      # Add flags to tell flextable exactly what to bold and indent
      mutate(
        is_var_header = (Category == Variable & Variable != "Overall"),
        is_category   = (Category != Variable & Variable != "Overall")
      )
    
    grouped_df <- 
      
      as_grouped_data(x = formatted_df, groups = "Outcome_Label")
    
    # Copy outcome title into category column so it doesn't get hidden
    grouped_df$Category <- ifelse(
      !is.na(grouped_df$Outcome_Label), 
      as.character(grouped_df$Outcome_Label), 
      as.character(grouped_df$Category)
    )
    
    # Build ft
    ft <- 
      
      flextable(
        grouped_df, 
        col_keys = c("Category", "M_Doc", "F_Doc", "All_Doc", "M_Nur", "F_Nur", "All_Nur")
      ) |> 
      set_header_labels(
        Category = "Subgroup",
        M_Doc = "Male", F_Doc = "Female", All_Doc = "Overall",
        M_Nur = "Male", F_Nur = "Female", All_Nur = "Overall"
      ) |>
      add_header_row(
        values = c("", "Doctor", "Nurse"),
        colwidths = c(1, 3, 3)
      ) |> 
      align(
        j = c("M_Doc", "F_Doc", "All_Doc", "M_Nur", "F_Nur", "All_Nur"),
        align = "center",
        part = "all"
      ) |> 
      align(
        j = "Category",
        align = "left",
        part = "all"
      ) |> 
      bold(i = ~!is.na(Outcome_Label), j = "Category") |> 
      bold(i = ~is_var_header == TRUE, j = "Category") |> 
      padding(i = ~is_category == TRUE, j = "Category", padding.left = 15) |> 
      
      # Display empty space instead of NA
      colformat_char(na_str = "") |> 
      bold(part = "header") |> 
      theme_booktabs() |> 
      add_header_lines(values = "Outcomes by Sociodemographics (Unweighted)") |> 
      add_footer_lines(values = "Unweighted Prevalence (95% CI). N represents total sample for that outcome.") |> 
      autofit()
    
    return(ft)
    
    
  }


# Unweighted prevalences --------------------------------------------------

# This function runs MA across outcomes using observed (i.e. unweighted) data


run_ma <- function(outcome) {
  ma <- ds |>
    drop_na(loc_2, all_of(outcome)) |>
    group_by(loc_2) |>
    summarise(
      event = sum(.data[[outcome]] == "Yes"),
      n = n(),
      .groups = "drop"
    ) |>
    metaprop(
      event = event,
      n = n,
      studlab = loc_2,
      sm = "PFT",
      backtransf = TRUE
    )
  
  #saves pdf
  
  pdf_filename <- paste0("out/figures/", outcome, "_forest.pdf")
  pdf(pdf_filename, width = 12, height = 8)
  forest(ma)
  dev.off()  # cerrar PDF
  
  # saves R object
  
  forest(ma, main = paste("Forest plot -", outcome))
  fp <- recordPlot()
  
  #output
  
  list(
    ma = ma,
    summary = summary(ma),
    forest = fp,
    pdf_file = pdf_filename
  )
}

run_mr <- function(ds, outcome, subgroup_var) {
  
  dat <- ds |>
    drop_na(loc_2, all_of(subgroup_var), all_of(outcome)) |>
    mutate(subgroup = .data[[subgroup_var]]) |>
    group_by(loc_2, subgroup) |>
    summarise(
      event = sum(.data[[outcome]] == "Yes"),
      n = n(),
      .groups = "drop"
    )
  
  metaprop(
    event = dat$event,
    n = dat$n,
    studlab = dat$loc_2,
    subgroup = dat$subgroup,
    sm = "PFT",
    backtransf = TRUE
  )
}


# Unweighted vs weighted --------------------------------------------------

get_country_prev_survey <- 
  function(group_var) {
  
  # Unweighted (Standard dplyr on the raw dataframe)
  unw <- 
    
    ds_w |>
    filter(!is.na(.data[[group_var]])) |>
    group_by(loc_2) |>
    summarise(Unweighted = mean(.data[[group_var]] == "Yes") * 100) |>
    select(Country = loc_2, Unweighted)
  
  # Weighted (Using survey::svyby)
  fmla <- as.formula(paste0("~", group_var))
  
  wht_raw <- 
    svyby(
      formula = fmla, 
      by = ~loc_2, 
      design = svy_design, 
      FUN = svymean, 
      na.rm = TRUE
      )
  
  # svyby outputs a column named "VariableLevel" (e.g., "work_33_dicYes")
  yes_col <- paste0(group_var, "Yes")
  
  wht <- 
    
    wht_raw |>
    select(Country = loc_2, Weighted = all_of(yes_col)) |>
    mutate(Weighted = Weighted * 100)
  
  # Join them together and tag the exposure name
  left_join(unw, wht, by = "Country") |>
    mutate(var = group_var)
}

# Unweighted country prevalences ------------------------------------------

# NOTE: This functions use the functions calc_unweighted_row_out and 
# calc_unweighted_row_exp, as it uses the same logic. 
# Be sure to run the full script or run such function before the following.

# Outcomes
# Modified Helper Function 
get_country_stats_out <- 
  function(subset_df, outcome) {
    
    # Overall Row (Global average for this subset)
    overall_val <- calc_unweighted_row_out(subset_df, outcome)
    
    overall_row <- tibble(
      Country = "Overall",
      Prevalence = overall_val
    )
    
    # Country Rows (Grouped by loc_2)
    country_rows <- 
      
      subset_df |>
      filter(!is.na(loc_2)) |> # Remove missing countries
      group_by(Country = loc_2) |>
      summarise(
        Prevalence = calc_unweighted_row_out(pick(everything()), outcome)
      ) |>
      select(Country, Prevalence)
    
    # Combine
    bind_rows(overall_row, country_rows)
  }

# Main Function to Build Flextable
generate_country_table_flex_out <- 
  function(data, outcome_var, outcome_title) {
    
    # Subsets (Same as before)
    d_m_doc <- 
      
      data |> 
      filter(scdm_2_rec == "Male", 
             work_2 == "Doctor")
    
    d_m_nur <- 
      
      data |> 
      filter(scdm_2_rec == "Male", 
             work_2 == "Nurse")
    
    d_f_doc <- 
      
      data |> 
      filter(scdm_2_rec == "Female", 
             work_2 == "Doctor")
    
    d_f_nur <-
      
      data |> 
      filter(scdm_2_rec == "Female", 
             work_2 == "Nurse")
    
    d_all_doc <- 
      
      data |> 
      filter(work_2 == "Doctor")
    
    d_all_nur <- 
      
      data |> 
      filter(work_2 == "Nurse")
    
    # Ns for headers (Same as before)  
    
    get_n <- 
      function(df) {
        
        df |> 
          filter(!is.na(.data[[outcome_var]])) |> 
          nrow() |> 
          style_number()
      }
    
    n_m_doc <- get_n(d_m_doc)
    n_m_nur <- get_n(d_m_nur)
    n_f_doc <- get_n(d_f_doc)
    n_f_nur <- get_n(d_f_nur)
    n_all_doc <- get_n(d_all_doc)
    n_all_nur <- get_n(d_all_nur)
    
    # Calculate Stats (Using NEW helper function)
    res_m_doc <- 
      get_country_stats_out(d_m_doc, outcome_var) |> 
      rename(M_Doc = Prevalence)
    
    res_m_nur <- 
      get_country_stats_out(d_m_nur, outcome_var) |> 
      rename(M_Nur = Prevalence)
    
    res_f_doc <- 
      get_country_stats_out(d_f_doc, outcome_var) |> 
      rename(F_Doc = Prevalence)
    
    res_f_nur <- 
      get_country_stats_out(d_f_nur, outcome_var) |> 
      rename(F_Nur = Prevalence)
    
    res_all_doc <- 
      get_country_stats_out(d_all_doc, outcome_var) |> 
      rename(All_Doc = Prevalence)
    
    res_all_nur <- 
      get_country_stats_out(d_all_nur, outcome_var) |> 
      rename(All_Nur = Prevalence)
    
    # Merge
    # Note: We join only on "Country" now, as "Variable" is not needed
    final_df <- 
      
      res_m_doc |>
      left_join(res_f_doc, by = "Country") |>
      left_join(res_all_doc, by = "Country") |>
      left_join(res_m_nur, by = "Country") |>
      left_join(res_f_nur, by = "Country") |> 
      left_join(res_all_nur, by = "Country") 
    
    # Flextable
    
    # Define labels with Ns for the second header row
    lbl_m_doc <- glue("Male\n(N = {n_m_doc})")
    lbl_f_doc <- glue("Female\n(N = {n_f_doc})")
    lbl_all_doc <- glue("Overall\n(N = {n_all_doc})")
    
    lbl_m_nur <- glue("Male\n(N = {n_m_nur})")
    lbl_f_nur <- glue("Female\n(N = {n_f_nur})")
    lbl_all_nur <- glue("Overall\n(N = {n_all_nur})")
    
    ft <- 
      
      final_df |> 
      flextable() |> 
      # Rename the columns
      set_header_labels(
        Country = "Country",
        M_Doc = lbl_m_doc,
        F_Doc = lbl_f_doc,
        All_Doc = lbl_all_doc,
        M_Nur = lbl_m_nur,
        F_Nur = lbl_f_nur,
        All_Nur = lbl_all_nur
      ) |> 
      # 2. Add the top header row (The Spanners)
      add_header_row(
        values = c("", "Doctor", "Doctor", "Doctor", "Nurse", "Nurse", "Nurse")
      ) |> 
      # 3. Merge the spanners horizontally
      merge_h(part = "header") |> 
      # 4. Styling
      theme_booktabs() |>  # A clean, professional look
      align(align = "center", part = "all") |> # Center everything
      align(j = 1, align = "left", part = "all") |> # Keep Country names left-aligned
      bold(part = "header") |> 
      fontsize(part = "all", size = 10) |> 
      autofit() |> 
      add_footer_lines(values = "Unweighted Prevalence (95% CI).")
    
    return(ft)
  }


# Exposures
# Modified Helper Function 
get_country_stats_exp <- 
  function(subset_df, exposure) {
    
    # Overall Row (Global average for this subset)
    overall_val <- calc_unweighted_row_exp(subset_df, exposure)
    
    overall_row <- tibble(
      Country = "Overall",
      Prevalence = overall_val
    )
    
    # Country Rows (Grouped by loc_2)
    country_rows <- 
      
      subset_df |>
      filter(!is.na(loc_2)) |> # Remove missing countries
      group_by(Country = loc_2) |>
      summarise(
        Prevalence = calc_unweighted_row_exp(pick(everything()), exposure)
      ) |>
      select(Country, Prevalence)
    
    # Combine
    bind_rows(overall_row, country_rows)
  }

# Main Function to Build Flextable
generate_country_table_flex_exp <- 
  function(data, exposure_var, exposure_title) {
    
    # Subsets (Same as before)
    d_m_doc <- 
      
      data |> 
      filter(scdm_2_rec == "Male", 
             work_2 == "Doctor")
    
    d_m_nur <- 
      
      data |> 
      filter(scdm_2_rec == "Male", 
             work_2 == "Nurse")
    
    d_f_doc <- 
      
      data |> 
      filter(scdm_2_rec == "Female", 
             work_2 == "Doctor")
    
    d_f_nur <-
      
      data |> 
      filter(scdm_2_rec == "Female", 
             work_2 == "Nurse")
    
    d_all_doc <- 
      
      data |> 
      filter(work_2 == "Doctor")
    
    d_all_nur <- 
      
      data |> 
      filter(work_2 == "Nurse")
    
    # Ns for headers (Same as before)  
    
    get_n <- 
      function(df) {
        
        df |> 
          filter(!is.na(.data[[exposure_var]])) |> 
          nrow() |> 
          style_number()
      }
    
    n_m_doc <- get_n(d_m_doc)
    n_m_nur <- get_n(d_m_nur)
    n_f_doc <- get_n(d_f_doc)
    n_f_nur <- get_n(d_f_nur)
    n_all_doc <- get_n(d_all_doc)
    n_all_nur <- get_n(d_all_nur)
    
    # --- Calculate Stats (Using NEW helper function) ---
    res_m_doc <- 
      get_country_stats_exp(d_m_doc, exposure_var) |> 
      rename(M_Doc = Prevalence)
    
    res_m_nur <- 
      get_country_stats_exp(d_m_nur, exposure_var) |> 
      rename(M_Nur = Prevalence)
    
    res_f_doc <- 
      get_country_stats_exp(d_f_doc, exposure_var) |> 
      rename(F_Doc = Prevalence)
    
    res_f_nur <- 
      get_country_stats_exp(d_f_nur, exposure_var) |> 
      rename(F_Nur = Prevalence)
    
    res_all_doc <- 
      get_country_stats_exp(d_all_doc, exposure_var) |> 
      rename(All_Doc = Prevalence)
    
    res_all_nur <- 
      get_country_stats_exp(d_all_nur, exposure_var) |> 
      rename(All_Nur = Prevalence)
    
    # Merge
    # Note: We join only on "Country" now, as "Variable" is not needed
    final_df <- 
      
      res_m_doc |>
      left_join(res_f_doc, by = "Country") |>
      left_join(res_all_doc, by = "Country") |>
      left_join(res_m_nur, by = "Country") |>
      left_join(res_f_nur, by = "Country") |> 
      left_join(res_all_nur, by = "Country") 
    
    # Flextable
    
    # Define labels with Ns for the second header row
    lbl_m_doc <- glue("Male\n(N = {n_m_doc})")
    lbl_f_doc <- glue("Female\n(N = {n_f_doc})")
    lbl_all_doc <- glue("Overall\n(N = {n_all_doc})")
    
    lbl_m_nur <- glue("Male\n(N = {n_m_nur})")
    lbl_f_nur <- glue("Female\n(N = {n_f_nur})")
    lbl_all_nur <- glue("Overall\n(N = {n_all_nur})")
    
    ft <- 
      
      final_df |> 
      flextable() |> 
      # Rename the columns
      set_header_labels(
        Country = "Country",
        M_Doc = lbl_m_doc,
        F_Doc = lbl_f_doc,
        All_Doc = lbl_all_doc,
        M_Nur = lbl_m_nur,
        F_Nur = lbl_f_nur,
        All_Nur = lbl_all_nur
      ) |> 
      # Add the top header row (The Spanners)
      add_header_row(
        values = c("", "Doctor", "Doctor", "Doctor", "Nurse", "Nurse", "Nurse")
      ) |> 
      # Merge the spanners horizontally
      merge_h(part = "header") |> 
      # Styling
      theme_booktabs() |>  # A clean, professional look
      align(align = "center", part = "all") |> # Center everything
      align(j = 1, align = "left", part = "all") |> # Keep Country names left-aligned
      bold(part = "header") |> 
      fontsize(part = "all", size = 10) |> 
      autofit() |> 
      add_footer_lines(values = "Unweighted Prevalence (95% CI).")
    
    return(ft)
  }

# All exposures/outcomes by country, unstratified


# Calculate Prevalence String for a single vector
# Returns a single string: "45.2% (40.1-50.5)"
calc_prev_string <- 
  function(vec, level = "Yes") {
    
    # Remove NAs from the vector
    vec <- vec[!is.na(vec)]
    n_total <- length(vec)
    
    # Handle empty data
    if (n_total == 0) return("-")
    
    # Count positives
    n_yes <- sum(vec == level)
    
    # CI calculation (using prop.test)
    ptest <- suppressWarnings(prop.test(n_yes, n_total, conf.level = 0.95))
    
    est <- ptest$estimate * 100
    ci_low <- ptest$conf.int[1] * 100
    ci_high <- ptest$conf.int[2] * 100
    
    # Format
    sprintf("%.1f%% (%.1f-%.1f)", est, ci_low, ci_high) 
  }

#  Main Function: Rows = Countries, Cols = Exposures
gen_country_tbl_flex <- 
  function(data, col_string) {
    
    # Prepare Column Labels
    # Create a named vector (Variable Name = Variable Label) for the table headers
    # This preserves your custom labels if they exist
    labels <- map_chr(col_string, function(v) {
      lbl <- tryCatch(attr(data[[v]], "label"), error = function(e) v)
      if (is.null(lbl)) return(v)
      return(lbl)
    })
    names(labels) <- col_string
    
    # Calculate "Overall" Row
    overall_row <- 
      
      data |> 
      summarise(
        Country = "Overall",
        across(
          all_of(col_string), 
          ~calc_prev_string(.x))
      )
    
    # Calculate Country Rows
    country_rows <- 
      
      data |>
      filter(!is.na(loc_2)) |> # Remove missing countries
      group_by(Country = loc_2) |>
      summarise(
        across(all_of(col_string), 
               ~calc_prev_string(.x))
      )
    
    # Combine  
    final_df <- bind_rows(overall_row, country_rows)
    
    # Create Flextable
    ft <- 
      
      final_df |>
      flextable() |>
      # Rename columns using the labels we extracted earlier
      set_header_labels(values = labels) |>
      # Formatting
      theme_booktabs() |>
      align(align = "center", part = "all") |>
      align(j = "Country", align = "left", part = "all") |> # Keep country names left-aligned
      bold(part = "header") |>
      fontsize(part = "all", size = 10) |>
      # 3. Add footer
      add_footer_lines("Values are Unweighted Prevalence (95% CI).") |> 
      autofit()
    
    return(ft)
  }

# Exposure or outcome prevalence for each country, stratified

prevalence_by_loc_tbl <- 
  function(var, data) {
    
    # Helper to calculate % and 95CI
    get_ci_str <- function(num, denom) {
      if (denom == 0) return("-")
    
    
    # binom.test gets exact CI safely for proportions
    res <- binom.test(num, denom)
    
    pct <- (num / denom) * 100
    ci_l <- res$conf.int[1] * 100
    ci_u <- res$conf.int[2] * 100
    
    sprintf("%.1f%% (%.1f-%.1f)", pct, ci_l, ci_u)
    }
    
    # Filter out NAs for the specific exposure so denominators are accurate
    ds_clean <- 
      
      data |>
      drop_na(loc_2, work_2, scdm_2_rec, all_of(var))
    
    # Calculate Subgroups (Female / Male)
    gender_stats <- 
      
      ds_clean |>
      group_by(loc_2, work_2, scdm_2_rec) |>
      summarise(
        n_denom = n(),
        n_num = sum(.data[[var]] == "Yes"),
        .groups = "drop"
      ) |>
      rename(Gender = scdm_2_rec) |>
      mutate(Gender = as.character(Gender))
    
    # Calculate Overall (By Profession only)
    overall_stats <- 
      
      ds_clean |>
      group_by(loc_2, work_2) |>
      summarise(
        n_denom = n(),
        n_num = sum(.data[[var]] == "Yes"),
        .groups = "drop"
      ) |>
      mutate(Gender = "Overall")
    
    # Combine, format the math, and pivot
    tbl_data <- 
      bind_rows(overall_stats, gender_stats) |>
      mutate(
        cell_val = purrr::map2_chr(n_num, n_denom, get_ci_str)
      ) |>
      # Create the column names that flextable will split later
      mutate(
        col_name = factor(
          paste(work_2, Gender, sep = "_"),
          levels = c("Doctor_Overall", "Doctor_Female", "Doctor_Male",
                     "Nurse_Overall", "Nurse_Female", "Nurse_Male")
        )
      ) |>
      select(loc_2, col_name, cell_val) |>
      # Pivot so Countries are rows
      pivot_wider(
        names_from = col_name,
        values_from = cell_val,
        values_fill = "-" # Drops a clean dash if a country has 0 obs
      ) |>
      rename(Country = loc_2) |>
      arrange(Country) |> 
      select(Country, Doctor_Overall, Doctor_Female, Doctor_Male, Nurse_Overall, Nurse_Female, Nurse_Male)
    
    # Build the flextable
    ft <- 
      
      flextable(tbl_data) |>
      # Magic function splits "Doctor_Overall" into grouped top headers!
      separate_header(split = "_") |> 
      theme_booktabs() |>
      align(j = 1, align = "left", part = "all") |>
      align(j = -1, align = "center", part = "all") |>
      # Add your footnote
      add_footer_lines(values = paste0("% (95% CI) for ", var)) |>
      autofit()
    
    return(ft)
  }


run_out_by_n_bin <- 
  
  function(data, outcomes) {
    
    # Helper: fit model and extract OR + CI
    fit_one <- function(outcome, predictor) {
      f <- as.formula(paste0(outcome, " ~ ", predictor))
      m <- glm(f, data = data, family = binomial(link = "log"))
      
      broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) |>
        filter(term != "(Intercept)") |> 
        mutate(outcome = outcome,
               predictor = predictor) |>
        select(outcome, predictor, estimate, conf.low, conf.high)
    }
    
    # All combinations: 4 outcomes × 3 predictors
    predictors <- c("n", "N", "rates_100")
    
    map_df(outcomes, \(o)
           map_df(predictors, \(p) fit_one(o, p))
    ) |>
      # reshape to wide format: rows = predictors, columns = outcomes
      mutate(label = paste0(
        sprintf("%.2f", estimate), 
        " (", sprintf("%.2f", conf.low), "–", sprintf("%.2f", conf.high), ")"
      )) |>
      select(outcome, predictor, label) |>
      tidyr::pivot_wider(
        names_from = outcome,
        values_from = label
      ) |>
      rename(model = predictor)
  }
# Sensitivity analyses ----------------------------------------------------

run_out_by_n_bin <- 
  
  function(data, outcomes) {
    
    # Helper: fit model and extract OR + CI
    fit_one <- function(outcome, predictor) {
      f <- as.formula(paste0(outcome, " ~ ", predictor))
      m <- glm(f, data = data, family = binomial(link = "log"))
      
      broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) |>
        filter(term != "(Intercept)") |> 
        mutate(outcome = outcome,
               predictor = predictor) |>
        select(outcome, predictor, estimate, conf.low, conf.high)
    }
    
    # All combinations: 4 outcomes × 3 predictors
    predictors <- c("n", "N", "rates_100")
    
    map_df(outcomes, \(o)
           map_df(predictors, \(p) fit_one(o, p))
    ) |>
      # reshape to wide format: rows = predictors, columns = outcomes
      mutate(label = paste0(
        sprintf("%.2f", estimate), 
        " (", sprintf("%.2f", conf.low), "–", sprintf("%.2f", conf.high), ")"
      )) |>
      select(outcome, predictor, label) |>
      tidyr::pivot_wider(
        names_from = outcome,
        values_from = label
      ) |>
      rename(model = predictor)
  }

run_out_by_n_poi <- function(data, outcomes){

  fit_one <- function(outcome, predictor) {
    f <- as.formula(paste0(outcome, " ~ ", predictor))
    m <- glm(f, data = data, family = poisson(link = "log"))
    
    # robust vcov
    V <- sandwich::vcovHC(m, type = "HC0")
    ct <- lmtest::coeftest(m, vcov = V)
    
    # extract coefficient row
    est  <- ct[2, "Estimate"]      # log(PR)
    se   <- ct[2, "Std. Error"]    # robust SE
    lcl  <- est - 1.96 * se        # CI on log scale
    ucl  <- est + 1.96 * se
    
    tibble(
      outcome = outcome,
      predictor = predictor,
      estimate = exp(est),
      conf.low = exp(lcl),
      conf.high = exp(ucl)
    )
  }
  
  predictors <- c("n", "N", "rates_100")
  
  map_df(outcomes, \(o)
         map_df(predictors, \(p) fit_one(o, p))
  ) |>
    mutate(label = paste0(
      sprintf("%.2f", estimate),
      " (", sprintf("%.2f", conf.low), "–", sprintf("%.2f", conf.high), ")"
    )) |>
    select(outcome, predictor, label) |>
    tidyr::pivot_wider(
      names_from = outcome,
      values_from = label,
      names_glue = "{.value}_{outcome}"
    ) |>
    rename(model = predictor)
}
