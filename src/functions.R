library(gt)
library(gtsummary)
library(glue)
library(meta)
library(lme4)
library(rlang)
library(emmeans)
library(broom.mixed)

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
  function(data, exposure_list, scdm) {

    # Loop through exposures and stack data
    long_df <- map_dfr(names(exposure_list), function(var_name) {
      title <- exposure_list[[var_name]]
      get_exposure_data_long(data, var_name, title, scdm)
    })


    # Prepare for gt: Insert "Header Rows" for Variables
    # We want a row that says "Age Group" followed by rows "<30", "30-50"
    formatted_df <- long_df |>
      group_by(Exposure_Label, Variable) |>
      group_split() |>
      map_dfr(function(chunk) {
        # Create a dummy header row
        header_row <- tibble(
          Exposure_Label = unique(chunk$Exposure_Label),
          Variable = unique(chunk$Variable),
          Category = unique(chunk$Variable), # The Category column becomes the Header text
          M_Doc = NA, F_Doc = NA, All_Doc = NA,
          M_Nur = NA, F_Nur = NA, All_Nur = NA
        )

        # Determine if this is the "Overall" row (optional: skip header for Overall if desired)
        if(unique(chunk$Variable) == "Overall") {
          return(chunk)
        } else {
          return(bind_rows(header_row, chunk))
        }
      })

    # Render with gt
    formatted_df |>

      gt(
        groupname_col = "Exposure_Label"
      ) |> # Group by the Exposure (Depression, etc)
      cols_label(
        Category = md("**Subgroup**"),
        M_Doc = md("**Male**"),
        F_Doc = md("**Female**"),
        All_Doc = md("**Overall**"),
        M_Nur = md("**Male**"),
        F_Nur = md("**Female**"),
        All_Nur = md("**Overall**")
      ) |>
      tab_spanner(
        label = md("**Doctor**"),
        columns = c(M_Doc, F_Doc, All_Doc)
      ) |>
      tab_spanner(
        label = md("**Nurse**"),
        columns = c(M_Nur, F_Nur, All_Nur)
      ) |>
      cols_align(
        align = "center",
        columns = contains("_")
      ) |>
      sub_missing(missing_text = "") |>

      # Make the inserted Variable Headers (where Category == Variable) Bold
      tab_style(
        style = list(cell_text(weight = "bold")),
        locations = cells_body(
          columns = Category,
          rows = Category == Variable & Variable != "Overall"
        )
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_row_groups()
      )|>

      # Indent the actual categories (rows where Category != Variable)
      tab_style(
        style = list(cell_text(indent = px(15))),
        locations = cells_body(
          columns = Category,
          rows = Category != Variable & Variable != "Overall"
        )
      ) |>

      #  Clean up: Hide the helper "Variable" column
      cols_hide(columns = Variable) |>

      tab_header(
        title = md("**Exposures by Sociodemographics (Unweighted)**")
      ) |>
      tab_footnote(
        footnote = "Unweighted Prevalence (95% CI). N represents total sample for that exposure."
      )
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
  function(data, outcome_list, scdm) {

    # Loop through exposures and stack data
    long_df <- map_dfr(names(outcome_list), function(var_name) {
      title <- outcome_list[[var_name]]
      get_outcome_data_long(data, var_name, title, scdm)
    })
    

    # Prepare for gt: Insert "Header Rows" for Variables
    # We want a row that says "Age Group" followed by rows "<30", "30-50"
    formatted_df <- long_df |>
      group_by(Outcome_Label, Variable) |>
      group_split() |>
      map_dfr(function(chunk) {
        # Create a dummy header row
        header_row <- tibble(
          Outcome_Label = unique(chunk$Outcome_Label),
          Variable = unique(chunk$Variable),
          Category = unique(chunk$Variable), # The Category column becomes the Header text
          M_Doc = NA, F_Doc = NA, All_Doc = NA,
          M_Nur = NA, F_Nur = NA, All_Nur = NA
        )

        # Determine if this is the "Overall" row (optional: skip header for Overall if desired)
        if(unique(chunk$Variable) == "Overall") {
          return(chunk)
        } else {
          return(bind_rows(header_row, chunk))
        }
      })

    # Render with gt
    formatted_df |>

      gt(
        groupname_col = "Outcome_Label"
        ) |> # Group by the Outcome (Depression, etc)
      cols_label(
        Category = md("**Subgroup**"),
        M_Doc = md("**Male**"),
        F_Doc = md("**Female**"),
        All_Doc = md("**Overall**"),
        M_Nur = md("**Male**"),
        F_Nur = md("**Female**"),
        All_Nur = md("**Overall**")
      ) |>
      tab_spanner(
        label = md("**Doctor**"),
        columns = c(M_Doc, F_Doc, All_Doc)
        ) |>
      tab_spanner(
        label = md("**Nurse**"),
        columns = c(M_Nur, F_Nur, All_Nur)
        ) |>
      cols_align(
        align = "center",
        columns = contains("_")
        ) |>
      sub_missing(missing_text = "") |>

      # Make the inserted Variable Headers (where Category == Variable) Bold
      tab_style(
        style = list(cell_text(weight = "bold")),
        locations = cells_body(
          columns = Category,
          rows = Category == Variable & Variable != "Overall"
        )
      ) |> tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    )|>

      # Indent the actual categories (rows where Category != Variable)
      tab_style(
        style = list(cell_text(indent = px(15))),
        locations = cells_body(
          columns = Category,
          rows = Category != Variable & Variable != "Overall"
        )
      ) |>

      #  Clean up: Hide the helper "Variable" column
      cols_hide(columns = Variable) |>

      tab_header(
        title = md("**Outcomes by Sociodemographics (Unweighted)**")
        ) |>
      tab_footnote(
        footnote = "Unweighted Prevalence (95% CI). N represents total sample for that outcome."
        )
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

# Unweighted associations -------------------------------------------------

# Function for GLM crude models

run_crude_glm <- 
  function(data, outcome_var, exposure_var, profession, 
           cluster_var, family_type = "poisson",
           output_dir = "out/mods") {
    
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
      
      message(glue("Loading ML model: {file_name_ml}"))
      ml <- readRDS(file_path_ml)
      
    } else {
      message(glue("Fitting ML model: {file_name_ml}"))
      
      
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
      
      message(glue("Loading Robust model: {file_name_rob}"))
      rob <- readRDS(file_path_rob)
      
    } else {
      message(glue("Fitting Robust model: {file_name_rob}"))
      
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
           output_dir = "out/mods") {
    
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
      
      message(glue("Loading ML model: {file_name_ml}"))
      ml <- tryCatch({
        readRDS(file_path_ml)
      }, error = function(e) {
        
        message(glue("Corrupted file detected! Deleting {file_name_ml} and refitting..."))
        file.remove(file_path_ml)
        return(NULL)
      })
    } 
    
    if(is.null(ml)) {
      message(glue("Fitting ML model: {file_name_ml}"))
      
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
      
      message(glue("Loading Robust SE model: {file_name_rob}"))
      rob <- tryCatch({
        readRDS(file_path_rob)
      }, error = function(e) {
        
        message(glue("Corrupted file detected! Deleting {file_name_rob} and refitting..."))
        file.remove(file_path_rob)
        return(NULL)
      })
    } 
    
    if(is.null(rob)) {
      message(glue("Fitting Robust model: {file_name_rob}"))
      
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
      Exposure = factor(Exposure, levels = exposure_order),
      Outcome = factor(Outcome, levels = c("Depression", "Anxiety", "Suicidal thoughts"))
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

run_glm_dose_emm <- function(
    data, outcome_var, exposure_var, confounder_list,
    cluster_var = "loc_2", family_type = "gaussian",
    model_label, output_dir = "out/mods_dose",
    force_fit = FALSE
) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  required_vars <- unique(c(outcome_var, exposure_var, cluster_var, confounder_list))
  
  # Clean dataset
  ds_clean <- data |> 
    select(all_of(required_vars)) |> 
    drop_na() |> 
    mutate(
      !!exposure_var := factor(
        .data[[exposure_var]],
        levels = c("No", "Yes, a few times", "Yes, monthly", "Yes, weekly", "Yes, daily"),
        ordered = TRUE
      )
    )
  
  n_size <- nrow(ds_clean)
  if (n_size < 10) return(NULL)
  
  # Short labels
  out_short <- case_when(
    outcome_var == "phq_sc"   ~ "dep",
    outcome_var == "gad_sc"   ~ "anx",
    outcome_var == "mh_phq_9" ~ "suic"
  )
  
  exp_short <- case_when(
    exposure_var == "work_30" ~ "har",
    exposure_var == "work_33" ~ "bul",
    exposure_var == "work_31" ~ "threats",
    exposure_var == "work_32" ~ "viol"
  )
  
  file_name <- paste0("dose_", out_short, "_", exp_short, "_",
                      tolower(gsub(" ", "_", model_label)), ".rds")
  file_path <- file.path(output_dir, file_name)
  
  # Fit or load model
  if(file.exists(file_path) & !force_fit){
    message(glue("Loading saved bundle: {file_name}"))
    bundle <- readRDS(file_path)
    model <- bundle$model
    balance_report <- bundle$balance
  } else {
    message(glue("Fitting Dose-response: {file_name}"))
    
    balance_report <- NULL
    
    model_formula <- as.formula(
      paste(
        outcome_var, "~",
        paste(c(exposure_var, confounder_list), collapse = " + "),
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
  emm_df <- as.data.frame(
    emmeans(model, specs = exposure_var, type = "response", df = NULL)
  ) %>%
    rename(estimate = emmean, std.error = SE) %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      Level = as.character(.data[[exposure_var]]),
      Outcome = outcome_var,
      Exposure = exposure_var,
      Model = model_label,
      N_Analysis = n_size
    ) %>%
    select(Outcome, Exposure, Level, estimate, conf.low, conf.high, Model, N_Analysis)
  
  
  if(is.null(emm_df) || nrow(emm_df) == 0){
    message(glue("Warning: No EMMs found for {exposure_var}"))
    return(NULL)
  }
  
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
      
      message(glue("Loading from disk: {name_crude}"))
      m_crude <- readRDS(file_crude)
      
    } else {
      
      message(glue("Fitting new model: {name_crude}"))
      
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
      
      message(glue("Loading from disk: {name_adj}"))
      m_adj <- readRDS(file_adj)
      
    } else {
      
      message(glue("Fitting new model: {name_adj}"))
      
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
    model_label, output_dir = "out/mods_moder",
    force_fit = FALSE
) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  required_vars <- 
    unique(
      c(outcome_var, exposure_var, cluster_var, confounder_list, moderator_var)
      )
  
  # Clean dataset
  ds_clean <- data |> 
    select(all_of(required_vars)) |> 
    drop_na() |> 
    mutate(
      !!exposure_var := factor(
        .data[[exposure_var]],
        levels = c("No", "Yes, a few times", "Yes, monthly", "Yes, weekly", "Yes, daily"),
        ordered = TRUE
      )
    )
  
  n_size <- nrow(ds_clean)
  if (n_size < 10) return(NULL)
  
  # Short labels
  out_short <- case_when(
    outcome_var == "phq_sc"   ~ "dep",
    outcome_var == "gad_sc"   ~ "anx",
    outcome_var == "mh_phq_9" ~ "suic"
  )
  
  exp_short <- case_when(
    exposure_var == "work_30" ~ "har",
    exposure_var == "work_33" ~ "bul",
    exposure_var == "work_31" ~ "threats",
    exposure_var == "work_32" ~ "viol"
  )
  
  moder_short <- case_when(
    moderator_var == "work_28_dic" ~ "bull_harass",
    moderator_var == "work_29_dic" ~ "threats_abuse_assault",
    
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
    message(glue("Loading saved bundle: {file_name}"))
    bundle <- readRDS(file_path)
    model <- bundle$model
    balance_report <- bundle$balance
  } else {
    message(glue("Fitting moderation: {file_name}"))
    
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
  emm_df <- as.data.frame(
    emmeans(model, 
            specs = c(exposure_var, moderator_var))
  ) %>%
    rename(estimate = emmean, std.error = SE) %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      Level = as.character(.data[[exposure_var]]),
      Moderator_Level = as.character(.data[[moderator_var]]),
      Outcome = outcome_var,
      Exposure = exposure_var,
      Moderator = moderator_var,
      Model = model_label,
      N_Analysis = n_size
    ) %>%
    select(Outcome, Exposure, Moderator, Level, Moderator_Level, estimate, conf.low, conf.high, Model, N_Analysis)
  
  
  if(is.null(emm_df) || nrow(emm_df) == 0){
    message(glue("Warning: No EMMs found for {exposure_var}"))
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
           output_dir = "out/mods") {
    
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
      
      message(glue("Loading ML model: {file_name_ml}"))
      ml <- readRDS(file_path_ml)
      
    } else {
      message(glue("Fitting ML model: {file_name_ml}"))
      
      
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
      
      message(glue("Loading Robust SE model: {file_name_rob}"))
      rob <- readRDS(file_path_rob)
      
    } else {
      message(glue("Fitting Robust SE model: {file_name_rob}"))
      
      
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

#

run_adj_moder_glm <- 
  
function(data, outcome_var, exposure_var, confounder_list, moderator, moder_val,
         cluster_var = "loc_2", family_type = "poisson", model_label,
         output_dir = "out/moder-dic/") {
  
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
    outcome_var == "suic_idea_poi"   ~ "suic"
  )
  
  # exposure abbreviation
  exp_short <- case_when( 
    exposure_var == "work_30_dic"   ~ "har",
    exposure_var == "work_33_dic"   ~ "bul",
    exposure_var == "work_31_dic"   ~ "threats",
    exposure_var == "work_32_dic"   ~ "viol"
  )
  
  
  # moderator abbreviation
  moder_short <- case_when(
    moderator == "work_28_dic" ~ "prot-har",
    moderator == "work_29_dic"  ~ "prot-vio"
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
    
    message(glue("Loading ML model: {file_name_ml}"))
    ml <- readRDS(file_path_ml)
    
  } else {
    message(glue("Fitting ML model: {file_name_ml}"))
    
    
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
    
    message(glue("Loading Robust SE model: {file_name_rob}"))
    rob <- readRDS(file_path_rob)
    
  } else {
    message(glue("Fitting Robust SE model: {file_name_rob}"))
    
    
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
   