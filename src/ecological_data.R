library(readr)
library(tidyverse)


# Load data ---------------------------------------------------------------


raw_doctors <- # n practising/active/licensed doctors by country
  
  readxl::read_xlsx("data/All_Employment_2025_v6 for prefills.xlsx",
                    sheet = "Physicians", col_names = FALSE)

raw_nurses <- # n practising/active/licensed nurses by country
  
  readxl::read_xlsx("data/All_Employment_2025_v6 for prefills.xlsx",
                    sheet = "Nurses", col_names = FALSE)

who_phys_age_sex_loc <- # doctors by sex, age, year, and country
  
  readxl::read_xlsx("data/All_Employment_2025_v6 for prefills.xlsx",
                    sheet = "Doctors_by_age")

who_nurs_age_sex_loc <- # nurses by sex, age, year, and country
  
  readxl::read_xlsx("data/All_Employment_2025_v6 for prefills.xlsx",
                    sheet = "prof_nurses_by_agegender")

# Cleaning data -----------------------------------------------------------


# Doctors
col_names <- c("loc_2", 
               "year", 
               "practising", 
               "active", 
               "licensed")

cleaned_data <- 
  
  raw_doctors |> 
  slice(10:n()) |> # data starts at row 9
  select(1, 2, 3, 5, 7) |> 
  setNames(col_names) |> 
  filter(!is.na(loc_2)) |> 
  mutate(across(!loc_2, as.integer),
         loc_2 = as.factor(loc_2),
         loc_2 = if_else(loc_2 == "Estland", "Estonia", loc_2),
         loc_2 = if_else(loc_2 == "Czech Republic", "Czechia", loc_2),
         loc_2 = if_else(loc_2 == "Slovak Republic", "Slovakia", loc_2)
         ) |> 
  filter(loc_2 %in% ds$loc_2)



cleaned_data <-   # gets most recent data on each category across countries
  
  cleaned_data |> 
  pivot_longer(practising:licensed,
               names_to = "category") |> 
  drop_na(value) |> 
  group_by(loc_2, category) |> 
  filter(year == max(year)) |> 
  ungroup()

who_phys_loc <- # prioritisises categories
  
  cleaned_data |> 
  mutate(category = factor(category, 
                           levels = # order matters
                             c("practising", "active", "licensed"))) |> 
  arrange(loc_2, category) |> 
  group_by(loc_2) |> 
  slice(1) |> 
  ungroup()

# Nurses
col_names <- c("loc_2", 
               "year", 
               "practising", 
               "prof_act_total",
               "active", 
               "licensed")

cleaned_data <- 
  
  raw_nurses |> 
  select(1, 2, 5, 9, 11, 17) |> # professional nurses only
  slice(9:n()) |> # data starts at row 9 
  setNames(col_names)  |> 
  filter(!is.na(loc_2)) |> 
  mutate(across(!loc_2, as.integer),
         loc_2 = if_else(loc_2 == "Estland", "Estonia", loc_2),
         loc_2 = if_else(loc_2 == "Czech Republic", "Czechia", loc_2),
         loc_2 = if_else(loc_2 == "Slovak Republic", "Slovakia", loc_2)
         ) |>  
  filter(loc_2 %in% ds$loc_2)

rm(col_names)

cleaned_data <-   # gets most recent data on each category across countries
  
  cleaned_data |> 
  pivot_longer(practising:licensed,
               names_to = "category") |> 
  drop_na(value) |> 
  group_by(loc_2, category) |> 
  filter(year == max(year)) |> 
  ungroup()

who_nurse_loc <- # prioritisises categories
  
  cleaned_data |> 
  mutate(category = factor(category, 
                           levels = # order matters
                             c("practising", "active", "licensed", "prof_act_total"))) |> 
  arrange(loc_2, category) |> 
  group_by(loc_2) |> 
  slice(1)

rm(cleaned_data)

# Merge

who_n_loc <- 
  
  bind_rows(who_phys_loc, who_nurse_loc, .id = "work_2") |> 
  mutate(work_2 = if_else(work_2 == "1", "Doctor", "Nurse")) |> 
  mutate(across(where(is.character), as.factor)) |> 
  tidyr::complete(
    loc_2,
    work_2,
    fill = list(N = NA)
  )

writexl::write_xlsx(
  who_n_loc,
  path = "ext/who_n_loc.xlsx"
  )

saveRDS(who_n_loc, "ext/who_n_loc.rds")

rm(raw_doctors, raw_nurses, who_phys_loc, who_nurse_loc)

# Disagreggated
# Doctors
col_names <- c("loc_2",
               "year",
               "fem_tot",
               "fem_less_35",
               "fem_35-44",
               "fem_45-54",
               "fem_55-64",
               "fem_65-74",
               "fem_75_over",
               "male_tot",
               "male_less_35",
               "male_35-44",
               "male_45-54",
               "male_55-64",
               "male_65-74",
               "male_75_over",
               "all_tot",
               "all_less_35",
               "all_35-44",
               "all_45-54",
               "all_55-64",
               "all_65-74",
               "all_75_over"
               )

who_phys_age_sex_loc <- # for IPW
  
  who_phys_age_sex_loc |>
  janitor::clean_names() |>
  select(1,2,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43) |> 
  slice(9:n()) |> 
  setNames(col_names) |> 
  mutate(loc_2 = if_else(loc_2 == "Estland", "Estonia", loc_2),
         loc_2 = if_else(loc_2 == "Czech Republic", "Czechia", loc_2),
         loc_2 = if_else(loc_2 == "Slovak Republic", "Slovakia", loc_2)
         ) |>
  filter(loc_2 %in% ds$loc_2) |> 
  pivot_longer(
    cols = !c(1, 2),
    names_to = c("scdm_2_rec", "age_eurostat"),
    names_pattern = "(fem|male|all)_(.*)",
    values_to = "N"
  ) |>
  mutate(
    scdm_2_rec = case_when(
      scdm_2_rec == "fem" ~ "Female",
      scdm_2_rec == "male" ~ "Male",
      scdm_2_rec == "all" ~ "All"),
    work_2 = "Doctor",
    age_eurostat = case_when(
      age_eurostat == "tot" ~ "All",
      age_eurostat == "less_35" ~ "Less than 35 years",
      age_eurostat == "35-44" ~ "From 35 to 44 years",
      age_eurostat == "45-54" ~ "From 45 to 54 years",
      age_eurostat == "55-64" ~ "From 55 to 64 years",
      age_eurostat == "65-74" ~ "From 65 to 74 years",
      age_eurostat == "75_over" ~ "75 years or over",
      TRUE ~ as.character(age_eurostat)
    )
  ) |>
  select(
    N_year = year,
    loc_2,
    work_2,
    scdm_2_rec,
    age_eurostat,
    N
  ) |>
  # keeps latest year with data available
  mutate(N = as.numeric(N)) |> 
  group_by(loc_2) |>
  filter(N_year == max(N_year)) |>
  ungroup()


# Nurses
col_names <- c("loc_2",
               "year",
               "fem_tot",
               "fem_less_25",
               "fem_25-34",
               "fem_35-44",
               "fem_45-54",
               "fem_55-64",
               "fem_65-74",
               "fem_75_over",
               "male_tot",
               "male_less_25",
               "male_25-34",
               "male_35-44",
               "male_45-54",
               "male_55-64",
               "male_65-74",
               "male_75_over",
               "all_tot",
               "all_less_25",
               "all_25-34",
               "all_35-44",
               "all_45-54",
               "all_55-64",
               "all_65-74",
               "all_75_over"
               )


who_nurs_age_sex_loc <-
  
  who_nurs_age_sex_loc |>
  janitor::clean_names() |>
  select(1,2,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49) |> 
  slice(10:n()) |> 
  setNames(col_names) |> 
  mutate(loc_2 = if_else(loc_2 == "Estland", "Estonia", loc_2),
           loc_2 = if_else(loc_2 == "Czech Republic", "Czechia", loc_2),
           loc_2 = if_else(loc_2 == "Slovak Republic", "Slovakia", loc_2)
         ) |>
  filter(loc_2 %in% ds$loc_2) |>
  mutate(
    across(
      all_of(3:26),
      \(x) as.numeric(x)
    )
  ) |> 
  mutate(
    # doctors and nurses have different age ranges in WHO and EUROSTAT
    
    fem_less_35 = fem_less_25 + `fem_25-34`,
    male_less_35 = male_less_25 + `male_25-34`,
    all_less_35 = all_less_25 + `all_25-34`
  ) |>
  select(-contains("_25")) |>
  relocate(fem_less_35, .after = fem_tot) |> 
  relocate(male_less_35, .after = male_tot) |> 
  relocate(all_less_35, .after = all_tot) |>  
  pivot_longer(
    cols = !c(1, 2),
    names_to = c("scdm_2_rec", "age_eurostat"),
    names_pattern = "(fem|male|all)_(.*)",
    values_to = "N"
  ) |>
  mutate(
    scdm_2_rec = case_when(
      scdm_2_rec == "fem" ~ "Female",
      scdm_2_rec == "male" ~ "Male",
      scdm_2_rec == "all" ~ "All"),
    work_2 = "Nurse",
    age_eurostat = case_when(
      age_eurostat == "tot" ~ "All",
      age_eurostat == "less_35" ~ "Less than 35 years",
      age_eurostat == "35-44" ~ "From 35 to 44 years",
      age_eurostat == "45-54" ~ "From 45 to 54 years",
      age_eurostat == "55-64" ~ "From 55 to 64 years",
      age_eurostat == "65-74" ~ "From 65 to 74 years",
      age_eurostat == "75_over" ~ "75 years or over",
      TRUE ~ as.character(age_eurostat)
    )
  ) |>
  select(
    N_year = year,
    loc_2,
    work_2,
    scdm_2_rec,
    age_eurostat,
    N
  ) |>
  # keeps latest year with data available
  
  group_by(loc_2) |>
  filter(N_year == max(N_year)) |>
  ungroup()

rm(col_names)

# Merge
who_ipw_population <- # for IPW
  
  bind_rows(who_phys_age_sex_loc, who_nurs_age_sex_loc) |>
  mutate(
    across(
      where(is.character),
      as_factor
    )
  ) |>
  # creates a dataset with all possible combinations
  
  tidyr::complete(
    loc_2,
    work_2,
    scdm_2_rec,
    age_eurostat,
    fill = list(N = NA)
  )



rm(who_nurs_age_sex_loc, who_phys_age_sex_loc)

who_ipw_population_red <- # CAUTION: percentages not very reliable due to NAs
  
  who_ipw_population |>
  filter(scdm_2_rec == "All" & age_eurostat == "All") |>
  mutate(perc_eu = round(N / sum(N, na.rm = FALSE) * 100, 2)) |>
  group_by(loc_2) |>
  mutate(perc_country = round(N / sum(N, na.rm = FALSE) * 100, 2))

who_ipw_population <-
  
  who_ipw_population |>
  filter(
    scdm_2_rec != "All" & age_eurostat != "All",
    age_eurostat != "75 years or over"
  ) |>
  mutate(age_eurostat = factor(age_eurostat,
                               levels = c(
                                 "Less than 35 years",
                                 "From 35 to 44 years",
                                 "From 45 to 54 years",
                                 "From 55 to 64 years",
                                 "From 65 to 74 years"
                               ),
                               ordered = TRUE
  )) |>
  # Overall percentage (EU level)
  
  mutate(perc_eu = round(N / sum(N) * 100, 2)) |>
  # Percentage by country
  
  group_by(loc_2) |>
  mutate(perc_country = round(N / sum(N) * 100, 2)) |>
  ungroup() |>
  # Percentage by country and profession
  
  group_by(loc_2, work_2) |>
  mutate(perc_country_work = round(N / sum(N) * 100, 2)) |>
  ungroup() |>
  # Percentage by country, profession, and sex
  
  group_by(loc_2, work_2, scdm_2_rec) |>
  mutate(perc_country_work_sex = round(N / sum(N) * 100, 2)) |>
  ungroup() |>
  # arrange for comparability
  
  arrange(loc_2, work_2, scdm_2_rec, age_eurostat)

writexl::write_xlsx(who_ipw_population,
                    "ext/who_ipw_population.xlsx")

saveRDS(who_ipw_population,
        "ext/who_ipw_population.rds")

# Compute weights for IPW -------------------------------------------------

# Create a dataset with Ns, percentages (by country), and percentages (EU),
# by country, type of job, sex, and age (recoded) (SAMPLE)

ipw_sample <-
  
  ds |>
  select(
    loc_2,
    work_2,
    scdm_2_rec,
    age_eurostat
  ) |>
  mutate(
    across(
      where(is.character), 
      as.factor),
    age_eurostat = factor(age_eurostat,
                          levels = c(
                            "Less than 35 years",
                            "From 35 to 44 years",
                            "From 45 to 54 years",
                            "From 55 to 64 years",
                            "From 65 to 74 years"
                            ),
                          ordered = TRUE
                          ),
    N_year = 2025
    ) |>
  group_by(
    N_year,
    loc_2,
    work_2,
    scdm_2_rec,
    age_eurostat
  ) |>
  summarise(N = n()) |>
  drop_na() |>
  ungroup() |>
  # Overall percentage (EU level)
  
  mutate(perc_eu = round(N / sum(N) * 100, 2)) |>
  # Percentage by country
  
  group_by(loc_2) |>
  mutate(perc_country = round(N / sum(N) * 100, 2)) |>
  ungroup() |>
  # Percentage by country and profession
  
  group_by(loc_2, work_2) |>
  mutate(perc_country_work = round(N / sum(N) * 100, 2)) |>
  ungroup() |>
  # Percentage by country, profession, and sex
  
  group_by(loc_2, work_2, scdm_2_rec) |>
  mutate(perc_country_work_sex = round(N / sum(N) * 100, 2)) |>
  ungroup() |>
  # arrange for comparability
  
  arrange(loc_2, work_2, scdm_2_rec, age_eurostat) |>
  # adds all possible combinations
  
  tidyr::complete(
    loc_2,
    work_2,
    scdm_2_rec,
    age_eurostat,
    fill = list(N = NA)
  ) 

saveRDS(ipw_sample,
        "ext/ipw_sample.rds")

who_weights <- # weights using who data
  
  ipw_sample |>
  left_join(who_ipw_population,
            by = c("loc_2", "work_2", "scdm_2_rec", "age_eurostat"),
            suffix = c("_sample", "_popul")
  ) |>
  mutate(
    w_eu =
      perc_eu_popul / perc_eu_sample,
    w_country =
      perc_country_popul / perc_country_sample,
    w_country_work =
      perc_country_work_popul / perc_country_work_sample,
    w_country_work_sex =
      perc_country_work_sex_popul / perc_country_work_sex_sample
  ) |> 
  select(
    N_year_sample,
    N_year_popul,
    loc_2,
    work_2,
    scdm_2_rec,
    age_eurostat,
    starts_with("w_")
  )

writexl::write_xlsx(who_weights, "ext/who_weights.xlsx")
saveRDS(who_weights,
        "ext/who_weights.rds")

# extract weights
w <- 
  
  ds |>
  left_join(
    who_weights |>
      select(loc_2, 
             work_2, 
             scdm_2_rec, 
             age_eurostat, 
             w_country_work_sex),
    by = c("loc_2", 
           "work_2", 
           "scdm_2_rec", 
           "age_eurostat")
  ) |> 
  select(user_id, w_country_work_sex)

saveRDS(w,
        "ext/weights.rds")
