

required_pkgs <- c("tidyverse")
future_model_pkgs <- c("lme4", "lmerTest", "geepack")

missing_required <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_required) > 0) {
  stop(
    "Install required package(s) before running this script: ",
    paste(missing_required, collapse = ", "),
    call. = FALSE
  )
}

missing_future <- setdiff(future_model_pkgs, rownames(installed.packages()))
if (length(missing_future) > 0) {
  message(
    "Note: these modeling packages will be needed in later phases but are not ",
    "required for Phase 1-2: ",
    paste(missing_future, collapse = ", ")
  )
}

library(tidyverse)

options(
  dplyr.summarise.inform = FALSE,
  readr.show_col_types = FALSE
)

set.seed(20260506)

raw_path <- "/Users/mahdibalotf/Desktop/project folder2/First step/df_women_menopause_clean_5cat_dynamic_v2.csv"
out_dir <- file.path(dirname(raw_path), "phase1_phase2_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Column map and helper functions -------------------------------------

# Keep all key variable names in one explicit map. If the raw extract changes,
# update this section first and the rest of the script should remain stable.
id_var <- "hhidpn"
wave_var <- "wave"
year_var <- "year"

age_var <- case_when(
  "agey_e" %in% names(read_csv(raw_path, n_max = 0, guess_max = Inf)) ~ "agey_e",
  "age_wave" %in% names(read_csv(raw_path, n_max = 0, guess_max = Inf)) ~ "age_wave",
  "agey_b" %in% names(read_csv(raw_path, n_max = 0, guess_max = Inf)) ~ "agey_b",
  TRUE ~ NA_character_
)

cog_var <- case_when(
  "cogtot27_imp" %in% names(read_csv(raw_path, n_max = 0, guess_max = Inf)) ~ "cogtot27_imp",
  "cog27" %in% names(read_csv(raw_path, n_max = 0, guess_max = Inf)) ~ "cog27",
  TRUE ~ NA_character_
)

cog_function_var <- "cogfunction"
female_var <- "Female"
gender_var <- "ragender"

meno_combo_var <- "meno_combo"
meno_type_var <- "meno_type"
meno_age5_var <- "meno_age5"
meno_status_var <- "meno_status_time"
meno_age_var <- "C253_age"

required_cols <- c(
  id_var, wave_var, year_var, age_var, cog_var, cog_function_var,
  female_var, gender_var, meno_combo_var, meno_type_var, meno_age5_var,
  meno_status_var, meno_age_var
)

if (any(is.na(required_cols))) {
  stop("At least one required variable could not be mapped. Check the column map.", call. = FALSE)
}

mode_nonmissing <- function(x) {
  # Deterministic mode: most frequent non-missing value; if tied, choose the
  # first tied value in the respondent's longitudinal order.
  x_nonmiss <- x[!is.na(x)]
  if (length(x_nonmiss) == 0) {
    return(x[NA_integer_][1])
  }

  tab <- table(as.character(x_nonmiss), useNA = "no")
  top_values <- names(tab)[tab == max(tab)]
  x_nonmiss[as.character(x_nonmiss) %in% top_values][1]
}

first_nonmissing <- function(x) {
  x_nonmiss <- x[!is.na(x)]
  if (length(x_nonmiss) == 0) {
    return(x[NA_integer_][1])
  }
  x_nonmiss[1]
}

safe_min <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

flag_female <- function(female, gender) {
  # HRS RAND convention: ragender == 2 is female. The derived Female variable
  # is expected to equal 1 for women when available.
  case_when(
    !is.na(female) ~ female == 1,
    is.na(female) & !is.na(gender) ~ gender == 2,
    TRUE ~ NA
  )
}

valid_meno_groups <- c(
  "1.Premature Natural",
  "2.Early Natural",
  "3.45-49 Natural",
  "4.50-55 Natural (ref)",
  "5.Late Natural",
  "6.Surgical"
)

natural_meno_groups <- c(
  "1.Premature Natural",
  "2.Early Natural",
  "3.45-49 Natural",
  "4.50-55 Natural (ref)",
  "5.Late Natural"
)

## ---- 2. Read raw data --------------------------------------------------------

# guess_max = Inf prevents readr from incorrectly guessing sparsely populated
# reproductive-history fields as logical.
df_raw <- read_csv(raw_path, guess_max = Inf)

missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop(
    "The dataset is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

parse_diagnostics <- problems(df_raw)
write_csv(parse_diagnostics, file.path(out_dir, "phase1_parse_diagnostics.csv"))

raw_inventory <- tibble(
  dataset = basename(raw_path),
  n_rows = nrow(df_raw),
  n_persons = n_distinct(df_raw[[id_var]]),
  min_wave = min(df_raw[[wave_var]], na.rm = TRUE),
  max_wave = max(df_raw[[wave_var]], na.rm = TRUE),
  min_year = min(df_raw[[year_var]], na.rm = TRUE),
  max_year = max(df_raw[[year_var]], na.rm = TRUE),
  cognition_variable = cog_var,
  age_variable = age_var,
  menopause_age_variable = meno_age_var
)

write_csv(raw_inventory, file.path(out_dir, "phase1_raw_inventory.csv"))

## ---- PHASE 1A. Temporal ordering check --------------------------------------

# The exposure must temporally precede the first available cognitive score.
# We compare resolved age at menopause against age at the first non-missing
# cognitive assessment. A respondent is excluded if menopause age is strictly
# greater than first cognitive-assessment age.

first_cognitive_assessment <- df_raw %>%
  filter(!is.na(.data[[cog_var]]), !is.na(.data[[age_var]])) %>%
  arrange(.data[[id_var]], .data[[age_var]], .data[[wave_var]]) %>%
  group_by(.data[[id_var]]) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    !!id_var := .data[[id_var]],
    first_cog_wave = .data[[wave_var]],
    first_cog_year = .data[[year_var]],
    first_cog_age = .data[[age_var]],
    first_cog_score = .data[[cog_var]]
  )

menopause_age_person <- df_raw %>%
  arrange(.data[[id_var]], .data[[wave_var]]) %>%
  group_by(.data[[id_var]]) %>%
  summarise(
    menopause_age_resolved = mode_nonmissing(.data[[meno_age_var]]),
    menopause_age_min = safe_min(.data[[meno_age_var]]),
    menopause_age_max = safe_max(.data[[meno_age_var]]),
    menopause_age_n_reports = sum(!is.na(.data[[meno_age_var]])),
    menopause_age_n_distinct = n_distinct(.data[[meno_age_var]], na.rm = TRUE),
    .groups = "drop"
  )

temporal_order_diagnostics <- first_cognitive_assessment %>%
  left_join(menopause_age_person, by = id_var) %>%
  mutate(
    temporal_order_violation = !is.na(menopause_age_resolved) &
      !is.na(first_cog_age) &
      menopause_age_resolved > first_cog_age,
    temporal_order_ambiguous = !is.na(menopause_age_min) &
      !is.na(menopause_age_max) &
      !is.na(first_cog_age) &
      menopause_age_min <= first_cog_age &
      menopause_age_max > first_cog_age
  )

write_csv(
  temporal_order_diagnostics,
  file.path(out_dir, "phase1_temporal_order_person_diagnostics.csv")
)

temporal_order_summary <- temporal_order_diagnostics %>%
  summarise(
    persons_with_valid_cognition = n(),
    persons_with_reported_menopause_age = sum(!is.na(menopause_age_resolved)),
    persons_temporal_order_violation = sum(temporal_order_violation, na.rm = TRUE),
    persons_temporal_order_ambiguous = sum(temporal_order_ambiguous, na.rm = TRUE)
  )

write_csv(temporal_order_summary, file.path(out_dir, "phase1_temporal_order_summary.csv"))


time_invariant_vars <- c(
  "rabyear", "raracem", "rahispan", "raeduc", "raedyrs",
  "Race", "Edu2", "Edu3", female_var, gender_var
) %>%
  intersect(names(df_raw))

time_invariant_long <- df_raw %>%
  select(all_of(c(id_var, wave_var, time_invariant_vars))) %>%
  arrange(.data[[id_var]], .data[[wave_var]]) %>%
  pivot_longer(
    cols = all_of(time_invariant_vars),
    names_to = "variable",
    values_to = "value",
    values_transform = list(value = as.character)
  )

time_invariant_person_diagnostics <- time_invariant_long %>%
  group_by(.data[[id_var]], variable) %>%
  summarise(
    n_nonmissing = sum(!is.na(value)),
    n_distinct_nonmissing = n_distinct(value, na.rm = TRUE),
    observed_values = paste(unique(value[!is.na(value)]), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(inconsistent = n_distinct_nonmissing > 1)

time_invariant_summary <- time_invariant_person_diagnostics %>%
  group_by(variable) %>%
  summarise(
    persons_checked = n(),
    persons_all_missing = sum(n_nonmissing == 0),
    persons_inconsistent = sum(inconsistent),
    pct_inconsistent = persons_inconsistent / persons_checked,
    max_distinct_values = max(n_distinct_nonmissing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(persons_inconsistent), variable)

time_invariant_inconsistent_cases <- time_invariant_person_diagnostics %>%
  filter(inconsistent) %>%
  arrange(variable, .data[[id_var]])

write_csv(
  time_invariant_summary,
  file.path(out_dir, "phase1_time_invariant_summary.csv")
)

write_csv(
  time_invariant_inconsistent_cases,
  file.path(out_dir, "phase1_time_invariant_inconsistent_cases.csv")
)

time_invariant_resolved <- df_raw %>%
  arrange(.data[[id_var]], .data[[wave_var]]) %>%
  group_by(.data[[id_var]]) %>%
  summarise(
    across(all_of(time_invariant_vars), mode_nonmissing, .names = "{.col}__resolved"),
    .groups = "drop"
  )

df_phase1_resolved <- df_raw %>%
  left_join(time_invariant_resolved, by = id_var)

for (v in time_invariant_vars) {
  resolved_v <- paste0(v, "__resolved")
  df_phase1_resolved[[v]] <- df_phase1_resolved[[resolved_v]]
}

df_phase1_resolved <- df_phase1_resolved %>%
  select(-ends_with("__resolved"))

## ---- PHASE 1C. Menopause transition logic -----------------------------------

# A respondent should not move backward from postmenopausal to premenopausal
# status across waves. Unknown or missing statuses are not treated as backward
# transitions; they are handled later as invalid exposure information.

menopause_transition_wave <- df_phase1_resolved %>%
  arrange(.data[[id_var]], .data[[wave_var]]) %>%
  mutate(
    meno_status_standard = case_when(
      str_detect(str_to_lower(.data[[meno_status_var]]), "pre") ~ "premenopausal",
      str_detect(str_to_lower(.data[[meno_status_var]]), "post") ~ "postmenopausal",
      TRUE ~ NA_character_
    ),
    meno_status_rank = case_when(
      meno_status_standard == "premenopausal" ~ 0L,
      meno_status_standard == "postmenopausal" ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  group_by(.data[[id_var]]) %>%
  mutate(
    post_observed = if_else(is.na(meno_status_rank), 0L, as.integer(meno_status_rank == 1L)),
    prior_post_observed = lag(cummax(post_observed), default = 0L),
    backward_meno_transition = meno_status_rank == 0L & prior_post_observed == 1L
  ) %>%
  ungroup()

menopause_transition_person_diagnostics <- menopause_transition_wave %>%
  group_by(.data[[id_var]]) %>%
  summarise(
    any_backward_meno_transition = any(backward_meno_transition, na.rm = TRUE),
    first_backward_wave = safe_min(if_else(backward_meno_transition, .data[[wave_var]], NA_real_)),
    status_sequence = paste(
      paste0(.data[[wave_var]], ":", coalesce(meno_status_standard, "missing_or_unknown")),
      collapse = " -> "
    ),
    .groups = "drop"
  )

write_csv(
  menopause_transition_person_diagnostics,
  file.path(out_dir, "phase1_menopause_transition_person_diagnostics.csv")
)

menopause_transition_summary <- menopause_transition_person_diagnostics %>%
  summarise(
    persons_checked = n(),
    persons_backward_meno_transition = sum(any_backward_meno_transition, na.rm = TRUE)
  )

write_csv(
  menopause_transition_summary,
  file.path(out_dir, "phase1_menopause_transition_summary.csv")
)

## ---- PHASE 1D. Cognitive missingness diagnostics ----------------------------

# Build a person-wave grid from each respondent's first observed wave through
# the last wave in the extract. This distinguishes an observed row with missing
# cognition from a missing person-wave row, which is attrition-compatible.

all_waves <- sort(unique(df_phase1_resolved[[wave_var]]))
max_wave <- max(all_waves, na.rm = TRUE)

person_wave_span <- df_phase1_resolved %>%
  group_by(.data[[id_var]]) %>%
  summarise(first_observed_wave = min(.data[[wave_var]], na.rm = TRUE), .groups = "drop")

person_wave_grid <- person_wave_span %>%
  tidyr::crossing(!!sym(wave_var) := all_waves) %>%
  filter(.data[[wave_var]] >= first_observed_wave) %>%
  select(all_of(c(id_var, wave_var)))

cognitive_by_wave <- df_phase1_resolved %>%
  group_by(.data[[id_var]], .data[[wave_var]]) %>%
  summarise(
    row_present = TRUE,
    year = first_nonmissing(.data[[year_var]]),
    age = first_nonmissing(.data[[age_var]]),
    cognition = first_nonmissing(.data[[cog_var]]),
    .groups = "drop"
  )

cognitive_missingness_panel <- person_wave_grid %>%
  left_join(cognitive_by_wave, by = c(id_var, wave_var)) %>%
  mutate(
    row_present = coalesce(row_present, FALSE),
    valid_cognition = row_present & !is.na(cognition),
    missing_state = case_when(
      valid_cognition ~ "valid_cognition",
      row_present & is.na(cognition) ~ "observed_row_missing_cognition",
      !row_present ~ "no_person_wave_row",
      TRUE ~ "unclassified"
    )
  )

cognitive_missingness_by_wave <- cognitive_missingness_panel %>%
  group_by(.data[[wave_var]]) %>%
  summarise(
    persons_in_risk_set = n(),
    person_wave_rows_present = sum(row_present),
    valid_cognition = sum(valid_cognition),
    observed_rows_missing_cognition = sum(missing_state == "observed_row_missing_cognition"),
    no_person_wave_row = sum(missing_state == "no_person_wave_row"),
    total_missing_cognition_or_row = observed_rows_missing_cognition + no_person_wave_row,
    pct_missing_cognition_or_row = total_missing_cognition_or_row / persons_in_risk_set,
    .groups = "drop"
  )

cognitive_missingness_person <- cognitive_missingness_panel %>%
  arrange(.data[[id_var]], .data[[wave_var]]) %>%
  group_by(.data[[id_var]]) %>%
  summarise(
    n_panel_waves = n(),
    n_valid_cognition = sum(valid_cognition),
    n_observed_rows_missing_cognition = sum(missing_state == "observed_row_missing_cognition"),
    n_no_person_wave_row = sum(missing_state == "no_person_wave_row"),
    first_valid_cog_wave = safe_min(if_else(valid_cognition, .data[[wave_var]], NA_real_)),
    last_valid_cog_wave = safe_max(if_else(valid_cognition, .data[[wave_var]], NA_real_)),
    any_missing_between_valid_cognition = any(
      !valid_cognition &
        .data[[wave_var]] > first_valid_cog_wave &
        .data[[wave_var]] < last_valid_cog_wave,
      na.rm = TRUE
    ),
    any_missing_after_last_valid_cognition = any(
      !valid_cognition & .data[[wave_var]] > last_valid_cog_wave,
      na.rm = TRUE
    ),
    any_missing_before_first_valid_cognition = any(
      !valid_cognition & .data[[wave_var]] < first_valid_cog_wave,
      na.rm = TRUE
    ),
    missingness_pattern = paste(
      case_when(
        valid_cognition ~ "O",
        missing_state == "observed_row_missing_cognition" ~ "M",
        missing_state == "no_person_wave_row" ~ "N",
        TRUE ~ "?"
      ),
      collapse = ""
    ),
    .groups = "drop"
  ) %>%
  mutate(
    missingness_class = case_when(
      n_valid_cognition == 0 ~ "zero_valid_cognitive_assessments",
      n_observed_rows_missing_cognition == 0 & n_no_person_wave_row == 0 ~ "complete_from_entry_to_last_wave",
      any_missing_between_valid_cognition ~ "intermittent_missingness",
      any_missing_after_last_valid_cognition ~ "monotone_missingness_attrition_compatible",
      any_missing_before_first_valid_cognition ~ "missing_before_first_valid_cognition_only",
      TRUE ~ "other_missingness_pattern"
    )
  )

cognitive_missingness_summary <- cognitive_missingness_person %>%
  count(missingness_class, sort = TRUE, name = "n_persons") %>%
  mutate(pct_persons = n_persons / sum(n_persons))

write_csv(
  cognitive_missingness_by_wave,
  file.path(out_dir, "phase1_cognitive_missingness_by_wave.csv")
)

write_csv(
  cognitive_missingness_person,
  file.path(out_dir, "phase1_cognitive_missingness_person_patterns.csv")
)

write_csv(
  cognitive_missingness_summary,
  file.path(out_dir, "phase1_cognitive_missingness_summary.csv")
)

## ---- PHASE 1E. Apply integrity exclusions -----------------------------------

# Strict Phase 1 exclusions:
# 1. Menopause age occurs after first available cognitive assessment.
# 2. Menopause status transitions backward from postmenopausal to premenopausal.
#
# The derived C253_major_inconsistency flag is exported for review but is not
# used as an automatic exclusion unless explicitly added below. This preserves
# the requested Phase 2 exposure-missing cascade.

major_menopause_age_inconsistency <- df_phase1_resolved %>%
  group_by(.data[[id_var]]) %>%
  summarise(
    any_c253_major_inconsistency = if (
      "C253_major_inconsistency" %in% names(df_phase1_resolved)
    ) {
      any(.data[["C253_major_inconsistency"]] == 1, na.rm = TRUE)
    } else {
      FALSE
    },
    .groups = "drop"
  )

integrity_person_flags <- df_phase1_resolved %>%
  distinct(.data[[id_var]]) %>%
  left_join(
    temporal_order_diagnostics %>%
      select(
        all_of(id_var),
        temporal_order_violation,
        temporal_order_ambiguous,
        first_cog_age,
        menopause_age_resolved
      ),
    by = id_var
  ) %>%
  left_join(
    menopause_transition_person_diagnostics %>%
      select(all_of(id_var), any_backward_meno_transition, first_backward_wave),
    by = id_var
  ) %>%
  left_join(major_menopause_age_inconsistency, by = id_var) %>%
  mutate(
    across(
      c(
        temporal_order_violation,
        temporal_order_ambiguous,
        any_backward_meno_transition,
        any_c253_major_inconsistency
      ),
      ~ replace_na(.x, FALSE)
    ),
    phase1_integrity_exclude = temporal_order_violation | any_backward_meno_transition,
    phase1_integrity_exclusion_reason = case_when(
      temporal_order_violation & any_backward_meno_transition ~
        "menopause_after_first_cognition_and_backward_transition",
      temporal_order_violation ~ "menopause_after_first_cognition",
      any_backward_meno_transition ~ "post_to_pre_menopause_backward_transition",
      TRUE ~ NA_character_
    )
  )

write_csv(integrity_person_flags, file.path(out_dir, "phase1_integrity_person_flags.csv"))

phase1_integrity_cascade <- tibble(
  step = c(
    "Raw person file",
    "Drop menopause-after-first-cognition temporal violations",
    "Drop post-to-pre menopause backward transitions"
  ),
  n_dropped = c(
    0L,
    sum(integrity_person_flags$temporal_order_violation, na.rm = TRUE),
    sum(
      integrity_person_flags$any_backward_meno_transition &
        !integrity_person_flags$temporal_order_violation,
      na.rm = TRUE
    )
  )
) %>%
  mutate(
    n_remaining = n_distinct(df_phase1_resolved[[id_var]]) - cumsum(n_dropped)
  )

write_csv(
  phase1_integrity_cascade,
  file.path(out_dir, "phase1_integrity_exclusion_cascade.csv")
)

df_phase1_clean <- df_phase1_resolved %>%
  left_join(
    integrity_person_flags %>%
      select(all_of(id_var), phase1_integrity_exclude, phase1_integrity_exclusion_reason),
    by = id_var
  ) %>%
  filter(!phase1_integrity_exclude) %>%
  select(-phase1_integrity_exclude)

write_rds(df_phase1_clean, file.path(out_dir, "hrs_phase1_integrity_clean_long.rds"))
write_csv(df_phase1_clean, file.path(out_dir, "hrs_phase1_integrity_clean_long.csv"))

## ---- PHASE 2A. Define eligibility fields ------------------------------------

# Valid exposure is an established menopause type/timing group. Premenopausal,
# unknown type, unknown/inconsistent, and missing menopause groups are not valid
# time-zero exposures for the main analytic cohort.

df_phase2_base <- df_phase1_clean %>%
  mutate(
    is_female = flag_female(.data[[female_var]], .data[[gender_var]]),
    valid_cognition = !is.na(.data[[cog_var]]),
    age_eligible = !is.na(.data[[age_var]]) & .data[[age_var]] >= 50,
    valid_menopause_exposure = .data[[meno_combo_var]] %in% valid_meno_groups,
    natural_menopause_group = .data[[meno_combo_var]] %in% natural_meno_groups,
    eligible_timezero_wave =
      is_female &
      age_eligible &
      valid_cognition &
      valid_menopause_exposure
  )

phase2_person_flags <- df_phase2_base %>%
  arrange(.data[[id_var]], .data[[wave_var]]) %>%
  group_by(.data[[id_var]]) %>%
  summarise(
    person_female = mode_nonmissing(is_female),
    any_valid_menopause_exposure = any(valid_menopause_exposure, na.rm = TRUE),
    any_valid_cognitive_assessment = any(valid_cognition, na.rm = TRUE),
    n_valid_cognitive_assessments = sum(valid_cognition, na.rm = TRUE),
    any_eligible_timezero_wave = any(eligible_timezero_wave, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(person_female = replace_na(person_female, FALSE))

exposure_availability_diagnostics <- df_phase2_base %>%
  group_by(.data[[id_var]]) %>%
  summarise(
    any_valid_menopause_exposure = any(valid_menopause_exposure, na.rm = TRUE),
    ever_premenopausal = any(.data[[meno_combo_var]] == "7.Premenopausal", na.rm = TRUE),
    ever_unknown_type = any(.data[[meno_combo_var]] == "8.Unknown type", na.rm = TRUE),
    ever_unknown_inconsistent = any(.data[[meno_combo_var]] == "9.Unknown/Inconsistent", na.rm = TRUE),
    all_meno_combo_missing = all(is.na(.data[[meno_combo_var]])),
    observed_meno_combo_values = paste(sort(unique(na.omit(.data[[meno_combo_var]]))), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    exposure_problem_class = case_when(
      any_valid_menopause_exposure ~ "has_valid_exposure",
      all_meno_combo_missing ~ "all_menopause_combo_missing",
      ever_premenopausal & !ever_unknown_type & !ever_unknown_inconsistent ~ "only_premenopausal_or_missing",
      ever_unknown_type | ever_unknown_inconsistent ~ "unknown_or_inconsistent_only",
      TRUE ~ "other_no_valid_exposure"
    )
  )

write_csv(
  exposure_availability_diagnostics,
  file.path(out_dir, "phase2_exposure_availability_person_diagnostics.csv")
)

write_csv(
  exposure_availability_diagnostics %>%
    count(exposure_problem_class, sort = TRUE, name = "n_persons"),
  file.path(out_dir, "phase2_exposure_availability_summary.csv")
)

## ---- PHASE 2B. Select Time Zero ---------------------------------------------

# Baseline/Time Zero is the first wave at which a respondent simultaneously
# satisfies: female, age >= 50, valid cognitive score, and established menopause
# type/timing. Rows before this wave are not part of the analytic risk period.

baseline_all_candidates <- df_phase2_base %>%
  filter(eligible_timezero_wave) %>%
  arrange(.data[[id_var]], .data[[wave_var]]) %>%
  group_by(.data[[id_var]]) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    !!id_var := .data[[id_var]],
    baseline_wave = .data[[wave_var]],
    baseline_year = .data[[year_var]],
    baseline_age = .data[[age_var]],
    baseline_cognition = .data[[cog_var]],
    baseline_cogfunction = .data[[cog_function_var]],
    baseline_dementia = .data[[cog_function_var]] == 3,
    baseline_normal = case_when(
      .data[[cog_function_var]] == 1 ~ 1L,
      .data[[cog_function_var]] == 2 ~ 0L,
      .data[[cog_function_var]] == 3 ~ 0L,
      TRUE ~ NA_integer_
    ),
    baseline_cind = .data[[cog_function_var]] == 2,
    menopause_group = .data[[meno_combo_var]],
    menopause_type = .data[[meno_type_var]],
    menopause_age_category = .data[[meno_age5_var]],
    menopause_age_resolved = .data[[meno_age_var]],
    natural_menopause_group = .data[[meno_combo_var]] %in% natural_meno_groups
  )

## ---- PHASE 2C. Strict exclusion cascade -------------------------------------

ids0 <- phase2_person_flags[[id_var]]
ids1 <- phase2_person_flags %>%
  filter(person_female) %>%
  pull(all_of(id_var))

ids2 <- phase2_person_flags %>%
  filter(.data[[id_var]] %in% ids1, any_valid_menopause_exposure) %>%
  pull(all_of(id_var))

ids3 <- phase2_person_flags %>%
  filter(.data[[id_var]] %in% ids2, any_valid_cognitive_assessment) %>%
  pull(all_of(id_var))

ids4 <- phase2_person_flags %>%
  filter(.data[[id_var]] %in% ids3, any_eligible_timezero_wave) %>%
  pull(all_of(id_var))

dementia_baseline_ids <- baseline_all_candidates %>%
  filter(replace_na(baseline_dementia, FALSE)) %>%
  pull(all_of(id_var))

ids5 <- setdiff(ids4, dementia_baseline_ids)

flow_row <- function(step, before_ids, after_ids, long_data) {
  tibble(
    step = step,
    n_persons_before = length(unique(before_ids)),
    n_persons_dropped = length(setdiff(unique(before_ids), unique(after_ids))),
    n_persons_remaining = length(unique(after_ids)),
    n_rows_remaining = long_data %>%
      filter(.data[[id_var]] %in% after_ids) %>%
      nrow()
  )
}

phase2_exclusion_cascade <- bind_rows(
  tibble(
    step = "Start after Phase 1 integrity exclusions",
    n_persons_before = NA_integer_,
    n_persons_dropped = 0L,
    n_persons_remaining = length(unique(ids0)),
    n_rows_remaining = nrow(df_phase2_base)
  ),
  flow_row("Drop men", ids0, ids1, df_phase2_base),
  flow_row("Drop women with missing/invalid menopause type-timing exposure", ids1, ids2, df_phase2_base),
  flow_row("Drop women with zero valid cognitive assessments", ids2, ids3, df_phase2_base),
  flow_row("Drop women with no eligible Time Zero wave", ids3, ids4, df_phase2_base),
  flow_row("Drop prevalent dementia at baseline", ids4, ids5, df_phase2_base)
)

write_csv(
  phase2_exclusion_cascade,
  file.path(out_dir, "phase2_strict_exclusion_cascade.csv")
)

baseline_person <- baseline_all_candidates %>%
  filter(.data[[id_var]] %in% ids5) %>%
  arrange(.data[[id_var]])

write_csv(baseline_person, file.path(out_dir, "hrs_phase2_baseline_person.csv"))
write_rds(baseline_person, file.path(out_dir, "hrs_phase2_baseline_person.rds"))

## ---- PHASE 2D. Build analytic longitudinal datasets -------------------------

analytic_long_allwaves <- df_phase2_base %>%
  inner_join(
    baseline_person %>%
      select(
        all_of(id_var),
        baseline_wave,
        baseline_year,
        baseline_age,
        baseline_cognition,
        baseline_cogfunction,
        baseline_normal,
        baseline_cind,
        menopause_group,
        menopause_type,
        menopause_age_category,
        menopause_age_resolved,
        natural_menopause_group
      ),
    by = id_var
  ) %>%
  filter(.data[[wave_var]] >= baseline_wave) %>%
  mutate(
    time_since_baseline_wave = .data[[wave_var]] - baseline_wave,
    time_since_baseline_year = .data[[year_var]] - baseline_year,
    age_centered_65 = .data[[age_var]] - 65,
    post_baseline_valid_cognition = !is.na(.data[[cog_var]])
  ) %>%
  arrange(.data[[id_var]], .data[[wave_var]])

analytic_long_modelready <- analytic_long_allwaves %>%
  filter(post_baseline_valid_cognition)

analytic_cohort_summary <- tibble(
  n_persons = n_distinct(analytic_long_allwaves[[id_var]]),
  n_person_waves_all_post_timezero = nrow(analytic_long_allwaves),
  n_person_waves_valid_cognition = nrow(analytic_long_modelready),
  mean_valid_cognitive_assessments_per_person =
    analytic_long_modelready %>%
    count(.data[[id_var]]) %>%
    summarise(mean_n = mean(n), .groups = "drop") %>%
    pull(mean_n),
  median_valid_cognitive_assessments_per_person =
    analytic_long_modelready %>%
    count(.data[[id_var]]) %>%
    summarise(median_n = median(n), .groups = "drop") %>%
    pull(median_n),
  n_baseline_normal = sum(baseline_person$baseline_normal == 1, na.rm = TRUE),
  n_baseline_cind = sum(baseline_person$baseline_cind, na.rm = TRUE),
  n_baseline_normal_missing = sum(is.na(baseline_person$baseline_normal))
)

write_csv(
  analytic_cohort_summary,
  file.path(out_dir, "phase2_analytic_cohort_summary.csv")
)

write_csv(
  baseline_person %>%
    count(menopause_group, sort = TRUE, name = "n_persons"),
  file.path(out_dir, "phase2_baseline_menopause_group_counts.csv")
)

write_rds(
  analytic_long_allwaves,
  file.path(out_dir, "hrs_phase2_analytic_long_allwaves.rds")
)

write_csv(
  analytic_long_allwaves,
  file.path(out_dir, "hrs_phase2_analytic_long_allwaves.csv")
)

write_rds(
  analytic_long_modelready,
  file.path(out_dir, "hrs_phase2_analytic_long_modelready.rds")
)

write_csv(
  analytic_long_modelready,
  file.path(out_dir, "hrs_phase2_analytic_long_modelready.csv")
)

## ---- PHASE 2E. Console review packet ----------------------------------------

cat("\n================ PHASE 1-2 REVIEW PACKET ================\n")
cat("\nRaw inventory:\n")
print(raw_inventory)

cat("\nPhase 1 integrity exclusion cascade:\n")
print(phase1_integrity_cascade)

cat("\nTemporal ordering summary:\n")
print(temporal_order_summary)

cat("\nTime-invariant variable summary:\n")
print(time_invariant_summary)

cat("\nMenopause transition summary:\n")
print(menopause_transition_summary)

cat("\nCognitive missingness summary:\n")
print(cognitive_missingness_summary)

cat("\nPhase 2 strict exclusion cascade:\n")
print(phase2_exclusion_cascade)

cat("\nAnalytic cohort summary:\n")
print(analytic_cohort_summary)

cat("\nBaseline menopause group counts:\n")
print(
  baseline_person %>%
    count(menopause_group, sort = TRUE, name = "n_persons")
)

cat("\nOutputs written to:\n", out_dir, "\n")
cat("Please review the printed tables and output CSV diagnostics before Phase 3.\n")

