#' Run all UPCR steps (without DQ investigation)
#'
#' UNDER DEVELOPMENT
#'
#' Please note that this function currently joins on date. A 24 hour window
#' should be considered for future iterations.
#'
#' This function generates a table of selected UPCRs, based on the following
#' approach:
#'
#' Exclude/include known urine protein mismappings
#' Exclude/include known urine creatinine mismappings
#' For a quantitative urine protein measurement on a given date for a patient:
#'    If UPCR directly reported, use
#'    If no UPCR directly reported, calculate if following criteria are met:
#'       Urine creatinine measurement for patient on same date
#'       Numeric value is populated
#'       Exclude measurements for which the source value indicates a range is
#'       present
#' Exception for “<” or “less than”: Study investigator informed us that lower
#' bounds are sometimes reported and these values should be used in UPCR
#' calculation
#'
#' The returned output should be further and analysed with input from the study
#' team:
#' Check for urine protein/urine creatinine mismappings for cohort(s) of
#' interest and add to known mismappings (scope-dependent)
#' If more than one urine creatinine available for a single urine protein,
#' calculate multiple and select via criteria agreed upon by study team,
#' potentially including the following:
#'
#' UPCR within sensible bounds:
#'    mg:mg/g:g: >= 0.01, <= 200
#'    mg:g: >= 10, <= 20,000
#'       Most should fall within 0.1 - 20 range (>2 is nephrotic range)
#' Selection hierarchy:
#'    Closest timestamp (if not 000000)
#'    Matching collection timing information (e.g. 12 hr or 24hr)
#'
#' @param proteinuria_urinalysis_codeset Codeset of of point-of-care urine
#' protein testing (qualitative / semi-quantitative)
#' @param urine_protein_codeset Codeset of quantitative urine protein
#' measurements
#' @param measurement_tbl Measurement table for study
#' @param urine_creatinine_codeset Codeset of quantitative urine creatinine
#' measurements
#' @param urine_prot_creat_ratio Codeset of directly reported UPCRs
#' @param prot_mismap_file csv with known mismappings associated with urine
#' protein measurements
#' @param creat_mismap_file csv with known mismappings associated with urine
#' creatinine measurements
#'
#' @return Table of selected UPCRs
#'
run_all_upcr_steps <- function(cohort,
                               proteinuria_urinalysis_codeset = load_codeset("proteinuria_urinalysis"),
                               urine_protein_codeset = load_codeset("urine_protein"),
                               measurement_tbl = cdm_tbl("measurement"),
                               urine_creatinine_codeset,
                               urine_prot_creat_ratio,
                               prot_mismap_file,
                               creat_mismap_file) {
  # qual/semi-quant codeset
  proteinuria_urinalysis_vctr <-
    proteinuria_urinalysis_codeset %>% select(concept_id) %>% pull()
  
  # quant codeset
  urine_protein_vctr <-
    urine_protein_codeset %>% select(concept_id) %>% pull()
  
  # urine creatinine codeset
  urine_creatinine_vctr <-
    urine_creatinine_codeset %>% select(concept_id) %>% pull()
  
  # upcr codeset
  urine_prot_creat_ratio_vctr <-
    urine_prot_creat_ratio %>% select(concept_id) %>% pull()
  
  urine_meas_vctr <-
    c(
      proteinuria_urinalysis_vctr,
      urine_protein_vctr,
      urine_creatinine_vctr,
      urine_prot_creat_ratio_vctr
    ) %>% unique()
  
  # urine protein measurements for cohort
  urine_meas_tbl <- measurement_tbl %>%
    filter(measurement_concept_id %in% urine_meas_vctr) %>%
    inner_join(select(cohort,
                      person_id), by = "person_id") %>%
    compute_new()
  
  # URINE PROTEIN ----
  
  # combine qual/semi-quant and quant codesets
  comb_urine_protein_vctr <-
    c(proteinuria_urinalysis_vctr, urine_protein_vctr) %>% unique()
  
  # collect all measurements for combined codesets
  # with flag for qual/semi_quant or quant codesets
  comb_urine_protein <- urine_meas_tbl %>%
    filter(measurement_concept_id %in% comb_urine_protein_vctr) %>%
    mutate(
      codeset = if_else(
        measurement_concept_id %in% urine_protein_vctr,
        "quant",
        "qual_semiquant"
      )
    ) %>%
    select(
      measurement_concept_id,
      measurement_date,
      measurement_datetime,
      measurement_id,
      measurement_source_concept_id,
      measurement_source_value,
      person_id,
      specimen_concept_id,
      specimen_source_value,
      unit_concept_id,
      unit_source_value,
      value_as_concept_id,
      value_as_number,
      value_source_value,
      visit_occurrence_id,
      measurement_concept_name,
      measurement_source_concept_name,
      specimen_concept_name,
      unit_concept_name,
      value_as_concept_name,
      site,
      codeset
    ) %>%
    collect_new()
  
  # URINE CREATININE ----
  
  # collect all measurements for urine creatinine codeset
  urine_creatinine <- urine_meas_tbl %>%
    filter(measurement_concept_id %in% urine_creatinine_vctr) %>%
    select(
      measurement_concept_id,
      measurement_date,
      measurement_datetime,
      measurement_id,
      measurement_source_concept_id,
      measurement_source_value,
      person_id,
      specimen_concept_id,
      specimen_source_value,
      unit_concept_id,
      unit_source_value,
      value_as_concept_id,
      value_as_number,
      value_source_value,
      visit_occurrence_id,
      measurement_concept_name,
      measurement_source_concept_name,
      specimen_concept_name,
      unit_concept_name,
      value_as_concept_name,
      site
    ) %>% # restrict to required fields
    collect_new()
  
  # COMBINE URINE PROTEIN AND URINE CREATININE ON SAME DATE ----
  
  # Prepare urine protein measurements
  # -- Restrict to required fields
  # -- Add value_as_number related variables: van_is_na
  # (whether value_as_number is NA) and parsed_van (extracting
  # value_as_number fields where "<" or "less than" in the source value and
  # value_as_number not extracted)
  # -- Convert to mg/dL as applicable, according to information
  # provided by unit_concept_id
  #  -- Add flag for whether useful unit information available
  # -- Add flag for whether source value contains range information
  # (with exclusion of "<" or "less than", which is expected)
  urine_protein_concs_prep <- comb_urine_protein %>%
    prepare_urine_conc()
  
  # Prepare urine creatinine measurements
  # -- Restrict to required fields
  # -- Add value_as_number related variables: van_is_na
  # (whether value_as_number is NA) and parsed_van (extracting
  # value_as_number fields where "<" or "less than" in the source value and
  # value_as_number not extracted)
  # -- Convert to mg/dL as applicable, according to information
  # provided by unit_concept_id
  #  -- Add flag for whether useful unit information available
  # -- Add flag for whether source value contains range information
  # (with exclusion of "<" or "less than", which is expected)
  urine_creat_concs_prep <- urine_creatinine %>%
    prepare_urine_conc()
  
  # Left join urine creatinine to urine protein on person_id and
  # measurement_date and add flags for whether there's a creatinine on the same
  # date and whether upcr criteria are met
  comb_urine_protein_creat <-
    join_prot_creat_concs(urine_protein_concs_prep,
                          urine_creat_concs_prep)
  
  # ADD KNOWN MISMAP FLAGS TO URINE PROTEIN AND URINE CREATININE ----
  
  # Add known mismap flags for both urine protein and urine creatine
  comb_urine_protein_creat_mismap_flag <-
    comb_urine_protein_creat %>%
    add_known_mismap_flags(prot_mismaps = prot_mismap_file,
                           creat_mismaps = creat_mismap_file)
  
  # CALCULATE UPCR ----
  # calculate upcr with parsed and converted (based on unit info - note unreliablity)
  upcr_calc <- comb_urine_protein_creat_mismap_flag %>%
    get_upcr_calcs(creat_value_name = "parsed_van_convert_creat",
                   prot_value_name = "parsed_van_convert")
  
  # CODESET-BASED UPCR ----
  
  # Get codeset-provided UPCR measurements
  upcr_codeset_meas <- urine_meas_tbl %>%
    filter(measurement_concept_id %in% urine_prot_creat_ratio_vctr) %>%
    select(
      person_id,
      site,
      measurement_id,
      measurement_date,
      measurement_datetime,
      measurement_source_value,
      measurement_concept_id,
      measurement_concept_name,
      specimen_concept_id,
      specimen_concept_name,
      unit_concept_id,
      unit_concept_name,
      unit_source_value,
      value_source_value,
      value_as_number,
      value_as_concept_id,
      value_as_concept_name
    ) %>%
    collect_new()
  
  # Prepare codeset upcr measurements by restricting to required fields,
  # adding van variables and upcr_codeset from parsed value as number
  # and converting to mg:mg based on provided unit_concept_id
  upcr_codeset_meas_prep <- prepare_codeset_upcrs(upcr_codeset_meas)
  
  # COMBINE CODESET PROVIDED UPCR AND CALCULATED UPCR ----
  
  # full join codeset provided and calculated upcrs on
  # person_id and measurement_date
  upcr_calc_or_codeset <-
    combine_codeset_calc_upcrs(upcr_codeset_meas_tbl = upcr_codeset_meas_prep,
                               upcr_calc_meas_tbl = upcr_calc)
  
  # select upcr from codeset and calculated upcrs (favor
  # codeset-provided where available) and add helper fields
  selected_upcr <- select_upcr(upcr_calc_or_codeset)
  
  return(selected_upcr)
  
}

#' Prepare urine concentration (for now, urine protein and urine creatinine)
#' measurements for UPCR calculations
#'
#' -- Restrict to required fields
#' -- Add value_as_number related variables: van_is_na
#' (whether value_as_number is NA) and parsed_van (extracting
#' value_as_number fields where "<" or "less than" in the source value and
#' value_as_number not extracted)
#' -- Convert to mg/dL as applicable, according to information
#' provided by unit_concept_id
#' -- Add flag for whether source value contains range information
#' (with exclusion of "<" or "less than", which is expected)
#'
#' @param urine_conc Collected urine concentration measurements
#'
#' @return Prepared urine concentration measurements
#'
prepare_urine_conc <- function(urine_conc) {
  urine_conc_prep <- urine_conc %>%
    select(
      starts_with("codeset"),
      # use starts_with() because only applicable to urine protein
      measurement_concept_id,
      measurement_concept_name,
      measurement_date,
      measurement_datetime,
      measurement_source_value,
      value_as_number,
      value_as_concept_id,
      value_as_concept_name,
      value_source_value,
      specimen_concept_id,
      specimen_concept_name,
      unit_concept_id,
      unit_concept_name,
      unit_source_value,
      site,
      measurement_id,
      person_id
    ) %>%
    add_van_variables() %>%
    convert_to_mg_per_dl() %>%
    add_src_val_is_range_flag()
  
  return(urine_conc_prep)
}

#' Left join urine creatinine to urine protein measurements
#' on person_id and measurement date, adding flags for whether
#' there's a creatinine measurement available on the same date
#' and whether upcr criteria are met
#'
#' Suffix _creat added for variables from creatinine measurement
#' table
#'
#' @param urine_protein_meas Urine protein measurements
#' @param urine_creatinine_meas Urine creatinine measurements
#'
#' @return Combined measurement table with additional flags
#'
join_prot_creat_concs <- function(urine_protein_meas,
                                  urine_creatinine_meas) {
  prot_creat_concs_join <- urine_protein_meas %>%
    left_join(
      urine_creatinine_meas,
      by = c("person_id", "measurement_date"),
      suffix = c("", "_creat")
    ) %>%
    mutate(creat_same_day = !is.na(measurement_id_creat)) %>%
    add_upcr_criteria_met_flag()
  
  return(prot_creat_concs_join)
}

#' Explore LOINC and CPT codes
#'
#' Extract characters after the | symbol from both protein and creatinine measurement source values
#' to get LOINC or CPT codes. Join these LOINC or CPT codes on the vocabulary. Return a table of
#' vocabulary information for these LOINC and CPT codes for manual review for mismappings.
#'
#' @param meas_table measurement table which must contain measurement_source_value (msv for protein)
#' and measurement_source_value_creat (msv for creatinine)
#'
#' @return Vocabulary information for LOINC and CPT extracted from measurement source values for
#' manual review
#'
explore_loinc_cpt <- function(meas_table) {
  meas_table_loinc_cpt_extract <- meas_table %>%
    mutate(
      prot_cpt_or_loinc = str_extract(measurement_source_value, "[^|]*$"),
      creat_cpt_or_loinc = str_extract(measurement_source_value_creat, "[^|]*$")
    ) %>%
    distinct(
      prot_cpt_or_loinc,
      measurement_concept_id,
      measurement_concept_name,
      creat_cpt_or_loinc,
      measurement_concept_id_creat,
      measurement_concept_name_creat
    )
  prot_cpt_loinc_check <- meas_table_loinc_cpt_extract %>%
    distinct(prot_cpt_or_loinc,
             measurement_concept_id,
             measurement_concept_name) %>%
    copy_to_new(
      dest = config('db_src'),
      df = .,
      name = "prot_cpt_loinc_check",
      temporary = TRUE
    ) %>%
    rename(cpt_or_loinc = prot_cpt_or_loinc) %>%
    mutate(component = "prot")
  
  creat_cpt_loinc_check <- meas_table_loinc_cpt_extract %>%
    distinct(
      creat_cpt_or_loinc,
      measurement_concept_id_creat,
      measurement_concept_name_creat
    ) %>%
    copy_to_new(
      dest = config('db_src'),
      df = .,
      name = "creat_cpt_loinc_check",
      temporary = TRUE
    ) %>%
    rename(cpt_or_loinc = creat_cpt_or_loinc) %>%
    mutate(component = "creat")
  
  loinc_cpt_summary <- prot_cpt_loinc_check %>%
    dplyr::union(creat_cpt_loinc_check) %>%
    inner_join(
      vocabulary_tbl("concept") %>%
        filter(vocabulary_id == "CPT4" |
                 vocabulary_id == "LOINC"),
      by = c("cpt_or_loinc" = "concept_code")
    ) %>%
    collect_new() %>%
    arrange(component, cpt_or_loinc)
  
  return(loinc_cpt_summary)
}

#' Add known mismap flags
#'
#' Add both urine protein and urine creatine mismaps flags via
#' add_prot_mismap_flags and add_creat_mismap_flags and add additional flags for
#' any urine protein mismapping (any_prot_mismap) and any urine creatinine
#' mismapping (any_creat_mismap) as well as any_mismap for either
#'
#' @param meas_tbl Measurement table (must contain measurement_source_value
#' field for urine protein and measurement_source_value_creat for urine
#' creatine)
#' @param prot_mismaps csv file with known potential LOINC/CPT urine protein
#' mismaps
#' @param creat_mismaps csv file with known potential LOINC/CPT urine protein
#' mismaps
#'
#' @return Measurement table with mismap flags
#'
add_known_mismap_flags <- function(meas_tbl,
                                   prot_mismaps,
                                   creat_mismaps) {
  meas_tbl_mismap_flags <- meas_tbl %>%
    add_prot_mismap_flags(mismaps = prot_mismaps) %>%
    add_creat_mismap_flags(mismaps = creat_mismaps) %>%
    mutate(
      any_creat_mismap = if_else(is.na(any_creat_mismap), FALSE, any_creat_mismap),
      any_prot_mismap = if_else(is.na(any_prot_mismap), FALSE, any_prot_mismap),
      any_mismap = any_creat_mismap | any_prot_mismap
    )
  
  return(meas_tbl_mismap_flags)
}

#' Calculate UPCR based on urine creatinine and urine protein measurements
#' Only include measurements where the following hold:
#' Urine protein from quantative codeset OR urine protein from qualitative/semi-quantitative codeset that we determined *should* be mapped to quantitative codeset
#' UPCR criteria are met
#' Numeric value populated
#' Urine creatinine is available on same date for patient
#' No range information in source value, with exception of < and less than
#' No mismappings, with exception of mismapping to a qualitative/semi-quantitative code
#'
#' @param comb_urine_protein_creat Table with urine protein and urine creatinine
#' measurements on same date for a patient
#' @param creat_value_name Name urine creatinine variable which should be used
#' for calculation
#' @param prot_value_name Name of urine protein variable which should be used
#' for calculation
#'
#' @return Table included calculated UPCR (upcr_calc) and other helpful fields
#'
get_upcr_calcs <- function(comb_urine_protein_creat,
                           creat_value_name,
                           prot_value_name) {
  upcr_calcs <- comb_urine_protein_creat %>%
    filter(
      (
        # use urine protein measurements if from the quantitative codeset AND
        # upcr criteria are met AND no mismaps identified for urine protein or
        # urine creatinine
        codeset == "quant" &
          upcr_criteria_met == TRUE &
          any_mismap == FALSE
      )
      |
        (
          # or if upcr criteria are met and the measurement is mismapped to a
          # qualitative code and there are no other mismaps identified for urine
          # protein or creatinine
          upcr_criteria_met == TRUE &
            prot_mismap_to_qual == TRUE &
            prot_urine_calcium == FALSE &
            prot_blood_specimen == FALSE &
            any_creat_mismap == FALSE
        )
    ) %>%
    mutate(
      creat_value = !!rlang::sym(creat_value_name),
      prot_value = !!rlang::sym(prot_value_name)
    ) %>%
    mutate(
      upcr_calc = ifelse(
        prot_value == 0 |
          is.na(prot_value) |
          creat_value == 0 |
          is.na(creat_value),
        NA,
        round(prot_value / creat_value, digits = 5)
      ),
      unit_info_useful_calc = unit_info_useful &
        unit_info_useful_creat
    ) %>%
    select(
      person_id,
      measurement_date,
      site,
      starts_with("measurement"),
      starts_with("specimen"),
      starts_with("unit"),
      starts_with("value"),
      starts_with("parsed"),
      prot_value,
      creat_value,
      upcr_calc
    )
  
  return(upcr_calcs)
}

#' Combine codeset-provided UPCRs and calculated UPCRs for comparison
#' #'
#' @param upcr_codeset_meas Codeset-provided UPCRs
#' @param upcr_calc_meas Calculated UPCRs
#'
#' @return Table with fields from calculated (no suffix) and codeset-based
#' UPCRs (suffix _codeset). Add "multiplier" field for insight as to whether
#' one UPCR is a multiple of another (indicating a unit issue, for example)
#'
combine_codeset_calc_upcrs <-
  function(upcr_codeset_meas_tbl,
           upcr_calc_meas_tbl) {
    combined_upcrs <- upcr_calc_meas_tbl %>%
      select(
        person_id,
        site,
        starts_with("measurement"),
        starts_with("specimen"),
        starts_with("unit"),
        starts_with("value"),
        starts_with("parsed"),
        starts_with("upcr"),
        starts_with("prot"),
        starts_with("creat")
      ) %>%
      full_join(
        select(
          upcr_codeset_meas_tbl,
          person_id,
          site,
          starts_with("measurement"),
          starts_with("specimen"),
          starts_with("unit"),
          starts_with("value"),
          starts_with("parsed"),
          starts_with("upcr")
        ),
        by = c("person_id", "site", "measurement_date"),
        suffix = c("", "_codeset")
      ) %>%
      mutate(
        multiplier = if_else(
          !is.na(upcr_calc) & !upcr_calc == 0 &
            !is.na(upcr_codeset) &
            !upcr_codeset == 0,
          upcr_codeset / # add "multiplier" field for insight about whether the
            # discrepancy is due to unit info
            upcr_calc,
          0
        )
      ) %>%
      distinct()
    
    return(combined_upcrs)
  }

#' Select codeset-provided UPCR over calculated and add fields for selection of
#' UPCR based on timestamp etc. Remove selected UPCR if NA
#'
#' @param upcr_calc_or_codeset Codeset with calculated and codeset-provided
#' UPCRs
#'
#' @return Table with selected UPCRs and additional fields:
#' - n_upcrs_on_date The number of distinct UPCRs available on the date for a
#' patient
#' - prot_creat_timestamp_match Whether the timestamp of the urine protein and
#' urine creatinine measurement match
#' - prot_creat_timestamp_diff_mins Time difference (in minutes) between the
#' urine protein and urine creatinine measurement (if applicable)
#'
select_upcr <- function(upcr_calc_or_codeset) {
  selected_upcr_tbl <- upcr_calc_or_codeset %>%
    mutate(
      selected_upcr = if_else(!is.na(upcr_codeset),
                              upcr_codeset,
                              upcr_calc),
      selected_meas_id = if_else(
        !is.na(upcr_codeset),
        measurement_id_codeset,
        measurement_id
      ),
      selected_flag = if_else(!is.na(upcr_codeset),
                              "codeset_provided",
                              "calculated"),
      selected_unit_info_useful = if_else(
        !is.na(upcr_codeset),
        unit_info_useful_codeset,
        unit_info_useful_calc
      )
    ) %>%
    select(
      person_id,
      site,
      person_id,
      measurement_date,
      site,
      starts_with("selected"),
      starts_with("measurement"),
      starts_with("specimen"),
      starts_with("unit"),
      starts_with("value"),
      starts_with("upcr"),
      starts_with("prot"),
      starts_with("creat")
    ) %>%
    collect_new() %>%
    group_by(person_id, measurement_date) %>%
    mutate(n_upcrs_on_date = n_distinct(measurement_id)) %>%
    ungroup() %>%
    mutate(
      prot_creat_timestamp_match = measurement_datetime ==
        measurement_datetime_creat,
      prot_creat_timestamp_diff_mins = abs(
        as.numeric(measurement_datetime -
                     measurement_datetime_creat)
      ) / 60
    ) %>%
    filter(!is.na(selected_upcr)) %>%
    distinct() %>%
    select(
      person_id,
      site,
      measurement_date,
      selected_upcr,
      selected_meas_id,
      selected_flag,
      selected_unit_info_useful,
      upcr_calc,
      upcr_codeset,
      n_upcrs_on_date,
      prot_creat_timestamp_match,
      prot_creat_timestamp_diff_mins,
      measurement_id_prot = measurement_id,
      measurement_datetime_prot = measurement_datetime,
      measurement_concept_id_prot = measurement_concept_id,
      measurement_concept_name_prot = measurement_concept_name,
      measurement_source_value_prot = measurement_source_value,
      unit_concept_id_prot = unit_concept_id,
      unit_concept_name_prot = unit_concept_name,
      unit_source_value_prot = unit_source_value,
      unit_info_useful_prot = unit_info_useful,
      prot_value,
      value_as_number_prot = value_as_number,
      value_as_concept_id_prot = value_as_concept_id,
      value_as_concept_name_prot = value_as_concept_name,
      value_source_value_prot = value_source_value,
      specimen_concept_id_prot = specimen_concept_id,
      specimen_concept_name_prot = specimen_concept_name,
      measurement_id_creat,
      measurement_datetime_creat,
      measurement_concept_id_creat,
      measurement_concept_name_creat,
      measurement_source_value_creat,
      unit_concept_id_creat,
      unit_concept_name_creat,
      unit_source_value_creat,
      unit_info_useful_creat,
      creat_value,
      value_as_number_creat,
      value_as_concept_id_creat,
      value_as_concept_name_creat,
      value_source_value_creat,
      specimen_concept_id_creat,
      specimen_concept_name_creat,
      measurement_id_codeset,
      measurement_datetime_codeset,
      measurement_concept_id_codeset,
      measurement_concept_name_codeset,
      measurement_source_value_codeset,
      unit_concept_id_codeset,
      unit_concept_name_codeset,
      unit_source_value_codeset,
      unit_info_useful_codeset,
      value_as_concept_id_codeset,
      value_as_concept_name_codeset,
      value_source_value_codeset,
      specimen_concept_id_codeset,
      specimen_concept_name_codeset
    )
  
  return(selected_upcr_tbl)
}

#' Add datetime workaround
#'
#' Due to issue whereby its not possible to collect and write to the database
#' some dates, created function to add measurement_datetimes via measurement_id
#' for table on the database (not local) without measurement_datetime
#'
#' @param meas_wo_datetimes Measurement table without datetime (must be on db,
#' not local)
#'
#' @return Measurement table with datetime
#'
add_datetime_workaround <-
  function(meas_wo_datetimes) {
    meas_w_datetimes <- meas_wo_datetimes %>%
      inner_join(select(
        cdm_tbl("measurement"),
        measurement_id,
        measurement_datetime
      ),
      by = "measurement_id") %>%
      inner_join(
        select(
          cdm_tbl("measurement"),
          measurement_id,
          measurement_datetime
        ),
        by = c("measurement_id_creat" = "measurement_id"),
        suffix = c("", "_creat")
      )
    
    return(meas_w_datetimes)
  }

#' Prepare codeset UPCR measurements
#'
#' Parse value_source_value where there's a "<" or "less than"
#' If the unit indicated by the unit_concept_id is mg per gram,
#' convert to gram per gram by dividing by 1000
#'
#' @param upcr_codeset_meas Measurement table with codeset-provided
#' UPCR measurements
#'
#' @return UPCR measurement table with additional fields
#' - parsed_van Value as number information with values extracted from
#' value_source_value where "<" or "less than" present
#' - van_is_na Boolean for whether value as number is NA
#' - upcr_codeset Standardized/parsed UPCR codeset value
#'
prepare_codeset_upcrs <- function(upcr_codeset_meas) {
  upcr_codeset_meas %>%
    add_van_variables() %>%
    mutate(
      # 9565: mg per mg, 9074 mg per mg of creatinine
      unit_info_useful = unit_concept_id %in% c(9017, 8723, 9565, 9074) &
        !is.na(unit_concept_id),
      upcr_codeset = if_else(
        unit_concept_id %in%
          c(9017, # milligram per gram of creatinine
            8723),
        # milligram per gram
        parsed_van / 1000,
        parsed_van
      )
    )
}

#' Add value_as_number variables
#'
#' Add 2 value_as_number variables
#' van_is_na is a boolean for whether value is number is NA (TRUE for NA)
#' parsed_van parses the value_source_value field if van_is_na == TRUE or
#' value_source_value contains a < symbol (as lower bounds are appropriate
#' to use in UPCR calculations)
#'
#' @param measurement_table measurement table with value_as_number field
#'
#' @return
add_van_variables <- function(measurement_table) {
  measurement_table %>%
    mutate(
      van_is_na = is.na(value_as_number),
      parsed_van = if_else(
        van_is_na == TRUE &
          !is.na(value_source_value) &
          (
            str_detect(value_source_value, "\\<") == TRUE |
              str_detect(tolower(value_source_value), "less than") == TRUE
          ),
        as.numeric(str_extract(
          value_source_value, "\\.?\\d+\\.?\\d*"
        )),
        value_as_number
      )
    )
}

#' Add source value is range flag (src_val_is_range) as TRUE if any of the
#' following strings are included in the value_source_value: "to", "-", ">"
#' and "greater than"
#' Note: "<" and "less than" are excluded because it expected that lower
#' bounds are reported
#'
#' @param meas_tbl Measurement table with value_source_value field
#'
#' @return Measurement table with additional src_val_is_range field
#'
add_src_val_is_range_flag <- function(meas_tbl) {
  meas_tbl %>%
    mutate(
      vsv_lower_case = tolower(value_source_value),
      src_val_is_range = str_detect(vsv_lower_case, 'to') |
        str_detect(vsv_lower_case, '-') |
        str_detect(vsv_lower_case, '>') |
        str_detect(vsv_lower_case, 'greater than')
    ) %>%
    select(-vsv_lower_case)
}

#' Add flag for whether UPCR criteria are met
#' - A numeric value is populated for urine protein in parsed_van (after parsing
#'  value_as_number)
#' - The source value does not contain range information (excluding "<" and
#' "less than")
#' - A numeric value is populated for urine creatinine in parsed_van_creat
#' (after parsing value_as_number)
#'
#' @param meas_tbl A measurement table with urine protein and urine creatinine
#' measurements
#'
#' @return Measurement table with additional upcr_criteria_met field
#'
add_upcr_criteria_met_flag <- function(meas_tbl) {
  meas_tbl %>%
    mutate(upcr_criteria_met = case_when(
      !is.na(parsed_van) &
        src_val_is_range == FALSE &
        !is.na(parsed_van_creat) ~ TRUE,
      TRUE ~ FALSE
    ))
}

#' Add urine protein mismap flags
#'
#' Add flags for mismapped LOINC or CPT codes (via add_prot_mismaps_from_file)
#' by extracting codes following the "|" symbol in measurement_source_value,
#' incorporating logic for flags for urine protein measurements erroneously
#' mapped to the quantitative or qualitative/semi-qualitative codesets
#' Add additional flags for
#' - presence of strings "blood" or "serum" in
#' measurement_source_value (decision not to rely on specimen as this info
#' sometimes conflicts with measurement_source_value)
#' - presence of strings "calcium" in measurement_source_value
#'
#' @param meas_tbl A measurement table with creatinine measurement variables,
#' must include "measurement_source_value"
#' @param mismaps A file of known LOINC/CPT mismappings associated with
#' urine protein measurements
#'
#' @return meas_tbl with additional flags:
#' - prot_mismap: urine protinine mismap based on LOINC/CPT
#' - prot_blood_specimen: blood/serum in measurement_source_value
#' - any_prot_mismap: flag if there's any urine protinine mismapping
#'
add_prot_mismap_flags <- function(df,
                                  mismaps) {
  df %>%
    add_prot_mismaps_from_file(mismaps = mismaps) %>%
    mutate(
      prot_mismap_to_qual = (
        prot_mismap_to_qual |
          str_detect(measurement_source_value, "2888-6") |
          str_detect(measurement_source_value, "35663-4")
      ) &
        codeset == "qual_semiqual",
      prot_mismap_to_quant = prot_mismap_to_quant &
        codeset == "quant",
      # check for these quantitative source values anywhere in
      # measurement_source_value (not restricted to after the '|')
      prot_urine_calcium = str_detect(tolower(measurement_source_value), "calcium"),
      prot_blood_specimen = str_detect(tolower(measurement_source_value), "blood") |
        str_detect(tolower(measurement_source_value), "serum"),
      # If a mismap is NA, set it to FALSE
      prot_mismap = if_else(is.na(prot_mismap), FALSE, prot_mismap),
      prot_mismap_to_qual = if_else(is.na(prot_mismap_to_qual), FALSE, prot_mismap_to_qual),
      prot_urine_calcium = if_else(is.na(prot_urine_calcium), FALSE, prot_urine_calcium),
      prot_blood_specimen = if_else(is.na(prot_blood_specimen), FALSE, prot_blood_specimen),
      any_prot_mismap = prot_mismap |
        prot_mismap_to_qual |
        prot_mismap_to_quant |
        prot_urine_calcium |
        prot_blood_specimen
    )
}

#' Add urine creatinine mismap flags
#'
#' Add flag for mismapped LOINC or CPT codes (via add_creat_mismaps_from_file)
#' by extracting codes following the "|" symbol in measurement_source_value
#' and add additional flag for presence of strings "blood" or "serum" in
#' measurement_source_value (decision not to rely on specimen as this info
#' sometimes conflicts with measurement_source_value)
#'
#' @param meas_tbl A measurement table with creatinine measurement variables
#' suffixed with "_creat", must include "measurement_source_value_creat"
#' @param mismaps A file of known LOINC/CPT mismappings associated with
#' urine creatinine measurements
#'
#' @return meas_tbl with additional flags:
#' - creat_mismap: urine creatinine mismap based on LOINC/CPT
#' - creat_blood_specimen: blood/serum in measurement_source_value
#' - any_creat_mismap: flag if there's any urine creatinine mismapping
#'
add_creat_mismap_flags <- function(meas_tbl,
                                   mismaps) {
  meas_tbl %>%
    add_creat_mismaps_from_file(mismaps = mismaps) %>%
    mutate(
      creat_blood_specimen = str_detect(tolower(measurement_source_value_creat), "blood") |
        str_detect(tolower(measurement_source_value_creat), "serum"),
      creat_mismap = if_else(is.na(creat_mismap), FALSE, creat_mismap),
      creat_blood_specimen = if_else(is.na(creat_blood_specimen), FALSE, creat_blood_specimen),
      any_creat_mismap = creat_mismap |
        creat_blood_specimen
    )
}

#' Add flags for mismapped CPT or LOINC codes for urine
#' creatinine extracted from
#' measurement_source_value (material following the '|' symbol)
#' according to file of potentially mismapped codes identified
#' through data quality exploration
#'
#' @param meas_tbl Measurement table with _creat suffix for creatinine
#' measurements
#' @param mismaps csv file of potentially mismapped fields
#'
#' @return Measurement table with flags
#'
add_creat_mismaps_from_file <- function(meas_tbl,
                                        mismaps) {
  meas_tbl %>%
    mutate(creat_cpt_or_loinc = str_extract(measurement_source_value_creat, "[^|]*$")) %>%
    left_join(
      select(mismaps, cpt_or_loinc, creat_mismap),
      by = c("creat_cpt_or_loinc" = "cpt_or_loinc"),
      suffix = c("", "_creat")
    ) %>%
    mutate(creat_mismap = if_else(creat_mismap == 0 |
                                    is.na(creat_mismap) , FALSE, TRUE))
}

#' Add flags for mismapped CPT or LOINC codes for urine
#' protein extracted from
#' measurement_source_value (material following the '|' symbol)
#' according to file of potentially mismapped codes identified
#' through data quality exploration
#'
#' @param meas_tbl Measurement table with no suffix for urine protein fields
#' @param mismaps csv file of potentially mismapped fields
#'
#' @return Measurement table with flags
#'
add_prot_mismaps_from_file <- function(meas_tbl,
                                       mismaps) {
  meas_tbl %>%
    mutate(prot_cpt_or_loinc = str_extract(measurement_source_value, "[^|]*$")) %>%
    left_join(
      select(
        mismaps,
        cpt_or_loinc,
        prot_mismap,
        prot_mismap_to_quant,
        prot_mismap_to_qual
      ),
      by = c("prot_cpt_or_loinc" = "cpt_or_loinc"),
      suffix = c("", "_prot")
    ) %>%
    mutate(
      prot_mismap = if_else(prot_mismap == 0 |
                              is.na(prot_mismap) , FALSE, TRUE),
      prot_mismap_to_quant = if_else(
        prot_mismap_to_quant == 0 |
          is.na(prot_mismap_to_quant) ,
        FALSE,
        TRUE
      ),
      prot_mismap_to_qual = if_else(
        prot_mismap_to_qual == 0 |
          is.na(prot_mismap_to_qual) ,
        FALSE,
        TRUE
      )
    )
}

#' Convert value_as_number and parsed_van to mg/dl
#' based on information in unit_concept_id
#'
#' @param meas_table Measurement table with
#' value_as_number and parsed_van fields
#'
#' @return Measurement table with additional
#' value_as_number_convert and parsed_van_convert fields
#' where value is converted to mg/dl based on
#' unit_concept_id
#'
convert_to_mg_per_dl <- function(meas_table) {
  meas_table %>%
    mutate(
      unit_info_useful = unit_concept_id %in% c(8840, 8713, 8751, 8861) &
        !is.na(unit_concept_id),
      value_as_number_convert = case_when(
        unit_concept_id == 8840 ~ value_as_number,
        # mg/dl
        unit_concept_id == 8713 ~ value_as_number * 1000,
        # g/dl
        unit_concept_id == 8751 ~ value_as_number / 10,
        # mg/l
        unit_concept_id == 8861 ~ value_as_number * 100,
        # mg/ml
        TRUE ~ value_as_number
      ),
      parsed_van_convert = case_when(
        unit_concept_id == 8840 ~ parsed_van,
        # mg/dl
        unit_concept_id == 8713 ~ parsed_van * 1000,
        # g/dl
        unit_concept_id == 8751 ~ parsed_van / 10,
        # mg/l
        unit_concept_id == 8861 ~ parsed_van * 100,
        # mg/ml
        TRUE ~ value_as_number
      )
    )
}

#' HSP Specific UPCR Clean-up
#' This function cleans and annotates UPCR for the HSP Project.
#' The following steps are done:
#' - derive units where determined not to be useful
#' - adjust the UPCR value to the mg/mg equivalent where not present
#' - remove any potential hour (12/24) mismapped calculations
#' - derive a unit category useful to the project
#' 
#' @param upcr_table UPCR Table generated from the run_all_upcr_steps function
#'
#' @return UCPR table with adjusted mg/mg equivalents and upcr categories
#'
clean_and_annotate_upcr<-function(upcr_table){
  upcr_table%>%
    mutate(upcr_calc=round(upcr_calc,2))%>%
    mutate(derived_unit= case_when(
      (selected_flag=='codeset_provided' & str_detect(unit_source_value_codeset,'mg/g'))  ~ 'mg/g',  
      
      (selected_flag=='calculated' & str_detect(tolower(unit_source_value_prot),'mg') & 
         str_detect(tolower(unit_source_value_creat),'mg')) ~ 'mg/mg',
      (selected_flag=='calculated' & str_detect(tolower(unit_source_value_prot),'mg') & 
         str_detect(tolower(unit_source_value_creat),'g')) ~ 'mg/g'))%>%
    mutate(adjusted_upcr = if_else(str_detect(derived_unit,'mg/g'),
                                   selected_upcr*0.001,selected_upcr
    ))%>%
    mutate(derived_unit=if_else(str_detect(derived_unit,'mg/g'),as.character('mg/mg'),derived_unit))%>% 
    mutate(selected_upcr=round(selected_upcr,2))%>%  
    mutate(final_upcr= round(case_when((selected_upcr==upcr_calc & str_detect(unit_source_value_codeset, 'mg/g cr')) ~ selected_upcr,
                                       (is.na(upcr_calc) & str_detect(unit_source_value_codeset, 'mg/g cr')) ~ selected_upcr,
                                       !is.na(adjusted_upcr) ~ adjusted_upcr,
                                       TRUE ~ selected_upcr
    ),2),
    final_unit=coalesce(derived_unit,unit_source_value_codeset))%>%
    mutate(hr_value_mismap= coalesce(if_else( selected_flag=='calculated' &
                                                ((str_detect(tolower(measurement_source_value_creat),'24') & !str_detect(tolower(measurement_source_value_prot),'24'))| (!str_detect(tolower(measurement_source_value_creat),'24') & str_detect(tolower(measurement_source_value_prot),'24')) |
                                                   (!str_detect(tolower(measurement_source_value_creat),'12') & str_detect(tolower(measurement_source_value_prot),'12'))|
                                                   (str_detect(tolower(measurement_source_value_creat),'12') & !str_detect(tolower(measurement_source_value_prot),'12')))
                                              ,TRUE,FALSE
    ),FALSE))%>%
    mutate(upcr_category=case_when(
      (final_upcr <0.5) ~ ' <0.5',
      (final_upcr>=0.5 & final_upcr<=2) ~ '0.5-2',
      (final_upcr>2) ~ '>2'
      
    ))%>%
    filter(!hr_value_mismap)
  
}
