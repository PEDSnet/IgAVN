#' Function to calculate blood pressure z score for age/sex/height using NHBPEP task force fourth report norms
#' https://www.nhlbi.nih.gov/files/docs/resources/heart/hbp_ped.pdf
#' restricts BP measurements to the correct unit_concept_id (mmHg)
#' @param meas_ht_tbl tbl in the format of the PEDSnet CDM table measurement 
#' containing height z-score measurements. Defaults to the PEDSnet CDM table measurement_anthro, where height z-scores are expected
#' with at least the cols: value_as_number, measurement_concept_id, measurement_date
#' @param min_ht_z the minimum height z score that will be acceptable to be used as the closest height for the derivation (greater than or equal to)
#' @param max_ht_z the maximum height z score that will be acceptable to be used as the closest height for the derivation (less than or equal to)
#' @param meas_bp_tbl tbl in the format of the PEDSnet CDM table measurement 
#' containing bp measurements. Defaults to the PEDSnet CDM table measurement_vitals, where BP values are expected
#' with at least the cols: value_as_number, unit_concept_id, measurement_concept_id, measurement_date
#' may be the same as the meas_ht_tbl
#' @param person_tbl tbl in the format of the PEDSnet CDM table person
#' with at least the cols: person_id, gender_concept_id, birth_date, birth_datetime
#' @param max_meas_dd integer value: maximum number of days between height z-score and BP measurement
#' @param bp_codes list of measurement_concept_ids for which to derive z-scores
#' defaults to all SBP and DBP concepts for all positions
#' @param min_age minimum age (in days) for which to calculate a BP
#' defaults to 365
#' @param max_age maximum age (in days) for which to calculate a BP
#' defaults to 365.25*18
#' @param insert_height boolean indicating whether to assert 50th percentile for height (z score of 0) if there is no height within the specified number of days around the BP measurement. if TRUE, the columns in the output for height_z_date, height_z_measurement_id, num_daily_height, diff_days, abs_diff will be NULL if a height value_as_number is inserted (i.e. not found within window and assumed to be 0). defaults to FALSE
#' @return tbl with a BP z-score for each blood pressure (i.e. one BP z-score row for each BP measurement_id) with the cols: 
#' person_id	
#' height_z_date: date of closest height z-score measurement chosen for the derivation	
#' height_z_value: value_as_number of the closest height z-score chosen for the derivation	
#' height_z_measurement_id: measurement_id of the closest height z-score chosen for the derivation
#' measurement_date: date of the blood pressure value
#' measurement_id: measurement_id of the blood pressure value  
#' measurement_concept_id: measurement_concept_id of the blood pressure value	
#' bp: value_as_number for the blood pressure value from the measurement table	
#' type: "Systolic" or "Diastolic" based on the measurement_concept_id	
#' diff_days: number of days between the height z-score and the blood pressure value. Negative indicates height z-score prior to the measurement date. Positive indicates height z-score after the measurement date	
#' abs_diff: absolute value of the number of days between the height z-score and the blood pressure value	
#' gender_concept_name: person's gender according to the person_tbl
#' measurement_age: person's age, in days, on the date of the blood pressure value	
#' bp_z: calculated blood pressure z-score
#' ht_match_type: string indicating whether or not a height was found within the window and used in the calculation ('real') or inserted as 0 if there was no height found within the window and the insert_height parameter is TRUE
#' num_daily_ht: the total number of heights that the patient had on the height_z_date. a central value_as_number will be chosen for the day so that there is only one height for the derivation of each blood pressure z-score
blood_pressure_z_scores <- function(meas_ht_tbl = cdm_tbl('measurement_anthro'),
                                    min_ht_z = -100L,
                                    max_ht_z = 100L,
                                    meas_bp_tbl = cdm_tbl('measurement_vitals'),
                                    person_tbl,
                                    max_meas_dd,
                                    bp_codes = c(3034703L,	#Diastolic Blood Pressure - Sitting
                                                 3019962L,	#Diastolic Blood Pressure - Standing
                                                 3013940L,	#Diastolic Blood Pressure - Supine
                                                 3012888L,	#Diastolic BP Unknown/Other
                                                 3018586L,	#Systolic Blood Pressure - Sitting
                                                 3035856L,	#Systolic Blood Pressure - Standing
                                                 3009395L,	#Systolic Blood Pressure - Supine
                                                 3004249L), #Systolic BP Unknown/Other
                                    min_age = 365,
                                    max_age = 365.25*18,
                                    insert_height = FALSE) {
  
  heights <- meas_ht_tbl %>%
    filter(measurement_concept_id == 2000000042L &
             value_as_number >= min_ht_z &
             value_as_number <= max_ht_z) %>%
    select(person_id, measurement_date, value_as_number, measurement_id) %>%
    group_by(measurement_date, person_id) %>%
    add_count(name = "num_daily_ht") %>% # find the total number of heights a person had on each day (and add, for each of those measurements i.e. for each row, a column with the person's total number for the day)
    arrange(value_as_number)%>% # order them in ascending order for each person+date
    filter(row_number()==ceiling(num_daily_ht/2)) %>% # Find the middle row - if there are an odd number, will choose the exact middle (and if 1, will choose the 1), if there is an even number, choose the higher? - median not supported against the db
    ungroup()%>%
    rename(height_z_date = measurement_date,
           height_z_value = value_as_number,
           height_z_measurement_id = measurement_id) %>%
    compute_new(indexes = list('person_id'))
      
  
  bps <- meas_bp_tbl %>%
    filter(measurement_concept_id %in% !!bp_codes) %>%
    filter(unit_concept_id == 8876L) %>%
    select(person_id, measurement_date, measurement_id, value_as_number, measurement_concept_id) %>%
    rename(bp = value_as_number) %>%
    mutate(type = case_when(measurement_concept_id %in% c(3034703L,	#Diastolic Blood Pressure - Sitting
                                                          3019962L,	#Diastolic Blood Pressure - Standing
                                                          3013940L,	#Diastolic Blood Pressure - Supine
                                                          3012888L) ~ 'Diastolic',
                            measurement_concept_id %in% c(3018586L,	#Systolic Blood Pressure - Sitting
                                                          3035856L,	#Systolic Blood Pressure - Standing
                                                          3009395L,	#Systolic Blood Pressure - Supine
                                                          3004249L) ~ 'Systolic')) %>%
    compute_new(indexes = list('person_id'))
  
  person_ages <- person_tbl %>%
    mutate(birth_dt = if_else(is.na(birth_date), sql('cast("birth_datetime" as date)'), birth_date)) %>%
    mutate(gender_concept_name = case_when(gender_concept_id == 8532L ~ "FEMALE",
                                           gender_concept_id == 8507L ~ "MALE")) %>%
    select(person_id, birth_dt, gender_concept_name)
  
  meas_closest <- heights %>%
    inner_join(bps, by = 'person_id') %>%
    mutate(diff_days = height_z_date - measurement_date) %>%
    mutate(abs_diff = abs(diff_days)) %>%
    filter(abs_diff <= max_meas_dd) %>%
    group_by(measurement_id) %>%
    filter(abs_diff==min(abs_diff, na.rm=TRUE))%>% # find closest height, in days
    filter(diff_days==min(diff_days, na.rm=TRUE))%>% # if there is a tie, favor height before BP
    ungroup() %>%
    mutate(ht_match_type = as.character('real')) %>%# is this neccessary? could also know this by looking at whether there is a height_z_date and height_z_measurement_id
    compute_new(indexes = list('person_id', 'measurement_id'))
  
  # insert a height z-score of 0 if specified to do so
  if(insert_height) {
    meas_all <- bps %>%
      anti_join(meas_closest, by = 'measurement_id') %>%
      mutate(height_z_value = 0,
             ht_match_type = as.character('inserted')) %>%
      dplyr::union(meas_closest) %>%
      inner_join(person_ages, by = 'person_id')
  }else {
    meas_all <- meas_closest %>%
      inner_join(person_ages, by = 'person_id')
  }
  
  bp_tbl <- meas_all %>%
    mutate(measurement_age = measurement_date - birth_dt) %>%
    filter(measurement_age >= min_age &
             measurement_age <= max_age) 
    
  # create lookup table
  BPcoeff <- tribble(
    ~gender_concept_name, ~SAexp, ~SB1, ~SB2, ~SB3, ~SB4, ~SG1, ~SG2, ~SG3, ~SG4, ~SSIG, ~DAexp, ~DB1, ~DB2, ~DB3, ~DB4, ~DG1, ~DG2, ~DG3, ~DG4, ~DSIG,
    "MALE", 102.19768, 1.82416, 0.12776, 0.00249, -0.00135, 2.73157, -0.19618, -0.04659, 0.00947, 10.7128, 61.01217, 0.68314, -0.09835, 0.01711, 0.00045, 1.46993, -0.07849, -0.03144, 0.00967, 11.6032,
    "FEMALE", 102.01027, 1.94397, 0.00598, -0.00789, -0.00059, 2.03526, 0.02534, -0.01884, 0.00121, 10.4855, 60.50510, 1.01301, 0.01157, 0.00424, -0.00137, 1.16641, 0.12795, -0.03869, -0.00079, 10.9573
  ) %>%
    pivot_longer(cols = -gender_concept_name) %>%
    mutate(type = case_when(grepl("^S", name) ~ "Systolic", TRUE ~ "Diastolic")) %>%
    mutate(name = substr(name, 2, 20)) %>%
    pivot_wider(id_cols = c(gender_concept_name, type)) %>%
    mutate(gender_concept_name = as.character(gender_concept_name))

  tmp <- bp_tbl %>%
    inner_join(BPcoeff, by = c('gender_concept_name', 'type'), copy = TRUE) %>%
    compute_new()

  tmp %>%
    mutate(ageyrs_m10 = measurement_age / 365.25 - 10.0) %>%
    # bp_exp: expected BP for the given age, height, sex
    mutate(bp_exp = Aexp + B1 * ageyrs_m10 + B2 * ageyrs_m10**2 + B3 * ageyrs_m10**3 + B4 * ageyrs_m10**4 +
      G1*height_z_value + G2 * height_z_value **2 + G3 * height_z_value **3 + G4 * height_z_value**4) %>%
    mutate(bp_z = (bp - bp_exp) / SIG) %>%
    dplyr::select(-Aexp, -B1, -B2, -B3, -B3, -B4, -G1, -G2, -G3, -G4, -SIG, -bp_exp, -ageyrs_m10, -birth_dt) %>%
    compute_new()
}
