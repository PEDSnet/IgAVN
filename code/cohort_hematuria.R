#' Function to derive positive versus negative hematuria urinalysis measurements
#'
#' This function parses value_as_concept_name, value_as_number, and
#' value_source_value fields to derive whether hematuria urinalysis measurements should
#' be classfied as positive or negative.
#' 
#' Assumptions:
#' 
#' Hematuria is typically classified as >=5 rbs per hpf.
#' 
#' Negative:
#'
#'- <5
#'- Negative
#'- Trace
#'- None
#'- None Seen
#'- 0-5
#'- 3-5
#'- Rare
#'
#'Positive:
#'  
#'  - >=5
#'- Positive
#'- Moderate (mod) (Includes hemolyzed results)
#'- Small (Includes hemolyzed results)
#'- Large (Includes hemolyzed results)
#'- TNTC (Too numerous to count)
#'- Trace Moderate/Small
#'
#' @param hematuria_urinalysis_measurements Collected (i.e. local) table with urine blood measurements
#'
#' @return Table of urine blood or hematuria urinalysis measurements with additional
#' derived_result variable which can be positive, negative or have an "unable to derive" indicator for 
#' ambiguous results
#'
derive_hematuria_urinalysis <- function(hematuria_urinalysis_measurements) {
  hematuria_urinalysis_measurements%>%
    collect()%>%
    mutate(vsv_lower_case = tolower(value_source_value),
           vsv_range_low = as.numeric(ifelse(str_detect(vsv_lower_case,'-|to'),(stringr::str_split_n(vsv_lower_case,'-|to',1)),NA)), 
           vsv_range_high= as.numeric(ifelse(str_detect(vsv_lower_case,'-|to'),(stringr::str_split_n(vsv_lower_case,'-|to',2)),NA)),
           vsv_less_than= as.numeric(ifelse(str_detect(vsv_lower_case,'<|less than'),(stringr::str_split_n(vsv_lower_case,'<|less than',2)),NA)),
           vsv_greater_than= as.numeric(ifelse(str_detect(vsv_lower_case,'>|greater than'),(stringr::str_split_n(vsv_lower_case,'>|greater than',2)),NA))
    )%>%
    mutate(derived_result= case_when(
      str_detect(value_as_concept_name,'\\+') ~ "POSITIVE",
      (str_detect(vsv_lower_case,"1\\+|2\\+|2\\+\\+|3\\+|4\\+|\\+|moderate|mod|med|large|larg|lg|lage|lrg|sm|small|smal|samll|positive|pos|sm heme|plus|tntc|too numerous|few|innumerable|many")
       & !str_detect(vsv_lower_case,'pos.*trace|trace \\+|^+.*trace|trace pos')) ~ 'POSITIVE',
      str_detect(value_as_concept_name,'Absent|Trace') ~ "NEGATIVE",
      str_detect(vsv_lower_case,'cancel|pending|insuff|clotted|unable|note|error|not sufficient|not done|specimen|not avail|clin lab') ~ "UNABLE TO DERIVE",
      str_detect(vsv_lower_case,'negative|none|trace|tr|tra|neg|not|rare|ng|absent|non seen|neagtive|neagative|occ|occasional') ~ "NEGATIVE",
      vsv_less_than <=5 ~ "NEGATIVE",
      vsv_greater_than >=5 ~ "POSITIVE",
      vsv_range_low >=5 ~ "POSITIVE", 
      vsv_range_high<=5 ~ "NEGATIVE",
      #DQA Issue observed where results reported as a range combines numbers. E.g. 2-4 becomes 24, 3 to 10 becomes 310.
      (!str_detect(vsv_lower_case,'-|to') & value_as_number >=5) ~ "POSITIVE",
      (!str_detect(vsv_lower_case,'-|to') & value_as_number <5) ~ "NEGATIVE",
      TRUE ~ "UNABLE TO DERIVE" )
    )
}
#' Function to annotate a set of urine blood measurements into dipstick or microscopy based on the result
#' @param urine_blood_measurements tbl with urine blood lab results
#' @return tbl with all columns in the urine_blood_measurements table plus the annotated result_type
annotate_urine_blood_labs <- function(urine_blood_measurements) {
  urine_blood_measurements%>%
    collect()%>%
    mutate(vsv_lower_case = tolower(value_source_value))%>%
    mutate(result_type= case_when(
      str_detect(vsv_lower_case,'negative|neg|trace|tr|rare|1\\+|2\\+|2\\+\\+|3\\+|4\\+|\\+|moderate|mod|med|large|lg|lage|larg|lrg|sm|small|smal|samll|positive|pos|sm heme|plus|neag|ng|neagtive|neagative') ~ "dipstick",
      str_detect(vsv_lower_case,'-|to|<|less than|>|greater than|tntc|none seen|none|too numerous|\\d|occasional|occ|many|seen|innumerable|few|absent') ~ "microscopy",
      TRUE ~ "UNABLE TO DETERMINE" )
    )
}

#' Function to apply a hierarchy to determine hematuria based on a combination of dipstick and microscopy test results.
#' 
#' Rules for the hierarchy are as follows:
#' 
#' * For patients with only dipstick results, use the derived_result to determine hematuria.
#' * For patients with micro only, use the derived_result to determine hematuria.
#' * For patients with a dipstick and micro:
#'  + If the dipstick is positive and (micro positive or unable to derive due to range) -> code positive
#'  + If the dipstick is negative and (micro negative or unable to derive due to range) -> code negative
#'  + If the dipstick is positive and micro negative -> code negative (e.g dipstick = "small" and micro = "0-2")
#'  + If the dipstick is negative and micro positive -> code positive (e.g dipstick = "trace" and micro = "5-10")
#' 
#' @param derived_hematuria_results tbl with urine blood lab results post annotation for the urine blood test type and hematuria derivation
#' @param urine_blood_timing variable to determine if the urine blood measurement is the FIRST or LAST. If it is of interest to determine
#' if a person EVER had hematuria within a given time period the option of EVER can be passed.
#' @return tbl with person id and anotation of positive or negative for hematuria
apply_hematuria_hierarchy <-function(derived_hematuria_results,
                                     urine_blood_timing){
  
  dipstick_results<-derived_hematuria_results%>%
    filter(result_type=='dipstick')
  
  microscopy_results<-derived_hematuria_results%>%
    filter(result_type=='microscopy')
  
  micro_no_dipstick<-microscopy_results%>%
    anti_join(select(dipstick_results,person_id),by='person_id')%>%
    filter(derived_result!='UNABLE TO DERIVE')%>%
    select(person_id,measurement_datetime,derived_result)%>%
    distinct()
  
  dipstick_comparison<-dipstick_results%>%
    left_join(microscopy_results,by='person_id',suffix=c('','_micro'))%>%
    mutate(in_24_hr= ifelse(between(measurement_datetime_micro,measurement_datetime,sql("(measurement_datetime + interval '24 hour')")),'Y','N'))%>%
    select(person_id,measurement_id,measurement_concept_id,measurement_source_value,measurement_datetime,
           result_type,vsv_lower_case,derived_result,
           measurement_id_micro,measurement_concept_id_micro,measurement_source_value_micro,measurement_datetime_micro,result_type_micro,vsv_lower_case_micro,derived_result_micro,in_24_hr)%>%
    distinct()%>%
    mutate(micro_avail = ifelse(!is.na(measurement_id_micro)&&in_24_hr=='Y','Y','N'))%>%
    group_by(measurement_id)%>%
    mutate(num_micro=row_number())%>%
    ungroup()
  
  dipstick_combined<-dipstick_comparison%>%
    filter(num_micro==1)%>%
    select(person_id,measurement_datetime,vsv_lower_case,derived_result,vsv_lower_case_micro,derived_result_micro,in_24_hr)%>%
    mutate(hierarchy_derived_result = case_when(
      derived_result=='POSITIVE' && (is.na(in_24_hr)|in_24_hr=='N') ~ 'POSITIVE',
      derived_result=='NEGATIVE' && (is.na(in_24_hr)|in_24_hr=='N') ~ 'NEGATIVE',
      derived_result=='POSITIVE' && (derived_result_micro=='POSITIVE' | derived_result_micro=='UNABLE TO DERIVE') && (in_24_hr=='Y')  ~ 'POSITIVE',
      derived_result=='NEGATIVE' && (derived_result_micro=='NEGATIVE' | derived_result_micro=='UNABLE TO DERIVE') && (in_24_hr=='Y') ~ 'NEGATIVE',
      derived_result=='POSITIVE' && derived_result_micro=='NEGATIVE' && (in_24_hr=='Y') ~ 'NEGATIVE',
      derived_result=='NEGATIVE' && derived_result_micro=='POSITIVE' && (in_24_hr=='Y') ~ 'POSITIVE',
      TRUE ~ 'UNABLE TO CATEGORIZE'
    ))%>%
    distinct()

  if(urine_blood_timing=='ever'){  
    micro_no_dipstick<-micro_no_dipstick%>%
      group_by(person_id)%>%
      summarise(derived_result_final=max(derived_result))%>%
      mutate(type='micro_only')
    
    dipstick_combined<-dipstick_combined%>%
      group_by(person_id)%>%
      summarise(derived_result_final=max(hierarchy_derived_result))%>%
      mutate(type='dipstick_combined')
    
    micro_no_dipstick%>%
      dplyr::union(dipstick_combined)%>%
      select(person_id,derived_result_final)
  }
  else if(urine_blood_timing=='first'){
    
    micro_no_dipstick<-micro_no_dipstick%>%
      group_by(person_id)%>%
      arrange(measurement_datetime)%>%
      filter(row_number()==1)%>%
      mutate(derived_result_final=derived_result)%>%
      mutate(type='micro_only')
    
    dipstick_combined<-dipstick_combined%>%
      group_by(person_id)%>%
      arrange(measurement_datetime)%>%
      filter(row_number()==1)%>%
      mutate(derived_result_final=hierarchy_derived_result)%>%
      mutate(type='dipstick_combined')
    
    micro_no_dipstick%>%
      dplyr::union(dipstick_combined)%>%
      select(person_id,derived_result_final)

  }
  else if(urine_blood_timing=='last'){
    
    micro_no_dipstick<-micro_no_dipstick%>%
      group_by(person_id)%>%
      arrange(measurement_datetime)%>%
      filter(row_number()==n())%>%
      mutate(derived_result_final=derived_result)%>%
      mutate(type='micro_only')
    
    dipstick_combined<-dipstick_combined%>%
      group_by(person_id)%>%
      arrange(measurement_datetime)%>%
      filter(row_number()==n())%>%
      mutate(derived_result_final=hierarchy_derived_result)%>%
      mutate(type='dipstick_combined')
    
    micro_no_dipstick%>%
      dplyr::union(dipstick_combined)%>%
      select(person_id,derived_result_final)
    
  }
}