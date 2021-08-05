#' Function to derive positive versus negative proteinuria urinalysis
#' measurements
#' NOTE: This function is a first pass and requires updates.
#'
#'
#'N.B. Change concept_names to concept_id (for now concept_name enhances readability)
#'
#' This function parses value_as_concept_name, value_as_number, and
#' value_source_value fields to derive a variable called derived_result
#' which classifies proteinuria urinalysis measurements as "POSITIVE",
#' "NEGATIVE", "UNAVAILABLE" (e.g. due to technical problem with testing), or
#' "UNABLE TO DERIVE" for measurements which cannot be classified with the
#' current approach. Trace is classified as "NEGATIVE".
#'
#' Information about measurements should be classified:
#' Negative – corresponds to 0 mg/dL
#' Trace – corresponds to 15-<30 mg/dL (classified as Negative)
#' Positive -
#'     1+ - corresponds to 30-<100 mg/dL
#'     2+  - corresponds to 100-<300 mg/dL
#'     3+ - corresponds to 300-<1000 mg/dL
#'     4+ - corresponds to 1000-<2000 mg/dL
#'
#' @param proteinuria_urinalysis_measurements Collected (i.e. local) table with
#' urine protein or proteinuria urinalysis measurements
#'
#' @return Table of urine protein or proteinuria urinalysis measurements with
#' an additional derived_result variable with the following categories:
#'
#' POSITIVE: Positive urine protein/proteinuria urinalysis
#' NEGATIVE: Negative urine protein/proteinuria urinalysis
#' UNAVAILABLE: Result unavailable (e.g. due to problem with testing)
#' UNABLE TO DERIVE: Current approach to classifying measurements cannot derive
#' result
#'
derive_proteinuria_urinalysis <-
  function(proteinuria_urinalysis_measurements) {
    proteinuria_urinalysis_measurements %>%
      collect() %>% 
      mutate(
        vsv_lower_case = tolower(value_source_value),
        derived_result_cont = as.numeric(str_extract(value_source_value, "\\.?\\d+\\.?\\d*")),
        derived_result_category = case_when(
          value_as_concept_id == 9190L | # "Not detected"
            value_as_concept_id == 4124457L | # "Normal range"
            value_as_concept_id == 9189L | # "Negative" 
            value_as_concept_id == 4125559L | # "Nothing"
            value_as_concept_id == 4188540L | # "No"
            value_as_concept_id == 4124462L ~ "NEGATIVE", # "None" 
          value_as_concept_id == 4123508L | value_source_value == "+" ~ "1+", # "+"
          value_as_concept_id == 4126673L | value_source_value == "++" ~ "2+", # "++"
          value_as_concept_id == 4125547L | value_source_value == "+++" ~ "3+", # "+++"
          value_as_concept_id == 4126674L | value_source_value == "++++" ~ "4+", # "++++"
          value_as_concept_id == 9192L  ~ "TRACE", # "Trace"
          derived_result_cont < 30 ~ "NEGATIVE",
          derived_result_cont >= 30 &
            derived_result_cont < 100 ~ "1+",
          derived_result_cont >= 100 &
            derived_result_cont < 300 ~ "2+",
          derived_result_cont >= 300 &
            derived_result_cont < 1000 ~ "3+",
          derived_result_cont >= 1000 & 
            derived_result_cont <= 2000 ~ "4+",
          derived_result_cont > 2000 ~ "UNABLE TO DERIVE", 
          str_detect(
            str_replace(vsv_lower_case, " ", ""),
            "1\\+") ~ '1+',
          str_detect(
            str_replace(vsv_lower_case, " ", ""),
            "2\\+|2\\+\\+") ~ '2+',
          str_detect(
            str_replace(vsv_lower_case, " ", ""),
            "3\\+") ~ '3+',
          str_detect(
            str_replace(vsv_lower_case, " ", ""),
            "4\\+") ~ '4+',
          str_detect(
            vsv_lower_case,
            "^tra|^tace$|^tr$|tarce|treace|^traqc$|^trace|ace$" 
          ) ~ "TRACE",
          str_detect(
            vsv_lower_case,
            "^ni$|can[']t|notified|duplicate|can not|leak|improper|pend|does not|incorr|error|not performed|contaminat|cancel|unable|refuse|insuff|not suff|wrong|omit|illegible|ilegible|not rec|recollect|not rep|no data|no result|not avail|no sample|no evid|no specimen|no suitable|not able|not done|not calc|poss|^qns.*|^[+][-]$" 
          ) | str_length(vsv_lower_case) > 30 ~ "UNABLE TO DERIVE",
          str_detect(
            vsv_lower_case,
            "neg|^ng$|^meg$|^beg$|^nrg$|^neag|^net$|^ne$|^nge$|^neh$|netative|ngative|ngeative|^neag|^hegative|^nag|^neative|^neb$|^nebg|^neeg|^nheg|gative$|^[-]$" 
          ) & !str_detect(vsv_lower_case, "ssa") ~ "NEGATIVE",
          str_detect(
            vsv_lower_case,
            "^plus$|^pos$|^positive$|^\\d[+]$|^large$|^present$|^positive|^postitive$|^postive$|^pos[.]$|^pos[*]$|^small$" 
          ) ~ "POSITIVE",
          str_detect(
            vsv_lower_case,
            "^[+]{1}$|^[\\*]{1}$|^[1][+]{1}|[1][t]|plus.*[1]|pos.*[1]|[+]{1}[1]$|grade.*[1]|[+].*[1]$|[1].*plus|[1].*positive|.* [+]{1} .*"  
          ) | vsv_lower_case == "1+" ~ "1+",
          str_detect(
            vsv_lower_case,
            "^[+]{2}$|^[\\*]{2}$|(100-200)|^[2][+]{1}|[2][t]|plus.*[2]|pos.*[2]|[+]{1}[2]$|grade.*[2]|[+].*[2]$|[2].*plus|[2].*positive|.* [+]{2} .*"
          )  ~ "2+",
          str_detect(
            vsv_lower_case,
            "^[+]{3}$|^[\\*]{3}$|^[3][+]|^[3][+]{1}|3t|plus.*[3]|pos.*[3]|[+]{1}[3]$|grade.*[3]|[+].*[3]$|[3].*plus|[3].*positive|.* [+]{3} .*"
          ) ~ "3+",
          str_detect(
            vsv_lower_case,
            "^[+]{4}$|^[\\*]{4}$|^[4][+]{1}|[4][t]|plus.*[4]|pos.*[4]|[+]{1}[4]$|grade.*[4]|[+].*[4]$|[4].*plus|[4].*positive|.*[+][4]$|.* [+]{4} .*"
          ) ~ "4+",
          value_source_value == "1+" ~ "1+",
          value_source_value == "2+" ~ "2+",
          value_source_value == "3+" ~ "3+",
          value_source_value == "4+" ~ "4+",
          TRUE ~ "UNABLE TO DERIVE"
        )) %>% 
      # stubborn ones 
      mutate(derived_result_category =
               case_when(
                 str_detect(
                   value_source_value,
                   "1[+]|[+]1|[+].1"
                 ) ~ "1+",
                 str_detect(
                   value_source_value,
                   "2[+]|[+]2|[+].2" 
                 ) ~ "2+",
                 str_detect(
                   value_source_value,
                   "3[+]|[+]3|[+].3"
                 ) ~ "3+",
                 str_detect(
                   value_source_value,
                   "4[+]|[+]4|[+].4"
                 ) ~ "4+",
                 TRUE ~ derived_result_category
               )) %>% 
      # the final stubborn ones 
      mutate(derived_result_category = 
               case_when(
                 is.na(as.numeric(vsv_lower_case)) & 
                   derived_result_cont < 30 ~ "NEGATIVE",
                 is.na(as.numeric(vsv_lower_case)) & 
                   derived_result_cont >= 30 &
                   derived_result_cont < 100  ~ "1+",
                 is.na(as.numeric(vsv_lower_case)) & 
                   derived_result_cont >= 100 &
                   derived_result_cont < 300  ~ "2+",
                 is.na(as.numeric(vsv_lower_case)) & 
                   derived_result_cont >= 300 &
                   derived_result_cont < 1000 ~ "3+",
                 is.na(as.numeric(vsv_lower_case)) & 
                   derived_result_cont >= 1000  & 
                   derived_result_cont <= 2000 ~ "4+",
                 derived_result_cont > 2000 ~ "UNABLE TO DERIVE",
                 str_detect(
                   vsv_lower_case, "qns"
                 ) ~ "UNABLE TO DERIVE",
                 TRUE ~ derived_result_category)) %>% 
      mutate(derived_result_category = 
               case_when(
                 str_detect(
                   vsv_lower_case, "plus one|plus 1|1[+]$") |
                   str_detect(
                     value_source_value, "[+]1$|[+].1$|1.[+]$" )~ "1+",
                 vsv_lower_case == "++" |
                   str_detect(
                     vsv_lower_case, "2[+]$"
                   ) |
                   str_detect(
                     value_source_value, "[+]2$|[+].2$|2.[+]$"
                   ) ~ "2+",
                 str_detect(
                   vsv_lower_case, "3[+]$"
                 ) |
                   str_detect(
                     value_source_value, "[+]3$|[+].3$|3.[+]$"
                   ) ~ "3+",
                 str_detect(
                   vsv_lower_case, "4[+]$"
                 ) | 
                   str_detect(
                     value_source_value, "[+]4$|[+].4$|4.[+]$"
                   ) ~ "4+",
                 TRUE ~ derived_result_category)) %>% 
      mutate(
        derived_result_pos_neg = case_when(
          derived_result_category == "TRACE" ~ "NEGATIVE",
          derived_result_category == "1+" | 
            derived_result_category == "2+" | 
            derived_result_category == "3+" |
            derived_result_category == "4+" ~ "POSITIVE",
          TRUE ~ derived_result_category
        )) %>% 
      mutate(
        derived_result_category_int = case_when(
          derived_result_category == "1+" ~ 1,
          derived_result_category == "2+" ~ 2,
          derived_result_category == "3+" ~ 3,
          derived_result_category == "4+" ~ 4
        ))
  }