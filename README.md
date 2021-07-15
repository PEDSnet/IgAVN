# IgAVN
This repository contains code and data used to create the cohort of children with IgA vasculitis with nephritis and assess their outcome.

## Cohort Development

#### IgAVN Diagnosis
- [SNOMED-CT codes](codesets/hsp_codeset.csv) were used to identify all patients with a diagnosis of IgAV/HSP on or after January 1, 2009. Patients diagnosed with another form of vasculitis were excluded. In order to focus on outcomes related to a patient’s IgAV and not attributable to other chronic medical conditions, we excluded patients with evidence of non-renal chronic disease. To identify chronic disease, we applied the taxonomy from the Pediatric Medical Complexity Algorithm (15, 16), which was developed and validated to identify chronic and complex medical conditions and aggregate related diagnoses according to body system. 

#### Nephrology Contact
- Contact with nephrology was defined as any in-person or administrative encounter either with a provider specializing in nephrology or at a nephrology specialty care site. Patients were divided into three groups based on frequency of contact with nephrology: group 0 included patients never seen by nephrology, group 1 were patients seen by nephrology a single time, and group 2 were those followed by nephrology two or more times.

## Outcomes

#### Blood Pressure Z-Score (Systolic/Diastolic Hypertension)
- Blood pressure z-score was derived using..
- Systolic hypertension was defined as SBP ≥ 95th percentile for children under 13 years of age or SBP ≥ 130 for patients 13 years of age or older. 
- Diastolic hypertension was defined as DBP ≥ 95th percentile for children under 13 years of age or DBP ≥ 80 for patients 13 years of age or older.

#### Estimated glomerular filtration rate (eGFR)
- eGFR was derived by applying the Bedside Schwartz equation to [height](codesets/height_codeset.csv) and [serum creatinine](codesets/serum_creatinine_codeset.csv) values.

#### Advanced Chronic Kidney Disease (CKD)
- Advanced CKD was defined as an eGFR less than 60 mL/min/1.73 m2 on two separate occasions at least three months apart, with no eGFR value above 90 mL/min/1.73 m2 between or after these values, except in the event that a patient received a kidney transplant.

#### Proteinuria
- Urine Protein was defined using a series of [LOINC Codes](codesets/urine_protein_codeset.csv)
- Proteinuria was defined as having ≥2+ on urinalysis.
#### UPCR
- Urine protein to creatinine ratio (UPCR) was reported as a ratio where available in the source data or otherwise calculated using a urine protein and urine creatinine value from the same date. We further categorized UPCR into three categories: <0.5, 0.5-2.0 and >2.0.
#### Hematuria
-  Urine Blood labs were defined using a series of [LOINC Codes](codesets/urine_blood_codeset.csv)
-  Hematuria was defined as having a qualitative result indicative of the presence of red blood cells (e.g positive, too numerous to count, moderate) or a quantitative result of >5 red blood cells observed in urinalysis dipstick and/or microscopic result. 
#### Serum Albumin/Hypoalbuminemia
- Serum Albumin labs were defined using a series of [LOINC Codes](codesets/serum_albumin_codeset.csv).
- Hypoalbuminemia was defined as serum albumin value of <2.5 g/dL.
#### Transplant
- A [broad approach](codesets/kidney_transplant_broad_codeset.csv) was applied for kidney transplant that takes into consideration the following:
  - Procedure codes for kidney transplant
  - Condition codes for kidney transplant
#### Biopsy
- A broad approach was applied for kidney biopsy that takes into consideration the following: 
    - Procedure code for kidney biopsy
   - Procedure code associated with biopsy (renal not specified), accompanied by source value string search for kidney/renal/kidney biopsy ICD code or kidney finding condition code on same date
   - Condition code for kidney biopsy
   - Condition code for biopsy result (renal not specified), accompanied by kidney finding condition code on same date
   - String search for kidney/renal biopsy in visit source value
#### Dialysis
- Dialysis was defined using a list of [CPT4, ICD9Proc, and HCPCS procedure codes](codesets/dialysis_broad.csv) which were designated by the study team as either chronic or not chronic. 
#### Chronic Dialysis
- Chronic dialysis was defined as presence of at least one chronic dialysis code or two non-chronic dialysis codes separated by at least 90 days. Initiation of chronic dialysis was defined as the first chronic dialysis procedure date for those with a chronic dialysis code or the first of any dialysis procedure date for those who met the 90 day separation criteria. 
#### End Stage Kidney Disease (ESKD)
- End stage kidney disease (ESKD) was defined as having outcomes of chronic dialysis, transplant and/or death. 

