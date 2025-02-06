# BIOSTAT696

**Data Dictionary for Recoded Variables**

This document describes the recoded variables in the cleaned dataset **birth_clean.RDS**.

---

### **Household and Demographic Variables**

| Variable Name         | Description |
|----------------------|-------------|
| clusterid           | Cluster ID from DHS data |
| household_number    | Household number within the cluster |
| household_rank      | Rank of individual within household |
| household_members   | Number of household members |
| household_under_5   | Number of children under 5 in household |
| household_head      | Gender of household head ("male", "female") |
| type_of_residence   | Urban or rural residence ("urban", "rural") |
| wealth_index        | Wealth index category |
| survey_weight       | Survey weighting variable |
| strata              | Stratification variable for survey design |
| primary_sampling_unit | PSU identifier for survey design |

---

### **Maternal Variables**

| Variable Name        | Description |
|---------------------|-------------|
| mother_current_age  | Age of the mother |
| mother_education    | Education level of the mother ("No education", "Primary", "Secondary or Higher") |
| mother_occupation   | Occupation type ("Manual labor", "Not working or non manual") |
| mother_height       | Mother's height (in cm) |
| mother_bmi         | Mother's BMI |
| mother_marital_status | Marital status of the mother ("married/partnered", "no partner") |
| antenatal_visits    | Number of antenatal visits attended |
| mother_hbp         | Whether mother was ever diagnosed with high blood pressure (Yes/No) |
| mother_diabetes    | Whether mother was ever diagnosed with diabetes (Yes/No) |

---

### **Child Variables**

| Variable Name       | Description |
|--------------------|-------------|
| birth_order       | Birth order of the child |
| sex_of_child      | Sex of child ("male", "female") |
| birth_weight_type | Source of birth weight measurement ("mother recall", "written card") |
| birth_weight      | Birth weight in grams (NA if missing) |
| birth_weight_cat  | Categorized birth weight ("low" < 2500g, "normal" 2500-4000g, "high" > 4000g) |
| birth_date        | Date of birth (YYYY-MM-DD) |

---

### **Household Infrastructure Variables**

| Variable Name           | Description |
|------------------------|-------------|
| household_water       | Type of drinking water source (binary: 1 = "Improved", 0 = "Unimproved") |
| household_toilet      | Type of toilet facility ("Flush toilet", "Pit latrine", "No facility", "Other") |
| household_electricity | Availability of electricity ("Yes", "No") |
| household_floor       | Floor material (binary: 1 = "Finished", 0 = "Unfinished") |
| household_cooking_fuel | Type of cooking fuel (binary: 1 = "Clean", 0 = "Non clean") |

---


