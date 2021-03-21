#### tidy functions for cleaning IWPC dataset
#### ** MATCHING ** IWPC ANALYSIS for reproduction

library(readxl)
library(tidyverse)
library(skimr)
library(tableone)
library(naniar)

'%!in%' <-Negate('%in%')


## 2.1 
iwpc <- read_excel("Q:/PharmPractice/Karnes Workgroup/WARFARIN PGX/ANALYSIS/IWPC/PS206767-553247439.xls",
                   sheet = "Subject Data")
iwpc = readxl::read_excel(path = path.expand("../PS206767-553247439.xls"),
                          sheet = 'Subject Data')
              

skim(iwpc)

iwpcdat = iwpc %>% 
  dplyr::select(id = `PharmGKB Subject ID`,
                sex = Gender,
                race_rep = `Race (Reported)`,
                race = `Race (OMB)`,
                ethnicity = `Ethnicity (OMB)`,
                eth_rep = `Ethnicity (Reported)`,
                age = Age,
                height = `Height (cm)`,
                weight = `Weight (kg)`,
                aspirin = `Aspirin`, 
                statin1 = `Simvastatin (Zocor)`,
                statin2 = `Atorvastatin (Lipitor)`, 
                statin3 = `Fluvastatin (Lescol)`,
                statin4 = `Lovastatin (Mevacor)`, 
                statin5 = `Pravastatin (Pravachol)`,
                statin6 =`Rosuvastatin (Crestor)`,
                statin7 = `Cerivastatin (Baycol)`,
                indication = `Indication for Warfarin Treatment`,
                carbamaz = `Carbamazepine (Tegretol)`,
                phenytoin = `Phenytoin (Dilantin)`,
                rifampin = `Rifampin or Rifampicin`,
                amio = `Amiodarone (Cordarone)`,
                diabetes = `Diabetes`, 
                target = `Target INR`,
                target_est = `Estimated Target INR Range Based on Indication`,
                cyp = `Cyp2C9 genotypes`,
                cyp_con = `CYP2C9 consensus`,
                vkor1639 = `VKORC1 -1639 consensus`, ## 9923231 
                vkor2255 = `VKORC1 2255 consensus`,  ## 2359612
                vkor1173 = `VKORC1 1173 consensus`, ## 9934438
                vkor1542 = `VKORC1 1542 consensus`, ## 8050894
                site = `Project Site`,
                dosewk = `Therapeutic Dose of Warfarin`,
                stable = `Subject Reached Stable Dose of Warfarin`,
                smoke = `Current Smoker`) %>% 
  na_if("NA") %>% 
  mutate(target = if_else(is.na(target),
                          target_est,
                          target)
  ) %>% 
  filter(stable == 1,
         !is.na(age),
         !is.na(dosewk),
         target %!in% c("3.5", "3.25", "3-4", "3", "2.5-3.5", "1.3", "1.75", "2")) %>% 
  mutate(
    target = if_else(is.na(target), "2.5", as.character(target)),
    target = factor(target,
                    levels = c("1.3",
                               "1.7-2.8",
                               "1.7-3.3",
                               "1.75",
                               "2",
                               "2-3",
                               "2-3.5",
                               "2.2000000000000002",
                               "2.2999999999999998", 
                               "2.5",
                               "2.5-3.5",
                               "2.6000000000000001",
                               "2.7000000000000002",
                               "2.7999999999999998",
                               "3",
                               "3-4",
                               "3.25",
                               "3.5"),
                    labels = c("1.3",
                               "2.25",
                               "2.5",
                               "1.75",
                               "2",
                               "2.5",
                               "2.75",
                               "2.2",
                               "2.3",
                               "2.5",
                               "3",
                               "2.6",
                               "2.7",
                               "2.8",
                               "3",
                               "3.5",
                               "3.25",
                               "3.5")),
    cyp = if_else(cyp_con == "*1/*11",
                  "*1/*2", 
                  if_else(cyp_con %in% 
                            c("*1/*5", 
                              "*1/*6", 
                              "*1/*13", 
                              "*1/*14"),
                          "*1/*3", 
                          cyp_con)),
    cyp = replace_na(cyp, "Missing"), 
    race = if_else(race == "Unknown",
                   race_rep,
                   race),
    race = if_else(race %in% 
                     c("Intermediate",
                       "other",
                       "NA",
                       "Other",
                       "Other Mixed Race"),
                   "Mixed or Missing",
                   if_else(grepl("Black*",
                                 race), 
                           "Black or African American",
                           race)),
    race = factor(race, 
                  levels = c("Asian",
                             "Black or African American",
                             "Mixed or Missing",
                             "White"),
                  labels = c("Asian",
                             "Black or African American",
                             "Mixed or Missing",
                             "white")) ,
    race = replace_na(race, "Mixed or Missing"),
    age = factor(age, 
                 levels = c("10 - 19", 
                            "20 - 29",
                            "30 - 39", 
                            "40 - 49",
                            "50 - 59", 
                            "60 - 69", 
                            "70 - 79", 
                            "80 - 89", 
                            "90+"),
                 labels = c("1",
                            "2", 
                            "3", 
                            "4", 
                            "5", 
                            "6",
                            "7", 
                            "8", 
                            "9")),
    age = as.numeric(age),
    vkor = vkor1639,
    statin = if_else(statin1 == "1" |
                       statin2 == "1" |
                       statin3 == "1" |
                       statin4 == "1" |
                       statin5 == "1" |
                       statin6 == "1" |
                       statin7 == "1",
                     "1",
                     "0"),
    statin = replace_na(statin, "0"),
    smoke = replace_na(smoke, "0"),
    ei = if_else(carbamaz == "1" |
                   phenytoin == "1" |
                   rifampin == "1",
                 "1",
                 "0"),
    ei = replace_na(ei, "0"),
    aspirin = replace_na(aspirin, "0"),
    diabetes = replace_na(diabetes, "0"),
    amio = replace_na(amio, "0"), 
    indication = if_else(indication %in% 
                           c("4", 
                             "4; 1", 
                             "4; 3",
                             "6; 4", 
                             "3; 2",
                             "3; 4", 
                             "3; 4; 6", 
                             "4; 6; 7; 8", 
                             "3; 4; 6; 8", 
                             "3; 4; 7", 
                             "3; 4; 7; 8",
                             "3; 4; 8", 
                             "4; 5",
                             "4; 7", 
                             "4; 7; 8",
                             "4; 8", 
                             "4;6"),
                         "MVR", 
                         if_else(indication %in% 
                                   c("1", 
                                     "1, 3",
                                     "1; 3", 
                                     "1; 6",
                                     "2",
                                     "9",
                                     "9,8",
                                     "9; 3", 
                                     "1 or 2",
                                     "1,2", 
                                     "1; 2",
                                     "1; 2; 3",
                                     "1; 2; 5; 8",
                                     "1; 2; 8",
                                     "1; 3"),
                                 "DVT/PE",
                                 if_else(indication %in% 
                                           c("3",
                                             "3; 1",
                                             "3; 1; 6", 
                                             "3; 1; 8",
                                             "3; 2",
                                             "3; 6",
                                             "3; 7",
                                             "3; 8"), 
                                         "AF",
                                         if_else(indication %in% 
                                                   c("6",
                                                     "6; 8", 
                                                     "6; 5"),
                                                 "TIA", 
                                                 "OTHER")))),
    height = as.numeric(height),
    weight = as.numeric(weight),
    sex = if_else(sex == "female", "1", "0")
    
  ) %>% 
  group_by(race, sex) %>% 
  mutate(height_imp = mean(height, na.rm = T),
         weight_imp = mean(weight, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(
    height = if_else(is.na(height), height_imp, as.numeric(height)),
    weight = if_else(is.na(weight), weight_imp, as.numeric(weight))) %>% 
  mutate_at(.vars = c("statin", 
                      "race",
                      "smoke", 
                      "sex",
                      "ei",
                      "ethnicity",
                      "aspirin",
                      "amio",
                      "diabetes",
                      "indication",
                      "site",
                      "cyp"), .funs = as.factor) %>% 
  dplyr::select(-cyp_con, -id, -race_rep,
                -statin1, -statin2, -statin3,
                -statin4, -statin5, -statin6,
                -statin7, -eth_rep, -phenytoin,
                -rifampin, -carbamaz, -target_est,
                -height_imp, -weight_imp) 

## impute with rs2359612
iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor1639) &
                             iwpcdat$race %!in% 
                             c("Black or African American", 
                               "Mixed or Missing") &
                             iwpcdat$vkor2255 == "C/C",
                           "G/G",
                           iwpcdat$vkor1639)

iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor_imp) &
                             iwpcdat$race %!in% 
                             c("Black or African American", 
                               "Mixed or Missing") &
                             iwpcdat$vkor2255 == "C/T",
                           "A/G",
                           iwpcdat$vkor_imp)

iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor_imp) &
                             iwpcdat$race %!in% 
                             c("Black or African American", 
                               "Mixed or Missing") &
                             iwpcdat$vkor2255 == "T/T",
                           "A/A",
                           iwpcdat$vkor_imp)

## impute with rs9934438
iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor_imp) &
                             iwpcdat$vkor1173 == "T/T",
                           "A/A",
                           iwpcdat$vkor_imp)

iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor_imp) &
                             iwpcdat$vkor1173 == "C/T",
                           "A/G",
                           iwpcdat$vkor_imp)

iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor_imp) &
                             iwpcdat$vkor1173 == "C/C",
                           "G/G",
                           iwpcdat$vkor_imp)

### impute with rs8050894
iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor_imp) &
                             iwpcdat$race %!in% 
                             c("Black or African American", 
                               "Mixed or Missing") &
                             iwpcdat$vkor1542 == "C/G",
                           "A/G",
                           iwpcdat$vkor_imp)

iwpcdat$vkor_imp = if_else(is.na(iwpcdat$vkor_imp) &
                             iwpcdat$race %!in% 
                             c("Black or African American",
                               "Mixed or Missing") &
                             iwpcdat$vkor1542 == "G/G",
                           "G/G",
                           iwpcdat$vkor_imp)

iwpcdat$vkor = if_else(is.na(iwpcdat$vkor_imp) &
                         iwpcdat$race %!in% 
                         c("Black or African American", 
                           "Mixed or Missing") &
                         iwpcdat$vkor1542 == "C/C",
                       "A/A",
                       iwpcdat$vkor_imp)


iwpcdat2 = iwpcdat %>% 
  mutate(vkor = replace_na(vkor, "Missing"),
         vkor = as.factor(vkor),
         race = replace_na(race, "Mixed or Missing"),
         sex = if_else(is.na(sex), "Missing",as.character(sex)),
         sex = as.factor(sex),
         sqrtdose = sqrt(dosewk),
         dosegroup = case_when(
           dosewk <= 21 ~  "low",
           dosewk >21 & dosewk < 49 ~ "intermediate",
           dosewk >= 49 ~ "high"
         ),
         dosegroup  =as.factor(dosegroup)) %>% 
  dplyr::select(-vkor1173,
                -vkor1639, 
                -vkor2255, 
                -vkor1542,
                -vkor_imp,
                -stable) %>% 
  droplevels()

### create testing training sets  ####
iwpcdat2$ID = seq(1:nrow(iwpcdat2))
train = iwpcdat2 %>% group_by(site) %>%  sample_frac(.8009) %>% ungroup()
test = anti_join(iwpcdat2, train)
test$train = '0'
train$train = "1"

iwpc_df = full_join(test, train)

iwpc_df$vkor = factor(iwpc_df$vkor, levels = c("G/G", "A/G", "A/A", "Missing"), ordered = T)
iwpc_df$race = factor(iwpc_df$race, levels = c("white", "Asian",
                                                "Black or African American",
                                                "Mixed or Missing"), ordered = T)


#### write data 
# write.csv(iwpc_df, file = "../iwpc_df.csv")

### table one ####
tableonedat = iwpc_df %>% 
  dplyr::select(dosewk, 
                vkor,
                cyp,
                age,
                height,
                weight,
                race,
                ei,
                amio,
                train)


table1 = CreateTableOne(strata = "train", data = tableonedat,
                        includeNA = F)
print(table1, nonnormal = c("dosewk", 
                            "weight",
                            "height"),
      quote =F, test = F)






