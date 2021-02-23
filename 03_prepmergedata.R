library(tidyverse)
library(readxl)
library(skimr)

#### prep for merging all hispanics + brazil data 




####  Chicago 
uic = read_xlsx("UIC.xlsx")
nrow(uic)
## race missing 
## antibiotics missing 

uicdat = uic %>% 
  dplyr::select(
    id = `PharmGKB ID`,
    race = Race,
    age = `Age (yr)`,
    sex = `Sex, M=0, F=1`,
    height = `Ht (cm)`,
    weight = `Wt (kg)`,
    target = `Goal INR`,
    dosewk = `Warfarin dose (mg/wk)`,
    diabetes = `Diabetes (Y=1, N=0)`,
    indication = `Indication for warfarin`,
    smoke = `Current Smoker (N=0, Y=1)`,
    amio = `Amiodarone N=0, Y=1`,
    aspirin = `AspirinN=0, Y=1`,
    statin = `Statin N=0, Y=1`,
    phenytoin = `Phenytoin N=0, Y=1`,
    carbamaz = `Carbamazepine N=0, Y=1`,
    cyp2 = `CYP2C9*2 3608C>T 144R>C`,
    cyp3 = `CYP2C9*3 42614A>C 359I>L`,
    vkor = `VKORC1 G3673A  -1639G>A`
  ) %>% 
  filter(is.na(id),
         dosewk <= 150) %>% 
  mutate(target = factor(target,
                         levels = c("2-3",
                                    "2.0-3.0",
                                    "2.0-3.0",
                                    "2.5-3.5"),
                         labels = c("2.5", 
                                    "2.5", 
                                    "2.5",
                                    "3")),
         indication = factor(indication,
                             levels = c("A fib",
                                        "A fib; PVD" ,
                                        "acute thrombotic microangiopathy" ,
                                        "Afib",
                                        "Afib,PVD",
                                        "cardiomyopathy",
                                        "DVT",
                                        "DVT/PE" ,
                                        "MVR" ,
                                        "MVR; A Fib",
                                        "Other (LV thrombus)",
                                        "OTHER:  brachial artery occlusion",
                                        "PE" ,
                                        "Prot C def" ,
                                        "PVD" ,
                                        "Stroke",
                                        "Stroke/TIA",
                                        "TIA; thrombosis aortic arch" ),
                             labels = c("AFIB", 
                                        "AFIB",
                                        "OTHER",
                                        "AFIB",
                                        "AFIB",
                                        "OTHER", 
                                        "DVT/PE",
                                        "DVT/PE",
                                        "MVR",
                                        "MVR",
                                        "DVT/PE",
                                        "OTHER",
                                        "DVT/PE",
                                        "OTHER",
                                        "OTHER",
                                        "TIA",
                                        "TIA",
                                        "TIA")),
         age = 
           cut(age, 
               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
               right = FALSE, 
               labels = c("10 - 19", "20 - 29", "30 - 39",
                          "40 - 49", "50 - 59", "60 - 69", 
                          "70 - 79", "80 - 89","90+")),
         age = as.numeric(age),
         cyp2 = if_else(cyp2 == "CT",
                        "*2", 
                        "*1"),
         cyp3 = if_else(cyp3 == "TT",
                        "*1", 
                        "*3"),
         ei = if_else(carbamaz == 1 | 
                      phenytoin == 1,
                      "1", "0"),
         ethnicity = "Hispanic or Latino") %>% 
  unite("cyp", c(cyp2, cyp3), sep = "/") %>% 
  mutate(cyp = factor(cyp, 
                      levels = c("*1/*1",
                                 "*1/*3",
                                 "*2/*1",
                                 "NA/NA"), 
                      labels = c("*1/*1", 
                                 "*1/*3", 
                                 "*1/*2", 
                                 "Missing")),
         vkor = if_else(is.na(vkor), "Missing", vkor),
         site = "22",
         race = "Mixed or Missing"
         ) %>% 
  dplyr::select(-carbamaz, -phenytoin, -id) %>% 
  mutate_at(.vars = c("sex", "diabetes", "smoke", "aspirin", "statin", "amio"),
            as.factor) %>% 
  mutate_if(is_character, as.factor) %>% 
  drop_na()


skim(uicdat)
nrow(uicdat
     )

#### Mt. Siani
mts = read_xlsx("MTS.xlsx", sheet = "Subject Data")

nrow(mts)
mtsdat = mts %>% 
  dplyr::select(
    id = `Subject ID`,
    race = `Race (OMB)`,
    age = `Age at Time of Consent`,
    sex = `Gender`,
    height = `Height (cm)`,
    weight = `Weight (kg)`,
    target = `Target INR`,
    dosewk = `Therapeutic Dose of Warfarin`,
    diabetes = `Diabetes`,
    indication = `Primary Indication for Warfarin Treatment`,
    smoke = `Current Smoker`,
    amio = `Amiodarone (Cordarone)`,
    aspirin = `Aspirin`,
    statin = `Statin`,
    phenytoin = `Phenytoin (Dilantin)`,
    carbamaz = `Carbamazepine (Tegretol)`,
    rifampin = `Rifampin or Rifampicin`,
    cyp = `Cyp2C9 genotypes`,
    vkor = `VKORC1 genotype:   -1639 G>A (3673); chr16:31015190(hg18); rs9923231`
  ) %>% 
  filter(dosewk <= 150) %>% 
  mutate(sex = if_else(sex == "male", "0", "1"),
         age = 
           cut(age, 
               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
               right = FALSE, 
               labels = c("10 - 19", "20 - 29", "30 - 39",
                          "40 - 49", "50 - 59", "60 - 69", 
                          "70 - 79", "80 - 89","90+")),
         age = as.numeric(age),
         ei = if_else(carbamaz == 1 | 
                        phenytoin == 1 | 
                        rifampin == 1,
                      "1", "0"),
         site = "23",
         ethnicity = "Hispanic or Latino",
         race = "Mixed or Missing",
         target = "2.5",
         indication = if_else(indication %in% 
                                c(1,2),
                              "DVT/PE", 
                              if_else(indication == 3,
                                      "AF",
                                      "OTHER")),
         cyp = factor(cyp,
                      levels = c("*1/*1",
                                 "*1/*2",
                                 "*1/*3",
                                 "*1/*8"),
                      labels = c("*1/*1", 
                                 "*1/*2",
                                 "*1/*3",
                                 "Missing")),
         diabetes = "0",
         ) %>% 
  mutate_at(.vars = c("smoke", 
                      "aspirin", 
                      "statin",
                      "amio",
                      "diabetes"),
            as.factor) %>% 
  mutate_if(is.character, as.factor) %>% 
  dplyr::select(-phenytoin, -carbamaz, -rifampin, -id)

skim(mtsdat)
nrow(mtsdat)
#### Arionza 

uaz = read_xlsx("UAZ.xlsx", sheet= "Subject Data")

uazdat = uaz %>% 
  dplyr::select(
    id = `Subject ID`,
    race = `Race (OMB)`,
    age = `Age at Time of Consent`,
    sex = `Gender`,
    height = `Height (cm)`,
    weight = `Weight (kg)`,
    target = `Target INR`,
    dosewk = `Therapeutic Dose of Warfarin`,
    diabetes = `Diabetes`,
    indication = `Primary Indication for Warfarin Treatment`,
    smoke = `Current Smoker`,
    amio = `Amiodarone (Cordarone)`,
    aspirin = `Aspirin`,
    statin = `Statin`,
    phenytoin = `Phenytoin (Dilantin)`,
    carbamaz = `Carbamazepine (Tegretol)`,
    rifampin = `Rifampin or Rifampicin`,
    cyp = `Cyp2C9 genotypes`,
    vkor = `VKORC1 genotype:   -1639 G>A (3673); chr16:31015190(hg18); rs9923231`
  ) %>% 
  filter(dosewk <= 150) %>% 
  mutate(sex = if_else(sex == "male", "0", "1"),
         age = 
           cut(age, 
               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
               right = FALSE, 
               labels = c("10 - 19", "20 - 29", "30 - 39",
                          "40 - 49", "50 - 59", "60 - 69", 
                          "70 - 79", "80 - 89","90+")),
         age = as.numeric(age),
         ei = if_else(carbamaz == 1 | 
                        phenytoin == 1 | 
                        rifampin == 1,
                      "1", "0"),
         site = "24",
         ethnicity = "Hispanic or Latino",
         indication = if_else(indication %in% c(1,2, 9), 
                              "DVT/PE", 
                              if_else(indication == 3,
                                      "AF",
                                      if_else(indication == 4,
                                              "MVR",
                                              if_else(indication == 6,
                                                      "TIA",
                                                      "OTHER"))))) %>% 
  mutate_at(.vars = c("diabetes",
                      "smoke", 
                      "aspirin", 
                      "statin",
                      "target",
                      "amio"),
            as.factor) %>% 
  mutate_if(is.character, as.factor) %>% 
  dplyr::select(-id, -phenytoin, -carbamaz, -rifampin)

skim(uazdat)




UPR <- read_excel("UPR.xlsx")


uprdat = UPR %>% 
  dplyr::select(
    id = `Subject ID`,
    race = `Race (OMB)`,
    age = `Age at Time of Consent`,
    sex = `Gender`,
    height = `Height (cm)`,
    weight = `Weight (kg)`,
    target = `Target INR`,
    dosewk = `Therapeutic Dose of Warfarin`,
    diabetes = `Diabetes`,
    indication = `Primary Indication for Warfarin Treatment`,
    smoke = `Current Smoker`,
    aspirin = `Aspirin`,
    amio = `Amiodarone (Cordarone)`,
    statin = `Statin`,
    phenytoin = `Phenytoin (Dilantin)`,
    carbamaz = `Carbamazepine (Tegretol)`,
    rifampin = `Rifampin or Rifampicin`,
    cyp = `Cyp2C9 genotypes`,
    vkor = `VKORC1 genotype:   -1639 G>A (3673); chr16:31015190(hg18); rs9923231`
  ) %>% 
  filter(dosewk <= 150) %>% 
  na_if("N/A") %>% 
  mutate(sex = if_else(sex == "male", "0", "1"),
         age = 
           cut(age, 
               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
               right = FALSE, 
               labels = c("10 - 19", "20 - 29", "30 - 39",
                          "40 - 49", "50 - 59", "60 - 69", 
                          "70 - 79", "80 - 89","90+")),
         age = as.numeric(age),
         ei = if_else(carbamaz == 1 | 
                        phenytoin == 1 | 
                        rifampin == 1,
                      "1", "0"),
         statin = replace_na(statin, "0"),
         aspirin = replace_na(aspirin, "0"),
         amio = replace_na(amio, "0"),
         ei = replace_na(ei, "0"),
         site = "25",
         ethnicity = "Hispanic or Latino",
         race = if_else(race %in% c("black", "white"),
                                   race, 
                                   "Mixed or Missing"),
         target = factor(target, 
                         levels = c("2-2.5",
                                    "2-3",
                                    "2.5-3.5"),
                         labels = c("2.25",
                                    "2.5",
                                    "3")),
         indication = if_else(indication %in% c("4", "4; 1", "4; 3", "6; 4", "3; 4"), "MVR", 
                              if_else(indication %in% c("1", "1, 3", "1; 3", "1; 6", "2", "9", "9,8", "9; 3"), "DVT/PE",
                                      if_else(indication %in% c("3", "3; 1", "3; 1; 6", "3; 1; 8", "3; 2", "3; 6", "3; 8"), "AF",
                                              if_else(indication %in% c("6", "6; 8"), "TIA", "OTHER")))),
         cyp = if_else(cyp %in% c("*1/*1", "*1/*2", 
                                  "*2/*2", "*1/*3",
                                  "*3/*3"), cyp, "Missing")) %>% 
  mutate_at(.vars = c("diabetes",
                      "smoke", 
                      "aspirin", 
                      "statin",
                      "target", 
                      "amio"),
            as.factor) %>%
  mutate_if(is.character, as.factor) %>% 
  dplyr::select(-id, -phenytoin, -carbamaz, -rifampin)


skim(uprdat)



### JOIN HISPANIC DATA TOGETHER 
latinos = full_join(uicdat, mtsdat)
skim(latinos)

latinos = full_join(latinos, uazdat)
skim(latinos)


latinos = full_join(latinos, uprdat)

skim(latinos)

#### prep for merging wth iwpc
latinosdat = 
  latinos %>% 
  filter(weight <=200) %>% 
  mutate(race = factor(race,
                       levels = c("Mixed or Missing",
                                  "American Indian or Alaska Native",
                                  "NA",
                                  "Other",
                                  "White",
                                  "black",
                                  "white"),
                       labels = c("Mixed or Missing",
                                  "Mixed or Missing",
                                  "Mixed or Missing",
                                  "Mixed or Missing",
                                  "white",
                                  "Black",
                                  "white")),
         ei = replace_na(ei, "0"),
         diabetes = if_else(diabetes == "1", "1", "0"),
         diabetes = as.factor(diabetes),
         vkor = factor(vkor,
                       levels = c("AA",
                                  "AG",
                                  "GG",
                                  "A/A",
                                  "G/A",
                                  "G/G",
                                  "NA"),
                       labels = c("AA",
                                  "AG",
                                  "GG",
                                  "AG",
                                  "AG",
                                  "GG",
                                  "Missing")))
  
  


skim(latinosdat)

#### BRAZIL 
brazil = read_excel("brazil.xlsx")

brazil_df = 
  brazil %>% 
  mutate(race = 
           if_else(White == 0, "Black or African American", "white"),
         ethnicity = "Hispanic or Latino", 
         eth = "South American",
         cyp = 
           if_else(CYP2C9_PM == 0, "*1/*1", 
                   if_else(CYP2C9_2 == 2, "*2/*2",
                           if_else(CYP2C9_2 == 1, "*1/*2",
                                   if_else(CYP2C9_3 == 2, "*3/*3",
                                           if_else(CYP2C9_3 == 1, "*1/*3", "Missing"))))),
         vkor = 
           if_else(VKORC1_G13673A == "#NULL!", 
                   "Missing",
                   if_else(VKORC1_G13673A == "2", "AA", 
                           if_else(VKORC1_G13673A == "1", "A/G", "GG"))),
         height = Height*100, 
         weight = Weight,
         sex = 
           if_else(Gender_Male == 0, 
                   "1", 
                   "0"),
         age = 
           cut(Age, 
               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
               right = FALSE, 
               labels = c("1", "2", "3", 
                          "4", "5", "6", 
                          "7", "8","9")),
         age = as.numeric(age),
         amio= 
           if_else(Amiodarone == "1",
                   "1", 
                   "0"),
         ei = 
           if_else(Carbamazepine == 1, 
                   "1",
                   "0"),
         dosewk = Dose_mg_week,
         site = "26",
         indication = if_else(indication_for_warfarin == "venous thrombosis" | indication_for_warfarin == "pulmonary thromboembolism", 
                              "DVT/PE",
                              if_else(indication_for_warfarin == "atrial fibrillation", "AFIB", 
                                      if_else(indication_for_warfarin == "cardiac valve replacement", "MVR", 
                                              if_else(indication_for_warfarin == "stroke", "TIA", "OTHER")))),
         target = "2.5", 
         smoke = if_else(smoking == "1", "1", "0"),
         diabetes = if_else(diabetic == "1", "1", "0"),
         aspirin = as.factor(aspirin),
         statin = if_else(statin == "1", "1", "0")) %>% 
  dplyr::select(site, dosewk, vkor, cyp,
                age, sex, 
                amio, ei, race, 
                ethnicity, 
                height, weight,
                indication, target, smoke,
                diabetes, statin, aspirin) %>% 
  mutate_if(is_character, 
            as_factor) %>% 
  filter(!is.na(age))



skim(brazil_df)
summary(brazil_df)

#### Second Brazil set
brazil2 = read_excel("brazil2.xls")

brazil_df2 = 
  brazil2 %>% 
  dplyr::select(
         race = `Color/Race`,
         age = `Age (years)`,
         height = `Height (m)`,
         weight = `Weight (kg)`,
         aspirin = Aspirin,
         statin = Statins,
         indication = Indication,
         ei = `Inducer status`,
         amio = Amiodarone,
         cyp = CYP2C9,
         vkor = `VKORC1 3673G>A`,
         dosewk = `Warfarin weekly  dose (mg)`,
         smoke = Smoking,
         diabetes = Diabetes,
         sex = Sex) %>% 
  mutate(race = if_else(race == "White", "white",
                   if_else(race == "Brown",
                   "Mixed or Missing", 
                   "Black")),
    ethnicity = "Hispanic or Latino",
    age = 
      cut(age, 
          breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
          right = FALSE, 
          labels = c("10 - 19", "20 - 29", "30 - 39",
                     "40 - 49", "50 - 59", "60 - 69", 
                     "70 - 79", "80 - 89","90+")),
    age = as.numeric(age),
    target = "2.5",
    aspirin = if_else(aspirin == "YES", "1", "0"),
    statin = if_else(statin == "YES", "1", "0"),
    indication = if_else(indication %in% c("Atrial Flutter", 
                                           "Atrial Fibrilation"),
                                           "AFIB",
                        if_else(indication == "Prosthetic valve",
                                "MVR",
                                if_else(indication == "Stroke",
                                        "TIA",
                                        if_else(indication %in% 
                                                  c("Throboembolism",
                                      "Thromboebolism"),
                                      "DVT/PE", 
                                      "OTHER")))),
    ei = "0",
    amio = if_else(amio == "YES", "1", '0'),
    cyp = factor(cyp,
                 levels = c("*1/*1", "*1/*2",
                            "*1/*3", "*1/*5",
                            "*2/*11", "*2/*2", 
                            "*2/*3", "*3/*3"),
                 labels = c("*1/*1", "*1/*2",
                            "*1/*3", "Missing", 
                            "Missing", "*2/*2",
                            "*2/*3", "*1/*2")),
    vkor = factor(vkor,
                  levels = c("AA","AG",
                             "GG"),
                  labels = c("AA","AG",
                             "GG")),
    site = "27",
    sex = if_else(sex == "F", "0", "1"),
    height = height*100,
    vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
    smoke = if_else(is.na(smoke) | smoke == "0", "0", "1"),
    diabetes = if_else(diabetes == "YES", "1", "0")) %>% 
  mutate_if(is.character, as.factor) %>% 
  drop_na()

skim(brazil_df2)

brazilmerged = full_join(brazil_df, brazil_df2)
skim(brazilmerged)


#### iwpc matching 
#### FIX INDICAATION, SMOKE, TARGET
iwpc = read_csv("iwpc_df.csv")

iwpc_df= iwpc %>% 
  mutate( vkor = factor(vkor, levels = c("G/G", "A/G", "A/A", "Missing"),
                        labels = c("GG", "AG", "AA", "NA")),
          target = if_else(target %in% c("1.3", "1.75", "2"), "<2",
                           if_else(target %in% c("2.2", "2.25", "2.3", "2.5", 
                                                 "2.6", "2.75", "2.8"), "2-3",
                                   if_else(target %in% c("3", "3.25", "3.5"), ">3", "Missing"))),
          race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing"),
                        labels = c("white",
                                   "Asian",
                                   "Black",
                                   "Mixed/NA"))) %>% 
  mutate_at(.vars = c("aspirin","amio","diabetes", "target", "site", "smoke","statin", "ei", "train","cyp", "dosegroup", "indication", "sex", "ethnicity"), .funs = as.factor) %>% 
  mutate_at(.vars = c("X1", "ID"),
            .funs = as.character) %>% 
  dplyr::select(-dosegroup, -train, -sqrtdose)





##### MERGE ALL AND INSPECT
data = full_join(brazilmerged, latinosdat)  %>% 
  mutate(vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
         vkor = factor(vkor, levels = c("A/G", "GG", "AA", "Missing", "AG", "NA"), labels = c("AG", "GG", "AA", "Missing", "AG", "Missing")))

data2 = full_join(data, iwpc_df) %>% 
  mutate(vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
         vkor = factor(vkor, levels = c("A/G", "GG", "AA", "Missing", "AG", "NA"), labels = c("AG", "GG", "AA", "Missing", "AG", "Missing")),
         target = factor(target, levels = c("2.5", "2-3", "2", "3", "2.25", "Missing"), labels = c("2.5", "2.5", "2", "3", "2.25", "2.5"))) %>% 
  filter(sex != "Missing") %>% 
  dplyr::select(-X1, -ID)

skim(data2)


# write.csv(data, file = "latinos.csv")
# write.csv(data2, file = "merged_iwpc_latinos_brazil.csv")

