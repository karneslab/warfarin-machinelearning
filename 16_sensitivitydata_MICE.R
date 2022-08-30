library(tidyverse)
library(readxl)
library(skimr)
library(naniar)
library(mice)
#### prep for merging all hispanics + brazil data 




####  Chicago 
####################
uic = read_xlsx("../UIC.xlsx")

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
         dosewk <= 175,
         dosewk >=7) 

# vis_miss(uicdat)

uicdat = uicdat  %>% 
  mutate(country = "US",
         target = factor(target,
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
  group_by(sex,race) %>% 
  mutate(height = if_else(is.na(height), mean(height, na.rm = T), height)) %>% 
  ungroup() %>% 
  dplyr::select(-carbamaz, -phenytoin, -id) %>% 
  mutate_at(.vars = c("sex", "diabetes", "smoke", "aspirin", "statin", "amio"),
            as.factor) %>% 
  mutate_if(is_character, as.factor) 


# skim(uicdat)
####################

#### Mt. Siani
####################
mts = read_xlsx("../MTS.xlsx", sheet = "Subject Data")

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
  na_if("NA")%>% 
  filter(dosewk <= 175,
         dosewk >=7) 

# vis_miss(mtsdat)

mtsdat = mtsdat %>% 
  mutate(country = "US",
         sex = if_else(sex == "male", "0", "1"),
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
                                 "*1/*5",
                                 "*1/*8"),
                      labels = c("*1/*1", 
                                 "*1/*2",
                                 "*1/*3",
                                 "Missing", 
                                 "Missing"))
         #diabetes = sample(c(0,1),nrow(mtsdat), prob = c(.5,.5), replace = T),
         # diabetes = if_else(is.na(diabetes), "0", as.character(diabetes))
  ) %>% 
  mutate_at(.vars = c("smoke", 
                      "aspirin", 
                      "statin",
                      "amio",
                      "diabetes"),
            as.factor) %>% 
  mutate_if(is.character, as.factor) %>% 
  dplyr::select(-phenytoin, -carbamaz, -rifampin, -id)

# skim(mtsdat)

####################


#### Arizona 
####################
uaz = read_xlsx("../UAZ.xlsx", sheet= "Subject Data")


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
  filter(dosewk <= 175,
         dosewk >=7,
         weight <= 200) %>% 
  na_if("NA")

# vis_miss(uazdat)

uazdat = uazdat %>% 
  mutate(country = "US",
         sex = if_else(sex == "male", "0", "1"),
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
         race = if_else(is.na(race), "Mixed or Missing", race),
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
                      "amio",
                      "race",
                      "sex",
                      "cyp",
                      "vkor",
                      "ethnicity",
                      "ei",
                      "site",
                      "indication"),
            as.factor) %>% 
  dplyr::select(-id, -phenytoin, -carbamaz, -rifampin)

# skim(uazdat)
####################

## UPR 
####################
UPR <- read_excel("../UPR.xlsx")


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
  filter(dosewk <= 175,
         dosewk >=7) %>% 
  na_if("N/A") 

# vis_miss(uprdat)


uprdat = uprdat %>% 
  mutate(country = "PR",
         sex = if_else(sex == "male", "0", "1"),
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
        #  statin = if_else(is.na(statin), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(statin)),
         # statin = if_else(is.na(statin), "0", as.character(statin)),
       #   aspirin = if_else(is.na(aspirin), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(aspirin)),
         # aspirin = if_else(is.na(aspirin), "0", as.character(aspirin)),
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


# skim(uprdat)

####################

### JOIN HISPANIC DATA TOGETHER 
####################
latinos = full_join(uicdat, mtsdat)
skim(latinos)

latinos = full_join(latinos, uazdat)
# skim(latinos)


latinos = full_join(latinos, uprdat)

# skim(latinos)

#### prep for merging wth iwpc
latinosdat = 
  latinos %>% 
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
         diabetes = as.factor(diabetes),
         vkor = factor(vkor,
                       levels = c("AA",
                                  "AG",
                                  "GG",
                                  "Missing",
                                  "A/A",
                                  "G/A",
                                  "G/G",
                                  "NA"),
                       labels = c("AA",
                                  "AG",
                                  "GG",
                                  "Missing",
                                  "AA",
                                  "AG",
                                  "GG",
                                  "Missing")))




# skim(latinosdat)
####################

#### University of Sao Paulo Brazil 
####################
brazil = read_excel("../brazil.xlsx")

brazil_df = 
  brazil %>% 
  na_if("#NULL!")

# vis_miss(brazil_df)

brazil_df = brazil_df %>% 
  mutate(country = "Brazil",
         race = 
           if_else(White == 0, "Black or African American", "white"),
         ethnicity = "Hispanic or Latino", 
         eth = "South American",
         cyp = 
           if_else(is.na(CYP2C9_PM), "Missing", 
                   if_else(CYP2C9_PM == 0, "*1/*1", 
                           if_else(CYP2C9_2 == 2, "*2/*2",
                                   if_else(CYP2C9_2 == 1, "*1/*2",
                                           if_else(CYP2C9_3 == 2, "*3/*3",
                                                   if_else(CYP2C9_3 == 1, "*1/*3", "Missing")))))),
         vkor = 
           if_else(is.na(VKORC1_G13673A), 
                   "Missing",
                   if_else(VKORC1_G13673A == "2", "AA", 
                           if_else(VKORC1_G13673A == "1", "AG", "GG"))),
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
         amio= if_else(is.na(amio), "0", "1"),
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
        #  smoke = if_else(is.na(smoking), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(smoking)),
         # smoke = if_else(is.na(smoking), "0", as.character(smoking)),
        #  diabetes = if_else(is.na(diabetic), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(diabetic)),
         # diabetes = if_else(is.na(diabetic), "0", as.character(diabetic)),
         diabetes = as.factor(diabetic),
         smoke = as.factor(smoking),
         aspirin = as.factor(aspirin),
        # statin = if_else(is.na(statin), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(statin)),
         # statin = if_else(is.na(statin), "0", as.character(statin)),
         statin = as.factor(statin)) %>% 
  dplyr::select(site, dosewk, vkor, cyp,
                age, sex, 
                amio, ei, race, 
                ethnicity, 
                height, weight,
                indication, target, smoke,
                diabetes, statin, aspirin,
                country) %>% 
  mutate_if(is_character, 
            as_factor) %>% 
  filter(!is.na(age), dosewk <= 175,
         dosewk >= 7)



# vis_miss(brazil_df)
####################

#### Instituto Nacional de Cardiologia Laranjeiras
####################
brazil2 = read_excel("../brazil2.xls")

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
  filter(!is.na(height),
         !is.na(weight),
         dosewk<=175,
         dosewk>=7)

# vis_miss(brazil_df2)

brazil_df2 = brazil_df2 %>% 
  mutate(country = "Brazil" ,
         race = if_else(race == "White", "white",
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
                       levels = c("AA","GA",
                                  "GG"),
                       labels = c("AA","AG",
                                  "GG")),
         site = "27",
         sex = if_else(sex == "F", "0", "1"),
         height = height*100,
         vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
         smoke = factor(smoke, levels = c("YES", "no", "No", "Ex"), labels = c("1","0","0", "0")),
         # smoke = if_else(is.na(smoke), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(smoke)),
         # smoke = if_else(is.na(smoke), "0", as.character(smoke)),
         smoke = as.factor(smoke),
         diabetes = if_else(diabetes == "YES", "1", "0")) %>% 
  mutate_if(is.character, as.factor) 

# skim(brazil_df2)
####################

### Porto Alegre
####################
brazil3  = read_excel("../brazil3.xlsx")

brazil_df3 = 
  brazil3 %>% 
  dplyr::select(
    race = `Race`,
    age = `Age`,
    height = `Height (m)`,
    weight = `Weight (kg)`,
    carbamaz = `Carbamazepine`,
    phenytoin = `Phenytoin`, 
    aspirin = Aspirin,
    statin = Statin,
    indication = `Anticoagulation indication`,
    amio = Amiodarone,
    cyp2 = `CYP2C9*2`,
    cyp3 = `CYP2C9*3`,
    vkor = `VKORC1 -1639`,
    dosewk = `Weekly dose (mg)`,
    smoke = Smoking,
    diabetes = Diabetes,
    sex = Gender) 

# vis_miss(brazil_df3)

brazil_df3 = brazil_df3 %>% 
  mutate(country = "Brazil", 
         race = if_else(race == "white", "white",
                        if_else(race == "black","Black or African American", 
                                "Mixed or Missing")),
         race = if_else(is.na(race), "Mixed or Missing", race),
         ethnicity = "Hispanic or Latino",
         ei = if_else(carbamaz == 1 | 
                        phenytoin == 1,
                      "1", "0"),
         age = 
           cut(age, 
               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
               right = FALSE, 
               labels = c("10 - 19", "20 - 29", "30 - 39",
                          "40 - 49", "50 - 59", "60 - 69", 
                          "70 - 79", "80 - 89","90+")),
         age = as.numeric(age),
         indication = if_else(str_detect(indication, "pros"), "MVR",
                              if_else(indication == "Atrial flutter" | str_detect(indication, "AF"), "AFIB",
                                      if_else(str_detect(indication, "stroke") | str_detect(indication, "TIA"), "TIA",
                                              if_else(str_detect(indication, "DVT") | str_detect(indication, "thrombus") | str_detect(indication, "PE"),
                                                      "DVT/PE", 
                                                      "OTHER")))),
         cyp = if_else(cyp3 == "AA" & cyp2 == "CC", "*1/*1", 
                       if_else(cyp3 == "AA" & cyp2 == "CT", "*1/*2",
                               if_else(cyp3 == "AA" & cyp2 == "TT", "*2/*2", 
                                       if_else(cyp3 == "AC" & cyp2 == "CC", "*1/*3", 
                                               if_else(cyp3 == "AC" & cyp2 == "CT", "*2/*3","*3/*3" ))))),
         vkor = factor(vkor,
                       levels = c("AA","GA",
                                  "GG"),
                       labels = c("AA","AG",
                                  "GG")),
         site = "28",
         sex = if_else(sex == "f", "0", "1"),
         height = as.numeric(height)*100,
         # smoke = if_else(is.na(smoke) | smoke == "n" | smoke == "ex", "0", "1"),
         smoke = if_else(smoke == "n",0,1) , 
         smoke = as.factor(smoke),
         diabetes = as.factor(diabetes),
         #smoke = if_else(is.na(smoke), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(smoke)),
         #aspirin = if_else(is.na(aspirin), sample(0:1, n(), replace =T, prob = c(0.5,0.5)), as.integer(aspirin)),
         # aspirin = if_else(is.na(aspirin), "0", as.character(aspirin)),
        #  statin = if_else(is.na(statin), sample(0:1, n(), replace =T, prob = c(0.5,0.5)), as.integer(statin)),
         # statin = if_else(is.na(statin), "0", as.character(statin)), 
        statin = as.factor(statin),
        aspirin = as.factor(aspirin),
         indication = if_else(is.na(indication),"OTHER", indication),
         amio = if_else(is.na(amio), "0", as.character(amio)),
         ei = if_else(is.na(ei), "0", as.character(ei)),
         # diabetes = if_else(is.na(diabetes), sample(0:1, n(), replace =T, prob = c(0.5,0.5)), as.integer(diabetes)),
         # diabetes = if_else(is.na(diabetes), "0", as.character(diabetes)),
         weight = as.numeric(weight), 
         dosewk = as.numeric(dosewk))%>% 
  mutate_if(is.character, as.factor) %>% 
  mutate_if(is.integer, as.factor) %>% 
  group_by(race, sex) %>% 
  mutate(height_imp = mean(height, na.rm = T),
         weight_imp = mean(weight, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(height = if_else(is.na(height), height_imp, as.numeric(height)),
         weight = if_else(is.na(weight), weight_imp, as.numeric(weight))) 

# vis_miss(brazil_df3)

brazilmerged = full_join(brazil_df, brazil_df2)
brazilmerged2 = full_join(brazilmerged, brazil_df3) %>% 
  dplyr::select(-height_imp, -weight_imp, -carbamaz, -phenytoin, -cyp2,-cyp3, -target) %>% 
  filter(dosewk <= 175,
         dosewk >=7)

# vis_miss(brazilmerged2)
####################

#### COLOMBIA DATA
####################
CO <- read_excel("../COL.xlsx") %>% 
  na_if("MD") %>% 
  filter(Dose >7,
         Dose <= 175,
         Stable == "1") %>% 
  mutate(Code = toupper(Code),
         Code = gsub(" ", "", Code))

CO1b =  read_excel(path = "../COL1b.xlsx",
                   skip = 2)  %>% 
  mutate(`...1` = toupper(`...1`),
         Code = gsub(" ", "", `...1`),
         diabetes = `MELLITUS DIABETES...2`)%>% 
  dplyr::select("Code", "diabetes") %>% 
  inner_join(CO)

# vis_miss(CO1b)

codat = CO1b %>% 
  mutate(
    race = "Mixed or Missing", 
    ethnicity = "Hispanic or Latino",
    country = "Colombia", 
    site = "29", 
    age = cut(Age, 
              breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
              right = FALSE, 
              labels = c("10 - 19", "20 - 29", "30 - 39",
                         "40 - 49", "50 - 59", "60 - 69", 
                         "70 - 79", "80 - 89","90+")),
    dosewk = Dose,
    age = as.numeric(age),
    height = if_else(is.na(Height), mean(as.numeric(Height), na.rm= T), as.numeric(Height)),
    height = height*100,
    weight = Weight,
    vkor = `VKORC1 rs9923231`,
    cyp = if_else(`CYP2C9 * 2 rs1799853` == "2 * / 2 *" & `CYP2C9 * 3 rs1057910` == "1 * / 1 *", "*2/*2",
                  if_else(`CYP2C9 * 2 rs1799853` == "1 * / 2 *" & `CYP2C9 * 3 rs1057910` == "1 * / 1 *", "*1/*2",
                          if_else(`CYP2C9 * 2 rs1799853` == "1 * / 2 *" & `CYP2C9 * 3 rs1057910` == "1 * / 3 *","*2/*3",
                                  if_else(`CYP2C9 * 2 rs1799853` == "1 * / 1 *" & `CYP2C9 * 3 rs1057910` == "1 * / 3 *", "*1/*3", "*1/*1")))),
    vkor = if_else(`VKORC1 rs9923231` == "A / A", "AA", 
                   if_else(`VKORC1 rs9923231` == "G / A", "GA", "GG")),
    ei= as_factor(Inducers),
    amio = as_factor(Amiodarone),
    indication = if_else(Indication == "Valvular replacement", "MVR",
                         if_else(Indication== "Atrial fibrillation", "AFIB", if_else(Indication == "Stroke", "TIA", if_else(Indication == "DVT" | Indication == "PE", "DVT/PE", "OTHER")))),
    statin = as_factor(Lovastatin),
    aspirin = if_else(grepl("aspirin", ignore.case = T, `Antiplatelet drugs`), "1", "0"),
    # smoke = "0",
   # smoke = sample(0:1, n(), replace =T, prob = c(0.5,0.5)),
    smoke = NA,
    diabetes = if_else(diabetes == 1, "1", diabetes),
 #  diabetes = as.integer(diabetes),
    # diabetes = if_else(is.na(diabetes),
    #                                  sample(0:1, n(),
    #                                         replace = T,
    #                                         prob = c(.5,.5)),
    #                                  as.integer(diabetes)),
    # diabetes = as.factor(diabetes),
    target = "2.5", 
    sex = if_else(Gender == "2", "0", "1")
  ) %>% 
  dplyr::select(site, dosewk, age, height, weight, cyp, vkor, race, ethnicity, country, amio, aspirin, statin, smoke, diabetes, ei, indication, sex, target)

# vis_miss(codat)



####################

### COL military 
####################
#################
# CO2 <- read_excel("../COL2.xlsx") %>% 
#   filter(WARMGSEM1 >7,
#          WARMGSEM1 <= 175) 
# 
# skim(CO2)
# 
# codat2 = CO2 %>% 
#   mutate(
#     race = "Mixed or Missing", 
#     ethnicity = "Hispanic or Latino",
#     country = "Colombia", 
#     site = "30", 
#     age = cut(EDAD_ANOS, 
#               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
#               right = FALSE, 
#               labels = c("10 - 19", "20 - 29", "30 - 39",
#                          "40 - 49", "50 - 59", "60 - 69", 
#                          "70 - 79", "80 - 89","90+")),
#     dosewk = WARMGSEM1,
#     age = as.numeric(age),
#     height = if_else(is.na(TALLA_CM), mean(as.numeric(TALLA_CM), na.rm= T), as.numeric(TALLA_CM)),
#     # height = height*100,
#     weight = PESO_KG,
#     vkor = VKORC1rs9923231,
#     cyp = if_else(CYP2C9rs1799853 == "22" & `CYP2CPrs1057910` == "11", "*2/*2",
#                   if_else(CYP2C9rs1799853 == "12*" & CYP2CPrs1057910 == "11", "*1/*2",
#                           if_else(CYP2C9rs1799853 == "12" & CYP2CPrs1057910 == "13*","*2/*3",
#                                   if_else(CYP2C9rs1799853 == "11" & CYP2CPrs1057910 == "13", "*1/*3", "*1/*1")))),
#     cyp = if_else(is.na(cyp), "Missing", cyp),
#     ei= if_else(FENITOINA == 1 | CARBAMACEPINA == 1 | RIFAMPICINA == 1, 1, 0),
#     ei = as_factor(ei),
#     amio = if_else(is.na(AMIODARONA), 0, if_else(AMIODARONA == 0, 0, 1)),
#     amio = as_factor(amio),
#     indication = if_else(VALVULOPATIACARDIACA == 1, "MVR",
#                          if_else(FIBRIAURICULAR == 1, "AFIB", 
#                                  if_else(ACVTROMBOTICO == 1, "TIA", 
#                                          if_else(TROMBOSISVENOSA == 1 | EMBPULMONAR==1, "DVT/PE", "OTHER")))),
#     statin = as_factor(ESTATINAS),
#     aspirin = as_factor(ASPIRINA),
#     # smoke = if_else(is.na(TABAQUISMO), 0,
#     #                 if_else(TABAQUISMO == 0, 0, 1)),
#     smoke = as.integer(TABAQUISMO),
#     smoke = if_else(is.na(smoke), sample(0:1, n(), replace =T, prob = c(0.5,0.5)), (smoke)),
#     smoke = as_factor(smoke),
#     diabetes = as.factor(DIABETESMELLITUS),
#     target = "2.5", 
#     sex = if_else(SEXO == "1", "0", "1"),
#     metformin = 0,
#     metformin= as_factor(metformin)
#   ) %>% 
#   dplyr::select(site, dosewk, age, height, weight, cyp, vkor, race, ethnicity, country, amio, aspirin, statin, smoke, diabetes, ei, indication, sex, target, metformin) %>% 
#   mutate_if(
#     is.character, as.factor
#   )

###########################

### MERGE SOUTH AMERICA
#############################
SA  = full_join(codat, brazilmerged2)


# vis_miss(SA)
# skim(SA)
#SA2 = full_join(SA, codat2)
  
#############################


#### iwpc matching 
####################
#### FIX INDICAATION, SMOKE, TARGET
'%!in%' <-Negate('%in%')


## 2.1 
## 2.1 
iwpc <- read_xls("../PS206767-553247439.xls", 
                 sheet = "Subject Data")




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
         # target %!in% c("3.5", "3.25", "3-4", "3", "2.5-3.5", "1.3", "1.75", "2")
  ) 

iwpcdat = iwpcdat %>% 
  mutate(
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
   # statin = if_else(is.na(statin), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(statin)),
    # statin = if_else(is.na(statin), "0", as.character(statin)),
    # smoke = if_else(is.na(smoke), "0", as.character(smoke)), 
   #  smoke = if_else(is.na(smoke), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(smoke)),
    ei = if_else(carbamaz == "1" |
                   phenytoin == "1" |
                   rifampin == "1",
                 "1",
                 "0"),
    ei = replace_na(ei, "0"),
     #aspirin = if_else(is.na(aspirin), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(aspirin)),
    # aspirin = if_else(is.na(aspirin),"0", as.character(aspirin)),
    # diabetes = if_else(is.na(diabetes), "0", as.character(diabetes)),
   # diabetes = if_else(is.na(diabetes), sample(0:1, n(), replace = T, prob = c(.5,.5)), as.integer(diabetes)),
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
                -height_imp, -weight_imp, -target) 

# vis_miss(iwpcdat)

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


iwpc_df = iwpcdat %>% 
  mutate(vkor = replace_na(vkor, "Missing"),
         vkor = as.factor(vkor),
         race = replace_na(race, "Mixed or Missing"),
         sex = if_else(is.na(sex), "Missing",as.character(sex)),
         sex = as.factor(sex),
         sqrtdose = sqrt(as.numeric(dosewk)),
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
iwpc_df$ID = seq(1:nrow(iwpc_df))
train = iwpc_df %>% group_by(site) %>%  sample_frac(.8009) %>% ungroup()
test = anti_join(iwpc_df, train)
test$train = '0'
train$train = "1"

iwpc_df = full_join(test, train)

iwpc_df$vkor = factor(iwpc_df$vkor, levels = c("G/G", "A/G", "A/A", "Missing"), ordered = T)
iwpc_df$race = factor(iwpc_df$race, levels = c("white", "Asian",
                                               "Black or African American",
                                               "Mixed or Missing"), ordered = F)



iwpc_df= iwpc_df %>% 
  mutate( vkor = factor(vkor, levels = c("G/G", "A/G", "A/A", "Missing"),
                        labels = c("GG", "AG", "AA", "NA")),
          # target = if_else(target %in% c("1.3", "1.75", "2"), "<2",
          #                  if_else(target %in% c("2.2", "2.25", "2.3", "2.5", 
          #                                        "2.6", "2.75", "2.8"), "2-3",
          #                          if_else(target %in% c("3", "3.25", "3.5"), ">3", "2.5"))),
          race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing"),
                        labels = c("white",
                                   "Asian",
                                   "Black",
                                   "Mixed/NA"),ordered = F),
          dosewk = as.numeric(dosewk)) %>% 
  mutate_at(.vars = c("aspirin","amio","diabetes", "site", "smoke","statin", "ei", "train","cyp", "dosegroup", "indication", "sex", "ethnicity"), .funs = as.factor) %>% 
  mutate_at(.vars = c( "ID"),
            .funs = as.character) %>% 
  dplyr::select(-dosegroup, -train, -sqrtdose) %>% 
  filter(
    dosewk <=175, dosewk >=7)


# vis_miss(iwpc_df)
###############################

##### MERGE ALL AND INSPECT
###################
data = full_join(SA, latinosdat)  %>% 
  mutate(vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
         dosewk = as.numeric(dosewk)) %>% 
  filter(site != "27",
         sex != "Missing",
         weight >=35,
         weight <= 150,
         height <= 200,
         height >= 130) %>% 
  select(-target)

data2 = full_join(data, iwpc_df) %>% 
  mutate(vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
         vkor = factor(vkor, levels = c("A/G", "GG","GA", "AA", "Missing", "AG", "NA"), labels = c("AG", "GG","AG", "AA", "Missing", "AG", "Missing")),
         ethnicity= factor(ethnicity, levels =c("Hispanic or Latino",   "not Hispanic or Latino", "Not Hispanic or Latino", "Unknown"), labels = c("Hispanic or Latino" , "not Hispanic or Latino", "not Hispanic or Latino", "Unknown")),
         race = factor(race, 
                       levels = c("white" , "Black or African American", "Black" , "Mixed or Missing","Asian", "Mixed/NA" ),
                       labels =c("white" , "Black or African American", "Black or African American" , "Mixed or Missing","Asian", "Mixed or Missing" )),
         smoke = factor(smoke, levels = c("0","1","2"), labels = c("0","1","1"))) %>% 
  filter(sex != "Missing",
         weight >=35,
         weight <= 150,
         height <= 200,
         height >= 130) %>% 
  dplyr::select( -ID) %>% 
  droplevels() %>% 
  mutate_if(is.character, as.factor)

vis_miss(data2)
visdat::vis_miss(data)

###################
# IMPUTATION # 
# !! ##
# MICE ### 
imp = mice(data2)
imp_data = complete(imp)
# Write data 

# write.csv(imp_data, file = "../merged_iwpc_ULLAsens_imp.csv")

