library(tidyverse)
library(readxl)
library(skimr)
library(naniar)
#### prep for merging all hispanics + brazil data 




####  Chicago 
uic = read_xlsx("UIC.xlsx")

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
  )

vis_miss(uicdat)

uicdat = uicdat  %>% 
  mutate(
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
         race = "Mixed or Missing") %>% 
  dplyr::select(-carbamaz, -phenytoin) %>% 
  mutate_at(.vars = c("sex", "diabetes", "smoke", "aspirin", "statin", "amio"),
            as.factor) %>% 
  filter(!is.na(height),is.na(id),
                  dosewk <= 150 )%>% 
  select(dosewk, age, height, weight, race, amio, ei, cyp,vkor)  %>%  mutate_if(is_character, as.factor) 


visdat::vis_miss(uicdat)



#### Mt. Siani
mts = read_xlsx("MTS.xlsx", sheet = "Subject Data")

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
  filter(dosewk <= 150) 

vis_miss(mtsdat)

mtsdat = mtsdat %>% 
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
                                 "Missing"))) %>% 
  mutate_at(.vars = c("smoke", 
                      "aspirin", 
                      "statin",
                      "amio",
                      "diabetes"),
            as.factor) %>% 
  mutate_if(is.character, as.factor) %>% 
  select(dosewk, age, height, weight, race, amio, ei, cyp,vkor)

vis_miss(mtsdat)



#### Arionza 

uaz = read_xlsx("UAZ.xlsx", sheet= "Subject Data")
sum(is.na(uaz))


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
  na_if("NA")


vis_miss(uazdat)


uazdat = uazdat %>% 
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
         race = if_else(is.na(race), "Mixed or Missing", race),
         ethnicity = "Hispanic or Latino") %>% 
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
  select(dosewk, age, height, weight, race, amio, ei,cyp,vkor)

vis_miss(uazdat)




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
  na_if("N/A") 

vis_miss(uprdat)


uprdat = uprdat %>% 
  mutate(sex = if_else(sex == "male", "0", "1"),
         age = 
           cut(age, 
               breaks=c(10, 20, 30, 40, 50, 60, 70,80, 90, 100), 
               right = FALSE, 
               labels = c("10 - 19", "20 - 29", "30 - 39",
                          "40 - 49", "50 - 59", "60 - 69", 
                          "70 - 79", "80 - 89","90+")),
         age = as.numeric(age),
         site = "25",
         ethnicity = "Hispanic or Latino",
         race = if_else(race %in% c("black", "white"),
                        race, 
                        "Mixed or Missing"),
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
  filter(!is.na(phenytoin),!is.na(carbamaz),!is.na(rifampin)) %>% 
  mutate(ei = if_else(carbamaz == 1 | 
                                 phenytoin == 1 | 
                                 rifampin == 1,
                               "1", "0")) %>% 
  select(dosewk, age, height, weight, race, amio, ei,cyp,vkor)


vis_miss(uprdat)



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
                                  "AG",
                                  "AG",
                                  "GG",
                                  "Missing")))




skim(latinosdat)

#### BRAZIL 
brazil = read_excel("brazil.xlsx")

brazil_df = 
  brazil %>% 
  na_if("#NULL!")

vis_miss(brazil_df)

brazil_df = brazil_df %>% 
  mutate(race = 
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
         target = "2.5") %>% 
  dplyr::select(site, dosewk, vkor, cyp,
                age, sex, 
                amio, ei, race, 
                ethnicity, 
                height, weight, target, statin, aspirin) %>% 
  mutate_if(is_character, 
            as_factor) %>% 
  filter(!is.na(age),!is.na(amio))%>% 
  select(dosewk, age, height, weight, race, amio, ei, cyp,vkor)

vis_miss(brazil_df)


#### iwpc matching 
#### FIX INDICAATION, SMOKE, TARGET
'%!in%' <-Negate('%in%')


## 2.1 
iwpc <- read_xls("PS206767-553247439.xls",
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
         target %!in% c("3.5", "3.25", "3-4", "3", "2.5-3.5", "1.3", "1.75", "2")) 

vis_miss(iwpcdat)

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
    ei = if_else(carbamaz == "1" |
                   phenytoin == "1" |
                   rifampin == "1",
                 "1",
                 "0"),
    height = as.numeric(height),
    weight = as.numeric(weight),
    sex = if_else(sex == "female", "1", "0")
  ) %>% 
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
                -rifampin, -carbamaz, -target_est) %>% 
  filter(!is.na(age),!is.na(height),!is.na(weight),!is.na(amio),
         !is.na(ei)) 

vis_miss(iwpcdat)


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
  droplevels() %>% 
  select(dosewk, age, height, weight, race, amio, ei,cyp,vkor)

### create testing training sets  ####
iwpc_df$ID = seq(1:nrow(iwpc_df))
train = iwpc_df %>%  sample_frac(.8009) %>% ungroup()
test = anti_join(iwpc_df, train)
test$train = '0'
train$train = "1"

iwpc_df = full_join(test, train)

iwpc_df$vkor = factor(iwpc_df$vkor, levels = c("G/G", "A/G", "A/A", "Missing"), ordered = T)
iwpc_df$race = factor(iwpc_df$race, levels = c("white", "Asian",
                                               "Black or African American",
                                               "Mixed or Missing"), ordered = T)



iwpc_df= iwpc_df %>% 
  mutate( vkor = factor(vkor, levels = c("G/G", "A/G", "A/A", "Missing"),
                        labels = c("GG", "AG", "AA", "NA")),
          
          race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing"),
                        labels = c("white",
                                   "Asian",
                                   "Black",
                                   "Mixed/NA"),ordered = F),
          dosewk = as.numeric(dosewk)) %>% 
  mutate_at(.vars = c("amio", "ei", "train","cyp"), .funs = as.factor) %>% 
  mutate_at(.vars = c( "ID"),
            .funs = as.character) %>% 
  select(dosewk, age, height, weight, race, amio, ei,cyp,vkor)





##### MERGE ALL AND INSPECT
data = full_join(brazil_df, latinosdat)  %>% 
  mutate(vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
         dosewk = as.numeric(dosewk))

data2 = full_join(data, iwpc_df) %>% 
  mutate(vkor = if_else(is.na(vkor), "Missing", as.character(vkor)),
         vkor = factor(vkor, levels = c("A/G", "GG", "AA", "Missing", "AG", "NA"), labels = c("AG", "GG", "AA", "Missing", "AG", "Missing")),
         race = factor(race, 
                       levels = c("white" , "Black or African American", "Black" , "Mixed or Missing","Asian", "Mixed/NA" ),
                       labels =c("white" , "Black or African American", "Black or African American" , "Mixed or Missing","Asian", "Mixed or Missing" ))) %>% 
  droplevels()

skim(data2)
vis_miss(data2)

# write.csv(data, file = "latinos_completecase.csv")
 write.csv(data2, file = "merged_iwpc_latinos_brazil_completecase.csv")

 tableone::CreateTableOne(data = data2)
