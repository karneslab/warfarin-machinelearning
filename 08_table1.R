library(readr)
library(tableone)
library(tidyverse)

set.seed(18)

mergeddat = read_csv("~/Documents/IWPC/merged_iwpc_latinos_brazil.csv") 

mergeddat = mergeddat %>% 
  mutate_if(is_character, as_factor) %>% 
  mutate_at(.vars = c("site", "sex", "amio", "ei", "target", "smoke", "diabetes", "statin", "aspirin"), .funs = as_factor) %>% 
  mutate(
    dosegroup = case_when(
      dosewk <= 21 ~  "low",
      dosewk >21 & dosewk < 49 ~ "intermediate",
      dosewk >= 49 ~ "high"
    )) %>% 
  mutate(dosegroup = factor(dosegroup),
         target= if_else(target == "2.25"| target == "2", "2.5", as.character(target)),
         race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing", "Mixed/NA", "Black"), labels = c("white", "Asian", "Black or African American", "Mixed or Missing", "Mixed or Missing", "Black or African American")),
         vkor = factor(vkor, levels = c("GG", "AG", "AA", "Missing"))) %>% 
  mutate(target = as.factor(target),
         sqrtdose = sqrt(dosewk),
         logdose = log(dosewk),
         cohort = if_else(site %in% c("22","23","24","25","26", "27"), "ULLA","IWPC"))

train = mergeddat %>% sample_frac(.7) %>% mutate(train = 1) 
mergeddat$train = if_else(mergeddat$X1 %in% train$X1, 1, 0)

tableonedat = mergeddat %>% 
  dplyr::select(age,height,weight,
                dosewk,dosegroup,
                cyp,vkor,
                race, ethnicity,
                amio, ei, cohort)

table1 = CreateTableOne(strata = "cohort", addOverall = F,data = tableonedat,
                        includeNA = F)

tab1mat =  print(table1, nonnormal = c("dosewk", "height", "weight"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(tab1mat, file = "~/Documents/IWPC/Table1.csv")


