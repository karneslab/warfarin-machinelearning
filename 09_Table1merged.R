library(readr)
library(tableone)
library(tidyverse)

set.seed(18)


###########################################################
#####################DATA##################################
###########################################################

mergeddat = read_csv("../merged_iwpc_ULLA.csv") %>% 
  mutate_at(.vars = c( "sex", "aspirin","vkor","cyp", "race", "amio", "diabetes", "site", "smoke", "statin", "ei", "ethnicity", "indication"), .funs = as.factor) %>% 
  mutate(
    dosegroup = case_when(
      dosewk <= 21 ~  "low",
      dosewk >21 & dosewk < 49 ~ "intermediate",
      dosewk >= 49 ~ "high"
    )) %>% 
  mutate(dosegroup = factor(dosegroup),
         sqrtdose = sqrt(dosewk),
         indication = factor(indication, levels = c("AF", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA"), labels = c("AFIB", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA")),
         race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing", "Mixed/NA", "Black"), labels = c("white", "Asian", "Black or African American", "Mixed or Missing", "Mixed or Missing", "Black or African American")),
         vkor = factor(vkor, levels = c("GG", "AG", "AA", "Missing")),
         ethnicity = factor(ethnicity, levels = c("not Hispanic or Latino", "Hispanic or Latino", "Unknown"))) %>% 
  dplyr::select( -site, -country)

tableonedat = mergeddat %>% 
  dplyr::select(age, height, weight, dosewk, cyp, vkor,race, ethnicity)

table1 = CreateTableOne(addOverall = F,data = tableonedat,
                        includeNA = F)

tab1mat =  print(table1, nonnormal = c("dosewk", "height", "weight"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(tab1mat, file = "~/Documents/IWPC/TableS2.csv")

