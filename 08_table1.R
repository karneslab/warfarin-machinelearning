library(readr)
library(tableone)
library(tidyverse)

set.seed(18)

## load and prep data
latinos = read_csv("../ULLA.csv") %>% 
  mutate(sqrtdose = sqrt(dosewk),
         race = factor(race, levels = c("white", "Black or African American", "Mixed or Missing", "Black"), labels = c("white",  "Black or African American", "Mixed or Missing", "Black or African American")),
         vkor = factor(vkor, levels = c("GG","GA", "AG", "AA", "Missing"),
                       labels = c("GG","AG", "AG", "AA" , "Missing")),
         indication = factor(indication, levels = c("AF", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA"), labels = c("AFIB", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA")),
         race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing"))) %>% 
  mutate_at(.vars = c("site", "sex", "amio", "ei", "smoke", "diabetes", "statin",
                      "aspirin", "cyp", "indication", "country", "ethnicity"), .funs = factor) %>% 
  dplyr::select(age,height,weight,
                dosewk,
                cyp,vkor,
                race, ethnicity, site)

iwpc = read_csv("../iwpc_df.csv") %>% 
  mutate( vkor = factor(vkor, levels = c("G/G", "A/G", "A/A", "Missing"),
                        labels = c("GG", "AG", "AA", "NA")),
          sqrtdose = sqrt(dosewk),
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
  filter(sex != "Missing")%>% 
  dplyr::select(age,height,weight,
                dosewk,
                cyp,vkor,
                race, ethnicity,
                site)

tableonedat = full_join(iwpc, latinos) %>% 
  mutate(cohort = if_else(site %in% c("22", "23", "24", "25", "26", "27", "28", "29"), "ULLA", "IWPC")) %>% 
  dplyr::select(-site)

table1 = CreateTableOne(strata = "cohort", addOverall = F,data = tableonedat,
                        includeNA = F)

tab1mat =  print(table1, nonnormal = c("dosewk", "height", "weight"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

## Save to a CSV file
write.csv(tab1mat, file = "~/Documents/IWPC/Table1.csv")


