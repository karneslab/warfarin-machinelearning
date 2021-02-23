library(readr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(bartMachine)
library(caret)
library(yardstick)
library(earth)
library(MASS)
library(e1071)

options(java.parameters = "-Xmx12g")
set_bart_machine_num_cores(4)

tgrid <- expand.grid(num_trees = 100,
                     k = 2,
                     alpha = 0.95, 
                     beta = 2,
                     nu =  3)

############################################################
################ADDITIONAL FUNCTIONS########################
################################################
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}
mae2 = function(actual, predicted, level=.95){
  
  d = actual-predicted^2
  mae = mean(abs(d))
  list('mae'=mae,
       'CI'=confidence_interval(abs(d), level)
  )
}

iwpc = readxl::read_xls("PS215192-569836982.xls", sheet = "Subject Data")[1:6256,] 

'%!in%' <-Negate('%in%')

###########################################
############### clean data #################
###########################################
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
                cyp_con = `CYP2C9 consensus`,
                vkor1639 = `VKORC1     -1639 consensus`, ## 9923231 
                vkor2255 = `VKORC1 2255 consensus`,  ## 2359612
                vkor1173 = `VKORC1 1173 consensus`, ## 9934438
                vkor1542 = `VKORC1 1542 consensus`, ## 8050894
                site = `Project Site`,
                dosewk = `Therapeutic Dose of Warfarin`,
                stable = `Subject Reached Stable Dose of Warfarin`,
                smoke = `Current Smoker`) %>% 
  na_if("NA") %>% 
  mutate(    cyp = if_else(cyp_con == "*1/*11",
                           "*1/*2", 
                           if_else(cyp_con %in% 
                                     c("*1/*5", 
                                       "*1/*6", 
                                       "*1/*13", 
                                       "*1/*14"),
                                   "*1/*3", 
                                   cyp_con)),
             height = as.numeric(height),
             weight = as.numeric(weight))



iwpc_df = iwpcdat %>% 
  mutate_at(.vars = c( "cyp", "race", "amio", "smoke"), .funs = as.factor) %>% 
  mutate(
    dosegroup = case_when(
      dosewk <= 21 ~  "low",
      dosewk >21 & dosewk < 49 ~ "intermediate",
      dosewk >= 49 ~ "high"
    )) %>% 
  mutate(dosegroup = factor(dosegroup),
         sqrtdose = sqrt(as.numeric(dosewk)),
         ei = if_else(carbamaz == "1" |
                        phenytoin == "1" |
                        rifampin == "1",
                      "1",
                      "0"),
         target= if_else(target == "2.25"| target == "2", "2.5", as.character(target)),
         race = factor(race, levels = c("White", "Asian", "Black or African American", "Unknown"), labels = c("white", "Asian", "Black or African American", "Mixed or Missing")),
         cyp= if_else(is.na(cyp), "Missing", as.character(cyp)),
         vkor1639 = if_else(is.na(vkor1639), "Missing", vkor1639),
         vkor1639 = factor(vkor1639, levels = c("G/G", "A/G", "A/A", "Missing"))) %>% 
  mutate(target = as.factor(target),
         ei = if_else(is.na(ei), "0", ei),
         amio = if_else(is.na(amio), "0", as.character(amio)),
         smoke = if_else(is.na(smoke), "0",as.character(smoke)),
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
         age = as.numeric(age)) %>% 
  mutate_at(.vars = c("amio", "cyp", "ei", "smoke"), .funs = factor) %>% 
  filter(
    !is.na(height),
    !is.na(weight),
    !is.na(age),
    !is.na(dosewk),
    !is.na(stable),
    stable != 0) 


iwpc_df2 = iwpc_df %>% 
  dplyr::select(age,
                height,
                weight,
                cyp,
                vkor1639,
                amio,
                ei,
                race,
                smoke,
                id,
                sqrtdose,
                dosegroup,
                dosewk) %>% 
  mutate(dosewk = as.numeric(dosewk))

###########################################
############### analysis #################
###########################################
### what percentage of the population is hispanic? 
iwpc_df %>% 
  group_by(ethnicity) %>% 
  summarize(n = n(),
            prop = n/nrow(iwpc_df))

### test/train data
train = iwpc_df2 %>% sample_frac(.8) %>% mutate(train = 1) 
iwpc_df2$train = if_else(iwpc_df2$id %in% train$id, 1, 0)
test = iwpc_df2 %>% filter(train == 0)


##Randomize testing and training
randTrainTesting<- function(data=iwpc_df2,colName="train"){
  data2<-as.data.frame(data)
  to_1<-sample(1: nrow(data2) ,length(which(data2[,colName] == 1)))
  data2[,colName]<- ifelse(as.numeric(row.names(data2)) %in% to_1, 1,0)
  data2
}


###############################################
#############FUNCTION FOR A SINGLE RUN#########
##############################################

fitOneSet <- function(){
  
  shu<-randTrainTesting()
  mnl<- min(unlist(lapply(shu, function(x) length(unique(x)))))
  
  train = shu[shu$train == 1,]
  test = shu[shu$train == 0,]
  
  itdoesntmatter = all(sort(as.character(unlist(lapply(Filter(is.factor, test), unique)))) %in% sort(as.character(unlist(lapply(Filter(is.factor, train), unique))) ) )
  
  itdoesntmatter
  
  if( mnl < 2 | itdoesntmatter == F){
    while (mnl < 2 | itdoesntmatter == F) {
      cat('\rTrying a new dataset')
      shu<-randTrainTesting()
      train = shu[shu$train == 1,]
      test = shu[shu$train == 0,]
      
      itdoesntmatter = all(sort(as.character(unlist(lapply(Filter(is.factor, test), unique)))) %in% sort(as.character(unlist(lapply(Filter(is.factor, train), unique))) ) )
      print(itdoesntmatter)
      mnl<- min(unlist(lapply(shu, function(x) length(unique(x)))))
    }
  }
  
  cat('\nOptimal dataset found\n')
  
#### IWPC LM 
fit = (lm(data = train, 
          sqrtdose ~ age + height + weight + vkor1639 + cyp + race + 
            ei + amio))

fit$coefficients <- c(5.6044, 
                      -0.2614,
                      0.0087,
                      0.0128,
                      -0.8677,
                      -1.6974,
                      -0.4854,
                      -0.5211,
                      -0.9357,
                      -1.0616,
                      -1.9206,
                      -2.3312,
                      -0.2188,
                      -0.1092,
                      -0.2760,
                      -0.1032,
                      1.1816,
                      -0.5503)

test$pred.iwpcmodel = predict(fit, test )

R2_iwpc = 1 - sum((test$pred.iwpcmodel - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2)


### by dose group, what percentage of patients was predicted w/in 20% of their actual dose
test %>% 
  mutate(dose = as.numeric(dosewk),
         pred = pred.iwpcmodel^2,
         change = abs(pred- dose)/dose,
         in20 = if_else(change > .2, 0, 1)) %>% 
  # group_by(dosegroup) %>% 
  summarise(n = n(),
            prop = n/nrow(test),
            nin = sum(in20),
            propin = nin/n)


#### Linear model 
newcovars = lm(sqrtdose ~ race + vkor1639 + cyp + age + height + weight + smoke + ei + amio, data = train)

test$pred_newcovars = predict(newcovars, test)

R2_newcovars = 1 - (sum((test$pred_newcovars - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2))


################ ML ######################
## SVM 

#Tune the SVM model
OptModelsvm=svm(sqrtdose ~ race + vkor1639 + cyp + age + height + weight + smoke + ei + amio, data = train)


#Predict Y using best model
test$pred_svm =predict(OptModelsvm,test)



test %>% 
  mutate(dose = sqrtdose^2,
         pred = pred_svm^2,
         change = abs(pred_svm^2 - dose)/dose,
         in20 = if_else(change > .2, 0, 1)) %>% 
  group_by(dosegroup) %>% 
  summarise(n = n(),
            prop = n/nrow(test),
            nin = sum(in20),
            propin = nin/n)

R2_svm = 1 - (sum((test$pred_svm^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))


## BART 
## Bayesian Additive Regression Trees (BART)


# trainbart = train %>% 
#   dplyr::select(sqrtdose,
#                 race,
#                 vkor1639,
#                 cyp,
#                 age,
#                 height,
#                 weight,
#                 smoke,
#                 ei,
#                 amio) %>% 
#   as.data.frame()
# 
# bart <- bartMachine(trainbart[,-1], trainbart$sqrtdose)
# 
# testbart = test %>% 
#   dplyr::select(sqrtdose,
#                 race,
#                 vkor1639,
#                 cyp,
#                 age,
#                 height,
#                 weight,
#                 smoke,
#                 ei,
#                 amio) %>% 
#   as.data.frame()
# 
# test$pred_bart = predict(bart, testbart)
# 
# 
# test %>% 
#   mutate(dose = dosewk,
#          pred = pred_bart^2,
#          change = abs(pred- dose)/dose,
#          in20 = if_else(change > .2, 0, 1)) %>% 
#   group_by(dosegroup) %>% 
#   summarise(n = n(),
#             prop = n/nrow(test),
#             nin = sum(in20),
#             propin = nin/n)
# 
# R2_bart = 1 - (sum((test$pred_bart^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))

## MARS 
# Fit a basic MARS model
mars =  train(sqrtdose ~ race + vkor1639 + cyp + age + height + weight + smoke + ei + amio, data = train,
              method = "earth",
              metric = "MAE"
)


test$pred_mars = predict(mars, test)


test %>% 
  mutate(dose = dosewk,
         pred = pred_mars^2,
         change = abs(pred- dose)/dose,
         in20 = if_else(change > .2, 0, 1)) %>% 
  group_by(dosegroup) %>% 
  summarise(n = n(),
            prop = n/nrow(test),
            nin = sum(in20),
            propin = nin/n)

R2_mars = 1 - (sum((test$pred_mars^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))


################ results #################
#########################################
resultsdat = test %>% 
  mutate(resid_svm = pred_svm - sqrtdose,
         # resid_bart = pred_bart - sqrtdose,
         resid_mars = pred_mars - sqrtdose,
         resid_newcovars= pred_newcovars - sqrtdose,
         off_svm = if_else((abs(pred_svm^2 - dosewk)/dosewk) > .2, 
                           1, 0),
         off_mars = if_else((abs(pred_mars^2 - dosewk)/dosewk) > .2,
                            1,0),
         # off_bart =if_else((abs(pred_bart^2 - dosewk)/dosewk) > .2,
         #                   1,0),
         off_iwpc = if_else((abs(pred.iwpcmodel^2 - dosewk)/dosewk) > .2,
                            1, 0),
         off_newlinear = if_else((abs(pred_newcovars^2 - dosewk )/dosewk) > .2,1,0)) %>% 
  pivot_longer(cols = starts_with("resid"), 
               names_to = "model",
               values_to = c("resid")) %>% 
  pivot_longer(cols = starts_with("pred"),
               names_to = "type", 
               values_to = "pred") %>% 
  dplyr::select(-type) %>% 
  distinct(resid, .keep_all = T) 


mae = cbind(iwpc = mae2(test$dosewk, test$pred.iwpcmodel)$mae,
            newlinear = mae2(test$dosewk, test$pred_newcovars)$mae,
            svm = mae2(test$dosewk, test$pred_svm)$mae,
            mars = mae2(test$dosewk, test$pred_mars)$mae)
            # bart = mae2(test$dosewk, test$pred_bart)$mae)

mae_low = cbind(iwpc = mae2(test$dosewk, test$pred.iwpcmodel[test$dosegroup == "low"])$mae,
                newlinear = mae2(test$dosewk[test$dosegroup == "low"], test$pred_newcovars[test$dosegroup == "low"])$mae,
                svm = mae2(test$dosewk[test$dosegroup == "low"], test$pred_svm[test$dosegroup == "low"])$mae,
                mars = mae2(test$dosewk[test$dosegroup == "low"], test$pred_mars[test$dosegroup == "low"])$mae)
                # bart = mae2(test$dosewk[test$dosegroup == "low"], test$pred_bart[test$dosegroup == "low"])$mae)


mae_int = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred.iwpcmodel[test$dosegroup == "intermediate"])$mae,
                newlinear = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_newcovars[test$dosegroup == "intermediate"])$mae,
                svm = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_svm[test$dosegroup == "intermediate"])$mae,
                mars = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_mars[test$dosegroup == "intermediate"])$mae)
                # bart = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_bart[test$dosegroup == "intermediate"])$mae)


mae_high = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "high"], test$pred.iwpcmodel[test$dosegroup == "high"])$mae,
                 newlinear = mae2(test$dosewk[test$dosegroup == "high"], test$pred_newcovars[test$dosegroup == "high"])$mae,
                 svm = mae2(test$dosewk[test$dosegroup == "high"], test$pred_svm[test$dosegroup == "high"])$mae,
                 mars = mae2(test$dosewk[test$dosegroup == "high"], test$pred_mars[test$dosegroup == "high"])$mae)
                 # bart = mae2(test$dosewk[test$dosegroup == "high"], test$pred_bart[test$dosegroup == "high"])$mae)


mae_ci = cbind(iwpc = mae2(test$dosewk, test$pred.iwpcmodel)$CI,               newlinear = mae2(test$dosewk, test$pred_newcovars)$CI,
               svm = mae2(test$dosewk, test$pred_svm)$CI,
               mars = mae2(test$dosewk, test$pred_mars)$CI)
               # bart = mae2(test$dosewk, test$pred_bart)$CI)


mae_ci_low = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred.iwpcmodel[test$dosegroup == "low"])$CI,
                   newlinear = mae2(test$dosewk[test$dosegroup == "low"], test$pred_newcovars[test$dosegroup == "low"])$CI,
                   svm = mae2(test$dosewk[test$dosegroup == "low"], test$pred_svm[test$dosegroup == "low"])$CI,
                   mars = mae2(test$dosewk[test$dosegroup == "low"], test$pred_mars[test$dosegroup == "low"])$CI)
                   # bart = mae2(test$dosewk[test$dosegroup == "low"], test$pred_bart[test$dosegroup == "low"])$CI)


mae_ci_int = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred.iwpcmodel[test$dosegroup == "intermediate"])$CI,
                   newlinear = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_newcovars[test$dosegroup == "intermediate"])$CI,
                   svm = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_svm[test$dosegroup == "intermediate"])$CI,
                   mars = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_mars[test$dosegroup == "intermediate"])$CI)
                   # bart = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_bart[test$dosegroup == "intermediate"])$CI)



mae_ci_high = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "high"], test$pred.iwpcmodel[test$dosegroup == "high"])$CI,
                    newlinear = mae2(test$dosewk[test$dosegroup == "high"], test$pred_newcovars[test$dosegroup == "high"])$CI,
                    svm = mae2(test$dosewk[test$dosegroup == "high"], test$pred_svm[test$dosegroup == "high"])$CI,
                    mars = mae2(test$dosewk[test$dosegroup == "high"], test$pred_mars[test$dosegroup == "high"])$CI)
                    # bart = mae2(test$dosewk[test$dosegroup == "high"], test$pred_bart[test$dosegroup == "high"])$CI)

mae_low_dat = t(rbind(mae_low, mae_ci_low))
mae_int_dat = t(rbind(mae_int, mae_ci_int))

mae_high_dat = t(rbind(mae_high, mae_ci_high))

mae_dosegroup_dat = cbind(low = mae_low_dat,
                          int = mae_int_dat,
                          high = mae_high_dat) %>% 
  as.data.frame() %>% 
  `colnames<-` (c("low","low_low","low_high","int", "int_low","int_high", "high", "high_low","high_high")) %>% 
  mutate(across(everything(), round, 2))

mae_dosegroup_dat = mae_dosegroup_dat %>% 
  as.data.frame() 

mae_dat = t(rbind(mae, mae_ci))

mae_dat = mae_dat %>% 
  as.data.frame() %>% 
  rownames_to_column("type") %>% 
  dplyr::select(type, MAE = V1, lower, upper)

resultsdat = test %>% 
  pivot_longer(cols = starts_with("pred")) %>% 
  mutate(name = gsub("pred*.", "", name)) %>% 
  group_by(name) %>% 
  mutate(preddose = value^2,
         change = abs((preddose - dosewk)/dosewk),
         in20 = if_else(change < .2,1,0),
         prop = sum(in20)/n(),
  ) %>%
  ungroup() %>% 
  mutate(a = cut(change, breaks = 30)) %>% 
  group_by(name, dosegroup) %>% 
  mutate(sum = n()) %>% 
  group_by(name, dosegroup,a ) %>% 
  mutate(n = n()/sum)


resultsdat2 = resultsdat %>% 
  group_by(name) %>% 
  summarise(prop = mean(prop)*100,
            n = sum(in20)) %>% 
  cbind(mae_dat) %>% 
  dplyr::select(-type)%>% 
  mutate_if(is.numeric,round,2) %>% 
  mutate(name = factor(name, levels = c("iwpcmodel", "newcovars", "svm", "mars")))

return(resultsdat2 %>% arrange(name) )

}


###############################################
#########FUNCTION FOR MULTIPLE RUNS###########
##############################################

resultsdat2 %>% 
  unite("CI", lower:upper, sep = "-") %>% 
  unite(`MAE (95% CI)`, MAE:CI, sep = " (") 

fullReplicates <- lapply(1:10, function(x){
  cat('\rreplicate=',x)
  cbind.data.frame('replicate'=x, fitOneSet())
})

dat = do.call(rbind, fullReplicates) %>% 
  group_by(name) %>% 
  mutate(name = factor(name, levels = c("iwpcmodel", "iwpcvars", "newcovars", "svm", "bart", "mars"), labels = c("IWPC", "IWPCV", "NLM", "SVR","BART","MARS")))


plot = ggboxplot(dat, x = "name", y = "prop", add = "jitter")

my_comparisons <- list( c("IWPC", "NLM"), c("IWPC", "SVR"), c("NLM", "SVR") )

plot2 = ggplot(dat, aes(name, prop, color = name)) +
  geom_boxplot() +
  geom_jitter(width = .2, alpha = .2) + 
  theme_minimal()+
  theme(legend.position = "none")+
  labs(x = "Model", y = "Proportion within 20%")+
  stat_compare_means(comparisons = my_comparisons)

pdf("merged100plot.pdf", width = 8, height = 5)
plot2
dev.off()

friedman_test(prop ~ name |replicate, data = as.data.frame(dat))

as.data.frame(dat) %>%
  wilcox_test(prop ~ name, paired = TRUE, p.adjust.method = "bonferroni")

kruskal.test(prop ~name, data = dat)

dat4 = dat %>% 
  filter(name == "NLM" |
           name == "SVR") %>% 
  as.data.frame()

t.test(data = dat4,prop ~ name)




library(tableone)
tableonedf = iwpc_df %>% 
  dplyr::select(dosewk,
                height,
                weight,
                race,
                vkor1639,
                cyp,
                age,
                sex,
                ei,
                amio,
                smoke) %>% 
  mutate(dosewk = as.numeric(dosewk))

table1 = CreateTableOne( data = tableonedf,
                        includeNA = F)

