options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()

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
library(viridis)

set.seed(18)

set_bart_machine_num_cores(4)



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

## load and prep data
latinos = read_csv("~/Documents/IWPC/latinos.csv") %>% 
  mutate(sqrtdose = sqrt(dosewk),
         dosegroup = case_when(
           dosewk <= 21 ~  "low",
           dosewk >21 & dosewk < 49 ~ "intermediate",
           dosewk >= 49 ~ "high"
         ),
         BSA = sqrt((height*weight)/3600),
         race = factor(race, levels = c("white", "Black or African American", "Mixed or Missing", "Black"), labels = c("white",  "Black or African American", "Mixed or Missing", "Black or African American")),
         vkor = factor(vkor, levels = c("GG", "AG", "AA", "Missing")),
         indication = factor(indication, levels = c("AF", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA"), labels = c("AFIB", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA")),
         race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing"))) %>% 
  mutate_at(.vars = c("site", "sex", "amio", "ei", "target", "smoke", "diabetes", "statin",
                      "aspirin", "cyp", "indication", "dosegroup", "ethnicity"), .funs = factor) 



train = latinos %>% sample_frac(.7) %>% mutate(train = 1) 
latinos$train = if_else(latinos$X1 %in% train$X1, 1, 0)
test = latinos %>% filter(train == 0 )


latinosdf = latinos %>% 
  dplyr::select(-X1, - ethnicity)

'%nin%' <-Negate('%in%')

##Randomize testing and training
randTrainTesting<- function(data=latinosdf,colName="train"){
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
            sqrtdose ~ age + height + weight + vkor + cyp + race + 
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
  
  
  
  
  #### IWPC VARIABLES
  
  iwpcvars = lm(sqrtdose ~ age + height + weight + vkor + cyp + 
                  race + ei + amio , 
                data = train)
  
  test$pred_iwpcvars = predict(iwpcvars, test)
  
  R2_iwpcvars = 1 - (sum((test$pred_iwpcvars - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2))
  
  
  ################ ML ######################
  
  
  OptModelsvm = svm(sqrtdose ~ age + height + weight + vkor + cyp + 
                      race + ei + amio , 
                    data = train)
  
  
  test$pred_svm_iwpc =stats::predict(OptModelsvm,test)
  
  
  R2_svm_iwpc = 1 - (sum((test$pred_svm_iwpc^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## BART 
  ## Bayesian Additive Regression Trees (BART)
  
  x1 = train %>%
    dplyr::select(c( age , height , weight , vkor , cyp , 
                     race , ei ,amio )) %>%
    as.data.frame()
  
  bart <- bartMachine(x1, train$sqrtdose, mem_cache_for_speed = F)
  
  
  testx = test %>%
    dplyr::select(c( age , height , weight , vkor , cyp , 
                     race , ei ,amio )) %>%
    as.data.frame()
  
  test$pred_bart_iwpc = predict(bart, testx)
  
  R2_bart_iwpc = 1 - (sum((test$pred_bart_iwpc^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## MARS 
  # Fit a basic MARS model
  
  
  mars = train(sqrtdose ~ age + height + weight + vkor + cyp + 
                 race + ei + amio , 
               data = train,
               method = "earth",
               metric = "MAE"
  )
  
  test$pred_mars_iwpc = predict(mars, test)
  
  
  
  R2_mars_iwpc = 1 - (sum((test$pred_mars_iwpc^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  
  
  #### Kitchen skin model 
  newcovars = lm(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei   +sex + smoke  + statin+ aspirin+diabetes,data = train)
  
  test$pred_NLM = predict(newcovars, test)
  
  R2_newcovars = 1 - (sum((test$pred_newcovars - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2))
  
  
  
  ################ ML ######################
  
  
  OptModelsvm = svm(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei   +sex + smoke  + statin+ aspirin+diabetes,data = train)
  
  
  test$pred_svm =stats::predict(OptModelsvm,test)
  
  
  R2_svm = 1 - (sum((test$pred_svm^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## BART 
  ## Bayesian Additive Regression Trees (BART)
  
  x1 = train %>%
    dplyr::select(c(vkor,cyp ,age  , weight  , target , amio  ,indication , height ,race  , ei  , smoke  ,statin,aspirin, diabetes, sex)) %>%
    as.data.frame()
  
  bart <- bartMachine(x1, train$sqrtdose, mem_cache_for_speed = F)

  
 testx = test %>% 
    dplyr::select(c(vkor,cyp ,age  , weight  , target , amio  ,indication , height ,race  , ei  , smoke  ,statin, aspirin, diabetes, sex)) %>%
    as.data.frame()
  
  test$pred_bart = predict(bart, testx)
  
  R2_bart = 1 - (sum((test$pred_bart^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## MARS 
  # Fit a basic MARS model

  
  mars = train(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei   +sex + smoke  + statin+aspirin+ diabetes,data = train,
               method = "earth",
               metric = "MAE"
  )
  
  test$pred_mars = predict(mars, test)
  
  
  
  R2_mars = 1 - (sum((test$pred_mars^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  
  #### CLINICAL 
  fitclin = (lm(data = train, 
            sqrtdose ~ age + height + weight + race + 
              ei + amio))
  
  fitclin$coefficients <- c(4.0376, 
                        -0.2546,
                        0.0118,
                        0.0134,
                        -0.6752,
                        0.4060,
                        0.0443,
                        1.2799,
                        -0.5695)
  
  test$pred_clin = predict(fitclin, test )
  
  R2_clin = 1 - sum((test$pred_clin - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2)
  
  ################ results #################
  #########################################
  resultsdat = test %>% 
    mutate(resid_svm = pred_svm - sqrtdose,
           resid_bart = pred_bart - sqrtdose,
           resid_mars = pred_mars - sqrtdose,
           resid_svm_iwpc = pred_svm_iwpc -sqrtdose,
           resid_bart_iwpc = pred_bart_iwpc - sqrtdose,
           resid_mars_iwpc = pred_mars_iwpc - sqrtdose,
           resid_iwpccovars = pred.iwpcmodel -sqrtdose,
           resid_iwpcvars = pred_iwpcvars - sqrtdose,
           resid_newcovars= pred_NLM - sqrtdose,
           resid_clinical = pred_clin - sqrtdose,
           off_svm = if_else((abs(pred_svm^2 - dosewk)/dosewk) > .2, 
                             1, 0),
           off_mars = if_else((abs(pred_mars^2 - dosewk)/dosewk) > .2,
                              1,0),
           off_bart =if_else((abs(pred_bart^2 - dosewk)/dosewk) > .2,
                             1,0),
           off_svm_iwpc = if_else((abs(pred_svm_iwpc^2 - dosewk)/dosewk) > .2, 
                                  1, 0),
           off_mars_iwpc = if_else((abs(pred_mars_iwpc^2 - dosewk)/dosewk) > .2, 
                                   1, 0),
           off_bart_iwpc = if_else((abs(pred_bart_iwpc^2 - dosewk)/dosewk) > .2, 
                                   1, 0),
           off_iwpc = if_else((abs(pred.iwpcmodel^2 - dosewk)/dosewk) > .2,
                              1, 0),
           off_iwpcvars = if_else((abs(pred_iwpcvars^2 - dosewk)/dosewk) > .2,
                                  1, 0),
           off_newlinear = if_else((abs(pred_NLM^2 - dosewk )/dosewk) > .2,1,0),
           off_clinical = if_else((abs(pred_clin^2-dosewk)/dosewk)>.2,1,0)) %>% 
    pivot_longer(cols = starts_with("resid"), 
                 names_to = "model",
                 values_to = c("resid")) %>% 
    pivot_longer(cols = starts_with("pred"),
                 names_to = "type", 
                 values_to = "pred") %>% 
    dplyr::select(-type) %>% 
    distinct(resid, .keep_all = T) 
  
  
  mae_white = cbind(iwpc = mae2(test$dosewk[test$race == "white"], test$pred.iwpcmodel[test$race == "white"])$mae,
              iwpcvars = mae2(test$dosewk[test$race == "white"], test$pred_iwpcvars[test$race == "white"])$mae,
              svm_iwpc = mae2(test$dosewk[test$race == "white"], test$pred_svm_iwpc[test$race == "white"])$mae,
              bart_iwpc = mae2(test$dosewk[test$race == "white"], test$pred_bart_iwpc[test$race == "white"])$mae,
              mars_iwpc = mae2(test$dosewk[test$race == "white"], test$pred_mars_iwpc[test$race == "white"])$mae,
              newlinear = mae2(test$dosewk[test$race == "white"], test$pred_NLM[test$race == "white"])$mae,
              svm = mae2(test$dosewk[test$race == "white"], test$pred_svm[test$race == "white"])$mae,
              mars = mae2(test$dosewk[test$race == "white"], test$pred_mars[test$race == "white"])$mae,
              bart = mae2(test$dosewk[test$race == "white"], test$pred_bart[test$race == "white"])$mae,
              clinical= mae2(test$dosewk[test$race == "white"], test$pred_clin[test$race == "white"])$mae)
  

  mae_black = cbind(iwpc = mae2(test$dosewk[test$race == "Black or African American"], test$pred.iwpcmodel[test$race == "Black or African American"])$mae,
                    iwpcvars = mae2(test$dosewk[test$race == "Black or African American"], test$pred_iwpcvars[test$race == "Black or African American"])$mae,
                    svm_iwpc = mae2(test$dosewk[test$race == "Black or African American"], test$pred_svm_iwpc[test$race == "Black or African American"])$mae,
                    bart_iwpc = mae2(test$dosewk[test$race == "Black or African American"], test$pred_bart_iwpc[test$race == "Black or African American"])$mae,
                    mars_iwpc = mae2(test$dosewk[test$race == "Black or African American"], test$pred_mars_iwpc[test$race == "Black or African American"])$mae,
                    newlinear = mae2(test$dosewk[test$race == "Black or African American"], test$pred_NLM[test$race == "Black or African American"])$mae,
                    svm = mae2(test$dosewk[test$race == "Black or African American"], test$pred_svm[test$race == "Black or African American"])$mae,
                    mars = mae2(test$dosewk[test$race == "Black or African American"], test$pred_mars[test$race == "Black or African American"])$mae,
                    bart = mae2(test$dosewk[test$race == "Black or African American"], test$pred_bart[test$race == "Black or African American"])$mae,
                    clinical= mae2(test$dosewk[test$race == "Black or African American"], test$pred_clin[test$race == "Black or African American"])$mae)
  
  mae_mixed = cbind(iwpc = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred.iwpcmodel[test$race == "Mixed or Missing"])$mae,
                    iwpcvars = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_iwpcvars[test$race == "Mixed or Missing"])$mae,
                    svm_iwpc = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_svm_iwpc[test$race == "Mixed or Missing"])$mae,
                    bart_iwpc = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_bart_iwpc[test$race == "Mixed or Missing"])$mae,
                    mars_iwpc = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_mars_iwpc[test$race == "Mixed or Missing"])$mae,
                    newlinear = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_NLM[test$race == "Mixed or Missing"])$mae,
                    svm = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_svm[test$race == "Mixed or Missing"])$mae,
                    mars = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_mars[test$race == "Mixed or Missing"])$mae,
                    bart = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_bart[test$race == "Mixed or Missing"])$mae,
                    clinical= mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_clin[test$race == "Mixed or Missing"])$mae) 
  
  
  mae_ci_white = cbind(iwpc = mae2(test$dosewk[test$race == "white"], test$pred.iwpcmodel[test$race == "white"])$CI,
                 iwpcvars = mae2(test$dosewk[test$race == "white"], test$pred_iwpcvars[test$race == "white"])$CI,
                 svm_iwpc =  mae2(test$dosewk[test$race == "white"], test$pred_svm_iwpc[test$race == "white"])$CI,
                 bart_iwpc =  mae2(test$dosewk[test$race == "white"], test$pred_bart_iwpc[test$race == "white"])$CI,
                 mars_iwpc =  mae2(test$dosewk[test$race == "white"], test$pred_mars_iwpc[test$race == "white"])$CI,
                 newlinear = mae2(test$dosewk[test$race == "white"], test$pred_NLM[test$race == "white"])$CI,
                 svm = mae2(test$dosewk[test$race == "white"], test$pred_svm[test$race == "white"])$CI,
                 mars = mae2(test$dosewk[test$race == "white"], test$pred_mars[test$race == "white"])$CI,
                 bart = mae2(test$dosewk[test$race == "white"], test$pred_bart[test$race == "white"])$CI,
                 clinical = mae2(test$dosewk[test$race == "white"], test$pred_clin[test$race == "white"])$CI)
  
  mae_ci_black = cbind(iwpc = mae2(test$dosewk[test$race == "Black or African American"], test$pred.iwpcmodel[test$race == "Black or African American"])$CI,
                       iwpcvars = mae2(test$dosewk[test$race == "Black or African American"], test$pred_iwpcvars[test$race == "Black or African American"])$CI,
                       svm_iwpc =  mae2(test$dosewk[test$race == "Black or African American"], test$pred_svm_iwpc[test$race == "Black or African American"])$CI,
                       bart_iwpc =  mae2(test$dosewk[test$race == "Black or African American"], test$pred_bart_iwpc[test$race == "Black or African American"])$CI,
                       mars_iwpc =  mae2(test$dosewk[test$race == "Black or African American"], test$pred_mars_iwpc[test$race == "Black or African American"])$CI,
                       newlinear = mae2(test$dosewk[test$race == "Black or African American"], test$pred_NLM[test$race == "Black or African American"])$CI,
                       svm = mae2(test$dosewk[test$race == "Black or African American"], test$pred_svm[test$race == "Black or African American"])$CI,
                       mars = mae2(test$dosewk[test$race == "Black or African American"], test$pred_mars[test$race == "Black or African American"])$CI,
                       bart = mae2(test$dosewk[test$race == "Black or African American"], test$pred_bart[test$race == "Black or African American"])$CI,
                       clinical = mae2(test$dosewk[test$race == "Black or African American"], test$pred_clin[test$race == "Black or African American"])$CI)
  
  mae_ci_mixed = cbind(iwpc = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred.iwpcmodel[test$race == "Mixed or Missing"])$CI,
                       iwpcvars = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_iwpcvars[test$race == "Mixed or Missing"])$CI,
                       svm_iwpc =  mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_svm_iwpc[test$race == "Mixed or Missing"])$CI,
                       bart_iwpc =  mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_bart_iwpc[test$race == "Mixed or Missing"])$CI,
                       mars_iwpc =  mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_mars_iwpc[test$race == "Mixed or Missing"])$CI,
                       newlinear = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_NLM[test$race == "Mixed or Missing"])$CI,
                       svm = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_svm[test$race == "Mixed or Missing"])$CI,
                       mars = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_mars[test$race == "Mixed or Missing"])$CI,
                       bart = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_bart[test$race == "Mixed or Missing"])$CI,
                       clinical = mae2(test$dosewk[test$race == "Mixed or Missing"], test$pred_clin[test$race == "Mixed or Missing"])$CI)
  
  # 
  mae_dat = cbind(mae_white, mae_black, mae_mixed)
  mae_ci_dat = cbind(mae_ci_white, mae_ci_black, mae_ci_mixed)
  
  mae = t(rbind(mae_dat, mae_ci_dat))
  
  
  mae = mae %>% 
    as.data.frame() %>% 
    rownames_to_column("name") %>% 
    dplyr::select(name, MAE = V1, lower, upper) %>% 
    mutate(name = gsub("\\..*", "", name),
           name = factor(name,levels = c("iwpc", "iwpcvars", "svm_iwpc", "bart_iwpc", "mars_iwpc", "newlinear", "svm", "mars", "bart", "clinical"), labels = c("IWPC", "IWPCV", "IWPC_SVM", "IWPC_BART", "IWPC_MARS", "NLM", "SVM", "MARS", "BART", "CLINICAL")),
    race = rep(1:3,each=10),
    race = factor(race, levels = c("1", "2", "3"), labels = c("white", "Black or African American", "Mixed or Missing")))
  
  resultsdat = test %>% 
    pivot_longer(cols = starts_with("pred")) %>% 
    mutate(name = gsub("pred*.", "", name)) %>% 
    group_by(name, race) %>% 
    mutate(preddose = value^2,
           change = abs((preddose - dosewk)/dosewk),
           in20 = if_else(change < .2,1,0),
           prop = sum(in20)/n(),
    ) %>%
    ungroup() %>% 
    mutate(a = cut(change, breaks = 30)) %>% 
    group_by(name,race) %>% 
    mutate(sum = n()) %>% 
    group_by(name,a ,race) %>% 
    mutate(n = n()/sum)
  
  
  
  resultsdat2 = resultsdat %>% 
    group_by(name,race) %>% 
    summarise(prop = mean(prop)*100,
              n = sum(in20)) %>% 
    mutate(name = factor(name, levels = c("iwpcmodel", "iwpcvars","svm_iwpc", "bart_iwpc", "mars_iwpc" ,"NLM", "svm", "bart", "mars", "clin"), labels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_BART", "IWPC_MARS", "NLM", "SVM", "BART", "MARS", "CLINICAL"))) %>% 
    arrange(name) %>% 
    inner_join(mae, by = c("name","race")) %>% 
    mutate_if(is.numeric,round,2) 
  
  return(resultsdat2)
  
}

###############################################
#########FUNCTION FOR MULTIPLE RUNS###########
##############################################

fullReplicates <- lapply(1:100, function(x){
  cat('\rreplicate=',x)
  cbind.data.frame('replicate'=x, fitOneSet())
})


dat = do.call(rbind, fullReplicates) %>% 
  group_by(name,race) %>% 
  mutate(name = factor(name, levels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_MARS", "IWPC_BART", "NLM","SVM", "MARS", "BART", "CLINICAL"),
                       labels = c("IWPC", "IWPCV", "IWPC SVR", "IWPC MARS", "IWPC BART", "NLM", "SVR","MARS", "BART", "CLINICAL")))

 
dat2 = dat %>% 
  group_by(race, name) %>% 
  summarise(prop = mean(prop),
            MAE = mean(MAE),
            lower = mean(lower),
            upper = mean(upper)) 
#write_csv(dat2, "dats_ULLAbyrace.csv" )

hlinedat = dat %>% 
  group_by(race) %>% 
  summarise(prop = median(prop))

race_n <- c(
  'white'="White (n = 895)",
  'Black or African American'="Black or African American (n = 342)",
  'Mixed or Missing'="Mixed or Missing (n = 184)"
)

p3 = ggplot(dat, aes(name, prop, color = name)) +
  geom_boxplot(width=.25, alpha = .2) + 
  geom_jitter(width = .2, alpha = .1) + 
  theme_minimal()+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.title = element_text(size = 22)) + 
  labs(x = "", y = " ", color = "Model") +
  facet_wrap(race~., labeller = as_labeller(race_n))+
  geom_hline(data = hlinedat, aes(yintercept=prop),
             color = "gray30") +
  scale_color_manual(values = c("#F8766D" ,"#D39200", "#93AA00", "#00BA38", "#00C19F" ,"#00B9E3" ,"#619CFF", "#DB72FB", "#FF61C3", "gray30"))

pdf("fig2ulla_byrace.pdf", width = 8, height = 5)

p3

dev.off()


friedman_test(prop ~ name |replicate, data = dat)

prop_wilcox = as.data.frame(dat) %>%
  group_by(race) %>% 
  wilcox_test(prop ~ name, paired = TRUE, p.adjust.method = "bonferroni")

MAE_wilcox = as.data.frame(dat) %>%
  group_by(race) %>% 
  wilcox_test(MAE ~ name, paired = TRUE, p.adjust.method = "bonferroni")

wilcox_results_ULLA = rbind(prop_wilcox, MAE_wilcox)
write_csv(wilcox_results_ULLA, "wilcox_results_ULLA_byrace.csv" )


### put the plots together
library(cowplot)

fig2 = plot_grid(p3, nrow = 1) 
save_plot("figure2.png", fig2, base_width = 18, base_height = 8, bg = "transparent")

