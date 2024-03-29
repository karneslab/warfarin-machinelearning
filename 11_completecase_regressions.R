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

set.seed(18)
options(java.parameters = "-Xmx20g")
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

mergeddat = read_csv("../merged_iwpc_ULLAcc.csv")

###########################################################
#####################DATA##################################
###########################################################

mergeddat = mergeddat %>% 
  mutate(
         sqrtdose = sqrt(dosewk),
         race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing", "Mixed/NA", "Black"), labels = c("white", "Asian", "Black or African American", "Mixed or Missing", "Mixed or Missing", "Black or African American")),
         vkor = factor(vkor, levels = c("GG", "AG", "AA", "Missing"))) %>% 
  mutate_at(.vars= c("cyp", "amio", "ei"),.funs = as_factor) %>% 
  dplyr::select(-country)


train = mergeddat %>% sample_frac(.7) %>% mutate(train = 1) 
mergeddat$train = if_else(mergeddat$X1 %in% train$X1, 1, 0)


'%nin%' <-Negate('%in%')

##Randomize testing and training
randTrainTesting<- function(data=mergeddat,colName="train"){
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
  
  bart <- bartMachine(x1, train$sqrtdose)
  
  
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
  
  

  ################ results #################
  #########################################
  resultsdat = test %>% 
    mutate(
           resid_svm_iwpc = pred_svm_iwpc -sqrtdose,
           resid_bart_iwpc = pred_bart_iwpc - sqrtdose,
           resid_mars_iwpc = pred_mars_iwpc - sqrtdose,
           resid_iwpccovars = pred.iwpcmodel -sqrtdose,
           resid_iwpcvars = pred_iwpcvars - sqrtdose,
        
           
           off_svm_iwpc = if_else((abs(pred_svm_iwpc^2 - dosewk)/dosewk) > .2, 
                                  1, 0),
           off_mars_iwpc = if_else((abs(pred_mars_iwpc^2 - dosewk)/dosewk) > .2, 
                                   1, 0),
           off_bart_iwpc = if_else((abs(pred_bart_iwpc^2 - dosewk)/dosewk) > .2, 
                                   1, 0),
           off_iwpc = if_else((abs(pred.iwpcmodel^2 - dosewk)/dosewk) > .2,
                              1, 0),
           off_iwpcvars = if_else((abs(pred_iwpcvars^2 - dosewk)/dosewk) > .2,
                                  1, 0)) %>% 
    pivot_longer(cols = starts_with("resid"), 
                 names_to = "model",
                 values_to = c("resid")) %>% 
    pivot_longer(cols = starts_with("pred"),
                 names_to = "type", 
                 values_to = "pred") %>% 
    dplyr::select(-type) %>% 
    distinct(resid, .keep_all = T) 
  
  maedat = test %>% 
    dplyr::select(starts_with("pred"), dosewk) %>% 
    pivot_longer(cols = pred.iwpcmodel:pred_mars_iwpc, values_to = "pred", names_to= "model") %>% 
    group_by( model) 
  
  mae = cbind(iwpc = mae2(test$dosewk, test$pred.iwpcmodel)$mae,
              iwpcvars = mae2(test$dosewk, test$pred_iwpcvars)$mae,
              svm_iwpc = mae2(test$dosewk, test$pred_svm_iwpc)$mae,
              bart_iwpc = mae2(test$dosewk, test$pred_bart_iwpc)$mae,
              mars_iwpc = mae2(test$dosewk, test$pred_mars_iwpc)$mae)
  
  
  mae_ci = cbind(iwpc = mae2(test$dosewk, test$pred.iwpcmodel)$CI,
                 iwpcvars = mae2(test$dosewk, test$pred_iwpcvars)$CI,
                 svm_iwpc =  mae2(test$dosewk, test$pred_svm_iwpc)$CI,
                 bart_iwpc =  mae2(test$dosewk, test$pred_bart_iwpc)$CI,
                 mars_iwpc =  mae2(test$dosewk, test$pred_mars_iwpc)$CI)
  

  mae_dat = t(rbind(mae, mae_ci))
  
  
  mae_dat = mae_dat %>% 
    as.data.frame() %>% 
    rownames_to_column("name") %>% 
    dplyr::select(name, MAE = V1, lower, upper) %>% 
    mutate(name = factor(name, levels = c("iwpc", "iwpcvars", "svm_iwpc", "bart_iwpc", "mars_iwpc"), labels = c("IWPC", "IWPCV", "IWPC_SVM", "IWPC_BART", "IWPC_MARS")))
  
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
    group_by(name) %>% 
    mutate(sum = n()) %>% 
    group_by(name,a ) %>% 
    mutate(n = n()/sum)
  
  
  
  resultsdat2 = resultsdat %>% 
    group_by(name) %>% 
    summarise(prop = mean(prop)*100,
              n = sum(in20)) %>% 
    mutate(name = factor(name, levels = c("iwpcmodel", "iwpcvars","svm_iwpc", "bart_iwpc", "mars_iwpc" ,"NLM", "svm", "bart", "mars"), labels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_BART", "IWPC_MARS", "NLM", "SVM", "BART", "MARS"))) %>% 
    arrange(name) %>% 
    # inner_join(mae_ethnicity_dat, by = c("ethnicity","name" )) %>% 
    inner_join(mae_dat, by = c("name")) %>% 
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
  group_by(name) %>% 
  mutate(name = factor(name, levels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_MARS", "IWPC_BART"),
                       labels = c("IWPC", "IWPCV",  "IWPC SVR","IWPC MARS", "IWPC BART")))

dat2= dat %>% 
  group_by(name) %>% 
  summarise(prop = median(prop),
            MAE = median(MAE),
            lower= median(lower),
            upper = median(upper))%>% 
  unite("CI", lower:upper , sep = "-") %>% 
  unite("MAE (95% CI)", MAE:CI, sep = "(")%>% 
  mutate(`MAE (95% CI)`= gsub("\\(", " (", `MAE (95% CI)`),
         `MAE (95% CI)` = paste(`MAE (95% CI)`,")", sep = ""))


#write_csv(dat2, "../dats_completecase.csv" )
### use same colors as full plot
library(scales)
show_col(hue_pal()(9))


p2 = ggplot(dat, aes(name, prop, color = name)) +
  geom_boxplot(width=.1, alpha = .1) + 
  geom_jitter(width = .2, alpha = .2) + 
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 15)) +
  labs(x = "", y = " ") +
  scale_color_manual(values = c("#F8766D", "#D39200", "#93AA00","#00BA38", "#00C19F"))

pdf("completecase.pdf", width = 8, height = 5)
p2

dev.off()

friedman_test(prop ~ name |replicate, data = as.data.frame(dat))

prop_wilcox = as.data.frame(dat) %>%
  wilcox_test(prop ~ name, paired = TRUE, p.adjust.method = "bonferroni")

MAE_wilcox = as.data.frame(dat) %>%
  wilcox_test(MAE ~ name, paired = TRUE, p.adjust.method = "bonferroni")



wilcox_results_Merged = rbind(prop_wilcox, MAE_wilcox)
  write_csv(wilcox_results_Merged, "../wilcox_results_Merged_completecase.csv" )

