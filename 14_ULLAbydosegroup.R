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
latinos = read_csv("../ULLA.csv") %>% 
  mutate(sqrtdose = sqrt(dosewk),
         dosegroup = case_when(
           dosewk <= 21 ~  "low",
           dosewk >21 & dosewk < 49 ~ "intermediate",
           dosewk >= 49 ~ "high"
         ),
         BSA = sqrt((height*weight)/3600),
         race = factor(race, levels = c("white", "Black or African American", "Mixed or Missing", "Black"), labels = c("white",  "Black or African American", "Mixed or Missing", "Black or African American")),
         vkor = factor(vkor, levels = c("GG", "AG","GA", "AA", "Missing"),
                       labels = c("GG", "AG", "AG","AA", "Missing")),
         indication = factor(indication, levels = c("AF", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA"), labels = c("AFIB", "AFIB", "DVT/PE", "MVR", "OTHER", "TIA")),
         race = factor(race, levels = c("white", "Asian", "Black or African American", "Mixed or Missing"))) %>% 
  mutate_at(.vars = c("site", "sex", "amio", "ei", "target", "smoke", "diabetes", "statin",
                      "aspirin", "cyp", "indication", "dosegroup", "ethnicity"), .funs = factor) %>% 
  dplyr::select(-target)



train = latinos %>% sample_frac(.7) %>% mutate(train = 1) 
latinos$train = if_else(latinos$X1 %in% train$X1, 1, 0)
test = latinos %>% filter(train == 0 )


latinosdf = latinos %>% 
  dplyr::select( - ethnicity, -site)

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
  newcovars = lm(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei   +sex + smoke  + statin+ aspirin+diabetes + country,data = train)
  
  test$pred_NLM = predict(newcovars, test)
  
  R2_newcovars = 1 - (sum((test$pred_NLM - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2))
  
  
  
  ################ ML ######################
  
  
  OptModelsvm = svm(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei   +sex + smoke  + statin+ aspirin+diabetes + country,data = train)
  
  
  test$pred_svm =stats::predict(OptModelsvm,test)
  
  
  R2_svm = 1 - (sum((test$pred_svm^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## BART 
  ## Bayesian Additive Regression Trees (BART)
  
  x1 = train %>%
    dplyr::select(c(vkor,cyp ,age  , weight   , amio  ,indication , height ,race  , ei  , smoke  ,statin,aspirin, diabetes, sex, country)) %>%
    as.data.frame()
  
  bart <- bartMachine(x1, train$sqrtdose, mem_cache_for_speed = F)
  
  
  testx = test %>% 
    dplyr::select(c(vkor,cyp ,age  , weight   , amio  ,indication , height ,race  , ei  , smoke  ,statin, aspirin, diabetes, sex,country)) %>%
    as.data.frame()
  
  test$pred_bart = predict(bart, testx)
  
  R2_bart = 1 - (sum((test$pred_bart^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## MARS 
  # Fit a basic MARS model
  
  
  mars = train(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei   +sex + smoke  + statin+aspirin+ diabetes+ country,data = train,
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
                            0.4060,
                            0.0443,
                            1.2799,
                            -0.5695)
  
  test$pred_clin = predict(fitclin, test )
  
  R2_clin = 1 - sum((test$pred_clin - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2)
  
  ################ results #################
  #########################################
  
  
  mae_high = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "high"], test$pred.iwpcmodel[test$dosegroup == "high"])$mae,
                    iwpcvars = mae2(test$dosewk[test$dosegroup == "high"], test$pred_iwpcvars[test$dosegroup == "high"])$mae,
                    svm_iwpc = mae2(test$dosewk[test$dosegroup == "high"], test$pred_svm_iwpc[test$dosegroup == "high"])$mae,
                    bart_iwpc = mae2(test$dosewk[test$dosegroup == "high"], test$pred_bart_iwpc[test$dosegroup == "high"])$mae,
                    mars_iwpc = mae2(test$dosewk[test$dosegroup == "high"], test$pred_mars_iwpc[test$dosegroup == "high"])$mae,
                    newlinear = mae2(test$dosewk[test$dosegroup == "high"], test$pred_NLM[test$dosegroup == "high"])$mae,
                    svm = mae2(test$dosewk[test$dosegroup == "high"], test$pred_svm[test$dosegroup == "high"])$mae,
                    mars = mae2(test$dosewk[test$dosegroup == "high"], test$pred_mars[test$dosegroup == "high"])$mae,
                    bart = mae2(test$dosewk[test$dosegroup == "high"], test$pred_bart[test$dosegroup == "high"])$mae,
                    clinical= mae2(test$dosewk[test$dosegroup == "high"], test$pred_clin[test$dosegroup == "high"])$mae)
  
  
  mae_int = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred.iwpcmodel[test$dosegroup == "intermediate"])$mae,
                    iwpcvars = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_iwpcvars[test$dosegroup == "intermediate"])$mae,
                    svm_iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_svm_iwpc[test$dosegroup == "intermediate"])$mae,
                    bart_iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_bart_iwpc[test$dosegroup == "intermediate"])$mae,
                    mars_iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_mars_iwpc[test$dosegroup == "intermediate"])$mae,
                    newlinear = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_NLM[test$dosegroup == "intermediate"])$mae,
                    svm = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_svm[test$dosegroup == "intermediate"])$mae,
                    mars = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_mars[test$dosegroup == "intermediate"])$mae,
                    bart = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_bart[test$dosegroup == "intermediate"])$mae,
                    clinical= mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_clin[test$dosegroup == "intermediate"])$mae)
  
  mae_low = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred.iwpcmodel[test$dosegroup == "low"])$mae,
                    iwpcvars = mae2(test$dosewk[test$dosegroup == "low"], test$pred_iwpcvars[test$dosegroup == "low"])$mae,
                    svm_iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred_svm_iwpc[test$dosegroup == "low"])$mae,
                    bart_iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred_bart_iwpc[test$dosegroup == "low"])$mae,
                    mars_iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred_mars_iwpc[test$dosegroup == "low"])$mae,
                    newlinear = mae2(test$dosewk[test$dosegroup == "low"], test$pred_NLM[test$dosegroup == "low"])$mae,
                    svm = mae2(test$dosewk[test$dosegroup == "low"], test$pred_svm[test$dosegroup == "low"])$mae,
                    mars = mae2(test$dosewk[test$dosegroup == "low"], test$pred_mars[test$dosegroup == "low"])$mae,
                    bart = mae2(test$dosewk[test$dosegroup == "low"], test$pred_bart[test$dosegroup == "low"])$mae,
                    clinical= mae2(test$dosewk[test$dosegroup == "low"], test$pred_clin[test$dosegroup == "low"])$mae) 
  
  
  #### MAE CIS 
  
  mae_ci_high = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "high"], test$pred.iwpcmodel[test$dosegroup == "high"])$CI,
                       iwpcvars = mae2(test$dosewk[test$dosegroup == "high"], test$pred_iwpcvars[test$dosegroup == "high"])$CI,
                       svm_iwpc =  mae2(test$dosewk[test$dosegroup == "high"], test$pred_svm_iwpc[test$dosegroup == "high"])$CI,
                       bart_iwpc =  mae2(test$dosewk[test$dosegroup == "high"], test$pred_bart_iwpc[test$dosegroup == "high"])$CI,
                       mars_iwpc =  mae2(test$dosewk[test$dosegroup == "high"], test$pred_mars_iwpc[test$dosegroup == "high"])$CI,
                       newlinear = mae2(test$dosewk[test$dosegroup == "high"], test$pred_NLM[test$dosegroup == "high"])$CI,
                       svm = mae2(test$dosewk[test$dosegroup == "high"], test$pred_svm[test$dosegroup == "high"])$CI,
                       mars = mae2(test$dosewk[test$dosegroup == "high"], test$pred_mars[test$dosegroup == "high"])$CI,
                       bart = mae2(test$dosewk[test$dosegroup == "high"], test$pred_bart[test$dosegroup == "high"])$CI,
                       clinical = mae2(test$dosewk[test$dosegroup == "high"], test$pred_clin[test$dosegroup == "high"])$CI)
  
  mae_ci_int = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred.iwpcmodel[test$dosegroup == "intermediate"])$CI,
                            iwpcvars = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_iwpcvars[test$dosegroup == "intermediate"])$CI,
                            svm_iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_svm_iwpc[test$dosegroup == "intermediate"])$CI,
                            bart_iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_bart_iwpc[test$dosegroup == "intermediate"])$CI,
                            mars_iwpc = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_mars_iwpc[test$dosegroup == "intermediate"])$CI,
                            newlinear = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_NLM[test$dosegroup == "intermediate"])$CI,
                            svm = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_svm[test$dosegroup == "intermediate"])$CI,
                            mars = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_mars[test$dosegroup == "intermediate"])$CI,
                            bart = mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_bart[test$dosegroup == "intermediate"])$CI,
                            clinical= mae2(test$dosewk[test$dosegroup == "intermediate"], test$pred_clin[test$dosegroup == "intermediate"])$CI)
  
  mae_ci_low = cbind(iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred.iwpcmodel[test$dosegroup == "low"])$CI,
                  iwpcvars = mae2(test$dosewk[test$dosegroup == "low"], test$pred_iwpcvars[test$dosegroup == "low"])$CI,
                  svm_iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred_svm_iwpc[test$dosegroup == "low"])$CI,
                  bart_iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred_bart_iwpc[test$dosegroup == "low"])$CI,
                  mars_iwpc = mae2(test$dosewk[test$dosegroup == "low"], test$pred_mars_iwpc[test$dosegroup == "low"])$CI,
                  newlinear = mae2(test$dosewk[test$dosegroup == "low"], test$pred_NLM[test$dosegroup == "low"])$CI,
                  svm = mae2(test$dosewk[test$dosegroup == "low"], test$pred_svm[test$dosegroup == "low"])$CI,
                  mars = mae2(test$dosewk[test$dosegroup == "low"], test$pred_mars[test$dosegroup == "low"])$CI,
                  bart = mae2(test$dosewk[test$dosegroup == "low"], test$pred_bart[test$dosegroup == "low"])$CI,
                  clinical= mae2(test$dosewk[test$dosegroup == "low"], test$pred_clin[test$dosegroup == "low"])$CI)
  
  # 
  mae_dat = cbind(mae_high, mae_int, mae_low)
  mae_ci_dat = cbind(mae_ci_high, mae_ci_int, mae_ci_low)
  
  mae = t(rbind(mae_dat, mae_ci_dat))
  
  
  mae = mae %>% 
    as.data.frame() %>% 
    rownames_to_column("name") %>% 
    dplyr::select(name, MAE = V1, lower, upper) %>% 
    mutate(name = gsub("\\..*", "", name),
           name = factor(name,levels = c("iwpc", "iwpcvars", "svm_iwpc", "bart_iwpc", "mars_iwpc", "newlinear", "svm", "mars", "bart", "clinical"), labels = c("IWPC", "IWPCV", "IWPC_SVM", "IWPC_BART", "IWPC_MARS", "NLM", "SVM", "MARS", "BART", "CLINICAL")),
           dosegroup = rep(1:3,each=10),
           dosegroup = factor(dosegroup, levels = c("1", "2", "3"), labels = c("High", "Intermediate", "Low")))
  
  resultsdat = test %>% 
    pivot_longer(cols = starts_with("pred")) %>% 
    mutate(name = gsub("pred*.", "", name)) %>% 
    group_by(name, dosegroup) %>% 
    mutate(preddose = value^2,
           change = abs((preddose - dosewk)/dosewk),
           in20 = if_else(change < .2,1,0),
           prop = sum(in20)/n(),
    ) %>%
    ungroup() %>% 
    mutate(a = cut(change, breaks = 30)) %>% 
    group_by(name,dosegroup) %>% 
    mutate(sum = n()) %>% 
    group_by(name,a ,dosegroup) %>% 
    mutate(n = n()/sum)
  
 
  
  
  resultsdat2 = resultsdat %>% 
    group_by(name,dosegroup) %>% 
    summarise(prop = mean(prop)*100,
              n = sum(in20)) %>% 
    mutate(name = factor(name, levels = c("iwpcmodel", "iwpcvars","svm_iwpc", "bart_iwpc", "mars_iwpc" ,"NLM", "svm", "bart", "mars", "clin"), labels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_BART", "IWPC_MARS", "NLM", "SVM", "BART", "MARS", "CLINICAL")),
           dosegroup = factor(dosegroup, levels = c("high", "intermediate", "low"),
                              labels = c("High","Intermediate" ,"Low"))) %>% 
    arrange(name) %>% 
    inner_join(mae, by = c("name","dosegroup")) %>% 
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
  group_by(name,dosegroup) %>% 
  mutate(name = factor(name, levels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_MARS", "IWPC_BART", "NLM","SVM", "MARS", "BART", "CLINICAL"),
                       labels = c("IWPC", "IWPCV", "IWPC SVR", "IWPC MARS", "IWPC BART", "NLM", "SVR","MARS", "BART", "CLINICAL"))) %>%   filter(name != "CLINICAL") %>% 
  as_data_frame()


dat2 = dat %>% 
  group_by(dosegroup, name) %>% 
  summarise(prop = mean(prop),
            MAE = mean(MAE),
            lower = mean(lower),
            upper = mean(upper)) %>% 
  mutate_if(is.numeric,funs(round(., digits = 2))) %>% 
  unite("CI", lower:upper , sep = "-") %>% 
  unite("MAE (95% CI)", MAE:CI, sep = "(")%>% 
  mutate(`MAE (95% CI)`= gsub("\\(", " (", `MAE (95% CI)`),
         `MAE (95% CI)` = paste(`MAE (95% CI)`,")", sep = ""))

#write_csv(dat2, "../dats_ULLAbydosegroup.csv" )

hlinedat = dat %>% 
  # group_by(dosegroup) %>% 
  summarise(prop = median(prop))

group_n <- c(
  'High'="High (n = 189)",
  'Intermediate'="Intermediate (n = 1,128)",
  'Low'="Low (n = 417)"
)


p3 = ggplot(dat, aes(name, prop, color = name)) +
  geom_hline(data = hlinedat, aes(yintercept=prop),
             color = "gray30") +
  geom_boxplot(width=.25, alpha = .2) + 
  geom_jitter(width = .2, alpha = .1) + 
  theme_minimal()+
  theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.y = element_text(size = 20),
        # axis.title.y = element_text(size = 22),
        # legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 5),
        # legend.title = element_text(size = 22),

        legend.position = "none",
        # legend.direction = 'horizontal'
        
        ) + 
  labs(x = "", y = " ", color = " ") +
  facet_wrap(.~dosegroup, nrow = 1, labeller = as_labeller(group_n))+
  scale_color_manual(values = c("#F8766D" ,"#D39200", "#93AA00", "#00BA38", "#00C19F" ,"#00B9E3" ,"#619CFF", "#DB72FB", "#FF61C3", "gray30")) + 
  guides(color = guide_legend(ncol = 5)) +
  ylim(9,75)

 #pdf("../fig4ulla_bydosegroup.pdf", width = 8, height = 5)

p3
p3_dosegroup = p3 
#dev.off()

dat3 = dat %>% 
  group_by(dosegroup) %>% 
  summarise(prop = mean(prop)) %>% 
  as_data_frame() %>% 
  mutate(dosegroup = factor(dosegroup))

friedman_test(prop ~ name |replicate, data = dat3)

prop_wilcox = as.data.frame(dat) %>%
  group_by(dosegroup) %>% 
  wilcox_test(prop ~ name, paired = TRUE, p.adjust.method = "bonferroni")

MAE_wilcox = as.data.frame(dat) %>%
  group_by(dosegroup) %>% 
  wilcox_test(MAE ~ name, paired = TRUE, p.adjust.method = "bonferroni")

wilcox_results_ULLA = rbind(prop_wilcox, MAE_wilcox)
#write_csv(wilcox_results_ULLA, "../wilcox_results_ULLA_bydosegroup.csv" )


### put the plots together
library(cowplot)

fig3 = plot_grid(p3_dosegroup,p3_race,p3_country, nrow = 3) 
ggsave("../fig2.tiff", fig3, height = 15, width = 12.7,units = "cm",  bg = "transparent")
export::graph2ppt(fig3, file="../fig3.pptx", width=8.4, height = 3.5, bg = "transparent")
