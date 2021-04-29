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
library(rsq)
library(broom)


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

mergeddat = read_csv("merged_iwpc_ULLA.csv")

###########################################################
#####################DATA##################################
###########################################################

mergeddat = mergeddat %>% 
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


train = mergeddat %>% sample_frac(.7) %>% mutate(train = 1) 
mergeddat$train = if_else(mergeddat$X1 %in% train$X1, 1, 0)
test = mergeddat %>% filter(train == 0)

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
  
  
  
  
  partialr2_iwpc =rsq.partial(fit)$partial.rsq
  summary_iwpc = summary(fit)$coefficients[,c(-3)]
  r2_iwpc = summary(fit)$r.squared

  #### IWPC VARIABLES
  
  iwpcvars = lm(sqrtdose ~ age + height + weight + cyp+vkor  + 
                  race + ei + amio , 
                data = train)
  
  test$pred_iwpcvars = predict(iwpcvars, test)
  
  R2_iwpcvars = 1 - (sum((test$pred_iwpcvars - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2))
  
  
  partialr2_iwpcv =rsq.partial(iwpcvars)$partial.rsq
  summary_iwpcv = summary(iwpcvars)$coefficients[,c(-3)]
  r2_iwpcv = summary(iwpcvars)$r.squared
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
  
  
  newcovars = lm(sqrtdose ~ age + height + weight +sex + cyp + vkor   +race  +ethnicity+ ei  + amio+ + statin + aspirin + indication + diabetes + smoke, data = train)
  
  test$pred_NLM = predict(newcovars, test)
  
  R2_newcovars = 1 - (sum((test$pred_NLM - test$sqrtdose) ^2)/sum((test$sqrtdose - mean(test$sqrtdose)) ^2))
  
  
  partialr2_nlm =rsq.partial(newcovars)$partial.rsq
  summary_nlm = summary(newcovars)$coefficients[,c(-3)]
  r2_nlm = summary(newcovars)$r.squared
  
  

  ################ ML ######################
  
  ### LASSO #### 
  lambdas <- 10^seq(2, -3, by = -.1)
  
  # Setting alpha = 1 implements lasso regression
  x1 = train %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate_if(is.factor, as.numeric) %>% 
    as.matrix()
  
  x1cols = train %>%
    dplyr::select(c( age , height , weight , vkor , cyp , 
                     race , ei ,amio,
                     ethnicity, diabetes, smoke, aspirin ,statin, sex , indication)) %>% 
    colnames()
  
  y1 = train$sqrtdose 
  
  lasso_reg <- cv.glmnet(x1[,x1cols],y1, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  
  # Best 
  lambda_best <- lasso_reg$lambda.min 
  lambda_best
  
  lasso_model <- glmnet(x1[,x1cols],y1, alpha = 1, lambda = lambda_best, standardize = TRUE)
  
  
  x2 = test %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate_if(is.factor, as.numeric) %>% 
    as.matrix()
  
  x2cols = test %>%
    dplyr::select(c( age , height , weight , vkor , cyp , 
                     race , ei ,amio,
                     ethnicity, diabetes, smoke, aspirin ,statin, sex , indication)) %>% 
    colnames()
  
  
  
  test$pred_lasso = predict(lasso_model, s = lambda_best, newx = x2[,x2cols])
  
  lambdas <- 10^seq(2, -3, by = -.1)
  
  # Setting alpha = 1 implements lasso regression
  x1 = train %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate_if(is.factor, as.numeric) %>% 
    as.matrix()
  
  x1cols = train %>%
    dplyr::select(c( age , height , weight , vkor , cyp , 
                     race , ei ,amio,
                     ethnicity, diabetes, smoke, aspirin ,statin, sex , indication)) %>% 
    colnames()
  
  y1 = train$sqrtdose 
  
  lasso_reg <- cv.glmnet(x1[,x1cols],y1, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  
  # Best 
  lambda_best <- lasso_reg$lambda.min 
  lambda_best
  
  lasso_model <- glmnet(x1[,x1cols],y1, alpha = 1, lambda = lambda_best, standardize = TRUE)
  
  
  x2 = test %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate_if(is.factor, as.numeric) %>% 
    as.matrix()
  
  x2cols = test %>%
    dplyr::select(c( age , height , weight , vkor , cyp , 
                     race , ei ,amio,
                     ethnicity, diabetes, smoke, aspirin ,statin, sex , indication)) %>% 
    colnames()
  
  
  
  test$pred_lasso = predict(lasso_model, s = lambda_best, newx = x2[,x2cols])
  
  #### random forest ##### 
  rf <- randomForest(
    sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei  + ethnicity +sex + smoke  + statin+ aspirin+ diabetes,data = train
  )
  
  test$pred_rf=stats::predict(rf,test)
  
  #### SVM ##### 
  
  OptModelsvm = svm(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei  + ethnicity +sex + smoke  + statin+ aspirin+ diabetes,data = train)
  
  
  test$pred_svm =stats::predict(OptModelsvm,test)
  
  
  R2_svm = 1 - (sum((test$pred_svm^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## BART ####
  ## Bayesian Additive Regression Trees (BART)
  
  x1 = train %>%
    dplyr::select(c(vkor,cyp ,age  , weight   , amio  ,indication , height ,race  , ei  , smoke  ,statin,aspirin, diabetes, ethnicity, sex)) %>%
    as.data.frame()
  
  bart <- bartMachine(x1, train$sqrtdose, mem_cache_for_speed = F)
  
  
  testx = test %>% 
    dplyr::select(c(vkor,cyp ,age  , weight   , amio  ,indication , height ,race  , ei  , smoke  ,statin,aspirin, diabetes, ethnicity, sex)) %>%
    as.data.frame()
  
  test$pred_bart = predict(bart, testx)
  
  R2_bart = 1 - (sum((test$pred_bart^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ## MARS #####
  # Fit a basic MARS model
  
  
  mars = train(sqrtdose ~ vkor + cyp + age  + weight   + amio  + indication + height +race  + ei  + ethnicity +sex + smoke  + statin+ aspirin+ diabetes,data = train,
               method = "earth",
               metric = "MAE"
  )
  
  test$pred_mars = predict(mars, test)
  
  
  
  R2_mars = 1 - (sum((test$pred_mars^2 - test$dosewk) ^2)/sum((test$dosewk - mean(test$dosewk)) ^2))
  
  ################ results #################
  #########################################

  
  mae = cbind(iwpc = mae2(test$dosewk, test$pred.iwpcmodel)$mae,
              iwpcvars = mae2(test$dosewk, test$pred_iwpcvars)$mae,
              svm_iwpc = mae2(test$dosewk, test$pred_svm_iwpc)$mae,
              bart_iwpc = mae2(test$dosewk, test$pred_bart_iwpc)$mae,
              mars_iwpc = mae2(test$dosewk, test$pred_mars_iwpc)$mae,
              newlinear = mae2(test$dosewk, test$pred_NLM)$mae,
              lasso = mae2(test$dosewk, test$pred_lasso)$mae,
              rfr = mae2(test$dosewk, test$pred_rf)$mae,
              svm = mae2(test$dosewk, test$pred_svm)$mae,
              mars = mae2(test$dosewk, test$pred_mars)$mae,
              bart = mae2(test$dosewk, test$pred_bart)$mae)
  
  
  
  mae_ci = cbind(iwpc = mae2(test$dosewk, test$pred.iwpcmodel)$CI,
                 iwpcvars = mae2(test$dosewk, test$pred_iwpcvars)$CI,
                 svm_iwpc =  mae2(test$dosewk, test$pred_svm_iwpc)$CI,
                 bart_iwpc =  mae2(test$dosewk, test$pred_bart_iwpc)$CI,
                 mars_iwpc =  mae2(test$dosewk, test$pred_mars_iwpc)$CI,
                 newlinear = mae2(test$dosewk, test$pred_NLM)$CI,
                 lasso = mae2(test$dosewk, test$pred_lasso)$CI,
                 rfr = mae2(test$dosewk, test$pred_rf)$CI,
                 svm = mae2(test$dosewk, test$pred_svm)$CI,
                 mars = mae2(test$dosewk, test$pred_mars)$CI,
                 bart = mae2(test$dosewk, test$pred_bart)$CI)
  
  
  mae_dat = t(rbind(mae, mae_ci))
  
  
  mae_dat = mae_dat %>% 
    as.data.frame() %>% 
    rownames_to_column("name") %>% 
    dplyr::select(name, MAE = V1, lower, upper) %>% 
    mutate(name = factor(name, levels = c("iwpc", "iwpcvars", "svm_iwpc", "bart_iwpc", "mars_iwpc", "newlinear","lasso", "rfr", "svm", "mars", "bart"), labels = c("IWPC", "IWPCV", "IWPC_SVM", "IWPC_BART", "IWPC_MARS", "NLM", "LASSO","RFR", "SVM", "MARS", "BART")))
  
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
    mutate(name = factor(name, levels = c("iwpcmodel", "iwpcvars","svm_iwpc", "bart_iwpc", "mars_iwpc" ,"NLM", "lasso","rf", "svm", "bart", "mars"), labels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_BART", "IWPC_MARS", "NLM", "LASSO","RFR", "SVM", "BART", "MARS"))) %>% 
    arrange(name) %>% 
    # inner_join(mae_ethnicity_dat, by = c("ethnicity","name" )) %>% 
    inner_join(mae_dat, by = c("name")) %>% 
    mutate_if(is.numeric,round,2) 
  
  partialr2 = list(iwpc = partialr2_iwpc, iwpcv = partialr2_iwpcv, nlm = partialr2_nlm)
  
  summary2 = list(iwpc = summary_iwpc, iwpcv = summary_iwpcv, nlm= summary_nlm)
  
  r2 = list(iwpc = r2_iwpc, iwpcv = r2_iwpcv, nlm = r2_nlm)
  
  attr(resultsdat2, "summary2") <- summary2
  attr(resultsdat2, "partialr2") <- partialr2
  attr(resultsdat2, "r2") <- r2
  return(resultsdat2)
  
}



###############################################
#########FUNCTION FOR MULTIPLE RUNS###########
##############################################
fullReplicates <- lapply(1:100, function(x){
  cat('\rreplicate=',x)
  k<-fitOneSet()
  list("summary2"=attr(k, "summary2"),
       "partialr2"=attr(k, "partialr2"),
       "r2" = attr(k, "r2"),
       'stats'=cbind.data.frame('replicate'=x, k))
})

dat = do.call(rbind, lapply(1:100, function(x) fullReplicates[[x]][[4]] )) %>% 
  group_by(name) %>% 
  mutate(name = factor(name, levels = c("IWPC", "IWPCV","IWPC_SVM", "IWPC_MARS", "IWPC_BART", "NLM","LASSO","RFR", "SVM", "MARS", "BART"),
                       labels = c("IWPC", "IWPCV", "IWPC SVR", "IWPC MARS", "IWPC BART", "NLM","LASSO","RFR", "SVR","MARS", "BART"))) %>% 
  as_data_frame()


partials = do.call(rbind, lapply(50, function(x) unlist(fullReplicates[[x]][[2]] ))) %>% as_data_frame()
partials2 = round(partials, digits = 2)

r2s = do.call(rbind,lapply(50, function(x) unlist(fullReplicates[[x]][[3]]))) %>% as.data.frame()

betas = do.call(rbind, lapply(50, function(x) do.call(rbind.data.frame,fullReplicates[[x]][[1]]))) %>% rownames_to_column()
betas$Estimate = round(betas$Estimate, digits = 2)
betas$`Std. Error` = round(betas$`Std. Error`, digits =2 )
betas$call = paste(betas$Estimate,"±",betas$`Std. Error`)
betas$call = str_replace_all(string=betas$call, pattern=" ", repl="")
betas$`Pr(>|t|)` = signif(betas$`Pr(>|t|)`, digits = 3)
betas$`Pr(>|t|)` = gsub("e", "x10", betas$`Pr(>|t|)`)
betas2 =betas %>% 
  dplyr::select(rowname, call, p = `Pr(>|t|)`)
# write_csv(partials2, "../partials_merged.csv" )
# write_csv(r2s, "../r2s_merged.csv" )
# write_csv(betas2, "../betas_merged.csv" )

dat2 = dat %>% 
  filter(name %in% c("IWPC", "IWPCV", "IWPC SVR", "IWPC MARS", "IWPC BART", "NLM","LASSO", "RFR","SVR","MARS", "BART")) %>% 
  group_by(name)

dats = dat2 %>% 
  group_by(name) %>% 
  summarise(prop = median(prop),
            MAE = median(MAE),
            lower = median(lower),
            upper = median(upper)) %>% 
  unite("CI", lower:upper , sep = "-") %>% 
  unite("MAE (95% CI)", MAE:CI, sep = "(")


# write_csv(dats, "../dats_merged.csv" )

p2 = ggplot(dat2, aes(name, prop, color = name)) +
  geom_boxplot(width=.1, alpha = .1) + 
  geom_jitter(width = .2, alpha = .1) + 
  theme_minimal()+
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 20, color = "gray30"),
        legend.title = element_text(size = 22, color = "gray30"))+
  labs(x = "", y = " ", color = "Model")+
  ylim(c(40,54))

## yellow, orange, navy, light orange, gray
# pdf("fig1merged.pdf", width = 8, height = 5)
p2

# dev.off()

friedman_test(prop ~ name |replicate, data = as.data.frame(dat))

prop_wilcox = as.data.frame(dat) %>%
  wilcox_test(prop ~ name, paired = TRUE, p.adjust.method = "bonferroni")

MAE_wilcox = as.data.frame(dat) %>%
  wilcox_test(MAE ~ name, paired = TRUE, p.adjust.method = "bonferroni")

# p5 = ggplot(dat, aes(name, MAE, color = name)) +
#   geom_boxplot(alpha =.1) +
#   geom_jitter(width = .2, alpha = .2) + 
#   theme_minimal()+
#   theme(legend.position = "none")+
#   labs(x = "", y = " ") +
#   ylim(7,11.5)
# p5

 wilcox_results_Merged = rbind(prop_wilcox, MAE_wilcox)
# write_csv(wilcox_results_Merged, "../wilcox_results_Merged.csv" )

 # save.image("~/Documents/IWPC/04_workspace.RData")
 
### put the plots together
library(cowplot)

fig1 = plot_grid(p3,p2, nrow = 1,rel_widths = 1:1, labels = c("A. ULLA (n = 1,591)", "B. Merged (n = 6,526)"), label_size = 20) 
save_plot("../fig1.png", fig1, base_width = 15, base_height = 8, bg = "transparent")


# fig2 = plot_grid(p4,p5,p6, nrow = 1, labels = c("A. IWPC", "B. Merged", "C. ULLA")) 
# save_plot("figure2.png", fig2, base_width = 15, base_height = 8, bg = "transparent")
