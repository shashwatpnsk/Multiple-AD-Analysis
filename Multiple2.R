# 000. Preliminaries -----

library(tidyverse)
library(hablar)
library(lfe)
library(readxl)
library(lubridate)
library(olsrr)
library(blorr)

# library(pastecs)
# library(rstatix)
# library(nortest)

#Functions to be used
###A.
logistic.regression.or.ci <- function(regress.out, level=0.95){
  ################################################################
  #                                                              #
  #  This function takes the output from a glm                   #
  #  (logistic model) command in R and provides not              #
  #  only the usual output from the summary command, but         #
  #  adds confidence intervals for all coefficients and OR's.    #
  #                                                              #
  #  This version accommodates multiple regression parameters    #
  #                                                              #
  ################################################################
  usual.output <- summary(regress.out)
  z.quantile <- qnorm(1-(1-level)/2)
  number.vars <- length(regress.out$coefficients)
  OR <- exp(regress.out$coefficients[-1])
  temp.store.result <- matrix(rep(NA, number.vars*2), nrow=number.vars)
  for(i in 1:number.vars)
  {
    temp.store.result[i,] <- summary(regress.out)$coefficients[i] +
      c(-1, 1) * z.quantile * summary(regress.out)$coefficients[i+number.vars]
  }
  intercept.ci <- temp.store.result[1,]
  slopes.ci <- temp.store.result[-1,]
  OR.ci <- exp(slopes.ci)
  output <- list(regression.table = usual.output, intercept.ci = intercept.ci,
                 slopes.ci = slopes.ci, OR=OR, OR.ci = OR.ci)
  return(output)
}

### B.
#Using logit to probability converter
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Multiple
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# I. Load the data -----

setwd("C:/Users/mailm/OneDrive/URW/Main Projects/Nandurbar CMAM/Multiple")

store <- read_csv("final_data.csv")

store <- store %>%
  convert(
    fct(admission_criterion),
    fct(treatment_group),
    fct(gender),
    fct(outcome),
    fct(block),
    fct(category)
  )
store$block[store$block == "nandurbar"] <- "Nandurbar"
store$block[store$block == "taloda"] <- "Taloda"

store$PHC[store$PHC == "Kusumsada"] <- "Kusumwada"
store$PHC[store$PHC == "pratappur (Taloda)"] <- "Pratappur (Taloda)"
store$PHC[store$PHC == "somawal"] <- "Somawal"


# II. Prepare for Analysis ------

## A. Variables that will be analysed ----
data <- store%>%
  # Calculate recovery time amongst those recovered
  mutate(rectime = ifelse(outcome == "cured", 
                          (date_of_discharge - date_of_admission) %>% as.numeric(), 
                          NA
                          )
         )%>%
  mutate(admittime = (date_of_discharge - date_of_admission) %>% as.numeric()
  )%>%
  # Create dummy for cured
  mutate(dcured = ifelse(outcome == 'cured', 1, 0)%>%
           as.factor()
         ) %>%
  #Dummy for treated or not
  mutate(
      dtreat = ifelse(
        outcome == 'cured' |
        outcome == 'not_recovered' |
        outcome == 'died',
        as.numeric(1),
        as.numeric(0)
      )%>%
       as.factor()
    )%>%
  #Calculate wt gain
  mutate(
    wtgain_abs = (discharge_weight - admission_weight)*1000
  )%>%
  mutate(
    wtgain_avg = (discharge_weight - admission_weight)*1000/(admission_weight*rectime)
  )%>%
  mutate(
    wtgain_avg_all = (discharge_weight - admission_weight)*1000/(admission_weight*admittime)
  )

## B. Create dummies ----

data <- data %>%
  mutate(dwfh = ifelse(admission_criterion == "WFH", 1, 0)%>%
           as.factor())%>%
  mutate(dwasting = ifelse(admission_wfh < -3, 1, 0)%>%
           as.factor())%>%
  mutate(dstunting = ifelse(admission_lfa < -3, 1, 0)%>%
           as.factor())%>%
  mutate(dunderweight = ifelse(admission_wfa < -3, 1, 0)%>%
           as.factor())%>%
  mutate(dage = ifelse(age_m >= 24, 1, 0)%>%
           as.factor())%>%
  mutate(dstunting_discharge = ifelse( discharge_lfa < -3, 1, 0)%>%
           as.factor())%>%
  mutate(dstunting_cure = ifelse(dstunting_discharge == dstunting, 0,1))%>%
  mutate(ismnt = ifelse(treatment_group == "MNT", 0,1))%>%
  mutate(isarf = ifelse(treatment_group == "ARF", 0,1))%>%
  mutate(issf = ifelse(treatment_group == "SF", 0,1))
  
## C. Age wise Recovery Days ----

result_days <- tibble(age = 0, N = 0, stunted = 0, non_stunted =  0)

for(i in 1:73) {
  
  result_days[i,1] <- i
  result_days[i,2] <- data %>% 
                        filter(dwfh == 1) %>% 
                        filter(is.na(rectime) == FALSE) %>%
                        filter(age_m == i) %>% 
                        nrow()
  result_days[i,3] <- NA
  result_days[i,4] <- NA
  
  tryCatch({
  reg <- data%>%
            filter(dwfh == 1)%>%
            filter(age_m == i)%>%
          lm(rectime ~ dstunting + dunderweight, data = .)

  result_days[i,3] <- coef(reg)[1] + coef(reg)[2]
  result_days[i,4] <- coef(reg)[1]
  } , error = function(e){})
}

result_days$non_stunted[3] <- NA

result_days$gap <- result_days$stunted - result_days$non_stunted

result_days_longer <- result_days %>% pivot_longer(names_to = "stunting", values_to = "recovery_days", cols = stunted:non_stunted)

## D. Age wise Recovery Rate ----

 result_rate <- data%>%
                  filter(dwfh == 1)%>%
                  filter(dtreat == 1)%>%
                  convert(num(dcured))%>%
                  group_by(age_m, dstunting)%>%
                  summarise(rec_rate = 100*mean(dcured, na.rm = TRUE))%>%
                  pivot_wider(names_from = "dstunting", values_from = "rec_rate")%>%
                  select(-`NA`)%>%
                  rename("stunted" = "1", "only_wasted" = "0")%>%
                  mutate(gap_stun = stunted - only_wasted)
    
    
    
# result_rate$gap_stun <- result_rate$stunted - result_rate$only_wasted
# 
# result_rate$gap_underweight <- result_rate$underweight - result_rate$only_wasted
# 
# result_rate$gap_stun_suw <- result_rate$stun_suw - result_rate$only_wasted
# 
result_rate_longer <- result_rate %>% pivot_longer(names_to = "category", values_to = "recovery_rate", cols = only_wasted:stunted)
# 
result_rate_longer_gap <- result_rate_longer %>% pivot_longer(names_to = "gap", values_to = "recovery_rate_gap", cols = gap_stun)


# III. Plots ----

#1. Recovery Days -----

## A. Stunting vs Non Stunted

  data%>%
  hablar::convert(fct(dstunting))%>%
  mutate(dstunting = ifelse(dstunting == 1, "Stunted", "Not Stunted"))%>%
  #   filter(dadmit == 1)%>%
  ggplot(aes(y=rectime, x=age_m, color=dstunting))+
  #   geom_point()+
  geom_smooth()+
  xlab("Age (Months)")+
  ylab("Recovery Time (days)")+
  theme_bw()+
  labs(color = "Stunting")

## A.2.

#G.  Stunting vs Non-Stunted gap
  data%>%
    hablar::convert(fct(dstunting))%>%
    mutate(dstunting = ifelse(dstunting == 1, "Stunted", "Not Stunted"))%>%
    #   filter(dadmit == 1)%>%
    ggplot(aes(y=rectime, x=age_m, color=dstunting))+
    #   geom_point()+
    geom_smooth(method = "lm")+
    xlab("Age (Months)")+
    ylab("Recovery Time (days)")+
    theme_bw()+
    labs(color = "Stunting")


## B. Stunted vs Non Stunted gap
  
    result_days%>%
      filter(is.na(gap) == FALSE)%>%
      ggplot(aes(y=gap, x=age, size = N))+
        geom_smooth(method = "lm")+
      #  geom_point(colour = "purple")+
        xlab("Age (Months)")+
        ylab("Recovery Time Gap (days)")+
      theme_bw()
    
## C. All Non-Stunted
        
    result_days%>%
          filter(is.na(non_stunted) == FALSE)%>%
          ggplot(aes(y=non_stunted, x=age, size = N))+
          geom_smooth(method = "lm")+
          geom_point(colour = "purple")+
          xlab("Age (Months)")+
          ylab("Recovery Time Non Stunted Kids (days)")+
          guides(size = guide_legend(title = "No. of \nChildren"))+
          scale_size(range = c(2,4.5))+
          theme_bw()
    
        
## D. All Stunted
        
    result_days%>%
          filter(is.na(stunted) == FALSE)%>%
          ggplot(aes(y=stunted, x=age, size = N))+
          geom_smooth(method = "lm")+
          geom_point(colour = "purple")+
          xlab("Age (Months)")+
          ylab("Recovery Time Stunted Kids (days)")+
          guides(size = guide_legend(title = "No. of \nChildren"))+
          theme_bw()
        labs(color = "Stunted?")
        
## E. Together
        
        result_days_longer%>%
          filter(is.na(recovery_days) == FALSE)%>%
          ggplot(aes(y=recovery_days, x=age, group=stunting, color = stunting, size = N))+
          geom_smooth(method = "lm", se = FALSE, size = 1.3)+
          geom_point(alpha = 0.4)+
          xlab("Age (Months)")+
          ylab("Recovery Time (days)")+
          labs(size = "No. of \nChildren")+
          labs(color = "Stunting")+
        scale_size(range = c(1,4))
        
## F. Relation to WFH Z-Scores
        
        data %>%
          filter(dwfh == 1)%>%
          filter(admission_wfh <= -3)%>%
          ggplot(aes(x=admission_wfh, y= rectime))+
          geom_point()+
          geom_smooth(method = "lm")

## G. Sex and Recovery Time
        
        data%>%
          hablar::convert(fct(dstunting))%>%
          mutate(dstunting = ifelse(dstunting == 1, "Stunted", "Not Stunted"))%>%
          #   filter(dadmit == 1)%>%
          ggplot(aes(y=rectime, x=age_m, color=gender))+
          #   geom_point()+
          geom_smooth(method = "lm")+
          xlab("Age (Months)")+
          ylab("Recovery Time (days)")+
          theme_bw()+
          labs(color = "Sex")
        
                  
#2. Recovery Rate (VERY EXPERIMENTAL) ------      
        
        ## A. None
        
        ## B. Stunted vs Non Stunted gap
        
        result_rate%>%
          filter(is.na(gap_stun) == FALSE)%>%
          ggplot(aes(y=gap_stun, x=age_m))+
          geom_smooth(method = "lm")+
          geom_point(colour = "purple")+
          xlab("Age (Months)")+
          ylab("Recovery Rate Gap \n for Stunded Kids (%)")+
          scale_size(range = c(2,4.5))+
          theme_bw()
        
        
        
        # ## C. SUW vs Only wasted gap
        # 
        # result_rate%>%
        #   filter(is.na(gap_underweight) == FALSE)%>%
        #   ggplot(aes(y=gap_underweight, x=age, size = N))+
        #   geom_smooth(method = "lm")+
        #   geom_point(colour = "purple")+
        #   xlab("Age (Months)")+
        #   ylab("Recovery Rate Gap \n for Underweight Kids (%)")+
        #   scale_size(range = c(2,4.5))
        # 
        # ## D. SUW + Stunted vs Only wasted gap
        # 
        # result_rate%>%
        #   filter(is.na(gap_underweight) == FALSE)%>%
        #   ggplot(aes(y=gap_underweight, x=age, size = N))+
        #   geom_smooth(method = "lm")+
        #   geom_point(colour = "purple")+
        #   xlab("Age (Months)")+
        #   ylab("Recovery Rate Gap \n for Underweight + Stunted Kids (%)")+
        #   scale_size(range = c(2,4.5))        
        
        
        ## E. All Stunted
        
        result_rate%>%
          filter(is.na(stunted) == FALSE)%>%
          ggplot(aes(y=stunted, x=age_m))+
          geom_smooth(method = "lm")+
          geom_point(colour = "purple")+
          xlab("Age (Months)")+
          ylab("Recovery Rate Stunted Kids (%)")+
          scale_size(range = c(2,4.5))+
          theme_bw()

        # ## E. All Underweight
        # 
        # result_rate%>%
        #   filter(is.na(underweight) == FALSE)%>%
        #   ggplot(aes(y=underweight, x=age, size = N))+
        #   geom_smooth(method = "lm")+
        #   geom_point(colour = "purple")+
        #   xlab("Age (Months)")+
        #   ylab("Recovery Rate \n Underweight Kids (%)")+
        #   scale_size(range = c(2,4.5))+
        #   theme_bw()
        # 
        # ## F. All Underweight + Stunted
        # 
        # result_rate%>%
        #   filter(is.na(stun_suw) == FALSE)%>%
        #   ggplot(aes(y=stun_suw, x=age, size = N))+
        #   geom_smooth(method = "lm")+
        #   geom_point(colour = "purple")+
        #   xlab("Age (Months)")+
        #   ylab("Recovery Rate Underweight Kids (%)")+
        #   scale_size(range = c(2,4.5))+
        #   theme_bw()
        # 
               
        ## E. Together All Category
        
        result_rate_longer%>%
          filter(is.na(recovery_rate) == FALSE)%>%
          ggplot(aes(y=recovery_rate, x=age_m, group=category, color = category))+
          geom_smooth(method = "lm", se = FALSE)+
          geom_point(alpha = 0.3)+
          xlab("Age (Months)")+
          ylab("Recovery Rate (%)")+
          scale_size(range = c(1,4))+
          theme_bw()
        

#3. Weight Gain

# A. Total        
               
# a. cured 
 data%>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      ggplot()+
      geom_density(aes(x = wtgain_avg, colour = dcured, fill = dcured), alpha = 0.25)+
      xlim(-2,3)+
      theme_bw()

#b. treatment    
    
    data%>%
    filter(dwfh == 1)%>%
      ggplot()+
      geom_density(aes(x = wtgain, colour = dtreat, fill = dtreat), alpha = 0.25)+
      xlim(-2,3)    
    
#c. admission criterion    
    
    data%>%
      ggplot()+
      geom_density(aes(x = wtgain, colour = admission_criterion, fill = admission_criterion), alpha = 0.25)+
      xlim(-2,3)    

#d. Stunting    
    
    data%>%
      filter(dtreat == 1,  dwfh == 1, is.na(dstunting) == FALSE)%>%
      ggplot()+
      geom_density(aes(x = wtgain, colour = dstunting, fill = dstunting), alpha = 0.25)+
      xlim(-2.5,2.5)        

#B. Velocity    
            
    #a. Stunting
    
    data%>%
      filter(dtreat == 1,  dwfh == 1, dcured == 1, is.na(dstunting) == FALSE)%>%
      ggplot()+
      geom_density(aes(x = wtgain_vel, colour = dstunting, fill = dstunting), alpha = 0.25)

#C. Discharge vs Admission
    x
    data%>%
      filter(dtreat == 1,  dwfh == 1, is.na(dstunting) == FALSE)%>%
      select(discharge_weight, admission_weight, dstunting)%>%
      hablar::convert(fct(dstunting))%>%
      pivot_longer(-dstunting, names_to = "Measured", values_to = "Weight")%>%
      ggplot()+
      geom_density(aes(x = Weight, colour = Measured, fill = dstunting), alpha = 0.55)+
      scale_color_grey()
      scale_fill_()

      
# IV. Tests on Data ----
        
#A Multi collinearity ----
       
#1. Block
  
#a. On WFH      
regA1a <- data %>%
          filter(dwfh == 1)%>%
          lm(admission_wfh ~ block, data = .)
        
summary(rega1)
        
#b. On LFA     
regA1b <- data %>%
  filter(dwfh == 1)%>%
  lm(admission_lfa ~ block, data = .)

summary(regA1b)

#c. On WFA      
regA1c <- data %>%
  filter(dwfh == 1)%>%
  lm(admission_wfa ~ block, data = .)

summary(regA1c)
       
#2. PHC

#a. On WFH      
regA2a <- data %>%
  filter(dwfh == 1)%>%
  lm(admission_wfh ~ PHC, data = .)

summary(regA2a)

#b. On LFA     
regA2b <- data %>%
  filter(dwfh == 1)%>%
  lm(admission_lfa ~ PHC, data = .)

summary(regA2b)

#c. On WFA      
regA2c <- data %>%
  filter(dwfh == 1)%>%
  lm(admission_wfa ~ PHC, data = .)

summary(regA2c)
 
#B. A check on whether being stunted is correlated with things ----

# 1. Taking up treatment

regB1 <- data %>%
  filter(dwfh == 1)%>%
  glm(dtreat ~ dstunting, family = "binomial", data = .)

summary(regB1)

#2. By treatment arm, gender and age

regB2 <- data%>%
  filter(dwfh ==1)%>%
  glm(dstunting ~ treatment_group + gender + age_m, family ="binomial", data = .)

summary(regB2)
# More stunted kids in SF and Males

#3.
table(data$dstunting, data$gender)


#C. A check on whether being underweight is correlated with things ----

# 1. Taking up treatment

regC1 <- data %>%
  filter(dwfh == 1)%>%
  glm(dtreat ~ dunderweight, family = "binomial", data = .)

summary(regC1)

#2. By treatment arm, gender and age

regC2 <- data%>%
  filter(dwfh ==1)%>%
  glm(dunderweight ~ treatment_group + gender + age_m, family ="binomial", data = .)

summary(regC2)

# Uncontrolled Group Means

#a. Recovery Rate

table(temp$dunderweight, temp$dcured)

table(temp$dstunting, temp$dcured)

table(temp$gender, temp$dcured)

#b. Recovery Times

data %>% 
  filter(dwfh == 1) %>% 
  select(rectime, dunderweight) %>% 
  group_by(dunderweight) %>%
  summarise(percentUW = mean(rectime, na.rm= TRUE))

data %>% 
  filter(dwfh == 1) %>% 
  select(rectime, dstunting) %>% 
  group_by(dstunting) %>%
  summarise(percentSt = mean(rectime, na.rm= TRUE))

data %>% 
  filter(dwfh == 1) %>% 
  select(rectime, gender) %>% 
  group_by(gender) %>%
  summarise(percentGen = mean(rectime, na.rm= TRUE))


temp <- data %>% filter(dwfh == 1)

table(temp$dunderweight, temp$dcured)

regaC1 <- data %>%
  filter(dwfh == 1)%>%
  lm(rectime ~ dstunting, data = .)

summary(regaC1)

#2. Underweight

regaC2 <- data %>%
  filter(dwfh == 1)%>%
  lm(rectime ~ dunderweight, data = .)

summary(regaC2)

data %>% group_by(dunderweight) %>% summary()

# D. Counts -----

#1. SUW vs Stunted numbers


limits <- c(-Inf, -3 , -2, Inf)

counts <- matrix(nrow = 3, ncol = 3, dimnames = (list(c("SUW", "MUW", "NUW") , c("SS","MS", "NS"))))

for (i in 1:3) {
  
  uw_u_bound <- limits[i+1]
  uw_l_bound <- limits[i]
  
  
  for(j in 1:3) {
    
    s_u_bound <- limits[j+1]
    s_l_bound <- limits[j]
    
    
    count <- data%>%
      filter(dwfh == 1)%>%
      filter( uw_l_bound < admission_wfa & admission_wfa <= uw_u_bound )%>%
      filter( s_l_bound < admission_lfa & admission_lfa <= s_u_bound )%>%
      nrow()
    
    counts[i,j] <- count 
    
    
  }
  
}

print(counts)


#2. Counts of gender across stunting


g_vals <- c("F", "M")
s_vals <- c(0,1)

counts <- matrix(nrow = 3, ncol = 5, dimnames = (list(c("SS","NS", "Total"), c("F", "M","Total","perc_F", "perc_M"))))

for (i in 1:2) {
  

  for(j in 1:2) {
    

    
    count <- data%>%
      filter(dwfh == 1)%>%
      filter( dstunting == s_vals[i] & gender == g_vals[j])%>%
      nrow()
    
    counts[i,j] <- count 
    
    
  }
  
}

counts[3,] <- counts[1,] + counts[2,]

counts[,3] <- counts[,1] + counts[,2]

sum(counts[1:2,], na.rm = TRUE)

counts[,4] <- counts[,1]/counts[,3]*100

counts[,5] <- counts[,2]/counts[,3]*100


#3. Counts of Age Categories across Stunting


age_vals <- c(0,1)
s_vals <- c(0,1)

counts <- matrix(nrow = 3, ncol = 5, dimnames = (list(c("SS","NS", "Total"), c("Below2", "Above2","Total","perc_B", "perc_A"))))

for (i in 1:2) {
  
  
  for(j in 1:2) {
    
    
    
    count <- data%>%
      filter(dwfh == 1)%>%
      filter( dstunting == s_vals[i] & dage == age_vals[j])%>%
      nrow()
    
    counts[i,j] <- count 
    
    
  }
  
}

counts[3,] <- counts[1,] + counts[2,]

counts[,3] <- counts[,1] + counts[,2]

counts[,4] <- counts[,1]/counts[,3]*100

counts[,5] <- counts[,2]/counts[,3]*100


print(counts)


#4. Weight gain across age groups






##V. Regressions-----

#A. Recovery days ------

#0. Rectime on WFH, LFA, and WFA scores. Controling for treatment arm

#a. WFH

rega0a <- data %>%
  filter(dwfh == 1)%>%
  filter(admission_wfh <= -3)%>% # There are some data that are high wfh and are distorting the results
  lm(rectime ~ admission_wfh + treatment_group, data = .)

summary(rega0a)
# Negative ** effect
x
#b. WFA

rega0b <- data %>%
  filter(dwfh == 1)%>%
  filter(admission_wfh <= -3)%>% # There are some data that are high wfh and are distorting the results
  lm(rectime ~ admission_wfa + treatment_group, data = .)

summary(rega0b)
# No effect
#a. HFA

rega0c <- data %>%
  filter(dwfh == 1)%>%
  filter(admission_wfh <= -3)%>% # There are some data that are high wfh and are distorting the results
  lm(rectime ~ admission_lfa + treatment_group, data = .)

summary(rega0c)

#Positive *** effect

#~~~~~~# Children with SUW or SS #~~~~~~#
        
#1. All children
  rega1 <- data %>%
    filter(dwfh == 1)%>%
    lm(rectime ~ dstunting + dunderweight, data = .)

  summary(rega1)
  
#2. Children less than 2 years
  rega2 <- data %>%
    filter(dwfh == 1)%>%
    filter(age_m < 24)%>%
    lm(rectime ~ dstunting + dunderweight, data = .)
  
  summary(rega2)
  
#3. Children more than 2 years
  rega3 <- data %>%
    filter(dwfh == 1)%>%
    filter(age_m >= 24)%>%
    lm(rectime ~ dstunting + dunderweight, data = .)
  
  summary(rega3)
    
#4a. Difference-in-Difference estimator (Only Stunted)
  
  rega4a <- data %>%
    filter(dwfh == 1)%>%
    lm(rectime ~ dstunting*dage, data = .)
  
  summary(rega4a)
  
# 5a. Recovery days across ages (DiD)
  
  rega5a <- data %>%
    filter(dwfh == 1)%>%
    lm(rectime ~ dstunting*age_m + dunderweight*age_m, data = .)
  
  summary(rega5a)
  
# Removing underweight as it is anyway not affecting rectime but biasing and reducing p value of dstunting
  
  #4b. Difference-in-Difference estimator
  
  rega4b <- data %>%
    filter(dwfh == 1)%>%
    lm(rectime ~ dstunting*dage, data = .)
  
  summary(rega4b)
  
  # 5b. Recovery days across ages
  
  rega5b <- data %>%
    filter(dwfh == 1)%>%
    lm(rectime ~ dstunting*age_m, data = .)
  
  summary(rega5b)

  # This tells us that the difference in recovery time is actually independent of age

  #5.II Gender controls
  
  rega5.II <- data %>%
    filter(dwfh == 1)%>%
    lm(rectime ~ dstunting + gender, data = .)
  
  summary(rega5.II)
  
  #5.II Gender DiD
  
  rega5.II <- data %>%
    filter(dwfh == 1)%>%
    lm(rectime ~ dstunting*gender, data = .)
  
    summary(rega5.II)
  
  
#6. With Controls, IV and FE

#Note LFE:FELM formula : Dependent ~ Independent + Control | Fixed Effects | IV | SE Clustering. 0 if that particular feature is not used.
  
  #0. Base - Nothing else
  rega10 <- data %>%
    filter(dwfh == 1)%>%
    felm(rectime ~ dstunting + dunderweight, data = .)
  
  summary(rega10)
    
  #A. Age controls
  rega6A <- data %>%
    filter(dwfh == 1)%>%
    felm(rectime ~ dstunting + dunderweight + age_m |0|0|0, data = .)
  
  summary(rega6A)
  
  # Age changes effect 
  
  #B. Age controls, Treat Control
  rega6B <- data %>%
    filter(dwfh == 1)%>%
    felm(rectime ~ dstunting + dunderweight + age_m + treatment_group  | 0 | 0 | 0, data = .)
  
  summary(rega6B)
  #Arm changes effect
  
  #C. Age controls, Treat Controls and PHC FE
  rega6C <- data %>%
    filter(dwfh == 1)%>%
    felm(rectime ~ dstunting + dunderweight + age_m + treatment_group | PHC | 0 |0 , data = .)
  
  summary(rega6C)
  
  #These results don't seem right. It cannot be true that recovery time in MNT is higher than ARF. Accounting for PHCs maybe increasing arbitrariness.
  #Won't use PHC FE
  
  #D. Age controls, Treat Control, Stunting and age dummy
  rega6D <- data %>%
    filter(dwfh == 1)%>%
    felm(rectime ~ dstunting + dunderweight + age_m + treatment_group + dstunting*age_m| 0 | 0 | 0, data = .)
  
  summary(rega6D)

  #E. Age controls, Treat, gender Control
  rega6E <- data %>%
    filter(dwfh == 1)%>%
    felm(rectime ~ dstunting + treatment_group + age_m + gender| 0 | 0 | 0, data = .)
  
  summary(rega6E)
  
    
  ### Hence final chosen models are 6B, 
  
  #7. Step wise 
  
  fitall <- lm(rectime ~ dstunting + treatment_group + age_m*dstunting + gender*dstunting + dunderweight  + age_m*dunderweight + gender*dunderweight, data =  data %>% filter(dwfh == 1))
  formula(fitall)
  
  fitstart <- lm(rectime ~ 1, data = data %>% filter(dwfh == 1))
  
  olsrr::ols_step_all_possible(fitall) %>% plot()
  
  step_fit <- olsrr::ols_step_both_p(fitall, details = TRUE)
  
  step_fit$model
  
#Final model therefore is:
# rectime ~ age_m + dstunting + dunderweight + gender + Treatment group fixed effects

  # 8. Final Model
  rega8 <- data %>%
    filter(dwfh == 1)%>%
    felm(rectime ~ dstunting  + dunderweight + age_m + gender+ treatment_group| 0 | 0 | 0, data = .)
  
  summary(rega8)

#B. Recovery percentage ------
  
  
  #0a. Difference-in-Difference estimator (Only Stunted)
  
  regB0a <- data %>%
    filter(dwfh == 1)%>%
    filter(dtreat == 1)%>%
    glm(rectime ~ dstunting*dage, data = .)
  
  summary(regB0a)
  
  logistic.regression.or.ci(regB0a)
  
  data %>%
    filter(dwfh == 1)%>%
    filter(dtreat == 1)%>%
    filter(dage == 1)%>%
    filter(dstunting == 0)%>%
  #  select(dstunting, dage, dcured)%>%
  #  group_by(dage, dstunting)%>%
    hablar::convert(num(dcured))%>%
    summarise(recR = mean(dcured, na.rm =TRUE))
  
  
  
  
  #Children with SUW or SS
  
  #1. 
  #a. All children ITT
  regB1a <- data %>%
    filter(dwfh == 1)%>%
    convert(num(dcured))%>%
    glm(dcured ~ dstunting + dunderweight, data = .)
  
  summary(regB1a)
  
  #b. All children followed protocol
  regB1b <- data %>%
    filter(dwfh == 1)%>%
    filter(dtreat == 1)%>%
    lm(dcured ~ dstunting + dunderweight, data = .)
  
  summary(regB1b)

  #2. All children below 2
  #a. All children ITT
  regB2a <- data %>%
    filter(age_m < 24)%>%
    filter(dwfh == 1)%>%
    lm(dcured ~ dstunting + dunderweight, data = .)
  
  summary(regB2a)
  
  #b. All children followed protocol
  regB2b <- data %>%
    filter(dwfh == 1)%>%
    filter(age_m < 24)%>%
    filter(dtreat == 1)%>%
    lm(dcured ~ dstunting + dunderweight, data = .)
  
  summary(regB2b)

  #3. All children above 2
    #a. All children ITT
    regB3a <- data %>%
      filter(age_m >= 10)%>%
      filter(dwfh == 1)%>%
      lm(dcured ~ dstunting + dunderweight, data = .)
    
    summary(regB3a)
  
    #b. All children followed protocol
    regB3b <- data %>%
      filter(dwfh == 1)%>%
      filter(age_m >= 24)%>%
      filter(dtreat == 1)%>%
      lm(dcured ~ dstunting + dunderweight, data = .)
    
    summary(regB3b)
    

  #4. Logistic regression on all children who took up treatment
    
    regB4 <- data %>%
      filter(dtreat == 1)%>%
      filter(dwfh == 1)%>%
      glm(dcured ~ dstunting + dunderweight, family = "binomial", data = .)
    
    summary(regB4)    
    
  #5. Logistic regression on all children who took up treatment with controls
    
    regB5 <- data %>%
      filter(dtreat == 1)%>%
      filter(dwfh == 1)%>%
      glm(dcured ~ dstunting + dunderweight + gender + treatment_group + age_m, family = "binomial", data = .)
    
    summary(regB5)    
    
    #5.2. Logistic regression stunting on gender
    
    regB5.2 <- data %>%
      filter(dtreat == 1)%>%
      filter(dwfh == 1)%>%
      glm(dcured ~ dstunting*gender, family = "binomial", data = .)
    
    summary(regB5.2)        
    
    #5.3. Logistic DiD on gender
    
    regB5.3 <- data %>%
      filter(dtreat == 1)%>%
      filter(dwfh == 1)%>%
      glm(dcured ~ dstunting*gender + dunderweight*gender, family = "binomial", data = .)
    
    summary(regB5.3)        
    #5.4. Logistic DiD on age_m
    
    regB5.4 <- data %>%
      filter(dtreat == 1)%>%
      filter(dwfh == 1)%>%
      glm(dcured ~ dstunting*age_m + dunderweight*age_m, family = "binomial", data = .)
    
    summary(regB5.4)        
    
    
    #6. Stepwise logistic regression  
    
    fitall <- glm(formula = dcured ~ dstunting + treatment_group + age_m*dstunting + gender*dstunting + dunderweight  + age_m*dunderweight + gender*dunderweight, 
                   data =  data %>% filter(dwfh == 1, dtreat == 1),
                   family = "binomial")
    formula(fitall)
    
    bestfit <- blorr::blr_step_aic_both(fitall)  
    
    bestfit$model
    
  #Best Fit formula : dcured ~ dstunting + treatment_group + age_m + gender + dunderweight + dstunting:gender
  #Modifiying this in the final regression  
    
  #7. Final regression
    
    regB7 <- data %>%
      filter(dtreat == 1)%>%
      filter(dwfh == 1)%>%
      glm(
          formula = dcured ~ dstunting + dunderweight + treatment_group + age_m + gender,
          family = "binomial", data = .
          )
    
    summary(regB7)    
    
    
  #8. Get the probabilities (INTERPRET VERY CAREFULLY)
    #a. The ORs
    
    regB7.OR.CI <- logistic.regression.or.ci(regB7)
    
    regB7.OR <- regB7.OR.CI$OR
    
    regB7.OR
    
    #b. Percentage Change in probabilities 
    #(One unit increase increases probabilty by X%)
    
    (regB7.OR - 1)*100 
    
    #b. The Probabilities (INTERPRET VERY CAREFULLY)
    
    logit2prob(regB7$coefficients[1])  
    
    #Using tanh (Preferred because it is exact)
    0.5 * (1 + tanh(regB7$coefficients[1])/2)
    
    #9. Running a lm also just to verify
    
    regB9 <- data %>%
      filter(dtreat == 1)%>%
      filter(dwfh == 1)%>%
      hablar::convert(num(dcured))%>%
      lm(
        formula = dcured ~ dstunting + dunderweight + treatment_group + age_m + gender,
        data = .
      )
    
    summary(regB9)    

    
        
    
#     
#   #Numbers
#   
#   filter(data, dstunting == 0) %>% nrow()
#   filter(data, dstunting == 1) %>% nrow()
#   
#   
# #Standard errors
#   dat_wasted <- filter(data, dwasting == 1)
#   dat_wasted_notstunt <- filter(data, dwasting == 1 & dstunting == 0)
#   dat_stunt_wast <- filter(data, dwasting == 1 & dstunting == 1)
#   
#   diff <- filter(data, dwasting == 1)
#   
#   var_stunt <- var(dat_stunt_wast$rectime, na.rm = TRUE)/nrow(dat_stunt_wast)
#   var_wasted <- var(dat_wasted$rectime, na.rm = TRUE)/nrow(dat_wasted)
#   
#   se <- sqrt(var_stunt + var_wasted)
#   
#   
    
#C. Weight Gain Abs ----

#1. Main Regressions
    
    #Note LFE::FELM formula : Dependent ~ Independent + Control | Fixed Effects | IV | SE Clustering. 0 if that particular feature is not used.
    
    #0. Base - Nothing else
    rega10 <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
     #filter(wtgain < 3 & wtgain > -3)%>%
      felm(wtgain ~ dstunting + dunderweight, data = .)
    
    summary(rega10)
    
    #A. Age controls
    rega1A <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
     #filter(wtgain < 3 & wtgain > -3)%>%
      felm(wtgain ~ dstunting + dunderweight + age_m |0|0|0, data = .)
    
    summary(rega1A)
    
    #B. Age controls, Treat Control
    rega1B <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      felm(wtgain ~ dstunting + dunderweight + age_m + treatment_group  | 0 | 0 | 0, data = .)
    
    summary(rega1B)
    
    ## SELECTED MODEL
    
    
    #C. Age controls, Treat Controls and PHC FE
    rega1C <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      felm(wtgain ~ dstunting + dunderweight + age_m + treatment_group | PHC | 0 |0 , data = .)
    
    summary(rega1C)
    
    #Won't use PHC FE. Not used it before also
    
    
#2. DiDs
      
    #A. UW
    
    #a. Age
    
    reg2Aa <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain ~ dunderweight*age_m, data = .)
    
    summary(reg2Aa)
    
    #b. Treatment Group
    
    reg2Ab <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain ~ dunderweight*treatment_group, data = .)
    
    summary(reg2Ab)      
    
    #B. Stunting
    
    #a. Age
    
    reg2Ba <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain ~ dstunting*age_m, data = .)
    
    summary(reg2Ba)
    
    #b. Treatment Group
    
    reg2Bb <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain ~ dstunting*treatment_group, data = .)
    
    summary(reg2Bb)      
    
#D. Average Weight Gain ----
    
    #1. Main Regressions
    
    #Note LFE::FELM formula : Dependent ~ Independent + Control | Fixed Effects | IV | SE Clustering. 0 if that particular feature is not used.
    
    #0. Base - Nothing else
    rega10 <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      #filter(wtgain < 3 & wtgain > -3)%>%
      felm(wtgain_avg ~ dstunting + dunderweight, data = .)
    
    summary(rega10)
    
    #A. Age controls
    rega1A <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      #filter(wtgain < 3 & wtgain > -3)%>%
      felm(wtgain_avg ~ dstunting + dunderweight + age_m |0|0|0, data = .)
    
    summary(rega1A)
    
    #B. Age controls, Treat Control, Gender Control for cured children
    rega1B <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      felm(wtgain_avg ~ dstunting + dunderweight + gender + age_m + treatment_group  | 0 | 0 | 0, data = .)
    
    summary(rega1B)
    
    
    #B.II Age controls, Treat Control, Gender Control for All children
    rega1B.II <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      felm(wtgain_avg_all ~ dstunting + dunderweight + gender + age_m + treatment_group  | 0 | 0 | 0, data = .)
    
    summary(rega1B.II)
    
    ## SELECTED MODEL
    
    
    #C. Age controls, Treat Controls and Gender Control
    rega1C <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      felm(wtgain_avg ~ dstunting + dunderweight + gender + age_m + treatment_group | 0 | 0 |0 , data = .)
    
    summary(rega1C)
    
    
    #2. DiDs
    
    #A. UW
    
    #a. Age
    
    reg2Aa <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain_avg ~ dunderweight*dage, data = .)
    
    summary(reg2Aa)
    
    #b. Treatment Group
    
    reg2Ab <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain_avg ~ dunderweight*treatment_group, data = .)
    
    summary(reg2Ab)      
    
    #c. Gender
    
    reg2Ac <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain_avg ~ dunderweight*gender, data = .)
    
    summary(reg2Ac)      
    
    

    #B. Stunting
    
    #a. Age
    
    reg2Ba <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain_avg ~ dstunting*dage, data = .)
    
    summary(reg2Ba)
    
    #b. Treatment Group
    
    reg2Bb <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain_avg ~ dstunting*treatment_group, data = .)
    
    summary(reg2Bb)      
    
    #c. Gender
    
    reg2Bc <- data %>%
      filter(dwfh == 1)%>%
      filter(dtreat == 1)%>%
      lm(wtgain_avg ~ dstunting*gender, data = .)
    
    summary(reg2Bc)      
    
    
##VI. Stunting Recovery ----

    #A. Overall percentage
    data%>% 
      filter(dwfh == 1)%>% 
      hablar::convert(num(dstunting))%>%
      summarise(mean(dstunting, na.rm = TRUE))

    data%>% 
      filter(dwfh == 1)%>% 
      hablar::convert(num(dstunting_discharge))%>%
      summarise(mean(dstunting_discharge, na.rm = TRUE))
    
    #B. Percentage of admission stunted that are stunted at discharge
    
    data%>% 
      filter(dwfh == 1)%>% 
      filter(dstunting == 1)%>%
      hablar::convert(num(dstunting_discharge))%>%
      summarise(mean(dstunting_discharge, na.rm = TRUE))
      
    #C. Stunting in admission correlated with stunting in discharge
    
    data%>% 
      filter(dwfh == 1)%>% 
      filter(is.na(dstunting_discharge) == FALSE)%>% 
      filter(is.na(dstunting) == FALSE)%>% 
      hablar::convert(num(dstunting_discharge),
                      num(dstunting))%>%
      select(dstunting_discharge, dstunting)%>%
      cor()
    
    #D. Percentage of admission stunted that are stunted at discharge
    
    data%>% 
      filter(dwfh == 1)%>% 
      filter(dstunting == 1)%>%
      hablar::convert(num(dstunting_discharge))%>%
      group_by(gender)%>%
      summarise((1 - mean(dstunting_discharge, na.rm = TRUE)))
    
    #E. Percentage across categories
    
    #1. Treatment groups
    data%>% 
      filter(dwfh == 1)%>% 
      filter(dstunting == 1)%>%
      hablar::convert(num(dstunting_discharge))%>%
      group_by(treatment_group)%>%
      summarise( (1 - mean(dstunting_discharge, na.rm = TRUE))*100)
    
    #2. Gender
    data%>% 
      filter(dwfh == 1)%>% 
      filter(dstunting == 1)%>%
      hablar::convert(num(dstunting_discharge))%>%
      group_by(gender)%>%
      summarise( (1 - mean(dstunting_discharge, na.rm = TRUE))*100)
    
    #1. Age Categories
    data%>% 
      filter(dwfh == 1)%>% 
      filter(dstunting == 1)%>%
      hablar::convert(num(dstunting_discharge))%>%
      group_by(dage)%>%
      summarise( (1 - mean(dstunting_discharge, na.rm = TRUE))*100)
    
  # F. Regressions
    #1. Overall
    
 regF1 <- data%>%
      filter(dwfh == 1)%>%
      filter(dstunting == 1)%>%
      filter(dtreat == 1)%>%
      hablar::convert(num(dstunting_discharge))%>%
      hablar::convert(num(dstunting))%>%
      lm(data = ., formula = dstunting_cure ~ age_m + treatment_group + gender)
      
summary(regF1)    

  #2. Dstunting_Cure

regF2 <-  data%>%
  filter(dwfh == 1)%>%
  filter(dstunting == 1)%>%
  hablar::convert(num(dstunting_discharge))%>%
  hablar::convert(num(dstunting))%>%
  lm(data = ., formula = dstunting_cure ~ 1)
    
summary(regF2)

  #3. Dstunting cure on a loop

results_agewise_stun_rec <- tibble(age_min = 0, age_max = 0, rec_rate = 0)


for (i in 1:6) {
  
  filt_data <- data%>%
    filter(admission_criterion == "WFH" & dtreat == 1)%>%
    filter(dstunting == 1)%>%
    filter(age_m < i*12 & age_m > (i-1)*12)
  
  
  results_agewise_stun_rec[i,1] <- (i-1)*12
  results_agewise_stun_rec[i,2] <- i*12
  
  results_agewise_stun_rec[i,3] <- filt_data %>% 
    summarise(mean(dstunting_cure, na.rm = TRUE))%>%
    as.numeric()
  
  lm(data = filt_data, formula = dstunting_cure ~ 1)%>%
    summary()%>%
    print()
  
  
  print(i)
  
}








    
#MMMM. Additional / Unsorted -----

      # I
    
    data %>% 
      filter(dwfh == 1) %>% 
      select(rectime, block, dstunting) %>% 
      group_by(block, dstunting) %>%
      summarise(rectime = mean(rectime, na.rm= TRUE))%>%
      View()

      # II
  SandUW  <-  data %>% 
      filter(dwfh == 1) %>% 
      select(dstunting, dwasting, dunderweight) %>%
        filter(dstunting == 1 & dunderweight == 1)%>%
      nrow()
      # 
      
  Sonly <-  data %>% 
        filter(dwfh == 1) %>% 
        select(dstunting, dwasting, dunderweight) %>%
        filter(dstunting == 1)%>%
        nrow()
      
      data %>% 
        filter(dwfh == 1) %>% 
        select(dstunting, dwasting, dunderweight) %>%
        filter(dstunting == 1)%>%
        nrow()
      
      # III
      
      #Rec time UW + Stunted
      data %>% 
        filter(dwfh == 1) %>% 
        select(rectime, dunderweight, dstunting) %>% 
        filter(dstunting == 1 & dunderweight == 1)%>%
        group_by(dstunting)%>%
        summarise(rectimes = sd(rectime, na.rm= TRUE))
     
      #Rec time UW
      data %>% 
        filter(dwfh == 1) %>% 
        select(rectime, dunderweight, dstunting) %>% 
        filter(dunderweight == 1)%>%
        group_by(dunderweight)%>%
        summarise(rectimes.sd = sd(rectime, na.rm= TRUE))
      
        
   t_val <-   (38.9 - 37.6)/(13.8/sqrt(2349) + 14.2/sqrt(2355))
      
      1- pnorm(t_val)
      
      
      #Rec time Stunted
      data %>% 
        filter(dwfh == 1) %>% 
        select(rectime, dunderweight, dstunting) %>% 
        filter(dstunting == 1)%>%
        group_by(dstunting)%>%
        summarise(rectimes = mean(rectime, na.rm= TRUE))
      
      #Rec rate stunted
      
      temp <- data %>% 
        filter(dwfh == 1) %>% 
        select(rectime, dunderweight, dstunting, dcured)
      
      # IV
      
      xtabs(~dcured + dunderweight + dstunting, data = temp)
      
      
      # V
      
      #Logit for UW and Stunting
      
      
      regt1 <- data %>%
        filter(dtreat == 1)%>%
        filter(dwfh == 1)%>%
        glm(dcured ~ dstunting*dunderweight + treatment_group + age_m + gender, family = "binomial", data = .)

      summary(regt1)      
      
      #a. The ORs
      
      regt1.OR.CI <- logistic.regression.or.ci(regt1)
      
      regt1.OR <- regt1.OR.CI$OR
      
      regt1.OR
      
      #b. Percentage Change in probabilities 
      #(One unit increase increases probabilty by X%)
      
      (regt1.OR - 1)*100 
      
      
      ### VI. dAge across stunting
      
      
      table(data$dage, data$dstunting)
      xtabs(~dage + dstunting, data = data)
      
      ### VII. Treatment group DiD
      
      regt2 <- data %>%
        filter(dtreat == 1)%>%
        filter(dwfh == 1)%>%
        glm(dcured ~ dstunting*isarf + dstunting*issf, family = "binomial", data = .)
      
      summary(regt2)      
      
      regt2 <- data %>%
        filter(dwfh == 1)%>%
        lm(rectime ~ dstunting*treatment_group, data = .)
      
      summary(regt2)  
      
      # VIII. 
      
      data %>%
        #filter(outcome == "died")%>%
        filter(dwfh == 1)%>%
        xtabs(~dstunting + outcome, data = .)
      
      
      
      # IX. SUW vs Stunted numbers
      
      
      limits <- c(-Inf, -3 , -2, Inf)
      
      counts <- matrix(nrow = 3, ncol = 3, dimnames = (list(c("SUW", "MUW", "NUW") , c("SS","MS", "NS"))))
      
      for (i in 1:3) {

        uw_u_bound <- limits[i+1]
        uw_l_bound <- limits[i]
        
        
        for(j in 1:3) {
          
          s_u_bound <- limits[j+1]
          s_l_bound <- limits[j]
          
          
          count <- data%>%
            filter(dwfh == 1)%>%
            filter( uw_l_bound < admission_wfa & admission_wfa <= uw_u_bound )%>%
            filter( s_l_bound < admission_lfa & admission_lfa <= s_u_bound )%>%
            nrow()
            
           counts[i,j] <- count 
            
          
        }
        
      }
      
      counts
      
      
      
            #2. a. SAM + Stunting + MUW
              # B. SAM + Stunting + SUW
      
      
      data%>%
        filter(admission_criterion == "WFH")%>%
        filter(admission_wfa >= -3)%>%
        filter(dstunting == 1)%>%
        View()

      data%>%
        filter(admission_criterion == "WFH" & dtreat == 1)%>%
        convert(num(dcured))%>%
        group_by(gender, dstunting)%>%
        summarise(percentage = mean(dcured, na.rm = TRUE))   
      
      
      
      #X. Gender DiD
      
      data %>%
        filter(dtreat == 1)%>%
        filter(dwfh == 1)%>%
        hablar::convert(num(dcured))%>%
        group_by(gender, dstunting)%>%
        summarise(rec_perc = mean(dcured, na.rm = TRUE))%>%
        pivot_wider(names_from = gender, values_from = rec_perc)
      
      
      regt2 <- data %>%
        filter(dtreat == 1)%>%
        filter(dwfh == 1)%>%
        glm(dcured ~ dstunting*gender, family = "binomial", data = .)
      
      summary(regt2)      
      
      regt2 <- data %>%
        filter(dwfh == 1)%>%
        lm(rectime ~ dstunting*gender, data = .)
      
      summary(regt2)  
      
      #XI. Treatment groups distribution
      
      data%>%
        filter(admission_criterion == "WFH" & dtreat == 1)%>%
        convert(num(dcured))%>%
        group_by(treatment_group, dstunting)%>%
        summarise(percentage = mean(dcured, na.rm = TRUE))%>%
        pivot_wider(names_from = treatment_group, values_from = percentage)
      
      
      ## Logistic to check significance
      
   regMMII <-  data%>%
        filter(admission_criterion == "WFH" & dtreat == 1)%>%
        convert(num(dcured))%>%
        glm(dcured ~ isarf + issf, data = .)
      
      summary(regMMII)
      
      
      #XII. Wasted + UW vs Wasted + No UW Recovery Percentage
      
      data%>%
        filter(admission_criterion == "WFH" & dtreat == 1)%>%
        convert(num(dcured))%>%
        group_by(dunderweight)%>%
        summarise(percentage = mean(dcured, na.rm = TRUE))
    #XIII.  
    regt13 <-  data%>%
        filter(admission_criterion == "WFH" & dtreat == 1)%>%
      glm(dcured ~ dunderweight, data = . , family ="binomial")

    summary(regt13)      
      
    #XIV. Gender Percentage recovery
    
    data%>%
      filter(admission_criterion == "WFH" & dtreat == 1)%>%
      convert(num(dcured))%>%
      group_by(dage)%>%
      summarise(percentage = mean(dcured, na.rm = TRUE))
    
    
    #XV. Stuff across age categories
    
    results_agewise_rate <- tibble(age_min = 0, age_max = 0, not_stunted = 0 , stunted = 0, all = 0)
    results_agewise_wtgain <- tibble(age_min = 0, age_max = 0, not_stunted = 0 , stunted = 0, all = 0)
    results_agewise_time <- tibble(age_min = 0, age_max = 0, not_stunted = 0 , stunted = 0, all = 0)
    results_agewise_stun_rec <- tibble(age_min = 0, age_max = 0, rec_rate = 0)
    
    
    for (i in 1:6) {
      
  filt_data <- data%>%
      filter(admission_criterion == "WFH" & dtreat == 1)%>%
      filter(age_m < i*12 & age_m > (i-1)*12)%>%
      convert(num(dcured),
              num(dstunting_cure)
              )
     
  results_agewise_rate[i,1] <- results_agewise_wtgain[i,1]  <- results_agewise_time[i,1] <- results_agewise_stun_rec[i,1] <- (i-1)*12
  results_agewise_rate[i,2] <- results_agewise_wtgain[i,2]  <- results_agewise_time[i,2] <-  results_agewise_stun_rec[i,2] <- i*12
  
results_agewise_rate[i, 3:4] <- 
  filt_data%>%
    convert(num(dcured))%>%
    group_by(dstunting)%>%
    summarise(rate = 100*mean(dcured, na.rm = TRUE)%>%
                     as.numeric(), 
              .groups = 'drop')%>%
    t()%>%
    as_tibble()%>%
    select(1:2)%>%
    slice(2)%>%
    convert(num(V1),
            num(V2)
    )%>%
    as.vector()

results_agewise_rate[i, 5] <- 
  filt_data%>%
  convert(num(dcured))%>%
  summarise(rate = 100*mean(dcured, na.rm = TRUE))%>%
              as.numeric()

results_agewise_time[i, 3:4] <- 
  filt_data%>%
  convert(num(dcured))%>%
  group_by(dstunting)%>%
  summarise(rate = mean(rectime, na.rm= TRUE)%>%
              as.numeric(), 
            .groups = 'drop')%>%
  t()%>%
  as_tibble()%>%
  select(1:2)%>%
  slice(2)%>%
  convert(num(V1),
          num(V2)
  )%>%
  as.vector()

results_agewise_time[i, 5] <- 
  filt_data%>%
  convert(num(dcured))%>%
  summarise(rate = mean(rectime, na.rm= TRUE))%>%
  as.numeric() 
        
results_agewise_wtgain[i, 3:4] <- 
  filt_data%>%
  convert(num(dcured))%>%
  group_by(dstunting)%>%
  summarise(rate = mean(wtgain_avg, na.rm = TRUE)%>%
              as.numeric(), 
            .groups = 'drop')%>%
  t()%>%
  as_tibble()%>%
  select(1:2)%>%
  slice(2)%>%
  convert(num(V1),
          num(V2)
  )%>%
  as.vector()

results_agewise_wtgain[i, 5] <- 
  filt_data%>%
  convert(num(dcured))%>%
  summarise(rate = mean(wtgain_avg, na.rm = TRUE))
              
results_agewise_stun_rec[i,3] <- filt_data %>% 
  summarise(mean(dstunting_cure, na.rm = TRUE))%>%
  as.numeric()

lm(data = filt_data, formula = dstunting_cure ~ 1)%>%
  summary()%>%
  print()


print(i)
    
    }
    
results_agewise_rate_wtgain_time <- bind_rows(results_agewise_rate, results_agewise_wtgain) %>% bind_rows(results_agewise_time)    
    
    
    #XVI. Composite Index


data%>%
  filter(admission_criterion == "WFH" & dtreat == 1)%>%
  group_by(dunderweight, dstunting)%>%
  count()
  


    
    
    
    
    
    
    # Preliminaries for the MUAC study ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    data %>% filter(admission_criterion == "WFH", admission_muac > 0)%>% 
      select(admission_muac)%>% 
      summary()
    
    
    data %>% filter(admission_criterion == "WFH", admission_muac > 0)%>% 
      select(admission_muac)%>% 
      ggplot(aes(admission_muac))+
      geom_histogram(bins = 0, colour = "dark green", fill = "green", alpha = 0.4)+
      scale_x_continuous(limits = c(10,15), breaks = seq(10,15,by = .5))+
      theme_bw()
    
    
    data %>% filter(admission_criterion == "WFH", admission_muac > 0)%>% 
      select(admission_muac)%>% 
      ggplot(aes(admission_muac))+
      stat_ecdf()+
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.1))+
      scale_x_continuous(limits = c(10,15), breaks = seq(10,15,by = .5))+
      ylab("Cumulative Proportion of Children")+
      theme_bw()
    
    
    data %>% filter(admission_criterion == "WFH", admission_wfh < -2, admission_muac > 0)%>% 
      ggplot(aes(admission_wfh, admission_muac))+
      geom_point(size = 0.3)+
      geom_smooth()+
      ggtitle("Children admitted on WFH")+
      theme_bw()
    
    
    
    data%>%
      convert(num(dcured))%>%
      group_by(admission_criterion, dstunting)%>%
      summarise(mean(dcured, na.rm = TRUE)*100)
    
    