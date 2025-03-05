library(heemod)
library(haven)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readxl)
library(stringr)
library(flextable)
library(ggsignif)
library(nlme)
library(cutpointr)
library(survival)
library(survminer)
library(pbapply) 
library(data.table)
library(flexsurv)
library(ConfidenceEllipse)



# load data
setwd("C:/Users/xinxia/OneDrive - Karolinska Institutet/CSF-registersamkörning/AD DMT HE model/Statistical analyses/New analyses_XX_202412/")

load("costdata_xx.RData")
load("models_xx.RData")
scb_mort<-as.data.frame(read_excel("life_table_scb.xlsx"))



# specify models for each transition,except for transitions from MCI
model_trans4<-final_model_trans1 # mild AD to moderate AD
model_trans5<-final_model_trans2 # mild AD to institutionalization
model_trans6<-final_model_trans3 # mild AD to death
model_trans7<-final_model_trans4 # moderate AD to severe AD
model_trans8<-final_model_trans5 # moderate AD to institutionalization
model_trans9<-final_model_trans6 # moderate AD to death
model_trans10<-final_model_trans7 # severe AD to institutionalization
model_trans11<-final_model_trans8 # severe AD to death

model_trans14<-final_model_trans9 # institutionalized mild AD to institutionalized moderate AD
model_trans15<-final_model_trans10 # institutionalized mild AD to death
model_trans16<-final_model_trans11 # institutionalized moderate AD to institutionalized severe AD
model_trans17<-final_model_trans12 # institutionalized moderate AD to death
model_trans18<-final_model_trans13 # institutionalized severe AD to death


# ************************************************************************************************
# Considerations:
# 1. The model for each transition was selected based on previous analyses specified in the 
#    R script "Flexible parametric multistate model_XX_202501.R".
# 2. Consider the time horizon for later HE modeling when predicting transition probabilities:
#    cycle_length=3 #in months (0.25 years)
#    cycles=12/cycle_length*10 #10 year time horizon (40 cycles)
# 3. Considering APOE e4 genotypes and adverse events
# ************************************************************************************************


# 1. Get transition-specific hazard ----
# get baseline hazard for each cycle for each transition
cycle_length<-0.25 #in years
cycles<-10/cycle_length #10 year time horizon


cumhaz_table <- do.call(rbind, lapply(c(4:11,14:18), function(i) {
  male_data <- predict(get(paste0("model_trans", i)), type = "cumhaz", 
                       newdata = data.frame(AGE = 70, SEX = "MALE"),
                       times = seq(0, 10, by = 0.25))[[1]][[1]] %>% mutate(SEX = 0)
  
  female_data <- predict(get(paste0("model_trans", i)), type = "cumhaz", 
                         newdata = data.frame(AGE = 70, SEX = "FEMALE"),
                         times = seq(0, 10, by = 0.25))[[1]][[1]] %>% mutate(SEX = 1)
  
  rbind(male_data, female_data) %>%
    mutate(trans = i)
})
)

cumhaz_table<-cumhaz_table %>% 
  group_by(trans,SEX) %>% 
  mutate(haz=.pred_cumhaz-lag(.pred_cumhaz),
         cycle=row_number()-1) %>% 
  ungroup %>% 
  mutate(haz=ifelse(.time==0,0,haz)) %>% 
  filter(cycle!=0)



# get hazard for each transition
# account for the effect of time-varying age (hazard ratio) on the hazard for each transition
get_haz<-function(transition_number, age_init, sex, hazards) {
  #Get hazard for each transition at each time point
  age_coef=(get(paste0("model_trans", transition_number)))$coefficients[["AGE"]]
  h=hazards %>% filter(trans==transition_number)   #Base hazard without covariates
  h=h %>% 
    mutate(hr_age=exp(age_coef*(age_init+.time-0.25-70)),#Hazard ratio changes over time due to increases in age (center around 70)
           haz=haz*hr_age,
           trans=transition_number) # calculate transition probability
  return(h %>% 
           filter(SEX==sex) %>% 
           select(haz,trans)) #return vector of transition probabilities by cycle
}


# 2. Extract costs and utilities by states ----
# create a function to extract costs
getcost<-function(df) {
  costtable<-df %>% group_by(state) %>%
    mutate(state_duration=as.numeric(next_date-DATE)) %>% 
    summarise(cost=sum(totalcost), 
              time=sum(state_duration)/365.25,
              annual_cost=cost/time, .groups='drop') %>%
    arrange(state) 
  return(costtable$annual_cost)
}



# extract utility
source("health utility R script_XX_20250103.R")
utilities<-as.vector(state_utility$mean)
utilities_lb<-as.vector(state_utility$lb)
utilities_ub<-as.vector(state_utility$ub)


# 3. Run HE model in heemod ----
statenames = c(
  "MCI",
  "Mild",
  "Moderate",
  "Severe",
  "MCI_inst",
  "Mild_inst",
  "Moderate_inst",
  "Severe_inst",
  "Death")



# create a customized "run_model" function to run the HE model
calimodel_cust<-function(rx_cycles, 
                        age_init, 
                        sex,
                        apoegeno,
                        start_state,
                        setting,
                        wane_fc, # treatment waning factor
                        cost_mci,
                        cost_ad
)
{
  haz=lapply(c(4:11,14:18), function(x) get_haz(x,age_init=age_init,sex=sex,cumhaz_table)) #get transition probabilities per transition and cycle
  haz<-do.call(rbind.data.frame, haz)
  
  param <- define_parameters(
    
    #utility of each state
    utilities_mci=utilities[1], 
    utilities_mild=utilities[2], 
    utilities_mod=utilities[3], 
    utilities_sev=utilities[4], 
    utilities_mci_inst=utilities[5], 
    utilities_mild_inst=utilities[6], 
    utilities_mod_inst=utilities[7], 
    utilities_sev_inst=utilities[8], 
    utilities_death=0, 
    
    #discount rate
    r=0.03,  
    disc=exp(-r*(model_time-1)*cycle_length), #continuously compounded discounting factor
    
    #disc=1/(1+r)^(model_time-1)*cycle_length
    
    rx_eff_log=log(0.69),
    rx_eff=exp(rx_eff_log),
    
    # consider waning effect after treatment ends
    rX=dispatch_strategy(
      standard=1,
      rx=case_when(model_time<=rx_cycles ~ rx_eff, 
                   model_time>rx_cycles ~ rx_eff^((1-wane_fc)^(model_time-rx_cycles)))
    ),
    
    
    #Variable to use to record time on treatment
    rxtime=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=rx_cycles ~ 1, T ~ 0)
    ),
    
    # additional costs and disutility associated with Lecanemab
    # infusion cost
    cost_infusion=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=6 ~ 3341*(12*cycle_length*2), # biweekly during the first 18 months
                   model_time>6&model_time<= rx_cycles~ 3341*(12*cycle_length), # every 4 weeks after 18 months
                   T ~ 0)
    ),
    
    # MRI cost
    cost_mri=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=1 ~ 2469*2, 
                   model_time<=4 ~ 2469,
                   T ~ 0) # 1, 5 (2.5 months), 7 (3.5 months), 14th infusion (7 months), and 52 week (13 months)
    ),
    
    # physician visit costs
    cost_phy_visit=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=1 ~ 6703*2, 
                   model_time<=4 ~ 6703,
                   T ~ 0) # same scheme as MRI
    ),
    
    
    cost_admin=cost_infusion+cost_mri+cost_phy_visit,
    
    
    # 18-month probability of non-severe ARIA-E among 
    # APOE e4 non-carrier in Lecanemab: 0.054*0.867, among heterozygote: 0.109*0.867
    
    prob_sev_ariae_18month=ifelse(apoegeno==0,0.054*0.133,
                                  0.109*0.133),
    
    prob_mild_ariae_18month=ifelse(apoegeno==0,0.054*(1-0.133),
                                   0.109*(1-0.133)),
    
    prob_ariae_soc=ifelse(apoegeno==0,0.003,0.019),
    
    
    haz_mild_ariae=dispatch_strategy(
      standard=0,
      rx=case_when(model_time==1~-log(1-prob_mild_ariae_18month*0.7), # hazard for first circle
                   model_time==2~-log(1-prob_mild_ariae_18month*0.2), # hazard for second circle
                   model_time>=3&model_time<=rx_cycles~-log(1-prob_mild_ariae_18month*0.1)/4, # constant hazard the next circles
                   T~-log(1-prob_ariae_soc)/6)),
    
    prob_mild_ariae=1-exp(-haz_mild_ariae), # transition probability
    
    
    # 18-month probability of severe ARIA-E among APOE e4 non-carrier in Lecanemab:  
    haz_severe_ariae=dispatch_strategy(
      standard=0,
      rx=case_when(model_time==1~-log(1-prob_sev_ariae_18month*0.7), # hazard for first circle
                   model_time==2~-log(1-prob_sev_ariae_18month*0.2), # hazard for second circle
                   model_time>=3&model_time<=rx_cycles~-log(1-prob_sev_ariae_18month*0.1)/4, # constant hazard the next circles
                   T~0)),
    
    prob_severe_ariae=1-exp(-haz_severe_ariae),
    
    # 18-month probability of non-severe ARIA-H among 
    # APOE e4 non-carrier in Lecanemab: 0.079*0.955, heterozygote: 0.081*0.955
    
    prob_sev_ariah_18month=ifelse(apoegeno==0,0.079*(1-0.955),
                                  0.081*(1-0.955)),
    
    prob_mild_ariah_18month=ifelse(apoegeno==0,0.079*0.955,
                                   0.081*0.955),
    
    prob_ariah_soc=ifelse(apoegeno==0,0.035,0.073),
    
    haz_mild_ariah=dispatch_strategy(
      standard=0, 
      rx=case_when(model_time<=rx_cycles~-log(1-prob_mild_ariah_18month)/6, # constant hazard 
                   T~-log(1-prob_ariah_soc)/6)),
    
    prob_mild_ariah=1-exp(-haz_mild_ariah), # transition probability per cycle
    
    # 18-month probability of severe ARIA-H among APOE e4 non-carrier in Lecanemab: 0.079*0.045    
    haz_severe_ariah=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=rx_cycles~-log(1-prob_sev_ariah_18month)/6, # constant hazard 
                   T~0)),
    
    prob_severe_ariah=1-exp(-haz_severe_ariah),
    
    # costs due to ARIA
    cost_mild_aria=13173,
    cost_severe_aria=498076.3223,
    
    
    prob_mild_aria=prob_mild_ariae+prob_mild_ariah,
    prob_severe_aria=prob_severe_ariae+prob_severe_ariah,
    
    
    # costs due to severe irr
    haz_irr=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=rx_cycles~-log(1-0.012)/6, # constant hazard 
                   T~0)),
    prob_irr=1-exp(-haz_irr),
    cost_irr=201.27,
    
    # only serious AE affect utility
    disutility_aria=-0.065,
    disutility_irr=-0.03,
    
    
    # transition probabilities for MCI from a previous study
    # hazard for APOE e4 non-carrier: 0.0969-0.1848
    # from a population of MCI with a mean age of 73
    mci_ad_haz=0.1338,
    mci_ad_age_loghr=log(1.06), # RR = 1.06 (1.03-1.1)
    mci_ad_female_loghr=log(1.33),
    mci_ad_sex_loghr=ifelse(sex==1,mci_ad_female_loghr,0), # RR = 1.33 (1.08-1.64)
    
    mci_ad_apoecarrier_loghr=log(1.82),
    mci_ad_apoe_loghr=ifelse(apoegeno==1,mci_ad_apoecarrier_loghr,0), # RR = 1.82 (1.49-2.22)
    
    # cumulative hazard for MCI per cycle, taking into account demographic factors
    mci_ad_haz_demo=mci_ad_haz*(exp(mci_ad_age_loghr*((model_time-1)*cycle_length-3)))*
      exp(mci_ad_sex_loghr)*exp(mci_ad_apoe_loghr)*cycle_length,
    
    # transition probability for MCI to institutionalization (from a population with mean age of 74.5)
    mci_inst_haz=-(log(1-0.043))/6,
    mci_inst_age_loghr=log(1.09), # 1.09 (1.05-1.13)
    mci_inst_haz_demo=mci_inst_haz*(exp(mci_inst_age_loghr*((model_time-1)*cycle_length-4.5)))*cycle_length,
    
    
    # transition probability for MCI to death, using life table from SCB
    age=round(age_init+model_time-50+1),
    tp_mci_death=ifelse(sex==0,scb_mort[age,2],
                        scb_mort[age,3]),
    
    
    t1=1-exp(-mci_ad_haz_demo*rX),
    t2=1-exp(-mci_inst_haz_demo),
    t3=tp_mci_death*cycle_length,
    
    t4=1-exp(-(haz %>% filter(trans==4))[model_time,1]*rX),
    t5=1-exp(-(haz %>% filter(trans==5))[model_time,1]),
    t6=1-exp(-(haz %>% filter(trans==6))[model_time,1]),
    t7=1-exp(-(haz %>% filter(trans==7))[model_time,1]),
    t8=1-exp(-(haz %>% filter(trans==8))[model_time,1]),
    t9=1-exp(-(haz %>% filter(trans==9))[model_time,1]),
    t10=1-exp(-(haz %>% filter(trans==10))[model_time,1]),
    t11=1-exp(-(haz %>% filter(trans==11))[model_time,1]),
    
    t12=1-exp(-mci_ad_haz_demo*rX),
    t13=tp_mci_death*cycle_length,
    
    t14=1-exp(-(haz %>% filter(trans==14))[model_time,1]*rX),
    t15=1-exp(-(haz %>% filter(trans==15))[model_time,1]),
    t16=1-exp(-(haz %>% filter(trans==16))[model_time,1]),
    t17=1-exp(-(haz %>% filter(trans==17))[model_time,1]),
    t18=1-exp(-(haz %>% filter(trans==18))[model_time,1]),
  )
  
  # Define states, costs and utilities----
  MCI = define_state(
    utility = (utilities_mci*(1-prob_severe_aria-prob_irr)+
                 (utilities_mci+disutility_aria)*prob_severe_aria+
                 (utilities_mci+disutility_irr)*prob_irr)
    *cycle_length*disc, 
    
    cost = ((cost_mci[1]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_mci[1]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_mci[1]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_mci[1]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length, # treatment time/cycle = 1/(no.cycles/year) - full treatment time every cycle
    rx.time.disc=rx.time*disc,
    cdr=0.5
  )
  Mild = define_state(
    utility = (utilities_mild*(1-prob_severe_aria-prob_irr)+
                 (utilities_mild+disutility_aria)*prob_severe_aria+
                 (utilities_mild+disutility_irr)*prob_irr)
    *cycle_length*disc, 
    
    cost = ((cost_ad[1]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_ad[1]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_ad[1]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_ad[1]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length,
    rx.time.disc=rx.time*disc,
    cdr=1
  )
  Moderate = define_state(
    utility = utilities_mod*cycle_length*disc,
    cost = cost_ad[2]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=2
  )
  Severe = define_state(
    utility = utilities_sev*cycle_length*disc,
    cost = cost_ad[3]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=3
  )
  MCI_inst = define_state(
    utility = (utilities_mci_inst*(1-prob_severe_aria-prob_irr)+
                 (utilities_mci_inst+disutility_aria)*prob_severe_aria+
                 (utilities_mci_inst+disutility_irr)*prob_irr)
    *cycle_length*disc, 
    
    cost = ((cost_mci[2]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_mci[2]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_mci[2]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_mci[2]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length,
    rx.time.disc=rx.time*disc,
    cdr=0.5
  )
  Mild_inst = define_state(
    utility = (utilities_mild_inst*(1-prob_severe_aria-prob_irr)+
                 (utilities_mild_inst+disutility_aria)*prob_severe_aria+
                 (utilities_mild_inst+disutility_irr)*prob_irr)
    *cycle_length*disc, 
    
    cost = ((cost_ad[4]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_ad[4]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_ad[4]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_ad[4]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length,
    rx.time.disc=rx.time*disc,
    cdr=1
  )
  Moderate_inst = define_state(
    utility = utilities_mod_inst*disc*cycle_length,
    cost = cost_ad[5]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=2
  )
  Severe_inst = define_state(
    utility = utilities_sev_inst*cycle_length*disc,
    cost = cost_ad[6]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=3
  )
  Death = define_state(
    utility = 0,
    cost = 0,
    rx.time=0,
    rx.time.disc=0,
    cdr=0
  )
  
  # Transitions----
  transmat<- define_transition(
    state_names = statenames,
    C, t1, 0,     0,     t2, 0,      0,      0,      t3,
    0, C,     t4, 0,     0,  t5,     0,      0,      t6, 
    0, 0,     C,     t7,    0,  0,      t8,     0,      t9, 
    0, 0,     0,     C,     0,  0,      0,      t10,    t11,
    0, 0,     0,     0,     C,  t12, 0,      0,      t13,
    0, 0,     0,     0,     0,  C,      t14, 0,      t15,
    0, 0,     0,     0,     0,  0,      C,      t16,    t17, 
    0, 0,     0,     0,     0,  0,      0,      C,      t18,
    0, 0,     0,     0,     0,  0,      0,      0,      C
  )
  
  # Strategies----
  strat_standard <- define_strategy(
    MCI=MCI,
    Mild=Mild,
    Moderate=Moderate,
    Severe=Severe,
    MCI_inst=MCI_inst,
    Mild_inst=Mild_inst,
    Moderate_inst=Moderate_inst,
    Severe_inst=Severe_inst,
    Death=Death,
    transition = transmat)
  
  init_state<-if (start_state == 0) {
    c(1, 0, 0, 0, 0, 0, 0, 0, 0)
  } else if(start_state == 1&setting==0) {
    c(0, 1, 0, 0, 0, 0, 0, 0, 0)
  } else{
    c(0, 0, 0, 0, 0, 1, 0, 0, 0)
  }
  
  res_mod <- run_model(
    standard = strat_standard,
    rx=strat_standard,
    parameters = param,
    init=init_state,
    cycles = cycles,
    cost = cost,
    effect = cdr,
    method='beginning'
  ) 
  return(res_mod)
}


# call starting population
load("pop_perc.RData")


cali_result<-lapply(1:nrow(pop_perc), function(i){
  age<-pop_perc[i,]$mean_age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  setting_i<-pop_perc[i,]$inst
  perc<-pop_perc[i,]$pop_perc
  
  temp_cali_model<-calimodel_cust(rx_cycles=6,age_init = age,sex=sex_i,
                             apoegeno=apoe_i,start_state=stage_i,
                             setting=setting_i,
                             wane_fc=1,
                             cost_mci=getcost(mci_svedem_costs),
                             cost_ad=getcost(costwide))
  
  dcdr6<-temp_cali_model$eval_strategy_list$standard$values$cdr[6]-
    temp_cali_model$eval_strategy_list$rx$values$cdr[6]
  
  dcdr14<-temp_cali_model$eval_strategy_list$standard$values$cdr[14]-
    temp_cali_model$eval_strategy_list$rx$values$cdr[14]
  
  return(data.frame(age=age,sex=sex_i,apoe=apoe_i,stage=stage_i,setting=setting_i, perc=perc,
              dcdr6=dcdr6,dcdr14=dcdr14))
}) %>% 
  bind_rows()


cali_result %>% 
  mutate(dcdr6=dcdr6*perc,
         dcdr14=dcdr14*perc) %>% 
  summarise(dcdr6=sum(dcdr6),
            dcdr14=sum(dcdr14))
