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



# specify models for each transition,except for transitions from MCI
model_trans1<-mci_ad_model # MCI to mild AD
model_trans3<-mci_death_model # MCI to death
model_trans4<-final_model_trans1 # mild AD to moderate AD
model_trans5<-final_model_trans2 # mild AD to institutionalization
model_trans6<-final_model_trans3 # mild AD to death
model_trans7<-final_model_trans4 # moderate AD to severe AD
model_trans8<-final_model_trans5 # moderate AD to institutionalization
model_trans9<-final_model_trans6 # moderate AD to death
model_trans10<-final_model_trans7 # severe AD to institutionalization
model_trans11<-final_model_trans8 # severe AD to death
model_trans12<-mci_ad_model # institutionalized MCI to institutionalized mild AD
model_trans13<-mci_death_model # institutionalized MCI to death
model_trans14<-final_model_trans9 # institutionalized mild AD to institutionalized moderate AD
model_trans15<-final_model_trans10 # institutionalized mild AD to death
model_trans16<-final_model_trans11 # institutionalized moderate AD to institutionalized severe AD
model_trans17<-final_model_trans12 # institutionalized moderate AD to death
model_trans18<-final_model_trans13 # institutionalized severe AD to death


# ************************************************************************************************
# Considerations:
# The objective of the calibration is to estimate the waning factor that makes the 
# difference in CDR at 18 months between the standard and the treatment arm equal to 
# the difference in CDR at 18 months + 2 years between the standard and the treatment arm.
# Reference (Eisai's presentation at AAIC2024): 
# Perspectives Session Lecanemab - 30July2024 FINAL (1)
# ************************************************************************************************


# 1. Get transition-specific hazard ----
# get baseline hazard for each cycle for each transition
cycle_length<-0.25 #in years
cycles<-10/cycle_length #10 year time horizon


# write a function to get hazard
get_haz<-function(trans, age_init, sex, apoe) {
  
  model<-get(paste0("model_trans", trans))
  
  sex<-ifelse(sex==0,"MALE","FEMALE")
  
  if (trans==1|trans==12){ # MCI transition to AD model includes APOE
    newdata<-data.frame(AGE = age_init, SEX = sex, NACCNE4S = apoe)
  } else {
    newdata<-data.frame(AGE = age_init, SEX = sex)
  }
  
  cumhaz_data<-predict(model, type = "cumhaz", 
                       newdata = newdata,
                       times = seq(0, 10, by = 0.25))[[1]][[1]] %>%
    rename(.time = .eval_time)
  
  cumhaz_data <- cumhaz_data %>% 
    mutate(haz=.pred_cumhaz-lag(.pred_cumhaz),
           cycle=row_number()-1) %>% 
    ungroup %>% 
    mutate(haz=ifelse(.time==0,0,haz)) %>% 
    filter(cycle!=0)
  
  return(cumhaz_data %>% select(haz))
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
                        cost_ad,
                        r # discount rate
)
{param <- define_parameters(
    
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
      rx=case_when(model_time<=6&model_time<= rx_cycles ~ 3626.581*(12*cycle_length*2), # biweekly during the first 18 months
                   model_time>6&model_time<= rx_cycles~ 3626.581*(12*cycle_length), # every 4 weeks after 18 months
                   T ~ 0)
    ),
    
    # MRI cost
    cost_mri=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=1 ~ 2401*2, 
                   model_time<=4&model_time<= rx_cycles ~ 2401,
                   T ~ 0) # 1, 5 (2.5 months), 7 (3.5 months), 14th infusion (7 months), and 52 week (13 months)
    ),
    
    # physician visit costs
    cost_phy_visit=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=1 ~ 3774*2, 
                   model_time<=4&model_time<= rx_cycles ~ 3774,
                   T ~ 0) # same scheme as MRI
    ),
    
    
    cost_admin=cost_infusion+cost_mri+cost_phy_visit,
    
    
    # Probability of non-severe ARIA-E
    prob_sev_ariae=ifelse(apoegeno==0,0.002,0.008),
    
    prob_mild_ariae=ifelse(apoegeno==0,0.014,0.013),
    
    haz_mild_ariae=dispatch_strategy(
      standard=0,
      rx=case_when(model_time==1~-log(1-prob_mild_ariae*0.7), # hazard for first circle
                   model_time==2~-log(1-prob_mild_ariae*0.2), # hazard for second circle
                   model_time>=3&model_time<=rx_cycles~-log(1-prob_mild_ariae*0.1)/5, # constant hazard the next circles
                   T~0)),
    
    prob_mild_ariae=1-exp(-haz_mild_ariae), # transition probability
    
    haz_severe_ariae=dispatch_strategy(
      standard=0,
      rx=case_when(model_time==1~-log(1-prob_sev_ariae*0.7), # hazard for first circle
                   model_time==2~-log(1-prob_sev_ariae*0.2), # hazard for second circle
                   model_time>=3&model_time<=rx_cycles~-log(1-prob_sev_ariae*0.1)/5, # constant hazard the next circles
                   T~0)),
    
    prob_severe_ariae=1-exp(-haz_severe_ariae),
    
    # Probability of non-severe ARIA-H
    prob_sev_ariah=ifelse(apoegeno==0,0.004,0.003),
    
    prob_mild_ariah=ifelse(apoegeno==0,0.003,0.004),
    
    haz_mild_ariah=dispatch_strategy(
      standard=0, 
      rx=case_when(model_time<=rx_cycles~-log(1-prob_mild_ariah)/7, # constant hazard 
                   T~0)),
    
    prob_mild_ariah=1-exp(-haz_mild_ariah), # transition probability per cycle
    
    haz_severe_ariah=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=rx_cycles~-log(1-prob_sev_ariah)/7, # constant hazard 
                   T~0)),
    
    prob_severe_ariah=1-exp(-haz_severe_ariah),
    
    # costs due to ARIA
    cost_mild_aria=13037,
    cost_severe_aria=105102,
    
    
    prob_mild_aria=prob_mild_ariae+prob_mild_ariah,
    prob_severe_aria=prob_severe_ariae+prob_severe_ariah,
    
    
    # costs due to severe irr
    haz_irr=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<=rx_cycles~-log(1-0.012)/7, # constant hazard 
                   T~0)),
    prob_irr=1-exp(-haz_irr),
    cost_irr=195,
    
    # only serious AE affect utility
    disutility_aria=-0.065,
    disutility_irr=-0.03,
    
    
    # transition probability for MCI to institutionalization (from a population with mean age of 74.5)
    mci_inst_haz=-(log(1-0.043))/6,
    mci_inst_age_loghr=log(1.09), # 1.09 (1.05-1.13)
    mci_inst_haz_demo=mci_inst_haz*cycle_length*
      (exp(mci_inst_age_loghr*(age_init-74.5))),
    
    
    t1=1-exp(-(get_haz(1,age_init,sex,apoegeno))[model_time,1]*rX),
    t2=1-exp(-mci_inst_haz_demo),
    t3=1-exp(-(get_haz(3,age_init,sex))[model_time,1]),
    t4=1-exp(-(get_haz(4,age_init,sex))[model_time,1]*rX),
    t5=1-exp(-(get_haz(5,age_init,sex))[model_time,1]),
    t6=1-exp(-(get_haz(6,age_init,sex))[model_time,1]),
    t7=1-exp(-(get_haz(7,age_init,sex))[model_time,1]),
    t8=1-exp(-(get_haz(8,age_init,sex))[model_time,1]),
    t9=1-exp(-(get_haz(9,age_init,sex))[model_time,1]),
    t10=1-exp(-(get_haz(10,age_init,sex))[model_time,1]),
    t11=1-exp(-(get_haz(11,age_init,sex))[model_time,1]),
    t12=1-exp(-(get_haz(12,age_init,sex,apoegeno))[model_time,1]*rX),
    t13=1-exp(-(get_haz(13,age_init,sex))[model_time,1]),
    t14=1-exp(-(get_haz(14,age_init,sex))[model_time,1]*rX),
    t15=1-exp(-(get_haz(15,age_init,sex))[model_time,1]),
    t16=1-exp(-(get_haz(16,age_init,sex))[model_time,1]),
    t17=1-exp(-(get_haz(17,age_init,sex))[model_time,1]),
    t18=1-exp(-(get_haz(18,age_init,sex))[model_time,1])
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
    cdr=0.5+(4.5-0.5)/2
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
    cdr=4.5+(9.5-4.5)/2
  )
  Moderate = define_state(
    utility = utilities_mod*cycle_length*disc,
    cost = cost_ad[2]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=9.5+(16-9.5)/2
  )
  Severe = define_state(
    utility = utilities_sev*cycle_length*disc,
    cost = cost_ad[3]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=16+(18-16)/2
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
    cdr=0.5+(4.5-0.5)/2
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
    cdr=4.5+(9.5-4.5)/2
  )
  Moderate_inst = define_state(
    utility = utilities_mod_inst*disc*cycle_length,
    cost = cost_ad[5]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=9.5+(16-9.5)/2
  )
  Severe_inst = define_state(
    utility = utilities_sev_inst*cycle_length*disc,
    cost = cost_ad[6]*cycle_length*disc,
    rx.time=0,
    rx.time.disc=0,
    cdr=16+(18-16)/2
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



load("pop_perc.RData")

# setting up parallel computing
ncpus = parallel::detectCores()-1
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod)))
clusterExport(cl, c("pop_perc","calimodel_cust","get_haz",
                    "statenames","cycles","cycle_length","utilities",
                    "getcost","mci_svedem_costs","costwide",
                    paste0("model_trans",c(1,3:18))))

cali_result_fun<-function(j){
  cali_result<-parLapply(cl,1:nrow(pop_perc), function(i){
    rx_cycles_i<-pop_perc[i,]$rx_cycles
    age<-pop_perc[i,]$age
    sex_i<-pop_perc[i,]$sex
    apoe_i<-pop_perc[i,]$apoe
    stage_i<-pop_perc[i,]$stage
    setting_i<-pop_perc[i,]$inst
    perc<-pop_perc[i,]$pop_perc
    
    temp_cali_model<-calimodel_cust(rx_cycles=rx_cycles_i,age_init = age,sex=sex_i,
                                    apoegeno=apoe_i,start_state=stage_i,
                                    setting=setting_i,
                                    wane_fc=j,
                                    cost_mci=getcost(mci_svedem_costs),
                                    cost_ad=getcost(costwide),
                                    r=0.03)
    
    cdr6_stand<-temp_cali_model$eval_strategy_list$standard$values$cdr[6]
    cdr6_rx<-temp_cali_model$eval_strategy_list$rx$values$cdr[6]
    
    prob_death6_stand<-temp_cali_model$eval_strategy_list$standard$counts$Death[6] 
    prob_death6_rx<-temp_cali_model$eval_strategy_list$rx$counts$Death[6]
    
    dcdr6<-cdr6_stand/(1-prob_death6_stand)-cdr6_rx/(1-prob_death6_rx)
    
    cdr14_stand<-temp_cali_model$eval_strategy_list$standard$values$cdr[14]
    cdr14_rx<-temp_cali_model$eval_strategy_list$rx$values$cdr[14]
    
    prob_death14_stand<-temp_cali_model$eval_strategy_list$standard$counts$Death[14] 
    prob_death14_rx<-temp_cali_model$eval_strategy_list$rx$counts$Death[14]
    
    dcdr14<-cdr14_stand/(1-prob_death14_stand)-cdr14_rx/(1-prob_death14_rx)
    
    return(data.frame(perc=perc,dcdr6=dcdr6,dcdr14=dcdr14))
  }) %>% 
    bind_rows()
  
  cali_result<-cali_result %>% 
    mutate(dcdr6=dcdr6*perc,
           dcdr14=dcdr14*perc) %>% 
    summarise(dcdr6=sum(dcdr6),
              dcdr14=sum(dcdr14))
  
  return(cali_result)
}
  
  

cali_result_fun(0.1)
cali_result_fun(0.5)
cali_result_fun(1)



stopCluster(cl)