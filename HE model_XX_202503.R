rm(list = ls())

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
library(parallel)
library(ggpubr)
library(boot)


# load data
load("costdata_xx.RData")
load("models_xx.RData")

# call starting population
load("profile_perc.RData")


# function for creating simulation population considering probability of discontinuing treatments
pop_perc_func<-function(i){
  prob_stop<-0.188/18*3
  
  rx_cycles_table<-data.frame(rx_cycles=seq(1,i,by=1))
  rx_cycles_table<-rx_cycles_table %>% 
    mutate(prob=(1-prob_stop)^(rx_cycles-1)*prob_stop,
           cum_prob=cumsum(prob),
           prop_stop=ifelse(rx_cycles==max(rx_cycles),1-cum_prob+prob,prob))
  
  profile_discontinue_perc<-profile_perc %>% 
    select(age_group,age,sex,apoe,stage,profile_perc) %>%
    rowwise() %>%
    mutate(rx_cycles = list(1:i)) %>% 
    unnest(rx_cycles) %>% 
    left_join(rx_cycles_table %>% 
                select(rx_cycles,prop_stop)) %>% 
    mutate(pop_perc=profile_perc*prop_stop) %>% 
    select(-prop_stop)
  
  return(profile_discontinue_perc)
}

# specify models for each transition, except for transition from MCI to institutionalization
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
                         times = seq(0, 10, by = cycle_length))[[1]][[1]] %>%
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


# cost summary
clustered_boot<-function(data,indices){
  sampled_ids<-unique(data$LOPNR)[indices] 
  resampled_data<-data %>% filter(LOPNR %in% sampled_ids)
  getcost(resampled_data)
}

set.seed(2025)  

getcost(mci_svedem_costs)
mci_costs_ci<-round(t(apply((boot(data=mci_svedem_costs, 
            statistic=clustered_boot, 
            R=10000,
            parallel="multicore", 
            ncpus=18))$t,2,function(x) {c(mean=mean(x),
                                          sd=sd(x),
                                          lb=quantile(x,probs=0.025),
                                          ub=quantile(x,probs=0.975))})),0)
mci_costs_ci<-mci_costs_ci[,3:4]

getcost(costwide)
ad_costs_ci<-round(t(apply((boot(data=costwide, 
            statistic=clustered_boot, 
            R=10000,
            parallel="multicore", 
            ncpus=18))$t,2,function(x) {c(mean=mean(x),
                                          sd=sd(x),
                                          lb=quantile(x,probs=0.025),
                                          ub=quantile(x,probs=0.975))})),0)

ad_costs_ci<-ad_costs_ci[,3:4]


# extract utility
source("health utility R script_XX_20250103.R")
utilities<-as.vector(state_utility$mean)
round(as.vector(state_utility$sd),digits = 3)
utilities_lb<-as.vector(state_utility$lb)
utilities_ub<-as.vector(state_utility$ub)

# utility summary
state_utility %>% 
  mutate(range=paste0(format(round(lb,2),nsmall=2),
                      "-",
                      format(round(ub,2),nsmall=2))) %>% 
  select(state,range)


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
runmodel_cust<-function(rx_cycles, 
                        rx_rr, # treatment effect in risk ratio
                        rx_inst, # whether to treat institutionalized people (yes/no)
                        wane_fc, # treatment waning factor
                        age_init, 
                        sex,
                        apoegeno,
                        start_state,
                        cost_mci,
                        cost_ad,
                        r_cost,
                        r_health) # discount rate)
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
  
  cost_mci_comm=cost_mci[1],
  cost_mci_inst=cost_mci[2],
  
  cost_mild=cost_ad[1],
  cost_mod=cost_ad[2],
  cost_sev=cost_ad[3],

  cost_mild_inst=cost_ad[4],
  cost_mod_inst=cost_ad[5],
  cost_sev_inst=cost_ad[6],
  
  disc_cost=exp(-r_cost*(model_time-1)*cycle_length), #continuously compounded discounting factor
  disc_health=exp(-r_health*(model_time-1)*cycle_length),
  
  #disc=1/(1+r)^(model_time-1)*cycle_length
  
  rx_eff_log=log(rx_rr), # specify log HR to facilitate later PSA
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
  
  
  # Incidence of ARIA-E in the trial
  inc_mild_ariae=ifelse(apoegeno==0,0.014,0.013),
  inc_sev_ariae=ifelse(apoegeno==0,0.002,0.008),
  
  haz_mild_ariae=dispatch_strategy(
    standard=0,
    rx=case_when(model_time==1~-log(1-inc_mild_ariae*0.7), # hazard for first circle
                 model_time==2~-log(1-inc_mild_ariae*0.2), # hazard for second circle
                 model_time>=3&model_time<=rx_cycles~-log(1-inc_mild_ariae*0.1)/5, # constant hazard the next circles
                 T~0)),
  
  prob_mild_ariae=1-exp(-haz_mild_ariae), # transition probability
  
  haz_severe_ariae=dispatch_strategy(
    standard=0,
    rx=case_when(model_time==1~-log(1-inc_sev_ariae*0.7), # hazard for first circle
                 model_time==2~-log(1-inc_sev_ariae*0.2), # hazard for second circle
                 model_time>=3&model_time<=rx_cycles~-log(1-inc_sev_ariae*0.1)/5, # constant hazard the next circles
                 T~0)),
  
  prob_severe_ariae=1-exp(-haz_severe_ariae),
  
  # Incidence of isolated ARIA-H in the trial
  inc_mild_ariah=ifelse(apoegeno==0,0.003,0.004),
  inc_sev_ariah=ifelse(apoegeno==0,0.004,0.003),
  
  haz_mild_ariah=dispatch_strategy(
    standard=0, 
    rx=case_when(model_time<=rx_cycles~-log(1-inc_mild_ariah)/7, # constant hazard 
                 T~0)),
  
  prob_mild_ariah=1-exp(-haz_mild_ariah), # transition probability per cycle
  
  haz_severe_ariah=dispatch_strategy(
    standard=0,
    rx=case_when(model_time<=rx_cycles~-log(1-inc_sev_ariah)/7, # constant hazard 
                 T~0)),
  
  prob_severe_ariah=1-exp(-haz_severe_ariah),
  
  # costs due to ARIA
  cost_mild_aria=8576,
  cost_severe_aria=107814,
  
  
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
  
  
  # transition probability for MCI to institutionalization 
  # from a 6-year study with a population with mean age of 74.5
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
  t12=1-exp(-(get_haz(12,age_init,sex,apoegeno))[model_time,1]*exp(log(rX)*rx_inst)),
  t13=1-exp(-(get_haz(13,age_init,sex))[model_time,1]),
  t14=1-exp(-(get_haz(14,age_init,sex))[model_time,1]*exp(log(rX)*rx_inst)),
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
    *cycle_length*disc_health, 
    
    cost = ((cost_mci_comm*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_mci_comm*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_mci_comm*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_mci_comm*cycle_length+cost_admin+cost_irr)*prob_irr)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length, # treatment time/cycle = 1/(no.cycles/year) - full treatment time every cycle
    
    cost_care=cost_mci_comm*cycle_length*disc_cost,
    
    cost_admin=cost_admin*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*disc_cost
)
  Mild = define_state(
    utility = (utilities_mild*(1-prob_severe_aria-prob_irr)+
                 (utilities_mild+disutility_aria)*prob_severe_aria+
                 (utilities_mild+disutility_irr)*prob_irr)
    *cycle_length*disc_health, 
    
    cost = ((cost_mild*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_mild*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_mild*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_mild*cycle_length+cost_admin+cost_irr)*prob_irr)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length,
    
    cost_care=cost_mild*cycle_length*disc_cost,
    
    cost_admin=cost_admin*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*disc_cost
    )
  Moderate = define_state(
    utility = utilities_mod*cycle_length*disc_health,
    cost = cost_mod*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_mod*cycle_length*disc_cost,
    cost_admin=0,
    cost_aria=0
    )
  Severe = define_state(
    utility = utilities_sev*cycle_length*disc_health,
    cost = cost_sev*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_sev*cycle_length*disc_cost,
    cost_admin=0,
    cost_aria=0
    )
  MCI_inst = define_state(
    utility = (utilities_mci_inst*(1-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
                 (utilities_mci_inst+disutility_aria)*prob_severe_aria*rx_inst+
                 (utilities_mci_inst+disutility_irr)*prob_irr*rx_inst)
    *cycle_length*disc_health, 
    
    cost = ((cost_mci_inst*cycle_length+cost_admin*rx_inst)*(1-prob_mild_aria*rx_inst-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
              (cost_mci_inst*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria*rx_inst+
              (cost_mci_inst*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria*rx_inst+
              (cost_mci_inst*cycle_length+cost_admin+cost_irr)*prob_irr*rx_inst)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length*rx_inst,
    
    cost_care=cost_mci_inst*cycle_length*disc_cost,
    
    cost_admin=cost_admin*rx_inst*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*rx_inst*disc_cost
    )
  Mild_inst = define_state(
    utility = (utilities_mild_inst*(1-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
                 (utilities_mild_inst+disutility_aria)*prob_severe_aria*rx_inst+
                 (utilities_mild_inst+disutility_irr)*prob_irr*rx_inst)
    *cycle_length*disc_health, 
    
    cost = ((cost_mild_inst*cycle_length+cost_admin*rx_inst)*(1-prob_mild_aria*rx_inst-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
              (cost_mild_inst*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria*rx_inst+
              (cost_mild_inst*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria*rx_inst+
              (cost_mild_inst*cycle_length+cost_admin+cost_irr)*prob_irr*rx_inst)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length*rx_inst,
    
    cost_care=cost_mild_inst*cycle_length*disc_cost,
    
    cost_admin=cost_admin*rx_inst*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*rx_inst*disc_cost
    )
  Moderate_inst = define_state(
    utility = utilities_mod_inst*cycle_length*disc_health,
    cost = cost_mod_inst*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_mod_inst*cycle_length*disc_cost,
    cost_admin=0,
    cost_aria=0
    )
  Severe_inst = define_state(
    utility = utilities_sev_inst*cycle_length*disc_health,
    cost = cost_sev_inst*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_sev_inst*cycle_length*disc_cost,
    cost_admin=0,
    cost_aria=0
    )
  Death = define_state(
    utility = 0,
    cost = 0,
    rx.time=0,
    cost_care=0,
    cost_admin=0,
    cost_aria=0
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
  } else{
    c(0, 1, 0, 0, 0, 0, 0, 0, 0)
  }
  
  res_mod <- run_model(
    standard = strat_standard,
    rx=strat_standard,
    parameters = param,
    init=init_state,
    cycles = cycles,
    cost = cost,
    effect = utility,
    method='beginning'
  ) 
  
  # summarise results from the hemmod model
  temp_ly_soc<-res_mod$eval_strategy_list$standard$counts %>% summarize_all(~sum(.))*cycle_length
  temp_ly_rx<-res_mod$eval_strategy_list$rx$counts %>% summarize_all(~sum(.))*cycle_length
  
  temp_qaly_soc<-sum(res_mod$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(res_mod$eval_strategy_list$rx$values$utility)
  
  temp_costs_soc<-sum(res_mod$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(res_mod$eval_strategy_list$rx$values$cost)
  
  temp_costs_care_soc<-sum(res_mod$eval_strategy_list$standard$values$cost_care)
  temp_costs_care_rx<-sum(res_mod$eval_strategy_list$rx$values$cost_care)
  
  temp_costs_admin_soc<-sum(res_mod$eval_strategy_list$standard$values$cost_admin)
  temp_costs_admin_rx<-sum(res_mod$eval_strategy_list$rx$values$cost_admin)
  
  temp_costs_aria_soc<-sum(res_mod$eval_strategy_list$standard$values$cost_aria)
  temp_costs_aria_rx<-sum(res_mod$eval_strategy_list$rx$values$cost_aria)
  
  temp_rxtime<-sum(res_mod$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  res_mod_summary<-data.frame(ly_mci_soc=temp_ly_soc[["MCI"]],ly_mci_rx=temp_ly_rx[["MCI"]],
                              ly_mild_soc=temp_ly_soc[["Mild"]],ly_mild_rx=temp_ly_rx[["Mild"]],
                              ly_mod_soc=temp_ly_soc[["Moderate"]],ly_mod_rx=temp_ly_rx[["Moderate"]],
                              ly_sev_soc=temp_ly_soc[["Severe"]],ly_sev_rx=temp_ly_rx[["Severe"]],
                              ly_mci_inst_soc=temp_ly_soc[["MCI_inst"]],ly_mci_inst_rx=temp_ly_rx[["MCI_inst"]],
                              ly_mild_inst_soc=temp_ly_soc[["Mild_inst"]],ly_mild_inst_rx=temp_ly_rx[["Mild_inst"]],
                              ly_mod_inst_soc=temp_ly_soc[["Moderate_inst"]],ly_mod_inst_rx=temp_ly_rx[["Moderate_inst"]],
                              ly_sev_inst_soc=temp_ly_soc[["Severe_inst"]],ly_sev_inst_rx=temp_ly_rx[["Severe_inst"]],
                              ly_death_soc=temp_ly_soc[["Death"]],ly_death_rx=temp_ly_rx[["Death"]],
                              qaly_soc=temp_qaly_soc,qaly_rx=temp_qaly_rx,
                              costs_soc=temp_costs_soc,costs_rx=temp_costs_rx,
                              costs_care_soc=temp_costs_care_soc,costs_care_rx=temp_costs_care_rx,
                              costs_admin_soc=temp_costs_admin_soc,costs_admin_rx=temp_costs_admin_rx,
                              costs_aria_soc=temp_costs_aria_soc,costs_aria_rx=temp_costs_aria_rx,
                              rxtime=temp_rxtime)
  
  # returns a list of a summary of the model and the model itself
  return(list(res_mod_summary,res_mod))
}


# 4. Create functions to facilitate summarising results ----
# create a function to summarise results
ce_summary_func<-function(results){
  results<-results %>%
    mutate(across(ly_mci_soc:rxtime, 
                  .fns = list(w = ~ . * results$perc),
                  .names = "{col}_{fn}" ) )
  
  return(results)
}


# create a function to organize results that can be exported to excel
ce_summary_excel_func<-function(summary){
  results_excel<-summary %>% 
    select(rx_cycles,rxtime,age_group,age,sex,apoe,stage,qaly_soc:costs_rx) %>% 
    mutate(sex=factor(sex,levels=0:1,labels=c("Male","Female")),
           apoe=factor(apoe,levels=0:1,labels=c("Non-carrier","Heterozygotes")),
           stage=factor(stage,levels=0:1,labels=c("MCI","Mild AD")),
           dqaly=qaly_rx-qaly_soc,
           dcost=costs_rx-costs_soc) %>% 
    set_names(c("Treatment cycle","Treatment duration (years)","Age group","Age","Sex","APOE genotype",
                "AD stage","QALY-SOC","QALY-treatment",
                "Costs-SOC","Costs-treatment",
                "QALY gain","Incremental costs"))
  return(results_excel)
}



# 5. Base-case model (waning factor = 1, i.e., no treatment effect after treatment stops)-----
# sex: 0 = male, 1 = female
# apoegeno: 0 = Non-carrier, 1 = heterozygote
# start_state: 0 = MCI, 1 = Mild AD


# 5.1. Calculate cost-effectiveness price at a WTP of 1000000 ----
QALYval=1000000
pop_perc_rx3yr<-pop_perc_func(12)

# setting up parallel computing
ncpus = parallel::detectCores()-2
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod)))
clusterExport(cl, c(ls(all.names = TRUE)))



base_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  
  pop_perc<-pop_perc_rx3yr
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  temp_df<-temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc)
  
  return(list(model=list(temp_model[[2]],
                         rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                         sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc),df=temp_df))
})

stopCluster(cl)

base_mod_summary<-ce_summary_func(bind_rows(lapply(base_mod,`[[`,"df")))
base_mod_excel<-ce_summary_excel_func(base_mod_summary)


10-sum(base_mod_summary$ly_death_soc_w)
10-sum(base_mod_summary$ly_death_rx_w)

sum(base_mod_summary$costs_care_soc_w)
sum(base_mod_summary$costs_care_rx_w)

sum(base_mod_summary$costs_admin_soc_w)
sum(base_mod_summary$costs_admin_rx_w)

sum(base_mod_summary$costs_aria_soc_w)
sum(base_mod_summary$costs_aria_rx_w)

qaly=sum(base_mod_summary$qaly_rx_w)-sum(base_mod_summary$qaly_soc_w)
cost=sum(base_mod_summary$costs_rx_w)-sum(base_mod_summary$costs_soc_w)
NMB=qaly*QALYval-cost
bc_price=NMB/sum(base_mod_summary$rxtime_w)
bc_price




# define factor labels and legend colors for later plotting
state_label<-c("MCI due to AD","Mild AD","Moderate AD","Severe AD",
         "Institutionalized MCI due to AD",
         "Institutionalized mild AD","Institutionalized moderate AD",
         "Institutionalized severe AD","Death")

cbp <- c("#999999","#CC79A7","#E69F00","#F0E442",
         "#D55E00","pink","#56B4E9","#0072B2",
         "#009E73")

# 5.2. Generate plots for sojourn time ----
# Plot - sojourn times by state
sojourn<-base_mod_summary %>% 
  select(matches("^ly_.*_w$")) %>% 
  summarize_all(~sum(.)) %>% 
  pivot_longer(
    cols = matches("^ly_.*_w$"),
    names_to = c("prefix", "state", "strategy", "suffix"),
    names_pattern = "^(ly)_(.*)_(.*)_(w)$", 
    values_to = "years"
  ) %>%
  select(-prefix, -suffix) %>% 
  mutate(state=factor(state, levels=rev(c("mci","mild","mod","sev",
                                      "mci_inst","mild_inst","mod_inst","sev_inst","death")),
                      labels = rev(state_label)),
         strategy=factor(strategy,levels=c("soc","rx"),
                         labels=c("SoC","Lecanemab+SoC")))

  
sojourn_total<-ggplot(data = sojourn %>% 
         group_by(strategy,state), 
       aes(y = years, x = strategy, fill = state)) +
  geom_col(position = 'stack') +
  geom_text(aes(label = ifelse(years >= 0.4, round(years, 2), "")), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3) + 
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Group", y = "Length of stay in each state (years)", fill = "State") +
  theme_minimal(base_size = 15) + 
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")  
  ) +
  scale_fill_manual(values = cbp) + 
  guides(fill = guide_legend(reverse = TRUE)) 
  

# Plot - differences in sojourn times between treatment groups by states
sojourn_diff<-sojourn %>% 
  pivot_wider(id_cols=state, names_from='strategy', values_from='years') %>%
  mutate(years=(`Lecanemab+SoC`-SoC), 
         strategy='difference') %>% select(state, years, strategy) %>% 
  mutate(month=years*12) 

sojourn_diff_plot<-ggplot(data=sojourn_diff, aes(y=month, x=state, fill=state)) + 
  geom_col(position='stack')+coord_flip() +
  coord_flip() +
  scale_y_continuous(limits = c(min(floor(sojourn_diff$month)),max(ceiling(sojourn_diff$month))),
                     breaks = seq(min(floor(sojourn_diff$month)),max(ceiling(sojourn_diff$month)), by = 1), expand = c(0.1, 0.1)) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  labs(x = "State", y = "Differences in lengths of stay in each state (months)", fill = "State") +
  theme_minimal(base_size = 15) + 
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = cbp) + 
  guides(fill = guide_legend(reverse = TRUE)) 

pdf(file = "sojourn time.pdf",height = 7,width = 8)
ggarrange(sojourn_total,sojourn_diff_plot,nrow = 2,common.legend = F,labels = c("A","B"),heights = c(0.7,1))
dev.off()


# Plot - transition probability by treatment groups
base_mod_models<-lapply(base_mod,`[[`,"model")
cohdist<-lapply(1:length(base_mod_models),function(i){
  
  temp_mod<-base_mod_models[[i]][[1]]
  temp_perc<-base_mod_models[[i]]$perc
  
  standard_counts<-temp_mod$eval_strategy_list$standard$counts*temp_perc
  rx_counts<-temp_mod$eval_strategy_list$rx$counts*temp_perc
  
  temp_cohdist <- bind_rows(
    standard_counts %>% mutate(group = "SoC"),
    rx_counts %>% mutate(group = "Treatment"))
  
  temp_cohdist$t=rep(1:cycles,2)
  
  return(temp_cohdist)
}) %>% 
  bind_rows() %>% 
  group_by(group,t) %>%
  summarize(across(everything(), ~ sum(.x))) %>% 
  ungroup
  
# export state trace to excel

  


coh.inst<-pivot_longer(cohdist,cols=MCI:Death,names_to='State') %>% 
  mutate(years=(t-1)*0.25,
         State=factor(State, levels=rev(statenames),labels = rev(state_label)))


coh.inst %>%
  mutate(value=round(value,3)) %>% 
  pivot_wider(
    names_from = State,
    values_from = value
  ) %>% 
  select(-t) %>% 
  rename(Group=group,
         Years=years)

ggplot(aes(x = years, y = value, group = interaction(State, group)), 
       data = coh.inst) +
  geom_line(aes(linetype = group, color = State)) +
  geom_point(aes(color = State)) +
  scale_x_continuous(breaks = seq(0, 10, by = 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
  scale_linetype_manual(values=c("solid","dashed"))+
  labs(x = "Time horizon (years)", y = "State occupation probability", color = "State",
       linetype="Group") +
  theme_minimal(base_size = 15) + 
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")  
  ) +
  scale_color_manual(values = cbp) + 
  guides(color = guide_legend(reverse = TRUE)) 



# 5.3. Generate plots for costs by treatment groups ----
costplot<-lapply(1:length(base_mod_models), function(i){
  
  temp_mod<-base_mod_models[[i]][[1]]
  temp_perc<-base_mod_models[[i]]$perc
  
  costs_soc<-data.frame(Markov.cycle=temp_mod$eval_strategy_list$standard$values$model_time, 
                             cost=(temp_mod$eval_strategy_list$standard$values$cost)*temp_perc,
                        cost_care=(temp_mod$eval_strategy_list$standard$values$cost_care)*temp_perc,
                        cost_admin=(temp_mod$eval_strategy_list$standard$values$cost_admin)*temp_perc,
                        cost_aria=(temp_mod$eval_strategy_list$standard$values$cost_aria)*temp_perc)
  
  costs_rx<-data.frame(Markov.cycle=temp_mod$eval_strategy_list$rx$values$model_time, 
                       cost=(temp_mod$eval_strategy_list$rx$values$cost)*temp_perc,
                       cost_care=(temp_mod$eval_strategy_list$rx$values$cost_care)*temp_perc,
                       cost_admin=(temp_mod$eval_strategy_list$rx$values$cost_admin)*temp_perc,
                       cost_aria=(temp_mod$eval_strategy_list$rx$values$cost_aria)*temp_perc)
  
  temp_costs <- bind_rows(
    costs_soc %>% mutate(group = "soc"),
    costs_rx %>% mutate(group = "rx")
  )
  return(temp_costs)
}) %>% 
  bind_rows() %>% 
  group_by(group,Markov.cycle) %>%
  summarize(across(everything(), ~ sum(.x))) %>% 
  ungroup %>% 
  mutate(years=Markov.cycle*0.25)


costplot_incr<-costplot %>% 
  pivot_wider(names_from=group,
              values_from=c(cost,cost_care,cost_admin,cost_aria),names_sep="_") %>% 
  mutate(incr_cost_total=cost_rx-cost_soc,
         incr_cost_care=cost_care_rx-cost_care_soc,
         incr_cost_admin=cost_admin_rx-cost_admin_soc,
         incr_cost_aria=cost_aria_rx-cost_aria_soc) %>% 
  select(years,contains("incr_")) %>% 
  pivot_longer(cols=starts_with("incr_"),
               names_to ="type",names_prefix = "incr_cost_") %>% 
  mutate(type=factor(type,levels=c("total","care","aria","admin"),
                     labels=c("Total","Background costs of care",
                              "Adverse event-related costs",
                              "Treatment administration and monitoring"))) %>% 
  mutate(dire=as.factor(ifelse(cost>=0,1,0)))


plot_costs<-ggplot(aes(y=cost, x=years, color=group), data=costplot) + 
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks = seq(0, 10, by = 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(floor(min(costplot$cost)/10000)*10000, 
                                  ceiling(max(costplot$cost)/10000)*10000, by = 10000),
                     limits = c(floor(min(costplot$cost)/10000)*10000,
                                ceiling(max(costplot$cost)/10000)*10000),
                     labels = scales::label_number(accuracy = 1, big.mark = ",")) +
  labs(x = "Time horizon (years)", y = "Costs (2023 SEK)", color = "Strategy",
       linetype="Group") +
  scale_color_manual(values = c("#E69F00","#56B4E9"),labels=c("Lecanemab+SoC","SoC")) +
  theme_minimal(base_size = 15) + 
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom")

plot_costs_incr<-ggplot(costplot_incr %>% 
         filter(type!="Total"),aes(x=years,y=value,fill=type)) +
  geom_bar(stat = "identity")+
  scale_x_continuous(breaks = seq(0, 10, by = 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(floor(min(costplot_incr$value)/10000)*10000, 
                                  ceiling(max(costplot_incr$value)/10000)*10000, by = 10000),
                     limits = c(floor(min(costplot_incr$value)/10000)*10000,
                                ceiling(max(costplot_incr$value)/10000)*10000),
                     labels = scales::label_number(accuracy = 1, big.mark = ",")) +
  labs(x = "Time horizon (years)", y = "Incremental costs (2023 SEK)", fill = "Type of costs",
       linetype="Group") +
  scale_fill_manual(breaks=c("Background costs of care",
                               "Treatment administration and monitoring",
                               "Adverse event-related costs"),
                      values = c("#009E73","#0072B2","pink"))+
  guides(fill=guide_legend(nrow=3))+
  theme_minimal(base_size = 15) + 
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom")


pdf(file = "costs_overtime.pdf",height = 6,width = 12)
ggarrange(plot_costs+theme(aspect.ratio = 1),
          plot_costs_incr+theme(aspect.ratio = 1),
          nrow = 1,labels = c("A","B"),
          heights = c(0.8, 1),widths = c(0.8,1))
dev.off()

# 6. One-way sensitivity analysis ----
# ************************************************************************************************
# 1. Varying treatment effects (varying HR)
# 2. Uncertainty around Utility (varying utility by state)
# 3. Varying discount rates for costs and QALYs: 0 and 0.05, same for costs and QALYs
# ************************************************************************************************
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod)))
clusterExport(cl, c("runmodel_cust",
                    ls(pattern = c("utiliti|model_trans|cycle|pop_perc")),
                    "get_haz","statenames","getcost","profile_perc",
                    "mci_svedem_costs","costwide","mci_costs_ci","ad_costs_ci"))

# 6.1. Define DSA parameters through define_dsa() ----
dsa_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)[[2]]
  
  dsa_para_mod <- define_dsa(
    rx_eff, 0.572, 0.833,
    
    utilities_mci, utilities_lb[1], utilities_ub[1],
    utilities_mild, utilities_lb[2], utilities_ub[2],
    utilities_mod, utilities_lb[3], utilities_ub[3],
    utilities_sev, utilities_lb[4], utilities_ub[4],
    utilities_mci_inst, utilities_lb[5], utilities_ub[5],
    utilities_mild_inst, utilities_lb[6], utilities_ub[6],
    utilities_mod_inst, utilities_lb[7], utilities_ub[7],
    utilities_sev_inst, utilities_lb[8], utilities_ub[8],
    
    cost_mci_comm, mci_costs_ci[1,1], mci_costs_ci[1,2],
    cost_mild, ad_costs_ci[1,1], ad_costs_ci[1,2],
    cost_mod, ad_costs_ci[2,1], ad_costs_ci[2,2],
    cost_sev, ad_costs_ci[3,1], ad_costs_ci[3,2],
    
    cost_mci_inst, mci_costs_ci[2,1], mci_costs_ci[2,2],
    cost_mild_inst, ad_costs_ci[4,1], ad_costs_ci[4,2],
    cost_mod_inst, ad_costs_ci[5,1], ad_costs_ci[5,2],
    cost_sev_inst, ad_costs_ci[6,1], ad_costs_ci[6,2]
  )
  
  temp_dsa<-run_dsa(
    model=temp_model,
    dsa=dsa_para_mod
  )
  
  temp_dsa_summary<-temp_dsa[["dsa"]] %>% 
    mutate(.par_value=format(round(.par_value_eval,2),nsmall=2),
           names=case_match(.par_names, "rx_eff"~"Treatment effect",
                            "utilities_mci"~"Utility of MCI",
                            "utilities_mild"~"Utility of mild AD",
                            "utilities_mod"~"Utility of moderate AD",
                            "utilities_sev"~"Utility of severe AD",
                            "utilities_mci_inst"~"Utility of institutionalized MCI",
                            "utilities_mild_inst"~"Utility of institutionalized mild AD",
                            "utilities_mod_inst"~"Utility of institutionalized moderate AD",
                            "utilities_sev_inst"~"Utility of institutionalized severe AD",
                            "cost_mci_comm"~"Annual costs of care of MCI",
                            "cost_mild"~"Annual costs of care of mild AD",
                            "cost_mod"~"Annual costs of care of moderate AD",
                            "cost_sev"~"Annual costs of care of severe AD",
                            "cost_mci_inst"~"Annual costs of care of institutionalized MCI",
                            "cost_mild_inst"~"Annual costs of care of institutionalized mild AD",
                            "cost_mod_inst"~"Annual costs of care of institutionalized moderate AD",
                            "cost_sev_inst"~"Annual costs of care of institutionalized severe AD"))
  
  temp_dsa_summary<-temp_dsa_summary %>% 
    filter(.strategy_names=="standard") %>% 
    dplyr::select(.cost,.effect,names,.par_value) %>% 
    rename(cost_soc=.cost,
           effect_soc=.effect) %>% 
    left_join(temp_dsa_summary %>% 
                filter(.strategy_names=="rx") %>% 
                dplyr::select(.cost,.effect,names,.par_value,rx.time) %>% 
                rename(cost_rx=.cost,
                       effect_rx=.effect)) %>% 
    mutate(dcost=cost_rx-cost_soc,
           effect=effect_rx-effect_soc) %>% 
    cbind(rx_cycles=rx_cycles_i,
          age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc)
  
  return(temp_dsa_summary)
}) %>% bind_rows()



dsa_mod_summary<-dsa_mod %>% 
  group_by(names,.par_value) %>% 
  summarise(cost_rx=sum(cost_rx*perc),
            cost_soc=sum(cost_soc*perc),
            effect_rx=sum(effect_rx*perc),
            effect_soc=sum(effect_soc*perc),
            rx.time=sum(rx.time*perc)) %>% 
  ungroup %>% 
  mutate(NMB=QALYval*(effect_rx-effect_soc)-(cost_rx-cost_soc),
         price=NMB/rx.time,
         price_diff=(price-bc_price)/bc_price)

# 6.2. Varying discounting rate ----
# discount rate = 0
disc0_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0,
                            r_health=0)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})

disc0_mod_summary<-ce_summary_func(bind_rows(disc0_mod)) %>% 
  summarise(cost_rx=sum(costs_rx*perc),
            cost_soc=sum(costs_soc*perc),
            effect_rx=sum(qaly_rx*perc),
            effect_soc=sum(qaly_soc*perc),
            rx.time=sum(rxtime*perc)) %>% 
  mutate(NMB=QALYval*(effect_rx-effect_soc)-(cost_rx-cost_soc),
         price=NMB/rx.time,
         price_diff=(price-bc_price)/bc_price)


# discount rate = 0.05
disc0.05_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.05,
                            r_health=0.05)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})

stopCluster(cl)

disc0.05_mod_summary<-ce_summary_func(bind_rows(disc0.05_mod)) %>% 
  summarise(cost_rx=sum(costs_rx*perc),
            cost_soc=sum(costs_soc*perc),
            effect_rx=sum(qaly_rx*perc),
            effect_soc=sum(qaly_soc*perc),
            rx.time=sum(rxtime*perc)) %>% 
  mutate(NMB=QALYval*(effect_rx-effect_soc)-(cost_rx-cost_soc),
         price=NMB/rx.time,
         price_diff=(price-bc_price)/bc_price)



# 6.3. Tornado plot ----  
dsa_plot_data<-dsa_mod_summary %>% 
  mutate(.par_value=as.numeric(.par_value)) %>% 
  rbind(disc0_mod_summary %>% 
          mutate(names="Discount rate for costs and health effects",
                 .par_value=0)) %>% 
  rbind(disc0.05_mod_summary %>% 
          mutate(names="Discount rate for costs and health effects",
                 .par_value=0.05)) %>% 
  group_by(names) %>% 
  mutate(level = as.factor(.par_value == max(.par_value)),
         max_absdiff = max(abs(price_diff))) %>%
  ungroup %>% 
  arrange(desc(max_absdiff))

factor_level <- dsa_plot_data %>%
  group_by(names) %>%
  slice_max(abs(price_diff)) %>%
  ungroup() %>%
  distinct(names, .keep_all = TRUE) %>%
  arrange(abs(price_diff)) %>%
  pull(names)

dsa_plot_data$names<-factor(dsa_plot_data$names,levels=factor_level)


pdf(file = "tonado.pdf",width = 12,height = 6)
ggplot(dsa_plot_data, aes(factor(names), price_diff, fill = level)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(y = price_diff, label = round(price_diff*100,0)),
            hjust = ifelse(dsa_plot_data$price_diff > 0, -0.2, 1.2),
            size = 3, color = "black") +
  coord_flip()+
  scale_y_continuous(labels = scales::label_number(scale = 100),
                     breaks = seq(-1.6,1.4,0.2)) +
  scale_fill_manual(values = c("#56B4E9","#E69F00"),
                    labels=c("Lower bound of the range","Upper bound of the range")) +
  labs(y = "Difference from the base-case threshold annual drug price of lecanemab (%)",
       x = "Parameter",
       fill = "Parameter value") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size=10),
        axis.text = element_text(size = 10),
        axis.ticks.length = unit(3, "pt"))
dev.off()



# 7. Bootstrapping Probabilistic sensitivity analysis -----
# ************************************************************************************************
# 1. Uncertainty around input parameters considered in the one-way sensitivity analyses
# 2. Uncertainty around transition probabilities estimated from individual-level data
# 3. Costs by disease states
# ************************************************************************************************
source("HE model_PSA_function_XX_202503.R")



n_sim<-10000 # add a sequence indicating number of simulations

psa_pop<-data.frame(rx_cycles=rep(12,4),
                    age=rep(70,4),
                    apoe=c(0,0,1,1),
                    sex=c(0,1,0,1),
                    stage=rep(0,4))

sim_prof<-expand.grid(sim=1:n_sim,prof=1:nrow(psa_pop))

ncpus = parallel::detectCores()-1
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod),
                   library(flexsurv),library(parallel)))
clusterExport(cl, c("psa_func",
                    ls(pattern = c("utiliti|model_trans|cycle|_sd")),
                    "mci_ad_model","mci_death_model","msm",
                    "transmat","timepoint_seq","psa_pop",
                    "boots_gethaz_mci","boots_gethaz_msm",
                    "get_haz","statenames","getcost",
                    "mci_svedem_costs","costwide","sim_prof"))




sim_psa<-parLapply(cl,1:nrow(sim_prof),function(n) {

  j<-sim_prof$sim[n]
  i<-sim_prof$prof[n]
  
    sim.seed<-2025 + j
    rx_cycles_i<-psa_pop[i,]$rx_cycles
    age_i<-psa_pop[i,]$age
    sex_i<-psa_pop[i,]$sex
    apoe_i<-psa_pop[i,]$apoe
    stage_i<-psa_pop[i,]$stage
    
    psa_summary<-tryCatch({ # catch error message instead of stopping interations
      psa_func(sim.seed,
               rx_cycles_psa=rx_cycles_i, rx_rr_psa=0.69,
               rx_inst_psa=0, wane_fc_psa=1,
               age_init_psa=age_i, sex_psa=sex_i,
               apoegeno_psa=apoe_i, start_state_psa=stage_i,
               cost_mci_psa=getcost(mci_svedem_costs),
               cost_ad_psa=getcost(costwide),
               r_cost_psa=0.03,
               r_health_psa=0.03) %>% 
        cbind(age=age_i,sex=sex_i, apoe=apoe_i, stage=stage_i,
              error_message=NA)
      
    }, error=function(e) {
      data.frame(age=age_i,sex=sex_i, apoe=apoe_i, stage=stage_i,
                 error_message=as.character(e$message)) 
    })
    return(psa_summary %>% cbind(sim=j))
}) %>% bind_rows()



stopCluster(cl)

save(sim_psa,file="sim_psa.RData")

sim_psa_summ<-sim_psa %>% 
  filter(!is.na(cost_soc)) %>% 
  mutate(dcost=cost_rx-cost_soc,
         deffect=effect_rx-effect_soc) %>% 
  select(sex,apoe,dcost,deffect,rx.time,sim) %>% 
  mutate(perc=ifelse(apoe==0,0.41,0.59),
         dcost=dcost*perc,
         deffect=deffect*perc,
         rx.time=rx.time*perc) %>% 
  group_by(sim,sex) %>% 
  summarise(dcost=sum(dcost),
            deffect=sum(deffect),
            rx.time=sum(rx.time)) %>% 
  ungroup %>% 
  group_by(sex) %>% 
  mutate(mean_cost=mean(dcost),
         mean_effect=mean(deffect)) %>% 
  ungroup

pdf(file = "psa plot.pdf",width = 8,height = 6)
ggplot(data=sim_psa_summ %>% 
                   mutate(Sex=factor(sex,levels=c(0,1),
                                     labels=c("Males","Females"))),
                 aes(x=deffect,y=dcost,color=Sex))+
  geom_vline(xintercept = 0,color="black",linetype="solid")+
  geom_hline(yintercept = 0,color="black",linetype="solid")+
  geom_point(alpha=0.4) +
  scale_x_continuous(breaks = seq(-0.1,0.6,by=0.1))+
  scale_y_continuous(breaks = seq(floor(min(sim_psa_summ$dcost)/100000)*100000-200000, 
                                  ceiling(max(sim_psa_summ$dcost)/100000)*100000, by = 100000),
                     limits = c(floor(min(sim_psa_summ$dcost)/100000)*10000-250000,
                                ceiling(max(sim_psa_summ$dcost)/100000)*100000),
                     labels = scales::label_number(accuracy = 1, big.mark = ",")) +
  scale_color_manual(values = c("#56B4E9","pink")) +
  labs(x = "Incremental QALY",
       y = "Incremental costs (2023 SEK)") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_point(data = sim_psa_summ %>% 
               mutate(Sex = factor(sex, levels = c(0, 1),
                                   labels = c("Males", "Females"))) %>%
               distinct(Sex, mean_effect, mean_cost),
             aes(x = mean_effect, y = mean_cost),
             fill=c("#56B4E9","pink"),
             color="black",size = 3,shape=23)
dev.off()



# get CI of incremental costs, QALYs, and prices
sim_psa_summ %>% 
  group_by(sex)%>%
  summarise(across(dcost:rx.time, 
                   list(mean = ~mean(. , na.rm = TRUE),
                        lb = ~quantile(., 0.025, na.rm = TRUE),
                        ub = ~quantile(., 0.975, na.rm = TRUE)), 
                   .names = "{.col}_{.fn}")) %>%
  ungroup %>% 
  mutate(price_mean=(QALYval*deffect_mean-dcost_mean)/rx.time_mean,
         price_lb=(QALYval*deffect_lb-dcost_lb)/rx.time_lb,
         price_ub=(QALYval*deffect_ub-dcost_ub)/rx.time_ub) %>% 
  mutate(dcost=paste0(round(dcost_mean, 0), 
                      " (", 
                      round(dcost_lb, 0), ",", 
                      round(dcost_ub, 0), ")"),
         deffect=paste0(round(deffect_mean, 2), 
                      " (", 
                      round(deffect_lb, 2), "-", 
                      round(deffect_ub, 2), ")"),
         price=paste0(round(price_mean, 0), 
                      " (", 
                      round(price_lb, 0), "-", 
                      round(price_ub, 0), ")")) %>% 
  select(sex,dcost,deffect,price) %>% 
  mutate(sex=factor(sex,levels=0:1,labels=c("Male","Female"))) %>% 
  set_names(c("Sex","Incremental costs","Incremental QALYs",
              "Threshold annual drug price"))



# 8. Scenario analysis ----
# ************************************************************************************************
# 1. Assuming effect modifications
# 2. Assuming treating for only 18 months or a sequence of 1:9 years
# 3. Assuming treating for whole time horizon, but with Waning factors = 0.5 when treatment stops
# 4. APOE heterozygotes vs. non-carriers
# 5. Discount rate of 3% for costs and 0 for health effects 
# ************************************************************************************************
# call functions for scenario analyses
source("HE model_scenario analyses_function_XX_202502.R")


# set cluster for parellel computing
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod)))
clusterExport(cl, c("scen_runmodel_cust","runmodel_cust",
                    ls(pattern = c("utiliti|model_trans|cycle|pop_perc")),
                    "get_haz","statenames","getcost","profile_perc",
                    "mci_svedem_costs","costwide"))

# 8.1. Assume effect modifications in subgroups ----
sub_rx_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-scen_runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,wane_fc=1,rx_inst=0,
                                 rx_sex_mod=1,rx_age_mod=1,rx_apoe_mod=1,rx_mci_mod=1,rx_ad_mod=1,
                                 age_init = age_i,sex=sex_i,
                                 apoegeno=apoe_i,start_state=stage_i,
                                 cost_mci=getcost(mci_svedem_costs),
                                 cost_ad=getcost(costwide),
                                 r_cost=0.03,
                                 r_health=0.03)
  
  return(temp_model %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})


sub_rx_mod_summary<-ce_summary_func(bind_rows(sub_rx_mod))
sub_rx_mod_excel<-ce_summary_excel_func(sub_rx_mod_summary)



# 8.2. Assuming treatment for varying durations of treatment ----
rx_duration_mod<-parLapply(cl,c(6,c(4:10)/cycle_length),function(j){
  pop_perc<-pop_perc_func(j)
  
  results<-lapply(1:nrow(pop_perc),function(i){
    rx_cycles_i<-pop_perc[i,]$rx_cycles
    age_group_i<-pop_perc[i,]$age_group
    age_i<-pop_perc[i,]$age
    sex_i<-pop_perc[i,]$sex
    apoe_i<-pop_perc[i,]$apoe
    stage_i<-pop_perc[i,]$stage
    perc<-pop_perc[i,]$pop_perc
    
    temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                              rx_inst=0,wane_fc=1,
                              age_init = age_i,sex=sex_i,
                              apoegeno=apoe_i,start_state=stage_i,
                              cost_mci=getcost(mci_svedem_costs),
                              cost_ad=getcost(costwide),
                              r_cost=0.03,
                              r_health=0.03)
    
    return(temp_model[1][[1]] %>% 
             cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                   sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
  }) %>% bind_rows() %>% cbind(rx_duration=j)
}) %>% bind_rows() 




rx_duration_mod_summary<-ce_summary_func(rx_duration_mod)
rx_duration_mod_excel<-lapply(c(6,c(4:10)/cycle_length), function(i){
  return(ce_summary_excel_func(rx_duration_mod_summary %>% 
                                 filter(rx_duration==i)) %>% 
           cbind(max_rx_cycles=i) %>% 
           select(max_rx_cycles,everything()))
}) %>% bind_rows()



# 8.3. Treatment with waning factor = 0.5, 0.1, or 0.25 ----
wf0.5_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=0.5,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})




wf0.5_mod_summary<-ce_summary_func(bind_rows(wf0.5_mod))
wf0.5_mod_excel<-ce_summary_excel_func(wf0.5_mod_summary)


wf0.1_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=0.1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})


wf0.25_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=0.25,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})


wf0.1_mod_summary<-ce_summary_func(bind_rows(wf0.1_mod))
wf0.1_mod_excel<-ce_summary_excel_func(wf0.1_mod_summary)

wf0.25_mod_summary<-ce_summary_func(bind_rows(wf0.25_mod))
wf0.25_mod_excel<-ce_summary_excel_func(wf0.25_mod_summary)


# 8.4. Assuming a starting population of only APOE non-carriers or only heterozygotes ----
# 8.4.1. Assuming a starting population of only APOE non-carriers ----
apoe0_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr %>% filter(apoe==0)),function(i){
  pop_perc<-pop_perc_rx3yr %>% 
    filter(apoe==0) %>% 
    mutate(pop_perc=pop_perc/sum(pop_perc))
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-scen_runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,wane_fc=1,
                            rx_inst=0,rx_sex_mod=0,rx_age_mod=0,rx_apoe_mod=1,
                            rx_mci_mod=0,rx_ad_mod=0,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})


apoe0_mod_summary<-ce_summary_func(bind_rows(apoe0_mod))
apoe0_mod_excel<-ce_summary_excel_func(apoe0_mod_summary)





# 8.4.2. Assuming a starting population of only APOE heterozygotes ----
apoe1_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr %>% filter(apoe==1)),function(i){
  pop_perc<-pop_perc_rx3yr %>% 
    filter(apoe==1) %>% 
    mutate(pop_perc=pop_perc/sum(pop_perc))
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-scen_runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,wane_fc=1,
                                 rx_inst=0,rx_sex_mod=0,rx_age_mod=0,rx_apoe_mod=1,
                                 rx_mci_mod=0,rx_ad_mod=0,
                                 age_init = age_i,sex=sex_i,
                                 apoegeno=apoe_i,start_state=stage_i,
                                 cost_mci=getcost(mci_svedem_costs),
                                 cost_ad=getcost(costwide),
                                 r_cost=0.03,
                                 r_health=0.03)
  
  return(temp_model %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})


apoe1_mod_summary<-ce_summary_func(bind_rows(apoe1_mod))
apoe1_mod_excel<-ce_summary_excel_func(apoe1_mod_summary)







# 8.5. Discounting rate 0.03 for costs and no discounting for QALY ----
disc_diff_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})




disc_diff_mod_summary<-ce_summary_func(bind_rows(disc_diff_mod))
disc_diff_mod_excel<-ce_summary_excel_func(disc_diff_mod_summary)


# 8.6. Assuming treatments for institutionalized patients ----
rx_inst_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=1,wane_fc=1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})




rx_inst_mod_summary<-ce_summary_func(bind_rows(rx_inst_mod))
rx_inst_mod_excel<-ce_summary_excel_func(rx_inst_mod_summary)



stopCluster(cl)


# 8.7. 4-year treatment with waning factor = 0.5, 0.1, or 0.25 ----
pop_perc_rx4yr<-pop_perc_func(16)


wf0.5_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=0.5,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})




wf0.5_mod_summary<-ce_summary_func(bind_rows(wf0.5_mod))
wf0.5_mod_excel<-ce_summary_excel_func(wf0.5_mod_summary)


wf0.1_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=0.1,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})


wf0.25_mod<-parLapply(cl,1:nrow(pop_perc_rx3yr),function(i){
  pop_perc<-pop_perc_rx3yr
  
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,rx_rr=0.69,
                            rx_inst=0,wane_fc=0.25,
                            age_init = age_i,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r_cost=0.03,
                            r_health=0.03)
  
  return(temp_model[1][[1]] %>% 
           cbind(rx_cycles=rx_cycles_i,age_group=age_group_i,age=age_i,
                 sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc))
})


wf0.1_mod_summary<-ce_summary_func(bind_rows(wf0.1_mod))
wf0.1_mod_excel<-ce_summary_excel_func(wf0.1_mod_summary)

wf0.25_mod_summary<-ce_summary_func(bind_rows(wf0.25_mod))
wf0.25_mod_excel<-ce_summary_excel_func(wf0.25_mod_summary)





# 8.8. Summarise results from scenario analyses ----
scen_sum_function<-function(mod_summary){
  price<-((sum(mod_summary$qaly_rx_w)-sum(mod_summary$qaly_soc_w))*QALYval-
            (sum(mod_summary$costs_rx_w)-sum(mod_summary$costs_soc_w)))/
    sum(mod_summary$rxtime_w)
  
  temp_sum<-mod_summary %>% 
    summarise(soc=10-sum(ly_death_soc_w),
              rx=10-sum(ly_death_rx_w),
              diff=rx-soc) %>% 
    mutate(out="ly") %>% 
    rbind(mod_summary %>% 
            summarise(soc=sum(qaly_soc_w),
                      rx=sum(qaly_rx_w),
                      diff=rx-soc) %>% 
            mutate(out="qaly")) %>% 
    rbind(mod_summary %>% 
            summarise(soc=sum(costs_soc_w),
                      rx=sum(costs_rx_w),
                      diff=rx-soc) %>% 
            mutate(out="costs")) %>% 
    rbind(mod_summary %>% 
            summarise(soc=sum(costs_care_soc_w),
                      rx=sum(costs_care_rx_w),
                      diff=rx-soc) %>% 
            mutate(out="costs_care")) %>% 
    rbind(mod_summary %>% 
            summarise(soc=sum(costs_admin_soc_w),
                      rx=sum(costs_admin_rx_w),
                      diff=rx-soc) %>% 
            mutate(out="costs_admin")) %>% 
    rbind(mod_summary %>% 
            summarise(soc=sum(costs_aria_soc_w),
                      rx=sum(costs_aria_rx_w),
                      diff=rx-soc) %>% 
            mutate(out="costs_aria")) %>% 
    rbind(data.frame(0, sum(mod_summary$rxtime_w),0,"rx.time") %>% 
            set_names(c("soc","rx","diff","out"))) %>% 
    rbind(data.frame(0, price,0,"price") %>% 
            set_names(c("soc","rx","diff","out"))) %>% 
    mutate(across(c(soc, rx, diff), 
                  ~ ifelse(grepl("costs|price",out), 
                           format(round(., 0),nsmall=0,big.mark = ",", scientific = FALSE), 
                           format(round(., 2),nsmall=2,big.mark = ",", scientific = FALSE)))) %>% 
    mutate(soc=ifelse(out %in% c("rx.time", "price"),"-",soc),
           diff=ifelse(out %in% c("rx.time", "price"),"-",diff)) %>% 
    rbind(c("-",round((price-bc_price)/bc_price*100,digits = 2),"-","price_dev"))
  
  return(temp_sum)
}


scen_sum_function(sub_rx_mod_summary)
scen_sum_function(wf0.5_mod_summary)
scen_sum_function(wf0.1_mod_summary)
scen_sum_function(wf0.25_mod_summary)
scen_sum_function(apoe0_mod_summary)
scen_sum_function(apoe1_mod_summary)
scen_sum_function(disc_diff_mod_summary)
scen_sum_function(rx_inst_mod_summary)

scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==6))
scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==16))
scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==20))
scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==24))
scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==28))
scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==32))
scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==36))
scen_sum_function(rx_duration_mod_summary %>% 
                    filter(rx_duration==40))






# export summary to excel
writexl::write_xlsx(list(base_mod=base_mod_excel,subgroup_mod=sub_rx_mod_excel,
                         rx_duration_mod=rx_duration_mod_excel,wf0.5_mod=wf0.5_mod_excel,
                         wf0.1_mod=wf0.1_mod_excel,wf0.25_mod=wf0.25_mod_excel,
                         apoe0_mod=apoe0_mod_excel,apoe1_mod=apoe1_mod_excel,
                         disc_diff_mod=disc_diff_mod_excel,rx_inst_mod=rx_inst_mod_excel),
                    "Model_summary.xlsx")
