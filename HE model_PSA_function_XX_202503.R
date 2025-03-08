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
library(mstate)

# load data
setwd("C:/Users/xinxia/OneDrive - Karolinska Institutet/CSF-registersamkörning/AD DMT HE model/Statistical analyses/New analyses_XX_202412/")

load("costdata_xx.RData")
load("models_xx.RData")



# 1. Construct a joint multistate model to avoid transition probability outside limits ----
transmat<-transMat(x=list(c(2,4,7), 
                          c(3,5,7),
                          c(6,7),
                          c(5,7),
                          c(6,7),
                          c(7),
                          c()),
                   names=c('Mild', 'Moderate', 'Severe', 
                           'Mild_inst', 'Moderate_inst', 
                           'Severe_inst', 'mort'))

msm<-fmsm(final_model_trans1,final_model_trans2,final_model_trans3,
          final_model_trans4,final_model_trans5,final_model_trans6,
          final_model_trans7,final_model_trans8,final_model_trans9,
          final_model_trans10,final_model_trans11,final_model_trans12,
          final_model_trans13,trans = transmat)


# specify time points to predict cumulative hazards
timepoint_seq<-seq(0, 10, by = 0.25)

# cumulative hazards for transitions from MCI to next stages
boots_gethaz_mci<-function(model, age_init, sex, apoe) {
  
  sex<-ifelse(sex==0,"MALE","FEMALE")
  
  newdata<-data.frame(AGE = age_init, SEX = sex, NACCNE4S = apoe)
  
  summfn<-function(x){
    cumhaz<-predict(x, type = "cumhaz", 
                    newdata = newdata,
                    times = timepoint_seq)
    list(cumhaz)
  }
  
  haz<-t(bootci.fmsm(x=model, B=1, fn=summfn) %>% 
           as.data.frame() %>% 
           select(V42:V82) %>% 
           slice_head(n=1)) %>% 
    as.data.frame() %>% 
    set_names("cumhaz") %>% 
    mutate(haz=cumhaz-lag(cumhaz)) %>% 
    filter(!is.na(haz)) %>% 
    select(haz)
  
  return(haz)
}

# cumulative hazards for transitions from AD dementia to next stages using the composed msm model
boots_gethaz_msm<-function(age_init, sex) {
  
  sex<-ifelse(sex==0,"MALE","FEMALE")
  newdata<-data.frame(AGE = age_init, SEX = sex)
  
  summfn<-function(x){
    cumhaz<-(msfit.flexsurvreg(x, t = timepoint_seq, trans = transmat,
                               newdata=newdata))$Haz
    list(cumhaz)
  }
  
  haz<-matrix(bootci.fmsm(x=msm, B=1, fn=summfn,sample = T)[1,], 
                     ncol = 3, byrow = FALSE) %>% 
    as.data.frame() %>% 
    set_names("time","cumhaz","trans") %>% 
    group_by(trans) %>% 
    mutate(haz=cumhaz-lag(cumhaz),
           cycle=row_number()-1) %>% 
    ungroup %>% 
    filter(cycle!=0) %>% 
    select(trans,haz)
  
  return(haz)
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



# 3. Specify parameters that do not need to be simulated ----
cycle_length<-0.25 #in years
cycles<-10/cycle_length #10 year time horizon

rx_eff_sd<-(log(0.833)-log(0.572))/(2*1.96)
utilities_sd<-state_utility$sd

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


# 4. Create a function to run PSA ----
psa_func<-function(sim.seed,
                   rx_cycles_psa, 
                   rx_rr_psa, # treatment effect in risk ratio
                   rx_inst_psa, # whether to treat institutionalized people (yes/no)
                   wane_fc_psa, # treatment waning factor
                   age_init_psa, 
                   sex_psa,
                   apoegeno_psa,
                   start_state_psa,
                   cost_mci_psa,
                   cost_ad_psa,
                   r_cost_psa,
                   r_health_psa){
  
  set.seed(sim.seed)
  
  # generate a random sample for costs
  lopnr_adlist<-data.frame(LOPNR=sample(unique(costwide$LOPNR), replace = TRUE))
  
  sample_adcost<-lopnr_adlist %>% 
    left_join(costwide,
              relationship = "many-to-many")
  
  lopnr_mcilist<-data.frame(LOPNR=sample(unique((mci_svedem_costs %>% 
                                                   filter(state=="MCI"))$LOPNR, 
                                                replace = TRUE))) %>% 
    rbind(data.frame(LOPNR=sample(unique((mci_svedem_costs %>% 
                                            filter(state=="MCI_inst"))$LOPNR, 
                                         replace = TRUE))))
  
  sample_mcicost<-lopnr_mcilist %>% 
    left_join(mci_svedem_costs,
              relationship = "many-to-many")
  
  # generate bootstrapping hazard tables
  boots_haz_mci_ad<-boots_gethaz_mci(mci_ad_model,age_init_psa,sex_psa,apoegeno_psa)
  boots_haz_mci_death<-boots_gethaz_mci(mci_death_model,age_init_psa,sex_psa,apoegeno_psa)
  boots_haz_dem<-boots_gethaz_msm(age_init_psa,sex_psa)
  
  psa_runmodel_cust<-function(rx_cycles, 
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
    utilities_death=0, 
    
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
    
    # transition probability for MCI to institutionalization 
    # from a 6-year study with a population with mean age of 74.5
    mci_inst_haz=-(log(1-0.043))/6,
    mci_inst_age_loghr=log(1.09), # 1.09 (1.05-1.13)
    mci_inst_haz_demo=mci_inst_haz*cycle_length*
      (exp(mci_inst_age_loghr*(age_init-74.5))),
    
    t1=1-exp(-(boots_haz_mci_ad$haz)[model_time]*rX),
    t2=1-exp(-mci_inst_haz_demo),
    t3=1-exp(-(boots_haz_mci_death$haz)[model_time]),
    t4=1-exp(-((boots_haz_dem %>% filter(trans==1))$haz)[model_time]*rX),
    t5=1-exp(-((boots_haz_dem %>% filter(trans==2))$haz)[model_time]),
    t6=1-exp(-((boots_haz_dem %>% filter(trans==3))$haz)[model_time]),
    t7=1-exp(-((boots_haz_dem %>% filter(trans==4))$haz)[model_time]),
    t8=1-exp(-((boots_haz_dem %>% filter(trans==5))$haz)[model_time]),
    t9=1-exp(-((boots_haz_dem %>% filter(trans==6))$haz)[model_time]),
    t10=1-exp(-((boots_haz_dem %>% filter(trans==7))$haz)[model_time]),
    t11=1-exp(-((boots_haz_dem %>% filter(trans==8))$haz)[model_time]),
    t12=1-exp(-(boots_haz_mci_ad$haz)[model_time]*exp(log(rX)*rx_inst)),
    t13=1-exp(-(boots_haz_mci_death$haz)[model_time]),
    t14=1-exp(-((boots_haz_dem %>% filter(trans==9))$haz)[model_time]*exp(log(rX)*rx_inst)),
    t15=1-exp(-((boots_haz_dem %>% filter(trans==10))$haz)[model_time]),
    t16=1-exp(-((boots_haz_dem %>% filter(trans==11))$haz)[model_time]),
    t17=1-exp(-((boots_haz_dem %>% filter(trans==12))$haz)[model_time]),
    t18=1-exp(-((boots_haz_dem %>% filter(trans==13))$haz)[model_time])
  )
  
  # Define states, costs and utilities----
  MCI = define_state(
    utility = (utilities_mci*(1-prob_severe_aria-prob_irr)+
                 (utilities_mci+disutility_aria)*prob_severe_aria+
                 (utilities_mci+disutility_irr)*prob_irr)
    *cycle_length*disc_health, 
    
    cost = ((cost_mci[1]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_mci[1]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_mci[1]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_mci[1]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length, # treatment time/cycle = 1/(no.cycles/year) - full treatment time every cycle
    
    cost_care=cost_mci[1]*cycle_length*disc_cost,
    
    cost_admin=cost_admin*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*disc_cost
  )
  Mild = define_state(
    utility = (utilities_mild*(1-prob_severe_aria-prob_irr)+
                 (utilities_mild+disutility_aria)*prob_severe_aria+
                 (utilities_mild+disutility_irr)*prob_irr)
    *cycle_length*disc_health, 
    
    cost = ((cost_ad[1]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_ad[1]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_ad[1]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_ad[1]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length,
    
    cost_care=cost_ad[1]*cycle_length*disc_cost,
    
    cost_admin=cost_admin*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*disc_cost
  )
  Moderate = define_state(
    utility = utilities_mod*cycle_length*disc_health,
    cost = cost_ad[2]*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_ad[2]*cycle_length*disc_cost,
    cost_admin=0,
    cost_aria=0
  )
  Severe = define_state(
    utility = utilities_sev*cycle_length*disc_health,
    cost = cost_ad[3]*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_ad[3]*cycle_length*disc_cost,
    cost_admin=0,
    cost_aria=0
  )
  MCI_inst = define_state(
    utility = (utilities_mci_inst*(1-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
                 (utilities_mci_inst+disutility_aria)*prob_severe_aria*rx_inst+
                 (utilities_mci_inst+disutility_irr)*prob_irr*rx_inst)
    *cycle_length*disc_health, 
    
    cost = ((cost_mci[2]*cycle_length+cost_admin*rx_inst)*(1-prob_mild_aria*rx_inst-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
              (cost_mci[2]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria*rx_inst+
              (cost_mci[2]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria*rx_inst+
              (cost_mci[2]*cycle_length+cost_admin+cost_irr)*prob_irr*rx_inst)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length*rx_inst,
    
    cost_care=cost_mci[2]*cycle_length*disc_cost,
    
    cost_admin=cost_admin*rx_inst*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*rx_inst*disc_cost
  )
  Mild_inst = define_state(
    utility = (utilities_mild_inst*(1-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
                 (utilities_mild_inst+disutility_aria)*prob_severe_aria*rx_inst+
                 (utilities_mild_inst+disutility_irr)*prob_irr*rx_inst)
    *cycle_length*disc_health, 
    
    cost = ((cost_ad[4]*cycle_length+cost_admin*rx_inst)*(1-prob_mild_aria*rx_inst-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
              (cost_ad[4]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria*rx_inst+
              (cost_ad[4]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria*rx_inst+
              (cost_ad[4]*cycle_length+cost_admin+cost_irr)*prob_irr*rx_inst)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length*rx_inst,
    
    cost_care=cost_ad[4]*cycle_length*disc_cost,
    
    cost_admin=cost_admin*rx_inst*disc_cost,
    
    cost_aria=(cost_mild_aria*prob_mild_aria+cost_severe_aria*prob_severe_aria+
                 cost_irr*prob_irr)*rx_inst*disc_cost
  )
  Moderate_inst = define_state(
    utility = utilities_mod_inst*cycle_length*disc_health,
    cost = cost_ad[5]*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_ad[5]*cycle_length*disc_cost,
    cost_admin=0,
    cost_aria=0
  )
  Severe_inst = define_state(
    utility = utilities_sev_inst*cycle_length*disc_health,
    cost = cost_ad[6]*cycle_length*disc_cost,
    rx.time=0,
    cost_care=cost_ad[6]*cycle_length*disc_cost,
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
  return(res_mod)
  }
  
  temp_model<-psa_runmodel_cust(rx_cycles=rx_cycles_psa,rx_rr=rx_rr_psa,
                                rx_inst=rx_inst_psa,wane_fc=wane_fc_psa,
                                age_init = age_init_psa,sex=sex_psa,
                                apoegeno=apoegeno_psa,start_state=start_state_psa,
                                cost_mci=cost_mci_psa,
                                cost_ad=cost_ad_psa,
                                r_cost=r_cost_psa,
                                r_health=r_health_psa)
  
  tryCatch({ # return error messages instead of stopping iterations
    psa_para_mod<-define_psa(rx_eff_log~normal(log(0.69),rx_eff_sd),
                             utilities_mci~normal(utilities[1], utilities_sd[1]),
                             utilities_mild~normal(utilities[2], utilities_sd[2]),
                             utilities_mod~normal(utilities[3], utilities_sd[3]),
                             utilities_sev~normal(utilities[4], utilities_sd[4]),
                             utilities_mci_inst~normal(utilities[5], utilities_sd[5]),
                             utilities_mild_inst~normal(utilities[6], utilities_sd[6]),
                             utilities_mod_inst~normal(utilities[7], utilities_sd[7]),
                             utilities_sev_inst~normal(utilities[8], utilities_sd[8])) 
    
    
    temp_model_psa<-run_psa(temp_model, psa = psa_para_mod, N = 1, keep = T)
    
    temp_ly_soc<-temp_model_psa$full$standard[[1]]$counts %>% summarize_all(~sum(.))*cycle_length
    temp_ly_rx<-temp_model_psa$full$rx[[1]]$counts %>% summarize_all(~sum(.))*cycle_length
    
    temp_model_psa_table<-data.frame(group=temp_model_psa$ps$.strategy_names,
                                       cost=temp_model_psa$ps$.cost,
                                     cost_care=temp_model_psa$ps$cost_care,
                                     cost_admin=temp_model_psa$ps$cost_admin,
                                     cost_aria=temp_model_psa$ps$cost_aria,
                                       qaly=temp_model_psa$ps$.effect,
                                       rx.time=temp_model_psa$ps$rx.time) %>% 
      pivot_wider(names_from=group, values_from=c(cost,cost_care,cost_admin,
                                                  cost_aria,qaly, rx.time)) %>% 
      select(-rx.time_standard) %>% 
      set_names(c("cost_soc","cost_rx","cost_care_soc","cost_care_rx",
                  "cost_admin_soc","cost_admin_rx",
                  "cost_aria_soc","cost_aria_rx",
                  "effect_soc","effect_rx","rx.time")) %>% 
      cbind(temp_ly_soc %>% 
              set_names(paste0("ly_",c("mci","mild","mod","sev",
                                       "mci_inst","mild_inst","mod_inst",
                                       "sev_inst","death"),"_soc"))) %>% 
      cbind(temp_ly_rx %>% 
              set_names(paste0("ly_",c("mci","mild","mod","sev",
                                       "mci_inst","mild_inst","mod_inst",
                                       "sev_inst","death"),"_rx")))
    
    return(temp_model_psa_table)
  }, error=function(e){
    return(paste("Error:", e$message))
  })
  }



