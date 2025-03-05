# ************************************************************************************************
# Script for senario analyses:
# In CLARITY trial, treatment favors MCI, APOE non-carrier, male, and age >=75
# NEJM paper, supplementary Figure S1 - percent slowing of decline
# ************************************************************************************************


# create a customized "run_model" function to run the HE model
scen_runmodel_cust<-function(rx_cycles, 
                        rx_rr, # treatment effect in risk ratio
                        wane_fc, # treatment waning factor
                        rx_inst, # whether to treat institutionalized people (yes/no)
                        rx_sex_mod, # whether treatment vary by sex (yes/no)
                        rx_age_mod, # whether treatment vary by age groups (yes/no)
                        rx_apoe_mod, # whether treatment vary by apoe (yes/no)
                        rx_mci_mod, # whether treatment vary by starting stage - mci (yes/no)
                        rx_ad_mod, # whether treatment vary by starting stage - mild AD (yes/no)
                        age_init, 
                        sex,
                        apoegeno,
                        start_state,
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
  
  rx_eff_log=log(rx_rr),
  rx_eff=exp(rx_eff_log),
  
  # consider effect modification by subgroups
  cd_sex_rr=ifelse(sex==0,1.59,0.44),
  rX_sex_rr=case_when(rx_sex_mod==0~1,
                      rx_sex_mod==1~(1+(rx_eff-1)*cd_sex_rr)/rx_eff),
  
  cd_age_rr=case_when(age_init<65~0.22,
                      age_init<75~0.85,
                      age_init>=75~1.48),
  rX_age_rr=case_when(rx_age_mod==0~1,
                      rx_age_mod==1~(1+(rx_eff-1)*cd_age_rr)/rx_eff),
  
  cd_apoe_rr=ifelse(apoegeno==0,1.52,1.11),
  rX_apoe_rr=case_when(rx_apoe_mod==0~1,
                       rx_apoe_mod==1~(1+(rx_eff-1)*cd_apoe_rr)/rx_eff),
  
  cd_mci_rr=ifelse(start_state==0,1.04,1),
  rX_mci_rr=case_when(rx_mci_mod==0~1,
                      rx_mci_mod==1~(1+(rx_eff-1)*cd_mci_rr)/rx_eff),
  
  cd_ad_rr=ifelse(start_state==1,1,1),
  rX_ad_rr=case_when(rx_ad_mod==0~1,
                     rx_ad_mod==1~(1+(rx_eff-1)*cd_ad_rr)/rx_eff),
  
  # final effect for MCI and mild AD
  rx_eff_mci=rx_eff*rX_sex_rr*rX_age_rr*rX_apoe_rr*rX_mci_rr,
  rx_eff_ad=rx_eff*rX_sex_rr*rX_age_rr*rX_apoe_rr*rX_ad_rr,

  
  # dispatch treatment effect - consider waning effect after treatment ends
  rX_mci=dispatch_strategy(
    standard=1,
    rx=case_when(model_time<=rx_cycles ~ rx_eff_mci, 
                 model_time>rx_cycles ~ rx_eff_mci^((1-wane_fc)^(model_time-rx_cycles)))
  ),
  
  rX_ad=dispatch_strategy(
    standard=1,
    rx=case_when(model_time<=rx_cycles ~ rx_eff_ad, 
                 model_time>rx_cycles ~ rx_eff_ad^((1-wane_fc)^(model_time-rx_cycles)))
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
  
  
  # Probability of ARIA-E
  prob_mild_ariae=ifelse(apoegeno==0,0.014,0.013),
  prob_sev_ariae=ifelse(apoegeno==0,0.002,0.008),
  
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
  
  # Probability of isolated ARIA-H
  prob_mild_ariah=ifelse(apoegeno==0,0.003,0.004),
  prob_sev_ariah=ifelse(apoegeno==0,0.004,0.003),
  
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
  
  
  t1=1-exp(-(get_haz(1,age_init,sex,apoegeno))[model_time,1]*rX_mci),
  t2=1-exp(-mci_inst_haz_demo),
  t3=1-exp(-(get_haz(3,age_init,sex))[model_time,1]),
  t4=1-exp(-(get_haz(4,age_init,sex))[model_time,1]*rX_ad),
  t5=1-exp(-(get_haz(5,age_init,sex))[model_time,1]),
  t6=1-exp(-(get_haz(6,age_init,sex))[model_time,1]),
  t7=1-exp(-(get_haz(7,age_init,sex))[model_time,1]),
  t8=1-exp(-(get_haz(8,age_init,sex))[model_time,1]),
  t9=1-exp(-(get_haz(9,age_init,sex))[model_time,1]),
  t10=1-exp(-(get_haz(10,age_init,sex))[model_time,1]),
  t11=1-exp(-(get_haz(11,age_init,sex))[model_time,1]),
  t12=1-exp(-(get_haz(12,age_init,sex,apoegeno))[model_time,1]*exp(log(rX_mci)*rx_inst)),
  t13=1-exp(-(get_haz(13,age_init,sex))[model_time,1]),
  t14=1-exp(-(get_haz(14,age_init,sex))[model_time,1]*exp(log(rX_ad)*rx_inst)),
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
  rx.time.disc=rx.time*disc
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
  rx.time.disc=rx.time*disc
)
Moderate = define_state(
  utility = utilities_mod*cycle_length*disc,
  cost = cost_ad[2]*cycle_length*disc,
  rx.time=0,
  rx.time.disc=0
)
Severe = define_state(
  utility = utilities_sev*cycle_length*disc,
  cost = cost_ad[3]*cycle_length*disc,
  rx.time=0,
  rx.time.disc=0
)
MCI_inst = define_state(
  utility = (utilities_mci_inst*(1-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
               (utilities_mci_inst+disutility_aria)*prob_severe_aria*rx_inst+
               (utilities_mci_inst+disutility_irr)*prob_irr*rx_inst)
  *cycle_length*disc, 
  
  cost = ((cost_mci[2]*cycle_length+cost_admin*rx_inst)*(1-prob_mild_aria*rx_inst-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
            (cost_mci[2]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria*rx_inst+
            (cost_mci[2]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria*rx_inst+
            (cost_mci[2]*cycle_length+cost_admin+cost_irr)*prob_irr*rx_inst)*disc, # cost/cycle = yearly cost/(no.cycles/year)
  
  rx.time=rxtime*cycle_length*rx_inst,
  rx.time.disc=rx.time*disc
)
Mild_inst = define_state(
  utility = (utilities_mild_inst*(1-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
               (utilities_mild_inst+disutility_aria)*prob_severe_aria*rx_inst+
               (utilities_mild_inst+disutility_irr)*prob_irr*rx_inst)
  *cycle_length*disc, 
  
  cost = ((cost_ad[4]*cycle_length+cost_admin*rx_inst)*(1-prob_mild_aria*rx_inst-prob_severe_aria*rx_inst-prob_irr*rx_inst)+
            (cost_ad[4]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria*rx_inst+
            (cost_ad[4]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria*rx_inst+
            (cost_ad[4]*cycle_length+cost_admin+cost_irr)*prob_irr*rx_inst)*disc, # cost/cycle = yearly cost/(no.cycles/year)
  
  rx.time=rxtime*cycle_length*rx_inst,
  rx.time.disc=rx.time*disc
)
Moderate_inst = define_state(
  utility = utilities_mod_inst*disc*cycle_length,
  cost = cost_ad[5]*cycle_length*disc,
  rx.time=0,
  rx.time.disc=0
)
Severe_inst = define_state(
  utility = utilities_sev_inst*cycle_length*disc,
  cost = cost_ad[6]*cycle_length*disc,
  rx.time=0,
  rx.time.disc=0
)
Death = define_state(
  utility = 0,
  cost = 0,
  rx.time=0,
  rx.time.disc=0
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


