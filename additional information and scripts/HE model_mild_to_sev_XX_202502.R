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
load("models_xx_mild_to_sev.RData")



# specify models for each transition,except for transitions from MCI
model_trans1<-mci_ad_model # MCI to mild AD
model_trans3<-mci_death_model # MCI to death
model_trans4<-final_model_trans1 # mild AD to moderate AD
model_trans5<-final_model_trans2 # mild AD to severe AD
model_trans6<-final_model_trans3 # mild AD to institutionalization
model_trans7<-final_model_trans4 # mild AD to death
model_trans8<-final_model_trans5 # moderate AD to severe AD
model_trans9<-final_model_trans6 # moderate AD to institutionalization
model_trans10<-final_model_trans7 # moderate AD to death
model_trans11<-final_model_trans8 # severe AD to institutionalization
model_trans12<-final_model_trans9 # severe AD to death
model_trans13<-mci_ad_model # institutionalized MCI to institutionalized mild AD
model_trans14<-mci_death_model # institutionalized MCI to death
model_trans15<-final_model_trans10 # institutionalized mild AD to institutionalized moderate AD
model_trans16<-final_model_trans11 # institutionalized mild AD to institutionalized severe AD
model_trans17<-final_model_trans12 # institutionalized mild AD to death
model_trans18<-final_model_trans13 # institutionalized moderate AD to institutionalized severe AD
model_trans19<-final_model_trans14 # institutionalized moderate AD to death
model_trans20<-final_model_trans15 # institutionalized severe AD to death


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
runmodel_cust<-function(rx_cycles, 
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
    t5=1-exp(-(get_haz(5,age_init,sex))[model_time,1]*rX),
    t6=1-exp(-(get_haz(6,age_init,sex))[model_time,1]),
    t7=1-exp(-(get_haz(7,age_init,sex))[model_time,1]),
    t8=1-exp(-(get_haz(8,age_init,sex))[model_time,1]),
    t9=1-exp(-(get_haz(9,age_init,sex))[model_time,1]),
    t10=1-exp(-(get_haz(10,age_init,sex))[model_time,1]),
    t11=1-exp(-(get_haz(11,age_init,sex))[model_time,1]),
    t12=1-exp(-(get_haz(12,age_init,sex))[model_time,1]),
    t13=1-exp(-(get_haz(13,age_init,sex,apoegeno))[model_time,1]*rX),
    t14=1-exp(-(get_haz(14,age_init,sex))[model_time,1]),
    t15=1-exp(-(get_haz(15,age_init,sex))[model_time,1]*rX),
    t16=1-exp(-(get_haz(16,age_init,sex))[model_time,1]*rX),
    t17=1-exp(-(get_haz(17,age_init,sex))[model_time,1]),
    t18=1-exp(-(get_haz(18,age_init,sex))[model_time,1]),
    t19=1-exp(-(get_haz(19,age_init,sex))[model_time,1]),
    t20=1-exp(-(get_haz(20,age_init,sex))[model_time,1])
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
    utility = (utilities_mci_inst*(1-prob_severe_aria-prob_irr)+
                 (utilities_mci_inst+disutility_aria)*prob_severe_aria+
                 (utilities_mci_inst+disutility_irr)*prob_irr)
    *cycle_length*disc, 
    
    cost = ((cost_mci[2]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_mci[2]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_mci[2]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_mci[2]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc, # cost/cycle = yearly cost/(no.cycles/year)

    rx.time=rxtime*cycle_length,
    rx.time.disc=rx.time*disc
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
    0, C,     t4, t5,     0,  t6,     0,      0,      t7, 
    0, 0,     C,     t8,    0,  0,      t9,     0,      t10, 
    0, 0,     0,     C,     0,  0,      0,      t11,    t12,
    0, 0,     0,     0,     C,  t13, 0,      0,      t14,
    0, 0,     0,     0,     0,  C,      t15, t16,      t17,
    0, 0,     0,     0,     0,  0,      C,      t18,    t19, 
    0, 0,     0,     0,     0,  0,      0,      C,      t20,
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
    effect = utility,
    method='beginning'
  ) 
  return(res_mod)
}


# 4. Base-case scenario-----
# sex: 0 = male, 1 = female
# apoegeno: 0 = Non-carrier, 1 = heterozygote
# start_state: 0 = MCI, 1 = Mild AD

# call starting population
load("pop_perc.RData")

# 4.1. Calculate cost-effectiveness price at a WTP of 1000000 ----
# setting up parallel computing
ncpus = parallel::detectCores()-1
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod)))
clusterExport(cl, c("pop_perc","runmodel_cust","get_haz",
                    "statenames","cycles","cycle_length","utilities",
                    "getcost","mci_svedem_costs","costwide",
                    paste0("model_trans",c(1,3:18))))

base_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  setting_i<-pop_perc[i,]$inst
  perc<-pop_perc[i,]$pop_perc
  
  temp_model<-runmodel_cust(rx_cycles=rx_cycles_i,age_init = age,sex=sex_i,
                            apoegeno=apoe_i,start_state=stage_i,
                            setting=setting_i,
                            wane_fc=1,
                            cost_mci=getcost(mci_svedem_costs),
                            cost_ad=getcost(costwide),
                            r=0.03)
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age=age,sex=sex_i,apoe=apoe_i,stage=stage_i,
              setting=setting_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})




base_mod_summary <- bind_rows(lapply(base_mod, function(x) {
  data.frame(
    age = x$age, sex = x$sex, apoe = x$apoe, stage = x$stage, 
    setting = x$setting, perc = x$perc, qaly_stand = x$qaly_stand,
    qaly_rx=x$qaly_rx,costs_stand = x$costs_stand,costs_rx=x$costs_rx,
    deffect = x$deffect, dcost = x$dcost, rxtime = x$rxtime
  )
}))


base_mod_summary<-base_mod_summary %>% 
  mutate(deffect_w=deffect*perc,
         dcost_w=dcost*perc,
         qaly_rx_w=qaly_rx*perc,
         qaly_stand_w=qaly_stand*perc,
         costs_rx_w=costs_rx*perc,
         costs_stand_w=costs_stand*perc,
         rxtime_w=rxtime*perc)

QALYval=1000000
qaly=sum(base_mod_summary$deffect_w)
cost=sum(base_mod_summary$dcost_w)
NMB=qaly*QALYval-cost
bc_price=NMB/sum(base_mod_summary$rxtime_w)
bc_price

sum(base_mod_summary$qaly_rx_w)
sum(base_mod_summary$qaly_stand_w)
sum(base_mod_summary$costs_rx_w)
sum(base_mod_summary$costs_stand_w)


# define factor labels and legend colors for later plotting
state_label<-c("MCI due to AD","Mild AD","Moderate AD","Severe AD",
               "Institutionalized MCI due to AD",
               "Institutionalized mild AD","Institutionalized moderate AD",
               "Institutionalized severe AD","Death")

cbp <- c("#999999","#CC79A7","#E69F00","#F0E442",
         "#D55E00","pink","#56B4E9","#0072B2",
         "#009E73")

# 4.2. Generate plots for sojourn time and transition probability comparisons ----
clusterExport(cl, c("base_mod","pop_perc","runmodel_cust","get_haz",
                    "statenames","cycles","cycle_length","utilities",
                    "getcost","mci_svedem_costs","costwide",
                    paste0("model_trans",c(1,3:18))))

# Plot - sojourn times by state
sojourn<-parLapply(cl,1:length(base_mod), function(i){
  
  temp_mod<-base_mod[[i]]$model[[1]]
  temp_model_info<-data.frame(model=i,perc=base_mod[[i]]$perc)
  temp_sojourn<-rbind(temp_mod$eval_strategy_list$standard$counts %>% 
                        summarize_all(~sum(.)) %>% 
                        mutate(strategy='standard'),
                      temp_mod$eval_strategy_list$rx$counts %>% 
                        summarize_all(~sum(.)) %>% 
                        mutate(strategy='rx')) %>% 
    pivot_longer(cols=-strategy, names_to='state', values_to = 'years') %>% 
    mutate(years=years*0.25,
           state=factor(state, levels=rev(statenames),labels = rev(state_label)),
           strategy=factor(strategy,levels=c("standard","rx"),
                           labels=c("Standard","Treatment")))
  
  return(bind_cols(temp_model_info, temp_sojourn))
}) %>% 
  bind_rows()


stopCluster(cl)

ggplot(data = sojourn %>% 
         mutate(years=years*perc) %>%  
         group_by(strategy,state) %>% 
         summarise(years=sum(years)) %>% 
         ungroup(), 
       aes(y = years, x = strategy, fill = state)) +
  geom_col(position = 'stack') +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Group", y = "Time horizon (years)", fill = "State") +
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
  mutate(years=years*perc) %>%  
  group_by(strategy,state) %>% 
  summarise(years=sum(years)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=state, names_from='strategy', values_from='years') %>%
  mutate(years=(Treatment-Standard), 
         strategy='difference') %>% select(state, years, strategy) %>% 
  mutate(month=years*12) 

ggplot(data=sojourn_diff, aes(y=month, x=state, fill=state)) + 
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



# Plot - transition probability by treatment groups
cohdist<-lapply(1:length(base_mod), function(i){
  
  temp_mod<-base_mod[[i]]$model[[1]]
  temp_perc<-base_mod[[i]]$perc
  
  standard_counts<-temp_mod$eval_strategy_list$standard$counts*temp_perc
  rx_counts<-temp_mod$eval_strategy_list$rx$counts*temp_perc
  
  temp_cohdist <- bind_rows(
    standard_counts %>% mutate(group = "SoC"),
    rx_counts %>% mutate(group = "Treatment")
  )
  
  temp_cohdist$t=rep(1:cycles,2)
  
  return(temp_cohdist)
}) %>% 
  bind_rows() %>% 
  group_by(group,t) %>%
  summarize(across(everything(), ~ sum(.x))) %>% 
  ungroup
  
  
  


coh.inst<-pivot_longer(cohdist,cols=MCI:Death,names_to='State') %>% 
  mutate(years=t*0.25,
         State=factor(State, levels=rev(statenames),labels = rev(state_label)))

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




# plot - total costs by treatment groups
costplot<-lapply(1:length(base_mod), function(i){
  
  temp_mod<-base_mod[[i]]$model[[1]]
  temp_perc<-base_mod[[i]]$perc
  
  standard_costs<-data.frame(Markov.cycle=temp_mod$eval_strategy_list$standard$values$model_time, 
                             cost=(temp_mod$eval_strategy_list$standard$values$cost)*temp_perc)
  
  rx_costs<-data.frame(Markov.cycle=temp_mod$eval_strategy_list$standard$values$model_time, 
                       cost=(temp_mod$eval_strategy_list$rx$values$cost)*temp_perc)
  
  temp_costs <- bind_rows(
    standard_costs %>% mutate(group = "SoC"),
    rx_costs %>% mutate(group = "Treatment")
  )
  return(temp_costs)
}) %>% 
  bind_rows() %>% 
  group_by(group,Markov.cycle) %>%
  summarize(across(everything(), ~ sum(.x))) %>% 
  ungroup %>% 
  mutate(years=Markov.cycle*0.25)


ggplot(aes(y=cost, x=years, color=group), data=costplot) + 
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks = seq(0, 10, by = 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(floor(min(costplot$cost)/10000)*10000, 
                                  ceiling(max(costplot$cost)/10000)*10000, by = 10000),
                     limits = c(floor(min(costplot$cost)/10000)*10000,
                                ceiling(max(costplot$cost)/10000)*10000)) +
  labs(x = "Time horizon (years)", y = "Costs (2023 SEK)", color = "Strategy",
       linetype="Group") +
  scale_color_manual(values = c("#E69F00","#56B4E9")) +
  theme_minimal(base_size = 15) + 
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")  
  )



# 5. One-way sensitivity analysis ----
# ************************************************************************************************
# Considerations (from NICE literature review):
# 1. Varying treatment effects (varying HR)
# 2. Uncertainty around Utility (varying utility by state)
# 3. Lecanemab - cost of administration
# # WTP reference from TLV (https://www.valueinhealthjournal.com/article/S1098-3015(22)03891-8/fulltext)
# ************************************************************************************************

# 5.1. Define DSA parameters ----
dsa_para_mod <- define_dsa(
  rx_eff_log, log(0.572), log(0.833),
  mci_ad_haz,0.0969,0.1848,
  mci_ad_age_loghr,log(1.03),log(1.1),
  mci_ad_female_loghr,log(1.08),log(1.64),
  
  utilities_mci, utilities_lb[1], utilities_ub[1],
  utilities_mild, utilities_lb[2], utilities_ub[2],
  utilities_mod, utilities_lb[3], utilities_ub[3],
  utilities_sev, utilities_lb[4], utilities_ub[4],
  utilities_mci_inst, utilities_lb[5], utilities_ub[5],
  utilities_mild_inst, utilities_lb[6], utilities_ub[6],
  utilities_mod_inst, utilities_lb[7], utilities_ub[7],
  utilities_sev_inst, utilities_lb[8], utilities_ub[8]
)

res_dsa<-run_dsa(
  model=base_mod,
  dsa=dsa_para_mod
)

res_dsa_summary<-res_dsa[["dsa"]] %>% 
  mutate(.par_value=ifelse(.par_names %in% c("rx_eff_log","mci_ad_haz","mci_ad_age_loghr","mci_ad_female_loghr"),
                           exp(.par_value_eval),res_dsa[["dsa"]]$.par_value_eval),
         .par_value=format(round(.par_value,2),nsmall=2),
         names=case_match(.par_names, "rx_eff_log"~"Treatment effect",
                          "mci_ad_haz"~"MCI hazard",
                          "mci_ad_age_loghr"~"MCI age hazard ratio",
                          "mci_ad_female_loghr"~"MCI sex hazard ratio",
                               "utilities_mci"~"Utility of MCI",
                               "utilities_mild"~"Utility of mild AD",
                               "utilities_mod"~"Utility of moderate AD",
                               "utilities_sev"~"Utility of severe AD",
                               "utilities_mci_inst"~"Utility of institutionalized MCI",
                               "utilities_mild_inst"~"Utility of institutionalized mild AD",
                               "utilities_mod_inst"~"Utility of institutionalized moderate AD",
                               "utilities_sev_inst"~"Utility of institutionalized severe AD"))

res_dsa_summary<-res_dsa_summary %>% 
  filter(.strategy_names=="standard") %>% 
  dplyr::select(.cost,.effect,names,.par_value) %>% 
  rename(cost_soc=.cost,
         effect_soc=.effect) %>% 
  left_join(res_dsa_summary %>% 
              filter(.strategy_names=="rx") %>% 
              dplyr::select(.cost,.effect,names,.par_value,rx.time) %>% 
              rename(cost_rx=.cost,
                     effect_rx=.effect)) %>% 
  mutate(dcost=cost_rx-cost_soc,
         effect=effect_rx-effect_soc,
         QALYval=1000000,
         NMB=effect*QALYval-dcost,
         rxtime=rx.time,
         price_diff=NMB/rxtime-bc_price,
         level=as.factor(ifelse(price_diff<0,0,1)))
  
factor_level<-(res_dsa_summary %>% 
    group_by(names) %>%
      slice_max(abs(price_diff)) %>% 
      ungroup %>% 
      distinct(names,.keep_all = T) %>%
    arrange(abs(price_diff)))$names

res_dsa_summary$names<-factor(res_dsa_summary$names,levels=factor_level)


# 5.2. Tornado plot ----  
ggplot(res_dsa_summary, aes(factor(names), price_diff, fill = level)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(y = price_diff, label = .par_value), # Place .par_value at the end of each bar
            hjust = ifelse(res_dsa_summary$price_diff > 0, -0.2, 1.2), # Adjust based on bar direction
            size = 3, color = "black") +
  coord_flip() +
  scale_y_continuous(breaks = seq(-70000, 60000, by = 10000),
                     labels = scales::label_number(accuracy = 1, big.mark = ",")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  labs(y = "Differences in cost-effective prices",
       x = "Parameter") +
  theme_minimal() +
  theme(legend.position = "none")




# 6. Bootstrapping Probabilistic sensitivity analysis -----
# 6.1. Define PSA parameters not estimated from SveDem data ----
rx_eff_sd<-(log(0.833)-log(0.572))/(2*1.96)
utilities_sd<-state_utility$sd

mci_ad_haz_sd<-(0.174-0.164)/1.96
mci_ad_age_loghr_sd<-(log(1.06)-log(1.03))/1.96
mci_ad_sex_loghr_sd<-(log(1.33)-log(1.08))/1.96

psa_para_mod<-define_psa(rx_eff_log~normal(-0.3710637,rx_eff_sd),
                         mci_ad_haz~normal(0.174,mci_ad_haz_sd),
                         mci_ad_age_loghr~normal(log(1.06),mci_ad_age_loghr_sd),
                         mci_ad_female_loghr~normal(log(1.33),mci_ad_female_loghr_sd),
                         
                         utilities_mci~normal(utilities[1], utilities_sd[1]),
                         utilities_mild~normal(utilities[2], utilities_sd[2]),
                         utilities_mod~normal(utilities[3], utilities_sd[3]),
                         utilities_sev~normal(utilities[4], utilities_sd[4]),
                         utilities_mci_inst~normal(utilities[5], utilities_sd[5]),
                         utilities_mild_inst~normal(utilities[6], utilities_sd[6]),
                         utilities_mod_inst~normal(utilities[7], utilities_sd[7]),
                         utilities_sev_inst~normal(utilities[8], utilities_sd[8])) 

res_psa<-run_psa(base_mod, psa = psa_para_mod, N = 100)
res_psa_table<-data.frame(group=res_psa$ps$.strategy_names,
                          iteration=res_psa$ps$.index,
                          cost=res_psa$ps$.cost,
                          qaly=res_psa$ps$.effect,
                          rx.time.disc=res_psa$ps$rx.time.disc)

res_psa_table_wide<-res_psa_table %>% 
  pivot_wider(id_cols=iteration, names_from=group, values_from=c(cost, qaly, rx.time.disc)) %>%
  mutate(incr_cost=cost_rx-cost_standard, 
         incr_qaly=qaly_rx-qaly_standard) 



# 6.2. Bootstrapping for transition probabilities and costs (estimated from SveDem data) ----
source("C:/Users/xinxia/OneDrive - Karolinska Institutet/CSF-registersamkörning/AD DMT HE model/Statistical analyses/New analyses_XX_202412/HE model_PSA_function_XX_20250107.R")


# 6.3. Summarize results and draw graphs ----
set.seed(2025)
n_sim<-200
# 1000 bootstrapping
psa<-pblapply(1:1000, function(i) {
  print(i)
  boots_func(age_init = 70, sex = "MALE", n_sim=n_sim) %>% 
    # n_sim: the number of simulations for psa based on parameter distributions
    mutate(iteration = iteration+n_sim*(i-1))}) %>%
  bind_rows()




# 6.3.1. Threshold cost-effective price ----
WTP=seq(100000, 3000000, 100000)

WTP_results<-cross_join(res_psa_table_wide, data.frame(WTP=WTP)) %>%
  mutate(NMB=incr_qaly*WTP-incr_cost,
         threshold_price=NMB/rx.time.disc_rx)

CEAC<-WTP_results %>% group_by(WTP) %>%  
  reframe(tibble::enframe(quantile(threshold_price,c(0.05, 0.5, 0.95)),'quantile', 'threshold_price' )) %>%
  pivot_wider(id_cols=WTP, names_from=quantile, values_from=threshold_price)


ggplot()+
  geom_ribbon(aes(ymin=`5%`, ymax=`95%`, x=WTP), data=CEAC, alpha=0.2, fill='red'#, outline.type = 'full'
  )+
  geom_line(aes(y=`50%`, x=WTP), data=CEAC)+
  scale_x_continuous(labels = scales::label_comma())+
  scale_y_continuous(labels = scales::label_comma())+
  ylab('Threshold cost-effective price')+
  xlab('Willingness to pay per QALY')



res_psa_table_wide<-res_psa_table_wide %>% 
  mutate(incr_cost_drug=incr_cost+rx.time.disc_rx*threshold.rxcost, # annual drug cost calculated from previous analysis
         cost_effective=if_else(incr_cost_drug/incr_qaly < 1000000, 
                                'Cost-effective', 'Not cost-effective') )

ellipse_95<-confidence_ellipse(res_psa_table_wide, 
                               x = incr_qaly, 
                               y = incr_cost_drug, 
                               conf_level = 0.95)



res_psa_table_wide %>% 
  ggplot(aes(y=incr_cost_drug, x=incr_qaly))+
  geom_point(aes(color=cost_effective), size=0.5
  )+
  xlim(pmin(0, min(res_psa_table_wide$incr_qaly)), 
       pmax(0,max(res_psa_table_wide$incr_qaly), 
            max(ellipse_95$x)))+
  ylim(pmin(0, min(res_psa_table_wide$incr_cost_drug)), 
       pmax(0,max(res_psa_table_wide$incr_cost_drug+res_psa_table_wide$rx.time.disc_rx), 
            max(ellipse_95$y)))+
  geom_path(data = ellipse_95, aes(x = x, y = y), color = "black", linewidth = 0.5) +
  geom_abline(intercept = 0, slope = 1000000, linetype = "dashed", color = "black") +
  coord_equal(ratio = 1/1000000) +
  annotate(
    "text",
    x = 0.01,
    y = 0.015 * 1000000,
    angle = atan(1000000 * 1/1000000) * 180/pi,
    label = "WTP=1000000"
  ) +
  labs(x="Incremental effectiveness",y="Incremental cost",color="Cost-effective")+
  theme_minimal() +
  scale_color_manual(values = c("#56B4E9","#E69F00"))









# 7. Scenario analysis ----
# ************************************************************************************************
# Considerations (from NICE literature review):
# 1. Varying discount rate for costs and QALYs
# 2. Varying baseline age
# 3. Lecanemab treatment effects
# 4. Baseline disease severity (proportion of MCI and mild AD)
# ************************************************************************************************





