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


# load data
setwd("C:/Users/xinxia/OneDrive - Karolinska Institutet/CSF-registersamkörning/AD DMT HE model/Statistical analyses/New analyses_XX_202412/")

load("costdata_xx.RData")
load("models_xx.RData")

# call starting population
load("pop_perc.RData")


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
    
    cost = ((cost_mci[1]*cycle_length+cost_admin)*(1-prob_mild_aria-prob_severe_aria-prob_irr)+
              (cost_mci[1]*cycle_length+cost_admin+cost_mild_aria)*prob_mild_aria+
              (cost_mci[1]*cycle_length+cost_admin+cost_severe_aria)*prob_severe_aria+
              (cost_mci[1]*cycle_length+cost_admin+cost_irr)*prob_irr)*disc_cost, # cost/cycle = yearly cost/(no.cycles/year)
    
    rx.time=rxtime*cycle_length # treatment time/cycle = 1/(no.cycles/year) - full treatment time every cycle
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
    
    rx.time=rxtime*cycle_length
    )
  Moderate = define_state(
    utility = utilities_mod*cycle_length*disc_health,
    cost = cost_ad[2]*cycle_length*disc_cost,
    rx.time=0
    )
  Severe = define_state(
    utility = utilities_sev*cycle_length*disc_health,
    cost = cost_ad[3]*cycle_length*disc_cost,
    rx.time=0
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
    
    rx.time=rxtime*cycle_length*rx_inst
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
    
    rx.time=rxtime*cycle_length*rx_inst
    )
  Moderate_inst = define_state(
    utility = utilities_mod_inst*cycle_length*disc_health,
    cost = cost_ad[5]*cycle_length*disc_cost,
    rx.time=0
    )
  Severe_inst = define_state(
    utility = utilities_sev_inst*cycle_length*disc_health,
    cost = cost_ad[6]*cycle_length*disc_cost,
    rx.time=0
    )
  Death = define_state(
    utility = 0,
    cost = 0,
    rx.time=0
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


# 4. Base-case model (waning factor = 1, i.e., no treatment effect after treatment stops)-----
# sex: 0 = male, 1 = female
# apoegeno: 0 = Non-carrier, 1 = heterozygote
# start_state: 0 = MCI, 1 = Mild AD


# 4.1. Calculate cost-effectiveness price at a WTP of 1000000 ----
QALYval=1000000

# setting up parallel computing
ncpus = parallel::detectCores()-1
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod)))
clusterExport(cl, c(ls(all.names = TRUE)))

base_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})

base_mod_summary <- bind_rows(lapply(base_mod, function(x) {
  data.frame(
    rx_cycles=x$rx_cycles,age_group = x$age_group, age = x$age, 
    sex = x$sex, apoe = x$apoe, stage = x$stage, 
    perc = x$perc, qaly_stand = x$qaly_stand,
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
         rxtime_w=rxtime*perc,
         NMB=deffect*QALYval-dcost,
         price=NMB/rxtime)


qaly=sum(base_mod_summary$deffect_w)
cost=sum(base_mod_summary$dcost_w)
NMB=qaly*QALYval-cost
bc_price=NMB/sum(base_mod_summary$rxtime_w)
bc_price

sum(base_mod_summary$qaly_rx_w)
sum(base_mod_summary$qaly_stand_w)
sum(base_mod_summary$costs_rx_w)
sum(base_mod_summary$costs_stand_w)

base_mod_excel<-base_mod_summary %>% 
  select(-c(age,perc,deffect_w:price)) %>% 
  select(rx_cycles,rxtime,age_group:dcost) %>% 
  mutate(sex=factor(sex,levels=0:1,labels=c("Male","Female")),
         apoe=factor(apoe,levels=0:1,labels=c("Non-carrier","Heterozygotes")),
         stage=factor(stage,levels=0:1,labels=c("MCI","Mild AD"))) %>% 
  set_names(c("Treatment cycle","Treatment duration (years)","Age","Sex","APOE genotype",
              "AD stage","QALY-SOC","QALY-treatment","Costs-SOC","Costs-treatment",
              "QALY gain","Incremental costs"))



# define factor labels and legend colors for later plotting
state_label<-c("MCI due to AD","Mild AD","Moderate AD","Severe AD",
         "Institutionalized MCI due to AD",
         "Institutionalized mild AD","Institutionalized moderate AD",
         "Institutionalized severe AD","Death")

cbp <- c("#999999","#CC79A7","#E69F00","#F0E442",
         "#D55E00","pink","#56B4E9","#0072B2",
         "#009E73")

# 4.2. Generate plots for sojourn time and transition probability comparisons ----
# Plot - sojourn times by state
sojourn<-lapply(1:length(base_mod), function(i){
  
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
  
# export state trace to excel

  


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
# 1. Varying treatment effects (varying HR)
# 2. Uncertainty around Utility (varying utility by state)
# 3. Varying discount rates for costs and QALYs: 0 and 0.05, same for costs and QALYs
# ************************************************************************************************
# 5.1. Define DSA parameters through define_dsa() ----
clusterExport(cl, c(ls(all.names = TRUE)))

dsa_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
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
  
  dsa_para_mod <- define_dsa(
    rx_eff, 0.572, 0.833,
    
    utilities_mci, utilities_lb[1], utilities_ub[1],
    utilities_mild, utilities_lb[2], utilities_ub[2],
    utilities_mod, utilities_lb[3], utilities_ub[3],
    utilities_sev, utilities_lb[4], utilities_ub[4],
    utilities_mci_inst, utilities_lb[5], utilities_ub[5],
    utilities_mild_inst, utilities_lb[6], utilities_ub[6],
    utilities_mod_inst, utilities_lb[7], utilities_ub[7],
    utilities_sev_inst, utilities_lb[8], utilities_ub[8]
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
                            "utilities_sev_inst"~"Utility of institutionalized severe AD"))
  
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
         price_diff=price-bc_price)


# 5.2. Varying discounting rate ----
# discount rate = 0
disc0_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})

disc0_mod_summary <- bind_rows(lapply(disc0_mod, function(x) {
  data.frame(
    rx_cycles=x$rx_cycles,age_group = x$age_group, age = x$age, 
    sex = x$sex, apoe = x$apoe, stage = x$stage, 
    perc = x$perc, qaly_stand = x$qaly_stand,
    qaly_rx=x$qaly_rx,costs_stand = x$costs_stand,costs_rx=x$costs_rx,
    deffect = x$deffect, dcost = x$dcost, rxtime = x$rxtime
  )
}))

disc0_mod_summary<-disc0_mod_summary %>% 
  summarise(effect_rx=sum(qaly_rx*perc),
            effect_soc=sum(qaly_stand*perc),
            cost_rx=sum(costs_rx*perc),
            cost_soc=sum(costs_stand*perc),
            rx.time=sum(rxtime*perc)) %>% 
  mutate(NMB=QALYval*(effect_rx-effect_soc)-(cost_rx-cost_soc),
         price=NMB/rx.time,
         price_diff=price-bc_price)


# discount rate = 0.05
disc0.05_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})

disc0.05_mod_summary <- bind_rows(lapply(disc0.05_mod, function(x) {
  data.frame(
    rx_cycles=x$rx_cycles,age_group = x$age_group, age = x$age, 
    sex = x$sex, apoe = x$apoe, stage = x$stage, 
    perc = x$perc, qaly_stand = x$qaly_stand,
    qaly_rx=x$qaly_rx,costs_stand = x$costs_stand,costs_rx=x$costs_rx,
    deffect = x$deffect, dcost = x$dcost, rxtime = x$rxtime
  )
}))

disc0.05_mod_summary<-disc0.05_mod_summary %>% 
  summarise(effect_rx=sum(qaly_rx*perc),
            effect_soc=sum(qaly_stand*perc),
            cost_rx=sum(costs_rx*perc),
            cost_soc=sum(costs_stand*perc),
            rx.time=sum(rxtime*perc)) %>% 
  mutate(NMB=QALYval*(effect_rx-effect_soc)-(cost_rx-cost_soc),
         price=NMB/rx.time,
         price_diff=price-bc_price)


dsa_mod_summary<-dsa_mod_summary %>% 
  rbind(disc0_mod_summary %>% 
          mutate(names="Discount rate for costs and health effects",
                 .par_value=0)) %>% 
  rbind(disc0.05_mod_summary %>% 
          mutate(names="Discount rate for costs and health effects",
                 .par_value=0.05)) %>% 
  group_by(names) %>% 
  mutate(level=as.factor(price_diff==max(price_diff))) %>% 
  mutate(max_absdiff=max(abs(price_diff))) %>% 
  ungroup %>% 
  arrange(desc(max_absdiff))


factor_level<-(dsa_mod_summary %>% 
                 group_by(names) %>%
                 slice_max(abs(price_diff)) %>% 
                 ungroup %>% 
                 distinct(names,.keep_all = T) %>%
                 arrange(abs(price_diff)))$names

dsa_mod_summary$names<-factor(dsa_mod_summary$names,levels=factor_level)


# 5.3. Tornado plot ----  
ggplot(dsa_mod_summary, aes(factor(names), price_diff, fill = level)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(y = price_diff, label = .par_value), # Place .par_value at the end of each bar
            hjust = ifelse(dsa_mod_summary$price_diff > 0, -0.2, 1.2), # Adjust based on bar direction
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
# ************************************************************************************************
# 1. Uncertainty around input parameters considered in the one-way sensitivity analyses
# 2. Uncertainty around transition probabilities estimated from individual-level data
# 3. Costs by disease states
# ************************************************************************************************
source("C:/Users/xinxia/OneDrive - Karolinska Institutet/CSF-registersamkörning/AD DMT HE model/Statistical analyses/New analyses_XX_202412/HE model_PSA_function_XX_202503.R")


# 6.1. Summarize results and draw graphs ----
ncpus = parallel::detectCores()-1
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod),
                   library(flexsurv),library(parallel)))
clusterExport(cl, c(ls(all.names = TRUE)))


n_sim<-1:5 # add a sequence indicating number of simulations

sim_psa<-parLapply(cl,n_sim,function(j) {
  results<-lapply(1,function(i){  # Changed vector range to 1:nrow(pop_perc)
    
    sim.seed<-2025 + j
    rx_cycles_i<-pop_perc[i,]$rx_cycles
    age_group_i<-pop_perc[i,]$age_group
    age_i<-pop_perc[i,]$age
    sex_i<-pop_perc[i,]$sex
    apoe_i<-pop_perc[i,]$apoe
    stage_i<-pop_perc[i,]$stage
    perc<-pop_perc[i,]$pop_perc
    
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
        cbind(rx_cycles=rx_cycles_i,
              age_group=age_group_i, age=age_i, 
              sex=sex_i, apoe=apoe_i, stage=stage_i, perc=perc,
              error_message=NA)
      
    }, error=function(e) {
      data.frame(rx_cycles=rx_cycles_i,
                 age_group=age_group_i, age=age_i, 
                 sex=sex_i, apoe=apoe_i, stage=stage_i, perc=perc,
                 error_message=as.character(e$message)) 
    })
    return(psa_summary %>% cbind(sim=j))
  })
  bind_rows(results)
}) %>% bind_rows()










# 6.2. Threshold cost-effective price ----
WTP=seq(100000, 3000000, 100000)

WTP_results<-cross_join(res_psa_table_wide, data.frame(WTP=WTP)) %>%
  mutate(NMB=incr_qaly*WTP-incr_cost,
         threshold_price=NMB/rx.time_rx)

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





# 7. Scenario analysis ----
# ************************************************************************************************
# 1. Assuming effect modifications
# 2. Assuming treating for only 18 months or a sequence of 1:9 years
# 3. Assuming treating for whole time horizon, but with Waning factors = 0.5 when treatment stops
# 4. APOE heterozygotes vs. non-carriers
# 5. Discount rate of 3% for costs and 0 for health effects 
# ************************************************************************************************

# create a function to summarise results
ce_summary_func<-function(model){
  results<-bind_rows(lapply(model, function(x) {
    data.frame(
      rx_cycles=x$rx_cycles,age_group = x$age_group, age = x$age, 
      sex = x$sex, apoe = x$apoe, stage = x$stage,perc = x$perc, 
      ly_mci_stand=x$ly_mci_stand,ly_mild_stand=x$ly_mild_stand,
      ly_mod_stand=x$ly_mod_stand,ly_sev_stand=x$ly_sev_stand,
      ly_mci_inst_stand=x$ly_mci_inst_stand,ly_mild_inst_stand=x$ly_mild_inst_stand,
      ly_mod_inst_stand=x$ly_mod_inst_stand,ly_sev_inst_stand=x$ly_sev_inst_stand,
      ly_death_stand=x$ly_death_stand,
      qaly_stand = x$qaly_stand,qaly_rx=x$qaly_rx,
      costs_stand = x$costs_stand,costs_rx=x$costs_rx,
      deffect = x$deffect, dcost = x$dcost, rxtime = x$rxtime
    )
  }))
 
  results<-results %>%
    mutate(across(ly_mci_stand:rxtime, 
                  .fns = list(w = ~ . * results$perc),
                  .names = "{col}_{fn}" ) ) %>%
    mutate(NMB=deffect*QALYval-dcost,
           price =NMB/rxtime)
  
  return(results)
}


# create a function to organize results that can be exported to excel
ce_summary_excel_func<-function(summary){
  results_excel<-summary %>% 
    select(-c(age,perc,ly_mci_stand_w:price)) %>% 
    select(rx_cycles,rxtime,age_group:dcost) %>% 
    mutate(sex=factor(sex,levels=0:1,labels=c("Male","Female")),
           apoe=factor(apoe,levels=0:1,labels=c("Non-carrier","Heterozygotes")),
           stage=factor(stage,levels=0:1,labels=c("MCI","Mild AD"))) %>% 
    set_names(c("Treatment cycle","Treatment duration (years)","Age","Sex","APOE genotype",
                "AD stage","QALY-SOC","QALY-treatment","Costs-SOC","Costs-treatment",
                "QALY gain","Incremental costs"))
}



# call functions for scenario analyses
source("HE model_scenario analyses_function_XX_202502.R")


# set cluster for parellel computing
clusterExport(cl, c(ls(all.names = TRUE)))

# 7.1. Assume effect modifications in subgroups ----
sub_rx_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_ly_stand<-temp_model$eval_strategy_list$standard$counts %>% summarize_all(~sum(.))*cycle_length
  temp_ly_rx<-temp_model$eval_strategy_list$rx$counts %>% summarize_all(~sum(.))*cycle_length
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              ly_mci_stand=temp_ly_stand[["MCI"]],ly_mci_rx=temp_ly_rx[["MCI"]],
              ly_mild_stand=temp_ly_stand[["Mild"]],ly_mild_rx=temp_ly_rx[["Mild"]],
              ly_mod_stand=temp_ly_stand[["Moderate"]],ly_mod_rx=temp_ly_rx[["Moderate"]],
              ly_sev_stand=temp_ly_stand[["Severe"]],ly_sev_rx=temp_ly_rx[["Severe"]],
              ly_mci_inst_stand=temp_ly_stand[["MCI_inst"]],ly_mci_inst_rx=temp_ly_rx[["MCI_inst"]],
              ly_mild_inst_stand=temp_ly_stand[["Mild_inst"]],ly_mild_inst_rx=temp_ly_rx[["Mild_inst"]],
              ly_mod_inst_stand=temp_ly_stand[["Moderate_inst"]],ly_mod_inst_rx=temp_ly_rx[["Moderate_inst"]],
              ly_sev_inst_stand=temp_ly_stand[["Severe_inst"]],ly_sev_inst_rx=temp_ly_rx[["Severe_inst"]],
              ly_death_stand=temp_ly_stand[["Death"]],ly_death_rx=temp_ly_rx[["Death"]],
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})


sub_rx_mod_summary<-ce_summary_func(sub_rx_mod)


qaly=sum(sub_rx_mod_summary$deffect_w)
cost=sum(sub_rx_mod_summary$dcost_w)
NMB=qaly*QALYval-cost
sub_rx_price=NMB/sum(sub_rx_mod_summary$rxtime_w)
sub_rx_price

sum(sub_rx_mod_summary$qaly_rx_w)
sum(sub_rx_mod_summary$qaly_stand_w)
sum(sub_rx_mod_summary$costs_rx_w)
sum(sub_rx_mod_summary$costs_stand_w)


sub_rx_mod_excel<-ce_summary_excel_func(sub_rx_mod_summary)
  

# 7.2. Assuming treatment for varying durations of treatment ----
rx_duration_mod<-parLapply(cl,c(6,seq(2:9)*12*cycle_length),function(j){
  pop_perc<-pop_perc %>% 
    filter(rx_cycles<=j) %>% 
    mutate(pop_perc=pop_perc/sum(pop_perc))
  
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
    
    temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
    temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
    
    temp_ly_stand<-temp_model$eval_strategy_list$standard$counts %>% summarize_all(~sum(.))*cycle_length
    temp_ly_rx<-temp_model$eval_strategy_list$rx$counts %>% summarize_all(~sum(.))*cycle_length
    
    temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
    temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
    
    temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
    temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
    
    temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
    
    return(list(rx_cycles=rx_cycles_i,
                age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
                deffect=temp_deffect,dcost=temp_dcost,
                ly_mci_stand=temp_ly_stand[["MCI"]],ly_mci_rx=temp_ly_rx[["MCI"]],
                ly_mild_stand=temp_ly_stand[["Mild"]],ly_mild_rx=temp_ly_rx[["Mild"]],
                ly_mod_stand=temp_ly_stand[["Moderate"]],ly_mod_rx=temp_ly_rx[["Moderate"]],
                ly_sev_stand=temp_ly_stand[["Severe"]],ly_sev_rx=temp_ly_rx[["Severe"]],
                ly_mci_inst_stand=temp_ly_stand[["MCI_inst"]],ly_mci_inst_rx=temp_ly_rx[["MCI_inst"]],
                ly_mild_inst_stand=temp_ly_stand[["Mild_inst"]],ly_mild_inst_rx=temp_ly_rx[["Mild_inst"]],
                ly_mod_inst_stand=temp_ly_stand[["Moderate_inst"]],ly_mod_inst_rx=temp_ly_rx[["Moderate_inst"]],
                ly_sev_inst_stand=temp_ly_stand[["Severe_inst"]],ly_sev_inst_rx=temp_ly_rx[["Severe_inst"]],
                ly_death_stand=temp_ly_stand[["Death"]],ly_death_rx=temp_ly_rx[["Death"]],
                qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
                costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
                rxtime=temp_rxtime,model=list(temp_model)))
  }) %>% bind_rows() %>% cbind(rx_duration=j)
}) %>% bind_rows() 




rx_duration_mod_summary<-ce_summary_func(rx_duration_mod)


rx_duration_mod_excel<-ce_summary_excel_func(rx_duration_mod_summary)


# 7.3. Treatment across time horizon with waning factor = 0.5 ----
wf0.5_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_ly_stand<-temp_model$eval_strategy_list$standard$counts %>% summarize_all(~sum(.))*cycle_length
  temp_ly_rx<-temp_model$eval_strategy_list$rx$counts %>% summarize_all(~sum(.))*cycle_length
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              ly_mci_stand=temp_ly_stand[["MCI"]],ly_mci_rx=temp_ly_rx[["MCI"]],
              ly_mild_stand=temp_ly_stand[["Mild"]],ly_mild_rx=temp_ly_rx[["Mild"]],
              ly_mod_stand=temp_ly_stand[["Moderate"]],ly_mod_rx=temp_ly_rx[["Moderate"]],
              ly_sev_stand=temp_ly_stand[["Severe"]],ly_sev_rx=temp_ly_rx[["Severe"]],
              ly_mci_inst_stand=temp_ly_stand[["MCI_inst"]],ly_mci_inst_rx=temp_ly_rx[["MCI_inst"]],
              ly_mild_inst_stand=temp_ly_stand[["Mild_inst"]],ly_mild_inst_rx=temp_ly_rx[["Mild_inst"]],
              ly_mod_inst_stand=temp_ly_stand[["Moderate_inst"]],ly_mod_inst_rx=temp_ly_rx[["Moderate_inst"]],
              ly_sev_inst_stand=temp_ly_stand[["Severe_inst"]],ly_sev_inst_rx=temp_ly_rx[["Severe_inst"]],
              ly_death_stand=temp_ly_stand[["Death"]],ly_death_rx=temp_ly_rx[["Death"]],
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})




wf0.5_mod_summary<-ce_summary_func(wf0.5_mod)


qaly=sum(wf0.5_mod_summary$deffect_w)
cost=sum(wf0.5_mod_summary$dcost_w)
NMB=qaly*QALYval-cost
wf0.5_price=NMB/sum(wf0.5_mod_summary$rxtime_w)
wf0.5_price

sum(wf0.5_mod_summary$qaly_rx_w)
sum(wf0.5_mod_summary$qaly_stand_w)
sum(wf0.5_mod_summary$costs_rx_w)
sum(wf0.5_mod_summary$costs_stand_w)


wf0.5_mod_excel<-ce_summary_excel_func(wf0.5_mod_summary)



# 7.4. Assuming a starting population of only APOE non-carriers or only heterozygotes ----
# 7.4.1. Assuming a starting population of only APOE non-carriers ----
apoe0_mod<-parLapply(cl,1:nrow(pop_perc %>% filter(apoe==0)),function(i){
  pop_perc<-pop_perc %>% 
    filter(apoe==0) %>% 
    mutate(pop_perc=pop_perc/sum(pop_perc))
  
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_ly_stand<-temp_model$eval_strategy_list$standard$counts %>% summarize_all(~sum(.))*cycle_length
  temp_ly_rx<-temp_model$eval_strategy_list$rx$counts %>% summarize_all(~sum(.))*cycle_length
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              ly_mci_stand=temp_ly_stand[["MCI"]],ly_mci_rx=temp_ly_rx[["MCI"]],
              ly_mild_stand=temp_ly_stand[["Mild"]],ly_mild_rx=temp_ly_rx[["Mild"]],
              ly_mod_stand=temp_ly_stand[["Moderate"]],ly_mod_rx=temp_ly_rx[["Moderate"]],
              ly_sev_stand=temp_ly_stand[["Severe"]],ly_sev_rx=temp_ly_rx[["Severe"]],
              ly_mci_inst_stand=temp_ly_stand[["MCI_inst"]],ly_mci_inst_rx=temp_ly_rx[["MCI_inst"]],
              ly_mild_inst_stand=temp_ly_stand[["Mild_inst"]],ly_mild_inst_rx=temp_ly_rx[["Mild_inst"]],
              ly_mod_inst_stand=temp_ly_stand[["Moderate_inst"]],ly_mod_inst_rx=temp_ly_rx[["Moderate_inst"]],
              ly_sev_inst_stand=temp_ly_stand[["Severe_inst"]],ly_sev_inst_rx=temp_ly_rx[["Severe_inst"]],
              ly_death_stand=temp_ly_stand[["Death"]],ly_death_rx=temp_ly_rx[["Death"]],
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})


apoe0_mod_summary<-ce_summary_func(apoe0_mod)


qaly=sum(apoe0_mod_summary$deffect_w)
cost=sum(apoe0_mod_summary$dcost_w)
NMB=qaly*QALYval-cost
apoe0_price=NMB/sum(apoe0_mod_summary$rxtime_w)
apoe0_price

sum(apoe0_mod_summary$qaly_rx_w)
sum(apoe0_mod_summary$qaly_stand_w)
sum(apoe0_mod_summary$costs_rx_w)
sum(apoe0_mod_summary$costs_stand_w)

apoe0_mod_excel<-ce_summary_excel_func(apoe0_mod_summary)

# 7.4.2. Assuming a starting population of only APOE heterozygotes ----
apoe1_mod<-parLapply(cl,1:nrow(pop_perc %>% filter(apoe==1)),function(i){
  pop_perc<-pop_perc %>% 
    filter(apoe==1) %>% 
    mutate(pop_perc=pop_perc/sum(pop_perc))
  
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_ly_stand<-temp_model$eval_strategy_list$standard$counts %>% summarize_all(~sum(.))*cycle_length
  temp_ly_rx<-temp_model$eval_strategy_list$rx$counts %>% summarize_all(~sum(.))*cycle_length
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              ly_mci_stand=temp_ly_stand[["MCI"]],ly_mci_rx=temp_ly_rx[["MCI"]],
              ly_mild_stand=temp_ly_stand[["Mild"]],ly_mild_rx=temp_ly_rx[["Mild"]],
              ly_mod_stand=temp_ly_stand[["Moderate"]],ly_mod_rx=temp_ly_rx[["Moderate"]],
              ly_sev_stand=temp_ly_stand[["Severe"]],ly_sev_rx=temp_ly_rx[["Severe"]],
              ly_mci_inst_stand=temp_ly_stand[["MCI_inst"]],ly_mci_inst_rx=temp_ly_rx[["MCI_inst"]],
              ly_mild_inst_stand=temp_ly_stand[["Mild_inst"]],ly_mild_inst_rx=temp_ly_rx[["Mild_inst"]],
              ly_mod_inst_stand=temp_ly_stand[["Moderate_inst"]],ly_mod_inst_rx=temp_ly_rx[["Moderate_inst"]],
              ly_sev_inst_stand=temp_ly_stand[["Severe_inst"]],ly_sev_inst_rx=temp_ly_rx[["Severe_inst"]],
              ly_death_stand=temp_ly_stand[["Death"]],ly_death_rx=temp_ly_rx[["Death"]],
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})




apoe1_mod_summary<-ce_summary_func(apoe1_mod)


qaly=sum(apoe1_mod_summary$deffect_w)
cost=sum(apoe1_mod_summary$dcost_w)
NMB=qaly*QALYval-cost
apoe1_price=NMB/sum(apoe1_mod_summary$rxtime_w)
apoe1_price

sum(apoe1_mod_summary$qaly_rx_w)
sum(apoe1_mod_summary$qaly_stand_w)
sum(apoe1_mod_summary$costs_rx_w)
sum(apoe1_mod_summary$costs_stand_w)

apoe1_mod_excel<-ce_summary_excel_func(apoe1_mod_summary)


# 7.5. Discounting rate 0.03 for costs and no discounting for QALY ----
disc_diff_mod<-parLapply(cl,1:nrow(pop_perc),function(i){
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
  
  temp_deffect<-summary(temp_model)$res_comp$.deffect[2]
  temp_dcost<-summary(temp_model)$res_comp$.dcost[2]
  
  temp_ly_stand<-temp_model$eval_strategy_list$standard$counts %>% summarize_all(~sum(.))*cycle_length
  temp_ly_rx<-temp_model$eval_strategy_list$rx$counts %>% summarize_all(~sum(.))*cycle_length
  
  temp_qaly_stand<-sum(temp_model$eval_strategy_list$standard$values$utility)
  temp_qaly_rx<-sum(temp_model$eval_strategy_list$rx$values$utility)
  
  temp_costs_stand<-sum(temp_model$eval_strategy_list$standard$values$cost)
  temp_costs_rx<-sum(temp_model$eval_strategy_list$rx$values$cost)
  
  temp_rxtime<-sum(temp_model$eval_strategy_list$rx$values$rx.time) # accumulated time on treatment
  
  return(list(rx_cycles=rx_cycles_i,
              age_group=age_group_i,age=age_i,sex=sex_i,apoe=apoe_i,stage=stage_i, perc=perc,
              deffect=temp_deffect,dcost=temp_dcost,
              ly_mci_stand=temp_ly_stand[["MCI"]],ly_mci_rx=temp_ly_rx[["MCI"]],
              ly_mild_stand=temp_ly_stand[["Mild"]],ly_mild_rx=temp_ly_rx[["Mild"]],
              ly_mod_stand=temp_ly_stand[["Moderate"]],ly_mod_rx=temp_ly_rx[["Moderate"]],
              ly_sev_stand=temp_ly_stand[["Severe"]],ly_sev_rx=temp_ly_rx[["Severe"]],
              ly_mci_inst_stand=temp_ly_stand[["MCI_inst"]],ly_mci_inst_rx=temp_ly_rx[["MCI_inst"]],
              ly_mild_inst_stand=temp_ly_stand[["Mild_inst"]],ly_mild_inst_rx=temp_ly_rx[["Mild_inst"]],
              ly_mod_inst_stand=temp_ly_stand[["Moderate_inst"]],ly_mod_inst_rx=temp_ly_rx[["Moderate_inst"]],
              ly_sev_inst_stand=temp_ly_stand[["Severe_inst"]],ly_sev_inst_rx=temp_ly_rx[["Severe_inst"]],
              ly_death_stand=temp_ly_stand[["Death"]],ly_death_rx=temp_ly_rx[["Death"]],
              qaly_stand=temp_qaly_stand,qaly_rx=temp_qaly_rx,
              costs_stand=temp_costs_stand,costs_rx=temp_costs_rx,
              rxtime=temp_rxtime,model=list(temp_model)))
})




disc_diff_mod_summary<-ce_summary_func(disc_diff_mod)


qaly=sum(disc_diff_mod_summary$deffect_w)
cost=sum(disc_diff_mod_summary$dcost_w)
NMB=qaly*QALYval-cost
disc_diff_price=NMB/sum(disc_diff_mod_summary$rxtime_w)
disc_diff_price

sum(disc_diff_mod_summary$qaly_rx_w)
sum(disc_diff_mod_summary$qaly_stand_w)
sum(disc_diff_mod_summary$costs_rx_w)
sum(disc_diff_mod_summary$costs_stand_w)

disc_diff_mod_excel<-ce_summary_excel_func(disc_diff_mod_summary)


stopCluster(cl)


# export summary to excel
writexl::write_xlsx(list(base_mod=base_mod_excel,subgroup_mod=sub_rx_mod_excel,
                         rx_duration_mod=rx_duration_mod_excel,wf0.5_mod=wf0.5_mod_excel,
                         apoe0_mod=apoe0_mod_excel,apoe1_mod=apoe1_mod_excel,
                         disc_diff_mod=disc_diff_mod_excel),
                    "Model_summary.xlsx")