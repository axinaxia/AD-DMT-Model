rm(list=ls())      # clear environment
cat("\014")        # clear console

# Load libraries----

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

#Load data ----

path<-'C:/Users/xinxia/OneDrive - Karolinska Institutet/Postdoc project plan/Postdoc project - obligation/AD DMT HE model_XX/Statistical analyses'
setwd(path)
#load('events.RData')
load('costtable.RData')
load('msdata.RData')

statenames = c(
  "Questionnable",
  "Mild",
  "Moderate",
  "Severe",
  "Questionnable_inst",
  "Mild_inst",
  "Moderate_inst",
  "Severe_inst",
  "Death")

#Functions----


gethaz<-function(df) {
  
  cx<-runcox(df)
  
  b<-as.data.frame(basehaz(fit=cx, centered=F)) %>%
    
    #Add cumhaz zero at time zero
    rbind(data.frame(time=0, hazard=0, strata=paste0('trans=',1:18))) %>%
    mutate(trans=as.numeric(strata)) %>% arrange(trans, time)
  
  #Set keys for rolling join
  b<-as.data.table(b)
  setkey(b, trans, time)
  
  
  #Create vector of time points of interest (x-month cycles), one for each transition
  
  timesvec<-cross_join(data.frame(trans=1:18), data.frame(time=round((365.25/12*cycle_length))*(0:cycles))) %>% as.data.table
  setkey(timesvec, trans, time)
  
  #Rolling join to get cumulative hazard at each half-year cycle. Calculate hazard over each cycle
  cumhaz<-b[timesvec, roll=T] %>%
    group_by(trans) %>%
    mutate(cumhaz=ifelse(is.na(hazard), 0, hazard),
           haz=cumhaz-lag(cumhaz),
           cycle=row_number()-1) %>% drop_na %>%
    ungroup %>%
    select(trans, cycle, haz) 
  r<-list()
  r$coef<-cx$coefficients
  r$haz=cumhaz
  return(r)
}

runcox<-function(df) {
  cx <- coxph(Surv(Tstart,Tstop,status)~
                AGE.1+AGE.2+AGE.4+AGE.5+AGE.6+AGE.8+AGE.9+AGE.10+AGE.15+AGE.17+
                SEXMALE.1+SEXMALE.2+SEXMALE.4+SEXMALE.5+SEXMALE.6+SEXMALE.8+SEXMALE.9+SEXMALE.10+SEXMALE.15+SEXMALE.17+
                strata(trans),
              data=df,method="breslow")
  return(cx)
}

get_tp<-function(transition_number, sex, age, hazards, coefficients) {
  #Get hazard for each transition at each time point
  if (paste0('AGE.', transition_number) %in% names(coefficients)) age_coef<-as.numeric(coefficients[names(coefficients)==paste0('AGE.', transition_number)]) else age_coef<-0
  if (paste0('SEXMALE.', transition_number) %in% names(coefficients)) sex_coef<-as.numeric(coefficients[names(coefficients)==paste0('SEXMALE.', transition_number)]) else sex_coef<-0
  h=hazards %>% filter(trans==transition_number)   #Base hazard without covariates
  hr=exp(age_coef*age+sex_coef*sex) #Hazard ratio
  h$p=1-exp(-h$haz*hr) #Probability of transition ??? I don't think this is correct
  return(h$p) #return vector of transition probabilities by cycle
}

getcost<-function(df) {
  costtable<-df %>% group_by(state, setting) %>%
    mutate(state=factor(state, levels=c('MMSE 26-30', 'MMSE 21-25', 'MMSE 10-20', 'MMSE 0-9')),
           setting=factor(setting, levels=c('Community', 'Institution'))
           ) %>%
    summarise(cost=sum(totalcost), 
              time=sum(sumtime),
              annual_cost=cost/time*365, .groups='drop') %>%
    arrange(setting, state) 
  return(costtable$annual_cost)
}

utilities<-c(0.8, 0.6, 0.4, 0.2, 0.8, 0.6, 0.4, 0.2)
cycle_length=3 #in months
cycles=12/cycle_length*10 #10 year time horizon

# Run model----


runmodel<-function(h, #Hazards
                   cost_df,#cost data
                   utility_df, #utility data
                   rx_dur=6, #treatment duration in cycles (18 month?)
                   rx_eff=0.7, #treatment effect
                   age_init=70, 
                   sex=0,
                   r=0.015 #discount rate
                   )
{
  #h<-gethaz(df) #retrieve hazards and coefficients
  
  tp=lapply(1:18, function(x) get_tp(x, sex, age_init, h$haz, h$coef)) #get transition probabilities per transition and cycle
  
  param <- define_parameters(
    age=age_init+model_time/(12/cycle_length),
    disc=exp(-r*model_time/(12/cycle_length)), #continuously compounded discounting factor
    

    #Variable to vary treatment effect by arm and over time
    rX=dispatch_strategy(
      standard=1,
      rx=case_when(model_time<rx_dur ~ rx_eff, T ~ 1)
    ),
    
    #Variable to use to record time on treatment
    rxtime=dispatch_strategy(
      standard=0,
      rx=case_when(model_time<rx_dur ~ 1, T ~ 0)
    ),
    
    # transition probabilities
    t1=tp[[1]][model_time],
    t2=tp[[2]][model_time],
    t3=tp[[3]][model_time],
    t4=tp[[4]][model_time],
    t5=tp[[5]][model_time],
    t6=tp[[6]][model_time],
    t7=tp[[7]][model_time],
    t8=tp[[8]][model_time],
    t9=tp[[9]][model_time],
    t10=tp[[10]][model_time],
    t11=tp[[11]][model_time],
    t12=tp[[12]][model_time],
    t13=tp[[13]][model_time],
    t14=tp[[14]][model_time],
    t15=tp[[15]][model_time],
    t16=tp[[16]][model_time],
    t17=tp[[17]][model_time],
    t18=tp[[18]][model_time]
  )
  
  # Define states, costs and utilities----
  Questionnable = define_state(
    utility = utility_df[1]/(12/cycle_length)*disc,
    cost = (cost_df[1]/(12/cycle_length))*disc,
    rx.time=rxtime/(12/cycle_length),
    rx.time.disc=rx.time*disc
  )
  Mild = define_state(
    utility = utility_df[2]/(12/cycle_length)*disc,
    cost = (cost_df[2]/(12/cycle_length))*disc,
    rx.time=rxtime/(12/cycle_length),
    rx.time.disc=rx.time*disc
  )
  Moderate = define_state(
    utility = utility_df[3]/(12/cycle_length)*disc,
    cost = (cost_df[3]/(12/cycle_length))*disc,
    rx.time=0,
    rx.time.disc=0
  )
  Severe = define_state(
    utility = utility_df[4]/(12/cycle_length)*disc,
    cost = cost_df[4]/(12/cycle_length)*disc,
    rx.time=0,
    rx.time.disc=0
  )
  Questionnable_inst = define_state(
    utility = utility_df[5]/(12/cycle_length)*disc,
    cost = (cost_df[5]/(12/cycle_length))*disc,
    rx.time=rxtime/(12/cycle_length),
    rx.time.disc=rx.time*disc
  )
  Mild_inst = define_state(
    utility = utility_df[6]/(12/cycle_length)*disc,
    cost = (cost_df[6]/(12/cycle_length))*disc,
    rx.time=rxtime/(12/cycle_length),
    rx.time.disc=rx.time*disc
  )
  Moderate_inst = define_state(
    utility = utility_df[7]/(12/cycle_length)*disc,
    cost = (cost_df[7]/(12/cycle_length))*disc,
    rx.time=0,
    rx.time.disc=0
  )
  Severe_inst = define_state(
    utility = utility_df[8]/(12/cycle_length)*disc,
    cost = cost_df[8]/(12/cycle_length)*disc,
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
    C, t1*rX, 0,     0,     t2, 0,      0,      0,      t3,
    0, C,     t4*rX, 0,     0,  t5,     0,      0,      t6, 
    0, 0,     C,     t7,    0,  0,      t8,     0,      t9, 
    0, 0,     0,     C,     0,  0,      0,      t10,    t11,
    0, 0,     0,     0,     C,  t12*rX, 0,      0,      t13,
    0, 0,     0,     0,     0,  C,      t14*rX, 0,      t15,
    0, 0,     0,     0,     0,  0,      C,      t16,    t17, 
    0, 0,     0,     0,     0,  0,      0,      C,      t18,
    0, 0,     0,     0,     0,  0,      0,      0,      C
  )
  
  # Strategies----
  strat_standard <- define_strategy(
    Questionnable=Questionnable,
    Mild=Mild,
    Moderate=Moderate,
    Severe=Severe,
    Questionnable_inst=Questionnable_inst,
    Mild_inst=Mild_inst,
    Moderate_inst=Moderate_inst,
    Severe_inst=Severe_inst,
    Death=Death,
    transition = transmat)
  
  init_state<-c(1,0,0,0,0,0,0,0,0)
    
  
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

#Base-case scenario-----
res_mod<-runmodel(gethaz(msdata), getcost(costtable), utilities)


QALYval=1000000
sum<-summary(res_mod)
NMB=sum$res_comp$.deffect[2]*QALYval-sum$res_comp$.dcost[2]
rxtime<-sum(res_mod$eval_strategy_list$rx$values$rx.time)
threshold.rxcost=NMB/rxtime

#Sojourn times
sojourn<-rbind(res_mod$eval_strategy_list$standard$counts %>% summarize_all(~sum(.)) %>% mutate(strategy='standard'),
               res_mod$eval_strategy_list$rx$counts %>% summarize_all(~sum(.)) %>% mutate(strategy='rx')
) %>% pivot_longer(cols=-strategy, names_to='state', values_to = 'years') %>% mutate(years=years/4,
                                                                                     state=factor(state, levels=rev(statenames))) 
ggplot(data=sojourn, aes(y=years, x=strategy, fill=state)) + geom_col(position='stack')+coord_flip()

sojourn_diff<-sojourn %>% 
  pivot_wider(id_cols=state, names_from='strategy', values_from='years') %>%
  mutate(years=rx-standard, 
         strategy='difference') %>% select(state, years, strategy) 

ggplot(data=sojourn_diff, aes(y=years, x=state, fill=state)) + 
  geom_col(position='stack')+coord_flip()


#c2<-get_counts(res_mod) %>% group_by(state_names) %>% summarize(sojourn_time=sum(count))
cohdist<-res_mod$eval_strategy_list$standard$counts
cohdist$t=1:40
cohdist$dist_Questionnable<-(cohdist$Questionnable+cohdist$Questionnable_inst)/(1-cohdist$Death)
cohdist$dist_Mild<-(cohdist$Mild+cohdist$Mild_inst)/(1-cohdist$Death)
cohdist$dist_Moderate<-(cohdist$Moderate+cohdist$Moderate_inst)/(1-cohdist$Death)
cohdist$dist_Severe<-(cohdist$Severe+cohdist$Severe_inst)/(1-cohdist$Death)

cohdist$dist_Inst<-(cohdist$Questionnable_inst+cohdist$Mild_inst+cohdist$Moderate_inst+cohdist$Severe_inst)/(1-cohdist$Death)
coh.inst<-pivot_longer(cohdist %>% select(t, dist_Inst),cols=starts_with('dist'),names_to='State'
)

ggplot(aes(x=t, y=value, color=State), data=coh.inst)+geom_line()+geom_point()

res_mod$eval_strategy_list$standard$parameters$t1
res_mod$eval_strategy_list$standard$counts

costplot<-rbind(data.frame(strategy='standard', 
                           Markov.cycle=res_mod$eval_strategy_list$standard$values$model_time, 
                           cost=res_mod$eval_strategy_list$standard$values$cost),
                data.frame(strategy='rx', 
                           Markov.cycle=res_mod$eval_strategy_list$standard$values$model_time, 
                           cost=res_mod$eval_strategy_list$rx$values$cost))

ggplot(aes(y=cost, x=Markov.cycle, color=strategy), data=costplot) + geom_line()+geom_point()                


#Sensitivity analysis----


se <- define_dsa(
  drugcost, 50000, 500000
)

res_dsa<-run_dsa(
  model=res_mod,
  dsa=se
)

res_mod

plot(res_dsa, 'difference')


#Model calibration

#Extract values

extract_values<-function(x) {
  d<-get_counts(x) %>% filter(state_names=='Death') %>% group_by(.strategy_names) %>%
    summarize(surv=sum(count)) 
  return(d[1,]$surv-d[2,]$surv)
}


#Calibrate overall survival
res_cal <- calibrate_model(
  res_mod,
  parameter_names = 'mortcal',
  fn_values = extract_values,
  target_values = 0,
  lower=1,
  upper=10
)
extract_values(res_mod)


#Bootstrapping PSA-----
patients<-msdata$LOPNR %>% unique

sample_patient<-function() {
  pass=FALSE
  attempts=0
  while(pass==FALSE) {
    s<-data.frame(LOPNR=sample(patients, replace=T))
    msdata_sample<-msdata %>% inner_join(s, relationship = 'many-to-many', by = join_by(LOPNR))
    cost_sample<-getcost(costtable %>% inner_join(s, relationship = 'many-to-many', by = join_by(LOPNR)))
    res<-try(gethaz(msdata_sample))
    if(length(cost_sample)==8 & sum(is.na(res$coef))==0 & !inherits(res, "try-error")) {
      pass<-TRUE
    }
    attempts=attempts+1
  }
  return(list(res, cost_sample))
}


run_bs<-function(iteration) {
  set.seed(iteration+1000)
  s<-sample_patient()
  res<-try(runmodel(s[[1]], s[[2]], utilities))
  if(!inherits(res, "try-error")) {
    out<-summary(res)$res_values
    out$iteration=iteration
    return(out)
  }
}



psa<-pblapply(1:100, run_bs) %>% bind_rows

results<-psa %>% 
  pivot_wider(id_cols=iteration, names_from=.strategy_names, values_from=c(cost, utility, rx.time.disc)) %>%
  mutate(incr_cost=cost_rx-cost_standard, 
         incr_utility=utility_rx-utility_standard) 

WTP=seq(100000, 3000000, 100000)

WTP_results<-cross_join(results, data.frame(WTP=WTP)) %>%
  mutate(NMB=incr_utility*WTP-incr_cost,
           threshold_price=NMB/rx.time.disc_rx)

CEAC<-WTP_results %>% group_by(WTP) %>%  
  reframe(tibble::enframe(quantile(threshold_price,c(0.05, 0.5, 0.95)),'quantile', 'threshold_price' )) %>%
  pivot_wider(id_cols=WTP, names_from=quantile, values_from=threshold_price)
  

ggplot()+
  #geom_point(aes(x=WTP, y=threshold_price), data=WTP_results,position='jitter', size=0.5)+
  geom_ribbon(aes(ymin=`5%`, ymax=`95%`, x=WTP), data=CEAC, alpha=0.2, fill='red'#, outline.type = 'full'
              )+
  geom_line(aes(y=`50%`, x=WTP), data=CEAC)+
  scale_x_continuous(labels = scales::label_comma())+
  scale_y_continuous(labels = scales::label_comma())+
  ylab('Threshold cost-effective price')+
  xlab('Willingness to pay per QALY')
#CE scatterplot

library(ConfidenceEllipse)

results<-results %>% 
  mutate(incr_cost=incr_cost+rx.time.disc_rx*78087, #???
         cost_effective=if_else(incr_cost/incr_utility < 1000000, 'Cost-effective', 'Not cost-effective') )

ellipse_95 <- confidence_ellipse(results, x = incr_utility, y = incr_cost, conf_level = 0.95)

icer_conf<-quantile(results$incr_cost/results$incr_utility, c(0.025, 0.0975))

results %>% 
  ggplot(aes(y=incr_cost, x=incr_utility))+
  geom_point(aes(color=cost_effective), size=0.5
             )+
  xlim(pmin(0, min(results$incr_utility)), pmax(0,max(results$incr_utility), max(ellipse_95$x)))+
  ylim(pmin(0, min(results$incr_cost)), pmax(0,max(results$incr_cost+results$rx.time.disc_rx), max(ellipse_95$y)))+
  geom_path(data = ellipse_95, aes(x = x, y = y), color = "black", linewidth = 0.5) +
  geom_abline(intercept = 0, slope = icer_conf[1], linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = icer_conf[2], linetype = "dashed", color = "black")+
    theme_minimal()




# library(dampack)
# 
# strategies <- c("standard", "rx")
# # Vector of WTP thresholds
# v.wtp <- seq(1000, 1500000, by = 100000)
# # Matrix of costs
# m.c <- psa %>% pivot_wider(id_cols=iteration, names_from=.strategy_names, values_from=cost) %>% select(-1)
# # Matrix of effectiveness
# m.e <- psa %>% pivot_wider(id_cols=iteration, names_from=.strategy_names, values_from=utility) %>% select(-1)
# # Compute CEAF
# out <- ceaf(v.wtp = v.wtp, strategies = strategies, 
#             m.e = m.e , m.c = m.c,
#             ceaf.out = TRUE)
# # Plot CEAF
# out$gg.ceaf

single_res<-summary(runmodel(msdata, c(0,0,getcost(costtable)), utilities))$res_values %>%
  mutate(
    strategy=factor(.strategy_names, levels=c('standard', 'rx')),
    cost=if_else(strategy=='rx', cost+rx.time.disc*150000, cost)
  ) %>% arrange(strategy) %>% 
  mutate(incr_cost=cost-lag(cost), 
         incr_utility=utility-lag(utility),
         ICER=incr_cost/incr_utility) 


single_res %>%   
  pivot_wider(id_cols=.n_indiv, 
              names_from=.strategy_names, 
              values_from=c(cost, utility, rx.time.disc)
              ) %>%
  mutate(incr_cost=cost_rx-cost_standard+150000, 
         incr_utility=utility_rx-utility_standard) 


  ggplot(aes(y=cost, x=utility), data=single_res)+
  geom_point()+
  geom_line()+
  geom_text(aes(label=.strategy_names), 
            nudge_x=0.001, 
            #nudge_y=0.01
            hjust=0
            )+
    geom_text(aes(x=xmean, y=ymin, label=paste('Incremental cost ', round(incr_cost))), data=
                single_res %>% 
                summarise(incr_cost=max(incr_cost, na.rm=T), 
                          xmean=mean(utility), 
                          ymean=mean(cost),
                          xmax=max(utility), 
                          ymin=min(cost))
    )+
  geom_text(aes(x=xmax, y=ymean, label=paste('Incremental utility ', round(incr_utility,3))), data=
              single_res %>% 
              summarise(incr_cost=max(incr_cost, na.rm=T), 
                        incr_utility=max(incr_utility, na.rm=T), 
                        xmean=mean(utility), 
                        ymean=mean(cost),
                        xmax=max(utility), 
                        ymin=min(cost)
                    ),
            angle=90
  )+
    geom_text(aes(x=xmean+0.01, y=ymean, label=paste('ICER ', round(ICER))), data=
                single_res %>% 
                summarise(ICER=max(ICER, na.rm=T), 
                          xmean=mean(utility), 
                          ymean=mean(cost),
                          xmax=max(utility), 
                          ymin=min(cost)
                )
    )
    
  



QALYval=1000000
sum<-summary(res_mod)
NMB=sum$res_comp$.deffect[2]*QALYval-sum$res_comp$.dcost[2]
rxtime<-sum(res_mod$eval_strategy_list$rx$values$rx.time)
threshold.rxcost=NMB/rxtime

# Process results-----


