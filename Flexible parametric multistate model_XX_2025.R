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
library(ggpubr)
library(rlang)
library(mstate)
library(table1)


# check baseline characteristics of SveDem
basechar_svedem<-mstatedata %>% 
  filter(Tstart==0) %>% 
  distinct(id,.keep_all = T) %>% 
  mutate(base_state=case_when(from==1|from==4~1,
                              from==2|from==5~2,
                              from==3|from==6~3))

table1::table1(~AGE+factor(base_state)|SEX,data = basechar_svedem)

# check baseline characteristics of NACC
table1::table1(~AGE+factor(NACCNE4S)|SEX,data = nacc_mci_select)


# ************************************************************************************************
# Steps:
# 1. Run proportional-hazard Royston-Parmar models for each transition
# 2. Compare cumulative hazard from the parametric models with the non-parametric model
# ************************************************************************************************
load("mstatedata_xx.RData")

mstatedata<-mstatedata
# check transitions in the data 
events(mstatedata)
table(mstatedata$from,mstatedata$to)


transmat<-table(eventdata$state,eventdata$next_state) %>% 
  as.data.frame() %>% 
  filter(Freq!=0) %>% 
  dplyr::select(-Freq) %>% 
  set_names(c("state","poss_next_state")) %>% 
  mutate(state=as.numeric(state),
         poss_next_state=as.numeric(poss_next_state)) %>% 
  arrange(state,poss_next_state) %>% 
  filter(state!=poss_next_state) %>% 
  mutate(trans=1:n())

# 1. Non-parametric model and Royston-Parmar model for each transition ----
rm(list = grep("_trans", ls(), value = TRUE))

for (i in 1:13){
  newdata<-mstatedata %>% 
    filter(trans==i)

  # run a Cox model
  model<-coxph(Surv(Tstart,Tstop,status) ~ 1, data=newdata, method="breslow")
  
  # estimate cumulative hazard
  cumhaz<-basehaz(model)
  
  # return one model for each transition
  assign(paste0("nonpara_model_trans",i),model)
  assign(paste0("nonpara_cumhaz_trans",i),cumhaz)
}



rm(newdata,model,cumhaz) # remove previous model outputs


for (i in 1:13){
  newdata<-mstatedata %>% 
    filter(trans==i)
  
for (j in 1:6){
  model<-tryCatch({
    flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                   data=newdata, k=j,scale="hazard",
                   timescale="log")
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with ",j," knots."))
  })
  
  cumhaz<-tryCatch({
    (predict(model, type = "cumhaz", newdata = data.frame(AGE = NA),
             times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
      rename(.time=.eval_time)
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with ",j," knots."))
  })  
  
  assign(paste0("RP_model_trans",i,"_with_",j,"knots"),model)
  assign(paste0("RP_cumhaz_trans",i,"_with_",j,"knots"),cumhaz)
}
}



# for transitions with difficulty converging, run simpler models first
for (i in 7:13){
  newdata<-mstatedata %>% 
    filter(trans==i)
  
  # model with 0 knots
  model0<-flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                         data=newdata, k=0,scale="hazard",
                         timescale="log")
  
  # model with 1 knots
  model1 <- tryCatch({
    flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                   data=newdata, k=1,scale="hazard",
                   inits=c(coef(model0), 0),
                   timescale="log")
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 1 knots."))
  })
  
  # estimate cumulative hazard
  cumhaz1 <- tryCatch({
    (predict(model1, type = "cumhaz", newdata = data.frame(AGE = NA),
             times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
      rename(.time=.eval_time)
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 1 knots."))
  })
  
  # model with 2 knots
  model2 <- tryCatch({
    flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                   data=newdata, k=2,scale="hazard",
                   inits=c(coef(model1), 0),
                   timescale="log")
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 2 knots."))
  })
  
  # estimate cumulative hazard
  cumhaz2 <- tryCatch({
    (predict(model2, type = "cumhaz", newdata = data.frame(AGE = NA),
             times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
      rename(.time=.eval_time)
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 2 knots."))
  })
  
  # model with 3 knots
  model3 <- tryCatch({
    flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                   data=newdata, k=3,scale="hazard",
                   inits=c(coef(model2), 0),
                   timescale="log")
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 3 knots."))
  })
  
  # estimate cumulative hazard
  cumhaz3 <- tryCatch({
    (predict(model3, type = "cumhaz", newdata = data.frame(AGE = NA),
             times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
      rename(.time=.eval_time)
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 3 knots."))
  })
  
  # model with 4 knots
  model4 <- tryCatch({
    flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                   data=newdata, k=4,scale="hazard",
                   inits=c(coef(model3), 0),
                   timescale="log")
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 4 knots."))
  })
  
  # estimate cumulative hazard
  cumhaz4 <- tryCatch({
    (predict(model4, type = "cumhaz", newdata = data.frame(AGE = NA),
             times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
      rename(.time=.eval_time)
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 4 knots."))
  })
  
  # model with 5 knots
  model5 <- tryCatch({
    flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                   data=newdata, k=5,scale="hazard",
                   inits=c(coef(model4), 0),
                   timescale="log")
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 5 knots."))
  })
  
  # estimate cumulative hazard
  cumhaz5 <- tryCatch({
    (predict(model5, type = "cumhaz", newdata = data.frame(AGE = NA),
             times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
      rename(.time=.eval_time)
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 5 knots."))
  })
  
  
  # model with 6 knots
  model6 <- tryCatch({
    flexsurvspline(Surv(Tstart,Tstop,status) ~ 1,
                   data=newdata, k=6,scale="hazard",
                   inits=c(coef(model4), 0),
                   timescale="log")
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 6 knots."))
  })
  
  # estimate cumulative hazard
  cumhaz6 <- tryCatch({
    (predict(model6, type = "cumhaz", newdata = data.frame(AGE = NA),
             times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
      rename(.time=.eval_time)
  }, error = function(e) {
    message(paste("Error for trans =", i, "with RP model with 6 knots."))
  })
  
  assign(paste0("RP_model_trans",i,"_with_1knots"),model1)
  assign(paste0("RP_cumhaz_trans",i,"_with_1knots"),cumhaz1)
  assign(paste0("RP_model_trans",i,"_with_2knots"),model2)
  assign(paste0("RP_cumhaz_trans",i,"_with_2knots"),cumhaz2)
  assign(paste0("RP_model_trans",i,"_with_3knots"),model3)
  assign(paste0("RP_cumhaz_trans",i,"_with_3knots"),cumhaz3)
  assign(paste0("RP_model_trans",i,"_with_4knots"),model4)
  assign(paste0("RP_cumhaz_trans",i,"_with_4knots"),cumhaz4)
  assign(paste0("RP_model_trans",i,"_with_5knots"),model5)
  assign(paste0("RP_cumhaz_trans",i,"_with_5knots"),cumhaz5)
  assign(paste0("RP_model_trans",i,"_with_6knots"),model6)
  assign(paste0("RP_cumhaz_trans",i,"_with_6knots"),cumhaz6)
  
}


rm(newdata,model0,model1,model2,model3,model4,model5,
   cumhaz1,cumhaz2,cumhaz3,cumhaz4,cumhaz5)





# 2. Assess model fitness ----
nonpara_cumhaz_combined<-do.call(rbind, lapply(1:13, function(i) {
  cumhaz_df<- get(paste0("nonpara_cumhaz_trans", i)) %>% 
    set_names(c("cumhaz","time"))
  cumhaz_df$trans <- i   
  cumhaz_df$model<-"Non-parametric"
  return(cumhaz_df)
}))


rp_cumhaz_combined <- do.call(rbind, lapply(1:13, function(i) {
  do.call(rbind, lapply(1:5, function(k) {
    cumhaz_df <- get(paste0("RP_cumhaz_trans", i, "_with_", k, "knots"))
    
    if (!is.null(cumhaz_df)) {
      cumhaz_df <- cumhaz_df %>%
        set_names(c("time", "cumhaz", "cumhaz_lb", "cumhaz_ub"))
    } else {
      cumhaz_df <- data.frame(time = NA, cumhaz = NA, cumhaz_lb = NA, cumhaz_ub = NA)
    }
    
    cumhaz_df$trans <- i
    cumhaz_df$model <- paste0("RP model with ", k, " knots")
    return(cumhaz_df)
  }))
}))



final_cumhaz_combined<-rbind(nonpara_cumhaz_combined %>% 
                               mutate(cumhaz_lb=NA,
                                      cumhaz_ub=NA), rp_cumhaz_combined) %>% 
  mutate(model=factor(model),
         trans=factor(trans))



model_fitness<-do.call(rbind, lapply(c(1:13), function(i) {
  do.call(rbind, lapply(1:5, function(k) {
  rpfit <- get(paste0("RP_model_trans", i, "_with_", k, "knots")) %>% glance()
  rpfit$model <- "Royston and Parmar" 
  rpfit$trans <- i           
  rpfit$knots<-k
  return(rpfit)
  }))
}))


model_fitness <- model_fitness %>%
  arrange(trans,BIC,AIC)



for (i in 1:13){
  assign(paste0("plot_trans",i),
         ggplot(final_cumhaz_combined %>% 
                  filter(trans==i), aes(x=time, y=cumhaz, color=model))+
           geom_line(aes(color=model))+
           scale_color_manual(values = c("black","#CC79A7",
                                "#D55E00","#56B4E9","#0072B2","#009E73"))+
           scale_x_continuous(breaks = seq(0,10,2),expand=c(0,0))+
           theme_minimal()+
           labs(title=paste("Transition",i),
                x='Time (years)',
                y='Cumulative hazard',
                linetype='Distribution',
                color="Transition")
  )
}

ggarrange(plot_trans1,plot_trans2,plot_trans3,ncol = 3,nrow = 1,common.legend = T)
ggarrange(plot_trans4,plot_trans5,plot_trans6,ncol = 3,nrow = 1,common.legend = T)
ggarrange(plot_trans7,plot_trans8,plot_trans9,ncol = 3,nrow = 1,common.legend = T)
ggarrange(plot_trans10,plot_trans11,plot_trans12,ncol = 3,nrow = 1,common.legend = T)
plot_trans13


# choose models with smallest BIC and plot estimated cumulative hazard with CI against non-parametric cumulative hazard
final_model<-model_fitness %>% 
  group_by(trans) %>% 
  slice_min(BIC) %>% 
  ungroup


for (i in 1:13){
  k<-(final_model %>% filter(trans==i))$knots
  state<-(transmat %>% filter(trans==i) %>% 
            mutate(state=factor(state,levels=c(1:6),
                                labels=c("mild AD","moderate AD","severe AD",
                                         "institutionalized mild AD",
                                         "institutionalized moderate AD",
                                         "institutionalized severe AD"))))$state
  next_state<-(transmat %>% filter(trans==i) %>% 
                 mutate(poss_next_state=factor(poss_next_state,levels=c(2:7),
                                     labels=c("moderate AD","severe AD",
                                              "institutionalized mild AD",
                                              "institutionalized moderate AD",
                                              "institutionalized severe AD","death"))))$poss_next_state
  
  assign(paste0("cumhaz_ci_trans",i),
         ggplot() +
           geom_line(data = nonpara_cumhaz_combined %>% 
                       filter(trans==i,
                              time<=10), 
                     aes(x = time, y = cumhaz), 
                     color = "black") +
           geom_line(data = get(paste0("RP_cumhaz_trans",i,"_with_",k,"knots")), 
                     aes(x = .time, y = .pred_cumhaz),color="#56B4E9") +
           geom_ribbon(data = get(paste0("RP_cumhaz_trans",i,"_with_",k,"knots")), 
                       aes(x = .time, 
                           ymin = .pred_lower, 
                           ymax = .pred_upper), 
                       alpha = .1,fill="#56B4E9") +
           scale_x_continuous(breaks = seq(0,10,2),expand=c(0,0))+
           theme_minimal() +
           theme(plot.title = element_text(size=10)) +
           labs(title=paste("From",state,"to",next_state),
                x = 'Time (years)',
                y = 'Cumulative hazard',
                linetype = 'Distribution',
                color = "Transition"))
}

ggarrange(cumhaz_ci_trans1,cumhaz_ci_trans2,cumhaz_ci_trans3,
          cumhaz_ci_trans4,cumhaz_ci_trans5,cumhaz_ci_trans6,
          cumhaz_ci_trans7,cumhaz_ci_trans8,cumhaz_ci_trans9,
          cumhaz_ci_trans10,cumhaz_ci_trans11,cumhaz_ci_trans12,cumhaz_ci_trans13,
          ncol = 3,nrow = 5,common.legend = T)

# transition 2, 8 does not fit very well

# 3. Deal with models that did not fit well ----
# for transition 2, choose 5 knots
cumhaz_ci_trans2<-ggplot() +
  geom_line(data = nonpara_cumhaz_combined %>% 
              filter(trans==2,
                     time<=10), 
            aes(x = time, y = cumhaz), 
            color = "black") +
  geom_line(data = get(paste0("RP_cumhaz_trans",2,"_with_",5,"knots")), 
            aes(x = .time, y = .pred_cumhaz),color="#56B4E9") +
  geom_ribbon(data = get(paste0("RP_cumhaz_trans",2,"_with_",5,"knots")), 
              aes(x = .time, 
                  ymin = .pred_lower, 
                  ymax = .pred_upper), 
              alpha = .1,fill="#56B4E9") +
  scale_x_continuous(breaks = seq(0,10,2),expand=c(0,0))+
  theme_minimal() +
  theme(plot.title = element_text(size=10)) +
  labs(title="From mild AD to institutionalized mild AD",
       x = 'Time (years)',
       y = 'Cumulative hazard',
       linetype = 'Distribution',
       color = "Transition")

final_model<-final_model %>% 
  mutate(knots=ifelse(trans==2,5,knots))



# for transition 8, choose 2 knots
cumhaz_ci_trans8<-ggplot() +
  geom_line(data = nonpara_cumhaz_combined %>% 
              filter(trans==8,
                     time<=10), 
            aes(x = time, y = cumhaz), 
            color = "black") +
  geom_line(data = get(paste0("RP_cumhaz_trans",8,"_with_",2,"knots")), 
            aes(x = .time, y = .pred_cumhaz),color="#56B4E9") +
  geom_ribbon(data = get(paste0("RP_cumhaz_trans",8,"_with_",2,"knots")), 
              aes(x = .time, 
                  ymin = .pred_lower, 
                  ymax = .pred_upper), 
              alpha = .1,fill="#56B4E9") +
  scale_x_continuous(breaks = seq(0,10,2),expand=c(0,0))+
  theme_minimal() +
  theme(plot.title = element_text(size=10)) +
  labs(title="From severe AD to death",
       x = 'Time (years)',
       y = 'Cumulative hazard',
       linetype = 'Distribution',
       color = "Transition")

final_model<-final_model %>% 
  mutate(knots=ifelse(trans==8,2,knots))





# generate cumulative hazard plots for the models
ggarrange(cumhaz_ci_trans1,cumhaz_ci_trans2,cumhaz_ci_trans3,
          cumhaz_ci_trans4,cumhaz_ci_trans5,cumhaz_ci_trans6,
          cumhaz_ci_trans7,cumhaz_ci_trans8,cumhaz_ci_trans9,
          cumhaz_ci_trans10,cumhaz_ci_trans11,cumhaz_ci_trans12,cumhaz_ci_trans13,
          ncol = 3,nrow = 5,common.legend = T)


# Add age and sex as predictors
# organize final models
for (i in c(1:8,10:13)){
  newdata<-mstatedata %>% 
    filter(trans==i)
  
  k<-(final_model %>% filter(trans==i))$knots
  
  assign(paste0("final_model_trans",i),
         flexsurvspline(formula = Surv(Tstart,Tstop,status) ~ AGE+SEX, data = newdata, 
                        k = k, scale = "hazard",
                        inits=c(coef(get(paste0("RP_model_trans",i,"_with_",k,"knots"))), 0)))
  print(i)
  print(k)
}

final_model_trans9<-flexsurvspline(Surv(Tstart,Tstop,status) ~ AGE+SEX,
                                   data=mstatedata %>% 
                                     filter(trans==9), knots = c(log(4),log(6),log(8)),scale="hazard",
                                   inits=c(coef(RP_model_trans9_with_2knots), 0))





# 4. Model for MCI progression to AD and death without AD ----
# 4.1. MCI progression to AD ----
nonpara_model_mci_ad<-coxph(Surv(fu_imp,prog) ~ 1, data=nacc_mci_select, method="breslow")

# estimate cumulative hazard
nonpara_model_cumhaz_mci_ad<-basehaz(nonpara_model_mci_ad)

for (k in 1:6){
  
 model <- flexsurvspline(Surv(fu_imp,prog) ~ 1,
                          data=nacc_mci_select, k=k,scale="hazard")
  
  # estimate cumulative hazard
  cumhaz <- (predict(model, type = "cumhaz", newdata = data.frame(AGE = NA),
                     times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
    rename(.time=.eval_time)
  
  # return one model for each transition
  assign(paste0("RP_model_mci_ad_with_",k,"knots"),model)
  assign(paste0("RP_cumhaz_mci_ad_with_",k,"knots"),cumhaz)
}


rm(newdata,model,cumhaz)




(do.call(rbind, lapply(1:6, function(k) {
  cumhaz_df <- get(paste0("RP_model_mci_ad_with_", k, "knots")) %>% glance()
  
  cumhaz_df$model <- paste0("RP model with ", k, " knots")
  return(cumhaz_df)
})) %>% slice_min(BIC))$model


mci_ad_plot<-ggplot() +
  geom_line(data = nonpara_model_cumhaz_mci_ad %>% 
              filter(time<=10), 
            aes(x = time, y = hazard), 
            color = "black") +
  geom_line(data = RP_cumhaz_mci_ad_with_1knots, 
            aes(x = .time, y = .pred_cumhaz),color="#56B4E9") +
  geom_ribbon(data = RP_cumhaz_mci_ad_with_1knots, 
              aes(x = .time, 
                  ymin = .pred_lower, 
                  ymax = .pred_upper), 
              alpha = .1,fill="#56B4E9") +
  scale_x_continuous(breaks = seq(0,10,2),expand=c(0,0))+
  theme_minimal() +
  theme(plot.title = element_text(size=10))+
  labs(title = "From MCI to mild AD",
       x = 'Time (years)',
       y = 'Cumulative hazard',
       linetype = 'Distribution',
       color = "Transition")



mci_ad_model<-flexsurvspline(Surv(fu_imp,prog) ~ AGE+factor(SEX)+factor(NACCNE4S),
                             data=nacc_mci_select,k=1,scale="hazard",
                             inits=c(coef(RP_model_mci_ad_with_1knots), 0))

mci_ad_model

mci_adsurv_data<-predict(mci_ad_model, type = "survival", 
                     newdata = data.frame(AGE = 77, SEX = "FEMALE", NACCNE4S = 0),
                     times = seq(0, 16, by = 0.5))[[1]][[1]] %>%
  rename(.time = .eval_time)

(ggsurvplot(
  survfit(Surv(fu_imp,prog) ~ 1, data = nacc_mci_select),
  break.time.by = 1, break.y.by=0.05, 
  ggtheme = theme_minimal() + theme(panel.grid.major = element_line(color = "gray80"))
))$plot+ geom_line(data = mci_adsurv_data, aes(y = .pred_survival, x = .time), color = "black")


# 4.2. MCI progression to death without AD ----
nonpara_model_mci_death<-coxph(Surv(fu_imp,death_new) ~ 1, data=nacc_mci_select, method="breslow")

# estimate cumulative hazard
nonpara_model_cumhaz_mci_death<-basehaz(nonpara_model_mci_death)


for (k in 1:5){# 6 knots do not work
  
  model <- flexsurvspline(Surv(fu_imp,death_new) ~ 1,
                          data=nacc_mci_select, k=k,scale="hazard")
  
  # estimate cumulative hazard
  cumhaz <- (predict(model, type = "cumhaz", newdata = data.frame(AGE = NA),
                     times = seq(0, 10, by = 0.1),conf.int=T))[[1]] %>% as.data.frame() %>% 
    rename(.time=.eval_time)
  
  # return one model for each transition
  assign(paste0("RP_model_mci_death_with_",k,"knots"),model)
  assign(paste0("RP_cumhaz_mci_death_with_",k,"knots"),cumhaz)
}


rm(newdata,model,cumhaz)




(do.call(rbind, lapply(1:5, function(k) {
  cumhaz_df <- get(paste0("RP_model_mci_death_with_", k, "knots")) %>% glance()
  
  cumhaz_df$model <- paste0("RP model with ", k, " knots")
  return(cumhaz_df)
})) %>% slice_min(BIC))$model



# 2 knots fit better
mci_death_plot<-ggplot() +
  geom_line(data = nonpara_model_cumhaz_mci_death %>% 
              filter(time<=10), 
            aes(x = time, y = hazard), 
            color = "black") +
  geom_line(data = RP_cumhaz_mci_death_with_2knots, 
            aes(x = .time, y = .pred_cumhaz),color="#56B4E9") +
  geom_ribbon(data = RP_cumhaz_mci_death_with_2knots, 
              aes(x = .time, 
                  ymin = .pred_lower, 
                  ymax = .pred_upper), 
              alpha = .1,fill="#56B4E9") +
  scale_x_continuous(breaks = seq(0,10,2),expand=c(0,0))+
  theme_minimal() +
  theme(plot.title = element_text(size=10))+
  labs(title = "From MCI to death",
       x = 'Time (years)',
       y = 'Cumulative hazard',
       linetype = 'Distribution',
       color = "Transition")



mci_death_model<-flexsurvspline(Surv(fu_imp,death_new) ~ AGE+factor(SEX)+factor(NACCNE4S),
                                data=nacc_mci_select, 
                                k = 2,scale="hazard")

mci_death_model

# apoe does not predict death

mci_death_model<-flexsurvspline(Surv(fu_imp,death_new) ~ AGE+factor(SEX),
                                data=nacc_mci_select, 
                                k = 2,scale="hazard")

mci_death_model


pdf(file="model fit.pdf",width = 12,height = 14)
ggarrange(mci_ad_plot,mci_death_plot,
          cumhaz_ci_trans1,cumhaz_ci_trans2,cumhaz_ci_trans3,
          cumhaz_ci_trans4,cumhaz_ci_trans5,cumhaz_ci_trans6,
          cumhaz_ci_trans7,cumhaz_ci_trans8,cumhaz_ci_trans9,
          cumhaz_ci_trans10,cumhaz_ci_trans11,cumhaz_ci_trans12,
          cumhaz_ci_trans13,
          ncol = 3,nrow = 5,common.legend = T)
dev.off()


save(list = c(paste0('final_model_trans', 1:13),"mci_ad_model","mci_death_model"), 
     file = 'models_xx.RData')
