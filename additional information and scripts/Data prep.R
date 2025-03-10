# LOAD LIBRARIES ----

library(RMySQL)
library(tidyverse)
library(consort)
library(ggplot2)
library(flextable)
library(labelled)
library(lme4)
library(janitor)
library(officer)
library(survival)
library(mstate)
library(data.table)
library(survminer)

#Conversion tables----

MOCA_to_MMSE<-data.frame(MOCA_VARDE=seq(0,30,1), MMSE_CONVERTED=c(5,7,9,10,11,12,13,14,15,16,17,18,19,20,20,21,22,23,24,25,25,26,27,27,28,28,29,29,29,30,30))

diag_convert<-data.frame(DIAGNOS=c(
  'DEMENS_UNS',
  'VASKULAR_DEMENS_INKL_SUBKORTIKAL_VASKULAR_DEMENS',
  'BLANDDEMENS_VID_ALZHEIMERS_SJUKDOM_VASKULAR_DEMENS',
  'DEMENS_VID_ALZHEIMERS_SJUKDOM_SEN_DEBUT' ,
  'DEMENS_VID_ALZHEIMERS_SJUKDOM_TIDIG_DEBUT',
  'DEMENS_VID_PARKINSONS_SJUKDOM',
  'FRONTOTEMPORAL_DEMENS' ,
  'LEWY_BODY_DEMENS' ,
  'LINDRIG_KOGNITIV_STORNING',
  'OVRIG_DEMENS' 
), diagnosis=c(
  'Unspecified',
  'VaD',
  'Mixed AD/VaD',
  'AD',
  'Early onset AD',
  'PD',
  'FTD',
  'LBD',
  'MCI',
  'Other')
)

# LOAD DATA AND PREPARE EVENT DATASET-----

# 
con <- dbConnect(MySQL(),
                 user='linus.jonsson',
                 password=rstudioapi::askForPassword(prompt = "Password"),
                 host='h1cogbase01.nvs.ki.se',
                 port=3306,
                 dbname='priv_linjon')


#Get tables from server----

tables<-c('INCLUSION', 'COG', 'DEMO', 'MORT', 'INST')
SVEDEM<-lapply(tables, function(x) dbReadTable(con, x))

SVEDEM$STATECOSTWIDE<-dbReadTable(con, 'SVEDEM_2021.STATECOSTWIDE')
names(SVEDEM)<-c(tables, 'STATECOSTWIDE')

incl<-SVEDEM$INCLUSION %>%  select(-n) %>%
  left_join(diag_convert) %>%
  mutate(excl_diag=ifelse(!DIAGNOS %in% c('DEMENS_VID_ALZHEIMERS_SJUKDOM_TIDIG_DEBUT', 'DEMENS_VID_ALZHEIMERS_SJUKDOM_SEN_DEBUT', 'LINDRIG_KOGNITIV_STORNING'
                                          #, 'BLANDDEMENS_VID_ALZHEIMERS_SJUKDOM_VASKULAR_DEMENS', 'DEMENS_UNS'
  ), diagnosis, NA),
  excl_date=ifelse(!DIAGNOSDATUM>='2014-01-01', year(DIAGNOSDATUM), NA),
  excl_age=factor(ifelse(!year(DIAGNOSDATUM)-d_YOB>=60, cut(year(DIAGNOSDATUM)-d_YOB, breaks=c(0,40,50,60)), NA), labels=c('<40', '40-49', '50-59')),
  excl_biomarker=case_when(amyloid==0 ~ 'Amyloid negative'),
  biomarker_status=case_when(amyloid==1 ~ 'Amyloid positive',
                             amyloid==0 | is.na(amyloid) ~ 'Amyloid status missing'),
  
  included=case_when(is.na(excl_date) & is.na(excl_age) & is.na(excl_diag) & amyloid==1 ~ 'Biomarker',
                     is.na(excl_date) & is.na(excl_age) & is.na(excl_diag) ~ 'Clinical',
                     T ~ 'Excluded')
  )

#Disease progression----
cog<-SVEDEM$COG %>%
  
  mutate(DATE=ifelse(substr(DATE,1,2)=='19', paste0('20',substr(DATE,3,10)), DATE),
         DATE=as.Date(DATE),
         MMSESR_VARDE=as.numeric(MMSESR_VARDE),
         untestable=is.na(MMSESR_VARDE) & is.na(MOCA_VARDE) & (MMSESR=='EJ_TESTBAR' | MOCA=='EJ_TESTBAR'),
  ) %>%
  left_join(MOCA_to_MMSE) %>%
  mutate(value=case_when(!is.na(MMSESR_VARDE) ~ MMSESR_VARDE,
                         !is.na(MOCA_VARDE) ~ MMSE_CONVERTED)) %>%
  select(LOPNR, DATE, value, untestable) %>% 
  
  #Impute untestable as Zero if the patient has previoustly tested <10 and no subsequent values
  group_by(LOPNR) %>%
  arrange(LOPNR, -as.numeric(DATE)) %>%
  mutate(nonmissing=!is.na(value),
         cumul_nonmissing=cumsum(nonmissing),
         below10=min(value, na.rm=T)<10,
         value=ifelse(untestable & below10 & cumul_nonmissing==0, 0, value)
  ) %>%
  drop_na %>%
  group_by(LOPNR, DATE) %>%
  summarise(value=max(value, na.rm=T)) %>%
  inner_join(incl %>% filter(included!='Excluded')) %>%
  group_by(LOPNR) %>%
  arrange(LOPNR,DATE) %>%
  mutate(time=as.numeric(as.Date(DATE)-as.Date(DIAGNOSDATUM))/365,
         basevalue=first(value),
         baseline_stage=factor(cut(basevalue, breaks=c(-Inf,10,20,25,30)), levels=c('(25,30]','(20,25]','(10,20]','(-Inf,10]'),  labels=c('MCI','Mild','Moderate','Severe')),
         stage=factor(cut(value, breaks=c(-Inf,10,20,25,30)), levels=c('(25,30]','(20,25]','(10,20]','(-Inf,10]'),  labels=c('MCI','Mild','Moderate','Severe')),
         n_values=n()
  ) %>%
  #  filter(n_values>1) %>%
  filter(time>=0)  

#Institutionalization and mortality - multistate model----

inst<-SVEDEM$INST %>%
  mutate(DATE=as.Date(paste0(PERIOD,'01'), format='%Y%m%d')) %>%
  group_by(LOPNR) %>% summarize(DATE=min(DATE)) %>% ungroup %>%
  inner_join(incl %>% select(LOPNR, DIAGNOSDATUM)) %>%
  mutate(DATE=pmin(DATE, DIAGNOSDATUM),
    insttime=as.numeric(DATE-as.Date(DIAGNOSDATUM)))

mort<-SVEDEM$MORT %>%
  mutate(DODSDAT_edit=ifelse(substr(DODSDAT,5,8)=='0000', paste0(substr(DODSDAT,1,4),'0101'), DODSDAT),
         DODSDAT_edit=ifelse(substr(DODSDAT_edit,7,8)=='00', paste0(substr(DODSDAT,1,6),'01'), DODSDAT_edit),
         DATE=as.Date(DODSDAT_edit, format='%Y%m%d'))

cens<-cog %>% 
  group_by(LOPNR) %>% 
  arrange(DATE) %>%
  slice_tail(n=1) %>%
  mutate(censtime=as.numeric(DATE-as.Date(DIAGNOSDATUM))+365*2) %>%
  select(LOPNR, censtime)

mmse<-cog %>% 
  mutate(mmsetime=as.numeric(DATE-as.Date(DIAGNOSDATUM))) %>% 
  select(LOPNR, DATE, MMSE=value, mmsetime)

trans<-cog %>% 
  select(LOPNR, DATE, value) %>%
  mutate(state=case_when(value<10 ~ 'Severe',
                         value<=20 ~ 'Moderate',
                         value<26 ~ 'Mild',
                         value>=26 ~ 'MCI')) %>%
  group_by(LOPNR) %>%
  arrange(LOPNR, DATE) %>%
  filter(is.na(lag(value)) | lag(value) > value) %>%
  filter(is.na(lag(state)) | state!=lag(state))

#For patients making transitions that 'jump' one or more states, add intermediate transition at half time
trans<-rbind(
  trans,
trans %>% group_by(LOPNR) %>%
  mutate(next_state=lead(state),
         next_date=lead(DATE),
         next_value=lead(value)) %>%
  filter(state=='MCI' & next_state=='Moderate') %>%
  mutate(DATE=DATE+(next_date-DATE)/2,
         value=round((value+next_value)/2),
         state='Mild'
         )%>%
  select(LOPNR, DATE, value, state),

trans %>% group_by(LOPNR) %>%
  mutate(next_state=lead(state),
         next_date=lead(DATE),
         next_value=lead(value)) %>%
  filter(state=='Mild' & next_state=='Severe') %>%
  mutate(state='Moderate',
         DATE=DATE+(next_date-DATE)/2,
         value=round((value+next_value)/2)) %>%
  select(LOPNR, DATE, value, state),

trans %>% group_by(LOPNR) %>%
  mutate(next_state=lead(state),
         next_date=lead(DATE),
         next_value=lead(value)) %>%
  filter(state=='MCI' & next_state=='Severe') %>%
  mutate(DATE=DATE+(next_date-DATE)/2,
         value=round((value+next_value)/2),
         state=case_when(value<10 ~ 'Severe',
                         value<=20 ~ 'Moderate',
                         value<26 ~ 'Mild',
                         value>=26 ~ 'MCI')) %>%
  select(LOPNR, DATE, value, state)
)

demo <- SVEDEM$DEMO  %>% 
  mutate(AGE=year(DIAGNOSDATUM)-d_YOB)%>% 
  select(LOPNR, AGE, SEX) %>%
  inner_join(incl)

eventdata<-incl %>% 
  filter(included !='Excluded') %>% select(LOPNR, DIAGNOSDATUM) %>%
  left_join(cens) %>%
  left_join(mort %>% select(LOPNR, mortdate=DATE)) %>%
  left_join(inst %>% select(LOPNR, insttime, instdate=DATE)) %>%
  left_join(trans %>% filter(state=='MCI') %>% group_by(LOPNR) %>% summarise(MCI_DATE=min(DATE)) ) %>%
  left_join(trans %>%  filter(state=='Mild')   %>% group_by(LOPNR) %>% summarise(Mild_DATE=min(DATE))) %>%
  left_join(trans %>% filter(state=='Moderate') %>% group_by(LOPNR) %>% summarise(Moderate_DATE=min(DATE))) %>%
  left_join(trans %>% filter(state=='Severe') %>% group_by(LOPNR) %>% summarise(Severe_DATE=min(DATE))) %>%
  left_join(trans %>% arrange(LOPNR, DATE) %>% group_by(LOPNR) %>% slice_head %>% select(LOPNR, START=state)) %>%
  mutate(START=as.numeric(factor(START, levels=c('MCI', 'Mild', 'Moderate', 'Severe'))),
         DIAGDATE=as.Date(DIAGNOSDATUM),
         morttime=ifelse(!is.na(mortdate), as.numeric(mortdate)-as.numeric(DIAGDATE), censtime),
         mort=as.numeric(!is.na(mortdate) & morttime<=censtime),
         insttime=ifelse(!is.na(instdate) & insttime<censtime, insttime, censtime),   
         insttime=ifelse(insttime==morttime, insttime-1, insttime), #deal with ties on institutionalization and mortality
         MCI_time=ifelse(!is.na(MCI_DATE), as.numeric(MCI_DATE)-as.numeric(DIAGDATE), censtime),
         Mild_time=ifelse(!is.na(Mild_DATE), as.numeric(Mild_DATE)-as.numeric(DIAGDATE), censtime),
         Moderate_time=ifelse(!is.na(Moderate_DATE), as.numeric(Moderate_DATE)-as.numeric(DIAGDATE), censtime),
         Severe_time=ifelse(!is.na(Severe_DATE), as.numeric(Severe_DATE)-as.numeric(DIAGDATE), censtime),
         MCI=as.numeric(!is.na(MCI_DATE)),
         Mild=as.numeric(!is.na(Mild_DATE)),
         Moderate=as.numeric(!is.na(Moderate_DATE)),
         Severe=as.numeric(!is.na(Severe_DATE)),
         MCI_inst_time=ifelse(!is.na(MCI_DATE) & !is.na(instdate) & insttime<=censtime, pmax(insttime, MCI_time), censtime),
         Mild_inst_time=ifelse(!is.na(Mild_DATE) & !is.na(instdate) & insttime<=censtime, pmax(insttime, Mild_time), censtime),
         Moderate_inst_time=ifelse(!is.na(Moderate_DATE) & !is.na(instdate) & insttime<=censtime, pmax(insttime, Moderate_time), censtime),
         Severe_inst_time=ifelse(!is.na(Severe_DATE) & !is.na(instdate) & insttime<=censtime, pmax(insttime, Severe_time), censtime),
         MCI_inst=as.numeric(!is.na(MCI_DATE) & !is.na(instdate) & insttime<=censtime),
         Mild_inst=as.numeric(!is.na(Mild_DATE) & !is.na(instdate) & insttime<=censtime),
         Moderate_inst=as.numeric(!is.na(Moderate_DATE) & !is.na(instdate) & insttime<=censtime),
         Severe_inst=as.numeric(!is.na(Severe_DATE) & !is.na(instdate) & insttime<=censtime)
  ) %>%
  left_join(demo)



msdata<-msprep(time=c('MCI_time', 'Mild_time', 'Moderate_time', 'Severe_time', 'MCI_inst_time','Mild_inst_time', 'Moderate_inst_time', 'Severe_inst_time', 'morttime'), 
       status=c('MCI', 'Mild', 'Moderate', 'Severe', 'MCI_inst', 'Mild_inst', 'Moderate_inst', 'Severe_inst', 'mort'), 
       data=eventdata, 
       id='LOPNR', 
        start=list(state=eventdata$START, 
                   time=rep(0, length(eventdata$START))),
       trans=transMat(x=list(c(2,5,9), 
                           c(3,6,9),
                           c(4,7,9),
                           c(8,9),
                           c(6,9),
                           c(7,9),
                           c(8,9),
                           c(9),
                           c()),
                      
                      names=c('MCI', 'Mild', 'Moderate', 'Severe', 'MCI_inst', 'Mild_inst', 'Moderate_inst', 'Severe_inst', 'mort')),
       #keep=c('AGE', 'd_YOB')
       ) 
msdata <- msdata %>% 
  inner_join(eventdata %>% 
               mutate(LOPNR=factor(LOPNR), 
                      AGEGRP=cut(AGE, breaks=c(0, 65, 75, 85, 120), labels=c('A', 'B', 'C', 'D' ))) %>%
  select(LOPNR, AGE, SEX) )

#OBS! check if inst < diagnosis is handled correctly
msdata <-msdata %>% filter(time>0 & time!=Inf) 

msdata<-expand.covs(msdata,covs=c("AGE","SEX"))

save(list=c('msdata', 'trans'), file='msdata.RData')  

library(tidyverse)
library(mstate)

load('msdata.RData')

runcox<-function(df) {
  cx <- coxph(Surv(Tstart,Tstop,status)~
                AGE.1+AGE.2+AGE.4+AGE.5+AGE.6+AGE.8+AGE.9+AGE.10+AGE.15+AGE.17+
                SEXMALE.1+SEXMALE.2+SEXMALE.4+SEXMALE.5+SEXMALE.6+SEXMALE.8+SEXMALE.9+SEXMALE.10+SEXMALE.15+SEXMALE.17+
                strata(trans),
              data=df,method="breslow")
  return(cx)
  }



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

cycle_length=3 #in months
cycles=12/cycle_length*10 #10 year time horizon

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



get_tp<-function(transition_number, sex, age, hazards, coefficients) {
  
  #Get hazard for each transition at each time point

  if (paste0('AGE.', transition_number) %in% names(coefficients)) age_coef<-as.numeric(coefficients[names(coefficients)==paste0('AGE.', transition_number)]) else age_coef<-0
  if (paste0('SEXMALE.', transition_number) %in% names(coefficients)) sex_coef<-as.numeric(coefficients[names(coefficients)==paste0('SEXMALE.', transition_number)]) else sex_coef<-0
  h=hazards %>% filter(trans==transition_number)   #Base hazard without covariates
  hr=exp(age_coef*age+sex_coef*sex) #Hazard ratio
  h$p=1-exp(-h$haz*hr) #Probability of transition
  return(h$p) #return vector of transition probabilities by cycle
}

#Draw survival curves for all transitions

ggsurvplot(survfit(Surv(Tstart/365, Tstop/365, status)~trans, 
                   data=msdata %>% 
                     filter(trans %in% c(1, 4, 7)), 
                   id=LOPNR), 
           risk.table=T, censor=F, conf.int = T)

ggsurvplot(survfit(Surv(Tstart/365, Tstop/365, status)~trans, 
                   data=msdata %>% 
                     filter(trans %in% c(12, 14, 16)), 
                   id=LOPNR), 
           risk.table=T, censor=F, conf.int = T)

ggsurvplot(survfit(Surv(Tstart/365, Tstop/365, status)~trans, 
                   data=msdata %>% 
                     filter(trans %in% c(2, 5, 8,10)), 
                   id=LOPNR), 
           risk.table=T, censor=F, conf.int = T)

#Mortality from community states
ggsurvplot(survfit(Surv(Tstart/365, Tstop/365, status)~trans, 
                   data=msdata %>% 
                     filter(trans %in% c(3, 6, 9,11)), 
                   id=LOPNR), 
           risk.table=T, censor=F, conf.int = T)

#Mortality from inst states
ggsurvplot(survfit(Surv(Tstart/365, Tstop/365, status)~trans, 
                   data=msdata %>% 
                     filter(trans %in% c(13, 15, 17, 18)), 
                   id=LOPNR), 
           risk.table=T, censor=F, conf.int = T)

#Get cost data from Sandar's table
costtable<-SVEDEM$STATECOSTWIDE %>% mutate(LOPNR=factor(LopNr))
save(costtable, file='costtable.RData')
