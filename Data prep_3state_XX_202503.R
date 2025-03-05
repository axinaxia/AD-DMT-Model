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




setwd("C:/Users/xinxia/OneDrive - Karolinska Institutet/CSF-registersamkörning/AD DMT HE model/Statistical analyses/New analyses_XX_202412")


# 1. Load data sets and clean the data ----
svedem_sql<-dbConnect(MySQL(),
                      user=rstudioapi::askForPassword(prompt = "Username"),
                      password=rstudioapi::askForPassword(prompt = "Password"),
                      host='h1cogbase01.nvs.ki.se',
                      port=3306,
                      dbname='SVEDEM_2021')
dbListTables(svedem_sql)

# 1.1. Data for dementia diagnosis and demographic factors ----
DEM_DIAG<-dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOS, d_YOB, DIAGNOSDATUM, SEX, LP
                      from SVEDEM_2021.GRUNDREG__2;') %>% 
  as.data.frame() %>% 
  mutate(setting="PV") %>% 
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOS, d_YOB, DIAGNOSDATUM, SEX, LP
                   from SVEDEM_2021.GRUNDREG__3;') %>% 
          as.data.frame() %>% 
          mutate(setting="SV")) %>% 
  mutate(DIAGNOSDATUM=as.Date(DIAGNOSDATUM)) %>% 
  group_by(LOPNR) %>% 
  filter(DIAGNOSDATUM==min(DIAGNOSDATUM)) %>% # keep the earliest diagnosis
  ungroup() 

DEM_DIAG<-DEM_DIAG %>% 
  mutate(AGE=year(DIAGNOSDATUM)-d_YOB)


# remove duplicated diagnostic records
DEM_DIAG<-DEM_DIAG %>% 
  distinct(LOPNR,.keep_all = TRUE) 


# 1.2. Cognition data (MMSE and MOCA converted to MMSE) ----
COG<-dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
                FROM SVEDEM_2021.GRUNDREG__2;') %>% 
  as.data.frame() %>%
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
          FROM SVEDEM_2021.GRUNDREG__3;') %>% 
          as.data.frame()) %>%
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
          FROM SVEDEM_2021.UPPF__2;') %>% 
          as.data.frame()) %>%
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
          FROM SVEDEM_2021.UPPF__3;') %>% 
          as.data.frame()) %>%
  rbind(dbGetQuery(svedem_sql, 
                   'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE
          FROM SVEDEM_2021.UPPF_HE_2;') %>%
          as.data.frame() %>% 
          mutate(MOCA=NA, MOCA_VARDE=NA)) %>%
  rbind(dbGetQuery(svedem_sql, 
                   'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE
          FROM SVEDEM_2021.UPPF_SA_2;') %>%
          as.data.frame() %>% 
          mutate(MOCA=NA, MOCA_VARDE=NA)) %>% 
  filter(DATE<=as.Date('2020-12-31'))



#https://www.mocatest.ch/umrechnung/Fasnacht_Wueest_2022.pdf ???
MOCA_to_MMSE<-data.frame(MOCA_VARDE=seq(0,30,1), 
                         MMSE_CONVERTED=c(5,7,9,10,11,12,13,14,15,16,17,18,19,20,20,
                                          21,22,23,24,25,25,26,27,27,28,28,29,29,29,30,30))



COG<-COG %>%
  # fix record errors in MMSE measurement dates
  mutate(DATE=ifelse(substr(DATE,1,2)=='19', paste0('20',substr(DATE,3,10)), DATE), 
         DATE=as.Date(DATE),
         MMSESR_VARDE=as.numeric(MMSESR_VARDE),
         MOCA_VARDE=as.numeric(MOCA_VARDE),
         untestable=is.na(MMSESR_VARDE) & is.na(MOCA_VARDE) & (MMSESR=='EJ_TESTBAR' | MOCA=='EJ_TESTBAR'),
  ) %>%
  
  # convert MoCA to MMSE for people missing MMSE scores
  left_join(MOCA_to_MMSE) %>%
  mutate(value=case_when(!is.na(MMSESR_VARDE) ~ MMSESR_VARDE,
                         !is.na(MOCA_VARDE) ~ MMSE_CONVERTED)) %>%
  select(LOPNR, DATE, value, untestable) %>%
  mutate(value=ifelse(untestable, 0, value)) %>%
  drop_na %>% 
  group_by(LOPNR, DATE) %>%
  summarise(value=max(value, na.rm=T)) %>%
  ungroup() %>% 
  inner_join(DEM_DIAG) %>%
  arrange(LOPNR,DATE) %>%
  group_by(LOPNR) %>%
  # generate MMSE measurement times, baseline cognitive states, and follow-up cognitive states
  mutate(time=as.numeric(as.Date(DATE)-as.Date(DIAGNOSDATUM)),
         n_values=n()) %>%
  ungroup %>% 
  filter(time>=0) 


# 1.3. Institutionalization data ----
INST<-dbGetQuery(svedem_sql, 'SELECT LOPNR, min(PERIOD) AS PERIOD 
            FROM SVEDEM_2021.SOL
            WHERE BOFORM=2
            GROUP BY LOPNR;') %>%
  as.data.frame() %>% 
  mutate(DATE=as.Date(paste0(PERIOD,'01'), format='%Y%m%d')) %>%
  rbind(dbGetQuery(svedem_sql, "SELECT LOPNR, min(VDAT) AS PERIOD 
            FROM SVEDEM_2021.SOL_20072012
            WHERE SBO=1 AND VDAT IS NOT NULL AND VDAT <> ''
            GROUP BY LOPNR;") %>%
          as.data.frame() %>% 
          mutate(DATE=as.Date(PERIOD, format='%Y%m%d'))) %>%
  group_by(LOPNR) %>% summarize(DATE=min(DATE))

# 1.4. CSF biomarker data ----
# define amyloid positivity according to Ab42/Ab40 ratio and Ab42/ptau ratio
sos_flexlab<-dbGetQuery(svedem_sql, 'SELECT * from SVEDEM_2021.FLEXLAB;')
source('biomarker_cutoff.R')



# 1.5. Mortality data ----
MORT<-dbGetQuery(svedem_sql, 'SELECT LOPNR, DODSDAT
            FROM SVEDEM_2021.DORS;') %>% 
  as.data.frame() %>%
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, DODSDAT
            FROM SVEDEM_2021.DORS_AVI;')) %>%
  mutate(DODSDAT_edit=ifelse(substr(DODSDAT,5,8)=='0000', paste0(substr(DODSDAT,1,4),'0101'), DODSDAT),
         DODSDAT_edit=ifelse(substr(DODSDAT_edit,7,8)=='00', paste0(substr(DODSDAT,1,6),'01'), DODSDAT_edit),
         DATE=as.Date(DODSDAT_edit, format='%Y%m%d'))


# 1.6. NPR, drug register, and SOL ----
# restric to before 2020-12-31 to load less data
NPR<-dbGetQuery(svedem_sql, "SELECT LopNr,INDATUM,MVO,DRG, DIAGNOS,ALDER, KON
                FROM SVEDEM_2021.PAR_SV
                WHERE INDATUM <= '2020-12-31';") %>% 
  as.data.frame() %>% 
  mutate(vtyp='sv') %>% 
  rbind(dbGetQuery(svedem_sql, "SELECT LopNr,INDATUM,MVO,DRG, DIAGNOS,ALDER, KON
                FROM SVEDEM_2021.PAR_OV
                WHERE INDATUM <= '2020-12-31';") %>% 
          as.data.frame() %>% 
          mutate(vtyp='ov'))

LMED<-dbGetQuery(svedem_sql, "SELECT LopNr,EDATUM,TKOST
                FROM SVEDEM_2021.LMED
                WHERE EDATUM >= '2014-01-01' AND EDATUM <= '2020-12-31' AND TKOST>0;")

SOL<-dbGetQuery(svedem_sql, "SELECT LopNr,PERIOD,BOFORM, HTJ, HTJTIM, MATD, 
                TRYGG, DAGV, BOSTO,KORTTID, KORTMAN
                FROM SVEDEM_2021.SOL
                WHERE PERIOD >= '2014';")


# 2. Define the study population based on inclusion criteria ----
# ************************************************************************************************
# Inclusion criteria (to be confirmed):
# 1. Diagnosed with early-onset AD, late-onset AD,or MCI
# 2. Age: 50-90 years when receiving a diagnosis
# 3. Abnormal CSF Ab42/Ab40 ratio and Ab42/ptau ratio
# ************************************************************************************************

incl<-DEM_DIAG %>%  
  left_join(bm.path,by=c('LOPNR'='LopNr')) %>%
  mutate(age_trial=as.numeric(AGE>=50&AGE<=90),
         clin_ad=as.numeric(DIAGNOS %in% c('DEMENS_VID_ALZHEIMERS_SJUKDOM_TIDIG_DEBUT',
                                           'DEMENS_VID_ALZHEIMERS_SJUKDOM_SEN_DEBUT')),
         clin_mci=as.numeric(DIAGNOS=='LINDRIG_KOGNITIV_STORNING'),
         no_csf=as.numeric(is.na(Amy42.40)&is.na(Amy42.pTau)),
         neg_csf=as.numeric((Amy42.40==0&Amy42.pTau==0)|(is.na(Amy42.40)&Amy42.pTau==0)|
                              (is.na(Amy42.pTau)&Amy42.40==0)),
         pos_csf=as.numeric(Amy42.40==1|Amy42.pTau==1))
         

# create a data set indicating whether the patient has stroke, bleeding disorders, or seizures
comorb<-NPR %>% 
  rename(LOPNR=LopNr) %>% 
  right_join(DEM_DIAG %>% dplyr::select(LOPNR,DIAGNOSDATUM)) %>% 
  filter(INDATUM<=DIAGNOSDATUM) %>% 
  mutate(comorb=ifelse(grepl("G46|I60|I61|I62|I63|I64|I67|I69|D65|D66|D67|D68|D69|G40|G41|R568",DIAGNOS),1,0)) %>% 
  group_by(LOPNR) %>% 
  slice_max(comorb) %>% 
  ungroup %>% 
  distinct(LOPNR,comorb)

table(comorb$comorb)



# 3. Compile and organize data sets for event outcomes ----
# 3.1. Prepare data for state transitions for multistate models ----
# ************************************************************************************************
# Data processing steps:
# 1. censoring rules:
# for people who do not enter severe dementia, censor at death or 2020-12-31 or 
# do not have MMSE measurements for a consecutive of 1 years
# for people who entered the state of severe dementia, censoring at death or 2020-12-31
# 2. Create cognitive states: Mild, Moderate, Severe
# 3. Address reverse transitions (e.g., Severe to Mild) - no reverse transition allowed
# 4. Add states for institutionalization
# 5. If the date of death equals the date of institutionalization, set the date to date of death,
# ignore institutionalization
# 6. Modify the next state for few people (<5) who transit from mild to institutionalized moderate
# to moderate dementia
# ************************************************************************************************
TRANS<-COG %>% 
  dplyr::select(LOPNR, DATE, value) %>% 
  left_join(MORT %>% dplyr::select(LOPNR,death_date=DATE)) %>% 
  # remove observations after 2020-12-31
  filter(DATE<=as.Date('2020-12-31')) %>%
  # remove observations after deaths 
  filter(!(!is.na(death_date)&DATE>=death_date)) %>% 
  dplyr::select(-death_date) %>%
  # categorize cognitive states
  mutate(state=case_when(value<10 ~ 3,
                         value<20 ~ 2,
                         value>=20 ~ 1)) %>% 
  # address reverse transitions, set the next state as the current state if the next state is lower
  group_by(LOPNR) %>%
  arrange(LOPNR, DATE) %>%
  mutate(state=cummax(state), 
         last_mmse_date=as.Date(max(DATE)),
         last_state=max(state)) %>%
  ungroup %>% 
  dplyr::select(-contains("value"))


# keep the earliest date of each state
TRANS<-TRANS %>% 
  group_by(LOPNR,state) %>%
  filter(DATE == min(DATE) | DATE == max(DATE)) %>%
  ungroup() %>% 
  group_by(LOPNR) %>% 
  mutate(next_state=lead(state),
         next_date=lead(DATE)) %>% 
  ungroup() %>% 
  filter(!is.na(next_state))

table(TRANS$state,TRANS$next_state)


TRANS<-TRANS %>% 
  group_by(LOPNR,state) %>%
  summarise(state=min(state),
            DATE=min(DATE),
            next_state=max(next_state),
            next_date=max(next_date)) %>% 
  ungroup %>% 
  group_by(LOPNR) %>% 
  mutate(last_state=max(next_state),
         last_state_date=max(next_date)) %>% 
  ungroup

table(TRANS$state,TRANS$next_state)


# combine death state for people who are not censored
CENS_STATE<-TRANS %>% 
  group_by(LOPNR) %>% 
  filter(next_state==max(next_state),
         next_date==max(next_date)) %>% 
  ungroup() %>%
  left_join(MORT %>% dplyr::select(LOPNR,death_date=DATE) %>% 
              filter(death_date<=as.Date('2020-12-31'))) %>% 
  mutate(censor_date=as.Date(pmin(death_date, as.Date('2020-12-31'),na.rm = T)),
         mmse_stop=as.numeric(censor_date-last_state_date>365.25&last_state!=3)) %>% 
  filter(mmse_stop==0|last_state==3,
         censor_date==death_date) %>% 
  mutate(state=last_state,
         DATE=last_state_date,
         next_state=4,
         next_date=as.Date(censor_date)) %>% 
  rbind(TRANS %>% 
          group_by(LOPNR) %>% 
          filter(next_state==max(next_state),
                 next_date==max(next_date)) %>% 
          ungroup() %>%
          left_join(MORT %>% dplyr::select(LOPNR,death_date=DATE) %>% 
                      filter(death_date<=as.Date('2020-12-31'))) %>% 
          mutate(censor_date=as.Date(pmin(death_date, as.Date('2020-12-31'),na.rm = T)),
                 mmse_stop=as.numeric(censor_date-last_state_date>365.25&last_state!=3)) %>% 
          filter(last_state==3,
                 is.na(death_date)) %>% 
          mutate(state=last_state,
                 DATE=last_state_date,
                 next_state=3,
                 next_date=as.Date(censor_date)))

TRANS<-TRANS %>% 
  left_join(CENS_STATE %>% dplyr::select(LOPNR,censor_date)) %>%
  dplyr::select(-c(censor_date,last_state,last_state_date)) %>% 
  rbind(CENS_STATE %>% 
          dplyr::select(LOPNR:next_date)) %>%
  arrange(LOPNR,DATE) %>% 
  filter(!is.na(next_state))

table(TRANS$state,TRANS$next_state)


TRANS<-TRANS %>% 
  group_by(LOPNR,state) %>%
  summarise(state=min(state),
            DATE=min(DATE),
            next_state=max(next_state),
            next_date=max(next_date)) %>% 
  ungroup

table(TRANS$state,TRANS$next_state)


# keep the data as one record per distinct state per person
TRANS<-TRANS %>% 
  group_by(LOPNR) %>% 
  mutate(dup_next_state=duplicated(next_state),
         dup_next_state=max(dup_next_state),
         no_trans=as.numeric(state==next_state),
         no_trans=max(no_trans)) %>% 
  ungroup

TRANS<-TRANS %>% 
  filter(dup_next_state==0) %>% 
  rbind(TRANS %>% 
          filter(dup_next_state==1) %>% 
          group_by(LOPNR,next_state) %>% 
          mutate(next_date=max(next_date)) %>% 
          ungroup %>% 
          filter(state!=next_state)) %>% 
  arrange(LOPNR,state) %>% 
  select(-c(dup_next_state,no_trans))

table(TRANS$state,TRANS$next_state)


# combine institutionalization data
TRANS_INST<-TRANS %>% 
  dplyr::select(LOPNR,state,DATE) %>% 
  rbind(TRANS %>% 
          dplyr::select(LOPNR,state=next_state,DATE=next_date)) %>%
  distinct_all %>% 
  arrange(LOPNR,DATE) %>% 
  left_join(INST %>% dplyr::select(LOPNR,inst_date=DATE)) %>% 
  mutate(inst=as.numeric(DATE>inst_date)) %>% 
  group_by(LOPNR) %>% 
  mutate(inst_before_study=ifelse(inst_date<=min(DATE),1,0),
         not_inst_across=as.numeric(max(inst)==0|is.na(inst))) %>%
  ungroup


# modify states: 1=Mild, 2=Moderate, 3=Severe, 4=Institutionalized mild, 
# 5=Institutionalized moderate, 6=Institutionalized severe, 7=Death
TRANS_INST<-TRANS_INST %>% 
  # for people who were already institutionalized before the start of the study
  filter(inst_before_study==1) %>% 
  mutate(state=state+3) %>% 
  # for people who were not institutionalized throughout the study
  rbind(TRANS_INST %>% 
          filter(not_inst_across==1) %>% 
          mutate(state=ifelse(state==4,7,state))) %>% # change the number for death from 4 to 7
  # for people who were not institutionalized before the study but were institutionalized during the study
  rbind(TRANS_INST %>% 
          filter(inst_before_study==0&not_inst_across==0) %>% 
          rbind(TRANS_INST %>% 
                  filter(inst_before_study==0&not_inst_across==0) %>% 
                  mutate(DATE=inst_date) %>% distinct(LOPNR,DATE,.keep_all = T)) %>% 
          arrange(LOPNR,DATE) %>% 
          group_by(LOPNR) %>% 
          mutate(state=cummax(state),
                 # for states after the dates of institutionalization,
                 # change the states to institutionalized states and 
                 # change death from 4 to 7
                 # the states before remain the same
                 state=ifelse(DATE>=inst_date,state+3,state))) %>% 
  group_by(LOPNR) %>% 
  mutate(next_state=lead(state),
         next_date=lead(DATE)) %>% 
  ungroup %>% 
  filter(!is.na(next_state)) %>% 
  arrange(LOPNR,DATE) %>% 
  dplyr::select(-c(inst_date,inst,inst_before_study,not_inst_across)) 


# check transitions
table(TRANS_INST$state,TRANS_INST$next_state)

# Do not allow jumping transitions
TRANS_INST<-TRANS_INST %>% 
  mutate(next_state=ifelse(state==1&(next_state==3|next_state==5|next_state==6),2,next_state),
         next_state=ifelse(state==2&next_state==6,3,next_state),
         next_state=ifelse(state==4&next_state==6,5,next_state))

TRANS_INST<-TRANS_INST %>% 
  filter(!DATE==next_date) %>% 
  group_by(LOPNR,state) %>% 
  summarise(DATE=min(DATE),
            next_state=max(next_state),
            next_date=max(next_date)) %>% 
  ungroup




# 3.2. Compile data for multistate models ----
eventdata<-incl %>% 
  left_join(comorb %>% dplyr::select(LOPNR,comorb)) %>% 
  filter(age_trial==1,
         clin_mci==0,
         pos_csf==1|(clin_ad==1&setting=="SV"&LP=="JA"),
         comorb==0) %>%
  dplyr::select(LOPNR,SEX,AGE,DIAGNOSDATUM) %>%
  inner_join(TRANS_INST) %>% 
  group_by(LOPNR) %>% 
  mutate(init_date=as.Date(min(DATE))) %>% 
  ungroup


# check sample sizes
length(unique(eventdata$LOPNR))
length(unique(TRANS_INST$LOPNR))


# create a table indicating possible transitions
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


# make the data format the same as msdata data (http://cran.nexr.com/web/packages/mstate/vignettes/Tutorial.pdf)
mstatedata<-eventdata %>% 
  inner_join(transmat,
             relationship = "many-to-many")


mstatedata<-mstatedata %>% 
  mutate(Tstart=as.numeric((DATE-init_date)/365.25),
         Tstop=as.numeric((next_date-init_date)/365.25),
         time=Tstop-Tstart,
         status=as.numeric(next_state==poss_next_state)) %>% 
  dplyr::select(id=LOPNR,from=state,to=poss_next_state,trans,Tstart,Tstop,time,status,AGE,SEX) %>% 
  arrange(id,from) %>% 
  #set time-varying age ???
  mutate(AGE=AGE+Tstart) %>% 
  mutate(trans=factor(trans))



# add class and attributes, so the output data set can use functions from mstate
class(mstatedata)<-c('msdata','data.frame')
attributes(mstatedata)$trans<-transMat(x=list(c(2,4,7), 
                                              c(3,5,7),
                                              c(6,7),
                                              c(5,7),
                                              c(6,7),
                                              c(7),
                                              c()),
                                       names=c('Mild', 'Moderate', 'Severe', 
                                               'Mild_inst', 'Moderate_inst', 
                                               'Severe_inst', 'mort'))


events(mstatedata)
table(mstatedata$from,mstatedata$to)


# check follow-up time per transition
mstatedata %>% 
  group_by(from,to) %>% 
  summarise(fu_mean=mean(time),
            fu_max=max(time),
            fu_min=min(time)) %>% 
  ungroup() %>% 
  arrange(from,to)

# adds type-specific covariates, type-specific covariates can be used to analyse 
# separate effects on all event types in a single analysis based on a stacked data set 
mstatedata<-expand.covs(mstatedata,covs=c("AGE","SEX")) # final data for multistate modeling

ggsurvplot(
  survfit(Surv(Tstart, Tstop,status) ~ 1, data = mstatedata %>% 
            filter(trans==1)),
  break.time.by = 1, break.y.by=0.05, 
  ggtheme = theme_minimal() + theme(panel.grid.major = element_line(color = "gray80"))
)


# 4. Compile data for costs ----
# 4.1. Import CPI and DRG weights and prices ----
cpi <- read_excel("CPI.xlsx") %>% mutate(year=as.numeric(year))

drgcost <- read_excel("Retrospektiva DRG vikter somatisk vard 2012-2021_uppdaterad.xlsx", 
                      n_max=1, col_types=c('skip', 'skip', rep('text',10 )))%>%
  pivot_longer(cols=as.character(2012:2021), names_to='year', values_to='drgcost') %>% 
  mutate(year=as.numeric(year),
         drgcost=as.numeric(drgcost)) 

drgcost <- drgcost %>% left_join(cpi, by='year') %>% 
  mutate(cpi2023= as.numeric('403.7'),
         factor= cpi2023/CPI, # convert drgcost to 2023 value
         drgcost= drgcost*factor) %>% select(year, drgcost)

drg<-read_excel('Retrospektiva DRG vikter somatisk vard 2012-2021_uppdaterad.xlsx', skip=2,  
                col_types=c('text', 'text', rep('numeric',10)))  %>% 
  pivot_longer(cols=as.character(2012:2021),names_to='year', values_to='weight') %>% 
  mutate(year=as.numeric(year)) %>%
  select(DRG, year, weight) %>% 
  drop_na #get rid of missing weights

# Create a table "drgnew" with all DRGs and all years
# Impute missing DRG weights by rolling forward and backward
alldrgyr<-crossing(DRG=unique(drg$DRG), year=2012:2021)

#set as data tables and set keys
drg<-data.table(drg)
alldrgyr<-data.table(alldrgyr)
setkey(drg, DRG, year )
setkey(alldrgyr, DRG, year)

#create new dataset rolling forward
temp<-drg[alldrgyr, roll=T]  %>% drop_na

#use new dataset to roll backward
setkey(temp, DRG, year)
drgnew<-temp[alldrgyr, roll=-Inf]  
rm(temp)



# 4.2. Outpatient and inpatient care costs ----
# Calculate costs using DRG codes and DRG weights*costs/weight
npr_cost<-NPR %>% 
  mutate(INDATUM=as.Date(INDATUM)) %>% 
  filter(INDATUM>=as.Date("2014-01-01")) %>% 
  rename(LOPNR=LopNr) %>% 
  mutate(year=year(as.Date(INDATUM))) %>% 
  left_join(drgnew) %>% 
  left_join(drgcost) %>%
  mutate(cost=weight*drgcost) 


# DRG codes in "drgnew" aren't exhaustive, impute missing costs with mvo costs when costs are missing
# generate costs per MVO (medicinska veksamhetsområde)
mvocost<-NPR %>% 
  rename(LOPNR=LopNr) %>% 
  mutate(year=year(as.Date(INDATUM))) %>% 
  distinct(year,DRG,MVO,vtyp) %>% 
  full_join(drgnew) %>% 
  full_join(drgcost) %>% 
  mutate(mvocost=weight*drgcost) %>%
  filter(!is.na(mvocost)) %>% 
  group_by(vtyp, MVO) %>%
  summarize(mvocost=mean(mvocost, na.rm=T)) %>% 
  filter(!MVO %in% c("-24"," 30"," 31"," 43"," 52","") & !is.na(MVO))


# impute missing costs with mvo costs
npr_cost<-npr_cost %>% 
  left_join(mvocost)
npr_cost[is.na(npr_cost$cost),]$cost<-npr_cost[is.na(npr_cost$cost),]$mvocost

npr_cost<-npr_cost %>%
  mutate(costdate=INDATUM,
         hcr=ifelse(vtyp=='ov','outpatient','inpatient')) %>% 
  dplyr::select(LOPNR, costdate, hcr, cost)

# 4.3. Drug costs ----
drugcost<-LMED %>% 
  rename(LOPNR=LopNr) %>% 
  mutate(EDATUM=as.Date(EDATUM),
         hcr="drug") %>% 
  rename(costdate=EDATUM,
         cost=TKOST)

# convert drugcost to 2023 value
drugcost <- drugcost %>% 
  mutate(year=as.numeric(year(costdate))) %>%
  left_join(cpi, by='year') %>% 
  mutate(cpi2023= as.numeric('403.7'),
         factor= cpi2023/CPI,
         cost= cost*factor) %>% 
  select(LOPNR, costdate, hcr, cost)


# 4.4. Social care costs ----
# ************************************************************************************************
# meaning of variables in SOL 
# (1) BOFORM: Brukarens permanentboende den sista dagen i månaden
#     1=ordinärt boende, 2=särskilt boende, 3=annat boende
# (2) HTJ: Brukaren har den sista dagen i månaden ett verkställt biståndsbeslut om hemtjänst
#     0=nej, 1=ja
# (3) HTJTIM: Antal beviljade/beräknade hemtjänsttimmar per månad för beslut som var verkställt 
#     den sista dagen i månaden
#     0-744 = antal timmar, 777 = omsorg dygnet runt, 888 = ej tillämpligt, 999 = uppgift saknas
# (4) KORTMAN: Antal dygn med korttidsvård/korttidsboende som brukaren erhöll under månaden
# (5) BOSTO: Brukaren har den sista dagen i månaden ett verkställt biståndsbeslut om boendestöd	
#     1=ja, 0=nej, 9=uppgift saknas
# (6) DAGV: Brukaren har den sista dagen i månaden ett verkställt biståndsbeslut om dagverksamhet
#     1=ja, 0=nej, 9=uppgift saknas
# ************************************************************************************************
comcare<-SOL %>% 
  rename(LOPNR=LopNr) %>% 
  mutate(costdate=as.Date(paste0(PERIOD,'01'), format='%Y%m%d'),
         inst=as.numeric(BOFORM==2),
         # set if BOFORM = 2, other support = zero
         hhtim=case_when(BOFORM == 2 ~ 0,
                         BOFORM!=2 & HTJ==1 & HTJTIM<777 ~ HTJTIM,
                         BOFORM!=2 &  HTJ==1 & HTJTIM==777 ~ 24*30,T ~0), 
         
         shortterm_days=as.numeric(KORTMAN),
         shortterm_days=case_when(BOFORM == 2 ~ 0, 
                                  shortterm_days>31|shortterm_days<0 ~ 0,
                                  T ~ as.numeric(shortterm_days)),
         support=case_when(BOFORM == 2 ~ 0,BOSTO=="1"~ 1, T ~ 0),
         dagv=case_when(BOFORM == 2 ~ 0, DAGV=='1'~ 1,T ~ 0)) %>% 
  # remove duplicated observations with same LopNr and costdate
  distinct(LOPNR, costdate,.keep_all = T)

# cannot find hemtjänst/timme kostnad and korttidsboende, use 2021 price*cpi2023/cpi2021
comcare_cost<-comcare %>%
  mutate(cost_homecare=ifelse(BOFORM==2,0,hhtim)*602.76, #people with inst cost does have homehelp cost and other costs
         cost_daytimeactivity=dagv*115000/12, # 2023 price
         cost_shortterm=shortterm_days*3859.29,
         cost_housingsupport=support*49100.05/12) # 2023 price


comcare_cost<-comcare_cost %>% 
  select(LOPNR, costdate, starts_with('cost_')) %>%
  pivot_longer(cols=starts_with('cost_'), names_prefix = 'cost_', values_to='cost') %>% 
  filter(cost >0) %>% 
  rename(hcr=name)


# calculate costs of institutionalization individually
# yearly cost of institutionalization: 1152118 in 2023
comcare_cost<-comcare_cost %>%
  rbind(eventdata %>%
          inner_join(eventdata %>% 
                       filter(DIAGNOSDATUM>=as.Date("2014-01-01")) %>% 
                       distinct(LOPNR)) %>% 
          filter(state>=4&state<=6) %>% 
          mutate(cost=(next_date-DATE)/365.25*1152118,
                 costdate=DATE,
                 hcr='inst') %>% 
          select(LOPNR, costdate, hcr, cost)) 


# 4.5. Compile all costs ----
allcost<-rbind(npr_cost,drugcost,comcare_cost)

costtypes<-allcost %>% select(hcr) %>% distinct %>% filter(!is.na(hcr))

# merging state with costs
statecost<- eventdata %>% 
  inner_join(eventdata %>% 
               filter(DIAGNOSDATUM>=as.Date("2014-01-01")) %>% 
               distinct(LOPNR)) %>% 
  select(LOPNR,state,DATE,next_state,next_date) %>% 
  left_join(allcost, relationship = "many-to-many") %>% 
  filter(costdate>=DATE & costdate<next_date) %>% 
  group_by(LOPNR,state, hcr) %>% 
  summarize(sumcost=sum(cost)) %>% 
  ungroup()

# set costs of other social care for institutionalized states to 0 to avoid cost duplication
statecost<-statecost %>% 
  mutate(sumcost=ifelse((state>=4&state<=6)& 
                          hcr %in% c('homecare', 'shortterm', 'daytimeactivity', 'housingsupport'), 0, sumcost))

# fill in all type of costs for each visit
statecost<-eventdata %>% 
  inner_join(eventdata %>% 
               filter(DIAGNOSDATUM>=as.Date("2014-01-01")) %>% 
               distinct(LOPNR)) %>% 
  select(LOPNR,state,DATE,next_state,next_date) %>% 
  group_by(LOPNR,state) %>% 
  cross_join(costtypes) %>% 
  ungroup() %>% 
  left_join(statecost) %>% 
  mutate(sumcost=ifelse(is.na(sumcost), 0, sumcost)) %>%
  mutate(hcr = factor(hcr, levels = c('drug','inpatient', 'outpatient',
                                      'inst', 'homecare', 'daytimeactivity', 
                                      'shortterm', 'housingsupport'))) 


costwide<-statecost %>% 
  pivot_wider(names_from=hcr, values_from=sumcost) %>% 
  mutate(totalcost=drug+inpatient+outpatient+inst+homecare+daytimeactivity+
           shortterm+housingsupport) %>% 
  select(LOPNR, state:next_date, totalcost) # final cost data for economic modeling



# costs for MCI
NPR_MCI_svedem<-NPR %>% 
  filter(grepl("F067",DIAGNOS),
         MVO %in% c(221, 241, 901, 944, 928),
         ALDER>=50&ALDER<=90) %>% 
  mutate(INDATUM=as.Date(INDATUM),
         LOPNR=LopNr) %>% 
  group_by(LOPNR) %>% 
  summarise(mci_date=min(INDATUM)) %>% 
  ungroup() %>% 
  distinct(LOPNR,.keep_all = T) %>% 
  inner_join(incl %>% 
               filter(clin_mci==0,
                      pos_csf==1|(clin_ad==1&setting=="SV"&LP=="JA")) %>%
               dplyr::select(LOPNR,DIAGNOSDATUM)) %>% 
  filter(mci_date<DIAGNOSDATUM&mci_date>=DIAGNOSDATUM-365.25*5,
         mci_date>=as.Date("2014-01-01")) %>% 
  inner_join(COG %>% 
               filter(DATE==DIAGNOSDATUM,value>=21) %>% 
               distinct(LOPNR))

# remove comorbidities that are exclusion criteria for Lecanemab
NPR_MCI_svedem<-NPR_MCI_svedem %>% 
  left_join(NPR %>% 
              rename(LOPNR=LopNr) %>% 
              right_join(NPR_MCI_svedem %>% dplyr::select(LOPNR,mci_date)) %>% 
              filter(INDATUM<=mci_date) %>% 
              mutate(comorb=ifelse(grepl("G46|I60|I61|I62|I63|I64|I67|I69|D65|D66|D67|D68|D69|G40|G41|R568",DIAGNOS),1,0)) %>% 
              group_by(LOPNR) %>% 
              slice_max(comorb) %>% 
              ungroup %>% 
              distinct(LOPNR,comorb)) %>% 
  filter(comorb==0)


NPR_MCI_svedem_inst<-NPR_MCI_svedem %>% 
  left_join(INST)

NPR_MCI_svedem_inst<-NPR_MCI_svedem_inst %>% 
  filter(DATE<=mci_date) %>% 
  mutate(state="MCI_inst",
         state_date=mci_date,
         next_state_date=DIAGNOSDATUM) %>% 
  rbind(NPR_MCI_svedem_inst %>% 
          filter(DATE>=DIAGNOSDATUM|is.na(DATE)) %>% 
          mutate(state="MCI",
                 state_date=mci_date,
                 next_state_date=DIAGNOSDATUM)) %>% 
  rbind(NPR_MCI_svedem_inst %>% 
          filter(DATE<DIAGNOSDATUM&DATE>mci_date) %>% 
          mutate(state="MCI",
                 state_date=mci_date,
                 next_state_date=DATE)) %>% 
  rbind(NPR_MCI_svedem_inst %>% 
          filter(DATE<DIAGNOSDATUM&DATE>mci_date) %>% 
          mutate(state="MCI_inst",
                 state_date=DATE,
                 next_state_date=DIAGNOSDATUM)) %>% 
  arrange(LOPNR,DATE)

# costs for MCI from SveDem
mci_svedem_costs<-NPR_MCI_svedem_inst %>% 
  inner_join(allcost, relationship = "many-to-many") %>% 
  filter(costdate>=state_date & costdate<next_state_date) %>% 
  group_by(LOPNR,hcr,state) %>% 
  summarize(sumcost=sum(cost)) %>% 
  ungroup()

mci_svedem_costs<-mci_svedem_costs %>% 
  distinct(LOPNR,state) %>% 
  cross_join(costtypes) %>% 
  left_join(mci_svedem_costs %>% dplyr::select(LOPNR,state,hcr,sumcost)) %>% 
  mutate(sumcost=ifelse(is.na(sumcost), 0, sumcost)) %>%
  mutate(hcr = factor(hcr, levels = c('drug','inpatient', 'outpatient',
                                      'inst', 'homecare', 'daytimeactivity', 
                                      'shortterm', 'housingsupport'))) %>% 
  pivot_wider(id_cols = c(LOPNR,state),names_from=hcr, values_from=sumcost) %>% 
  left_join(NPR_MCI_svedem_inst %>% dplyr::select(LOPNR,state,next_state_date,state_date)) %>% 
  mutate(duration=as.numeric(next_state_date-state_date)/365.25,
         # add institutionalization costs
         inst=ifelse(state=="MCI_inst",duration*1152118,0),
         # set other social costs as 0 if institutionalized
         homecare=ifelse(state=="MCI_inst",0,homecare), 
         shortterm=ifelse(state=="MCI_inst",0,shortterm), 
         daytimeactivity=ifelse(state=="MCI_inst",0,daytimeactivity), 
         housingsupport=ifelse(state=="MCI_inst",0,housingsupport), 
         totalcost=drug+inpatient+outpatient+inst+homecare+daytimeactivity+
           shortterm+housingsupport)

mci_svedem_costs<-mci_svedem_costs %>% 
  mutate(DATE=state_date,
         next_date=next_state_date) %>% 
  dplyr::select(LOPNR,state,DATE,next_date,totalcost)



save(list=c('costwide','mci_svedem_costs'), file='costdata_xx.RData')




# 5. Clean NACC data for modeling MCI transitions ----
source("NACC_mci_progression_data_prep.R")

save(list=c('mstatedata', 'eventdata','nacc_mci_select'), file='mstatedata_xx.RData')  