rm(list = ls())

options(scipen=999)

#LOAD LIBRARIES----
library(tidyverse)
library(stringr)
library(ggplot2)
library(readxl)
library(stringr)
library(survival)
library(data.table)
library(nlme)
library(consort)
library(lubridate)

#LOAD DATA----
path<- '~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/CSF-registersamkörning/Data'
setwd(path)

fallkontroll<-read_excel('Leverans/AW_Lev_fall_och_kontroller.xlsx') #case and control
#load('matched_data.Rdata') #case and control matched on age, sex, lan and no. of comorbidities
#load('comorb_summary.RData') #comorbidity summary
load('svedem221114.Rdata')
load('bmpath.RData') #biomarker data #cut point defined at 14.9 amy42.ptau ratio
load('lmed.RData')
cpi <- read_excel("CPI.xlsx") %>% mutate(year=as.numeric(year))


svedem$grundreg<-rbind(svedem$sos_dat_669_grundreg__2.sas7bd %>% select(LopNr, MMSESR_VARDE, DIAGNOS, DIAGNOSDATUM), 
                       svedem$sos_dat_669_grundreg__3.sas7bd%>% select(LopNr, MMSESR_VARDE, DIAGNOS, DIAGNOSDATUM)) %>%
  group_by(LopNr) %>% summarize(MMSE=max(MMSESR_VARDE),
                                DIAGNOS=max(DIAGNOS),
                                DIAGNOSDATUM=min(DIAGNOSDATUM)) %>% ungroup()


deaths<-rbind(svedem$ut_r_dors_11408_2021.sas7bdat %>% select(LopNr, DODSDAT),
              svedem$ut_r_dors_avi_11408_2021.sas7bdat%>% select(LopNr, DODSDAT)) %>%
  mutate(deathdate=as.Date(DODSDAT, format='%Y%m%d')) %>% select(-DODSDAT)

diags<-c(
  'DEMENS_UNS',
  'VASKULAR_DEMENS_INKL_SUBKORTIKAL_VASKULAR_DEMENS',
  'BLANDDEMENS_VID_ALZHEIMERS_SJUKDOM_VASKULAR_DEMENS',
  'DEMENS_VID_ALZHEIMERS_SJUKDOM_SEN_DEBUT' ,
  'DEMENS_VID_ALZHEIMERS_SJUKDOM_TIDIG_DEBUT',
  'DEMENS_VID_PARKINSONS_SJUKDOM',
  'FRONTOTEMPORAL_DEMENS' ,
  'LEWY_BODY_DEMENS' ,
  'LINDRIG_KOGNITIV_STORNING',
  'MCI',
  'OVRIG_DEMENS' 
)
diaggrps<-c(
  'Unspecified',
  'VaD',
  'Mixed AD/VaD',
  'AD',
  'Early onset AD',
  'PD',
  'FTD',
  'LBD',
  'MCI',
  'MCI',
  'Other')

# DEMOGRAPHICS-----

cutdate=as.Date('2013-01-01') #SOL started in 2012JUL

### MCI diagnosis in medical register, include only observations which are within 5 years before the diagnosis of mild dementia
mci <- rbind(svedem$ut_r_par_sv_11408_2021.sas7bdat %>% select(LopNr,INDATUM,HDIA) %>% mutate(vtyp='sv'),
             svedem$ut_r_par_ov_11408_2021.sas7bdat %>% select(LopNr,INDATUM,HDIA) %>% mutate(vtyp='ov'))%>% 
  filter(HDIA =='F067') %>% 
  left_join(svedem$grundreg %>% select(LopNr, DIAGNOSDATUM, MMSE), by='LopNr') %>%
  filter(MMSE>20 & (INDATUM<DIAGNOSDATUM) & 
           (DIAGNOSDATUM-INDATUM < 365*5)) %>% # sensitivity analysis here 7 years, 3 years
  arrange(LopNr, INDATUM) %>%
  group_by(LopNr) %>%
  mutate(mcidate = min(INDATUM), #first date of MCI diagnosis
         nextdate = pmin(lead(INDATUM,1), DIAGNOSDATUM, na.rm=TRUE),
         time=as.numeric(DIAGNOSDATUM)-as.numeric(INDATUM),
         nexttime=lead(time),
         duration= ifelse(!is.na(nexttime), time-nexttime, as.numeric(nextdate)-as.numeric(INDATUM))) %>% 
  ungroup() %>% 
  filter(mcidate > cutdate)

#check number of unique patients in mci
n_distinct(mci$LopNr) # 5889 patients with mci due to all causes

# check duration of MCI
mci %>% mutate(
         duration = (as.numeric(DIAGNOSDATUM)-as.numeric(mcidate))/365.25) %>% 
  #include only the first observation of each lopnr
  filter(!duplicated(LopNr)) %>% 
  summarise(mean_time=mean(duration),sd_time=sd(duration),
            median_time=median(duration),IQR_U=quantile(duration,3/4),IQR_L=quantile(duration,1/4),
            min_time=min(duration),max_time=max(duration))


# create a new diagnosis date including MCI diagnosis date
grundreg <- svedem$grundreg %>% select(LopNr, DIAGNOSDATUM, MMSE, DIAGNOS) %>% 
  left_join(mci %>% select(LopNr, mcidate) %>% distinct, by='LopNr') %>%
  mutate(newdiagdate=ifelse(!is.na(mcidate),mcidate,DIAGNOSDATUM)) 

grundreg$newdiagdate<-as.Date(grundreg$newdiagdate)


#create base population
pop<-bind_rows(fallkontroll %>% rename(LopNr=lopnr_fall) %>% mutate(ref_case=LopNr) %>% select(LopNr, ref_case, fodelsear, kon, lan) %>% distinct %>% mutate(group='Case'),
               fallkontroll %>% rename(LopNr=lopnr_kontroll) %>% mutate(ref_case=lopnr_fall) %>% select(LopNr, ref_case, fodelsear, kon, lan) %>% distinct %>% mutate(group='Control'),) %>% 
#pop<-matched_df %>%   #sensitivity analysis of population matched with number of multimorbidity
  left_join(grundreg %>% select(LopNr, DIAGNOSDATUM, newdiagdate,MMSE, DIAGNOS), by=c('ref_case'='LopNr')) %>% 
  filter(!is.na(DIAGNOSDATUM) & DIAGNOSDATUM>='2007-01-01') %>%
  left_join(deaths) %>%
  mutate(censdate=as.Date('20201231', format='%Y%m%d'),
         leftcensdate=as.Date('20070101', format='%Y%m%d'),
         lastdate=case_when(!is.na(deathdate) & deathdate<=censdate ~ deathdate,T ~ censdate),
         death=ifelse(!is.na(deathdate) & deathdate<=lastdate,1,0),
         baseage=as.numeric(substr(DIAGNOSDATUM,1,4))-as.numeric(fodelsear), 
         baseagegrp=cut(baseage, breaks=c(0,65,70,75,80,85,105)), #there are too few observations in 90+ amy+AD, combine the agegrp as85+
         sex=factor(kon, levels=c('1', '2'), labels=c('Male', 'Female')),
         state=case_when(DIAGNOS=='LINDRIG_KOGNITIV_STORNING' ~ 'MCI', 
                         MMSE<10 ~ 'MMSE 0-9',
                         MMSE<=20 ~ 'MMSE 10-20',
                         #MMSE<26 ~ 'MMSE 21-25',
                         MMSE>=21 ~'MMSE 21-30'),
         diaggrp=factor(DIAGNOS, levels=diags, labels=diaggrps),
         diag_year=(substr(DIAGNOSDATUM,1,4)))


# rename the levels in agegrp
levels(pop$baseagegrp)<-c('0-65', '66-70', '71-75', '76-80', '81-85', '>85')

# combine with comorbidity data
#load('multimorb.RData')
#pop<-pop %>% left_join(multimorb %>% select(LopNr, n_comorb), by='LopNr')

#Select population
p <- pop  %>% left_join(bm.path) %>%
  mutate(flag_diag=as.numeric(DIAGNOS %in% c('BLANDDEMENS_VID_ALZHEIMERS_SJUKDOM_VASKULAR_DEMENS', 
                                             'DEMENS_VID_ALZHEIMERS_SJUKDOM_SEN_DEBUT', 
                                             'DEMENS_VID_ALZHEIMERS_SJUKDOM_TIDIG_DEBUT',
                                             'LINDRIG_KOGNITIV_STORNING')),
         flag_amy=ifelse(is.na(Amy42.pTau),99,
                         ifelse(!is.na(Amy42.pTau) & Amy42.pTau==1,1,0)),
         flag_date=as.numeric(DIAGNOSDATUM>=cutdate & DIAGNOSDATUM<=as.Date('2020-12-31') & (is.na(deathdate) | DIAGNOSDATUM<deathdate)),
         follow = (as.numeric(lastdate)-as.numeric(newdiagdate))/365.25) %>% 
  filter(flag_diag==1 & flag_date==1) %>%  #only include those with clinical AD diagnosis and diagnosed between 2013-2020
  filter(!is.na(state)) %>% 
  select(-Amy42.40, -Amy42.pTau,-date_Amy42.40, -date_Amy42.pTau)

# remove ref_case that doesn't have a control
p <- p %>% group_by(ref_case) %>%
  filter(n() == 3 & any(group == "Case") & any(group == "Control")) %>% ungroup()

table(p$group)


#Combine with MMSE data
mmse<-rbind(svedem$sos_dat_669_grundreg__2.sas7bd %>% select(LopNr, MMSESR, MMSE=MMSESR_VARDE, DATE=DIAGNOSDATUM),
            svedem$sos_dat_669_grundreg__3.sas7bd %>% select(LopNr, MMSESR, MMSE=MMSESR_VARDE, DATE=DIAGNOSDATUM),
            svedem$sos_dat_669_uppf__2.sas7bdat %>% select(LopNr, MMSESR, MMSE=MMSESR_VARDE, DATE=UPPFOLJNINGSDATUM),
            svedem$sos_dat_669_uppf__3.sas7bdat %>% select(LopNr,MMSESR, MMSE=MMSESR_VARDE, DATE=UPPFOLJNINGSDATUM),
            svedem$sos_dat_669_uppf_sa_2.sas7bdat %>% select(LopNr, MMSESR, MMSE=MMSESR_VARDE, DATE=UPPFOLJNINGSDATUM),
            svedem$sos_dat_669_uppf_he_2.sas7bdat %>% select(LopNr, MMSESR,MMSE=MMSESR_VARDE, DATE=UPPFOLJNINGSDATUM)) %>% 
  mutate(MMSE=case_when(MMSESR=='EJ_TESTBAR' ~0,T ~MMSE)) %>% drop_na %>% 
  group_by(LopNr, DATE) %>% 
  summarize(MMSE=max(MMSE)) %>% ungroup 

#Comcare data (both case and control)
sol1<-svedem$ut_r_sol_11408_2021.sas7bdat %>%  #Socialtjänstinsatser till äldre och personer med funktionsnedsättning
  select(LopNr, PERIOD, BOFORM, HTJ, HTJTIM, MATD, TRYGG, DAGV, BOSTO,KORTTID, KORTMAN) %>%
  # HTJ - homecare, HTJTIM - home care hour during the month, MATD - home service for food dsitritution, TRYGG - assistance for security alarm
  # DAGV - assistance for day activity, BOSTO - housing support for daily life, 
  # KORTTID - short-term care, KORTMAN - no of KORTTID day received during the month
  mutate(costdate=as.Date(paste0(PERIOD,'01'), format='%Y%m%d'),
         inst=as.numeric(BOFORM==2),
         # set if BOFORM = 2, other support = zero
         hhtim=case_when(BOFORM == 2 ~ 0,
                         BOFORM!=2 & HTJ==1 & HTJTIM<777 ~ HTJTIM,
                         BOFORM!=2 &  HTJ==1 & HTJTIM==777 ~ 24*30,T ~0), #
         shortterm_days=case_when(BOFORM == 2 ~ 0,
                                  KORTMAN %in% c('88', '99','32','60','62',' .','-1', '-2', '-3', '-7', '-9', '0', '') ~ 0,
                                  T ~ as.numeric(KORTMAN)),
         support=case_when(BOFORM == 2 ~ 0,BOSTO=="1"~ 1, T ~ 0),
         dagv=case_when(BOFORM == 2 ~ 0, DAGV=='1'~ 1,T ~ 0)) %>% 
  filter(costdate>=cutdate & costdate<as.Date('2021-01-01')) %>% 
  # remove duplicate observations with same LopNr and costdate
  group_by(LopNr, costdate) %>%
  filter(row_number()==1) %>% ungroup()

#take out inst date
instdate<-sol1 %>% 
  filter(inst==1) %>% 
  group_by(LopNr) %>%
  summarize(instdate=min(costdate))


## combine longitudinal data
MMSE <- rbind(mci %>% select(LopNr, DATE = INDATUM) %>% mutate(MMSE = 31), #arbitary number to represent MCI
              mmse) %>% 
  #create a variale called mmse valide date = last visit date + 2 years
  group_by(LopNr) %>%
  arrange(LopNr, DATE) %>%
  mutate(mmsevaliddate=ifelse(!is.na(lead(DATE)), pmin(lead(DATE), DATE+730), DATE+730),
         mmsevaliddate = as.Date(mmsevaliddate)) %>% ungroup()

#combine with institution
m <- MMSE %>% select(LopNr, date=DATE, MMSE,mmsevaliddate) %>% mutate(inst=NA) %>% 
  rbind(instdate %>% select(LopNr, date=instdate) %>% mutate(inst=1, MMSE=NA, mmsevaliddate=NA)) %>%
  group_by(LopNr) %>% 
  arrange(LopNr, date, inst) %>% 
  fill(inst, .direction='down') %>%
  arrange(LopNr, date, MMSE) %>%
  fill(mmsevaliddate, .direction='down') %>% 
  #impute MMSE value if inst=1 and date < mmsevaliddate
  mutate (MMSE = case_when(is.na(MMSE) & inst==1 & date<mmsevaliddate ~ lag(MMSE),T ~ MMSE)) %>% 
  #select MMSE for the study population
  inner_join(p %>% select(LopNr, DIAGNOS,newdiagdate,baseagegrp, flag_amy, group) %>% filter(group=='Case')) %>%
  mutate(inst=ifelse(is.na(inst),0,inst),
         state=case_when(MMSE ==31 | DIAGNOS =='LINDRIG_KOGNITIV_STORNING' ~ 'MCI',
                         MMSE<10 ~ 'MMSE 0-9',MMSE<=20 ~ 'MMSE 10-20', MMSE<26 ~'MMSE 21-25', MMSE>=26 ~'MMSE 26-30'), #for HE model
                         #MMSE<10 ~ 'Severe',MMSE<=20 ~ 'Moderate', MMSE>=21 ~'Mild'),
         setting=ifelse(inst==1,'Institution', 'Community')) %>% 
  filter(!is.na(MMSE)) %>% ungroup() %>% select(-group)


m2 <- m %>% 
  group_by(LopNr) %>%
  arrange(LopNr, date) %>%
  # make state ordinal 
  mutate(state = factor(state, levels=c('MCI', 'MMSE 26-30', 'MMSE 21-25', 'MMSE 10-20', 'MMSE 0-9'), ordered=TRUE)) %>%
  #mutate(state = factor(state, levels=c('MCI', 'Mild', 'Moderate', 'Severe'), ordered=TRUE)) %>%
  # exlude the obervations that have lower/milder state than previous state
  filter(is.na(lag(state)) | lag(state) <= state) %>% ungroup()

table(m2$state)

#combine longitudinal MMSE with population data
m3 <- p %>% select(LopNr, group, ref_case, fodelsear,sex, deathdate) %>%
  inner_join(m2, by=c('ref_case'= 'LopNr'), relationship='many-to-many') %>%
  mutate(age=as.numeric(substr(date,1,4))-as.numeric(fodelsear),
         agegrp=cut(age, breaks=c(0,65,70,75,80,85,105)),
         censdate = ifelse(!is.na(deathdate),pmin(deathdate, mmsevaliddate, as.Date('2020-12-31')), pmin(mmsevaliddate, as.Date('2020-12-31'))),
         censdate = as.Date(censdate),
         follow = (as.numeric(censdate) - as.numeric(newdiagdate))/365) %>% 
  filter(date< censdate & !is.na(agegrp)) %>%
  group_by(LopNr) %>%
  arrange(LopNr, date) %>%
  mutate(time=as.numeric(date)-as.numeric(newdiagdate),
         nexttime=ifelse(!is.na(lead(time)), lead(time), as.numeric(censdate)-as.numeric(newdiagdate)),
         state_duration = pmin(730,nexttime - time))%>% ungroup() %>% 
  filter(time>=0 ,state_duration>0) %>% 
  select(-fodelsear) 


#DRG weights and prices-----
drgcost <- read_excel("KPP/Retrospektiva DRG vikter somatisk vard 2012-2021_uppdaterad.xlsx", 
                      n_max=1, col_types=c('skip', 'skip', rep('text',10 )))%>%
  pivot_longer(cols=as.character(2012:2021), names_to='year', values_to='drgcost') %>% 
  mutate(year=as.numeric(year),
         drgcost=as.numeric(drgcost)) 

# convert drgcost to 2021 value
drgcost <- drgcost %>% left_join(cpi, by='year') %>% 
  mutate(cpi2021= as.numeric('343.19'),
        factor= cpi2021/CPI,
        drgcost= drgcost*factor) %>% select(year, drgcost)

drg<-read_excel('KPP/Retrospektiva DRG vikter somatisk vard 2012-2021_uppdaterad.xlsx', skip=2,  
                col_types=c('text', 'text', rep('numeric',10)))  %>% 
  pivot_longer(cols=as.character(2012:2021),names_to='year', values_to='weight') %>% 
  mutate(year=as.numeric(year)) %>%
  select(DRG, year, weight)

#get rid of missing weights
drg<-drg %>% drop_na 

#Create table with all DRGs and all years
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

#Add costs per DRG point

#drgnew<-drgnew %>% left_join(drgcost) %>% mutate(cost=drgcost*weight)

#Group remaining encounters by MVO (medical acitivity area), assign mean cost within MVO (for each medical activity)
mvodata<-rbind(
  svedem$ut_r_par_sv_11408_2021.sas7bdat %>% mutate(year=as.numeric(substr(INDATUM,1,4))) %>% select(year,DRG, MVO) %>% distinct %>% mutate(vtyp='sv'),
  svedem$ut_r_par_ov_11408_2021.sas7bdat %>% mutate(year=as.numeric(substr(INDATUM,1,4))) %>% select(year,DRG, MVO) %>% distinct %>% mutate(vtyp='ov') )

mvocost<-mvodata %>% full_join(drgnew) %>% full_join(drgcost) %>% mutate(mvocost=weight*drgcost) %>%
  filter(!is.na(mvocost)) %>% group_by(vtyp, MVO) %>%
  summarize(mvocost=mean(mvocost, na.rm=T))

#DIRECT MEDICAL COSTS----
#Assign costs to all outpatient and inpatient care
medcost<-rbind(
  svedem$ut_r_par_sv_11408_2021.sas7bdat %>% mutate(year=as.numeric(substr(INDATUM,1,4))) %>% select(LopNr,INDATUM,year,DRG, MVO) %>% mutate(vtyp='sv'),
  svedem$ut_r_par_ov_11408_2021.sas7bdat %>% mutate(year=as.numeric(substr(INDATUM,1,4))) %>% select(LopNr,INDATUM,year,DRG, MVO)%>% mutate(vtyp='ov')) %>% 
  filter(INDATUM>=cutdate & INDATUM<as.Date('2021-01-01'))

medcost<-medcost %>% left_join(drgnew) %>% left_join(drgcost) %>% left_join(mvocost) %>%
  mutate(cost=weight*drgcost) 
 

# 8271016 observations have missing drgweight in medcost which is stragne because all the DRG weight was filled in already, check
n_distinct(medcost$DRG) #1428
n_distinct(drgnew$DRG) #1525
# which was because the DRG code in drgnew was not exhaustive, so we need to impute the missing cost with mvo cost
#mStats::summ(medcost,cost)

#impute missing medcost with mvocost
medcost[is.na(medcost$cost),]$cost<-medcost[is.na(medcost$cost),]$mvocost

Medcost<-medcost %>% 
  inner_join(p %>% select(LopNr, newdiagdate)) %>% 
  filter(INDATUM>=newdiagdate) %>%
  mutate(time=as.numeric(INDATUM-newdiagdate),
         vtyp= ifelse(vtyp == 'ov', 'outpatient', 'inpatient')) %>%
  select(LopNr, name=vtyp, time, costdate=INDATUM,cost) 


#DIRECT NON-MEDICAL COSTS-----
comcare<-sol1 %>% 
  inner_join(p %>% select(LopNr, newdiagdate)) %>% 
  filter(costdate>=newdiagdate) %>%
  mutate(time=as.numeric(costdate-newdiagdate),
         cost_inst=inst*1056142/12,
         cost_homecare=ifelse(BOFORM==2,0,hhtim)*592, #people with inst cost does have homehelp cost and other costs
         cost_daytimeactivity=dagv*113800/12,
         cost_shortterm=shortterm_days*3300,
         cost_housingsupport=support*62874/12) %>% 
  select(LopNr, costdate, time,starts_with('cost_')) %>%
  pivot_longer(cols=starts_with('cost_'), names_prefix = 'cost_', values_to='cost') %>% 
  filter(cost >0)

instdaycost <- 1056142/365

#drug cost
drugcost <- lmed %>% 
  inner_join(p %>% select(LopNr, newdiagdate)) %>%
  mutate(costdate=as.Date(EDATUM, format='%Y%m%d'),
         time=as.numeric(costdate-newdiagdate),
         name='drugs') %>% 
  filter(costdate>=newdiagdate) %>%
  select(LopNr, costdate, time, name, TKOST) %>% 
  filter(costdate>=cutdate & costdate<as.Date('2021-01-01'))

# convert drugcost to 2021 value
drugcost <- drugcost %>% 
  mutate(year=as.numeric(substr(costdate,1,4))) %>%
  left_join(cpi, by='year') %>% 
  mutate(cpi2021= as.numeric('343.19'),
         factor= cpi2021/CPI,
         cost= TKOST*factor) %>% select(LopNr, costdate, time, name, cost) %>% 
  filter(cost>0)

#TOTAL COSTS-----
allcost<-rbind(Medcost ,drugcost,comcare) %>% rename(costtime=time) 

costtypes<-allcost %>% select(name) %>% distinct

# merginng longitudinal mmse with all costs
statecost<- m3 %>% select(LopNr, group, agegrp, sex, flag_amy, state,setting,date,time, nexttime,state_duration) %>% 
  left_join(allcost, by='LopNr', relationship = "many-to-many") %>% 
  filter(costtime>=time & costtime<(time+state_duration)) %>% 
  group_by(LopNr, group, time, nexttime,state_duration,state, flag_amy,name,setting) %>% 
  summarize(sumcost=sum(cost)) %>% ungroup()

# set if setting = institution & name=homecare or shortterm or daytimeactivity or housingsupport, sumcost =0 to avoid cost duplication
statecost <- statecost %>% 
  mutate(sumcost=ifelse(setting=='Institution' & 
                           name %in% c('homecare', 'shortterm', 'daytimeactivity', 'housingsupport'), 0, sumcost))

# fill in all type of costs for each visit
statecosts<-m3 %>% 
  group_by(LopNr, group, agegrp,sex, state, setting, flag_amy, time, nexttime,state_duration) %>% 
  summarize(sumtime=sum(state_duration)) %>% 
  cross_join(costtypes) %>% 
  left_join(statecost) %>% 
  mutate(sumcost=ifelse(is.na(sumcost), 0, sumcost), 
         # adjust cost of institution due to censoring (L&R) and missing costs
         sumcost=ifelse(group=='Case' & setting=='Institution' & name =='inst', instdaycost*state_duration, sumcost)) %>% 
  ungroup() %>%
  mutate(name = factor(name, levels = c('drugs','inpatient', 'outpatient','inst', 'homecare', 'daytimeactivity', 'shortterm', 'housingsupport'))) 
         
#check case, setting = institution, name=inst and sumcost =0
statecosts %>% filter(group=='Case' & setting=='Institution' & name=='homecare' & sumcost!=0) %>% nrow() 

levels(statecosts$agegrp)<-c('0-65', '66-70', '71-75', '76-80', '81-85', '>85')

#############################################################################
#make data for HE model
costwide <- statecosts %>% filter(group=='Case' &state!='MCI' & flag_amy==1) %>% 
  pivot_wider(names_from=name, values_from=sumcost) %>% 
  mutate(totalcost=drugs+inpatient+outpatient+inst+homecare+daytimeactivity+shortterm+housingsupport)

# save costwide table R
save(costwide, file='costwide.RData')


# function to calculate meancost
meancosttotal<-function(data) {
  data %>% summarize(
    sumtime=sum(sumtime),
    sumcost=sum(totalcost)) %>%
    mutate(cost=round(sumcost/(sumtime/365),0))
}
meancosttotal(costwide %>% group_by(state))

# upload costwide to server
library(DBI)
library(RMySQL)
library(rstudioapi)

con <- dbConnect(MySQL(),
                 user='sandar.aye',
                 password=rstudioapi::askForPassword(prompt = "Password"),
                 host='h1cogbase01.nvs.ki.se',
                 port=3306,
                 dbname='SVEDEM_2021')

dbWriteTable(con, 'STATECOSTWIDE', costwide, overwrite = TRUE)


rs<-dbGetQuery(con, 'show databases')
rs


#Select a database
dbExecute(con, 'USE SVEDEM_2021;')

#List tables in database
dbListTables(con)

#Retrieve new table
tbl<-dbGetQuery(con, 'SELECT * FROM SVEDEM_2021.STATECOSTWIDE;')




