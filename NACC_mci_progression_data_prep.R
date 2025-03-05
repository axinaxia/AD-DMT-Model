library(haven)
library(ggplot2)
library(survival)
library(survminer)
library(pbapply) 
library(RMySQL)
library(tidyverse)
library(icenReg)

nacc2024<-dbConnect(MySQL(),
                    user=rstudioapi::askForPassword(prompt = "Username"),
                    password=rstudioapi::askForPassword(prompt = "Password"),
                    host='h1cogbase01.nvs.ki.se',
                    port=3306,dbname="NACC_2024")

dbListTables(nacc2024)

NACC_A2024<-dbGetQuery(nacc2024,"SELECT * FROM NACC_2024.NACC_A;")
NACC_B2024<-dbGetQuery(nacc2024,"SELECT * FROM NACC_2024.NACC_B;")

# ************************************************************************************************
# Information in NACC datasets
# Variables for dates of visits: NACCVNUM, VISITYR, VISITMO, VISITDAY
# demographics: BIRTHYR, BIRTHMO, SEX, RACE (1=white)
# Mortality: NACCDIED, NACCYOD, NACCMOD
# Attrition: NACCACTV
# 0 = Died, discontinued, lost to follow-up
# 1 = Annual follow-up (no discontinuation/loss to follow-up or minimal contact)
# 2 = Minimal contact with Center, no annual follow-up
# Cognitive states: CDRSUM (CDR-SB), CDRGLOB (global CDR)
# NACCUDSD (Cognitive status at UDS visit NACC, 3=MCI, 4=dementia), 
# Ab pathology: AMYLPET, AMYLCSF
# Comorbidity: CBSTROKE, HXSTROKE, STKIMAG
# Genetic: NACCNE4S
# 0 = no e4 allele
# 1 = 1 copy of e4 allele
# 2 = 2 copies of e4 allele
# institutionalization: NACCNURP,NACCNRMO,NACCNRDY,NACCNRYR
# ************************************************************************************************

nacc_com<-NACC_A2024 %>% 
  select(NACCID,NACCVNUM,VISITYR,VISITMO,VISITDAY,
         BIRTHYR, BIRTHMO, SEX, RACE,CBSTROKE, HXSTROKE) %>% 
  left_join(NACC_B2024 %>% 
              select(NACCID,NACCVNUM,NACCNE4S,NACCALZD,NACCALZP,
                     NACCUDSD,NACCACTV,NACCDIED, NACCYOD, NACCMOD,
                     AMYLPET, AMYLCSF, STKIMAG,
                     NACCNURP,NACCNRMO,NACCNRDY,NACCNRYR)) %>% 
  arrange(NACCID,NACCVNUM) %>% 
  mutate(date=as.Date(paste(VISITYR,VISITMO,VISITDAY, sep='-')),
         death_date=ifelse(NACCDIED==1,
                           ymd(paste(NACCYOD,ifelse(NACCMOD==88,01,NACCMOD),'01',sep='-')),
                           NA),
         death_date=as.Date(death_date),
         inst_date=ifelse(NACCNURP==1,
                          ymd(paste(NACCNRYR,NACCNRMO,NACCNRDY,sep='-')),
                          NA),
         inst_date=as.Date(inst_date),
         AGE=VISITYR-BIRTHYR,
         SEX=factor(SEX,levels = c(2,1),labels = c("FEMALE","MALE"))) %>% 
  group_by(NACCID) %>% 
  mutate(AMYLPET=max(AMYLPET,na.rm = T), 
         AMYLCSF=max(AMYLCSF,na.rm = T)) %>% 
  ungroup

nacc_mci<-nacc_com %>% 
  filter(NACCUDSD==3) %>%
  mutate(mci_date=as.Date(paste(VISITYR,VISITMO,VISITDAY, sep='-'))) %>% 
  group_by(NACCID) %>% 
  slice_min(mci_date) %>% 
  ungroup %>% 
  select(NACCID,mci_date) %>%
  left_join(nacc_com,by=c("NACCID"="NACCID")) %>% 
  arrange(NACCID,NACCVNUM) %>% 
  filter(date>=mci_date) %>% 
  rename(naccid=NACCID,
         naccnum=NACCVNUM) %>% 
  select(naccid,naccnum,date,mci_date,NACCACTV,NACCUDSD,death_date,AGE,SEX,NACCNE4S,
         AMYLPET,AMYLCSF,NACCALZD,NACCALZP,RACE,HXSTROKE,inst_date) %>% 
  mutate(inst=as.numeric(inst_date<=mci_date),
         inst=ifelse(!is.na(inst),0,inst)) %>% 
  arrange(naccid,date) %>%
  group_by(naccid) %>%
  mutate(n_obs = n(),
         prop_ad_etio=sum(NACCALZD==1)/n_obs,
         prop_ad_ci=sum(NACCALZP==1)/n_obs) %>%
  ungroup() %>%
  filter(n_obs > 1) 


nacc_mci<-nacc_mci %>% 
  group_by(naccid) %>% 
  mutate(prog=ifelse(NACCUDSD==4,1,0),
         
         dem_date=ifelse(NACCUDSD==4,date,NA),
         dem_date=as.Date(min(dem_date,na.rm = T)),
         dem_date=as.Date(ifelse(is.infinite(dem_date),NA,dem_date)),
         
         end_date=pmin(dem_date,death_date,na.rm = T),
         end_date=ifelse(is.na(end_date),max(date),end_date),
         end_date=as.Date(end_date),
         fu=as.numeric((end_date-mci_date)/365.25),
         
         death=as.numeric(death_date==end_date),
         death=ifelse(is.na(death),0,1)) %>% 
  ungroup %>% 
  distinct(naccid,prog,end_date,.keep_all = T) %>% 
  group_by(naccid) %>% 
  slice_max(prog) %>% 
  ungroup %>% 
  left_join(nacc_mci %>% 
              group_by(naccid) %>% 
              filter(NACCUDSD==3) %>% 
              summarise(last_mci_date=max(date)) %>% 
              ungroup) %>% 
  mutate(last_mci_date=pmin(last_mci_date,dem_date,end_date),
         last_mci_interval=as.numeric((end_date-last_mci_date)/365.25))

table(nacc_mci$prog)
table(nacc_mci$death)

nacc_mci_select<-nacc_mci %>% 
  filter(AMYLPET==1|AMYLCSF==1|prop_ad_etio>0.5|prop_ad_ci>0.5,
         RACE==1,
         HXSTROKE!=2,
         NACCNE4S<=1,
         AGE>=50&AGE<=90)

n_distinct(nacc_mci$naccid)
n_distinct(nacc_mci_select$naccid)

table(nacc_mci_select$prog)
table(nacc_mci_select$death)


ggsurvplot(
  survfit(Surv(fu,prog) ~ 1, data = nacc_mci),
  break.time.by = 1, break.y.by=0.05, # 1-year intervals (assuming 'fu' is in days)
  ggtheme = theme_minimal() + theme(panel.grid.major = element_line(color = "gray80"))
)


ggsurvplot(
  survfit(Surv(fu,prog) ~ 1, data = nacc_mci_select),
  break.time.by = 1, break.y.by=0.05, # 1-year intervals (assuming 'fu' is in days)
  ggtheme = theme_minimal() + theme(panel.grid.major = element_line(color = "gray80"))
)

# impute event time to deal with interval censoring
nacc_mci_imp<-nacc_mci_select %>% 
  filter(prog==1) %>% 
  mutate(left=ifelse(last_mci_interval<=1.5,as.numeric(last_mci_date-as.Date("2005-01-01")),
                     as.numeric(end_date-as.Date("2005-01-01"))-1.5*365.25),
         right=as.numeric(end_date-as.Date("2005-01-01"))) %>% 
  select(naccid,left,right,mci_date,end_date,SEX,AGE,NACCNE4S)


par_fit <- ic_par(cbind(left,right) ~ SEX+AGE+NACCNE4S, model = "ph", 
                  data = nacc_mci_imp,
                  dist = "loglogistic")

imputedValues<-imputeCens(par_fit, imputeType = 'fullSample',samples=2000)


nacc_mci_imp<-nacc_mci_imp %>% 
  cbind(imputedValues=as.numeric(rowMeans(imputedValues)))


nacc_mci_select<-nacc_mci_select %>% 
  left_join(nacc_mci_imp %>% select(naccid,imputedValues)) %>%
  mutate(end_date_imp=ifelse(is.na(imputedValues),as.Date(end_date),
                             as.Date(imputedValues+as.Date("2005-01-01"))),
         end_date_imp=as.Date(end_date_imp),
         fu_imp=as.numeric((end_date_imp-mci_date)/365.25),
         death_new=ifelse(death==1&fu_imp<fu,0,death))


ggsurvplot(
  survfit(Surv(fu_imp,prog) ~ 1, data = nacc_mci_select),
  break.time.by = 1, break.y.by=0.05, 
  ggtheme = theme_minimal() + theme(panel.grid.major = element_line(color = "gray80"))
)

ggsurvplot(
  survfit(Surv(fu_imp,death_new) ~ 1, data = nacc_mci_select),
  break.time.by = 1, break.y.by=0.05, 
  ggtheme = theme_minimal() + theme(panel.grid.major = element_line(color = "gray80"))
)

table(nacc_mci_select$prog)
table(nacc_mci_select$death_new)