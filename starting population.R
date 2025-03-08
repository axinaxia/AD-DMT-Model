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

# 2. Data for dementia diagnosis and demographic factors ----
sos_flexlab<-dbGetQuery(svedem_sql, 'SELECT * from SVEDEM_2021.FLEXLAB;')
source('biomarker_cutoff.R')



dem2020<-dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOS, d_YOB, DIAGNOSDATUM, SEX, MMSESR_VARDE,BOENDEFORM 
                      from SVEDEM_2021.GRUNDREG__2;') %>% 
  as.data.frame() %>% 
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOS, d_YOB, DIAGNOSDATUM, SEX, MMSESR_VARDE,BOENDEFORM
                   from SVEDEM_2021.GRUNDREG__3;') %>% 
          as.data.frame()) %>% 
  mutate(DIAGNOSDATUM=as.Date(DIAGNOSDATUM)) %>% 
  group_by(LOPNR) %>% 
  filter(DIAGNOSDATUM==min(DIAGNOSDATUM)) %>% # keep the earliest diagnosis
  ungroup() 

dem2020<-dem2020 %>% 
  mutate(AGE=year(DIAGNOSDATUM)-d_YOB)

dem2020<-dem2020 %>% 
  filter(year(DIAGNOSDATUM) %in% c(2019,2020))


dem2020<-dem2020 %>% 
  left_join(bm.path %>% rename(LOPNR=LopNr)) %>% 
  mutate(age_group=cut(AGE,breaks = seq(50,90,by=5)),
         age_group=ifelse(AGE==50,1,age_group),
         mild=as.numeric(MMSESR_VARDE>=21),
         ad=as.numeric(DIAGNOS %in% c("DEMENS_VID_ALZHEIMERS_SJUKDOM_SEN_DEBUT",
                                      "DEMENS_VID_ALZHEIMERS_SJUKDOM_TIDIG_DEBUT")),
         bio_ad=as.numeric(Amy42.40==1|Amy42.pTau==1),
         mild_ad=as.numeric(mild==1&(ad==1|bio_ad==1))) %>% 
  filter(!is.na(age_group),
         !is.na(mild_ad))



prev_mildad<-dem2020 %>% 
  group_by(age_group,SEX) %>% 
  summarise(n_total=n()) %>% 
  ungroup %>% 
  left_join(dem2020 %>% 
          filter(mild_ad==1) %>% 
          group_by(age_group,SEX) %>% 
          summarise(no_mildad=n()) %>% 
          ungroup) %>% 
  mutate(prop_mildad=no_mildad/n_total,
         sex=as.numeric(SEX=="FEMALE")) %>% 
  select(age_group,sex,prop_mildad)


# calculate proportion of institutionalized mild AD patients
n_mildad_inst<-n_distinct((dem2020 %>% 
         filter(mild_ad==1,
                BOENDEFORM %in% 
                  c("SARSKILT_BOENDE_PERMANENT",
                    "SARSKILT_BOENDE_PERMANENT_ANPASSAT_FOR_PERSONER_MED_DEMENSSJUKDOM")))$LOPNR)

n_mildad_inst/table(dem2020$mild_ad)[2]

# calculate proportion of institutionalized MCI patients (0)
n_mci_inst<-n_distinct((dem2020 %>% 
                             filter(DIAGNOS=="LINDRIG_KOGNITIV_STORNING",
                                    BOENDEFORM %in% 
                                      c("SARSKILT_BOENDE_PERMANENT",
                                        "SARSKILT_BOENDE_PERMANENT_ANPASSAT_FOR_PERSONER_MED_DEMENSSJUKDOM")))$LOPNR)



scb_demo<-as.data.frame(read_delim("SCB_demographics.csv", 
                                   delim = ";", escape_double = FALSE, trim_ws = TRUE))

scb_demo<-scb_demo %>% 
  mutate(prop_age=Number/sum(Number),
         age_group=case_match(Age,
                              "50-54 years"~1,
                              "55-59 years"~2,
                              "60-64 years"~3,
                              "65-69 years"~4,
                              "70-74 years"~5,
                              "75-79 years"~6,
                              "80-84 years"~7,
                              "85-89 years"~8),
         sex=as.numeric(Sex=="women")) %>% 
  select(age_group,sex,prop_age)
  
prev_mildad<-prev_mildad %>% 
  left_join(scb_demo) %>% 
  mutate(prop_mci=ifelse(age_group<=6,0.13,0.19),
         prop_dem=case_when(age_group<=4~0.006,
                            age_group %in% c(5,6)&sex==0~0.026,
                            age_group %in% c(5,6)&sex==1~0.028,
                            age_group %in% c(7,8)&sex==0~0.083,
                            age_group %in% c(7,8)&sex==1~0.098),
         prev_mildad=prop_age*prop_dem*prop_mildad,
         prev_mci=prop_age*prop_mci*0.18)

prev_mci_mildad<-prev_mildad %>% 
  select(age_group,sex,prev_mci) %>% 
  rename(prop=prev_mci) %>% 
  mutate(stage=0) %>% 
  rbind(prev_mildad %>% 
          select(age_group,sex,prev_mildad) %>% 
          rename(prop=prev_mildad) %>% 
          mutate(stage=1)) %>% 
  arrange(stage,sex,age_group) %>% 
  mutate(profile_perc=prop/sum(prop),
         profile_perc_form=round(profile_perc*100,digits = 2))


prev_mci_mildad_apoe<-prev_mci_mildad %>% 
  select(age_group,sex,stage,prop) %>% 
  mutate(apoe=0) %>% 
  rbind(prev_mci_mildad %>% 
          select(age_group,sex,stage,prop) %>% 
          mutate(apoe=1)) %>% 
  mutate(prop=ifelse(apoe==0,prop*0.41,prop*0.59))


# create a table indicate proportion of different patient profiles, using average age per age group
profile_perc<-prev_mci_mildad_apoe %>% 
  mutate(profile_perc=prop/sum(prop)) %>% 
  arrange(stage,apoe,sex,age_group) %>% 
  mutate(age=(age_group-1)*5+52.5,
         age_group=factor(age_group,levels=1:8,
                          labels=c("50-54 years","55-59 years","60-64 years",
                                   "65-69 years","70-74 years","75-79 years",
                                   "80-84 years","85-89 years"))) %>% 
  select(-prop)

sum((profile_perc %>% filter(sex==0))$profile_perc)
sum((profile_perc %>% filter(sex==1))$profile_perc)

sum((profile_perc %>% mutate(age=age*profile_perc))$age)

sum((profile_perc %>% filter(stage==0))$profile_perc)
sum((profile_perc %>% filter(stage==1))$profile_perc)

sum((profile_perc %>% filter(apoe==0))$profile_perc)
sum((profile_perc %>% filter(apoe==1))$profile_perc)


# create a table indicate proportion of different patient profiles, using detailed age
profile_perc_detailed<-prev_mci_mildad_apoe %>% 
  mutate(profile_perc=prop/sum(prop)) %>% 
  arrange(stage,apoe,sex,age_group) %>% 
  mutate(age=case_match(age_group,
                             1~50,
                             2~55,
                             3~60,
                             4~65,
                             5~70,
                             6~75,
                             7~80,
                             8~85)) %>%
  rowwise() %>%
  mutate(expanded_rows = list(0:5)) %>%
  unnest(expanded_rows) %>% 
  filter(!(age_group!=8&expanded_rows==5)) %>% 
  mutate(age=age+expanded_rows,
         profile_perc=ifelse(age_group!=8,profile_perc/5,profile_perc/6),
         age_group=factor(age_group,levels=1:8,
                          labels=c("50-54 years","55-59 years","60-64 years",
                                   "65-69 years","70-74 years","75-79 years",
                                   "80-84 years","85-89 years"))) %>% 
  select(-prop)



sum((profile_perc_detailed %>% filter(sex==0))$profile_perc)
sum((profile_perc_detailed %>% filter(sex==1))$profile_perc)

sum((profile_perc_detailed %>% mutate(age=age*profile_perc))$age)

sum((profile_perc_detailed %>% filter(stage==0))$profile_perc)
sum((profile_perc_detailed %>% filter(stage==1))$profile_perc)

sum((profile_perc_detailed %>% filter(apoe==0))$profile_perc)
sum((profile_perc_detailed %>% filter(apoe==1))$profile_perc)





# 3. Consider dropout ----
# the probability of discontinuing treatment per cycle: 0.188/18*3

pop_perc_func<-function(i){
  prob_stop<-0.188/18*3
  
  rx_cycles_table<-data.frame(rx_cycles=seq(1,i,by=1))
  rx_cycles_table<-rx_cycles_table %>% 
    mutate(prob=(1-prob_stop)^(rx_cycles-1)*prob_stop,
           cum_prob=cumsum(prob),
           prop_stop=ifelse(rx_cycles==max(rx_cycles),1-cum_prob+prob,prob))
  
  profile_discontinue_perc<-profile_perc %>% 
    select(age_group,age,sex,apoe,stage,profile_perc) %>%
    rowwise() %>%
    mutate(rx_cycles = list(1:i)) %>% 
    unnest(rx_cycles) %>% 
    left_join(rx_cycles_table %>% 
                select(rx_cycles,prop_stop)) %>% 
    mutate(pop_perc=profile_perc*prop_stop) %>% 
    select(-prop_stop)
  
  return(profile_discontinue_perc)
}



  

save(profile_perc,file="profile_perc.RData")
save(profile_perc_detailed,file="profile_perc_detailed.RData")
