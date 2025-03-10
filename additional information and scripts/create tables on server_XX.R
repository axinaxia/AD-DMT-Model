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


svedem_sql<-dbConnect(MySQL(),
                      user=rstudioapi::askForPassword(prompt = "Username"),
                      password=rstudioapi::askForPassword(prompt = "Password"),
                      host='h1cogbase01.nvs.ki.se',
                      port=3306,
                      dbname='SVEDEM_2021')

INCLUSION<-dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOS, d_YOB, DIAGNOSDATUM, SEX
                      from SVEDEM_2021.GRUNDREG__2;') %>% 
  as.data.frame() %>% 
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOS, d_YOB, DIAGNOSDATUM, SEX 
                   from SVEDEM_2021.GRUNDREG__3;') %>% 
          as.data.frame()) %>% 
  mutate(DIAGNOSDATUM=as.Date(DIAGNOSDATUM)) %>% 
  group_by(LOPNR) %>% 
  filter(DIAGNOSDATUM==min(DIAGNOSDATUM)) %>% # keep the earliest diagnosis
  ungroup() 


DEMO<-COG<-dbGetQuery(svedem_sql, 'SELECT LOPNR, DIAGNOSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
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
                   'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, NULL AS MOCA, NULL AS MOCA_VARDE
          FROM SVEDEM_2021.UPPF_HE_2;') %>%
          as.data.frame()) %>%
  rbind(dbGetQuery(svedem_sql, 
                   'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, NULL AS MOCA, NULL AS MOCA_VARDE
          FROM SVEDEM_2021.UPPF_SA_2;') %>%
          as.data.frame())




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
                   'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, NULL AS MOCA, NULL AS MOCA_VARDE
          FROM SVEDEM_2021.UPPF_HE_2;') %>%
          as.data.frame()) %>%
  rbind(dbGetQuery(svedem_sql, 
                   'SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, NULL AS MOCA, NULL AS MOCA_VARDE
          FROM SVEDEM_2021.UPPF_SA_2;') %>%
          as.data.frame())



MORT<-dbGetQuery(svedem_sql, 'SELECT LOPNR, DODSDAT
            FROM SVEDEM_2021.DORS;') %>% 
  as.data.frame() %>%
  rbind(dbGetQuery(svedem_sql, 'SELECT LOPNR, DODSDAT
            FROM SVEDEM_2021.DORS_AVI;')) %>%
  mutate(DODSDAT_edit=ifelse(substr(DODSDAT,5,8)=='0000', paste0(substr(DODSDAT,1,4),'0101'), DODSDAT),
         DODSDAT_edit=ifelse(substr(DODSDAT_edit,7,8)=='00', paste0(substr(DODSDAT,1,6),'01'), DODSDAT_edit),
         DATE=as.Date(DODSDAT_edit, format='%Y%m%d'))




INST<-dbGetQuery(svedem_sql, 'SELECT LOPNR, min(PERIOD) AS PERIOD 
            FROM SVEDEM_2021.SOL
            WHERE BOFORM=2
            GROUP BY LOPNR;') %>%
  as.data.frame() %>% 
  mutate(DATE=as.Date(paste0(PERIOD,'01'), format='%Y%m%d')) %>%
  group_by(LOPNR) %>% summarize(DATE=min(DATE))

sos_flexlab<-dbGetQuery(svedem_sql, 'SELECT * from SVEDEM_2021.FLEXLAB;')
source('biomarker_cutoff.R')

SVEDEM<-list()

SVEDEM$INCLUSION<-INCLUSION %>% 
  left_join(bm.path,by=c('LOPNR'='LopNr')) %>% 
  mutate(amyloid=case_when((Amy42.40==0&Amy42.pTau==0)|(is.na(Amy42.40)&Amy42.pTau==0)|
                             (is.na(Amy42.pTau)&Amy42.40==0)~0,
                           Amy42.40==1|Amy42.pTau==1~1,
                           is.na(Amy42.40)&is.na(Amy42.pTau)~NA))

SVEDEM$COG<-COG
SVEDEM$MORT<-MORT
SVEDEM$INST<-INST

