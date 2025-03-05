library(haven)
library(tidyverse)
library(ggplot2)
library(readxl)
library(flextable)
library(ggsignif)
library(nlme)
library(cutpointr)
library(survival)
library(survminer)
library(flexsurv)
library(data.table)
library(scales)
library(officer)





# generate biomarker cutoff ----
sos_flexlab$res.num<-as.numeric(gsub(',','.',gsub('\\*','',sos_flexlab$Resultat)))


#Adjust for switching methods
bm.means<-sos_flexlab %>% group_by(Kortnamn) %>% summarize(mean=mean(res.num, na.rm=T))
conv.40<-bm.means[bm.means$Kortnamn=='C-AMY40',]$mean/bm.means[bm.means$Kortnamn=='qAm40xE',]$mean
conv.42<-bm.means[bm.means$Kortnamn=='C-AMY42',]$mean/bm.means[bm.means$Kortnamn=='Am42xE',]$mean
conv.nfl<-bm.means[bm.means$Kortnamn=='C-NFL',]$mean/bm.means[bm.means$Kortnamn=='qNFL/qqNFL',]$mean
conv.tau<-bm.means[bm.means$Kortnamn=='C-TAU',]$mean/bm.means[bm.means$Kortnamn=='qC-TAU',]$mean

#Group tests 
sos_flexlab$biomarker<-case_when(sos_flexlab$Kortnamn=='Am42xE' ~ 'Amy42',
                                        sos_flexlab$Kortnamn=='C-AMY40' ~ 'Amy40',
                                        sos_flexlab$Kortnamn=='C-AMY42' ~ 'Amy42',
                                        sos_flexlab$Kortnamn=='C-FOSF' ~ 'pTau',
                                        sos_flexlab$Kortnamn=='C-TAU' ~ 'Tau',
                                        sos_flexlab$Kortnamn=='qAm40xE' ~ 'Amy40',
                                        sos_flexlab$Kortnamn=='qC-AMY' ~ 'Amy42',  #Question: is this Amy 42??
                                        sos_flexlab$Kortnamn=='qC-FOSF' ~ 'pTau',
                                        sos_flexlab$Kortnamn=='qC-TAU' ~ 'Tau',
                                        sos_flexlab$Kortnamn=='qNFL/qqNFL' ~ 'NfL',
                                        sos_flexlab$Kortnamn=='C-NFL' ~ 'NfL',
                                        sos_flexlab$Kortnamn=='C-GFA' ~ 'GFAP',
                                        sos_flexlab$Kortnamn=='P-NFL' ~ 'Plasma NfL'
)

sos_flexlab <- sos_flexlab %>% 
  mutate(res.num.adj=case_when(
    Kortnamn=='qAm40xE' ~ res.num*conv.40,
    Kortnamn=='Am42xE' ~ res.num*conv.42,
    T ~res.num))


#go wide again, this time grouped by biomarker
bm.wide2<-sos_flexlab %>% 
  filter(!is.na(biomarker)) %>% 
  group_by(LopNr, biomarker, Provdatum) %>% 
  mutate(seqnum=row_number()) %>%
  pivot_wider(id_cols=c('LopNr', 'Provdatum', 'seqnum'), names_from=biomarker, values_from=res.num.adj) %>%
  #Create ratios
  mutate(Amy42.40=Amy42/Amy40,
         Amy42.pTau=Amy42/pTau)

bm.long2<-bm.wide2 %>% 
  pivot_longer(-c('LopNr', 'Provdatum', 'seqnum'), names_to='biomarker', values_to='res.num.adj') %>% 
  filter(!is.na(res.num.adj))


#Calculate youden index and optimal cutpoint
yoden<-bm.wide2 %>% 
  ungroup %>% 
  dplyr::select(LopNr, Amy42.40, Amy42.pTau) %>% 
  ungroup %>% drop_na %>%
  mutate(Amy42.40.path=as.numeric(Amy42.40<0.072))

cp <- cutpointr(data=yoden, x=Amy42.pTau, class=Amy42.40.path, 
                method = oc_youden_kernel,
                pos_class=1, direction='<=')

#Mark pathological
bm.long2 <-bm.long2 %>%
  mutate(pathological=case_when(
    biomarker=='Amy42.40' & res.num.adj<0.072~ TRUE,
    biomarker=='Amy42.40' & res.num.adj>=0.072 ~ FALSE,
    biomarker=='Amy42.pTau' & res.num.adj<cp$optimal_cutpoint ~ TRUE,
    biomarker=='Amy42.pTau' & res.num.adj>=cp$optimal_cutpoint ~ FALSE))


bm.path<-bm.long2 %>% 
  filter(biomarker=="Amy42.40"|biomarker=="Amy42.pTau") %>% 
  filter(!is.na(pathological)) %>% 
  ungroup() %>% 
  dplyr::select(-c(res.num.adj,seqnum,Provdatum)) %>% 
  mutate(pathological=ifelse(pathological=="TRUE",1,0)) %>% 
  distinct_all() %>% 
  group_by(LopNr,biomarker) %>% 
  mutate(amy_pos=max(pathological)) %>% 
  filter(pathological==amy_pos) %>% 
  ungroup() %>% 
  dplyr::select(-amy_pos) %>% 
  pivot_wider(id_cols='LopNr', names_from=biomarker, values_from=pathological)
  


