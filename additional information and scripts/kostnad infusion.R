#Set Options----
Sys.setenv(LANG = "eng")
options(dplyr.summarise.inform=FALSE)
options(scipen = 1)

#Load Libraries----

library(tidyverse)
library(haven)
library(labelled)
library(readxl)
#library(plyr)
library(tidyr)
library(haven)
library(flextable)
library(eq5d)
library(ggplot2)
library(robustbase)
#library(summarytools)
library(ggplot2)


# LOAD DATA-----
path<-'~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Projekt/CSF-registersamkörning/Kostnadsdata infusion mm/Data KPP 231002/'
setwd(path)
d<-read_excel('immunoterapier vid Alzheimers sjukdom.xlsx')


d %>% filter(!is.na(`Kostnad för läkemedel`))  %>% summarize(n())

d %>% filter(is.na(`Kostnad för läkemedel`)) %>% summarize(n())
d2<-d %>% filter(!is.na(`Kostnad för läkemedel`)) %>%
  mutate(kostnad=ifelse(is.na(`Grundkostnad mott.`),0,`Grundkostnad mott.`)+
                         ifelse(is.na(`Kostnad för behandlade personal`),0,`Kostnad för behandlade personal`)+
                         ifelse(is.na(`Kostnad för material` ),0,`Kostnad för material`)     
  ) %>% filter(kostnad>0)

ggplot(aes(x=kostnad),data=d2)+geom_histogram(fill='blue', bins=50,aes(y=after_stat(density)))+geom_density(adjust=4)
table(d2$kostnad)

d2 %>% summarize(mean(kostnad),
                  median(kostnad),
                  min(kostnad),
                  max(kostnad),
                  sd(kostnad))
  

  
 # select(År, kostnad, Region, mvotxt, sjkh)
                        