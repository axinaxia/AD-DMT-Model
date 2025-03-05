library(readxl)
library(tidyverse)


# Read in the data
# setwd('~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Projekt/CSFregister/ADmodel/')
setwd("C:/★Work files in KI/R environment/health utility")


df<-read_excel("Eq-5D utility lit review (4 Sep 2024).xlsx", skip=2, col_names=c('Study', 
                                                                    'Timepoint',
                                                                    'Country', 
                                                                    'Study_design',
                                                                    'Setting',
                                                                    'n_self',
                                                                    'mci_self',
                                                                    'mild_self',
                                                                    'mod_self',
                                                                    'mild_mod_self',
                                                                    'sev_self',
                                                                    'NS_self',
                                                                    'n_proxy',
                                                                    'mci_proxy',
                                                                    'mild_proxy',
                                                                    'mod_proxy',
                                                                    'mild_mod_proxy',
                                                                    'sev_proxy',
                                                                    'NS_proxy'
                                                                    ))

#parse means and stdev

parse_mean_sd<-function(x){
  r=data.frame(mean=rep(NA, length(x)),
               sd=rep(NA, length(x)))
  for (i in 1:length(x)) {
    split<-strsplit(x[i], '\\(') %>% unlist
    if (!is.na(split[1]) & !is.null(split[1]))
       r$mean[i]<-as.numeric(split[1])
    if (!is.na(split[2]) & !is.null(split[2]))
      split[2]<-str_replace(split[2], '\\)', '') 
      r$sd[i]<-as.numeric(split[2])
  }
    return(r)
}

df<-df %>%
  pivot_longer(cols=c( 'mci_self',
                   'mild_self',
                   'mod_self',
                   'mild_mod_self',
                   'sev_self',
                   'NS_self',
                   'mci_proxy',
                   'mild_proxy',
                   'mod_proxy',
                   'mild_mod_proxy',
                   'sev_proxy',
                   'NS_proxy'), names_to='group', values_to='value') %>%
  filter(!is.na(value)) %>%
  mutate(parse_mean_sd(value)) 

# df %>% group_by(group) %>%
#   select(mean.value=mean) %>%
#   drop_na %>%
#   summarise(mean=mean(mean.value))

# df %>% group_by(group) %>%
#   rename(mean.value=mean) %>%
#   # drop_na %>% 
#   summarise(mean = mean(mean.value, na.rm = TRUE),
#             sd = mean(sd, na.rm = TRUE),
#             n = sum(n_self, na.rm = TRUE) + sum(n_proxy, na.rm = TRUE)
#   )


#meta/analytic estimate
library(meta)
metares<-metamean(n_self, mean, sd, Study
                  , data=df %>% filter(group=='mild_self'))
forest(metares
       )




            