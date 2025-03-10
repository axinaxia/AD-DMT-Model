library(readxl)
library(tidyverse)


# Read in the data
# setwd('~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Projekt/CSFregister/ADmodel/')
setwd('C:/Users/xinxia/OneDrive - Karolinska Institutet/Postdoc project plan/Postdoc project - obligation/AD DMT HE model_XX/Statistical analyses')

df<-read_excel("Eq-5D utility lit review (18 Sep 2024).xlsx", skip=2, col_names=c('Study', 
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
                                                                                 'n_proxy_tot',
                                                                                 'n_mci', 
                                                                                 'mean_mci',
                                                                                 'n_mild',
                                                                                 'mean_mild',
                                                                                 'n_mod',
                                                                                 'mean_mod',
                                                                                 'n_mildmod',
                                                                                 'mean_mildmod',
                                                                                 'n_sev',
                                                                                 'mean_sev',
                                                                                 'n_ns',
                                                                                 'mean_ns',
                                                                                 'comments'
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

df <- df %>%
  pivot_longer(cols = n_mci:mean_ns, 
               names_to=c(".value", "group"), 
               names_pattern = "(.*)_(.*)"
               ) %>%
  filter(!is.na(n)) %>%
  mutate(parse_mean_sd(mean)) %>% 
  filter(!is.na(sd)) 
  

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
estimate_sd <- function(level) {
  
  metares <- metamean(n, mean, sd, Study, 
                      data = df %>% filter(group==level))
  data.frame(n=sum(metares$n), mean=metares$TE.common, sd=metares$seTE.common * sqrt(sum(metares$n)))
}

estimate_sd("mci")
estimate_sd("mild")
estimate_sd("mod")
estimate_sd("mildmod")
estimate_sd("sev")
estimate_sd("ns")

group_level <- c("mci", "mild", "mod", "mildmod", "sev", "ns") 
group_name <- c("mci_proxy", "mild_proxy", "moderate_proxy", "mild_moderate_proxy", "severe_proxy", "Ns_proxy") 

group_level %>% 
  map_df(estimate_sd) %>% 
  mutate(seq = 1:6, 
         group = group_name, 
         mean = round(mean, 2),
         sd = round(sd, 2),
         .before = 1) %>% 
  set_names(c("", "group", "sample size", "mean", "SD"))
  



