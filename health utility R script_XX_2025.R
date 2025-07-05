library(readxl)
library(tidyverse)



# Read in the data
utility_df_raw<-read_excel("literature review of utility.xlsx", sheet = "included studies"
                           ,skip=2, col_names=c('Study', 'Timepoint','Country', 
                                                'Study_design','Study_population','Setting',
                                                'n_total_self','n_mci_self','mean_mci_self',
                                                'n_mild_self','mean_mild_self',
                                                'n_mod_self','mean_mod_self',
                                                'n_sev_self','mean_sev_self',
                                                'n_tot_proxy','n_mci_proxy', 'mean_mci_proxy',
                                                'n_mild_proxy','mean_mild_proxy',
                                                'n_mod_proxy','mean_mod_proxy',
                                                'n_sev_proxy','mean_sev_proxy')) %>%
  filter(!is.na(Study)) %>%
  mutate(row_num=row_number()) %>%
  filter(row_num<which(Study=="Glossary of terms")[1]) %>%
  select(-row_num)


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

utility_df_long <- utility_df_raw %>%
  pivot_longer(cols = n_mci_self:mean_sev_proxy, 
               names_to=c(".value", "group", "rater"), 
               names_pattern = "(.*)_(.*)_(.*)") %>%
  mutate(parse_mean_sd(mean)) %>% 
  filter(!is.na(mean))


#meta/analytic estimate
library(meta)
estimate_sd <- function(data) {
  
  metares <- metamean(n, mean, sd, Study, 
                      data = data, random = T)
  data.frame(n=sum(metares$n), mean=metares$TE.random, sd=metares$seTE.random,
             lb=metares$lower.random,ub=metares$upper.random, i2=metares$I2)
}


table(utility_df_long$Setting,utility_df_long$group)



# compile utility by states and settings
state_utility<-estimate_sd(data=utility_df_long %>% 
                             filter(group=="mci",
                                    rater=="self",
                                    Setting %in% c("Comm","Mixed"))) %>% 
  mutate(state="MCI") %>% 
  
  rbind(estimate_sd(data=utility_df_long %>% 
                  filter(group=="mild",
                         rater=="proxy",
                         Setting %in% c("Comm","Mixed"))) %>% 
          mutate(state="Mild")) %>% 
  
  rbind(estimate_sd(data=utility_df_long %>% 
                          filter(group=="mild",
                                 rater=="proxy",
                                 Setting %in% c("Institute","Mixed"))) %>% 
          mutate(state="Mild_inst")) %>% 
  
  rbind(estimate_sd(data=utility_df_long %>% 
                      filter(group=="mod",
                             rater=="proxy",
                             Setting %in% c("Comm","Mixed"))) %>% 
          mutate(state="Moderate")) %>% 
  
  rbind(estimate_sd(data=utility_df_long %>% 
                      filter(group=="mod",
                             rater=="proxy",
                             Setting %in% c("Institute","Mixed"))) %>% 
          mutate(state="Moderate_inst")) %>% 
  
  rbind(estimate_sd(data=utility_df_long %>%
                      filter(group=="sev",
                             rater=="proxy",
                             Setting %in% c("Comm","Mixed"))) %>% 
          mutate(state="Severe")) %>% 
  
  rbind(estimate_sd(data=utility_df_long %>% 
                      filter(group=="sev",
                             rater=="proxy",
                             Setting %in% c("Institute","Mixed"))) %>% 
          mutate(state="Severe_inst"))


state_utility<-state_utility %>% 
  rbind(data=utility_df_long %>%
          filter(group=="mci",
                 rater=="self",
                 Setting %in% c("Institute","Mixed")) %>% 
          mutate(state="MCI_inst",
                 lb=mean-1.98*(sd/sqrt(n)),
                 ub=mean+1.98*(sd/sqrt(n)),
                 i2=NA) %>% 
          select(n,mean,sd,lb,ub,i2,state))




state_utility<-state_utility %>% 
  mutate(state=factor(state, levels=c("MCI","Mild","Moderate","Severe",
                                      "MCI_inst","Mild_inst","Moderate_inst","Severe_inst"),
                      labels=c("MCI","Mild","Moderate","Severe",
                               "MCI_inst","Mild_inst","Moderate_inst","Severe_inst"))) %>% 
  arrange(state)



n_study<-c("MCI",n_distinct(utility_df_long %>% 
                       filter(group=="mci",
                              rater=="self",
                              Setting %in% c("Comm","Mixed")) %>% 
                      select(Study))) %>% 
  
  rbind(c("MCI_inst",n_distinct(utility_df_long %>% 
                                   filter(group=="mci",
                                          rater=="self",
                                          Setting %in% c("Institute","Mixed")) %>% 
                                   select(Study)))) %>% 
  
  rbind(c("Mild",n_distinct(utility_df_long %>% 
                      filter(group=="mild",
                             rater=="proxy",
                             Setting %in% c("Comm","Mixed")) %>% 
                     select(Study)))) %>% 
  
  rbind(c("Mild_inst",n_distinct(utility_df_long %>% 
                      filter(group=="mild",
                             rater=="proxy",
                             Setting %in% c("Institute","Mixed")) %>% 
                     select(Study)))) %>% 
  
  rbind(c("Moderate",n_distinct(utility_df_long %>% 
                      filter(group=="mod",
                             rater=="proxy",
                             Setting %in% c("Comm","Mixed")) %>% 
                     select(Study)))) %>% 
  
  rbind(c("Moderate_inst",n_distinct(utility_df_long %>% 
                      filter(group=="mod",
                             rater=="proxy",
                             Setting %in% c("Institute","Mixed")) %>% 
                     select(Study)))) %>% 
  
  rbind(c("Severe",n_distinct(utility_df_long %>%
                      filter(group=="sev",
                             rater=="proxy",
                             Setting %in% c("Comm","Mixed")) %>% 
                     select(Study)))) %>% 
  
  rbind(c("Severe_inst",n_distinct(utility_df_long %>% 
                      filter(group=="sev",
                             rater=="proxy",
                             Setting %in% c("Institute","Mixed")) %>% 
                      select(Study)))) %>% 
  as.data.frame() %>% 
  set_names(c("state", "n_study"))


state_utility<-state_utility %>% 
  left_join(n_study)


state_utility %>% 
  mutate(ci=paste0(format(round(mean, 2), nsmall=2), " (", 
                  format(round(lb, 2), nsmall=2), "-", 
                  format(round(ub, 2), nsmall=2), ")"),
         sd_r=format(round(sd, 3), nsmall=3),
         mean_r=format(round(mean, 2), nsmall=2))
