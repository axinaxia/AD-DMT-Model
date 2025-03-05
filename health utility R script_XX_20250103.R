library(readxl)
library(tidyverse)


# Read in the data
utility_df_raw<-read_excel("Eq-5D utility lit review (22 Jan 2025) Yunfei.xlsx", skip=2, col_names=c('Study', 
                                                                                 'Timepoint',
                                                                                 'Country', 
                                                                                 'Study_design',
                                                                                 'Setting',
                                                                                 'n_total_self',
                                                                                 'n_mci_self',
                                                                                 'mean_mci_self',
                                                                                 'n_mild_self',
                                                                                 'mean_mild_self',
                                                                                 'n_mod_self',
                                                                                 'mean_mod_self',
                                                                                 'n_mildmod_self',
                                                                                 'mean_mildmod_self',
                                                                                 'n_severe_self',
                                                                                 'mean_sev_self',
                                                                                 'n_ns_self',
                                                                                 'mean_ns_self',
                                                                                 'n_tot_proxy',
                                                                                 'n_mci_proxy', 
                                                                                 'mean_mci_proxy',
                                                                                 'n_mild_proxy',
                                                                                 'mean_mild_proxy',
                                                                                 'n_mod_proxy',
                                                                                 'mean_mod_proxy',
                                                                                 'n_mildmod_proxy',
                                                                                 'mean_mildmod_proxy',
                                                                                 'n_sev_proxy',
                                                                                 'mean_sev_proxy',
                                                                                 'n_ns_proxy',
                                                                                 'mean_ns_proxy',
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

utility_df_long <- utility_df_raw %>%
  pivot_longer(cols = n_mci_self:mean_ns_proxy, 
               names_to=c(".value", "group", "rater"), 
               names_pattern = "(.*)_(.*)_(.*)"
               ) %>%
  filter(!is.na(n)) %>% # remove studies without sample sizes for each disease state
  mutate(parse_mean_sd(mean)) %>% 
  filter(!is.na(sd)) %>% # remove studies without standard deviations for each disease state
  filter(group!="ns") # remove studies without information on disease state
  

table((utility_df_long %>% 
         filter(rater=="self"|rater=="proxy",
                group %in% c("mci")))$Setting)

table((utility_df_long %>% 
         filter(rater=="proxy",
                group %in% c("mild")))$Setting)

table((utility_df_long %>% 
         filter(rater=="proxy",
                group %in% c("mod")))$Setting)

table((utility_df_long %>% 
         filter(rater=="proxy",
                group %in% c("sev")))$Setting)


#meta/analytic estimate
library(meta)
estimate_sd <- function(data) {
  
  metares <- metamean(n, mean, sd, Study, 
                      data = data)
  data.frame(n=sum(metares$n), mean=metares$TE.common, sd=metares$seTE.common,
             lb=metares$lower.common,ub=metares$upper.common)
}


table(utility_df_long$Setting,utility_df_long$group)


# compile utility by states and settings
state_utility<-estimate_sd(data=utility_df_long %>% 
                             filter(group=="mci",
                                    rater=="self"|rater=="proxy",
                                    Setting %in% c("Comm","Mixed"))) %>% 
  mutate(state="MCI") %>% 
  rbind(estimate_sd(data=utility_df_long %>% 
                      filter(group=="mci",
                             rater=="self"|rater=="proxy",
                             Setting %in% c("Institute","Mixed"))) %>% 
          mutate(state="MCI_inst")) %>% 
  
  rbind(estimate_sd(data=utility_df_long %>% 
                  filter(group=="mild",
                         rater=="proxy",
                         Setting %in% c("Comm","Gen pub","Mixed"))) %>% 
      mutate(state="Mild") %>% 
      rbind(estimate_sd(data=utility_df_long %>% 
                          filter(group=="mild",
                                 rater=="proxy",
                                 Setting %in% c("Institute","Mixed"))) %>% 
              mutate(state="Mild_inst"))) %>% 
  rbind(estimate_sd(data=utility_df_long %>% 
                  filter(group=="mod",
                         rater=="proxy",
                         Setting %in% c("Comm","Gen pub","Mixed"))) %>% 
      mutate(state="Moderate") %>% 
      rbind(estimate_sd(data=utility_df_long %>% 
                          filter(group=="mod",
                                 rater=="proxy",
                                 Setting %in% c("Institute","Mixed"))) %>% 
              mutate(state="Moderate_inst"))) %>% 
  rbind(estimate_sd(data=utility_df_long %>% 
                  filter(group=="sev",
                         rater=="proxy",
                         Setting %in% c("Comm","Gen pub","Mixed"))) %>% 
      mutate(state="Severe") %>% 
      rbind(estimate_sd(data=utility_df_long %>% 
                          filter(group=="sev",
                                 rater=="proxy",
                                 Setting %in% c("Institute","Mixed"))) %>% 
              mutate(state="Severe_inst"))
  )




state_utility<-state_utility %>% 
  mutate(state=factor(state, levels=c("MCI", "Mild", "Moderate", "Severe", 
                                      "MCI_inst", "Mild_inst", "Moderate_inst", "Severe_inst"),
                      labels=c("MCI", "Mild", "Moderate", "Severe",
                               "MCI_inst", "Mild_inst", "Moderate_inst", "Severe_inst"))) %>% 
  arrange(state)

