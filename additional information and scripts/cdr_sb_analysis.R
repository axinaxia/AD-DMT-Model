library(readr)
library(tidyr)
library(ggplot2)


CDR_SB_plot_data <- read_delim("additional information/CDR_SB_plot-data.csv", 
                               delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(CDR_SB_plot_data)


CDR_SB_plot_data<-CDR_SB_plot_data %>% 
  cbind(n_ar=c(859,875,801,792,765,757,677,649,642,603,554,524,299,265)) %>% 
  group_by(group) %>% 
  mutate(n_drop=n-lead(n),
         prob_part=1-n_drop/n) %>% 
  ungroup




# simulate a dataset with individual-level CDR scores
CDR_SB_plot_data <- CDR_SB_plot_data %>% 
  mutate(across(c(cdr_sb,cdr_sb_lb,cdr_sb_ub), ~ as.numeric(gsub(",", ".", .))),
         sd=(cdr_sb_ub-cdr_sb_lb)/(2*1.96))




# at time 0, CDR-SB for Lecanemab group: 3.17±1.34, for placebo group: 3.22±1.34
set.seed(2025)
n_sim <- 10000

cdr_lc_sim <- data.frame(
  id = 1:n_sim,
  cdr0 = rnorm(n_sim, mean = 3.17, sd = 1.34)  
)

for (t in seq(3, 18, by = 3)) {
  cdr_lc_sim <- cdr_lc_sim %>% 
    left_join(cdr_lc_sim %>%
                filter(!is.na(!!sym(paste0("cdr", t-3)))) %>% 
                sample_frac(
                  (CDR_SB_plot_data %>% filter(time == t - 3, group == 1))$prob_part
                ) %>% 
                mutate(
                  !!sym(paste0("cdr", t)) := cdr0 + rnorm(
                    n(),
                    mean = (CDR_SB_plot_data %>% filter(group == 1 & time == t))$cdr_sb,
                    sd = (CDR_SB_plot_data %>% filter(group == 1 & time == t))$sd
                  )
                ))
}






# placebo group
set.seed(2025)
cdr_pl_sim <- data.frame(
  id = 1:n_sim,
  cdr0 = rnorm(n_sim, mean = 3.22, sd = 1.34)
)

for (t in seq(3, 18, by = 3)) {
  cdr_pl_sim <- cdr_pl_sim %>% 
    left_join(cdr_pl_sim %>%
                filter(!is.na(!!sym(paste0("cdr", t-3)))) %>% 
                sample_frac(
                  (CDR_SB_plot_data %>% filter(time == t - 3, group == 0))$prob_part
                ) %>% 
                mutate(
                  !!sym(paste0("cdr", t)) := cdr0 + rnorm(
                    n(),
                    mean = (CDR_SB_plot_data %>% filter(group == 0 & time == t))$cdr_sb,
                    sd = (CDR_SB_plot_data %>% filter(group == 0 & time == t))$sd
                  )
                ))
}


round_to_half <- function(x) {
  round(x * 2) / 2
}

cdr_sim<-cdr_lc_sim %>% 
  mutate(group=1) %>% 
  bind_rows(cdr_pl_sim %>% 
              mutate(group=0)) %>% 
  mutate_at(vars(cdr0:cdr18), round_to_half) %>% 
  mutate(state=ifelse(cdr0<3,1,2), # CDR: MCI=[0,3), mild=[3-9.5), moderate=[9.5-16),severe=[16-18]
         prog3=ifelse((state==1&cdr3>=3)|(state==2&cdr3>=9.5),1,0),
         prog6=ifelse((state==1&cdr6>=3)|(state==2&cdr6>=9.5),1,0),
         prog9=ifelse((state==1&cdr9>=3)|(state==2&cdr9>=9.5),1,0),
         prog12=ifelse((state==1&cdr12>=3)|(state==2&cdr12>=9.5),1,0),
         prog15=ifelse((state==1&cdr15>=3)|(state==2&cdr15>=9.5),1,0),
         prog18=ifelse((state==1&cdr18>=3)|(state==2&cdr18>=9.5),1,0),
         prog=ifelse(prog3==1|prog6==1|prog9==1|prog12==1|prog15==1|prog18==1,1,0)) %>% 
  filter(cdr0>=0&cdr0<=8)

table(cdr_sim$group,cdr_sim$prog)





