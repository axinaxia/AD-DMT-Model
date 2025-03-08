n_sim<-10000 # add a sequence indicating number of simulations
sim_prof<-expand.grid(sim=1:n_sim,prof=1:nrow(pop_perc_func(12)))

ncpus = parallel::detectCores()-1
cl = makeCluster(ncpus, type="PSOCK")
clusterEvalQ(cl, c(library(tidyverse),library(heemod),
                   library(flexsurv),library(parallel)))
clusterExport(cl, c("psa_func",
                    ls(pattern = c("utiliti|model_trans|cycle|_sd|pop_perc")),
                    "mci_ad_model","mci_death_model","msm",
                    "transmat","timepoint_seq","profile_perc",
                    "boots_gethaz_mci","boots_gethaz_msm",
                    "get_haz","statenames","getcost",
                    "mci_svedem_costs","costwide","sim_prof"))



sim_psa<-parLapply(cl,1:nrow(sim_prof),function(n) {
  pop_perc<-pop_perc_rx3yr
  
  j<-sim_prof$sim[n]
  i<-sim_prof$prof[n]
  
  sim.seed<-2025 + j
  rx_cycles_i<-pop_perc[i,]$rx_cycles
  age_group_i<-pop_perc[i,]$age_group
  age_i<-pop_perc[i,]$age
  sex_i<-pop_perc[i,]$sex
  apoe_i<-pop_perc[i,]$apoe
  stage_i<-pop_perc[i,]$stage
  perc<-pop_perc[i,]$pop_perc
  
  psa_summary<-tryCatch({ # catch error message instead of stopping interations
    psa_func(sim.seed,
             rx_cycles_psa=rx_cycles_i, rx_rr_psa=0.69,
             rx_inst_psa=0, wane_fc_psa=1,
             age_init_psa=age_i, sex_psa=sex_i,
             apoegeno_psa=apoe_i, start_state_psa=stage_i,
             cost_mci_psa=getcost(mci_svedem_costs),
             cost_ad_psa=getcost(costwide),
             r_cost_psa=0.03,
             r_health_psa=0.03) %>% 
      cbind(rx_cycles=rx_cycles_i,
            age_group=age_group_i, age=age_i, 
            sex=sex_i, apoe=apoe_i, stage=stage_i, perc=perc,
            error_message=NA)
    
  }, error=function(e) {
    data.frame(rx_cycles=rx_cycles_i,
               age_group=age_group_i, age=age_i, 
               sex=sex_i, apoe=apoe_i, stage=stage_i, perc=perc,
               error_message=as.character(e$message)) 
  })
  return(psa_summary %>% cbind(sim=j))
}) %>% bind_rows()



stopCluster(cl)
