require(dplyr)
mpi.global.mpi <- function(period_eff_s, period_eff_first, p_techeff_1t, p_techeff_2t){
  
  period_merge <- 
    left_join(period_eff_s, period_eff_first) %>% 
    left_join(p_techeff_1t) %>% 
    left_join(p_techeff_2t)
  
  global_mpi <- 
    period_merge %>% 
    mutate(eff_second_t = eff_s_t/eff_first_t) %>% 
    mutate(sec_first_t = eff_first_t/tech_first_t,
           sec_second_t = eff_second_t/tech_second_t) %>% 
    group_by(idx.var) %>% 
    mutate(mpi_s = eff_s_t*(1/dplyr::lag(eff_s_t)),
           mpi_1 = eff_first_t*(1/dplyr::lag(eff_first_t)),
           mpi_2 = eff_second_t*(1/dplyr::lag(eff_second_t)),
           mpi_tech_1 = tech_first_t*(1/dplyr::lag(tech_first_t)),
           mpi_tech_2 = tech_second_t*(1/dplyr::lag(tech_second_t)),
           mpi_sec_1 = sec_first_t*(1/dplyr::lag(sec_first_t)),
           mpi_sec_2 = sec_second_t*(1/dplyr::lag(sec_second_t)))
  
  return(global_mpi)
  
}
