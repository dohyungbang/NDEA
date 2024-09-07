mpi.period.first.eff <- function(dea_data, epsilon, agg.sys.eff, period.sys.eff, idx.var, time.var, x.vars, z.vars, y.vars){
  
  n_x <- length(x.vars)
  n_y <- length(y.vars)
  n_z <- length(z.vars)
  
  target_df <- 
    dea_data %>% 
    select(idx.var, time.var, x.vars, z.vars, y.vars) %>% 
    rename(idx.var = idx.var,
           time.var = time.var)
  
  dmu_list <- unique(dea_data[,idx.var] %>% as.character())
  year_list <- unique(dea_data[,time.var] %>% as.numeric())
  
  agg_data <- 
    target_df %>% 
    select(-time.var) %>% 
    group_by(idx.var) %>% 
    summarise(across(everything(), list(sum)))
  
  names(agg_data) <- c("idx.var", x.vars, z.vars, y.vars)
  
  
  first_period_constr_mat <- c()
  second_period_constr_mat <- c()
  for (dmu_j in dmu_list){
    
    for (year_p in year_list){
      
      period_sub_target_y <- target_df %>% filter(idx.var == dmu_j & time.var == year_p) %>% select(y.vars) %>% as.numeric()
      period_sub_target_x <- target_df %>% filter(idx.var == dmu_j & time.var == year_p) %>% select(x.vars) %>% as.numeric()
      period_sub_target_z <- target_df %>% filter(idx.var == dmu_j & time.var == year_p) %>% select(z.vars) %>% as.numeric()
      
      first_period_constr_mat <- rbind(first_period_constr_mat, 
                                       c(-period_sub_target_x, period_sub_target_z, rep(0, n_y)) )
      
      second_period_constr_mat <- rbind(second_period_constr_mat, 
                                        c(rep(0, n_x), -period_sub_target_z, period_sub_target_y) )
    }
  }
  
  # non-zero constraint
  nonzero_constr_mat <- diag(1,  n_x + n_z + n_y)
  
  eff_1t_list <- c()
  for (dmu_k in dmu_list){
    
    agg_sub_target <- agg_data %>% filter(idx.var == dmu_k)
    agg_eff_s <- agg.sys.eff %>% filter(idx.var == dmu_k) %>% select(eff_s) %>% as.numeric()
    
    fixed_eff_s_constr_mat <- c(-agg_eff_s*(agg_sub_target[,x.vars] %>% as.numeric()), 
                                rep(0, n_z), 
                                agg_sub_target[,y.vars] %>% as.numeric())
    
    for (year_t in year_list){
      
      period_sub_target <- target_df %>% filter(idx.var == dmu_k & time.var == year_t)
      period_eff_s_t <- period.sys.eff %>% filter(idx.var == dmu_k & time.var == year_t) %>% select(eff_s_t) %>% as.numeric()

      obj_func <- c(rep(0, n_x),
                    period_sub_target[, c(z.vars)] %>% as.numeric(),
                    rep(0, n_y))
      
      linear_constr_mat <- c(period_sub_target[,x.vars] %>% as.numeric(), 
                             rep(0, n_z), 
                             rep(0, n_y))
      
      fixed_eff_s_t_constr_mat <- c(-period_eff_s_t*(period_sub_target[,x.vars] %>% as.numeric()), 
                                    rep(0, n_z), 
                                    period_sub_target[,y.vars] %>% as.numeric())
      
      constr_mat <- rbind(linear_constr_mat, # 1
                          fixed_eff_s_constr_mat, # 1
                          fixed_eff_s_t_constr_mat, # 1
                          first_period_constr_mat, # length of dmu * length of year
                          second_period_constr_mat, # length of dmu * length of year
                          nonzero_constr_mat)
      
      constr_dir <- c("=", 
                      "=", 
                      "=", 
                      rep("<=", length(dmu_list)*length(year_list)),
                      rep("<=", length(dmu_list)*length(year_list)),
                      rep(">=", nrow(nonzero_constr_mat))
      )
      
      constr_rhs <- c(1, 
                      0,
                      0,
                      rep(0, length(dmu_list)*length(year_list)),
                      rep(0, length(dmu_list)*length(year_list)),
                      rep(epsilon, nrow(nonzero_constr_mat))
      )
      
      mod <- lp("max",
                obj_func,
                constr_mat,
                constr_dir,
                constr_rhs)
      
      eff_first_t <- mod$objval
      
      result_list <- data.frame(idx.var = dmu_k,
                                time.var = year_t,
                                eff_first_t = eff_first_t)
      
      eff_1t_list <- rbind.data.frame(eff_1t_list, result_list)
      
    }
    
  }

  
  
  return(eff_1t_list)
  
}