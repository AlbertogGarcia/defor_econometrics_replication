library(fixest)
library(survival)
library(msm)

source(here::here('unbiased_dgp', 'pixloc_DGP.R'))

all_specification_run <- function(nobs, years, b0, b1, b2_0, b2_1, b3, std_a = std_a, std_v = std_v, std_p = std_p, ppoints, cpoints,
                                  cellsize_list, pixloc){
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### set up results table
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  n_models <- 8 + length(cellsize_list)*2
  
  results <- data.frame('b0'= rep(b0, n_models), 'b1'= rep(b1, n_models), 'b2_0'= rep(b2_0, n_models), 'b2_1'= rep(b2_1, n_models), 'b3'= rep(b3, n_models), 
                             'std_a'= rep(std_a, n_models), 'std_v'= rep(std_v, n_models), 'std_p'= rep(std_p, n_models), "years" = years,
                             'gridsize' = rep(NA, n_models), 
                             'pixel'=rep(NA, n_models),'grid'=rep(NA, n_models),'property'=rep(NA, n_models),'county'=rep(NA, n_models),
                             'pixel fe'=rep(NA, n_models),'grid fe'=rep(NA, n_models),'property fe'=rep(NA, n_models),'county fe'=rep(NA, n_models),'treatment fe'=rep(NA, n_models),
                             'weights'=rep(NA, n_models), 'cox'=rep(NA, n_models), 'HE estimator'=rep(NA, n_models), 
                             'se_pixel'=rep(NA, n_models), 'se_grid'=rep(NA, n_models), 'se_property'=rep(NA, n_models), 'se_county'=rep(NA, n_models),
                             'bias'=rep(NA, n_models), 'cover'=rep(NA, n_models), 'power'=rep(NA, n_models), 
                             'notes'=rep(NA, n_models),
                             stringsAsFactors=FALSE)
  
  firstcol = which(colnames(results)=="gridsize")
  lastcol = which(colnames(results)=="power")
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Generate data based on DGP
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  dataframe_sim <- pixloc_DGP(nobs, years, b0, b1, b2_0, b2_1, b3, std_a = std_a, std_v = std_v, std_p = std_p, ppoints, cpoints,
                         cellsize_list, pixloc)
  
  countylevel_df = dataframe_sim$countylevel_df
  proplevel_df = dataframe_sim$proplevel_df
  grid_df_list = dataframe_sim$grid_df_list
  panels = dataframe_sim$panels
  surv_df = dataframe_sim$surv_df
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Pixel TWFE
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  pix_twfe <- suppressMessages(feols(y_it ~  post*treat|year + pixels, data = panels))
  
  estimate <- tail(pix_twfe$coefficients, n = 1)
  bias <- estimate - ATT
  se <- tail(summary(pix_twfe, cluster = ~pixels)$se, n=1)
  cover <- dplyr::between(ATT, estimate - 1.96 * se, estimate + 1.96 * se)*1
  power <- dplyr::between(0, estimate - 1.96 * se, estimate + 1.96 * se)*1
  
  results[1, c(firstcol:lastcol)] <- c(
    NA,
    1,0,0,0, # unit of analysis:pixel, grid, property, county
    1,0,0,0,0, # FE: pixel, grid, property, county, treatment
    0,0,0, # weights, cox, HE estimator
    1,0,0,0, # se cluster level: pixel, grid, property, county
    bias, cover, power
  )
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### DID
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  DID <- suppressMessages(feols(y_it ~  post*treat, data = panels))
  
  estimate <- tail(DID$coefficients, n = 1)
  bias <- estimate - ATT
  se <- tail(summary(DID, cluster = ~pixels)$se, n=1)
  cover <- dplyr::between(ATT, estimate - 1.96 * se, estimate + 1.96 * se)*1
  power <- dplyr::between(0, estimate - 1.96 * se, estimate + 1.96 * se)*1
  
  results[2, c(firstcol:lastcol)] <- c(
    NA,
    1,0,0,0, # unit of analysis:pixel, grid, property, county
    0,0,0,0,1, # FE: pixel, grid, property, county, treatment
    0,0,0, # weights, cox, HE estimator
    1,0,0,0, # se cluster level: pixel, grid, property, county
    bias, cover, power
  )
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Aggregated FE
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # County fe model
  
  countyfe_DID <- suppressMessages(feols(y_it ~  post:treat|year + county, data = panels))
  
  estimate <- tail(countyfe_DID$coefficients, n = 1)
  bias <- estimate - ATT
  se <- tail(summary(countyfe_DID, cluster = ~county)$se, n=1)
  cover <- dplyr::between(ATT, estimate - 1.96 * se, estimate + 1.96 * se)*1
  power <- dplyr::between(0, estimate - 1.96 * se, estimate + 1.96 * se)*1
  
  results[3, c(firstcol:lastcol)] <- c(
    NA,
    1,0,0,0, # unit of analysis:pixel, grid, property, county
    0,0,0,1,0, # FE: pixel, grid, property, county, treatment
    0,0,0, # weights, cox, HE estimator
    0,0,0,1, # se cluster level: pixel, grid, property, county
    bias, cover, power
  )
  
  # Property fe model
  proptreatfe_DID <- suppressMessages(feols(y_it ~  post:treat|year + treat + property, data = panels))
  
  estimate <- tail(proptreatfe_DID$coefficients, n = 1)
  bias <- estimate - ATT
  se <- tail(summary(proptreatfe_DID, cluster = ~property)$se, n=1)
  cover <- dplyr::between(ATT, estimate - 1.96 * se, estimate + 1.96 * se)*1
  power <- dplyr::between(0, estimate - 1.96 * se, estimate + 1.96 * se)*1
  
  results[4, c(firstcol:lastcol)] <- c(
    NA,
    1,0,0,0, # unit of analysis:pixel, grid, property, county
    0,0,1,0,1, # FE: pixel, grid, property, county, treatment
    0,0,0, # weights, cox, HE estimator
    0,0,1,0, # se cluster level: pixel, grid, property, county
    bias, cover, power
  )
  
  
  for(k in cellsize_list){
    
    pos <- as.numeric(which(cellsize_list == k))
    
    this_grid_panels <- panels %>%
      rename(this_grid = paste0("grid_", k))
    
    gridfe_DID <- suppressMessages(feols(y_it ~  post:treat|year + treat + this_grid, data = this_grid_panels))
    
    estimate <- tail(gridfe_DID$coefficients, n = 1)
    bias <- estimate - ATT
    se <- tail(summary(gridfe_DID, cluster = ~this_grid)$se, n=1)
    cover <- dplyr::between(ATT, estimate - 1.96 * se, estimate + 1.96 * se)*1
    power <- dplyr::between(0, estimate - 1.96 * se, estimate + 1.96 * se)*1
    
    results[4 + pos, c(firstcol:lastcol)] <- c(
      k,
      1,0,0,0, # unit of analysis:pixel, grid, property, county
      0,1,0,0,1, # FE: pixel, grid, property, county, treatment
      0,0,0, # weights, cox, HE estimator
      0,1,0,0, # se cluster level: pixel, grid, property, county
      bias, cover, power
    )
    
  }
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Aggregated unit of analysis
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Property UOA model
  prop_DID <- feols(deforrate ~  post:treat|year+property, data = proplevel_df)
  
  estimate <- tail(prop_DID$coefficients, n = 1)
  bias <- estimate - ATT
  se <- tail(summary(prop_DID, cluster = ~property)$se, n=1)
  cover <- dplyr::between(ATT, estimate - 1.96 * se, estimate + 1.96 * se)*1
  power <- dplyr::between(0, estimate - 1.96 * se, estimate + 1.96 * se)*1
  
  results[5 + length(cellsize_list), c(firstcol:lastcol)] <- c(
    NA,
    0,0,1,0, # unit of analysis:pixel, grid, property, county
    0,0,1,0,0, # FE: pixel, grid, property, county, treatment
    0,0,0, # weights, cox, HE estimator
    0,0,1,0, # se cluster level: pixel, grid, property, county
    bias, cover, power
  )
  
  # County UOA model
  county_DID <- feols(deforrate ~  post:treat|year+county, data = countylevel_df)
  
  estimate <- tail(county_DID$coefficients, n = 1)
  bias <- estimate - ATT
  se <- tail(summary(county_DID, cluster = ~county)$se, n=1)
  cover <- dplyr::between(ATT, estimate - 1.96 * se, estimate + 1.96 * se)*1
  power <- dplyr::between(0, estimate - 1.96 * se, estimate + 1.96 * se)*1
  
  
  results[6 + length(cellsize_list), c(firstcol:lastcol)] <- c(
    NA,
    0,0,0,1, # unit of analysis:pixel, grid, property, county
    0,0,0,1,0, # FE: pixel, grid, property, county, treatment
    0,0,0, # weights, cox, HE estimator
    0,0,0,1, # se cluster level: pixel, grid, property, county
    bias, cover, power
  )
  
  
  
  # Grid unit of analysis
  grid_DID <- purrr::map(grid_df_list, ~feols(deforrate ~ post:treat | year + treat + this_grid, data = .))
  
  grid_estimate = purrr::map_dbl(grid_DID, ~ tail(.$coefficients, n = 1))
  grid_bias = grid_estimate - ATT
  grid_se = purrr::map_dbl(grid_DID, ~tail(summary(., cluster = ~ this_grid)$se, n=1))
  grid_cover = purrr::map2_dbl(grid_estimate, grid_se, ~ dplyr::between(ATT, .x - 1.96 * .y, .x + 1.96 * .y)*1)
  grid_power = purrr::map2_dbl(grid_estimate, grid_se, ~ dplyr::between(0, .x - 1.96 * .y, .x + 1.96 * .y)*1)
  
  for(k in cellsize_list){
    
    pos <- as.numeric(which(cellsize_list == k))
    
    results[6 + length(cellsize_list) + pos, c(firstcol:lastcol)] <- c(
      k,
      0,1,0,0, # unit of analysis:pixel, grid, property, county
      0,1,0,0,0, # FE: pixel, grid, property, county, treatment
      0,0,0, # weights, cox, HE estimator
      0,1,0,0, # se cluster level: pixel, grid, property, county
      grid_bias[pos], grid_cover[pos], grid_power[pos]
    )
    
  }
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Survival analysis
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  ## this model is more similar to the traditional DID setup but does not actually recover the desired HRR
  cox_did <- coxph(Surv(t_start, t_end, outcome) ~ post*treat
                   , data = surv_df )
  #summary(cox_interaction)
  coefs <- cox_did$coefficients
  hr_did <- coefs[[3]] %>% exp()
  
  vcoeff <- cox_did$var[3,3]
  se_cox_did <- deltamethod(~ exp(x1), coefs[[3]], vcoeff)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### now we'll try to estimate separately the components of the desired hazard ratio
  
  cox_post <- coxph(Surv(t_start, t_end, outcome) ~ treat
                    , data = surv_df %>% filter(post==1))
  #summary(cox_post)
  hr_11_01 <- cox_post$coefficients %>% exp()
  
  cox_treat <- coxph(Surv(t_start, t_end, outcome) ~ post
                     , data = surv_df %>% filter(treat==1))
  #summary(cox_treat)
  hr_11_10 <- cox_treat$coefficients %>% exp()
  
  cox_notreat <- coxph(Surv(t_start, t_end, outcome) ~ post
                       , data = surv_df %>% filter(treat==0))
  #summary(cox_notreat)
  hr_01_00 <- cox_notreat$coefficients %>% exp()
  
  hr_11_cf <- 1/(1/hr_11_10 + 1/hr_11_01 - (1/(hr_11_01*hr_01_00)))
  
  vcov_hr <- matrix(0, nrow = 3, ncol = 3)
  vcov_hr[1,1] <- cox_treat$var
  vcov_hr[2,2] <- cox_post$var
  vcov_hr[3,3] <- cox_notreat$var
  
  se_hr_11_cf <- deltamethod(
    ~ 1/(1/exp(x1) + 1/exp(x2) - (1/(exp(x2)*exp(x3)))),
    c(cox_treat$coefficients, cox_post$coefficients, cox_notreat$coefficients),
    vcov_hr
  )
  
  #### calculating observed deforestation rate to transition to ATT estimate
  defor_summary <- panels %>%
    group_by(treat, post) %>%
    summarise(mean_y_it = mean(y_it, na.rm = TRUE))%>%
    ungroup()
  
  d_obs = defor_summary[4,3]$mean_y_it
  
  ATT_11_cf <- d_obs - d_obs/hr_11_cf
  ATT_cox <- d_obs - d_obs/hr_did
  
  
  cox_cover <- dplyr::between(haz_rat, hr_did - 1.96 * se_cox_did, hr_did + 1.96 * se_cox_did)*1
  cox_power <- dplyr::between(0, hr_did - 1.96 * se_cox_did, hr_did + 1.96 * se_cox_did)*1
  cox_bias <- ATT_cox - ATT
  
  results[7 + 2*length(cellsize_list), c(firstcol:lastcol)] <- c(
    NA,
    1,0,0,0, # unit of analysis:pixel, grid, property, county
    0,0,0,0,1, # FE: pixel, grid, property, county, treatment
    0,1,0, # weights, cox, HE estimator
    1,0,0,0, # se cluster level: pixel, grid, property, county
    cox_bias, cox_cover, cox_power
  )
  
  hr_11_cf_bias <- ATT_11_cf- ATT
  hr_11_cf_cover <- dplyr::between(haz_rat, hr_11_cf - 1.96 * se_hr_11_cf, hr_11_cf + 1.96 * se_hr_11_cf)*1
  hr_11_cf_power <- dplyr::between(0, hr_11_cf - 1.96 * se_hr_11_cf, hr_11_cf + 1.96 * se_hr_11_cf)*1
  
  
  results[8 + 2*length(cellsize_list), c(firstcol:lastcol)] <- c(
    NA,
    1,0,0,0, # unit of analysis:pixel, grid, property, county
    0,0,0,0,0, # FE: pixel, grid, property, county, treatment
    0,1,1, # weights, cox, HE estimator
    1,0,0,0, # se cluster level: pixel, grid, property, county
    hr_11_cf_bias, hr_11_cf_cover, hr_11_cf_power
  )
  
  return(results)
}
