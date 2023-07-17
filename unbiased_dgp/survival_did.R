library(tidyverse)
library(tictoc)
library(here)
library(DeclareDesign)
library(fixest)
source(here::here('unbiased_dgp', 'full_landscape.R'))

library(survival)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(msm)

survival_did <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, std_p = 0.0, cellsize, ppoints, cpoints){
  
  countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  pixloc <- pixloc_df
  
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  
  n_mod = 6
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'b3'= rep(b3, summ_row), 
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row),
                             'iteration' = rep(NA, summ_row), "model" = rep(NA, summ_row),
                             'HRR' = rep(NA, summ_row), "d_obs" = rep(NA, summ_row),
                             'bias'=rep(NA, summ_row), 'cover'=rep(NA, summ_row),'notes'=rep(NA, summ_row),
                             stringsAsFactors=FALSE)

for(i in 1:n){
  tic("loop")
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # simulating data ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Nobs <- length(pixloc$treat)  
  panels <- fabricate(
    pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
    year = add_level(N = (years*2), nest = FALSE),
    obs = cross_levels(
      by = join(pixels, year),
      post = ifelse(year > years, 1, 0),
      v_it = rnorm(N, 0, std_v),
      ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
      ystar_cf = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
    )
  )
  
  #generate random 
  errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
  
  panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
  
  panels <- panels %>%
    inner_join(pixloc, by = c("pixels", "treat")) %>%
    inner_join(errortable_property, by = "property") %>%
    #inner_join(errortable_county, by = "county") %>%
    mutate(ystar = ystar + p_err,
           ystar_cf = ystar_cf + p_err,
           y = (ystar > 0)*1 ,
           y_cf = (ystar_cf > 0)*1 )
  
  
  panels <- panels %>%
    mutate(pixels = as.numeric(pixels),
           year = as.numeric(year),
           defor_indic = ifelse(y==1, year, 99),
           defor_indic_cf = ifelse(y_cf==1, year, 99))%>%
    group_by(pixels)%>%
    mutate(defor_year = min(defor_indic),
           defor_year_cf = min(defor_indic_cf),
           defor = (year>=defor_year)*1,
           y_it = ifelse(year>defor_year, NA, defor),
           defor_cf = (year>=defor_year_cf)*1,
           y_it_counter = ifelse(year>defor_year_cf, NA, defor_cf)
    )%>%
    ungroup()
  
  DID <- feols(y_it ~  post*treat, data = panels) 
  DID_se <- summary(DID, se = "hetero")$se[4]
  DID_cover <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_se, DID$coefficients[4] + 1.96 * DID_se)*1
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Create counterfactual panel dataframe ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # this dataframe only contains treated units, both with and without treatment
  
  panels_long_cf <- subset(panels, treat == 1)%>%
    gather(key = "group", value = "outcome", c(y_it, y_it_counter))%>%
    mutate(observed = ifelse(group == "y_it", 1, 0)) # observed = 1 means that treatment did actually happen
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Set up survival dataframe ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Ref: https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
  surv_df <- panels %>% 
    mutate(t_start = year - 1,
           t_end = year,
           # t_end = ifelse((t_end==20) & (y_it==0), Inf, t_end),
           outcome = y_it,
           treat_now = treat * post) %>% 
    select(pixels, t_start, t_end, outcome, treat, treat_now, post, year) %>% 
    drop_na()
  
  
  ### counterfactual survival dataframe
  surv_cf_df <- panels_long_cf %>% 
    mutate(t_start = year - 1,
           t_end = year,
           # t_end = ifelse((t_end==20) & (y_it==0), Inf, t_end),
           outcome = outcome,
           treat_now = observed * post) %>% 
    select(pixels, t_start, t_end, outcome, treat, group, treat_now, post, year, observed) %>% 
    drop_na()
  
  ### corrected dataframe so we can include post dummy 
  surv_df_correction <- surv_df %>% 
    mutate(t_start = ifelse(post==1, t_start -10, t_start),
           t_end = ifelse(post==1, t_end -10, t_end))
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Run Cox proportional hazards models ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # running main cox model ----
  
  # km_AG_fit <- survfit(Surv(tstart, tend, outcome) ~ treat + treat_now, data=surv_df)
  # autoplot(km_AG_fit)
  cox <- coxph(Surv(t_start, t_end, outcome) ~ treat + treat_now 
               , data = surv_df)
  
  
  #ggforest(cox)
  # print(summary(cox))
  
  coefs <- cox$coefficients
  hr <- coefs[[2]] %>% exp()
  vcoeff <- cox$var[2,2]
  se_cox <- deltamethod(~ exp(x1), coefs[[2]], vcoeff)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # counterfactual cox proportional hazards model 
  
  # we use the counterfactual survival dataframe, which only includes the treated observations and either their treated or untreated outcomes
  cox_counterfactual <- coxph(Surv(t_start, t_end, outcome) ~ treat_now 
                              , data = surv_cf_df)
  
  #ggforest(cox)
  #print(summary(cox_counterfactual))
  coefs <- cox_counterfactual$coefficients
  hr_cf <- coefs[[1]] %>% exp()
  
  vcoeff <- cox_counterfactual$var
  se_cox_cf <- deltamethod(~ exp(x1), coefs[[1]], vcoeff)
  
  # this hr corresponds really well to what we'd expect if the hazard ratio we want is 
  haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # running main cox model ----
  
  ## this model is more similar to the traditional DID setup but does not actually recover the desired HRR
  cox_interaction <- coxph(Surv(t_start, t_end, outcome) ~ post*treat
                           , data = surv_df_correction )
  #summary(cox_interaction)
  coefs <- cox_interaction$coefficients
  hr_int <- coefs[[3]] %>% exp()
  
  vcoeff <- cox_interaction$var[3,3]
  se_cox_int <- deltamethod(~ exp(x1), coefs[[3]], vcoeff)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### now we'll try to estimate separately the components of the desired hazard ratio
  
  cox_post <- coxph(Surv(t_start, t_end, outcome) ~ treat
                    , data = surv_df_correction %>% filter(post==1))
  #summary(cox_post)
  hr_11_01 <- cox_post$coefficients %>% exp()
  
  cox_treat <- coxph(Surv(t_start, t_end, outcome) ~ post
                     , data = surv_df_correction %>% filter(treat==1))
  #summary(cox_treat)
  hr_11_10 <- cox_treat$coefficients %>% exp()
  
  cox_notreat <- coxph(Surv(t_start, t_end, outcome) ~ post
                       , data = surv_df_correction %>% filter(treat==0))
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
  
  hr_2 <- hr_11_10/hr_01_00
  
  vcov_hr <- matrix(0, nrow = 2, ncol = 2)
  vcov_hr[1,1] <- cox_treat$var
  vcov_hr[2,2] <- cox_notreat$var
  se_hr_rat <- deltamethod(
    ~ exp(x1) / exp(x2) ,
    c(cox_treat$coefficients, cox_notreat$coefficients),
    vcov_hr
  )
  
  #### calculating observed deforestation rate to transition to ATT estimate
  defor_summary <- panels %>%
    group_by(treat, post) %>%
    summarise(mean_y_it = mean(y_it, na.rm = TRUE),
              mean_y_it_counter = mean(y_it_counter, na.rm = TRUE))%>%
    ungroup()
  
  d_obs = defor_summary[4,3]$mean_y_it
  
  
  ATT_11_cf <- d_obs - d_obs/hr_11_cf
  ATT_cox <- d_obs - d_obs/hr
  ATT_stratcox <- d_obs - d_obs/hr_2
  ATT_cf <- d_obs - d_obs/hr_cf
  ATT_int <- d_obs - d_obs/hr_int
  
  hr_11_cf_cover <- dplyr::between(haz_rat, hr_11_cf - 1.96 * se_hr_11_cf, hr_11_cf + 1.96 * se_hr_11_cf)*1
  
  hr_cover <- dplyr::between(haz_rat, hr - 1.96 * se_cox, hr + 1.96 * se_cox)*1
  
  hr_int_cover <- dplyr::between(haz_rat, hr_int - 1.96 * se_cox_int, hr_int + 1.96 * se_cox_int)*1
  
  hr_rat_cover <- dplyr::between(haz_rat, hr_2 - 1.96 * se_hr_rat, hr_2 + 1.96 * se_hr_rat)*1
  
  
  
  
  firstcol = which(colnames(summary_long)=="iteration")
  lastcol = which(colnames(summary_long)=="cover")
  
  summary_long[i,c(firstcol:lastcol)] <- c(
    i,
    "cox_counterfactual",#model
    hr_cf, #HRR
    d_obs, #d_obs
    ATT_cf - ATT, #bias
    NA#cover
  )
  
  summary_long[i+n,c(firstcol:lastcol)] <- c(
    i,
    "cox",#model
    hr, #HRR
    d_obs, #d_obs
    ATT_cox - ATT, #bias
    hr_cover#cover
  )
  
  summary_long[i+n*2,c(firstcol:lastcol)] <- c(
    i,
    "cox interaction",#model
    hr_int, #HRR
    d_obs, #d_obs
    ATT_int - ATT, #bias
    hr_int_cover#cover
  )
  
  summary_long[i+n*3,c(firstcol:lastcol)] <- c(
    i,
    "cox did",#model
    hr_11_cf, #HRR
    d_obs, #d_obs
    ATT_11_cf - ATT, #bias
    hr_11_cf_cover#cover
  )
  
  summary_long[i+n*4,c(firstcol:lastcol)] <- c(
    i,
    "cox HR1/HR2",#model
    hr_2, #HRR
    d_obs, #d_obs
    ATT_stratcox - ATT, #bias
    hr_rat_cover#cover
  )
  
  summary_long[i+n*5,c(firstcol:lastcol)] <- c(
    i,
    "pixel DID",#model
    NA, #HRR
    d_obs, #d_obs
    DID$coefficients[4] - ATT, #bias
    DID_cover#cover
  )
  
  print(i)
  toc()
  
}
  
  outputs = list(
    "summary_long" = summary_long
    )
  return(outputs)
  
  
}



