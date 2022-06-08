#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level

library(ggplot2)
library(clubSandwich)
library(reshape2)
library(matrixStats)
library(ggplot2)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tidyverse)
library(tictoc)
library(fixest)
library(here)
library(DeclareDesign)
library(msm)
library(survival)

source(here::here('unbiased_dgp', 'multigrid_landscape.R'))
source(here::here('unbiased_dgp', 'nested_landscape.R'))
source(here::here('unbiased_dgp', 'proptreat_landscape.R'))

#begin function
specifications <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.5, std_p = 0.0, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE){
  
  cellsize_list <- list(cellsize_small, cellsize_med, cellsize_large)
  
  minpropsize = 0
  while(minpropsize < 10){
    
    if (proptreatassign) {
      countyscape = proptreat_landscape(nobs, cellsize_list, ppoints, cpoints)
      pixloc_df = countyscape$pixloc_df
    } else {
      
      if (nestedprops) {
        countyscape = nested_landscape(nobs, cellsize_list, ppoints, cpoints)
        pixloc_df = countyscape$pixloc_df
      
      } else {
        countyscape = multigrid_landscape(nobs, cellsize_list, ppoints, cpoints)
        pixloc_df = countyscape$pixloc_df
        
      }
      
    }
    
    minpropsize <- min(pixloc_df$parea)
    
    }
  
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  pixloc <- pixloc_df
  
  
  grid_covermat <- matrix(nrow = n, ncol = 3)
  grid_coeffmatrix <- matrix(nrow = n, ncol = 3)
  gridfe_covermat <- matrix(nrow = n, ncol = 3)
  gridfe_coeffmatrix <- matrix(nrow = n, ncol = 3)
  
  n_mod = 19
  # 3 grid sizes plus prop and county, aggregated 
  # and at fixed effects gives 10
  # weighted county and property models gives 2
  # pixel DID, pixel fe and pixel keeping pixels gives 3
  # two cox ph models gives 2
  
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'b3'= rep(b3, summ_row), 
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row), "years" = years,
                             'gridsize' = rep(NA, summ_row), 'iteration' = rep(NA, summ_row), 
                             'pixel'=rep(NA, summ_row),'grid'=rep(NA, summ_row),'property'=rep(NA, summ_row),'county'=rep(NA, summ_row),
                             'pixel fe'=rep(NA, summ_row),'grid fe'=rep(NA, summ_row),'property fe'=rep(NA, summ_row),'county fe'=rep(NA, summ_row),'treatment fe'=rep(NA, summ_row),
                             'weights'=rep(NA, summ_row), 'cox'=rep(NA, summ_row), 'HE estimator'=rep(NA, summ_row), 
                             'se_pixel'=rep(NA, summ_row), 'se_grid'=rep(NA, summ_row), 'se_property'=rep(NA, summ_row), 'se_county'=rep(NA, summ_row),
                             'bias'=rep(NA, summ_row), 'cover'=rep(NA, summ_row), 'prop_concern' = rep(0, summ_row), 'notes'=rep(NA, summ_row),
                             stringsAsFactors=FALSE)
  
  options(dplyr.summarise.inform = FALSE)
  print(min(pixloc$parea))
  
  for(i in 1:n){
    tic("loop")
    Nobs <- length(pixloc$treat)  
    panels <- fabricate(
      pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
      year = add_level(N = (years*2), nest = FALSE),
      obs = cross_levels(
        by = join(pixels, year),
        post = ifelse(as.numeric(year) > years, 1, 0),
        v_it = rnorm(N, 0, std_v),
        ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
        ystar_cf = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
      )
    )
    
    #generate random 
    errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    errortable_county <- data.frame(county = as.character(unique(pixloc$county)), c_err = rnorm(length(unique(pixloc$county)), 0, std_c))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(errortable_property, by = "property") %>%
      inner_join(errortable_county, by = "county") %>%
      mutate(year = as.numeric(year),
             ystar = ystar + p_err,
             ystar_cf = ystar_cf + p_err,
             y = (ystar > 0)*1 ,
             y_cf = (ystar_cf > 0)*1 )
    
    #need to determine which year deforestation occurred
    year_df <- panels %>%
      dplyr::select(pixels, year, y) %>%
      dcast(pixels ~ year , value.var = "y")
    
    rownames(year_df) <- year_df$pixels
    
    year_df <- year_df %>%
      dplyr::select(- pixels)
    
    #creating variable for the year a pixel is deforested
    not_defor <- rowSums(year_df)<1 *1
    defor_year <- max.col(year_df, ties.method = "first") 
    defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
    defor_df <- tibble::rownames_to_column(defor_df)
    names(defor_df)[1] <- paste("pixels")
    
    panels <- defor_df %>%
      dplyr::select(pixels, defor_year) %>%
      inner_join(panels, by = "pixels") %>%
      mutate_at(vars(pixels, county, year, property, defor_year), as.numeric)
    
    panels <- panels %>%
      mutate(indic = year - defor_year,
             defor = ifelse(indic > 0, 1, y),
             y_it = ifelse(indic > 0, NA, y))
    
    ppanels <- panels %>%
      group_by(property, year)%>%
      mutate(ptreat = mean(treat))
    
     
    for(k in cellsize_list){
      
      this_panel <- panels %>%
        rename(this_grid = paste0("grid_", k),
               this_garea = paste0("garea_", k))%>%
        select(this_grid, treat, post, year, defor, this_garea, y_it, property, county)%>%
        group_by(this_grid, year) %>%
        mutate(gtreat = mean(treat))
      
      
      this_grid_df <- as.data.frame(this_panel) %>%
        dplyr::group_by(this_grid, year, post) %>%
        dplyr::summarise(this_garea = mean(this_garea),
                         defor = mean(defor),
                         treat = mean(treat))%>%
        ungroup()
      
      
      this_grid_df <- as.data.frame(this_grid_df[order(this_grid_df$this_grid, this_grid_df$year),])
      this_grid_df <- slide(this_grid_df, Var = "defor", GroupVar = "this_grid", NewVar = "deforlag",
                            slideBy = -1, reminder = FALSE)
      
      ##### creating forested share variable #####
      this_grid_df$forshare <- 1 -  this_grid_df$defor
      this_grid_df$forsharelag <- 1 -  this_grid_df$deforlag
      
      #generate outcome var
      this_grid_df$deforrate <- ((this_grid_df$forsharelag- this_grid_df$forshare) / this_grid_df$forsharelag)
      #remove any infinite values
      #gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
      this_grid_df <- 
        this_grid_df %>% 
        filter_all(all_vars(!is.infinite(.)))%>%
        drop_na(deforrate)
      
      pos <- as.numeric(which(cellsize_list == k))
      
      # grid unit of analysis
      DID_grid <- feols(deforrate ~  post:treat|year+this_grid, data = this_grid_df)
      grid_coeffmatrix[i, pos] <- DID_grid$coefficients - ATT
      #clustering at grid level
      clse    <- tail(summary(DID_grid, cluster = ~this_grid)$se, n=1)
      grid_covermat[i,pos] <- dplyr::between(ATT, tail(DID_grid$coefficients, n=1) - 1.96 * clse, tail(DID_grid$coefficients, n=1) + 1.96 * clse)*1
      
      # aggregated fixed effects
      pix_twfe <- suppressMessages(feols(y_it ~  post:treat| year + this_grid + county, data = this_panel))
      gridfe_coeffmatrix[i,pos] <- tail(pix_twfe$coefficients, n=1) - ATT
      #clustering at grid level
      clse    <- tail(summary(pix_twfe, cluster = ~this_grid)$se, n=1)
      gridfe_covermat[i,pos] <- dplyr::between(ATT, tail(pix_twfe$coefficients, n=1) - 1.96 * clse, tail(pix_twfe$coefficients, n=1) + 1.96 * clse)*1
      
    }
    
    # aggregate up to property in each year 
    proplevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(property, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat),
                       parea = mean(parea))%>%
      ungroup()
    
    proplevel_df <- as.data.frame(proplevel_df[order(proplevel_df$property, proplevel_df$year),])
    proplevel_df <- slide(proplevel_df, Var = "defor", GroupVar = "property", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    proplevel_df$forshare <- 1 -  proplevel_df$defor
    proplevel_df$forsharelag <- 1 -  proplevel_df$deforlag
    
    #generate outcome var
    proplevel_df$deforrate <- ((proplevel_df$forsharelag- proplevel_df$forshare) / proplevel_df$forsharelag)
    #remove any infinite values
    
    proplevel_df <- 
      proplevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))%>%
      na.omit(deforrate)%>%
      drop_na(deforrate)
    
    # aggregate up to county in each year 
    countylevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(county, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat),
                       carea = mean(carea))%>%
      ungroup()
    
    countylevel_df <- as.data.frame(countylevel_df[order(countylevel_df$county, countylevel_df$year),])
    countylevel_df <- slide(countylevel_df, Var = "defor", GroupVar = "county", NewVar = "deforlag",
                            slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    countylevel_df$forshare <- 1 -  countylevel_df$defor
    countylevel_df$forsharelag <- 1 -  countylevel_df$deforlag
    
    #generate outcome var
    countylevel_df$deforrate <- ((countylevel_df$forsharelag- countylevel_df$forshare) / countylevel_df$forsharelag)
    #remove any infinite values
    
    countylevel_df <- 
      countylevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))%>%
      drop_na(deforrate)
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ##### RUNNING MODELS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # simple DID
    DID <- suppressMessages(feols(y_it ~  post*treat, data = panels))
    
    # aggregated units of analysis
    prop_DID <- feols(deforrate ~  post:treat|year+property, data = proplevel_df)
    county_DID <- feols(deforrate ~  post:treat|year+county, data = countylevel_df)
    
    ### TWFE regressions with aggregated fixed effects
    countyfe_DID <- suppressMessages(feols(y_it ~  post:treat|year + county, data = panels))
    propcountyfe_DID <- suppressMessages(feols(y_it ~  post:treat|year + property + county, data = panels))
    proptreatfe_DID <- suppressMessages(feols(y_it ~  post:treat|year + treat + property, data = panels))
    
    
    # weighting by county and property area
    weight_DIDp <- suppressMessages(feols(deforrate ~  post:treat|year+property, data = proplevel_df, weights = proplevel_df$parea))
    weight_DIDc <- suppressMessages(feols(deforrate ~  post:treat|year+county, weights = countylevel_df$carea, data = countylevel_df))
    
    # problematic specifications
    bad_DID1 <- feols(defor ~  post*treat, data = panels)
    bad_DID2 <- suppressMessages(feols(y_it ~  post*treat|year + pixels, data = panels))
    
    # regular DID clustered at pixel
    DID_se  <- tail(summary(DID, cluster = ~pixels)$se, n=1)
    
    #clustering at group level for aggregated analyses
    prop_se    <- tail(summary(prop_DID, cluster = ~property)$se, n=1)
    county_se    <- tail(summary(county_DID, cluster = ~county)$se, n=1)
    
    # same for weighted
    weightp_se    <- tail(summary(weight_DIDp, cluster = ~property)$se, n=1)
    weightc_se   <- tail(summary(weight_DIDc, cluster = ~county)$se, n=1)
    
    # aggregated fixed effects
    propcountyfe_pse    <- tail(summary(propcountyfe_DID, cluster = ~property)$se, n=1)
    propcountyfe_cse    <- tail(summary(propcountyfe_DID, cluster = ~county)$se, n=1)
    proptreatfe_se    <- tail(summary(proptreatfe_DID, cluster = ~property)$se, n=1)
    countyfe_se    <- tail(summary(countyfe_DID, cluster = ~county)$se, n=1)
    
    # se for bad specifications clustered at pixel level
    bad_se1    <- tail(summary(bad_DID1, cluster = ~pixels)$se, n=1)
    bad_se2    <- tail(summary(bad_DID2, cluster = ~pixels)$se, n=1)
    
    
    DID_cover <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_se, DID$coefficients[4] + 1.96 * DID_se)*1
    
    prop_cover <- dplyr::between(ATT, prop_DID$coefficients - 1.96 * prop_se, prop_DID$coefficients + 1.96 * prop_se)*1
    county_cover <- dplyr::between(ATT, county_DID$coefficients - 1.96 * county_se, county_DID$coefficients + 1.96 * county_se)*1
    
    weightp_cover <- dplyr::between(ATT, weight_DIDp$coefficients - 1.96 * weightp_se, weight_DIDp$coefficients + 1.96 * weightp_se)*1
    weightc_cover <- dplyr::between(ATT, weight_DIDc$coefficients - 1.96 * weightc_se, weight_DIDc$coefficients + 1.96 * weightc_se)*1
    
    propcountyfe_pcover <- dplyr::between(ATT, tail(propcountyfe_DID$coefficients, n=1) - 1.96 * propcountyfe_pse, tail(propcountyfe_DID$coefficients, n=1) + 1.96 * propcountyfe_pse)*1
    propcountyfe_ccover <- dplyr::between(ATT, tail(propcountyfe_DID$coefficients, n=1) - 1.96 * propcountyfe_cse, tail(propcountyfe_DID$coefficients, n=1) + 1.96 * propcountyfe_cse)*1
    
    countyfe_cover <- dplyr::between(ATT, tail(countyfe_DID$coefficients, n=1) - 1.96 * countyfe_se, tail(countyfe_DID$coefficients, n=1) + 1.96 * countyfe_se)*1
    proptreatfe_cover <- dplyr::between(ATT, tail(proptreatfe_DID$coefficients, n=1) - 1.96 * proptreatfe_se, tail(proptreatfe_DID$coefficients, n=1) + 1.96 * proptreatfe_se)*1
    
    
    bad_cover1 <- dplyr::between(ATT, tail(bad_DID1$coefficients, n=1) - 1.96 * bad_se1, tail(bad_DID1$coefficients, n=1) + 1.96 * bad_se1)*1
    bad_cover2 <- dplyr::between(ATT, bad_DID2$coefficients - 1.96 * bad_se2, bad_DID2$coefficients + 1.96 * bad_se2)*1
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### survival analysis
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    haz_rat <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) / pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
    
    surv_df <- panels %>% 
      mutate(t_start = year - 1,
             t_end = year,
             # t_end = ifelse((t_end==20) & (y_it==0), Inf, t_end),
             outcome = y_it,
             treat_now = treat * post) %>% 
      select(pixels, t_start, t_end, outcome, treat, treat_now, post, year) %>% 
      drop_na()%>%
      mutate(t_start = ifelse(post==1, t_start -10, t_start),
             t_end = ifelse(post==1, t_end -10, t_end))
    
    
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
    
    hr_11_cf_cover <- dplyr::between(haz_rat, hr_11_cf - 1.96 * se_hr_11_cf, hr_11_cf + 1.96 * se_hr_11_cf)*1
    cox_cover <- dplyr::between(haz_rat, hr_did - 1.96 * se_cox_did, hr_did + 1.96 * se_cox_did)*1
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### filling in summary long dataframe output
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="notes")
    
    # keeping deforested pixels
    summary_long[i,c(firstcol:lastcol)] <- c(
      i, 
      1,0,0,0,
      0,0,0,0,1,
      0,0,0,
      1,0,0,0,
      bad_DID1$coefficients[4] - ATT,
      bad_cover1, 0, 
      "keeping pixels after deforestation event"
    )
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="prop_concern")
    
    # Pixel fixed effects
    summary_long[i+n,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      1,0,0,0,0,
      0,0,0,
      1,0,0,0,
      bad_DID2$coefficients - ATT,
      bad_cover2,
      0
    )
    
    # traditional DID
    summary_long[i+n*2,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,0,0,
      1,0,0,0,
      DID$coefficients[4] - ATT,
      DID_cover,
      0
    )
    
    # property uoa
    summary_long[i+n*3,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0,
      0,0,1,0,0,
      0,0,0,
      0,0,1,0,
      prop_DID$coefficients - ATT,
      prop_cover,
      0
    )
    
    # county uoa
    summary_long[i+n*4,c(firstcol:lastcol)] <- c(
      i,
      0,0,0,1,
      0,0,0,1,0,
      0,0,0,
      0,0,0,1,
      county_DID$coefficients - ATT,
      county_cover,
      0
    )
    
    proptest <- panels %>%
      filter(year == max(year))%>%
      group_by(property)%>%
      dplyr::summarise(test = mean(y_it, na.rm = T))%>%
      filter(is.na(test))
    
    if(length(proptest$test) > 0){
      
      print(paste0(length(proptest$test), " properties fully deforested before end of study period"))
      
    }
    
    
    # property county fe
    summary_long[i+n*5,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,1,1,0,
      0,0,0,
      0,0,1,0,
      tail(propcountyfe_DID$coefficients, n=1) - ATT,
      propcountyfe_pcover,
      ifelse(length(proptest$test) > 0, 1, 0)
    )
    
    # property county fe
    summary_long[i+n*6,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,1,1,0,
      0,0,0,
      0,0,0,1,
      tail(propcountyfe_DID$coefficients, n=1) - ATT,
      propcountyfe_ccover,
      ifelse(length(proptest$test) > 0, 1, 0)
    )
    
    # property treat fe
    summary_long[i+n*7,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,1,0,1,
      0,0,0,
      0,0,1,0,
      tail(proptreatfe_DID$coefficients, n=1) - ATT,
      proptreatfe_cover,
      ifelse(length(proptest$test) > 0, 1, 0)
    )
    
    # county fe
    summary_long[i+n*8,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,1,0,
      0,0,0,
      0,0,0,1,
      tail(countyfe_DID$coefficients, n=1) - ATT,
      countyfe_cover,
      0
    )
    
    firstcol = which(colnames(summary_long)=="gridsize")
    lastcol = which(colnames(summary_long)=="cover")
    
    ### GRID specifications
    # AS UNIT OF ANALYSIS
    summary_long[i+n*9,c(firstcol:lastcol)] <- c(
      cellsize_small, i,
      0,1,0,0,
      0,1,0,0,0,
      0,0,0,
      0,1,0,0,
      grid_coeffmatrix[i, 1],
      grid_covermat[i,1]
    )
    
    summary_long[i+n*10,c(firstcol:lastcol)] <- c(
      cellsize_med, i,
      0,1,0,0,
      0,1,0,0,0,
      0,0,0,
      0,1,0,0,
      grid_coeffmatrix[i, 2],
      grid_covermat[i,2]
    )
    
    summary_long[i+n*11,c(firstcol:lastcol)] <- c(
      cellsize_large, i,
      0,1,0,0,
      0,1,0,0,0,
      0,0,0,
      0,1,0,0,
      grid_coeffmatrix[i, 3],
      grid_covermat[i,3]
    )
    
    # AS FE
    
    summary_long[i+n*12,c(firstcol:lastcol)] <- c(
      cellsize_small, i,
      1,0,0,0,
      0,1,0,1,0,
      0,0,0,
      0,1,0,0,
      gridfe_coeffmatrix[i, 1],
      gridfe_covermat[i,1]
    )
    
    summary_long[i+n*13,c(firstcol:lastcol)] <- c(
      cellsize_med, i,
      1,0,0,0,
      0,1,0,1,0,
      0,0,0,
      0,1,0,0,
      gridfe_coeffmatrix[i, 2],
      gridfe_covermat[i,2]
    )
    
    summary_long[i+n*14,c(firstcol:lastcol)] <- c(
      cellsize_large, i,
      1,0,0,0,
      0,1,0,1,0,
      0,0,0,
      0,1,0,0,
      gridfe_coeffmatrix[i, 3],
      gridfe_covermat[i,3]
    )
    
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="cover")
    
    ### WEIGHTED specifications
    # weighted property
    summary_long[i+n*15,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0,
      0,0,1,0,0,
      1,0,0,
      0,0,1,0,
      weight_DIDp$coefficients - ATT,
      weightp_cover
    )
    
    #weighted county
    summary_long[i+n*16,c(firstcol:lastcol)] <- c(
      i,
      0,0,0,1,
      0,0,0,1,0,
      1,0,0,
      0,0,0,1,
      weight_DIDc$coefficients - ATT,
      weightc_cover
    )
    
    ### SURVIVAL specifications
    # Cox DID
    summary_long[i+n*17,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,1,0,
      1,0,0,0,
      ATT_cox- ATT,
      cox_cover
    )
    
    # our proposed estimator
    summary_long[i+n*18,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,0,
      0,1,1,
      1,0,0,0,
      ATT_11_cf- ATT,
      hr_11_cf_cover
    )
    
    
    
    
    print(i)
    toc()
  }
  
  outputs = list(
    "summary_long" = summary_long)
  
  return(outputs)
  
  #end function  
}  