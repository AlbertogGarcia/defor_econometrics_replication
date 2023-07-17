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
source(here::here('unbiased_dgp', 'full_landscape.R'))
source(here::here('unbiased_dgp', 'prop_landscape.R'))

join <- fabricatr::join_using

options(dplyr.summarise.inform = FALSE)


#begin function
heterogeneous_propertyarea <- function(n, nobs, years, b0, b1, b2_0, b2_1, std_a = 0, std_v = 0.25, std_p = 0.0, std_b3 = .1, given_ATT, cellsize, ppoints, cpoints){
  
  # if (proptreatassign) {
  #   countyscape = prop_landscape(nobs, cellsize, ppoints, cpoints)
  #   pixloc_df = countyscape$pixloc_df
  # } else {
  #   countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
  #   pixloc_df = countyscape$pixloc_df
  # }
  # 
  
  countyscape = prop_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  
  
  pixloc <- pixloc_df
  
  unit_area <- data.frame("county area" = pixloc_df$carea, "property area" = pixloc_df$parea, "grid area" = pixloc_df$garea)
  
  
  n_mod = 4
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'std_b3'= rep(std_b3, summ_row), 
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row),
                             'iteration' = rep(NA, summ_row), 
                             'pixel'=rep(NA, summ_row),'grid'=rep(NA, summ_row),'property'=rep(NA, summ_row),'county'=rep(NA, summ_row),
                             'pixel fe'=rep(NA, summ_row),'grid fe'=rep(NA, summ_row),'property fe'=rep(NA, summ_row),'county fe'=rep(NA, summ_row),'treatment fe'=rep(NA, summ_row),
                             'weights'=rep(NA, summ_row),
                             'se_pixel'=rep(NA, summ_row), 'se_grid'=rep(NA, summ_row), 'se_property'=rep(NA, summ_row), 'se_county'=rep(NA, summ_row),
                             'estimate'=rep(NA, summ_row), 'cover'=rep(NA, summ_row),'notes'=rep(NA, summ_row), 
                            "p_ATT" = rep(NA, summ_row), "ls_ATT" = rep(NA, summ_row), "given_ATT" = rep(NA, summ_row),
                             stringsAsFactors=FALSE)
  
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
        ystar_cf = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
      )
    )
    
    pixloc <- pixloc %>%
      mutate(z_parea = (parea - mean(parea)) / sd(parea),
             mu_adjust = z_parea / (1/std_b3))
    
    std_avpt = (std_a^2+std_v^2 + std_p^2 + sd(pixloc$mu_adjust)^2)^.5
    
    #generate random 
    errortable_property <- data.frame(property = as.character(unique(pixloc$property)), 
                                      p_err = rnorm(length(unique(pixloc$property)), 0, std_p)
    )
    
    b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)
    
    
    this_ATT <- pnorm(b0+b1+b2_1+b3_mu, 0, sd = std_avpt) - pnorm(b0+b1+b2_1, 0, sd = std_avp)
    
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(errortable_property, by = "property") %>%
      #inner_join(errortable_county, by = "county") %>%
      mutate(year = as.numeric(year),
             ystar = ystar_cf + (b3_mu+mu_adjust)*post*treat + p_err,
             ystar_cf = ystar_cf + p_err,
             y = (ystar > 0)*1 ,
             y_cf = (ystar_cf > 0)*1 )
    
    #need to determine which year deforestation occurred
    year_df <- panels %>%
      dplyr::select(pixels, year, y) %>%
     # dcast(pixels ~ year , value.var = "y")
      pivot_wider(names_from = year, values_from = y)
    
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
      inner_join(panels, by = "pixels")
    
    
    cols.num <- c("pixels", "grid", "property", "county", "year")
    panels[cols.num] <- sapply(panels[cols.num],as.numeric)
    
    panels <- panels %>%
      mutate(indic = year - defor_year,
             defor = ifelse(indic > 0, 1, y),
             y_it = ifelse(indic > 0, NA, y))
    
    treat_post <- subset(panels, treat==1 & post == 1)
    
    
    prop_treat_post <- treat_post %>%
      group_by(property)%>%
      summarise(y = mean(y),
                y_it = mean(y_it, na.rm = TRUE),
                y_cf = mean(y_cf))
      
      
    ls_ATT = mean( subset(panels, treat==1&post==1)$y)-mean( subset(panels, treat==1&post==1)$y_cf)
    p_ATT = mean( prop_treat_post$y)-mean( prop_treat_post$y_cf)
    
    
    
    # Aggregate to grid level
    gridlevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(grid, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat),
                       garea = mean(garea)) %>%
      ungroup()
    
    gridlevel_df <- as.data.frame(gridlevel_df[order(gridlevel_df$grid, gridlevel_df$year),])%>%
      group_by(grid)%>%
      mutate(deforlag = lag(defor))%>%
      ungroup%>%
      mutate(forshare = 1 - defor,
             forsharelag = 1 - deforlag,
             deforrate = (forsharelag - forshare) / forsharelag)%>% 
      filter_all(all_vars(!is.infinite(.)))%>%
      drop_na(deforrate)
    
    # aggregate up to property in each year 
    proplevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(property, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat),
                       parea = mean(parea)) %>%
      ungroup()
    
    proplevel_df <- as.data.frame(proplevel_df[order(proplevel_df$property, proplevel_df$year),])%>%
      group_by(property)%>%
      mutate(deforlag = lag(defor))%>%
      ungroup%>%
      mutate(forshare = 1 - defor,
             forsharelag = 1 - deforlag,
             deforrate = (forsharelag - forshare) / forsharelag)%>% 
      filter_all(all_vars(!is.infinite(.)))%>%
      drop_na(deforrate)
    
    
    # simple DID
    
    DID <- feols(y_it ~  post*treat, data = panels)
    
    # regular did coefficient
    DID_coef <- DID$coefficients[4] 
    
    #coverage of simple DID w/ standard errors clustered at pixel
    DID_clse_pixel  <- tail(summary(DID, cluster = ~pixels)$se, n=1)
    DID_cover <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_clse_pixel, DID$coefficients[4] + 1.96 * DID_clse_pixel)*1
    
    
    # aggregated units of analysis
    # aggregated units of analysis
    
    agg_DID_grid <- feols(deforrate ~  post*treat|year+grid, data = gridlevel_df)
    agg_DID_prop <- feols(deforrate ~  post*treat|year+property, data = proplevel_df)
    
    # coefficients
    DID_grid_coef <- agg_DID_grid$coefficients 
    DID_prop_coef <- agg_DID_prop$coefficients 
    
    #clustering at group level for aggregated analyses
    agg_clse_grid    <- tail(summary(agg_DID_grid, cluster = ~grid)$se, n=1)
    agg_clse_prop    <- tail(summary(agg_DID_prop, cluster = ~property)$se, n=1)
    
    # coverage with clustered standard errors at group level
    DID_grid_cover <- dplyr::between(ATT, agg_DID_grid$coefficients - 1.96 * agg_clse_grid, agg_DID_grid$coefficients + 1.96 * agg_clse_grid)*1
    DID_prop_cover <- dplyr::between(ATT, agg_DID_prop$coefficients - 1.96 * agg_clse_prop, agg_DID_prop$coefficients + 1.96 * agg_clse_prop)*1
    
    
    # weighting by property area
    weight_DIDp <- feols(deforrate ~  post*treat|year+property, data = proplevel_df, weights = proplevel_df$parea)
    
    # with weights coefficients
    weight_DIDp_coef <- weight_DIDp$coefficients 
    
    # same for weighted standard errors (group)
    weight_clsep    <- tail(summary(weight_DIDp, cluster = ~property)$se, n=1)
    
    #weighted coverage
    weight_DIDp_cover <- dplyr::between(ATT, weight_DIDp$coefficients - 1.96 * weight_clsep, weight_DIDp$coefficients + 1.96 * weight_clsep)*1
    
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="given_ATT")
    
    
    summary_long[i,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0, #pixel, grid, property, county
      0,0,0,0,1, # FE: pixel, grid, property, county, treatment
      0, # weights
      1,0,0,0, # SE: pixel, grid, prop, county
      DID_coef, # coef
      DID_cover, # coverage
      NA, #notes
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n,c(firstcol:lastcol)] <- c(
      i,
      0,1,0,0, #pixel, grid, property, county
      0,1,0,0,0, # FE: pixel, grid, property, county, treatment
      0, # weights
      0,1,0,0, # SE: pixel, grid, prop, county
      DID_grid_coef, # coef
      DID_grid_cover, # coverage
      NA, #notes
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*2,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0, #pixel, grid, property, county
      0,0,1,0,0, # FE: pixel, grid, property, county, treatment
      0, # weights
      0,0,1,0, # SE: pixel, grid, prop, county
      DID_prop_coef, # coef
      DID_prop_cover, # coverage
      NA, #notes
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*3,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0, #pixel, grid, property, county
      0,0,1,0,0, # FE: pixel, grid, property, county, treatment
      1, # weights
      0,0,1,0, # SE: pixel, grid, prop, county
      weight_DIDp_coef, # coef
      weight_DIDp_cover, # coverage
      NA, #notes
      p_ATT, ls_ATT, ATT
    )
    
    
    
    print(i)
    toc()
  }
  
  
  summary_long <- summary_long %>%
    mutate_at(vars(p_ATT, ls_ATT, estimate, cover), as.numeric)%>%
    select(iteration, everything())
  
  
  outputs = list("summary_long" = summary_long, "unit_area" = unit_area)
  return(outputs)
  
  #end function  
}  
