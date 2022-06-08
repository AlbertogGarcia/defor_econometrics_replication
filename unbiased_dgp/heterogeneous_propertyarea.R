#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level

library(ggplot2)
library(clubSandwich)
library(reshape2)
library(matrixStats)
library(ggplot2)
# library(plm)
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
  
  did_covermat <- matrix(nrow = n, ncol = 1)
  agg_covermat <- matrix(nrow = n, ncol = 3)
  weight_covermat <- matrix(nrow = n, ncol = 2)
  fix_covermat <- matrix(nrow = n, ncol = 3)
  
  
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  weight_coeffmatrix <- matrix(nrow = n, ncol = 3)
  fix_coeffmatrix <- matrix(nrow = n, ncol = 4)
  

  n_mod = 9
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
    
    #ATT = ls_ATT
    # aggregate up to county in each year 
    gridlevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(grid, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat),
                       garea = mean(garea))%>%
      ungroup()
    
    
    gridlevel_df <- gridlevel_df[order(gridlevel_df$grid, gridlevel_df$year),]
    gridlevel_df <- slide(gridlevel_df, Var = "defor", GroupVar = "grid", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    gridlevel_df$forshare <- 1 -  gridlevel_df$defor
    gridlevel_df$forsharelag <- 1 -  gridlevel_df$deforlag
    
    #generate outcome var
    gridlevel_df$deforrate <- ((gridlevel_df$forsharelag- gridlevel_df$forshare) / gridlevel_df$forsharelag)
    #remove any infinite values
    #gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
    gridlevel_df <- 
      gridlevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregate up to property in each year 
    proplevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(property, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat),
                       parea = mean(parea))%>%
      ungroup()
    
    
    proplevel_df <- proplevel_df[order(proplevel_df$property, proplevel_df$year),]
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
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregate up to county in each year 
    countylevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(county, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat),
                       carea = mean(carea))%>%
      ungroup()
    
    countylevel_df <- countylevel_df[order(countylevel_df$county, countylevel_df$year),]
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
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregated units of analysis
    
    agg_DID1 <- feols(deforrate ~  post*treat|year+grid, data = gridlevel_df)
    
    
    agg_DID2 <- feols(deforrate ~  post*treat|year+property, data = proplevel_df)
    
    
    agg_DID3 <- feols(deforrate ~  post*treat|year+county, data = countylevel_df)
    
    # weighting by county and property area
    
    weight_DIDp <- feols(deforrate ~  post*treat|year+property, data = proplevel_df, weights = proplevel_df$parea)
    
    
    weight_DIDc <- feols(deforrate ~  post*treat|year+county, weights = countylevel_df$carea, data = countylevel_df)
    
    
    ### TWFE regressions with aggregated fixed effects
    
    fix_DID1 <- feols(y_it ~  post*treat|grid, data = panels)
    
    fix_DID2 <- feols(y_it ~  post:treat|property, data = panels)
    
    fix_DID3 <- feols(y_it ~  post*treat|year+county, data = panels)
    
    # simple DID
    
    DID <- feols(y_it ~  post*treat, data = panels)
    
    
    #calculating bias from each method
    
    # aggregated units of analysis
    coeffmatrix[i,1] <- agg_DID1$coefficients 
    coeffmatrix[i,2] <- agg_DID2$coefficients 
    coeffmatrix[i,3] <- agg_DID3$coefficients 
    
    # with weights
    weight_coeffmatrix[i,1] <- weight_DIDp$coefficients 
    weight_coeffmatrix[i,2] <- weight_DIDc$coefficients 
    
    # aggregated fixed effects
    fix_coeffmatrix[i,1] <- tail(fix_DID1$coefficients, n=1) 
    fix_coeffmatrix[i,2] <- tail(fix_DID2$coefficients, n=1) 
    fix_coeffmatrix[i,3] <- tail(fix_DID3$coefficients, n=1) 
    
    # regular did
    fix_coeffmatrix[i,4] <- DID$coefficients[4] 
    
    #regular DID standard errors
    # clustered at pixel
    DID_clse_pixel  <- tail(summary(DID, cluster = ~pixels)$se, n=1)
   
    did_covermat[i,1] <- dplyr::between(ATT, DID$coefficients[4] - 1.96 * DID_clse_pixel, DID$coefficients[4] + 1.96 * DID_clse_pixel)*1
    
    #clustering at group level for aggregated analyses
    agg_clse1    <- tail(summary(agg_DID1, cluster = ~grid)$se, n=1)
    agg_clse2    <- tail(summary(agg_DID2, cluster = ~property)$se, n=1)
    agg_clse3    <- tail(summary(agg_DID3, cluster = ~county)$se, n=1)
    
    # same for weighted
    weight_clsep    <- tail(summary(weight_DIDp, cluster = ~property)$se, n=1)
    weight_clsec   <- tail(summary(weight_DIDc, cluster = ~county)$se, n=1)
    
    # coverage with clustered standard errors at group level
    agg_covermat[i,1] <- dplyr::between(ATT, agg_DID1$coefficients - 1.96 * agg_clse1, agg_DID1$coefficients + 1.96 * agg_clse1)*1
    agg_covermat[i,2] <- dplyr::between(ATT, agg_DID2$coefficients - 1.96 * agg_clse2, agg_DID2$coefficients + 1.96 * agg_clse2)*1
    agg_covermat[i,3] <- dplyr::between(ATT, agg_DID3$coefficients - 1.96 * agg_clse3, agg_DID3$coefficients + 1.96 * agg_clse3)*1
    
    #weighted
    weight_covermat[i,1] <- dplyr::between(ATT, weight_DIDp$coefficients - 1.96 * weight_clsep, weight_DIDp$coefficients + 1.96 * weight_clsep)*1
    weight_covermat[i,2] <- dplyr::between(ATT, weight_DIDc$coefficients - 1.96 * weight_clsec, weight_DIDc$coefficients + 1.96 * weight_clsec)*1
    
    ### clustered standard errors for aggregated fixed effects
    fix_clse1    <- tail(summary(fix_DID1, cluster = ~grid)$se, n=1)
    fix_clse2    <- tail(summary(fix_DID2, cluster = ~property)$se, n=1)
    fix_clse3    <- tail(summary(fix_DID3, cluster = ~county)$se, n=1)
    
    #whether att is within CI
    
    
    fix_covermat[i,1] <- dplyr::between(ATT, tail(fix_DID1$coefficients, n=1) - 1.96 * fix_clse1, tail(fix_DID1$coefficients, n=1) + 1.96 * fix_clse1)*1
    fix_covermat[i,2] <- dplyr::between(ATT, tail(fix_DID2$coefficients, n=1) - 1.96 * fix_clse2, tail(fix_DID2$coefficients, n=1) + 1.96 * fix_clse2)*1
    fix_covermat[i,3] <- dplyr::between(ATT, tail(fix_DID3$coefficients, n=1) - 1.96 * fix_clse3, tail(fix_DID3$coefficients, n=1) + 1.96 * fix_clse3)*1
    
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="given_ATT")
    
    
    summary_long[i,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,0,1,
      0,
      1,0,0,0,
      fix_coeffmatrix[i,4],
      did_covermat[i,1],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n,c(firstcol:lastcol)] <- c(
      i,
      0,1,0,0,
      0,1,0,0,0,
      0,
      0,1,0,0,
      coeffmatrix[i,1],
      agg_covermat[i,1],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*2,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0,
      0,0,1,0,0,
      0,
      0,0,1,0,
      coeffmatrix[i,2],
      agg_covermat[i,2],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*3,c(firstcol:lastcol)] <- c(
      i,
      0,0,0,1,
      0,0,0,1,0,
      0,
      0,0,0,1,
      coeffmatrix[i,3],
      agg_covermat[i,3],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*4,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,1,0,0,1,
      0,
      0,1,0,0,
      fix_coeffmatrix[i,1],
      fix_covermat[i,1],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*5,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,1,0,0,
      0,
      0,0,1,0,
      fix_coeffmatrix[i,2],
      fix_covermat[i,2],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*6,c(firstcol:lastcol)] <- c(
      i,
      1,0,0,0,
      0,0,0,1,1,
      0,
      0,0,0,1,
      fix_coeffmatrix[i,3],
      fix_covermat[i,3],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*7,c(firstcol:lastcol)] <- c(
      i,
      0,0,1,0,
      0,0,1,0,0,
      1,
      0,0,1,0,
      weight_coeffmatrix[i,1],
      weight_covermat[i,1],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    summary_long[i+n*8,c(firstcol:lastcol)] <- c(
      i,
      0,0,0,1,
      0,0,0,1,0,
      1,
      0,0,0,1,
      weight_coeffmatrix[i,2],
      weight_covermat[i,2],
      NA,
      p_ATT, ls_ATT, ATT
    )
    
    
    print(i)
    toc()
  }
  
  
  summary_long <- summary_long %>%
    mutate_at(vars(p_ATT, ls_ATT, estimate, cover), as.numeric)%>%
    mutate(ATT_diff = p_ATT - ls_ATT)%>%
    select(iteration, everything())
  
  
  outputs = list("summary_long" = summary_long, "unit_area" = unit_area)
  return(outputs)
  
  #end function  
}  
