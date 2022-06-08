library(tidyverse)
library(tictoc)
library(fixest)
library(here)
library(DeclareDesign)
library(DataCombine)
source(here::here('multigroup_dgp', 'multi_group_landscape.R'))
source(here::here('multigroup_dgp', 'my_event_study.R'))

#begin function
trends_fcn <- function(nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a = 0.0, std_v = 0.25, std_p = 0.0, cellsize, ppoints, cpoints){
  
  countyscape = multi_group_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  pixloc <- pixloc_df
  
  #create data frame with 0 rows and 5 columns
  pixel_es_long <- data.frame(matrix(ncol = 5, nrow = 0))
  county_es_long <- pixel_es_long
  
  #provide column names
  cnames <- c('term', 'estimate', 'std.error', 'estimator', 'iteration')
  colnames(pixel_es_long) <- cnames
  colnames(county_es_long) <- cnames
  
  
  std_avp = (std_a^2+std_v^2+std_p^2)^.5
  
  
  ################################################################################################
  ####### early group
  ################################################################################################
  b0a = qnorm(base_a, mean = 0, sd = std_avp)
  
  b1a = qnorm(trend1 + base_a, mean = 0, sd = std_avp) -b0a
  b2a = qnorm(trend2 + pnorm(b0a+b1a, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1a - b0a
  b3a = qnorm(trend3 + pnorm(b0a+b1a+b2a, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1a - b0a - b2a
  b4a = qnorm(trend3 + pnorm(b0a+b1a+b2a+b3a, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1a - b0a - b2a - b3a
  
  tau_a = qnorm( pnorm(b0a+b1a+b2a, mean = 0, sd = std_avp) + ATT_a , mean = 0, sd = std_avp) - (b0a+b1a+b2a)
  tau_a2 = qnorm( pnorm(b0a+b1a+b2a+b3a+tau_a, mean = 0, sd = std_avp) + dyn_ATT_a , mean = 0, sd = std_avp) - (b0a+b1a+b2a+b3a+tau_a)
  tau_a3 = qnorm( pnorm(b0a+b1a+b2a+b3a+b4a+tau_a2, mean = 0, sd = std_avp) + dyn_ATT_a , mean = 0, sd = std_avp) - (b0a+b1a+b2a+b3a++b4a+tau_a2)
  
  ################################################################################################
  ####### late group
  ################################################################################################
  b0b = qnorm(base_b, mean = 0, sd = std_avp)
  
  b1b = qnorm(trend1 + base_b, mean = 0, sd = std_avp) -b0b
  b2b = qnorm(trend2 + pnorm(b0b+b1b, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1b - b0b
  b3b = qnorm(trend3 + pnorm(b0b+b1b+b2b, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1b - b0b - b2b
  b4b = qnorm(trend3 + pnorm(b0b+b1b+b2b+b3b, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1b - b0b - b2b - b3b
  
  tau_b = qnorm( pnorm(b0b+b1b+b2b+b3b, mean = 0, sd = std_avp) + ATT_b , mean = 0, sd = std_avp) - (b0b+b1b+b2b+b3b)
  tau_b2 = qnorm( pnorm(b0b+b1b+b2b+b3b+b4b+tau_b, mean = 0, sd = std_avp) + dyn_ATT_b , mean = 0, sd = std_avp) - (b0b+b1b+b2b+b3b+b4b+tau_b)
  ################################################################################################
  ###########      never treated group
  ################################################################################################
  b0c = qnorm(base_c, mean = 0, sd = std_avp)
  
  b1c = qnorm(trend1 + base_c, mean = 0, sd = std_avp) -b0c
  b2c = qnorm(trend2 + pnorm(b0c+b1c, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1c - b0c
  b3c = qnorm(trend3 + pnorm(b0c+b1c+b2c, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1c - b0c - b2c
  b4c = qnorm(trend3 + pnorm(b0c+b1c+b2c+b3c, mean = 0, sd = std_avp), mean = 0, sd = std_avp) -b1c - b0c - b2c - b3c
  
  
  
  
  Nobs <- length(pixloc$G)  
    
  panels <- fabricate(
    pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), G = pixloc$G),
    year = add_level(N = 5, nest = FALSE),
    obs = cross_levels(
      by = join(pixels, year),
      year = as.numeric(year),
      post = (year >= G)*1,
      GU = (G==0)*1,
      G3 = (G==3)*1,
      G4 = (G==4)*1,
      v_it = rnorm(N, 0, std_v),
      b0 = b0a*G3 + b0b*G4 + b0c*GU,
      b1 = b1a*G3 + b1b*G4 + b1c*GU,
      b2 = b2a*G3 + b2b*G4 + b2c*GU,
      b3 = b3a*G3 + b3b*G4 + b3c*GU,
      b4 = b4a*G3 + b4b*G4 + b4c*GU,
      tau = tau_a*G3*post + tau_b*G4*post ,
      tau_2 = tau_a2*G3*((year >=G+1)*1) + tau_b2*G4*((year >=G+1)*1),
      tau_3 = tau_a3*G3*((year >=G+2)*1),
      ystar = b0 + b1*((year >=2)*1) + b2*((year >=3)*1) + b3*((year >=4)*1) + b4*((year >=5)*1) + tau +tau_2 + tau_3 + a_i + v_it,
      y = (ystar > 0)*1
    )
  )
    
    #generate random 
  errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    
  panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
  panels <- panels %>%
    inner_join(pixloc, by = c("pixels", "G")) %>%
    inner_join(errortable_property, by = "property") %>%
      #inner_join(errortable_county, by = "county") %>%
    mutate(ystar = ystar + p_err) %>%
    mutate(y = (ystar > 0)*1,
            pixels = as.numeric(pixels),
            year = as.numeric(year),
            defor_indic = ifelse(y==1, year, 99))%>%
    group_by(pixels)%>%
    mutate(defor_year = min(defor_indic),
            defor = (year>=defor_year)*1,
            y_it = ifelse(year>defor_year, NA, defor)
    )
    
    
  trends_df <- panels %>%
      group_by(G, year)%>%
      summarise(defor = mean(y, na.rm=TRUE))%>%
      mutate(Group = ifelse(G==0, "never treated", "early group\n(treated in year 3)"),
             Group = ifelse(G==4, "late group\n(treated in year 4)", Group))
    
    
  outputs = list(
    "trends_df" = trends_df
  )
  
  return(outputs)
  
}
