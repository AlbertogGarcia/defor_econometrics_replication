# this function's intention is to provide distributional parameters for the various specifications that have the typical binary outcome variable
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(Metrics)
library(fixest)
source(here::here('unbiased_dgp', 'deforestation_DGP.R'))
#source('deforestation_DGP.R')
#begin function

TWFE_expost <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25){
  
  #preallocate n x 4 matrix
  
  n_mod = 10
  coeffmatrix <- matrix(nrow = n, ncol = n_mod)
  
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'b3'= rep(b3, summ_row), 
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(0, summ_row),
                             'iteration' = rep(NA, summ_row), 
                             'model'=rep(NA, summ_row),
                             'bias'=rep(NA, summ_row), 
                             stringsAsFactors=FALSE)
  
  for(i in 1:n){
    
    # call defor_sim function to simulate dataframe, returned as panels  
    dgp_results <- deforestation_DGP(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
    panels = dgp_results$panels
    
    final_panels <- panels %>%
      filter(year == max(year))
    
    dropped_panels <- panels %>%
      drop_na(y_it)
    
    # drop first period
    ex_post_panels <- panels %>%
      filter(post == 1)
      
    survivors_panel_period1 <- panels %>%
      mutate(drop = ifelse(post == 0 & y_it == 1, 1, 0))%>%
      group_by(pixels)%>%
      filter(max(drop) == 0)
    
    survivors_panel <- panels %>%
      mutate(drop = ifelse(year == max(year) & is.na(y_it)
                           , 1, 0))%>%
      group_by(pixels)%>%
      filter(max(drop) == 0)
      
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
    
    # run two-way fixed effects    
    coeffmatrix[i, 2] <- feols(y_it ~  post*treat|year+pixels, data = panels
    )$coefficients - ATT
    
    coeffmatrix[i,3]  <- lm(y_it ~  post*treat, 
                            data = survivors_panel
    )$coefficients[4] - ATT
    
    coeffmatrix[i,4]  <- lm(y_it ~  post*treat, 
                            data = survivors_panel_period1
    )$coefficients[4] - ATT
    
    coeffmatrix[i,5] <- mean(subset(ex_post_panels, treat == 1)$y_it, na.rm = TRUE) - mean(subset(ex_post_panels, treat == 0)$y_it, na.rm = TRUE) - ATT
    
    coeffmatrix[i,6] <- mean(subset(final_panels, treat == 1)$y_it, na.rm = TRUE) - mean(subset(final_panels, treat == 0)$y_it, na.rm = TRUE) - ATT
    
    coeffmatrix[i, 7] <- feols(y_it ~  post*treat|year+pixels, data = survivors_panel
    )$coefficients - ATT
    
    coeffmatrix[i, 8] <- feols(y_it ~  post*treat|year+pixels, data = survivors_panel_period1
    )$coefficients - ATT
    
    coeffmatrix[i, 9] <- tail(feols(y_it ~  post*treat, data = survivors_panel_period1
    )$coefficients, n=1) - ATT
    
    coeffmatrix[i, 10] <- tail(feols(y ~  post*treat, data = dropped_panels
    )$coefficients, n=1) - ATT
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="bias")
    
    summary_long[i,c(firstcol:lastcol)] <- c(
      i,
      "DID",
      coeffmatrix[i,1]
    )
    
    summary_long[i+n,c(firstcol:lastcol)] <- c(
      i,
      "TWFE",
      coeffmatrix[i,2]
    )
    
    summary_long[i+2*n,c(firstcol:lastcol)] <- c(
      i,
      "DID on dataset dropping deforested pixels prior last period",
      coeffmatrix[i,3]
    )
    
    summary_long[i+3*n,c(firstcol:lastcol)] <- c(
      i,
      "DID on dataset dropping deforested pixels prior to treatment",
      coeffmatrix[i,4]
    )
    
    summary_long[i+4*n,c(firstcol:lastcol)] <- c(
      i,
      "ex-post difference in means",
      coeffmatrix[i,5]
    )
    
    summary_long[i+5*n,c(firstcol:lastcol)] <- c(
      i,
      "final period difference in means",
      coeffmatrix[i,6]
    )
    
    summary_long[i+6*n,c(firstcol:lastcol)] <- c(
      i,
      "TWFE on dataset dropping deforested pixels prior last period",
      coeffmatrix[i,7]
    )
    
    summary_long[i+7*n,c(firstcol:lastcol)] <- c(
      i,
      "TWFE on dataset dropping deforested pixels prior to treatment",
      coeffmatrix[i,8]
    )
    
    summary_long[i+8*n,c(firstcol:lastcol)] <- c(
      i,
      "TWFE on only post-treatment period",
      coeffmatrix[i,9]
    )
    
    summary_long[i+9*n,c(firstcol:lastcol)] <- c(
      i,
      "TWFE on latent outcome y",
      coeffmatrix[i,10]
    )
    
    #end for loop
    print(i)
  }  
  # 
  
  
  summary_long <- summary_long %>%
    dplyr::select(iteration, everything())
  
  
  outputs = list("summary_long" = summary_long)
  return(outputs)
  
  #end function  
}  




