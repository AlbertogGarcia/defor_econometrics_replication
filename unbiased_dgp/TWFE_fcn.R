# this function's intention is to provide distributional parameters for the various specifications that have the typical binary outcome variable
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(fixest)
library(dplyr)
source(here::here('unbiased_dgp', 'deforestation_DGP.R'))
#source('deforestation_DGP.R')
#begin function

TWFE_fcn <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 2)
  
  n_mod = 2
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
    ATT = dgp_results$ATT
    
    
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
    
    # run two-way fixed effects    
    coeffmatrix[i, 2] <- feols(y_it ~  post*treat|year+pixels, data = panels
    )$coefficients - ATT
    
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
   
    
    #end for loop
    print(i)
  }  
  
  # get distribution information from matrix  
  
  summary_long <- summary_long %>%
    dplyr::select(iteration, everything())
  
  outputs = list("summary_long" = summary_long)
  
  return(outputs)
  
  #end function  
}  




