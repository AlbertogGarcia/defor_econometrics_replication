library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
source('deforestation_DGP.R')

#begin function
quickmonte <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, outcome_var = "y"){
  
  #preallocate n x 1 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 1)
  
  for(i in 1:n){
    # call defor_dgp function to simulate dataframe, returned as defor_df  
    dgp_results <- deforestation_DGP(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
    panels = dgp_results$panels
    ATT = dgp_results$ATT
    panels['outcome'] = panels[outcome_var]
    
    mod <- lm(outcome ~  post*treat, data = panels)
    coeffmatrix[i,1] <- mod$coefficients[4] - ATT
    print(i)
    
    #end for loop  
  }  
  
  #assign('coeff_didbias',coeffmatrix, envir=.GlobalEnv)
  
  did_coeff <- as.data.frame(coeffmatrix)
  actual <- rep(0, times = n)
  
  plot = ggplot() +
    geom_density(data = did_coeff , aes(x = V1), alpha = .2, fill="#29CD44")+
    #geom_vline(aes(xintercept= (DID_estimand - ATT), color="DID estimand - ATT"), linetype="dashed")+
    geom_vline(data = did_coeff, xintercept = 0, linetype = "dashed")+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Mean:", round(mean(did_coeff$V1), digits = 4),
                                    ", RMSE:", round(rmse(actual, did_coeff$V1), digits = 4)) 
    )
  
  
  
  outputs = list("plot" = plot, "did_biases" = did_coeff)
  return(outputs)
  
  #end function  
}  




