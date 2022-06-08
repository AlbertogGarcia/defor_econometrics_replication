# this function's intention is to provide distributional parameters for the various specifications that have the typical binary outcome variable
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(fixest)

source(here::here('unbiased_dgp', 'deforestation_DGP.R'))
#begin function

DID_keep <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.5){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 2)
  
  for(i in 1:n){
    
    # call defor_sim function to simulate dataframe, returned as panels  
    dgp_results <- deforestation_DGP(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
    panels = dgp_results$panels
    ATT = dgp_results$ATT
    
    P_a = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2)^.5)
    P_b = pnorm(b0+b1, 0, (std_a^2+std_v^2)^.5)
    P_c = pnorm(b0+b2_0, 0, (std_a^2+std_v^2)^.5)
    P_d = pnorm(b0, 0, (std_a^2+std_v^2)^.5)
    
    DID_keep_estimand = P_a + P_b - P_a*P_b - P_b - (P_c + P_d - P_c*P_d - P_d)
    
    bias = P_b - P_d - P_a*P_b + P_c*P_d
    
    coeffmatrix[i,1]  <- tail(feols(y_it ~  post*treat, 
                            data = panels
    )$coefficients, n = 1) - ATT
    
    # DID keeping variables

    coeffmatrix[i, 2] <- tail(feols(defor ~  post*treat,
                           data = panels
    )$coefficients, n=1) - ATT
    
    
    #end for loop
    print(i)
  }  
  
  # get distribution information from matrix  
  
  b_coeff <- as.data.frame(coeffmatrix)
  actual <- rep(0, times = n)
  
  names(b_coeff)[1] <- paste("DID drop")
  names(b_coeff)[2] <- paste("DID keep")
  
  suppressWarnings(cbias <- melt(b_coeff, value.name = "bias"))
  
  plot = ggplot(data = cbias, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("DID1", "DID2"), labels=c("DID dropping deforested pixels", "DID keeping deforested pixels"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = (DID_keep_estimand-ATT), color = "red", linetype = "dashed")+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("DID dropping pixels mean:", round(colMeans(coeffmatrix)[1], digits = 4),"RMSE:",round(rmse(actual, coeffmatrix[1]), digits = 5), "\n", 
                                    "DID keeping pixels mean::", round(colMeans(coeffmatrix)[2], digits = 4), "RMSE:",round(rmse(actual, coeffmatrix[2]), digits = 5)
                                    )
    )+
    theme_minimal()
  
  
  outputs = list("plot" = plot, "did_keeps" = b_coeff)
  return(outputs)
  
  #end function  
}  




