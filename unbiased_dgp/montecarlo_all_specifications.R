library(purrr)

montecarlo_all_specifications <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.0, std_v = 0.5, std_p = 0.0, 
                                          cellsize_list, ppoints, cpoints, 
                                          min_psize = 1,
                                          nestedprops = FALSE, proptreatassign = FALSE){
  
  source(here::here('unbiased_dgp', 'all_specification_run.R'))
  source(here::here('unbiased_dgp', 'multigrid_landscape.R'))
  source(here::here('unbiased_dgp', 'nested_landscape.R'))
  source(here::here('unbiased_dgp', 'proptreat_landscape.R'))
  
  options(dplyr.summarise.inform = FALSE)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##### Generate landscape, which will be used in all iterations of monte carlo
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  minpropsize = 0
  while(minpropsize < max(min_psize, 1)){
    
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
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##### Run monte carlo n times
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  all_spec_results <- purrr::map(seq_len(n), 
                                 ~all_specification_run(
                                   nobs, years, b0, b1, b2_0, b2_1, b3, std_a = std_a, std_v = std_v, std_p = std_p, ppoints, cpoints,
                                   cellsize_list, pixloc
                                 )
                     , .progress = TRUE)
  
  results_long <- rbindlist(all_spec_results, idcol = 'iteration')
  
  return(results_long)
  
  # end function
}