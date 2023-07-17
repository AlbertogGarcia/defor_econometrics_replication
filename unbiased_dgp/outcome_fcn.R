library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(Metrics)
library(DataCombine)
library(tictoc)
library(fixest)
source(here::here('unbiased_dgp', 'grid_landscape.R'))

join <- fabricatr::join_using

#begin function
outcome_fcn <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.5, cellsize){
  
  gridscape = grid_landscape(nobs, cellsize)
  pixloc_df = gridscape$pixloc_df
  gridcoords = gridscape$gridcoords
  #N_treat = gridscape$N_treat    
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2)^.5)
  
  pixloc <- pixloc_df %>%
    mutate_at(vars(pixels), as.character)
  
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  
  for(i in 1:n){
    tic("loop")
    
    Nobs <- length(pixloc$treat)   
    panels <- fabricate(
      pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
      year = add_level(N = (years*2), nest = FALSE),
      obs = cross_levels(
        by = join(pixels, year),
        year = as.numeric(year),
        post = ifelse(year > years, 1, 0),
        v_it = rnorm(N, 0, std_v),
        ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
        y = (ystar > 0)*1
      )
    )
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      mutate(year = as.numeric(year),
             y = (ystar > 0)*1 
             )
    
    #need to determine which year deforestation occurred
    year_df <- panels %>%
      dplyr::select(pixels, year, y) %>%
      #dcast(pixels ~ year , value.var = "y")
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
      inner_join(panels, by = "pixels") %>%
      mutate_at(vars(pixels, year, defor_year), as.numeric)
    
    panels <- panels %>%
      mutate(indic = year - defor_year,
             defor = ifelse(indic > 0, 1, y),
             y_it = ifelse(indic > 0, NA, y))
    
    
    # Aggregate to grid level
    gridlevel_df <- as.data.frame(panels) %>%
      dplyr::group_by(grid, year, post) %>%
      dplyr::summarise(defor = mean(defor),
                       treat = mean(treat)) %>%
      ungroup()
    
    
    gridlevel_df <- as.data.frame(gridlevel_df[order(gridlevel_df$grid, gridlevel_df$year),])%>%
      group_by(grid)%>%
      mutate(deforlag = lag(defor))%>%
      ungroup%>%
      mutate(forshare = 1 - defor,
             forsharelag = 1 - deforlag)
    
    baseline_forshare <- gridlevel_df %>%
      filter(year == min(year)) %>%
      mutate(forshare0 = forshare)%>%
      select(grid, forshare0)
    
    #generate outcome vars
    gridlevel_df <- gridlevel_df %>%
      left_join(baseline_forshare, by = "grid")%>%
      mutate(deforrate1 = (forsharelag - forshare) / forsharelag) %>%
      mutate(deforrate2 = (forsharelag - forshare) / forshare0) %>%
      mutate(deforrate3 = log(forsharelag / forshare))
    
    #remove any infinite values
    #gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
    gridlevel_df <- gridlevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # run two-way fixed effects with outcomes 
    coeffmatrix[i,1] <- feols(deforrate1 ~  post*treat|year+grid, data = gridlevel_df)$coefficients - ATT
    
    coeffmatrix[i,2] <- feols(deforrate2 ~  post*treat|year+grid, data = gridlevel_df)$coefficients - ATT
    
    coeffmatrix[i,3] <- feols(deforrate3 ~  post*treat|year+grid, data = gridlevel_df)$coefficients - ATT
    
    print(i)
    toc()
  }
  
  
  coeff_bias <- as.data.frame(coeffmatrix)
  names(coeff_bias)[1] <- paste("outcome1")
  names(coeff_bias)[2] <- paste("outcome2")
  names(coeff_bias)[3] <- paste("outcome3")
  
  outputs = list("coeff_bias" = coeff_bias)
  return(outputs)
  
  #end function  
}  
