library(tidyverse)
library(DeclareDesign)

join <- fabricatr::join_using

pixloc_DGP <- function(nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.0, std_v = 0.5, std_p = 0.0, ppoints, cpoints,
                                  cellsize_list, pixloc){
  
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
  
  panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
  
  panels <- panels %>%
    inner_join(pixloc, by = c("pixels", "treat")) %>%
    inner_join(errortable_property, by = "property") %>%
    mutate(year = as.numeric(year),
           ystar = ystar + p_err,
           ystar_cf = ystar_cf + p_err,
           y = (ystar > 0)*1 ,
           y_cf = (ystar_cf > 0)*1 )
  
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
    mutate_at(vars(pixels, county, year, property, defor_year), as.numeric)
  
  panels <- panels %>%
    mutate(indic = year - defor_year,
           defor = ifelse(indic > 0, 1, y),
           y_it = ifelse(indic > 0, NA, y))
  
  ppanels <- panels %>%
    group_by(property, year)%>%
    mutate(ptreat = mean(treat))
  
  grid_df_list <- list()
  
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
      ungroup
    
    
    this_grid_df <- as.data.frame(this_grid_df[order(this_grid_df$this_grid, this_grid_df$year),])%>%
      group_by(this_grid)%>%
      mutate(deforlag = lag(defor))%>%
      ungroup%>%
      mutate(forshare = 1 - defor,
             forsharelag = 1 - deforlag,
             deforrate = (forsharelag - forshare) / forsharelag,
             cellsize = k)%>% 
      filter_all(all_vars(!is.infinite(.)))%>%
      drop_na(deforrate)
    
    # assign(paste0("grid_", k, "_df"), this_grid_df)
    
    pos <- as.numeric(which(cellsize_list == k))
    
    grid_df_list[[pos]] <- this_grid_df
    
  }
  
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
  
  
  # aggregate up to county in each year 
  countylevel_df <- as.data.frame(panels) %>%
    dplyr::group_by(county, year, post) %>%
    dplyr::summarise(defor = mean(defor),
                     treat = mean(treat),
                     carea = mean(carea))%>%
    ungroup()
  
  countylevel_df <- as.data.frame(countylevel_df[order(countylevel_df$county, countylevel_df$year),])%>%
    group_by(county)%>%
    mutate(deforlag = lag(defor))%>%
    ungroup%>%
    mutate(forshare = 1 - defor,
           forsharelag = 1 - deforlag,
           deforrate = (forsharelag - forshare) / forsharelag)%>% 
    filter_all(all_vars(!is.infinite(.)))%>%
    drop_na(deforrate)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### survival data frame
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
  
  
  outputs = list(
    "countylevel_df" = countylevel_df,
    "proplevel_df" = proplevel_df,
    "panels" = panels,
    "surv_df" = surv_df,
    "grid_df_list" = grid_df_list
  )
  
  return(outputs)
  
}