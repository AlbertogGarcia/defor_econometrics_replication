# function to generate generic landscape
# pixel value realizations can then be simulated afterward
library(sf)
library(rlist)
library(tidyverse)
library(spatstat)
#library(ggpattern)

multigrid_landscape <- function(nobs, cellsize_list, ppoints, cpoints){
  
  
  rootn <- ceiling(sqrt(nobs))
  landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))
  
  #generate pixels within landscape
  DT <- data.frame(
    pixels= seq(from = 1, to = rootn^2),
    longitude= (rep(seq(from = 1, to = rootn), rootn) - .33),
    latitude= (rep(1:rootn, each= rootn) - .33)
  ) 
  
  #create sf object by assigning long lat coordinates as geometry
  pixloc_df <- st_as_sf(DT, coords = c("longitude", "latitude"))
  
  #generate voronoi pts for 
  vorpts_prop <- st_sample(landscape, type = "unifpoint", n = ppoints)
  vorpts_county <- st_sample(landscape, type = "unifpoint", n = cpoints)
  
  v_county <- vorpts_county %>%  # consider the sampled points
    st_geometry() %>% #  as geometry only 
    st_union() %>% # unite them 
    st_voronoi() %>% # perform the voronoi tessellation
    st_collection_extract(type = "POLYGON") %>%
    st_intersection(landscape)
  
  v_property <- vorpts_prop %>%  # consider the sampled points
    st_geometry() %>% #  as geometry only 
    st_union() %>% # unite them 
    st_voronoi() %>% # perform the voronoi tessellation
    st_collection_extract(type = "POLYGON") %>% # select the polygons
    st_intersection(landscape)  # limit to within landscape boundaries
  
  pixloc_grid <- pixloc_df %>%
    st_set_geometry(NULL)
  
  for(i in cellsize_list){
    
    
    overgrid <- st_make_grid(landscape, i, square = TRUE)
    
    #trim grid to landscape
    overgrid <- st_intersection(overgrid, landscape)
    
    wgrid <- st_within(pixloc_df, overgrid, sparse = FALSE, prepared = TRUE)*1
    
    pixloc_df[ , ncol(pixloc_df) + 1] <- as.character(max.col(wgrid))                  # Append new column
    colnames(pixloc_df)[ncol(pixloc_df)] <- paste0("grid_", i)  # Rename column name
    
    gareas <- tibble::rownames_to_column(
        data.frame(matrix(unlist(st_area(overgrid))))
        )
    
    names(gareas)[1] <- paste0("grid_", i)
    names(gareas)[2] <- paste0("garea_", i)
    
    
    pixloc_df <- pixloc_df %>%
      inner_join(gareas, by = paste0("grid_", i)
      )
    
  }
  
  wprop <- st_within(pixloc_df, v_property, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$property <- max.col(wprop)
  
  wcounty <- st_within(pixloc_df, v_county, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$county <- max.col(wcounty)
  
  #determine treated vs. untreated counties
  treat_counties <- sample(1:length(v_county), (length(v_county)/2) )
  county= seq(from = 1, to = length(v_county))
  
  treatcounty <- data.frame(
    county= seq(from = 1, to = length(v_county)),
    treat = ifelse(county %in% treat_counties, 1, 0)
  ) 
  
  
  # determine which pixels are treated vs. untreated
  pixloc_df <- pixloc_df %>%
    inner_join(treatcounty, by = "county")
  
  #getting areas for properies and counties
  pareas <- data.frame(matrix(unlist(st_area(v_property))))
  pareas <- tibble::rownames_to_column(pareas)
  names(pareas)[1] <- paste("property")
  names(pareas)[2] <- paste("parea")
  careas <- data.frame(matrix(unlist(st_area(v_county))))
  careas <- tibble::rownames_to_column(careas)
  names(careas)[1] <- paste("county")
  names(careas)[2] <- paste("carea")
  
  
  
  #sapply(pixloc_df, class)
  pixloc_df <- pixloc_df %>%
    mutate_at(vars(pixels, county, property), as.character)
  
  #merging back areas
  
  pixloc_df <- pixloc_df %>%
    inner_join(careas, by = "county") %>%
    inner_join(pareas, by = "property")
  
  intervention_plot <- st_as_sf(v_county) %>%
    tibble::rownames_to_column("county") 
  intervention_plot$county <- as.integer(intervention_plot$county)
  intervention_plot <- intervention_plot %>%
    inner_join(treatcounty, by = "county") 
  intervention_plot$treat <- as.factor(intervention_plot$treat)
  
  intervention_area <- intervention_plot %>%
    filter( treat == 1)
  control_area <- intervention_plot %>%
    filter( treat == 0)
  
  outputs = list('pixloc_df' = pixloc_df, 'control_area' = control_area, 'intervention_area' = intervention_area, "p_bounds" = v_property, "c_bounds" = v_county)
  return(outputs)
  
}
