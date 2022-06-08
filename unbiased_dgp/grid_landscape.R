# function to generate generic landscape
# pixel value realizations can then be simulated afterward
library(sf)
library(ggplot2)
library(rlist)
library(tidyverse)

grid_landscape <- function(nobs, cellsize){
  
  
  rootn <- ceiling(sqrt(nobs))
  landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))
  
  overgrid <- st_make_grid(landscape, cellsize, square = TRUE)
  
  #trim grid to landscape
  overgrid <- st_intersection(overgrid, landscape)
  
  
  #generate pixels within landscape
  DT <- data.frame(
    pixels= seq(from = 1, to = rootn^2),
    longitude= (rep(seq(from = 1, to = rootn), rootn) - .33),
    latitude= (rep(1:rootn, each= rootn) - .33)
  ) 
  
  pixloc_df <- st_as_sf(DT, coords = c("longitude", "latitude"))
  
  
  #determine which pixels are in each grid 
  wgrid <- st_within(pixloc_df, overgrid, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$grid <- max.col(wgrid)
  
  #determine treated vs. untreated grids
  treat_grids <- sample(1:length(overgrid), (length(overgrid)/2) )
  grid= seq(from = 1, to = length(overgrid))
  treatgrid <- data.frame(
    grid= seq(from = 1, to = length(overgrid)),
    treat = ifelse(grid %in% treat_grids, 1, 0)
  ) 
  
  # determine which pixels are treated vs. untreated
  pixloc_df <- merge(pixloc_df, treatgrid, by = "grid")
  gridcoords_df <- data.frame(overgrid)
  
  outputs = list('pixloc_df' = pixloc_df, 'gridcoords' = gridcoords_df)
  return(outputs)
  
}
