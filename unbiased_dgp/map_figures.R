#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Alberto Garcia and Robert Heilmayr
# Paper: Conservation Impact Evaluation Using Remotely Sensed Data
# Date: 6/12/22
# Purpose: Create map of example landscape (Figure 1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import packages --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(patchwork)
library(tidyverse)
library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(DeclareDesign)
library(data.table)
source(here::here("unbiased_dgp", "full_landscape.R"))
library(ggpubr)

set.seed(6226)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set parameters --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format_fig <- function(figure){
  figure <- figure +
    theme(panel.background = element_rect(fill = "white",colour = "white"),
          plot.background = element_rect(fill = "white",colour = "white"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.background=element_rect(fill = "white"),
          legend.key=element_rect(fill = "white"),
          legend.text=element_text(size=12),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0,  # Left margin
                               unit = "in")
          )
  return(figure)
}

palette <- list("white" = "white",
                "light_grey" = "#d9d9d9",
                "dark" = "#0c2230",
                "red" = "#ed195a",
                "blue" = "#1c86ee",
                "green" = "#7CAE7A",
                "dark_green" = "#496F5D",
                "light_blue" = "#56B4E9",
                "gold" = "#DAA520")

base_0 = .04
base_1 = .08
trend = 0.02
ATT = -0.07

std_a = 0.1
std_v = 0.5
years = 1
std_p = 0.2

main_nobs = 150^2
main_ppoints = 225
main_cpoints = 25
csize = main_nobs/main_cpoints
psize = main_nobs/main_ppoints

nobs = 75^2
ppoints = floor(nobs/psize)+1
cpoints = floor(nobs/csize)+1

cellsize = 12

nobs_ratio = main_nobs/ nobs

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)


countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
pixloc_df = countyscape$pixloc_df
control_area = countyscape$control_area
intervention_area = countyscape$intervention_area
intervention_area_merge = intervention_area %>% st_union()
p_bounds = countyscape$p_bounds
c_bounds = countyscape$c_bounds


ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)

pixloc <- pixloc_df

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
    ystar_counterfactual = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
  )
)


#generate random 
error_table <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, sd = std_p))

panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)

panels <- panels %>%
  inner_join(pixloc, by = c("pixels", "treat")) %>%
  inner_join(error_table, by = "property") %>%
  mutate(ystar = ystar + p_err, 
         y = (ystar > 0)*1, 
         ystar_counterfactual = ystar_counterfactual+p_err,
         y_counterfactual = (ystar_counterfactual > 0)*1) 

panels_counterfactual <- panels

#need to determine which year deforestation occurred
year_df <- panels %>%
  dplyr::select(pixels, year, y) %>%
  dcast(pixels ~ year , value.var = "y")

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
  inner_join(panels, by = "pixels")

cols.num <- c("pixels", "grid", "property", "county", "year")
panels[cols.num] <- sapply(panels[cols.num],as.numeric)

panels <- panels %>%
  mutate(indic = year - defor_year,
         defor = ifelse(indic > 0, 1, y),
         y_it = ifelse(indic > 0, NA, y) )



data_df <- panels %>%
  mutate(y_it = replace_na(y_it, -1),
         y_plot = ifelse(y_counterfactual != y, 2, y_it),
         y_it = y_it %>% recode("-1" = "Previously deforested", 
                                "0" = "Not deforested", 
                                "1" = "Deforested"),
         y_plot = y_plot %>% recode("-1" = "Deforested", 
                                "0" = "Not deforested", 
                                "1" = "Deforested",
                                "2" = "Counterfactual deforestation"),
         y_counterfactual = y_counterfactual %>% recode("-1" = "Previously deforested", 
                                                        "0" = "Not deforested", 
                                                        "1" = "Deforested"),
         y = y %>% recode("-1" = "Previously deforested", 
                          "0" = "Not deforested", 
                          "1" = "Deforested"),
         treat = as.factor(treat),
         treat = treat %>% recode("0" = "Stable forest - not treated",
                                  "1" = "Stable forest - treated")) %>% 
  st_as_sf() %>%
  dplyr::select(pixels, year, treat, defor, y_it, defor_year, y, y_counterfactual, y_plot) 





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Landscape plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fills = c(#"Previously deforested" = palette$light_grey,
          #"Not deforested" = palette$dark, 
          "Deforested" = palette$red,
          "Stable forest - treated" = palette$dark,
          "Stable forest - not treated" = palette$dark_green,
          "Counterfactual deforestation" = palette$light_blue,
          "County boundaries" = palette$gold,
          "Property boundaries" = "grey65",
          "Grid cell boundaries" = "white")

colors = c("County boundaries" = palette$green,
  "Property boundaries" = "grey65",
  "Grid cell boundaries" = "white")




plot_df <- data_df %>% 
  filter(year==2) %>% 
  mutate(plot_var = ifelse(y_plot %in% c("Deforested", "Counterfactual deforestation"), as.character(y_plot), as.character(treat)))


# County-level analysis
county_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 1.6) +
  geom_sf(data = p_bounds, size = .75, fill = "NA", color = "grey60") +
  geom_sf(data = c_bounds, size = 1.1, fill = "NA", color = palette$gold) +
  geom_vline(xintercept = c(seq(from = cellsize, to = sqrt(nobs), by = cellsize)), color="white", size=1.1) + 
  geom_hline(yintercept = c(seq(from = cellsize, to = sqrt(nobs), by = cellsize)), color="white", size=1.1) +
  scale_fill_manual(values = fills)+
  guides(fill = guide_legend(override.aes = list(size=5)))+
  xlab("")+ylab("")
county_plot %>% 
  format_fig()

ggsave(county_plot %>% 
         format_fig(),
       path = "paper/figs",
       filename = "landscape_map.png", 
       width = 7,
       height = 4)

