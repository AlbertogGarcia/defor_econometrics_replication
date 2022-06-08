library(DeclareDesign)
library(reshape2)  
deforestation_DGP <- function(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v){
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2)^.5)
  
  panels <- fabricate(
    pixels = add_level(N = nobs, a_i = rnorm(N, 0, std_a), treat = rbinom(N, 1, 0.5)),
    year = add_level(N = (years*2), nest = FALSE),
    obs = cross_levels(
      by = join(pixels, year),
      v_it = rnorm(N, 0, std_v)
      )
    )%>%
    mutate(pixels = as.numeric(pixels),
           year = as.numeric(year),
           post = ifelse(year > years, 1, 0),
           ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
           y = (ystar > 0)*1,
           ystar_cf = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it,
           y_cf = (ystar_cf > 0)*1,
           defor_indic = ifelse(y==1, year, 99),
           defor_indic_cf = ifelse(y_cf==1, year, 99))%>%
      group_by(pixels)%>%
      mutate(defor_year = min(defor_indic),
             defor_year_cf = min(defor_indic_cf),
             defor = (year>=defor_year)*1,
             y_it = ifelse(year>defor_year, NA, defor),
             defor_cf = (year>=defor_year_cf)*1,
             y_it_counter = ifelse(year>defor_year_cf, NA, defor_cf)
      )
    
    
    outputs = list("panels" = panels, "ATT" = ATT)
    # assign('panels',panels, envir=.GlobalEnv)
    # assign('ATT',ATT, envir=.GlobalEnv)
    return(outputs)
  }



