library(tidyverse)
library(here)
library(data.table)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameterization 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

years = 6 # years in each period (total years = years*2)
nobs = 150^2
n = 500 # number of monte carlo iterations


cellsize_p = 10
cellsize_c = 30
cellsize_list <- c(cellsize_p, cellsize_c)
ppoints = 225
cpoints = 25

avg_parea  = nobs/ppoints
avg_carea = nobs/cpoints

std_v = 0.5
std_a = 0
std_p = 0

# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

###############################################################################################################
######## Baseline showing aggregation resolves pixel fixed effects issue
###############################################################################################################

source(here::here('unbiased_dgp', 'montecarlo_all_specifications.R'))

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)


set.seed(0930)
# summary function that estimates all of the different specifications
aggregation <- montecarlo_all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, 
                                             cellsize_list, ppoints, cpoints, 
                                             min_psize = years*2,
                                             nestedprops = FALSE, proptreatassign = FALSE)

library(rio)
export(aggregation, here::here("paper", "results", "results_aggregation.rds"))


###############################################################################################################
######## Introduction pixel level unobservables, which impact non-random selection
###############################################################################################################
std_a = 0.1
# we'll need to recompute the parameters if we change the value of sigma_p

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

# summary function that estimates all of the different specifications
aggregation_0 <- montecarlo_all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p,
                                               cellsize_list, ppoints, cpoints, 
                                               min_psize = years*2,
                                               nestedprops = FALSE, proptreatassign = FALSE)

export(aggregation_0, here::here("paper", "results", "results_selection.rds"))
###############################################################################################################
######## Adding in property level disturbances
###############################################################################################################

#### 0.1
std_p = 0.1

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_1 <- montecarlo_all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p,
                                               cellsize_list, ppoints, cpoints, 
                                               min_psize = years*2,
                                               nestedprops = FALSE, proptreatassign = FALSE)

##### 0.2
std_p = 0.2

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

#ATT = pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 )^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 )^.5)

aggregation_2 <- montecarlo_all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p,
                                               cellsize_list, ppoints, cpoints, 
                                               min_psize = years*2,
                                               nestedprops = FALSE, proptreatassign = FALSE)


#### 0.3
std_p = 0.3

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_3 <- montecarlo_all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p,
                                               cellsize_list, ppoints, cpoints, 
                                               min_psize = years*2,
                                               nestedprops = FALSE, proptreatassign = FALSE)

results_full <- rbind(aggregation_0, aggregation_1, aggregation_2, aggregation_3)

export(results_full, here::here("paper", "results", "results_full.rds"))


###############################################################################################################
######## alternative parameterization
###############################################################################################################

std_p = 0.3
std_a = 0.1

base_0 = .05
base_1 = .02
trend = -.005
ATT = -.01

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)


aggregation_alt <- montecarlo_all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p,
                                                 cellsize_list, ppoints, cpoints, 
                                                 min_psize = years*2,
                                                 nestedprops = FALSE, proptreatassign = FALSE)

export(aggregation_alt, here::here("paper", "results", "results_alt.rds"))

base_0 = .02
base_1 = .05
trend = .005
ATT = .01

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_alt2 <- montecarlo_all_specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p,
                                                  cellsize_list, ppoints, cpoints, 
                                                  min_psize = years*2,
                                                  nestedprops = FALSE, proptreatassign = FALSE)

export(aggregation_alt2, here::here("paper", "results", "results_alt2.rds"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### weighting analysis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(0930) # re-setting seed, so one can split runs up

source(here::here('unbiased_dgp', 'heterogeneous_propertyarea.R'))

# we start with our base parameterization without property level perturbations
std_a = 0
std_v = 0.5
std_p = 0.0
std_b3 = .05
ppoints = 225

# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01


# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1


cellsize = cellsize_p

weights <- heterogeneous_propertyarea(n, nobs, years, b0, b1, b2_0, b2_1, std_a, std_v, std_p, std_b3, given_ATT = ATT, cellsize, ppoints, cpoints)
results_pweights <- weights$summary_long

library(rio)
export(results_pweights, here::here("paper", "results", "results_pweights.rds"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### outcome analysis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source(here::here('unbiased_dgp', 'outcome_fcn.R'))

# we start with our base parameterization without property level perturbations
std_a = 0
std_v = 0.5
std_p = 0
# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

set.seed(0930)
outcomes <- outcome_fcn(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, cellsize)
outcome <- outcomes$coeff_bias

library(rio)
export(outcome, here::here("paper", "results", "results_outcomeFormula.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### DID keeping vs. dropping obs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source(here::here('unbiased_dgp', 'DID_keep.R'))

set.seed(0930)
keeps <- DID_keep(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
keeps <- keeps$did_keeps

library(rio)
export(keeps, here::here("paper", "results", "results_pixKeep.rds"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######## show TWFE is equivalent to dropping all pixels deforested in first period
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source(here::here('unbiased_dgp', 'TWFE_expost.R'))

estimator_comp <- TWFE_expost(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)

summary_coeff <- estimator_comp$summary_long %>%
  mutate_at(vars(bias), as.numeric)

summary_wide  <- summary_coeff %>%
  group_by(model)%>%
  summarise(RMSE = rmse(bias, 0),
            q25 = quantile(bias, probs = .25),
            q75 = quantile(bias, probs = .75),
            Bias = mean(bias))

export(summary_wide, here::here("paper", "results", "twfe_comp_summary.rds"))
