#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Alberto Garcia and Robert Heilmayr
# Paper: Conservation Impact Evaluation Using Remotely Sensed Data
# Date: 6/12/22
# Purpose: Workflow to generate primary paper results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import packages --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(tictoc)
library(here)
library(DeclareDesign)
library(survival)
library(ggplot2)
library(dplyr)
library(ggfortify)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameterization 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

years = 6
nobs = 150^2
n = 500

cellsize_small = 5
cellsize_med = 10
cellsize_large = 30
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

source(here::here('unbiased_dgp', 'specifications.R'))

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

set.seed(0930)
# summary function that estimates all of the different specifications
aggregation <- specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE)

summary_long <- aggregation$summary_long

library(rio)
export(summary_long, "unbiased_dgp/results/summary_long.rds")

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
aggregation_0 <- specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE)
summary_long_0 <- aggregation_0$summary_long 

export(summary_long_0, "unbiased_dgp/results/summary_selection.rds")
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

aggregation_1 <- specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE)
summary_long_1 <- aggregation_1$summary_long

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

aggregation_2 <- specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE)
summary_long_2 <- aggregation_2$summary_long


#### 0.3
std_p = 0.3

std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

aggregation_3 <- specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE)
summary_long_3 <- aggregation_3$summary_long

summary_full <- rbind(summary_long_0, summary_long_1, summary_long_2, summary_long_3)

export(summary_full, "unbiased_dgp/results/summary_full.rds")


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

set.seed(0930)
cellsize = cellsize_med

aggregation_alt <- specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE)

summary_long_alt <- aggregation_alt$summary_long
export(summary_long_alt, "unbiased_dgp/results/summary_long_alt.rds")

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

aggregation_alt2 <- specifications(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, std_c = 0.0, cellsize_small, cellsize_med, cellsize_large, ppoints, cpoints, nestedprops = FALSE, proptreatassign = FALSE)

summary_long_alt2 <- aggregation_alt2$summary_long
export(summary_long_alt2, "unbiased_dgp/results/summary_long_alt2.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### weighting analysis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

set.seed(0930)
weights <- heterogeneous_propertyarea(n, nobs, years, b0, b1, b2_0, b2_1, std_a, std_v, std_p, std_b3, given_ATT = ATT, cellsize = 10, ppoints, cpoints)
summary_pweights <- weights$summary_long

library(rio)
export(summary_pweights, "unbiased_dgp/results/summary_pweights.rds")

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
outcomes <- outcome_fcn(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, cellsize = 10)
outcome <- outcomes$coeff_bias

library(rio)
export(outcome, "unbiased_dgp/results/outcomes.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### DID keeping vs. dropping obs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source(here::here('unbiased_dgp', 'DID_keep.R'))

set.seed(0930)
keeps <- DID_keep(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
keeps <- keeps$did_keeps

library(rio)
export(keeps, "unbiased_dgp/results/keeps.rds")



source(here::here('unbiased_dgp', 'TWFE_expost.R'))

set.seed(0930)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######## show TWFE is equivalent to dropping all pixels deforested in first period
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimator_comp <- TWFE_expost(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)

summary_coeff <- estimator_comp$summary_long %>%
  mutate_at(vars(bias), as.numeric)

summary_wide  <- summary_coeff %>%
  group_by(model)%>%
  summarise(RMSE = rmse(bias, 0),
            q25 = quantile(bias, probs = .25),
            q75 = quantile(bias, probs = .75),
            Bias = mean(bias))

export(summary_wide, "unbiased_dgp/results/twfe_comp.rds")
