source(here::here('multigroup_dgp', 'multipleGT.R'))
source(here::here('multigroup_dgp', 'multipleGT_agg.R'))
source(here::here('multigroup_dgp', 'multipleGT_pix.R'))
source(here::here('multigroup_dgp', 'trends_fcn.R'))
library(rio)
base_a = .07
base_b = .05
base_c = .02
trend1 = .00
trend2 = .00
trend3 = .00

dyn_ATT_a = 0
dyn_ATT_b = 0
ATT_a = -0.02
ATT_b =  -0.02

std_a = 0.0
std_v = 0.5
std_p = 0.0

cpoints = 30
cellsize=30 
ppoints=100

nobs = 100^2
n=200

set.seed(0930)

multiGT_agg <- multipleGT_agg(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize, ppoints, cpoints)

export(multiGT_agg$es_long, "multigroup_dgp/results_multi/county_long.rds")

county_es <- multiGT_agg$es_long %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

export(county_es, "multigroup_dgp/results_multi/county_es.rds")


multiGT_pix <- multipleGT_pix(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize=10, ppoints=50, cpoints)

export(multiGT_pix$es_long, "multigroup_dgp/results_multi/pixel_long.rds")

pixel_es <- multiGT_pix$es_long %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05, na.rm = T),
            q95 = quantile(estimate, probs = 0.95, na.rm = T),
            estimate = mean(estimate, na.rm = T))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

export(pixel_es, "multigroup_dgp/results_multi/pixel_es.rds")



set.seed(0930)
landscape <- trends_fcn(500000, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize, ppoints, cpoints)

export(landscape$trends_df, "multigroup_dgp/results_multi/landscape.rds")
  
##########################################################################
##########################################################################
dyn_ATT_a = 0.01
dyn_ATT_b =  -0.02
ATT_a = -0.03
ATT_b =  0.02

set.seed(0930)

multiGT_agg <- multipleGT_agg(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize=10, ppoints=50, cpoints)
library(rio)
export(multiGT_pix$es_long, "multigroup_dgp/results_multi/county_long_hetTE.rds")

county_es <- multiGT_agg$es_long %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

export(county_es, "multigroup_dgp/results_multi/county_es_hetTE.rds")


my_event_study_plot(county_es, seperate = FALSE)+
  ggtitle("estimates with aggregated unit of analysis (county)")+
  geom_segment(aes(x = -2.5, y = 0, xend = -0.5, yend = 0), color = "limegreen")
  #geom_segment(aes(x = -0.5, y = -0.02, xend = 2.5, yend = -0.02), color = "limegreen")


multiGT_pix <- multipleGT_pix(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize=10, ppoints=50, cpoints)

export(multiGT_pix$es_long, "multigroup_dgp/results_multi/pixel_long_hetTE.rds")

pixel_es <- multiGT_pix$es_long %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05, na.rm = T),
            q95 = quantile(estimate, probs = 0.95, na.rm = T),
            estimate = mean(estimate, na.rm = T))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

export(pixel_es, "multigroup_dgp/results_multi/pixel_es_heterogTE.rds")


set.seed(0930)
landscape_het <- trends_fcn(500000, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize, ppoints, cpoints)

export(landscape_het$trends_df, "multigroup_dgp/results_multi/het_landscape.rds")