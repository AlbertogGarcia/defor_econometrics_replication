Overview
--------

The code in this replication package constructs the analysis from Garcia and Heilmayr, 2022. Three main files run all of the code to generate the data for the 1 figure and 16 results files needed to run the R markdown file to generate the paper. 

Data Availability and Provenance Statements
----------------------------

This paper does not involve analysis of external data. The only data are generated via simulation).


Instructions to Replicators
---------------------------

- Run `unbiased_dgp/map_figures.R`. Image file `landscape_map.png` will be output to `unbiased_dgp/figs`. Copy `unbiased_dgp/figs` into `paper/figs`. This is the lone figure used by `paper/defor_metrics_draft.Rmd`.

- Run `unbiased_dgp/analysis_main.R`. Results files will be output to `unbiased_dgp/results`. Copy `unbiased_dgp/results` to `paper/results`.

- Run `multigroup_dgp/analysis_multiple_gt.R`. Results files will be output to `unbiased_dgp/results_multi`. Copy `unbiased_dgp/results_multi` to `paper/results_multi`.

- Replicators should now be able to knit `paper/defor_metrics_draft.Rmd` into the paper pdf. 