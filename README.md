Overview
--------

The code in this replication package reproduces all results from Garcia and Heilmayr, 2022. Three main scripts generate all input files (1 figure and 16 results files) needed to knit the paper (`paper/defor_metrics_draft.Rmd`). 


Data Availability and Provenance Statements
----------------------------

This paper does not involve analysis of external data. All data are generated via simulation.



Instructions for replication
---------------------------

- Set up compute environment: Open `defor_econometrics_replication.Rproj`. Renv will be automatically installed. Run renv::restore() to pull package versions used for original paper production.

- Run `unbiased_dgp/map_figures.R`. Image file `landscape_map.png` will be output to `unbiased_dgp/figs`. This is the lone figure used by `paper/defor_metrics_draft.Rmd`.

- Run `unbiased_dgp/analysis_main.R`. Results files will be output to `paper/results`.

- Run `multigroup_dgp/analysis_multiple_gt.R`. Results files will be output to `paper/results_multi`.

- You should now be able to knit `paper/defor_metrics_draft.Rmd` into the paper pdf. 