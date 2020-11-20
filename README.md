# GLMM-small-sample-adjusted-CRT

This repo contains electronic materials for "Evaluating tests for cluster-randomized trials with few clusters under generalized linear mixed models with covariate adjustment: a simulation study".

## Electronic tables of simulation results

The electronic tables of simulation results are held at https://QIU-Hongxiang-David.github.io/small-sample-adjusted-GLMM-CRT. The raw data are available at [simulation1_typeIerror.csv](https://github.com/QIU-Hongxiang-David/small-sample-adjusted-GLMM-CRT/simulation1_typeIerror.csv) and [simulation2_typeIerror.csv](https://github.com/QIU-Hongxiang-David/small-sample-adjusted-GLMM-CRT/simulation2_typeIerror.csv).


## R code for simulations

The R code for the two simulations are under `simulation1/` and `simulation2/`. Under each folder, `generate_data.R` contains functions to generate data; `simulation.R` contains code to call functions in `generate_data.R` and run the simulation.
