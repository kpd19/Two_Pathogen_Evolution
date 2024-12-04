# Summary

This library contains code used in the manuscript titled Synthesizing Selection-Mosaic Theory and Host-Pathogen Theory to Explain Large-Scale Pathogen Coexistence in an Insect-Virus System. This paper uses data collected from field experiments and observational studies to explore the mechanisms supporting the coexistence of two viral morphotypes across populations of the Douglas-fir tussock moth. The model and analysis are implemented in R and Julia. 

# Requirements and Setup

The code was built using R version 4.3.2 and Julia version 1.7.2.  R can be downloaded [here](https://www.r-project.org) and Julia can be downloaded [here](https://julialang.org/downloads/). The code requires several packages that are not part of the base R installations. After installing R, navigate to the main repository directory and run the installation.R script to install necessary packages. After installing Julia, run the installation.jl script to install necessary packages. 

# Fitting Bayesian model to field transmission experiment results

In the `TransmissionExperiment` directory, the script `fit_transmission_models.R` uses the Stan programming language to fit Bayesian hierarchical models to the data from our field experiments, contained in `STAN_input_data.csv`. Examples of the output can be found in `output` and `figures`. Each Stan model should take approximately 1 minute to run. The script calculates the leave-one-out cross-validation (LOO-CV) criterion and sum of squared errors (SSE) for each model. The script also plots the model and the data for the best model determined by LOO-CV and SSEs at the morphotype and isolate level. 

The script `thinned_posterior.R` draws 225 parameters from the posterior distribution for the best model. The script `noevo_posteriors` calculates the fraction infected at the end of a single epizootic over a range of initial host densities for the single-pathogen model without host evolution. The fraction infected results are analyzed in `thinned_posterior.R`.
 
# Morphotype frequency and forest composition from Forest Inventory Analysis

In the `MorphotypeFrequency` directory, the spatial dataset on morphotype frequency data from the literature and our field collections was combined in the `aggregating_morphotype_data.R` script in R. We tabulated the total number of isolates identified as the multi-capsid morphotype or the single-capsid morphotype as a result of our PCR analyses. The code includes analysis of coinfections. 

To calculate the percent of the forest that was made up of each tree species, we used the National Forest Type Dataset credated by the Forest Service Inventory and Analysis (FIA) program. The raster data can be downloaded [here](https://data.fs.usda.gov/geodata/rastergateway/forest_type/). In the script `forest_composition.R`, for the United States sites, we calculated the percent of each tree species in a 5 kilometer radius around each field site. For sites that had no trees identified in the FIA dataset, we extended the radius to 10 kilometers. The script then bins the percent that is composed of Douglas-fir to the size of spatial grid used in our mechanistic model, which has 37 patches. 

The script `morph_douglas_glm.R` compares a generalized linear model there the percent of the multi-capsid morphotype is a function of percent Douglas-fir at each site to a model where there is no relationship to forest composition using AIC analysis. 

The script `map_figure.R` visualizes the morphotype distribution on a map with distributions of Douglas-fir and *Abies* species using shapefiles.

# Running the parallelized fitting routine with Message Passing

The `ModelFitting` directory contains the scripts for running the line search routine and increased realizations for the evolution model with a selection mosaic, the model without the selection mosaic, and the model without host evolution. These scripts are highly parallelized using the Message Passing Interface (MPI) and are only designed to be executed on a computing cluster. The realization scripts have been adapted to run only one value for percent Douglas-fir. To parallelize and run as a batch on slurm, the line that says `idx = 18` can be changed to `idx = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])`. The scripts have also been adapted to run only one parameter set for 10 realizations as an example, instead of the top 30 parameter sets for 2000 realizations.  Examples of the line search output for evolution model with a selection mosaic are in `linesearch_op/op1`: the output for all values tried are in `tried_morep1.csv` and the best recorded likelihood score as each round progresses is `best_morep1.csv`. Examples for 10 realizations from the realization output are in `realization_op/op1/`: the likelihood values for all the data points that correspond to forests of that percent Douglas-fir are `ll/ndoug19_p1_sig0.5.csv` and the population values over time are in `all/ndoug19_p1sig0.5.csv`. The naming convention for this example indicates this is the first parameter set for the top 30 parameter sets from the line search, epsilon = 0.5 for environmental stochasitcity, and the model is performed on a hexagonal grid where is 19 out of 37 trees are Douglas-fir. The distances for all the grid coordinates for hexagonal grids with a radius of 2-8 are included in `data`. Examples for the evolution model without a selection mosaic are in `realization_op/op_nt1/` and examples for the model without host evolution are in `realization_op/ne1/` The top 30 parameter sets for each model as a result of the line search routines are in the `data` directory. 

# Model comparison

The `ModelComparison` directory contains the script `model_comparison.R`, which compares the results of all models tested using AIC analysis, which are reported in `data/aic_table.csv`. The likelihood values for the parameter sets for each model are contained in the `data` directory. The `data` directory also contains the results of the fraction infected by the multi-capsid morphotype for each parameter set for each model. The script plots the percent multi-capsid morphotype as a function of percent Douglas-fir, which are shown in `figures` directory.

