# Summary

This library contains code used in the manuscript titled Synthesizing Selection-Mosaic Theory and Host-Pathogen Theory to Explain Large-Scale Pathogen Coexistence in an Insect-Virus System. This paper uses data collected from field experiments and observational studies to explore the mechanisms supporting the coexistence of two viral morphotypes across populations of the Douglas-fir tussock moth. The model is implement in R and Julia. 

# Requirements and Setup

The code was built using R version 4.3.2 and Julia version 1.7.2.  R can be downloaded [here](https://www.r-project.org) and Julia can be downloaded [here](https://julialang.org/downloads/). The code requires several packages that are not part of the base R installations. After installing R, navigate to the main repository directory and run the installation.R script. To run this script from the command line, simply navigate to the directory and execute:

# Fitting Bayesian model to field transmission experiment results

The script `fit_transmission_models.R` uses the Stan programming language to fit Bayesian hierarchical models to the data from our field experiments, contained in `STAN_input_data.csv`. Examples of the output can be found in `output` and `figures`. Each Stan model should take approximately 1 minute to run. The script calculates the leave-one-out cross-validation (LOO-CV) criterion and sum of squared errors (SSE) for each model. The script also plots the model and the data for the best model determined by LOO-CV and SSEs at the morphotype and isolate level. 

The script `thinned_posterior.R` draws 225 parameters from the posterior distribution for the best model. The script `noevo_posteriors` calculates the fraction infected at the end of a single epizootic over a range of initial host densities for the single-pathogen model without host evolution. The fraction infected results are analyzed in `thinned_posterior.R`.
 
# Morphotype frequency and forest composition from Forest Inventory Analysis

The spatial dataset on morphotype frequency data from the literature and our field collections was combined in the `aggregating_morphotype_data.R` script in R. We tabulated the total number of isolates identified as the multi-capsid morphotype or the single-capsid morphotype as a result of our PCR analyses. The code includes analysis of coinfections. 

To calculate the percent of the forest that was made up of each tree species, we used the National Forest Type Dataset credated by the Forest Service Inventory and Analysis (FIA) program. The raster data can be downloaded [here](https://data.fs.usda.gov/geodata/rastergateway/forest_type/). In the script `forest_composition.R`, for the United States sites, we calculated the percent of each tree species in a 5 kilometer radius around each field site. For sites that had no trees identified in the FIA dataset, we extended the radius to 10 kilometers. The script then bins the percent that is composed of Douglas-fir to the size of spatial grid used in our mechanistic model, which has 37 patches. 

The script `morph_douglas_glm.R` compares a generalized linear model there the percent of the multi-capsid morphotype is a function of percent Douglas-fir at each site to a model where there is no relationship to forest composition using AIC analysis. 

The script `map_figure.R` visualizes the morphotype distribution 

# Running the parallelized fitting routine with Message Passing

# Model comparison


