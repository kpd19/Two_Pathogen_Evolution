# Summary

This library contains code used in the manuscript titled Synthesizing Selection-Mosaic Theory and Host-Pathogen Theory to Explain Large-Scale Pathogen Coexistence in an Insect-Virus System. This paper uses data collected from field experiments and observational studies to explore the mechanisms supporting the coexistence of two viral morphotypes across populations of the Douglas-fir tussock moth. The model is implement in R and Julia

# Requirements and Setup

The code was built using R version 4.3.2 and Julia version 1.7.2.  R can be downloaded [here](https://www.r-project.org) and Julia can be downloaded [here](https://julialang.org/downloads/). The code requires several packages that are not part of the base R installations. After installing R, navigate to the main repository directory and run the installation.R script. To run this script from the command line, simply navigate to the directory and execute:

# Fitting Bayesian model to field transmission experiment results
 
# Morphotype frequency and forest composition from Forest Inventory Analysis

The spatial dataset on morphotype frequency data from the literature and our field collections was combined in the `aggregating_morphotype_data.R` script in R. We tabulated the total number of isolates identified as the multi-capsid morphotype or the single-capsid morphotype as a result of our PCR analyses. The code includes analysis of coinfections. 

To calculate the percent of the forest that was made up of each tree species, we used the National Forest Type Dataset credated by the Forest Service Inventory and Analysis (FIA) program. The raster data can be downloaded [here](https://data.fs.usda.gov/geodata/rastergateway/forest_type/). For the United States sites, we calculated the percent of each tree species in a 5 kilometer radius around each field site. For sites that had no trees identified in the FIA dataset, we extended the radius to 10 kilometers. The script then bins the percent that is composed of Douglas-fir to the size of spatial grid used in our mechanistic model, which has 37 patches.  

# Running the parallelized fitting routine with Message Passing

# Model comparison


