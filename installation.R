################################################
## Install packages for the analysis 
################################################

## Install multiple packages
install.packages(c("tidyverse","gridExtra","raster","sp",'rstan','loo','bayesplot','posterior','sf',
                   'geosphere','geodata','ggnewscale','scatterpie','ggtext','WaveletComp', 'gg3D',
                   'mgcv'),
                dependencies = T)
