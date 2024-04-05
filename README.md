# BAM Landbird Models

## A generalized modeling framework for spatially extensive species abundance prediction and population estimation

BAM's landbird models bridge the gap between local studies and large-scale management needs by compiling and harmonizing data from many soruces to predict avian abundance at a fine resolution and broad extent. The modelling workflow is broken up into eleven steps:

1. Download from [WildTrax](wildtrax.ca) and harmonization
2. Calculate of statistical offsets with the [`QPAD`](github.com/borealbirds/QPAD) package
3. Extract covariates
4. Stratify dataset by region and year
5. Gap analysis to identify areas of extrapolation
6. Tune models for learning rate
7. Bootstrap models
8. Spatial prediction
9. Mosaic and average predictions
10. Validate models
11. Interpret models and predictions

## Version 4

Version 4 of BAM's landbird models can be reproduced using the data object stored in the data repository and the generalized model script available in the 'V4' folder of this repository. 

## Version 5

This repository contains code for Version 5 of the workflow, including GIS (gis) and the modelling workflow (analysis), and is currently under active development.

New developments to BAM's Landbird Models for Version 5 include:

* Expanded study area to include hemi-boreal regions of the continental United States using natural biographic boundaries rather than political boundaries to delineate regions
* Expanded list of environmental covariates including time-matched variables for vegetation biomass, human disturbance, and annual climate
* Acquisition and inclusion of new datsets for data-sparse regions, including appropriate eBird data
* Predictions for five year intervals from 1985 to 2020
* More rigorous model tuning and covariate thinning
