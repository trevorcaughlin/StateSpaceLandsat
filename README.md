# StateSpaceLandsat
This repository contains the code necessary to run the state-space models from the paper titled, "Monitoring tropical forest succession at landscape scales despite uncertainty in Landsat time series."
The .STAN files correspond to the four state-space models mentioned in the text.
"BASIC.STAN" is the model with no covariates in the measurement error term.
"MM.STAN" is the intra-annual phenology model, with date-of-acquisition as a covariate for measurement error.
"YEAR.STAN" is the inter- and intra-annual phenology model with date-of-acquisition and year as covariates for measurement error.
"HT.STAN" is the canopy height model with date-of-acquisition, year, and Landsat-lidar fusion as components of measurement error.

The file titled "running models.R" is the code necessary to run the Stan models, and "landsat_lidar_ndvi.csv," "TEN_squares.Rdata," and "prior_list.Rdata,"
represent inputs to the Stan models.

For more information, email trevorcaughlin@boisestate.edu.
