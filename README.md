# Individual herbivores take risks based on resource quality: stoichiometric distribution models with snowshoe hares

### Isabella C. Richmond; Juliana Balluffi-Fry; Eric Vander Wal; Shawn J. Leroux; Matteo Rizzuto; Travis R. Heckford; Joanie L. Kennah; Gabrielle R. Riefesel; Yolanda F. Wiersma

#### Manuscript is being prepared for submission at the Journal of Mammalogy. All data and code for the manuscript are available here and published on Zenodo.

This repository contains code, data, and results for Richmond et al. Individual herbivores take risks based on resource quality: stoichiometric distribution models with snowshoe hares. 

Zenodo Repository:

**script**  
This folder contains all the R scripts necessary to reproduce our analyses.  

* 1 - DataTriangulationRazimuth.R takes raw VHF relocation data and triangulates them with error ellipses for each collar and date using the package ```razimuth```.
* 2 - kUDRazimuth.R uses triangulated data to calculate kernel utilization distributions for each individual using ```adehabitatHR```. These kernel utilization distributions are NOT used in future analysis.
* 3 - aKDERazimuth.R uses triangulated data to calculate autocorrelated kernel density estimates for each individual using ```ctmm```. These kernels are used for all future analysis.
* 4 - HorizontalComplexity.R processes raw horizontal complexity data collected in the field and produces finalized versions for analysis.  
* 5 - RiskOrdination.R takes all habitat complexity data and analyzes understory and overstory habitat complexity using Principal Components Analyses. PCA axes are then extracted for future analysis.
* 6 - HomeRangeExtraction.R extracts kernel utilization distribution data for each individual hare, food quality data, and habitat complexity data at each habitat complexity sampling point.
* 7a - lme4Models.R models our data with intensity of space use (kernel utilization distribution) as a response variable and habitat complexity and food quality as explanatory variables using linear mixed effects models in the ```lme4``` package.
* 7b - MCMCglmmModels.R models our data using a Bayesian approach with intensity of space use (kernel utilization distribution) as a response variable and habitat complexity and food quality as explanatory variables using a Markov chain Monte Carlo sampler for linear mixed effects models in the ```MCMCglmm``` package. 
* 7c - MCMCglmmPlots.R plots of random slopes and correlations used in the manuscript from the Bayesian models produced using ```MCMCglmm``` in 7b - MCMCglmmModels.R.
* 8 - Maps.R - plots maps used in the manuscript.
* function-plotOutliers.R -a function used to plot outliers in 1 - DataTriangulationRazimuth.R 
* function-plotVariograms.R - a function used to plot variograms in 3 - aKDERazimuth.R.
* function-residPlots.R - a function used to produce diagnostic plots for Frequentist models in 7a- lme4Models.R.

**graphics**  
This folder contains all the graphics and figures produced during the analysis.  

 
**input**  
This folder and all folders within contain the raw data used in analysis.

**output**  
This folder contains all the output from our analysis.  

* outliers
  + Contains outlier plots from triangulation of VHF data in 1 - DataTriangulationRazimuth.R 
* variograms
  + Contains variograms from movement models produced using ctmm in 3 - aKDERazimuth.R