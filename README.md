# Structural nested cumulative failure time models
Author: Joy Shi

Last updated: June 21, 2021

G-estimation of structural nested cumulative failure time models (SNCFTMs) is a method to estimate the causal effect of a time-varying treatment on a failure-time outcome in the presence of informative right-censoring and treatment-confounder feedback. Under certain "no interaction" conditions, g-estimates of SNCFTMs can be used to calculate marginal cumulative risks of the outcome under static treatment regimes. 

R code is provided to g-estimate SNCFTMs based on (1) adjustment of exposure-outcome confounders, or (2) instrumental variable analysis. Additional information about the methods are available at:

1. [Picciotto S, Hernán MA, Page JH, Young JG, Robins JM. Structural nested cumulative failure time models to estimate the effects of interventions. Journal of the American Statistical Association. 2012 Sep 1;107(499):886-900.](https://pubmed.ncbi.nlm.nih.gov/24347749/ "Picciotto S, Hernán MA, Page JH, Young JG, Robins JM. Structural nested cumulative failure time models to estimate the effects of interventions. Journal of the American Statistical Association. 2012 Sep 1;107(499):886-900.")
2. Shi J, Swanson SA, Kraft P, Rosner B, De Vivo I, Hernán MA. Instrumental variable estimation for a time-varying treatment and a time-to-event outcome via structural nested cumulative failure time models. In preparation. 2021. 

The SNCFTM functions are based on the SNCFTM SAS macro by Sally Picciotto (more information is available [here](https://www.hsph.harvard.edu/causal/software/ "here")). 

Please see the data example as an illustration of how the data must be set up to run these models.

Simulations to assess the performance of these models are also available. Data were generated according to the following DAGs:

| Data-generating model | DAG |
|------------|-------------|
| Time-fixed instrument, exposure, confounder and outcome | <img src="/dags/dag1.png" width=40%> |
| Time-fixed instrument, exposure and confounder; time-varying outcome | <img src="/dags/dag2.png" width=65%> |
| Time-fixed instrument; time-varying exposure, confounder and outcome | <img src="/dags/dag3.png" width=85%> |
