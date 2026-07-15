## Survival analysis of adult puffin capture and resight data from Great Island, NL, 2021-2023

### Data belongs to Dr. Pierre-Paul Bitton, Memorial University of Newfoundland

<h3>Contents</h3>

Relatively standard capture-mark-recapture data set spanning 3 years of adult puffins on Great Island. With the one addition is that daily observations were recorded, providing a data set of number of times an individual is seen per season as opposed to the usual seen/not seen. This script leverages those multiple observations per season and applies the [Riecke et al. (2022)](https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.9274) method, with the number of observations treated as an over-dispersed and zero-inflated Poisson process. All `R` and `JAGS` code is adapted from example code provided in Riecke et al. (2022).

<h3>Parameters estimated</h3>

* $\phi$: apparent survival
* $\epsilon$: mean number of times seen per detection occasion
* $\gamma$: transition probabilities of moving between availability states
* $\theta$: dispersion parameter to model individual heterogenity in detection probabilities

Gregory J. Robertson

last edit: 14 July 2026
