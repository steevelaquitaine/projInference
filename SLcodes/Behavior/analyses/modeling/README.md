
##README
_________

### Models of motion direction estimation

#### The basic Bayesian observer (assumed statistically optimal)

**Predictions**: The basic Bayesian observer makes very flexible predictions about bias and variability when both are studied in isolation. The model can predict any bias amounts, ranging from zero when estimates match displayed directions, to full when all estimates lie at prior mean). The model can also flexibly produce a wide range of variabilities ranging from minimal when likelihoods and priors are infinitely strong to maximal when both are flat. However and critically, the model makes narrow predictions about the ratio between bias and variability that strongly constrains what is expected of the fitted data. It prescribes that variability must decrease when bias increases. This makes any plausible outcome incompatible with that prediction inconsistent with the theory. Thus that the basic Bayesian observer simultaneously fits both mean and variability well provides impressive evidence in favor of Bayesian integration.

**fitting**: fminsearch, Nelder Mead, 200 iterations/per parameter (200*9), 10*200*9 function evaluations



#### Likelihood observer (fixed uniform prior)
**fitting**: fminsearch, Nelder Mead, 1000 iterations, 10000 function evaluations


#### The basic learning Bayesian observer

We evaluated how fast subjects learned prior strengths kp with the power law function:
```matlab
kp(t) = (1-exp(-(t-1)/tau))*kp;
```
Such that prior strengths are updated at each trial by 
```matlab
exp(-(t-1)/tau
```
Tau quantifies how fast subject learn (learning rate, see figure) and can range from 0 (1 trial learning) to inf (no learning). We initialized tau at 20, assuming slow learning converging at 100 trials:

<p align="center">
  <img src="figures/BayesianLearning.png?raw=true" alt="Sublime's custom image"/>
</p>

