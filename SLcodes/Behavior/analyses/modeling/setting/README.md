
##README
________


##### Setting the parameters of fminsearch - Nelder mead fit


######  Objective function : loglikelihood

- We minimized -logl (i.e., - sum(logl over trial data))
- -logl was rounded to obtain integer precision which speeded up optimization convergence (red: rounded vs gray: floats; example shown for 100 iterations and 1000 function evaluations)
<p align="center">
  <img src="ObjFunPrecision.png?raw=true" alt="Objective function convergence for different precisions"/>
</p>

######  Number of iterations and function evaluations

- iterations : 200 * number of fit parameters
- function evaluations : 10 per iteration (set at 10 * 200 * number of fit parameters)


######  Tolerance functions for x and the objective value
