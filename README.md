
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PROjection-based Functional Invariance Testing (PROFIT)

<!-- badges: start -->
<!-- badges: end -->

The goal of PROFIT is to test whether the mean function on a
longitudinally collected functional data setting is time-invariant or
not. Mathematically speaking, suppose we observe replications of
functional data $\{Y_j(s), t_{j}\}_{j=1}^m$ observed longitudinally over
different time points in a compact time domain. Suppose
$\mu(s,t_j) = \mathbb{E}(Y_j(s))$ is the mean function of the response
at time $t_j$. We want to test whether

``` math
\mu(s,t) = \mu(s,t^\prime) \;\; \text{for all } s \text{ and } t \neq t^\prime \qquad \text{vs} \qquad \mu(s,t) \neq \mu(s,t^\prime) \;\; \text{for some } s \text{ and } t \neq t^\prime
```

In this GitHub repository, we present the R code for implementing
PROFIT. Specifically, the
[simulation](https://github.com/SalilKoner/PROFIT/tree/main/Simulation)
folder contains all the necessary code to reproduce the results
presented in the [manuscript](https://arxiv.org/abs/2104.11355).

To obtain the empirical power and size properties of PROFIT for various
simulation settings presented in the manuscript, we need to run the
[profit_sim.R](https://github.com/SalilKoner/PROFIT/blob/main/Simulation/profit_sim.R)
file. This file takes seven user input to run from a unix terminal using
Rscript command. The input required to run this file are as follows

- `factor` : factor to determine sample size and sparsity level; see
  [utilities/factor.R](https://github.com/SalilKoner/PROFIT/blob/main/Simulation/utilities/factor.R)
  for meaning of different values of this variable. It must be an
  integer between 11 and 21.
- `smooth_mean` : must be logical. TRUE for smooth mean, FALSE for
  non-smooth mean; see
  [utilities/sim_settings.R](https://github.com/SalilKoner/PROFIT/blob/main/Simulation/utilities/sim_settings.R)
  and Section S4.4 of Supplementary material for more details for the
  structure of the mean.
- `delta` : a number that controls departure from null; `delta` = 0
  indicates null, otherwise alternate hypothesis. This parameter is the
  same as the $\delta$ in Section 4.1 of the manuscript.
- `PVE` : Percentage of variation explained set to extract the
  eigenfunctions, should be at least 0.9, and less than 1.
- `N_Knots` : Number of knots, $Q$, used to represent eta_k(t) using
  truncated linear basis, should be at least 10;
- `n_rep` : Number of replication used to empirically compute the power,
  should be $\geq 5000$ for $\delta=0$ to get reliable results, for
  $\delta > 0$, it could be $1000$ or $2000$.
- `ncores` : Number of cores allocated to run the `n_rep` replications.
  Must be a positive integer $\geq 1$.

The following R script needs to be submitted in an unix terminal to
compute the the probability of rejecting the null hypothesis for PROFIT
under the setting: `factor=11`, `smooth_mean=TRUE`, `delta=0`,
`PVE=0.9`, `N_Knots=20`, `n_rep=100` (please increase it for real
experiment), and `ncores=4`

``` bash
cd ./Simulation
Rscript profit_sim.R 11 TRUE 0 0.9 20 100 4
```

The result will be saved inside the `./Results/` folder in the
appropriate sub-directory specified in the argument list.
