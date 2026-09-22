# Estimation of Cpm* under PMITICS

**Authors:** Dr. Riyadh Al-Mosawi and Prof. Sanku Dey

## Overview

This repository contains the R and Rcpp code for the paper:

> *Estimation of Cpm* under Progressive Multiple-Interval Type-I Censoring:
> Bayesian and Classical Approaches with a Log-Logistic Distribution.*

The code computes maximum likelihood estimates (MLEs), maximum
product-spacing estimates (MPSEs), bootstrap confidence intervals, Lindley and
Tierney–Kadane Bayesian approximations, and MCMC estimates with HPD credible
intervals under SEL and LINEX losses.

## Files

- `UserFunLogLogistPMITICS_cpp.cpp`: main Rcpp functions.
- `Simulation.R`: Monte Carlo simulation code.
- `README.md`: GitHub documentation.

## Required packages

```r
packages <- c(
  "Rcpp", "RcppArmadillo", "RcppEigen", "RcppNumerical",
  "RcppProgress", "roptim", "pracma", "coda", "progress"
)

missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0L) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
```

Compile the C++ functions using:

```r
Rcpp::sourceCpp("UserFunLogLogistPMITICS_cpp.cpp")
```

## Main functions

| Function | Description |
|:--|:--|
| `cdf_cpp()`, `pdf_cpp()`, `sur_cpp()`, `qnt_cpp()` | Log-logistic distribution functions |
| `Cpm_star_fun_cpp()` | Computes the process capability index Cpm* |
| `simulate_PMITICS_cpp()` | Generates a PMITICS sample |
| `loglik_cpp()`, `logmps_cpp()` | Log-likelihood and log-product-spacing functions |
| `Estim_profile_cpp()` | Computes MLE or MPSE, ESEs, and ACIs |
| `Boot_fun_cpp()` | Computes bootstrap confidence intervals |
| `Lindley_GEL_cpp()`, `Lindley_LINEX_cpp()` | Lindley approximations |
| `TK_GEL_cpp()`, `TK_LINEX_cpp()` | Tierney–Kadane approximations |
| `MCMC_cpp()` | MCMC estimates and HPD credible intervals |

Use `type = "LK"` for likelihood-based inference and `type = "PS"` for
product-spacing-based inference.

## Simulation settings

```r
index_set <- 8
index_scheme <- 3

Sim_no <- 500
B <- 1000

MC_size <- 20500
MC_burn <- 500

c_val1 <- -1.5
c_val2 <- 1.5
```

## Running the simulation

Place the R and C++ files in the same folder. Remove or modify any
machine-specific `setwd()` command, and then run:

```r
source("Simulation.R")
```

The program estimates beta, gamma, and Cpm*. It reports the estimate, bias,
squared error, estimated standard error, interval limits, interval length, and
coverage indicator. A checkpoint is saved after every 50 successful
replications.
