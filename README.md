eiv
================

The R package **eiv** provides a main function `eivgmm` to correct
measurement error in a longitudinal setting.

## Installation

You can install the development version at
[GitHub](https://github.com/ytalau/eiv) using the follwing code,

    if (! require(devtools)) install.packages("devtools")
    devtools::install_github("https://github.com/ytalau/eiv")

## Supported families

-   Gaussian (identity link)
-   Binomial (logit link)

## Examples

Linear regression:

    library(eiv)
    set.seed(8593)
    n <- 200 # number of independent clusters
    m <- 5 # cluster size
    X <- rnorm(n * m, mean = 0, sd = 1)
    W <- X + rnorm(n * m, mean = 0, sd = 0.5)
    Y <- X * 1 + rnorm(n * m, 0, 1)
    dat <- data.frame(X = X, W = W, Y = Y, time = rep(1:m, n), 
                          id = rep(1:n, each = m))
    init <- glm(Y ~ W, data = dat)
    coef(init)

    (Intercept)           W 
      -0.027616    0.805806 

    cres <- eiv::eivgmm(Y ~ W ,
                        data = dat, me.var = "W",
                        mcov = lapply(1:m, function(i) as.matrix(0.25)),
                        time.var = "time", id.var = "id",
                        family = "gaussian",
                        init.beta = init$coefficients,
                        modify_inv = 0,
                        finsam_cor = 0)
    coef(cres)

    (Intercept)           W 
      -0.029009    0.987741 
