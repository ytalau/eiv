#' @useDynLib eiv
#' @importFrom Rcpp evalCpp
#' @importFrom stats as.formula model.matrix model.response model.extract
#' model.frame
#'
#'
#' @title Generate and Correct Generalized Estimating Equations
#'
#' @description The \code{mcesteqn} function generates the corrected generalized
#' estimating equations for each subject and each time point. Then the corrected
#' estimating equations are combined using General Method of Moments (GMM) for
#' each time point.
#'
#' @name mcesteqn
#'
#' @details This function is used in \code{mcgmm} function to get the estimates.
#'
#' @param beta The value of covariates
#' @param lb the dimension of covariates
#' @param n the total number of observations
#' @param m the number of observations for each subject
#' @param X model matrix for the covariates for each ID
#' @param ind the index of the surrogate covariates
#' @param mcov the covariance matrix for the surrogate variables
#' @param Y the response variable vector for each ID
#' @export mcesteqn

mcesteqn <- function(beta, lb, n, m, X, mcov,
                       ind, Y) {
    d <- matrix(0, nrow = lb, ncol = lb*m)
    v <- matrix(0, nrow = lb*m, ncol = lb*m)
    us <- numeric(lb*m)
    for (i in 1:n) {
        ui <- numeric(lb*m)
        di <- matrix(0, nrow = lb, ncol = lb*m)
        m <- nrow(X[[i]])
        for (j in 1:m) {
            u <- X[[i]][j, ]
            a <- u %*% beta
            tmp <- drop(mcov %*% beta[ind])
            b <- tmp %*% beta[ind]
            a <- a + b/2
            a <-  exp(-a) + 1
            ay <- Y[[i]][j]
            ui[((j - 1)*lb + 1):(j*lb)] <- c(a*ay-1L)*u
            ui[(j - 1)*lb + ind] <- ui[(j - 1)*lb + ind] +
                c((a-1L)*ay)*tmp
            u1 <- u
            u1[ind] <- u1[ind] + tmp
            d1 <- u1 %o% u1
            d2 <- matrix(0, nrow = lb, ncol = lb)
            d2[ind, ind] <- mcov
            di[, ((j - 1)*lb + 1):(j*lb)] = c(ay*(a-1))*(d2-d1)

        }
        vi <- ui %o% ui
        d <- d + di
        v <- v + vi
        us <- us + ui
    }
    us <- us/n
    vi <- us %o% us
    v <- v/n - vi
    v <- v/n
    d <- d/n
    vinv <- solve(v)
    dold <- d %*% vinv
    out <- dold %*% us
    return(out)
}


#' @title Estimation via Generalized Method of Moments Approach with
#' Measurement Errors Corrected
#'
#' @description The \code{mcgmm} function returns the estimation by combining
#' and solving measurement errors corrected estimating equations via generalized
#' method of moments.
#'
#' @name mcgmm
#'
#' @details The input of data must be a data.frame.
#'
#' @param init.beta The initial guess of the estimates
#' @param formula a symbolic description of the model
#' @param data a data frame that contains variables in the model
#' and corresponding a variables identifying subjects and time points
#' @param me.var names of variables with measurement errors in the data
#' @param mcov the covariance matrix for the surrogate variables
#' @param time.var name of variable that identifies different time points in
#' the data
#' @param id.var name of variable that identifies clusters in the data
#' @param control a list of parameters that pass into the estimating process
#' @param family the family of response variable
#' @importFrom nleqslv nleqslv
#' @export mcgmm

## need to add the control element
mcgmm <- function(formula, data, me.var, mcov,
                  time.var, id.var, init.beta, family = "binomial",
                  control = list()) {
    ## step one order the data set
    idx_id <- which(colnames(data) %in% id.var)
    idx_tm <- which(colnames(data) %in% time.var)
    ## cannot work with data.table; need to figure this out in the future
    idx_id.time <- order(data[id.var], data[time.var])
    dat <- data[idx_id.time, ]
    dat1 <- split(dat, dat$id)
    formula <- as.formula(formula)
    mf <- lapply(1:length(dat1), function(i) model.frame(formula, dat1[[i]]))
    X <- lapply(1:length(dat1), function(i) model.matrix(formula, dat1[[i]]))
    Y <-  lapply(1:length(mf), function(i) model.response(mf[[i]]))
    ind <- which(colnames(X[[1]]) %in% me.var)
    n <- length(mf)
    lb <- length(init.beta)
    m <- nrow(X[[1]])
    control <- do.call("mcgmm.control", control)
    res <- nleqslv(init.beta, fn = mcesteqn,
                   lb = lb, n = n, m = m, X = X, mcov = mcov,
                   ind = ind, Y = Y,
                   control = control)
    coeffs <- res$x
    convergence <- res$termcd
    return(list(coefficients = coeffs, convergence = convergence))

}


#' @title Estimation via Generalized Method of Moments Approach with
#' Measurement Errors Corrected
#'
#' @description Auxiliary function for \code{mcgmm} function.
#'
#' @name mcgmm.control
#'
#' @details Only used internally in \code{mcgmm}
#'
#' @param epsilon The relative steplength tolerance; the default value is 1e-8
#' @param maxit The maximum iternation number of major iternations.
#' See corresponding documentation to \code{\link{nleqslv}{nleqslv}}
#' @param trace A logical variable indicating if detailed report of the progress
#' of iternation is given
#'
#' @export mcgmm.control


mcgmm.control <- function(epsilon = 1e-8 ,
                          trace = FALSE, maxit = 150) {
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of epsilon must be > 0")
    if(!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iternations must be > 0")
    list(xtol = epsilon,
         trace = as.integer(trace), maxit = maxit)
}



