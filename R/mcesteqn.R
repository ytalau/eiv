#' @importFrom stats as.formula model.matrix model.response model.extract
#' model.frame terms pchisq printCoefmat rnorm glm coef
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
#' @details This function is used in \code{eivgmm} function to get the estimates.
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


#' @title Data processing for \code{eivgmm} and \code{eivgmm_R} function
#'
#' @description This function is built to process data in a way that all
#' the obversations are grouped by the time variable.
#'
#' @name process_eivgmm
#'
#'
#' @param formula a symbolic description of the model
#' @param data a data frame that contains variables in the model
#' and corresponding a variables identifying subjects and time points
#' @param me.var names of variables with measurement errors in the data
#' @param time.var name of variable that identifies different time points in
#' the data
#' @param id.var name of variable that identifies clusters in the data
#' @export process_eivgmm


process_eivgmm <- function(formula, data, me.var,
                          time.var, id.var) {
    ## step one order the data set
    idx_id <- which(colnames(data) %in% id.var)
    idx_tm <- which(colnames(data) %in% time.var)
    ## cannot work with data.table; need to figure this out in the future
    idx_id.time <- order(data[id.var], data[time.var])
    dat <- data[idx_id.time, ]
    dat1 <- split(dat, dat$id)
    #formula <- as.formula(formula)
    lb <- ncol(model.matrix(formula, data))
    mf <- lapply(1:length(dat1), function(i) model.frame(formula, dat1[[i]]))
    X <- lapply(1:length(dat1), function(i) model.matrix(formula, dat1[[i]]))
    xnames <- dimnames(X[[1]])[[2]]
    Y <-  lapply(1:length(mf), function(i) model.response(mf[[i]]))
    ind <- which(colnames(X[[1]]) %in% me.var)
    n <- length(mf)
    m <- nrow(X[[1]])
    return(list(X = X, Y = Y, ind = ind, n = n, m = m, lb = lb,
                xnames = xnames))
}



#' @title Estimation via Generalized Method of Moments Approach with
#' Errors in Variables Corrected
#'
#' @description The \code{eivgmm} function returns the estimation by combining
#' and solving measurement errors corrected estimating equations via generalized
#' method of moments.
#'
#' @name eivgmm
#'
#' @details The input of data must be a data.frame.
#'
#' @param start The starting values for the parameters in the regression
#' @param formula a symbolic description of the model
#' @param data a data frame that contains variables in the model
#' and corresponding a variables identifying subjects and time points
#' @param me.var names of variables with measurement errors in the data
#' @param mcov list of the covariance matrices for the surrogate variables
#' @param time.var name of variable that identifies different time points in
#' the data
#' @param id.var name of variable that identifies clusters in the data
#' @param control a list of parameters that pass into the estimating process
#' @param family the family of response variable
#' @param modify_inv a logical variable specifying whether a not invertible
#' matrix should be fixed
#' @param finsam_cor a logical variable specifying whether or not the finite
#' sample bias should be corrected
#' @importFrom nleqslv nleqslv
#' @export

## need to add the control element
eivgmm <- function(formula, data, me.var, mcov = list(),
                  time.var, id.var, start = NULL, family = "binomial",
                  control = list(), modify_inv = FALSE, finsam_cor = TRUE) {
    call <- match.call()
    formula <- as.formula(formula)

    stopifnot("The family is not considered in the function" =
                  family %in% c("gaussian", "binomial"))

    dat_out <- process_eivgmm(formula, data, me.var,
                             time.var, id.var)
    control <- do.call("eivgmm.control", control)
    if(is.null(start)) {
       start <- coef(glm(formula = formula, family = family, data = data))
    }

    family <- which(c("gaussian", "binomial") %in% family)
    res <- nleqslv(start, fn = rcpp_mcesteqn,
                   lb = dat_out$lb, n = dat_out$n, m = dat_out$m,
                   X = dat_out$X, mcov = mcov,
                   family = family,
                   ind = dat_out$ind, Y = dat_out$Y,
                   modify_inv = as.numeric(modify_inv))
    coeffs <- res$x
    convergence_code <- res$termcd
    names(coeffs) <- dat_out$xnames
    inf <- rcpp_inference(beta = coeffs,
                          lb = dat_out$lb, n = dat_out$n, m = dat_out$m,
                          X = dat_out$X, mcov = mcov,
                          family = family,
                          ind = dat_out$ind, Y = dat_out$Y,
                          finsam_cor = as.numeric(finsam_cor),
                          modify_inv = as.numeric(modify_inv))
    convergence_message <-
        switch(convergence_code,
               "Convergence of function values has been achieved",
               paste("The relative distance between consecutive solutions",
                     "is smaller than specified xtol value without convergence",
                     "of function values"),
               "No better point found without convergence of function values",
               "Iteration limit maxit exceeded without convergence")
    fit <- list(call = call,
                mcov = mcov,
                inf = inf,
                formula = formula,
                X = dat_out$X,
                Y = dat_out$Y,
                coefficients = coeffs,
                convergence_code = convergence_code,
                convergence_message = convergence_message,
                me.var = me.var)

    class(fit) <- "eivgmm"
    fit

}

#' @title Estimation via Generalized Method of Moments Approach with
#' Measurement Errors Corrected (R implementation)
#'
#' @description The \code{eivgmm_R} function returns the estimation by combining
#' and solving measurement errors corrected estimating equations via generalized
#' method of moments. This function is a pure R implementation.
#'
#' @name eivgmm_R
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
#' @export eivgmm_R


eivgmm_R <- function(formula, data, me.var, mcov,
                  time.var, id.var, init.beta, family = "binomial",
                  control = list()) {
    dat_out <- process_eivgmm(formula, data, me.var,
                             time.var, id.var)
    control <- do.call("eivgmm.control", control)
    res <- nleqslv(init.beta, fn = mcesteqn,
                   lb = dat_out$lb, n = dat_out$n, m = dat_out$m,
                   X = dat_out$X, mcov = mcov,
                   ind = dat_out$ind, Y = dat_out$Y,
                   control = control)
    coeffs <- res$x
    convergence <- res$termcd
    return(list(coefficients = coeffs, convergence = convergence))

}




#' @title Estimation via Generalized Method of Moments Approach with
#' Measurement Errors Corrected
#'
#' @description Auxiliary function for \code{eivgmm} function.
#'
#' @name eivgmm.control
#'
#' @details Only used internally in \code{eivgmm}
#'
#' @param epsilon The relative steplength tolerance; the default value is 1e-8
#' @param maxit The maximum iternation number of major iternations.
#' See corresponding documentation to \code{\link{nleqslv}{nleqslv}}
#' @param trace A logical variable indicating if detailed report of the progress
#' of iternation is given
#'
#' @export eivgmm.control


eivgmm.control <- function(epsilon = 1e-8 ,
                          trace = FALSE, maxit = 150) {
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of epsilon must be > 0")
    if(!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iternations must be > 0")
    list(xtol = epsilon,
         trace = as.integer(trace), maxit = maxit)
}


#' @export
summary.eivgmm <- function(object, ...) {

    inf <- object$inf
    coef.matrix <- data.frame(Estimate = unname(object$coefficients),
                              Std.Err = sqrt(diag(inf$acov)))
    coef.matrix$wald <- (coef.matrix$Estimate / coef.matrix$Std.Err)^2
    coef.matrix$chi.squared <- 1 - pchisq(coef.matrix$wald, df = 1)
    colnames(coef.matrix) <- c("Estimate", "Std.Err", "Wald", "Pr(>z)")
    rownames(coef.matrix) <- names(object$coefficients)
    # think about the design since mcov is no longer a matrix but a list
    # mcov <- as.matrix(object$mcov)
    # dimnames(mcov)[[1]] <- object$me.var
    # dimnames(mcov)[[2]] <- object$me.var
    out <- list(formula = object$formula,
                call = object$call,
                mcov = object$mcov,
                X = object$X,
                Y = object$Y,
                convergence_code = object$convergence_code,
                convergence_message = object$convergence_message,
                D_matrix = inf$d,
                vcov = inf$acov,
                modify_vinv_inf = inf$modify_vinv,
                modify_acovinv = inf$modify_acovinv,
                V_inv = inf$vinv,
                esteqn = inf$us,
                coef.matrix = coef.matrix)
    class(out) <- "summary.eivgmm"
    out
}



