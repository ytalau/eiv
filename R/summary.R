
#' @export
print.mcgmm <- function(x, ...) {
    ## set up the framework like glm/geeglm/geese
    cat("\nCall:\n")
    dput(x$call)
    ## coefficients
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\n")
    cat(x$convergence_message)
    invisible(x)
}


#' @export
print.summary.mcgmm <- function(x, ...) {
    cat("\nFormula: ")
    print(x$formula)
    cat("\nVariance matrix for measurement errors\n")
    print(x$mcov)
    cat("\n")
    printCoefmat(x$coef.matrix, digits = 3, P.values = TRUE)
    cat("\n")
    cat(x$convergence_message)
    invisible(x)
}
