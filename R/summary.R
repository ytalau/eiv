
#' @export
print.mcgmm <- function(x, ...) {
    ## set up the framework like glm/geeglm/geese
    cat("\nCall:\n")
    dput(x$call)
    ## coefficients
    cat("\n Coefficients:\n")
    print(x$coefficients)
    invisible(x)
}
