function f = wz_wrapped_gaussian(theta, mu, rho, sd = 1, acc = 1e-05, tol = acc) 
% based on the CircStats library for R: https://cran.r-project.org/web/packages/CircStats
{
    if (missing(rho)) {
        rho <- exp(-sd^2/2)
    }
    if (rho < 0 | rho > 1) 
        stop("rho must be between 0 and 1")
    var <- -2 * log(rho)
    term <- function(theta, mu, var, k) {
        1/sqrt(var * 2 * pi) * exp(-((theta - mu + 2 * pi * k)^2)/(2 * 
            var))
    }
    k <- 0
    Next <- term(theta, mu, var, k)
    delta <- 1
    while (delta > tol) {
        k <- k + 1
        Last <- Next
        Next <- Last + term(theta, mu, var, k) + term(theta, 
            mu, var, -k)
        delta <- abs(Next - Last)
    }
    Next
}