#' Box-Tidwell Power Transformation for GLMs
#'
#' Estimates the optimal Box-Tidwell power transformation for one or more
#' predictors in a GLM via maximum likelihood.  The transformation power
#' \eqn{\lambda} is found by maximising the log-likelihood of the GLM over
#' \eqn{\lambda}, with \eqn{\lambda} reparameterised through a scaled probit
#' link to keep the optimisation unconstrained.
#'
#' @section Formula specification:
#' The formula must contain the symbol \code{lambda} wherever the
#' transformation is to be applied.  For example, to estimate the optimal
#' power transformation of \code{x}:
#'
#' \preformatted{y ~ I(x^lambda) + z}
#'
#' During optimisation, \code{lambda} is inserted into the data frame at each
#' iteration and the formula is re-evaluated with that value, so any
#' expression involving \code{lambda} that is valid inside \code{I()} or
#' other formula terms will work.
#'
#' @section Parameterisation:
#' To allow unconstrained optimisation, the transformation power is
#' reparameterised as
#' \deqn{\lambda = w \cdot \Phi(\theta) - w/2}
#' where \eqn{\Phi} is the standard-normal CDF and \eqn{w} is \code{width}.
#' This maps \eqn{\theta \in (-\infty, \infty)} to
#' \eqn{\lambda \in (-w/2,\, w/2)}.  The default \code{width = 6} gives
#' \eqn{\lambda \in (-3, 3)}, which covers all practically relevant power
#' transformations.  Results are reported on the \eqn{\lambda} scale.
#'
#' @param formula A formula containing \code{lambda} as described above.
#'   The formula is passed directly to \code{\link[stats]{glm}}, so it
#'   retains its original environment.
#' @param data A data frame containing all variables referenced in
#'   \code{formula} (except \code{lambda}, which is supplied internally).
#' @param width Positive scalar controlling the range of \eqn{\lambda}:
#'   \eqn{\lambda \in (-\text{width}/2,\, \text{width}/2)}.  Default
#'   \code{6}, giving \eqn{\lambda \in (-3, 3)}.
#' @param glm_args A named list of additional arguments passed to
#'   \code{\link[stats]{glm}} at each iteration (e.g.
#'   \code{list(family = binomial)}).  Defaults to an empty list, which
#'   uses \code{glm}'s defaults (Gaussian family with identity link).
#' @param ... Additional arguments passed to \code{\link[maxLik]{maxLik}}
#'   (e.g. \code{method}, \code{iterlim}, \code{tol}).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{lambda}}{A named numeric vector of length 3 giving the MLE
#'     and 95\% confidence interval for \eqn{\lambda} on the original
#'     transformation scale: \code{c(estimate, lower, upper)}.}
#'   \item{\code{mle_res}}{The \code{maxLik} result object, which contains
#'     the estimates and standard errors on the \eqn{\theta} scale along
#'     with convergence information.}
#' }
#'
#' @examples
#' \dontrun{
#' data(Prestige, package = "carData")
#' Prestige$income <- Prestige$income / 1000
#'
#' ## Estimate optimal power transformation of income in a Gaussian GLM
#' fit <- glm_bt(
#'   formula  = prestige ~ I(income^lambda) + education + women,
#'   data     = Prestige,
#'   glm_args = list(family = gaussian)
#' )
#' fit$lambda   # point estimate and 95% CI on the lambda scale
#' }
#'
#' @importFrom stats glm logLik confint
#' @importFrom maxLik maxLik
#' @export
glm_bt <- function(formula, data, width = 6, glm_args = list(), ...) {

  llfun <- function(theta, .formula, .data, .width, .glm_args) {
    lambda <- .width * pnorm(theta) - (.width / 2)
    .data$lambda <- lambda
    m <- do.call(glm, c(list(formula = .formula, data = .data), .glm_args))
    as.numeric(logLik(m))
  }

  fit <- maxLik::maxLik(
    logLik   = llfun,
    start    = 1,
    .formula = formula,
    .data    = data,
    .width   = width,
    .glm_args = glm_args,
    ...
  )

  out    <- c(fit$estimate, confint(fit))
  lambda <- setNames(width * pnorm(out) - (width / 2),
                     c("estimate", "lower", "upper"))
  res <- list(lambda = lambda, mle_res = fit)
  class(res) <- "glm_bt"
  res
}

#' @rdname glm_bt
#' @param x A \code{glm_bt} object.
#' @export
print.glm_bt <- function(x, ...) {
  cat("Box-Tidwell Power Transformation\n")
  cat("---------------------------------\n")
  print(x$lambda, digits = 4)
  invisible(x)
}
