#' Component-Plus-Residual Data for a Model Term
#'
#' Extracts the data needed to draw a component-plus-residual (partial
#' residual) plot for a single term in a fitted model.  The function returns
#' the fitted component \eqn{\hat{f}(x)}, its confidence band, the residuals,
#' the C+R value \eqn{\hat{f}(x) + e}, and the untransformed predictor
#' variable \eqn{x}.
#'
#' When a term is specified as a transformation of a variable (e.g.\
#' \code{"log(income)"}), \code{f_hat} reflects the transformed contribution
#' to the linear predictor while \code{x} contains the original, untransformed
#' variable extracted from \code{data}.  This makes it straightforward to plot
#' the C+R values against the predictor on its natural scale.
#'
#' For \code{glm} objects the component \eqn{\hat{f}(x)} is on the link scale
#' (as returned by \code{predict(..., type = "terms")}), and the residuals are
#' working residuals (\code{residuals(obj, type = "working")}), so the C+R
#' values are also on the link scale.  For \code{lm} objects raw residuals are
#' used.
#'
#' @param obj A fitted model object.  Methods are provided for \code{lm} and
#'   \code{glm}.
#' @param term Character string naming the model term for which C+R data are
#'   required.  Must match one of the term labels returned by
#'   \code{attr(terms(obj), "term.labels")} exactly (e.g. \code{"log(income)"},
#'   not \code{"income"}).
#' @param data The data frame used to fit \code{obj}.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return A data frame with one row per observation and the following columns:
#' \describe{
#'   \item{f_hat}{The fitted component for \code{term} (i.e. the term's
#'     contribution to the linear predictor, centred as returned by
#'     \code{predict(..., type = "terms")}).}
#'   \item{lwr}{Lower bound of the pointwise confidence interval for
#'     \code{f_hat}.}
#'   \item{upr}{Upper bound of the pointwise confidence interval for
#'     \code{f_hat}.}
#'   \item{e}{Residuals: raw residuals for \code{lm}, working residuals for
#'     \code{glm}.}
#'   \item{cr}{Component-plus-residual value: \code{f_hat + e}.}
#'   \item{x}{The untransformed predictor variable extracted from \code{data}.}
#' }
#'
#' @examples
#' data(Prestige, package = "carData")
#' Prestige$income <- Prestige$income / 1000
#' Prestige <- Prestige[complete.cases(
#'   Prestige[, c("prestige", "income", "education", "women")]), ]
#'
#' ## lm method
#' pmod <- lm(prestige ~ education + women + log(income), data = Prestige)
#' cr <- get_cr_data(pmod, "log(income)", Prestige)
#' plot(cr$x, cr$cr, xlab = "income (thousands)", ylab = "Component + Residual")
#' lines(cr$x[order(cr$x)], cr$f_hat[order(cr$x)])
#'
#' @importFrom stats predict terms setNames residuals
#' @export
get_cr_data <- function(obj, term, data, ...) {
  UseMethod("get_cr_data")
}

#' @rdname get_cr_data
#' @export
get_cr_data.lm <- function(obj, term, data, ...) {
  preds <- predict(obj,
                   type = "terms",
                   interval = "confidence",
                   terms = term)
  pred_df <- do.call(data.frame, preds[c(1, 3:4)]) |>
    setNames(c("f_hat", "lwr", "upr"))
  pred_df$e <- obj$residuals
  pred_df$cr <- pred_df$f_hat + pred_df$e
  var_name <- all.vars(str2lang(term))
  pred_df$x <- data[[var_name]]
  pred_df
}

#' @rdname get_cr_data
#' @export
get_cr_data.glm <- function(obj, term, data, ...) {
  preds <- predict(obj,
                   type = "terms",
                   interval = "confidence",
                   terms = term)
  pred_df <- do.call(data.frame, preds[c(1, 3:4)]) |>
    setNames(c("f_hat", "lwr", "upr"))
  pred_df$e <- residuals(obj, type = "working")
  pred_df$cr <- pred_df$f_hat + pred_df$e
  var_name <- all.vars(str2lang(term))
  pred_df$x <- data[[var_name]]
  pred_df
}
