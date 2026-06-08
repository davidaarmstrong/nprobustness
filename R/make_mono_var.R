#' Impose Monotonicity on a Factor Variable via Order-Restricted Regression
#'
#' Fits an order-restricted linear model using \code{ic.infer::orlm} and
#' collapses factor levels of \code{mono_var} that receive identical
#' restricted fitted values, returning a new ordered factor appended to
#' \code{data}.
#'
#' @param obj A fitted model object (e.g. from \code{lm}).
#' @param data A data frame containing \code{mono_var}.
#' @param mono_var Character string naming the factor variable in \code{data}
#'   whose levels should be constrained to be monotone.
#' @param direction One of \code{"up"} (non-decreasing) or \code{"down"}
#'   (non-increasing). Determines the direction of the monotonicity constraint.
#'
#' @return \code{data} with an additional column \code{<mono_var>_ord}, an
#'   ordered factor whose integer codes respect the specified monotone ordering.
#'   Levels that collapse to the same restricted fitted value receive the same
#'   code.
#'
#' @details
#' The function constructs the monotonicity constraint matrix via
#' \code{ic.infer::make.mon.ui}, identifies the columns of the model matrix
#' that correspond to \code{mono_var}, and passes them to
#' \code{ic.infer::orlm}. Predictions are obtained on a grid spanning all
#' factor levels (held at their observed order) using
#' \code{marginaleffects::datagrid}. Adjacent levels with identical restricted
#' fitted values (compared after rounding to 12 decimal places) are pooled into
#' a single ordinal category via \code{base::rle}.
#'
#' @seealso \code{\link[ic.infer]{orlm}}, \code{\link[ic.infer]{make.mon.ui}},
#'   \code{\link[marginaleffects]{datagrid}}
#'
#' @examples
#' \dontrun{
#' data(mtcars)
#' mtcars$cyl <- factor(mtcars$cyl)
#' fit <- lm(mpg ~ cyl + wt, data = mtcars)
#' mtcars2 <- make_mono_var(fit, mtcars, "cyl", direction = "down")
#' table(mtcars2$cyl, mtcars2$cyl_ord)
#' }
#'
#' @importFrom ic.infer make.mon.ui orlm
#' @importFrom marginaleffects datagrid
#' @importFrom dplyr arrange
#' @importFrom rlang .data
#' @export
make_mono_var <- function(obj, data, mono_var, direction = c("up", "down")) {
  direction <- match.arg(direction)

  stopifnot(mono_var %in% names(data))
  stopifnot(is.factor(data[[mono_var]]))

  sign <- ifelse(direction == "up", 1, -1)

  # Constraint matrix
  ui <- ic.infer::make.mon.ui(data[[mono_var]]) * sign

  # Columns in model matrix associated with mono_var
  trm <- match(mono_var, attr(terms(obj), "term.labels"))
  if (is.na(trm)) stop("mono_var is not a term in the model.")

  X <- model.matrix(obj)
  idx <- which(attr(X, "assign") == trm)

  # Order-restricted model
  out <- ic.infer::orlm(obj, ui = ui, index = idx)

  # Build grid over original factor levels
  levs <- levels(data[[mono_var]])

  grid_args <- list(model = obj, factor(levs, levels = levs))
  names(grid_args)[2] <- mono_var
  d <- do.call(marginaleffects::datagrid, grid_args) |>
    dplyr::arrange(.data[[mono_var]])

  # Predictions from restricted model on the grid
  Xg <- model.matrix(formula(obj), data = d)
  yhat <- as.vector(Xg %*% coef(out))

  # Pool adjacent levels with identical restricted fitted values
  pooled <- rle(round(yhat, 12))
  ord_vals <- rep(seq_along(pooled$lengths), pooled$lengths)

  lookup <- data.frame(
    old = levs,
    new = ord_vals
  )

  new_var <- paste0(mono_var, "_ord")

  data[[new_var]] <- lookup$new[
    match(as.character(data[[mono_var]]), lookup$old)
  ]

  data[[new_var]] <- factor(data[[new_var]], ordered = TRUE)

  message("Ordinal variable ", new_var, " created in returned data frame.")

  data
}
