#' Truncated Power Basis for Regression Splines
#'
#' Constructs the truncated power basis (TPB) design matrix for fitting a
#' piecewise polynomial (regression spline) of a given degree with interior
#' knots.  The matrix contains global polynomial terms
#' \eqn{x, x^2, \ldots, x^p} followed by one truncated power term
#' \eqn{(x - \xi_i)_+^p = \max(x - \xi_i, 0)^p} for each interior knot
#' \eqn{\xi_i}.  Binding an intercept column and passing the result to
#' \code{lm()} or \code{glm()} fits a natural regression spline with
#' continuous derivatives up to order \eqn{p-1} at each knot.
#'
#' @param x Numeric vector of predictor values.
#' @param degree Non-negative integer giving the polynomial degree.
#'   The default (\code{3}) produces a cubic spline.
#' @param nknots Positive integer giving the number of interior knots to place
#'   at equally-spaced quantiles of \code{x} when \code{knot_loc = NULL}.
#'   Ignored if \code{knot_loc} is supplied.  Default \code{3}.
#' @param knot_loc Optional numeric vector of knot locations.  When supplied,
#'   these values are used directly and \code{nknots} is ignored.  Must
#'   contain unique values.
#'
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{degree + length(knots)} columns, where \code{knots} is either
#'   \code{knot_loc} (if supplied) or the \code{nknots} interior quantiles of
#'   \code{x}.  Columns are named \code{tpb1}, \code{tpb2}, \ldots
#'
#' @details
#' The basis is constructed in two steps:
#' \enumerate{
#'   \item \strong{Global polynomial} — columns 1 through \code{degree} hold
#'     \eqn{x^1, x^2, \ldots, x^p}.  Together with an intercept (not included
#'     here) these span all polynomials of degree \eqn{\le p}.
#'   \item \strong{Truncated power terms} — one column per knot \eqn{\xi_i}
#'     holding \eqn{(x - \xi_i)^p \cdot \mathbf{1}(x \ge \xi_i)}.  Each term
#'     adds a change in the \eqn{p}-th derivative at \eqn{\xi_i} while
#'     preserving continuity of all lower derivatives.
#' }
#' Auto-placed knots are located at equally-spaced probability quantiles
#' (\eqn{1/(k+1), 2/(k+1), \ldots, k/(k+1)} for \eqn{k} = \code{nknots}),
#' so they adapt to the marginal distribution of \code{x} rather than its
#' range.
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#'
#' # Default: cubic spline with 3 interior knots -> 6-column basis
#' B <- tpb(x)
#' dim(B)     # 100 x 6
#'
#' # Linear spline (degree = 1) with 4 knots
#' B2 <- tpb(x, degree = 1, nknots = 4)
#'
#' # Explicit knot placement
#' B3 <- tpb(x, knot_loc = c(0.25, 0.5, 0.75))
#'
#' # Use in a regression
#' y  <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' df <- as.data.frame(tpb(x))
#' fit <- lm(y ~ ., data = df)
#'
#' @importFrom stats quantile
#' @export
tpb <- function(x, degree = 3, nknots = 3, knot_loc = NULL) {

  # --- Global polynomial terms: x^1, x^2, ..., x^degree ---
  out <- sapply(seq_len(degree), function(d) x^d)
  out <- matrix(out, nrow = length(x))   # stays a matrix when length(x) == 1

  # --- Knot locations ---
  if (is.null(knot_loc)) {
    # Place knots at equally-spaced quantiles, excluding the endpoints
    q <- seq(0, 1, length.out = nknots + 2L)
    q <- q[-c(1L, length(q))]
    s <- quantile(x, q, na.rm = TRUE)
  } else {
    s <- knot_loc
  }

  # Knot uniqueness is required for the basis to be full rank
  if (length(s) != length(unique(s))) {
    stop("The knot locations are not unique.\n")
  }

  # --- Truncated power terms: (x - s_i)^degree * I(x >= s_i) ---
  for (i in seq_along(s)) {
    out <- cbind(out, (x - s[i])^degree * (x >= s[i]))
  }

  colnames(out) <- paste0("tpb", seq_len(ncol(out)))
  out
}
