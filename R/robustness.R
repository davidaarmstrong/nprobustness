#' Distribution Similarity and Robustness Measures
#'
#' Computes one or more measures comparing two distributions, dispatching on
#' the class of \code{target}: a density function (parametric) or a numeric
#' vector of posterior draws (non-parametric).  Both arguments must be the
#' same class.
#'
#' @param target The target (baseline) distribution: a density function such as
#'   \code{dnorm}, a numeric vector of posterior/bootstrap draws, or a numeric
#'   \emph{matrix} of draws (columns are independent comparisons).
#' @param alternative The alternative distribution.  For the function method
#'   must be the same density function as \code{target}.  For the numeric
#'   method may be a matrix even when \code{target} is a vector, in which case
#'   each column of \code{alternative} is compared to \code{target}.  For the
#'   matrix method, \code{alternative} may be a matrix (column-wise pairing)
#'   or a vector (each column of \code{target} vs the same vector).
#' @param type Character string or vector of strings naming the measure(s) to
#'   compute.  Any combination of \code{"ovl"} (Overlap Coefficient),
#'   \code{"js"} (Jensen-Shannon Divergence), \code{"kl"}
#'   (Kullback-Leibler Divergence), and \code{"np"} (Neumayer-Plümper
#'   robustness) is accepted.  Partial matching is supported.
#' @param target_args Named list of parameters for \code{target} when it is a
#'   density function (e.g. \code{list(mean = 1, sd = 0.3)} for \code{dnorm}).
#'   Ignored for the numeric method.
#' @param alt_args Named list of parameters for \code{alternative} when it is
#'   a density function (e.g. \code{list(mean = 1.5, sd = 0.5)} for
#'   \code{dnorm}).  To evaluate multiple alternatives against a common target,
#'   supply vectors for any parameter that varies across alternatives
#'   (e.g. \code{list(mean = c(1, 2, 3), sd = c(0.3, 0.5, 0.4))}); scalars
#'   are recycled across alternatives.  The result is then a tibble with one
#'   row per alternative.  Ignored for the numeric and matrix methods.
#' @param target_q The quantile function corresponding to \code{target} (e.g.
#'   \code{qnorm} when \code{target = dnorm}).  Used only by the
#'   \strong{function} method with \code{"np"} in \code{type} to obtain the
#'   2.5th and 97.5th percentiles of the target distribution exactly.  When
#'   \code{NULL} (the default), the function method attempts to infer the
#'   quantile function automatically by replacing a leading \code{"d"} in the
#'   name of \code{target} with \code{"q"} (e.g. \code{dnorm}
#'   \eqn{\to} \code{qnorm}, \code{dt} \eqn{\to} \code{qt}) and checking
#'   whether that function exists on the search path.  If inference succeeds
#'   the found function is used; if it fails a warning is issued and the
#'   fallback numerical CDF inversion is used instead.  Ignored for all
#'   \code{type} values other than \code{"np"} and for the numeric method.
#' @param support Numeric vector of length 2 giving the integration bounds for
#'   the function method.  Defaults to \code{c(-Inf, Inf)}.  Distributions
#'   with bounded support (e.g. \code{dbeta}) are handled automatically:
#'   their density functions return \code{NaN} outside the support, which is
#'   replaced with 0 internally, so explicit bounds are not required (though
#'   supplying them, e.g. \code{support = c(0, 1)} for Beta, improves
#'   efficiency).  Also used as the search window for numerical CDF inversion
#'   when \code{"np"} is requested and \code{target_q = NULL}.
#' @param n For the \strong{numeric} method: number of equally-spaced grid
#'   points passed to \code{stats::density()}.  For the \strong{function}
#'   method with \code{"np"} in \code{type} and \code{target_q = NULL}:
#'   resolution of the numerical CDF grid used to invert quantiles.
#'   Default \code{512L}.
#' @param ... For the \strong{numeric} method, additional arguments forwarded
#'   to \code{stats::density()} (e.g. \code{bw}, \code{kernel}).  Unused for
#'   the function method.
#'
#' @details
#' \strong{Asymmetric measures:} \code{"kl"} computes
#' \eqn{KL(\mathrm{target} \| \mathrm{alternative})}.  \code{"np"} computes
#' the proportion of \code{alternative}'s probability mass that falls within
#' the 95\,\% interval of \code{target}.  Both treat \code{target} as the
#' reference (baseline) distribution.
#'
#' \strong{Symmetric measures:} \code{"ovl"} (\eqn{\int \min(f,g)\,dx},
#' range \eqn{[0,1]}) and \code{"js"} (\eqn{\frac{1}{2}KL(f\|m) +
#' \frac{1}{2}KL(g\|m)} where \eqn{m=(f+g)/2}) are symmetric in
#' \code{target} and \code{alternative}.  The JS result is divided by
#' \eqn{\log(2)}, its theoretical maximum, so it is returned on \eqn{[0,1]}.
#' Note that \code{"js"} and \code{"kl"} are \emph{divergences} (0 = identical
#' distributions) while \code{"ovl"} and \code{"np"} are \emph{similarities}
#' (1 = identical distributions).
#'
#' \strong{Function method:} All integrals are evaluated with
#' \code{stats::integrate()}.  For \code{"np"}, \code{target_q} is used to
#' obtain the 95\,\% interval exactly.  When \code{target_q = NULL}, the
#' quantile function is inferred automatically from the name of \code{target}
#' (replacing a leading \code{"d"} with \code{"q"}) if a matching function
#' exists; otherwise a warning is issued and quantiles fall back to numerical
#' CDF inversion via a two-stage scan (\eqn{\pm 100} then \eqn{\pm 10^4}).
#' When multiple measures are requested, \code{f} and \code{g} are
#' constructed once and reused across all measures.
#'
#' \strong{Numeric method:} Kernel density estimates are built on a shared
#' grid via \code{stats::density()} once, then all requested measures are
#' evaluated using the trapezoidal rule.  For \code{"np"} the 95\,\% interval
#' is the empirical \code{stats::quantile(target, c(0.025, 0.975))}.
#'
#' @return For a single target–alternative pair (the default), a named numeric
#'   vector in canonical order: \code{"ovl"}, \code{"np"}, \code{"js"},
#'   \code{"kl"}.  For vectorised calls — multiple \code{alt_args} lists
#'   (function method) or a matrix \code{alternative} (numeric/matrix methods)
#'   — a \link[tibble]{tibble} with a leading \code{model} column and one row
#'   per alternative.
#'
#' @seealso \code{\link{np_robust}}, \code{\link{ind_robust}}
#'
#' @examples
#' ## ---- Function method (parametric) ----
#'
#' # All four measures at once; qnorm is inferred automatically from "dnorm"
#' robustness(dnorm, dnorm,
#'            type        = c("ovl", "js", "kl", "np"),
#'            target_args = list(mean = 1.0, sd = 0.3),
#'            alt_args    = list(mean = 1.5, sd = 0.5))
#'
#' # Single measure
#' robustness(dnorm, dnorm, type = "np",
#'            target_args = list(mean = 1.0, sd = 0.3),
#'            alt_args    = list(mean = 1.5, sd = 0.5))
#'
#' # Supply target_q explicitly when auto-detection cannot apply
#' robustness(dnorm, dnorm, type = "np",
#'            target_args = list(mean = 1.0, sd = 0.3),
#'            alt_args    = list(mean = 1.5, sd = 0.5),
#'            target_q    = qnorm)
#'
#' ## ---- Numeric method (posterior draws) ----
#' set.seed(42)
#' base_draws <- rnorm(2000, mean = 1.0, sd = 0.3)
#' alt_draws  <- rnorm(2000, mean = 1.5, sd = 0.5)
#'
#' robustness(base_draws, alt_draws, type = c("ovl", "js", "kl", "np"))
#'
#' @importFrom stats integrate density approxfun quantile
#' @importFrom dplyr as_tibble
#' @export
robustness <- function(target, alternative,
                       type        = c("ovl", "js", "kl", "np"),
                       target_args = list(),
                       alt_args    = list(),
                       target_q    = NULL,
                       support     = c(-Inf, Inf),
                       n           = 512L,
                       ...) {
  type <- match.arg(type, several.ok = TRUE)
  UseMethod("robustness")
}


# ---- function method --------------------------------------------------------

#' @export
#' @method robustness function
robustness.function <- function(target, alternative,
                                type        = c("ovl", "js", "kl", "np"),
                                target_args = list(),
                                alt_args    = list(),
                                target_q    = NULL,
                                support     = c(-Inf, Inf),
                                n           = 512L,
                                ...) {
  types <- match.arg(type, several.ok = TRUE)

  # Auto-detect the matching quantile function from the density function's name
  # (d* -> q*) when target_q is not supplied and "np" is requested.
  # substitute() captures the caller's original expression because UseMethod
  # forwards promises unchanged.
  if ("np" %in% types && is.null(target_q)) {
    d_name      <- deparse(substitute(target))
    d_name_bare <- sub("^.*::", "", d_name)       # strip any "pkg::" prefix
    q_name      <- sub("^d", "q", d_name_bare)
    if (q_name != d_name_bare) {                  # name started with "d"
      if (exists(q_name, mode = "function", envir = .GlobalEnv, inherits = TRUE)) {
        target_q <- get(q_name, mode = "function", envir = .GlobalEnv, inherits = TRUE)
      } else {
        warning("Quantile function '", q_name, "' not found on the search path; ",
                "falling back to numerical CDF inversion. ",
                "Supply `target_q` explicitly to silence this warning.")
      }
    }
    # Name does not start with "d" (custom density) — fall through silently.
  }

  # --- Vectorised path: any alt_args element is a vector of length > 1 ---
  # API: list(mean = c(1, 2, 3), sd = c(0.3, 0.5, 0.4)).
  # Scalars are recycled; names on the first multi-valued element become row
  # labels; falls back to "M1", "M2", ... when unnamed.
  if (is.list(alt_args) && length(alt_args) > 0L &&
      max(lengths(alt_args)) > 1L) {
    max_len   <- max(lengths(alt_args))
    first_vec <- which(lengths(alt_args) == max_len)[1L]
    nms       <- names(alt_args[[first_vec]])
    if (is.null(nms)) nms <- paste0("M", seq_len(max_len))
    res_list  <- lapply(seq_len(max_len), function(i) {
      aa <- lapply(alt_args, function(v) if (length(v) == 1L) v else v[[i]])
      robustness(target, alternative,
                 type        = types,
                 target_args = target_args,
                 alt_args    = aa,
                 target_q    = target_q,
                 support     = support,
                 n           = n, ...)
    })
    df           <- as.data.frame(do.call(rbind, res_list))
    df           <- cbind(model = nms, df)
    rownames(df) <- NULL
    return(dplyr::as_tibble(df))
  }

  # Build density closures once; reused across all requested measures.
  # Replace NaN with 0: distributions with bounded support (e.g. dbeta) return
  # NaN rather than 0 outside their domain, which breaks integrate().
  # Inf is left as-is; integrable singularities at boundaries (e.g. Beta(0.5,0.5)
  # at x = 0, 1) do not affect pmin() and are handled by the KL floor.
  f <- function(x) { v <- do.call(target,      c(list(x), target_args)); v[is.nan(v)] <- 0; v }
  g <- function(x) { v <- do.call(alternative, c(list(x), alt_args));    v[is.nan(v)] <- 0; v }

  # Determine effective finite integration bounds.
  # When support is infinite, integrate() maps the real line to (-1,1) and
  # places quadrature points that are far too spread out to resolve narrow
  # distributions (e.g. posteriors with small SEs).  Scanning for the actual
  # support of f+g and integrating over that finite interval instead is both
  # more stable and faster.  The user-supplied support acts as a hard
  # constraint if finite.
  eff_lo <- support[1L]
  eff_hi <- support[2L]
  if (is.infinite(eff_lo) || is.infinite(eff_hi)) {
    # Primary path: locate each density's mode, then walk outward from it in
    # each direction until the density drops below tol * peak.  This adapts
    # automatically to any location and scale without needing a fixed search
    # window.  Each density gets its own threshold so distributions with very
    # different peak heights (e.g. narrow vs. wide) are handled fairly.
    pf <- .find_peak(f)
    pg <- .find_peak(g)
    if (!is.null(pf) && !is.null(pg)) {
      tol  <- 1e-8
      sd_f <- 1 / (pf$peak * sqrt(2 * pi))
      sd_g <- 1 / (pg$peak * sqrt(2 * pi))
      lo_f <- .bound_search(f, pf$mode, sd_f, -1L, pf$peak * tol)
      hi_f <- .bound_search(f, pf$mode, sd_f, +1L, pf$peak * tol)
      lo_g <- .bound_search(g, pg$mode, sd_g, -1L, pg$peak * tol)
      hi_g <- .bound_search(g, pg$mode, sd_g, +1L, pg$peak * tol)
      if (is.infinite(eff_lo)) eff_lo <- min(lo_f, lo_g)
      if (is.infinite(eff_hi)) eff_hi <- max(hi_f, hi_g)
    }
    # Fallback if .find_peak() returns NULL for either distribution (true peak
    # lies beyond ±10^6 or the density is zero everywhere in those ranges).
    if (is.infinite(eff_lo) || is.infinite(eff_hi)) {
      h <- function(x) f(x) + g(x)
      for (win in list(c(-10, 10), c(-100, 100), c(-1e4, 1e4))) {
        b <- .locate_support(h, win, 1000L)
        if (!is.null(b)) {
          mass <- tryCatch(integrate(h, b[1L], b[2L])$value, error = function(e) 0)
          if (mass >= 0.5) {
            if (is.infinite(eff_lo)) eff_lo <- b[1L]
            if (is.infinite(eff_hi)) eff_hi <- b[2L]
            break
          }
        }
      }
    }
  }

  vapply(types, function(t) switch(t,

    ovl = {
      integrand <- function(x) pmin(f(x), g(x))
      integrate(integrand, eff_lo, eff_hi)$value
    },

    kl = {
      integrand <- function(x) {
        fx <- f(x); gx <- g(x)
        # Floor gx at the smallest representable positive double to prevent
        # log(f/0) = Inf when g underflows while f is still representable.
        # Normal distributions have full support, so true KL is always finite;
        # this restores numerical stability without distorting the result.
        ifelse(fx <= 0, 0, fx * log(fx / pmax(gx, .Machine$double.xmin)))
      }
      integrate(integrand, eff_lo, eff_hi)$value
    },

    js = {
      m     <- function(x) pmax((f(x) + g(x)) / 2, .Machine$double.xmin)
      kl_fm <- function(x) { fx <- f(x); mx <- m(x); ifelse(fx <= 0, 0, fx * log(fx / mx)) }
      kl_gm <- function(x) { gx <- g(x); mx <- m(x); ifelse(gx <= 0, 0, gx * log(gx / mx)) }
      (integrate(kl_fm, eff_lo, eff_hi)$value +
         integrate(kl_gm, eff_lo, eff_hi)$value) / (2 * log(2))
    },

    np = {
      # Prefer exact quantiles from the matching q-function; fall back to
      # numerical CDF inversion when target_q is not supplied.
      ci <- if (!is.null(target_q)) {
        do.call(target_q, c(list(c(0.025, 0.975)), target_args))
      } else {
        .fn_quantile(f, c(0.025, 0.975), support, as.integer(n))
      }
      integrate(g, ci[1L], ci[2L])$value
    }

  ), numeric(1L))[intersect(.measure_order, types)]
}


# ---- numeric method ---------------------------------------------------------

#' @export
#' @method robustness numeric
robustness.numeric <- function(target, alternative,
                               type        = c("ovl", "js", "kl", "np"),
                               target_args = list(),
                               alt_args    = list(),
                               target_q    = NULL,
                               support     = c(-Inf, Inf),
                               n           = 512L,
                               ...) {
  types <- match.arg(type, several.ok = TRUE)

  # --- Vectorised path: alternative is a matrix ---
  # Each column of alternative is compared independently to target.
  if (is.matrix(alternative)) {
    nms      <- if (!is.null(colnames(alternative))) colnames(alternative)
                else paste0("M", seq_len(ncol(alternative)))
    res_list <- lapply(seq_len(ncol(alternative)), function(j) {
      robustness.numeric(target, alternative[, j],
                         type = types, target_args = target_args,
                         alt_args = alt_args, target_q = target_q,
                         support = support, n = n, ...)
    })
    df           <- as.data.frame(do.call(rbind, res_list))
    df           <- cbind(model = nms, df)
    rownames(df) <- NULL
    return(dplyr::as_tibble(df))
  }

  # Build the shared KDE grid once; reused across all requested measures.
  x_range <- range(c(target, alternative))
  pad     <- diff(x_range) * 0.15
  x_from  <- x_range[1L] - pad
  x_to    <- x_range[2L] + pad

  if (!is.infinite(support[1L])) x_from <- max(x_from, support[1L])
  if (!is.infinite(support[2L])) x_to   <- min(x_to,   support[2L])

  kd1 <- density(target,      from = x_from, to = x_to, n = n, ...)
  kd2 <- density(alternative, from = x_from, to = x_to, n = n, ...)

  # kd1$x == kd2$x because from / to / n all match
  x   <- kd1$x
  f   <- kd1$y
  g   <- kd2$y
  eps <- 1e-10        # floor to avoid log(0)

  vapply(types, function(t) switch(t,

    ovl = .trapz(x, pmin(f, g)),

    kl = {
      fs <- pmax(f, eps); gs <- pmax(g, eps)
      .trapz(x, fs * log(fs / gs))
    },

    js = {
      fs <- pmax(f, eps); gs <- pmax(g, eps)
      m  <- (fs + gs) / 2
      (.trapz(x, fs * log(fs / m)) + .trapz(x, gs * log(gs / m))) / (2 * log(2))
    },

    np = {
      # Empirical quantiles of the target draws give the 95 % interval exactly
      ci    <- quantile(target, c(0.025, 0.975))
      in_ci <- x >= ci[1L] & x <= ci[2L]
      if (sum(in_ci) < 2L) return(0)
      .trapz(x[in_ci], g[in_ci])
    }

  ), numeric(1L))[intersect(.measure_order, types)]
}


# ---- matrix method ----------------------------------------------------------

#' @export
#' @method robustness matrix
robustness.matrix <- function(target, alternative,
                               type        = c("ovl", "js", "kl", "np"),
                               target_args = list(),
                               alt_args    = list(),
                               target_q    = NULL,
                               support     = c(-Inf, Inf),
                               n           = 512L,
                               ...) {
  types <- match.arg(type, several.ok = TRUE)

  if (is.matrix(alternative)) {
    # Column-wise pairing: target[,j] vs alternative[,j]
    if (ncol(target) != ncol(alternative))
      stop("`target` and `alternative` matrices must have the same number of columns.")
    nms      <- if (!is.null(colnames(alternative))) colnames(alternative)
                else paste0("M", seq_len(ncol(alternative)))
    res_list <- lapply(seq_len(ncol(target)), function(j) {
      robustness.numeric(target[, j], alternative[, j],
                         type = types, support = support, n = n, ...)
    })
  } else {
    # Each column of target vs the single alternative vector
    nms      <- if (!is.null(colnames(target))) colnames(target)
                else paste0("M", seq_len(ncol(target)))
    res_list <- lapply(seq_len(ncol(target)), function(j) {
      robustness.numeric(target[, j], alternative,
                         type = types, support = support, n = n, ...)
    })
  }

  df           <- as.data.frame(do.call(rbind, res_list))
  df           <- cbind(model = nms, df)
  rownames(df) <- NULL
  dplyr::as_tibble(df)
}


# ---- shared internal helpers ------------------------------------------------

## Canonical output order: similarities first ([0,1]), then divergences.
.measure_order <- c("ovl", "np", "js", "kl")

## Core NP score: proportion of N(est, se^2) that falls within [conf_low, conf_high].
## Used by np_robust(), ind_robust(), and robustness(..., type = "np") (function method).
.np_score <- function(conf_low, conf_high, est, se) {
  pnorm(conf_high, est, se) - pnorm(conf_low, est, se)
}

## Trapezoidal integration over an arbitrary grid.
.trapz <- function(x, y) {
  dx <- diff(x)
  sum(dx * (y[-1L] + y[-length(y)]) / 2)
}

## Fallback: numerically invert a density function's CDF to obtain quantiles.
## Only called by robustness.function() for type = "np" when target_q = NULL.
.fn_quantile <- function(f, probs, support, n) {

  # Step 1 — locate finite bounds for the numerical CDF
  if (!is.infinite(support[1L]) && !is.infinite(support[2L])) {
    lo <- support[1L]
    hi <- support[2L]
  } else {
    # Two-stage scan: try ±100 first (covers most regression-scale posteriors),
    # fall back to ±1e4 for heavier-tailed or shifted distributions.
    bounds <- .locate_support(f, c(-100,  100),  500L)
    if (is.null(bounds))
      bounds <- .locate_support(f, c(-1e4,  1e4), 500L)
    if (is.null(bounds))
      stop("Cannot locate support of target density on [-1e4, 1e4]. ",
           "Supply `target_q` or an explicit `support` argument.")
    lo <- if (!is.infinite(support[1L])) max(bounds[1L], support[1L]) else bounds[1L]
    hi <- if (!is.infinite(support[2L])) min(bounds[2L], support[2L]) else bounds[2L]
  }

  # Step 2 — build fine CDF on [lo, hi] via trapezoidal rule
  x   <- seq(lo, hi, length.out = n)
  fx  <- vapply(x, f, numeric(1L))
  dx  <- x[2L] - x[1L]
  cdf <- c(0, cumsum((fx[-length(fx)] + fx[-1L]) / 2 * dx))
  cdf <- cdf / cdf[length(cdf)]     # normalise to [0, 1]

  # Step 3 — invert by linear interpolation
  vapply(probs, function(p) approxfun(cdf, x)(p), numeric(1L))
}

## Walk from `start` in the given direction (+1 or -1) until f(x) falls to or
## below `threshold`, using exponential step-doubling to bracket the crossing,
## then bisect to pinpoint it.  `step` is the initial stride (≈ 1 estimated
## SD of the distribution); it need not be exact — the doubling loop corrects
## for under- or over-estimates.
.bound_search <- function(f, start, step, direction, threshold) {
  step  <- abs(step)
  x     <- start                        # at the mode: f(x) > threshold
  x_new <- start + direction * step
  for (i in seq_len(64L)) {             # exponential expansion to bracket
    if (f(x_new) <= threshold) break
    x    <- x_new
    step <- step * 2
    x_new <- start + direction * step
  }
  for (i in seq_len(60L)) {             # bisect to refine the crossing point
    mid <- (x + x_new) / 2
    if (f(mid) > threshold) x <- mid else x_new <- mid
    if (abs(x_new - x) < 1e-12 * (1 + abs(mid))) break
  }
  x_new
}

## Locate the approximate mode and peak density value of a density function.
## Scans at progressively coarser resolutions (±100, ±10^4, ±10^6) so that
## distributions centred anywhere on the real line are found.  Returns
## list(mode, peak) where peak is the maximum density observed, or NULL if
## f appears zero throughout all three scan ranges.
.find_peak <- function(f) {
  for (rng in list(c(-100, 100), c(-1e4, 1e4), c(-1e6, 1e6))) {
    x    <- seq(rng[1L], rng[2L], length.out = 1000L)
    y    <- vapply(x, f, numeric(1L))
    peak <- max(y, na.rm = TRUE)
    if (is.finite(peak) && peak > .Machine$double.eps)
      return(list(mode = x[which.max(y)], peak = peak))
  }
  NULL
}

## Identify the region of x where f(x) > tol * max(f).
## Returns c(lo, hi) padded by one grid step, or NULL if f is zero throughout.
## An absolute floor (abs_tol) guards against false peaks from denormal
## floating-point values when the true support lies entirely outside the
## scan window (e.g. distributions centred far from zero).
.locate_support <- function(f, range, n, tol = 1e-8, abs_tol = .Machine$double.eps) {
  x     <- seq(range[1L], range[2L], length.out = n)
  y     <- vapply(x, f, numeric(1L))
  max_y <- max(y)
  if (max_y < abs_tol) return(NULL)
  idx   <- which(y > max_y * tol)
  if (length(idx) == 0L) return(NULL)
  dx    <- x[2L] - x[1L]
  c(x[max(1L, min(idx) - 1L)] - dx,
    x[min(n,  max(idx) + 1L)] + dx)
}
