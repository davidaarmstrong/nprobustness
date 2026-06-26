#' Robustness for marginaleffects Objects
#'
#' Computes distribution-overlap robustness measures row-wise between two
#' \code{marginaleffects} objects (the output of \code{comparisons},
#' \code{avg_comparisons}, \code{slopes}, \code{avg_slopes},
#' \code{predictions}, or \code{avg_predictions}).  Each row's
#' \code{estimate} and \code{std.error} are treated as defining a normal
#' sampling distribution; the requested measures compare the baseline and
#' alternative sampling distributions for corresponding effects.
#'
#' @param target A \code{marginaleffects} object representing the baseline
#'   model effects.
#' @param alternative A \code{marginaleffects} object representing the
#'   alternative model effects.  Must have the same number of rows as
#'   \code{target}; rows are matched positionally.
#' @param type Character vector of measures to compute.  Same choices as
#'   \code{\link{robustness}}: any combination of \code{"ovl"}, \code{"js"},
#'   \code{"kl"}, \code{"np"}.
#' @param target_args,alt_args,target_q,support,n Unused; present only so the
#'   method signature matches the generic.
#' @param ... Unused.
#'
#' @details
#' \strong{OVL / KL / JS} are computed by numerical integration of two normal
#' densities \eqn{N(\hat\beta_{\text{base}},\, \text{SE}_{\text{base}}^2)}
#' and \eqn{N(\hat\beta_{\text{alt}},\, \text{SE}_{\text{alt}}^2)}, using
#' finite bounds wide enough to cover both distributions (no expensive support
#' scan).
#'
#' \strong{NP} uses \code{target$conf.low} and \code{target$conf.high}
#' directly (as \code{np_robust()} does), falling back to \eqn{\pm 1.96
#' \times \text{SE}_{\text{base}}} only when those columns are absent.
#'
#' @return A tibble with all columns of \code{alternative} plus
#'   \code{base_est}, \code{base_lwr}, \code{base_upr} (from \code{target})
#'   and one column per requested measure, in canonical order (\code{ovl},
#'   \code{np}, \code{js}, \code{kl}).
#'
#' @seealso \code{\link{robustness}}, \code{\link{np_robust}}
#'
#' @importFrom dplyr as_tibble
#' @importFrom pbapply pblapply
#' @importFrom stats dnorm qnorm
#' @export
#' @method robustness marginaleffects
robustness.marginaleffects <- function(target, alternative,

                                       type        = c("ovl", "js", "kl", "np"),
                                       target_args = list(),
                                       alt_args    = list(),
                                       target_q    = NULL,
                                       support     = c(-Inf, Inf),
                                       n           = 512L,
                                       ...) {
  types <- match.arg(type, several.ok = TRUE)

  n_rows <- nrow(target)
  if (nrow(alternative) != n_rows)
    stop("`target` and `alternative` must have the same number of rows.")

  types_fn <- setdiff(types, "np")
  do_np    <- "np" %in% types

  measures <- do.call(rbind, pbapply::pblapply(seq_len(n_rows), function(i) {

    est_t <- target$estimate[i]
    se_t  <- target$std.error[i]
    est_a <- alternative$estimate[i]
    se_a  <- alternative$std.error[i]

    res <- numeric(0)

    if (length(types_fn) > 0L) {
      # Finite bounds covering both normals at ±8 SDs — avoids the full
      # .find_peak() / .bound_search() support scan for each row.
      pad <- 8 * max(se_t, se_a, na.rm = TRUE)
      lo  <- min(est_t, est_a, na.rm = TRUE) - pad
      hi  <- max(est_t, est_a, na.rm = TRUE) + pad
      r <- robustness(dnorm, dnorm,
                      type        = types_fn,
                      target_args = list(mean = est_t, sd = se_t),
                      alt_args    = list(mean = est_a, sd = se_a),
                      target_q    = qnorm,
                      support     = c(lo, hi))
      res <- c(res, r)
    }

    if (do_np) {
      cl <- if (!is.null(target$conf.low)  && !is.na(target$conf.low[i]))
              target$conf.low[i]
            else qnorm(0.025, est_t, se_t)
      cu <- if (!is.null(target$conf.high) && !is.na(target$conf.high[i]))
              target$conf.high[i]
            else qnorm(0.975, est_t, se_t)
      res <- c(res, np = .np_score(cl, cu, est_a, se_a))
    }

    res[intersect(.measure_order, types)]
  }))

  # measures is an n_rows × length(types_used) matrix after rbind
  if (!is.matrix(measures))
    measures <- matrix(measures, nrow = 1L,
                       dimnames = list(NULL, intersect(.measure_order, types)))

  result          <- as.data.frame(alternative)
  result$base_est <- target$estimate
  result$base_lwr <- if (!is.null(target$conf.low))  target$conf.low  else NA_real_
  result$base_upr <- if (!is.null(target$conf.high)) target$conf.high else NA_real_
  for (m in colnames(measures)) result[[m]] <- measures[, m]

  dplyr::as_tibble(result)
}

# In marginaleffects >= 0.19, comparisons() and predictions() dropped
# "marginaleffects" from their class vector, so robustness.marginaleffects
# no longer dispatches for them.  These two aliases restore coverage.

#' @export
#' @method robustness comparisons
robustness.comparisons <- robustness.marginaleffects

#' @export
#' @method robustness predictions
robustness.predictions <- robustness.marginaleffects
