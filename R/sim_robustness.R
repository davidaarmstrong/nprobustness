#' Simulate the Distribution of Robustness Under the Baseline Specification
#'
#' Repeatedly simulates the response from the baseline model, refits all
#' models on the same simulated response, computes marginal effects via a
#' user-supplied function, and evaluates \code{\link{robustness}} for each
#' replicate.  The resulting distribution contextualises the observed
#' robustness statistic: values far below the simulation bulk indicate that
#' the alternative model produces effects meaningfully different from what the
#' baseline specification would predict.
#'
#' @param base_mod The baseline model object (class \code{lm} or \code{glm}).
#' @param alt_mod The alternative model object, or a \emph{named} list of
#'   alternative model objects.
#' @param orig_data The data frame used to fit \code{base_mod}.
#' @param me_fn A function \code{function(model) -> marginaleffects_object}
#'   applied to the (simulated) \strong{baseline} model at each iteration.
#'   Capture any fixed arguments — \code{variables}, \code{newdata},
#'   \code{by}, etc. — in the closure.
#' @param alt_me_fn A function (or list of functions, one per element of
#'   \code{alt_mod}) applied to the (simulated) \strong{alternative} model at
#'   each iteration.  Defaults to \code{me_fn}.  When \code{alt_mod} is a
#'   list, \code{alt_me_fn} may be either a single function (reused for every
#'   alternative) or a list of the same length as \code{alt_mod}.  Mixing a
#'   list \code{alt_me_fn} with a scalar \code{alt_mod} (or vice versa) is an
#'   error.
#' @param type Character vector of measures to compute; see
#'   \code{\link{robustness}}.
#' @param R Number of simulation replicates.  Default \code{100}.
#' @param show_progress If \code{TRUE} (the default), a progress bar is shown
#'   for each alternative model's simulation loop via \pkg{pbapply}.
#' @param ... Unused; reserved for future arguments.
#'
#' @details
#' All models — baseline and every alternative — are refit to the \emph{same}
#' simulated response at each replicate (via \code{\link{sim_models}}).
#'
#' \code{me_fn} and \code{alt_me_fn} must produce \code{marginaleffects}
#' objects with the same number of rows in the same order.  A clear error is
#' raised if they do not, either during the actual-robustness computation or
#' inside the simulation loop (the iteration number is reported to assist
#' debugging).
#'
#' Iterations where either \code{me_fn} or \code{alt_me_fn} throws an error
#' (e.g. a singular simulated model) are silently dropped.  If many
#' iterations fail, inspect the simulated models manually.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{simulated}}{A tibble with one row per replicate × effect.
#'     Contains all columns returned by \code{\link{robustness.marginaleffects}}
#'     (identifying ME columns, \code{base_est}/\code{base_lwr}/\code{base_upr},
#'     and the requested measure columns), plus an \code{iter} column.  When
#'     \code{alt_mod} is a list, a leading \code{model} column is also
#'     present.}
#'   \item{\code{actual}}{The robustness computed from the original
#'     (un-simulated) models — same structure as one iteration of
#'     \code{simulated} but without \code{iter}.}
#' }
#'
#' @seealso \code{\link{robustness}}, \code{\link{robustness.marginaleffects}},
#'   \code{\link{sim_models}}, \code{\link{sim_robust}}
#'
#' @importFrom dplyr bind_rows
#' @importFrom pbapply pblapply
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' \dontrun{
#' data(mtcars)
#' base_mod <- lm(qsec ~ wt + hp + disp, data = mtcars)
#' alt_mod  <- lm(qsec ~ wt + hp + disp + cyl + vs, data = mtcars)
#'
#' me_fn <- function(m) marginaleffects::avg_slopes(m, variables = c("wt", "hp", "disp"))
#'
#' out <- sim_robustness(base_mod, alt_mod, mtcars, me_fn = me_fn, R = 200)
#' head(out$simulated)
#' out$actual
#' }
sim_robustness <- function(base_mod,
                           alt_mod,
                           orig_data,
                           me_fn,
                           alt_me_fn     = me_fn,
                           type          = c("ovl", "js", "kl", "np"),
                           R             = 100,
                           show_progress = TRUE,
                           ...) {
  type <- match.arg(type, several.ok = TRUE)

  # ---- Normalise alt_mod and alt_me_fn to parallel lists -------------------
  multi_alt <- inherits(alt_mod, "list")
  alt_mods  <- if (multi_alt) alt_mod else list(alt_mod)
  n_alt     <- length(alt_mods)

  if (is.list(alt_me_fn)) {
    if (!multi_alt)
      stop("`alt_me_fn` is a list but `alt_mod` is a single model; ",
           "supply a single function for `alt_me_fn` or wrap `alt_mod` in a list.")
    if (length(alt_me_fn) != n_alt)
      stop("`alt_me_fn` has ", length(alt_me_fn), " element(s) but `alt_mod` has ",
           n_alt, "; lengths must match.")
    alt_fns <- alt_me_fn
  } else {
    alt_fns <- rep(list(alt_me_fn), n_alt)
  }

  alt_names <- if (!is.null(names(alt_mods))) names(alt_mods)
               else paste0("M", seq_len(n_alt))

  appfun <- if (show_progress) pbapply::pblapply else lapply

  # ---- Compute actual robustness on the real models ------------------------
  act_list <- lapply(seq_len(n_alt), function(k) {
    b_me <- me_fn(base_mod)
    a_me <- alt_fns[[k]](alt_mods[[k]])
    if (nrow(b_me) != nrow(a_me) && nrow(b_me) != 1L)
      stop("me_fn(base_mod) returned ", nrow(b_me), " row(s) but alt_me_fn ",
           "returned ", nrow(a_me), " row(s). Both functions must produce ",
           "the same set of effects, or me_fn may return 1 row (recycled).")
    robustness(b_me, a_me, type = type)
  })
  actual <- if (n_alt == 1L) {
    act_list[[1L]]
  } else {
    dplyr::bind_rows(stats::setNames(act_list, alt_names), .id = "model")
  }

  # ---- Simulate all models from the same responses -------------------------
  # sim_models() fits base + all alternatives to the same R simulated datasets.
  # When alt_mods is a list, sm$robust[[k]][[i]] is iteration i of alternative k.
  sm <- sim_models(base_mod,
                   if (n_alt == 1L) alt_mods[[1L]] else alt_mods,
                   orig_data, R)

  base_sims    <- sm$base
  alt_sims_all <- if (n_alt == 1L) list(sm$robust) else sm$robust

  # ---- Simulation loop: one progress bar per alternative -------------------
  sim_list <- lapply(seq_len(n_alt), function(k) {
    afn      <- alt_fns[[k]]
    alt_sims <- alt_sims_all[[k]]

    rows <- appfun(seq_along(base_sims), function(i) {
      b_me <- tryCatch(suppressWarnings(me_fn(base_sims[[i]])),
                       error = function(e) NULL)
      a_me <- tryCatch(suppressWarnings(afn(alt_sims[[i]])),
                       error = function(e) NULL)
      if (is.null(b_me) || is.null(a_me)) return(NULL)
      if (nrow(b_me) != nrow(a_me) && nrow(b_me) != 1L)
        stop("me_fn(base) returned ", nrow(b_me), " row(s) but alt_me_fn(alt) ",
             "returned ", nrow(a_me), " row(s) in simulation ", i,
             ". Both functions must produce the same set of effects, ",
             "or me_fn may return 1 row (recycled).")
      robustness(b_me, a_me, type = type)
    })

    dplyr::bind_rows(rows, .id = "iter")
  })

  simulated <- if (n_alt == 1L) {
    sim_list[[1L]]
  } else {
    dplyr::bind_rows(stats::setNames(sim_list, alt_names), .id = "model")
  }

  result <- list(simulated = simulated, actual = actual)
  attr(result, "type") <- type
  class(result) <- c("sim_robustness", "list")
  result
}

#' Summarise a sim_robustness Object
#'
#' For each measure requested in \code{\link{sim_robustness}}, computes the
#' proportion of simulated robustness values that are greater than or equal to
#' the observed value (an empirical p-value under the baseline specification),
#' the mean of the simulated distribution, and the observed (actual) value.
#'
#' @param object A \code{sim_robustness} object returned by
#'   \code{\link{sim_robustness}}.
#' @param .by Character vector of column names to group by when summarising
#'   (e.g. \code{c("term", "contrast")}).  These columns are also used as the
#'   join key between the simulated and actual tibbles.  When \code{NULL}
#'   (the default), no grouping is applied; this is valid only when
#'   \code{object$actual} has a single row — an error is raised otherwise.
#' @param ... Unused.
#'
#' @return A tibble with one row per group.  For each measure \code{m}
#'   (e.g. \code{ovl}, \code{np}) three columns are produced:
#'   \describe{
#'     \item{\code{m_pct}}{Proportion of simulated values \eqn{\ge} the
#'       observed value.}
#'     \item{\code{m_mean_sim}}{Mean of the simulated distribution.}
#'     \item{\code{m_obs}}{The observed value from the real data.}
#'   }
#'
#' @seealso \code{\link{sim_robustness}}
#'
#' @importFrom dplyr across all_of cur_column left_join rename_with select summarise
#' @export
#' @method summary sim_robustness
summary.sim_robustness <- function(object, .by = NULL, ...) {
  types <- attr(object, "type")

  if (is.null(.by) && nrow(object$actual) > 1L)
    stop("`.by` is NULL but `object$actual` has ", nrow(object$actual),
         " rows. Specify `.by` with the columns that identify each effect ",
         "(e.g., `.by = \"term\"`).")

  by_cols <- if (is.null(.by)) character(0) else .by

  object$simulated |>
    dplyr::select(dplyr::all_of(c("iter", .by, types))) |>
    dplyr::left_join(
      object$actual |>
        dplyr::select(dplyr::all_of(c(.by, types))) |>
        dplyr::rename_with(\(x) paste0(x, "_obs"), dplyr::all_of(types)),
      by = by_cols
    ) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(types),
        list(
          pct      = \(x) mean(get(paste0(dplyr::cur_column(), "_obs")) >= x, na.rm = TRUE),
          mean_sim = \(x) mean(x, na.rm = TRUE),
          obs      = \(x) mean(get(paste0(dplyr::cur_column(), "_obs")), na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .by = .by
    )
}
