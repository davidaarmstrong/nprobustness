globalVariables(c("std.error", "gamlss"))

#' Neumayer-Plumper Robustness
#'
#' Calculate the Neumayer and Pluemper (2017) measure
#' of robustness for a robustness model relative to
#' a baseline model
#' @param base_mod A baseline model object
#' @param robust_mod An alternative robustness object
#' @param vbl A character string giving the variable
#' whose robustness is being tested.
#' @param type The quantity being tested - first difference "fd",  marginal effect "slope", or prediction "pred".
#' @param base_args Arguments to be passed to `avg_slopes()` or `avg_comparisons()` that govern the baseline model effects.
#' @param robust_args Arguments to be passed to `avg_slopes()` or `avg_comparisons()` that govern the baseline model effects; these should produce the same number of
#' effect estimates as the base_args specification does.
#' @param ... Other arguments to be passed down, not implemented.
#'
#' @returns A data frame that contains the robust model estimate and standard error along with the 95% CI from the baseline model and the robustness calculation.
#'
#' @importFrom marginaleffects avg_comparisons avg_predictions avg_slopes datagrid
#' @importFrom stats pnorm
#' @importFrom dplyr as_tibble select mutate %>%
#' @export
np_robust <- function(base_mod,
                      robust_mod,
                      vbl,
                      base_args = list(),
                      robust_args = list(),
                      type = c("fd", "slope", "pred"),
                      ...){
  type <- match.arg(type)
  base_args$model <- base_mod
  robust_args$model <- robust_mod
  base_args$variables <- vbl
  robust_args$variables <- vbl
  if(type == "fd"){
    base_args$variables <- vbl
    robust_args$variables <- vbl
    b_comps <- suppressWarnings(do.call(avg_comparisons, base_args))
    b_rob <- suppressWarnings(do.call(avg_comparisons, robust_args))
  }
  if(type == "slope"){
    base_args$variables <- vbl
    robust_args$variables <- vbl
    b_comps <- suppressWarnings(do.call(avg_slopes, base_args))
    b_rob <- suppressWarnings(do.call(avg_slopes, robust_args))
  }
  if(type == "pred"){
    base_args$variables <- vbl
    robust_args$variables <- vbl
    b_comps <- suppressWarnings(do.call(avg_predictions, base_args))
    b_rob <- suppressWarnings(do.call(avg_predictions, robust_args))
  }
  rob <- pnorm(b_comps$conf.high, b_rob$estimate, b_rob$std.error) -
    pnorm(b_comps$conf.low, b_rob$estimate, b_rob$std.error)
  res <- select(b_rob, 1:std.error) %>%
    mutate(conf.low = b_comps$conf.low,
           conf.high = b_comps$conf.high)
  res$robust<- rob
  as_tibble(res)
}


globalVariables("std.error")

#' Neumayer-Plumper Individual Robustness
#'
#' Calculate the Neumayer and Pluemper (2017) measure
#' of robustness for each observation in a robustness model relative to
#' a baseline model
#' @param base_mod A baseline model object
#' @param robust_mod An alternative robustness object
#' @param vbl A character string giving the variable
#' whose robustness is being tested.  If `type = "pred"`, no variable is needed as predictions are simply made for all observations. 
#' @param type The quantity being tested - first difference "fd",  marginal effect "slope", or prediction "pred".
#' @param base_args Arguments to be passed to `avg_slopes()` or `avg_comparisons()` that govern the baseline model effects.
#' @param robust_args Arguments to be passed to `avg_slopes()` or `avg_comparisons()` that govern the baseline model effects; these should produce the same number of
#' effect estimates as the base_args specification does.
#' @param ... Other arguments to be passed down, not implemented.
#'
#' @returns A data frame that contains the robust model estimate and standard error along with the 95% CI from the baseline model and the robustness calculation.
#'
#' @importFrom marginaleffects comparisons predictions slopes 
#' @importFrom insight get_data
#' @export
ind_robust <- function(base_mod,
                           robust_mod,
                           vbl = NULL, 
                           base_args = list(),
                           robust_args = list(),
                           type = c("fd", "slope", "pred"),
                           ...){
  type <- match.arg(type)
  base_args$model <- base_mod
  robust_args$model <- robust_mod
  bdat <- get_data(base_mod)
  rdat <- get_data(robust_mod)
  base_args$newdata <- bdat
  robust_args$newdata <- rdat
  if(type == "fd"){
    if(is.null(vbl)){
      stop("Must specify variable for comparisons.\n")
    }else{
      base_args$variables <- vbl
      robust_args$variables <- vbl
    } 
    b_comps <- suppressWarnings(do.call(comparisons, base_args))
    b_rob <- suppressWarnings(do.call(comparisons, robust_args))
  }
  if(type == "slope"){
    if(is.null(vbl)){
      stop("Must specify variable for comparisons.\n")
    }else{
      base_args$variables <- vbl
      robust_args$variables <- vbl
    } 
    b_comps <- suppressWarnings(do.call(slopes, base_args))
    b_rob <- suppressWarnings(do.call(slopes, robust_args))
  }
  if(type == "pred"){
    b_comps <- suppressWarnings(do.call(predictions, base_args))
    b_rob <- suppressWarnings(do.call(predictions, robust_args))
  }
  rob <- pnorm(b_comps$conf.high, b_rob$estimate, b_rob$std.error) -
    pnorm(b_comps$conf.low, b_rob$estimate, b_rob$std.error)
  res <- b_rob %>%
    mutate(conf.low = b_comps$conf.low,
           conf.high = b_comps$conf.high)
  res$robust<- rob
  as_tibble(res)
}  

##' Estimate GAMLSS with Robustness Weights
##' 
##' Estimates the GAMLSS model with an iterative step that
##' identifies outliers and downweights them according
##' to Huber's Bisquare function.  
##' 
##' @param formula The formula for the mean in a gamlss model
##' @param data A data frame that can be used to fit the model.  As it is being used in a gamlss model, it cannot contain any missing data. 
##' @param ... Other arguments that will be passed down to `gamlss()`. 
##' @details
##' The function creates a variable called `weight` in the dataset that holds the robust regression weights.  If there is an existing
##' variable in your data called `weight`, it will be overwritten.  
##' @returns A list with elements `mod` which holds the `gamlss` model object and `weights` which give the robust regression weights
##' 
##' @importFrom MASS psi.bisquare
##' @importFrom stats deviance residuals get_all_vars na.omit
##' 
##' @export
rob_gamlss <- function(formula, data, ...){
  if(!requireNamespace("gamlss"))stop("Please install and load the gamlss package to use this function.\n")
  tmp <- get_all_vars(formula, data=data)
  tmp$weight <- 1
  tmp <- na.omit(tmp)
  mod1 <- gamlss(formula, data=tmp, weights=weight, ...)
  devDiff <- 1
  prevDev <- deviance(mod1)
  maxit <- 30
  k <- 1
  while(k < maxit && devDiff > 0){
    e <- residuals(mod1, type="simple")
    S2e <- sum(e^2)/mod1$df.residual
    se <- e/sqrt(S2e)
    w <- psi.bisquare(se)
    tmp$weight <- w
    mod1 <- gamlss(formula, data=tmp, weights=weight, ...)
    devDiff <- abs(deviance(mod1) - prevDev)
    k <- k+1
  }
  invisible(list(mod = mod1, weights = w))
}  
  


