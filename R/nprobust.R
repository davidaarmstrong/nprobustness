globalVariables(c("std.error", "gamlss"))

#' Neumayer-Plumper Robustness
#'
#' Calculate the Neumayer and Pluemper (2017) measure
#' of robustness for a robustness model relative to
#' a baseline model
#' @param base_mod A baseline model object
#' @param robust_mod An alternative model object or a list of model objects against which robustness of the baseline model will be evaluated.
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
#' @importFrom dplyr as_tibble select mutate %>% bind_rows
#' @export
np_robust <- function(base_mod,
                      robust_mod,
                      vbl,
                      base_args = list(),
                      robust_args = list(),
                      type = c("fd", "slope", "pred"),
                      ...){
  type <- match.arg(type)
  if(!"variables" %in% names(base_args)){
    base_args$variables <- vbl
  }
  base_args$model <- base_mod
  if(!"variables" %in% names(robust_args)){
    robust_args$variables <- vbl
  }
  if(!inherits(robust_mod, "list")){
    robust_args$model <- robust_mod
  }else{
    robust_args <- lapply(robust_mod, \(m){
      args <- robust_args
      args$model <- m
      args
    })
  }
  if(type == "fd"){
    b_comps <- suppressWarnings(do.call(avg_comparisons, base_args))
    if(!inherits(robust_mod, "list")){
      b_rob <- suppressWarnings(do.call(avg_comparisons, robust_args))
    }else{
      b_rob <- lapply(robust_args, \(a){
        suppressWarnings(do.call(avg_comparisons, a))})
    }
  }
  if(type == "slope"){
    b_comps <- suppressWarnings(do.call(avg_slopes, base_args))
    if(!inherits(robust_mod, "list")){
      b_rob <- suppressWarnings(do.call(avg_slopes, robust_args))
    }else{
      b_rob <- lapply(robust_args, \(a){
        suppressWarnings(do.call(avg_slopes, a))
    })
    }    
  }
  if(type == "pred"){
    b_comps <- suppressWarnings(do.call(avg_predictions, base_args))
    if(!inherits(robust_mod, "list")){
      b_rob <- suppressWarnings(do.call(avg_predictions, robust_args))
    }else{
      b_rob <- lapply(robust_args, \(a){
        suppressWarnings(do.call(avg_predictions, a))
      })
    }    
  }
  if(!inherits(b_rob, "list")){
    rob <- pnorm(b_comps$conf.high, b_rob$estimate, b_rob$std.error) -
      pnorm(b_comps$conf.low, b_rob$estimate, b_rob$std.error)
    res <- select(b_rob, 1:std.error) %>%
      mutate(base_est = b_comps$estimate, 
             base_lwr = b_comps$conf.low,
             base_upr = b_comps$conf.high)
    res$robust<- rob
  }else{
    rob <- lapply(b_rob, \(b){
      rob_score <- pnorm(b_comps$conf.high, b$estimate, b$std.error) -
      pnorm(b_comps$conf.low, b$estimate, b$std.error)
      res <- select(b, 1:std.error) %>%
        mutate(base_est = b_comps$estimate, 
               base_lwr = b_comps$conf.low,
               base_upr = b_comps$conf.high)
      res$robust <- rob_score
      res
    })
    if(is.null(names(robust_mod))){
      names(rob) <- paste("M", seq_along(rob), sep="")  
    }else{
      names(rob) <- names(robust_mod)
    }
    res <- bind_rows(rob, .id="model")
  }
  as_tibble(res)
}


globalVariables(c("std.error", "weight"))

#' Neumayer-Plumper Individual Robustness
#'
#' Calculate the Neumayer and Pluemper (2017) measure
#' of robustness for each observation in a robustness model relative to
#' a baseline model
#' @param base_mod A baseline model object
#' @param robust_mod An alternative model object or a list of model objects against which robustness of the baseline model will be evaluated
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
  if(!"variables" %in% names(base_args)){
    base_args$variables <- vbl
  }
  base_args$model <- base_mod
  if(!"newdata" %in% names(base_args)){
    bdat <- get_data(base_mod)
    base_args$newdata <- bdat
  }
  
  if(!"variables" %in% names(robust_args)){
    robust_args$variables <- vbl
  }
  if(!inherits(robust_mod, "list")){
    robust_args$model <- robust_mod
    robust_args$newdata <- get_data(robust_mod)
  }else{
    robust_args <- lapply(robust_mod, \(m){
      args <- robust_args
      args$model <- m
      if(!"newdata" %in% names(a)){
        args$newdata <- get_data(m)
      }
      args
    })
  }
  if(type == "fd"){
    b_comps <- suppressWarnings(do.call(comparisons, base_args))
    if(!inherits(robust_mod, "list")){
      b_rob <- suppressWarnings(do.call(comparisons, robust_args))
    }else{
      b_rob <- lapply(robust_args, \(a){
        suppressWarnings(do.call(comparisons, a))
      })
    }
  }
  if(type == "slope"){
    b_comps <- suppressWarnings(do.call(slopes, base_args))
    if(!inherits(robust_mod, "list")){
      b_rob <- suppressWarnings(do.call(slopes, robust_args))
    }else{
      b_rob <- lapply(robust_args, \(a){
        suppressWarnings(do.call(slopes, a))
      })
    }
  }
  if(type == "pred"){
    b_comps <- suppressWarnings(do.call(predictions, base_args))
    if(!inherits(robust_mod, "list")){
      b_rob <- suppressWarnings(do.call(predictions, robust_args))
    }else{
      b_rob <- lapply(robust_args, \(a){
        suppressWarnings(do.call(predictions, a))
      })
    }
  }
  if(!inherits(b_rob, "list")){
    rob <- pnorm(b_comps$conf.high, b_rob$estimate, b_rob$std.error) -
      pnorm(b_comps$conf.low, b_rob$estimate, b_rob$std.error)
    res <- select(b_rob, 1:std.error) %>%
      mutate(base_est = b_comps$estimate, 
             base_lwr = b_comps$conf.low,
             base_upr = b_comps$conf.high)
    res$robust<- rob
  }else{
    rob <- lapply(b_rob, \(b){
      rob_score <- pnorm(b_comps$conf.high, b$estimate, b$std.error) -
        pnorm(b_comps$conf.low, b$estimate, b$std.error)
      res <- select(b, 1:std.error) %>%
        mutate(base_est = b_comps$estimate, 
               base_lwr = b_comps$conf.low,
               base_upr = b_comps$conf.high)
      res$robust <- rob_score
      res
    })
    if(is.null(names(robust_mod))){
      names(rob) <- paste("M", seq_along(rob), sep="")  
    }else{
      names(rob) <- names(robust_mod)
    }
    res <- bind_rows(rob, .id="model")
  }
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
##' @importFrom utils combn
##' @importFrom MASS psi.bisquare
##' @importFrom stats deviance residuals get_all_vars na.omit
##' 
##' @export
rob_gamlss <- function(formula, data, ...){
  if(!requireNamespace("gamlss"))stop("Please install and load the gamlss package to use this function.\n")
  tmp <- get_all_vars(formula, data=data)
  tmp$weight <- 1
  tmp <- na.omit(tmp)
  mod1 <- gamlss(formula, data=tmp, weights=tmp$weight, ...)
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

#' Exclude all permutations of variables
#' 
#' Update model by excluding all permutations and combinations of indicated variables. 
#' 
#' @param vec A vector of strings giving variable names to be excluded from the model. 
#' @param always_include A vector of strings giving variable names that should always be included in the model.  Defaults to `NULL`.
#' @param baseline_model A model object that will be updated with new formula.  As such, it must be compatible with `update()`. 
#' @param ... Other arguments, currently not used. 
#' 
#' @export
#' @importFrom stats update reformulate
#' @importFrom insight get_data
exclusion_mods <- function(vec, always_include = NULL, baseline_model, ...){
  combs <- lapply(1:length(vec), \(n)combn(vec, n))
  combs <- lapply(combs, \(x)c(lapply(1:ncol(x), \(i)as.vector(x[,i]))))
  combs <- do.call(c, combs)
  incl <- lapply(combs, \(x)c(x, always_include))
  forms <- sapply(incl, \(x)reformulate(x))
  forms <- c(reformulate(always_include), forms)
  all_vars <- names(insight::get_data(baseline_model))
  included <- t(sapply(incl, \(x)all_vars %in% x))
  colnames(included) <- all_vars
  included <- included[,-1]
  mods <- lapply(forms, \(f)update(baseline_model, formula=f))
  included <- rbind(rep(FALSE, ncol(included)), included)
  attr(mods, "included") <- included
  mods
}   
