globalVariables(c("std.error", "gamlss", "a", "base_mod", "robust_args", "robust", "label", "model", "term"))

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
#' 
#' @examples
#' data(mtcars)
#' # Fit a baseline model
#' base_mod <- lm(qsec ~ wt + hp + disp, data = mtcars)
#' # Fit an alternative model to evaluate robustness
#' rob_mod <- lm(qsec ~ wt + hp + disp + cyl + vs + carb, data = mtcars)
#' # calculate robustness
#' np_robust(base_mod, rob_mod, type = "fd", 
#'          vbl = list(wt = "2sd", hp="2sd", disp="2sd"))
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
#' 
#' @examples
#' data(mtcars)
#' # Fit a baseline model
#' base_mod <- lm(qsec ~ wt + hp + disp, data = mtcars)
#' # Fit an alternative model to evaluate robustness
#' rob_mod <- lm(qsec ~ wt + hp + disp + cyl + vs + carb, data = mtcars)
#' # calculate robustness
#' ind_robust(base_mod, rob_mod, type = "pred")
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
##' 
##' @examples
##' \dontrun{
##' data(mtcars)
##' baseg <- gamlss::gamlss(qsec ~ wt + hp + disp + vs + carb, data = mtcars)
##' robg <- rob_gamlss(qsec ~ wt + hp + disp + vs + carb, data = mtcars)
##' summary(robg$mod)
##' summary(robg$weights)
##' }
rob_gamlss <- function(formula, data, ...){
  if(!requireNamespace("gamlss"))stop("Please install and load the gamlss package to use this function.\n")
  tmp <- get_all_vars(formula, data=data)
  tmp$weight <- 1
  tmp <- na.omit(tmp)
#  mod1 <- gamlss::gamlss(formula, data=tmp, weights=tmp$weight)
  mod1 <- gamlss::gamlss(formula, data=tmp, weights=tmp$weight, ...)
  devDiff <- 1
  prevDev <- deviance(mod1)
  maxit <- 30
  k <- 1
  while(k < maxit && devDiff > 0){
    e <- residuals(mod1, type="simple")
    S2e <- sum(e^2)/mod1$df.residual
    se <- e/sqrt(S2e)
    w <- MASS::psi.bisquare(se)
    tmp$weight <- w
    mod1 <- gamlss::gamlss(formula, data=tmp, weights=weight, ...)
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
#' @importFrom dplyr across everything
#' @importFrom stats update reformulate confint model.matrix family formula lm runif rnorm model.response model.frame
#' @importFrom insight get_data
#' 
#' @examples 
#' data(mtcars)
#' base_mod <- lm(qsec ~ wt + hp + disp, data = mtcars)
#' # build models
#' mods <- exclusion_mods(c("cyl", "vs", "carb"),
#'                        always_include = c("wt", "hp", "disp"),
#'                        base_mod)
#' np_robust(base_mod, mods, vbl=c("wt", "hp", "disp"), type="fd")
exclusion_mods <- function(vec, always_include = NULL, baseline_model, ...){
  combs <- lapply(1:length(vec), \(n)combn(vec, n))
  combs <- lapply(combs, \(x)c(lapply(1:ncol(x), \(i)as.vector(x[,i]))))
  combs <- do.call(c, combs)
  incl <- lapply(combs, \(x)c(x, always_include))
  forms <- sapply(incl, \(x)reformulate(x))
  forms <- c(reformulate(always_include), forms)
  all_vars <- names(insight::get_data(baseline_model))
  included <- t(sapply(incl, \(x)vec %in% x))
  colnames(included) <- vec
  mods <- lapply(forms, \(f)update(baseline_model, formula=f))
  included <- rbind(rep(FALSE, ncol(included)), included)
  included <- as_tibble(included)
  included <- included %>% mutate(across(everything(), ~ifelse(.x, "Yes", "No")))
  attr(mods, "included") <- included
  mods
}   

#' Simulate model results
#' 
#' Simulates model results under the assumption that the input model 
#' is correct.  
#' 
#' @param base_model A model object of class `lm` or `glm`. 
#' @param robust_models A list of model objects of class `lm` or `glm` that will be used to simulate calculate robustness.
#' @param orig_data A data frame with the original data. 
#' @param R Integer giving the number of times the model should be simulated. 
#' @param return_data If `TRUE`, the simulated data are returned otherwise they are not (the default). 
#' @param ... Other arguments, currently not used. 
#' 
#' @importFrom stats alias complete.cases rbinom reorder rpois
#' @export
#' @examples
#' data(mtcars)
#' lmod <- lm(qsec ~ hp + wt + vs, data=mtcars)
#' rmod <- lm(qsec ~ log(hp) + wt + vs, data=mtcars)
#' mods <- sim_models(lmod, rmod, orig_data=mtcars, R=10)
sim_models <- function(base_model, robust_models, orig_data, R=100, 
                       return_data = FALSE, ...){
  if(!inherits(base_model, "lm") & !inherits(base_model, "glm")){
    stop("Model must be an lm or glm.\n")
  }
  a <- alias(base_model)
  if(!is.null(a$Complete)){
    avars <- rownames(a$Complete)
    stop(paste0("Base model contains aliased coefficients, please remove ", 
                paste(avars, collapse=", "), "and re-estimate the model.\n"))
  }
  if(inherits(robust_models, "list")){
    for(i in 1:length(robust_models)){
      a <- alias(robust_models[[i]])
      if(!is.null(a$Complete)){
        avars <- rownames(a$Complete)
        stop(paste0("Robust model ", i, " contains aliased coefficients, please remove ", 
                    paste(avars, collapse=", "), "and re-estimate the model.\n"))
      }
    }
  }else{
    a <- alias(robust_models)
    if(!is.null(a$Complete)){
      avars <- rownames(a$Complete)
      stop(paste0("Robust model contains aliased coefficients, please remove ", 
                  paste(avars, collapse=", "), "and re-estimate the model.\n"))
    }
  }
  dats <- sim_data(base_model, orig_data, R=R, ...)
  base_mods <- lapply(dats, \(d){
    cl <- base_model$call
    cl$data <- d
    eval(cl)
  })
  if(inherits(robust_models, "list")){
    rob_mods <- lapply(robust_models, \(x){
      lapply(dats, \(d){
        cl <- x$call
        cl$data <- d
        eval(cl)
      })
    })
  }else{
    rob_mods <- lapply(dats, \(d){
      cl <- robust_models$call
      cl$data <- d
      eval(cl)
    })
  }
  res <- list(base = base_mods, robust=rob_mods)
  if(return_data){
    res$data <- dats
  }
  return(res)
}

#' Simulate data under the base model specification
#' 
#' Simulates data for lm or glm (binomial or poisson families) assuming that the
#' baseline model is correctly specified.  
#' 
#' @param model The baseline model for a robustness analysis. 
#' @param orig_data The data frame used to fit the model. 
#' @param R Number of simulated datasets to make. 
#' @param ... Other arguments, currently not implemented. 
#' 
#' @export
#' @examples
#' data(mtcars)
#' mod <- lm(qsec ~ wt + disp + hp, data=mtcars)
#' sdat <- sim_data(mod, mtcars, R=1)
#' head(sdat[[1]])
sim_data <- function(model, orig_data, R=100, ...){
  UseMethod("sim_data")
}

#' @export
#' @method sim_data lm
#' @rdname sim_data
sim_data.lm <- function(model, orig_data, R=100, ...){
  R <- as.integer(R)
  d <- get_all_vars(formula(model), data=orig_data)
  dat <- orig_data[complete.cases(d), ]
  X <- model.matrix(model)
  ci_b <- confint(model)
  sd_e <- sqrt(sum(model$residuals^2)/model$df.residual)
  lapply(1:R, \(x){
    b <- apply(ci_b, 1, \(x)runif(1, x[1], x[2]))    
    d[[1]] <- X %*%b + rnorm(nrow(X), 0, sd_e)
    cbind(d, dat[, setdiff(names(orig_data), names(d))])
  })
}

#' @export
#' @method sim_data glm
#' @rdname sim_data
sim_data.glm <- function(model, orig_data, R=100, ...){
  R <- as.integer(R)
  d <- insight::get_data(model)
  dat <- orig_data[complete.cases(d), ]
  X <- model.matrix(model)
  ci_b <- confint(model)
  lapply(1:R, \(x){
    b <- apply(ci_b, 1, \(x)runif(1, x[1], x[2]))
    eta <- X %*%b
    ystar <- family(model)$linkinv(eta)
    if(family(model)$family == "binomial"){
      if(!inherits(model.response(model.frame(model)), "matrix")){
        d[[1]] <- rbinom(nrow(X), 1, ystar)
      }else{
        n <- rowSums(model.response(model.frame(model)))
        d[[1]] <- rbinom(nrow(X), n, ystar)
        d[[2]] <- n-d[[1]]
      }
    }
    if(family(model)$family == "poisson"){
      d[[1]] <- rpois(nrow(X), ystar)
    }
    cbind(d, dat[, setdiff(names(orig_data), names(d))])
  })
}

#' Simulate robustness scores under baseline model 
#' 
#' Simulates robustness scores assuming that the baseline model is correctly
#' specified. 
#' 
#' @param base_model A model, currently lm and glm (binomial and poisson families) are supported. 
#' @param robust_models A model or list of models against which robustness will be evaluated
#' @param orig_data The data frame from which the data for `base_model` and `robust_models` are obtained. 
#' @param R Number of simulations
#' @param rob_type Type of robustness analysis to do `"avg"` uses the `np_robust()` function and `"ind"` uses the `ind_robust()` function. 
#' @param vbl Variable(s) for which robustness is to be evaluated. 
#' @param type Quantity for which robustness will be evaluated - first difference (`"fd"`), slope (`"slope"`) or predictions (`"pred"`).  
#' @param base_args Arguments to be passed to the robustness analysis for the baseline model.
#' @param robust_args Arguments to be passed to the robustness analysis for the robustness models. 
#' @param show_progress If `TRUE`, `pbapply::pblapply()` will be used instead of `lapply()`. 
#' @param arrange_robust If `TRUE`, the data will be re-arranged according to average simulated robustness value and a new factor variable
#' `label` will be created that has the appropriate levels. 
#' @param ... Other arguments passed down to the robustness functions. 
#' 
#' @export
#' @importFrom dplyr bind_rows mutate
#' @importFrom pbapply pblapply
sim_robust <- function(base_model,
                       robust_models, 
                       orig_data, 
                       R = 100, 
                       rob_type = c("avg", "ind"), 
                       vbl = NULL, 
                       type = c("fd", "slope", "pred"), 
                       base_args = NULL, 
                       robust_args = NULL, 
                       show_progress=TRUE, 
                       arrange_robust = TRUE, 
                       ...){
  rob_type <- match.arg(rob_type)
  appfun <- ifelse(show_progress, pblapply, lapply)
  robfun <- switch(rob_type,
    avg = np_robust,
    ind = ind_robust
  )
  type <- match.arg(type)
    sim_mods <- sim_models(base_model, robust_models, 
                         orig_data, R)

  base_mods <- sim_mods$base
  robust_mods <- sim_mods$robust
  n_robust <- ifelse(inherits(robust_mods[[1]], "list"), length(robust_mods), 1)
  if(n_robust == 1){
    robs <- appfun(seq_along(base_mods), \(i){
      robfun(base_mods[[i]], robust_mods[[i]], vbl = vbl, base_args = base_args, robust_args=robust_args, type=type, ...)  
    })
   robs <- bind_rows(robs, .id="iter") %>% mutate(model = "M1")
  }else{
    robs <- lapply(seq_along(robust_mods), \(m){
      tmp <- appfun(seq_along(base_mods), \(i){
        robfun(base_mods[[i]], robust_mods[[m]][[i]], vbl = vbl, base_args = base_args, robust_args=robust_args, type, ...)  
      })
      tmp <- bind_rows(tmp, .id="iter")
      tmp$model <- paste0("M", m)
      tmp
    })
    robs <- bind_rows(robs)
  }
  act <- robfun(base_model, robust_models, vbl = vbl, base_args = base_args, robust_args=robust_args, type=type, ...) 
  if(!"model" %in% names(act)){
    act$model <- "M1"
  }
  if(arrange_robust){
    robs <- robs  %>% mutate(label = factor(model), label = reorder(label, robust, mean))
    act <- act %>% mutate(label = factor(model, levels=levels(robs$label)))
  }  
  res <- list(simulated = robs, actual = act )
  return(res)
}


#' Get F-test Results
#' 
#' Get the results of an F-test for dropping each model term independently. 
#' 
#' @param x A model object
#' @param vbl A string giving the variable name associated with the effect of interest
#' @param ... other arguments, not currently implemented. 
#' 
#' @importFrom stats drop1
#' 
#' @export
#' 
#' @returns A scalar giving the p-value for the variable of interest. 
#' @examples 
#' data(mtcars)
#' mod <- lm(mpg ~ poly(wt, 2) + hp + vs, data=mtcars)
#' get_anova(mod, "wt")
get_anova <- function(x, vbl, ...){
  a <- try(drop1(x, test="F"))
  if(!inherits(a, "try-error")){
    p <- a[grep(vbl, rownames(a)), 6]  
  }else{
    p <- NA
  }
  p
}

## Make Confounding and Error Variable Models
#' 
#' Makes models where the confounding variable's correlation with the 
#' variable of interest is controlled.  
#' 
#' @param model A model object of class - `lm` and `glm` should work, though `glm` will 
#' likely require `scale_dv=FALSE`.  
#' @param data A data frame from which model variables can be retrieved. 
#' @param vbl A string giving the name of the variable to be evaluated. 
#' @param scale_dv If `scale_dv=TRUE` and `scale_vbl=FALSE`, then the dependent variable
#' in the model will be scaled to have the same variance as the variable of interest. If
#' `scale_dv=TRUE` and `scale_vbl=TRUE`, then both variables will be scaled to have variance of 1. 
#' @param scale_vbl If `scale_dv=FALSE` and `scale_vbl=TRUE`, then the variable of interest
#' in the model will be scaled to have the same variance as the dependent variable. If
#' `scale_dv=TRUE` and `scale_vbl=TRUE`, then both variables will be scaled to have variance of 1. 
#' @param n_levels The number of levels used in the weights between the common variance and error. 
#' @param ... other arguments, not currently implemented. 
#' 
#' @export 
#' @importFrom stats princomp cor sd terms
#' @importFrom broom tidy
#' @importFrom dplyr filter
#' @returns A list with three elements: 
#' * params: a data frame with the variables returned by `broom::tidy()` along with `rzx` and `rzy` #' - the correlation between the confounding variable and the independent and dependent variables as #' well as `p_anova` the p-value for the F-test of the joint significance of the regressors 
#' belonging to the variable of interest. 
#' * models: a list of the models estimated. 
#' * data: A data frame containing the variables for the models. 
#' @examples
#' data(mtcars)
#' base_mod <- lm(mpg ~ disp + wt + hp + vs, data=mtcars)
#' m_hp <- make_ME_mods(base_mod, mtcars, "hp", scale_dv=TRUE, scale_vbl=TRUE)
make_ME_mods <- function(model, data, vbl, scale_dv = TRUE, scale_vbl=TRUE, n_levels=21, ...){
  D <- get_all_vars(formula(model), data)
  D <- na.omit(D)
  dv_name <- setdiff(names(attr(terms(model), "dataClasses")), attr(terms(model), "term.labels"))
  if(scale_dv & !scale_vbl){
    D[[dv_name]] <- D[[dv_name]]*sd(D[[vbl]])/sd(D[[dv_name]])
  }
  if(!scale_dv & scale_vbl){
    D[[vbl]] <- D[[vbl]]*sd(D[[dv_name]])/sd(D[[vbl]])
  }
  if(scale_dv & scale_vbl){
    D[[vbl]] <- c(scale(D[[vbl]]))
    D[[dv_name]] <- c(scale(D[[dv_name]]))
  }
  if(!scale_dv & !scale_vbl){
    stop("For PCA to control correlation among variables, either scale_dv or scale_vbl (or both) should be true.\n")
  }
  M <- princomp(D[,c(vbl, dv_name)])$scores[,1]
  E <- rnorm(length(M))
  E <- lm(E ~ ., data=D)$residuals
  E <- E * sd(M)/sd(E)
  p <- seq(0,1, length=n_levels)
  D$M <- M; D$E <- E
  mods <- lapply(p, \(w){
    d <- D
    d$C <- w*M + (1-w)*E
    args <- list(data= d, formula = . ~ .+C, object=base_mod)
    do.call(update, args)
  })
  mod_res <- lapply(mods, \(x){
    D <- insight::get_data(x)         
    out <- broom::tidy(x) %>% filter(grepl(vbl, term))%>% 
      mutate(rzy = cor(D$C, D[[dv_name]]), 
             rzx = cor(D$C, D[[vbl]]))
    out$p_anova = get_anova(x, vbl)[1]
    out}) %>% 
    bind_rows()
  return(list(params = mod_res, models = mods, data= D))
}




