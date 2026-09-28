# Copyright (C) 2013-2026 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Adapted from jaspAcceptanceSampling/R/bayesianSSMforAS.R at
# 7c7183f96578f26b148cff5affce1b3e7df6c986.

binomialStateSpaceInternal <- function(jaspResults, dataset = NULL, options) {
  ready <- .binomialSSMReady(options)
  jaspResults[["placeholder"]] <- NULL

  if (ready) {
    dataset <- .binomialSSM_readData(dataset, options)
    .binomialSSM_fitModel(jaspResults, dataset, options)
  } else {
    jaspResults[["fit"]] <- NULL
  }

  if (!is.null(jaspResults[["placeholder"]]))
    return()

  .binomialSSM_createSummaryTable(jaspResults, options)
  .binomialSSM_createMcmcSummaryTable(jaspResults, options)
  .binomialSSM_plotState(jaspResults, options)
  .binomialSSM_plotPosteriorDist(jaspResults, options)
  .binomialSSM_plotBetaCoefficient(jaspResults, options)
}

.binomialSSM_setAnalysisError <- function(jaspResults, message, options) {
  errorTable <- createJaspTable(title = gettext("Binomial State Space Models"))
  errorTable$dependOn(.binomialSSM_fitDependencies(options))
  errorTable$position <- 0
  jaspResults[["placeholder"]] <- errorTable
  errorTable$setError(message)
}

.binomialSSM_covariates <- function(options) {
  variables <- unlist(options[["covariates"]], use.names = FALSE)
  variables[!is.na(variables) & nzchar(variables)]
}

.binomialSSMReady <- function(options) {
  assigned <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (!assigned(options[["successes"]]))
    return(FALSE)
  if (identical(options[["trialsMode"]], "variable"))
    return(assigned(options[["trials"]]))
  # Invalid trial modes are reported by validation, not treated as missing data.
  TRUE
}

.binomialSSM_requiredColumns <- function(options) {
  columns <- c(options[["successes"]], options[["time"]], .binomialSSM_covariates(options))
  if (identical(options[["trialsMode"]], "variable"))
    columns <- c(columns, options[["trials"]])
  unique(columns[!is.na(columns) & nzchar(columns)])
}

.binomialSSM_readData <- function(dataset, options) {
  if (is.null(dataset))
    dataset <- .readDataSetToEnd(columns.as.numeric = .binomialSSM_requiredColumns(options),
                               exclude.na.listwise = NULL)

  time <- options[["time"]]
  if (length(time) == 1L && !is.na(time) && nzchar(time) && time %in% names(dataset))
    dataset <- dataset[order(dataset[[time]], na.last = TRUE), , drop = FALSE]
  dataset
}

.binomialSSM_isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

.binomialSSM_trials <- function(dataset, options) {
  if (identical(options[["trialsMode"]], "constant"))
    rep(options[["trialsConstant"]], nrow(dataset))
  else
    dataset[[options[["trials"]]]]
}

.binomialSSM_time <- function(dataset, options) {
  if (length(options[["time"]]) == 1L && nzchar(options[["time"]]))
    dataset[[options[["time"]]]]
  else
    seq_len(nrow(dataset))
}

.binomialSSM_validateOptions <- function(options) {
  mode <- options[["trialsMode"]]
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) || !mode %in% c("constant", "variable"))
    return(gettext("Select a trials variable or a constant number of trials."))
  time <- options[["time"]]
  if (!is.null(time) && (!is.character(time) || length(time) != 1L || is.na(time)))
    return(gettext("Assign a single numeric time variable or use row order."))

  scalar <- function(value, label, lower, whole = FALSE, strict = FALSE) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
      return(gettextf("%1$s must be a finite number.", label))
    if (whole && (!.binomialSSM_isWholeNumber(value) || value > .Machine$integer.max))
      return(gettextf("%1$s must be a whole number no greater than %2$s.", label, .Machine$integer.max))
    if ((strict && value <= lower) || (!strict && value < lower))
      return(if (strict) gettextf("%1$s must be greater than %2$s.", label, lower)
             else gettextf("%1$s must be at least %2$s.", label, lower))
    NULL
  }

  if (identical(options[["trialsMode"]], "constant")) {
    error <- scalar(options[["trialsConstant"]], gettext("Number of trials"), 0, whole = TRUE)
    if (!is.null(error)) return(error)
  }

  positive <- c(priorTheta1Shape1 = gettext("Initial probability prior alpha"),
                priorTheta1Shape2 = gettext("Initial probability prior beta"),
                priorSigmaSD = gettext("State-noise prior SD"))
  if (length(.binomialSSM_covariates(options)) > 0L) {
    distribution <- options[["priorBetaDistribution"]]
    if (!is.character(distribution) || length(distribution) != 1L || is.na(distribution) ||
        !distribution %in% c("normal", "studentT", "cauchy"))
      return(gettext("Select a Normal, Student t, or Cauchy prior for the covariate coefficient."))
    positive <- c(positive, switch(distribution,
      normal = c(priorBetaNormalSD = gettext("Normal coefficient prior SD")),
      studentT = c(priorBetaScale = gettext("Student t coefficient prior scale"),
                   priorBetaDf = gettext("Student t coefficient prior degrees of freedom")),
      cauchy = c(priorBetaCauchyScale = gettext("Cauchy coefficient prior scale"))))
  }
  for (name in names(positive)) {
    error <- scalar(options[[name]], positive[[name]], 0, strict = TRUE)
    if (!is.null(error)) return(error)
  }

  integers <- c(advancedMcmcBurnin = gettext("Warmup iterations"),
                advancedMcmcSamples = gettext("Sampling iterations"),
                advancedMcmcChains = gettext("Chains"),
                advancedMcmcThin = gettext("Thinning"),
                advancedMcmcSeed = gettext("Seed"),
                advancedMcmcMaxTreeDepth = gettext("Maximum tree depth"))
  for (name in names(integers)) {
    lower <- if (name == "advancedMcmcBurnin") 0 else 1
    error <- scalar(options[[name]], integers[[name]], lower, whole = TRUE)
    if (!is.null(error)) return(error)
  }
  if (options[["advancedMcmcBurnin"]] + options[["advancedMcmcSamples"]] > .Machine$integer.max)
    return(gettext("The total number of iterations exceeds the supported integer range."))
  if (options[["advancedMcmcSamples"]] < 2 * options[["advancedMcmcThin"]])
    return(gettext("Request at least two retained sampling iterations per chain after thinning."))
  error <- scalar(options[["advancedMcmcAdaptDelta"]], gettext("Target acceptance probability"), 0, strict = TRUE)
  if (!is.null(error)) return(error)
  if (options[["advancedMcmcAdaptDelta"]] >= 1)
    return(gettext("Target acceptance probability must be less than one."))
  NULL
}

.binomialSSM_validateInputs <- function(dataset, options) {
  error <- .binomialSSM_validateOptions(options)
  if (!is.null(error)) return(error)
  if (!is.data.frame(dataset) || nrow(dataset) == 0L)
    return(gettext("The dataset is empty."))
  if (nrow(dataset) < 2L)
    return(gettext("At least two time points are required."))

  covariates <- .binomialSSM_covariates(options)
  if (anyDuplicated(covariates))
    return(gettext("Each covariate must be assigned only once."))
  missingColumns <- setdiff(.binomialSSM_requiredColumns(options), names(dataset))
  if (length(missingColumns) > 0L)
    return(gettextf("The following assigned variables are missing from the dataset: %1$s.",
                    paste(missingColumns, collapse = ", ")))

  numericInput <- function(x, label, count = FALSE) {
    if (!is.numeric(x))
      return(gettextf("%1$s must be numeric.", label))
    if (any(!is.finite(x)))
      return(gettextf("%1$s cannot contain missing or non-finite values.", label))
    if (count && any(!.binomialSSM_isWholeNumber(x)))
      return(gettextf("%1$s must contain whole numbers only.", label))
    if (count && any(x < 0))
      return(gettextf("%1$s cannot contain negative values.", label))
    if (count && any(x > .Machine$integer.max))
      return(gettextf("%1$s exceeds the supported integer range.", label))
    NULL
  }

  successes <- dataset[[options[["successes"]]]]
  trials <- .binomialSSM_trials(dataset, options)
  time <- .binomialSSM_time(dataset, options)
  for (field in list(list(successes, gettext("Successes"), TRUE),
                     list(trials, gettext("Trials"), TRUE),
                     list(time, gettext("Time"), FALSE))) {
    error <- numericInput(field[[1]], field[[2]], field[[3]])
    if (!is.null(error)) return(error)
  }
  if (any(round(successes) > round(trials)))
    return(gettext("Successes cannot exceed the number of trials at any time point."))
  if (!any(round(trials) > 0))
    return(gettext("At least one time point must contain observed trials."))

  steps <- diff(sort(time))
  if (any(!is.finite(steps)))
    return(gettext("The range of time values is too large."))
  if (any(steps <= 0))
    return(gettext("Time values must be unique."))
  if (any(abs(steps - steps[1]) > sqrt(.Machine$double.eps) * max(abs(steps))))
    return(gettext("Time points must be equally spaced. Represent unobserved time points with zero trials and zero successes."))

  for (covariate in covariates) {
    predictor <- dataset[[covariate]]
    label <- gettextf("Covariate '%1$s'", jaspBase::decodeColNames(covariate))
    error <- numericInput(predictor, label)
    if (!is.null(error)) return(error)
    predictorSD <- stats::sd(predictor)
    if (!is.finite(predictorSD) || predictorSD <= 0)
      return(gettextf("%1$s must vary across time points.", label))
  }
  if (length(covariates) > 1L) {
    predictor <- base::scale(as.matrix(dataset[, covariates, drop = FALSE]))
    tolerance <- .Machine$double.eps * max(dim(predictor))
    if (qr(predictor, tol = tolerance)$rank < length(covariates))
      return(gettext("The covariates are linearly dependent. Remove redundant covariates; their separate effects cannot be distinguished by the data."))
  }
  NULL
}

.binomialSSM_fitDependencies <- function(options) {
  dependencies <- c("successes", "trialsMode", "time", "covariates",
                    "priorTheta1Shape1", "priorTheta1Shape2", "priorSigmaSD",
                    "advancedMcmcBurnin", "advancedMcmcSamples", "advancedMcmcChains",
                    "advancedMcmcThin", "advancedMcmcSeed", "advancedMcmcAdaptDelta",
                    "advancedMcmcMaxTreeDepth")
  dependencies <- c(dependencies, if (identical(options[["trialsMode"]], "constant")) "trialsConstant" else "trials")
  if (length(.binomialSSM_covariates(options)) > 0L) {
    dependencies <- c(dependencies, "priorBetaDistribution")
    distribution <- options[["priorBetaDistribution"]]
    if (identical(distribution, "normal"))
      dependencies <- c(dependencies, "priorBetaNormalSD")
    else if (identical(distribution, "studentT"))
      dependencies <- c(dependencies, "priorBetaScale", "priorBetaDf")
    else if (identical(distribution, "cauchy"))
      dependencies <- c(dependencies, "priorBetaCauchyScale")
  }
  dependencies
}

.binomialSSM_prepareStanData <- function(dataset, options) {
  stanData <- list(
    T = nrow(dataset),
    n = as.integer(round(.binomialSSM_trials(dataset, options))),
    y = as.integer(round(dataset[[options[["successes"]]]])),
    prior_theta1_a = options[["priorTheta1Shape1"]],
    prior_theta1_b = options[["priorTheta1Shape2"]],
    prior_sigma_sd = options[["priorSigmaSD"]]
  )
  covariates <- .binomialSSM_covariates(options)
  predictorCenter <- predictorScale <- numeric(0)
  modelName <- "binomialLocalLevel"
  if (length(covariates) > 0L) {
    predictor <- base::scale(as.matrix(dataset[, covariates, drop = FALSE]))
    predictorCenter <- stats::setNames(as.numeric(attr(predictor, "scaled:center")), covariates)
    predictorScale <- stats::setNames(as.numeric(attr(predictor, "scaled:scale")), covariates)
    stanData$K <- as.integer(length(covariates))
    stanData$predictor <- unname(predictor)
    stanData$prior_beta_distribution <- match(options[["priorBetaDistribution"]],
                                             c("normal", "studentT", "cauchy"))
    scaleOption <- switch(options[["priorBetaDistribution"]],
                          normal = "priorBetaNormalSD", studentT = "priorBetaScale",
                          cauchy = "priorBetaCauchyScale")
    stanData$prior_beta_scale <- options[[scaleOption]]
    stanData$prior_beta_df <- if (identical(options[["priorBetaDistribution"]], "studentT"))
      options[["priorBetaDf"]] else 3
    modelName <- "binomialLocalLevelRegression"
  }
  list(stanData = stanData, modelName = modelName, time = .binomialSSM_time(dataset, options),
       covariates = covariates, predictorCenter = predictorCenter, predictorScale = predictorScale)
}

.binomialSSM_runModel <- function(prepared, options) {
  model <- .binomialSSM_getStanModel(prepared$modelName)
  chains <- as.integer(options[["advancedMcmcChains"]])
  availableCores <- suppressWarnings(parallel::detectCores())
  if (!is.finite(availableCores) || availableCores < 1)
    availableCores <- 1L

  stanFit <- rstan::sampling(
    object = model, data = prepared$stanData, chains = chains,
    iter = options[["advancedMcmcBurnin"]] + options[["advancedMcmcSamples"]],
    warmup = options[["advancedMcmcBurnin"]], thin = options[["advancedMcmcThin"]],
    cores = min(chains, as.integer(availableCores)), refresh = 0,
    seed = options[["advancedMcmcSeed"]],
    control = list(adapt_delta = options[["advancedMcmcAdaptDelta"]],
                   max_treedepth = as.integer(options[["advancedMcmcMaxTreeDepth"]]))
  )
  sampler <- rstan::get_sampler_params(stanFit, inc_warmup = FALSE)
  diagnostics <- list(
    divergences = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
    maxTreedepth = sum(vapply(sampler, function(x)
      sum(x[, "treedepth__"] >= options[["advancedMcmcMaxTreeDepth"]]), numeric(1)))
  )
  c(list(samples = rstan::extract(stanFit), summary = rstan::summary(stanFit)$summary,
         samplerDiagnostics = diagnostics), prepared)
}

.binomialSSM_fitModel <- function(jaspResults, dataset, options) {
  if (is.null(jaspResults[["fit"]])) {
    fit <- createJaspState()
    fit$dependOn(.binomialSSM_fitDependencies(options))
    jaspResults[["fit"]] <- fit
  }
  cached <- jaspResults[["fit"]]$object
  if (!is.null(cached)) {
    if (!is.null(cached$error))
      .binomialSSM_setAnalysisError(jaspResults, cached$error, options)
    return()
  }

  validationError <- .binomialSSM_validateInputs(dataset, options)
  if (!is.null(validationError)) {
    jaspResults[["fit"]]$object <- list(error = validationError)
    .binomialSSM_setAnalysisError(jaspResults, validationError, options)
    return()
  }

  jaspBase::startProgressbar(1, gettext("Sampling..."))
  on.exit(jaspBase::progressbarTick())
  tryCatch({
    prepared <- .binomialSSM_prepareStanData(dataset, options)
    jaspResults[["fit"]]$object <- .binomialSSM_runModel(prepared, options)
  }, error = function(e) {
    message <- gettextf("Estimation failed: %1$s", conditionMessage(e))
    jaspResults[["fit"]]$object <- list(error = message)
    .binomialSSM_setAnalysisError(jaspResults, message, options)
  })
}
