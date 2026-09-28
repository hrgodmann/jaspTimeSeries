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
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Adapted from jaspAcceptanceSampling/R/bayesianSSMforAS.R,
# source commit 7c7183f96578f26b148cff5affce1b3e7df6c986.

.binomialSSM_createSummaryTable <- function(jaspResults, options) {
  if (!is.null(jaspResults[["summaryTable"]])) return()

  table <- createJaspTable(title = gettext("Posterior Summary"))
  table$dependOn(.binomialSSM_fitDependencies(options))
  table$position <- 0
  jaspResults[["summaryTable"]] <- table

  table$addColumnInfo(name = "parameter", title = gettext("Parameter"), type = "string")
  table$addColumnInfo(name = "mean", title = gettext("Mean"), type = "number")
  table$addColumnInfo(name = "sd", title = gettext("SD"), type = "number")
  table$addColumnInfo(name = "median", title = gettext("Median"), type = "number")
  table$addColumnInfo(name = "lower", title = gettext("Lower"), type = "number")
  table$addColumnInfo(name = "upper", title = gettext("Upper"), type = "number")

  state <- jaspResults[["fit"]]
  if (is.null(state) || is.null(state$object)) return()
  fit <- state$object

  table$setData(.binomialSSM_summaryRows(fit))
  table$addFootnote(gettext("Intervals are 95% equal-tail credible intervals. The latest probability is on the 0 to 1 scale and refers to the final time point in the supplied series, not a forecast beyond it. State probabilities are estimated using all observations."))
  table$addFootnote(gettext("The state-noise SD is on the log-odds scale. Probability plots display percentages."))
  if (length(fit$covariates) == 1L)
    table$addFootnote(gettext("The covariate is centered and scaled to unit SD. Its coefficient is the log-odds change for a one-SD increase in the covariate, holding the latent level fixed."))
  else if (length(fit$covariates) > 1L)
    table$addFootnote(gettext("Each covariate is centered and scaled to unit SD. Its coefficient is the log-odds change for a one-SD increase in that covariate, holding the other covariates and the latent level fixed."))

  diagnostics <- fit$samplerDiagnostics
  if (isTRUE(diagnostics$divergences > 0)) {
    table$addFootnote(gettextf("Sampling produced %1$i divergent transitions. Posterior estimates may be unreliable.",
                             as.integer(diagnostics$divergences)))
  }
  if (isTRUE(diagnostics$maxTreedepth > 0)) {
    table$addFootnote(gettextf("Sampling reached the maximum tree depth in %1$i transitions.",
                             as.integer(diagnostics$maxTreedepth)))
  }
}

.binomialSSM_summaryRows <- function(fit) {
  samples <- fit$samples
  draws <- list(samples$theta[, ncol(samples$theta)], samples$sigma_eta)
  labels <- c(gettext("Latest probability"), gettext("State-noise SD"))
  if (length(fit$covariates) > 0L) {
    beta <- .binomialSSM_coefficientSamples(fit)
    draws <- c(draws, lapply(seq_along(fit$covariates), function(i) beta[, i]))
    labels <- c(labels, gettextf("Standardized coefficient: %1$s", .binomialSSM_covariateLabels(fit)))
  }

  values <- t(vapply(draws, function(x) {
    c(mean(x), stats::sd(x), stats::median(x),
      stats::quantile(x, probs = c(0.025, 0.975), names = FALSE))
  }, numeric(5)))

  data.frame(parameter = labels, mean = values[, 1], sd = values[, 2],
             median = values[, 3], lower = values[, 4], upper = values[, 5])
}

.binomialSSM_covariateLabels <- function(fit) {
  jaspBase::decodeColNames(fit$covariates)
}

.binomialSSM_coefficientSamples <- function(fit) {
  beta <- fit$samples$beta
  if (is.null(dim(beta)))
    beta <- matrix(beta, ncol = 1L)
  beta
}

.binomialSSM_plotState <- function(jaspResults, options) {
  if (!isTRUE(options[["statePlot"]])) {
    jaspResults[["statePlot"]] <- NULL
    return()
  }
  if (!is.null(jaspResults[["statePlot"]])) return()

  plot <- createJaspPlot(title = gettext("Success Probability Over Time"), width = 600, height = 400)
  plot$dependOn(c(.binomialSSM_fitDependencies(options), "statePlot"))
  plot$position <- 1
  jaspResults[["statePlot"]] <- plot

  state <- jaspResults[["fit"]]
  if (is.null(state) || is.null(state$object)) return()
  df <- .binomialSSM_statePlotData(state$object)

  xBreaks <- jaspGraphs::getPrettyAxisBreaks(df$time)
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(df$lower, df$upper))

  plot$plotObject <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "#c2c2c2", alpha = 0.6) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = gettext("Time"), y = gettext("Probability (%)")) +
    ggplot2::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks)) +
    ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()
}

.binomialSSM_statePlotData <- function(fit) {
  theta <- fit$samples$theta * 100
  data.frame(time = fit$time, mean = colMeans(theta),
             lower = apply(theta, 2, stats::quantile, 0.025),
             upper = apply(theta, 2, stats::quantile, 0.975))
}

.binomialSSM_plotPosteriorDist <- function(jaspResults, options) {
  if (!isTRUE(options[["posteriorDistPlot"]])) {
    jaspResults[["postDist"]] <- NULL
    return()
  }
  if (!is.null(jaspResults[["postDist"]])) return()

  plot <- createJaspPlot(title = gettext("Posterior Probability at the Latest Time Point"), width = 500, height = 400)
  plot$dependOn(c(.binomialSSM_fitDependencies(options), "posteriorDistPlot"))
  plot$position <- 2
  jaspResults[["postDist"]] <- plot

  state <- jaspResults[["fit"]]
  if (is.null(state) || is.null(state$object)) return()
  theta <- state$object$samples$theta
  df <- data.frame(x = theta[, ncol(theta)] * 100)

  xBreaks <- jaspGraphs::getPrettyAxisBreaks(df$x)
  density <- stats::density(df$x)
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(density$y)))

  plot$plotObject <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_density(fill = "#b3b3b3") +
    ggplot2::labs(x = gettext("Probability (%)"), y = gettext("Density")) +
    ggplot2::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks)) +
    ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()
}

.binomialSSM_plotBetaCoefficient <- function(jaspResults, options) {
  state <- jaspResults[["fit"]]
  if (!isTRUE(options[["plotBeta"]]) || is.null(state) ||
      is.null(state$object) || length(state$object$covariates) == 0L) {
    jaspResults[["betaPlot"]] <- NULL
    return()
  }
  fit <- state$object
  if (length(fit$covariates) == 1L) {
    .binomialSSM_createCoefficientPlot(jaspResults, "betaPlot", fit, options, 1L, 3L)
    return()
  }

  if (is.null(jaspResults[["betaPlot"]])) {
    container <- createJaspContainer(title = gettext("Covariate Coefficients"))
    container$dependOn(c("covariates", "plotBeta"))
    container$position <- 3
    jaspResults[["betaPlot"]] <- container
  }

  for (i in seq_along(fit$covariates))
    .binomialSSM_createCoefficientPlot(jaspResults[["betaPlot"]], paste0("coefficient", i),
                                     fit, options, i, i)
}

.binomialSSM_createCoefficientPlot <- function(parent, key, fit, options, index, position) {
  if (!is.null(parent[[key]])) return()

  plot <- createJaspPlot(title = gettext("Covariate Coefficient"), width = 500, height = 400)
  plot$dependOn(c(.binomialSSM_fitDependencies(options), "plotBeta"))
  plot$position <- position
  parent[[key]] <- plot

  plot$title <- gettextf("Covariate Coefficient: %1$s", .binomialSSM_covariateLabels(fit)[index])
  beta <- .binomialSSM_coefficientSamples(fit)[, index]
  df <- data.frame(x = beta)
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(c(beta, 0))
  density <- stats::density(beta)
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(density$y)))

  plot$plotObject <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_density(fill = "#b3b3b3") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
    ggplot2::labs(x = gettext("Standardized Covariate Coefficient"), y = gettext("Density")) +
    ggplot2::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks)) +
    ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()
}

.binomialSSM_createMcmcSummaryTable <- function(jaspResults, options) {
  if (!isTRUE(options[["showMcmcSummary"]])) {
    jaspResults[["mcmcSummary"]] <- NULL
    return()
  }
  if (!is.null(jaspResults[["mcmcSummary"]])) return()

  table <- createJaspTable(title = gettext("MCMC Diagnostics"))
  table$dependOn(c(.binomialSSM_fitDependencies(options), "showMcmcSummary",
                  "showTheta", "showSigma", "showBeta"))
  table$position <- 4
  jaspResults[["mcmcSummary"]] <- table

  table$addColumnInfo(name = "parameter", title = gettext("Parameter"), type = "string")
  table$addColumnInfo(name = "mean", title = gettext("Mean"), type = "number")
  table$addColumnInfo(name = "sd", title = gettext("SD"), type = "number")
  table$addColumnInfo(name = "n_eff", title = gettext("ESS"), type = "number")
  table$addColumnInfo(name = "Rhat", title = gettext("R-hat"), type = "number")

  state <- jaspResults[["fit"]]
  if (is.null(state) || is.null(state$object)) return()
  fit <- state$object
  rows <- .binomialSSM_mcmcRows(fit, options)

  if (nrow(rows) == 0L) {
    table$setError(gettext("No parameters selected for the MCMC diagnostics table."))
    return()
  }

  table$setData(rows)
}

.binomialSSM_mcmcRows <- function(fit, options) {
  sm <- fit$summary
  allNames <- rownames(sm)
  pars <- character(0)

  if (isTRUE(options[["showTheta"]]))
    pars <- c(pars, grep("^theta\\[", allNames, value = TRUE))
  if (isTRUE(options[["showSigma"]]) && "sigma_eta" %in% allNames)
    pars <- c(pars, "sigma_eta")
  labels <- pars
  if (length(fit$covariates) > 0L && isTRUE(options[["showBeta"]])) {
    betaPars <- paste0("beta[", seq_along(fit$covariates), "]")
    if (length(fit$covariates) == 1L && !betaPars %in% allNames && "beta" %in% allNames)
      betaPars <- "beta"
    present <- betaPars %in% allNames
    pars <- c(pars, betaPars[present])
    labels <- c(labels, gettextf("Standardized coefficient: %1$s", .binomialSSM_covariateLabels(fit)[present]))
  }

  selected <- sm[pars, , drop = FALSE]
  data.frame(parameter = labels, mean = selected[, "mean"], sd = selected[, "sd"],
             n_eff = selected[, "n_eff"], Rhat = selected[, "Rhat"], row.names = NULL)
}
