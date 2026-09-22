#
# Copyright (C) 2013-2022 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Main function ----

bayesianStateSpaceInternal <- function(jaspResults, dataset, options) {


  # check if results can be computed
  ready <- (options$dependent != "" && any(options[c("autoregressiveComponent","localLevelComponent","localLinearTrendComponent","semiLocalLinearTrendComponent")]==TRUE,length(options$seasonalities)>0))
  # Init options: add variables to options to be used in the remainder of the analysis

  # ensure the prediction horizon is 0 whenever there are covariates or fixedFactors
  if (length(options[["covariates"]]) > 0L || length(options[["fixedFactors"]]) > 0L)
    options[["predictionHorizon"]] <- 0L

  # read dataset
  dataset <- .bstsReadData(options,ready)
  # error checking
  .bstsErrorHandling(dataset, options)

  .bstsCreateContainerMain(jaspResults,options,ready)

  # Compute (a list of) results from which tables and plots can be created
  .bstsComputeResults(jaspResults, dataset, options,ready)

  # Compute burn amount and pass it to/create options$burn
  options <- .bstsBurnHelper(jaspResults,options)
  options <- .bstsTimeHelper(jaspResults,dataset,options)

  # Output containers, tables, and plots based on the results. These functions should not return anything!
  .bstsCreateModelSummaryTable(jaspResults,options,ready)
  .bstsCreateCoefficientTable(jaspResults,options,ready)
  .bstsCreateStatePlots(jaspResults,dataset,options,ready)
  .bstsCreatePredictionPlot(jaspResults,options,ready)
  .bstsCreateControlPlots(jaspResults,options,ready)
  .bstsCreateErrorPlots(jaspResults,options,ready)

  return()
}

# Init functions ----
.bstsInitOptions <- function(jaspResults, options) {
  # Determine if analysis can be run with user input
  # Calculate any options common to multiple parts of the analysis

  return(options)
}

.bstsReadData <- function(options,ready) {
  # Read in the dataset using the built-in functions
  if(!ready) return()
  numericVariables  <- c(options$dependent,unlist(options$covariates))
  numericVariables  <- numericVariables[numericVariables != ""]
  nominalVars       <- unlist(options$fixedFactors)
  timeVar           <- unlist(options$time)
  timeVar           <- timeVar[timeVar != ""]
  dataset <- .readDataSetToEnd(columns.as.numeric  = numericVariables,
                               columns.as.factor = nominalVars,
                               columns = timeVar)

  return(dataset)


}

.bstsErrorHandling <- function(dataset, options) {
  # Custom function to check whether we missing values in predictors


  .hasErrors(dataset = dataset,
             type = 'missingValues',
             missingValues.target = options$covariates,
             exitAnalysisIfErrors = TRUE)



}

.bstsModelDependencies <- function() {
  return(c("dependent",
           "covariates",
           "posteriorSummaryTable",
           "expectedPredictors",
           "posteriorSummaryCiLevel",
           "distFam",
           "samples",
           "modelTerms",
           "seasonalities",
           "autoregressiveComponent","lagSelectionMethod","lags","maxLags","arSdPrior","arSigmaGuess","arSigmaWeight",
           "localLevelComponent",'localLevelSdPrior','localLevelSigmaGuess','localLevelSigmaWeight',
           "localLinearTrendComponent",'lltLevelPrior','lltLevelSigmaGuess','lltLevelSigmaWeight','lltSlopePrior',
           'lltSlopeSigmaGuess','lltSlopeSigmaWeight',
           "semiLocalLinearTrendComponent",
           "dynamicRegregressionComponent"
  ))
}
.bstsStatePlotDependencies <- function(){
  return(c('aggregatedStatesPlot','aggregatedStatesPlotCiLevel',"aggregatedStatesPlotObservationsShown",'componentStatesPlot'))
}

.bstsPredictionDependencies <- function(){
  return(c("predictionHorizon"))
}


.bstsControlDependencies <- function(){
  return(c('controlChartPlot',"controlPeriod","controlSigma","probalisticControlPlot"))
}


# Results functions "----",

.bstsComputeResults <- function(jaspResults, dataset, options,ready) {
  if (!ready) return()


  if (is.null(jaspResults[["bstsMainContainer"]][["bstsModelResults"]])) {
    bstsModelResultsState <- createJaspState()

    bstsModelResults <- .bstsResultsHelper(dataset,options)
    bstsModelResultsState$object <- bstsModelResults
    jaspResults[["bstsMainContainer"]][["bstsModelResults"]] <- bstsModelResultsState
  }

  if (is.null(jaspResults[["bstsMainContainer"]][["bstsModelPredictions"]]) && options$predictionHorizon > 0) {
    bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object
    bstsModelPredictionsState <- createJaspState()
    bstsModelPredictionsState$dependOn(.bstsPredictionDependencies())
    bstsPredictionResults <- bsts::predict.bsts(object = bstsResults, horizon = options$predictionHorizon, seed = options$seed)
    bstsModelPredictionsState$object <- bstsPredictionResults
    jaspResults[["bstsMainContainer"]][["bstsModelPredictions"]] <- bstsModelPredictionsState
  }

  return()
}

.bstsResultsHelper <- function(dataset,options) {

  predictors = NULL
  if (length(options$covariates)>0|length(options$fixedFactors) >0)
    predictors <- .bstsGetPredictors(options$modelTerms)
  formula = .bstsGetFormula(dependent=dataset[,options[["dependent"]]],predictors = predictors,options)

  ss   <- list()
  #AddAr
  if(options$autoregressiveComponent){
    if(options$lagSelectionMethod == "manual")
      ss <- bsts::AddAr(ss,y = dataset[,options[["dependent"]]],lags =options$lags)

    if(options$lagSelectionMethod == "auto")
      ss <- bsts::AddAutoAr(ss,y=dataset[,options[["dependent"]]],lags=options$maxLags)
  }

  #Add local level Component

  if(options$localLevelComponent)
    ss <- bsts::AddLocalLevel(ss,y=dataset[,options[["dependent"]]])

  # Add Local Linear trend component

  if(options$localLinearTrendComponent)
    ss <- bsts::AddLocalLinearTrend(ss,y=dataset[,options[["dependent"]]])

  # Add semi-local linear trend

  if(options$semiLocalLinearTrendComponent)
    ss <- bsts::AddSemilocalLinearTrend(ss,y=dataset[,options[["dependent"]]])


  if(options$dynamicRegregressionComponent)# & options$dynamicRegregressionLags==0)
    ss <- bsts::AddDynamicRegression(ss,formula=formula,data=dataset)


  if (!is.null(options$seasonalities)) {

    for (seas in options$seasonalities) {

      normalPriorSd.prior <- if(seas$inverseGammaPriorSd=="") NULL else{Boom::SdPrior(as.numeric(seas$inverseGammaPriorSd),seas$inverseGammaPriorN)}
      normal.prior <- if(seas$normalPriorSd=="") NULL else Boom::NormalPrior(seas$normalPriorMean,as.numeric(seas$normalPriorSd))
      ss <- bsts::AddSeasonal(ss,
                              y = dataset[,options[["dependent"]]],
                              nseasons = seas$number,
                              season.duration = seas$duration,
                              sigma.prior = normalPriorSd.prior,
                              initial.state.prior = normal.prior
      )
      if(!seas$name == "")
        ss[[length(ss)]]$name <- seas$name

    }
  }


  model <- bsts::bsts(formula = formula,
                      data=dataset,
                      state.specification = ss,
                      niter = options$samples,
                      timestamps=NULL,
                      seed = options$seed,
                      expected.model.size = options$expectedPredictors,
                      model.options = bsts::BstsOptions(timeout.seconds = options$timeout )
  )

  return(model)
}





# Helper function to calculate burn-in amount and pass to options$burn-> needed for summary,plots,prediction
.bstsBurnHelper <- function(jaspResults,options) {

  bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object
  if(options$burninMethod == "auto")
    options$burn <- bsts::SuggestBurn(options$automaticBurninProportion,bstsResults)

  if(options$burninMethod == "manual")
    options$burn <- options$manualBurninAmount


  return(options)
}


.bstsTimeHelper <- function(jaspResults,dataset,options){
  bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object

  actualValues <- as.numeric(bstsResults$original.series)
  options$time0 <- options$time
  if(options$time != "")
    options$time <- as.POSIXct(dataset[,options[["time"]]], tz = "UTC")
  else
    options$time <- 1:length(actualValues)

  return(options)

}

# Helper function to create regression formula that is passed into bsts functions

.bstsGetPredictors <- function(modelTerms, modelType = "alternative") {

  if (!is.character(modelType) || !modelType %in% c("alternative", "null"))
    stop(gettext("Unknown value provided for modelType, possible values: `alternative`, `null`"))

  predictors <- NULL

  for (i in seq_along(modelTerms)) {
    components <- unlist(modelTerms[[i]]$components)

    predictor <- paste0(components, collapse = ":")

    if (modelType == "alternative") {
      predictors <- c(predictors, predictor)
    } else if (modelType == "null") {
      isNuisance <- modelTerms[[i]]$isNuisance
      if (isNuisance)
        predictors <- c(predictors, predictor)
    }
  }

  return(predictors)
}

.bstsGetFormula <- function(dependent,options, predictors = NULL, includeConstant) {


  if (is.null(predictors))
    # if we don't have any predictors, the bsts function takes a vector as input
    return(dependent)
  else
    # if bsts has regression then we need an actual formula
    dependent = options[["dependent"]]
  formula <- paste(dependent, "~", paste(predictors, collapse = "+"))

  return(as.formula(formula, env = parent.frame(1)))
}



# helper function that determines which percentile rank the 2*sigma threshold is in the credible interval for each state
quantInv <- function(distr, value){

  if(sum(distr < value) == length(distr))
    return(1) # 1 if valuelarger than all MCMC draws
  else
    ecdf(distr)(value)
}





# Output functions ----
.bstsCreateContainerMain <- function(jaspResults, options,ready) {

  if (!is.null(jaspResults[["bstsMainContainer"]])) return()

  bstsMainContainer <- createJaspContainer()
  jaspResults[["bstsMainContainer"]] <- bstsMainContainer

  jaspResults[["bstsMainContainer"]]$dependOn(.bstsModelDependencies())

  return()
}


.bstsCreateModelSummaryTable <- function(jaspResults,options,ready){
  if(!is.null(jaspResults[["bstsMainContainer"]][["bstsModelSummaryTable"]])) return()

  bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object

  bstsTable <- createJaspTable(title = gettext("Model Summary"))
  bstsTable$dependOn(c("burninMethod",'automaticBurninProportion','manualBurninAmount'))
  bstsTable$position <- 1

  bstsTable$addColumnInfo(name="resSd",   title=gettext("Residual SD"),               type= "number")
  bstsTable$addColumnInfo(name="predSd",  title=gettext("Prediction SD"),             type= "number")
  bstsTable$addColumnInfo(name="R2",      title =gettextf("R%s", "\u00B2"),           type= "number")
  bstsTable$addColumnInfo(name="relGof",  title=gettext("Harvey's goodness of fit"),  type= "number")

  if(!is.null(bstsResults)){
    message <- if (bstsResults$niter < options$samples )
      gettextf("Only %1$s draws were sampled out of the desired %2$s. Additionally, %3$s MCMC draws out of %1$s are discarded as burn in.", bstsResults$niter, options$samples, options$burn)
    else
      gettextf("%1$d MCMC draws out of %2$d are discarded as burn in.", options$burn, bstsResults$niter)

    if (options$"localLevelComponent")
      message <- paste(message, gettext("The local level component has been selected (by default). If you wish to adjust the model, you can do so under 'Model Components'."), sep = "\n")
    bstsTable$addFootnote(message)
  }

  if(!is.null(bstsResults)){
    .bstsFillModelSummaryTable(bstsTable,bstsResults,ready)
  }

  jaspResults[["bstsMainContainer"]][["bstsModelSummaryTable"]] <- bstsTable

  return()
}


.bstsFillModelSummaryTable <- function(bstsTable,bstsResults,ready) {
  if(!ready) return()

  res <- summary(bstsResults)

  bstsTable$addRows(list(
    resSd   = res$residual.sd,
    predSd   = res$prediction.sd,
    R2      = res$rsquare,
    relGof  = res$relative.gof

  ))

}

# table for regression coefficients"
.bstsCreateCoefficientTable <- function(jaspResults,options,ready){
  if(!is.null(jaspResults[["bstsMainContainer"]][["bstsCoefficientSummaryTable"]]) | !length(options$modelTerms) >0 | !ready | !options$posteriorSummaryTable) return()


  bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object
  bstsCoefficientTable <- createJaspTable(title = gettext("Posterior Summary of Coefficients"))
  bstsCoefficientTable$dependOn(c("posteriorSummaryTable","showCoefMeanInc","posteriorSummaryCiLevel"))
  bstsCoefficientTable$position <- 2

  overtitle <- gettextf("%s%% Credible Interval", format(100*options[["posteriorSummaryCiLevel"]], digits = 3))
  bstsCoefficientTable$addColumnInfo(name = "coef",       title = gettext("Coefficients"),           type = "string")
  bstsCoefficientTable$addColumnInfo(name = "priorIncP",  title = gettext("P(incl)"),                type = "number")
  bstsCoefficientTable$addColumnInfo(name = "postIncP",   title = gettext("P(incl|data)"),           type = "number")
  bstsCoefficientTable$addColumnInfo(name = "BFinc",      title = gettext("BF<sub>inclusion</sub>"), type = "number")
  bstsCoefficientTable$addColumnInfo(name = 'mean',       title = gettext("Mean"),                   type = "number", format= "dp:3")
  bstsCoefficientTable$addColumnInfo(name = "sd",         title = gettext("SD"),                     type = "number",format= "dp:3")
  if (options$showCoefMeanInc){
    bstsCoefficientTable$addColumnInfo(name = 'meanInc',  title = gettext("Mean<sub>inclusion</sub>"), type = "number")
    bstsCoefficientTable$addColumnInfo(name = "sdInc",    title = gettext("SD<sub>inclusion</sub>"), type = "number")
  }
  bstsCoefficientTable$addColumnInfo(name = "lowerCri",   title = gettext("Lower"),         type = "number", overtitle = overtitle)
  bstsCoefficientTable$addColumnInfo(name = "upperCri",   title = gettext("Upper"),         type = "number", overtitle = overtitle)

  .bstsFillCoefficientTable(bstsResults,bstsCoefficientTable,options,ready)

  jaspResults[["bstsMainContainer"]][["bstsCoefficientSummaryTable"]] <- bstsCoefficientTable

}

.bstsFillCoefficientTable <- function(bstsResults,bstsCoefficientTable,options,ready) {
  res <- as.data.frame(summary(bstsResults,order = F)$coefficients)
  res$priorInc <- bstsResults$prior$prior.inclusion.probabilities
  res$BFinc <- res$inc.prob/(1-res$inc.prob)



  condQuantile <- function(beta,ci){
    if (length(beta)>0)
      return(quantile(beta,ci))
    return(0)
  }
  ci <- options$posteriorSummaryCiLevel
  res$lo_ci <- apply(bstsResults$coefficients, 2,condQuantile,((1- ci)/2))
  res$hi_ci <- apply(bstsResults$coefficients, 2,condQuantile,1-((1- ci)/2))
  res <- res[order(res$inc.prob,decreasing = TRUE),]


  for (i in 1:nrow(res)) {
    row <- list(
      coef = rownames(res[i,]),
      priorIncP = res$priorInc[i],
      postIncP = res$inc.prob[i],
      BFinc = res$BFinc[i],
      mean = res$mean[i],
      sd = res$sd[i],
      lowerCri = res$lo_ci[i],
      upperCri = res$hi_ci[i]
    )

    if(options$showCoefMeanInc){
      row <- append(row,c(meanInc = res$mean.inc[i],
                          sdInc = res$sd.inc[i]) )
    }

    bstsCoefficientTable$addRows(row)

  }


}


.bstsCreateStatePlots <- function(jaspResults,dataset,options,ready) {


  if (!is.null(jaspResults[["bstsStatePlots"]])) return()


  bstsStatePlots <- createJaspContainer(title = gettext("State Plots"))


  bstsStatePlots$dependOn(.bstsStatePlotDependencies())


  bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object

  if(options$aggregatedStatesPlot) .bstsAggregatedStatePlot(bstsStatePlots,bstsResults,dataset,options,ready)
  if(options$componentStatesPlot)
    .bstsComponentStatePlot(bstsStatePlots,bstsResults,options,ready)


  jaspResults[["bstsMainContainer"]][["bstsStatePlots"]] <- bstsStatePlots
  return()
}

.bstsAggregatedStatePlot <- function(bstsStatePlots,bstsResults,dataset,options,ready) {
  if (!ready | !options$aggregatedStatesPlot) return()


  bstsAggregatedStatePlot <- createJaspPlot(title= gettext("Aggregated State"), height = 320, width = 480)


  # get all states
  state <- bstsResults$state.contribution
  # discard burn ins
  state <- state[-(1:options$burn), , , drop = FALSE]
  # sum to final state
  state <- rowSums(aperm(state, c(1, 3, 2)), dims = 2)
  actualValues <- as.numeric(bstsResults$original.series)



  ymin <- apply(state,2,quantile,probs= ((1- options$aggregatedStatesPlotCiLevel)/2))
  ymax <- apply(state,2,quantile,probs= 1-((1- options$aggregatedStatesPlotCiLevel)/2))

  time <- options$time

  mean <-colMeans(state)


  p <-  ggplot2::ggplot(NULL,ggplot2::aes(x=time,y=mean)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=ymin,ymax =ymax),
                         fill ="blue",alpha=0.5) + ggplot2::xlab(gettext("Time")) +
    ggplot2::ylab(gettext("Distribution")) +
    ggplot2::geom_line(size=0.7)


  if(options$aggregatedStatesPlotObservationsShown)
    p <- p + ggplot2::geom_point(ggplot2::aes(y=actualValues))

  p <- jaspGraphs::themeJasp(p)

  bstsAggregatedStatePlot$plotObject <- p

  bstsStatePlots[["bstsAggregatedStatePlot"]] <- bstsAggregatedStatePlot

  return()
}

.bstsComponentStatePlot <- function(bstsStatePlots,bstsResults,options,ready){

  bstsComponentStatePlot <- createJaspPlot(title= gettext("Component States"), height = 320, width = 480)

  means <- as.data.frame(apply(bstsResults$state.contribution, 2, colMeans))

  means$time <- options$time

  ymin <- as.data.frame(apply(bstsResults$state.contribution, 2,matrixStats::colQuantiles,probs=c(0.025)))
  ymax <- as.data.frame(apply(bstsResults$state.contribution, 2,matrixStats::colQuantiles,probs=c(0.975)))
  ymin$time <- options$time
  ymax$time <- options$time

  ymin <- reshape2::melt(ymin,id.vars='time',variable.name="component")
  ymax <- reshape2::melt(ymax,id.vars='time',variable.name="component")
  means2 <- reshape2::melt(means,id.vars='time',variable.name="component")

  p <-  ggplot2::ggplot(data=means2, ggplot2::aes(x=time, y=value)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ymin$value,ymax=ymax$value),
                         fill ="blue",alpha=0.5) +
    ggplot2::theme_bw() + ggplot2::theme(legend.title = ggplot2::element_blank()) + ggplot2::ylab("") + ggplot2::xlab("")   +
    ggplot2::facet_grid(component ~ ., scales="free") +
    ggplot2::guides(colour=FALSE) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = -90, hjust = 0))


  bstsComponentStatePlot$plotObject <- p

  bstsStatePlots[["bstsComponentStatePlot"]] <- bstsComponentStatePlot

  return()
}

.bstsCreatePredictionPlot <- function(jaspResults,options,ready) {
  if(!ready | !options$predictionHorizon>0) return()

  bstsPredictionPlot <- createJaspPlot(title=gettext("Prediction plot"), height = 320, width = 480)
  bstsPredictionPlot$dependOn(.bstsPredictionDependencies())

  bstsModelResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object

  bstsPredictionResults <- jaspResults[["bstsMainContainer"]][["bstsModelPredictions"]]$object

  predictionRange <- (length(bstsPredictionResults$original.series)+1):(length(bstsPredictionResults$original.series)+options$predictionHorizon)

  predDF <- data.frame(time=seq_along(bstsPredictionResults$original.series),mean=bstsPredictionResults$original.series)
  predDF[predictionRange,"mean"] <-bstsPredictionResults$mean
  predDF[predictionRange,"LL"] <- bstsPredictionResults$interval[1,]
  predDF[predictionRange,"UL"] <- bstsPredictionResults$interval[2,]
  predDF[predictionRange,"time"] <- predictionRange

  p <- ggplot2::ggplot(data=predDF,ggplot2::aes(x=time,y=mean)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=LL,ymax=UL),
                         fill ="blue",alpha=0.5) +
    ggplot2::xlab(gettext("Time")) +
    ggplot2::geom_line(size=0.7) +
    ggplot2::geom_vline(xintercept=length(bstsPredictionResults$original.series),linetype=2) +
    ggplot2::theme_classic()

  p <- jaspGraphs::themeJasp(p)

  bstsPredictionPlot$plotObject <- p
  jaspResults[["bstsMainContainer"]][["bstsPredictionPlot"]] <- bstsPredictionPlot

  return()
}
# control plot

.bstsCreateControlPlots <- function(jaspResults,options,ready){
  if (!is.null(jaspResults[["bstsControlPlots"]]) | !ready) return()


  bstsControlPlots <- createJaspContainer(title = gettext("Control Plots"))

  bstsControlPlots$dependOn(.bstsControlDependencies())
  bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object


  if(options$controlChartPlot)
    .bstsFillControlPlotThreshold(bstsControlPlots,bstsResults,options,ready)

  if(options$probalisticControlPlot)
    .bstsFillControlPlotProbability(bstsControlPlots,bstsResults,options,ready)


  jaspResults[["bstsMainContainer"]][["bstsControlPlots"]] <- bstsControlPlots
  return()

}

.bstsFillControlPlotThreshold <- function(bstsControlPlots,bstsResults,options,ready){

  control = c(1:options$controlPeriod)
  L = options$controlSigma
  CI=.95
  bstsControlPlotThreshold <- createJaspPlot(title=gettext("Threshold Plot"), height = 320, width = 480)

  # extract model components
  state <- bstsResults$state.contribution
  # discard burn ins
  state <- state[-(1:options$burn), , , drop = FALSE]
  # sum to final state
  state <- rowSums(aperm(state, c(1, 3, 2)), dims = 2)

  # determine 0.95% credible interval
  ymin <- apply(state,2,quantile,probs= ((1- CI)/2))
  ymax <- apply(state,2,quantile,probs= 1-((1- CI)/2))

  time <- options$time
  mean <-colMeans(state)

  # compute limits based on the estimated state
  state_mean <- mean(mean[control])
  state_sd <- sd(mean[control])


  ul_state <- state_mean + L* state_sd
  ll_state <- state_mean - L* state_sd

  # ymin is lower 95 cI bound of state
  # ul is the threshold that has to be exceeded according to 2 sigmas
  # thus ymin must be larger  than the ul threshold
  CI_reached <- ymin > ul_state # cases where 95% confident that above threshold

  p <- ggplot2::ggplot(NULL,ggplot2::aes(x=time,y=mean))+
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=ymin,ymax =ymax ),
                         alpha=0.5,fill ="blue") + ggplot2::xlab(gettext("Time")) +
    ggplot2::ylab(gettext("Distribution")) +
    ggplot2::geom_line(size=0.7) +
    ggplot2::theme_classic() +
    ggplot2::geom_hline(yintercept=ul_state, linetype="dashed", color = "red") +
    ggplot2::geom_hline(yintercept=ll_state, linetype="dashed", color = "red")

  first_threshold <- time[which(ymin > ul_state)[1]]

    p <- jaspGraphs::themeJasp(p)

  if (!is.na(first_threshold)) {
    captionText <- if (options$time0 != "")
      gettextf("First date where %1$.3f*sigma is exceeded by %2$.1f%% state credible interval: %3$s", L, 100 * CI, first_threshold)
    else
      gettextf("First time point where %1$.3f*sigma is exceeded by %2$.1f%% state credible interval: %3$s", L, 100 * CI, first_threshold)
  } else {
    captionText <- gettext("Threshold never reached")
  }

  p <- p + ggplot2::labs(caption = captionText)



  p <- p+ ggplot2::theme(plot.caption = ggplot2::element_text(face="italic",hjust = 0,size=12)) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 15, 10, 15))

  bstsControlPlotThreshold$plotObject <- p



  bstsControlPlots[["bstsControlPlotThreshold"]] <- bstsControlPlotThreshold
  return()
}

.bstsFillControlPlotProbability <- function(bstsControlPlots,bstsResults,options,ready){

  bstsControlPlotProbability <- createJaspPlot(title=gettext("Probability Plot"), height = 320, width = 480)


  control = c(1:options$controlPeriod)
  L = options$controlSigma
  CI=.95 # add option in JASP later
  bstsControlPlotThreshold <- createJaspPlot(title=gettext("Threshold Plot"))

  # extract model components
  state <- bstsResults$state.contribution
  # discard burn ins
  state <- state[-(1:options$burn), , , drop = FALSE]
  # sum to final state
  state <- rowSums(aperm(state, c(1, 3, 2)), dims = 2)

  # determine 0.95% credible interval
  ymin <- apply(state,2,quantile,probs= ((1- CI)/2))
  ymax <- apply(state,2,quantile,probs= 1-((1- CI)/2))

  time <- options$time
  mean <-colMeans(state)

  # compute limits based on the estimated state
  state_mean <- mean(mean[control])
  state_sd <- sd(mean[control])


  ul_state <- state_mean + L* state_sd
  ll_state <- state_mean - L* state_sd

  # ymin is lower 95 cI bound of state
  # ul is the threshold that has to be exceeded according to 2 sigmas
  # thus ymin must be larger  than the ul threshold
  CI_reached <- ymin > ul_state # cases where 95% confident that above threshold

  # now let's compute probability of exeeding the threshold by getting the percentile rank from the MCMC draws

  threshold_prob <- 1-apply(state, 2, quantInv,ul_state)



  p <- ggplot2::ggplot(NULL,ggplot2::aes(time,threshold_prob)) + ggplot2::xlab(gettext("Time")) +
    ggplot2::ylab(gettext("Probability")) + ggplot2::geom_line() + ggplot2::theme_classic() + ggplot2::ylim(0,1)


  p <- jaspGraphs::themeJasp(p) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 15, 10, 15))

  bstsControlPlotProbability$plotObject <- p

  bstsControlPlots[["bstsControlPlotProbability"]] <- bstsControlPlotProbability
  return()
}

.bstsCreateErrorPlots <- function(jaspResults,options,ready) {
  if (!is.null(jaspResults[["bstsErrorPlots"]])) return()


  bstsErrorPlots <- createJaspContainer(title = gettext("Error Plots"))

  bstsResults <- jaspResults[["bstsMainContainer"]][["bstsModelResults"]]$object

  if(options$residualPlot) .bstsResidualPlot(bstsErrorPlots,bstsResults,options,ready)
  if(options$forecastErrorPlot)
    .bstsForecastErrorPlot(bstsErrorPlots,bstsResults,options,ready)


  jaspResults[["bstsMainContainer"]][["bstsErrorPlots"]] <- bstsErrorPlots
  return()
}


.bstsResidualPlot <- function(bstsErrorPlots,bstsResults,options,ready){
  if (!ready | !options$residualPlot) return()


  bstsResidualPlot <- createJaspPlot(title= gettext("Residual Plot"), height = 320, width = 480)

  residuals <- bsts::residuals.bsts(bstsResults,options$burn)

  mean <- colMeans(residuals)
  ymin <- apply(residuals,2,quantile,probs= ((1- options$aggregatedStatesPlotCiLevel)/2))
  ymax <- apply(residuals,2,quantile,probs= 1-((1- options$aggregatedStatesPlotCiLevel)/2))

  time <- options$time

  p <-  ggplot2::ggplot(NULL,ggplot2::aes(x=time,y=mean)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=ymin,ymax =ymax),
                         fill ="blue",alpha=0.5) + ggplot2::xlab(gettext("Time")) +
    ggplot2::ylab(gettext("Distribution")) +
    ggplot2::geom_line(size=0.7)

  p <- jaspGraphs::themeJasp(p)

  bstsResidualPlot$plotObject <- p

  bstsErrorPlots[["bstsResidualPlot"]] <- bstsResidualPlot

  return()
}

.bstsForecastErrorPlot <- function(bstsErrorPlots,bstsResults,options,ready){
  if (!ready | !options$forecastErrorPlot) return()


  bstsForecastErrorPlot <- createJaspPlot(title= gettext("Forecast Error Plot"), height = 320, width = 480)

  errors <- bsts::bsts.prediction.errors(bstsResults,burn=options$burn)$in.sample

  mean <- colMeans(errors)
  ymin <- apply(errors,2,quantile,probs= ((1- options$aggregatedStatesPlotCiLevel)/2))
  ymax <- apply(errors,2,quantile,probs= 1-((1- options$aggregatedStatesPlotCiLevel)/2))

  time <- options$time

  p <-  ggplot2::ggplot(NULL,ggplot2::aes(x=time,y=mean)) +
    ggplot2::geom_ribbon(mapping=ggplot2::aes(ymin=ymin,ymax =ymax),
                         fill ="blue",alpha=0.5) + ggplot2::xlab(gettext("Time")) +
    ggplot2::ylab(gettext("Distribution")) +
    ggplot2::geom_line(size=0.7)

  p <- jaspGraphs::themeJasp(p)

  bstsForecastErrorPlot$plotObject <- p

  bstsErrorPlots[["bstsForecastErrorPlot"]] <- bstsForecastErrorPlot

  return()
}

