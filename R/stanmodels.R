# Copyright (C) 2013-2026 University of Amsterdam
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Adapted from jaspAcceptanceSampling/R/asSSMStanLoader.R.
# Stan models are compiled when the package is installed from source.
# Analysis execution only loads those precompiled modules; it never invokes
# a native compiler. Model metadata is constructed lazily and cached.
# This is the module-owned loader, not rstantools-generated code. Keeping it
# at R/stanmodels.R lets rstan_config() preserve it during source installation.

.binomialSSM_getStanModel <- local({
  cache <- new.env(parent = emptyenv())

  .stanFile <- function(model_name) {
    path <- system.file("stan", paste0(model_name, ".stan"), package = "jaspTimeSeries")

    # Allow source-tree development only from this package's project root.
    if (!nzchar(path) && file.exists("DESCRIPTION")) {
      package_name <- unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L])
      if (identical(package_name, "jaspTimeSeries"))
        path <- file.path("inst", "stan", paste0(model_name, ".stan"))
    }

    if (!nzchar(path) || !file.exists(path))
      stop(gettextf("Stan model file not found for model: %1$s.", model_name), call. = FALSE)

    path
  }

  .precompiled <- function(model_name, stan_file) {
    expected_suffix <- paste0("_", model_name, "_mod")
    module_candidates <- paste0("stan_fit4", model_name, "_mod")

    reg <- tryCatch(getDLLRegisteredRoutines("jaspTimeSeries"), error = function(e) NULL)
    if (!is.null(reg) && !is.null(reg[[".Call"]])) {
      call_names <- names(reg[[".Call"]])
      boot_names <- grep("^_rcpp_module_boot_.*_mod$", call_names, value = TRUE)
      if (length(boot_names) > 0L) {
        discovered <- sub("^_rcpp_module_boot_", "", boot_names)
        discovered <- discovered[endsWith(discovered, expected_suffix)]
        if (length(discovered) > 0L)
          module_candidates <- unique(c(discovered, module_candidates))
      }
    }

    module_obj <- NULL
    module_class <- NULL

    for (module_name in module_candidates) {
      class_name <- NULL
      module_obj <- tryCatch(
        Rcpp::Module(module_name, PACKAGE = "jaspTimeSeries", mustStart = TRUE),
        error = function(e) NULL
      )
      if (is.null(module_obj))
        next

      module_print <- capture.output(print(module_obj))
      class_line <- grep("^\\s*\\d+\\s+classes\\s*:\\s*$", module_print)
      if (length(class_line) > 0L) {
        after <- trimws(module_print[(class_line[1L] + 1L):length(module_print)])
        after <- after[nzchar(after)]
        if (length(after) > 0L)
          class_name <- after[1L]
      }
      if (is.null(class_name) || !nzchar(class_name))
        class_name <- paste0("rstantools_model_", model_name)

      module_class <- tryCatch(
        eval(substitute(M$X, list(M = module_obj, X = as.name(class_name)))),
        error = function(e) NULL
      )
      if (!is.null(module_class))
        break
    }

    if (is.null(module_class))
      stop(gettextf("Precompiled Stan module class not available for model: %1$s.", model_name), call. = FALSE)

    # stanc_builder supplies metadata, not native compilation. The sampler
    # comes exclusively from the package's precompiled Rcpp module above.
    stanfit <- rstan::stanc_builder(
      stan_file,
      allow_undefined = TRUE,
      obfuscate_model_name = FALSE
    )
    stanfit$model_cpp <- list(
      model_cppname = stanfit$model_name,
      model_cppcode = stanfit$cppcode
    )

    methods::new(
      Class = "stanmodel",
      model_name = stanfit$model_name,
      model_code = stanfit$model_code,
      model_cpp = stanfit$model_cpp,
      mk_cppmodule = function(x) module_class
    )
  }

  function(model_name) {
    if (!model_name %in% c("binomialLocalLevel", "binomialLocalLevelRegression"))
      stop(gettextf("Unknown Stan model: %1$s.", model_name), call. = FALSE)

    if (exists(model_name, envir = cache, inherits = FALSE))
      return(get(model_name, envir = cache, inherits = FALSE))

    stan_file <- .stanFile(model_name)
    model <- tryCatch(
      .precompiled(model_name, stan_file),
      error = function(e) {
        stop(gettextf("Precompiled Stan model could not be loaded (%1$s): %2$s", model_name, conditionMessage(e)), call. = FALSE)
      }
    )

    assign(model_name, model, envir = cache)
    model
  }
})
