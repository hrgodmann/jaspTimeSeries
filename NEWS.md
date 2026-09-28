# jaspTimeSeries Changelog

> **HOW TO READ AND UPDATE THIS CHANGELOG:**
>
> This document follows a modified [Keep a Changelog](https://keepachangelog.com/) format adapted for the R/JASP ecosystem. Releases are listed in reverse chronological order (newest first).
> As an example see [jaspModuleTemplate](https://github.com/jasp-stats/jaspModuleTemplate/blob/master/NEWS.md).
> * **Adding New Changes (For Contributors):** Add new entries under `# jaspTimeSeries (development version)` and the appropriate category (`## Added`, `## Changed`, `## Fixed`, etc.).
> * **Issue References:** Reference the relevant GitHub issue when one exists.
> * **Format Categories:**
>   * **Added:** New analyses, features, options, or output.
>   * **Changed:** Updates to existing analyses, defaults, dependencies, or output.
>   * **Fixed:** Bug fixes in analyses, plots, tables, help, QML layouts, or module infrastructure.
>   * **Deprecated / Removed:** Outdated analyses, options, or legacy code.

---
# jaspTimeSeries (development version)

## Added
* Integrated the Bayesian State Space Models analysis from the standalone `jaspBsts` module, originally developed by Fridtjof Petersen.
* Added Binomial State Space Models for success/trial counts, with optional continuous covariates, probability plots, and MCMC diagnostics. Covariate effects are additive and constant over time, with named summaries, diagnostics, and separate coefficient plots. Stan models are precompiled during package installation.

## Changed
* Increased the Binomial State Space Models default target acceptance probability to 0.97 and maximum tree depth to 15.
* Added Normal, Student t, and Cauchy coefficient-prior choices for Binomial State Space Models as radio buttons with separate scale/SD settings. Student t is the default, with scale 0.5 and editable degrees of freedom defaulting to 3; Normal SD defaults to 1 and Cauchy scale to 0.707. Earlier unreleased prototype analyses using the shared scale for Normal or Cauchy need their settings reviewed.
* Updated module metadata to use the project website and a consistent package title.
* Declared GPL-3-or-later for the combined module including its new RStan native infrastructure; existing source notices are retained.

## Fixed
* Synchronized the module version in `inst/Description.qml` with `DESCRIPTION`.
* Corrected the translation workflow to target the jaspTimeSeries Weblate components.
* Declared the directly used `tseries` package dependency.
