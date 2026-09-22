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

## Changed
* Updated module metadata to use the project website and a consistent package title.

## Fixed
* Synchronized the module version in `inst/Description.qml` with `DESCRIPTION`.
* Corrected the translation workflow to target the jaspTimeSeries Weblate components.
* Declared the directly used `tseries` package dependency.
