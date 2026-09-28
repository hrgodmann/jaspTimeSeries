# Binomial state-space native models

The two `.stan` files are the model sources. They are compiled into the
package's shared library during source installation, not when an analysis runs.

`configure` and `configure.win` call `rstantools::rstan_config()` to regenerate
the native wrappers and platform Makevars. Do not hand-edit generated files in
`src/`. After editing a Stan source, run `rstantools::rstan_config(".")` from
the package root with its development dependencies available and review the
resulting generated-source changes as well as the model source.

`R/stanmodels.R` is deliberately module-owned, adapted from Acceptance Sampling's
lazy loader. It has no rstantools generated-file marker, so `rstan_config()`
preserves it and reports an expected "not overwritten" warning for this file.
Do not replace it with the standard eager loader or add a second loader. Its
model names must agree with the two registered native modules. Model loading
may translate Stan source for metadata but never invokes a native compiler.

The generated infrastructure retains its GPL-3-or-later notice. The combined
package declares GPL-3-or-later; pre-existing source notices are preserved.
