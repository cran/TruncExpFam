# TruncExpFam 1.2.1

* Increased tolerance on unit test (issue #114)

# TruncExpFam 1.2.0

* Implemented `ptrunc()` and `qtrunc()` for all distributions (issue #54)
* Refactoring (issue #104, #112)
* Fixed bugs related to using the Negative Binomial with `mu` instead of `prob` (issue #107)
* Fixed domain validation on Negative Binomial and Inverse Gamma
* Added domain validation to `rtrunc(..., faster = TRUE)` (issue #109)
* Added `faster` argument to `rtrunc()` aliases (issue #110)
* Improved calculation of cumulative densities (issue #113)

# TruncExpFam 1.1.1

* Fixed `UseMethod()` no longer forwarding local variables from the generic on R-devel (issue #103)

# TruncExpFam 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added support for simple numeric-class vectors on `mlEstimationTruncDist()` (issue #95)
* Improved linting (issue #83)
* Improved unit test coverage (issues #66, #82)
* Improved output format (issues #89, #94)
* Improved documentation (issue #86)
* Improved ML estimation (issues #84, #85, #90, #92)
* Improved parameter parsing (issue #74)
* Simplified codebase (issues #91, #77, #78)
* Added a vignette (issue #53)
* Enabled usage of `parameters2natural()` and `natural2parameters()` with numeric vectors (issue #88)
* Faster calculation of `rtrunc()` for continuous distributions (issue #78)
* Faster calculation of `rtrunc()` for discrete distributions (issue #77)
* Created package page on https://ocbe-uio.github.io/TruncExpFam/ (issues #100, #101)
* Reclassified `eta` and `parms` arguments as `parms_*` instead of `trunc_*` (issue #97)

# TruncExpFam 1.0.1

* Hotfix for package documentation (https://github.com/r-lib/roxygen2/issues/1491)
* Improved argument consistency between generic and methods

# TruncExpFam 1.0.0

* CRAN debut
