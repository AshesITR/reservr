## Test environments

* GitHub Actions (ubuntu-latest): devel, release, oldrel, oldrel-1, oldrel-2
* GitHub Actions (windows): release
* GitHub Actions (macOS): release
* win-builder: devel

## R CMD check results

    checking for GNU extensions in Makefiles ... NOTE
      GNU make is a SystemRequirements.

This is due to the usage of `RcppParallel`.

0 errors | 0 warnings | 1 note

* This is a bugfix release, fixing some problems with the tensorflow integration when trying to train a `reservr_keras_model` from arbitrary distributions and observations. Regression tests were added for all distributions.
