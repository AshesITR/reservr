## Test environments

* GitHub Actions (ubuntu-latest): devel, release, oldrel, 3.6
* GitHub Actions (windows): release, 3.6
* GitHub Actions (macOS): release
* win-builder: devel

## R CMD check results

    checking for GNU extensions in Makefiles ... NOTE
      GNU make is a SystemRequirements.

This is due to the usage of `RcppParallel`.

0 errors | 0 warnings | 1 note

* This is a new release.
