## Test environments

* GitHub Actions (ubuntu-latest): devel, release, oldrel, oldrel-1, oldrel-2
* GitHub Actions (windows): release
* GitHub Actions (macOS): release
* win-builder: devel

## R CMD check results

    checking installed package size ... NOTE
      installed size is 13.7Mb
      sub-directories of 1Mb or more:
        R      1.5Mb
        libs  11.7Mb

The R source code is 409,1kB; I am not sure why installed size is almost triple that amount.
`libs` is large due to the size of the built C++ library (`reservr.so` is 8,6MB on my machine). 

    checking for GNU extensions in Makefiles ... NOTE
      GNU make is a SystemRequirements.

This is due to the usage of `RcppParallel`.

0 errors | 0 warnings | 1 note

* This is a bugfix release, fixing some problems with the tensorflow integration when trying to train a `reservr_keras_model` from arbitrary distributions and observations. Regression tests were added for all distributions.
