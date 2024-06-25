# reservr (development version)

# reservr 0.0.3

* Fixed segfaults on r-devel caused by zero-length input to C++ routines.
* Migrated to `{keras3}` for keras support.

# reservr 0.0.2

* Fixed tensorflow log-density implementation for `dist_erlangmix()` and `dist_exponential()` to work with censored data.
* Multiple bug fixes related to tensorflow training integration, where input tensor shapes can be unknown.
* Improved testing of tensorflow integration.

# reservr 0.0.1

* Initial CRAN release
