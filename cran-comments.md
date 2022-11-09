## Test environments

* GitHub Actions (ubuntu-latest): devel, release, oldrel, 3.6
* GitHub Actions (windows): release, 3.6
* GitHub Actions (macOS): release
* win-builder: devel

## R CMD check results

    checking R code for possible problems ... NOTE
    Found the following possibly unsafe calls:
    File ‘reservr/R/tf_compile.R’:
      assignInNamespace("tf_function", function(f, ...) f, ns = "tensorflow")
      assignInNamespace("tf_function", old_tf_function, ns = "tensorflow")

This is intentional behaviour of the function `tf_with_disable_graph()`.
It's self resetting via `on.exit()` hooks so should be safe.

    checking for GNU extensions in Makefiles ... NOTE
      GNU make is a SystemRequirements.

This is due to the usage of `RcppParallel`.

0 errors | 0 warnings | 2 notes

* This is a new release.
