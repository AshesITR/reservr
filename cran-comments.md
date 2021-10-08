## Test environments

* GitHub Actions (ubuntu-20.04): devel, release, oldrel, 3.6
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

    checking compilation flags used ... NOTE
      Compilation used the following non-portable flag(s):
        ‘-march=x86-64’

This seems to be a local configuration issue and should not affect installation.

    checking for detritus in the temp directory ... NOTE
      Found the following files/directories:
        ‘__pycache__’ ‘tmp1cmt78c8.py’ ‘tmp95ck51wz.py’ ‘tmpmyliaqr1.py’
        ‘tmpo73o7yd8.py’ ‘tmpvpukjl38.py’ ‘tmpyt71tn9h.py’

This seems to be caused by running the tensorflow integration tests.

0 errors | 0 warnings | 4 notes

* This is a new release.
