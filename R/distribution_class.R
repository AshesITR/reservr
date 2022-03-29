#' Base class for Distributions
#'
#' Represents a modifiable Distribution family
#'
#' @examples
#' # Example for param_bounds:
#'
#' # Create an Exponential Distribution with rate constrained to (0, 2)
#' # instead of (0, Inf)
#' my_exp <- dist_exponential()
#' my_exp$param_bounds$rate <- interval(c(0, 2))
#' my_exp$get_param_bounds()
#'
#' fit_dist(my_exp, rexp(100, rate = 3), start = list(rate = 1))$params$rate
#'
#' @family Distributions
Distribution <- R6Class( # nolint: cyclocomp_linter.
  classname = "Distribution",
  public = list(
    #' @details Construct a Distribution instance
    #'
    #' Used internally by the `dist_*` functions.
    #'
    #' @param type Type of distribution. This is a string constant for the
    #' default implementation. Distributions with non-constant type must
    #' override the `get_type()` function.
    #' @param caps Character vector of capabilities to fuel the default
    #' implementations of `has_capability()` and `require_capability()`.
    #' Distributions with dynamic capabilities must override the
    #' `has_capability()` function.
    #' @param params Initial parameter bounds structure, backing the
    #' `param_bounds` active binding (usually a list of intervals).
    #' @param name Name of the Distribution class. Should be `CamelCase` and end
    #' with `"Distribution"`.
    #' @param default_params Initial fixed parameters backing the
    #' `default_params` active binding (usually a list of numeric / NULLs).
    initialize = function(type, caps, params, name, default_params) {
      private$.type <- type
      private$.caps <- caps
      private$.params <- params
      private$.default_params <- default_params
      private$.name <- name
    },

    # nocov start (stubs with stop() as implementation)

    #' @details Sample from a Distribution
    #'
    #' @param n number of samples to draw.
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `n`. In that
    #' case the `i`-th sample will use the `i`-th parameters.
    #'
    #' @return A length `n` vector of i.i.d. random samples from the
    #' Distribution with the specified parameters.
    #'
    #' @examples
    #' dist_exponential(rate = 2.0)$sample(10)
    sample = function(n, with_params = list()) {
      stop("sample() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Density of a Distribution
    #'
    #' @param x Vector of points to evaluate the density at.
    #' @param log Flag. If `TRUE`, return the log-density instead.
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(x)`.
    #' In that case, the `i`-th density point will use the `i`-th parameters.
    #'
    #' @return A numeric vector of (log-)densities
    #'
    #' @examples
    #' dist_exponential()$density(c(1.0, 2.0), with_params = list(rate = 2.0))
    density = function(x, log = FALSE, with_params = list()) {
      stop("density() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Compile a TensorFlow function for log-density evaluation
    #'
    #' @return A `tf_function` taking arguments `x` and `args` returning the
    #' log-density of the Distribution evaluated at `x` with parameters `args`.
    tf_logdensity = function() {
      stop("tf_logdensity() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Cumulative probability of a Distribution
    #'
    #' @param q Vector of points to evaluate the probability function at.
    #' @param lower.tail If `TRUE`, return P(X <= q). Otherwise return P(X > q).
    #' @param log.p If `TRUE`, probabilities are returned as `log(p)`.
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(q)`.
    #' In that case, the `i`-th probability point will use the `i`-th
    #' parameters.
    #'
    #' @return A numeric vector of (log-)probabilities
    #'
    #' @examples
    #' dist_exponential()$probability(
    #'   c(1.0, 2.0),
    #'   with_params = list(rate = 2.0)
    #' )
    probability = function(q, lower.tail = TRUE, log.p = FALSE,
                           with_params = list()) {
      stop("probability() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Compile a TensorFlow function for log-probability evaluation
    #'
    #' @return A `tf_function` taking arguments `qmin`, `qmax` and `args`
    #' returning the log-probability of the Distribution evaluated over the
    #' closed interval \[`qmin`, `qmax`\] with parameters `args`.
    tf_logprobability = function() {
      stop("tf_logprobability() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Quantile function of a Distribution
    #'
    #' @param p Vector of probabilities.
    #' @param lower.tail If `TRUE`, return P(X <= q). Otherwise return P(X > q).
    #' @param log.p If `TRUE`, probabilities are returned as `log(p)`.
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(p)`.
    #' In that case, the `i`-th quantile will use the `i`-th parameters.
    #'
    #' @return A numeric vector of quantiles
    #'
    #' @examples
    #' dist_exponential()$quantile(c(0.1, 0.5), with_params = list(rate = 2.0))
    quantile = function(p, lower.tail = TRUE, log.p = FALSE,
                        with_params = list()) {
      stop("quantile() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Hazard function of a Distribution
    #'
    #' @param x Vector of points.
    #' @param log Flag. If `TRUE`, return the log-hazard instead.
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(x)`.
    #' In that case, the `i`-th hazard point will use the `i`-th parameters.
    #'
    #' @return A numeric vector of (log-)hazards
    #'
    #' @examples
    #' dist_exponential(rate = 2.0)$hazard(c(1.0, 2.0))
    hazard = function(x, log = FALSE, with_params = list()) {
      stop("hazard() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Gradients of the density of a Distribution
    #'
    #' @param x Vector of points.
    #' @param log Flag. If `TRUE`, return the gradient of the log-density
    #' instead.
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(x)`.
    #' In that case, the `i`-th density point will use the `i`-th parameters.
    #'
    #' @return A list structure containing the (log-)density gradients of all
    #' free parameters of the Distribution evaluated at `x`.
    #'
    #' @examples
    #' dist_exponential()$diff_density(
    #'   c(1.0, 2.0),
    #'   with_params = list(rate = 2.0)
    #' )
    diff_density = function(x, log = FALSE, with_params = list()) {
      stop("diff_density() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Gradients of the cumulative probability of a Distribution
    #'
    #' @param q Vector of points to evaluate the probability function at.
    #' @param lower.tail If `TRUE`, return P(X <= q). Otherwise return P(X > q).
    #' @param log.p If `TRUE`, probabilities are returned as `log(p)`.
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(q)`.
    #' In that case, the `i`-th probability point will use the `i`-th
    #' parameters.
    #'
    #' @return A list structure containing the cumulative (log-)probability
    #' gradients of all free parameters of the Distribution evaluated at `q`.
    #'
    #' @examples
    #' dist_exponential()$diff_probability(
    #'   c(1.0, 2.0),
    #'   with_params = list(rate = 2.0)
    #' )
    diff_probability = function(q, lower.tail = TRUE, log.p = FALSE,
                                with_params = list()) {
      stop("diff_probability() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Determine if a value is in the support of a Distribution
    #'
    #' @param x Vector of points
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(x)`.
    #' In that case, the `i`-th point will use the `i`-th parameters.
    #'
    #' @return A logical vector with the same length as `x` indicating whether
    #' `x` is part of the support of the distribution given its parameters.
    #'
    #' @examples
    #' dist_exponential(rate = 1.0)$is_in_support(c(-1.0, 0.0, 1.0))
    is_in_support = function(x, with_params = list()) {
      stop("is_in_support() is not implemented for ", class(self)[1L], ".")
    },

    #' @details Determine the support of a Distribution
    #'
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of constant length `n`.
    #' In that case, the `i`-th result will use the `i`-th parameters.
    #'
    #' @return A length `n` list of [interval_union()] describing the support.
    #' Discrete values will be given as integer intervals whereas continuous intervals describe the continuous support.
    #'
    get_support = function(with_params = list()) {
      stop("get_support() is not implemented for ", class(self)[1L], ".")
    },

    # nocov end

    #' @details Determine if a value has positive probability
    #'
    #' @param x Vector of points
    #' @param with_params Distribution parameters to use.
    #' Each parameter value can also be a numeric vector of length `length(x)`.
    #' In that case, the `i`-th point will use the `i`-th parameters.
    #'
    #' @return A logical vector with the same length as `x` indicating whether
    #' there is a positive probability mass at `x` given the Distribution
    #' parameters.
    #'
    #' @examples
    #' dist_dirac(point = 0.0)$is_discrete_at(c(0.0, 1.0))
    is_discrete_at = function(x, with_params = list()) {
      if (self$is_continuous()) {
        rep_len(FALSE, length(x))
      } else {
        stop("is_discrete_at() is not implemented for ", class(self)[1L], ".") # nocov # nolint: line_length_linter.
      }
    },

    #' @details Compile a TensorFlow function for discrete support checking
    #'
    #' @return A `tf_function` taking arguments `x` and `args` returning whether
    #' the Distribution has a point mass at `x` given parameters `args`.
    tf_is_discrete_at = function() {
      check_installed("tensorflow")

      if (self$is_continuous()) {
        private$.tf_retrieve_or_call(
          "c_false",
          function() function(x, args) {
            tensorflow::tf$broadcast_to(FALSE, shape = tensorflow::tf$shape(x))
          }
        )
      } else {
        # discrete and mixed case needs to check support as well
        # nocov start
        stop("tf_is_discrete_at() is not implemented for ",
              class(self)[1L], ".")
        # nocov end
      }
    },

    #' @details
    #' Check if a capability is present
    #'
    #' @param caps Character vector of capabilities
    #' @return A logical vector the same length as `caps`.
    #'
    #' @examples
    #' dist_exponential()$has_capability("density")
    has_capability = function(caps) {
      caps %in% private$.caps
    },

    #' @details
    #' Get the type of a Distribution. Type can be one of `discrete`,
    #' `continuous` or `mixed`.
    #'
    #' @return A string representing the type of the Distribution.
    #'
    #' @examples
    #' dist_exponential()$get_type()
    #' dist_dirac()$get_type()
    #'
    #' dist_mixture(list(dist_dirac(), dist_exponential()))$get_type()
    #' dist_mixture(list(dist_dirac(), dist_binomial()))$get_type()
    get_type = function() {
      private$.type
    },

    #' @details Get the component Distributions of a transformed Distribution.
    #'
    #' @return A possibly empty list of Distributions
    #'
    #' @examples
    #' dist_trunc(dist_exponential())$get_components()
    #' dist_dirac()$get_components()
    #' dist_mixture(list(dist_exponential(), dist_gamma()))$get_components()
    get_components = function() {
      list()
    },

    #' @details
    #' Check if a Distribution is discrete, i.e. it has a density with respect
    #' to the counting measure.
    #'
    #' @return `TRUE` if the Distribution is discrete, `FALSE` otherwise.
    #' Note that mixed distributions are not discrete but can have point masses.
    #'
    #' @examples
    #' dist_exponential()$is_discrete()
    #' dist_dirac()$is_discrete()
    is_discrete = function() {
      identical(self$get_type(), "discrete")
    },

    #' @details
    #' Check if a Distribution is continuous, i.e. it has a density with respect
    #' to the Lebesgue measure.
    #'
    #' @return `TRUE` if the Distribution is continuous, `FALSE` otherwise.
    #' Note that mixed distributions are not continuous.
    #'
    #' @examples
    #' dist_exponential()$is_continuous()
    #' dist_dirac()$is_continuous()
    is_continuous = function() {
      identical(self$get_type(), "continuous")
    },

    #' @details
    #' Ensure that a Distribution has all required capabilities.
    #' Will throw an error if any capability is missing.
    #'
    #' @param caps Character vector of Capabilities to require
    #' @param fun_name Frienly text to use for generating the error message in
    #' case of failure.
    #'
    #' @return Invisibly `TRUE`.
    #'
    #' @examples
    #' dist_exponential()$require_capability("diff_density")
    require_capability = function(caps,
                                  fun_name = paste0(sys.call(-1)[[1]], "()")) {
      if (!all(self$has_capability(caps))) {
        missing_caps <- setdiff(caps, private$.caps)
        if (length(missing_caps) > 1) {
          missing_caps <- paste0(
            paste(head(missing_caps, -1), collapse = ", "),
            " and ",
            tail(missing_caps, 1)
          )
        }

        stop(
          fun_name, " requires missing capabilites.\n",
          private$.name, " doesn't provide ", missing_caps, "."
        )
      }
      invisible(TRUE)
    },

    #' @details
    #' Get the number of degrees of freedom of a Distribution family.
    #' Only parameters without a fixed default are considered free.
    #'
    #' @return An integer representing the degrees of freedom suitable e.g. for
    #' AIC calculations.
    #'
    #' @examples
    #' dist_exponential()$get_dof()
    #' dist_exponential(rate = 1.0)$get_dof()
    get_dof = function() {
      res <- private$.default_params

      sum_dof <- function(proto) {
        list_elems <- vapply(proto, is.list, logical(1))
        distr_elems <- vapply(proto, is.Distribution, logical(1))
        null_elems <- vapply(proto, is.null, logical(1))
        sum(null_elems) +
          sum(vapply(proto[list_elems], sum_dof, numeric(1))) +
          sum(vapply(
            proto[distr_elems], function(d) d$get_dof(), numeric(1)
          ))
      }

      sum_dof(res)
    },

    #' @details
    #' Get Placeholders of a Distribution family.
    #' Returns a list of free parameters of the family.
    #' Their values will be `NULL`.
    #'
    #' If the Distribution has Distributions as parameters, placeholders will be
    #' computed recursively.
    #'
    #' @return A named list containing any combination of (named or unnamed)
    #' lists and `NULL`s.
    #'
    #' @examples
    #' dist_exponential()$get_placeholders()
    #' dist_mixture(list(dist_dirac(), dist_exponential()))$get_placeholders()
    get_placeholders = function() {
      res <- private$.default_params

      prune <- function(proto) {
        null_elems <- vapply(proto, is.null, logical(1))
        list_elems <- vapply(proto, is.list, logical(1))
        distr_elems <- vapply(proto, is.Distribution, logical(1))
        keep <- null_elems | list_elems | distr_elems
        proto[list_elems] <- lapply(proto[list_elems], prune)
        proto[distr_elems] <- lapply(
          proto[distr_elems],
          function(d) d$get_placeholders()
        )
        proto[keep]
      }

      prune(res)
    },

    #' @details
    #' Get a full list of parameters, possibly including placeholders.
    #'
    #' @param with_params Optional parameter overrides with the same structure
    #' as `dist$get_params()`. Given Parameter values are expected to be length
    #' 1.
    #'
    #' @return A list representing the (recursive) parameter structure of the
    #' Distribution with values for specified parameters and `NULL` for free
    #' parameters that are missing both in the Distributions parameters and in
    #' `with_params`.
    #'
    #' @examples
    #' dist_mixture(list(dist_dirac(), dist_exponential()))$get_params(
    #'   with_params = list(probs = list(0.5, 0.5))
    #' )
    get_params = function(with_params = list()) {
      my_params <- private$.make_params(with_params, 1L)

      get_distr_params <- function(elem) {
        if (is.list(elem) &&
          length(elem) == 2L &&
          hasName(elem, "dist") &&
          hasName(elem, "params") &&
          is.Distribution(elem$dist)) {
          elem$dist$get_params(elem$params)
        } else if (is.list(elem)) {
          lapply(elem, get_distr_params)
        } else {
          elem
        }
      }

      lapply(my_params, get_distr_params)
    },

    #' @details Get a list of constant TensorFlow parameters
    #'
    #' @param with_params Optional parameter overrides with the same structure
    #' as `dist$tf_make_constants()`. Given Parameter values are expected to be
    #' length 1.
    #'
    #' @return A list representing the (recursive) constant parameters of the
    #' Distribution with values sprecified by parameters. Each constant is a
    #' TensorFlow Tensor of dtype `floatx`.
    tf_make_constants = function(with_params = list()) {
      check_installed(c("keras", "tensorflow"))
      my_params <- private$.make_params(with_params, 1L)

      get_tf_consts <- function(elem) {
        if (is.list(elem) &&
          length(elem) == 2L &&
          hasName(elem, "dist") &&
          hasName(elem, "params") &&
          is.Distribution(elem$dist)) {
          elem$dist$tf_make_constants(elem$params)
        } else if (is.list(elem)) {
          are_null <- vapply(elem, is.null, logical(1L))
          lapply(elem[!are_null], get_tf_consts)
        } else if (is.numeric(elem)) {
          keras::k_constant(elem)
        } else {
          # nocov start
          stop(
            "Unsupported parameter class ", class(params)[1L],
            " for tf_make_constants of ", class(self)[1L]
          )
          # nocov end
        }
      }

      get_tf_consts(my_params)
    },

    #' @details Compile distribution parameters into tensorflow outputs
    #'
    #' @param input A keras layer to bind all outputs to
    #' @param name_prefix Prefix to use for layer names
    #'
    #' @return A list with two elements
    #'
    #' * `outputs` a flat list of keras output layers, one for each parameter.
    #' * `output_inflater` a function taking keras output layers and
    #'   transforming them into a list structure suitable for passing to the
    #'   loss function returned by
    #'   [tf_compile_trunc_loss](tf_compile_trunc_loss(dist))
    tf_compile_params = function(input, name_prefix = "") {
      bounds <- self$get_param_bounds()
      out <- self$get_placeholders()

      for (ph in names(out)) {
        b <- bounds[[ph]]
        if (!is.Interval(b)) {
          stop("Unsupported parameter ", ph, " for ", class(self)[1L], ".") # nocov # nolint: line_length_linter.
        }

        layer_name <- paste0(name_prefix, ph)
        out[[ph]] <- b$tf_make_layer(input, layer_name)
      }

      ph_names <- names(out)

      list(
        outputs = out,
        output_inflater = eval(bquote(function(outputs) {
          if (!is.list(outputs)) outputs <- list(outputs)
          names(outputs) <- .(ph_names)
          outputs
        }))
      )
    },

    #' @details Get Interval bounds on all Distribution parameters
    #'
    #' @return A list representing the free (recursive) parameter structure of
    #' the Distribution with `Interval` objects as values representing the
    #' bounds of the respective free parameters.
    #'
    #' @examples
    #' dist_mixture(
    #'   list(dist_dirac(), dist_exponential()),
    #'   probs = list(0.5, 0.5)
    #' )$get_param_bounds()
    #'
    #' dist_mixture(
    #'   list(dist_dirac(), dist_exponential())
    #' )$get_param_bounds()
    #'
    #' dist_genpareto()$get_param_bounds()
    #' dist_genpareto1()$get_param_bounds()
    get_param_bounds = function() {
      proto <- private$.make_params(list(), 1L)

      get_bounds <- function(elem, ranges) {
        if (is.null(elem)) {
          ranges
        } else if (is.list(elem) &&
          length(elem) == 2L &&
          hasName(elem, "dist") &&
          hasName(elem, "params") &&
          is.Distribution(elem$dist)) {
          elem$dist$get_param_bounds()
        } else if (is.list(elem) && rlang::is_named(elem)) {
          mapply(get_bounds, elem,
                 ranges[names(elem)], SIMPLIFY = FALSE)
        } else if (is.list(elem)) {
          mapply(get_bounds, elem,
                 rep_len(ranges, length(elem)), SIMPLIFY = FALSE)
        } else {
          NULL
        }
      }

      prune <- function(elem) {
        if (is.list(elem)) {
          non_null <- !vapply(elem, is.null, logical(1L))
          lapply(elem[non_null], prune)
        } else if (is.Interval(elem)) {
          elem
        } else {
          stop("non-list and non-Interval encountered during pruning.") # nocov
        }
      }

      res <- get_bounds(elem = proto, ranges = private$.params[names(proto)])

      prune(res)
    },

    #' @details Get additional (non-linear) equality constraints on Distribution
    #' parameters
    #'
    #' @return `NULL` if the box constraints specified by
    #' `dist$get_param_bounds()` are sufficient, or a function taking full
    #' Distribution parameters and returning either a numeric vector
    #' (which must be 0 for valid parameter combinations) or a list with
    #' elements
    #'
    #'  * `constraints`: The numeric vector of constraints
    #'  * `jacobian`: The Jacobi matrix of the constraints with respect to the
    #'    parameters
    #'
    #' @examples
    #' dist_mixture(
    #'   list(dist_dirac(), dist_exponential())
    #' )$get_param_constraints()
    get_param_constraints = function() {
      NULL
    },

    #' @details Export sampling, density, probability and quantile functions
    #' to plain R functions
    #'
    #' Creates new functions in `envir` named `{r,d,p,q}<name>` which implement
    #' `dist$sample`, `dist$density`, `dist$probability` and `dist$quantile` as
    #' plain functions with default arguments specified by `with_params` or the
    #' fixed parameters.
    #'
    #' The resulting functions will have signatures taking all parameters as
    #' separate arguments.
    #'
    #' @param name common suffix of the exported functions
    #' @param envir Environment to export the functions to
    #' @param with_params Optional list of parameters to use as default values
    #' for the exported functions
    #'
    #' @return Invisibly `NULL`.
    #'
    #' @examples
    #' tmp_env <- new.env(parent = globalenv())
    #' dist_exponential()$export_functions(
    #'   name = "exp",
    #'   envir = tmp_env,
    #'   with_params = list(rate = 2.0)
    #' )
    #' evalq(
    #'   fitdistrplus::fitdist(rexp(100), "exp"),
    #'   envir = tmp_env
    #' )
    export_functions = function(name, envir = parent.frame(),
                                with_params = list()) {
      params <- self$get_params(with_params = with_params)
      params_flat <- tryCatch(
        flatten_params(params),
        # If flattening is impossible, allow no arguments
        error = function(e) numeric()
      )

      make_fun <- function(prefix, first_param, func_name,
                           general_params = list()) {
        if (length(params_flat)) {
          param_list <- lapply(
            names(params_flat),
            function(nm) substitute(as.name(nm), list(nm = nm))
          )
          names(param_list) <- names(params_flat)
          param_list <- do.call(call, c("c", param_list))
          param_list <- as.call(list(
            quote(reservr::inflate_params), param_list
          ))
        } else {
          param_list <- list()
        }

        ffmls <- c(alist(x = ), params_flat, general_params) # nolint lintr/#532
        names(ffmls)[1L] <- first_param

        general_params <- c(list(x = as.name(first_param)), general_params)
        names(general_params)[1L] <- first_param

        # Can't use bquote(..., splice = TRUE) for backward compatibility with R < 4.0
        # Instead we construct the spliced call manually
        spliced_call <- as.call(c(
          substitute(self$func_name, list(func_name = func_name)),
          general_params,
          alist(with_params = params)
        ))

        fbody <- bquote({
          params <- .(param_list)
          tryCatch(
            .(spliced_call),
            error = function(e) {
              warning("Error during evaluation; returning NaNs.\n", e)
              rep_len(NaN, length(first_param))
            }
          )
        })

        fun <- as.function(c(ffmls, fbody))
        fun_name <- paste0(prefix, name)
        message("Exported `", fun_name, "()`.")
        assign(fun_name, fun, envir = envir)
      }

      if (self$has_capability("density")) {
        make_fun("d", "x", "density", list(log = FALSE))
      }
      if (self$has_capability("sample")) {
        make_fun("r", "n", "sample")
      }
      if (self$has_capability("probability")) {
        make_fun("p", "q", "probability",
                 list(lower.tail = TRUE, log.p = FALSE))
      }
      if (self$has_capability("quantile")) {
        make_fun("q", "p", "quantile",
                 list(lower.tail = TRUE, log.p = FALSE))
      }

      invisible(NULL)
    }
  ),
  private = list(
    .make_params = function(with_params, n) {
      # Make a list of objects
      # 1. list(...) => list()
      # 2. numeric => numeric
      # 3. distribution => distribution + params

      make_params_elem <- function(proto, value) {
        if (is.list(proto)) {
          if (is.null(value)) {
            value <- vector("list", length(proto))
          } else if (length(value) != length(proto)) {
            value <- rep_len(value, length(proto))
          }
          mapply(make_params_elem, proto, value, SIMPLIFY = FALSE)
        } else if (is.Distribution(proto)) {
          if (is.list(value) && setequal(names(value), c("dist", "params")) && is.Distribution(value[["dist"]])) {
            # Prevent re-wrapping when called twice e.g. via hazard()
            value
          } else {
            list(
              dist = proto,
              params = value
            )
          }
        } else {
          if (!is.null(value)) {
            if (is.vector(value)) {
              rep_len(value, length.out = n)
            } else {
              value
            }
          } else {
            if (is.vector(proto)) {
              rep_len(proto, length.out = n)
            } else {
              proto
            }
          }
        }
      }

      res <- lapply(names(private$.params), function(nm) {
        make_params_elem(private$.default_params[[nm]], with_params[[nm]])
      })
      names(res) <- names(private$.params)

      res
    },
    .caps = character(),
    .type = character(),
    .params = list(),
    .default_params = list(),
    .name = "",
    .tf_functions = list(
      nograph = list(),
      float32 = list(),
      float64 = list()
    ),
    .tf_retrieve_or_call = function(name, impl) {
        check_installed(c("keras", "tensorflow"))

        cache_key <- if (tf_is_graph_disabled()) {
          "nograph"
        } else {
          keras::k_floatx()
        }

        res <- private$.tf_functions[[cache_key]][[name]]
        use_cache <- getOption("reservr.cache_tf_function", default = TRUE)

        if (is.null(res) ||
          (cache_key != "nograph" && reticulate::py_is_null_xptr(res)) ||
          !use_cache) {
          fn <- impl()
          assign("tf", tensorflow::tf, environment(fn))
          res <- maybe_tf_function(fn)

          private$.tf_functions[[cache_key]][[name]] <- res
        }

        res
      }
  ),
  active = list(
    #' @field default_params Get or set (non-recursive) default parameters of a
    #' Distribution
    default_params = function(value) {
      if (missing(value)) {
        private$.default_params
      } else {
        assert_that(
          is.list(value),
          setequal(names(private$.default_params), names(value)),
          msg = paste0(
            "`default_params` must be a list with names ",
            enumerate_strings(names(private$.default_params), quote = "'")
          )
        )
        # Perform bounds checking
        bounds <- private$.params
        for (param in names(private$.default_params)) {
          curr_bounds <- bounds[[param]]
          if (is.Interval(curr_bounds)) {
            assert_that(
              is.null(value[[param]]) || is.numeric(value[[param]]) &&
                all(curr_bounds$contains(value[[param]])),
              msg = sprintf(
                "`default_params$%s` must be numeric in %s, or NULL.",
                param, curr_bounds
              )
            )
          } else if (is.list(curr_bounds) &&
            length(curr_bounds) &&
            is.Interval(curr_bounds[[1L]])) {
            assert_that(
              is.list(value[[param]]),
              all(vapply(
                value[[param]],
                function(val) is.null(val) || (is.numeric(val) &&
                  all(curr_bounds[[1L]]$contains(val))),
                logical(1L))
              ),
              msg = sprintf(
                paste(
                  "`default_params$%s` must be a list of",
                  "numeric values in %s, or NULLs"
                ),
                param,
                curr_bounds
              )
            )
          }
        }
        # New defaults are assumed valid -> accept.
        # Also ensure the order of elements stays the same just in case.
        private$.default_params <- value[names(private$.default_params)]
      }
    },

    #' @field param_bounds Get or set (non-recursive) parameter bounds
    #' (box constraints) of a Distribution
    param_bounds = function(value) {
      if (missing(value)) {
        private$.params
      } else {
        assert_that(
          is.list(value),
          setequal(names(private$.params), names(value)),
          msg = paste0("`param_bounds` must be a list with names ",
                       enumerate_strings(names(private$.params), quote = "'"))
        )
        # Check type safety
        bounds <- private$.params

        check_types <- function(old, new, path = "") {
          expected <- if (is.list(old) &&
            length(old) == 1L &&
            is.null(names(old)) &&
            is.Interval(old[[1L]])) {
            "a list containing an Interval"
          } else if (is.list(old) && !length(old)) {
            "an empty list"
          } else if (is.list(old)) {
            paste0("a list with ", length(old), "elements")
          } else if (is.Interval(old)) {
            "an Interval"
          }

          ok <- if (is.list(old)) {
            if (length(old) == 1L &&
              is.null(names(old)) &&
              is.Interval(old[[1L]])) {
              length(new) == 1L &&
                is.null(names(new)) &&
                is.Interval(new[[1L]])
            } else {
              is.list(new) && length(new) == length(old)
            }
          } else if (is.Interval(old)) {
            is.Interval(new)
          } else {
            FALSE
          }

          if (!ok) {
            stop(
              "param_bounds", path, " must be ", expected, ".",
              call. = FALSE
            )
          }
          if (is.list(old) &&
            length(old) &&
            !is.Interval(old[[1L]])) {
            sub_names <- if (!is.null(names(old))) {
              paste(path, names(old), sep = "$")
            } else {
              paste0(path, "[[", seq_along(old), "]]")
            }
            if (!is.null(names(old))) {
              if (!setequal(names(old), names(new))) {
                stop(
                  "param_bounds", path, " must be a list with names ",
                  enumerate_strings(names(old), quote = "'"), ".",
                  call. = FALSE
                )
              }
              new <- new[names(old)]
            } else {
              names(new) <- NULL
            }
            new <- mapply(
              check_types,
              old = old,
              new = new,
              path = sub_names,
              SIMPLIFY = FALSE
            )
          }
          # With fixed names / ordering
          new
        }

        private$.params <- check_types(bounds, value)
      }
    }
  )
)

# nocov start
distribution_class <- function(
  name,
  type = c("continuous", "discrete"),
  params = list(),
  sample = NULL,
  density = NULL,
  probability = NULL,
  quantile = NULL,
  diff_density = NULL,
  diff_probability = NULL,
  hazard = function(x, log, params) {
    if (log) {
      self$density(x, log = TRUE, with_params = params) -
        self$probability(x, lower.tail = FALSE, log.p = TRUE,
                         with_params = params)
    } else {
      self$density(x, with_params = params) /
        self$probability(x, lower.tail = FALSE, with_params = params)
    }
  },
  support = I_REALS,
  is_discrete = if (type == "discrete") {
    function(x, params) {
      self$is_in_support(x, params)
    }
  } else if (type == "continuous") {
    function(x, params) {
      logical(length(x))
    }
  },
  tf_logdensity = NULL,
  tf_logprobability = NULL,
  tf_is_discrete_at = NULL,
  ...,
  active = list()
) {
  clsname <- paste0(
    toupper(substr(name, 1, 1)),
    substr(name, 2, nchar(name)),
    "Distribution"
  )
  type <- match.arg(type)
  caps <- c(
    character(),
    if (!missing(sample)) "sample" else NULL,
    if (!missing(density)) "density" else NULL,
    if (!missing(probability)) "probability" else NULL,
    if (!missing(quantile)) "quantile" else NULL,
    if (!missing(diff_density)) "diff_density" else NULL,
    if (!missing(diff_probability)) "diff_probability" else NULL,
    if (!missing(tf_logdensity)) "tf_logdensity" else NULL,
    if (!missing(tf_logprobability)) "tf_logprobability" else NULL,
    if (!missing(support) && is.Interval(support)) "support" else NULL
  )

  if (is.Interval(support)) {
    support_interval <- support

    is_in_support <- function(x, params) {
      support_interval$contains(x)
    }

    get_support <- if (type == "discrete") {
      function(params) list(discrete = list(support_interval))
    } else if (type == "continuous") {
      function(params) list(continuous = list(support_interval))
    }
  } else {
    is_in_support <- support

    get_support <- function(params) {
      super$get_support(params)
    }
  }

  R6Class(
    classname = clsname,
    inherit = Distribution,
    public = list(
      initialize = function(...) {
        super$initialize(
          type = type,
          caps = caps,
          params = params,
          name = name,
          default_params = list(...)
        )
      },
      sample = function(n, with_params = list()) {
        self$require_capability("sample")
        params <- private$.make_params(with_params, n)
        if (n == 0L) return(numeric())
        private$.sample_impl(n = n, params = params)
      },
      density = function(x, log = FALSE, with_params = list()) {
        self$require_capability("density")
        params <- private$.make_params(with_params, length(x))
        if (!length(x)) return(numeric())
        private$.density_impl(x = x, log = log, params = params)
      },
      probability = function(q, lower.tail = TRUE, log.p = FALSE,
                             with_params = list()) {
        self$require_capability("probability")
        params <- private$.make_params(with_params, length(q))
        if (!length(q)) return(numeric())
        private$.probability_impl(q = q, lower.tail = lower.tail, log.p = log.p,
                                  params = params)
      },
      quantile = function(p, lower.tail = TRUE, log.p = FALSE,
                          with_params = list()) {
        self$require_capability("quantile")
        params <- private$.make_params(with_params, length(p))
        if (!length(p)) return(numeric())
        private$.quantile_impl(p = p, lower.tail = lower.tail, log.p = log.p,
                               params = params)
      },
      hazard = function(x, log = FALSE, with_params = list()) {
        self$require_capability(c("density", "probability"))
        params <- private$.make_params(with_params, length(x))
        if (!length(x)) return(numeric())
        private$.hazard_impl(x = x, log = log, params = params)
      },
      diff_density = function(x, log = FALSE, with_params = list()) {
        self$require_capability("diff_density")
        params <- private$.make_params(with_params, length(x))
        if (!length(x)) return(empty_derivative(self$get_placeholders()))
        private$.diff_density_impl(x = x, vars = self$get_placeholders(),
                                   log = log, params = params)
      },
      diff_probability = function(q, lower.tail = TRUE, log.p = FALSE,
                                  with_params = list()) {
        self$require_capability("diff_probability")
        params <- private$.make_params(with_params, length(q))
        if (!length(q)) return(empty_derivative(self$get_placeholders()))
        private$.diff_probability_impl(q = q, vars = self$get_placeholders(),
                                       lower.tail = lower.tail, log.p = log.p,
                                       params = params)
      },
      is_in_support = function(x, with_params = list()) {
        params <- private$.make_params(with_params, length(x))
        if (!length(x)) return(logical())
        private$.is_in_support_impl(x = x, params = params)
      },
      get_support = function(with_params = list()) {
        if (length(with_params)) {
          max_length <- function(ll) {
            is_list <- vapply(ll, is.list, logical(1L))
            n_here <- max(lengths(ll)[!is_list])
            max(vapply(ll[is_list], max_length, logical(1L)), n_here)
          }
          n <- max(max_length(with_params), 1L)
        } else {
          n <- 1L
        }
        params <- private$.make_params(with_params, n)
        private$.get_support_impl(params = params)
      },
      is_discrete_at = function(x, with_params = list()) {
        params <- private$.make_params(with_params, length(x))
        if (!length(x)) return(logical())
        private$.is_discrete_impl(x = x, params = params)
      },
      tf_logdensity = function() {
        if (is.null(private$.tf_logdensity_impl)) {
          super$tf_logdensity()
        } else {
          private$.tf_retrieve_or_call(
            "tf_logdensity",
            private$.tf_logdensity_impl
          )
        }
      },
      tf_logprobability = function() {
        if (is.null(private$.tf_logprobability_impl)) {
          super$tf_logprobability()
        } else {
          private$.tf_retrieve_or_call(
            "tf_logprobability",
            private$.tf_logprobability_impl
          )
        }
      },
      tf_is_discrete_at = function() {
        if (is.null(private$.tf_is_discrete_at_impl)) {
          super$tf_is_discrete_at()
        } else {
          private$.tf_retrieve_or_call(
            "tf_is_discrete_at",
            private$.tf_is_discrete_at_impl
          )
        }
      },
      ...
    ),
    private = list(
      .sample_impl = sample,
      .density_impl = density,
      .probability_impl = probability,
      .quantile_impl = quantile,
      .hazard_impl = hazard,
      .diff_density_impl = diff_density,
      .diff_probability_impl = diff_probability,
      .is_in_support_impl = is_in_support,
      .get_support_impl = get_support,
      .is_discrete_impl = is_discrete,
      .tf_logdensity_impl = tf_logdensity,
      .tf_logprobability_impl = tf_logprobability,
      .tf_is_discrete_at_impl = tf_is_discrete_at
    ),
    active = active
  )
}

distribution_class_simple <- function(name,
                                      fun_name,
                                      type = c("continuous", "discrete"),
                                      params = list(),
                                      support = I_REALS,
                                      envir = parent.frame(),
                                      ...) {
  type <- match.arg(type)

  sample_fun <- get(paste0("r", fun_name), mode = "function", envir = envir)
  density_fun <- get(paste0("d", fun_name), mode = "function", envir = envir)
  probability_fun <- get(
    paste0("p", fun_name), mode = "function", envir = envir
  )
  quantile_fun <- get(paste0("q", fun_name), mode = "function", envir = envir)

  eval(substitute(
    distribution_class(
      name = name,
      type = type,
      params = params,
      sample = function(n, params) {
        do.call(sample_fun, c(list(n = n), params))
      },
      density = function(x, log, params) {
        do.call(density_fun, c(list(x = x, log = log), params))
      },
      probability = function(q, lower.tail, log.p, params) {
        do.call(probability_fun,
                c(list(q = q, lower.tail = lower.tail, log.p = log.p), params))
      },
      quantile = function(p, lower.tail, log.p, params) {
        do.call(quantile_fun,
                c(list(p = p, lower.tail = lower.tail, log.p = log.p), params))
      },
      support = support,
      ...
    ),
    env = list(
      sample_fun = sample_fun,
      density_fun = density_fun,
      probability_fun = probability_fun,
      quantile_fun = quantile_fun
    )
  ))
}
# nocov end

#' Test if object is a Distribution
#'
#' @param object An R object.
#'
#' @examples
#' is.Distribution(dist_dirac())
#'
#' @export
#' @return `TRUE` if `object` is a Distribution, `FALSE` otherwise.
is.Distribution <- function(object) {
  inherits(object, "Distribution")
}
