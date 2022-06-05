compile_simple_function <- function(simple_fun, dist) {
  fmls <- formals(simple_fun)
  fcall <- call("simple_fun")
  fcall[1 + seq_along(fmls)] <- fmls
  names(fcall) <- c("", names(fmls))
  for (gen_arg in intersect(names(fmls), c("x", "n", "p", "q", "lower.tail", "log", "log.p"))) {
    fcall[[gen_arg]] <- as.name(gen_arg)
  }
  i <- 1L
  for (ph in names(dist$get_placeholders())) {
    fcall[[ph]] <- substitute(param_matrix[, i], list(i = i))
    i <- i + 1L
  }
  for (para in names(dist$default_params)) {
    if (!is.null(dist$default_params[[para]])) {
      fcall[[para]] <- dist$default_params[[para]]
    }
  }
  fmls_outer <- c(
    fmls[names(fmls) %in% c("x", "n", "p", "q")],
    alist(param_matrix = ),
    fmls[names(fmls) %in% c("lower.tail", "log", "log.p")]
  )
  as_compiled_distribution_function(
    as.function(c(fmls_outer, fcall)),
    i - 1L
  )
}

compile_simple_prob_interval_continuous <- function(fun, dist) {
  fmls <- formals(fun)
  fcall <- call("fun")
  fcall[1 + seq_along(fmls)] <- fmls
  names(fcall) <- c("", names(fmls))

  i <- 1L
  for (ph in names(dist$get_placeholders())) {
    fcall[[ph]] <- substitute(param_matrix[, i], list(i = i))
    i <- i + 1L
  }
  for (para in names(dist$default_params)) {
    if (!is.null(dist$default_params[[para]])) {
      fcall[[para]] <- dist$default_params[[para]]
    }
  }
  fcall[["lower.tail"]] <- TRUE
  fcall[["log.p"]] <- FALSE

  fcall_upper <- fcall
  fcall_upper[["q"]] <- as.name("qmax")
  fcall_lower <- fcall
  fcall_lower[["q"]] <- as.name("qmin")
  as_compiled_distribution_function(
    eval(bquote(function(qmin, qmax, param_matrix, log.p = FALSE) {
      prob <- .(fcall_upper) - .(fcall_lower)
      if (log.p) log(prob) else prob
    })),
    i - 1L
  )
}

compile_simple_prob_interval_discrete <- function(pfun, dfun, dist) {
  fmls <- formals(pfun)
  fcall <- call("pfun")
  fcall[1 + seq_along(fmls)] <- fmls
  names(fcall) <- c("", names(fmls))

  i <- 1L
  for (ph in names(dist$get_placeholders())) {
    fcall[[ph]] <- substitute(param_matrix[, i], list(i = i))
    i <- i + 1L
  }
  for (para in names(dist$default_params)) {
    if (!is.null(dist$default_params[[para]])) {
      fcall[[para]] <- dist$default_params[[para]]
    }
  }
  fcall[["lower.tail"]] <- TRUE
  fcall[["log.p"]] <- FALSE

  fmls_outer <- alist(qmin =, qmax =, param_matrix = )
  fcall_upper <- fcall
  fcall_upper[["q"]] <- as.name("qmax")
  fcall_lower <- fcall
  fcall_lower[["q"]] <- as.name("qmin")
  fcall_lower_d <- fcall
  fcall_lower_d[[1L]] <- as.name("dfun")
  names(fcall_lower_d)[2L] <- "x"
  fcall_lower_d[["x"]] <- as.name("qmin")
  fcall_lower_d[["lower.tail"]] <- NULL
  fcall_lower_d[["log.p"]] <- NULL
  fcall_lower_d[["log"]] <- FALSE

  as_compiled_distribution_function(
    eval(bquote(function(qmin, qmax, param_matrix, log.p = FALSE) {
      prob <- .(fcall_upper) - .(fcall_lower) + .(fcall_lower_d)
      if (log.p) log(prob) else prob
    })),
    i - 1L
  )
}

as_compiled_distribution_function <- function(fun, n_params) {
  fun <- as.function(fun)
  fun <- compiler::cmpfun(fun, options = list(optimize = 3L))
  n_params <- as.integer(n_params)
  class(fun) <- "compiled_distribution_function"
  attr(fun, "n_params") <- n_params
  fun
}
