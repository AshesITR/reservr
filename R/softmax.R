#' Soft-Max function
#'
#' Softmax for a vector x is defined as
#'
#' \eqn{s_i = \exp(x_i) / \sum_k \exp(x_k)}
#'
#' It satisfies `sum(s) == 1.0` and can be used to smoothly enforce a sum
#' constraint.
#'
#' @param x A numeric vector or matrix
#'
#' @return `softmax` returns the softmax of `x`; rowwise if `x` is a matrix.
#' @export
#'
#' @examples
#' softmax(c(5, 5))
#' softmax(diag(nrow = 5, ncol = 6))
softmax <- function(x) {
  if (is.matrix(x)) {
    softmax_mat(x)
  } else {
    softmax_vec(x)
  }
}

#' @export
#' @rdname softmax
#' @return `dsoftmax` returns the Jacobi-matrix of `softmax(x)` at `x`. `x` must be a vector.
dsoftmax <- function(x) {
  dsoftmax_vec(x)
}
