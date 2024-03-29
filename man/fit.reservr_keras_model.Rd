% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tf_fit.R
\name{fit.reservr_keras_model}
\alias{fit.reservr_keras_model}
\title{Fit a neural network based distribution model to data}
\usage{
\method{fit}{reservr_keras_model}(
  object,
  x,
  y,
  batch_size = NULL,
  epochs = 10,
  verbose = getOption("keras.fit_verbose", default = 1),
  callbacks = NULL,
  view_metrics = getOption("keras.view_metrics", default = "auto"),
  validation_split = 0,
  validation_data = NULL,
  shuffle = TRUE,
  class_weight = NULL,
  sample_weight = NULL,
  initial_epoch = 0,
  steps_per_epoch = NULL,
  validation_steps = NULL,
  ...
)
}
\arguments{
\item{object}{A compiled \code{reservr_keras_model} as obtained by \code{\link[=tf_compile_model]{tf_compile_model()}}.}

\item{x}{A list of input tensors (predictors)}

\item{y}{A \code{trunc_obs} tibble of observed outcomes, or something convertible via \code{\link[=as_trunc_obs]{as_trunc_obs()}}.}

\item{batch_size}{Integer or \code{NULL}. Number of samples per gradient update.
If unspecified, \code{batch_size} will default to 32.}

\item{epochs}{Number of epochs to train the model. Note that in conjunction
with \code{initial_epoch}, \code{epochs} is to be understood as "final epoch". The
model is not trained for a number of iterations given by \code{epochs}, but
merely until the epoch of index \code{epochs} is reached.}

\item{verbose}{Verbosity mode (0 = silent, 1 = progress bar, 2 = one line
per epoch). Defaults to 1 in most contexts, 2 if in knitr render or running
on a distributed training server.}

\item{callbacks}{List of callbacks to be called during training.}

\item{view_metrics}{View realtime plot of training metrics (by epoch). The
default (\code{"auto"}) will display the plot when running within RStudio,
\code{metrics} were specified during model \code{\link[keras:compile]{compile()}}, \code{epochs > 1} and
\code{verbose > 0}. Use the global \code{keras.view_metrics} option to establish a
different default.}

\item{validation_split}{Float between 0 and 1. Fraction of the training data
to be used as validation data. The model will set apart this fraction of
the training data, will not train on it, and will evaluate the loss and any
model metrics on this data at the end of each epoch. The validation data is
selected from the last samples in the \code{x} and \code{y} data provided, before
shuffling.}

\item{validation_data}{Data on which to evaluate the loss and any model
metrics at the end of each epoch. The model will not be trained on this
data. This could be a list (x_val, y_val) or a list (x_val, y_val,
val_sample_weights). \code{validation_data} will override \code{validation_split}.}

\item{shuffle}{shuffle: Logical (whether to shuffle the training data before
each epoch) or string (for "batch"). "batch" is a special option for
dealing with the limitations of HDF5 data; it shuffles in batch-sized
chunks. Has no effect when \code{steps_per_epoch} is not \code{NULL}.}

\item{class_weight}{Optional named list mapping indices (integers) to a
weight (float) value, used for weighting the loss function (during training
only). This can be useful to tell the model to "pay more attention" to
samples from an under-represented class.}

\item{sample_weight}{Optional array of the same length as x, containing
weights to apply to the model's loss for each sample. In the case of
temporal data, you can pass a 2D array with shape (samples,
sequence_length), to apply a different weight to every timestep of every
sample. In this case you should make sure to specify
\code{sample_weight_mode="temporal"} in \code{\link[keras:compile]{compile()}}.}

\item{initial_epoch}{Integer, Epoch at which to start training (useful for
resuming a previous training run).}

\item{steps_per_epoch}{Total number of steps (batches of samples) before
declaring one epoch finished and starting the next epoch. When training
with input tensors such as TensorFlow data tensors, the default \code{NULL} is
equal to the number of samples in your dataset divided by the batch size,
or 1 if that cannot be determined.}

\item{validation_steps}{Only relevant if \code{steps_per_epoch} is specified.
Total number of steps (batches of samples) to validate before stopping.}

\item{...}{Unused. If old arguments are supplied, an error message will be raised informing how to fix the issue.}
}
\value{
A \code{history} object that contains all information collected during training.
The model object will be updated in-place as a side-effect.
}
\description{
This function delegates most work to \code{\link[keras:fit.keras.engine.training.Model]{keras::fit.keras.engine.training.Model()}} and performs additional consistency
checks to make sure \code{\link[=tf_compile_model]{tf_compile_model()}} was called with the appropriate options to support fitting the observations
\code{y} as well as automatically converting \code{y} to a n x 6 matrix needed by the compiled loss function.
}
\details{
Additionally, the default \code{batch_size} is \code{min(nrow(y), 10000)} instead of keras default of \code{32} because the latter
is a very bad choice for fitting most distributions since the involved loss is much less stable than typical losses
used in machine learning, leading to divergence for small batch sizes.
}
\examples{
dist <- dist_exponential()
params <- list(rate = 1.0)
N <- 100L
rand_input <- runif(N)
x <- dist$sample(N, with_params = params)

if (interactive() && keras::is_keras_available()) {
  tf_in <- keras::layer_input(1L)
  mod <- tf_compile_model(
    inputs = list(tf_in),
    intermediate_output = tf_in,
    dist = dist,
    optimizer = keras::optimizer_adam(),
    censoring = FALSE,
    truncation = FALSE
  )

  tf_fit <- fit(
    object = mod,
    x = k_matrix(rand_input),
    y = x,
    epochs = 10L,
    callbacks = list(
      callback_debug_dist_gradients(mod, k_matrix(rand_input), x, keep_grads = TRUE)
    )
  )
}

}
\seealso{
predict.reservr_keras_model tf_compile_model keras::fit.keras.engine.training.Model
}
