create_qq_dataframe <- function(res, uniform = TRUE) {
  if (!inherits(res, "DHARMa")) {
    stop("Input must be a DHARMa simulationOutput object")
  }

  # Extract scaled residuals
  residuals <- res$scaledResiduals
  n <- length(residuals)

  # Create expected quantiles
  expected <- (1:n)/(n + 1)

  # Create dataframe
  qq_df <- data.frame(
    expected = expected,
    observed = residuals
  )

  # Sort values for proper QQ plotting
  sx <- sort(qq_df$expected)
  sy <- sort(qq_df$observed)

  # Handle different lengths if necessary
  lenx <- length(sx)
  leny <- length(sy)

  if (leny < lenx)
    sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx)
    sy <- approx(1L:leny, sy, n = lenx)$y

  # Convert to normal quantiles if requested
  if (!uniform) {
    sx <- qnorm(sx)
    sy <- qnorm(sy)
  }

  # Return final dataframe
  data.frame(expected = sx, observed = sy)
}

#' Create residual plots using ggplot from DHARMa objects
#'
#' @param res A DHARMa simulationOutput object from simulateResiduals()
#' @param uniform Whether to use uniform or normal quantiles in QQ plot (default: TRUE)
#' @param combine Whether to combine plots with patchwork (default: TRUE)
#' @return Either a combined plot object or a list of plots
#' @export
gg_dharma <- function(res, uniform = TRUE, combine = TRUE) {
  if (!inherits(res, "DHARMa")) {
    stop("Input must be a DHARMa simulationOutput object")
  }

  # Create residuals vs predicted plot data
  residuals <- res$scaledResiduals
  predicted <- res$fittedPredictedResponse
  predicted_rank <- rank(predicted, ties.method = "average")
  predicted_rank <- predicted_rank/max(predicted_rank)

  residual_df <- data.frame(
    residuals = residuals,
    predicted = predicted_rank
  )

  # Create QQ plot data
  qq_df <- create_qq_dataframe(res, uniform = uniform)

  # Create residuals vs predicted plot
  p1 <- ggplot(residual_df, aes(x = predicted, y = residuals)) +
    geom_point(alpha = .5, size = 0.5) +
    xlab("Predicted Value") +
    ylab("Residuals") +
    geom_hline(yintercept = .25, color = "black", alpha = 0.5) +
    geom_hline(yintercept = .50, color = "black", alpha = 0.5) +
    geom_hline(yintercept = .75, color = "black", alpha = 0.5) +
    geom_quantile(color = "tomato3", linetype = "dashed")

  # Create QQ plot
  p2 <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(alpha = .5, size = 0.5) +
    xlab("Expected") +
    ylab("Observed") +
    geom_abline(intercept = 0, slope = 1)

  # Add axis limits based on whether using uniform or normal quantiles
  if (uniform) {
    p2 <- p2 + xlim(0, 1) + ylim(0, 1)
  }

  # Return combined or separate plots
  if (combine) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      warning("patchwork package not available, returning plots as a list")
      return(list(qq_plot = p2, resid_plot = p1))
    }
    return(p2 + p1)
  } else {
    return(list(qq_plot = p2, resid_plot = p1))
  }
}
