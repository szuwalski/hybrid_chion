# Survey selectivity utilities for the two-fishery sensitivity model.
stacked_survey_curve <- function(sizes, theta) {
  z <- sizes / 30
  midpoint1 <- theta[1]
  midpoint2 <- midpoint1 + exp(theta[2])
  weight <- plogis(theta[5])
  weight * plogis(exp(theta[3]) * (z - midpoint1)) +
    (1 - weight) * plogis(exp(theta[4]) * (z - midpoint2))
}

survey_curve <- function(b, theta) {
  if (!is.null(b$data$simple_survey_selectivity) &&
      b$data$simple_survey_selectivity == 1) {
    plogis(exp(theta[3]) * (b$data$sizes / 30 - theta[1]))
  } else {
    stacked_survey_curve(b$data$sizes, theta)
  }
}

stacked_survey_matrix <- function(b, theta) {
  pre <- survey_curve(b, theta[, 1])
  post <- survey_curve(b, theta[, 2])
  outer(as.numeric(b$data$years < 1988), pre) +
    outer(as.numeric(b$data$years >= 1988), post)
}

stacked_survey_prior_nll <- function(b, theta) {
  prior <- b$data$stacked_survey_prior
  sd <- b$data$stacked_survey_prior_sd
  use <- is.finite(prior) & is.finite(sd) & sd > 0
  -sum(dnorm(survey_curve(b, theta[, 1])[use], prior[use], sd[use], log = TRUE)) -
    sum(dnorm(survey_curve(b, theta[, 2])[use], prior[use], sd[use], log = TRUE))
}

stacked_survey_start <- function(b) {
  prior <- b$data$stacked_survey_prior
  sd <- b$data$stacked_survey_prior_sd
  use <- is.finite(prior) & is.finite(sd) & sd > 0
  if (!is.null(b$data$simple_survey_selectivity) &&
      b$data$simple_survey_selectivity == 1) {
    loss <- function(p) {
      pred <- plogis(exp(p[2]) * (b$data$sizes / 30 - p[1]))
      sum(((pred[use] - prior[use]) / sd[use])^2)
    }
    best <- optim(c(100 / 30, log(3)), loss, method = "BFGS")$par
    one <- c(best[1], 0, best[2], 0, 0)
  } else {
    loss <- function(p) sum(((stacked_survey_curve(b$data$sizes, p)[use] - prior[use]) / sd[use])^2)
    starts <- lapply(c(.3, .5, .7), function(w)
      c(40 / 30, log(60 / 30), log(3), log(2), qlogis(w)))
    fits <- lapply(starts, function(p) optim(p, loss, method = "BFGS", control = list(maxit = 1000)))
    one <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]$par
  }
  matrix(rep(one, 2), 5, 2,
         dimnames = list(c("midpoint1_30mm", "log_gap_30mm", "log_slope1",
                           "log_slope2", "logit_weight"), c("pre1988", "post1988")))
}

