library(RTMB)

source("models/model/R/joint_classification_model.R")

args <- commandArgs(trailingOnly = TRUE)
argval <- function(key, default) {
  z <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(z)) sub(paste0("^", key, "="), "", z[1]) else default
}
fit_path <- argval("--fit", "models/model/fits/model_fit.rds")
diagnostic_path <- argval("--diagnostics", "models/model/fits/hessian_diagnostics.csv")
weak_path <- argval("--weak", "models/model/fits/weak_hessian_direction.csv")

x <- readRDS(fit_path)
d <- x$data
p <- J1_parameters(d)
for (nm in intersect(names(p), names(x$par_list))) {
  if (length(p[[nm]]) == length(x$par_list[[nm]])) p[[nm]] <- x$par_list[[nm]]
}
obj <- MakeADFun(function(z) do.call(J1_objective(d), z), p,
                 map = x$map, silent = TRUE)

if (!identical(names(obj$par), names(x$fit$par))) {
  stop("Saved parameter layout does not match the current objective source")
}
theta <- x$fit$par
objective <- obj$fn(theta)
gradient <- obj$gr(theta)
cat("Computing", length(theta), "by", length(theta), "numerical Hessian\n")
H <- optimHess(theta, obj$fn, obj$gr)
H <- (H + t(H)) / 2
finite_hessian <- all(is.finite(H))
if (!finite_hessian) stop("Hessian contains non-finite entries")
eig <- eigen(H, symmetric = TRUE)
ev <- eig$values
chol_result <- try(chol(H), silent = TRUE)
chol_ok <- !inherits(chol_result, "try-error")
inverse_result <- try(solve(H), silent = TRUE)
inverse_ok <- !inherits(inverse_result, "try-error") && all(is.finite(inverse_result))

diagnostics <- data.frame(
  metric = c("parameters", "objective", "max_gradient", "finite_hessian",
             "chol_positive_definite", "direct_inverse_succeeded",
             "eigenvalues_below_1e-8", "eigenvalues_below_1e-6",
             "minimum_eigenvalue", "maximum_eigenvalue",
             "scaled_minimum_eigenvalue", "condition_number"),
  value = c(length(theta), objective, max(abs(gradient)), finite_hessian, chol_ok,
            inverse_ok, sum(ev < 1e-8), sum(ev < 1e-6),
            min(ev), max(ev), min(ev) / max(abs(ev)),
            if (min(ev) > 0) max(ev) / min(ev) else Inf)
)
x$report <- obj$report(theta)
x$hessian_diagnostics <- list(
  checked = Sys.time(), method = "optimHess plus symmetric eigenvalues/chol",
  positive_definite = chol_ok && min(ev) > 0,
  minimum_eigenvalue = min(ev), maximum_eigenvalue = max(ev),
  scaled_minimum_eigenvalue = min(ev) / max(abs(ev)),
  condition_number = if (min(ev) > 0) max(ev) / min(ev) else Inf
)
if (x$hessian_diagnostics$positive_definite) {
  x$sdr <- sdreport(obj, par.fixed = theta, getJointPrecision = FALSE)
  x$uncertainty_type <- "MAP inverse-Hessian / delta-method uncertainty"
}
saveRDS(x, fit_path)
write.csv(diagnostics, diagnostic_path, row.names = FALSE)
weak <- data.frame(parameter = make.unique(names(theta)),
                   loading = eig$vectors[, which.min(ev)])
weak <- weak[order(abs(weak$loading), decreasing = TRUE), ]
write.csv(head(weak, 50), weak_path, row.names = FALSE)
print(diagnostics)
cat("Largest loadings in weakest Hessian direction:\n")
print(head(weak, 15), row.names = FALSE)
