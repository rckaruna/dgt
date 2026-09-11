# vanbelle_design.R -- reproduce the reviewer's subjects x raters probit
# design with a known population phi, and compare dgt_icc() (v0.5.0) with
# the sample intraclass kappa, the mean pairwise phi, and the value the
# pre-0.5.0 code returned.
#
# Generating model (probit threshold representation):
#   eta*_ij = alpha + S_i + R_j + e_ij,  S ~ N(0, var_sub), R ~ N(0, var_rat),
#   e ~ N(0, 1),  Y_ij = 1{eta*_ij > 0}.
# Population quantities:
#   tetrachoric rho = var_sub / (1 + var_sub + var_rat)     (= ICC_eta, probit)
#   P = Phi(alpha / sqrt(1 + var_sub + var_rat))
#   phi_abs = (P11 - P^2) / (P (1 - P)),  P11 = P(both latent > 0)
# phi_abs is the population intraclass kappa (Fleiss) and the absolute
# response-scale ICC of Definition 1.
#
# Run time: about one minute per brms fit on an M-series Mac with
# cmdstanr; set R_REPS to taste.

suppressPackageStartupMessages({
  library(brms); library(dgt); library(mvtnorm)
})

set.seed(20260910)
R_REPS  <- 5          # replicate datasets
n_sub   <- 50;  n_rat <- 5           # the reviewer's design
alpha   <- -0.9; var_sub <- 8.8; var_rat <- 1.15   # near her JAGS posterior means

# ---- population truth --------------------------------------------------
V    <- 1 + var_sub + var_rat
rho  <- var_sub / V
P    <- pnorm(alpha / sqrt(V))
P11  <- as.numeric(pmvnorm(lower = c(0, 0), upper = c(Inf, Inf),
                           mean = c(alpha, alpha),
                           sigma = matrix(c(V, var_sub, var_sub, V), 2, 2)))
phi_pop <- (P11 - P^2) / (P * (1 - P))
cat(sprintf("Population: P = %.3f, tetrachoric rho = %.3f, phi (intraclass kappa) = %.3f\n",
            P, rho, phi_pop))

# ---- sample estimators ---------------------------------------------------
fleiss_kappa_binary <- function(Y) {           # Y: subjects x raters, 0/1
  n <- ncol(Y); x <- rowSums(Y); pbar <- mean(Y)
  1 - sum(x * (n - x)) / (nrow(Y) * n * (n - 1) * pbar * (1 - pbar))
}
mean_pairwise_phi <- function(Y) {
  m <- suppressWarnings(cor(Y)); mean(m[upper.tri(m)], na.rm = TRUE)
}

simulate_design <- function(n_sub, n_rat, alpha, var_sub, var_rat) {
  S <- rnorm(n_sub, 0, sqrt(var_sub)); Rr <- rnorm(n_rat, 0, sqrt(var_rat))
  E <- matrix(rnorm(n_sub * n_rat), n_sub, n_rat)
  Y <- (alpha + outer(S, Rr, "+") + E > 0) * 1L
  long <- data.frame(y = as.vector(Y),
                     sub = factor(rep(seq_len(n_sub), n_rat)),
                     met = factor(rep(seq_len(n_rat), each = n_sub)))
  list(wide = Y, long = long)
}

fit_and_icc <- function(long, link) {
  fit <- brm(y ~ 1 + (1 | sub) + (1 | met), data = long,
             family = bernoulli(link = link),
             chains = 2, iter = 2000, warmup = 1000, cores = 2,
             backend = "cmdstanr", refresh = 0, silent = 2)
  res <- dgt_icc(fit, person_group = "sub", seed = 1)
  s <- res$summary
  # pre-0.5.0 value: logistic inverse link applied to the fitted parameters
  post <- posterior::as_draws_df(fit)
  old <- mean(sapply(seq_len(nrow(post)), function(i)
    dgt:::.binomial_icc_core(post$b_Intercept[i], post$sd_sub__Intercept[i],
                             post$sd_met__Intercept[i], 1,
                             invlink = stats::plogis)[["icc_abs"]]))
  c(icc_Y_abs = s$estimate[s$measure == "ICC_Y (response-scale), absolute"],
    icc_Y_rel = s$estimate[s$measure == "ICC_Y (response-scale), relative"],
    icc_eta   = s$estimate[s$measure == "ICC_eta (link-scale)"],
    icc_I     = s$estimate[s$measure == "ICC_I (information)"],
    pre_0_5_0 = old)
}

out <- vector("list", R_REPS)
for (r in seq_len(R_REPS)) {
  d <- simulate_design(n_sub, n_rat, alpha, var_sub, var_rat)
  kap <- fleiss_kappa_binary(d$wide); phi_s <- mean_pairwise_phi(d$wide)
  pr <- fit_and_icc(d$long, "probit"); lg <- fit_and_icc(d$long, "logit")
  out[[r]] <- data.frame(rep = r, kappa_sample = kap, phi_pairwise = phi_s,
                         t(setNames(pr, paste0("probit_", names(pr)))),
                         t(setNames(lg, paste0("logit_",  names(lg)))))
  cat(sprintf("rep %d: kappa = %.3f  probit ICC_Y = %.3f (pre-0.5.0 %.3f)  logit ICC_Y = %.3f  probit ICC_eta = %.3f\n",
              r, kap, pr["icc_Y_abs"], pr["pre_0_5_0"], lg["icc_Y_abs"], pr["icc_eta"]))
}
res <- do.call(rbind, out)
print(round(colMeans(res[, -1]), 3))
cat(sprintf("\nPopulation phi = %.3f; population tetrachoric = %.3f\n", phi_pop, rho))
cat("Expect: probit ICC_Y_abs and logit ICC_Y_abs both near the population phi\n",
    "and the sample kappa; probit ICC_eta near the tetrachoric rho;\n",
    "pre_0_5_0 markedly below phi (the reported bug).\n")
write.csv(res, "vanbelle_design_results.csv", row.names = FALSE)
