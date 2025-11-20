
suppressPackageStartupMessages({
  library(gamlss)
  library(VGAM)
})

# Simulate Beta-Binomial Data
set.seed(123)
n <- 1000
bd <- 50
mu_true <- 0.3
theta_true <- 0.1 # ASPEN theta = 1/(alpha+beta)
# alpha+beta = 1/theta = 10
# alpha = mu * 10 = 3
# beta = (1-mu) * 10 = 7

alpha <- 3
beta <- 7
prob <- rbeta(n, alpha, beta)
y <- rbinom(n, bd, prob)

# Fit GAMLSS
df <- data.frame(y=y, bd=bd)
m <- gamlss(cbind(y, bd-y) ~ 1, family=BB, data=df, trace=FALSE)

mu_est <- plogis(coef(m, "mu"))
sigma_est <- exp(coef(m, "sigma"))

cat("True Mu:", mu_true, "\n")
cat("Est Mu:", mu_est, "\n")
cat("True Theta (ASPEN):", theta_true, "\n")
cat("Est Sigma (GAMLSS):", sigma_est, "\n")

# Check relationship
# If sigma = theta, they should be close.
# If sigma = rho = 1/(alpha+beta+1), then rho = 1/11 = 0.0909.
# If sigma = theta, then sigma = 0.1.

cat("Hypothesis 1: Sigma = Theta (1/(alpha+beta)) -> Expect 0.1\n")
cat("Hypothesis 2: Sigma = Rho (1/(alpha+beta+1)) -> Expect 0.0909\n")
