library(GIGrvg)

data <- readRDS("../data/fit_data.Rds")
attach(data)
source("./sampler_sketch.R")
source("./utilities.R")
max_iter <- 2000
# Initialize the working parameters...
N <- dim(Z_s)[1]; P <- dim(Z_s)[2]; Q <- dim(Y_s)[2]; J <- max(ids)
rho <- .99
l <- 15 / 768 * (2*pi) 
prec_y <- rho * Lap + diag(1-rho, Q)
prec_z <- solve(wendland_c2(distMat, l))
r_y <- chol(prec_y); r_z <- chol(prec_z)
beta <- matrix(0, P*Q, max_iter)
beta[,1] <- c(coef(lm(Y_s ~ Z_s-1))) # will increase in length dynamically
sig2inv <- numeric(max_iter); sig2inv[1] <- 1/.7
tau2inv <- numeric(max_iter); tau2inv[1] <- 10
Y_s[mis_inds] <- rnorm(dim(mis_inds)[1]) * 3 # Initialize the missing values

c0 <- 1/3 
c2 <- matrix(1, length(small_inds), max_iter)
c2[,1] <- rep(.25^2, length(small_inds))
beta_diags <- rep(100, length(beta)); beta_diags[small_inds] <- c2[,1] * 100

gamma <- array(0, dim = c(P, Q, J))
for (j in 1:J) gamma[,,j] <- t(backsolve(r_y, t(backsolve(r_z, matrix(rnorm(P*Q), P))))) / tau2inv[1]
gamma_mean_a <- matrix(0, N, Q) # array of Z %*% gamma's
beta_m <- matrix(beta[,1], P) # Matricized form of beta (at each iter)
for (j in 1:J) gamma_mean_a[id_list[[j]],] <- matrix(Z_s[id_list[[j]],], ncol=P) %*% gamma[,,j]
theta <- (Y_s - Z_s %*% beta_m - gamma_mean_a)
theta <- t(scale(t(theta), T, F)) # sum-to-zero constraint
lpd <- matrix(0, dim(Y_s)[1], max_iter)

# Let's stick with a fixed effect across patients
# A loop for computing block diagonal cholesky
beta_pprec_r_a <- matrix(0, P, P*Q) # container of choleskys
umean_a <- numeric(P*Q) # container of vectors needed for sampling sig2inv

for (iter in 2:max_iter) {
  sig_tr_total <- 0 # keeping sum of all traces
  beta_umean <- c(crossprod(Z_s, Y_s - gamma_mean_a - theta))
  beta_uz <- rnorm(P*Q)
  for (q in 1:Q) {
    cholInd <- ((q-1)*P+1):(q*P)
    beta_pprec_r_a[,cholInd] <- chol(diag(1/beta_diags[cholInd]) + crossprod(Z_s)*(1-rho))
    umean_a[cholInd] <- forwardsolve(t(beta_pprec_r_a[,cholInd]), beta_umean[cholInd])
    sig_tr_total <- sig_tr_total + (1-rho)^2*sum(umean_a[cholInd]^2)
    beta[cholInd,iter] <- (1-rho)*backsolve(beta_pprec_r_a[,cholInd], umean_a[cholInd])
    beta_uz[cholInd] <- backsolve(beta_pprec_r_a[,cholInd], beta_uz[cholInd])
  }
  sig2inv[iter] <- rgamma(1, .5+.5*N*Q, .5+.5*
                          ((1-rho)*sum(diag(crossprod(Y_s-gamma_mean_a-theta)))-sig_tr_total) )
  beta[,iter] <- beta[,iter] + beta_uz / sig2inv[iter]
  beta_m <- matrix(beta[,iter], P)
  
  for (j in 1:J) {
    Z_slice <- matrix(Z_s[id_list[[j]],], ncol=P)
    Y_slice <- matrix(Y_s[id_list[[j]],], ncol=Q)
    prec_gamma_r <- chol(tau2inv[iter-1]*prec_z + crossprod(Z_slice))
    gamma[,,j] <- t(backsolve(r_y, t(backsolve(prec_gamma_r, matrix(rnorm(P*Q), P))))) / (sig2inv[iter] * tau2inv[iter-1]) + 
      forbacksolve(prec_gamma_r, c(crossprod(Z_slice, Y_slice - Z_slice %*% beta_m)))
    gamma_mean_a[id_list[[j]],] <- Z_slice %*% gamma[,,j]
  }

  tau2_sse <- sum(apply(gamma, 3, function(g) sum(diag(prec_y %*% t(g) %*% prec_z %*% g))))*sig2inv[iter]
  tau2inv[iter] <- rgamma(1, .5+.5*J*P*Q, .5+.5*tau2_sse)
  
  c2[,iter] <- sapply(beta[small_inds,iter], function(b) rgig(1, 0, sum(b^2)*sig2inv[iter], c0^-2))
  beta_diags[small_inds] <- c2[,iter] * 100
  
  theta <- t(backsolve(r_y, matrix(rnorm(N*Q), Q))) / sig2inv[iter] + (1-rho) * 
      t(forbacksolve(r_y, t(Y_s - Z_s %*% beta_m - gamma_mean_a)))
  theta <- t(scale(t(theta), T, F)) # sum-to-zero constraint
  
  Y_s <- impute_car(Y_s, Z_s %*% beta_m + gamma_mean_a + theta, 
                    sig2inv[iter], rho, mis_inds, mis_adjs)
  
  lpd[,iter] <- rowSums(dnorm(Y_s - Z_s %*% beta_m - gamma_mean_a - theta, 
                              sd = (sig2inv[iter]*(1-rho))^-.5, log = T))
}
