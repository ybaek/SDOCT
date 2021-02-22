library(Rcpp)
tf <- tempfile(fileext = "Module.cpp")
exposeClass("constrained_nnet",
            constructors = list(c("Eigen::VectorXd", "Eigen::VectorXd",
                                  "Eigen::MatrixXd", "int")),
            fields = c("y", "x", "K", "L"),
            methods = c("log_prob<>", "gradient<>"),
            rename = c(log_prob = "log_prob<>",
                       gradient = "gradient<>"),
            header = c("// [[Rcpp::depends(BH)]]",
                       "// [[Rcpp::depends(RcppEigen)]]",
                       "// [[Rcpp::depends(StanHeaders)]]",
                       "// [[Rcpp::plugins(cpp14)]]",
                       paste0("#include <", file.path(getwd(),
                              "constrained_nnet.hpp"), ">")),
      file = tf,
      Rfile = FALSE)
Sys.setenv(PKG_CXXFLAGS = paste0(Sys.getenv("PKG_CXXFLAGS"), " -I",
                                 system.file("include", "src", 
                                 package = "StanHeaders", mustWork = TRUE)))
sourceCpp(tf)

y <- rnorm(100)
x <- seq(1, 100, length.out = 10)
K <- matrix(0, 100, 100)
K[lower.tri(K)] <- dist(1:100)
K[upper.tri(K)] <- t(K)[upper.tri(K)]
K <- K[, x]
L <- 2
cn <- new(constrained_nnet, y = y, x = x, K = K, L = L)

xi <- rnorm(2)
vecW <- rnorm(2 * length(x))
alpha <- 1.0
beta <- 1.0
l_ou <- 1.0
tau <- 1.0
sigma <- 1.0
cn$log_prob(c(xi = xi, vecW = vecW, alpha = alpha, beta = beta, l_ou = l_ou, tau = tau, sigma = sigma))
cn$gradient(c(xi = xi, vecW = vecW, alpha = alpha, beta = beta, l_ou = l_ou, tau = tau, sigma = sigma))
