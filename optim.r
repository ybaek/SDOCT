library(Rcpp)
# library(StanHeaders)
# Sys.setenv(PKG_CXXFLAGS = StanHeaders:::CxxFlags(as_character = TRUE))
# SH <- system.file("lib", .Platform$r_arch,
#                   package = "StanHeaders", mustWork = TRUE)
# Sys.setenv(PKG_LIBS = paste0(StanHeaders:::LdFlags(as_character = TRUE),
#            " -L", shQuote(SH), " -lStanHeaders"))
# Rcpp::sourceCpp("autodiff_test.cpp")
# X <- matrix(rnorm(8), 2)
# b <- rnorm(4)
# g(X, b)

tf <- tempfile(fileext = "Module.cpp")
exposeClass("gp_mean",
            constructors = list(c("Eigen::VectorXd", "Eigen::VectorXd")),
            fields = c("x", "y"),
            methods = c("log_prob<>", "gradient<>"),
            rename = c(log_prob = "log_prob<>", "gradient<>"),
            header = c("// [[Rcpp::depends(BH)]]",
                       "// [[Rcpp::depends(RcppEigen)]]",
                       "// [[Rcpp::depends(RcppParallel)]",
                       "// [[Rcpp::depends(StanHeaders)]]",
                       "// [[Rcpp::plugins(cpp14)]]",
                       paste0("#include <", file.path(getwd(), "gp_mean.hpp"), ">")),
      file = tf,
      Rfile = FALSE)
Sys.setenv(PKG_CXXFLAGS = paste0(Sys.getenv("PKG_CXXFLAGS"), " -I",
                                 system.file("include", "src", 
                                             package = "StanHeaders", mustWork = TRUE)))
sourceCpp(tf)
y <- rnorm(2)
x <- 1:2
D <- matrix(0, 2, 2)
D[lower.tri(D)] <- dist(x)
D[upper.tri(D)] <- t(D)[upper.tri(D)]

K <- exp(-D/2.) * 1.^2 + diag(.5^2, 2)
-log(2*pi)-.5*log(det(K))-.5*c(t(y) %*% solve(K) %*% y)

gp <- new(gp_mean, x = x, y = y)
gp$log_prob(c(tau = 0., rho = log(2.), sigma = log(.5)))
round(gp$gradient(c(tau = 0., rho = log(2.), sigma = log(.5))), digits = 4)
