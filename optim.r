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
            constructors = list(c("Eigen::MatrixXd", "Eigen::VectorXd")),
            fields = c("D", "y"),
            methods = c("log_prob<>", "gradient<>"),
            rename = c(log_prob = "log_prob<>", gradient = "gradient<>"),
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
y <- rnorm(10)
D <- matrix(0, 10, 10)
D[lower.tri(D)] <- dist(1:10)
D[upper.tri(D)] <- t(D)[upper.tri(D)]
gp <- new(gp_mean, D = D, y = y)
gp$log_prob(c(tau = 1.0, sigma = 1.0, rho = 2.0))
round(gp$gradient(c(tau = 1.0, sigma = 1.0, rho = 2.0)), digits = 4)