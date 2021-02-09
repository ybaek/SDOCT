library(Rcpp)
library(StanHeaders)
Sys.setenv(PKG_CXXFLAGS = StanHeaders:::CxxFlags(as_character = TRUE))
SH <- system.file("lib", .Platform$r_arch,
                  package = "StanHeaders", mustWork = TRUE)
Sys.setenv(PKG_LIBS = paste0(StanHeaders:::LdFlags(as_character = TRUE),
           " -L", shQuote(SH), " -lStanHeaders"))
Rcpp::sourceCpp("autodiff_test.cpp")

# Evaluate the gradient
X <- matrix(rnorm(8), 2)
b <- rnorm(4)
g(X, b)
