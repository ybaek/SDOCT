// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math.hpp>                         // Includes prim/ and rev/ (reverse-mode)
#include <RcppEigen.h>                           // After including stan/math

// [[Rcpp::plugins(cpp14)]]

using Eigen::Dynamic;
using Eigen::Matrix;

/* Objective function as a C++11 functor                 */
/* To be optimized with parameters passed in operator () */
/* Example: f(B) = tr(tanh(XB)'tanh(XB))                 */
/* Reshaping vector <-> matrix (Eigen::Map) is tricky    */
/* Analytically phrase f to operate ONLY with vector     */

struct foo {
  const Matrix<double, Dynamic, Dynamic> X_;

  foo(const Matrix<double, Dynamic, Dynamic>& X) : X_(X) {}

  template <typename T>
  T operator()(const Matrix<T, Dynamic, 1>& b) const {
    using namespace stan::math;
    Eigen::Index D = X_.cols();
    Eigen::Index N = b.size() / D;
    // Returns meaningless result when b.size() != N * D
    T sum = 0.;
    for (Eigen::Index i = 0; i < N; ++i) {
      /* For correct usage of Eigen::Ref, see the question */
      /* https://stackoverflow.com/questions/21132538/correct-usage-of-the-eigenref-class */
      const Eigen::Ref<const Matrix<T, Dynamic, 1>>& bsub = b.segment(i * D, D);
      sum += dot_self(tanh((X_ * bsub).eval()));
    }
    return sum;
  }
};


/* Functional returns the gradient of objective WRT b */

// [[Rcpp::export]]
auto g(Matrix<double, Dynamic, Dynamic> X, Matrix<double, Dynamic, 1> b) {
  foo f(X);
  double fb;
  Matrix<double, Dynamic, 1> grad_fb;
  stan::math::gradient(f, b, fb, grad_fb);
  return grad_fb;
}
