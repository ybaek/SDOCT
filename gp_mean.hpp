#include <stan/model/model_header.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

class gp_mean {

public: // these would ordinarily be private in the C++ code generated by Stan
  Eigen::MatrixXd D;
  Eigen::VectorXd y;

  gp_mean(Eigen::MatrixXd D, Eigen::VectorXd y) :
    D(D), y(y) {}

  template <bool propto__ = false, bool jacobian__ = false, typename T__ = double>
  // propto__ is usually true but causes log_prob() to return 0 when called from R
  // jacobian__ is usually true for MCMC but typically is false for optimization
  T__ log_prob(std::vector<T__>& params_r__) const {
    using namespace stan::math;
    T__ lp__(0.0);
    accumulator<T__> lp_accum__;

    // set up model parameters
    std::vector<int> params_i__;
    stan::io::reader<T__> in__(params_r__, params_i__);
    T__ tau;
    T__ rho;
    T__ sigma;
    if (jacobian__) {
        tau   = in__.scalar_lb_constrain(0, lp__);
        sigma = in__.scalar_lb_constrain(0, lp__);
        rho   = in__.scalar_lb_constrain(0, lp__);
    }
    else {
        tau   = in__.scalar_lb_constrain(0);
        sigma = in__.scalar_lb_constrain(0);
        rho   = in__.scalar_lb_constrain(0);
    }
    // FIXME: this operation yields the wrong log density!
    Eigen::Matrix<T__, Eigen::Dynamic, Eigen::Dynamic> K =
      tau * (exp(-D.array() / rho)).matrix() +
      sigma *
      Eigen::MatrixXd::Identity(y.size(), y.size());

    // log-likelihood (should add priors)
    lp_accum__.add(lp__);
    lp_accum__.add(multi_normal_lpdf<propto__>(y,
      Eigen::VectorXd::Zero(y.size()).eval(),
      K));
    return lp_accum__.sum();
  }

  template <bool propto__ = false, bool jacobian__ = false>
  std::vector<double> gradient(std::vector<double>& params_r__) const {
    // Calculate gradients using reverse-mode autodiff
    // although you could do them analytically in this case

    using std::vector;
    using stan::math::var;
    double lp;
    std::vector<double> gradient;
    try {
      vector<var> ad_params_r(params_r__.size());
      for (size_t i = 0; i < params_r__.size(); ++i) {
        var var_i(params_r__[i]);
        ad_params_r[i] = var_i;
      }
      var adLogProb
        = this->log_prob<propto__, jacobian__>(ad_params_r);
      lp = adLogProb.val();
      adLogProb.grad(ad_params_r, gradient);
    } catch (const std::exception &ex) {
      stan::math::recover_memory();
      throw;
    }
    stan::math::recover_memory();
    return gradient;
  }
};