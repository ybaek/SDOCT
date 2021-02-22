#include <stan/model/model_header.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

class constrained_nnet {

public: 
  Eigen::VectorXd y; // Observed data
  Eigen::VectorXd x; // Set of knot locations
  Eigen::MatrixXd K; // Kernel matrix (obs. domain)
  size_t L; // Deepest layer dim.

  constrained_nnet(Eigen::VectorXd y, Eigen::VectorXd x, Eigen::MatrixXd K, int L) :
    y(y), x(x), K(K), L(L) {}

  // TODO: raise error after checking K and x dimensions
  // TODO 2: now make this a loop across multiple y's (but same x, K, L)

  template <bool propto__ = false, bool jacobian__ = false, typename T__ = double>
  T__ log_prob(std::vector<T__>& params_r__) const {
    using namespace stan::math;
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::Ref;

    T__ lp__(0.0);
    accumulator<T__> lp_accum__;

    // set up model parameters
    std::vector<int> params_i__;
    stan::io::reader<T__> in__(params_r__, params_i__);
 
    size_t L_mid = x.size();
    size_t LmL      = L_mid * L;
    auto xi   = in__.vector_constrain(L);
    auto vecW = in__.vector_constrain(LmL);
    T__ alpha;
    T__ beta;
    T__ l_ou;
    T__ tau;
    T__ sigma;
    alpha = in__.scalar_lb_constrain(0);
    beta  = in__.scalar_lb_constrain(0);
  	l_ou  = in__.scalar_lb_constrain(0);
    tau   = in__.scalar_lb_constrain(0);
    sigma = in__.scalar_lb_constrain(0);

    // 1. Forming the covariance matrix
    // TODO: this code will work only for scalar x's (in R1)
    std::vector<double> x_s(x.data(), x.data() + L_mid);
    auto C = gp_exponential_cov(x_s, tau, l_ou);
    auto Quad = quad_form_sym(C, transpose(K).eval());
    size_t K_rows = rows(K);
    for (size_t n = 0; n < K_rows; ++n) {
      Quad(n,n) += square(sigma);
    }

    // 2. Forming the mean vector
     T__ alpha_inv = 1.0 / alpha;
     Matrix<T__, Dynamic, 1> eta;
     eta.resize(L_mid, 1);
     for (size_t n = 0; n < L_mid; ++n) {
       const Ref<const Matrix<T__, 1, Dynamic>>& vecWsub = vecW.segment(n * L, L).transpose();
       eta(n) = tanh(alpha_inv * (vecWsub * xi).value());
     }
     Matrix<T__, Dynamic, 1> meanvec = beta * K * eta;

    // log-likelihood (should add priors)
    lp_accum__.add(lp__);
    lp_accum__.add(multi_normal_lpdf<propto__>(y, meanvec, Quad));
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
