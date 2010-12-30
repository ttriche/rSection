library(inline)
library(Rcpp)

inc <- '
#include <iostream>
#include <armadillo>
#include <cppbugs/cppbugs.hpp>

using namespace arma;
using namespace cppbugs;

class TestModel: public MCModel {
public:
  const mat& y; // given
  const mat& X; // given

  Normal<vec> b;
  Uniform<double> tau_y;
  Deterministic<mat> y_hat;
  Normal<mat> likelihood;
  Deterministic<double> rsq;

  TestModel(const mat& y_,const mat& X_): y(y_), X(X_),
                                          b(randn<vec>(X_.n_cols)),
                                          tau_y(1),
                                          y_hat(X*b.value),
                                          likelihood(y_,true),
                                          rsq(0)
  {
    add(b);
    add(tau_y);
    add(y_hat);
    add(likelihood);
    add(rsq);
  }

  void update() {
    y_hat.value = X*b.value;
    rsq.value = as_scalar(1 - var(y - y_hat.value) / var(y));
  }
  double logp() const {
    return b.logp(0.0, 0.0001) + tau_y.logp(0,100) +
likelihood.logp(y_hat.value,tau_y.value);
  }
};
'

src <- '
mat X(REAL(XR),100,2);
mat y(REAL(yr),100,1);

TestModel m(y,X);
int iterations = 1e5;
m.sample(iterations, 1e4, 10);
return Rcpp::List::create(Rcpp::Named("b", m.b.mean()),
Rcpp::Named("ar", m.acceptance_ratio()));
'
fun <- cxxfunction(signature(XR="numeric", yr="numeric"), body=src,
include=inc, plugin="Rcpp")
NR <- 100
NC <- 2
x <- matrix(rnorm(NR*NC),NR,NC)
y <- rnorm(NR)
print(system.time(ans <- fun(x,y)))
print(ans)
