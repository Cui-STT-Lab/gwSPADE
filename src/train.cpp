#include <Rcpp.h>

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>

#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <string>

// Sampler
#include "sampler.h"

// keyATM models
#include "keyATM_base.h"

// Weighted LDA models
#include "LDA_weight.h"


// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Eigen;
using namespace Rcpp;
using namespace std;



//' Run the Collapsed Gibbs sampler for weighted LDA
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_LDA(List model, bool resume = false)
{
  LDAweight LDAweight_model(model);
  if (resume) {
    LDAweight_model.resume_fit();
  } else {
    LDAweight_model.fit();
  }
  model = LDAweight_model.return_model();
  return model;
}


