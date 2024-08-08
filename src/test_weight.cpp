#include <Rcpp.h>
#include <RcppEigen.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <iostream>
#include <algorithm>
#include <string>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
void readWeightMat(List model) {
    // Convert Weight_Mat from R matrix to Eigen::MatrixXd
    Eigen::MatrixXd Weight_Mat = Rcpp::as<Eigen::MatrixXd>(model["Weight_Mat"]);

    // Print the matrix
    Rcpp::Rcout << "Weight_Mat:" << std::endl << Weight_Mat << std::endl;
}

