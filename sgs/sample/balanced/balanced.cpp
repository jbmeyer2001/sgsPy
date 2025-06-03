/******************************************************************************
 *
 * Project: sgs
 * Purpose: Integrate BalancedSampling package into sgs
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 * Adapted from code authored by Wilmer Prentius in the following files:
 * https://github.com/envisim/BalancedSampling/blob/2.1.1/src/cube.cc
 * https://github.com/envisim/BalancedSampling/blob/2.1.1/src/cube_stratified.cc
 * https://github.com/envisim/BalancedSampling/blob/2.1.1/src/hlpm2.cc
 * License: GPL (>=2)
 *
 ******************************************************************************/


//#include "CubeClass.h"
//#include "CubeStratifiedClass.h"
//#include "IndexListClass.h"
//#include "KDTreeClass.h"
//#include "LpmClass.h"


// [[Rcpp::export(.lcube_cpp)]]
Rcpp::IntegerVector lcube_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &xbal,
  Rcpp::NumericMatrix &xspread,
  size_t treeBucketSize,
  int treeMethod,
  double eps
) {
  size_t N = xbal.nrow();
  size_t pbal = xbal.ncol();
  size_t pspread = xspread.nrow();

  if (N != (size_t)xspread.ncol())
    throw std::invalid_argument("xbal and xspread does not match");
  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob and x does not match");

  Cube cube(
    REAL(prob),
    REAL(xbal),
    N,
    pbal,
    eps,
    REAL(xspread),
    pspread,
    treeBucketSize,
    treeMethod
  );

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample.begin(), cube.sample.end());

  return sample;
}


// [[Rcpp::export(.lcube_stratified_cpp)]]
Rcpp::IntegerVector lcube_stratified_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &xbalance,
  Rcpp::NumericMatrix &xspread,
  Rcpp::IntegerVector &strata,
  const size_t bucketSize,
  const int method,
  const double eps
) {
  size_t N = xbalance.nrow();
  size_t p = xbalance.ncol();
  size_t pxs = xspread.nrow();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob and x does not match");
  if (N != (size_t)strata.length())
    throw std::range_error("strata and x does not match");
  if (N != (size_t)xspread.ncol())
    throw std::range_error("xspread and xbal does not match");

  CubeStratified cube(
    INTEGER(strata),
    REAL(prob),
    REAL(xbalance),
    N,
    p,
    eps,
    REAL(xspread),
    pxs,
    bucketSize,
    method
  );

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample_.begin(), cube.sample_.end());

  return sample;
}


// [[Rcpp::export(.hlpm2_cpp)]]
Rcpp::IntegerMatrix hlpm2_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &sizes,
  size_t treeBucketSize,
  int treeMethod,
  double eps
) {
  size_t sn = sizes.length();
  size_t N = x.ncol();
  size_t p = x.nrow();

  if (N != (size_t)prob.length())
    throw std::invalid_argument("prob an x does not match");

  Lpm lpm(
    LpmMethod::lpm2,
    REAL(prob),
    REAL(x),
    N,
    p,
    eps,
    treeBucketSize,
    treeMethod
  );

  // Make a copy of the tree for reuse later
  KDTree* orgTree = lpm.tree->Copy();
  IndexList* orgIdx = lpm.idx;

  // Run the algorithm to get a base sample
  lpm.Run();

  orgIdx->Fill();
  delete lpm.tree;

  // Set the return matrix
  size_t orgSampleSize = lpm.sample.size();
  Rcpp::IntegerMatrix sample(orgSampleSize, 2);

  for (size_t i = 0, j = 0; i < N; i++) {
    if (j < orgSampleSize && i == lpm.sample[j] - 1) {
      sample(j, 0) = lpm.sample[j];
      sample(j, 1) = sn;
      j += 1;
    } else {
      orgTree->RemoveUnit(i);
      orgIdx->Erase(i);
    }
  }

  size_t remainingSize = orgSampleSize;

  for (size_t i = 0; i < sn - 1; i++) {
    double subprob = (double)sizes[i] / (double)remainingSize;
    lpm.tree = orgTree->Copy();
    lpm.idx = orgIdx->CopyLen();
    lpm.sample.resize(0);

    for (size_t j = 0; j < orgIdx->Length(); j++) {
      lpm.probabilities[orgIdx->Get(j)] = subprob;
    }

    lpm.Run();

    // Remove all selected unit from orgTree and orgIdx, and set their subsample
    for (size_t j = 0, k = 0; j < orgSampleSize && k < lpm.sample.size(); j++) {
      if ((size_t)sample(j, 0) != lpm.sample[k])
        continue;

      size_t id = lpm.sample[k] - 1;
      orgTree->RemoveUnit(id);
      orgIdx->Erase(id);
      sample(j, 1) = i + 1;
      k += 1;
    }

    remainingSize -= lpm.sample.size();
    delete lpm.tree;
    delete lpm.idx;
  }

  lpm.tree = nullptr;
  lpm.idx = nullptr;
  delete orgTree;
  delete orgIdx;

  return sample;
}

