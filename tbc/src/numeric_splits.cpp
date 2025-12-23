// src/numeric_splits.cpp
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace Rcpp;

// -------- optional RcppParallel (compile-time detect) -------------------------
#ifndef TBC_HAS_RCPP_PARALLEL
# if defined(__has_include)
#  if __has_include(<RcppParallel.h>)
#   include <RcppParallel.h>
#   define TBC_HAS_RCPP_PARALLEL 1
#  else
#   define TBC_HAS_RCPP_PARALLEL 0
#  endif
# else
#  define TBC_HAS_RCPP_PARALLEL 0
# endif
#endif
#if TBC_HAS_RCPP_PARALLEL
#include <RcppParallel.h>
using namespace RcppParallel;
#endif

// Simpson (= Gini) mit Priors
static inline double simpson_risk(const std::vector<double>& cnt,
                                  const std::vector<double>& prior) {
  const int K = (int)cnt.size();
  double m = 0.0, sumsq = 0.0;
  for (int k = 0; k < K; ++k) m += prior[k] * cnt[k];
  if (m <= 0.0) return 0.0;
  for (int k = 0; k < K; ++k) {
    const double pk = (prior[k] * cnt[k]) / m;
    sumsq += pk * pk;
  }
  return m * (1.0 - sumsq);
}

struct ColBest {
  double crit   = -std::numeric_limits<double>::infinity();
  double cut    = std::numeric_limits<double>::quiet_NaN();
  int    var    = 0;    // 1-based for R
  bool   naLeft = true; // preferred NA side (true = left)
};

// ---------- Kern (seriell, Rcpp) ---------------------------------------------
static inline ColBest eval_one_feature_serial(
    const NumericMatrix& X,
    const IntegerVector& y,              // Klassen 1..K
    const std::vector<double>& prior,
    const int j                          // 0-based
) {
  const int n = X.nrow();
  const int K = (int)prior.size();

  std::vector<double> v_known; v_known.reserve(n);
  std::vector<int>    cls_known; cls_known.reserve(n);
  std::vector<int>    naCounts(K, 0);

  for (int i = 0; i < n; ++i) {
    const double xij = X(i, j);
    const int lab = y[i];
    if (R_finite(xij)) {
      v_known.push_back(xij);
      cls_known.push_back((lab >= 1 && lab <= K) ? lab : 1);
    } else if (lab >= 1 && lab <= K) {
      naCounts[lab - 1] += 1;
    }
  }
  if ((int)v_known.size() < 2) return ColBest{};

  std::vector<int> ord(v_known.size());
  for (int i = 0; i < (int)ord.size(); ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return v_known[a] < v_known[b]; });

  const int T = (int)ord.size();
  std::vector<double> v_sorted(T);
  std::vector<int>    cls_sorted(T);
  for (int t = 0; t < T; ++t) {
    v_sorted[t]   = v_known[ord[t]];
    cls_sorted[t] = cls_known[ord[t]];
  }

  std::vector<double> parentCnt(K, 0.0);
  for (int i = 0; i < y.size(); ++i) {
    const int lab = y[i];
    if (lab >= 1 && lab <= K) parentCnt[lab - 1] += 1.0;
  }
  const double parentRisk = simpson_risk(parentCnt, prior);

  std::vector<double> leftKnown(K, 0.0), rightKnown(K, 0.0);
  for (int t = 0; t < T; ++t) {
    const int lab = cls_sorted[t];
    if (lab >= 1 && lab <= K) rightKnown[lab - 1] += 1.0;
  }

  ColBest best;
  for (int t = 0; t < T - 1; ++t) {
    const int lab = cls_sorted[t];
    leftKnown[lab - 1]  += 1.0;
    rightKnown[lab - 1] -= 1.0;

    if (!(v_sorted[t] < v_sorted[t + 1])) continue;

    double nL = 0.0, nR = 0.0;
    for (int k = 0; k < K; ++k) { nL += leftKnown[k]; nR += rightKnown[k]; }
    const bool route_left = (nL >= nR);

    std::vector<double> leftAll(leftKnown), rightAll(rightKnown);
    if (route_left) { for (int k = 0; k < K; ++k) leftAll[k]  += (double)naCounts[k]; }
    else            { for (int k = 0; k < K; ++k) rightAll[k] += (double)naCounts[k]; }

    double massL = 0.0, massR = 0.0;
    for (int k = 0; k < K; ++k) { massL += leftAll[k]; massR += rightAll[k]; }
    if (massL <= 0.0 || massR <= 0.0) continue;

    const double riskL = simpson_risk(leftAll,  prior);
    const double riskR = simpson_risk(rightAll, prior);
    const double crit  = parentRisk - (riskL + riskR);

    if (crit > best.crit) {
      best.crit = crit;
      best.cut  = 0.5 * (v_sorted[t] + v_sorted[t + 1]);
      best.var  = j + 1;
      best.naLeft = route_left;
    }
  }
  return best;
}

#if TBC_HAS_RCPP_PARALLEL
// ---------- Kern (parallel, RcppParallel) ------------------------------------
static inline ColBest eval_one_feature_parallel(
    const RcppParallel::RMatrix<double>& X,
    const RcppParallel::RVector<int>& y,
    const std::vector<double>& prior,
    const int j
) {
  const int n = X.nrow();
  const int K = (int)prior.size();

  std::vector<double> v_known; v_known.reserve(n);
  std::vector<int>    cls_known; cls_known.reserve(n);
  std::vector<int>    naCounts(K, 0);

  for (int i = 0; i < n; ++i) {
    const double xij = X(i, j);
    const int lab = y[i];
    if (R_finite(xij)) {
      v_known.push_back(xij);
      cls_known.push_back((lab >= 1 && lab <= K) ? lab : 1);
    } else if (lab >= 1 && lab <= K) {
      naCounts[lab - 1] += 1;
    }
  }
  if ((int)v_known.size() < 2) return ColBest{};

  std::vector<int> ord(v_known.size());
  for (int i = 0; i < (int)ord.size(); ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return v_known[a] < v_known[b]; });

  const int T = (int)ord.size();
  std::vector<double> v_sorted(T);
  std::vector<int>    cls_sorted(T);
  for (int t = 0; t < T; ++t) {
    v_sorted[t]   = v_known[ord[t]];
    cls_sorted[t] = cls_known[ord[t]];
  }

  std::vector<double> parentCnt(K, 0.0);
  for (std::size_t i = 0, N = (std::size_t)y.length(); i < N; ++i) {
    const int lab = y[(int)i];
    if (lab >= 1 && lab <= K) parentCnt[lab - 1] += 1.0;
  }
  const double parentRisk = simpson_risk(parentCnt, prior);

  std::vector<double> leftKnown(K, 0.0), rightKnown(K, 0.0);
  for (int t = 0; t < T; ++t) {
    const int lab = cls_sorted[t];
    if (lab >= 1 && lab <= K) rightKnown[lab - 1] += 1.0;
  }

  ColBest best;
  for (int t = 0; t < T - 1; ++t) {
    const int lab = cls_sorted[t];
    leftKnown[lab - 1]  += 1.0;
    rightKnown[lab - 1] -= 1.0;

    if (!(v_sorted[t] < v_sorted[t + 1])) continue;

    double nL = 0.0, nR = 0.0;
    for (int k = 0; k < K; ++k) { nL += leftKnown[k]; nR += rightKnown[k]; }
    const bool route_left = (nL >= nR);

    std::vector<double> leftAll(leftKnown), rightAll(rightKnown);
    if (route_left) { for (int k = 0; k < K; ++k) leftAll[k]  += (double)naCounts[k]; }
    else            { for (int k = 0; k < K; ++k) rightAll[k] += (double)naCounts[k]; }

    double massL = 0.0, massR = 0.0;
    for (int k = 0; k < K; ++k) { massL += leftAll[k]; massR += rightAll[k]; }
    if (massL <= 0.0 || massR <= 0.0) continue;

    const double riskL = simpson_risk(leftAll,  prior);
    const double riskR = simpson_risk(rightAll, prior);
    const double crit  = parentRisk - (riskL + riskR);

    if (crit > best.crit) {
      best.crit = crit;
      best.cut  = 0.5 * (v_sorted[t] + v_sorted[t + 1]);
      best.var  = j + 1;
      best.naLeft = route_left;
    }
  }
  return best;
}

struct BestReducer : public RcppParallel::Worker {
  const RcppParallel::RMatrix<double> X;
  const RcppParallel::RVector<int>    y;
  const std::vector<double>           prior;
  const std::vector<int>              order0; // 0-based
  ColBest localBest;

  BestReducer(const NumericMatrix& X_,
              const IntegerVector& y_,
              const std::vector<double>& prior_,
              const std::vector<int>& order_)
    : X(X_), y(y_), prior(prior_), order0(order_), localBest() {}

  BestReducer(BestReducer& other, RcppParallel::Split)
    : X(other.X), y(other.y), prior(other.prior), order0(other.order0), localBest() {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t idx = begin; idx < end; ++idx) {
      const int j = order0[idx];
      ColBest cb = eval_one_feature_parallel(X, y, prior, j);
      if (cb.crit > localBest.crit) localBest = cb;
    }
  }
  void join(const BestReducer& rhs) {
    if (rhs.localBest.crit > localBest.crit) localBest = rhs.localBest;
  }
};
#endif

// [[Rcpp::export]]
Rcpp::List tbc_best_numeric_splits(
    const Rcpp::NumericMatrix& X,
    const Rcpp::IntegerVector& y,
    const Rcpp::NumericVector& priorRatio,
    Rcpp::Nullable<Rcpp::IntegerVector> feature_order = R_NilValue
) {
  const std::size_t p = (std::size_t)X.ncol();
  if (p == 0) {
    return Rcpp::List::create(_["bestCrit"]=R_NilValue, _["bestVar"]=0, _["bestCut"]=R_NilValue);
  }

  const int K = priorRatio.size();
  std::vector<double> prior(K, 1.0);
  for (int k = 0; k < K; ++k) prior[k] = (double)priorRatio[k];

  // 0-based order
  std::vector<int> order0; order0.reserve(p);
  if (feature_order.isNotNull()) {
    IntegerVector fo(feature_order.get());
    if ((std::size_t)fo.size() == p) {
      for (std::size_t i = 0; i < p; ++i) {
        const int j1 = fo[(int)i];
        const int j0 = (j1 >= 1 && j1 <= (int)p) ? (j1 - 1) : (int)i;
        order0.push_back(j0);
      }
    }
  }
  if (order0.empty()) for (std::size_t j = 0; j < p; ++j) order0.push_back((int)j);

  ColBest best;

#if TBC_HAS_RCPP_PARALLEL
  {
    BestReducer reducer(X, y, prior, order0);
    parallelReduce((std::size_t)0, (std::size_t)order0.size(), reducer);
    best = reducer.localBest;
  }
#else
  for (std::size_t idx = 0; idx < order0.size(); ++idx) {
    const int j = order0[idx];
    ColBest cb = eval_one_feature_serial(X, y, prior, j);
    if (cb.crit > best.crit) best = cb;
  }
#endif

  return Rcpp::List::create(
    _["bestCrit"] = best.crit,
    _["bestVar"]  = best.var,
    _["bestCut"]  = best.cut,
    _["naLeft"]   = best.naLeft
  );
}
