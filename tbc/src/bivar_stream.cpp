// src/bivar_stream.cpp
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

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
using namespace RcppParallel;
#endif

// -----------------------------------------------------------------------------
// Simpson/Gini risk with prior weighting
static inline double simpson_risk(const std::vector<double>& cnt,
                                  const std::vector<double>& prior) {
  const int K = static_cast<int>(cnt.size());
  double m = 0.0, sumsq = 0.0;
  for (int k = 0; k < K; ++k) m += prior[k] * cnt[k];
  if (m <= 0.0) return 0.0;
  for (int k = 0; k < K; ++k) {
    const double pk = (prior[k] * cnt[k]) / m;
    sumsq += pk * pk;
  }
  return m * (1.0 - sumsq);
}

struct SplitResult {
  double crit   = R_NegInf;
  double cut    = NA_REAL;
  bool   naLeft = true;
};

// best split on projected vector v (length n), majority NA routing
static SplitResult best_split_for_vector(const std::vector<double>& v,
                                         const int* y1,                // classes 1..K
                                         int n,
                                         const std::vector<double>& prior,
                                         double parentRisk) {
  const int K = static_cast<int>(prior.size());

  std::vector<double> vals; vals.reserve(n);
  std::vector<int>    lab;  lab.reserve(n);
  std::vector<double> naCounts(K, 0.0);

  for (int r = 0; r < n; ++r) {
    const double vr = v[r];
    if (R_finite(vr)) {
      vals.push_back(vr);
      int c = y1[r];
      lab.push_back((c >= 1 && c <= K) ? c : 1);
    } else {
      int c = y1[r] - 1;
      if (c >= 0 && c < K) naCounts[c] += 1.0;
    }
  }

  const int nKnown = (int)vals.size();
  if (nKnown < 2) return SplitResult{};

  std::vector<int> ord(nKnown);
  for (int i = 0; i < nKnown; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return vals[a] < vals[b]; });

  std::vector<double> sv(nKnown);
  std::vector<int>    sc(nKnown);
  for (int i = 0; i < nKnown; ++i) { sv[i] = vals[ord[i]]; sc[i] = lab[ord[i]]; }

  std::vector<double> leftKnown(K, 0.0), rightKnown(K, 0.0);
  for (int t = 0; t < nKnown; ++t) {
    const int cls = sc[t] - 1;
    if (cls >= 0 && cls < K) rightKnown[cls] += 1.0;
  }

  SplitResult best;
  std::vector<double> leftAll(K), rightAll(K);

  for (int t = 0; t < nKnown - 1; ++t) {
    const int cls = sc[t] - 1;
    if (cls >= 0 && cls < K) { leftKnown[cls] += 1.0; rightKnown[cls] -= 1.0; }

    if (!(sv[t] < sv[t + 1])) continue;

    double nL = 0.0, nR = 0.0;
    for (int j = 0; j < K; ++j) { nL += leftKnown[j]; nR += rightKnown[j]; }

    const bool route_left = (nL >= nR); // tie -> left
    for (int j = 0; j < K; ++j) {
      leftAll[j]  = leftKnown[j]  + (route_left ? naCounts[j] : 0.0);
      rightAll[j] = rightKnown[j] + (route_left ? 0.0 : naCounts[j]);
    }

    // einzige Bedingung: beide Kinder > 0 nach Routing
    double massL = 0.0, massR = 0.0;
    for (int j = 0; j < K; ++j) { massL += leftAll[j]; massR += rightAll[j]; }
    if (massL <= 0.0 || massR <= 0.0) continue;

    const double riskL = simpson_risk(leftAll,  prior);
    const double riskR = simpson_risk(rightAll, prior);
    const double gain  = parentRisk - (riskL + riskR);

    if (gain > best.crit) {
      best.crit   = gain;
      best.cut    = 0.5 * (sv[t] + sv[t + 1]);
      best.naLeft = route_left;
    }
  }

  return best;
}

static inline int remap_label(int v, const std::vector<int>& uniq) {
  auto it = std::lower_bound(uniq.begin(), uniq.end(), v);
  if (it == uniq.end() || *it != v) return -1;
  return static_cast<int>(std::distance(uniq.begin(), it));
}

// ------------------------- SERIAL PATH ---------------------------------------
#if !TBC_HAS_RCPP_PARALLEL
static void bivar_serial(const NumericMatrix& X,
                         const IntegerVector& y1,         // 1..K
                         const std::vector<double>& prior,
                         double parentRisk,
                         const NumericVector& coefs,
                         double& outCrit, double& outCut,
                         int& outI, int& outJ, double& outCoef, bool& outNaLeft) {
  const int n = X.nrow(), p = X.ncol();
  std::vector<double> v(n);

  outCrit = R_NegInf; outCut = NA_REAL; outI = 0; outJ = 0; outCoef = NA_REAL; outNaLeft = true;

  for (int i = 0; i < p - 1; ++i) {
    for (int j = i + 1; j < p; ++j) {
      for (R_xlen_t c = 0; c < coefs.size(); ++c) {
        const double a = coefs[c];
        for (int r = 0; r < n; ++r) {
          const double xi = X(r, i), xj = X(r, j);
          v[r] = xi + a * xj; // NA/Inf propagate
        }
        SplitResult res = best_split_for_vector(v, INTEGER(y1), n, prior, parentRisk);
        if (res.crit > outCrit) {
          outCrit = res.crit; outCut = res.cut;
          outI = i + 1; outJ = j + 1; outCoef = a; outNaLeft = res.naLeft;
        }
      }
    }
  }
}
#endif

#if TBC_HAS_RCPP_PARALLEL
// ------------------------- PARALLEL PATH -------------------------------------
struct PairWorker : public RcppParallel::Worker {
  const RMatrix<double> X;
  const int*            y1;     // 1..K
  const int             n, p;
  const std::vector<double>& prior;
  const double          parentRisk;
  const std::vector<double> coefs;

  // per-pair outputs
  RVector<double> bestCrit;
  RVector<double> bestCut;
  RVector<int>    bestI;
  RVector<int>    bestJ;
  RVector<double> bestCoef;
  RVector<int>    bestNaLeft;

  PairWorker(const NumericMatrix& X_,
             const IntegerVector& y1_,
             const std::vector<double>& prior_,
             double parentRisk_,
             const NumericVector& coefs_,
             NumericVector& oCrit,
             NumericVector& oCut,
             IntegerVector& oI,
             IntegerVector& oJ,
             NumericVector& oCoef,
             IntegerVector& oNaLeft)
    : X(X_), y1(INTEGER(y1_)), n(X_.nrow()), p(X_.ncol()),
      prior(prior_), parentRisk(parentRisk_),
      coefs(as<std::vector<double>>(coefs_)),
      bestCrit(oCrit), bestCut(oCut), bestI(oI), bestJ(oJ), bestCoef(oCoef), bestNaLeft(oNaLeft) {}

  inline void pid_to_pair(int pid, int& i, int& j) const {
    int rem = pid;
    for (i = 0; i < p - 1; ++i) {
      const int cnt = p - i - 1;
      if (rem < cnt) { j = i + 1 + rem; return; }
      rem -= cnt;
    }
    i = 0; j = 1;
  }

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> v(n);
    for (std::size_t pid = begin; pid < end; ++pid) {
      int i, j;
      pid_to_pair(static_cast<int>(pid), i, j);

      double locCrit = R_NegInf, locCut = NA_REAL, locCoef = NA_REAL;
      int    locNaL  = 1;

      for (std::size_t c = 0; c < coefs.size(); ++c) {
        const double a = coefs[c];
        for (int r = 0; r < n; ++r) {
          const double xi = X(r, i), xj = X(r, j);
          v[r] = xi + a * xj;
        }
        SplitResult res = best_split_for_vector(v, y1, n, prior, parentRisk);
        if (res.crit > locCrit) {
          locCrit = res.crit; locCut = res.cut; locCoef = a; locNaL = res.naLeft ? 1 : 0;
        }
      }

      bestCrit[pid]   = locCrit;
      bestCut[pid]    = locCut;
      bestI[pid]      = i + 1;
      bestJ[pid]      = j + 1;
      bestCoef[pid]   = locCoef;
      bestNaLeft[pid] = locNaL;
    }
  }
};
#endif

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List tbc_best_bivariate_split(Rcpp::NumericMatrix X,
                                    Rcpp::IntegerVector y,
                                    Rcpp::NumericVector priorRatio,
                                    Rcpp::NumericVector coefs) {
  const int n = X.nrow();
  const int p = X.ncol();
  if (n <= 0 || p < 2) {
    return Rcpp::List::create(
      _["bestCrit"] = R_NegInf, _["bestCut"] = R_NaReal,
      _["i"] = 0, _["j"] = 0, _["coef"] = R_NaReal, _["naLeft"] = true
    );
  }

  // unique sorted labels -> map to 0..K-1, then to 1..K
  std::vector<int> uniq; uniq.reserve(n);
  for (int r = 0; r < n; ++r) uniq.push_back(y[r]);
  std::sort(uniq.begin(), uniq.end());
  uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
  const int K = static_cast<int>(uniq.size());

  IntegerVector y1(n);
  for (int r = 0; r < n; ++r) {
    int cls0 = -1;
    auto it = std::lower_bound(uniq.begin(), uniq.end(), y[r]);
    if (it != uniq.end() && *it == y[r]) cls0 = (int)std::distance(uniq.begin(), it);
    y1[r] = (cls0 >= 0 ? cls0 + 1 : 1);
  }

  std::vector<double> prior(K, 1.0);
  for (int k = 0; k < K && k < priorRatio.size(); ++k) prior[k] = priorRatio[k];

  // parent risk over ALL rows in node
  std::vector<double> parentCnt(K, 0.0);
  for (int r = 0; r < n; ++r) {
    const int cls = y1[r] - 1;
    if (cls >= 0) parentCnt[cls] += 1.0;
  }
  const double parentRisk = simpson_risk(parentCnt, prior);

#if TBC_HAS_RCPP_PARALLEL
  {
    const int npairs = p * (p - 1) / 2;
    NumericVector bestCrit(npairs, R_NegInf);
    NumericVector bestCut(npairs,  NA_REAL);
    IntegerVector bestI(npairs), bestJ(npairs), bestNaLeft(npairs);
    NumericVector bestCoef(npairs);

    PairWorker worker(X, y1, prior, parentRisk, coefs,
                      bestCrit, bestCut, bestI, bestJ, bestCoef, bestNaLeft);

    parallelFor(0, npairs, worker);

    // reduce to global best
    double gCrit = R_NegInf, gCut = NA_REAL, gCoef = NA_REAL;
    int gI = 0, gJ = 0; int gNaL = 1;
    for (int pid = 0; pid < npairs; ++pid) {
      const double c = bestCrit[pid];
      if (NumericVector::is_na(c) || !(c > gCrit)) continue;
      gCrit = c; gCut = bestCut[pid]; gI = bestI[pid]; gJ = bestJ[pid];
      gCoef = bestCoef[pid]; gNaL = bestNaLeft[pid];
    }

    return Rcpp::List::create(
      _["bestCrit"] = gCrit, _["bestCut"]  = gCut,
      _["i"] = gI, _["j"] = gJ, _["coef"] = gCoef,
      _["naLeft"] = static_cast<bool>(gNaL != 0)
    );
  }
#else
  // serial fallback
  double outCrit, outCut, outCoef; int outI, outJ; bool outNaLeft;
  bivar_serial(X, y1, prior, parentRisk, coefs,
               outCrit, outCut, outI, outJ, outCoef, outNaLeft);

  return Rcpp::List::create(
    _["bestCrit"] = outCrit, _["bestCut"]  = outCut,
    _["i"] = outI, _["j"] = outJ, _["coef"] = outCoef,
    _["naLeft"] = outNaLeft
  );
#endif
}
