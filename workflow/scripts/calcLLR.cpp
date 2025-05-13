#include <Rcpp.h>

double predq(const Rcpp::NumericVector &qest, const double gc, const double gcscale) {
  return(qest[0]+gc*qest[1]+gc*gc*qest[2]);
}

double expcPerBin(const Rcpp::IntegerVector &bincounts, const int l, const int r, const Rcpp::NumericVector &gc, const Rcpp::NumericVector &mapnorm, const Rcpp::NumericVector &qest, const double gcscale) {
  double slr = std::accumulate(bincounts.begin()+l, bincounts.begin()+r, 0.0);
  double wlr = 0;
  for(int i=l; i != r; ++i) {
    wlr += predq(qest, gc[i], gcscale)*mapnorm[i];
  }
  double explr = slr/wlr;
  return(explr);
}

void calcp(Rcpp::NumericVector &vec, const Rcpp::IntegerVector &bincounts, const double expli, const int l, const int r) {
  //  vec.resize(r-l);
  int k=0;
  for(int i=l; i != r; ++i) {
    if(k < vec.length()) { // for testing
      vec[k++] = R::dpois(bincounts[i], expli, 1);
    }
  }
}

// [[Rcpp::export]]
double calcLLR(Rcpp::IntegerVector &bincounts, int i, Rcpp::IntegerVector w, const Rcpp::NumericVector &mapnorm, const Rcpp::NumericVector &gc, double gcscale, const Rcpp::NumericVector &qest) {
  int l = i-w[0];
  int r = i+w[w.length()-1];
  if(l < 1 || r < 1) {
    return 0;
  }
  double expli = expcPerBin(bincounts, l, i, gc, mapnorm, qest,  gcscale);
  double expir = expcPerBin(bincounts, i, r, gc, mapnorm, qest,  gcscale);
  double explr = expcPerBin(bincounts, l, r, gc, mapnorm, qest,  gcscale);

  Rcpp::NumericVector pli(i-l);
  Rcpp::NumericVector pir(r-i);
  Rcpp::NumericVector plrli(i-l);
  Rcpp::NumericVector plrir(r-i);

  calcp(pli, bincounts, expli, l, i);
  calcp(plrli, bincounts, explr, l, i);

  calcp(pir, bincounts, expir, i, r);
  calcp(plrir, bincounts, explr, i, r);

  double retLLR = .0;
  for(int j = 0; j != pli.length(); ++j) {
    retLLR += pli[j]-plrli[j];
  }
  for(int j = 0; j != pir.length(); ++j) {
    retLLR += pir[j]-plrir[j];
  }
  
  return(retLLR);
}
