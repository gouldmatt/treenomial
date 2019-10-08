#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double logL1(const NumericMatrix coeffMatA, const NumericMatrix coeffMatB) {
  return(sum((log(1+abs(coeffMatA - coeffMatB)))));
}

// [[Rcpp::export]]
double logL1Complex(const ComplexVector coeffMatA, const ComplexVector coeffMatB) {
  return(sum(log(1+Mod(coeffMatA - coeffMatB))));
}

