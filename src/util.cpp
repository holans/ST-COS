#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector seq_int_ordered(double from, double to)
{
  if (from > to) {
    return Rcpp::IntegerVector(0);
  }

  int a = ceil(from);
  int b = floor(to);
  int n = b - a + 1;
  return Rcpp::seq_len(n) + a - 1;
}
