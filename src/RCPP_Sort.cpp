#include <math.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
List InsertSorted(NumericVector x, double x_new) {

  x.insert(0,x_new);

  x.sort();

  NumericVector v = {x_new};

  List L = List::create(_["Vector"] =  x, _["Position"] = match(v,x));


  return(L);
}
