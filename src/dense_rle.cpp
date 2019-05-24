#include <Rcpp.h>
using namespace Rcpp;

// Calculate geometric means
double calcGeoMeanVal(NumericVector inputCounts){
  // Collect counts here
  NumericVector raw_counts;
  double sizeFactor;
  raw_counts = inputCounts[inputCounts > 0];
  sizeFactor = exp((mean(log(raw_counts))));
  return sizeFactor;
}

//' Calculate geometric means for dense matrices
//'
//' @param x A NumericMatrix
//' @export
//' 
// [[Rcpp::export]]
NumericVector calcGeoMeansDense(NumericMatrix x){
  int nrow = x.nrow();
  NumericVector out(nrow);
  
  for (int i = 0; i < nrow; i++){
    double sizeFactor = calcGeoMeanVal(x(i, _));
    out(i) = sizeFactor;
  }
  return(out);
}

// Calculate sizeFactors from dense matrix.
double calcNormFactor(NumericVector x, 
                      NumericVector y){
  
  // Filter zero values out of matrix
  NumericVector tempList;
  
  // Loop through rows, collect nonzero indices
  for (int i = 0; i < x.size(); i++){
    if (x(i) > 0){
      double val = x(i) / y(i);
      tempList.push_back(val);
    }
  }
  
  return median(tempList);
}

//' Calculate size factors for dense matrices
//'
//'
//' @param x matrix
//' @param y geometric means
//' @export
// [[Rcpp::export]]
NumericVector calcNormFactorDense(NumericMatrix x, NumericVector y){
  int ncol = x.ncol();
  NumericVector out(ncol);
  
  for (int i = 0; i < ncol; i++){
    // Turn row into a vector
    NumericVector counts = x(_, i);
    
    // Create a matrix from indices, counts, geoMeans
    NumericVector geoMeans = y;
    
    double normFactor = calcNormFactor(counts, geoMeans);
    out(i) = normFactor;
  }
  return(out);
}

//' Apply size factors to dense matrix
//' 
//' @param x A dense matrix
//' @param y Size factors
//' @export
//' 
// [[Rcpp::export]]
NumericMatrix applyRLE(NumericMatrix x, NumericVector y){
  // Number of columns
  int ncol = x.ncol();
  
  // Output matrix
  NumericMatrix out(x.size());
  
  // Iterate over columns, divide all values by sizeFactor
  for (int i = 0; i < ncol; i++){
    out(_, i) = x(_, i) / y(i);
  }
  return out;
}

