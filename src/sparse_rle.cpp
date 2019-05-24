#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Convert R sparse matrix to RcppArmadillo sparse matrix
//' 
//' @param mat Sparse matrix object from R
//' @export
//' 
// [[Rcpp::export]]
arma::sp_mat convertSparse(S4 mat) {
  // obtain dim, i, p. x from S4 object
  IntegerVector dims = mat.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));
  
  int nrow = dims[0], ncol = dims[1];
  
  // use Armadillo sparse matrix constructor
  arma::sp_mat res(i, p, x, nrow, ncol);
  return res;
}

//' Calculate geometric means of sparse matrices
//' 
//' @param sparseMatrix Sparse matrix from R
//' @export
//' 
// [[Rcpp::export]]
NumericVector calcGeoMeansSparse(S4 sparseMatrix){
  // Conversion of sparse matrices is automatic
  arma::sp_mat counts = convertSparse(sparseMatrix);
  counts = counts.t();
  
  // Get number of genes - columns
  int ncols = counts.n_cols;
  
  // Populate vector with this many columns
  NumericVector geoMeans(ncols);
  
  // Iterate through each column to calculate values
  for (int i = 0; i < ncols; ++i){
    // Retrieve column
    arma::sp_mat col_counts = counts.col(i);
    
    // Get nonzeros
    arma::vec nonzeros = arma::nonzeros(col_counts);
    
    // Calculate geoMeans
    double geoMean = exp(mean(log(nonzeros)));
    geoMeans(i) = geoMean;
  }
  return(geoMeans);
}

//' Calculate size factors from sparse matrices
//' 
//' @param sparseMatrix A sparse matrix from R
//' @param geoMeans List of size factors
//' @export
//' 
// [[Rcpp::export]]
NumericVector calcNormFactorSparse(S4 sparseMatrix, NumericVector geoMeans){
  // Conversion of sparse matrices is automatic
  arma::sp_mat counts = convertSparse(sparseMatrix);
  
  // Get number of cells - columns
  int ncols = counts.n_cols;
  
  // Generate sizeFactors to output
  NumericVector sizeFactors(ncols);
  
  // Loop over columns to calculate sizeFactors
  for (int i = 0; i < ncols; i++){
    // Nonzero counts
    NumericVector values;
    
    // Get column
    arma::sp_mat cellCounts = counts.col(i);
    
    // Get start and end point
    for (arma::sp_mat::const_iterator j = cellCounts.begin(); j != cellCounts.end(); ++j){
      // Extract the corresponding geoMean
      double geoMean = geoMeans(j.row());
      
      if (*j > 0){
        values.push_back(*j / geoMean);
      }
    }
    
    sizeFactors(i) = median(values);
  }
  return(sizeFactors);
}
