#include <Rcpp.h>
using namespace Rcpp;
//' This function returns sum of minimum peptide/spectra count across AP-MS
//' runs.
//' @title GetPPN
//' 
//' @param mat A matrix of peptide/spectra count in the format of 
//' proteins (column) by AP-MS runs (row)
//' @return A matrix of proteins by proteins indicating sum of minimum 
//' peptide/spectra count across all AP-MS runs.
//' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
//' @references Guruharsha, K. G., et al. "A protein complex network of 
//' Drosophila melanogaster." Cell 147.3 (2011): 690-703.
//' \url{https://doi.org/10.1016/j.cell.2011.08.047}
//' @examples
//' GetPPN(mat)
//' 
// [[Rcpp::export]]
NumericMatrix GetPPN(NumericMatrix mat) {
  int nc = mat.ncol();
  int rstart = 0;
  int rend = mat.nrow();
  NumericMatrix rmat(nc, nc);
  for (int c1 = 0; c1 < nc; c1 ++){
    for(int c2 = 0; c2 < c1; c2++){
      int sMinXY = 0;
      for(int r = rstart; r < rend; r++){
        int minXY = 0;
        if(mat(r, c1) <= mat(r, c2)){
          minXY = mat(r, c1);
        }else{
          minXY = mat(r, c2);
        }
        sMinXY += minXY;
      }
      rmat(c1, c2) = sMinXY;
    }
  }
  return(rmat);
}
