// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
void wedgeFill(NumericMatrix baseMat, NumericMatrix shiftsMat, NumericMatrix resMat) {

  int baseRows = baseMat.nrow();
  int shiftsRows = shiftsMat.nrow();
  int rowRes = 0;
  int colRes = 0;
  long double coeffRes = 0;

  // add the extra + y term
  resMat(1,0) = 1;
  // loop through all combinations of each coordinate list (nonzero elements of each matrix)
  for(int bRow = 0; bRow<baseRows; bRow++){
    for(int sRow = 0; sRow<shiftsRows; sRow++){
        rowRes = baseMat(bRow,0) + shiftsMat(sRow,0) - 2;
        colRes = baseMat(bRow,1) + shiftsMat(sRow,1) - 2;
        coeffRes = baseMat(bRow,2)*shiftsMat(sRow,2);
        //std::cout << "brow|srow" << bRow << "|" << sRow << std::endl <<  rowRes << " " << colRes << " " <<  Result(rowRes,colRes) + coeffRes << std::endl;
        resMat(rowRes,colRes) = resMat(rowRes,colRes) + coeffRes;
      }
  }
}


// [[Rcpp::export]]
void wedgeFillComplex(ComplexMatrix baseMat, ComplexMatrix shiftsMat, ComplexMatrix resMat) {

  int baseRows = baseMat.nrow();
  int shiftsRows = shiftsMat.nrow();

  int rowRes; // = Rcomplex(0);
  int colRes; //  = Rcomplex(0);
  Rcomplex coeffRes;

  // add the extra + 1 + i term
  resMat(0,0).i = resMat(0,0).i + 1;
  resMat(0,0).r = resMat(0,0).r + 1;

  // loop through all combinations of each coordinate list (nonzero elements of each matrix)
  for(int bRow = 0; bRow<baseRows; bRow++){
    for(int sRow = 0; sRow<shiftsRows; sRow++){
      rowRes = baseMat(bRow,0).r + shiftsMat(sRow,0).r - 2;
      colRes = baseMat(bRow,1).r + shiftsMat(sRow,1).r - 2;
      coeffRes = baseMat(bRow,2)*shiftsMat(sRow,2);
      resMat(rowRes,colRes) = resMat(rowRes,colRes) + coeffRes;
    }
  }
}
