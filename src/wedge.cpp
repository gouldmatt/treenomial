// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "wedge.h"

using namespace arma;

// convert the sparse matrix to coordinate form matrix of [row, col, val]
mat coordReal(SpMat<double> &x) {

  mat  ans(x.n_nonzero, 3);
  int k = 0;
  sp_mat::const_iterator i;

  for (i = x.begin(); i != x.end(); ++i){
    ans(k, 0) = i.row(); // Row position
    ans(k, 1) = i.col(); // Col position
    ans(k++,2) = *i;
  }


  return(ans);
}

// performs the polynomial multiplication of A,B with the extra +y term
SpMat<double> wedge(SpMat<double> &A, SpMat<double> &B){

  Mat<double> baseMat = coordReal(A);
  Mat<double> shiftsMat = coordReal(B);

  int baseRows = baseMat.n_rows;
  int shiftsRows = shiftsMat.n_rows;
  int rowRes = 0;
  int colRes = 0;
  double coeffRes = 0;

  mat resMat(std::ceil(((double)A.n_cols + (double)B.n_cols)/2), A.n_cols + B.n_cols-1,fill::zeros);

  // add the extra + y term
  resMat(1,0) = 1;
  // loop through all combinations of each coordinate list (nonzero elements of each matrix)
  for(int bRow = 0; bRow<baseRows; bRow++){

    for(int sRow = 0; sRow<shiftsRows; sRow++){
      rowRes = baseMat(bRow,0) + shiftsMat(sRow,0);
      colRes = baseMat(bRow,1) + shiftsMat(sRow,1);
      coeffRes = baseMat(bRow,2)*shiftsMat(sRow,2);
      resMat(rowRes,colRes) = resMat(rowRes,colRes) + coeffRes;

    }
  }

  return(sp_mat(resMat));
}

// sp_cx_mat version of coordReal
cx_mat coordComplex(sp_cx_mat &x) {
  cx_mat res(x.n_nonzero, 3,fill::zeros);

  sp_cx_mat::const_iterator it     = x.begin();
  sp_cx_mat::const_iterator it_end = x.end();
  //
  int k = 0;
  for(; it != it_end; ++it)
  {
    res.at(k, 0) = (it.row());
    res.at(k, 1) = (it.col());
    res.at(k, 2) = (*it);
    k++;
  }

  return(res);
}

sp_cx_mat wedgeTipLabel(sp_cx_mat &A, sp_cx_mat &B){

  cx_mat baseMat = coordComplex(A);
  cx_mat shiftsMat = coordComplex(B);

  int baseRows = baseMat.n_rows;
  int shiftsRows = shiftsMat.n_rows;
  int rowRes = 0;
  int colRes = 0;
  cx_double coeffRes = 0;

  sp_cx_mat resMat(A.n_rows + B.n_rows - 1, A.n_cols + B.n_cols - 1);

  // add the extra + y term
  resMat(0,0) = cx_double(1,1);
  // loop through all combinations of each coordinate list (nonzero elements of each matrix)
  for(int bRow = 0; bRow<baseRows; bRow++){

    for(int sRow = 0; sRow<shiftsRows; sRow++){

      rowRes = real(baseMat(bRow,0)) + real(shiftsMat(sRow,0));
      colRes = real(baseMat(bRow,1)) + real(shiftsMat(sRow,1));
      coeffRes = baseMat(bRow,2)*shiftsMat(sRow,2);
      resMat(rowRes,colRes) += coeffRes;

    }
  }

  return(resMat);
}

cx_rowvec wedgeConv(cx_rowvec &A, cx_rowvec &B, cx_double y){
  cx_rowvec res = conv(A, B);
  res[0] += y;
  return(res);
}
