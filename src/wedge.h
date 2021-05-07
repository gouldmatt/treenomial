#ifndef WEDGE_H
#define WEDGE_H

// convert the sparse matrix to coordinate form matrix of [row, col, val]
arma::mat coordReal(arma::SpMat<double> &x); 

// performs the polynomial multiplication of A,B with the extra +y term
arma::SpMat<double> wedge(arma::SpMat<double> &A, arma::SpMat<double> &B); 

// sp_cx_mat version of coordReal
arma::cx_mat coordComplex(arma::sp_cx_mat &x);

arma::sp_cx_mat wedgeTipLabel(arma::sp_cx_mat &A, arma::sp_cx_mat &B);

arma::cx_rowvec wedgeConv(arma::cx_rowvec &A, arma::cx_rowvec &B, arma::cx_double y); 

#endif  // WEDGE_H