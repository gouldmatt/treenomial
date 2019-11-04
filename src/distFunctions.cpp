// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


using namespace arma;

double logDiff(const mat coeffMatA, const mat coeffMatB) {
  return(accu(log(1+abs(coeffMatA - coeffMatB))));
}

double wLogDiff(const mat coeffMatA, const mat coeffMatB) {

  vec w(coeffMatA.n_rows+1,fill::ones);
  w = regspace(1,(coeffMatA.n_rows));
  w = shift(w, 1);
  w = pow(w,-2);

  mat logDiff(log(1+abs(coeffMatA - coeffMatB)));

  return(accu(logDiff.each_col()%w));
}

double logDiffComplex(const cx_rowvec coeffMatA, const cx_rowvec coeffMatB) {
  return(accu(log(1+abs(coeffMatA - coeffMatB))));
}


double binB(const mat coeffMatA, const mat coeffMatB) {
  return(accu((coeffMatA == true) && (coeffMatB == false)));
}

double binC(const mat coeffMatA, const mat coeffMatB) {
  return(accu((coeffMatA == false) && (coeffMatB == true)));
}

double tipLab(const cx_mat coeffMatA, const cx_mat coeffMatB) {
  return(accu(log(1+abs(coeffMatA - coeffMatB))));
}

// [[Rcpp::export]]
std::vector<double> compareCoeffRcpp(Rcpp::List coeffsList, std::string method){
  int coeffsLength = coeffsList.length();

  if(method == "logDiff"){
    std::vector<mat> Y = Rcpp::as<std::vector<mat>>(coeffsList);
    std::vector<double> distVect(coeffsLength-1);

    RcppThread::parallelFor(1, coeffsLength, [&distVect,&Y] (unsigned int i) {
      distVect[i-1] = logDiff(Y[0],Y[i]);
    });

    return(distVect);

  } else if (method == "wLogDiff"){
    std::vector<mat> Y = Rcpp::as<std::vector<mat>>(coeffsList);
    std::vector<double> distVect(coeffsLength-1);

    RcppThread::parallelFor(1, coeffsLength, [&distVect,&Y] (unsigned int i) {
      distVect[i-1] = wLogDiff(Y[0],Y[i]);
    });

    return(distVect);

  } else if (method == "pa"){
    std::vector<mat> Y = Rcpp::as<std::vector<mat>>(coeffsList);
    std::vector<double> distVect(coeffsLength-1);

    RcppThread::parallelFor(1, coeffsLength, [&distVect,&Y] (unsigned int i) {
      distVect[i-1] = binB(Y[0],Y[i]);
    });

    return(distVect);

  } else if (method == "ap"){
    std::vector<mat> Y = Rcpp::as<std::vector<mat>>(coeffsList);
    std::vector<double> distVect(coeffsLength-1);

    RcppThread::parallelFor(1, coeffsLength, [&distVect,&Y] (unsigned int i) {
      distVect[i-1] = binC(Y[0],Y[i]);
    });

    return(distVect);

  } else if (method == "logDiffComplex"){
    std::vector<cx_rowvec> Y = Rcpp::as<std::vector<cx_rowvec>>(coeffsList);
    std::vector<double> distVect(coeffsLength-1);

    RcppThread::parallelFor(1, coeffsLength, [&distVect,&Y] (unsigned int i) {
      distVect[i-1] = logDiffComplex(Y[0],Y[i]);
    });

    return(distVect);

  } else if (method == "tipLab"){
    std::vector<cx_mat> Y = Rcpp::as<std::vector<cx_mat>>(coeffsList);
    std::vector<double> distVect(coeffsLength-1);

    RcppThread::parallelFor(1, coeffsLength, [&distVect,&Y] (unsigned int i) {
      distVect[i-1] = logDiffComplex(Y[0],Y[i]);
    });


    return(distVect);

  } else {
    throw std::invalid_argument("invalid method");
  }
}


// [[Rcpp::export]]
Rcpp::NumericMatrix coeffDistRcpp(Rcpp::List coeffsList, std::string method){

  if(method == "logDiff"){
    std::vector<mat> coeffs = Rcpp::as<std::vector<mat>>(coeffsList);
    int numCoeffs = coeffs.size();
    mat distMat(numCoeffs, numCoeffs, fill::zeros);

    RcppThread::parallelFor(0, numCoeffs, [&distMat,&numCoeffs,&coeffs] (unsigned int i) {
      for(int j = i+1; j < numCoeffs; j++){
        distMat(i,j) = logDiff(coeffs[i],coeffs[j]);
      }
    });

    distMat = distMat.t() + distMat;

    return(Rcpp::wrap(distMat));

  } else if (method == "wLogDiff"){
    std::vector<mat> coeffs = Rcpp::as<std::vector<mat>>(coeffsList);
    int numCoeffs = coeffs.size();
    mat distMat(numCoeffs, numCoeffs, fill::zeros);

    RcppThread::parallelFor(0, numCoeffs, [&distMat,&numCoeffs,&coeffs] (unsigned int i) {
      for(int j = i+1; j < numCoeffs; j++){
        distMat(i,j) = wLogDiff(coeffs[i],coeffs[j]);
      }
    });

    distMat = distMat.t() + distMat;

    return(Rcpp::wrap(distMat));

  } else if (method == "pa"){
    std::vector<mat> coeffs = Rcpp::as<std::vector<mat>>(coeffsList);
    int numCoeffs = coeffs.size();
    mat distMat(numCoeffs, numCoeffs, fill::zeros);

    RcppThread::parallelFor(0, numCoeffs, [&distMat,&numCoeffs,&coeffs] (unsigned int i) {
      for(int j = i+1; j < numCoeffs; j++){
        distMat(i,j) = binB(coeffs[i],coeffs[j]);
      }
    });

    distMat = distMat.t() + distMat;

    return(Rcpp::wrap(distMat));

  } else if (method == "ap"){
    std::vector<mat> coeffs = Rcpp::as<std::vector<mat>>(coeffsList);
    int numCoeffs = coeffs.size();
    mat distMat(numCoeffs, numCoeffs, fill::zeros);

    RcppThread::parallelFor(0, numCoeffs, [&distMat,&numCoeffs,&coeffs] (unsigned int i) {
      for(int j = i+1; j < numCoeffs; j++){
        distMat(i,j) = binC(coeffs[i],coeffs[j]);
      }
    });

    distMat = distMat.t() + distMat;

    return(Rcpp::wrap(distMat));

  } else if (method == "logDiffComplex"){
    std::vector<cx_rowvec> coeffs = Rcpp::as<std::vector<cx_rowvec>>(coeffsList);
    int numCoeffs = coeffs.size();
    mat distMat(numCoeffs, numCoeffs, fill::zeros);

    RcppThread::parallelFor(0, numCoeffs, [&distMat,&numCoeffs,&coeffs] (unsigned int i) {
      for(int j = i+1; j < numCoeffs; j++){
        distMat(i,j) = logDiffComplex(coeffs[i],coeffs[j]);
      }
    });

    distMat = distMat.t() + distMat;

    return(Rcpp::wrap(distMat));

  } else if (method == "tipLab"){
    std::vector<cx_mat> coeffs = Rcpp::as<std::vector<cx_mat>>(coeffsList);
    int numCoeffs = coeffs.size();
    mat distMat(numCoeffs, numCoeffs, fill::zeros);

    RcppThread::parallelFor(0, numCoeffs, [&distMat,&numCoeffs,&coeffs] (unsigned int i) {
      for(int j = i+1; j < numCoeffs; j++){
        distMat(i,j) = tipLab(coeffs[i],coeffs[j]);
      }
    });

    distMat = distMat.t() + distMat;

    return(Rcpp::wrap(distMat));

  } else {
    throw std::invalid_argument("invalid method");
  }
}