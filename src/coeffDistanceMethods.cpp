#include <Rcpp.h>

using namespace Rcpp;


double logL1(const NumericMatrix coeffMatA, const NumericMatrix coeffMatB) {
  return(sum((log(1+abs(coeffMatA - coeffMatB)))));
}

double wLogL1(const NumericMatrix coeffMatA, const NumericMatrix coeffMatB) {
  return(sum((log(1+abs(coeffMatA - coeffMatB)))));
}

double logL1Complex(const ComplexVector coeffMatA, const ComplexVector coeffMatB) {
  return(sum(log(1+Mod(coeffMatA - coeffMatB))));
}

double binB(const LogicalMatrix coeffMatA, const LogicalMatrix coeffMatB) {
  return(sum((coeffMatA == true) & (coeffMatB == false)));
}

double binC(const LogicalMatrix coeffMatA, const LogicalMatrix coeffMatB) {
  return(sum((coeffMatA == false) & (coeffMatB == true)));
}

double binTipLab(const ComplexMatrix coeffMatA, const ComplexMatrix coeffMatB) {

  NumericVector w (coeffMatA.nrow()); // as<NumericVector>(seq(1,coeffMatA.nrow()));
  w = Range(1,coeffMatA.nrow());
  w.push_front(1);
  w = pow(w,-2);

  return(sum(w*log(1+Mod(coeffMatA - coeffMatB))));
}



// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::export]]
NumericMatrix coeffDistRcpp(const List coeffs, String method,  bool progressBar){

  double numCoeffs = coeffs.length();
  NumericMatrix distMat(numCoeffs, numCoeffs);
  Progress p((numCoeffs*numCoeffs-numCoeffs)/2, progressBar);


  if(method == "logL1"){
    for(int i = 0; i < numCoeffs; i++){
      for(int j = i+1; j < numCoeffs; j++){

        if (Progress::check_abort() ){
          return(-1);
        }

        distMat(i,j) = logL1(as<NumericMatrix>(coeffs[i]),as<NumericMatrix>(coeffs[j]));


        distMat(j,i) = distMat(i,j);

        p.increment();
      }
    }
  } else if(method == "wLogL1"){
    for(int i = 0; i < numCoeffs; i++){
      for(int j = i+1; j < numCoeffs; j++){

        if (Progress::check_abort() ){
          return(-1);
        }

        distMat(i,j) = wLogL1(as<NumericMatrix>(coeffs[i]),as<NumericMatrix>(coeffs[j]));


        distMat(j,i) = distMat(i,j);

        p.increment();
      }
    }


  } else if(method == "logL1Complex"){

    for(int i = 0; i < numCoeffs; i++){
      for(int j = i+1; j < numCoeffs; j++){

        if (Progress::check_abort() ){
          return(-1);
        }

        distMat(i,j) = logL1Complex(as<ComplexVector>(coeffs[i]),as<ComplexVector>(coeffs[j]));


        distMat(j,i) = distMat(i,j);

        p.increment();
      }
    }
  } else if(method == "b"){
    for(int i = 0; i < numCoeffs; i++){
      for(int j = i+1; j < numCoeffs; j++){

        if (Progress::check_abort() ){
          return(-1);
        }

        distMat(i,j) = binB(as<LogicalMatrix>(coeffs[i]),as<LogicalMatrix>(coeffs[j]));


        distMat(j,i) = distMat(i,j);

        p.increment();
      }
    }


  } else if(method == "c"){

    for(int i = 0; i < numCoeffs; i++){
      for(int j = i+1; j < numCoeffs; j++){

        if (Progress::check_abort() ){
          return(-1);
        }

        distMat(i,j) = binC(as<LogicalMatrix>(coeffs[i]),as<LogicalMatrix>(coeffs[j]));


        distMat(j,i) = distMat(i,j);

        p.increment();
      }
    }

  } else if(method == "tipLabel"){

    for(int i = 0; i < numCoeffs; i++){
      for(int j = i+1; j < numCoeffs; j++){

        if (Progress::check_abort() ){
          return(-1);
        }

        distMat(i,j) = binTipLab(as<ComplexMatrix>(coeffs[i]),as<ComplexMatrix>(coeffs[j]));


        distMat(j,i) = distMat(i,j);

        p.increment();
      }
    }

  } else {
    throw std::invalid_argument("invalid method");
  }


  return(distMat);
}
