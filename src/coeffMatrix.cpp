// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "wedge.h"

using namespace arma;


// main function to convert a wedge order to a coefficient matrix
// the complex and tip label version use similiar strategy using the different wedge versions
mat coeffMatrixReal(std::vector<std::string> wedgeOrder){

  long unsigned int j = 0;
  int subTreeNum = 2;
  std::string op1;
  std::string op2;
  std::string subPattern;

  // this will store the unique coefficient matrices as the total matrix is built
  // intialize with the leaf and cherry matrices
  std::vector<SpMat<double>> subCoeffMats(2);

  SpMat<double> leaf(1,2);
  leaf(0,1) = 1;

  SpMat<double> cherry(2,3);
  cherry(1,0) = 1;
  cherry(0,2) = 1;

  subCoeffMats[0] = leaf;
  subCoeffMats[1] = cherry;

  // each valid string in the wedge order maps to a matrix
  std::map<std::string, int> subCoeffOrder;
  subCoeffOrder.insert(std::make_pair("0", 0));
  subCoeffOrder.insert(std::make_pair("001", 1));

  // loop until the entire wedgeorder is one element
  while(wedgeOrder.size() != 1){
    j = 2;
    // determine the wedge operands
    while(true){
      // RcppThread::checkUserInterrupt();
      if(wedgeOrder[j] == "1"){
        op1 = wedgeOrder[j-2];
        op2 = wedgeOrder[j-1];
        break;
      }
      j++;
    }
    // RcppThread::checkUserInterrupt();

    // insert the new wedge in the map
    subPattern = op1 + op2 + "1";
    subCoeffOrder.insert(std::make_pair(subPattern,subTreeNum));

    // access the two operand matrices and perform the wedge
    SpMat<double> op1Mat = subCoeffMats[subCoeffOrder.find(op1)->second];
    SpMat<double> op2Mat = subCoeffMats[subCoeffOrder.find(op2)->second];
    subCoeffMats.push_back(wedge(op1Mat,op2Mat));

    subTreeNum++;


    // go through and flag any repeats of the current wedge operation
    for(j = 0; j < wedgeOrder.size()-2; j++){

      // account for order swapped option of wedge operands
      if(wedgeOrder[j] == op1 && wedgeOrder[j+1] == op2 && wedgeOrder[j+2] == "1"){

        wedgeOrder[j] = subPattern;
        wedgeOrder[j+1] = " ";
        wedgeOrder[j+2] = " ";

      } else if(wedgeOrder[j] == op2 && wedgeOrder[j+1] == op1 && wedgeOrder[j+2] == "1"){

        wedgeOrder[j] = subPattern;
        wedgeOrder[j+1] = " ";
        wedgeOrder[j+2] = " ";

      }
    }

    // erase the repeated locations
    wedgeOrder.erase(std::remove_if(wedgeOrder.begin(),
                                    wedgeOrder.end(),
                                    [](std::string x){return x == " ";}),
                                    wedgeOrder.end());

  }

  return(mat(subCoeffMats.back()));
}


// [[Rcpp::export]]
arma::mat wedgeExport( const arma::mat& A,  const arma::mat& B){

 SpMat<double> a(A);
 SpMat<double> b(B);
 return(arma::mat(wedge(a,b)));

}



// [[Rcpp::export]]
Rcpp::List coeffMatListDefault(std::vector<std::vector<std::string>> wedgeOrders, arma::cx_double y, std::string tipLabA = " ", std::string tipLabB = " ", int nThreads = -1){
  int numCoeffs = wedgeOrders.size();

  size_t numThreads = std::thread::hardware_concurrency();
  if(nThreads != -1){
    numThreads = nThreads;
  }

  Rcpp::List output(numCoeffs);


  arma::field<arma::mat> coeffs(numCoeffs);

  RcppThread::parallelFor(0, numCoeffs, [&coeffs, &wedgeOrders] (unsigned int i) {
    coeffs[i] = coeffMatrixReal(wedgeOrders[i]);
  },numThreads,0);

  output = Rcpp::wrap(coeffs);


  return(output);
}



