// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "wedge.h"

using namespace arma;


cx_rowvec coeffMatrixComplex(std::vector<std::string> wedgeOrder, cx_double y){
  long unsigned int j = 0;
  int subTreeNum = 2;
  std::string op1;
  std::string op2;
  std::string subPattern;


  std::vector<cx_rowvec> subCoeffMats(2);

  cx_rowvec  leaf(2,fill::zeros);
  leaf[1] = cx_double(1,0);


  cx_rowvec cherry(3,fill::zeros);
  cherry[0] = y;
  cherry[2] = cx_double(1, 0);

  subCoeffMats[0] = leaf;
  subCoeffMats[1] = cherry;


  std::map<std::string, int> subCoeffOrder;
  subCoeffOrder.insert(std::make_pair("0", 0));
  subCoeffOrder.insert(std::make_pair("001", 1));


  while(wedgeOrder.size() != 1){

    j = 2;
    while(true){
      // RcppThread::checkUserInterrupt();
      if(wedgeOrder[j] == "1"){
        op1 = wedgeOrder[j-2];
        op2 = wedgeOrder[j-1];
        break;
      }
      j++;
    }


    subPattern = op1 + op2 + "1";

    subCoeffOrder.insert(std::make_pair(subPattern,subTreeNum));




    cx_rowvec  op1Mat = subCoeffMats[subCoeffOrder.find(op1)->second];
    cx_rowvec  op2Mat = subCoeffMats[subCoeffOrder.find(op2)->second];
    subCoeffMats.push_back(wedgeConv(op1Mat,op2Mat,y));



    subTreeNum++;


    for(j = 0; j < wedgeOrder.size()-2; j++){

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

    wedgeOrder.erase(std::remove_if(wedgeOrder.begin(),
                                    wedgeOrder.end(),
                                    [](std::string x){return x == " ";}),
                                    wedgeOrder.end());

  }

  return(subCoeffMats.back());
}

sp_cx_mat coeffMatrixTipLabel(std::vector<std::string> wedgeOrder, std::string tipLabA, std::string tipLabB){
  long unsigned int j = 0;
  int subTreeNum = 2;
  std::string op1;
  std::string op2;
  std::string subPattern;


  std::vector<sp_cx_mat> subCoeffMats(2);

  sp_cx_mat  leafA(1,2);
  leafA(0,1) = cx_double(1,0);

  sp_cx_mat  leafB(2,1);
  leafB(1,0) = cx_double(1,0);


  subCoeffMats[0] = leafA;
  subCoeffMats[1] = leafB;

  std::map<std::string, int> subCoeffOrder;
  // use the tip labels to define the two types of leafs
  subCoeffOrder.insert(std::make_pair(tipLabA, 0));
  subCoeffOrder.insert(std::make_pair(tipLabB, 1));


  while(wedgeOrder.size() != 1){

    j = 2;
    while(true){
      // RcppThread::checkUserInterrupt();
      if(wedgeOrder[j] == "1"){
        op1 = wedgeOrder[j-2];
        op2 = wedgeOrder[j-1];
        break;
      }
      j++;
    }


    subPattern = op1 + op2 + "1";

    subCoeffOrder.insert(std::make_pair(subPattern,subTreeNum));




    sp_cx_mat  op1Mat = subCoeffMats[subCoeffOrder.find(op1)->second];
    sp_cx_mat  op2Mat = subCoeffMats[subCoeffOrder.find(op2)->second];
    subCoeffMats.push_back(wedgeTipLabel(op1Mat,op2Mat));



    subTreeNum++;


    for(j = 0; j < wedgeOrder.size()-2; j++){

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

    wedgeOrder.erase(std::remove_if(wedgeOrder.begin(),
                                    wedgeOrder.end(),
                                    [](std::string x){return x == " ";}),
                                    wedgeOrder.end());

  }

  return(subCoeffMats.back());
}

// [[Rcpp::export]]
arma::cx_rowvec wedgeExportConv(arma::cx_rowvec A, arma::cx_rowvec B, arma::cx_double y){

  cx_rowvec res = wedgeConv(A,B,y);

  return(res);
}


// [[Rcpp::export]]
Rcpp::List coeffMatListEvalY(std::vector<std::vector<std::string>> wedgeOrders, arma::cx_double y, std::string tipLabA = " ", std::string tipLabB = " ", int nThreads = -1){
  int numCoeffs = wedgeOrders.size();

  size_t numThreads = std::thread::hardware_concurrency();
  if(nThreads != -1){
    numThreads = nThreads;
  }

  Rcpp::List output(numCoeffs);

  arma::field<arma::cx_rowvec> coeffs(numCoeffs);

  RcppThread::parallelFor(0, numCoeffs, [&coeffs, &wedgeOrders, &y] (unsigned int i) {
    coeffs[i] = coeffMatrixComplex(wedgeOrders[i],y);
  },numThreads,0);

  output = Rcpp::wrap(coeffs);


  return(output);
}

// [[Rcpp::export]]
Rcpp::List coeffMatListTipLabel(std::vector<std::vector<std::string>> wedgeOrders, arma::cx_double y, std::string tipLabA = " ", std::string tipLabB = " ", int nThreads = -1){
  int numCoeffs = wedgeOrders.size();

  size_t numThreads = std::thread::hardware_concurrency();
  if(nThreads != -1){
    numThreads = nThreads;
  }

  Rcpp::List output(numCoeffs);


  arma::field<arma::cx_mat> coeffs(numCoeffs);

  RcppThread::parallelFor(0, numCoeffs, [&coeffs, &wedgeOrders, &tipLabA, &tipLabB] (unsigned int i) {
    coeffs[i] = coeffMatrixTipLabel(wedgeOrders[i], tipLabA, tipLabB);
  },numThreads,0);

  output = Rcpp::wrap(coeffs);

  return(output);
}
