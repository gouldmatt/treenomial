// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "wedge.h"

using namespace arma;

// [[Rcpp::export]]
arma::field<arma::mat> alignCoeffs_(arma::field<arma::mat> coeffs){
  unsigned int numCoeffs = coeffs.n_elem;

  uword maxSizeR = 0;
  uword maxSizeC = 0;

  /*
   *  find max row/col present, resize smaller matrices and then
   *  align so that max x exponent is aligned with larger matrix
   */

  Mat<double> currCoeff = coeffs[0];
  vec coeffRowSizes(numCoeffs);
  vec coeffColSizes(numCoeffs);

  for(unsigned int i = 0; i < numCoeffs; i++){
    Mat<double> currCoeff = coeffs[i];
    coeffRowSizes[i] = currCoeff.n_rows;
    coeffColSizes[i] = currCoeff.n_cols;
  }

  maxSizeR = coeffRowSizes.max();
  maxSizeC = coeffColSizes.max();


  for(unsigned int i = 0; i < numCoeffs; i++){
    Mat<double> currCoeff = coeffs[i];
    if(currCoeff.n_rows != maxSizeR || currCoeff.n_cols != maxSizeC){
      int coeffs_cols = currCoeff.n_cols;
      currCoeff.resize( maxSizeR, maxSizeC);
      currCoeff = shift(currCoeff,maxSizeC - coeffs_cols, 1);
      coeffs[i] = currCoeff;
    }
  }



  return(coeffs);
}


arma::mat latticePositions(std::vector<std::string> wedgeOrder, std::vector<std::string> wedgeOrderNodes, arma::mat L, int numTips){


  vec latticePos = L.col(2);

  vec dpt = L.col(3);
  vec pl = L.col(1);

  long unsigned int j = 0;
  int subTreeNum = 2;
  std::string op1;
  std::string op2;
  int pnode;
  int A;
  int B;
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


  std::map<int, int> subTreeNode;


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

        pnode = stoi(wedgeOrderNodes[j]);
        A = stoi(wedgeOrderNodes[j-2]);
        B = stoi(wedgeOrderNodes[j-1]);

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


        pnode = stoi(wedgeOrderNodes[j+2]);
        subTreeNode.insert(std::make_pair(pnode,subTreeNum-1));
        // std::cout << "making pair " << pnode << " " << subTreeNum << std::endl;


        wedgeOrderNodes[j] = " ";
        wedgeOrderNodes[j+1] = " ";


        wedgeOrder[j] = subPattern;
        wedgeOrder[j+1] = " ";
        wedgeOrder[j+2] = " ";

      } else if(wedgeOrder[j] == op2 && wedgeOrder[j+1] == op1 && wedgeOrder[j+2] == "1"){


        pnode = stoi(wedgeOrderNodes[j+2]);
        subTreeNode.insert(std::make_pair(pnode,subTreeNum-1));
        // std::cout << "making pair " << pnode << " " << subTreeNum << std::endl;


        wedgeOrderNodes[j] = " ";
        wedgeOrderNodes[j+1] = " ";

        wedgeOrder[j] = subPattern;
        wedgeOrder[j+1] = " ";
        wedgeOrder[j+2] = " ";

      }
    }

    wedgeOrderNodes.erase(std::remove_if(wedgeOrderNodes.begin(),
                                         wedgeOrderNodes.end(),
                                    [](std::string x){return x == " ";}),
                                    wedgeOrderNodes.end());

    // erase the repeated locations
    wedgeOrder.erase(std::remove_if(wedgeOrder.begin(),
                                    wedgeOrder.end(),
                                    [](std::string x){return x == " ";}),
                                    wedgeOrder.end());

  }



  A = 0;
  B = 0;
  mat A_poly;
  mat B_poly;


  uword h = max(dpt);

  for(uword i = 0; i <= h; i++){
    uvec CLNodes = find(dpt == i);
    for(uword j = 0; j < CLNodes.n_elem; j++){

      int pnode = CLNodes[j] + 1;
      // std::cout << "pnode " << pnode << std::endl;
        uvec children = find(pl == pnode);


      if(children.n_elem != 0){

        rowvec binX = L( pnode-1, span(4, 4+h));
        // std::cout << "h " << h << std::endl;
        // std::cout << "binx n_elem " << binX.n_elem << std::endl;
        rowvec binA = shift(binX,1);
        rowvec binB = shift(binX,1);
        binB[0] = 1;


        int x = latticePos[pnode-1];

        int a = 2*x;
        int b = 2*x+1;

        A = children[0] + 1;

        B = children[1] + 1;

        if(A <= numTips && B <= numTips){

          if(L(A-1,2) >=  L(B-1,2)){
            // latticePos[A-1] = a;
            // latticePos[B-1] = b;

            L(A-1, span(4, 4+h)) = binA;
            L(B-1, span(4, 4+h)) = binB;


          } else {


            L(A-1, span(4, 4+h)) = binB;
            L(B-1, span(4, 4+h)) = binA;

//
//             latticePos[A-1] = b;
//             latticePos[B-1] = a;

          }

        } else {
          // latticePos[A-1] = a;
          // latticePos[B-1] = b;

          L(A-1, span(4, 4+h)) = binA;
          L(B-1, span(4, 4+h)) = binB;



          int A_poly_index = subTreeNode.find(A)->second;
          int B_poly_index = subTreeNode.find(B)->second;
          // std::cout << "A " << A << std::endl;
          // std::cout << "B " << B << std::endl;
          // std::cout << A_poly_index << std::endl;
          // std::cout << B_poly_index << std::endl;

          if(A <= numTips && B > numTips){
            A_poly = mat(leaf);
            B_poly = mat(subCoeffMats[B_poly_index]);

          } else if(A > numTips && B <= numTips){
            A_poly = mat(subCoeffMats[A_poly_index]);
            B_poly = mat(leaf);

          } else {
            A_poly = mat(subCoeffMats[A_poly_index]);
            B_poly = mat(subCoeffMats[B_poly_index]);

            // std::cout << A_poly << std::endl;
            // std::cout << B_poly << std::endl;
          }

          arma::field<arma::mat> polys = alignCoeffs_( arma::field<arma::mat>({A_poly,B_poly}));

          A_poly = polys[0];
          B_poly = polys[1];


          // std::cout << A_poly - B_poly << std::endl;
          if(approx_equal(A_poly, B_poly, "absdiff", 0.1)){
            // std::cout << "A" << std::endl;
            // latticePos[A-1] = a;
            // latticePos[B-1] = b;
            L(A-1, span(4, 4+h)) = binA;
            L(B-1, span(4, 4+h)) = binB;
          } else {

            if(accu(A_poly - B_poly) >= 0){
              // std::cout << "B" << std::endl;
              // latticePos[A-1] = a;
              // latticePos[B-1] = b;
              L(A-1, span(4, 4+h)) = binA;
              L(B-1, span(4, 4+h)) = binB;
            } else {
              // std::cout << "C" << std::endl;
              // latticePos[A-1] = b;
              // latticePos[B-1] = a;
              L(A-1, span(4, 4+h)) = binB;
              L(B-1, span(4, 4+h)) = binA;
            }

          }


        }

      }


    }



  }

  return L;
}

vec vecUnion(vec x, vec y){
  return(unique(join_cols(x,y)));
}

// [[Rcpp::export]]
arma::vec setDiff(arma::vec x, arma::vec y){
  if(x.n_elem != 0 || y.n_elem != 0){
    uvec pos;
    for(uword i = 0; i < y.n_elem; i++){
        pos = find(x == y[i]);
        x.elem(pos) = zeros(pos.n_elem);
    }

    x = x.rows( find(x > 0));
    return(x);
  } else {
    return(x);
  }

}



int getNumSetDiff(vec &x, vec &y){
  // std::cout << vecUnion(setDiff(x,y),setDiff(y,x)) << std::endl;
  return(vecUnion(setDiff(x,y),setDiff(y,x)).n_elem);
}


// [[Rcpp::export]]
double lattDistance(const arma::mat &L1, const arma::mat &L2, const double w){
  // std::cout << "L1" << L1 << std::endl;
  // std::cout << "L2" << L2 << std::endl;

  vec CC;
  uvec iA;
  uvec iB;

  // vec lPos1 = L1.col(2);
  // vec lPos2 = L2.col(2);


  mat lPos1 = L1.cols(4, L1.n_cols-1);
  mat lPos2 = L2.cols(4, L2.n_cols-1);


  intersect(CC, iA, iB,lPos1,lPos2);

  std::cout << "iA" << iA << std::endl;
  std::cout << "iB" << iB << std::endl;
  std::cout << "CC" << CC << std::endl;


  // intersect(CC, iA, iB, L1.col(2) ,L2.col(2) );

  vec bl1 = L1.col(2);
  vec bl2 = L2.col(2);
  //
  // bl1 = bl1.elem(iA);
  // bl2 = bl2.elem(iB);
  //
  // bl1 = bl1.rows( find(bl1 > 0));
  // bl2 = bl2.rows( find(bl2 > 0));
  //
  // vec flat1 = lPos1.as_col();
  // vec flat2 = lPos2.as_col();
  //
  // int numSetDiff = getNumSetDiff(flat1, flat2);
  // // std::cout << "numsetdiff " << numSetDiff << std::endl;
  // return(w*numSetDiff + accu(abs(bl1 - bl2)/(bl1 + bl2) ));
  return 0;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix lattDistMat(Rcpp::List lattList, double w, int nThreads = -1){

  size_t numThreads = std::thread::hardware_concurrency();
  if(nThreads != -1){
    numThreads = nThreads;
  }

  std::vector<mat> latts = Rcpp::as<std::vector<mat>>(lattList);
  int listLength = latts.size();
  mat distMat(listLength, listLength, fill::zeros);

  // std::cout << distMat << std::endl;

  RcppThread::parallelFor(0, listLength, [&distMat,&listLength,&latts,&w] (unsigned int i) {
    for(int j = i+1; j < listLength; j++){
      distMat(i,j) = lattDistance(latts[i],latts[j], w);
    }
  }, numThreads,0);

  // std::cout << distMat << std::endl;

  distMat = distMat.t() + distMat;

  return(Rcpp::wrap(distMat));

}


// [[Rcpp::export]]
Rcpp::List latticeList(std::vector<std::vector<std::string>> wedgeOrders, std::vector<std::vector<std::string>> wedgeOrdersNodes, std::vector<arma::mat> latList,  std::vector<int> numTips, int nThreads = -1){
  int listLength = wedgeOrders.size();

  size_t numThreads = std::thread::hardware_concurrency();
  if(nThreads != -1){
    numThreads = nThreads;
  }

  Rcpp::List output(listLength);


  arma::field<arma::mat> lattices(listLength);

  RcppThread::parallelFor(0, listLength, [&lattices, &wedgeOrders, &wedgeOrdersNodes, &latList, &numTips] (unsigned int i) {
    lattices[i] = latticePositions(wedgeOrders[i], wedgeOrdersNodes[i], latList[i], numTips[i]);
  },numThreads,0);

  output = Rcpp::wrap(lattices);


  return(output);
}
