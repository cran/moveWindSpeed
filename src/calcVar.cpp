#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector calcVar(NumericMatrix groundSpeeds, NumericVector windSpeed, NumericVector phi){
NumericVector As= sqrt(pow((groundSpeeds(_,0) - windSpeed(0)),2)+pow((groundSpeeds(_,1) - windSpeed(1)),2));
int n =As.size();
NumericVector x=rep(1-phi,n-2);
x.push_back(1.0);
x.push_front(1.0);

NumericVector u= As*x;
NumericVector res= NumericVector::create(1.0);
NumericVector meanA= NumericVector::create(1.0);

meanA(0)=std::accumulate(u.begin(),u.end(), 0.0)/std::accumulate(x.begin(),x.end(),0.0);

res(0)=pow(As(0)-meanA(0),2)*(1-pow(phi(0),2));

for(int i=1; i<n; i++){
  res(0)=res(0)+pow(As(i)-((1-phi(0)) * meanA(0)+phi(0)*As(i-1)),2);
}
return(wrap(res(0)/(n-1)));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
