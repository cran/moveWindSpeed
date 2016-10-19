// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calcVar
NumericVector calcVar(NumericMatrix groundSpeeds, NumericVector windSpeed, NumericVector phi);
RcppExport SEXP moveWindSpeed_calcVar(SEXP groundSpeedsSEXP, SEXP windSpeedSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type groundSpeeds(groundSpeedsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type windSpeed(windSpeedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(calcVar(groundSpeeds, windSpeed, phi));
    return rcpp_result_gen;
END_RCPP
}