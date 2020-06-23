// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mcesteqn
arma::mat mcesteqn(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, arma::mat mcov, arma::uvec ind);
RcppExport SEXP _sme_mcesteqn(SEXP lbSEXP, SEXP mSEXP, SEXP nSEXP, SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP, SEXP mcovSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mcov(mcovSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(mcesteqn(lb, m, n, X, Y, beta, mcov, ind));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sme_mcesteqn", (DL_FUNC) &_sme_mcesteqn, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_sme(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}