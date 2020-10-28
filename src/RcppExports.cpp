// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// shrink_est
Rcpp::List shrink_est(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, arma::mat mcov, arma::uvec pos, arma::mat v);
RcppExport SEXP _eiv_shrink_est(SEXP lbSEXP, SEXP mSEXP, SEXP nSEXP, SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP, SEXP mcovSEXP, SEXP posSEXP, SEXP vSEXP) {
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
    Rcpp::traits::input_parameter< arma::uvec >::type pos(posSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(shrink_est(lb, m, n, X, Y, beta, mcov, pos, v));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mcesteqn
arma::vec rcpp_mcesteqn(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, arma::mat mcov, arma::uvec ind, bool modify_inv);
RcppExport SEXP _eiv_rcpp_mcesteqn(SEXP lbSEXP, SEXP mSEXP, SEXP nSEXP, SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP, SEXP mcovSEXP, SEXP indSEXP, SEXP modify_invSEXP) {
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
    Rcpp::traits::input_parameter< bool >::type modify_inv(modify_invSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mcesteqn(lb, m, n, X, Y, beta, mcov, ind, modify_inv));
    return rcpp_result_gen;
END_RCPP
}
// calculate_G
arma::mat calculate_G(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, arma::mat mcov, arma::uvec ind, arma::mat acov, arma::mat vinv, arma::vec us, arma::mat d, double rho);
RcppExport SEXP _eiv_calculate_G(SEXP lbSEXP, SEXP mSEXP, SEXP nSEXP, SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP, SEXP mcovSEXP, SEXP indSEXP, SEXP acovSEXP, SEXP vinvSEXP, SEXP usSEXP, SEXP dSEXP, SEXP rhoSEXP) {
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
    Rcpp::traits::input_parameter< arma::mat >::type acov(acovSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vinv(vinvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type us(usSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_G(lb, m, n, X, Y, beta, mcov, ind, acov, vinv, us, d, rho));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_inference
Rcpp::List rcpp_inference(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, arma::mat mcov, arma::uvec ind, bool modify_inv);
RcppExport SEXP _eiv_rcpp_inference(SEXP lbSEXP, SEXP mSEXP, SEXP nSEXP, SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP, SEXP mcovSEXP, SEXP indSEXP, SEXP modify_invSEXP) {
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
    Rcpp::traits::input_parameter< bool >::type modify_inv(modify_invSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_inference(lb, m, n, X, Y, beta, mcov, ind, modify_inv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eiv_shrink_est", (DL_FUNC) &_eiv_shrink_est, 9},
    {"_eiv_rcpp_mcesteqn", (DL_FUNC) &_eiv_rcpp_mcesteqn, 9},
    {"_eiv_calculate_G", (DL_FUNC) &_eiv_calculate_G, 13},
    {"_eiv_rcpp_inference", (DL_FUNC) &_eiv_rcpp_inference, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_eiv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
