// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// binomial_fun
void binomial_fun(double a, double b, double ay, arma::rowvec& tmp, arma::rowvec& u, arma::mat mcov_j, arma::uvec pos, arma::mat& d2);
RcppExport SEXP _eiv_binomial_fun(SEXP aSEXP, SEXP bSEXP, SEXP aySEXP, SEXP tmpSEXP, SEXP uSEXP, SEXP mcov_jSEXP, SEXP posSEXP, SEXP d2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type ay(aySEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type tmp(tmpSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mcov_j(mcov_jSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos(posSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type d2(d2SEXP);
    binomial_fun(a, b, ay, tmp, u, mcov_j, pos, d2);
    return R_NilValue;
END_RCPP
}
// shrink_est
Rcpp::List shrink_est(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, Rcpp::List mcov, arma::uvec pos, arma::mat v);
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
    Rcpp::traits::input_parameter< Rcpp::List >::type mcov(mcovSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos(posSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(shrink_est(lb, m, n, X, Y, beta, mcov, pos, v));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mcesteqn
arma::vec rcpp_mcesteqn(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, Rcpp::List mcov, arma::uvec ind, bool modify_inv);
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
    Rcpp::traits::input_parameter< Rcpp::List >::type mcov(mcovSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< bool >::type modify_inv(modify_invSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mcesteqn(lb, m, n, X, Y, beta, mcov, ind, modify_inv));
    return rcpp_result_gen;
END_RCPP
}
// calculate_G
arma::mat calculate_G(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, Rcpp::List mcov, arma::uvec ind, arma::mat acov, arma::mat vinv, arma::vec us, arma::mat d);
RcppExport SEXP _eiv_calculate_G(SEXP lbSEXP, SEXP mSEXP, SEXP nSEXP, SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP, SEXP mcovSEXP, SEXP indSEXP, SEXP acovSEXP, SEXP vinvSEXP, SEXP usSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mcov(mcovSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type acov(acovSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vinv(vinvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type us(usSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_G(lb, m, n, X, Y, beta, mcov, ind, acov, vinv, us, d));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_inference
Rcpp::List rcpp_inference(int lb, int m, int n, Rcpp::List X, Rcpp::List Y, arma::vec beta, Rcpp::List mcov, arma::uvec ind, bool finsam_cor, bool modify_inv);
RcppExport SEXP _eiv_rcpp_inference(SEXP lbSEXP, SEXP mSEXP, SEXP nSEXP, SEXP XSEXP, SEXP YSEXP, SEXP betaSEXP, SEXP mcovSEXP, SEXP indSEXP, SEXP finsam_corSEXP, SEXP modify_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mcov(mcovSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< bool >::type finsam_cor(finsam_corSEXP);
    Rcpp::traits::input_parameter< bool >::type modify_inv(modify_invSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_inference(lb, m, n, X, Y, beta, mcov, ind, finsam_cor, modify_inv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eiv_binomial_fun", (DL_FUNC) &_eiv_binomial_fun, 8},
    {"_eiv_shrink_est", (DL_FUNC) &_eiv_shrink_est, 9},
    {"_eiv_rcpp_mcesteqn", (DL_FUNC) &_eiv_rcpp_mcesteqn, 9},
    {"_eiv_calculate_G", (DL_FUNC) &_eiv_calculate_G, 12},
    {"_eiv_rcpp_inference", (DL_FUNC) &_eiv_rcpp_inference, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_eiv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
