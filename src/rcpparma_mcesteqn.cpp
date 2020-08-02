// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec rcpp_mcesteqn(int lb, int m, int n, Rcpp::List X, Rcpp::List Y,
                        arma::vec beta,
                        arma::mat mcov,
                        arma::uvec ind,
                        bool modify_vinv) {
    // maybe I should move this to other functions
    arma::uvec pos = ind - 1;
    arma::mat d = arma::zeros(lb, lb*m);
    arma::mat v = arma::zeros(lb*m, lb*m);
    arma::vec us = arma::zeros(lb*m);
    for (int i = 0; i < n; i++) {
        arma::mat ui = arma::zeros(lb*m);
        arma::mat di = arma::zeros(lb, lb*m);
        arma::mat X1 = X[i];
        arma::vec Y1 = Y[i];
        for (int j = 0; j < X1.n_rows; j++) {
            // do I need to use clone here?
            arma::rowvec u = X1.row(j);
            // figure out the dimension here what is double what is not?
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov;
            double b = as_scalar(tmp * beta.elem(pos));
            // not necessarily weights; think about other families for future design
            double weight = exp(-a - b/2) + 1;
            double ay = as_scalar(Y1.row(j));
            // think about the dimension
            ui.rows(j*lb, (j + 1L)*lb - 1L) = (weight*ay-1L)*u.t();
            ui.elem(j*lb + pos) += (weight-1L)*ay*tmp.t();
            arma::rowvec u1 = u;
            u1.elem(pos) += tmp;
            arma::mat d1 = u1.t() * u1;
            arma::mat d2 = arma::zeros(lb, lb);
            d2(pos, pos) = mcov;
            di.cols(j*lb, (j + 1L)*lb - 1L) = ay*(weight-1L)*(d2 - d1);
        }
        arma::mat vi = ui * ui.t();
        d += di;
        v += vi;
        us += ui;
    }
    us = us/n;
    arma::mat vi = us * us.t();
    v = v/n - vi;
    v = v/n;
    double  k = 1/n;
    if (modify_vinv) v = v + k * arma::eye(lb*m, lb*m);
    d = d/n;
    arma::mat vinv = inv(v);
    arma::mat dold = d * vinv;
    arma::vec out = dold * us;
    return out;

}


// [[Rcpp::export]]
Rcpp::List rcpp_inference(int lb, int m, int n, Rcpp::List X, Rcpp::List Y,
                          arma::vec beta,
                          arma::mat mcov,
                          arma::uvec ind) {
    arma::uvec pos = ind - 1;
    arma::mat d = arma::zeros(lb, lb*m);
    arma::mat v = arma::zeros(lb*m, lb*m);
    arma::vec us = arma::zeros(lb*m);
    for (int i = 0; i < n; i++) {
        arma::mat ui = arma::zeros(lb*m);
        arma::mat di = arma::zeros(lb, lb*m);
        arma::mat X1 = X[i];
        arma::vec Y1 = Y[i];
        for (int j = 0; j < X1.n_rows; j++) {
            arma::rowvec u = X1.row(j);
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov;
            double b = as_scalar(tmp * beta.elem(pos));
            double weight = exp(-a - b/2) + 1;
            double ay = as_scalar(Y1.row(j));
            // think about the dimension
            ui.rows(j*lb, (j + 1L)*lb - 1L) = (weight*ay-1L)*u.t();
            ui.elem(j*lb + pos) += (weight-1L)*ay*tmp.t();
            arma::rowvec u1 = u;
            u1.elem(pos) += tmp;
            arma::mat d1 = u1.t() * u1;
            arma::mat d2 = arma::zeros(lb, lb);
            d2(pos, pos) = mcov;
            di.cols(j*lb, (j + 1L)*lb - 1L) = ay*(weight-1L)*(d2 - d1);
        }
        arma::mat vi = ui * ui.t();
        d += di;
        v += vi;
        us += ui;
    }
    us = us/n;
    arma::mat vi = us * us.t();
    v = v/n - vi;
    v = v/n;
    double  k = 1/n;
    arma::vec s = svd(v);
    double cond = s.max()/s.min();
    bool modify_vinv = 0;
    if (cond == arma::datum::inf) {
        v = v + k * arma::eye(lb*m, lb*m);
        modify_vinv = 1;
    }
    d = d/n;
    arma::mat vinv = inv(v);
    arma::mat dold = d * vinv;
    arma::mat acovinv = dold * d.t();
    arma::vec c = svd(acovinv);
    double cond_c = c.max()/c.min();
    bool modify_acovinv = 0;
    if (cond_c == arma::datum::inf) {
        acovinv = acovinv + k * arma::eye(lb, lb);
        modify_acovinv = 1;
    }
    arma::mat acov = inv(acovinv);
    arma::vec out = dold * us;
    return Rcpp::List::create(
        Rcpp::Named("v") = v,
        Rcpp::Named("vinv") = vinv,
        Rcpp::Named("d") = d,
        Rcpp::Named("us") = us,
        Rcpp::Named("acov") = acov,
        Rcpp::Named("modify_vinv") = modify_vinv,
        Rcpp::Named("modify_acovinv") = modify_acovinv,
        Rcpp::Named("mcesteqn") = out
    );

}
