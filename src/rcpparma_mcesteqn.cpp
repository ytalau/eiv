// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mcesteqn(int lb, int m, int n, Rcpp::List X, Rcpp::List Y,
                   arma::vec beta,
                   arma::mat mcov,
                   arma::uvec ind) {
    // maybe I should move this to other functions
    arma::uvec pos = ind - 1;
    arma::mat d = arma::zeros(lb, lb*m);
    arma::mat v = arma::zeros(lb*m, lb*m);
    arma::mat us = arma::zeros(lb*m);
    for (int i = 0; i < n; i++) {
        arma::mat ui = arma::zeros(lb*m);
        arma::mat di = arma::zeros(lb, lb*m);
        arma::mat X1 = X[i];
        arma::vec Y1 = Y[i];
        for (int j = 0; j < X1.n_rows; j++) {
            // do I need to use clone here?
            arma::mat u = X1.row(j);
            // figure out the dimension here what is double what is not?
            double a = as_scalar(u * beta);
            arma::mat tmp = beta.elem(pos).t() * mcov;
            arma::mat b = tmp * beta.elem(pos);
            // not necessarily weights; think about other families for future design
            double weight = as_scalar(exp(-a - b/2) + 1L);
            double ay = as_scalar(Y1.row(j));
            // think about the dimension
            ui.rows(j*lb, (j + 1L)*lb - 1L) = (weight*ay-1L)*u.t();
            ui.elem(j*lb + pos) += (weight-1L)*ay*tmp.t();
            arma::mat u1 = u;
            u1.elem(pos) += tmp;
            arma::mat d1 = u1.t() * u1;
            arma::mat d2 = arma::zeros(lb, lb);
            d2(pos, pos) = mcov;
            di.cols(j*lb, (j + 1L)*lb - 1L) = ay*(a-1L)*(d2 - d1);
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
    d = d/n;
    arma::mat vinv = inv(v);
    arma::mat dold = d * vinv;
    arma::mat out = dold * us;
    return dold;

}

