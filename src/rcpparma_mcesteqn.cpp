// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void binomial_fun(double a, double b, double ay, arma::rowvec& tmp,
                        arma::rowvec& u, arma::mat mcov_j, arma::uvec pos,
                        arma::mat& d2) {

    double p1 = exp(-a/2 - b/8) * ay;
    double p2 = exp(a/2 - b/8) * (ay-1);
    arma::rowvec u1 = u;
    u1.elem(pos) += -tmp/2;
    arma::mat d1 = u1.t() * u1;
    u = (p2 + p1) * u;
    tmp = (p1 - p2) * tmp/2;
    d2(pos, pos) = mcov_j;
    d2 = (p2 - p1)/2 * (d2 - d1);
}

// [[Rcpp::export]]
void gaussian_fun(double a, double b, double ay, arma::rowvec& tmp,
                  arma::rowvec& u, arma::mat mcov_j, arma::uvec pos,
                  arma::mat& d2) {
    arma::mat d1 = u.t() * u;
    u = (ay - a) * u.t();
    tmp = b;
    d2(pos, pos) = mcov_j;
    d2 = d2 - d1;
}

// [[Rcpp::export]]
Rcpp::List shrink_est(int lb, int m, int n, Rcpp::List X, Rcpp::List Y,
                     arma::vec beta, Rcpp::List mcov, arma::uvec pos,
                     arma::mat v) {
    double b_bar2 = 0;
    for (int i = 0; i < n; i++) {
        arma::mat ui = arma::zeros(lb*m);
        arma::mat X1 = X[i];
        arma::vec Y1 = Y[i];
        for (int j = 0; j < X1.n_rows; j++) {
            arma::rowvec u = X1.row(j);
            arma::mat mcov_j = mcov[j];
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov_j;
            double b = as_scalar(tmp * beta.elem(pos));
            double ay = as_scalar(Y1.row(j));
            double weight1 = exp(-a/2 - b/8);
            double weight2 = exp(a/2 - b/8);
            ui.rows(j*lb, (j + 1L)*lb - 1L) = (weight2*(ay-1) +
                weight1*ay)*u.t();
            ui.elem(j*lb + pos) += (weight1*ay - weight2*(ay-1)) * tmp.t()/2;

        }
        //scale it first
        arma::mat vi = ui * ui.t();
        // scaled it check this...
        double b_bar2i = pow(arma::norm(vi - v, "fro") / n, 2.0);
        b_bar2 += b_bar2i;
    }
    arma::mat I;
    I.eye(size(v));
    double m_n = arma::mean(arma::trace(v * I));
    double d2 = pow(arma::norm((v - m_n * I), "fro"), 2.0);
    double b2 = std::min(d2, b_bar2);
    double a2 = d2 - b2;
    arma::mat v_m = b2/d2 * m_n * I + a2/d2 * v;
    double rho = b2/d2;
    return Rcpp::List::create(
        Rcpp::Named("v_m") = v_m,
        Rcpp::Named("rho") = rho
    );

}
// [[Rcpp::export]]
arma::vec rcpp_mcesteqn(int lb, int m, int n, Rcpp::List X,
                        Rcpp::List Y,
                        arma::vec beta,
                        Rcpp::List mcov,
                        arma::uvec ind,
                        bool modify_inv) {
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
            arma::mat mcov_j = mcov[j];
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov_j;
            double b = as_scalar(tmp * beta.elem(pos));
            double ay = as_scalar(Y1.row(j));
            arma::mat d2 = arma::zeros(lb, lb);
            binomial_fun(a, b, ay, tmp, u, mcov_j, pos, d2);
            ui.rows(j*lb, (j + 1L)*lb - 1L) = u.t();
            ui.elem(j*lb + pos) += tmp.t();
            di.cols(j*lb, (j + 1L)*lb - 1L) = d2;
        }
        arma::mat vi = ui * ui.t();
        d += di;
        v += vi;
        us += ui;
    }
    us = us/n;
    arma::mat vi = us * us.t();
    v = v/n;
    //insert the calculation here; follow that lnshrink
    if (modify_inv) {
        v = v - vi ;
        Rcpp::List outlist = shrink_est(lb, m, n, X, Y, beta, mcov, pos, v);
        arma::mat v_m = outlist["v_m"];
        v = v_m;
    }
    v = v/n;
    d = d/n;
    arma::mat vinv = inv(v);
    arma::mat dold = d * vinv;
    arma::vec out = dold * us;
    return out;

}

// [[Rcpp::export]]
arma::mat calculate_G(int lb, int m, int n,
                      Rcpp::List X, Rcpp::List Y,
                      arma::vec beta,
                      Rcpp::List mcov,
                      arma::uvec ind,
                      arma::mat acov,
                      arma::mat vinv,
                      arma::vec us,
                      arma::mat d) {
    arma::uvec pos = ind - 1;
    arma::mat g = arma::zeros(lb, lb);
    for (int i = 0; i < n; i++) {
        arma::mat ui = arma::zeros(lb*m);
        arma::mat di = arma::zeros(lb, lb*m);
        arma::mat X1 = X[i];
        arma::vec Y1 = Y[i];
        arma::mat gi = arma::zeros(lb, lb);
        for (int j = 0; j < X1.n_rows; j++) {
            arma::rowvec u = X1.row(j);
            arma::mat mcov_j = mcov[j];
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov_j;
            double b = as_scalar(tmp * beta.elem(pos));
            double ay = as_scalar(Y1.row(j));
            arma::mat d2 = arma::zeros(lb, lb);
            binomial_fun(a, b, ay, tmp, u, mcov_j, pos, d2);
            ui.rows(j*lb, (j + 1L)*lb - 1L) = u.t();
            ui.elem(j*lb + pos) += tmp.t();
            di.cols(j*lb, (j + 1L)*lb - 1L) = d2;
        }
        for (arma::uword k = 0; k < lb; k++) {
            arma::rowvec dk = di.rows(k, k);
            gi.cols(k, k) = acov * d * vinv * (ui * dk + dk.t() * ui.t()) *
                vinv * us;
        }
        g += gi / n;
    }
    g = g/n;
    arma::mat I_lb = arma::eye(lb, lb);
    arma::mat term = I_lb + g;
    arma::mat acov_c = term * acov * term.t();
    return acov_c;
}

// [[Rcpp::export]]
Rcpp::List rcpp_inference(int lb, int m, int n,
                          Rcpp::List X, Rcpp::List Y,
                          arma::vec beta,
                          Rcpp::List mcov,
                          arma::uvec ind,
                          bool finsam_cor,
                          bool modify_inv) {
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
            arma::mat mcov_j = mcov[j];
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov_j;
            double b = as_scalar(tmp * beta.elem(pos));
            double ay = as_scalar(Y1.row(j));
            arma::mat d2 = arma::zeros(lb, lb);
            binomial_fun(a, b, ay, tmp, u, mcov_j, pos, d2);
            ui.rows(j*lb, (j + 1L)*lb - 1L) = u.t();
            ui.elem(j*lb + pos) += tmp.t();
            di.cols(j*lb, (j + 1L)*lb - 1L) = d2;
        }
        arma::mat vi = ui * ui.t();
        d += di;
        v += vi;
        us += ui;
    }
    us = us/n;
    arma::mat vi = us * us.t();
    v = v/n;
    if (modify_inv) {
        v = v - vi;
        Rcpp::List outlist = shrink_est(lb, m, n, X, Y, beta, mcov, pos, v);
        arma::mat v_m = outlist["v_m"];
        v = v_m;
    }
    v = v/n;
    d = d/n;
    arma::mat vinv = inv(v);
    arma::mat dold = d * vinv;
    arma::mat acovinv = dold * d.t();
    arma::mat acov;
    acov = inv(acovinv);
    arma::vec out = dold * us;
    if (finsam_cor) {
        arma::mat acov_c = calculate_G(lb, m, n,
                                     X, Y,
                                     beta,
                                     mcov,
                                     ind,
                                     acov,
                                     vinv,
                                     us, d);
        acov = acov_c;
    }

    return Rcpp::List::create(
        Rcpp::Named("v") = v,
        Rcpp::Named("vinv") = vinv,
        Rcpp::Named("d") = d,
        Rcpp::Named("us") = us,
        Rcpp::Named("acov") = acov,
        Rcpp::Named("mcesteqn") = out
    );

}




