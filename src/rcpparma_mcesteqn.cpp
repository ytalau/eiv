// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat nearPD(arma::mat m,
                 int maxit,
                 double eig_tol,
                 double conv_tol) {
    // check whether the matrix is symmetric
    int iter = 0;
    bool converged = false;
    arma::vec eigval;
    arma::mat eigvec;
    arma::mat X = m;
    while (iter < maxit && ! converged) {
        arma::mat Y = X;
        eig_sym(eigval, eigvec, Y);
        arma::uvec ids = find(eigval >= eig_tol * eigval.min());
        arma::mat Q = eigvec.cols(ids);
        arma::mat vector = eigval.elem(ids);
        int n = ids.n_elem;
        arma::vec tmp = arma::vectorise(arma::repelem(vector,
                                                      n, 1));
        arma::vec Q_v = arma::vectorise(Q) % tmp;
        arma::mat Q_tmp(Q_v);
        Q_tmp.reshape(Q.n_rows, Q.n_cols);
        arma::mat X = Q_tmp * Q.t();
        double conv = arma::norm(Y-X, "inf") / arma::norm(Y, "inf");
        iter = iter + 1L;
        // do I need a check for negative semi-definite?
        converged = (conv <= conv_tol);
    }
    return X;
}

arma::mat neaRPD(arma::mat m, int maxit,
                 double eig_tol,
                 double conv_tol){

    Rcpp::Environment Matrix("package:Matrix"); // Load the Matrix package in R!
    Rcpp::Function nearPD = Matrix["nearPD"];   // Extract nearPD() R function

    // Compute with R function an S4 object
    Rcpp::List PD = nearPD(m, Rcpp::Named("maxit", maxit),
                           Rcpp::Named("eig.tol", eig_tol),
                           Rcpp::Named("conv.tol", conv_tol));
    Rcpp::S4 D_s4 = PD["mat"];

    // Convert the S4 object to an Armadillo matrix
    Rcpp::NumericVector temp = Rcpp::NumericVector(D_s4.slot("x"));
    Rcpp::NumericVector dims = D_s4.slot("Dim");

    // Advanced armadillo matrix ctor that reuses memory
    arma::mat D(temp.begin(), // pointer to NumericVector
                dims[0],      // Number of Rows
                    dims[1],      // Number of Columns
                        false,        // Avoid copying by disabling `copy_aux_mem`
                        true);
    return D;
}


arma::mat inv_mod(arma::mat m,
                  int maxit,
                  double eig_tol,
                  double conv_tol) {
    // first check whether it is invertible
    // double c = arma::cond(m);
    // if (c > 1e6) {
    //     m = neaRPD(m, maxit, eig_tol, conv_tol);
    // }
    m = neaRPD(m, maxit, eig_tol, conv_tol);
    arma::mat inv = arma::inv(m);
    return inv;
}



// [[Rcpp::export]]
arma::vec rcpp_mcesteqn(int lb, int m, int n, Rcpp::List X, Rcpp::List Y,
                        arma::vec beta,
                        arma::mat mcov,
                        arma::uvec ind,
                        int maxit,
                        double eig_tol,
                        double conv_tol,
                        bool modify_inv) {
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
            double weight1 = exp(-a/2 - b/8);
            double weight2 = exp(a/2 - b/8);
            double ay = as_scalar(Y1.row(j));
            // think about the dimension
            ui.rows(j*lb, (j + 1L)*lb - 1L) = (weight2*(ay-1) +
                weight1*ay)*u.t();
            ui.elem(j*lb + pos) += (weight1*ay - weight2*(ay-1))*tmp.t()/2;
            arma::rowvec u1 = u;
            u1.elem(pos) += -tmp/2;
            arma::mat d1 = u1.t() * u1;
            arma::mat d2 = arma::zeros(lb, lb);
            d2(pos, pos) = mcov;
            di.cols(j*lb, (j + 1L)*lb - 1L) = (weight2*(ay-1) -
                weight1*ay)/2*(d2 - d1);

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
    arma::mat vinv;
    if (modify_inv) {
        vinv = inv_mod(v, maxit, eig_tol, conv_tol);
    } else {
        vinv = inv(v);
    }
    arma::mat dold = d * vinv;
    arma::vec out = dold * us;
    return out;

}

// [[Rcpp::export]]
arma::mat calculate_G(int lb, int m, int n,
                      Rcpp::List X, Rcpp::List Y,
                      arma::vec beta,
                      arma::mat mcov,
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
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov;
            double b = as_scalar(tmp * beta.elem(pos));
            double weight1 = exp(-a/2 - b/8);
            double weight2 = exp(a/2 - b/8);
            double ay = as_scalar(Y1.row(j));
            // think about the dimension
            ui.rows(j*lb, (j + 1L)*lb - 1L) = (weight2*(ay-1) +
                weight1*ay)*u.t();
            ui.elem(j*lb + pos) += (weight1*ay - weight2*(ay-1))*tmp.t()/2;
            arma::rowvec u1 = u;
            u1.elem(pos) += -tmp/2;
            arma::mat d1 = u1.t() * u1;
            arma::mat d2 = arma::zeros(lb, lb);
            d2(pos, pos) = mcov;
            di.cols(j*lb, (j + 1L)*lb - 1L) = (weight2*(ay-1) -
                weight1*ay)/2*(d2 - d1);
        }
        for (arma::uword k = 0; k < lb; k++) {
            arma::rowvec dk = di.rows(k, k);
            gi.cols(k, k) = acov * d * vinv * (ui * dk + dk.t() * ui.t()) * vinv * us;
        }
        g += gi;
    }
    g = g/n;
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
                          arma::mat mcov,
                          arma::uvec ind,
                          bool bootstrap,
                          arma::mat meat,
                          int maxit,
                          double eig_tol,
                          double conv_tol,
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
            double a = as_scalar(u * beta);
            arma::rowvec tmp = beta.elem(pos).t() * mcov;
            double b = as_scalar(tmp * beta.elem(pos));
            double weight1 = exp(-a/2 - b/8);
            double weight2 = exp(a/2 - b/8);
            double ay = as_scalar(Y1.row(j));
            // think about the dimension
            ui.rows(j*lb, (j + 1L)*lb - 1L) = (weight2*(ay-1) +
                weight1*ay)*u.t();
            ui.elem(j*lb + pos) += (weight1*ay - weight2*(ay-1))*tmp.t()/2;
            arma::rowvec u1 = u;
            u1.elem(pos) += -tmp/2;
            arma::mat d1 = u1.t() * u1;
            arma::mat d2 = arma::zeros(lb, lb);
            d2(pos, pos) = mcov;
            di.cols(j*lb, (j + 1L)*lb - 1L) = (weight2*(ay-1) -
                weight1*ay)/2*(d2 - d1);
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
    arma::mat I;
    I.eye(size(v));
    double m_n = arma::accu(arma::trace(vi * I))/(lb * m);
    double d2 = pow(arma::norm((vi - m_n*I), "fro"), 2.0)/(lb * m);
    // arma::mat tmp = arma::zeros(lb*m, lb*m);
    // for (int k = 0; k < lb*m; k++) {
    //     arma::colvec usk = us.row(k);
    //     arma::mat tmpk = usk * usk.t() - v;
    //     tmp += tmpk;
    // }
    double b_bar2 = pow(arma::norm(v, "fro"), 2.0)/pow(lb * m, 2.0);
    double b2 = std::min(b2, b_bar2);
    double a2 = d2 - b2;
    arma::mat v_m = b2/d2 * m_n * I + a2/d2 * v;
    arma::mat vinv;
    if (modify_inv) {
        vinv = inv(v_m);
    } else {
        vinv = inv(v);
    }
    arma::mat dold = d * vinv;
    arma::mat acovinv = dold * d.t();
    arma::mat acov;
    // if (modify_inv) {
    //     acov = inv_mod(acovinv, maxit, eig_tol, conv_tol);
    // } else {
    //     acov = inv(acovinv);
    // }
    // acov = inv(acovinv);
    acov = inv(acovinv);
    if (bootstrap) {
        acov = acov * meat * acov * n;
    }
    arma::vec out = dold * us;
    arma::mat acov_c = calculate_G(lb, m, n,
                                   X, Y,
                                   beta,
                                   mcov,
                                   ind,
                                   acov,
                                   vinv,
                                   us, d);
    return Rcpp::List::create(
        Rcpp::Named("v") = v,
        Rcpp::Named("vinv") = vinv,
        Rcpp::Named("d") = d,
        Rcpp::Named("us") = us,
        Rcpp::Named("acov") = acov_c,
        Rcpp::Named("mcesteqn") = out
    );

}




