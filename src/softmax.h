#ifndef SOFTMAX_H
#define SOFTMAX_H
#include<RcppArmadillo.h>

arma::mat dsoftmax_vec(arma::vec x);
arma::mat softmax_mat(arma::mat x);
std::vector<double> softmax_vec(arma::vec x);

#endif
