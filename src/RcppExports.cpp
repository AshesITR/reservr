// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dgamma_matrix
NumericMatrix dgamma_matrix(NumericVector x, NumericVector shape, double scale);
RcppExport SEXP _reservr_dgamma_matrix(SEXP xSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(dgamma_matrix(x, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// pgamma_diff_matrix
NumericMatrix pgamma_diff_matrix(NumericVector lower, NumericVector upper, NumericVector shape, double scale);
RcppExport SEXP _reservr_pgamma_diff_matrix(SEXP lowerSEXP, SEXP upperSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(pgamma_diff_matrix(lower, upper, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// trunc_erlangmix_ellik
Rcpp::List trunc_erlangmix_ellik(const arma::vec& xmin, const arma::vec& xmax, const arma::vec& tmin, const arma::vec& tmax, const arma::vec& weight, const arma::vec& shapes, double scale, const arma::mat& zadj, bool parallel);
RcppExport SEXP _reservr_trunc_erlangmix_ellik(SEXP xminSEXP, SEXP xmaxSEXP, SEXP tminSEXP, SEXP tmaxSEXP, SEXP weightSEXP, SEXP shapesSEXP, SEXP scaleSEXP, SEXP zadjSEXP, SEXP parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tmin(tminSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type shapes(shapesSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type zadj(zadjSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(trunc_erlangmix_ellik(xmin, xmax, tmin, tmax, weight, shapes, scale, zadj, parallel));
    return rcpp_result_gen;
END_RCPP
}
// do_integrate_gk_lst
Rcpp::List do_integrate_gk_lst(const Rcpp::Function& fun, const arma::vec& lower, const arma::vec& upper, const Rcpp::List& params, const double tolerance, const int max_iter, bool debug);
RcppExport SEXP _reservr_do_integrate_gk_lst(SEXP funSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP paramsSEXP, SEXP toleranceSEXP, SEXP max_iterSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type fun(funSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(do_integrate_gk_lst(fun, lower, upper, params, tolerance, max_iter, debug));
    return rcpp_result_gen;
END_RCPP
}
// do_integrate_gk_mat
Rcpp::List do_integrate_gk_mat(const Rcpp::Function& fun, const arma::vec& lower, const arma::vec& upper, const arma::mat& params, const double tolerance, const int max_iter, bool debug);
RcppExport SEXP _reservr_do_integrate_gk_mat(SEXP funSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP paramsSEXP, SEXP toleranceSEXP, SEXP max_iterSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type fun(funSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(do_integrate_gk_mat(fun, lower, upper, params, tolerance, max_iter, debug));
    return rcpp_result_gen;
END_RCPP
}
// softmax_mat
arma::mat softmax_mat(arma::mat x);
RcppExport SEXP _reservr_softmax_mat(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(softmax_mat(x));
    return rcpp_result_gen;
END_RCPP
}
// softmax_vec
std::vector<double> softmax_vec(arma::vec x);
RcppExport SEXP _reservr_softmax_vec(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(softmax_vec(x));
    return rcpp_result_gen;
END_RCPP
}
// dsoftmax_vec
arma::mat dsoftmax_vec(arma::vec x);
RcppExport SEXP _reservr_dsoftmax_vec(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(dsoftmax_vec(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_reservr_dgamma_matrix", (DL_FUNC) &_reservr_dgamma_matrix, 3},
    {"_reservr_pgamma_diff_matrix", (DL_FUNC) &_reservr_pgamma_diff_matrix, 4},
    {"_reservr_trunc_erlangmix_ellik", (DL_FUNC) &_reservr_trunc_erlangmix_ellik, 9},
    {"_reservr_do_integrate_gk_lst", (DL_FUNC) &_reservr_do_integrate_gk_lst, 7},
    {"_reservr_do_integrate_gk_mat", (DL_FUNC) &_reservr_do_integrate_gk_mat, 7},
    {"_reservr_softmax_mat", (DL_FUNC) &_reservr_softmax_mat, 1},
    {"_reservr_softmax_vec", (DL_FUNC) &_reservr_softmax_vec, 1},
    {"_reservr_dsoftmax_vec", (DL_FUNC) &_reservr_dsoftmax_vec, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_reservr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
