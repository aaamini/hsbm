
#ifndef __UTILS__
#define __UTILS__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

int find_tunc(arma::vec beta, double threshold);

arma::uvec get_freq(arma::uvec z, int K);
arma::uvec get_up_freq(arma::uvec freq);

arma::umat get_freq_minus_self(arma::uvec z, int K);

arma::vec fast_agg(arma::vec x, arma::uvec z, int K);

arma::uvec fast_agg_u(arma::uvec x, arma::uvec z, int K) ;

arma::mat comp_blk_sums(arma::sp_mat At, arma::uvec z, int Kcap);

arma::mat sp_compress_col(arma::sp_mat At, arma::uvec z, int Kcap);

List comp_blk_sums_and_sizes(arma::sp_mat At, arma::uvec z, int Kcap, bool div_diag = true);

arma::vec sp_single_col_compress(arma::sp_mat A, int col_idx, arma::uvec z, int Kcap);

arma::mat comp_blk_sums_diff(arma::sp_mat& A, int s, int zs_new, arma::uvec& z, int Kcap);

void print_progress(int itr, int itr_max);                    
#endif /* __UTILS__ */
