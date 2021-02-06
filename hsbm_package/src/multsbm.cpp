// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
#include <Rcpp/Benchmark/Timer.h>

#include "sampling.h" 
#include "utils.h" 

using namespace Rcpp;

// [[Rcpp::export]]
arma::umat multsbm_gibbs_sampler_fast(arma::sp_mat A, const int K, 
                    const double alpha = 1, const double beta = 1,
                    const int niter = 100) {

    int n = A.n_rows;

    arma::umat z_hist(n, niter, arma::fill::zeros);
    arma::uvec z = sample_int_vec(K, n);
    arma::mat B(K, K, arma::fill::zeros);   
    arma::vec pri(K); 

    for (int iter = 0; iter < niter; iter++) {
        List out = comp_blk_sums_and_sizes(A, z, K);

        arma::mat lambda = out["lambda"];
        arma::umat NN = out["NN"]; 
        B = symmat_rbeta(lambda + alpha, NN - lambda + beta);
        arma::mat uu = log(B/(1-B) + 1e-11);
        arma::mat vv = log(1-B + 1e-11);

        arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
        // Rcout << nn << "\n" << nn+1;
        pri = rdirichlet(nn + 1);

        for (int i = 0; i < n; i++) {
            arma::vec taui = sp_single_col_compress(A, i, z, K);
            arma::uvec mmi = get_freq(z, K);
            mmi(z(i))--;
            arma::vec Hi = uu * taui + vv * mmi +  log(pri);

            z(i) = sample_index(safe_exp(Hi));
        }
        z_hist.col(iter) = z;
    }
    return z_hist;
}

