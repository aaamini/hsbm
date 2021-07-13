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


// [[Rcpp::export]]
arma::mat beta_fun_symmat(arma::mat a, arma::mat b) {
    // Computes beta function over symmetric matrix arguments
    // only computes the lower triangular part + diagonal; set the rest to 1
    int K = a.n_rows;
    arma::mat out(K, K, arma::fill::ones);
    for (int k = 0; k < K; k++) {
        for (int el = 0; el <= k; el++) {
            out(k, el) = R::beta(a(k, el), b(k, el));
        }
    }
    return out;
}

// [[Rcpp::export]]
arma::mat comp_beta_matrix(arma::sp_mat& A, arma::uvec& z, const int K, double alpha, double beta) {
    List out = comp_blk_sums_and_sizes(A, z, K);
    arma::mat lambda = out["lambda"];
    arma::umat NN = out["NN"]; 
            
    return beta_fun_symmat(lambda + alpha, NN - lambda + beta);
}


// [[Rcpp::export]]
arma::umat multsbm_collapsed_gibbs_sampler(arma::sp_mat& A, const int K, 
                    const double alpha = 1, const double beta = 1,
                    const int niter = 10) {

    int n = A.n_rows;
    arma::umat z_hist(n, niter, arma::fill::zeros);
    
    // Randomly initialize labels
    arma::uvec z = sample_int_vec(K, n);
    arma::uvec z_new(n, arma::fill::zeros);
    arma::vec pri(K, arma::fill::zeros);
    arma::mat Bet(K, K);

    for (int iter = 0; iter < niter; iter++) {
        for (int s = 0; s < n; s++) {
            Bet = comp_beta_matrix(A, z, K, alpha, beta);
            arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
            pri = rdirichlet(nn + 1);

            arma::vec prob(K, arma::fill::ones);
            for (int rp = 0; rp < K; rp++) { // rp is the potential new value of z(s)
                z_new = z;
                z_new(s) = rp;
                arma::mat Bet_new =  comp_beta_matrix(A, z_new, K, alpha, beta);;
                prob(rp) *= arma::prod(arma::prod(Bet_new / Bet)) * pri(rp);
                // prob(rp) *= static_cast<double>((nn(rp) + 1)) / nn(z(s)); 
            }

            z(s) = sample_index(prob);; // update z
        }
        z_hist.col(iter) = z;
    }

    return z_hist;    
    
}
