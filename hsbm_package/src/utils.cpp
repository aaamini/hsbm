// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "sampling.h"

using namespace Rcpp;


void print_progress(int itr, int itr_max) {
  int width = ceil(log10(itr_max));
  if (itr % 10 == 0) Rcout << "("
                            << std::setw(width) << itr 
                            << " / " 
                            << std::setw(width) << itr_max << ")\r";
}


// List gem_posterior_counts(arma::uvec z, int K) {
//     arma::uvec count1(K, arma::fill::zeros);
//     arma::uvec count2(K, arma::fill::zeros);
//     int n = z.n_elem;

    
//     for (int k = 0; k < K; k++){
//       for (int i = 0; i < n; i++){
//         if (z(i) == k) {
//           count1(k)++;
//         } else if (z(i) > k) {
//           count2(k)++;
//         }
//       } //i
//     }      

//     return List::create(Named("count1") = count1, Named("count2") = count2);
// }

// [[Rcpp::export]]
arma::uvec get_up_freq(arma::uvec freq) {
    
    const int K = freq.n_elem;
    arma::uvec up_freq(K, arma::fill::zeros);

    for (int k = K-2; k >= 0; k--) {
        up_freq(k) = up_freq(k+1) + freq(k+1);
    }
    return up_freq;
}


// [[Rcpp::export]]
int find_tunc(arma::vec beta, double threshold) {
  int n = beta.n_elem;
  int idx = n;
  double cumsum = 0.0;
  for (idx = 0; idx < n - 1; idx++) {
    cumsum += beta(idx);
    if (cumsum > 1 - threshold) break;
  }
  return idx;
}


// [[Rcpp::export]]
arma::vec fast_agg(arma::vec x, arma::uvec z, int K) {
    int n = z.n_elem;

    arma::vec S(K, arma::fill::zeros);

    for (int i = 0; i < n; i++) {
        S(z(i)) += x(i);
    }
    return S;
}

arma::uvec fast_agg_u(arma::uvec x, arma::uvec z, int K) {
    int n = z.n_elem;

    arma::uvec S(K, arma::fill::zeros);

    for (int i = 0; i < n; i++) {
        S(z(i)) += x(i);
    }
    return S;
}


// counts the occurence frequency of integers 0, 1, 2, ..., K-1 in z
// the numbers in z should be in [0, K-1], otherwise an error occurs
// [[Rcpp::export]]
arma::uvec get_freq(arma::uvec z, int K) { 
  int n = z.n_elem;
  // int K = max(z)+1;
  arma::uvec freq(K,  arma::fill::zeros);
  // Rcout << K << freq; 
  for (int i = 0; i < n; i++) {
    freq(z(i))++;
  }
  return freq;
}

// [[Rcpp::export]]
arma::umat get_freq_minus_self(arma::uvec z, int K) {
    int n = z.n_elem;
    arma::uvec zcounts = get_freq(z, K);
    arma::umat freq = arma::repmat(zcounts.t(), n, 1);

    for (int i = 0; i < n; i++) {
        freq(i, z(i))--;
    }
    return freq;
}


// [[Rcpp::export]]
arma::mat comp_blk_sums(arma::sp_mat At, arma::uvec z, int Kcap) {
    // Compute block sums of a sparse matrix At w.r.t. labels in "z". The labels
    // in z are in the interval [0,1,...,Kcap]

    arma::sp_mat::const_iterator it     = At.begin();
    arma::sp_mat::const_iterator it_end = At.end();

    arma::mat lambda(Kcap, Kcap, arma::fill::zeros); 
    for(; it != it_end; ++it) {
        lambda(z(it.row()), z(it.col())) += (*it);
    }

    return lambda;
}

// [[Rcpp::export]]
arma::mat sp_compress_col(arma::sp_mat At, arma::uvec z, int Kcap) {
    int n = At.n_rows;
    arma::mat B(n, Kcap, arma::fill::zeros);

    arma::sp_mat::const_iterator it     = At.begin();
    arma::sp_mat::const_iterator it_end = At.end();

    for(; it != it_end; ++it) {
        B(it.row(), z(it.col())) += (*it);
    }

    return B;
}

// [[Rcpp::export]]
arma::vec sp_single_col_compress(arma::sp_mat A, int col_idx, arma::uvec z, int Kcap) {

    arma::sp_mat Acol(A.col(col_idx));
    arma::vec b(Kcap, arma::fill::zeros);

    // Rcpp::print(wrap(Acol));
    for (arma::sp_mat::iterator it = Acol.begin(); it != Acol.end(); ++it) {
        b(z(it.row())) += (*it);
        //Rcout << it.row() << " ";
    }
    return b;
}

// TODO: This should be named comp_blk_sums_and_sizes
// [[Rcpp::export]]
List comp_blk_sums_and_sizes(arma::sp_mat At, arma::uvec z, int Kcap, bool div_diag = true) {
    // z is a Kcap x 1 vector
    // At is a sparse n x n matrix

    int n = At.n_rows;
    // arma::mat lambda(Kcap, Kcap, arma::fill::zeros);

    arma::uvec zc = get_freq(z, Kcap); // zcounts

    arma::mat lambda = comp_blk_sums(At, z, Kcap);
    arma::umat NN = zc * zc.t() - arma::diagmat(zc);

    if (div_diag) {
        lambda.diag() /= 2; 
        NN.diag() /= 2;     // assumes that the diagonal of At is zero, 
                            // otherwise have to remove diag. first
    }
    
    return List::create( 
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("NN") = NN
    );
}



// An illustration of how to access only the zero elements of a sparse matrix
// [[Rcpp::export]]
void iter_over_sp_mat(arma::sp_mat At) {
    arma::sp_mat::const_iterator it     = At.begin();
    arma::sp_mat::const_iterator it_end = At.end();

    Rcout << std::setw(3) << "i"
        << std::setw(3) << "j"
        << std::setw(3) << "x"
        << std::endl;
    for(; it != it_end; ++it) {
        Rcout << std::setw(3) << it.row() + 1 
            << std::setw(3) << it.col() + 1
            << std::setw(3) << (*it) 
            << std::endl;
    }

}



