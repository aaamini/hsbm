// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <random>

#include "utils.h"

using namespace Rcpp;
// using namespace arma;   // cuases Reference to RcppArmadillo is ambiguous

// [[Rcpp::export]]
int sample_int(int N) {  // This can perhaps be removed later
  Rcpp::IntegerVector pool = Rcpp::seq(0, N-1);
  return RcppArmadillo::sample(pool,1,TRUE)[0];
} 

// [[Rcpp::export]]
arma::uvec sample_int_vec(int N, int size) {
  Rcpp::IntegerVector pool = Rcpp::seq(0, N-1);
  return Rcpp::as<arma::uvec>(RcppArmadillo::sample(pool, size, TRUE));
} 


// [[Rcpp::export]]
arma::vec stick_break(arma::vec x) {
  int n = x.n_elem;
  arma::vec out = arma::ones(n+1);
  out(arma::span(1,n)) = arma::cumprod( 1 - x );
  return( out(arma::span(0,n-1)) %  x);
}

// // [[Rcpp::export]]
// arma::vec stick_break(arma::vec x) {
//   int n = x.n_elem;
//   arma::vec out = arma::ones(n+1);
//   out(arma::span(1,n)) = arma::cumprod( arma::ones<arma::vec>(n)-x );
//   return( out(arma::span(0,n-1)) %  x);
// }


// [[Rcpp::export]]
arma::mat symmat_rbeta(arma::mat alpha, arma::mat beta) {
    // Samples form a symmetric beta-distributed matrix, with independent entries
    // drawn as ~ Beta(alpha(i,j), beta(i,j))
    int K = alpha.n_rows;
    arma::mat samp(K, K, arma::fill::zeros);

    for (int i = 0; i < K; i++) {
        for (int j = i; j < K; j++) {
            samp(i,j) = R::rbeta(alpha(i,j), beta(i,j));
            samp(j,i) = samp(i,j);
        }
    }
    return samp;
}



// // [[Rcpp::export]]
// int find_tunc(arma::vec beta, double threshold) {
//   
//   int n = beta.n_elem;
//   arma::vec temp = arma::zeros(n+1);
//   temp(arma::span(1,n)) = cumsum(beta);
//   
//   int idx = 0;
//   for (int j = 0; j < n+1; j++) {
//     idx = n-j;
//     if (temp(idx) < 1-threshold) break;
//   }
//   if (idx > n-1) idx--; 
//   return ( idx );
// }

// [[Rcpp::export]]
arma::vec safe_exp(arma::vec log_prob) {
  log_prob -= arma::max(log_prob);
  return exp(log_prob) + 1e-11;
}


// [[Rcpp::export]]
int my_sampler(arma::vec prob_vec) {
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::vector<double> weights = arma::conv_to< std::vector<double>>::from(prob_vec);
  
  std::discrete_distribution<> d(weights.begin(), weights.end());
  return( d(gen) );
}

// [[Rcpp::export]]
int sample_index(arma::vec prob) {
  prob = prob / sum(prob);
  return( arma::as_scalar(arma::find(arma::randu() < arma::cumsum(prob),1)) );
}


// std::vector<arma::uvec> sample_markov_labels(const arma::uvec &zinit, const int nlayers, const int K) {
//   int n = zinit.n_elem;
//   std::vector<arma::uvec> z_list(nlayers);
//   z_list[0] = zinit;
//   for (int t = 1; t < nlayers; t++) {
//     for ()
//     z_list[t]
//   }
// }


// int sample_labels(arma::mat a, arma::mat b, arma::uvec zeta, arma::uvec R) {
//   arma::vec temp = a * zeta + b * R; // computes sum_l a_{k l} zeta_l + b_{k l} R_l
//   temp = temp - arma::max(temp)*arma::ones(arma::size(temp));
//   arma::vec prob_vec = exp(temp) + 1e-11;

//   return sample_index( prob_vec );
//   //return my_sampler( prob_vec );
// }



// int sample_labels2(arma::mat a, arma::mat b, 
//                   arma::uvec zeta, arma::uvec R, int xi, int OO) {
//   arma::vec temp = a * zeta + b * R; // computes sum_l a_{k l} zeta_l + b_{k l} R_l
//   arma::vec allones = arma::ones(arma::size(temp));
//   temp += arma::diagvec(a) % (xi*allones) + arma::diagvec(b) % (OO*allones); // compute a_{kk} xi + b_{kk} OO
//   temp = temp - arma::max(temp)*allones;
//   arma::vec prob_vec = exp(temp) + 1e-11;

//   return sample_index( prob_vec );
//   //return my_sampler( prob_vec );
// }

// [[Rcpp::export]]
 arma::vec rgamma_vec(arma::vec shape) {
   // scale = 1
  int n = shape.n_elem;
  arma::vec out(n);

  for (int i = 0; i < n; i++) {
    out(i) = R::rgamma(shape(i), 1);
  }

  return out;
 }

 // [[Rcpp::export]]
 arma::vec rbeta_vec(arma::vec alpha, arma::vec beta) {
   // scale = 1
  int n = alpha.n_elem;
  arma::vec out(n);

  for (int i = 0; i < n; i++) {
    out(i) = R::rbeta(alpha(i),beta(i));
  }

  return out;
 }

// [[Rcpp::export]]
arma::vec rdirichlet(arma::vec theta) {
  arma::vec out = rgamma_vec(theta);
  return out / arma::sum(out);
 }

// [[Rcpp::export]]
arma::vec gem_gibbs_update(const arma::uvec z, const int Zcap, 
                        const double concent_param) {
  arma::uvec count1 = get_freq(z, Zcap);
  arma::uvec count2 = get_up_freq(count1);

  return(
    rbeta_vec(arma::conv_to<arma::vec>::from(count1) + 1, 
             arma::conv_to<arma::vec>::from(count2) + concent_param)
  );
}



// int sample_klabels(arma::mat a, arma::mat b, 
//                   arma::vec zeta, arma::uvec R, double xi, int OO) {
//  // computes a_{kk} xi + b_{kk} OO + sum_l a_{k l} zeta_l + b_{k l} R_l
//   arma::vec prob_vec = safe_exp( xi*arma::diagvec(a) + OO*arma::diagvec(b) + a*zeta + b*R );
//   //return prob_vec;
//   return sample_index( prob_vec );
// }



// int sample_glabels(arma::mat a, arma::mat b, 
//                     arma::vec zeta, arma::uvec R) {
//   // computes sum_l a_{k l} zeta_l + b_{k l} R_l
//   arma::vec prob_vec = safe_exp( a * zeta + b * R ); 
//   return sample_index( prob_vec );
//   //return my_sampler( prob_vec );
// }


