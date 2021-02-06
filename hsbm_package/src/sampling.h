#ifndef __SAMPLING__
#define __SAMPLING__

#include <random>

arma::vec safe_exp(arma::vec log_prob);

int sample_int(int N);
arma::uvec sample_int_vec(int N, int size);

arma::vec rgamma_vec(arma::vec shape);  // scale = 1
arma::vec rdirichlet(arma::vec theta);
arma::vec rbeta_vec(arma::vec alpha, arma::vec beta);


arma::vec gem_gibbs_update(const arma::uvec z, const int Zcap, const double concent_param);
arma::vec stick_break(arma::vec x);

int find_tunc(arma::vec beta, double threshold);

int my_sampler(arma::vec prob_vec);

int sample_index(arma::vec prob);

arma::mat symmat_rbeta(arma::mat alpha, arma::mat beta);

int sample_labels(arma::mat a, arma::mat b, arma::uvec zeta, arma::uvec R);


int sample_labels2(arma::mat a, arma::mat b, 
                    arma::uvec zeta, arma::uvec R, int xi, int OO);

int sample_glabels(arma::mat a, arma::mat b, 
                    arma::vec zeta, arma::uvec R);
int sample_klabels(arma::mat a, arma::mat b, 
                    arma::vec zeta, arma::uvec R, double xi, int OO);

#endif // __SAMPLING__