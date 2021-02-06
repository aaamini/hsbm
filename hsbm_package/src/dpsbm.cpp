// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

#include "sampling.h" 
#include "utils.h" 

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::umat fit_dpsbm(arma::sp_mat & A, 
                const double gam0 = 1,  
                const double alpha = 1, const double beta = 1,
                const int niter = 50, const int Zcap = 20, 
                const bool verb = true, const bool slice = false) {

    int n = A.n_rows;

    arma::umat z_hist(n, niter, arma::fill::zeros);
    // arma::uvec z(n, arma::fill::ones); 
    arma::uvec z = sample_int_vec(Zcap, n);
    arma::mat B(Zcap, Zcap, arma::fill::zeros);   
    // arma::vec pri(Zcap); 
    arma::vec gamp(Zcap), gam(Zcap);
    arma::vec unif(n, arma::fill::randu);

    for (int iter = 0; iter < niter; iter++) {
        List out = comp_blk_sums_and_sizes(A, z, Zcap);

        arma::mat lambda = out["lambda"];
        arma::umat NN = out["NN"]; 
        B = symmat_rbeta(lambda + alpha, NN - lambda + beta);
        arma::mat uu = log(B/(1-B) + 1e-11);
        arma::mat vv = log(1-B + 1e-11);

        gamp = gem_gibbs_update(z, Zcap, gam0);

        // arma::uvec count1 = get_freq(z, Zcap);
        // arma::uvec count2 = get_up_freq(count1);

        // gamp = rbeta_vec(arma::conv_to<arma::vec>::from(count1) + 1, 
        //                  arma::conv_to<arma::vec>::from(count2) + gam0);
        gam = stick_break(gamp);
        
        arma::vec Hi(Zcap, arma::fill::zeros);

        for (int i = 0; i < n; i++) {
            arma::vec taui = sp_single_col_compress(A, i, z, Zcap); // Zcap x 1 vector
            arma::uvec mmi = get_freq(z, Zcap);
            mmi(z(i))--;

            if (slice) {
              int zi_support = find_tunc( gam, R::runif(0, gam(z(i))) );
              Hi = uu.rows(0, zi_support) * taui + vv.rows(0, zi_support) * mmi;
              // Rcpp::print(wrap(zi_support));
            } else {
              int zi_support = Zcap-1;
              Hi = uu.rows(0, zi_support) * taui + vv.rows(0, zi_support) * mmi + log(gam);            
            }
            
            z(i) = sample_index(safe_exp(Hi));
        }
        z_hist.col(iter) = z + 1; // +1 is to put the labels on 1-based indexing
    }
    return z_hist;

}

//' @export
// [[Rcpp::export]]
List fit_mult_dpsbm(List A, 
                const double gam0 = 1,  
                const double alpha = 1, const double beta = 1,
                const int niter = 50, const int Zcap = 20, 
                const bool verb = true, const bool slice = false) {

  std::vector< std::vector<arma::uvec> > zb_list(niter);

  int tmax = A.length();
  // arma::uvec n(tmax);

  for (int t=0; t < tmax; t++) {
    // mat temp = A[t];
    arma::sp_mat temp = A[t];
    arma::umat zb_mat = fit_dpsbm(temp, gam0, alpha, beta, niter, Zcap, verb, slice);
    for (int it=0; it < niter; it++) {
      zb_list[it].push_back(zb_mat.col(it));
    }
  }
  return  List::create(Named("zb_list") = zb_list);
}

/* Older implementation

// Single layer slice sampler
// [[Rcpp::export]]
arma::umat dpsbm_slice_sampler_fast(sp_mat A, double gam0=1,  
                double alpha_eta = 1, double beta_eta = 1,
                int ITRmax=50, int Zcap=20, bool verb = true) {
  
  
  int n = A.n_rows;
  
  // lists to store the history of the chain
  // std::vector<IntegerVector> zb_list;
  //std::vector<mat> gam_list;
  //std::vector<mat> eta_list;
  List gam_list(ITRmax);
  List eta_list(ITRmax);
  
  // variables storing the current state of the chain
  arma::vec u(n);
  arma::uvec zb(n);
  uvec Z_all(n);
  arma::umat zb_hist(n, ITRmax);
  
  arma::vec gamp(Zcap), gam(Zcap);
  arma::mat eta(Zcap, Zcap), aa(Zcap,Zcap), bb(Zcap,Zcap);
  
  for (int i=0; i < n; i++){
    zb(i) = 1;
    u(i) = R::runif(0,1);
  }
    
  int itr = 0;
  bool CONVERGED = false;
  
  while (!CONVERGED && itr < ITRmax) {
    
    //update gamma
    for (int z = 0; z < Zcap; z++){
      int count1=0, count2=0;
      for (int i = 0; i < n; i++){
        if (zb(i) == z) {
          count1++;
        } else if (zb(i) > z) {
          count2++;
        }
      } //i
      gamp(z) = R::rbeta(count1 + 1, count2 + gam0);
      
    } // z
    gam = stick_break( gamp );

    // Rcout << 1;
    
    for (int i=0; i < n; i++){
      Z_all(i) = find_tunc( gam, u(i) );
    } // i 
   
    // Rcout << 2;

    // update eta
    List out = comp_blk_sums_and_sizes(A, zb, Zcap);

    arma::mat lambda = out["lambda"];
    arma::umat NN = out["NN"];    
    
    eta = symmat_rbeta(lambda + alpha_eta, NN - lambda + beta_eta);
    aa = log(eta/(1-eta) + 1e-11);
    bb = log(1-eta + 1e-11);

    // Rcout << 3;

    // update zb
    arma::mat tau = sp_compress_col(A, zb , Zcap);
    arma::umat mm = get_freq_minus_self(zb, Zcap);

    // Rcout << 4;       

    // Rcpp::print(wrap(mm));

    for (int i = 0; i < n; i++){
      int Z = Z_all(i); 
      zb(i) = sample_glabels( aa.rows(0,Z), bb.rows(0,Z), tau.row(i).t(), mm.row(i).t());
    } //i

    // Rcout << 5;
     
    // update u
    for (int i = 0; i < n; i++){
      u(i) = R::runif(0, gam(zb(i)) );
    }
    

    // Rcout << 6;

    // zb_list.push_back(clone(wrap(zb)));
    zb_hist.col(itr) = zb;

    // gam_list[itr] = gam;
    // eta_list[itr] = eta;
      

    if (verb) print_progress(itr, ITRmax);
    itr++;
  
  }
  if (verb) print_progress(itr, ITRmax);
  
  // return Rcpp::List::create( Rcpp::Named("zb_list") = zb_list );
  // return List::create(
  //      _["zb_list"] = zb_list,
  //      _["gam"] = gam_list,
  //      _["eta"] = eta_list
  // );
  return zb_hist;
  
} // dpsbm_slice_sampler_fast


*/