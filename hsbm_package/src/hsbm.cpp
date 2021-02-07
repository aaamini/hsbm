// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

#include "sampling.h"
#include "utils.h"

using namespace Rcpp;


List update_eta(List Alist, NumericVector n, std::vector<arma::uvec> zb, // arma::umat zb,
                    double alpha_eta, double beta_eta, int Kcap);

// --- could export like this using Roxygen2:
//' @export
// [[Rcpp::export]]
List fit_hsbm(List A, double beta0=0.1, double gam0=0.5,
                double alpha_eta = 1, double beta_eta = 5,
                int niter=50, int Kcap=10, int Gcap=10, bool verb = true,
                const bool rand_init = true, const bool seq_g_update = true) {

  int T = A.length();
  int tmax = A.length();
  NumericVector n(T);

  for (int t=0; t < T; t++) {
    // mat temp = A[t];
    arma::sp_mat temp = A[t];
    n(t) = temp.n_rows;
  }

  int nmax = max(n);

  // lists to store the history of the chain
  std::vector< std::vector<arma::uvec> > zb_list;

  std::vector<arma::uvec> zb(tmax), zb_1based(tmax), gb(tmax);
  std::vector<arma::vec> u(tmax);
  std::vector<arma::vec> gamp(tmax), gam(tmax), v(tmax);
  arma::umat kb(tmax, Gcap);
  //std::vector <arma::uvec> kb;
  arma::vec pip(Kcap), pi(Kcap);
  arma::mat eta;

  for (int t=0; t < T; t++) {
    gb[t] = sample_int_vec(Gcap, n[t]);
    zb[t] = sample_int_vec(Kcap, n[t]);
    u[t] = Rcpp::runif(n[t]);
    // kb[t] = arma::regspace<arma::uvec>(0, Gcap-1);
    // kb.row(t) = arma::regspace<arma::uvec>(0, Gcap-1).t();
    for (int g = 0; g < Gcap; g++){
      if (rand_init) {
        kb(t,g) = g;
      } else {
        kb(t,g) = 1;
      }
    }
  }

  int itr = 0;
  bool CONVERGED = false;

  while (!CONVERGED && itr < niter) {
    //update gamma
    for (int t = 0; t < tmax; t++) {
      // comment 1
      gamp[t] =  gem_gibbs_update(gb[t], Gcap, gam0);
      gam[t] = stick_break(gamp[t]);
    } // t

    // update pi
    pip = gem_gibbs_update(kb.as_col(), Kcap, beta0);
    pi = stick_break( pip );

    // update eta
    List out = update_eta(A, n, zb, alpha_eta, beta_eta, Kcap);
    arma::mat aa = out["aa"];
    arma::mat bb = out["bb"];
    arma::mat temp_eta = out["eta"];
    eta = temp_eta;

    // update gb
    for (int t=0; t < T; t++) {
      arma::sp_mat At = A[t];
      // compute aok(g,l) = a(k_{tg},l) and bok(g,l) = b(k_{tg},l) which also depend on "t"
      arma::mat aok(Gcap, Kcap), bok(Gcap, Kcap);
      for (int g = 0; g < Gcap; g++){
        for (int el = 0; el < Kcap; el++) {
          aok(g,el) = aa(kb(t,g),el);
          bok(g,el) = bb(kb(t,g),el);
        } //el
      } //g

      if (seq_g_update) {
        for (int i = 0; i < n[t]; i++) {
          // comment 3
          // arma::uvec zz = zb.row(t).t();
          arma::uvec zz = zb[t];
          arma::vec taui = sp_single_col_compress(At, i, zz, Kcap);
          arma::uvec mmi = get_freq(zz, Kcap);
          mmi(zz(i))--;

          //comment 12
          gb[t](i) = sample_index(safe_exp(aok * taui +  bok * mmi + log(gam[t])));
          zb[t](i) = kb(t, gb[t](i)); // update z
        }//i

      } else {
        arma::mat tau = sp_compress_col(At, zb[t], Kcap);
        arma::umat mm = get_freq_minus_self(zb[t], Kcap);
        for (int i = 0; i < n[t]; i++){
          // comment 3
          // comment 14
          gb[t](i) = sample_index(safe_exp(aok *  tau.row(i).t() +  bok * mm.row(i).t() + log(gam[t])));
          zb[t](i) = kb(t, gb[t](i)); // update z
        }//i
      }

    } //t

    // update kb
    for (int t=0; t < T; t++) {
        arma::sp_mat At = A[t];

        // List out = comp_blk_sums_and_sizes(At, gb.row(t).t(), Gcap, false);
        List out = comp_blk_sums_and_sizes(At, gb[t], Gcap, false);
        arma::mat  xit = out["lambda"]; // Gcap x Gcap
        arma::umat OOt = out["NN"];  // Gcap x Gcap

        // compute xi_{tgg} and O_{tgg}
        arma::vec diag_xi = xit.diag()/2; // Gcap x 1   old code did not divide by 2
        arma::uvec diag_OO = OOt.diag()/2; // Gcap x 1  old code did not divide by 2

        for (int g = 0; g < Gcap; g++) {
            // compute zeta_{tg} = (zeta_{tgl}) and RR_{tg} = (R_{tgl})
            arma::vec zeta(Kcap);
            arma::uvec RR(Kcap);

            zeta = fast_agg(xit.row(g).t(), kb.row(t).t(), Kcap); // zeta[k_{g'}] += xit_{tgg'} over g'
            zeta(kb(t,g)) -= xit(g,g); // zeta[k_{g}] -= xit_{tgg} effectivly excluding g' == g from the sum

            RR = fast_agg_u(OOt.row(g).t(), kb.row(t).t(), Kcap);
            RR(kb(t,g)) -= OOt(g,g);

            kb(t,g) = sample_index(
              safe_exp(diag_xi(g)*arma::diagvec(aa) +
                      diag_OO(g)*arma::diagvec(bb) +
                      aa*zeta + bb*RR + log(pi)) );
        } // g
    } //t

    //udpdate zb
    for (int t = 0; t < T; t++) {
      for (int i = 0; i < n[t]; i++){
        zb[t](i) = kb(t, gb[t](i));
      }
      zb_1based[t] = zb[t] + 1;  // +1 is to put the labels on 1-based indexing
    }

    // zb_list.push_back(clone(wrap(zb)));
    zb_list.push_back(zb_1based);

    if (verb) print_progress(itr, niter);

    itr++;

  }
  if (verb) print_progress(itr, niter);

  return Rcpp::List::create( Rcpp::Named("zb_list") = zb_list,
                             Rcpp::Named("gb") = gb,
                             Rcpp::Named("kb") = kb,
                             Rcpp::Named("eta") = eta);

} // hsbm_infer


// update eta
// [[Rcpp::export]]
List update_eta(List Alist, NumericVector n, std::vector<arma::uvec> zb, // arma::umat zb,
                    double alpha_eta, double beta_eta, int Kcap) {
    int T = Alist.length();
    arma::mat lambda(Kcap, Kcap, arma::fill::zeros);
    arma::umat NN(Kcap, Kcap, arma::fill::zeros);

    for (int tt = 0; tt < T; tt++) {
        // List out = comp_blk_sums_and_sizes(Alist[tt], zb.row(tt).t(), Kcap);
        List out = comp_blk_sums_and_sizes(Alist[tt], zb[tt], Kcap);
        arma::mat temp_lambda = out["lambda"];
        arma::umat temp_NN = out["NN"];
        lambda += temp_lambda;
        NN += temp_NN;
    }
    // List out = update_lambda_N(Alist[1], zb.row(1).t(), Kcap);
    // arma::mat lambda = out["lambda"];
    // arma::umat NN = out["NN"];

    arma::mat eta = symmat_rbeta(lambda + alpha_eta, NN - lambda + beta_eta);
    arma::mat aa = log(eta/(1-eta) + 1e-11);
    arma::mat bb = log(1-eta + 1e-11);

    return List::create(
        Rcpp::Named("eta") = eta,
        Rcpp::Named("aa") = aa,
        Rcpp::Named("bb") = bb,
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("NN") = NN
    );
}
