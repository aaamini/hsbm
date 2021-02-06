// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

#include "sampling.h"
#include "utils.h"
// #include "old_updates.h"

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

  // comment 15
  std::vector<arma::uvec> zb(tmax), zb_1based(tmax), gb(tmax);
  std::vector<arma::vec> u(tmax);
  std::vector<arma::vec> gamp(tmax), gam(tmax), v(tmax);
  arma::umat kb(tmax, Gcap);
  //std::vector <arma::uvec> kb;
  arma::vec pip(Kcap), pi(Kcap);
  arma::mat eta;

  // mat eta(Kcap, Kcap), aa(Kcap,Kcap), bb(Kcap,Kcap);

  for (int t=0; t < T; t++) {
    gb[t] = sample_int_vec(Gcap, n[t]);
    zb[t] = sample_int_vec(Kcap, n[t]);
    u[t] = Rcpp::runif(n[t]);

    // comment 17
    // kb[t] = arma::regspace<arma::uvec>(0, Gcap-1);
    // kb.row(t) = arma::regspace<arma::uvec>(0, Gcap-1).t();
    for (int g = 0; g < Gcap; g++){
      if (rand_init) {
        kb(t,g) = g;
      } else {
        kb(t,g) = 1;
      }
    }
    // v[t] = Rcpp::runif(Gcap);
  }

  // Rcout << "--- 1 --- \n";

  int itr = 0;
  bool CONVERGED = false;

  //NumericMatrix time_res(niter,10);
  NumericVector timer_res(10);


  while (!CONVERGED && itr < niter) {

    Timer timer;
    timer.step("start");

    //update gamma
    for (int t = 0; t < tmax; t++) {
      // comment 1
      gamp[t] =  gem_gibbs_update(gb[t], Gcap, gam0);
      gam[t] = stick_break(gamp[t]);

      //comment 10
    } // t

    timer.step("gamma");

    //Rcout << "--- 2 --- \n";

    // update pi
    // comment 2
    pip = gem_gibbs_update(kb.as_col(), Kcap, beta0);
    pi = stick_break( pip );

    timer.step("pi");

    //Rcout << "--- 3 --- \n";

    // comment 18

    // update eta
    List out = update_eta(A, n, zb, alpha_eta, beta_eta, Kcap);
    arma::mat aa = out["aa"];
    arma::mat bb = out["bb"];
    arma::mat temp_eta = out["eta"];
    eta = temp_eta;

    timer.step("eta");

    //Rcout << "--- 4 --- \n";


    // update gb
    for (int t=0; t < T; t++) {
      // mat At = A[t];
      arma::sp_mat At = A[t];
      // compute aok(g,l) = a(k_{tg},l) and bok(g,l) = b(k_{tg},l) which also depend on "t"
      arma::mat aok(Gcap, Kcap), bok(Gcap, Kcap);
      for (int g = 0; g < Gcap; g++){
        for (int el = 0; el < Kcap; el++) {
          aok(g,el) = aa(kb(t,g),el);
          bok(g,el) = bb(kb(t,g),el);
        } //el
      } //g

      //Rcout << "--- 5 --- \n";

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

    //Rcout << "---6 --- \n";
    timer.step("gb");

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


    // Rcpp::print( wrap(kb) );

    timer.step("kb");

    // comment 13
    //udpdate zb
    for (int t = 0; t < T; t++) {
      for (int i = 0; i < n[t]; i++){
        //zb(t,i) = kb(t, gb(t,i));
        zb[t](i) = kb(t, gb[t](i));
      }
      zb_1based[t] = zb[t] + 1;  // +1 is to put the labels on 1-based indexing
    }

    //Rcout << "--- 8 --- \n";

    timer.step("zb");

    // zb_list.push_back(clone(wrap(zb)));
    zb_list.push_back(zb_1based); 


    NumericVector temp_res(timer);
    timer_res += temp_res;

    if (verb) print_progress(itr, niter);

    itr++;

  }
  if (verb) print_progress(itr, niter);

  return Rcpp::List::create( Rcpp::Named("zb_list") = zb_list,
                             Rcpp::Named("gb") = gb,
                             Rcpp::Named("kb") = kb,
                             Rcpp::Named("eta") = eta,
                             Rcpp::Named("times") = timer_res / itr);
  // comment X

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





// // update kb
// arma::umat update_kb(arma::umat kb, std::vector<arma::uvec> gb, // arma::umat gb,
//                 List A, NumericVector n, arma::vec pi, std::vector<arma::vec> v, // arma::mat v, // IntegerMatrix K_all,
//                 arma::mat aa, arma::mat bb, int Kcap, int Gcap) {

//     int T = A.length();
//     for (int t=0; t < T; t++) {
//         arma::sp_mat At = A[t];

//         // List out = comp_blk_sums_and_sizes(At, gb.row(t).t(), Gcap, false);
//         List out = comp_blk_sums_and_sizes(At, gb[t], Gcap, false);
//         arma::mat  xit = out["lambda"]; // Gcap x Gcap
//         arma::umat OOt = out["NN"];  // Gcap x Gcap

//         // compute xi_{tgg} and O_{tgg}
//         arma::vec diag_xi = xit.diag()/2; // Gcap x 1   old code did not divide by 2
//         arma::uvec diag_OO = OOt.diag()/2; // Gcap x 1  old code did not divide by 2

//         for (int g = 0; g < Gcap; g++) {
//             // compute zeta_{tg} = (zeta_{tgl}) and RR_{tg} = (R_{tgl})
//             arma::vec zeta(Kcap);
//             arma::uvec RR(Kcap);

//             zeta = fast_agg(xit.row(g).t(), kb.row(t).t(), Kcap); // zeta[k_{g'}] += xit_{tgg'} over g'
//             zeta(kb(t,g)) -= xit(g,g); // zeta[k_{g}] -= xit_{tgg} effectivly excluding g' == g from the sum

//             RR = fast_agg_u(OOt.row(g).t(), kb.row(t).t(), Kcap);
//             RR(kb(t,g)) -= OOt(g,g);

//             // Sample k_{tg}
//             //int K = find_tunc( pi, v(t,g) );
//             int K = find_tunc( pi, v[t](g) );
//             // int K = K_all(t,g); // need to truncate "a" and "b" to K_{tg} number of rows;
//                                 // this is where slice sampling comes into play


//             kb(t,g) = sample_index(  safe_exp(diag_xi(g)*arma::diagvec(aa) +
//                                               diag_OO(g)*arma::diagvec(bb) +
//                                               aa*zeta + bb*RR + log(pi)) );

//             // kb(t,g) = sample_index(  safe_exp(diag_xi(g)*arma::diagvec(aa.rows(0,K)) +
//             //                                   diag_OO(g)*arma::diagvec(bb.rows(0,K)) +
//             //                                   aa.rows(0,K)*zeta + bb.rows(0,K)*RR) );


//             // kb(t,g) = sample_klabels(aa.rows(0,K), bb.rows(0,K), zeta, RR,  diag_xi(g), diag_OO(g));
//         } // g
//     } //t

//     return kb;
//     // return List::create(
//     //             Named("xit") = xit, Named("OOt") = OOt,
//     //             Named("diag_xi") = diag_xi, Named("diag_OO") = diag_OO,
//     //             Named("zeta") = zeta, Named("RR") = RR
//     //             );

// }

// arma::umat update_gb(arma::umat kb, arma::umat gb, arma::umat zb,
//                 List A, NumericVector n, IntegerMatrix G_all,
//                 arma::mat aa, arma::mat bb, int Kcap, int Gcap) {
//     // update gb
//     int T = A.length();

//     for (int t=0; t < T; t++) {
//       // mat At = A[t];
//       arma::sp_mat At = A[t];
//       // compute aok(g,l) = a(k_{tg},l) and bok(g,l) = b(k_{tg},l) which also depend on "t"
//       arma::mat aok(Gcap, Kcap), bok(Gcap, Kcap);
//       for (int g = 0; g < Gcap; g++){
//         for (int el = 0; el < Kcap; el++) {
//           aok(g,el) = aa(kb(t,g),el);
//           bok(g,el) = bb(kb(t,g),el);
//         } //el
//       } //g

//       arma::mat tau = sp_compress_col(At, zb.row(t).t(), Kcap);
//       arma::umat mm = get_freq_minus_self(zb.row(t).t(), Kcap);

//       for (int i = 0; i < n[t]; i++){

//         // comment 3

//         int G = G_all(t,i); // will truncate to G_{ti}; this is where slice sampling comes into play
//         // if (G+1  >= Gcap) Rcout << G;
//         // Rcout << aok.n_rows << " " << G << "\n";
//         gb(t,i) = sample_glabels( aok.rows(0, G), bok.rows(0, G), tau.row(i).t(), mm.row(i).t() );
//       }//i
//     } //t

//     return gb;
// }


// gb = update_gb(kb, gb, zb, A, n, G_all, aa, bb, Kcap, Gcap);
// kb_list.push_back(kb);
// gb_list.push_back(gb);
// pi_list.push_back(pi);
// gam_list.push_back(gam);
// aa_list.push_back(aa);
// bb_list.push_back(bb);

 //  Rcpp::Named("gb_list") = gb_list,
//  Rcpp::Named("kb_list") = kb_list,
//  Rcpp::Named("gam_list") = gam_list,
//  Rcpp::Named("pi_list") = pi_list,
//  Rcpp::Named("aa_list") = aa_list,
//  Rcpp::Named("bb_list") = bb_list,

//kb = update_kb_old(kb, gb, A, n, K_all, aa, bb, Kcap, Gcap);
// kb = update_kb( kb, Rcpp::as<arma::umat>(gb), A, n, K_all, aa, bb, Kcap, Gcap);
// kb = update_kb(kb, gb, A, n, K_all, aa, bb, Kcap, Gcap);

// if (itr > 100) {
//   kb = update_kb(kb, gb, A, n, pi, v, aa, bb, Kcap, Gcap);
// }

// List out = update_eta_old(A, n, zb, alpha_eta, beta_eta, Kcap);
// List out = update_eta(A, n, Rcpp::as<arma::umat>(zb), alpha_eta, beta_eta, Kcap);


// comment 18
    // for (int t=0; t < T; t++) {
    //   for (int g = 0; g < Gcap; g++){
    //     K_all(t,g) = find_tunc( pi, v(t,g) );
    //   } //g
    // } //t

// comment 17
 // for (int i=0; i < n[t]; i++){
    //   if (rand_init) {
    //     // gb(t, i) = sample_int(Gcap);
    //     gb[t](i) = sample_int(Gcap);
    //     // zb(t, i) = sample_int(Kcap);
    //     zb[t](i) = sample_int(Kcap);
    //   } else {
    //     // gb(t,i) = 2; // sample_int(Gcap); // sample_int(10);
    //     gb[t](i) = 2; // sample_int(Gcap); // sample_int(10);
    //     // zb(t,i) = 1;
    //     zb[t](i) = 1;
    //   }
    //   // u(t,i) = R::runif(0,1);
    //   u[t](i) = R::runif(0,1);
    // }

// comment 1
// arma::uvec count1 = get_freq(gb.row(t).t(), Gcap);
      // arma::uvec count2 = get_up_freq(count1);
      // // for (int g = 0; g < Gcap; g++){
      // //   int count1=0, count2=0;
      // //   for (int i = 0; i < nmax; i++){
      // //     if (gb(t,i) == g) {
      // //       count1++;
      // //     } else if (gb(t,i) > g) {
      // //       count2++;
      // //     }
      // //   }
      // //   gamp(t,g) = R::rbeta(count1 + 1, count2 + gam0);

      // // } // g

      // gamp.row(t) = rbeta_vec(arma::conv_to<arma::vec>::from(count1) + 1,
      //           arma::conv_to<arma::vec>::from(count2) + gam0).t();


// comment 2
    // for (int k=0; k < Kcap; k++) {
    //   int count1=0, count2=0;
    //   for (int t=0; t < T; t++) {
    //     for (int g = 0; g < Gcap; g++){
    //       if (kb(t,g) == k) {
    //         count1++;
    //       } else if (kb(t,g) > k) {
    //         count2++;
    //       }
    //     } // g
    //   } // t
    //   pip(k) = R::rbeta(count1 + 1, count2 + beta0);
    // } // k

//comment 3
    // uvec tau(Kcap), mm(Kcap); // these are in fact tau_{ti} and m_{ti}
    // // compute coordinates tau_{til} and m_{til} for all "l"
    // for (int el = 0; el < Kcap; el++) {
    //   tau(el) = 0;
    //   mm(el) = 0;
    //   for (int j = 0; j < n[t]; j++) {
    //     if ( j != i && zb(t,j) == el ) {
    //       tau(el) += At(i,j);
    //       mm(el)++;
    //     }
    //   } // j
    // } // el

    // if (t == 0 && i == 10) { // print diagnostics
    //   //Rcpp::print(wrap(tau.t()));
    //   // Rcpp::print(wrap(arma::join_cols(arma::conv_to<arma::rowvec>::from(tau), tau2.row(i))));
    //   Rcpp::print(wrap(arma::join_cols(arma::conv_to<arma::urowvec>::from(mm), mm2.row(i))));
    //   // Rcout << At.row(i);
    //   // Rcpp::print(wrap(zb.row(t)));
    //   Rcout << std::endl;
    // }

// comment X
// return List::create(
  //      _["gb_list"]= gb_list,
  //      _["kb_list"] = kb_list,
  //      _["zb_list"] = zb_list,
  //      _["pi"] = pi_list,
  //      _["gam"] = gam_list,
  //      _["eta"] = eta_list
  // );


// comment 10
 // gamp.row(t) = gem_gibbs_update(gb.row(t).t(), Gcap, gam0).t();
      // gam.row(t) = stick_break(gamp.row(t).t()).t();

      // // gam(t,span::all) = stick_break( gamp(t,span::all).t() ).t();

      // // for (int i=0; i < n[t]; i++){
      // //   // G_all(t,i) = find_tunc( gam(t,span::all).t(), u(t,i) );
      // //   G_all(t,i) = find_tunc( gam.row(t).t(), u(t,i) );
      // // } // i

// comment 12
// int gti_supp = find_tunc( gam[t], u[t](i) );
          // int gti_supp = find_tunc( gam.row(t).t(), u(t,i) );
          // //double uti = (itr == 0) ?  R::runif(0,1) : R::runif(0, gam(t, gb(t,i)));
          // //int gti_supp = find_tunc(gam.row(t).t(), uti);

          // gb(t,i) = sample_glabels( aok.rows(0, gti_supp), bok.rows(0, gti_supp), taui, mmi ); //this is where slice sampling comes into play
          // zb(t,i) = kb(t, gb(t,i)); // update z


          // gb[t](i) = sample_index( safe_exp( aok.rows(0, gti_supp) * taui +  bok.rows(0, gti_supp) * mmi ) );


// comment 13
// // update u
    // for (int t = 0; t < T; t++) {
    //   for (int i = 0; i < n[t]; i++){
    //     //u(t,i) = R::runif(0, gam(t, gb(t,i)) );
    //     u[t](i) = R::runif(0, gam[t](gb[t](i)) );
    //   }
    // }

    //Rcout << "--- 7 --- \n";

    // timer.step("u");

    // //update v
    // for (int t = 0; t < T; t++) {
    //   for (int g = 0; g < Gcap; g++){
    //     //v(t,g) = R::runif(0, pi( kb(t,g) ) );
    //     v[t](g) = R::runif(0, pi( kb(t,g) ) );
    //   }
    // }

    // timer.step("v");

// comment 14
// int gti_supp = find_tunc( gam.row(t).t(), u(t,i) );
          // // double uti = (itr == 0) ?  R::runif(0,1) : R::runif(0, gam(t, gb(t,i)));
          // // int gti_supp = find_tunc(gam.row(t).t(), uti);
          // gb(t,i) = sample_glabels( aok.rows(0, gti_supp), bok.rows(0, gti_supp), tau.row(i).t(), mm.row(i).t() );
          // zb(t,i) = kb(t, gb(t,i)); // update z

          // int gti_supp = find_tunc( gam[t], u[t](i) );
          // double uti = (itr == 0) ?  R::runif(0,1) : R::runif(0, gam(t, gb(t,i)));
          // int gti_supp = find_tunc(gam.row(t).t(), uti);

// comment 15
  // std::vector< std::vector<arma::uvec> > zb_list, gb_list;
  // std::vector< arma::umat > kb_list;
  // std::vector<IntegerMatrix> zb_list;
  // // std::vector<IntegerMatrix> kb_list;
  // // std::vector<IntegerMatrix> gb_list;
  // std::vector<arma::vec> pi_list;
  // std::vector<std::vector<arma::vec>> gam_list;
  // std::vector<arma::mat> aa_list, bb_list;

  // variables storing the current state of the chain
  // NumericMatrix u(T,nmax);
  // IntegerMatrix gb(T,nmax), zb(T,nmax);
  // arma::umat gb(T,nmax), zb(T,nmax);



  // IntegerMatrix G_all(T,nmax);
  // IntegerMatrix K_all(T,Kcap);
  // IntegerMatrix kb(T,Gcap);
  // arma::umat kb(T,Gcap);


// List update_eta_half(arma::sp_mat At, arma::uvec z, int Kcap) {
//     // update eta
//     // k != el
//     int n = At.n_rows;
//     arma::mat lambda(Kcap, Kcap, arma::fill::zeros);
//     arma::mat NN(Kcap, Kcap, arma::fill::zeros);
//     for (int k = 0; k < Kcap; k++){
//       for (int el = k+1; el < Kcap; el++) {

//           for (int i = 0; i < n; i++){
//             for (int j = 0; j < n; j++){
//               if (i != j && z(i) == k && z(j) == el) {
//                 lambda(k, el) += At(i,j);
//                 NN(k, el)++;
//               }
//             } // j
//           } // i

//         lambda(el, k) = lambda(k, el);
//         NN(el, k) = NN(k, el);
//         // eta(k,el) = R::rbeta(lambda + alpha_eta, NN - lambda + beta_eta);
//         // eta(el,k) = eta(k,el);
//         // aa(k,el) = log( eta(k,el)/(1-eta(k,el)) + 1e-11);
//         // aa(el,k) = aa(k,el);
//         // bb(k,el) = log( 1-eta(k,el) + 1e-11);
//         // bb(el,k) = bb(k,el);
//       } // el
//     } // k

//     for (int k = 0; k < Kcap; k++) {
//         for (int i = 0; i < n; i++){
//             for (int j = i+1; j < n; j++){
//                 if (z(i) == k && z(j) == k) {
//                     lambda(k, k) += At(i,j);
//                     NN(k, k)++;
//                 } //if
//             } // j
//         } // i
//     }

//     return List::create( Rcpp::Named("lambda") = lambda,
//                          Rcpp::Named("NN") = NN);
// }
