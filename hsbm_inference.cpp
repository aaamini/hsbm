// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include<random>
#include <Rcpp/Benchmark/Timer.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
//using namespace arma;
using arma::sp_mat;
using arma::randu;
using arma::mat;
using arma::vec;
using arma::uvec;
using arma::cube;
using arma::rowvec;
using arma::zeros;
using arma::ones;
using arma::span;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// // [[Rcpp::export]]
// int sample_int(int N) {
//   Rcpp::IntegerVector pool = Rcpp::seq(0, N-1);
//   std::random_shuffle(pool.begin(), pool.end());
//   return pool[0];
// } 

// [[Rcpp::export]]
int sample_int(int N) {
  Rcpp::IntegerVector pool = Rcpp::seq(0, N-1);
  return RcppArmadillo::sample(pool,1,TRUE)[0];
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
int my_sampler(arma::vec prob_vec) {
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::vector<double> weights = arma::conv_to< std::vector<double>>::from(prob_vec);
  
  std::discrete_distribution<> d(weights.begin(), weights.end());
  return( d(gen) );
}

// [[Rcpp::export]]
int sample_index(vec prob) {
  prob = prob / sum(prob);
  return( as_scalar(find(randu() < cumsum(prob),1)) );
}

// [[Rcpp::export]]
int sample_labels(mat a, mat b, uvec zeta, uvec R) {
  vec temp = a * zeta + b * R; // computes sum_l a_{k l} zeta_l + b_{k l} R_l
  temp = temp - max(temp)*ones(size(temp));
  vec prob_vec = exp(temp) + 1e-11;

  return sample_index( prob_vec );
  //return my_sampler( prob_vec );
}

// [[Rcpp::export]]
int sample_labels2(mat a, mat b, uvec zeta, uvec R, int xi, int OO) {
  vec temp = a * zeta + b * R; // computes sum_l a_{k l} zeta_l + b_{k l} R_l
  vec allones = ones(size(temp));
  temp += diagvec(a) % (xi*allones) + diagvec(b) % (OO*allones); // compute a_{kk} xi + b_{kk} OO
  temp = temp - max(temp)*allones;
  vec prob_vec = exp(temp) + 1e-11;

  return sample_index( prob_vec );
  //return my_sampler( prob_vec );
}


// [[Rcpp::export]]
List hsbm_inferC(List A, double beta0=3, double gam0=1,  
                double alpha_eta = 1, double beta_eta = 1,
                int ITRmax=50, int Kcap=20, int Gcap=20) {
  
  int T = A.length();
  NumericVector n(T);
  
  for (int t=0; t < T; t++) {
    // mat temp = A[t];
    sp_mat temp = A[t];
    n(t) = temp.n_rows; 
  }
  
  int nmax = max(n);
  
  // lists to store the history of the chain
  std::vector<IntegerMatrix> zb_list;
  std::vector<IntegerMatrix> kb_list;
  std::vector<IntegerMatrix> gb_list;
  std::vector<vec> pi_list;
  std::vector<mat> gam_list;
  std::vector<mat> eta_list;
  
  
  // variables storing the current state of the chain
  NumericMatrix u(T,nmax); 
  IntegerMatrix gb(T,nmax), zb(T,nmax);
  IntegerMatrix G_all(T,nmax);
  IntegerMatrix K_all(T,Kcap);
  
  mat v(T,Gcap), gamp(T,Gcap), gam(T,Gcap);
  IntegerMatrix kb(T,Gcap);
  vec pip(Kcap), pi(Kcap);
  
  mat eta(Kcap, Kcap), aa(Kcap,Kcap), bb(Kcap,Kcap);
  
  for (int t=0; t < T; t++) {
    for (int i=0; i < n[t]; i++){
      
      gb(t,i) = 2;// sample_int(10);
      
      zb(t,i) = 1;
      u(t,i) = R::runif(0,1);
    }
    for (int g = 0; g < Gcap; g++){
      kb(t,g) = 1;
      v(t,g) = R::runif(0,1);
    }
    //u(j, _) = runif(nmax);
  }
  
  int itr = 0;
  bool CONVERGED = false;
  
  //NumericMatrix time_res(ITRmax,10);
  // NumericVector time_res(10);
  vec time_res(10);
  
  while (!CONVERGED && itr < ITRmax) {
    
    Timer timer;
    
    timer.step("start");
    
    //update gamma
    for (int t = 0; t < T; t++) {
      for (int g = 0; g < Gcap; g++){
        
        int count1=0, count2=0;
        for (int i = 0; i < nmax; i++){
          if (gb(t,i) == g) {
            count1++;
          } else if (gb(t,i) > g) {
            count2++;
          }
        } 
        gamp(t,g) = R::rbeta(count1 + 1, count2 + gam0);
        
      } // g
      gam(t,span::all) = stick_break( gamp(t,span::all).t() ).t();
      
      for (int i=0; i < n[t]; i++){
        G_all(t,i) = find_tunc( gam(t,span::all).t(), u(t,i) );
      } // i 
    } // t
    
    timer.step("gamma");
    
    // update pi
    for (int k=0; k < Kcap; k++) {
      int count1=0, count2=0;
      for (int t=0; t < T; t++) {
        for (int g = 0; g < Gcap; g++){
          if (kb(t,g) == k) {
            count1++;
          } else if (kb(t,g) > k) {
            count2++;
          }
        } // g
      } // t
      pip(k) = R::rbeta(count1 + 1, count2 + beta0);
    } // k
    
    pi = stick_break( pip );
    
    timer.step("pi");
    
    for (int t=0; t < T; t++) {
      for (int g = 0; g < Gcap; g++){
        K_all(t,g) = find_tunc( pi, v(t,g) );
      } //g
    } //t

    timer.step("K_all");
    
    // update eta
    // k != el
    for (int k = 0; k < Kcap; k++){
      for (int el = k+1; el < Kcap; el++) {
        int lambda = 0;
        int NN = 0;
        for (int t = 0; t < T; t++) {
          //mat At = A[t];
          sp_mat At = A[t];
          for (int i = 0; i < n[t]; i++){
            for (int j = 0; j < n[t]; j++){
              if (i != j && zb(t,i) == k && zb(t,j) == el) {
                lambda += At(i,j);
                NN++;
              }
            } // j 
          } // i
        } // t
        
        eta(k,el) = R::rbeta(lambda + alpha_eta, NN - lambda + beta_eta);
        eta(el,k) = eta(k,el);
        aa(k,el) = log( eta(k,el)/(1-eta(k,el)) + 1e-11);
        aa(el,k) = aa(k,el);
        bb(k,el) = log( 1-eta(k,el) + 1e-11);
        bb(el,k) = bb(k,el);
      } // el
    } // k

    // k == el
    for (int k = 0; k < Kcap; k++) {
      int lambda = 0;
      int NN = 0;
      for (int t = 0; t < T; t++) {
        // mat At = A[t];
        sp_mat At = A[t];
        for (int i = 0; i < n[t]; i++){
          for (int j = i+1; j < n[t]; j++){
            if (zb(t,i) == k && zb(t,j) == k) {
              lambda += At(i,j);
              NN++;
            } //if 
          } // j
        } // i
      } // t
      eta(k,k) = R::rbeta(lambda + alpha_eta, NN - lambda + beta_eta);
      aa(k,k) = log( eta(k,k)/(1-eta(k,k)) + 1e-11);
      bb(k,k) = log( 1-eta(k,k) + 1e-11);
    } // k

    timer.step("eta");
    
    // update kb 
    for (int t=0; t < T; t++) {
      // mat At = A[t];
      sp_mat At = A[t];
      for (int g = 0; g < Gcap; g++) {
        // compute xi_{tgg} and O_{tgg}
        int xi = 0;
        int OO = 0;
        for (int i = 0; i < n[t]; i++) {
          for (int j = i+1; j < n[t]; j++) {
            if (i != j && gb(t,i) == g && gb(t,j) == g) {
              xi += At(i,j);
              OO++;
            }
          } //j   
        } //i  

        // compute zeta_{tg} = (zeta_{tgl}) and RR_{tg} = (R_{tgl})
        uvec zeta(Kcap), RR(Kcap);
        for (int el = 0; el < Kcap; el++) {
          zeta(el) = 0;
          RR(el) = 0;
          for (int gp = 0; gp < Gcap; gp++) {
            for (int i = 0; i < n[t]; i++) {
              for (int j = 0; j < n[t]; j++) {
                if (gp != g && i != j && gb(t,i) == g && gb(t,j) == gp && kb(t,gp) == el ) {
                  zeta(el) += At(i,j);
                  RR(el)++;
                }
              } // j 
            } // i
          } // gp  
        } // el
        
        // Sample k_{tg}  
        int K = K_all(t,g); // need to truncate "a" and "b" to K_{tg} number of rows; ; this is where slice sampling comes into play
        //Rcout << aa.n_rows << " " << K << "\n";
        kb(t,g) = sample_labels2( aa.rows(0,K), bb.rows(0,K), zeta, RR,  xi, OO);
        //kb(t,g) = sample_labels( aa.rows(0,K), bb.rows(0,K), zeta, RR);
      } // g
    } //t
   
    timer.step("kb");
   
    // update gb
    for (int t=0; t < T; t++) {
      // mat At = A[t];
      sp_mat At = A[t];
      // compute aok(g,l) = a(k_{tg},l) and bok(g,l) = b(k_{tg},l) which also depend on "t"
      mat aok(Gcap,Kcap), bok(Gcap,Kcap);
      for (int g = 0; g < Gcap; g++){
        for (int el = 0; el < Kcap; el++) {
          aok(g,el) = aa(kb(t,g),el);
          bok(g,el) = bb(kb(t,g),el);
        } //el
      } //g

      for (int i = 0; i < n[t]; i++){
        uvec tau(Kcap), mm(Kcap); // these are in fact tau_{ti} and m_{ti}
        // compute coordinates tau_{til} and m_{til} for all "l"
        for (int el = 0; el < Kcap; el++) {
          tau(el) = 0;
          mm(el) = 0;
          for (int j = 0; j < n[t]; j++) {
            if ( j != i && zb(t,j) == el ) {
              tau(el) += At(i,j);
              mm(el)++;
            }        
          } // j
        } // el
      
        int G = G_all(t,i); // will truncate to G_{ti}; this is where slice sampling comes into play
        // if (G+1  >= Gcap) Rcout << G;
        // Rcout << aok.n_rows << " " << G << "\n";
        gb(t,i) = sample_labels( aok.rows(0,G), bok.rows(0,G), tau, mm );
      }//i
    } //t
    
    timer.step("gb");
   
    // update u
    for (int t = 0; t < T; t++) {
      for (int i = 0; i < n[t]; i++){
      u(t,i) = R::runif(0, gam(t, gb(t,i)) );
      }
    }
    
    timer.step("u");
   
    //update v
    for (int t = 0; t < T; t++) {
      for (int g = 0; g < Gcap; g++){
        v(t,g) = R::runif(0, pi( kb(t,g) ) );
      }
    }
    
    timer.step("v");
   
    //udpdate zb
    for (int t = 0; t < T; t++) {
      for (int i = 0; i < n[t]; i++){
        zb(t,i) = kb(t, gb(t,i));
      }
    }
    
    timer.step("zb");
   
    //Zb.slice(itr) = zb;
    zb_list.push_back(clone(zb));
    // gb_list.push_back(clone(gb));
    // kb_list.push_back(clone(kb));
    // pi_list.push_back(clone(pi));
    // gam_list.push_back(clone(gam));
    // eta_list.push_back(clone(eta));
    
    NumericVector temp_res(timer);
    // time_res(itr, _) = temp_res;
    time_res += temp_res;
      
    if (itr % 10 == 0) Rcout << '*';
    if (itr % 100 == 0 && itr/100 > 1) Rcout << itr << '\n';
    
    itr++;
  
  }
  
  return Rcpp::List::create( Rcpp::Named("zb_list") = zb_list,
                             Rcpp::Named("times") = time_res);
  // return List::create(
  //      _["gb_list"]= gb_list,
  //      _["kb_list"] = kb_list,
  //      _["zb_list"] = zb_list,
  //      _["pi"] = pi_list,
  //      _["gam"] = gam_list,
  //      _["eta"] = eta_list
  // );
  
} // hsbm_infer



// Single layer slice sampler
// [[Rcpp::export]]
List dpsbm_inferC(sp_mat A, double gam0=1,  
                double alpha_eta = 1, double beta_eta = 1,
                int ITRmax=50, int Zcap=20) {
  
  
  int n = A.n_rows;
  
  // lists to store the history of the chain
  std::vector<IntegerVector> zb_list;
  //std::vector<mat> gam_list;
  //std::vector<mat> eta_list;
  List gam_list(ITRmax);
  List eta_list(ITRmax);
  
  // variables storing the current state of the chain
  vec u(n);
  IntegerVector zb(n);
  uvec Z_all(n);
  
  vec gamp(Zcap), gam(Zcap);
  mat eta(Zcap, Zcap), aa(Zcap,Zcap), bb(Zcap,Zcap);
  
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
    
    for (int i=0; i < n; i++){
      Z_all(i) = find_tunc( gam, u(i) );
    } // i 
   
    // update eta
    // k != el
    for (int k = 0; k < Zcap; k++){
      for (int el = k+1; el < Zcap; el++){

        int lambda = 0;
        int NN = 0;
        for (int i = 0; i < n; i++){
          for (int j = 0; j < n; j++){
            if (i != j && zb(i) == k && zb(j) == el) {
              lambda += A(i,j);
              NN++;
            }
          } // j 
        } // i
      
      
      eta(k,el) = R::rbeta(lambda + alpha_eta, NN - lambda + beta_eta);
      eta(el,k) = eta(k,el);
      aa(k,el) = log( eta(k,el)/(1-eta(k,el)) + 1e-11);
      aa(el,k) = aa(k,el);
      bb(k,el) = log( 1-eta(k,el) + 1e-11);
      bb(el,k) = bb(k,el);
      } // l
    } // k

    // k == el
    for (int k = 0; k < Zcap; k++){
        int lambda = 0;
        int NN = 0;
      
        for (int i = 0; i < n; i++){
          for (int j = i+1; j < n; j++){
            if (zb(i) == k && zb(j) == k) {
              lambda += A(i,j);
              NN++;
            }
          }
        }
        
        eta(k,k) = R::rbeta(lambda + alpha_eta, NN - lambda + beta_eta);
        aa(k,k) = log( eta(k,k)/(1-eta(k,k)) + 1e-11);
        bb(k,k) = log( 1-eta(k,k) + 1e-11);
    } // k

    // update zb
    for (int i = 0; i < n; i++) {
      uvec tau(Zcap), mm(Zcap); // these are in fact tau_{i} and m_{i}
      // compute coordinates tau_{til} and m_{til} for all "l"
      for (int el = 0; el < Zcap; el++) {
        tau(el) = 0;
        mm(el) = 0;
        for (int j = 0; j < n; j++) {
          if ( j != i && zb(j) == el ) {
            tau(el) += A(i,j);
            mm(el)++;
          }        
        } // j
      } // el
     
      int Z = Z_all(i); // will truncate to Z_{i}; this is where slice sampling comes into play

      zb(i) = sample_labels( aa.rows(0,Z), bb.rows(0,Z), tau, mm );
    } //i
     
    // update u
    for (int i = 0; i < n; i++){
      u(i) = R::runif(0, gam(zb(i)) );
    }
    
    zb_list.push_back(clone(zb));
    // gam_list.push_back(gam);
    // eta_list.push_back(eta);
    gam_list[itr] = gam;
    eta_list[itr] = eta;
      
    if (itr % 10 == 0) Rcout << '*';
    if (itr % 100 == 0 && itr/100 > 1) Rcout << itr << '\n';
    
    itr++;
  
  }
  
  // return Rcpp::List::create( Rcpp::Named("zb_list") = zb_list );
  return List::create(
       _["zb_list"] = zb_list,
       _["gam"] = gam_list,
       _["eta"] = eta_list
  );
  
} // dpsbm_infer





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# n <- 10000
# temp <- numeric(n)
# prob <- c(1,3,4,4)
# for (j in 1:n) {
#   temp[j] <- sample_index(prob)
# }
# table(temp)/n
# prob/sum(prob)

# res <- hsbm_infer(Amat)


  
*/
