//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include<RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;




// function to store correlations from matrix to vector
// [[Rcpp::export]]

colvec corvec(mat V){
  int l = V.n_rows;
  colvec x(1);
  
  for(int j=0; j<l-1; j++){
    
    colvec temp = trans(V(j, span(j+1,l-1)));
    x = join_cols(x,temp);
    
  }
  
  x.shed_row(0);
  
  return x;
}



// [[Rcpp::export]]
NumericVector realRoot(double a, double b, double c){
  
  NumericVector out(2);
  
  double D = pow(b,2)-4*a*c;
  
  if(D<=0){
    cout << "Error: no real roots";
    return false;
  }
  
  out[0] = (-b-sqrt(D))/(2*a);
  out[1] = (-b+sqrt(D))/(2*a);
  
  return out; 
}


// draw from multivariate normal distribution

// [[Rcpp::export]]
mat mvrnorm(int n, colvec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return (repmat(mu, 1, n).t() + Y * chol(sigma));
}


//=============================================================
// check positive definite matrix
//=============================================================


// [[Rcpp::export]]
bool ispd(mat X, double tol=1e-8){
  colvec eigenval = eig_sym(X);
  
  int n = X.n_rows;
  //double maxeig = max(eigenval);
  //if(!all(eigenval>=-tol*abs(maxeig))){
  //  return false;
  //}
  colvec abseig = abs(eigenval);
  
  for(int i=0;i<n;i++){
    if(abseig(i) < tol){
      eigenval(i) = 0;
    }
  }
  
  if(any(eigenval<=0)){
    return false;
  }else{
    return true;
  }
  
}




//[[Rcpp::export]]
mat sigMatrix(colvec sig2draw, int k, std::string sigstructure){
  
  //int nsig = sig2draw.n_elem;
  mat Sigma(k,k);
  colvec sig2(k);
  
  if(sigstructure=="equalvar"){
    sig2.fill(as_scalar(sqrt(sig2draw))); 
  }
  
  if(sigstructure=="unequalvar"){
    sig2 = sqrt(sig2draw);   
  }
  
  Sigma = diagmat(sig2);
  
  return Sigma;
}



// Generate Initail pars

//[[Rcpp::export]]
List initialGen(colvec& y, mat& X, int k, 
                std::string sigstructure,
                std::string Rstructure){
  
  //RNGScope scope;
  int n = X.n_rows, p = X.n_cols;
  
  //initial value for betadraw
  colvec betadraw = inv(trans(X)*X)*(trans(X)*y);
  
  colvec res = y-X*betadraw;
  
  double sig2init = as_scalar((trans(res)*res)/(n-p));
  
  colvec sig2draw, rhodraw;
  
  
  
  if(sigstructure=="equalvar"){
    sig2draw.resize(1);
  }
  
  if(sigstructure=="unequalvar"){
    sig2draw.resize(k);
  }
  
  sig2draw.fill(sig2init);
  
  mat Sdraw = sigMatrix(sig2draw, k, sigstructure); 
  
  //initail rhodraw and Rdraw
  //double rhoinit = R::runif(0.0,1.0);
  double rhoinit = 0.9;
  
  mat Rdraw(k,k);
  Rdraw.fill(rhoinit);
  Rdraw.diag().ones();
  
  if(Rstructure=="CS" || Rstructure=="AR"){
    rhodraw.resize(1);
  }
  
  if(Rstructure=="UN"){
    rhodraw.resize(k*(k-1)/2);
  }
  
  rhodraw.fill(rhoinit);
  
  return List::create(Named("betadraw")=betadraw,
                      Named("sig2draw")=sig2draw,
                      Named("rhodraw")=rhodraw,
                      Named("Sdraw")=Sdraw, 
                      Named("Rdraw")=Rdraw);
  
}


//  Function to Update beta

// [[Rcpp::export]]
colvec betaUpdate(colvec& y, mat& X, int k, mat& Sdraw, 
                  mat& Rdraw, colvec& B0, mat& sigB){
  
  int n = X.n_rows, p = X.n_cols;
  mat VB = zeros<mat>(p,p);
  colvec muB = zeros<colvec>(p);
  
  mat invR = inv(Rdraw);
  mat invS = inv(Sdraw);
  mat invSig = invS*invR*invS;
  
  for(int j=0; j<n; j+=k){
    VB += trans(X.rows(j,j+k-1))*invSig*X.rows(j,j+k-1);
    muB += trans(X.rows(j,j+k-1))*invSig*y.rows(j,j+k-1);
  }
  
  mat V1 = inv(inv(sigB)+VB);
  colvec mu1 = V1*(inv(sigB)*B0+muB);
  
  colvec betadraw = vectorise(mvrnorm(1, mu1,V1));
  
  return betadraw;
  
}



// function to update variance

// [[Rcpp::export]]
double sigmapost(colvec& y, mat& X, int k, colvec& betadraw, mat& Sdraw, 
                 mat& Rdraw, colvec s0, mat sigS0){
  
  int n = X.n_rows, l = n/k;
  
  
  mat invR = inv(Rdraw);
  mat invS = inv(Sdraw);
  mat invSig = invS*invR*invS;
  
  double empSS = 0, logprior, logval;
  
  
  for (int j=0; j<n; j+=k){
    
    colvec res = y.rows(j,j+k-1)-X.rows(j,j+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  logprior = as_scalar(-0.5*trans(log(Sdraw.diag())-s0)*inv(sigS0)*(log(Sdraw.diag())-s0));
  // unequal variance, posterior distribution
  
  logval = -0.5*l*log(det(Sdraw*Rdraw*Sdraw))-0.5*empSS+logprior;
  
  return logval;
  
}



// [[Rcpp::export]]
List sigpostpar(colvec& y, mat& X, int k, colvec& betadraw,
                mat& Rdraw, double nu0, double tau0){
  
  int n = X.n_rows, l = n/k;
  double nu1, tau1;
  
  mat invR = inv(Rdraw);
  
  double empSS = 0;
  
  for (int j=0; j<n; j+=k){
    
    colvec res = y.rows(j,j+k-1)-X.rows(j,j+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invR*res);
    
  }
  
  // equal variance, posterior distribution
  
  nu1 = k*l/2+nu0;
  tau1 = 0.5*empSS+tau0; // rate for gamma distribution
  
  return List::create(Named("nu1")=nu1,
                      Named("tau1")=tau1);
  
}


// [[Rcpp::export]]
List sigmaUpdate(colvec& y, mat& X, int k,  std::string sigstructure, 
                 colvec& betadraw, mat& Rdraw, mat& Sdraw, colvec sig2draw, 
                 double nu0, double tau0, colvec& s0, mat& sigS0,
                 colvec sigma_sd, colvec sigcount){
  
  //int n = X.n_rows;
  double nu1, tau1;
  
  colvec sig2new = sig2draw;
  mat Snew = Sdraw;
  
  double u, signewpost, sigdrawpost, accprob;
  
  List sigpostout;
  
  if(sigstructure== "equalvar"){
    
    sigpostout = sigpostpar(y,X,k,betadraw,Rdraw,nu0,tau0);
    nu1 = as<double>(sigpostout["nu1"]);
    tau1 = as<double>(sigpostout["tau1"]); // scale for inverse-gamma distribution
    
    //RNGScope scope;
    sig2draw = 1/(R::rgamma(nu1,1/tau1)); //rate(gamma)=scale(inversegamma),scale(gamma)=1/tau1
    Sdraw = sigMatrix(sig2draw,k,sigstructure);
    
  }
  
  if(sigstructure== "unequalvar"){
    
    
    NumericVector srow, srandrow;
    int srowsize;
    
    NumericVector::iterator srowit;
    
    srow = seq(0,k-1);
    srowsize = srow.size();
    srandrow = Rcpp::RcppArmadillo::sample(srow,srowsize,FALSE);
    
    
    u = log(R::runif(0,1));
    
    for(srowit=srandrow.begin();srowit!=srandrow.end();srowit++){
      
      sig2new = sig2draw;
      
      sig2new(*srowit) = sig2draw(*srowit)+R::rnorm(0,sigma_sd(*srowit));
      
      if(sig2new(*srowit)>0){
        Snew = diagmat(sqrt(sig2new));
        signewpost = sigmapost(y,X,k,betadraw,Snew,Rdraw,s0,sigS0);
      }else{
        signewpost = R_NegInf;
      }
      
      sigdrawpost = sigmapost(y,X,k,betadraw,Sdraw,Rdraw,s0,sigS0);
      
      accprob = signewpost - sigdrawpost;
      
      
      if(u<accprob){
        
        sig2draw(*srowit) = sig2new(*srowit);
        Sdraw = Snew;
        sigcount(*srowit) += 1;    
      }     
    }
    
  }
  
  
  
  return List::create(Named("sig2draw")=sig2draw,
                      Named("Sdraw")=Sdraw,
                      Named("sigcount")=sigcount);
  
}




// function to update correlation matrix

// [[Rcpp::export]]
double detR(mat& Rdraw, int row, int col, int change){
  
  mat Rtemp = Rdraw;
  
  Rtemp(row,col) = change;
  Rtemp(col,row) = change;
  
  return (det(Rtemp));
}


// [[Rcpp::export]]
double rhoUpdate(mat Rdraw, double rhodraw, 
                 int row, int col, double division){
  
  
  //int nrow = Rdraw.n_rows, ncol = Rdraw.n_cols;
  
  double d1 = detR(Rdraw, row, col, 1),
    d2 = detR(Rdraw, row, col,-1),
    d0 = detR(Rdraw, row, col, 0);
  
  double a = 0.5*(d1+d2-2*d0), b = 0.5*(d1-d2), c = d0;
  
  NumericVector sol = realRoot(a,b,c);
  
  NumericVector diff = abs(rhodraw-sol);
  
  double mindiff = min(diff)/division;
  
  double rhonew = rhodraw + R::runif(-mindiff,mindiff);
  
  
  return rhonew;
  
}




// [[Rcpp::export]]
double Rpost(colvec& y, mat& X, int k, 
             colvec& betadraw, mat& Sdraw, mat& R){
  
  int n = X.n_rows, l = n/k;
  
  mat invR = inv(R);
  mat invS = inv(Sdraw);
  mat invSig = invS*invR*invS;
  
  double empSS = 0;
  
  for (int j=0; j<n; j+=k){
    
    colvec res = y.rows(j,j+k-1)-X.rows(j,j+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  return (-0.5*l*log(det(R))-0.5*empSS);
  
}


// function to update correlation matrix for CS, AR,UN

// [[Rcpp::export]]
List RUpdate(colvec& y, mat& X, int k, std::string Rstructure,
             colvec& betadraw, mat& Sdraw, colvec rhodraw, mat& Rdraw,  
             mat& division, mat& rhocount){
  
  
  mat Rnew = Rdraw;
  
  double rhonew, Rnewpost, Rdrawpost, accprob;
  
  //colvec rhodraw;
  
  double u = log(R::runif(0,1));
  
  //cout<<u<<endl;
  
  if(Rstructure=="CS"){
    
    rhonew = rhoUpdate(Rdraw, Rdraw(0,1), 0, 1, as_scalar(division)); 
    
    Rnew.fill(rhonew);
    Rnew.diag().ones();
    
    if(ispd(Rnew)){
      Rnewpost = Rpost(y, X, k, betadraw, Sdraw, Rnew);
    }else{
      Rnewpost = R_NegInf;
    }
    
    Rdrawpost = Rpost(y, X, k, betadraw, Sdraw, Rdraw);
    
    accprob = Rnewpost-Rdrawpost;
    //cout<<accprob<<endl;
    
    if(u < accprob){
      
      rhodraw = rhonew;
      Rdraw = Rnew;    
      rhocount += 1;
      
    }
    
    
  }
  
  
  if(Rstructure=="AR"){
    
    //rhodraw.resize(1);
    
    rhonew = rhoUpdate(Rdraw, Rdraw(0,1), 0, 1, as_scalar(division)); 
    
    for(int i=0; i<k-1; i++){
      for(int j=1; j<k; j++){
        
        Rnew(i,j) = pow(rhonew,abs(j-i));
        Rnew(j,i) = Rnew(i,j);
        
      }
    }   
    
    Rnew.diag().ones();
    
    if(ispd(Rnew)){
      Rnewpost = Rpost(y, X, k, betadraw, Sdraw, Rnew);
    }else{
      Rnewpost = R_NegInf;
    }
    
    Rdrawpost = Rpost(y, X, k, betadraw, Sdraw, Rdraw);
    
    accprob = Rnewpost-Rdrawpost;
    
    if(u < accprob){
      
      rhodraw = rhonew;
      Rdraw = Rnew;    
      rhocount += 1;
      
    }
  }
  
  
  
  if(Rstructure=="UN"){
    
    NumericVector row, col, randrow, randcol;
    int rowsize, colsize;
    
    NumericVector::iterator rowit, colit;
    
    row = seq(0,k-2);
    rowsize = row.size();
    randrow = Rcpp::RcppArmadillo::sample(row,rowsize,FALSE);
    
    
    for(rowit=randrow.begin();rowit!=randrow.end();rowit++){
      // cout<< *rowit<<endl;
      
      col = seq(*rowit+1,k-1);
      colsize = col.size();
      randcol = Rcpp::RcppArmadillo::sample(col,colsize,FALSE);
      
      for(colit=randcol.begin();colit!=randcol.end();colit++){
        
        //rhodraw = Rdraw(i,j);
        //cout<< *colit<<endl;
        
        rhonew = rhoUpdate(Rdraw,Rdraw(*rowit,*colit),*rowit,*colit,division(*rowit,*colit));
        
        Rnew(*rowit,*colit) = rhonew;
        Rnew(*colit,*rowit) = rhonew;
        Rnew.diag().ones();
        
        if(ispd(Rnew)){
          Rnewpost = Rpost(y, X, k, betadraw, Sdraw, Rnew);
        }else{
          Rnewpost = R_NegInf;
        }
        
        Rdrawpost = Rpost(y, X, k, betadraw, Sdraw, Rdraw);
        
        accprob = Rnewpost-Rdrawpost;
        
        if(u < accprob){
          //rhodraw = rhonew;
          Rdraw = Rnew;    
          rhocount(*rowit,*colit) += 1;
        }
        // cout<< Rdraw<<endl;
      }
    }
    
    rhodraw = corvec(Rdraw);
  }
  
  
  
  return List::create(Named("rhodraw")=rhodraw,
                      Named("Rdraw")=Rdraw,   
                      Named("rhocount")=rhocount);
  
}


// function of likelihood for all pars


// function of likelihood for all pars


// [[Rcpp::export]]
double logLik(colvec& y, mat& X, int k, std::string sigstructure,
              colvec& betadraw, mat& Sdraw, colvec& sig2draw, mat& Rdraw,
              colvec& B0, mat& sigB, double nu0, double tau0,
              colvec& s0, mat& sigS0){
  
  int n = X.n_rows, l = n/k, p = sigB.n_rows;
  
  double logy, logbeta, logsigma = 0, logfull;

  mat invR = inv(Rdraw);
  mat invS = inv(Sdraw);
  mat invSig = invS*invR*invS;
  
  double empSS = 0;
  
  for (int i=0; i<n; i+=k){
    
    colvec res = y.rows(i,i+k-1)-X.rows(i,i+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  // log likelihood of data y, pi = M_PI
  //logy = as_scalar(-(k/4)*log(2*M_PI)-l*log(det(Sdraw))-(l/2)*log(det(Rdraw))-0.5*empSS);
  logy = as_scalar(-0.5*k*l*log(2*M_PI)-0.5*l*log(det(Sdraw*Rdraw*Sdraw))-0.5*empSS);
  
  
  
  // log prior for beta
  //logbeta = as_scalar(-(p/4)*log(2*M_PI)-0.5*log(det(sigB))-0.5*(trans(betadraw-B0)*inv(sigB)*(betadraw-B0)));
  logbeta = as_scalar(-0.5*p*log(2*M_PI)-0.5*log(det(sigB))-0.5*(trans(betadraw-B0)*inv(sigB)*(betadraw-B0)));
  
  if(sigstructure=="equalvar"){
    
    //logsigma.resize(1);
    logsigma = as_scalar(nu0*log(tau0)-log(Rf_gammafn(nu0))-(nu0+1)*log(sig2draw)-tau0/sig2draw);
    
  }
  
  if(sigstructure=="unequalvar"){
    
    //logsigma.resize(k);
    
    logsigma = as_scalar(-0.5*k*log(2*M_PI)-0.5*log(det(sigS0))-0.5*(trans(log(Sdraw.diag())-s0))*inv(sigS0)*(log(Sdraw.diag())-s0));
  }

  logfull = logy+logbeta+logsigma;
  
  return logfull;
  
}









// [[Rcpp::export]]
List mcmcUpdate(colvec& y, mat& X, int k, std::string sigstructure,
                std::string Rstructure, colvec& B0, mat& sigB, 
                double nu0, double tau0, colvec& s0, mat& sigS0, 
                colvec sigma_sd, mat division, int niter){
  
  int p = X.n_cols;
  //int nkeep = niter-nburn;
  
  colvec betadraw, sig2draw, rhodraw;
  
  mat Sdraw, Rdraw;
  
  double deviance;    
  
  colvec sigcount;
  mat rhocount;
  
  List initialpars, Sout, Rout;
  
  
  // get initial values
  initialpars = initialGen(y, X, k,sigstructure, Rstructure);
  
  betadraw = as<arma::colvec>(initialpars["betadraw"]);
  sig2draw = as<arma::colvec>(initialpars["sig2draw"]);
  rhodraw = as<arma::colvec>(initialpars["rhodraw"]);
  Sdraw = as<arma::mat>(initialpars["Sdraw"]);
  Rdraw = as<arma::mat>(initialpars["Rdraw"]);
  
  int nsig = sig2draw.n_elem;
  
  
  if(sigstructure=="equalvar"){
    sigcount.resize(1);
  }
  
  if(sigstructure=="unequalvar"){
    sigcount.resize(k);
  }
  
  sigcount.fill(0);
  
  
  if(Rstructure=="CS" || Rstructure=="AR"){
    rhodraw.resize(1);
    rhocount.resize(1,1);
  }
  
  if(Rstructure=="UN"){
    rhodraw.resize(k*(k-1)/2);
    rhocount.resize(k,k);
  }
  
  rhocount.fill(0);
  
  int ncor = rhodraw.n_elem;
  
  
  mat betamat = zeros<mat>(niter,p),
    sig2mat = zeros<mat>(niter,nsig),
    rhomat = zeros<mat>(niter,ncor);
  colvec devmat = zeros<colvec>(niter);
  
  // Update iterations
  
  for(int iter=0; iter<niter; iter++){
    
    betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);
    betamat.row(iter) = trans(betadraw);
    
    Sout = sigmaUpdate(y,X,k,sigstructure,betadraw,Rdraw,Sdraw, sig2draw, 
                       nu0,tau0,s0,sigS0,sigma_sd,sigcount);
    sig2draw = as<colvec>(Sout["sig2draw"]);
    Sdraw = as<mat>(Sout["Sdraw"]);
    sigcount = as<colvec>(Sout["sigcount"]);
    sig2mat.row(iter) = trans(sig2draw);  
    
    Rout = RUpdate(y,X,k,Rstructure,betadraw,Sdraw,
                   rhodraw,Rdraw, division, rhocount);
    rhodraw = as<colvec>(Rout["rhodraw"]);
    Rdraw = as<mat>(Rout["Rdraw"]);
    rhocount = as<mat>(Rout["rhocount"]);
    rhomat.row(iter) = trans(rhodraw);
    
    deviance = -2*logLik(y,X,k,sigstructure,betadraw,Sdraw,sig2draw, 
                         Rdraw,B0,sigB,nu0,tau0,s0,sigS0);
    
    devmat(iter) = deviance;
          
  }
  
  return List::create(Named("sigstructure")=sigstructure,
                      Named("Rstructure")=Rstructure,
                      Named("sigcount")=sigcount,
                      Named("rhocount")=rhocount,
                      Named("betamat")=betamat,
                      Named("sig2mat")=sig2mat,
                      Named("rhomat")=rhomat,
                      Named("devmat")=devmat);
  
  
}


// Adaptive MCMC function for all pars

// [[Rcpp::export]]
List adaptiveUpdate(colvec& y, mat& X, int k, std::string sigstructure,
                    std::string Rstructure, colvec& B0, mat& sigB, 
                    double nu0, double tau0, colvec& s0, mat& sigS0,
                    colvec sigma_sd, mat division, int nadap, int niter){
  
  int p = X.n_cols;
  //int nkeep = niter-nburn;
  const double lowaccrate = 0.2, upaccrate = 0.5;
  
  //lowaccrate = 0.234; 
  //upaccrate = 0.44;
  
  
  double lowacc_count = lowaccrate*nadap, 
    upacc_count = upaccrate*nadap;
  
  
  colvec betadraw, sig2draw, rhodraw;
  
  mat Sdraw, Rdraw;
  
  double deviance;    
  
  colvec sigcount;
  mat rhocount;
  NumericVector count;
  
  List initialpars, Sout, Rout;
  
  
  // get initial values
  initialpars = initialGen(y, X, k,sigstructure, Rstructure);
  
  betadraw = as<arma::colvec>(initialpars["betadraw"]);
  sig2draw = as<arma::colvec>(initialpars["sig2draw"]);
  rhodraw = as<arma::colvec>(initialpars["rhodraw"]);
  Sdraw = as<arma::mat>(initialpars["Sdraw"]);
  Rdraw = as<arma::mat>(initialpars["Rdraw"]);
  
  int nsig = sig2draw.n_elem;
  
  
  if(sigstructure=="equalvar"){
    sigcount.resize(1);
  }
  
  if(sigstructure=="unequalvar"){
    sigcount.resize(k);
  }
  
  sigcount.fill(0);
  
  
  if(Rstructure=="CS" || Rstructure=="AR"){
    rhodraw.resize(1);
    rhocount.resize(1,1);
  }
  
  if(Rstructure=="UN"){
    rhodraw.resize(k*(k-1)/2);
    rhocount.resize(k,k);
  }
  
  rhocount.fill(0);
  
  // set up initial count
  if(sigstructure=="equalvar"){
    
    if(Rstructure=="CS" || Rstructure=="AR"){
      count = wrap(conv_to<colvec>::from(rhocount));
    }
    
    if(Rstructure=="UN"){
      count = wrap(corvec(rhocount));
    }
  }
  
  
  if(sigstructure=="unequalvar"){
    
    if(Rstructure=="CS" || Rstructure=="AR"){
      count = wrap(join_cols(sigcount,conv_to<colvec>::from(rhocount)));
    }
    
    if(Rstructure=="UN"){
      
      count = wrap(join_cols(sigcount,corvec(rhocount)));
    }
  }
  
  
  
  int ncor = rhodraw.n_elem;
  
  
  mat betamat = zeros<mat>(niter,p),
    sig2mat = zeros<mat>(niter,nsig),
    rhomat = zeros<mat>(niter,ncor);
  colvec devmat = zeros<colvec>(niter);
  
  // Update iterations
  double testno = 0;
  
  while(is_true(any(count<lowacc_count)) || is_true(any(count>upacc_count))){
    
    testno += 1;
    cout<<"test "<< testno<<endl;
    
    for(int iadap=0; iadap<nadap; iadap++){
      
      betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);
      
      Sout = sigmaUpdate(y,X,k,sigstructure,betadraw,Rdraw,Sdraw, sig2draw, 
                         nu0,tau0,s0,sigS0,sigma_sd,sigcount);
      sig2draw = as<colvec>(Sout["sig2draw"]);
      Sdraw = as<mat>(Sout["Sdraw"]);
      sigcount = as<colvec>(Sout["sigcount"]);
        
        
      Rout = RUpdate(y,X,k,Rstructure,betadraw,Sdraw,
                     rhodraw,Rdraw,division,rhocount);
      rhodraw = as<colvec>(Rout["rhodraw"]);
      Rdraw = as<mat>(Rout["Rdraw"]);
      rhocount = as<mat>(Rout["rhocount"]);
      
    }
    
    
    // adaptive MH phase
    if(sigstructure=="equalvar"){
      
      if(Rstructure=="CS" || Rstructure=="AR"){
        count = wrap(conv_to<colvec>::from(rhocount));
      }
      
      if(Rstructure=="UN"){
        count = wrap(corvec(rhocount));
      }
    }
    
    
    if(sigstructure=="unequalvar"){
      
      if(Rstructure=="CS" || Rstructure=="AR"){
        count = wrap(join_cols(sigcount,conv_to<colvec>::from(rhocount)));
      }
      
      if(Rstructure=="UN"){
        
        count = wrap(join_cols(sigcount,corvec(rhocount)));
      }
    }
    
    
    
    if(sigstructure=="unequalvar"){
      
      for(int i=0; i<k; i++){
        
        if(sigcount(i)<lowacc_count && sigma_sd(i)>(0.05*sigma_sd(i))){
          sigma_sd(i) -= 0.05*sigma_sd(i);
        }
        
        if(sigcount(i)>upacc_count){
          sigma_sd(i) += 0.05*sigma_sd(i);
        }
        
      }
      
    }
    
    
    if(Rstructure=="CS" || Rstructure=="AR"){
      
      if(as_scalar(rhocount)<lowacc_count){
        division += 0.1*division;
      }
      
      if(as_scalar(rhocount)>upacc_count && as_scalar(division)>as_scalar(0.1*division)){
        division -= 0.1*division;
      }
      
    }
    
    
    
    if(Rstructure=="UN"){
      
      for(int i=0; i<k-1; i++){
        for(int j=i+1; j<k; j++){
          
          if(rhocount(i,j)<lowacc_count){
            division(i,j) += 0.1*division(i,j);
          }
          
          if(rhocount(i,j)>upacc_count && division(i,j)>(0.1*division(i,j))){
            division(i,j) -= 0.1*division(i,j);
          }
          
        }
      } 
    }
    
    
    
    
    cout<<sigcount<<endl;
    cout<<sigma_sd<<endl;
    cout<<rhocount<<endl;
    cout<<division<<endl;
    
    sigcount.zeros();
    rhocount.zeros();
    //count.zeros();
    
  }
  
  cout<<"pass" <<endl;
  
  // Update iterations
  
  for(int iter=0; iter<niter; iter++){
    
    betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);
    betamat.row(iter) = trans(betadraw);
    
    Sout = sigmaUpdate(y,X,k,sigstructure,betadraw,Rdraw,Sdraw, sig2draw, 
                       nu0,tau0,s0,sigS0,sigma_sd,sigcount);
    sig2draw = as<colvec>(Sout["sig2draw"]);
    Sdraw = as<mat>(Sout["Sdraw"]);
    sigcount = as<colvec>(Sout["sigcount"]);
    sig2mat.row(iter) = trans(sig2draw);  
      
    Rout = RUpdate(y,X,k,Rstructure,betadraw,Sdraw,
                   rhodraw,Rdraw, division, rhocount);
    rhodraw = as<colvec>(Rout["rhodraw"]);
    Rdraw = as<mat>(Rout["Rdraw"]);
    rhocount = as<mat>(Rout["rhocount"]);
    rhomat.row(iter) = trans(rhodraw);
    
    deviance = -2*logLik(y,X,k,sigstructure,betadraw,Sdraw,sig2draw, 
                         Rdraw,B0,sigB,nu0,tau0,s0,sigS0);
    
    devmat(iter) = deviance;
    
  }
  
  return List::create(Named("sigcount")=sigcount,
                      Named("rhocount")=rhocount,
                      Named("sigma_sd")=sigma_sd,
                      Named("division")=division,
                      Named("betamat")=betamat,
                      Named("sig2mat")=sig2mat,
                      Named("rhomat")=rhomat,
                      Named("devmat")=devmat);
  
}
