//[[Rcpp::depends(RcppArmadillo)]]

#include<RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat mergeMatrix(mat B1, mat B2, mat B3){
  
  // B1 and B3 has equal rows, B4,
  
  mat M1 = join_rows(B1,B3);
  mat M2 = join_rows(trans(B3),B2);
  mat BM = join_cols(M1,M2);
  
  return(BM);
  
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


// [[Rcpp::export]]
mat mvrnorm(int n, colvec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return (repmat(mu, 1, n).t() + Y * chol(sigma));
}



//[[Rcpp::export]]
mat sigMatrix(colvec& gamma, colvec& alpha, List& t){
  
  int nsize = t.size();
  
  colvec sig = zeros<colvec>(1), tempsig;
  
  //colvec sig = gamma(0)+alpha(0)*as<arma::colvec>(t[0]);
  
  for(int i=0;i<nsize;i++){
    colvec tempsig = gamma[i]+alpha[i]*as<arma::colvec>(t[i]);
    //cout<<tempsig<<endl;
    
    sig = join_cols(sig,tempsig);
  }
  
  
  sig.shed_row(0);
  
  //cout<<sig<<endl;
  //colvec sqrts= sqrt(sig);                                           
  
  mat Sigma = diagmat(sqrt(sig));
  
  return Sigma;
}



// Generate Initail pars

//[[Rcpp::export]]
List initialGen(colvec& y, mat& X, List& t){
  
  //RNGScope scope;
  
  int nblock = t.size();
  
  int k = 0;
  
  for(int nb=0; nb<nblock; nb++){
    k += as<NumericVector>(t[nb]).size();
  }
  
  IntegerVector ntemp = seq_len(nblock); 
  int ncor = sum(ntemp);
  
  int bsize ;
  
  
  //initial value for betadraw
  colvec betadraw = inv(trans(X)*X)*(trans(X)*y);
  
  // generate initial values of alpha and gamma, Sdraw
  colvec alphadraw = zeros<colvec>(nblock);
  colvec gammadraw = zeros<colvec>(nblock);
  alphadraw.fill(R::runif(0.0,1.0));
  
  for(int nb=0; nb<nblock; nb++){
    gammadraw(nb) = max(-alphadraw(nb)*as<arma::colvec>(t[nb]))+1;
  }
  
  mat Sdraw = sigMatrix(gammadraw,alphadraw,t); 
  
  //initail rhodraw and Rdraw
  
  colvec rhodraw = zeros<colvec>(ncor);
  rhodraw.fill(R::runif(0.0,1.0));
  
  List BM(ncor);
  
  mat Btemp;
  // mat Btemp=zeros<mat>(brow,brow);
  
  for(int i=0; i<nblock;i++){
    bsize = as<NumericVector>(t[i]).size();
    Btemp.resize(bsize,bsize);
    
    Btemp.fill(rhodraw(i));
    BM[i] = Btemp;
    
  }
  
  Btemp = zeros<mat>(as<NumericMatrix>(BM[0]).nrow(),as<NumericMatrix>(BM[1]).ncol());
  Btemp.fill(rhodraw(2));
  BM[2] = Btemp;
  
  mat Rdraw = mergeMatrix(BM[0],BM[1],BM[2]);
  
  Rdraw.diag().ones();
  
  return List::create(Named("betadraw")=betadraw,
                      Named("alphadraw")=alphadraw,
                      Named("gammadraw")=gammadraw,
                      Named("rhodraw")=rhodraw,
                      Named("Sdraw")=Sdraw, 
                      Named("Rdraw")=Rdraw,
                      Named("BM")=BM);
  
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



// Functions to Update gamma and alpha that give sigma^2 = gamma+alpha*t

// [[Rcpp::export]]
double gammapost(colvec& y, mat& X, int k, colvec t, colvec& betadraw, 
                 mat& Rdraw, mat& S, double gamma, double alphadraw, 
                 double delta, double eta){
  
  
  int n = X.n_rows, l = n/k;
  
  mat invR = inv(Rdraw);
  mat invS = inv(S);
  mat invSig = invS*invR*invS;
  
  double empSS = 0;
  
  for (int j=0; j<n; j+=k){
    
    colvec res = y.rows(j,j+k-1)-X.rows(j,j+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  //colvec alphat = -alphadraw*t, gammat = - gamma/t;
  //double maxalphat = max(alphat), maxgammat = max(gammat);
  double zalpha = max(-alphadraw*t)/delta, 
    zgamma = max(-gamma/t)/eta;
  
  //NumericVector valpha = NumericVector::create(zalpha);
  //NumericVector vgamma = NumericVector::create(zgamma);
  
  //RNGScope scope;
  
  double dalpha = R::dnorm(zalpha,0.0,1.0,0);
  double dgamma = R::dnorm(zgamma,0.0,1.0,0);
  
  double logval = -l*log(det(S))-0.5*empSS-0.5*pow(gamma/delta,2)-log(1-dalpha)-log(1-dgamma);
  
  return logval;
}



// [[Rcpp::export]]
double alphapost(colvec& y,mat& X, int k, colvec t, colvec& betadraw, 
                 mat& Rdraw, mat& S, double gammadraw, double alpha, 
                 double delta, double eta){
  
  int n = X.n_rows, l = n/k;
  
  mat invR = inv(Rdraw);
  mat invS = inv(S);
  mat invSig = invS*invR*invS;
  
  double empSS = 0;
  
  for (int j=0; j<n; j+=k){
    
    colvec res = y.rows(j,j+k-1)-X.rows(j,j+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  //colvec alphat = -alpha*t, gammat = - gammadraw/t;
  //double maxalphat = max(alphat), maxgammat = max(gammat);
  double zalpha = max(-alpha*t)/delta, 
    zgamma = max(-gammadraw/t)/eta;
  
  double dalpha = R::dnorm(zalpha,0.0,1.0,0);
  double dgamma = R::dnorm(zgamma,0.0,1.0,0);
  
  double logval =  -0.5*l*log(det(S*Rdraw*S))-0.5*empSS-0.5*pow((alpha/eta),2)-log(1-dalpha)-log(1-dgamma);
  
  return logval;
  
}


// [[Rcpp::export]]
List sigmaUpdate(colvec& y, mat& X,  List& t, 
                 colvec& betadraw, mat& Rdraw, mat& Sdraw, 
                 colvec gammadraw, colvec alphadraw, 
                 double delta, double eta,
                 colvec gamma_sd, colvec alpha_sd,
                 colvec gammacount, colvec alphacount){
  
  int nblock = t.size();
  int k  = 0;
  for(int nb=0; nb<nblock; nb++){
    k += as<NumericVector>(t[nb]).size();
  }
  
  
  colvec gammanew = gammadraw, 
    alphanew = alphadraw;
  
  mat Snew_ga = Sdraw, Snew_al = Sdraw;
  
  //int tlen = k/nblock;
  double gammanewpost, gammadrawpost,
  alphanewpost, alphadrawpost,
  accgamma, accalpha;
  
  double alphaj, maxalphat, gammaj, maxgammat;
  NumericVector tj ;
  
  //colvec alphat = zeros<colvec>(tlen),
  //       gammat = zeros<colvec>(tlen);
  
  colvec alphat, gammat;
  
  double u = log(R::runif(0,1));
  
  //int j=0;
  for(int j=0; j<nblock; j++){
    
    tj = t[j];
    
    // MH Update for gamma
    gammanew = gammadraw;
    gammanew(j) = gammadraw(j)+R::rnorm(0, gamma_sd(j));
    
    alphaj = alphadraw(j);
    alphat = -alphaj*tj;
    maxalphat = max(alphat);
    
    if(gammanew(j) > maxalphat){
      
      Snew_ga = sigMatrix(gammanew, alphadraw, t);
      
      gammanewpost = gammapost(y, X, k, t[j], betadraw, Rdraw, Snew_ga,
                               gammanew(j), alphadraw(j), delta, eta);
      
      gammadrawpost = gammapost(y,X, k, t[j], betadraw, Rdraw, Sdraw,
                                gammadraw(j), alphadraw(j), delta, eta);
      
      accgamma = gammanewpost-gammadrawpost;
      
      if(u < accgamma){
        gammadraw(j) = gammanew(j);
        Sdraw = Snew_ga;
        gammacount(j) += 1;
      } 
    }
    
    else{
      gammadraw(j) = gammadraw(j);
      Sdraw = Sdraw;
      gammacount(j) = gammacount(j); 
    } 
    
    
    
    // MH Update for alpha
    
    alphanew = alphadraw;
    alphanew(j) = alphadraw(j)+R::rnorm(0,alpha_sd(j));
    
    gammaj = gammadraw(j);
    gammat = - gammaj/tj;
    maxgammat = max(gammat);
    
    
    if(alphanew(j)>maxgammat){
      
      Snew_al = sigMatrix(gammadraw, alphanew, t);
      
      alphanewpost = alphapost(y, X, k, t[j], betadraw, Rdraw, Snew_al,
                               gammadraw(j), alphanew(j), delta, eta);
      
      alphadrawpost = alphapost(y, X, k, t[j], betadraw, Rdraw, Sdraw,
                                gammadraw(j), alphadraw(j), delta, eta);
      
      accalpha = alphanewpost-alphadrawpost;
      
      if(u < accalpha){
        alphadraw(j) = alphanew(j);
        Sdraw = Snew_al;
        alphacount(j) +=1;
      }      
    }
    
    else{
      alphadraw(j) = alphadraw(j);
      Sdraw = Sdraw;
      alphacount(j) = alphacount(j);
    }
    
  }
  
  
  return List::create(Named("gammadraw")=gammadraw,
                      Named("alphadraw")=alphadraw,
                      Named("Sdraw")=Sdraw,
                      Named("gammacount")=gammacount,
                      Named("alphacount")=alphacount);
  
}



// Functions to update correlations for block matrix with CS-CS structure

// [[Rcpp::export]]
double detR(mat& Rdraw, int row, int col, int change){
  
  mat Rtemp = Rdraw;
  
  Rtemp(row,col) = change;
  Rtemp(col,row) = change;
  
  return (det(Rtemp));
}


// [[Rcpp::export]]
List blockUpdate(mat& Rdraw, double rhodraw, mat B, 
                 int row, int col, double division){
  
  double d1 = detR(Rdraw, row, col, 1),
    d2 = detR(Rdraw, row, col,-1),
    d0 = detR(Rdraw, row, col, 0);
  
  double a = 0.5*(d1+d2-2*d0), b = 0.5*(d1-d2), c = d0;
  
  NumericVector sol = realRoot(a,b,c);
  
  NumericVector diff = abs(rhodraw-sol);
  
  double mindiff = min(diff)/division;
  
  double rhonew = rhodraw + R::runif(-mindiff,mindiff);
  
  B.fill(rhonew);
  
  return List::create(Named("rhonew")=rhonew, Named("B")=B);
  
}


// [[Rcpp::export]]
double Rpost(colvec& y, mat& X, int k, colvec& betadraw, 
             mat& Sdraw, mat& R){
  
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


// [[Rcpp::export]]
List RUpdate(colvec& y, mat& X, int k, int k1,colvec& betadraw, mat& Sdraw, 
               mat& Rdraw, colvec& rhodraw, List& BM, 
               colvec division, colvec rhocount){
  
  int nblock = BM.size(); // number of blocks in block matrix
  
  colvec rowindex = as<arma::colvec>(NumericVector::create(0,k1,0));
  colvec colindex = as<arma::colvec>(NumericVector::create(1,k1+1,k1));
  
  mat Rnew = Rdraw;
  
  List out;
  
  double rhonew, Rnewpost, Rdrawpost, accprob;
  
  double u = log(R::runif(0,1));
  
  //colvec rhocount = zeros<colvec>(nblock);
  
  for(int i=0; i<nblock; i++){
    
    mat Bi = as<arma::mat>(BM[i]);
    
    out = blockUpdate(Rdraw, rhodraw(i), Bi, rowindex(i), colindex(i), division(i));
    
    rhonew = as_scalar(as<arma::colvec>(out["rhonew"]));
    
    mat Bnew = as<arma::mat>(out["B"]);
    
    if(i==0){      
      Rnew = mergeMatrix(Bnew,as<arma::mat>(BM[1]), as<arma::mat>(BM[2]));
      //Rnew.diag().ones();
    }
    
    else if(i==1){
      Rnew = mergeMatrix(as<arma::mat>(BM[0]), Bnew, as<arma::mat>(BM[2]));
      //Rnew.diag().ones(); 
    }
    
    else{
      Rnew = mergeMatrix(as<arma::mat>(BM[0]), as<arma::mat>(BM[1]), Bnew);
      // Rnew.diag().ones();
    }
    
    Rnew.diag().ones();
    
    Rnewpost = Rpost(y, X, k, betadraw, Sdraw, Rnew);
    Rdrawpost = Rpost(y, X, k, betadraw, Sdraw, Rdraw);
    
    accprob = Rnewpost-Rdrawpost;
    
    if(u < accprob){
      
      rhodraw(i) = rhonew;
      Rdraw = Rnew;
      BM[i] = Bnew;
      rhocount(i) += 1;
    }
  }
  
  return List::create(Named("rhodraw")=rhodraw,
                      Named("Rdraw")=Rdraw,
                      Named("BM")=BM,
                      Named("rhocount")=rhocount);
  
}


// Function of log likelihood function for data and parameters, use for DIC

// [[Rcpp::export]]
double logLik(colvec& y, mat& X, int k, int nblock, List t,
              colvec& betadraw, mat& Sdraw, mat& Rdraw, 
              colvec& B0, mat& sigB, double delta, double eta,
              colvec& gammadraw, colvec& alphadraw){
  
  int n = X.n_rows, p = X.n_cols, l = n/k;
  
  double logy, logbeta, logprior,logfull;
  
  
  int tlen = k/nblock;
  double alphaj, gammaj, maxalphat, maxgammat; 
  
  colvec alphat = zeros<colvec>(tlen),
    gammat = zeros<colvec>(tlen);
  
  colvec loggamma(nblock), logalpha(nblock);
  colvec valpha(nblock), vgamma(nblock), 
  dalpha(nblock), dgamma(nblock);
  
  NumericVector tj(tlen);
  
  mat invR = inv(Rdraw);
  mat invS = inv(Sdraw);
  mat invSig = invS*invR*invS;
  
  double empSS = 0;
  
  for (int j=0; j<n; j+=k){
    
    colvec res = y.rows(j,j+k-1)-X.rows(j,j+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  // log likelihood of data y, pi = M_PI
  logy = as_scalar(-0.5*k*l*log(2*M_PI)-0.5*l*log(det(Sdraw*Rdraw*Sdraw))-0.5*empSS);
  
  // log prior for beta
  logbeta = as_scalar(-0.5*p*log(2*M_PI)-0.5*log(det(sigB))-0.5*(trans(betadraw-B0)*inv(sigB)*(betadraw-B0)));
  
  
  //RNGScope scope;
  
  for(int j=0; j<nblock; j++){
    
    tj = t[j];
    
    alphaj = alphadraw(j);
    alphat = -alphaj*tj;
    maxalphat = max(alphat);
    
    gammaj = gammadraw(j);
    gammat = - gammaj/tj;
    maxgammat = max(gammat);
    
    valpha(j) = maxalphat/delta;
    vgamma(j) = maxgammat/eta;
    dalpha(j) = R::dnorm(valpha(j),0.0,1.0,0);
    dgamma(j) = R::dnorm(vgamma(j),0.0,1.0,0);
    
    loggamma(j) = -0.5*log(2*M_PI)-log(delta)-0.5*pow(gammadraw(j)/delta,2)-log(1-dalpha(j));
    logalpha(j) = -0.5*log(2*M_PI)-log(eta)-0.5*pow(alphadraw(j)/eta,2)-log(1-dgamma(j));
    
  }
  
  logprior = logbeta + sum(loggamma) + sum(logalpha);
  
  logfull = logy + logprior;
  
  return logfull;
} 




// [[Rcpp::export]]
List mcmcUpdate(colvec& y, mat& X,  List& t, 
                colvec& B0, mat& sigB, double delta, double eta, 
                colvec gamma_sd, colvec alpha_sd, colvec division,               
                int niter){
  
  int nblock = t.size();
  
  int k = 0;
  
  for(int nb=0; nb<nblock; nb++){
    k += as<NumericVector>(t[nb]).size();
  }

  int k1 = as<NumericVector>(t[0]).size();
  
  //IntegerVector ntemp = seq_len(nblock); 
  //int ncor = sum(ntemp);
  int nrho = nblock*(nblock+1)/2;
 
  int p = X.n_cols;
  //int nkeep = niter-nburn;
  
  
  colvec betadraw = zeros<vec>(p),
    gammadraw = zeros<vec>(nblock),
    alphadraw = zeros<vec>(nblock),
    rhodraw = zeros<vec>(nrho);
  
  mat Sdraw = zeros<mat>(k,k),
    Rdraw = zeros<mat>(k,k);
  
  double deviance;    
  
  List BM(nrho);
  
  
  mat betamat = zeros<mat>(niter,p),
    gammamat = zeros<mat>(niter,nblock),
    alphamat = zeros<mat>(niter,nblock),
    rhomat = zeros<mat>(niter,nrho);
 colvec devmat = zeros<colvec>(niter);
  
  colvec gammacount = zeros<vec>(nblock),
    alphacount = zeros<vec>(nblock),
    rhocount = zeros<vec>(nrho);
  
  List initialpars, Sout, Rout;

  // get initial values
  initialpars = initialGen(y, X, t);
  
  betadraw = as<arma::colvec>(initialpars["betadraw"]);
  gammadraw = as<arma::colvec>(initialpars["gammadraw"]);
  alphadraw = as<arma::colvec>(initialpars["alphadraw"]);
  rhodraw = as<arma::colvec>(initialpars["rhodraw"]);
  Sdraw = as<arma::mat>(initialpars["Sdraw"]);
  Rdraw = as<arma::mat>(initialpars["Rdraw"]);
  BM = initialpars["BM"];
  
  
  // Update iterations
  
  for(int iter=0; iter<niter; iter++){
    
    betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);
    betamat.row(iter) = trans(betadraw);
    
    Sout = sigmaUpdate(y,X,t,betadraw,Rdraw,Sdraw,gammadraw,alphadraw, 
                       delta,eta,gamma_sd,alpha_sd,gammacount,alphacount);
    gammadraw = as<colvec>(Sout["gammadraw"]);
    alphadraw = as<colvec>(Sout["alphadraw"]);
    Sdraw = as<mat>(Sout["Sdraw"]);
    gammacount = as<mat>(Sout["gammacount"]);
    alphacount = as<mat>(Sout["alphacount"]);
    
    gammamat.row(iter) = trans(gammadraw);
    alphamat.row(iter) = trans(alphadraw);
    
    
    Rout = RUpdate(y,X,k,k1,betadraw,Sdraw,Rdraw,
                   rhodraw,BM,division, rhocount);
    rhodraw = as<colvec>(Rout["rhodraw"]);
    Rdraw = as<mat>(Rout["Rdraw"]);
    BM = Rout["BM"];
    rhocount = as<colvec>(Rout["rhocount"]);
    rhomat.row(iter) = trans(rhodraw);
    
    deviance = logLik(y,X,k,nblock,t,betadraw,Sdraw,Rdraw,
                      B0,sigB,delta,eta,gammadraw,alphadraw);
    devmat(iter) = deviance;
    
  }
  
  return List::create(Named("gammacount")=gammacount,
                      Named("alphacount")=alphacount,
                      Named("rhocount")=rhocount,
                      Named("betamat")=betamat,
                      Named("gammamat")=gammamat,
                      Named("alphamat")=alphamat,
                      Named("rhomat")=rhomat,
                      Named("devmat")=devmat);
  
}




// Adaptive MCMC function for all pars

// [[Rcpp::export]]
List adaptiveUpdate(colvec& y, mat& X, List& t, colvec& B0, mat& sigB, 
                    double delta, double eta, colvec gamma_sd,
                    colvec alpha_sd, colvec division, 
                    int nadap, int niter){
  
  int nblock = t.size();
  
  int k = 0;
  
  for(int nb=0; nb<nblock; nb++){
    k += as<NumericVector>(t[nb]).size();
  }
  
  int k1 = as<NumericVector>(t[0]).size();
  
  //IntegerVector ntemp = seq_len(nblock); 
  //int ncor = sum(ntemp);
  int nrho = nblock*(nblock+1)/2;
 
  
  
  int p = X.n_cols;
  //int nkeep = niter-nburn;
  const double lowaccrate = 0.2, 
    upaccrate = 0.5;
  double lowacc_count = lowaccrate*nadap, 
    upacc_count = upaccrate*nadap;
  
  
  
  colvec betadraw = zeros<vec>(p),
    gammadraw = zeros<vec>(nblock),
    alphadraw = zeros<vec>(nblock),
    rhodraw = zeros<vec>(nrho);
  colvec sig2draw;
  
  mat Sdraw = zeros<mat>(k,k),
    Rdraw = zeros<mat>(k,k);
  
  double deviance;    
  
  List BM(nrho);
  
  
  mat betamat = zeros<mat>(niter,p),
    gammamat = zeros<mat>(niter,nblock),
    alphamat = zeros<mat>(niter,nblock),
    rhomat = zeros<mat>(niter,nrho),
    sig2mat = zeros<mat>(niter,k);
 colvec devmat = zeros<colvec>(niter);
  
  colvec gammacount = zeros<vec>(nblock),
    alphacount = zeros<vec>(nblock),
    rhocount = zeros<vec>(nrho);
  
  NumericVector count = wrap(join_cols(join_cols(gammacount,alphacount),rhocount));
  
  List initialpars, Sout, Rout;
  
  //double lowgamma, lowalpha, lowdivision; 
  
  
  // get initial values
  initialpars = initialGen(y, X, t);
  
  betadraw = as<arma::colvec>(initialpars["betadraw"]);
  gammadraw = as<arma::colvec>(initialpars["gammadraw"]);
  alphadraw = as<arma::colvec>(initialpars["alphadraw"]);
  rhodraw = as<arma::colvec>(initialpars["rhodraw"]);
  Sdraw = as<arma::mat>(initialpars["Sdraw"]);
  Rdraw = as<arma::mat>(initialpars["Rdraw"]);
  BM = initialpars["BM"];
  
  
  // adaptive MH phase
  
  double testno = 0;
  
  while(is_true(any(count<lowacc_count)) || is_true(any(count>upacc_count))){  
    
    testno += 1;
    cout<<"test"<<testno<<endl;
    
    for(int iadap=0; iadap<nadap; iadap++){
      
      betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);

      Sout = sigmaUpdate(y,X,t,betadraw,Rdraw,Sdraw,gammadraw,alphadraw, 
                         delta,eta,gamma_sd,alpha_sd,gammacount,alphacount);
      gammadraw = as<colvec>(Sout["gammadraw"]);
      alphadraw = as<colvec>(Sout["alphadraw"]);
      Sdraw = as<mat>(Sout["Sdraw"]);
      gammacount = as<mat>(Sout["gammacount"]);
      alphacount = as<mat>(Sout["alphacount"]);

      Rout = RUpdate(y,X,k,k1,betadraw,Sdraw,Rdraw,
                     rhodraw,BM,division, rhocount);
      rhodraw = as<colvec>(Rout["rhodraw"]);
      Rdraw = as<mat>(Rout["Rdraw"]);
      BM = Rout["BM"];
      rhocount = as<colvec>(Rout["rhocount"]);

    }
    
    count = wrap(join_cols(join_cols(gammacount,alphacount),rhocount));
    
    
    for(int j=0; j<nblock; j++){
      
      
      
      if(gammacount(j)<lowacc_count && gamma_sd(j)>(0.05*gamma_sd(j))){
        gamma_sd(j) -= 0.05*gamma_sd(j);
      }
      
      
      if(gammacount(j)>upacc_count){
        gamma_sd(j) += 0.05*gamma_sd(j);
      }
      
      
      if(alphacount(j)<lowacc_count && alpha_sd(j)>(0.05*alpha_sd(j))){
        alpha_sd(j) -= 0.05*alpha_sd(j);
      }
      
      if(alphacount(j)>upacc_count){
        alpha_sd(j) += 0.05*alpha_sd(j);
      }
      
    }
    
    for(int i=0; i<nrho;i++){
      
      if(rhocount(i)<lowacc_count){
        division(i) += 0.1*division(i);
      }
      
      if(rhocount(i)>upacc_count && division(i)>(0.1*division(i))){
        division(i) -= 0.1*division(i);
      }
      
    } 
    
    
    cout<<gammacount<<endl;
    cout<<gamma_sd<<endl;
    cout<<alphacount<<endl;
    cout<<alpha_sd<<endl;
    cout<<rhocount<<endl;
    cout<<division<<endl;
    
    gammacount.zeros();
    alphacount.zeros();
    rhocount.zeros();
    //count.zeros();
    
  }
  
  cout<<"pass" <<endl;
  
  // Update iterations
  
  for(int iter=0; iter<niter; iter++){
    
    betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);
    betamat.row(iter) = trans(betadraw);
    
    Sout = sigmaUpdate(y,X,t,betadraw,Rdraw,Sdraw,gammadraw,alphadraw, 
                       delta,eta,gamma_sd,alpha_sd,gammacount,alphacount);
    gammadraw = as<colvec>(Sout["gammadraw"]);
    alphadraw = as<colvec>(Sout["alphadraw"]);
    Sdraw = as<mat>(Sout["Sdraw"]);
    gammacount = as<mat>(Sout["gammacount"]);
    alphacount = as<mat>(Sout["alphacount"]);
    
    gammamat.row(iter) = trans(gammadraw);
    alphamat.row(iter) = trans(alphadraw);
    
    sig2mat.row(iter) = trans(pow(Sdraw.diag(),2));
    
    Rout = RUpdate(y,X,k,k1,betadraw,Sdraw,Rdraw,
                   rhodraw,BM,division, rhocount);
    rhodraw = as<colvec>(Rout["rhodraw"]);
    Rdraw = as<mat>(Rout["Rdraw"]);
    BM = Rout["BM"];
    rhocount = as<colvec>(Rout["rhocount"]);
    rhomat.row(iter) = trans(rhodraw);
    
    deviance = -2*logLik(y,X,k,nblock,t,betadraw,Sdraw,Rdraw,
                         B0,sigB,delta,eta,gammadraw,alphadraw);
    devmat(iter) = deviance;
    
    
  }
  
  return List::create(Named("gammacount")=gammacount,
                      Named("alphacount")=alphacount,
                      Named("rhocount")=rhocount,
                      Named("gamma_sd")=gamma_sd,
                      Named("alpha_sd")=alpha_sd,
                      Named("division")=division,
                      Named("betamat")=betamat,
                      Named("gammamat")=gammamat,
                      Named("alphamat")=alphamat,
                      Named("sig2mat")=sig2mat,
                      Named("rhomat")=rhomat,
                      Named("devmat")=devmat);
  
}


