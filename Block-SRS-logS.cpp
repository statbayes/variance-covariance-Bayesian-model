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
mat sigMatrix(colvec sig2draw, List t){
  
  
  int nsize = t.size(), nsig;
  
  colvec sig = zeros<colvec>(1), tempsig;
  
  for(int i=0; i<nsize; i++){
    nsig = as<NumericVector>(t[i]).size();
    
    tempsig.set_size(nsig);
    tempsig.fill(sqrt(sig2draw(i)));
    sig = join_cols(sig,tempsig);
    
  }
  
  sig.shed_row(0);
  mat Sigma = diagmat(sig);
  
  return Sigma;
}

// Generate Initail pars

//[[Rcpp::export]]
List initialGen(colvec& y, mat& X, List t){
  
  //RNGScope scope;
  
  int nblock = t.size();
  
  IntegerVector ntemp = seq_len(nblock); 
  int ncor = sum(ntemp);
  
  int bsize ;
  
  //initial value for betadraw
  colvec betadraw = inv(trans(X)*X)*(trans(X)*y);
  
  colvec sig2draw = ones<colvec>(nblock);
  
  mat Sdraw = sigMatrix(sig2draw, t); 
  
  //initail rhodraw and Rdraw
  
  colvec rhodraw = zeros<colvec>(ncor);
  rhodraw.fill(0.95);
  //R::runif(0,1)
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
                      Named("sig2draw")=sig2draw,
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




// Functions to Update variances: sigma

// [[Rcpp::export]]
double sigmapost(colvec& y, mat& X, int k, colvec& betadraw, mat& Sdraw, 
                 mat& Rdraw, colvec& sig2draw, colvec& s0, mat& sigS0){
  
  int n = X.n_rows, l = n/k;
  
  mat invR = inv(Rdraw);
  mat invS = inv(Sdraw);
  mat invSig = invS*invR*invS;
  
  double empSS = 0;
  
  for (int j=0; j<n; j+=k){
    
    colvec res = y.rows(j,j+k-1)-X.rows(j,j+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  double logprior = as_scalar(-0.5*trans(log(sig2draw)-s0)*inv(sigS0)*(log(sig2draw)-s0));
  double logval = as_scalar(-0.5*l*log(det(Sdraw*Rdraw*Sdraw))-0.5*empSS)+logprior;
  
  return logval;
  
}


// [[Rcpp::export]]
List sigmaUpdate(colvec& y, mat& X, List t,
                 colvec& betadraw, colvec sig2draw, mat& Sdraw, 
                 mat& Rdraw, colvec s0, mat sigS0,
                 colvec sigma_sd, colvec sigcount){
  
  // sig2draw is nblock*1 vector
  // s0 (nblock*1)
  // sigS0 (nblock*nblock)
  
  int nblock = t.size();
  int k  = 0;
  for(int nb=0; nb<nblock; nb++){
    k += as<NumericVector>(t[nb]).size();
  }
  
  
  colvec sig2new = sig2draw;
  mat Snew = Sdraw;
  
  
  double signewpost, sigdrawpost, accprob;
  
  
  double u = log(R::runif(0,1));
  
  
  for(int j=0; j<nblock; j++){
    
    sig2new = sig2draw;
    
    sig2new(j) = sig2draw(j)+R::rnorm(0,sigma_sd(j));
    
    if(sig2new(j)>0){
      Snew = sigMatrix(sig2new, t);
      
      signewpost = sigmapost(y,X,k,betadraw,Snew,Rdraw,sig2new,s0,sigS0);
    }else{
      signewpost = R_NegInf;
    }
    
    sigdrawpost = sigmapost(y,X,k,betadraw,Sdraw,Rdraw,sig2draw,s0,sigS0);
    
    accprob = signewpost - sigdrawpost;
    
    
    if(u<accprob){
      
      sig2draw(j) = sig2new(j);
      Sdraw = Snew;
      sigcount(j) += 1;    
    }     
  }
  
  
  
  return List::create(Named("sig2draw")=sig2draw,
                      Named("Sdraw")=Sdraw,
                      Named("sigcount")=sigcount);
  
}



// Fuctions to Update correlations with different block structures: 
// AR-CS, AR-AR,CS-AR,CS-CS

// [[Rcpp::export]]
double detR(mat& Rdraw, int row, int col, int change){
  
  mat Rtemp = Rdraw;
  
  Rtemp(row,col) = change;
  Rtemp(col,row) = change;
  
  return (det(Rtemp));
}


// [[Rcpp::export]]
List blockUpdate(mat Rdraw, double rhodraw, mat B, 
                 std::string structure,
                 int row, int col, double division){
  int nrow = B.n_rows, ncol = B.n_cols;
  
  double d1 = detR(Rdraw, row, col, 1),
    d2 = detR(Rdraw, row, col,-1),
    d0 = detR(Rdraw, row, col, 0);
  
  double a = 0.5*(d1+d2-2*d0), b = 0.5*(d1-d2), c = d0;
  
  NumericVector sol = realRoot(a,b,c);
  
  NumericVector diff = abs(rhodraw-sol);
  
  double mindiff = min(diff)/division;
  
  double rhonew = rhodraw + R::runif(-mindiff,mindiff);
  
  if(structure=="CS"){
    
    B.fill(rhonew);
    
  }
  
  else if(structure=="AR"){
    
    for(int i=0; i<nrow; i++){
      for(int j=i+1; j<ncol;j++){
        B(i,j) = pow(rhonew,abs(j-i));
        B(j,i) = B(i,j);
      }
    }   
  }
  
  else { 
    B.fill(rhonew); 
  }
  
  return List::create(Named("rhonew")=rhonew, Named("B")=B);
  
}


// [[Rcpp::export]]
double Rpost(colvec& y, mat& X, int k, colvec& betadraw, mat& Sdraw, mat& R){
  
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
List RUpdate(colvec& y, mat& X, int k, int k1, colvec& betadraw, 
             mat& Sdraw, mat& Rdraw, colvec& rhodraw, 
             List& BM, std::vector< std::string > structure,
             colvec division, colvec rhocount){
  
  int ncor = BM.size(); // number of blocks in block matrix
  
  colvec rowindex = as<arma::colvec>(NumericVector::create(0,k1,0));
  colvec colindex = as<arma::colvec>(NumericVector::create(1,k1+1,k1));
  
  mat Rnew = Rdraw;
  
  List out;
  
  double rhonew, Rnewpost, Rdrawpost, accprob;
  
  double u = log(R::runif(0,1));
  
  //colvec rhocount = zeros<colvec>(nblock);
  
  for(int i=0; i<ncor; i++){
    
    mat Bi = as<arma::mat>(BM[i]);
    
    out = blockUpdate(Rdraw, rhodraw(i), Bi, structure[i], 
                      rowindex(i), colindex(i), division(i));
    
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
    
    if(ispd(Rnew)){
      Rnewpost = Rpost(y, X, k, betadraw, Sdraw, Rnew);
    }else{
      Rnewpost = R_NegInf;
    }
    
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






// [[Rcpp::export]]
double logLik(colvec& y, mat& X, int k, int nblock, colvec& betadraw,
              colvec& sig2draw, mat& Sdraw, mat& Rdraw,
              colvec& B0, mat& sigB, colvec& s0, mat& sigS0){
  
  int n = X.n_rows, p = X.n_cols, l = n/k;
  
  double logy, logbeta, logsigma, logprior, logfull;
  
  
  mat invR = inv(Rdraw);
  mat invS = inv(Sdraw);
  mat invSig = invS*invR*invS;
  
  double empSS = 0;
  
  for (int i=0; i<n; i+=k){
    
    colvec res = y.rows(i,i+k-1)-X.rows(i,i+k-1)*betadraw;
    empSS += as_scalar(trans(res)*invSig*res);
    
  }
  
  // log likelihood of data y, pi = M_PI
  logy = as_scalar(-0.5*k*l*log(2*M_PI)-0.5*l*log(det(Sdraw*Rdraw*Sdraw))-0.5*empSS);
  //logy = as_scalar(-l*log(det(Sdraw))-(l/2)*log(det(Rdraw))-0.5*empSS);
  
 
  // log prior for beta
  logbeta = as_scalar(-0.5*p*log(2*M_PI)-0.5*log(det(sigB))-0.5*(trans(betadraw-B0)*inv(sigB)*(betadraw-B0)));
  //logbeta = as_scalar(-0.5*log(det(sigB))-0.5*(trans(betadraw-B0)*inv(sigB)*(betadraw-B0)));
  
  logsigma = as_scalar(-0.5*k*log(2*M_PI)-0.5*log(det(sigS0))-0.5*(trans(log(sig2draw)-s0))*inv(sigS0)*(log(sig2draw)-s0));
  logprior = logbeta + logsigma;
  
  logfull = logy + logprior;
  
  return logfull;
  
}




//[[Rcpp::export]]

colvec predUpdate(mat X, int k, colvec betadraw, mat Sdraw, mat Rdraw){
  
  int n = X.n_rows, l = n/k;
  
  //colvec predy = zeros<colvec>(n);
  
  colvec muhat = X*betadraw;
  
  mat cov = Sdraw*Rdraw*Sdraw;
  
  mat DiagM = zeros<mat>(l,l);
  DiagM.diag().ones();
  
  mat Omega = kron(DiagM, cov);
  
  
  
  return (trans(mvrnorm(1,muhat,Omega)));
  
}





// [[Rcpp::export]]
List mcmcUpdate(colvec& y, mat& X, List t, 
                std::vector< std::string > structure,
                colvec& B0, mat& sigB, colvec& s0, mat& sigS0, 
                colvec sigma_sd, colvec division, int niter){
  
  int  p = X.n_cols; //n = X.n_rows,
  int nblock = t.size();
  
  int k = 0;
  
  for(int nb=0; nb<nblock; nb++){
    k += as<NumericVector>(t[nb]).size();
  }
  
  int k1 = as<NumericVector>(t[0]).size();
  
  //IntegerVector ntemp = seq_len(nblock); 
  //int ncor = sum(ntemp);
  int nrho = nblock*(nblock+1)/2;
  
  colvec betadraw = zeros<vec>(p),
    sig2draw = zeros<vec>(nblock),
    rhodraw = zeros<vec>(nrho);

  
  mat Sdraw = zeros<mat>(k,k),
    Rdraw = zeros<mat>(k,k);
  
  
  double deviance;    
  
  List BM(nrho);
  
  
  
  mat betamat = zeros<mat>(niter,p),
      sig2mat = zeros<mat>(niter,nblock),
      rhomat = zeros<mat>(niter,nrho);
  colvec devmat = zeros<colvec>(niter);
  
 
  colvec sigcount = zeros<vec>(nblock),
    rhocount = zeros<vec>(nrho);
  

  List initialpars, Sout, Rout;
  
  
  // get initial values
  initialpars = initialGen(y, X, t);
  
  betadraw = as<arma::colvec>(initialpars["betadraw"]);
  sig2draw = as<arma::colvec>(initialpars["sig2draw"]);
  rhodraw = as<arma::colvec>(initialpars["rhodraw"]);
  Sdraw = as<arma::mat>(initialpars["Sdraw"]);
  Rdraw = as<arma::mat>(initialpars["Rdraw"]);
  BM = initialpars["BM"];
  
  
  // Update iterations
  
  for(int iter=0; iter<niter; iter++){
    
    betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);
    betamat.row(iter) = trans(betadraw);
    
    Sout = sigmaUpdate(y,X,t,betadraw,sig2draw,Sdraw,Rdraw, 
                       s0,sigS0,sigma_sd,sigcount);  
    sig2draw = as<colvec>(Sout["sig2draw"]);
    Sdraw = as<mat>(Sout["Sdraw"]);
    sigcount = as<mat>(Sout["sigcount"]);
    sig2mat.row(iter) = trans(sig2draw);
    
    Rout = RUpdate(y,X,k,k1,betadraw,Sdraw,Rdraw,rhodraw,
                   BM,structure,division,rhocount);
    rhodraw = as<colvec>(Rout["rhodraw"]);
    Rdraw = as<mat>(Rout["Rdraw"]);
    BM = Rout["BM"];
    rhocount = as<arma::colvec>(Rout["rhocount"]);
    rhomat.row(iter) = trans(rhodraw);
    
    
    deviance = -2*logLik(y,X,k,nblock,betadraw,sig2draw,Sdraw,
                         Rdraw,B0,sigB,s0,sigS0);
    devmat(iter) = deviance;

  }
  
  return List::create(Named("structure")=structure,
                      Named("sigcount")=sigcount,
                      Named("rhocount")=rhocount,
                      Named("betamat")=betamat,
                      Named("sig2mat")=sig2mat,
                      Named("rhomat")=rhomat,
                      Named("devmat")=devmat);
  
}




// Adaptive MCMC function for all pars

// [[Rcpp::export]]
List adaptiveUpdate(colvec& y, mat& X, List t,
                    std::vector< std::string > structure,
                    colvec& B0, mat& sigB, colvec& s0, mat& sigS0, 
                    colvec sigma_sd, colvec division,               
                    int nadap, int niter){
  //n = X.n_rows,
  int  p = X.n_cols;
  //int nkeep = niter-nburn;
  
  int nblock = t.size();
  
  int k = 0;
  
  for(int nb=0; nb<nblock; nb++){
    k += as<NumericVector>(t[nb]).size();
  }
  
  int k1 = as<NumericVector>(t[0]).size();
  
  //IntegerVector ntemp = seq_len(nblock); 
  //int ncor = sum(ntemp);
  int nrho = nblock*(nblock+1)/2;
  
  const double lowaccrate = 0.2, 
    upaccrate = 0.45;
  double lowacc_count = lowaccrate*nadap, 
    upacc_count = upaccrate*nadap;
  
  
  //colvec betadraw, sig2draw, rhodraw;
  colvec betadraw = zeros<vec>(p),
    sig2draw = zeros<vec>(nblock),
    rhodraw = zeros<vec>(nrho),
    predydraw;
  
  //mat Sdraw, Rdraw;
  mat Sdraw = zeros<mat>(k,k),
    Rdraw = zeros<mat>(k,k);
  
  double deviance;    
  
  //List BM;
  List BM(nrho);
  
  
  mat betamat = zeros<mat>(niter,p),
    sig2mat = zeros<mat>(niter,nblock),
    rhomat = zeros<mat>(niter,nrho);
  colvec devmat = zeros<mat>(niter);

  
  colvec sigcount = zeros<vec>(nblock),
         rhocount = zeros<vec>(nrho);
  
  //colvec sigcount, rhocount;
  
  NumericVector count = wrap(join_cols(sigcount,rhocount));
  
  List initialpars, Sout, Rout;
  
  //double lowgamma, lowalpha, lowdivision; 
  
  
  // get initial values
  initialpars = initialGen(y, X, t);
  
  betadraw = as<arma::colvec>(initialpars["betadraw"]);
  sig2draw = as<arma::colvec>(initialpars["sig2draw"]);
  rhodraw = as<arma::colvec>(initialpars["rhodraw"]);
  Sdraw = as<arma::mat>(initialpars["Sdraw"]);
  Rdraw = as<arma::mat>(initialpars["Rdraw"]);
  BM = initialpars["BM"];
  
  
  // adaptive MH phase
  
  double testno = 0;
  while(is_true(any(count<lowacc_count)) || is_true(any(count>upacc_count))){  
    
    testno += 1;
    cout<<"test "<<testno<<endl;
    
    for(int iadap=0; iadap<nadap; iadap++){
      
      betadraw = betaUpdate(y,X,k,Sdraw,Rdraw,B0,sigB);
      
      Sout = sigmaUpdate(y,X,t,betadraw,sig2draw,Sdraw,Rdraw, 
                         s0,sigS0,sigma_sd,sigcount);  
      sig2draw = as<colvec>(Sout["sig2draw"]);
      Sdraw = as<mat>(Sout["Sdraw"]);
      sigcount = as<mat>(Sout["sigcount"]);
     
      Rout = RUpdate(y,X,k,k1,betadraw,Sdraw,Rdraw,rhodraw,
                     BM,structure,division,rhocount);
      rhodraw = as<colvec>(Rout["rhodraw"]);
      Rdraw = as<mat>(Rout["Rdraw"]);
      BM = Rout["BM"];
      rhocount = as<arma::colvec>(Rout["rhocount"]);
    
    }
    
    count = wrap(join_cols(sigcount,rhocount));
    
    
    for(int j=0; j<nblock; j++){
      
      
      
      if(sigcount(j)<lowacc_count && sigma_sd(j)>(0.05*sigma_sd(j))){
        sigma_sd(j) -= 0.05*sigma_sd(j);
      }
      
      
      if(sigcount(j)>upacc_count){
        sigma_sd(j) += 0.05*sigma_sd(j);
      }
      
    }
    
    for(int k=0; k<nrho;k++){
      
      if(rhocount(k)<lowacc_count){
        division(k) += 0.1*division(k);
      }
      
      if(rhocount(k)>upacc_count && division(k)>(0.1*division(k))){
        division(k) -= 0.1*division(k);
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
    
    Sout = sigmaUpdate(y,X,t,betadraw,sig2draw,Sdraw,Rdraw, 
                       s0,sigS0,sigma_sd,sigcount);  
    sig2draw = as<colvec>(Sout["sig2draw"]);
    Sdraw = as<mat>(Sout["Sdraw"]);
    sigcount = as<mat>(Sout["sigcount"]);
    sig2mat.row(iter) = trans(sig2draw);
    
    Rout = RUpdate(y,X,k,k1,betadraw,Sdraw,Rdraw,rhodraw,
                   BM,structure,division,rhocount);
    rhodraw = as<colvec>(Rout["rhodraw"]);
    Rdraw = as<mat>(Rout["Rdraw"]);
    BM = Rout["BM"];
    rhocount = as<arma::colvec>(Rout["rhocount"]);
    rhomat.row(iter) = trans(rhodraw);
    
    
    deviance = -2*logLik(y,X,k,nblock,betadraw,sig2draw,Sdraw,
                         Rdraw,B0,sigB,s0,sigS0);
    devmat(iter) = deviance;
  }
  
  return List::create(Named("structure")=structure,
                      Named("sigcount")=sigcount,
                      Named("rhocount")=rhocount,
                      Named("sigma_sd")=sigma_sd,
                      Named("division")=division,
                      Named("betamat")=betamat,
                      Named("sig2mat")=sig2mat,
                      Named("rhomat")=rhomat,
                      Named("devmat")=devmat);
  
}
