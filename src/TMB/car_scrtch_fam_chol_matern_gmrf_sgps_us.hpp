#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

//namespace car_gmrf {
//using namespace Eigen;
//using namespace tmbutils;
//
//
//template<class Type>
//struct car_gmrf_t{  
//  SparseMatrix<Type> diagI;        // G0 eqn (10) in Lindgren 
//  SparseMatrix<Type> N;        // G1 eqn (10) in Lindgren 
//  SparseMatrix<Type> Mpm5;        // G2 eqn (10) in Lindgren 
//  car_gmrf_t(SEXP x){  /* x = List passed from R */
//  diagI = asSparseMatrix<Type>(getListElement(x,"diagI"));
//  N = asSparseMatrix<Type>(getListElement(x,"N"));
//  Mpm5 = asSparseMatrix<Type>(getListElement(x,"Mpm5"));
//}
//};
//
//template<class Type> 
//  SparseMatrix<Type> Q_car(car_gmrf_t<Type> car_mats,Type rho, Type s2_s){
//	  //See Proposition 1 J.M.VerHoefetal./SpatialStatistics25(2018)68–8
//    Type tau2=Type(1.0) / s2_s;
//  return tau2*(car_mats.Mpm5*(car_mats.diagI - rho*car_mats.N)*car_mats.Mpm5); 
//  }
//};

template<class Type>
  Type car_scrtch_fam_chol_matern_gmrf_sgps_us(objective_function<Type>* obj)
{
//  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace car_gmrf; // 
  using namespace density; 
//  using namespace distributions_R;
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_INTEGER( binary_ind );
  DATA_VECTOR( y );                
  //DATA_VECTOR( cens );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SCALAR(rhomax);
  DATA_SPARSE_MATRIX(Zc_par);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zc_sib);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zs);         //Design matrix for spatial random effects 
  DATA_SPARSE_MATRIX(Qa);         //Design matrix for genetic random effects 
  DATA_STRUCT(car_mats,car_gmrf_t); //
  //DATA_INTEGER(flag);  
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_s );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_c_par );
  PARAMETER( log_sdvc_c_sib );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( uasc );
  PARAMETER_VECTOR( uc_par );
  PARAMETER_VECTOR( uc_sib );
  
  PARAMETER_VECTOR(beta);  
  PARAMETER(logit_rho);
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  Type rho = exp(logit_rho)/(1+exp(logit_rho))*rhomax;
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_s = pow(exp(log_sdvc_s),2);
  Type vc_c_par = pow(exp(log_sdvc_c_par),2);
  Type vc_c_sib = pow(exp(log_sdvc_c_sib),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  // Spatial interpolation
  //matrix<Type> AA = Zs*A; 
  //vector<Type> delta = (AA*x)/tau;
  vector<Type> delta = Zs*x;
  vector<Type> ua = uasc*exp(log_sdvc_a) ;
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q = Q_car(car_mats,rho,vc_s);
  //---------------------------------------------
  
  vector<Type> uc_par_l = Zc_par*uc_par ;
  vector<Type> uc_sib_l = Zc_sib*uc_sib ;
  vector<Type> eta = X*beta + delta + ua + uc_par_l + uc_sib_l;

  //for( int j=0; j< Lt_Ga.cols(); j++){
  //  nll -= dnorm( ua(j) , Type(0.0), exp(log_sdvc_a), true );
  //} 
  
  //Calculates nll-------------------------------
  Type nll = 0.0;

  //nll = GMRF(Q, false)(x);
  //if (flag == 0)
  //  return nll; 
  nll+= GMRF(Q)(x);
  nll+=GMRF(Qa)(uasc); 
 
	  //nll = GMRF(Q, false)(x);  
  // Return un-normalized density on request
  //if (flag == 0) return nll;
  
  for( int j=0; j< Zc_par.cols(); j++){
    nll -= dnorm( uc_par(j) , Type(0.0), exp(log_sdvc_c_par), true );
  } 
  for( int j=0; j< Zc_sib.cols(); j++){
    nll -= dnorm( uc_sib(j) , Type(0.0), exp(log_sdvc_c_sib), true );
  } 
  for( int j=0; j< X.rows(); j++){
    nll -= dnorm( y(j) , eta(j), exp(log_sdvc_res), true );
  } 
  //car_mats,rho,vc_s
//  for( int j=0; j< X.rows(); j++){
//    nll -= cens(j) * -log(1 + pow(eta(j)*y(j),shp)) + //(S(y)
//	    (eta(j)*shp)*pow(eta(j)*y(j),shp-1.0) - log(1 + pow(eta(j)*y(j),shp)) ;//h(y)
//  } 
//
//  for( int j=0; j< X.rows(); j++){
//    nll -=  -log(1 + pow(eta(j)*y(j),shp)) + //(S(y)
//	    cens(j) * (eta(j)*shp)*pow(eta(j)*y(j),shp-1.0) - 
//	    cens(j)* log(1 + pow(eta(j)*y(j),shp)) ;//h(y)
//  } 
  //---------------------------------------------
  
  //Report what we want to report----------------
  //Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011);
  //Type vc_s = pow(exp(.5*log(Type(1.0)/(Type(4.0)*M_PI)) - log(kappa) - log_tau),2);
  //ADREPORT(range);
  //

  if( binary_ind == 1){

  Type vc_pp = vc_c_par; 
  Type vc_ps = .5*vc_a;
  Type vc_ss = .5*vc_a + vc_c_sib; 
  Type khat = y.mean();
  Type t=qnorm(1-khat);
  Type z=dnorm(t, Type(0.0), Type(1.0));
  Type i=z/khat;
  
  Type Eb_pp=khat + vc_pp/khat;
  Type Eb_ps=khat + vc_ps/khat;
  Type Eb_ss=khat + vc_ss/khat;

  Type tvc_pp=qnorm(1-Eb_pp);
  Type tvc_ps=qnorm(1-Eb_ps);
  Type tvc_ss=qnorm(1-Eb_ss);

  Type vc_pp_adj = (t - tvc_pp*sqrt(1-(pow(t,2.0) - pow(tvc_pp,2.0))*(1-t/i))) / (i + pow(tvc_pp,2.0)*(i-t));
  Type vc_ps_adj = (t - tvc_ps*sqrt(1-(pow(t,2.0) - pow(tvc_ps,2.0))*(1-t/i))) / (i + pow(tvc_ps,2.0)*(i-t));
  Type vc_ss_adj = (t - tvc_ss*sqrt(1-(pow(t,2.0) - pow(tvc_ss,2.0))*(1-t/i))) / (i + pow(tvc_ss,2.0)*(i-t));

  Type vc_a_lia=2*vc_ps_adj; 
  Type vc_c_par_lia=vc_pp_adj; 
  Type vc_c_sib_lia=2*(vc_ss_adj - vc_ps_adj); 

  ADREPORT(vc_a_lia);
  ADREPORT(vc_c_par_lia);
  ADREPORT(vc_c_sib_lia);
  ADREPORT(khat);
  }
  
  ADREPORT(vc_a);
  ADREPORT(vc_c_par);
  ADREPORT(vc_c_sib);
  ADREPORT(vc_res);
  ADREPORT(vc_s);
  ADREPORT(rho);
  //ADREPORT(shp);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this



//
/* List of sparse matrices */
//template<class Type>
//struct losm_t : vector<SparseMatrix<Type> > {
//  losm_t(SEXP x){  /* x = List passed from R */
//    (*this).resize(LENGTH(x));
//    for(int i=0; i<LENGTH(x); i++){
//      SEXP sm = VECTOR_ELT(x, i);
//      if(!isValidSparseMatrix(sm))
//            error("Not a sparse matrix");
//      (*this)(i) = asSparseMatrix<Type>(sm);
//    }
//  }
//};
