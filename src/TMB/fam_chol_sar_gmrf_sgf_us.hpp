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
  Type fam_chol_sar_gmrf_sgf_us(objective_function<Type>* obj)
{
  using namespace car_gmrf; //
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_INTEGER( binary_ind );
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SCALAR(rhomax);
  DATA_SPARSE_MATRIX(Zc_fam);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zs);         //Design matrix for spatial random effects 
  DATA_SPARSE_MATRIX(Lt_Ga);         //Design matrix for genetic random effects 
  DATA_STRUCT(car_mats,car_gmrf_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_s );
  PARAMETER( log_sdvc_c_fam );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( ua );
  PARAMETER_VECTOR( uc_fam );
  
  PARAMETER_VECTOR(beta);  
  PARAMETER(logit_rho);
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lua = Lt_Ga * ua; 
  Type rho = exp(logit_rho)/(1+exp(logit_rho))*rhomax;
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_s = pow(exp(log_sdvc_s),2);
  Type vc_c_fam = pow(exp(log_sdvc_c_fam),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  // Spatial interpolation
  vector<Type> delta = Zs*x*sqrt(vc_s);
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q = Q_sar(car_mats,rho,vc_s);
  //---------------------------------------------
  
  //Calculates nll-------------------------------
  Type nll = 0.0;
  nll = GMRF(Q)(x);                              

  //nll = GMRF(Q, false)(x);  
  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> uc_fam_l = Zc_fam*uc_fam ;
  vector<Type> eta = X*beta + delta + Lua + uc_fam_l ;

  for( int j=0; j< Lt_Ga.cols(); j++){
    nll -= dnorm( ua(j) , Type(0.0), exp(log_sdvc_a), true );
  } 
  for( int j=0; j< Zc_fam.cols(); j++){
    nll -= dnorm( uc_fam(j) , Type(0.0), exp(log_sdvc_c_fam), true );
  } 
  for( int j=0; j< X.rows(); j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------

  if( binary_ind == 1){

  Type vc_pp = vc_c_fam; 
  Type vc_ps = .5*vc_a + vc_c_fam;

  Type khat = y.mean();
  Type t=qnorm(1-khat);
  Type z=dnorm(t, Type(0.0), Type(1.0));
  Type i=z/khat;
  
  Type Eb_pp=khat + vc_pp/khat;
  Type Eb_ps=khat + vc_ps/khat;

  Type tvc_pp=qnorm(1-Eb_pp);
  Type tvc_ps=qnorm(1-Eb_ps);

  Type vc_pp_adj = (t - tvc_pp*sqrt(1-(pow(t,2.0) - pow(tvc_pp,2.0))*(1-t/i))) / (i + pow(tvc_pp,2.0)*(i-t));
  Type vc_ps_adj = (t - tvc_ps*sqrt(1-(pow(t,2.0) - pow(tvc_ps,2.0))*(1-t/i))) / (i + pow(tvc_ps,2.0)*(i-t));
  
  Type vc_a_lia=2*(vc_ps_adj - vc_pp_adj); 
  Type vc_c_fam_lia=vc_pp_adj; 
  
  ADREPORT(vc_a_lia);
  ADREPORT(vc_c_fam_lia);
  ADREPORT(khat);
  
  }

  //Report what we want to report----------------
  ADREPORT(vc_a);
  ADREPORT(vc_c_fam);
  ADREPORT(vc_res);
  ADREPORT(vc_s);
  ADREPORT(rho);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
