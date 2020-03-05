#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

//namespace R_inla_generalized {
//using namespace Eigen;
//using namespace tmbutils;
//using namespace R_inla;
//
//template<class Type>
//  SparseMatrix<Type> Q_spde_generalized(spde_t<Type> spde, Type kappa, int alpha=2){
//  Type kappa_pow2 = kappa*kappa;
//  Type kappa_pow4 = kappa_pow2*kappa_pow2;
//  	
//  if( alpha==1 ) return kappa_pow2*spde.M0 + spde.M1;
//  if( alpha==2 ) return kappa_pow4*spde.M0 + Type(2.0)*kappa_pow2*spde.M1 + spde.M2;
//}
//}

template<class Type>
  Type ratescrtch_fam_chol_car_gmrf_sgps_us(objective_function<Type>* obj)
{
  using namespace car_gmrf; //
//  using namespace R_inla_generalized; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  DATA_VECTOR( offset );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SCALAR(rhomax);
  DATA_SPARSE_MATRIX(Zc_par);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zc_sib);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zs);         //Design matrix for spatial random effects 
  DATA_SPARSE_MATRIX(Lt_Ga);         //Design matrix for genetic random effects 
  DATA_STRUCT(car_mats,car_gmrf_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating points witin triangles 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_phi );
  PARAMETER( logit_rho );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_s );
  PARAMETER( log_sdvc_c_par );
  PARAMETER( log_sdvc_c_sib );
//  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( ua );
  PARAMETER_VECTOR( uc_par );
  PARAMETER_VECTOR( uc_sib );
  
  PARAMETER_VECTOR(beta);  
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lua = Lt_Ga * ua; 
  Type rho = exp(logit_rho)/(1+exp(logit_rho))*rhomax;
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_s = pow(exp(log_sdvc_s),2);
  Type vc_c_par = pow(exp(log_sdvc_c_par),2);
  Type vc_c_sib = pow(exp(log_sdvc_c_sib),2);
//  Type vc_res = pow(exp(log_sdvc_res),2);
  Type phi = exp(log_phi);
  //------------------------------------------

  // Spatial interpolation
  vector<Type> delta = Zs*x; 
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q = Q_car(car_mats,rho,vc_s);
  //---------------------------------------------
  
  //Calculates nll-------------------------------
  Type nll = 0.0;
  nll = GMRF(Q)(x);                              

  //nll = GMRF(Q, false)(x);  
  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> uc_par_l = Zc_par*uc_par ;
  vector<Type> uc_sib_l = Zc_sib*uc_sib ;
  vector<Type> eta = X*beta + delta + Lua + uc_par_l + uc_sib_l;
  vector<Type> muhat = eta + offset;

  for( int j=0; j< Lt_Ga.cols(); j++){
    nll -= dnorm( ua(j) , Type(0.0), exp(log_sdvc_a), true );
  } 
  for( int j=0; j< Zc_par.cols(); j++){
    nll -= dnorm( uc_par(j) , Type(0.0), exp(log_sdvc_c_par), true );
  } 
  for( int j=0; j< Zc_sib.cols(); j++){
    nll -= dnorm( uc_sib(j) , Type(0.0), exp(log_sdvc_c_sib), true );
  } 

  for( int j=0; j< X.rows(); j++){
    nll -=  dnbinom_robust(y(j),
		    muhat(j),
		    Type(2.0)*muhat(j) + log_phi,
		    true) ;//h(y)
  } 
  //---------------------------------------------
  
  //Report what we want to report----------------
  ADREPORT(vc_a);
  ADREPORT(vc_c_par);
  ADREPORT(vc_c_sib);
  ADREPORT(vc_s);
  ADREPORT(phi);
  ADREPORT(rho);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
