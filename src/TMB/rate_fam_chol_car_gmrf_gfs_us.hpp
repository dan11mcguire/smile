#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj


template<class Type>
  Type ratescrtch_fam_chol_car_gmrf_gfs_us(objective_function<Type>* obj)
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
  DATA_SPARSE_MATRIX(Zc_fam);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zc_sib);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Lt_Ga);         //Design matrix for genetic random effects 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_phi );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_c_fam );
  PARAMETER( log_sdvc_c_sib );
  PARAMETER_VECTOR( ua );
  PARAMETER_VECTOR( uc_fam );
  PARAMETER_VECTOR( uc_sib );
  
  PARAMETER_VECTOR(beta);  
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lua = Lt_Ga * ua; 
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_c_fam = pow(exp(log_sdvc_c_fam),2);
  Type vc_c_sib = pow(exp(log_sdvc_c_sib),2);
  Type phi = exp(log_phi);
  //------------------------------------------

  
  //Calculates nll-------------------------------
  Type nll = 0.0;

  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> uc_fam_l = Zc_fam*uc_fam ;
  vector<Type> uc_sib_l = Zc_sib*uc_sib ;
  vector<Type> eta = X*beta + Lua + uc_fam_l + uc_sib_l;
  vector<Type> muhat = eta + offset;

  for( int j=0; j< Lt_Ga.cols(); j++){
    nll -= dnorm( ua(j) , Type(0.0), exp(log_sdvc_a), true );
  } 
  for( int j=0; j< Zc_fam.cols(); j++){
    nll -= dnorm( uc_fam(j) , Type(0.0), exp(log_sdvc_c_fam), true );
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
  ADREPORT(vc_c_fam);
  ADREPORT(vc_c_sib);
  ADREPORT(phi);
  ADREPORT(rho);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
