#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
  Type fam_chol_car_sf_us(objective_function<Type>* obj)
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         // Design matrix for fixed effects
  DATA_SPARSE_MATRIX(Vw);        // eigenvectors of neighborhood (Weight) matrix W
  DATA_VECTOR(wj);        // eigenvalues of neighborhood (Weight) matrix W
  DATA_SCALAR(minrho);    // 1/min(eigenvalue)
  DATA_SCALAR(maxrho);    // 1/max(eigenvalue)
  DATA_SPARSE_MATRIX(Zc_fam);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zs);         //Design matrix for spatial random effects 

  // Parameters
  //PARAMETER( mu );
  PARAMETER( unscaled_rhocar );
  PARAMETER( log_sdvc_s );
  PARAMETER( log_sdvc_c_fam );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( us );
  PARAMETER_VECTOR( uc_fam );
  
  PARAMETER_VECTOR(beta);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lus = Vw * us; 
  Type sd_s = exp(log_sdvc_s);
  Type vc_s = pow(sd_s,2);
  Type vc_c_fam = pow(exp(log_sdvc_c_fam),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  //Calculates nll-------------------------------
  Type nll = 0.0;
  
  vector<Type> usl = Zs*Lus ;
  vector<Type> uc_fam_l = Zc_fam*uc_fam ;
  vector<Type> eta = X*beta + usl  + uc_fam_l ;
  
  Type rhocar = (maxrho-minrho)*(exp(unscaled_rhocar)/(1+exp(unscaled_rhocar))) + minrho;
  
  for( int j=0; j< Zs.cols(); j++){
    nll -= dnorm( us(j) , Type(0.0), sd_s/sqrt(1-rhocar*wj(j)), true );
  } 
  for( int j=0; j< Zc_fam.cols(); j++){
    nll -= dnorm( uc_fam(j) , Type(0.0), exp(log_sdvc_c_fam), true );
  } 
  for( int j=0; j< X.rows(); j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------
  
  //Report what we want to report----------------
  ADREPORT(rhocar);
  ADREPORT(vc_c_fam);
  ADREPORT(vc_res);
  ADREPORT(vc_s);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
