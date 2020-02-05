#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
  Type fam_chol_car_sgps_us(objective_function<Type>* obj)
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         // Design matrix for fixed effects
  DATA_MATRIX(Vw);        // eigenvectors of neighborhood (Weight) matrix W
  DATA_VECTOR(wj);        // eigenvalues of neighborhood (Weight) matrix W
  DATA_SCALAR(minrho);    // 1/min(eigenvalue)
  DATA_SCALAR(maxrho);    // 1/max(eigenvalue)
  DATA_SPARSE_MATRIX(Zc_par);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zc_sib);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zs);         //Design matrix for spatial random effects 
  DATA_SPARSE_MATRIX(Lt_Ga);         //Design matrix for genetic random effects 

  // Parameters
  //PARAMETER( mu );
  PARAMETER( unscaled_rhocar );
  PARAMETER( log_sdvc_s );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_c_par );
  PARAMETER( log_sdvc_c_sib );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( us );
  PARAMETER_VECTOR( ua );
  PARAMETER_VECTOR( uc_par );
  PARAMETER_VECTOR( uc_sib );
  
  PARAMETER_VECTOR(beta);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lus = Vw * us; 
  vector<Type> Lua = Lt_Ga * ua; 
  Type sd_s = exp(log_sdvc_s);
  Type vc_s = pow(sd_s,2);
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_c_par = pow(exp(log_sdvc_c_par),2);
  Type vc_c_sib = pow(exp(log_sdvc_c_sib),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  //Calculates nll-------------------------------
  Type nll = 0.0;
  
  vector<Type> usl = Zs*Lus ;
  vector<Type> uc_par_l = Zc_par*uc_par ;
  vector<Type> uc_sib_l = Zc_sib*uc_sib ;
  vector<Type> eta = X*beta + usl + Lua + uc_par_l + uc_sib_l ;
  
  Type rhocar = (maxrho-minrho)*(exp(unscaled_rhocar)/(1+exp(unscaled_rhocar))) + minrho;
  
  for( int j=0; j< Zs.cols(); j++){
    nll -= dnorm( us(j) , Type(0.0), sd_s/sqrt(1-rhocar*wj(j)), true );
  } 
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
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------
  
  //Report what we want to report----------------
  ADREPORT(rhocar);
  ADREPORT(vc_a);
  ADREPORT(vc_c_par);
  ADREPORT(vc_c_sib);
  ADREPORT(vc_res);
  ADREPORT(vc_s);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
