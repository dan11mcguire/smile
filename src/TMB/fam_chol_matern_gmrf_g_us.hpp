#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
  Type fam_chol_matern_gmrf_g_us(objective_function<Type>* obj)
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SPARSE_MATRIX(Lt_Ga);         //Design matrix for genetic random effects 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( ua );
  
  PARAMETER_VECTOR(beta);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lua = Lt_Ga * ua; 
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  
  
  //Calculates nll-------------------------------
  Type nll = 0.0;

  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> eta = X*beta + Lua ;

  for( int j=0; j< Lt_Ga.cols(); j++){
    nll -= dnorm( ua(j) , Type(0.0), exp(log_sdvc_a), true );
  } 
  for( int j=0; j< X.rows(); j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------
  
  //Report what we want to report----------------
  ADREPORT(vc_a);
  ADREPORT(vc_res);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
