#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
  Type fam_chol_matern_gmrf_sfs_us(objective_function<Type>* obj)
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SPARSE_MATRIX(Zc_fam);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zc_sib);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zs);         //Design matrix for spatial random effects 
  DATA_SPARSE_MATRIX(Lt_Gc_fam);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Lt_Gc_sib);         //Design matrix for genetic random effects 
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating points witin triangles 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_c_fam );
  PARAMETER( log_sdvc_c_sib );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( uc_fam );
  PARAMETER_VECTOR( uc_sib );
  
  PARAMETER_VECTOR(beta);  
  PARAMETER(log_tau);
  PARAMETER(log_kappa);
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lucfam = Lt_Gc_fam * uc_fam; 
  vector<Type> Lucsib = Lt_Gc_sib * uc_sib; 
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type vc_c_fam = pow(exp(log_sdvc_c_fam),2);
  Type vc_c_sib = pow(exp(log_sdvc_c_sib),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  // Spatial interpolation
  matrix<Type> AA = Zs*A; 
  vector<Type> delta = (AA*x)/tau;
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q = Q_spde(spdeMatrices,kappa);
  //---------------------------------------------
  
  //Calculates nll-------------------------------
  Type nll = 0.0;
  nll = GMRF(Q)(x);                              

  //nll = GMRF(Q, false)(x);  
  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> uc_fam_l = Zc_fam*Lucfam ;
  vector<Type> uc_sib_l = Zc_sib*Lucsib ;
  vector<Type> eta = X*beta + delta + uc_fam_l + uc_sib_l;

  for( int j=0; j< Zc_fam.cols(); j++){
    nll -= dnorm( uc_fam(j) , Type(0.0), exp(log_sdvc_c_fam), true );
  } 
  for( int j=0; j< Zc_sib.cols(); j++){
    nll -= dnorm( uc_sib(j) , Type(0.0), exp(log_sdvc_c_sib), true );
  } 
  for( int j=0; j< X.rows(); j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------
  
  //Report what we want to report----------------
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011);
  Type vc_s = pow(exp(.5*log(Type(1.0)/(Type(4.0)*M_PI)) - log_kappa - log_tau),2);
  ADREPORT(range);
  ADREPORT(vc_c_fam);
  ADREPORT(vc_c_sib);
  ADREPORT(vc_res);
  ADREPORT(vc_s);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
