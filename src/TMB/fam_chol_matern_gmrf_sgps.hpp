#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
	
template<class Type>
  Type fam_chol_matern_gmrf_sgps(objective_function<Type>* obj)
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
//  DATA_SPARSE_MATRIX(Za);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zc_fam);         //Design matrix for family-shared random effects 
  DATA_SPARSE_MATRIX(Zc_par);         //Design mat parent-within-family shared random effects 
  DATA_SPARSE_MATRIX(Zc_sib);         //Design mat sydep-within-family shared random effects 
  DATA_SPARSE_MATRIX(Lt_Ga);         // cholesky decomposition of relationship matrix 
 // DATA_SPARSE_MATRIX(Lt_Gc_par);         // not realy needed (identity matrix)
 // DATA_SPARSE_MATRIX(Lt_Gc_sib);         // not realy needed (identity matrix) 
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating points witin triangles 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_c_par );
  PARAMETER( log_sdvc_c_sib );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( ua );
  PARAMETER_VECTOR( uc_par );
  PARAMETER_VECTOR( uc_sib );
  
  PARAMETER_VECTOR(beta);  
  PARAMETER(log_tau);
  PARAMETER(log_kappa);
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lua = Lt_Ga * ua; 
 // vector<Type> Lucpar = Lt_Gc_par * uc_par; 
 // vector<Type> Lucsib = Lt_Gc_sib * uc_sib; 
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_c_par = pow(exp(log_sdvc_c_par),2);
  Type vc_c_sib = pow(exp(log_sdvc_c_sib),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  // Spatial interpolation
  matrix<Type> AA = Zc_fam*A; 
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

//  vector<Type> ual = Za*Lua ;
  vector<Type> uc_par_l = Zc_par*uc_par ;
  vector<Type> uc_sib_l = Zc_sib*uc_sib ;
  vector<Type> eta = X*beta + delta + Lua + uc_par_l + uc_sib_l;

  for( int j=0; j< Lua.cols(); j++){
    nll -= dnorm( ua(j) , Type(0.0), exp(log_sdvc_a), true );
  } 
  for( int j=0; j< Zc_par.cols(); j++){
    nll -= dnorm( uc_par(j) , Type(0.0), exp(log_sdvc_c_par), true );
  } 
  for( int j=0; j< Zc_sib.cols(); j++){
    nll -= dnorm( uc_sib(j) , Type(0.0), exp(log_sdvc_c_sib), true );
  } 
  for( int j=0; j< Lua.rows(); j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------
  
  //Report what we want to report----------------
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011);
  Type vc_s = pow(exp(.5*log(Type(1.0)/(Type(4.0)*M_PI)) - log_kappa - log_tau),2);
  ADREPORT(vc_s);

//  Type twin_vc_a = 2*vc_m;
//  Type twin_vc_c = vc_pair - vc_m ;
//  Type vc_res = pow(exp(log_sdvc_res),2);
//  Type f_n_mz = Type(npair_mz);
//  Type f_n_dz = Type(npair_dz);
//  Type p_dz = Type(2.0)*f_n_dz / (f_n_dz + f_n_mz);    // if dz => os and mz => ss, calculate estimated vcs
//  Type p_mz_given_ss = (1 - p_dz) / (f_n_mz / (f_n_mz + f_n_dz));
//  Type sex_vc_a = 2/p_mz_given_ss*vc_m;
//  //Type sex_vc_c = ((p_mz_given_ss + 1)*vc_pair - vc_m) / p_mz_given_ss ;
//  Type sex_vc_c = vc_pair - vc_m / p_mz_given_ss ;
//  Type vc_a=0.0 ;
//  Type vc_c=0.0 ;
// 
//   if(is_twin_model == 1){
//    vc_a=twin_vc_a;
//    vc_c=twin_vc_c;
//  } else {
//    vc_a=sex_vc_a;
//    vc_c=sex_vc_c;
//  }
  ADREPORT(range);
  ADREPORT(vc_a);
  ADREPORT(vc_c_par);
  ADREPORT(vc_c_sib);
  ADREPORT(vc_res);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
