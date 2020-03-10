#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
  Type fam_chol_matern_gmrf_gfs_us(objective_function<Type>* obj)
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  DATA_INTEGER( binary_ind );

  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SPARSE_MATRIX(Zc_fam);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Zc_sib);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Lt_Ga);         //Design matrix for genetic random effects 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_c_fam );
  PARAMETER( log_sdvc_c_sib );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( ua );
  PARAMETER_VECTOR( uc_fam );
  PARAMETER_VECTOR( uc_sib );
  
  PARAMETER_VECTOR(beta);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lua = Lt_Ga * ua; 
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_c_fam = pow(exp(log_sdvc_c_fam),2);
  Type vc_c_sib = pow(exp(log_sdvc_c_sib),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  
  
  //Calculates nll-------------------------------
  Type nll = 0.0;

  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> uc_fam_l = Zc_fam*uc_fam ;
  vector<Type> uc_sib_l = Zc_sib*uc_sib ;
  vector<Type> eta = X*beta + Lua + uc_fam_l + uc_sib_l;

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
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------
  if( binary_ind == 1){

  Type vc_pp = vc_c_fam; 
  Type vc_ps = .5*vc_a + vc_c_fam;
  Type vc_ss = .5*vc_a + vc_c_fam + vc_c_sib; 

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
  
  Type vc_a_lia=2*(vc_ps_adj - vc_pp_adj); 
  Type vc_c_fam_lia=vc_pp_adj; 
  Type vc_c_sib_lia=vc_ss_adj - vc_ps_adj; 
  
  ADREPORT(vc_a_lia);
  ADREPORT(vc_c_fam_lia);
  ADREPORT(vc_c_sib_lia);
  ADREPORT(khat);
  }

  //Report what we want to report----------------
  ADREPORT(vc_a);
  ADREPORT(vc_c_fam);
  ADREPORT(vc_c_sib);
  ADREPORT(vc_res);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
