#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
  Type fam_chol_matern_gmrf_gp_us(objective_function<Type>* obj)
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SPARSE_MATRIX(Zc_par);         //Design matrix for genetic random effects 
  DATA_SPARSE_MATRIX(Lt_Ga);         //Design matrix for genetic random effects 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_a );
  PARAMETER( log_sdvc_c_par );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( ua );
  PARAMETER_VECTOR( uc_par );
  
  PARAMETER_VECTOR(beta);  
  //---------------------------------------
  
  //Transform parameters-------------------
  vector<Type> Lua = Lt_Ga * ua; 
  Type vc_a = pow(exp(log_sdvc_a),2);
  Type vc_c_par = pow(exp(log_sdvc_c_par),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  
  
  //Calculates nll-------------------------------
  Type nll = 0.0;

  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> uc_par_l = Zc_par*uc_par ;
  vector<Type> eta = X*beta + Lua + uc_par_l ;

  for( int j=0; j< Lt_Ga.cols(); j++){
    nll -= dnorm( ua(j) , Type(0.0), exp(log_sdvc_a), true );
  } 
  for( int j=0; j< Zc_par.cols(); j++){
    nll -= dnorm( uc_par(j) , Type(0.0), exp(log_sdvc_c_par), true );
  } 
  for( int j=0; j< X.rows(); j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------

  Type vc_pp = vc_c_par; 
  Type vc_ps = .5*vc_a;

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
  
  Type vc_a_lia=2*vc_ps_adj; 
  Type vc_c_par_lia=vc_pp_adj; 
  
  ADREPORT(vc_a_lia);
  ADREPORT(vc_c_par_lia);
  ADREPORT(khat);


  //Report what we want to report----------------
  ADREPORT(vc_a);
  ADREPORT(vc_c_par);
  ADREPORT(vc_res);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
