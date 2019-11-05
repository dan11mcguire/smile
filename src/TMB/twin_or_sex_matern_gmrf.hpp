#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_INTEGER( is_twin_model);
  DATA_INTEGER( npair_dz );
  DATA_INTEGER( npair_mz );
  //DATA_FACTOR( pair_id );
  //DATA_FACTOR( m_id );                
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SPARSE_MATRIX(Zp);         //Design matrix for sy-pair random effects 
  DATA_SPARSE_MATRIX(Zm);         //Design matrix for sy-sex or mz random effects 
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating points witin triangles 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_pair );
  PARAMETER( log_sdvc_m );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( up );
  PARAMETER_VECTOR( um );
  
  
  PARAMETER_VECTOR(beta);  
  PARAMETER(log_tau);
  PARAMETER(log_kappa);
  //  PARAMETER(log_omega);  
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type vc_m = pow(exp(log_sdvc_m),2);
  Type vc_pair = pow(exp(log_sdvc_pair),2);
  // Type omega = exp(log_omega);  // Parameter of Weibull distribution
  //------------------------------------------

  // Spatial interpolation
  matrix<Type> AA = Zp*A; //Added here
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

  vector<Type> upl = Zp*up;
  vector<Type> uml = Zm*um;
  vector<Type> eta = X*beta + delta + upl + uml;

  for( int j=0; j<npair_mz + npair_dz*2; j++){
    nll -= dnorm( um(j), Type(0.0), exp(log_sdvc_m), true );
  }

  for( int j=0; j<npair_dz + npair_mz; j++){
    nll -= dnorm( up(j), Type(0.0), exp(log_sdvc_pair), true );
  }

  for( int j=0; j<npair_dz*2 + npair_mz*2; j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

//  for(int i=0; i<time.size(); i++){    
//    Type lambda = exp(eta(i));
//    Type t_omega = pow(time(i),omega);
//    Type S = exp(-lambda*t_omega);        // Survival function
//    Type f = lambda*omega*t_omega/time(i)*S;// Weibull density
//    if(notcens(i)){
//      nll -= log(f);
//    }else{
//      nll -= log(S); //The pasient survived until cencoring
//    }
//  }
  //---------------------------------------------
  
  //Report what we want to report----------------
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011);
  Type twin_vc_a = 2*vc_m;
  Type twin_vc_c = vc_pair - vc_m ;
  Type vc_res = pow(exp(log_sdvc_res),2);
  Type f_n_mz = Type(npair_mz);
  Type f_n_dz = Type(npair_dz);
  Type p_dz = Type(2.0)*f_n_dz / (f_n_dz + f_n_mz);    // if dz => os and mz => ss, calculate estimated vcs
  Type p_mz_given_ss = (1 - p_dz) / (f_n_mz / (f_n_mz + f_n_dz));
  Type sex_vc_a = 2/p_mz_given_ss*vc_m;
  //Type sex_vc_c = ((p_mz_given_ss + 1)*vc_pair - vc_m) / p_mz_given_ss ;
  Type sex_vc_c = vc_pair - vc_m / p_mz_given_ss ;
  Type vc_a=0.0 ;
  Type vc_c=0.0 ;
 
   if(is_twin_model == 1){
    vc_a=twin_vc_a;
    vc_c=twin_vc_c;
  } else {
    vc_a=sex_vc_a;
    vc_c=sex_vc_c;
  }
  ADREPORT(range);
  ADREPORT(vc_a);
  ADREPORT(vc_c);
  ADREPORT(vc_res);
  //ADREPORT(f_n_mz);
  //ADREPORT(f_n_dz);
  //ADREPORT(p_mz_given_ss);
  //ADREPORT(p_dz);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
