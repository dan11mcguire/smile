#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

//namespace car_gmrf {
//using namespace Eigen;
//using namespace tmbutils;
//
//
//template<class Type>
//struct car_gmrf_t{  
//  SparseMatrix<Type> diagI;        // G0 eqn (10) in Lindgren 
//  SparseMatrix<Type> N;        // G1 eqn (10) in Lindgren 
//  SparseMatrix<Type> Mpm5;        // G2 eqn (10) in Lindgren 
//  car_gmrf_t(SEXP x){  /* x = List passed from R */
//  diagI = asSparseMatrix<Type>(getListElement(x,"diagI"));
//  N = asSparseMatrix<Type>(getListElement(x,"N"));
//  Mpm5 = asSparseMatrix<Type>(getListElement(x,"Mpm5"));
//}
//};
//
//template<class Type> 
//  SparseMatrix<Type> Q_car(car_gmrf_t<Type> car_mats,Type rho, Type s2_s){
//	  //See Proposition 1 J.M.VerHoefetal./SpatialStatistics25(2018)68–8
//    Type tau2=Type(1.0) / s2_s;
//  return tau2*(car_mats.Mpm5*(car_mats.diagI - rho*car_mats.N)*car_mats.Mpm5); 
//  }
//};


template<class Type>
  Type fam_chol_iid_gmrf_s_us(objective_function<Type>* obj)
{
  using namespace car_gmrf; //
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR( y );                
  //Load data and parameters----------------
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_SCALAR(rhomax);
  DATA_SPARSE_MATRIX(Zs);         //Design matrix for spatial random effects 
  DATA_STRUCT(car_mats,car_gmrf_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated

  // Parameters
  //PARAMETER( mu );
  PARAMETER( log_sdvc_s );
  PARAMETER( log_sdvc_res );
  
  PARAMETER_VECTOR(beta);  
  PARAMETER(logit_rho);
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  Type rho = exp(logit_rho)/(1+exp(logit_rho))*rhomax;
  Type vc_s = pow(exp(log_sdvc_s),2);
  Type vc_res = pow(exp(log_sdvc_res),2);
  //------------------------------------------

  // Spatial interpolation
  vector<Type> delta = Zs*x;
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q = Q_iid(car_mats,rho,vc_s);
  //---------------------------------------------
  
  //Calculates nll-------------------------------
  Type nll = 0.0;
  nll = GMRF(Q)(x);                              

  //nll = GMRF(Q, false)(x);  
  // Return un-normalized density on request
  //if (flag == 0) return nll;

  vector<Type> eta = X*beta + delta ;

  for( int j=0; j< X.rows(); j++){
    nll -= dnorm( y(j), eta(j), exp(log_sdvc_res), true );
  } 

  //---------------------------------------------
  
  //Report what we want to report----------------
  ADREPORT(vc_res);
  ADREPORT(vc_s);
  ADREPORT(rho);

  //---------------------------------------------
  
  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
