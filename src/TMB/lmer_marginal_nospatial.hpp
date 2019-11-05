#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
	
template<class Type>
Type objective_function<Type>::operator() ()
{  
  using namespace density;
  // Data
  DATA_VECTOR( y );               
  DATA_SPARSE_MATRIX(Zp);         //Design matrix for sy-pair random effects 
  DATA_SPARSE_MATRIX(Zm);         //Design matrix for sy-m random effects 

  // Parameters
  PARAMETER( mu );
  PARAMETER_VECTOR(log_sd);
  
  
  int i,j;
  int n=y.size();
  
  //Transformed parameters 
  vector<Type> sd = exp(log_sd);  
  vector<Type> muvec(y.size()) ;  
  for(i=0; i < n; i++)
  {
    muvec(i) = mu;
  }

  matrix<Type> Gp(Zp.cols(), Zp.cols());
  for(i=0; i < Zp.cols(); i++)
  { 
    Gp(i,i)=pow(sd[0], 2.0); 
    for(j=0; j < i; j++){
      Gp(i,j)=Type(0.0);
      Gp(j,i)=Type(0.0);
    } 
  } 

  matrix<Type> Gm(Zm.cols(),Zm.cols());
  for(i=0; i < Zm.cols(); i++)
  { 
    Gm(i,i)=pow(sd[1],2.0); 
    for(j=0; j < i; j++){
      Gm(i,j)=Type(0.0);
      Gm(j,i)=Type(0.0);
    } 
  } 
  matrix<Type> R(n,n);
  for(i=0; i < n; i++)
  {
    //R(i,i)=Type(1.0);
    R(i,i)=pow(sd[2],2.0);
    for(j=0; j < i; j++)
    {
      R(i,j)=Type(0.0);
      R(j,i)=Type(0.0);
    }
  }
  
  Type mnll = 0;

  //MVNORM_t<Type> mvnorm(pow(sd[2],2.0)*R + Zp*Gp*Zp.transpose() + Zm*Gm*Zm.transpose());
  MVNORM_t<Type> mvnorm(R + Zp*Gp*Zp.transpose() + Zm*Gm*Zm.transpose());
  vector<Type> res = y - muvec ; 
  
  mnll=mvnorm(res) ;
  

  Type vc_a = 2*pow(exp(sd[1]),2);
  Type vc_c = pow(exp(sd[0]),2) - pow(exp(sd[1]),2) ;
  Type vc_res = pow(exp(sd[2]),2);
  ADREPORT(vc_a);
  ADREPORT(vc_c);
  ADREPORT(vc_res);
  
 
  return mnll;
}

 
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
