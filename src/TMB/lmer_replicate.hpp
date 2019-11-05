#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type lmer_replicate(objective_function<Type>* obj)
{
  // Data
  DATA_INTEGER( npair_dz );
  DATA_INTEGER( npair_mz );
  DATA_FACTOR( pair_id );
  DATA_FACTOR( m_id );                
  DATA_VECTOR( y );                
//  DATA_SCALAR( sdvc_pair ); 
  // Parameters
  PARAMETER( mu );
  PARAMETER( log_sdvc_pair );
  PARAMETER( log_sdvc_m );
  PARAMETER( log_sdvc_res );
  PARAMETER_VECTOR( up );
  PARAMETER_VECTOR( um );
  
  Type jnll = 0;
  
  for( int i=0; i<npair_dz*2 + npair_mz*2; i++){
    jnll -= dnorm( y(i), mu + up(pair_id(i)) + um(m_id(i)), exp(log_sdvc_res), true );
  } 

  for( int j=0; j<npair_mz + npair_dz*2; j++){
    jnll -= dnorm( um(j), Type(0.0), exp(log_sdvc_m), true );
  }

  for( int j=0; j<npair_dz + npair_mz; j++){
    jnll -= dnorm( up(j), Type(0.0), exp(log_sdvc_pair), true );
  }
  
  // Reporting
//Type vc_res = pow(exp(log_sdvc_res),2);
//Type vc_m = pow(exp(log_sdvc_m),2);
//ADREPORT( vc_res );
//REPORT( vc_res );
//ADREPORT( vc_m );
//REPORT( vc_m );
ADREPORT( um );
REPORT( um );
ADREPORT( up );
REPORT( up ); 
ADREPORT( mu );
REPORT( mu );
//  
//  // bjas-correction testing
//  Type MeanZ = Z.sum() / Z.sjze();
//  Type SampleVarZ = ( (Z-MeanZ) * (Z-MeanZ) ).sum();
//  Type SampleSDZ = pow( SampleVarZ + 1e-20, 0.5);
//  REPORT( SampleVarZ );
//  REPORT( SampleSDZ );  
//  ADREPORT( SampleVarZ );
//  ADREPORT( SampleSDZ );  
//  
//
  return jnll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
