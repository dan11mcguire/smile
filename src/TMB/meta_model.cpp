// Generated by TMB: do not edit by hand

#define TMB_LIB_INIT R_init_meta_model
#include <TMB.hpp>
#include "fam_chol_matern_gmrf_sgps.hpp"
#include "fam_chol_matern_gmrf_s.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "fam_chol_matern_gmrf_sgps") {
    return fam_chol_matern_gmrf_sgps(this);
  } else if(model == "fam_chol_matern_gmrf_s") {
    return fam_chol_matern_gmrf_s(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
