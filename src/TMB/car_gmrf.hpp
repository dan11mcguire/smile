

namespace car_gmrf {
using namespace Eigen;
using namespace tmbutils;


template<class Type>
struct car_gmrf_t{  
  SparseMatrix<Type> diagI;        // G0 eqn (10) in Lindgren 
  SparseMatrix<Type> N;        // G1 eqn (10) in Lindgren 
  SparseMatrix<Type> Mpm5;        // G2 eqn (10) in Lindgren 
  car_gmrf_t(SEXP x){  /* x = List passed from R */
  diagI = asSparseMatrix<Type>(getListElement(x,"diagI"));
  N = asSparseMatrix<Type>(getListElement(x,"N"));
  Mpm5 = asSparseMatrix<Type>(getListElement(x,"Mpm5"));
}
};

template<class Type> 
  SparseMatrix<Type> Q_car(car_gmrf_t<Type> car_mats,Type rho, Type s2_s){
	  //See Proposition 1 J.M.VerHoefetal./SpatialStatistics25(2018)68–8
    Type tau2=Type(1.0) / s2_s;
  return tau2*(car_mats.Mpm5*(car_mats.diagI - rho*car_mats.N)*car_mats.Mpm5); 
  }
};
