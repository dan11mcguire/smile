to_liability<-function(vc,K){
  t<-qnorm(1-K)
  z<-dnorm(t)
  i<-z/K
  Eb<-K + vc/K
  tvc<-qnorm(1-Eb)
  rvc <- (t - tvc*sqrt(1-(t^2 - tvc^2)*(1-t/i))) / (i + tvc^2*(i-t))
  return(rvc)
}


to_observed<-function(r,K){
  t<-qnorm(1-K)
  z<-dnorm(t)
  i<-z/K
  tvc<-(t - i*r) / sqrt(1 - r^2*i*(i - t))
  Eb<-1-pnorm(tvc) 
  vc<-K*(Eb - K)
  return(vc)
 }

