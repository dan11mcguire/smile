#n_dz=500
#	n_mz=500
#	p_ssdz=0.5
#	A=0.7 
#	C=0.1
#	mu=10

simTwinsSex<-
  function(n_dz=500,
           n_mz=500,
           p_ssdz=0.5,
           A=0.7, 
           C=0.1,
           mu=10){
  
  E=1-A-C
  s2Pair<-.5*A + C
  s2M<-.5*A
  s2MZ<-s2Pair + s2M
  s2DZ<-s2Pair
  u_pair<-rep(rnorm((n_dz + n_mz), 0, sqrt(s2Pair)), each=2)
  u_m<-c(rnorm(n_dz*2, 0, sqrt(s2M)), rep(rnorm(n_mz, 0, sqrt(s2M)), each=2))
  eps<-rnorm((n_dz + n_mz)*2, 0, sqrt(E))
  
  d<-data.table(u_pair, u_m, eps, 
                pair=rep(1:(n_dz + n_mz),each=2),
                m=c((1:(n_dz*2)), 
                    rep((n_dz*2+1):(n_dz*2 + n_mz),
                        each=2)),
                sib=rep(c(1,2), times=n_mz + n_dz),
                ttype=rep(c("dz", "mz"),times=c(n_dz*2, n_mz*2)),
                mu=mu)

  ssdz.pairs<-d[ttype=="dz",unique(pair)]
  ssdz.pairs<-ssdz.pairs[which(rbinom(n_dz, 1,0.5)==1)]
  d[ttype=="mz", ss:=1L]
  d[pair %in% ssdz.pairs,ss:=1L]
  d[is.na(ss),ss:=0L]
  d[,y:=u_pair + u_m + eps + mu]
  n_ss<-d[ss==1,length(unique(pair))]
  setkey(d,pair)
  d[ss==1,m_ss:=rep(1:n_ss, each=2)]
  d[ss==0,m_ss:=(n_ss+1):(n_ss + d[,sum(ss==0)])]
  d[,pid:=1:.N] 
}

#mllz<-function(pair_i,twin_type="mz",par){
# ui_pair<-par[1];
#if(twin_type=="mz"){
# uij_mzextrass<-par[2]
# ll<-d[pair==pair_i,
#      dmvnorm(y,
#      mean=rep(mu.hat + ui_pair + uij_mzextrass,2),
#      sigma=diag(resvar.hat,2),log=TRUE)]
#  llu<-dmvnorm(matrix(rep(c(ui_pair, uij_mzextrass),times=c(2,2)),ncol=2),
#              mean=c(0,0),
#              sigma=diag(c(pair_ss.hat, mz_extrass.hat)),log=TRUE)
#} else {
#  uij_mzextrass<-par[2:3]
#  ll<-d[pair==pair_i,
#      dmvnorm(y,
#      mean=c(mu.hat + ui_pair + uij_mzextrass),
#      sigma=diag(resvar.hat,2))]
#  llu<-dmvnorm(matrix(c(ui_pair, ui_pair, uij_mzextrass),ncol=2),
#	        mean=c(0,0),
#	        sigma=diag(c(pair_ss.hat, mz_extrass.hat)),log=TRUE)
#}
#  
# ll<-sum(c(ll,llu))
#     
# return(-1*ll)
#}
#


extract_vc<-function(model,name){
  vc<-as.data.table(VarCorr(model))[grp==name,vcov]
  return(vc)
}

