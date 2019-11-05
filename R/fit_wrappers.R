fit_by_chol<-function(ss,ycol="y"){
   
  ### ss$d must be sorted by "ttype", "pair", "sib" 
  y=ss$d[,ycol,with=F][[1]]
  X<-matrix(rep(1, nrow(ss$d)))
  spde<-inla.spde2.matern(ss$sp$mesh, alpha=2)
  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
  A<-ss$sp$A
  
  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
  vcd<-ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps", names(ss$d),value=T)]
  vcs<-sapply(unlist(vcd), function(x) x/sum(vcd))
  npair_dz<-ss$d[ttype=="dz",uniqueN(pair)]
  npair_mz<-ss$d[ttype=="mz",uniqueN(pair)]
  
  print("making sparse design matrices.")
  uftype<-ss$d[,unique(ttype)]
  Za_ref<-list("mz"=matrix(1,nrow=2,ncol=1),"dz"=diag(2))
  Za<-c()
  for(ff in 1:length(uftype)){
    Za<-c(Za,rep(list(Za_ref[uftype[ff]][[1]]), 
                 ss$d[ttype==uftype[ff], uniqueN(pair)])
         )
  }
  Za<-.bdiag(Za)


  Zc<-.bdiag(rep(list(matrix(1, nrow=2,ncol=1)), npair_dz + npair_mz)) 
  
  Ga_ref<-list("mz"=matrix(1,nrow=1,ncol=1),"dz"=matrix(c(1, .5, .5, 1),ncol=2))
  Ga<-c()
  for(ff in 1:length(uftype)){
    Ga<-c(Ga,rep(list(Ga_ref[uftype[ff]][[1]]), 
                 ss$d[ttype==uftype[ff], uniqueN(pair)])
         )
  }
  Ga<-.bdiag(Ga)

  Gc<-.bdiag(rep(list(matrix(1, nrow=1,ncol=1)),npair_dz + npair_mz)) 

  #GG<-.bdiag(list(Ga,Gc))
  #Lt<-chol(GG)
  #all.equal(GG, as(Lt%*%t(Lt), "dgTMatrix"))
  #all.equal(GG, as(t(Lt)%*%(Lt), "dgTMatrix"))

  print("finished making sparse design matrices.")
  Lta<-Matrix::chol(Ga)
  Ltc<-Matrix::chol(Gc)
  print("finished cholesky decomposition") 
  #Lta%*%matrix(rnorm(ncol(Za)),ncol=1) 
  data<-list(
             y=y,
             X=X,
             Za=Za,
             Zc=Zc,
             Lta=as(t(Lta),"dgTMatrix"),
             Ltc=as(t(Ltc),"dgTMatrix"),
             spdeMatrices=spdeMatrices,
             A=A
             )
  parameters <- list(
                     log_sdvc_a = 0.01,
                     log_sdvc_c = 0.01,
                     log_sdvc_res    = 0,
                     ua            = rep(0.0, ncol(Za)),
                     uc            = rep(0.0, ncol(Zc)),
                     beta          = rep(0, ncol(X)),
                     log_tau       = 0,
                     log_kappa     = 0,
                     x             = rep(0.0, nrow(data$spdeMatrices$M0))
                    )
#  compile("twin_or_sex_chol_direct_matern_gmrf.cpp")
#  dyn.load(dynlib("twin_or_sex_chol_direct_matern_gmrf"))
  print("making ad function")
  obj_chol <- MakeADFun(data, parameters, random=c("x", "ua", "uc"), DLL="twin_or_sex_chol_direct_matern_gmrf")
  print("fitting model.")   
  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
  print("finished fitting model.")   
  if(!inherits(opt_chol, "try-error")){
    rep_chol <- try(sdreport(obj_chol))
  }
  if(inherits(rep_chol, "try-error")){
    rt<-data.table(par=c("range", "vc_a", "vc_c", "vc_res"), est=NA, sd=NA, how="chol")
  } else {
    rt<- summ(rep_chol)[,how:="chol"]
  } 
  return(list(rep=rep_chol,opt=opt_chol))
}

fit_by_pair<-function(ss,ycol="y"){

  y=ss$d[,ycol,with=F][[1]]
  X<-matrix(rep(1, nrow(ss$d)))
  spde<-inla.spde2.matern(ss$sp$mesh, alpha=2)
  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
  A<-ss$sp$A
  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
  vcd<-ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps", names(ss$d),value=T)]
  vcs<-sapply(unlist(vcd), function(x) x/sum(vcd))
  
  npair_dz<-ss$d[ttype=="dz",uniqueN(pair)]
  npair_mz<-ss$d[ttype=="mz",uniqueN(pair)]
  #Zp<-as(model.matrix(~factor(pair) - 1, data=ss$d), "dgTMatrix")
  #Zm<-as(model.matrix(~factor(m) - 1, data=ss$d), "dgTMatrix")
  irow<-1:ss$d[,.N]
  jpair<-rep(1:ss$d[,uniqueN(pair)],times=ss$d[,uniqueN(pid),by=pair][,V1])  
  jm<-rep(1:ss$d[,uniqueN(m)],times=ss$d[,uniqueN(pid),by=m][,V1])  
  Zp<-sparseMatrix(i=irow, j=jpair, x=1,giveCsparse=F)
  Zm<-sparseMatrix(i=irow, j=jm, x=1,giveCsparse=F)
  #all.equal(Zp, Zp2)
  #all.equal(Zm, Zm2)
   
  data<-list(
             is_twin_model=1, 
             npair_dz=npair_dz,
             npair_mz=npair_mz,
             y=y,
             X=X,
             Zp=Zp,
             Zm=Zm,
             spdeMatrices=spdeMatrices,
             A=A
             )
  parameters <- list(
                     log_sdvc_pair = 0,
                     log_sdvc_m    = 0,
                     log_sdvc_res  = 1,
                     up            = rep(0.0, npair_mz + npair_dz),
                     um            = rep(0.0, npair_mz + npair_dz*2),
                     beta          = rep(0, ncol(X)),
                     log_tau       = 0,
                     log_kappa     = 0,
                     x             = rep(0.0, nrow(data$spdeMatrices$M0))
                    )

  obj_twin <- MakeADFun(data, parameters, random=c("x", "um", "up"), DLL="twin_or_sex_matern_gmrf")
   
  opt_twin <- try(nlminb(obj_twin$par, obj_twin$fn, obj_twin$gr))
  if(!inherits(opt_twin, "try-error")){
    rep_twin <- try(sdreport(obj_twin))
  }
  return(list(rep=rep_twin,
              opt=opt_twin))
}

fit_by_sex<-function(ss,ycol="y"){

  y=ss$d[,ycol,with=F][[1]]
  X<-matrix(rep(1, nrow(ss$d)))
  spde<-inla.spde2.matern(ss$sp$mesh, alpha=2)
  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
  A<-ss$sp$A
  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
  vcd<-ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps", names(ss$d),value=T)]
  vcs<-sapply(unlist(vcd), function(x) x/sum(vcd))
  
  npair_os<-ss$d[ss==0,uniqueN(pair)]
  npair_ss<-ss$d[ss==1,uniqueN(pair)]
  irow<-1:ss$d[,.N]
  jpair<-rep(1:ss$d[,uniqueN(pair)],times=ss$d[,uniqueN(pid),by=pair][,V1])  
  jm<-rep(1:ss$d[,uniqueN(m_ss)],times=ss$d[,uniqueN(pid),by=m_ss][,V1])  
  Zp<-sparseMatrix(i=irow, j=jpair, x=1,giveCsparse=F)
  Zm<-sparseMatrix(i=irow, j=jm, x=1,giveCsparse=F)
  #Zp<-model.matrix(~factor(pair) - 1, data=ss$d)
  #Zm<-model.matrix(~factor(m_ss) - 1,data=ss$d)
  
  data<-list(
             is_twin_model=0, 
             npair_dz=npair_os,
             npair_mz=npair_ss,
             y=y,
             X=X,
             Zp=Zp,
             Zm=Zm,
             spdeMatrices=spdeMatrices,
             A=A
             )
  parameters <- list(
                     log_sdvc_pair = 0,
                     log_sdvc_m    = 0,
                     log_sdvc_res  = 1,
                     up            = rep(0.0, npair_ss + npair_os),
                     um            = rep(0.0, npair_ss + npair_os*2),
                     beta          = rep(0, ncol(X)),
                     log_tau       = 0,
                     log_kappa     = 0,
                     x             = rep(0.0, nrow(data$spdeMatrices$M0))
                    )
  data$Zp<-as(data$Zp, "dgTMatrix")   
  data$Zm<-as(data$Zm, "dgTMatrix")  
 
  obj_sex <- MakeADFun(data, parameters, random=c("x", "um", "up"), DLL="twin_or_sex_matern_gmrf")
  opt_sex <-try(nlminb(obj_sex$par, obj_sex$fn, obj_sex$gr))
  if(!inherits(opt_sex, "try-error")){
    rep_sex <- try(sdreport(obj_sex)) 
  }
  if(inherits(rep_sex, "try-error")){
    rs<-data.table(par=c("range", "vc_a", "vc_c", "vc_res"), est=NA, sd=NA, how="sex")
  } else {
    rs<- summ(rep_sex)[,how:="sex"]
  } 
  return(list(opt=opt_sex,rep=rep_sex))
}
