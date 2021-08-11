invlogit<-function(x){ exp(x) / (1 + exp(x))}
logit<-function(x){ log(x/(1-x))}

sim_twins_spatial<-function(n_dz = 5, 
                            n_mz = 5, 
                            X=NULL,loc=NULL,
                            p_ssdz = 0.5, 
                            A = 0.2, C = 0.1, S=0.2, 
                            mu = 10, ...){   
    E = 1 - A - C - S
    s2Pair <- 0.5 * A + C
    s2M <- 0.5 * A
    s2MZ <- s2Pair + s2M
    s2DZ <- s2Pair
    d<-gen_twin_ids(n_dz, n_mz) 
    d[,u_pair:=rnorm(1, 0, sqrt(s2Pair)),by=pair]
    d[,u_m:=rnorm(1, 0, sqrt(s2M)),by=m]
    d[,eps:=rnorm(1, 0, sqrt(E)),by=pid] 
    if(!missing(X)){
      beta<-runif(ncol(X))
      d[,fixed:=as.vector(X%*%beta)]
      d<-cbind(d, X)
    } else {
      d[,fixed:=mu]
    } 

    loc<-sim_spatial(n=d[,uniqueN(pair)],...) 
    locdt<-data.table(loc$loc)[,pair:=1:.N]
    locdt[,u_s:=as.vector(loc$A%*%loc$u)]
    d<-merge(d, locdt, by="pair",sort=F)
    d[,y:=fixed + u_pair + u_m + u_s + eps]
   return(list(
      d=d,
      loc=loc,
      s2M=s2M,
      s2Pair=s2Pair,
      s2res=E
   )) 
}

gen_twin_ids<-function(n_dz = 500, n_mz = 500){

    id_dt<-data.table(pair = rep(1:(n_dz + n_mz), each = 2),    
                      m = c((1:(n_dz * 2)), rep((n_dz * 2 + 1):(n_dz * 2 + n_mz), each = 2)), 
                      sib = rep(c(1, 2), times = n_mz + n_dz), 
                      ttype = rep(c("dz", "mz"), times = c(n_dz * 2, n_mz * 2)))

    ssdz.pairs <- id_dt[ttype == "dz", unique(pair)]
    ssdz.pairs <- ssdz.pairs[which(rbinom(n_dz, 1, 0.5) == 1)]
    id_dt[ttype == "mz", `:=`(ss, 1L)]
    id_dt[pair %in% ssdz.pairs, `:=`(ss, 1L)]
    id_dt[is.na(ss), `:=`(ss, 0L)]
    n_ss <- id_dt[ss == 1, length(unique(pair))]
    setkey(id_dt, pair)
    id_dt[ss == 1, `:=`(m_ss, rep(1:n_ss, each = 2))]
    id_dt[ss == 0, `:=`(m_ss, (n_ss + 1):(n_ss + id_dt[, sum(ss == 0)]))]
    id_dt[, `:=`(pid, 1:.N)]

    return(id_dt)
}


sim_spatial<-function(n,inla.seed ,use_my_sample=T, sigma.u=1, rng=2,bmax=NULL){

  if(missing(bmax)) bmax<-10 
  fake.locations = matrix(c(0,0,bmax,bmax, 0, bmax, bmax, 0), 
			  nrow = 4, byrow = T)
  mesh = inla.mesh.2d(loc = fake.locations, max.edge=c(0.5, 1))

  #if(mesh$n > n){
  #  mesh = inla.mesh.2d(loc = fake.locations, max.edge=c(0.5, 1),max.n=n/2)
  #}
  if(sigma.u > 0){
    spde = inla.spde2.pcmatern(mesh, prior.range=c(rng, 0.5), prior.sigma=c(sigma.u,0.5))
    spdeMatrices = spde$param.inla[c("M0","M1","M2")]
    
    Qu = inla.spde.precision(spde, theta=c(log(rng), log(sigma.u)))
    
    set.seed(inla.seed)
    if(use_my_sample){
    #  cholQu<-Matrix::chol(Qu) 
    #  cholQuinv<-t(chol(chol2inv(cholQu)))  ## triple check this line later (t)?
    #  u<-cholQuinv%*%matrix(rnorm(dim(Qu)[1],sd=sigma.u),ncol=1)
      cholQu<-t(Matrix::chol(Qu)) 
     # all.equal(cholQu%*%t(cholQu), Qu)
      cholQuinv<-chol(chol2inv(cholQu))
     # all.equal(cholQuinv%*%t(cholQuinv), solve(Qu))
      u<-cholQuinv%*%matrix(rnorm(dim(Qu)[1]),ncol=1)
    } else {
      u = inla.qsample(n=1, Q=Qu, seed = inla.seed)
    }
      # cma <- chol(ma  <- cbind(1, 1:3, c(1,3,7)))
      # ma %*% chol2inv(cma)
    
    u = u[ ,1]
    
    set.seed(inla.seed)
  } else {
    u<-rep(0, mesh$n) 
    Qu<-diag(mesh$n)
  }
  loc = matrix(runif(2*n), n)*bmax # coordinates
  colnames(loc)<-c("lat", "lon")
  A = inla.spde.make.A(mesh,loc)
  #plot(mesh)
  #points(loc)
  #var(as.vector(A%*%u))
  return(list(
    loc=loc,
    Qu=Qu,
    mesh=mesh, 
    A=A,
    u=u
  ))
}



#corr<-function(x, y){
#  n<-length(x)
# (sum(x*y) - n*mean(x)*mean(y)) / ((n-1)*sd(x)*sd(y))
#}
#mse<-function(x,y){
#  mean((x-y)^2)
#}
#for(i in 1:1000){
# o[i,]<-c(mse(x)
#
#}

sim_fam_spatial<-function(  n_trio = 5,
                            n_quad = 5,
			    create_ustar=TRUE,
                            inla.seed=NULL,
			    mesh=NULL,
			    ustar=NULL,
                            X=NULL,
                            loc=NULL,
                            Ag = 0.2,
                            C_par = 0.1,
                            C_sib = 0.1,
                            C_fam = 0.1,
                            S=0.2,
                            mu = 0,
			    sigma.u=1,
			    rng=.1,
			    bmax=NULL,
			    rescale_to_S=FALSE,
			    use_my_sample=TRUE,
                            ...){

     
    if(missing(inla.seed)) inla.seed<-sample.int(10^3,1)

    set.seed(inla.seed) 

    E = 1 - Ag - C_par - C_sib - C_fam - S
    id_dt<-gen_fam_ids(n_trio, n_quad) 
    id_dt[emprel %in% c(1:2),empreltype:="parent"]
    id_dt[emprel %in% 3,empreltype:="sib"]
    mlists<-make_matrix_lists(id_dt) 
     
    u_a<-as.vector(mlists$Ltlist$Lt_Ga %*%matrix(rnorm(ncol(mlists$Ltlist$Lt_Ga), 0, sqrt(Ag)),ncol=1))
    u_sib<-as.vector(mlists$Zlist$Zc_sib %*% matrix(rnorm(ncol(mlists$Zlist$Zc_sib), 0, sqrt(C_sib)),ncol=1))
    u_par<-as.vector(mlists$Zlist$Zc_par %*% matrix(rnorm(ncol(mlists$Zlist$Zc_par), 0, sqrt(C_par)),ncol=1))
    u_fam<-as.vector(mlists$Zlist$Zc_fam %*% matrix(rnorm(ncol(mlists$Zlist$Zc_fam), 0, sqrt(C_fam)),ncol=1))
    
    id_dt[,u_a:=u_a]
    id_dt[,u_sib:=u_sib]
    id_dt[,u_par:=u_par]
    id_dt[,u_fam:=u_fam]
    id_dt[,eps:=rnorm(1, 0, sqrt(E)),by=pid]

    if(!missing(X)){
      beta<-runif(ncol(X))
      id_dt[,fixed:=as.vector(X%*%beta)]
      id_dt<-cbind(d, X)
    } else {
      id_dt[,fixed:=mu]
    } 
   # if(all(!missing(mesh, ustar))) 
    if(create_ustar){
      loc <- sim_spatial(n = id_dt[, uniqueN(set)], inla.seed = inla.seed, 
        sigma.u = sigma.u, rng = rng, use_my_sample = use_my_sample) 
      locdt <- data.table(loc$loc)[, `:=`(set, 1:.N)]
      locdt[, `:=`(u_s, as.vector(loc$A %*% loc$u))] 

      sp<-list(loc=loc$loc,
              mesh=loc$mesh,
              A=loc$A, 
              u=loc$u,
              Qu=loc$Qu)

    } else {
      if(missing(mesh) | missing(ustar)) stop("need to provide the 'mesh' and generated random effect 'ustar' to project at new locations")
      if(missing(bmax)) bmax<-10
      n<-id_dt[,uniqueN(set)] 
      loc <- matrix(runif(2 * n), n) * bmax
      colnames(loc)<-c("lat", "lon")
      A <- inla.spde.make.A(mesh, loc)
      locdt<-data.table(loc)[,set:=1:.N]
      locdt[,u_s:=as.vector(A %*% as.matrix(ustar,ncol=1))]

      sp<-list(loc=loc,
               mesh=mesh,
               A=A, 
               u=ustar)
    }
    
    if(rescale_to_S){
      vus<-locdt[,var(u_s)*(.N-1)/.N]
      scaleus<-sqrt(S/vus)
      locdt[,u_s:=scaleus*u_s]
    }
    id_dt<-merge(id_dt, locdt, by="set",sort=F)
    id_dt[,y:=fixed + u_a + u_sib + u_par + u_fam + u_s + eps]
    
    return(list(
      d=id_dt,
      sp=sp,
      Zlist=mlists$Zlist,
      Glist=mlists$Glist,
      Ltlist=mlists$Ltlist,
      inla.seed=inla.seed,
      Ag=Ag,
      C_par=C_par,
      C_fam=C_fam,
      C_sib=C_sib,
      S=S,
      E=E
   )) 
}

    

make_spatial_mats<-function(family_data, sigma.u=1, rng=2,bmax=NULL,mesh=NULL){
  
  if(missing(mesh)){
    if(missing(bmax)) bmax<-10 
    fake.locations = matrix(c(0,0,bmax,bmax, 0, bmax, bmax, 0), 
          		  nrow = 4, byrow = T)
    mesh = inla.mesh.2d(loc = fake.locations, max.edge=c(0.5, 1))
    #if(mesh$n > n){
    #  mesh = inla.mesh.2d(loc = fake.locations, max.edge=c(0.5, 1),max.n=n/2)
    #}
  } 
  spde = inla.spde2.pcmatern(mesh, prior.range=c(rng, 0.5), prior.sigma=c(sigma.u,0.5))
  spdeMatrices = spde$param.inla[c("M0","M1","M2")]
  Qu = inla.spde.precision(spde, theta=c(log(rng), log(sigma.u)))
  loc = as.matrix(family_data[,.(lat,lon)][,unique(.SD)])
  colnames(loc) <- c("lat", "lon")
  A = inla.spde.make.A(mesh,loc)
  return(list(
    loc=loc,
    Qu=Qu,
    mesh=mesh, 
    A=A
  ))
}

#sim_fam_spatial<-function(  n_trio = 5,
#                            n_quad = 5,
#                            inla.seed=NULL,
#                            X=NULL,
#                            loc=NULL,
#                            Ag = 0.2,
#                            C_par = 0.1,
#                            C_sib = 0.1,
#                            C_fam = 0.1,
#                            S=0.2,
#                            mu = 0,
#			    sigma.u=1,
#			    rng=.1,
#			    rescale_to_S=FALSE,
#			    use_my_sample=T,
#                            ...){
#
#    Za_ref<-list("mz"=matrix(1,nrow=2,ncol=1),
#                 "dz"=diag(2), "",
#                 "trio" = diag(3),
#                 "quad" = diag(4)) 
#  
#  
#  
#  
#    Ga_ref<-list("mz"=matrix(1,nrow=1,ncol=1),
#                 "dz"=matrix(c(1, .5, .5, 1),ncol=2),
#                 "trio"=matrix(c(1, .5, .5, 
#                                 .5, 1, .5, 
#                                 .5, .5, 1),ncol=3),
#                 "quad"=matrix(c(1, 0, .5, .5, 
#                                 0, 1, .5, .5, 
#                                 .5, .5, 1, .5, 
#                                 .5, .5, .5, 1),ncol=4, byrow=T))
#  
#  
#    Zc_sib_ref<-list("mz"=matrix(1,nrow=2,ncol=1),
#                     "dz"=matrix(1,nrow=2,ncol=1),
#                     "trio"=matrix(c(1, 0, 
#                                     0, 1, 
#                                     0, 1),nrow=3,byrow=T),
#                     "quad"=matrix(c(1, 0, 0, 
#                                     0, 1, 0, 
#                                     0, 0, 1, 
#                                     0, 0, 1),nrow=4,byrow=T)
#                     )
#  
#    Gc_sib_ref<-list("mz"=matrix(1),
#                     "dz"=matrix(1),
#                     "trio"=diag(2),
#                     "quad"=diag(3))
#  
#    Zc_par_ref<-list("mz"=diag(2),
#                     "dz"=diag(2),
#                     "trio"=diag(3),
#                     "quad"=matrix(c(1, 0, 0, 
#                                     1, 0, 0, 
#                                     0, 1, 0, 
#                                     0, 0, 1),nrow=4,byrow=T)
#                     )
#  
#    Gc_par_ref<-list("mz"=diag(2),
#                     "dz"=diag(2),
#                     "trio"=diag(3),
#                     "quad"=diag(3))
#  
#    Zc_fam_ref<-list("mz"=matrix(1, nrow=2),
#                     "dz"=matrix(1, nrow=2),
#                     "trio"=matrix(1, nrow=3),
#                     "quad"=matrix(1, nrow=4)
#                     )
#    
#  
#  
#    Gc_fam_ref<-list("mz"=diag(1),
#                     "dz"=diag(1),
#                     "trio"=diag(1),
#                     "quad"=diag(1))
#      
#     
#     
#    if(missing(inla.seed)) inla.seed<-sample.int(10^3,1)
#
#    set.seed(inla.seed) 
#
#    E = 1 - Ag - C_par - C_sib - C_fam - S
#    id_dt<-gen_fam_ids(n_trio, n_quad) 
#    id_dt[emprel %in% c(1:2),empreltype:="parent"]
#    id_dt[emprel %in% 3,empreltype:="sib"]
#     
#    uftype<-id_dt[,unique(ttype)] 
#    Za<-Ga<-c()
#    Zc_sib<-Gc_sib<-c()
#    Zc_par<-Gc_par<-c()
#    Zc_fam<-Gc_fam<-c()
#
#
#    for(ff in 1:length(uftype)){
#
#      Za<-c(Za,rep(list(Za_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#      Zc_sib<-c(Zc_sib,rep(list(Zc_sib_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#      Zc_par<-c(Zc_par,rep(list(Zc_par_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#      Zc_fam<-c(Zc_fam,rep(list(Zc_fam_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#     
#      Ga<-c(Ga,rep(list(Ga_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#      Gc_sib<-c(Gc_sib,rep(list(Gc_sib_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#      Gc_par<-c(Gc_par,rep(list(Gc_par_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#      Gc_fam<-c(Gc_fam,rep(list(Gc_fam_ref[uftype[ff]][[1]]),id_dt[ttype==uftype[ff], uniqueN(set)]))
#
#    }
#    
#    Za<-.bdiag(Za)
#    Zc_sib<-.bdiag(Zc_sib)
#    Zc_par<-.bdiag(Zc_par) 
#    Zc_fam<-.bdiag(Zc_fam)
#
#    Ga<-.bdiag(Ga)
#    Gc_sib<-.bdiag(Gc_sib)
#    Gc_par<-.bdiag(Gc_par) 
#    Gc_fam<-.bdiag(Gc_fam)
#
#    Lt_Ga<-t(Matrix::chol(Ga))
#    Lt_Gc_sib<-t(Matrix::chol(Gc_sib))
#    Lt_Gc_par<-t(Matrix::chol(Gc_par))
#    Lt_Gc_fam<-t(Matrix::chol(Gc_fam))
#
#    u_a<-as.vector(Za %*% Lt_Ga %*%matrix(rnorm(ncol(Za), 0, sqrt(Ag)),ncol=1))
#    u_sib<-as.vector(Zc_sib %*% Lt_Gc_sib%*%matrix(rnorm(ncol(Zc_sib), 0, sqrt(C_sib)),ncol=1))
#    u_par<-as.vector(Zc_par %*% Lt_Gc_par%*%matrix(rnorm(ncol(Zc_par), 0, sqrt(C_par)),ncol=1))
#    u_fam<-as.vector(Zc_fam %*% Lt_Gc_fam%*%matrix(rnorm(ncol(Zc_fam), 0, sqrt(C_fam)),ncol=1))
#
#    id_dt[,u_a:=u_a]
#    id_dt[,u_sib:=u_sib]
#    id_dt[,u_par:=u_par]
#    id_dt[,u_fam:=u_fam]
#    id_dt[,eps:=rnorm(1, 0, sqrt(E)),by=pid]
#
#    if(!missing(X)){
#      beta<-runif(ncol(X))
#      id_dt[,fixed:=as.vector(X%*%beta)]
#      id_dt<-cbind(d, X)
#    } else {
#      id_dt[,fixed:=mu]
#    } 
#
#    #loc<-sim_spatial(n=id_dt[,uniqueN(set)],inla.seed=seed+1L,...) 
#    loc<-sim_spatial(n=id_dt[,uniqueN(set)],
#		     inla.seed=inla.seed,
#		     sigma.u=sigma.u,
#		     rng=rng,
#		     use_my_sample=use_my_sample) 
#
#    locdt<-data.table(loc$loc)[,set:=1:.N]  
#    locdt[,u_s:=as.vector(loc$A%*%loc$u)]
#    
#    if(rescale_to_S){
#      vus<-locdt[,var(u_s)*(.N-1)/.N]
#      scaleus<-sqrt(S/vus)
#      locdt[,u_s:=scaleus*u_s]
#    }
#    
#    id_dt<-merge(id_dt, locdt, by="set",sort=F)
#    id_dt[,y:=fixed + u_a + u_sib + u_par + u_fam + u_s + eps]
#   
#  return(list(
#      d=id_dt,
#      sp=list(loc=loc$loc,
#                  mesh=loc$mesh,
#                  A=loc$A, 
#                  u=loc$u,
#                  Qu=loc$Qu),
#      Zlist=list(Za=Za,Zc_sib=Zc_sib,Zc_par=Zc_par,Zc_fam=Zc_fam),
#      Glist=list(Ga=Ga,Gc_sib=Gc_sib,Gc_par=Gc_par,Gc_fam=Gc_fam),
#      Ltlist=list(Lt_Ga=Lt_Ga,Lt_Gc_sib=Lt_Gc_sib,Lt_Gc_par=Lt_Gc_par,Lt_Gc_fam=Lt_Gc_fam),
#      inla.seed=inla.seed,
#      Ag=Ag,
#      C_par=C_par,
#      C_fam=C_fam,
#      C_sib=C_sib,
#      S=S,
#      E=E
#   )) 
#}





plot_fitted_spatial<-function(sim_twin_ob=ss, tmb_rep=rep){
  idx<-names(rep$par.random)=="x"
  idlogtau<-which(names(rep$par.fixed)=="log_tau")
  us_est<-rep$par.random[idx]/exp(rep$par.fixed[idlogtau])
  quilt.plot(x=ss$loc$loc[,1],
             y=ss$loc$loc[,2],
             z=as.vector(ss$loc$A%*%us_est),nx=80,ny=80, 
             col = plasma(101), 
             main="Field projected to data locations", 
             zlim = range(us_est))

}

plot_true_spatial<-function(sim_twin_ob=ss){

  quilt.plot(x=ss$loc$loc[, 1],y=ss$loc$loc[, 2],z=as.vector(ss$loc$A%*%ss$loc$u),nx=80,ny=80, 
             col = plasma(101), main="True spatial effects", 
             zlim = range(ss$loc$u))

}



