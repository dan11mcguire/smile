summ<-function(report){
  dt<-data.table(par=names(report$value),
             est=report$value, 
  	   sd=sqrt(diag(report$cov)))
  return(dt)
}

###################################################################
##  fit variance component models with TMB and nlminb, wrappers ###
###################################################################
###################################################################
##                                                               ##
## s = spatial  (first s)                                        ##
## g = additive genetics                                         ##
## f = within-family shared random effect                        ##
## p = parent-within-family shared random effect                 ##
## s = sibling within family shared random effect (last s)       ##
##                                                               ##
###################################################################
###################################################################

fit_fam<-function(ss, vc_struct, ycol="y",family="gaussian", silent=TRUE){

  vc_str<-vc_f(vc_struct)
  link<-link_f(family) 
  tmplt<-paste0(vc_str, link)
  
  dl<-data_f(ss, tmplt, ycol)
  strt.fixed<-coef(glm(dl$y~dl$X-1,family=gsub("over","",family)))
  params<-param_f(ss, dl, tmplt,strt.fixed)
  re_params<-grep("^x|^u", names(params),value=T) 
  if(grepl("over",family)) re_params<-c("res", re_params) 
  
  obj_chol <- MakeADFun(data=c(model=tmplt, dl), 
                        parameters=params, 
                        random=re_params, 
                        silent=silent,
                        DLL="meta_model")
   
  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
  if(!inherits(opt_chol, "try-error")){
    rep_chol <- try(sdreport(obj_chol))
  }
  if(inherits(rep_chol, "try-error")){
    rt<-data.table(par=NA, est=NA, sd=NA, how="chol")
  } else {
    rt<- summ(rep_chol)[,how:=vc_struct]
  } 
  return(list(rep=rep_chol,opt=opt_chol,rt=rt))
 # return(list(re_params=re_params,dl=dl, tmplt=tmplt, params=params))
}

link_f<-function(link){
  switch(link,
            gaussian = "", 
            binomial  = "_logit", 
            overbinomial  = "_overlogit", 
            poisson   = "_pois",
            overpoisson = "_overpois")
}

vc_f<-function(vc_struct){
  if(!vc_struct %in% c("sgfps",
                       "sgfp",
                       "sgfs",
                       "sgps",
                       "gfps",
                       "sgf",
                       "sgp",
                       "sgs",
                       "gfp",
                       "gfs",
                       "gps",
                       "sg",
                       "gf",
                       "gp",
                       "gs",
                       "g",
		       "s")) stop(paste0("template not found. did you mean to specify '", vc_struct, "'?"))
#  switch(vc_struct,
#            sgfps = "fam_chol_matern_gmrf_sgfps", 
#            gfps  = "fam_chol_matern_gmrf_gfps", 
#            sgps  = "fam_chol_matern_gmrf_sgps", 
#            gps   = "fam_chol_matern_gmrf_gps")
  paste0("fam_chol_matern_gmrf_",vc_struct)
}




data_f<-function(ss,tmplt,ycol){
  
  y=ss$d[,ycol,with=F][[1]]
  X<-matrix(rep(1, nrow(ss$d)))
#  spde<-inla.spde2.pcmatern(ss$sp$mesh, prior.range=c(0.5,0.5),prior.sigma=c(0.5,0.5))
  spde<-inla.spde2.matern(ss$sp$mesh, alpha=2)
  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
  A<-ss$sp$A
  
  d<-switch(tmplt,
            fam_chol_matern_gmrf_sgfps =  list(
                          y=y,
                          X=X,
                          Za=ss$Zlist[["Za"]],
                          Zc_fam=ss$Zlist[["Zc_fam"]],
                          Zc_par=ss$Zlist[["Zc_par"]],
                          Zc_sib=ss$Zlist[["Zc_sib"]],
                          Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
                          Lt_Gc_fam=as(ss$Ltlist[["Lt_Gc_fam"]],"dgTMatrix"),
                          Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
                          Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
                          spdeMatrices=spdeMatrices,
                          A=A
                          ),

            fam_chol_matern_gmrf_gfps =  list(
                          y=y,
                          X=X,
                          Za=ss$Zlist[["Za"]],
                          Zc_fam=ss$Zlist[["Zc_fam"]],
                          Zc_par=ss$Zlist[["Zc_par"]],
                          Zc_sib=ss$Zlist[["Zc_sib"]],
                          Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
                          Lt_Gc_fam=as(ss$Ltlist[["Lt_Gc_fam"]],"dgTMatrix"),
                          Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
                          Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix")
                          ),

            fam_chol_matern_gmrf_sgps =  list(
                          y=y,
                          X=X,
                          Za=ss$Zlist[["Za"]],
                          Zc_fam=ss$Zlist[["Zc_fam"]],
                          Zc_par=ss$Zlist[["Zc_par"]],
                          Zc_sib=ss$Zlist[["Zc_sib"]],
                          Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
                          Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
                          Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
                          spdeMatrices=spdeMatrices,
                          A=A
                          ),

            fam_chol_matern_gmrf_gps =  list(
                          y=y,
                          X=X,
                          Za=ss$Zlist[["Za"]],
                          Zc_par=ss$Zlist[["Zc_par"]],
                          Zc_sib=ss$Zlist[["Zc_sib"]],
                          Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
                          Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
                          Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix")
                          ),

              fam_chol_matern_gmrf_sgps_logit= list(
                            y=y,
                            X=X,
                            Za=ss$Zlist[["Za"]],
                            Zc_fam=ss$Zlist[["Zc_fam"]],
                            Zc_par=ss$Zlist[["Zc_par"]],
                            Zc_sib=ss$Zlist[["Zc_sib"]],
                            Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
                            Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
                            Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
                            spdeMatrices=spdeMatrices,
                            A=A
                            ),

               fam_chol_matern_gmrf_sgps_overlogit=list(
                           y=y,
                           X=X,
                           Za=ss$Zlist[["Za"]],
                           Zc_fam=ss$Zlist[["Zc_fam"]],
                           Zc_par=ss$Zlist[["Zc_par"]],
                           Zc_sib=ss$Zlist[["Zc_sib"]],
                           Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
                           Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
                           Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
                           spdeMatrices=spdeMatrices,
                           A=A
                           ),

               fam_chol_matern_gmrf_sgps_overpois=list(
                          y=y,
                          X=X,
                          Za=ss$Zlist[["Za"]],
                          Zc_fam=ss$Zlist[["Zc_fam"]],
                          Zc_par=ss$Zlist[["Zc_par"]],
                          Zc_sib=ss$Zlist[["Zc_sib"]],
                          Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
                          Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
                          Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
                          spdeMatrices=spdeMatrices,
                          A=A
                          )
        ) 
   return(d)
} 
#data<-data_f(ss, "sgfps", "y")

param_f<-function(ss, data, tmplt,strt.fixed){
  switch(tmplt,
             fam_chol_matern_gmrf_sgfps=    list(
                            log_sdvc_a = 0.1,
                            log_sdvc_c_fam = log(sqrt(0.1)),
                            log_sdvc_c_par = log(sqrt(0.1)),
                            log_sdvc_c_sib = 2,
                            log_sdvc_res   = log(sqrt(0.1)),
                            ua             = rep(0.0, ncol(ss$Zlist$Za)),
                            uc_fam         = rep(0.0, ncol(ss$Zlist$Zc_fam)),
                            uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                            uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                            beta           = strt.fixed,
                            log_tau        = 0,
                            log_kappa      = 0,
                            x              = rep(0.0, nrow(data$spdeMatrices$M0))
                           ),

             fam_chol_matern_gmrf_gfps=    list(
                            log_sdvc_a = 0.1,
                            log_sdvc_c_fam = log(sqrt(0.1)),
                            log_sdvc_c_par = log(sqrt(0.1)),
                            log_sdvc_c_sib = 2,
                            log_sdvc_res   = log(sqrt(0.1)),
                            ua             = rep(0.0, ncol(ss$Zlist$Za)),
                            uc_fam         = rep(0.0, ncol(ss$Zlist$Zc_fam)),
                            uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                            uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                            beta           = strt.fixed 
                           ),

             fam_chol_matern_gmrf_sgps=    list(
                            log_sdvc_a = 0.1,
                            log_sdvc_c_par = log(sqrt(0.1)),
                            log_sdvc_c_sib = 2,
                            log_sdvc_res   = log(sqrt(0.1)),
                            ua             = rep(0.0, ncol(ss$Zlist$Za)),
                            uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                            uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                            beta           = strt.fixed,
                            log_tau        = 0,
                            log_kappa      = 0,
                            x              = rep(0.0, nrow(data$spdeMatrices$M0))
                           ),

             fam_chol_matern_gmrf_gps=    list(
                            log_sdvc_a = 0.1,
                            log_sdvc_c_par = log(sqrt(0.1)),
                            log_sdvc_c_sib = 2,
                            log_sdvc_res   = log(sqrt(0.1)),
                            ua             = rep(0.0, ncol(ss$Zlist$Za)),
                            uc_fam         = rep(0.0, ncol(ss$Zlist$Zc_fam)),
                            uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                            uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                            beta           = strt.fixed 
                           ),

               fam_chol_matern_gmrf_sgps_logit=list(
                                  log_sdvc_a = log(sqrt(0.1)),
                                  log_sdvc_c_par = log(sqrt(0.2)),
                                  log_sdvc_c_sib = log(sqrt(0.2)),
                                  #log_sdvc_res   = log(sqrt(0.2)),
                                  ua             = rep(0.0, ncol(ss$Zlist$Za)),
                                  uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                                  uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                                  #res            = rep(0.0, nrow(ss$Zlist$Zc_sib)),
                                  beta           = strt.fixed,
                                  log_tau        = 0,
                                  log_kappa      = 0,
                                  x              = rep(0.0, nrow(data$spdeMatrices$M0))
                                 ),

               fam_chol_matern_gmrf_sgps_overlogit = list(
                                  log_sdvc_a = log(sqrt(.35)),
                                  log_sdvc_c_par = log(sqrt(.1)),
                                  log_sdvc_c_sib = log(sqrt(.2)),
                                  log_sdvc_res   = log(sqrt(.3)),
                                  ua             = runif(ncol(ss$Zlist$Za)),    #rep(0.0, ncol(ss$Zlist$Za)),
                                  uc_par         = runif(ncol(ss$Zlist$Zc_par)),#rep(0.0, ncol(ss$Zlist$Zc_par)),
                                  uc_sib         = runif(ncol(ss$Zlist$Zc_sib)),#rep(0.0, ncol(ss$Zlist$Zc_sib)),
                                  res            = runif(nrow(ss$Zlist$Zc_sib)),#rep(0.0, nrow(ss$Zlist$Zc_sib)),
                                  beta           = strt.fixed,                #rep(0, ncol(X)),
                                  log_tau        = 0,
                                  log_kappa      = 0,
                                    x              = rep(0.0, nrow(data$spdeMatrices$M0))
                                   ),

               fam_chol_matern_gmrf_sgps_overpois=list(
                                  log_sdvc_a = log(sqrt(0.1*500)),
                                  log_sdvc_c_par = log(sqrt(0.1*500)),
                                  log_sdvc_c_sib = log(sqrt(0.1*500)),
                                  log_sdvc_res   = log(sqrt(0.8*500)),
                                  ua             = rep(0.0, ncol(ss$Zlist$Za)),
                                  uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                                  uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                                  res            = rep(0.0, nrow(ss$Zlist$Zc_sib)),
                                  beta           = strt.fixed,
                                  log_tau        = -1.58,
                                  log_kappa      = 1.4,
                                  x              = rep(0.0, nrow(data$spdeMatrices$M0))
                                 )
            )
}


## testing

#silent<-TRUE
#vc_struct<-"sgps"
#ycol<-"y"
#chk<-fit_fam(ss=ss, vc_struct="sgps", ycol="y",silent=T)
#chk2<-fit_fam(ss=ss, vc_struct="sgfps", ycol="y",silent=T)
#chk3<-fit_fam(ss=ss, vc_struct="sgps", ycol="y",silent=T)
#chk4<-fit_fam(ss=ss, vc_struct="sgps", ycol="y",silent=T)
#class(data_f(ss, vc_struct))
#rel<-data.table(re=names(chk$rep$par.random), val=chk$rep$par.random)
#rel[,var(val),by=re]
#chk$opt$p
#
#1-sum(chk$rep$value[-1])
#ss$d[,.(u_a=unique(u_a),uahat=rel[re=="ua",val])][,plot(u_a,uahat)]
#ss$d[,.(u_sib=unique(u_sib),uahat=rel[re=="uc_sib",val])][,plot(u_sib,uahat)]
#ss$d[,.(u_par=unique(u_par),uahat=rel[re=="uc_par",val])][,plot(u_par,uahat)]
##ss$d[,.(u_s=unique(u_s),uahat=rel[re=="x",val])][,plot(u_par,uahat)]
#
#tau.hat<-exp(chk$opt$par['log_tau'])
#estdt<-data.table(
#  ss$d[,.(set,emprel, ttype,pid,empreltype,lat,lon,tru_l=y-eps)],
#  ua=as.vector(ss$Zlist$Za %*% rel[re=="ua",val]),
#  uc_par=as.vector(ss$Zlist$Zc_par %*% rel[re=="uc_par",val]),
#  uc_sib=as.vector(ss$Zlist$Zc_sib %*% rel[re=="uc_sib",val]),
#  u_s=as.vector(ss$Zlist$Zc_fam %*% ss$sp$A %*% rel[re=="x",val]/tau.hat)
#)[,mypred:=ua + uc_par + uc_sib + u_s]
#
#
#unsc<-as.vector(ss$sp$A%*%ss$sp$u)
#xhat<-rel[re=="x",val]
#
#plot(ss$d[,u_s],as.vector(ss$Zlist$Zc_fam%*%ss$sp$A%*%xhat)/tau.hat)

  #strt.fixed<-unname(coef(glm(y~X-1,family=binomial("logit"))))




fit_fam_sgps_logit<-function(ss,ycol="yb"){

  y=ss$d[,ycol,with=F][[1]]
  X<-matrix(rep(1, nrow(ss$d)))
  spde<-inla.spde2.matern(ss$sp$mesh, alpha=0.2)
  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
  A<-ss$sp$A
  
  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
   
  data<-list(
             y=y,
             X=X,
             Za=ss$Zlist[["Za"]],
             Zc_fam=ss$Zlist[["Zc_fam"]],
             Zc_par=ss$Zlist[["Zc_par"]],
             Zc_sib=ss$Zlist[["Zc_sib"]],
             Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
             Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
             Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
             spdeMatrices=spdeMatrices,
             A=A
             )
  parameters <- list(
                     log_sdvc_a = log(sqrt(0.1)),
                     log_sdvc_c_par = log(sqrt(0.2)),
                     log_sdvc_c_sib = log(sqrt(0.2)),
                     #log_sdvc_res   = log(sqrt(0.2)),
                     ua             = rep(0.0, ncol(ss$Zlist$Za)),
                     uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                     uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                     #res            = rep(0.0, nrow(ss$Zlist$Zc_sib)),
                     beta           = rep(0, ncol(X)),
                     log_tau        = 0,
                     log_kappa      = 0,
                     x              = rep(0.0, nrow(data$spdeMatrices$M0))
                    )
  #compile("fam_chol_matern_gmrf_sgps_logit.cpp")
  #dyn.load(dynlib("fam_chol_matern_gmrf_sgps_logit"))
  print("fitting logit model") 
  obj_chol <- MakeADFun(data, 
                        parameters, 
                        random=c("x", "ua", "uc_par", "uc_sib"),  # , "res"
                        DLL="fam_chol_matern_gmrf_sgps_logit",
                        silent=T)
   
  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
  if(!inherits(opt_chol, "try-error")){
    rep_chol <- try(sdreport(obj_chol))
  }

  if(inherits(rep_chol, "try-error")){
    rt<-data.table(par=NA, est=NA, sd=NA, how="logit")
  } else {
    rt<- summ(rep_chol)[,how:="logit"]
  } 
  return(list(rep=rep_chol,opt=opt_chol,rt=rt))
}

fit_fam_sgps_overlogit<-function(ss,ycol="yb",silent=T){

  y=ss$d[,ycol,with=F][[1]]
  X<-matrix(rep(1, nrow(ss$d)))
  spde<-inla.spde2.matern(ss$sp$mesh, alpha=0.2)
  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
  A<-ss$sp$A
  
  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
   
  data<-list(
             y=y,
             X=X,
             Za=ss$Zlist[["Za"]],
             Zc_fam=ss$Zlist[["Zc_fam"]],
             Zc_par=ss$Zlist[["Zc_par"]],
             Zc_sib=ss$Zlist[["Zc_sib"]],
             Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
             Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
             Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
             spdeMatrices=spdeMatrices,
             A=A
             )
  strt.fixed<-unname(coef(glm(y~X-1,family=binomial("logit"))))
  parameters <- list(
                     log_sdvc_a = log(sqrt(.35)),
                     log_sdvc_c_par = log(sqrt(.1)),
                     log_sdvc_c_sib = log(sqrt(.2)),
                     log_sdvc_res   = log(sqrt(.3)),
                     ua             = runif(ncol(ss$Zlist$Za)),    #rep(0.0, ncol(ss$Zlist$Za)),
                     uc_par         = runif(ncol(ss$Zlist$Zc_par)),#rep(0.0, ncol(ss$Zlist$Zc_par)),
                     uc_sib         = runif(ncol(ss$Zlist$Zc_sib)),#rep(0.0, ncol(ss$Zlist$Zc_sib)),
                     res            = runif(nrow(ss$Zlist$Zc_sib)),#rep(0.0, nrow(ss$Zlist$Zc_sib)),
                     beta           = strt.fixed,                #rep(0, ncol(X)),
                     log_tau        = 0,
                     log_kappa      = 0,
                     x              = rep(0.0, nrow(data$spdeMatrices$M0))
                    )
  #compile("fam_chol_matern_gmrf_sgps_overlogit.cpp")
  #dyn.load(dynlib("fam_chol_matern_gmrf_sgps_overlogit"))

  print("fitting over-dispersed logit model") 
  obj_chol <- MakeADFun(data, 
                        parameters, 
                        random=c("x", "ua", "uc_par", "uc_sib", "res"),  
                        DLL="fam_chol_matern_gmrf_sgps_overlogit",
                        silent=silent)
   
  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
  if(!inherits(opt_chol, "try-error")){
    rep_chol <- try(sdreport(obj_chol))
  }

  if(inherits(rep_chol, "try-error")){
    rt<-data.table(par=NA, est=NA, sd=NA, how="overlogit")
  } else {
    rt<- summ(rep_chol)[,how:="overlogit"]
  } 
  return(list(rep=rep_chol,opt=opt_chol,rt=rt))
}



fit_fam_sgps_pois<-function(ss,ycol="yb"){

  y=ss$d[,ycol,with=F][[1]]
  X<-matrix(rep(1, nrow(ss$d)))
  spde<-inla.spde2.matern(ss$sp$mesh, alpha=0.2)
  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
  A<-ss$sp$A
  
  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
   
  data<-list(
             y=y,
             X=X,
             Za=ss$Zlist[["Za"]],
             Zc_fam=ss$Zlist[["Zc_fam"]],
             Zc_par=ss$Zlist[["Zc_par"]],
             Zc_sib=ss$Zlist[["Zc_sib"]],
             Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
             Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
             Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
             spdeMatrices=spdeMatrices,
             A=A
             )
  parameters <- list(
                     log_sdvc_a = log(sqrt(0.1)),
                     log_sdvc_c_par = log(sqrt(0.2)),
                     log_sdvc_c_sib = log(sqrt(0.2)),
                     log_sdvc_res   = log(sqrt(0.2)),
                     ua             = rep(0.0, ncol(ss$Zlist$Za)),
                     uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
                     uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
                     res            = rep(0.0, nrow(ss$Zlist$Zc_sib)),
                     beta           = rep(0, ncol(X)),
                     log_tau        = 0,
                     log_kappa      = 0,
                     x              = rep(0.0, nrow(data$spdeMatrices$M0))
                    )

  print("fitting poisson-approx model") 
  #compile("fam_chol_matern_gmrf_sgps_pois.cpp")
  #dyn.load(dynlib("fam_chol_matern_gmrf_sgps_pois"))
  obj_chol <- MakeADFun(data, 
                        parameters, 
                        random=c("x", "ua", "uc_par", "uc_sib", "res"), 
                        DLL="fam_chol_matern_gmrf_sgps_pois",
                        silent=T)
   
  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
  if(!inherits(opt_chol, "try-error")){
    rep_chol <- try(sdreport(obj_chol))
  }


  if(inherits(rep_chol, "try-error")){
    rt<-data.table(par=NA, est=NA, sd=NA, how="pois")
  } else {
    rt<- summ(rep_chol)[,how:="pois"]
  } 
  return(list(rep=rep_chol,opt=opt_chol,rt=rt))
}






#fit_fam_sgfps<-function(ss,ycol="y"){
#
#  y=ss$d[,ycol,with=F][[1]]
#  X<-matrix(rep(1, nrow(ss$d)))
#  spde<-inla.spde2.pcmatern(ss$sp$mesh, prior.range=c(0.5,0.5),prior.sigma=c(0.5,0.5))
#  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
#  A<-ss$sp$A
#  
#  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
#   
#  data<-list(
#             y=y,
#             X=X,
#             Za=ss$Zlist[["Za"]],
#             Zc_fam=ss$Zlist[["Zc_fam"]],
#             Zc_par=ss$Zlist[["Zc_par"]],
#             Zc_sib=ss$Zlist[["Zc_sib"]],
#             Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
#             Lt_Gc_fam=as(ss$Ltlist[["Lt_Gc_fam"]],"dgTMatrix"),
#             Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
#             Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
#             spdeMatrices=spdeMatrices,
#             A=A
#             )
#  parameters <- list(
#                     log_sdvc_a = 0.1,
#                     log_sdvc_c_fam = log(sqrt(0.1)),
#                     log_sdvc_c_par = log(sqrt(0.1)),
#                     log_sdvc_c_sib = 2,
#                     log_sdvc_res   = log(sqrt(0.1)),
#                     ua             = rep(0.0, ncol(ss$Zlist$Za)),
#                     uc_fam         = rep(0.0, ncol(ss$Zlist$Zc_fam)),
#                     uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
#                     uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
#                     beta           = rep(0, ncol(X)),
#                     log_tau        = 0,
#                     log_kappa      = 0,
#                     x              = rep(0.0, nrow(data$spdeMatrices$M0))
#                    )
#  obj_chol <- MakeADFun(data, 
#                        parameters, 
#                        random=c("x", "ua", "uc_fam", "uc_par", "uc_sib"), 
#                        DLL="fam_chol_matern_gmrf_sgfps")
#   
#  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
#  if(!inherits(opt_chol, "try-error")){
#    rep_chol <- try(sdreport(obj_chol))
#  }
#  if(inherits(rep_chol, "try-error")){
#    rt<-data.table(par=c("range", "vc_a", "vc_c", "vc_res"), est=NA, sd=NA, how="chol")
#  } else {
#    rt<- summ(rep_chol)[,how:="chol"]
#  } 
#  return(list(rep=rep_chol,opt=opt_chol))
#}
#
#fit_fam_gfps<-function(ss,ycol="y"){
#
#  y=ss$d[,ycol,with=F][[1]]
#  X<-matrix(rep(1, nrow(ss$d)))
#  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
#   
#  data<-list(
#             y=y,
#             X=X,
#             Za=ss$Zlist[["Za"]],
#             Zc_fam=ss$Zlist[["Zc_fam"]],
#             Zc_par=ss$Zlist[["Zc_par"]],
#             Zc_sib=ss$Zlist[["Zc_sib"]],
#             Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
#             Lt_Gc_fam=as(ss$Ltlist[["Lt_Gc_fam"]],"dgTMatrix"),
#             Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
#             Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix")
#             )
#  parameters <- list(
#                     log_sdvc_a = 0.1,
#                     log_sdvc_c_fam = log(sqrt(0.1)),
#                     log_sdvc_c_par = log(sqrt(0.1)),
#                     log_sdvc_c_sib = 2,
#                     log_sdvc_res   = log(sqrt(0.1)),
#                     ua             = rep(0.0, ncol(ss$Zlist$Za)),
#                     uc_fam         = rep(0.0, ncol(ss$Zlist$Zc_fam)),
#                     uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
#                     uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
#                     beta           = rep(0, ncol(X)),
#                    )
#  obj_chol <- MakeADFun(data, 
#                        parameters, 
#                        random=c("x", "ua", "uc_fam", "uc_par", "uc_sib"), 
#                        DLL="fam_chol_matern_gmrf_gfps")
#   
#  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
#  if(!inherits(opt_chol, "try-error")){
#    rep_chol <- try(sdreport(obj_chol))
#  }
#  if(inherits(rep_chol, "try-error")){
#    rt<-data.table(par=c("range", "vc_a", "vc_c", "vc_res"), est=NA, sd=NA, how="chol")
#  } else {
#    rt<- summ(rep_chol)[,how:="chol"]
#  } 
#  return(list(rep=rep_chol,opt=opt_chol))
#}
#
#
#
#
#fit_fam_sgps<-function(ss,ycol="y"){
#
#  y=ss$d[,ycol,with=F][[1]]
#  X<-matrix(rep(1, nrow(ss$d)))
#  spde<-inla.spde2.matern(ss$sp$mesh, alpha=0.2)
#  spdeMatrices<-spde$param.inla[c("M0", "M1", "M2")]
#  A<-ss$sp$A
#  
#  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
#   
#  data<-list(
#             y=y,
#             X=X,
#             Za=ss$Zlist[["Za"]],
#             Zc_fam=ss$Zlist[["Zc_fam"]],
#             Zc_par=ss$Zlist[["Zc_par"]],
#             Zc_sib=ss$Zlist[["Zc_sib"]],
#             Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
#             Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
#             Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix"),
#             spdeMatrices=spdeMatrices,
#             A=A
#             )
#  parameters <- list(
#                     log_sdvc_a = log(sqrt(0.3)),
#                     log_sdvc_c_par = log(sqrt(0.1)),
#                     log_sdvc_c_sib = log(sqrt(0.1)),
#                     log_sdvc_res   = log(sqrt(0.3)),
#                     ua             = rep(0.0, ncol(ss$Zlist$Za)),
#                     uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
#                     uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
#                     beta           = rep(0, ncol(X)),
#                     log_tau        = 0,
#                     log_kappa      = 0,
#                     x              = rep(0.0, nrow(data$spdeMatrices$M0))
#                    )
#  obj_chol <- MakeADFun(data, 
#                        parameters, 
#                        random=c("x", "ua", "uc_par", "uc_sib"), 
#                        DLL="fam_chol_matern_gmrf_sgps")
#   
#  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
#  if(!inherits(opt_chol, "try-error")){
#    rep_chol <- try(sdreport(obj_chol))
#  }
#  if(inherits(rep_chol, "try-error")){
#    rt<-data.table(par=c("range", "vc_a", "vc_c", "vc_res"), est=NA, sd=NA, how="chol")
#  } else {
#    rt<- summ(rep_chol)[,how:="chol"]
#  } 
#  return(list(rep=rep_chol,opt=opt_chol))
#}
#
#
#
#fit_fam_gps<-function(ss,ycol="y"){
#
#  y=ss$d[,ycol,with=F][[1]]
#  X<-matrix(rep(1, nrow(ss$d)))
#  ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps|^y$", names(ss$d),value=T)]
#   
#  data<-list(
#             y=y,
#             X=X,
#             Za=ss$Zlist[["Za"]],
#             Zc_par=ss$Zlist[["Zc_par"]],
#             Zc_sib=ss$Zlist[["Zc_sib"]],
#             Lt_Ga=as(ss$Ltlist[["Lt_Ga"]],"dgTMatrix"),
#             Lt_Gc_par=as(ss$Ltlist[["Lt_Gc_par"]],"dgTMatrix"),
#             Lt_Gc_sib=as(ss$Ltlist[["Lt_Gc_sib"]],"dgTMatrix")
#             )
#  parameters <- list(
#                     log_sdvc_a = log(sqrt(0.3)),
#                     log_sdvc_c_par = log(sqrt(0.1)),
#                     log_sdvc_c_sib = log(sqrt(0.1)),
#                     log_sdvc_res   = log(sqrt(0.3)),
#                     ua             = rep(0.0, ncol(ss$Zlist$Za)),
#                     uc_par         = rep(0.0, ncol(ss$Zlist$Zc_par)),
#                     uc_sib         = rep(0.0, ncol(ss$Zlist$Zc_sib)),
#                     beta           = rep(0, ncol(X))
#                    )
#  obj_chol <- MakeADFun(data, 
#                        parameters, 
#                        random=c("ua", "uc_par", "uc_sib"), 
#                        DLL="fam_chol_matern_gmrf_gps")
#   
#  opt_chol <- try(nlminb(obj_chol$par, obj_chol$fn, obj_chol$gr))
#  if(!inherits(opt_chol, "try-error")){
#    rep_chol <- try(sdreport(obj_chol))
#  }
#  if(inherits(rep_chol, "try-error")){
#    rt<-data.table(par=c("range", "vc_a", "vc_c", "vc_res"), est=NA, sd=NA, how="chol")
#  } else {
#    rt<- summ(rep_chol)[,how:="chol"]
#  } 
#  return(list(rep=rep_chol,opt=opt_chol))
#}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
