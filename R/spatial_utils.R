
invlogit<-function(x){ exp(x) / (1 + exp(x))}
logit<-function(x){ log(x/(1-x))}


#'Calculate a single summary variance component for covariance matrix with non-constant variances/covariances.
#'
#'@param vc variance component
#'@returns the summary statistic
gowfac<-function(vc=NULL, Z=NULL, G=NULL){
  if(!missing(G) & !missing(Z)){
    gows<-sum(rep(diag(G) - 
		  rowMeans(G),times=colSums(Z)))
  } else if(!missing(Z)){
    gows<-sum(diag(Z%*%t(Z)) - 
      	     rowMeans(Z%*%t(Z)))
  } else {
    gows<-sum(rep(diag(G) - 
	         	  rowMeans(G)))
  }
  if(!missing(vc)) gows<-vc*gows
  return(gows)
}


fillin_nearest_neighb<-function(lgeo_opt, w){
  geo<-copy(lgeo_opt)
  geo[,rn:=1:.N]
  whch<-rowSums(w)
  fillin<-which(whch==0)
  if(length(fillin)==0) break("no rows/columns need filling")
  whch<-which(whch!=-9)
  ivec<-unlist(lapply(whch, function(x) which(w[,x]!=0)))
  jvec<-rep(whch, times=rowSums(w))
  bst<-c()
  for(jj in 1:length(fillin)){
    dchk<-merge(geo[rn==fillin[jj],.(rn,X,Y,tmp=1)],
          geo[rn!=fillin[jj],.(rn,X,Y,tmp=1)],
          allow.cartesian=T,
          by="tmp",
          suffix=c("_a", "_b"))
    dchk[,dist:=sqrt((X_a - X_b)^2 + (Y_a - Y_b)^2)]
    bst[jj]<-dchk[order(dist)][1,rn_b]
  
  }
  ivec<-c(ivec, fillin, bst)
  jvec<-c(jvec, bst, fillin)
  wnew<-sparseMatrix(i=ivec,j=jvec)
  return(list(wnew=wnew,
	      fillin=fillin,
	      bstneighb=bst))
}

gen_fam_ids<-function(n_quad = 5){

    id_dt<-data.table(set = rep(1:(n_quad), times = rep(c(4), times=c(n_quad))),    
                      emprel = c(rep(c(1,2,3,3),times=n_quad)), 
                      ttype = rep(c("quad"), times = c(n_quad*4)))

    id_dt[, `:=`(pid, 1:.N)]
    setorderv(id_dt, c("ttype", "set", "emprel","pid"))

    return(id_dt)
}



sim_fam_nonspatial<-function(
                             n_quad = 5                
                            ,X=NULL                    
                            ,Ag = 0.2                  
                            ,C_par = 0.1               
                            ,C_sib = 0.1               
                            ,C_fam = 0.1               
                            ,S=0.2                     
                            ,mu = 0                    
                             ){
    E<-1-Ag-C_fam - C_sib - C_par -S
    id_dt<-gen_fam_ids(n_quad) 
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
    
    return(list(
      d=id_dt,
      Zlist=mlists$Zlist,
      Glist=mlists$Glist,
      Ltlist=mlists$Ltlist,
      Ag=Ag,
      C_par=C_par,
      C_fam=C_fam,
      C_sib=C_sib,
      S=S,
      E=E
   )) 
}



make_matrix_lists<-function(family_data,
                            ...){
    
  
    Ga_ref<-list("mz"=matrix(1,nrow=1,ncol=1),
                 "dz"=matrix(c(1, .5, .5, 1),ncol=2),
                 "trio"=matrix(c(1, .5, .5, 
                                 .5, 1, .5, 
                                 .5, .5, 1),ncol=3),
                 "quad"=matrix(c(1, 0, .5, .5, 
                                 0, 1, .5, .5, 
                                 .5, .5, 1, .5, 
                                 .5, .5, .5, 1),ncol=4, byrow=T))
  
    Zc_sib_ref<-list("mz"=matrix(1,nrow=2,ncol=1),
                     "dz"=matrix(1,nrow=2,ncol=1),
                     "trio"=matrix(c(1, 0, 
                                     0, 1, 
                                     0, 1),nrow=3,byrow=T),
                     "quad"=matrix(c(1, 0, 0, 
                                     0, 1, 0, 
                                     0, 0, 1, 
                                     0, 0, 1),nrow=4,byrow=T)
                     )
  
  
    Zc_par_ref<-list("mz"=diag(2),
                     "dz"=diag(2),
                     "trio"=diag(3),
                     "quad"=matrix(c(1, 0, 0, 
                                     1, 0, 0, 
                                     0, 1, 0, 
                                     0, 0, 1),nrow=4,byrow=T)
                     )
  
  
    Zc_fam_ref<-list("mz"=matrix(1, nrow=2),
                     "dz"=matrix(1, nrow=2),
                     "trio"=matrix(1, nrow=3),
                     "quad"=matrix(1, nrow=4)
                     )
    
    id_dt<-copy(family_data)
    id_dt[emprel %in% c(1:2),empreltype:="parent"]
    id_dt[emprel %in% 3,empreltype:="sib"]
    if(!"cc_group" %in% names(id_dt)){
      id_dt[,cc_group:=0]
    }
    ordervars<-c("cc_group", "ttype", "set", "emprel","pid")
    ordervars<-ordervars[ordervars %in% names(id_dt)] 
    setorderv(id_dt, ordervars)
    uftype<-id_dt[,unique(ttype)] 
    
    Ga<-c()
    Zc_sib<-c()
    Zc_par<-c()
    Zc_fam<-c()
    
    uftype<-id_dt[,.N,by=.(set,ttype,cc_group)][,.N,by=.(ttype,cc_group)] 
    
    for(ff in 1:uftype[,.N]){
      ffn<-uftype[ff,N]
      fftype<-uftype[ff,ttype]
      Zc_sib<-c(Zc_sib,
      	  rep(list(Zc_sib_ref[fftype][[1]]),
      	      ffn))
      Zc_par<-c(Zc_par,
      	  rep(list(Zc_par_ref[fftype][[1]]),
      	      ffn))
      Zc_fam<-c(Zc_fam,
      	  rep(list(Zc_fam_ref[fftype][[1]]),
      	      ffn))
      
      Ga<-c(Ga,rep(list(Ga_ref[fftype][[1]]),
      	     ffn))
      
    }
    
    Zc_sib<-.bdiag(Zc_sib)
    Zc_par<-.bdiag(Zc_par) 
    Zc_fam<-.bdiag(Zc_fam)
    
    Ga<-.bdiag(Ga)
    Lt_Ga<-t(Matrix::chol(Ga))
    
  return(list(
      d=id_dt,
      Zlist=list(Zc_sib=Zc_sib,Zc_par=Zc_par,Zc_fam=Zc_fam),
      Glist=list(Ga=Ga),
      Ltlist=list(Lt_Ga=Lt_Ga)
   )) 
}





vc_model<-function(
                   X
                  ,y
                  ,spatial_structure = c("car", "sar", "iid"
                  ,mlists
                  ,mat   
                  ){
   binary_ind<-ifelse(all(sort(unique(y))==c(0,1)), 1, 0)
}


fitmod<-function(s_type,
		 ycol, 
		 vc_struct, 
		 y,
		 X,
		 Zs,
		 mlists,
                 offset=NULL,
		 outcome="binary",
		 binary_ind=1,
		 A=NULL,
		 strt.fixed=0,
                 minrng=1e-3,
		 car_mats=NULL,
		 rhomax=NULL,
		 tmplt=NULL,
		 re.hess=FALSE,
		 enrolid=NULL,
		 vc_inits=NULL,
		 beta_inits=NULL){
  
  if(is.null(tmplt)){
    tmplt<-paste0("fam_chol_", s_type, "_",vc_struct, "_us")
  }
  vcl<-strsplit(vc_struct, "")[[1]]
  #vary<-d[,var(get(ycol))]
  vary<-var(y)
  if(missing(vc_inits)){
    vc_s_init<-.5*(vary)/length(vcl)   
    vc_c_sib_init<-.5*(vary)/length(vcl)   
    vc_c_par_init<-.5*(vary)/length(vcl)   
    vc_c_fam_init<-.5*(vary)/length(vcl)   
    vc_a_init<-.5*(vary)/length(vcl)   
  } else {
    for(item in 1:length(vc_inits)){
       assign(names(vc_inits[item]), unname(vc_inits[item]))
    } 
  
  }
  if(outcome %in% c("binary", "gaussian")){
    dl<-list(y = y, 
            X = X,
            binary_ind=binary_ind)
    fmod<-lm(dl$y~X-1)
    print(summary(fmod))
    strt.fixed<-coef(fmod)
    if(!missing(beta_inits)){
      strt.fixed<-beta_inits 
    }
    params<-list(beta=strt.fixed,
                 log_sdvc_res = log(sqrt(vary/2))
    	     )
    
  } else if (outcome=="count"){
    
    dl<-list(y = y, 
             X = X,
             offset=offset)
#    lmod<-lm(dl$y~X-1)
    glmd<-as.data.frame(cbind(y=dl$y, cbind(offsetx=offset, X)))
    fmod<-MASS::glm.nb(y ~ .-offsetx-1 + offset(offsetx),data=glmd,trace=4)
    fmod<-list(summary=summary(fmod), coef=fmod$coefficients,theta=fmod$theta)
    print(fmod$summary)
    strt.fixed<-fmod$coef
    params<-list(beta=strt.fixed,log_phi=log(fmod$theta))
#    rm(fmod)
    rm(glmd)
  
  } else {
    stop(paste0("outcome type: ", outcome, " not recognized.")) 
  }
  
  if(vcl[1]=="s"){
    if(s_type=="car_gmrf"){ 
      dl$Zs=as(Zs,"dgTMatrix");
      dl$car_mats=car_mats;
      dl$rhomax=rhomax
      params$logit_rho<-0  #logit((.5*maxrho - minrho)/(maxrho-minrho))
      params$log_sdvc_s<-log(sqrt(vc_s_init))
      params$x<-rep(0, nrow(car_mats$N))
    } else {
      dl$Zs=as(Zs,"dgTMatrix");
      dl$minrho<-minrho 
      dl$maxrho<-maxrho
      dl$Vw<-as(eigw$vectors, "dgTMatrix")
      dl$wj<-eigw$values
      vc_s_init<-.5*(vary)/length(vcl)   
      params$unscaled_rhocar<-0  
      params$log_sdvc_s<-log(sqrt(vc_s_init))
      params$us<-rep(0, nrow(eigw$vectors))
    }
  }
  if("g" %in%  vcl){
  #  dl$Za = mlists$Zlist[["Za"]];
    dl$Lt_Ga = as(mlists$Ltlist[["Lt_Ga"]], "dgTMatrix");
    params$log_sdvc_a = log(sqrt(vc_a_init));
  #  params$ua=rep(0,ncol(dl$Za));
    params$ua=rep(0,ncol(dl$Lt_Ga));
    #params$ua=rnorm(ncol(dl$Lt_Ga),0, sqrt(vc_g_init)/5);
  }
  if("f" %in%  vcl){
  #  dl$Lt_Gc_fam = as(ss$Ltlist[["Lt_Gc_fam"]],  "dgTMatrix");
    dl$Zc_fam = mlists$Zlist[["Zc_fam"]];
    params$log_sdvc_c_fam = log(sqrt(vc_c_fam_init));
    params$uc_fam = rep(0, ncol(dl$Zc_fam));
    #params$uc_fam = rnorm(ncol(dl$Zc_fam),0,sqrt(vc_fam_init)/5);
  }
  if("p" %in%  vcl){
  #  dl$Lt_Gc_par = as(ss$Ltlist[["Lt_Gc_par"]],  "dgTMatrix");
    dl$Zc_par = mlists$Zlist[["Zc_par"]];
    params$log_sdvc_c_par = log(sqrt(vc_c_par_init)); 
    params$uc_par = rep(0, ncol(dl$Zc_par));
    #params$uc_par = rnorm(ncol(dl$Zc_par),0,sqrt(vc_par_init)/5);
  }
  if(length(vcl) > 1 & vcl[length(vcl)]=="s"){
  #  dl$Lt_Gc_sib = as(ss$Ltlist[["Lt_Gc_sib"]], "dgTMatrix");
    dl$Zc_sib = mlists$Zlist[["Zc_sib"]];
    params$log_sdvc_c_sib = log(sqrt(vc_c_sib_init)); 
    params$uc_sib = rep(0, ncol(dl$Zc_sib));
    #params$uc_sib = rnorm(ncol(dl$Zc_sib),0,sqrt(vc_sib_init)/5);
  }
  
  re_params <- setdiff(grep("^x|^u", names(params), value = T), "unscaled_rhocar")
  
   
  #library(pryr)
#  rm(A)
#  rm(X)
#  rm(Zs); rm(dw); rm(dicdi); rm(bigA);  rm(di); rm(case_fam); 
#  rm(dsmk);  rm(i); rm(ii); rm(jj); rm(master); rm(ms); rm(tbco);
#  rm(mlists)
  #rm(ff);
  if(!grepl("_iid_",tmplt)){
    obj_chol <- MakeADFun(data = c(model = tmplt, dl), parameters = params, 
        random = re_params, silent = FALSE, DLL = "meta_model")
  }
  else{
    obj_chol <- MakeADFun(data = c(model = tmplt, dl), parameters = params, 
                 random = re_params, silent = FALSE, DLL = "meta_model",
		 map=list(logit_rho=factor(NA)))
  }
  obj_chol$env$tracepar<-TRUE
  
  fittime<-system.time(
        opt_chol <- try(nlminb(
      		   obj_chol$par, 
            		   obj_chol$fn, 
            		   obj_chol$gr,
            		   control=list(eval.max=5000,
            				iter.max=5000)
      		   ))
      	    )
   
  
  sdtime<-system.time(
    rep_chol <- try(sdreport(obj_chol,getJointPrecision=re.hess))
		      ) 
#  sdtime<-system.time(
#    rep_chol <- try(sdreport(obj_chol))
#		      ) 
  re<-data.table(eff=rep_chol$par.random, re=names(rep_chol$par.random))    
  
  odti<-as.data.table(summary(rep_chol),rownames(summary(rep_chol)))
  setnames(odti, old=names(odti)[1], new="param") 
  odti[,ll:=Estimate - qnorm(.975)*`Std.\ Error`]
  odti[,ul:=Estimate + qnorm(.975)*`Std.\ Error`]
  odti[,objective:=opt_chol$objective]
  odti[,mess:=opt_chol$message]
  odti[,fit_time:=fittime[3]]
  odti[,sd_time:=sdtime[3]]
  odti[,modeltype:=s_type] 
  
  fe<-X%*%odti[grep("beta",param),Estimate]

  indlev<-data.table(X,y, fe=as.vector(fe))
  if(!missing(enrolid)){
    indlev[,enrolid:=enrolid] 
  }

  Vnl<-list()

  if(vcl[1]=="s"){
    us<-Zs %*% odti[param=="x",Estimate] 
    indlev[,us:=as.vector(us)]
  }
  if("g" %in%  vcl){
    ua<-odti[param=="ua",Estimate]
    indlev[,ua:=as.vector(ua)]

    Vnl[["g"]]<-odti[param=="vc_a",Estimate] * (mlists$Glist$Ga)
  }
  if("f" %in%  vcl){
    ufam<-mlists$Zlist$Zc_fam %*% odti[param=="uc_fam",Estimate] 
    indlev[,ufam:=as.vector(ufam)]
    
    Vnl[["f"]]<-odti[param=="vc_c_fam",Estimate] * 
	    (mlists$Zlist$Zc_fam %*% t(mlists$Zlist$Zc_fam))
  }
  if("p" %in%  vcl){
    upar<-mlists$Zlist$Zc_par %*% odti[param=="uc_par",Estimate] 
    indlev[,upar:=as.vector(upar)]
    
    Vnl[["p"]]<-odti[param=="vc_c_par",Estimate] * 
	    (mlists$Zlist$Zc_par %*% t(mlists$Zlist$Zc_par))
  }
  if(length(vcl) > 1 & vcl[length(vcl)]=="s"){
    usib<-mlists$Zlist$Zc_sib %*% odti[param=="uc_sib",Estimate] 
    indlev[,usib:=as.vector(usib)]
    Vnl[["sib"]]<-odti[param=="vc_c_sib",Estimate] * 
	    (mlists$Zlist$Zc_sib %*% t(mlists$Zlist$Zc_sib))
  }
  Vnl[["res"]]<-sparseMatrix(i=1:nrow(X), j=1:nrow(X),x=odti[param=="vc_res",Estimate]) 
  VnlS<-Vnl[["res"]]
  for(i in setdiff(names(Vnl),"res")){
    VnlS<-VnlS + Vnl[[i]] 
  }
  rm(Vnl) 
  recols<-grep("^fe$",names(indlev))
  recols<-grep("^u|^x$",names(indlev)[(recols+1):length(names(indlev))],value=T)
  indlev[,yhat:=rowSums(.SD),.SDcols=c("fe", recols),by=.I]
  #brd<-solve(t(as.matrix(X)) %*% as.matrix(X))
  #hc1<-sparseMatrix(i=1:nrow(X),j=1:nrow(X),x=fmod$residuals^2)
  #chz<-t(as.matrix(X)) %*% hc1 %*% as.matrix(X)
  #vcovhc<-brd %*% chz %*% brd
  mean_yhat<-indlev[,mean(yhat)] 
  var_fe<-var(X %*% odti[grep("beta",param),Estimate])
  VnlS.inv<-chol2inv(chol(VnlS)) 
  txtvinv<-t(t(X) %*% VnlS.inv)
  res.fe<-indlev[,y - fe]
  scoresi<-apply(txtvinv, 2, function(x) x * res.fe)
  vc.bread<-t(X) %*% VnlS.inv %*% X 
  olsrprt<-as.data.table(tidy(fmod))

  robustse<-sqrt(diag(vcovHC(fmod,type="HC0")))
  olsrprt[,hc.se:=robustse]
  olsrprt[,hc.p:=2*pnorm(-abs(estimate/hc.se))]
  olsrprt[,param:=paste0("beta_",term)]

  return(list(odti=odti,
              rep_chol=rep_chol,
              opt_chol=opt_chol,
	      obj_chol=obj_chol,
              minrng=minrng,
#	      fmod=fmod,
	      indlev=indlev,
	      template=tmplt,
              mean_yhat=mean_yhat,
              var_fe=var_fe,
	      vc.bread=vc.bread,
	      scoresi=scoresi,
              olsrprt=olsrprt))
}



