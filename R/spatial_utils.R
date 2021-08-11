
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



chk<-function(X){
  if(missing(X)) X<-3
  return(X)
}

smile<-function(
                 X
                ,y
                ,Zs
                ,wmatd
                ,vc_structure=c("gpcs")
                ,spatial_structure = c("car", "sar", "iid")
                ,mlists
                ,enrolid=NULL
                ,vc_inits=NULL
                ,beta_inits=NULL
                ){
   
   binary_ind<-ifelse(all(sort(unique(y))==c(0,1)), 1, 0)
   if(missing(X)) X<-matrix(rep(1,nrow(ss$d),ncol=1))
   if(missing(y)) stop("response vector 'y' required") 
   
   car_mats<-NULL
   if(!missing(spatial_structure) & !missing(wmatd)){
    
       rs<-rowSums(wmatd)
       wmatplus<-as(do.call(rbind, 
           		 lapply(1:nrow(wmatd), function(i) wmatd[i,]/rs[i])), "dgTMatrix")
       
       Mpm5<-as(sparseMatrix(i=1:length(rs), 
       		      j=1:length(rs),
       		      x=1/rs^-.5),
       	 "dgTMatrix")
       
       Mpp5<-as(sparseMatrix(i=1:length(rs), 
       		      j=1:length(rs),
       		      x=1/rs^.5),
       	 "dgTMatrix")
       
       diagI<-as(sparseMatrix(i=1:length(rs), 
       		      j=1:length(rs),
       		      x=1),
       	 "dgTMatrix")
       N<-Mpm5%*% wmatplus %*% Mpp5
       car_mats<-list(diagI=diagI,N=N,Mpm5=Mpm5)
  } 
   
   s_type<-ifelse(!missing(car_mats), "car_gmrf", "matern_gmrf")
   vcl<-strsplit(vc_structure, "")[[1]]
   
   mtchvcl<-c("g", "p", "c", "s") 
   vcl<-mtchvcl[mtchvcl %in% vcl]
   vc_struct<-paste0(vcl,collapse="") 
   tmplt<-paste0("fam_chol_car_gmrf_",vc_structure,"_us")
   
   if(grepl("iid",spatial_structure)){
     tmplt<-gsub("_car_","_iid_",tmplt)
   } 
   if(grepl("sar",spatial_structure)){
     tmplt<-gsub("_car_","_sar_",tmplt)
   } 
   if(grepl("matern",spatial_structure)){
     tmplt<-gsub("_car_","_matern_",tmplt)
   } 
   
   ### temp fix until template filenames are updated
   tmplt<-gsub("_gpcs_","_sgps_",tmplt) 
   tmplt<-gsub("_gpc_" ,"_gps_",tmplt) 
   tmplt<-gsub("_gps_" ,"_sgp_",tmplt) 
   tmplt<-gsub("_gcs_" ,"_sgs_",tmplt) 
   tmplt<-gsub("_gp_"  ,"_gp_",tmplt) 
   tmplt<-gsub("_gs_"  ,"_sg_",tmplt) 
   tmplt<-gsub("_gc_"  ,"_gs_",tmplt) 
   tmplt<-gsub("_ps_"  ,"_ps_",tmplt) 
   tmplt<-gsub("_cs_"  ,"_ss_",tmplt) 
   tmplt<-gsub("_pc_"  ,"_pc_",tmplt) 
   tmplt<-gsub("_g_"   ,"_g_",tmplt) 
   tmplt<-gsub("_p_"   ,"_p_",tmplt) 
   tmplt<-gsub("_s_"   ,"_s_",tmplt) 
    
   ### fit the model ###
   mqc<-fitmod(s_type=s_type, 
     	         vcl=vcl, 
     	         y=y,
                 X=X,
                 Zs=Zs,
                 mlists=mlists,
     	         car_mats=car_mats,
                 rhomax=1,
                 tmplt=tmplt,
   	         outcome="gaussian",
                 binary_ind=binaryind)
     
     rhohat<-mqc$odti[param=="rho",Estimate]
     vc_s_hat<-mqc$odti[param=="vc_s",Estimate]
     sigma_hat<-vc_s_hat*chol2inv(chol(Mpm5%*%(diagI - rhohat*N)%*%Mpm5))
     gvc_s<-gowfac(Z=Zs,G=sigma_hat)  / nrow(Zs)
             
     ##### summarize output ########
     cv<-mqc$rep_chol$cov
     rownames(cv)<-colnames(cv)<-names(mqc$rep_chol$value)
     params<-c("vc_a", "vc_c_par", "vc_c_sib", "vc_s","rho")
     inclparams<-match(rownames(cv),params)
     params<-params[inclparams[!is.na(inclparams)]]
     vc.cov<-cv[params,params]
            
     vc.vec<-mqc$odti[match(params,param),Estimate]
     names(vc.vec)<-params
     vc.vec<-t(t(vc.vec))
     print("calculating liability scale variance components")
     if(!any(is.na(vc.cov)) & K < 1){
       lia_est_se<-calculate_se_lia(vc_s_hat,vc.vec, vc.cov, gvc_s, sigma_hat, N,Mpm5, diagI,Zs, khat=mean(ss$d[,mean(yb)]),
                                    s_type="car",binary=TRUE)
       
       setnames(lia_est_se$outrprt,old=c("est", "se"), new=c("Estimate", "Std.Error")) 
       #mqc$odti<-merge(mqc$odti, 
       #                lia_est_se$outrprt[,unique(.SD)][!is.na(Estimate)],
       #                all=T,suffix=c("", "2",sort=F),by="param")
       #mqc$odti[!is.na(Estimate2) & is.na(Estimate),`:=`(Estimate=Estimate2,`Std. Error`=`Std. Error2`)]
       mqc$lia_est_se<-lia_est_se
     }
      mqc$odti[param=="x",param:=paste0("location_",1:nrow(wmatd))]
      mqc$odti[param=="beta",param:=mqc$olsrprt[,param]]
      #mqc$odti[grep("beta",param)][,pval:=2*pnorm(-abs(Estimate/`Std. Error`))][]
     
     output<-list(tmb_rpt=mqc$odti[!grepl("lia",param)]
                  ,vc_rpt=mqc$lia_est_se
                  ,indlevel_dt=mqc$indlev
                  ,ols_rpt=mqc$olsrpt
                  ,tmb_mod=mqc$opt_chol
                  )
     return(output)
   
}


fitmod<-function(s_type,
		 ycol, 
		 vcl, 
		 y,
		 X,
		 Zs,
		 mlists,
                 offset=NULL,
		 outcome="binary",
		 binary_ind=binary_ind,
		 car_mats=NULL,
		 rhomax=1,
		 tmplt=NULL,
		 re.hess=FALSE,
		 enrolid=NULL,
		 vc_inits=NULL,
		 beta_inits=NULL){
  
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
    
  
  if("s" %in% vcl){
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
    dl$Lt_Ga = as(mlists$Ltlist[["Lt_Ga"]], "dgTMatrix");
    params$log_sdvc_a = log(sqrt(vc_a_init));
    params$ua=rep(0,ncol(dl$Lt_Ga));
  }
  if("f" %in%  vcl){
    dl$Zc_fam = mlists$Zlist[["Zc_fam"]];
    params$log_sdvc_c_fam = log(sqrt(vc_c_fam_init));
    params$uc_fam = rep(0, ncol(dl$Zc_fam));
  }
  if("p" %in%  vcl){
    dl$Zc_par = mlists$Zlist[["Zc_par"]];
    params$log_sdvc_c_par = log(sqrt(vc_c_par_init)); 
    params$uc_par = rep(0, ncol(dl$Zc_par));
  }
  if("c" %in% vcl){
    dl$Zc_sib = mlists$Zlist[["Zc_sib"]];
    params$log_sdvc_c_sib = log(sqrt(vc_c_sib_init)); 
    params$uc_sib = rep(0, ncol(dl$Zc_sib));
  }
  
  re_params <- setdiff(grep("^x|^u", names(params), value = T), "unscaled_rhocar")
  
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

  if("s" %in% vcl){
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
  if("c" %in% vcl){
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
	      indlev=indlev,
	      template=tmplt,
              mean_yhat=mean_yhat,
              var_fe=var_fe,
	      vc.bread=vc.bread,
	      scoresi=scoresi,
              olsrprt=olsrprt))
}

calculate_se_lia<-function(vc_s_hat,vc.vec, vc.cov,gvc_s,sigma_hat,N,Mpm5,diagI,Zs,khat,s_type,binary=TRUE){

  #### GVC_S SE !!! ###
 if(s_type %in% c("sar", "car")){
   if(s_type=="sar"){

     pSigmasprho<- -vc_s_hat*sigma_hat %*% (-2*(Mpm5%*%Mpm5)%*%(diagI - rhohat*N)%*%N) %*% sigma_hat

   } else if(s_type=="car"){

     pSigmasprho<- -vc_s_hat*sigma_hat %*%  (-Mpm5 %*% Mpm5 %*% N ) %*% sigma_hat

   }
   pSigmaspvcs<- (1/vc_s_hat)*sigma_hat
   pmat<-diag(ncol(sigma_hat)) - 1/nrow(sigma_hat) * matrix(1, ncol(sigma_hat), nrow(sigma_hat)) 
   ZtZ<-t(Zs) %*% Zs

   gpvcs<-sum(diag((t(pmat) %*% ZtZ ) %*% pSigmaspvcs))/nrow(Zs)
   gprho<-sum(diag((t(pmat) %*% ZtZ) %*% pSigmasprho))/nrow(Zs)
   gprime.gow<-matrix(c(gpvcs,gprho 
                         ), ncol=1,nrow=2)
  
    
   g0<-matrix(c(
                1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, gprime.gow[1,1], gprime.gow[2,1]
                ),nrow=4,ncol=5,byrow=T)
   rownames(g0)<-c("vc_a", "vc_c_par", "vc_c_sib", "gvc_s")
   colnames(g0)<-c("vc_a", "vc_c_par", "vc_c_sib", "vc_s", "rho")
   g0<-g0[,rownames(vc.vec)] 

   vc.cov2<-g0 %*% vc.cov %*% t(g0)
   vc.cov2<-as.matrix(vc.cov2[which(diag(vc.cov2)!=0),which(diag(vc.cov2)!=0)])
   vc.vec2<-matrix(c(vc.vec[-which(rownames(vc.vec) %in% c("vc_s", "rho")),1],gvc_s))
   rownames(vc.vec2)<-c(rownames(vc.vec)[-which(rownames(vc.vec) %in% c("vc_s", "rho"))], "gvc_s")
   rownames(vc.cov2)<-colnames(vc.cov2)<-rownames(vc.vec2)

 } else if (s_type=="iid"){

   vc.cov2<-vc.cov
   vc.vec2<-vc.vec
   whchchng<-which(rownames(vc.vec2)=="vc_s")
   rownames(vc.vec2)[whchchng]<-"gvc_s"
   rownames(vc.cov2)<-colnames(vc.cov2)<-rownames(vc.vec2)

 } else {
   stop("s_type not recognized.  should be sar,car, or iid")
 }
  
  g1<-matrix(c(0, 1, 0, 1,
               0.5, 0, 0, 1,
               0.5, 0, 1, 1,
               0, 0, 0, 1),
             byrow=T,nrow=4,ncol=4)
  
  rownames(g1)<-c("vc_pp", "vc_ps", "vc_ss", "gvc_s")
  colnames(g1)<-c("vc_a", "vc_c_par", "vc_c_sib", "gvc_s")
  g1<-as.matrix(g1[,rownames(vc.vec2)] )
  robs.vec<-g1 %*% vc.vec2
  robs.cov<-g1 %*% vc.cov2 %*% t(g1)
  rownames(robs.vec)<-c("vc_pp", "vc_ps", "vc_ss", "gvc_s")
  rownames(robs.cov)<-colnames(robs.cov)<-c("vc_pp", "vc_ps", "vc_ss", "gvc_s")
  
  
  if(binary){
    rlia.vec<-sapply(lapply(1:nrow(robs.vec), function(i){
      to_lia2(robs.vec[i,1], robs.cov[i,i], khat)
    }), "[[", 1)
    rlia.vec<-t(t(rlia.vec))
    rownames(rlia.vec)<-rownames(robs.vec)
    gp2<-sapply(lapply(1:nrow(robs.vec), function(i){
      to_lia2(robs.vec[i,1], robs.cov[i,i], khat)
    }), "[[", 4)
    gp2<-diag(gp2)
    
    rlia.cov<-t(gp2) %*% robs.cov %*% (gp2)
    rownames(rlia.cov)<-colnames(rlia.cov)<-c("vc_pp", "vc_ps", "vc_ss", "gvc_s")
    
    
    
    g3<-matrix(c(
                 0, 2, 0, -2, 
                 1, 0, 0, -1, 
                 0, -1, 1, 0, 
                 0, 0, 0, 1
                 ),ncol=4,nrow=4,byrow=T)

    colnames(g3)<-c("vc_pp", "vc_ps", "vc_ss", "gvc_s")
    rownames(g3)<-c("vc_a", "vc_c_par", "vc_c_sib", "gvc_s")
    g3<-g3[rownames(vc.vec2),]
    if(!is.matrix(g3)){
      cn<-names(g3)
      g3<-matrix(g3,nrow=1,byrow=T)
      colnames(g3)<-cn
    }
    vclia.vec<-g3 %*% rlia.vec 
    vclia.cov<-g3 %*% rlia.cov %*% t(g3)
    rownames(vclia.cov)<-rownames(vclia.vec)<-colnames(vclia.cov)<-rownames(vc.vec2)
    outrprt<-data.table(param=c(paste0(rownames(vclia.cov),"_lia"), 
                                paste0(rownames(vc.vec2)), 
                                paste0(rownames(rlia.vec),"_lia"), 
                                paste0(rownames(robs.vec))),
                        est=c(vclia.vec[,1], vc.vec2[,1], rlia.vec[,1], robs.vec[,1]),
                        se=c(sqrt(diag(vclia.cov)), sqrt(diag(vc.cov2)), sqrt(diag(rlia.cov)), sqrt(diag(robs.cov)))
               )
  returnl<-list(outrprt=outrprt, 
                vclia.cov=vclia.cov,
                vcobs.cov=vc.cov2,
                rlia.cov=rlia.cov,
                robs.cov=robs.cov)

  } else {

  outrprt<-data.table(param=c(paste0(rownames(vc.vec2)), 
                              paste0(rownames(robs.vec))),
                      est=c(vc.vec2[,1], robs.vec[,1]),
                      se=c(sqrt(diag(vc.cov2)), sqrt(diag(robs.cov)))
             )
  returnl<-list(outrprt=outrprt, 
                vcobs.cov=vc.cov2,
                robs.cov=robs.cov)

  } 
  return(returnl)
}
 


to_lia2<-function (vc, varvc,K) 
{
    t <- qnorm(1 - K)
    z <- dnorm(t)
    i <- z/K
    Kr <- K + vc/K
    tr <- qnorm(1 - Kr)
    rvc <- (t - tr * sqrt(1 - (t^2 - tr^2) * (1 - t/i)))/(i + 
        tr^2 * (i - t))
   
 
    zr<-dnorm(tr)

    vtr<-varvc/K^2/zr^2

    gprim2_num1<- -1*sqrt(1-(t^2 - tr^2)*(1 - t/i)) - 
               0.5*tr*(1-(t^2 - tr^2)*(1-t/i))^(-.5)*2*tr*(1-t/i)
    gprim2_num2<- -(t - tr*sqrt(1-(t^2 - tr^2)*(1-t/i)))*2*tr*(i-t)
    
    gprime2_denom<-i + tr^2*(i-t)
    gprime<-(gprim2_num1/gprime2_denom + gprim2_num2 / gprime2_denom^2)
    gprime.wrt.robs<- gprime * (1/(dnorm(qnorm(1-Kr))) * -1 * (1/K))
    gprime2<-gprime^2

    
    vrlia<-vtr*gprime2 
    return(c(rlia=rvc,vrlia=vrlia,gprime=gprime,gprime.wrt.robs=gprime.wrt.robs))
}

