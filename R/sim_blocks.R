#inla.seed=tg[ind,seed]
sim_blocks<-function(n_dz,n_mz,Ag=0.2,C=0.1,S=0.2,X=NULL,inla.seed=NULL,mu=10){

  if(missing(inla.seed)) inla.seed<-sample.int(10^3,1)
  fake.locations = matrix(c(0,0,5,5, 0, 5, 5, 0), nrow = 4, byrow = T)
  #mesh.sim = inla.mesh.2d(loc = fake.locations, max.edge=c(0.4, 1))
  mesh.sim = inla.mesh.2d(loc = fake.locations, max.edge=c(1, 1))
  
  
  E = 1 - Ag - C - S
  s2Pair <- 0.5 * Ag + C
  s2M <- 0.5 * Ag
  s2MZ <- s2Pair + s2M
  s2DZ <- s2Pair
  
  
  spde = inla.spde2.pcmatern(mesh.sim, prior.range = c(.5, .5), prior.sigma = c(.5, .5))
  
  rng<-0.1
  sigma.u<-1
  Qu = inla.spde.precision(spde, theta=c(log(rng), log(sigma.u)))
  u = inla.qsample(n=1, Q=Qu, seed = inla.seed)
  u = u[ ,1]
  
  n<-n_dz + n_mz

  set.seed(inla.seed) ## added
  loc = matrix(runif(2*n), n)*5 # coordinates
  A = inla.spde.make.A(mesh=mesh.sim, loc=loc)
  
 # plot(mesh.sim)
 # points(loc,cex=0.5)
 # points(A[1,]%*%mesh.sim$loc[,1],
 #        A[1,]%*%mesh.sim$loc[,2],pch=19,col="red")
 
 #  xy<-loc[1,c(1,2)]
 
 # points(xy[1], xy[2],pch=4)
 # points(mesh.sim$loc[which(A[1,]!=0),1],mesh.sim$loc[which(A[1,]!=0),2],pch=18, col="darkgreen")

  ### which triangle does each spatial location belong to?  Each row is a spatial location, the columns
  ### represent the 3 row indices of the mesh.sim$loc matrix corresponding to the triangle.     
  trip<-data.table(matrix(unlist(lapply(1:A@Dim[1],
  		                      function(x) which(as.matrix(A[x,]!=0)))), 
  		         ncol=3,byrow=T))
  names(trip)<-tolower(names(trip))

  ### unique triangles
  utrip<-trip[,unique(.SD)]
  ### poly data.table has utrip[,.N] * 3 rows, 
  ### with id corresponding to the unique triangle each point belongs to 
  
  for(i in 1:nrow(utrip)){
    polyi<-data.table(mesh.sim$loc[unlist(utrip[i,]),], id=i)
    if(i==1){
      poly<-polyi 
    }else{
      poly<-rbindlist(list(poly,polyi)) 
    }
  }
  names(poly)[1:2]<-c("lat", "lon")
  ### pol groups points by id to make the polygon (triangle) object of length utrip[,.N] 
  pol<-st_as_sf(poly, coords=c("lat", "lon")) %>%
   group_by(id) %>%
   summarize(do_union=FALSE) %>%
   st_cast("POLYGON")
  
  ## random centroid point within triangle (unrelated to the spatial location data of the twins) 
  #txtpts<-st_centroid_within_poly(pol) %>% st_coordinates() %>% as.data.table()
  #txtpts[,id:=pol$id]
 
  ### reference for which triangle each spatial row for sydeps / twins belongs to 
  reftrip<-merge(trip, utrip[,id:=1:.N],sort=F)
  blocklist<-list()
  locdt<-data.table(loc)[,pair:=1:.N]
  setnames(locdt, old=c("V1", "V2"),new=c("lat", "lon"))
  locdt[,blockid:=reftrip[,id]]
  locdt[,u_s:=as.vector(A%*%u)] 
  vus<-locdt[,var(u_s)*(.N-1)/.N]
  scaleus<-sqrt(S/vus)
  locdt[,u_s:=scaleus*u_s]
  


  for(i in 1:nrow(pol)){
    alln<-pol[pol[i,],] 
    alln<-alln[alln$id>=i,]
    allnid<-alln$id
    blocklist[[i]]<-allnid
   # g2<-g1 + geom_sf(data=alln[alln$id!=i,],fill="blue",alpha=0.2) + 
   #          geom_sf(data=alln[alln$id==i,],fill="red",alpha=0.2)# +
   #         # geom_text(data=txtpts[id %in% allnid], aes(x=X,y=Y,label=id))
   # xyi<-loc[reftrip[id %in% blocklist[[i]],which=T],]
   # pd<-data.table(x=xyi[,1], y=xyi[,2])
   # g3<-g2 + geom_point(data=pd, aes(x=x, y=y))
   #  readline(prompt="Press [enter] to continue")
   #  print(g3) 
  }
  
  d<-gen_twin_ids(n_dz, n_mz)
  d<-merge(d, locdt, by="pair", sort=F)
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

  d[,y:=fixed + u_pair + u_m + u_s + eps]
  return(list(
    d=d,
    blocklist=blocklist,
    sp=list(loc=loc,
            mesh=mesh.sim,
            A=A, 
            u=u,
            Qu=Qu),
            s2M=s2M,
            s2Pair=s2Pair,
            s2S = S,
            s2res=E,
            seed=inla.seed
  ))
}
