Example SMILE model Analysis
================
Dan McGuire
2024-04-28

# SMILE


[![DOI](https://zenodo.org/badge/394834563.svg)](https://zenodo.org/doi/10.5281/zenodo.11081927)


A **S**patial-**MI**xed **L**inear **E**ffects model for estimating the
contributions of genetic heritability, family-, and community-level
environment to phenotype variation.

## Installation

    library(devtools)
    install_github("dan11mcguire/smile")

## Detailed example

1.  Simulate genetic, spatial community-level, and family-level
    environment variance components.
2.  Then, try to recapture the simulated variances with the estimated
    SMILE model.

##### Part 1.

``` r
library(smile)
#> Loading required package: Rcpp
#> Loading required package: TMB
library(data.table)
library(Matrix)
library(parallel)
library(sf)
#> Linking to GEOS 3.7.2, GDAL 2.4.2, PROJ 5.2.0
library(ggplot2)
library(broom)
library(spdep)
#> Loading required package: sp
#> Loading required package: spData
#> To access larger datasets in this package, install the spDataLarge
#> package with: `install.packages('spDataLarge',
#> repos='https://nowosad.github.io/drat/', type='source')`
library(sandwich)


data(lgeo_opt)
lgeo_opt[1:20]
#>             X        Y   key opt_map    N   state   geoname  gid  nto geometry
#>  1: -86.64274 32.53492 01001  empcty  699 Alabama   Autauga 1924  699  <XY[1]>
#>  2: -87.72257 30.72748 01003  empcty 1239 Alabama   Baldwin 1826 1467  <XY[1]>
#>  3: -85.39321 31.86958 01005  empcty  431 Alabama   Barbour 2073  431  <XY[1]>
#>  4: -87.12648 32.99863 01007  empcty  260 Alabama      Bibb 1876  260  <XY[1]>
#>  5: -86.56738 33.98087 01009  empcty  420 Alabama    Blount 1930  420  <XY[1]>
#>  6: -85.71568 32.10053 01011  empcty   88 Alabama   Bullock 2028   88  <XY[1]>
#>  7: -86.68028 31.75243 01013  empcty  132 Alabama    Butler 1919  132  <XY[1]>
#>  8: -85.82604 33.77143 01015  empcty  817 Alabama   Calhoun 2012  926  <XY[1]>
#>  9: -85.39203 32.91435 01017  empcty  180 Alabama  Chambers 2074  180  <XY[1]>
#> 10: -85.60379 34.17592 01019  empcty  252 Alabama  Cherokee 2043  252  <XY[1]>
#> 11: -86.71880 32.84786 01021  empcty  489 Alabama   Chilton 1916  489  <XY[1]>
#> 12: -88.26318 32.01977 01023  empcty  248 Alabama   Choctaw 1761  248  <XY[1]>
#> 13: -87.83081 31.67669 01025  empcty  304 Alabama    Clarke 1811  304  <XY[1]>
#> 14: -85.86058 33.26902 01027  empcty   68 Alabama      Clay 2008   68  <XY[1]>
#> 15: -85.51881 33.67452 01029  empcty  212 Alabama  Cleburne 2052  212  <XY[1]>
#> 16: -85.98815 31.40265 01031  empcty  314 Alabama    Coffee 1997  314  <XY[1]>
#> 17: -87.80493 34.70047 01033  empcty  604 Alabama   Colbert 1815  604  <XY[1]>
#> 18: -86.99367 31.42923 01035  empcty   80 Alabama   Conecuh 1890   80  <XY[1]>
#> 19: -86.24765 32.93624 01037  empcty  116 Alabama     Coosa 1963  116  <XY[1]>
#> 20: -86.45127 31.24849 01039  empcty  315 Alabama Covington 1942  315  <XY[1]>
#>     rn         prop
#>  1:  1 2.087205e-04
#>  2:  2 3.699638e-04
#>  3:  3 1.286961e-04
#>  4:  4 7.763567e-05
#>  5:  5 1.254115e-04
#>  6:  6 2.627669e-05
#>  7:  7 3.941503e-05
#>  8:  8 2.439552e-04
#>  9:  9 5.374777e-05
#> 10: 10 7.524688e-05
#> 11: 11 1.460148e-04
#> 12: 12 7.405249e-05
#> 13: 13 9.077402e-05
#> 14: 14 2.030471e-05
#> 15: 15 6.330293e-05
#> 16: 16 9.376000e-05
#> 17: 17 1.803536e-04
#> 18: 18 2.388790e-05
#> 19: 19 3.463745e-05
#> 20: 20 9.405860e-05
seed<-808
n_quad=10*10**3
```

Calculate spatial neighborhood matrix First take a weighted random
sample of locations to assign a community to each of the 10,000
families.

``` r

sampleloc<-as.character(sample(x=lgeo_opt[,key], 
                   size=n_quad, 
                   prob=lgeo_opt[,prop],replace=TRUE))

ngeo_opt<-lgeo_opt[key %in% sampleloc,head(.SD,1),by=key]
ngeo_opt<-as.data.table(ngeo_opt)[,rn:=1:.N]
wsp<-poly2nb(as(st_as_sf(ngeo_opt), "Spatial"),snap=0.01)
whch<-sapply(1:length(wsp), function(x) ifelse(wsp[[x]][1]!=0,x,-9))
fillin<-which(whch==-9)
whch<-which(whch!=-9)
ivec<-rep(whch, 
      times=sapply(wsp[whch], length))
jvec<-unlist(lapply(1:length(whch), 
            function(x){
              wsp[[whch[x]]] 
            }))

wspm<-sparseMatrix(i=ivec, j=jvec)
isSymmetric(wspm)
#> [1] TRUE
wmatd<-fillin_nearest_neighb(as.data.table(ngeo_opt),wspm)$wnew
```

`wmatd` is a spatial neighborhood matrix for each of the locations in
the `ngeo_opt` data.table.

``` r
ngeo_opt[1:10]
#>       key         X        Y opt_map    N   state   geoname  gid  nto geometry
#>  1: 01003 -87.72257 30.72748  empcty 1239 Alabama   Baldwin 1826 1467  <XY[1]>
#>  2: 01011 -85.71568 32.10053  empcty   88 Alabama   Bullock 2028   88  <XY[1]>
#>  3: 01013 -86.68028 31.75243  empcty  132 Alabama    Butler 1919  132  <XY[1]>
#>  4: 01015 -85.82604 33.77143  empcty  817 Alabama   Calhoun 2012  926  <XY[1]>
#>  5: 01019 -85.60379 34.17592  empcty  252 Alabama  Cherokee 2043  252  <XY[1]>
#>  6: 01021 -86.71880 32.84786  empcty  489 Alabama   Chilton 1916  489  <XY[1]>
#>  7: 01029 -85.51881 33.67452  empcty  212 Alabama  Cleburne 2052  212  <XY[1]>
#>  8: 01031 -85.98815 31.40265  empcty  314 Alabama    Coffee 1997  314  <XY[1]>
#>  9: 01033 -87.80493 34.70047  empcty  604 Alabama   Colbert 1815  604  <XY[1]>
#> 10: 01039 -86.45127 31.24849  empcty  315 Alabama Covington 1942  315  <XY[1]>
#>     rn         prop
#>  1:  1 3.699638e-04
#>  2:  2 2.627669e-05
#>  3:  3 3.941503e-05
#>  4:  4 2.439552e-04
#>  5:  5 7.524688e-05
#>  6:  6 1.460148e-04
#>  7:  7 6.330293e-05
#>  8:  8 9.376000e-05
#>  9:  9 1.803536e-04
#> 10: 10 9.405860e-05
wmatd[1:10,1:10]
#> 10 x 10 sparse Matrix of class "ngCMatrix"
#>                          
#>  [1,] . . . . . . . . . .
#>  [2,] . . . . . . . . . .
#>  [3,] . . . . . . . . . |
#>  [4,] . . . . | . | . . .
#>  [5,] . . . | . . | . . .
#>  [6,] . . . . . . . . . .
#>  [7,] . . . | | . . . . .
#>  [8,] . . . . . . . . . |
#>  [9,] . . . . . . . . . .
#> [10,] . . | . . . . | . .
```

``` r
rs<-rowSums(wmatd)
print(mean(1/rs))
#> [1] 0.2753396
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
      

### Zs:  loading matrix mapping individuals to their location-specific spatial random effects  ###
ds<-data.table(kid=rep(sampleloc,each=4))
ds<-merge(ds,ngeo_opt[,.(gid,kid=key,rn)],by=c("kid"),sort=F)    
ivec<-lapply(1:ds[,max(rn)], function(x) ds[rn==x,which=TRUE]) 
jvec<-lapply(1:length(ivec), function(x) rep(x, length(ivec[[x]])))  
Zs<-sparseMatrix(i=unlist(ivec), j=unlist(jvec), x=1,giveCsparse=F)
#> Warning in sparseMatrix(i = unlist(ivec), j = unlist(jvec), x = 1, giveCsparse =
#> F): 'giveCsparse' has been deprecated; setting 'repr = "T"' for you
```

Simulate (S) spatial random effects with CAR covariance structure, and
simulate (G) genetic ,(P) parental-, and (C) child- level environment
random effects.

Let G=0.4, P=0.05, C=0.2, S=0.03

``` r

S<-0.03
rho<-0.9
sigma_inv<-(Mpm5%*%(diagI - rho*N)%*%Mpm5) ### CAR specification of covariance matrix
sigma_s<-chol2inv(chol(sigma_inv))
sets2<-S/(gowfac(Z=Zs,G=sigma_s) / nrow(Zs))
all.equal(t(chol(sigma_inv)) %*% (chol(sigma_inv)), sigma_inv)
#> [1] TRUE
summary(apply(sigma_inv, 2, function(x) sum(x!=0)))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   2.000   4.000   6.000   6.506   8.000  36.000
summary(apply(sigma_s, 2, function(x) sum(x!=0)))   
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     2.0   147.0  1095.0   766.5  1095.0  1095.0
chsigma_s<-t(chol(sigma_s))

ustar<-rnorm(nrow(wmatd), sd=1)*sqrt(sets2)
var(ustar)
#> [1] 0.0921445
ustar<-chsigma_s %*% ustar
ustar<-as.vector(Zs %*% ustar)
ds[["u_s"]]<-ustar

G=0.4
Fam=0
P=0.05
C=0.2
E<- 1 - G - Fam - P - C - S

n_trio<-0
mu<-0
print("simulating data.")
#> [1] "simulating data."


ss<-smile::sim_fam_nonspatial(
                        n_quad=n_quad
                       ,Ag=G
                       ,C_fam=Fam
                       ,C_par=P
                       ,C_sib=C
                       ,S=S
                       ,mu=0
                          )


ss$d<-cbind(ss$d, ds)
ss$d[,y:=u_a + u_sib + u_par + u_fam + eps + fixed + u_s]
ss$d[1:8]
#>    set emprel ttype pid empreltype         u_a      u_sib        u_par u_fam
#> 1:   1      1  quad   1     parent -0.64204155  0.4022357  0.008452374     0
#> 2:   1      2  quad   2     parent -0.06847872 -0.8762082  0.008452374     0
#> 3:   1      3  quad   3        sib -1.48573593  0.3123996  0.065532458     0
#> 4:   1      3  quad   4        sib -0.68617384  0.3123996 -0.022099407     0
#> 5:   2      1  quad   5     parent -0.84183545  0.1092574 -0.039303675     0
#> 6:   2      2  quad   6     parent  1.25537106  0.4995305 -0.039303675     0
#> 7:   2      3  quad   7        sib  0.58248898  0.1856743  0.048942750     0
#> 8:   2      3  quad   8        sib  0.70423850  0.1856743 -0.050728489     0
#>           eps fixed   kid  gid  rn           u_s          y
#> 1: -0.3597001     0 06059  142  97 -0.0008791048 -0.5919327
#> 2: -0.3685423     0 06059  142  97 -0.0008791048 -1.3056559
#> 3:  0.7036684     0 06059  142  97 -0.0008791048 -0.4050146
#> 4:  0.5712522     0 06059  142  97 -0.0008791048  0.1744995
#> 5:  0.9354838     0 29189 1518 698 -0.0266756550  0.1369264
#> 6:  0.1425446     0 29189 1518 698 -0.0266756550  1.8314668
#> 7: -0.3202865     0 29189 1518 698 -0.0266756550  0.4701439
#> 8:  0.8179898     0 29189 1518 698 -0.0266756550  1.6304984

ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps",names(ss$d),value=T)]
#>         u_a     u_sib      u_par u_fam       eps        u_s
#> 1: 0.407533 0.2016964 0.04980313     0 0.3168906 0.02974938


ss$d[,lapply(.SD, var),.SDcols=c(grep("^[yu]", names(ss$d),value=T))] ## Actual variance of the outcome and the random effects
#>         u_a     u_sib      u_par u_fam        u_s         y
#> 1: 0.407533 0.2016964 0.04980313     0 0.02974938 0.9917908
mlists<-list(Ltlist=ss$Ltlist, Zlist=ss$Zlist,Glist=ss$Glist)
ss$Zlist<-NULL
ss$Glist<-NULL
ss$Ltlist<-NULL


mlists$Glist$Ga[1:12,1:12]
#> 12 x 12 sparse Matrix of class "dgTMatrix"
#>                                                      
#>  [1,] 1.0 .   0.5 0.5 .   .   .   .   .   .   .   .  
#>  [2,] .   1.0 0.5 0.5 .   .   .   .   .   .   .   .  
#>  [3,] 0.5 0.5 1.0 0.5 .   .   .   .   .   .   .   .  
#>  [4,] 0.5 0.5 0.5 1.0 .   .   .   .   .   .   .   .  
#>  [5,] .   .   .   .   1.0 .   0.5 0.5 .   .   .   .  
#>  [6,] .   .   .   .   .   1.0 0.5 0.5 .   .   .   .  
#>  [7,] .   .   .   .   0.5 0.5 1.0 0.5 .   .   .   .  
#>  [8,] .   .   .   .   0.5 0.5 0.5 1.0 .   .   .   .  
#>  [9,] .   .   .   .   .   .   .   .   1.0 .   0.5 0.5
#> [10,] .   .   .   .   .   .   .   .   .   1.0 0.5 0.5
#> [11,] .   .   .   .   .   .   .   .   0.5 0.5 1.0 0.5
#> [12,] .   .   .   .   .   .   .   .   0.5 0.5 0.5 1.0
mlists$Zlist$Zc_par[1:12,1:9]
#> 12 x 9 sparse Matrix of class "dgTMatrix"
#>                        
#>  [1,] 1 . . . . . . . .
#>  [2,] 1 . . . . . . . .
#>  [3,] . 1 . . . . . . .
#>  [4,] . . 1 . . . . . .
#>  [5,] . . . 1 . . . . .
#>  [6,] . . . 1 . . . . .
#>  [7,] . . . . 1 . . . .
#>  [8,] . . . . . 1 . . .
#>  [9,] . . . . . . 1 . .
#> [10,] . . . . . . 1 . .
#> [11,] . . . . . . . 1 .
#> [12,] . . . . . . . . 1
mlists$Zlist$Zc_sib[1:12,1:9]
#> 12 x 9 sparse Matrix of class "dgTMatrix"
#>                        
#>  [1,] 1 . . . . . . . .
#>  [2,] . 1 . . . . . . .
#>  [3,] . . 1 . . . . . .
#>  [4,] . . 1 . . . . . .
#>  [5,] . . . 1 . . . . .
#>  [6,] . . . . 1 . . . .
#>  [7,] . . . . . 1 . . .
#>  [8,] . . . . . 1 . . .
#>  [9,] . . . . . . 1 . .
#> [10,] . . . . . . . 1 .
#> [11,] . . . . . . . . 1
#> [12,] . . . . . . . . 1
mlists$Zlist$Zc_fam[1:12,1:3]
#> 12 x 3 sparse Matrix of class "dgTMatrix"
#>            
#>  [1,] 1 . .
#>  [2,] 1 . .
#>  [3,] 1 . .
#>  [4,] 1 . .
#>  [5,] . 1 .
#>  [6,] . 1 .
#>  [7,] . 1 .
#>  [8,] . 1 .
#>  [9,] . . 1
#> [10,] . . 1
#> [11,] . . 1
#> [12,] . . 1


vc_struct<-"sgps"
vcl<-strsplit(vc_struct,"")[[1]]
tmplt<-paste0("fam_chol_car_gmrf_",vc_struct,"_us")
#tmplt<-paste0("fam_chol_sar_gmrf_",vc_struct,"_us")
if(vcl[1]!="s"){
  tmplt<-gsub("car","matern",tmplt)
}


K<-0.05
ss$d[,yb:=as.numeric(y > qnorm(1-K))]

 
binaryind<-ifelse(K==1,0,1) 
khat<-ss$d[,mean(yb)]
z<-dnorm(qnorm(1-khat))
yvec<-ss$d[,yb]
```

##### Part 2.

Fit the smile model using the `smile()` function like below.

``` r
model_output<-
     smile(
       X=matrix(rep(1,nrow(ss$d),ncol=1))
      ,y=yvec
      ,Zs=Zs
      ,wmatd=wmatd
      ,vc_structure="gpcs"
      ,spatial_structure="car"
      ,mlists=mlists
        )
#> 
#> Call:
#> lm(formula = dl$y ~ X - 1)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -0.04902 -0.04902 -0.04902 -0.04902  0.95098 
#> 
#> Coefficients:
#>   Estimate Std. Error t value Pr(>|t|)    
#> X  0.04902    0.00108   45.41   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2159 on 39999 degrees of freedom
#> Multiple R-squared:  0.04902,    Adjusted R-squared:  0.049 
#> F-statistic:  2062 on 1 and 39999 DF,  p-value: < 2.2e-16
#> 
#> Order of parameters:
#>  [1] "beta"           "log_sdvc_res"   "logit_rho"      "log_sdvc_s"    
....

....

#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  52.1227 
#> iter: 1  value: -439342 mgc: 0.08531226 ustep: 1 
#> iter: 2  mgc: 3.463896e-14 
#> outer mgc:  52.1551 
#> iter: 1  value: -439333.8 mgc: 0.7442273 ustep: 1 
#> iter: 2  mgc: 5.107026e-14 
#> outer mgc:  352.0399 
#> iter: 1  value: -439333.8 mgc: 0.7442273 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  352.0401 
#> iter: 1  value: -439334 mgc: 0.001293073 ustep: 1 
#> iter: 2  mgc: 3.907985e-14 
#> outer mgc:  0.0375477 
#> iter: 1  value: -439334.1 mgc: 0.001293921 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> outer mgc:  0.03776781 
#> outer mgc:  0.7791964 
#> [1] "calculating liability scale variance components"
```

`vc_a_lia`, `vc_c_par_lia`, `vc_c_sib_lia`, `gvc_s_lia`, represent the
G, P, C, and S variance compoenents respectively. The model estimates
are close to the true simulated variances.

``` r
model_output$vc_rpt
#> $outrprt
#>            param     Estimate    Std.Error
#>  1:     vc_a_lia 4.773216e-01 2.857057e-02
#>  2: vc_c_par_lia 1.849986e-08 3.260371e-05
#>  3: vc_c_sib_lia 1.202324e-01 2.119843e-02
#>  4:    gvc_s_lia 3.536662e-02 9.745047e-03
#>  5:         vc_a 7.226481e-03 5.061615e-04
#>  6:     vc_c_par 2.095046e-10 3.692258e-07
#>  7:     vc_c_sib 2.618706e-03 4.905807e-04
#>  8:        gvc_s 3.824204e-04 1.103593e-04
#>  9:    vc_pp_lia 3.536664e-02 9.745099e-03
#> 10:    vc_ps_lia 2.740274e-01 1.346810e-02
#> 11:    vc_ss_lia 3.942598e-01 1.937179e-02
#> 12:    gvc_s_lia 3.536662e-02 9.745047e-03
#> 13:        vc_pp 3.824206e-04 1.103599e-04
#> 14:        vc_ps 3.995661e-03 2.601611e-04
#> 15:        vc_ss 6.614367e-03 4.728298e-04
#> 16:        gvc_s 3.824204e-04 1.103593e-04
#> 
#> $vclia.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      8.162775e-04  1.427130e-10 -2.371067e-04 -1.176456e-04
#> vc_c_par  1.427130e-10  1.063002e-09 -4.374916e-11 -2.774337e-11
#> vc_c_sib -2.371067e-04 -4.374916e-11  4.493736e-04 -9.195274e-06
#> gvc_s    -1.176456e-04 -2.774337e-11 -9.195274e-06  9.496594e-05
#> 
#> $vcobs.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      2.561994e-07  2.578843e-14 -8.386966e-08 -8.545263e-09
#> vc_c_par  2.578843e-14  1.363277e-13 -9.600165e-15 -2.975186e-15
#> vc_c_sib -8.386966e-08 -9.600165e-15  2.406695e-07 -4.577766e-10
#> gvc_s    -8.545263e-09 -2.975186e-15 -4.577766e-10  1.217918e-08
#> 
#> $rlia.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 9.496695e-05 3.614318e-05 2.694787e-05 9.496591e-05
#> vc_ps 3.614318e-05 1.813897e-04 5.364111e-05 3.614314e-05
#> vc_ss 2.694787e-05 5.364111e-05 3.752661e-04 2.694787e-05
#> gvc_s 9.496591e-05 3.614314e-05 2.694787e-05 9.496594e-05
#> 
#> $robs.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 1.217931e-08 7.906558e-09 7.448772e-09 1.217918e-08
#> vc_ps 7.906558e-09 6.768378e-08 2.529117e-08 7.906548e-09
#> vc_ss 7.448772e-09 2.529117e-08 2.235680e-07 7.448772e-09
#> gvc_s 1.217918e-08 7.906548e-09 7.448772e-09 1.217918e-08
```

The full model report from TMB with all fitted parameters are found in
the `tmb_rpt` object.

``` r
model_output$tmb_rpt
#>                  param      Estimate   Std. Error            ll            ul
#>      1:     log_sdvc_a -2.465002e+00 3.502130e-02 -2.533642e+00 -2.396361e+00
#>      2:     log_sdvc_s -3.252116e+00 1.442907e-01 -3.534921e+00 -2.969311e+00
#>      3: log_sdvc_c_par -1.114314e+01 8.811880e+02 -1.738240e+03  1.715954e+03
#>      4: log_sdvc_c_sib -2.972537e+00 9.366853e-02 -3.156124e+00 -2.788950e+00
#>      5:   log_sdvc_res -1.656691e+00 7.922043e-03 -1.672218e+00 -1.641164e+00
#>     ---                                                                      
#> 101630:       vc_c_par  2.095046e-10 3.692258e-07 -7.234599e-07  7.238789e-07
#> 101631:       vc_c_sib  2.618706e-03 4.905807e-04  1.657186e-03  3.580227e-03
#> 101632:         vc_res  3.639288e-02 5.766119e-04  3.526274e-02  3.752302e-02
#> 101633:           vc_s  1.497090e-03 4.320323e-04  6.503222e-04  2.343858e-03
#> 101634:            rho  8.277991e-01 9.894636e-02  6.338678e-01  1.021730e+00
#>         objective                     mess fit_time sd_time modeltype
#>      1: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#>      2: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#>      3: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#>      4: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#>      5: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#>     ---                                                              
#> 101630: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#> 101631: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#> 101632: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#> 101633: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
#> 101634: -4815.588 relative convergence (4)    32.96   5.582  car_gmrf
```

The `indlevel_dt` object is a data-table showing all fitted effects and
the response variable.

``` r
model_output$indlevel_dt
#>        V1 y         fe          us           ua          upar         usib
#>     1:  1 0 0.04866724 -0.01496494 -0.012333923 -3.840237e-10 -0.002400056
#>     2:  1 0 0.04866724 -0.01496494 -0.012333923 -3.840237e-10 -0.002400056
#>     3:  1 0 0.04866724 -0.01496494 -0.004038160 -1.655637e-10 -0.004138933
#>     4:  1 0 0.04866724 -0.01496494 -0.004038160 -1.655637e-10 -0.004138933
#>     5:  1 0 0.04866724 -0.13458885 -0.020426972  4.283612e-09 -0.001546104
#>    ---                                                                    
#> 39996:  1 0 0.04866724 -0.39071304 -0.002817301 -1.155087e-10 -0.002887607
#> 39997:  1 0 0.04866724 -0.17347930  0.048617766 -9.728082e-10 -0.006079817
#> 39998:  1 0 0.04866724 -0.17347930  0.048617766 -9.728082e-10 -0.006079817
#> 39999:  1 0 0.04866724 -0.17347930 -0.017622375 -7.225137e-10  0.047395496
#> 40000:  1 1 0.04866724 -0.17347930  0.110105430  4.514300e-09  0.047395496
#>               yhat
#>     1:  0.01896832
#>     2:  0.01896832
#>     3:  0.02552521
#>     4:  0.02552521
#>     5: -0.10789468
#>    ---            
#> 39996: -0.34775071
#> 39997: -0.08227411
#> 39998: -0.08227411
#> 39999: -0.09503893
#> 40000:  0.03268888
```
