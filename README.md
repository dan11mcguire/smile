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
seed<-240428
n_quad=5*10**4
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
#>       key         X        Y opt_map    N   state  geoname  gid  nto geometry
#>  1: 01001 -86.64274 32.53492  empcty  699 Alabama  Autauga 1924  699  <XY[1]>
#>  2: 01003 -87.72257 30.72748  empcty 1239 Alabama  Baldwin 1826 1467  <XY[1]>
#>  3: 01005 -85.39321 31.86958  empcty  431 Alabama  Barbour 2073  431  <XY[1]>
#>  4: 01007 -87.12648 32.99863  empcty  260 Alabama     Bibb 1876  260  <XY[1]>
#>  5: 01009 -86.56738 33.98087  empcty  420 Alabama   Blount 1930  420  <XY[1]>
#>  6: 01011 -85.71568 32.10053  empcty   88 Alabama  Bullock 2028   88  <XY[1]>
#>  7: 01013 -86.68028 31.75243  empcty  132 Alabama   Butler 1919  132  <XY[1]>
#>  8: 01015 -85.82604 33.77143  empcty  817 Alabama  Calhoun 2012  926  <XY[1]>
#>  9: 01017 -85.39203 32.91435  empcty  180 Alabama Chambers 2074  180  <XY[1]>
#> 10: 01019 -85.60379 34.17592  empcty  252 Alabama Cherokee 2043  252  <XY[1]>
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
wmatd[1:10,1:10]
#> 10 x 10 sparse Matrix of class "ngCMatrix"
#>                          
#>  [1,] . . . . . . . . . .
#>  [2,] . . . . . . . . . .
#>  [3,] . . . . . | . . . .
#>  [4,] . . . . . . . . . .
#>  [5,] . . . . . . . . . .
#>  [6,] . . | . . . . . . .
#>  [7,] . . . . . . . . . .
#>  [8,] . . . . . . . . . |
#>  [9,] . . . . . . . . . .
#> [10,] . . . . . . . | . .
```

``` r
rs<-rowSums(wmatd)
print(mean(1/rs))
#> [1] 0.1895747
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
#>   2.000   6.000   8.000   7.991   9.000  43.000
summary(apply(sigma_s, 2, function(x) sum(x!=0)))   
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>       2    2253    2253    2002    2253    2253
chsigma_s<-t(chol(sigma_s))

ustar<-rnorm(nrow(wmatd), sd=1)*sqrt(sets2)
var(ustar)
#> [1] 0.1494072
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
#>    set emprel ttype pid empreltype         u_a      u_sib       u_par u_fam
#> 1:   1      1  quad   1     parent  0.50166875  0.1620568 -0.24587356     0
#> 2:   1      2  quad   2     parent  0.71797562 -0.2160529 -0.24587356     0
#> 3:   1      3  quad   3        sib -0.08317736  0.2716949  0.10990682     0
#> 4:   1      3  quad   4        sib  1.28227941  0.2716949  0.14598332     0
#> 5:   2      1  quad   5     parent  0.35128758  0.9183841  0.26587638     0
#> 6:   2      2  quad   6     parent -1.09811046  0.6610450  0.26587638     0
#> 7:   2      3  quad   7        sib -0.19649152 -0.1438240  0.02954308     0
#> 8:   2      3  quad   8        sib -0.66115767 -0.1438240  0.17750730     0
#>            eps fixed   kid  gid   rn         u_s          y
#> 1:  0.42294199     0 36081 3162 1308 -0.11621849  0.7245755
#> 2: -0.35392075     0 36081 3162 1308 -0.11621849 -0.2140901
#> 3:  0.05691046     0 36081 3162 1308 -0.11621849  0.2391163
#> 4: -0.10553498     0 36081 3162 1308 -0.11621849  1.4782041
#> 5:  0.81495917     0 35840 2524 2429 -0.04499378  2.3055135
#> 6: -0.26226425     0 35840 2524 2429 -0.04499378 -0.4784472
#> 7: -0.84607883     0 35840 2524 2429 -0.04499378 -1.2018451
#> 8:  0.20221374     0 35840 2524 2429 -0.04499378 -0.4702544

ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps",names(ss$d),value=T)]
#>          u_a     u_sib      u_par u_fam       eps        u_s
#> 1: 0.3995635 0.2009005 0.05008297     0 0.3239125 0.02618131


ss$d[,lapply(.SD, var),.SDcols=c(grep("^[yu]", names(ss$d),value=T))] ## Actual variance of the outcome and the random effects
#>          u_a     u_sib      u_par u_fam        u_s        y
#> 1: 0.3995635 0.2009005 0.05008297     0 0.02618131 1.003025
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
#> -0.05057 -0.05057 -0.05057 -0.05057  0.94943 
#> 
#> Coefficients:
#>   Estimate Std. Error t value Pr(>|t|)    
#> X  0.05057    0.00049   103.2   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2191 on 199999 degrees of freedom
#> Multiple R-squared:  0.05057,    Adjusted R-squared:  0.05057 
#> F-statistic: 1.065e+04 on 1 and 199999 DF,  p-value: < 2.2e-16
#> 
#> Order of parameters:
#>  [1] "beta"           "log_sdvc_res"   "logit_rho"      "log_sdvc_s"    
....

....

#> iter: 2  mgc: 3.164136e-13 
#> outer mgc:  247.4443 
#> iter: 1  value: -1042481 mgc: 0.08116664 ustep: 1 
#> iter: 2  mgc: 1.181277e-13 
#> outer mgc:  247.5641 
#> iter: 1  value: -1042437 mgc: 3.621699 ustep: 1 
#> iter: 2  mgc: 3.766987e-13 
#> outer mgc:  760.4196 
#> iter: 1  value: -1042437 mgc: 3.621699 ustep: 1 
#> iter: 2  mgc: 3.08642e-13 
#> outer mgc:  760.2132 
#> iter: 1  value: -1042437 mgc: 0.0006558111 ustep: 1 
#> iter: 2  mgc: 3.43725e-13 
#> outer mgc:  0.2231608 
#> iter: 1  value: -1042437 mgc: 0.0006563234 ustep: 1 
#> iter: 2  mgc: 1.998401e-13 
#> outer mgc:  0.02336281 
#> outer mgc:  0.671621 
#> [1] "calculating liability scale variance components"
```

`vc_a_lia`, `vc_c_par_lia`, `vc_c_sib_lia`, `gvc_s_lia`, represent the
G, P, C, and S variance compoenents respectively. The model estimates
are close to the true simulated variances.

``` r
model_output$vc_rpt
#> $outrprt
#>            param     Estimate    Std.Error
#>  1:     vc_a_lia 0.4009703941 1.325774e-02
#>  2: vc_c_par_lia 0.0629779581 1.604484e-02
#>  3: vc_c_sib_lia 0.1986157174 9.432908e-03
#>  4:    gvc_s_lia 0.0252976462 3.284498e-03
#>  5:         vc_a 0.0059297453 2.344287e-04
#>  6:     vc_c_par 0.0007896177 2.165804e-04
#>  7:     vc_c_sib 0.0044262157 2.256567e-04
#>  8:        gvc_s 0.0002833920 3.802443e-05
#>  9:    vc_pp_lia 0.0882756043 1.615656e-02
#> 10:    vc_ps_lia 0.2257828433 6.600923e-03
#> 11:    vc_ss_lia 0.4243985606 7.893954e-03
#> 12:    gvc_s_lia 0.0252976462 3.284498e-03
#> 13:        vc_pp 0.0010730097 2.184326e-04
#> 14:        vc_ps 0.0032482647 1.203870e-04
#> 15:        vc_ss 0.0076744803 2.113882e-04
#> 16:        gvc_s 0.0002833920 3.802443e-05
#> 
#> $vclia.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      1.757677e-04  2.824167e-05 -6.679698e-05 -1.115768e-05
#> vc_c_par  2.824167e-05  2.574369e-04 -1.010014e-05 -3.595112e-06
#> vc_c_sib -6.679698e-05 -1.010014e-05  8.897975e-05 -1.720228e-06
#> gvc_s    -1.115768e-05 -3.595112e-06 -1.720228e-06  1.078793e-05
#> 
#> $vcobs.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      5.495683e-08  6.199954e-09 -2.069250e-08 -6.920302e-10
#> vc_c_par  6.199954e-09  4.690706e-08 -2.444331e-09 -3.200578e-10
#> vc_c_sib -2.069250e-08 -2.444331e-09  5.092093e-08 -1.825046e-11
#> gvc_s    -6.920302e-10 -3.200578e-10 -1.825046e-11  1.445857e-09
#> 
#> $rlia.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 2.610346e-04 1.573481e-05 3.914446e-06 7.192815e-06
#> vc_ps 1.573481e-05 4.357219e-05 8.453469e-06 5.209088e-06
#> vc_ss 3.914446e-06 8.453469e-06 6.231450e-05 3.488860e-06
#> gvc_s 7.192815e-06 5.209088e-06 3.488860e-06 1.078793e-05
#> 
#> $robs.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 4.771280e-08 3.879761e-09 1.417180e-09 1.125800e-09
#> vc_ps 3.879761e-09 1.449303e-08 4.128536e-09 1.099842e-09
#> vc_ss 1.417180e-09 4.128536e-09 4.468496e-08 1.081592e-09
#> gvc_s 1.125800e-09 1.099842e-09 1.081592e-09 1.445857e-09
```

The full model report from TMB with all fitted parameters are found in
the `tmb_rpt` object.

``` r
model_output$tmb_rpt
#>                  param      Estimate   Std. Error           ll           ul
#>      1:     log_sdvc_a -2.5638870102 0.0197671841 -2.602629979 -2.525144041
#>      2:     log_sdvc_s -3.2526713482 0.0670880832 -3.384161575 -3.121181121
#>      3: log_sdvc_c_par -3.5719808125 0.1371425555 -3.840775282 -3.303186343
#>      4: log_sdvc_c_sib -2.7101051560 0.0254909241 -2.760066449 -2.660143863
#>      5:   log_sdvc_res -1.6538187706 0.0047114505 -1.663053044 -1.644584497
#>     ---                                                                    
#> 502592:       vc_c_par  0.0007896177 0.0002165804  0.000365128  0.001214107
#> 502593:       vc_c_sib  0.0044262157 0.0002256567  0.003983937  0.004868495
#> 502594:         vc_res  0.0366025437 0.0003449021  0.035926548  0.037278539
#> 502595:           vc_s  0.0014954282 0.0002006508  0.001102160  0.001888697
#> 502596:            rho  0.8904948061 0.0393929948  0.813285955  0.967703657
#>         objective                     mess fit_time sd_time modeltype
#>      1: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#>      2: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#>      3: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#>      4: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#>      5: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#>     ---                                                              
#> 502592: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#> 502593: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#> 502594: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#> 502595: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
#> 502596: -21051.31 relative convergence (4)   145.44  23.025  car_gmrf
```

The `indlevel_dt` object is a data-table showing all fitted effects and
the response variable.

``` r
model_output$indlevel_dt
#>         V1 y         fe            us           ua          upar         usib
#>      1:  1 0 0.04895031 -0.2460843097 -0.008128981 -0.0011603001 -0.003252041
#>      2:  1 0 0.04895031 -0.2460843097 -0.008128981 -0.0011603001 -0.003252041
#>      3:  1 0 0.04895031 -0.2460843097 -0.002667393 -0.0005023227 -0.005631557
#>      4:  1 0 0.04895031 -0.2460843097 -0.002667393 -0.0005023227 -0.005631557
#>      5:  1 1 0.04895031  0.0001537616  0.108286814  0.0150442160  0.089294315
#>     ---                                                                      
#> 199996:  1 0 0.04895031 -0.0252876017 -0.003244945 -0.0006110871 -0.006850919
#> 199997:  1 0 0.04895031 -0.1352654626 -0.009012389 -0.0012863944 -0.003605453
#> 199998:  1 0 0.04895031 -0.1352654626 -0.009012389 -0.0012863944 -0.003605453
#> 199999:  1 0 0.04895031 -0.1352654626 -0.002957269 -0.0005569120 -0.006243560
#> 200000:  1 0 0.04895031 -0.1352654626 -0.002957269 -0.0005569120 -0.006243560
#>                yhat
#>      1: -0.20967532
#>      2: -0.20967532
#>      3: -0.20593527
#>      4: -0.20593527
#>      5:  0.26172942
#>     ---            
#> 199996:  0.01295576
#> 199997: -0.10021939
#> 199998: -0.10021939
#> 199999: -0.09607289
#> 200000: -0.09607289
```
