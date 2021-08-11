Example SMILE model Analysis
================
Dan McGuire
2021-08-11

# SMILE

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
#>  1: -86.64274 32.53492 01001  empcty  699 Alabama   Autauga 1924  699     <XY>
#>  2: -87.72257 30.72748 01003  empcty 1239 Alabama   Baldwin 1826 1467     <XY>
#>  3: -85.39321 31.86958 01005  empcty  431 Alabama   Barbour 2073  431     <XY>
#>  4: -87.12648 32.99863 01007  empcty  260 Alabama      Bibb 1876  260     <XY>
#>  5: -86.56738 33.98087 01009  empcty  420 Alabama    Blount 1930  420     <XY>
#>  6: -85.71568 32.10053 01011  empcty   88 Alabama   Bullock 2028   88     <XY>
#>  7: -86.68028 31.75243 01013  empcty  132 Alabama    Butler 1919  132     <XY>
#>  8: -85.82604 33.77143 01015  empcty  817 Alabama   Calhoun 2012  926     <XY>
#>  9: -85.39203 32.91435 01017  empcty  180 Alabama  Chambers 2074  180     <XY>
#> 10: -85.60379 34.17592 01019  empcty  252 Alabama  Cherokee 2043  252     <XY>
#> 11: -86.71880 32.84786 01021  empcty  489 Alabama   Chilton 1916  489     <XY>
#> 12: -88.26318 32.01977 01023  empcty  248 Alabama   Choctaw 1761  248     <XY>
#> 13: -87.83081 31.67669 01025  empcty  304 Alabama    Clarke 1811  304     <XY>
#> 14: -85.86058 33.26902 01027  empcty   68 Alabama      Clay 2008   68     <XY>
#> 15: -85.51881 33.67452 01029  empcty  212 Alabama  Cleburne 2052  212     <XY>
#> 16: -85.98815 31.40265 01031  empcty  314 Alabama    Coffee 1997  314     <XY>
#> 17: -87.80493 34.70047 01033  empcty  604 Alabama   Colbert 1815  604     <XY>
#> 18: -86.99367 31.42923 01035  empcty   80 Alabama   Conecuh 1890   80     <XY>
#> 19: -86.24765 32.93624 01037  empcty  116 Alabama     Coosa 1963  116     <XY>
#> 20: -86.45127 31.24849 01039  empcty  315 Alabama Covington 1942  315     <XY>
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

Calculate spatial neighborhood matrix

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


rs<-rowSums(wmatd)
print(mean(1/rs))
#> [1] 0.2793672
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
#>   2.000   4.000   6.000   6.446   8.000  36.000
summary(apply(sigma_s, 2, function(x) sum(x!=0)))   
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     2.0   176.0  1074.0   744.8  1074.0  1074.0
chsigma_s<-t(chol(sigma_s))

ustar<-rnorm(nrow(wmatd), sd=1)*sqrt(sets2)
var(ustar)
#> [1] 0.08990551
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
#>    set emprel ttype pid empreltype        u_a      u_sib      u_par u_fam
#> 1:   1      1  quad   1     parent -0.4547499 -0.7180906 0.05526169     0
#> 2:   1      2  quad   2     parent -0.1465048 -0.2057012 0.05526169     0
#> 3:   1      3  quad   3        sib -0.6589747  0.4520573 0.27486008     0
#> 4:   1      3  quad   4        sib -0.1542205  0.4520573 0.20155806     0
#> 5:   2      1  quad   5     parent  0.9949360  0.3110330 0.12758908     0
#> 6:   2      2  quad   6     parent  0.4148473  0.3524012 0.12758908     0
#> 7:   2      3  quad   7        sib  1.0403881  0.2292141 0.36512113     0
#> 8:   2      3  quad   8        sib  1.4617490  0.2292141 0.24731909     0
#>              eps fixed   kid  gid   rn         u_s          y
#> 1:  3.458450e-01     0 06059  142   91 -0.14677277 -0.9185066
#> 2: -8.266989e-05     0 06059  142   91 -0.14677277 -0.4437997
#> 3: -2.050654e-01     0 06059  142   91 -0.14677277 -0.2838955
#> 4:  2.260212e-01     0 06059  142   91 -0.14677277  0.5786433
#> 5: -4.399024e-01     0 47093 2298 1117  0.06761274  1.0612683
#> 6:  4.957666e-01     0 47093 2298 1117  0.06761274  1.4582169
#> 7: -4.111011e-02     0 47093 2298 1117  0.06761274  1.6612260
#> 8:  2.811920e-01     0 47093 2298 1117  0.06761274  2.2870869

ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps",names(ss$d),value=T)]
#>          u_a     u_sib      u_par u_fam       eps        u_s
#> 1: 0.4001115 0.2031868 0.04983674     0 0.3232322 0.02684367


ss$d[,lapply(.SD, var),.SDcols=c(grep("^[yu]", names(ss$d),value=T))] ## Actual variance of the outcome and the random effects
#>          u_a     u_sib      u_par u_fam        u_s     y
#> 1: 0.4001115 0.2031868 0.04983674     0 0.02684367 1.011
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
#> -0.05088 -0.05088 -0.05088 -0.05088  0.94912 
#> 
#> Coefficients:
#>   Estimate Std. Error t value Pr(>|t|)    
#> X 0.050875   0.001099    46.3   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2197 on 39999 degrees of freedom
#> Multiple R-squared:  0.05088,    Adjusted R-squared:  0.05085 
#> F-statistic:  2144 on 1 and 39999 DF,  p-value: < 2.2e-16
#> 
#> Order of parameters:
#>  [1] "beta"           "log_sdvc_res"   "logit_rho"      "log_sdvc_s"    
....

....

#> iter: 2  mgc: 3.697043e-14 
#> outer mgc:  46.90366 
#> iter: 1  value: -204883.7 mgc: 0.08051005 ustep: 1 
#> iter: 2  mgc: 2.220446e-14 
#> outer mgc:  46.94764 
#> iter: 1  value: -204873.9 mgc: 0.6733089 ustep: 1 
#> iter: 2  mgc: 2.442491e-14 
#> outer mgc:  389.1064 
#> iter: 1  value: -204873.9 mgc: 0.6733089 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  389.0942 
#> iter: 1  value: -204874 mgc: 0.0008267623 ustep: 1 
#> iter: 2  mgc: 3.197442e-14 
#> outer mgc:  0.04804158 
#> iter: 1  value: -204874.1 mgc: 0.0008272814 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  0.03586956 
#> outer mgc:  0.7559082 
#> [1] "calculating liability scale variance components"
```

`vc_a_lia`, `vc_c_par_lia`, `vc_c_sib_lia`, `gvc_s_lia`, represent the
G, P, C, and S variance compoenents respectively. The model estimates
are close to the true simulated variances.

``` r
model_output$vc_rpt
#> $outrprt
#>            param     Estimate    Std.Error
#>  1:     vc_a_lia 0.4645948577 0.0295112262
#>  2: vc_c_par_lia 0.0627896126 0.0360400403
#>  3: vc_c_sib_lia 0.1863731181 0.0198047245
#>  4:    gvc_s_lia 0.0254107764 0.0087622823
#>  5:         vc_a 0.0071897123 0.0005393348
#>  6:     vc_c_par 0.0007944552 0.0004896795
#>  7:     vc_c_sib 0.0043928880 0.0005008102
#>  8:        gvc_s 0.0002874069 0.0001024222
#>  9:    vc_pp_lia 0.0882003890 0.0360341514
#> 10:    vc_ps_lia 0.2577082052 0.0139037130
#> 11:    vc_ss_lia 0.4440813233 0.0169588679
#> 12:    gvc_s_lia 0.0254107764 0.0087622823
#> 13:        vc_pp 0.0010818621 0.0004913923
#> 14:        vc_ps 0.0038822631 0.0002724399
#> 15:        vc_ss 0.0082751511 0.0004745721
#> 16:        gvc_s 0.0002874069 0.0001024222
#> 
#> $vclia.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      0.0008709125  1.931494e-04 -2.818659e-04 -1.011925e-04
#> vc_c_par  0.0001931494  1.298885e-03 -5.627765e-05 -3.860101e-05
#> vc_c_sib -0.0002818659 -5.627765e-05  3.922271e-04 -8.035605e-06
#> gvc_s    -0.0001011925 -3.860101e-05 -8.035605e-06  7.677759e-05
#> 
#> $vcobs.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      2.908821e-07  4.179080e-08 -9.969333e-08 -8.987305e-09
#> vc_c_par  4.179080e-08  2.397860e-07 -1.485417e-08 -4.404922e-09
#> vc_c_sib -9.969333e-08 -1.485417e-08  2.508108e-07 -6.115553e-11
#> gvc_s    -8.987305e-09 -4.404922e-09 -6.115553e-11  1.049030e-08
#> 
#> $rlia.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 1.298460e-03 8.415502e-05 1.984176e-05 3.817658e-05
#> vc_ps 8.415502e-05 1.933132e-04 4.434466e-05 2.618135e-05
#> vc_ss 1.984176e-05 4.434466e-05 2.876032e-04 1.814575e-05
#> gvc_s 3.817658e-05 2.618135e-05 1.814575e-05 7.677759e-05
#> 
#> $robs.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 2.414664e-07 2.248713e-08 7.571811e-09 6.085383e-09
#> vc_ps 2.248713e-08 7.422351e-08 2.431570e-08 5.996652e-09
#> vc_ss 7.571811e-09 2.431570e-08 2.252187e-07 5.935497e-09
#> gvc_s 6.085383e-09 5.996652e-09 5.935497e-09 1.049030e-08
```

The full model report from TMB with all fitted parameters are found in
the `tmb_rpt` object.

``` r
model_output$tmb_rpt
#>                  param      Estimate   Std. Error            ll           ul
#>      1:     log_sdvc_a -2.4675520600 0.0375073995 -2.5410652121 -2.394038908
#>      2:     log_sdvc_s -3.3964258547 0.1781833203 -3.7456587452 -3.047192964
#>      3: log_sdvc_c_par -3.5689270048 0.3081857107 -4.1729598983 -2.964894111
#>      4: log_sdvc_c_sib -2.7138841996 0.0570023834 -2.8256068182 -2.602161581
#>      5:   log_sdvc_res -1.6674749684 0.0109850991 -1.6890053669 -1.645944570
#>     ---                                                                     
#> 101634:       vc_c_par  0.0007944552 0.0004896795 -0.0001652989  0.001754209
#> 101635:       vc_c_sib  0.0043928880 0.0005008102  0.0034113181  0.005374458
#> 101636:         vc_res  0.0356163692 0.0007824987  0.0340827000  0.037150038
#> 101637:           vc_s  0.0011217653 0.0003997597  0.0003382506  0.001905280
#> 101638:            rho  0.8138778224 0.1163323459  0.5858706141  1.041885031
#>         objective                     mess fit_time sd_time modeltype
#>      1:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#>      2:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#>      3:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#>      4:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#>      5:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#>     ---                                                              
#> 101634:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#> 101635:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#> 101636:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#> 101637:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
#> 101638:  -4134.37 relative convergence (4)   25.572   5.113  car_gmrf
```

The `indlevel_dt` object is a data-table showing all fitted effects and
the response variable.

``` r
model_output$indlevel_dt
#>        V1 y        fe          us           ua          upar         usib
#>     1:  1 0 0.0504759  0.04759013 -0.012628762 -0.0015065164 -0.004165092
#>     2:  1 0 0.0504759  0.04759013 -0.012628762 -0.0015065164 -0.004165092
#>     3:  1 0 0.0504759  0.04759013 -0.004109617 -0.0006422059 -0.007102071
#>     4:  1 0 0.0504759  0.04759013 -0.004109617 -0.0006422059 -0.007102071
#>     5:  1 0 0.0504759 -0.08025647  0.100026730 -0.0056460223 -0.015609656
#>    ---                                                                   
#> 39996:  1 0 0.0504759  0.27146140 -0.021957900 -0.0034313395  0.071859871
#> 39997:  1 0 0.0504759 -0.07964166 -0.011595236 -0.0013832245 -0.003824225
#> 39998:  1 0 0.0504759 -0.07964166 -0.011595236 -0.0013832245 -0.003824225
#> 39999:  1 0 0.0504759 -0.07964166 -0.003773290 -0.0005896484 -0.006520844
#> 40000:  1 0 0.0504759 -0.07964166 -0.003773290 -0.0005896484 -0.006520844
#>               yhat
#>     1:  0.07976567
#>     2:  0.07976567
#>     3:  0.08621214
#>     4:  0.08621214
#>     5:  0.04899048
#>    ---            
#> 39996:  0.36840793
#> 39997: -0.04596844
#> 39998: -0.04596844
#> 39999: -0.04004954
#> 40000: -0.04004954
```
