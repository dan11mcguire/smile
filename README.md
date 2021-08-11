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
#>  1: 01001 -86.64274 32.53492  empcty  699 Alabama  Autauga 1924  699     <XY>
#>  2: 01003 -87.72257 30.72748  empcty 1239 Alabama  Baldwin 1826 1467     <XY>
#>  3: 01005 -85.39321 31.86958  empcty  431 Alabama  Barbour 2073  431     <XY>
#>  4: 01007 -87.12648 32.99863  empcty  260 Alabama     Bibb 1876  260     <XY>
#>  5: 01009 -86.56738 33.98087  empcty  420 Alabama   Blount 1930  420     <XY>
#>  6: 01015 -85.82604 33.77143  empcty  817 Alabama  Calhoun 2012  926     <XY>
#>  7: 01017 -85.39203 32.91435  empcty  180 Alabama Chambers 2074  180     <XY>
#>  8: 01021 -86.71880 32.84786  empcty  489 Alabama  Chilton 1916  489     <XY>
#>  9: 01031 -85.98815 31.40265  empcty  314 Alabama   Coffee 1997  314     <XY>
#> 10: 01047 -87.10647 32.32598  empcty  268 Alabama   Dallas 1879  268     <XY>
#>     rn         prop
#>  1:  1 2.087205e-04
#>  2:  2 3.699638e-04
#>  3:  3 1.286961e-04
#>  4:  4 7.763567e-05
#>  5:  5 1.254115e-04
#>  6:  6 2.439552e-04
#>  7:  7 5.374777e-05
#>  8:  8 1.460148e-04
#>  9:  9 9.376000e-05
#> 10: 10 8.002446e-05
wmatd[1:10,1:10]
#> 10 x 10 sparse Matrix of class "ngCMatrix"
#>                          
#>  [1,] . . . . . . . | . |
#>  [2,] . . . . . . . . . .
#>  [3,] . . . . . . . . . .
#>  [4,] . . . . . . . | . .
#>  [5,] . . . . . . . . . .
#>  [6,] . . . . . . . . . .
#>  [7,] . . . . . . . . . .
#>  [8,] | . . | . . . . . |
#>  [9,] . . . . . . . . . .
#> [10,] | . . . . . . | . .
```

``` r
rs<-rowSums(wmatd)
print(mean(1/rs))
#> [1] 0.2746493
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
#>    2.00    4.00    6.00    6.42    8.00   32.00
summary(apply(sigma_s, 2, function(x) sum(x!=0)))   
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     2.0  1218.0  1218.0   935.3  1218.0  1218.0
chsigma_s<-t(chol(sigma_s))

ustar<-rnorm(nrow(wmatd), sd=1)*sqrt(sets2)
var(ustar)
#> [1] 0.09364819
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
#>    set emprel ttype pid empreltype        u_a      u_sib       u_par u_fam
#> 1:   1      1  quad   1     parent  0.2856087 -0.1759811  0.05023624     0
#> 2:   1      2  quad   2     parent  0.3521049  0.7804724  0.05023624     0
#> 3:   1      3  quad   3        sib -0.2596390 -0.5636323 -0.13923269     0
#> 4:   1      3  quad   4        sib  0.2847946 -0.5636323  0.29238055     0
#> 5:   2      1  quad   5     parent  0.8362206  0.6912616 -0.20792723     0
#> 6:   2      2  quad   6     parent  0.3829977 -0.1380916 -0.20792723     0
#> 7:   2      3  quad   7        sib  0.6736118  0.3535593 -0.17374926     0
#> 8:   2      3  quad   8        sib  0.3489379  0.3535593 -0.04854736     0
#>            eps fixed   kid  gid   rn        u_s          y
#> 1: -0.42298122     0 12057 2528  164  0.2005934 -0.0625241
#> 2:  0.07482450     0 12057 2528  164  0.2005934  1.4582314
#> 3: -0.88686273     0 12057 2528  164  0.2005934 -1.6487733
#> 4:  0.12518786     0 12057 2528  164  0.2005934  0.3393240
#> 5: -0.14872827     0 34580   61 1501 -0.1215483  1.0492784
#> 6:  1.07584283     0 34580   61 1501 -0.1215483  0.9912734
#> 7: -0.05439081     0 34580   61 1501 -0.1215483  0.6774828
#> 8: -0.16374868     0 34580   61 1501 -0.1215483  0.3686529

ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps",names(ss$d),value=T)]
#>          u_a     u_sib      u_par u_fam       eps        u_s
#> 1: 0.4013279 0.2026454 0.05032043     0 0.3186449 0.02759281


ss$d[,lapply(.SD, var),.SDcols=c(grep("^[yu]", names(ss$d),value=T))] ## Actual variance of the outcome and the random effects
#>          u_a     u_sib      u_par u_fam        u_s        y
#> 1: 0.4013279 0.2026454 0.05032043     0 0.02759281 0.996959
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
#> -0.05188 -0.05188 -0.05188 -0.05188  0.94812 
#> 
#> Coefficients:
#>   Estimate Std. Error t value Pr(>|t|)    
#> X 0.051875   0.001109   46.78   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2218 on 39999 degrees of freedom
#> Multiple R-squared:  0.05188,    Adjusted R-squared:  0.05185 
#> F-statistic:  2188 on 1 and 39999 DF,  p-value: < 2.2e-16
#> 
#> Order of parameters:
#>  [1] "beta"           "log_sdvc_res"   "logit_rho"      "log_sdvc_s"    
....

....

#> iter: 2  mgc: 3.197442e-14 
#> outer mgc:  49.5361 
#> iter: 1  value: -216167.5 mgc: 0.07881897 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  49.53198 
#> iter: 1  value: -216158.5 mgc: 0.3821762 ustep: 1 
#> iter: 2  mgc: 2.264855e-14 
#> outer mgc:  369.3994 
#> iter: 1  value: -216158.5 mgc: 0.3821762 ustep: 1 
#> iter: 2  mgc: 3.130829e-14 
#> outer mgc:  369.3748 
#> iter: 1  value: -216158.6 mgc: 0.0003169883 ustep: 1 
#> iter: 2  mgc: 4.307665e-14 
#> outer mgc:  0.03245477 
#> iter: 1  value: -216158.7 mgc: 0.0003172778 ustep: 1 
#> iter: 2  mgc: 2.176037e-14 
#> outer mgc:  0.009659882 
#> outer mgc:  0.723282 
#> [1] "calculating liability scale variance components"
```

`vc_a_lia`, `vc_c_par_lia`, `vc_c_sib_lia`, `gvc_s_lia`, represent the
G, P, C, and S variance compoenents respectively. The model estimates
are close to the true simulated variances.

``` r
model_output$vc_rpt
#> $outrprt
#>            param     Estimate    Std.Error
#>  1:     vc_a_lia 0.4472942939 2.891835e-02
#>  2: vc_c_par_lia 0.0310615839 3.998070e-02
#>  3: vc_c_sib_lia 0.1856823683 2.044139e-02
#>  4:    gvc_s_lia 0.0137619851 6.586203e-03
#>  5:         vc_a 0.0068637168 5.385110e-04
#>  6:     vc_c_par 0.0003779470 5.052097e-04
#>  7:     vc_c_sib 0.0043059296 5.079936e-04
#>  8:        gvc_s 0.0001580622 7.701079e-05
#>  9:    vc_pp_lia 0.0448235690 4.001669e-02
#> 10:    vc_ps_lia 0.2374091320 1.410631e-02
#> 11:    vc_ss_lia 0.4230915004 1.736471e-02
#> 12:    gvc_s_lia 0.0137619851 6.586203e-03
#> 13:        vc_pp 0.0005360092 5.061095e-04
#> 14:        vc_ps 0.0035899206 2.721769e-04
#> 15:        vc_ss 0.0078958502 4.773074e-04
#> 16:        gvc_s 0.0001580622 7.701079e-05
#> 
#> $vclia.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      8.362707e-04  1.731992e-04 -3.067454e-04 -5.345775e-05
#> vc_c_par  1.731992e-04  1.598457e-03 -6.362209e-05 -2.024965e-05
#> vc_c_sib -3.067454e-04 -6.362209e-05  4.178503e-04 -4.279795e-06
#> gvc_s    -5.345775e-05 -2.024965e-05 -4.279795e-06  4.337807e-05
#> 
#> $vcobs.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      2.899941e-07  3.801667e-08 -1.047540e-07 -4.348948e-09
#> vc_c_par  3.801667e-08  2.552368e-07 -1.522507e-08 -2.510351e-09
#> vc_c_sib -1.047540e-07 -1.522507e-08  2.580575e-07  2.193458e-10
#> gvc_s    -4.348948e-09 -2.510351e-09  2.193458e-10  5.930662e-09
#> 
#> $rlia.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 1.601335e-03 8.299912e-05 1.509723e-05 2.312841e-05
#> vc_ps 8.299912e-05 1.989880e-04 4.133548e-05 1.664919e-05
#> vc_ss 1.509723e-05 4.133548e-05 3.015333e-04 1.236939e-05
#> gvc_s 2.312841e-05 1.664919e-05 1.236939e-05 4.337807e-05
#> 
#> $robs.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 2.561468e-07 2.025417e-08 5.248451e-09 3.420311e-09
#> vc_ps 2.025417e-08 7.408024e-08 2.192256e-08 3.756188e-09
#> vc_ss 5.248451e-09 2.192256e-08 2.278224e-07 3.975534e-09
#> gvc_s 3.420311e-09 3.756188e-09 3.975534e-09 5.930662e-09
```

The full model report from TMB with all fitted parameters are found in
the `tmb_rpt` object.

``` r
model_output$tmb_rpt
#>                  param      Estimate   Std. Error            ll            ul
#>      1:     log_sdvc_a -2.4907530913 0.0392288190 -2.567640e+00 -2.4138660190
#>      2:     log_sdvc_s -3.9667241273 0.2436090941 -4.444189e+00 -3.4892590765
#>      3: log_sdvc_c_par -3.9403783023 0.6683604923 -5.250341e+00 -2.6304158087
#>      4: log_sdvc_c_sib -2.7238811133 0.0589876779 -2.839495e+00 -2.6082673891
#>      5:   log_sdvc_res -1.6422822191 0.0106473167 -1.663151e+00 -1.6214138619
#>     ---                                                                      
#> 101618:       vc_c_par  0.0003779470 0.0005052097 -6.122458e-04  0.0013681398
#> 101619:       vc_c_sib  0.0043059296 0.0005079936  3.310280e-03  0.0053015787
#> 101620:         vc_res  0.0374568963 0.0007976309  3.589357e-02  0.0390202241
#> 101621:           vc_s  0.0003585479 0.0001746911  1.615972e-05  0.0007009361
#> 101622:            rho  0.9563912052 0.0293073324  8.989499e-01  1.0138325212
#>         objective                     mess fit_time sd_time modeltype
#>      1: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#>      2: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#>      3: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#>      4: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#>      5: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#>     ---                                                              
#> 101618: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#> 101619: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#> 101620: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#> 101621: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
#> 101622: -3726.723 relative convergence (4)    29.15   4.966  car_gmrf
```

The `indlevel_dt` object is a data-table showing all fitted effects and
the response variable.

``` r
model_output$indlevel_dt
#>        V1 y         fe          us           ua          upar         usib
#>     1:  1 0 0.05148101  0.47622747 -0.013926692 -0.0008279503 -0.004716396
#>     2:  1 0 0.05148101  0.47622747 -0.013926692 -0.0008279503 -0.004716396
#>     3:  1 0 0.05148101  0.47622747 -0.004531625 -0.0003528909 -0.008040935
#>     4:  1 0 0.05148101  0.47622747 -0.004531625 -0.0003528909 -0.008040935
#>     5:  1 0 0.05148101 -0.06614711 -0.011562537 -0.0006873998 -0.003915754
#>    ---                                                                    
#> 39996:  1 0 0.05148101  0.17344462 -0.004102174 -0.0003194483 -0.007278915
#> 39997:  1 0 0.05148101 -0.55024405 -0.009452408 -0.0005619514 -0.003201141
#> 39998:  1 0 0.05148101 -0.55024405 -0.009452408 -0.0005619514 -0.003201141
#> 39999:  1 0 0.05148101 -0.55024405 -0.003075732 -0.0002395163 -0.005457592
#> 40000:  1 0 0.05148101 -0.55024405 -0.003075732 -0.0002395163 -0.005457592
#>               yhat
#>     1:  0.50823744
#>     2:  0.50823744
#>     3:  0.51478303
#>     4:  0.51478303
#>     5: -0.03083179
#>    ---            
#> 39996:  0.21322510
#> 39997: -0.51197854
#> 39998: -0.51197854
#> 39999: -0.50753588
#> 40000: -0.50753588
```
