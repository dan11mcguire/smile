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

``` r
library(data.table)
library(smile)
#> Loading required package: Rcpp
#> Loading required package: TMB
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
#source("mod_fit_func_sim.R")

### 
#lgeo_opt<-st_read("../pipeline_scripts/fam2_geo_opt_key.shp") %>% as.data.table()
#lgeo_opt[,rn:=1:.N]
#lgeo_opt[,prop:=N/sum(N)]
#lgeo_opt[1:10]
data(lgeo_opt)
```

Calculate spatial neighborhood matrix

``` r
seed<-808
n_quad=50*10**3
n_quad=5*10**3
#n_quad=10*10**3

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
#> [1] 0.32758
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
random effects. \\ Let G=0.4, P=0.05, C=0.2, S=0.03

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
#>   2.000   4.000   5.000   5.849   7.000  36.000
summary(apply(sigma_s, 2, function(x) sum(x!=0)))   
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     2.0    25.0   205.0   300.8   560.0   560.0
chsigma_s<-t(chol(sigma_s))

ustar<-rnorm(nrow(wmatd), sd=1)*sqrt(sets2)
var(ustar)
#> [1] 0.07064532
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
#> 1:   1      1  quad   1     parent -0.06541632  0.5868310  0.222964576     0
#> 2:   1      2  quad   2     parent -0.50730795  0.2749536  0.222964576     0
#> 3:   1      3  quad   3        sib -1.05807568 -0.2301593 -0.149106044     0
#> 4:   1      3  quad   4        sib  0.41517717 -0.2301593 -0.323246357     0
#> 5:   2      1  quad   5     parent -0.25286521 -0.2885624  0.026776329     0
#> 6:   2      2  quad   6     parent  0.15088998 -0.2890426  0.026776329     0
#> 7:   2      3  quad   7        sib  0.32982146 -0.6688153  0.286072104     0
#> 8:   2      3  quad   8        sib -0.31943207 -0.6688153 -0.001287953     0
#>            eps fixed   kid  gid  rn          u_s          y
#> 1: -0.27104256     0 22097 1353 395 -0.009263186  0.4640735
#> 2:  0.59977828     0 22097 1353 395 -0.009263186  0.5811253
#> 3:  0.50681966     0 22097 1353 395 -0.009263186 -0.9397845
#> 4: -0.52861327     0 22097 1353 395 -0.009263186 -0.6761049
#> 5: -0.22266276     0 48085  884 889  0.044714687 -0.6925993
#> 6: -0.02548822     0 48085  884 889  0.044714687 -0.0921498
#> 7:  1.20247902     0 48085  884 889  0.044714687  1.1942720
#> 8: -0.30814861     0 48085  884 889  0.044714687 -1.2529692

ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps",names(ss$d),value=T)]
#>          u_a     u_sib      u_par u_fam       eps        u_s
#> 1: 0.4022313 0.1997663 0.04939735     0 0.3211561 0.02808586


ss$d[,lapply(.SD, var),.SDcols=c(grep("^[yu]", names(ss$d),value=T))] ## Actual variance of the outcome and the random effects
#>          u_a     u_sib      u_par u_fam        u_s        y
#> 1: 0.4022313 0.1997663 0.04939735     0 0.02808586 1.007999
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
#> -0.05195 -0.05195 -0.05195 -0.05195  0.94805 
#> 
#> Coefficients:
#>   Estimate Std. Error t value Pr(>|t|)    
#> X 0.051950   0.001569    33.1   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2219 on 19999 degrees of freedom
#> Multiple R-squared:  0.05195,    Adjusted R-squared:  0.0519 
#> F-statistic:  1096 on 1 and 19999 DF,  p-value: < 2.2e-16
#> 
#> Order of parameters:
#>  [1] "beta"           "log_sdvc_res"   "logit_rho"      "log_sdvc_s"    
#>  [5] "x"              "log_sdvc_a"     "ua"             "log_sdvc_c_par"
#>  [9] "uc_par"         "log_sdvc_c_sib" "uc_sib"        
#> Not matching template order:
#>  [1] "log_sdvc_a"     "log_sdvc_s"     "log_sdvc_c_par" "log_sdvc_c_sib"
#>  [5] "log_sdvc_res"   "ua"             "uc_par"         "uc_sib"        
#>  [9] "beta"           "logit_rho"      "x"             
#> Your parameter list has been re-ordered.
#> (Disable this warning with checkParameterOrder=FALSE)
#> Constructing atomic qnorm1
#> Constructing atomic qnorm1
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>      -2.545107      -2.545107      -2.545107      -2.545107      -1.851959 
#>           beta      logit_rho 
#>       0.051950       0.000000 
#> Constructing atomic qnorm1
#> Optimizing tape... Done
#> iter: 1  value: -88356.25 mgc: 76.99326 ustep: 1 
#> iter: 2  mgc: 1.025846e-13 
#> iter: 1  mgc: 1.025846e-13 
#> Matching hessian patterns... Done
#> outer mgc:  1798.806 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.366666174   -2.562918730   -2.442059415   -2.378294552   -0.887985649 
#>           beta      logit_rho 
#>    0.040046589    0.003118651 
#> iter: 1  value: -70266.07 mgc: 26.59381 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5247416509  -2.5471394223  -2.5333460769  -2.5260687661  -1.7419436825 
#>           beta      logit_rho 
#>   0.0505914961   0.0003559232 
#> iter: 1  value: -86849.7 mgc: 6.88383 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  186.8681 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.537854957   -2.569241833   -2.606362293   -2.545031410   -1.739367467 
#>           beta      logit_rho 
#>    0.132206276    0.002805981 
#> iter: 1  value: -87719.6 mgc: 85.50044 ustep: 1 
#> iter: 2  mgc: 1.948219e-12 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5268713812  -2.5487887945  -2.5397808442  -2.5285928852  -1.7476099065 
#>           beta      logit_rho 
#>   0.0571550135   0.0005317879 
#> iter: 1  value: -87041.57 mgc: 7.104461 ustep: 1 
#> iter: 2  mgc: 1.065814e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5257732345  -2.5479383334  -2.5364629035  -2.5272913809  -1.7446882476 
#>           beta      logit_rho 
#>   0.0537706858   0.0004411073 
#> iter: 1  value: -86943.58 mgc: 3.424029 ustep: 1 
#> iter: 2  mgc: 5.018208e-14 
#> iter: 1  mgc: 5.018208e-14 
#> outer mgc:  306.2894 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5266485094  -2.5487532910  -2.5394923054  -2.5283671643  -1.7464343018 
#>           beta      logit_rho 
#>   0.0498036458   0.0005305546 
#> iter: 1  value: -87024.25 mgc: 4.34731 ustep: 1 
#> iter: 2  mgc: 5.551115e-14 
#> iter: 1  mgc: 5.551115e-14 
#> outer mgc:  310.6953 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5289886748  -2.5512061285  -2.5483261936  -2.5313157691  -1.7499571764 
#>           beta      logit_rho 
#>   0.0531985958   0.0008041783 
#> iter: 1  value: -87241.36 mgc: 3.64426 ustep: 1 
#> iter: 2  mgc: 2.697842e-14 
#> iter: 1  mgc: 2.697842e-14 
#> outer mgc:  219.397 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.530985072   -2.554087934   -2.557843826   -2.534011981   -1.748620278 
#>           beta      logit_rho 
#>    0.050006792    0.001144238 
#> iter: 1  value: -87421.77 mgc: 3.521929 ustep: 1 
#> iter: 2  mgc: 6.794565e-14 
#> iter: 1  mgc: 6.794565e-14 
#> outer mgc:  282.0266 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.532600832   -2.556977157   -2.566735905   -2.536284608   -1.744070283 
#>           beta      logit_rho 
#>    0.052497745    0.001513509 
#> iter: 1  value: -87559.5 mgc: 2.644535 ustep: 1 
#> iter: 2  mgc: 4.263256e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.534477304   -2.560298478   -2.576937569   -2.538908311   -1.738694103 
#>           beta      logit_rho 
#>    0.052520328    0.001940523 
#> iter: 1  value: -87717.55 mgc: 1.124152 ustep: 1 
#> iter: 2  mgc: 4.352074e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.539698250   -2.569539455   -2.605321859   -2.546208291   -1.723735854 
#>           beta      logit_rho 
#>    0.052583163    0.003128613 
#> iter: 1  value: -88161.7 mgc: 3.152317 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> iter: 1  mgc: 3.552714e-14 
#> outer mgc:  199.307 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.55389324    -2.58659968    -2.65522379    -2.56308700    -1.73176643 
#>           beta      logit_rho 
#>     0.04964765     0.00570596 
#> iter: 1  value: -89332.37 mgc: 3.330059 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  352.2222 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56910244    -2.60629401    -2.70067455    -2.57836830    -1.72077141 
#>           beta      logit_rho 
#>     0.06653291     0.01049285 
#> iter: 1  value: -90302.76 mgc: 16.19822 ustep: 1 
#> iter: 2  mgc: 2.133849e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.555051028   -2.588998774   -2.660845412   -2.564399623   -1.726405721 
#>           beta      logit_rho 
#>    0.061247749    0.006238685 
#> iter: 1  value: -89400.32 mgc: 11.65141 ustep: 1 
#> iter: 2  mgc: 7.904788e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.553906445   -2.586885577   -2.655908968   -2.563144878   -1.730405627 
#>           beta      logit_rho 
#>    0.052743365    0.005760396 
#> iter: 1  value: -89333.53 mgc: 3.147542 ustep: 1 
#> iter: 2  mgc: 3.463896e-14 
#> iter: 1  mgc: 3.463896e-14 
#> outer mgc:  156.4138 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.554286508   -2.587677346   -2.657682289   -2.563554236   -1.728544689 
#>           beta      logit_rho 
#>    0.050643051    0.005947243 
#> iter: 1  value: -89356.18 mgc: 2.127996 ustep: 1 
#> iter: 2  mgc: 3.064216e-14 
#> iter: 1  mgc: 3.064216e-14 
#> outer mgc:  188.4621 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.556261188   -2.590191847   -2.662887434   -2.565375085   -1.726842879 
#>           beta      logit_rho 
#>    0.052632301    0.006645281 
#> iter: 1  value: -89471.01 mgc: 1.962786 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  138.672 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.560531119   -2.595414289   -2.673840592   -2.569249648   -1.724405276 
#>           beta      logit_rho 
#>    0.050876863    0.008185707 
#> iter: 1  value: -89723.24 mgc: 1.822159 ustep: 1 
#> iter: 2  mgc: 4.13003e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56943213    -2.60578693    -2.69517827    -2.57718329    -1.72169166 
#>           beta      logit_rho 
#>     0.05077471     0.01133206 
#> iter: 1  value: -90244.78 mgc: 1.813557 ustep: 1 
#> iter: 2  mgc: 3.774758e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60311671    -2.64504067    -2.77592759    -2.60720702    -1.71142237 
#>           beta      logit_rho 
#>     0.05038813     0.02323899 
#> iter: 1  value: -92235.83 mgc: 7.290917 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  310.2447 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.65701879    -2.62411103    -2.86807973    -2.58821204    -1.66854702 
#>           beta      logit_rho 
#>     0.05818349     0.10169043 
#> iter: 1  value: -93973.22 mgc: 11.28259 ustep: 1 
#> iter: 2  mgc: 9.85878e-14 
#> iter: 1  mgc: 9.85878e-14 
#> outer mgc:  1081.918 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.5847190     -2.6293883     -2.9389885     -2.5448372     -1.6795585 
#>           beta      logit_rho 
#>      0.0410567      0.1892740 
#> iter: 1  value: -93146.97 mgc: 15.14097 ustep: 1 
#> iter: 2  mgc: 1.163514e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.62814689    -2.62642511    -2.89735659    -2.57116179    -1.67587514 
#>           beta      logit_rho 
#>     0.04178252     0.13717098 
#> iter: 1  value: -93670.22 mgc: 14.40922 ustep: 1 
#> iter: 2  mgc: 1.205702e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.64735123    -2.62511474    -2.87894647    -2.58280284    -1.67424630 
#>           beta      logit_rho 
#>     0.04210348     0.11413038 
#> iter: 1  value: -93905.53 mgc: 14.08761 ustep: 1 
#> iter: 2  mgc: 4.156675e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.65743968    -2.62431981    -2.86877998    -2.58877844    -1.67187933 
#>           beta      logit_rho 
#>     0.04720265     0.10176612 
#> iter: 1  value: -94018.73 mgc: 9.577916 ustep: 1 
#> iter: 2  mgc: 1.071365e-13 
#> iter: 1  mgc: 1.071365e-13 
#> outer mgc:  774.8201 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.65047905    -2.62581264    -2.87221602    -2.58594816    -1.67490434 
#>           beta      logit_rho 
#>     0.05252279     0.10674310 
#> iter: 1  value: -93921.63 mgc: 4.608033 ustep: 1 
#> iter: 2  mgc: 4.92939e-14 
#> iter: 1  mgc: 4.92939e-14 
#> outer mgc:  209.4053 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63670160    -2.63266120    -2.87885438    -2.58213552    -1.68312542 
#>           beta      logit_rho 
#>     0.04889541     0.11918556 
#> iter: 1  value: -93758.78 mgc: 3.34398 ustep: 1 
#> iter: 2  mgc: 2.875478e-14 
#> iter: 1  mgc: 2.875478e-14 
#> outer mgc:  490.9611 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.6134772     -2.6509368     -2.8846957     -2.5832991     -1.6854521 
#>           beta      logit_rho 
#>      0.0532004      0.1536927 
#> iter: 1  value: -93444.31 mgc: 3.442785 ustep: 1 
#> iter: 2  mgc: 4.751755e-14 
#> iter: 1  mgc: 4.751755e-14 
#> outer mgc:  245.7672 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61501301    -2.66215653    -2.88272964    -2.60067291    -1.68292743 
#>           beta      logit_rho 
#>     0.05157901     0.19468755 
#> iter: 1  value: -93651.66 mgc: 1.511387 ustep: 1 
#> iter: 2  mgc: 3.952394e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61647121    -2.67379558    -2.88034684    -2.61870601    -1.67981961 
#>           beta      logit_rho 
#>     0.05164857     0.23764442 
#> iter: 1  value: -93858.15 mgc: 1.479119 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61956677    -2.69850366    -2.87528847    -2.65698780    -1.67322212 
#>           beta      logit_rho 
#>     0.05179626     0.32883595 
#> iter: 1  value: -94301.08 mgc: 3.202904 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  122.4633 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60907768    -2.74577934    -2.96264324    -2.62337482    -1.69304681 
#>           beta      logit_rho 
#>     0.04709975     0.70706713 
#> iter: 1  value: -94971.84 mgc: 5.442731 ustep: 1 
#> iter: 2  mgc: 5.195844e-14 
#> iter: 1  mgc: 5.195844e-14 
#> outer mgc:  774.8047 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60625777    -2.77408568    -3.04891386    -2.58837998    -1.68003739 
#>           beta      logit_rho 
#>     0.05270959     1.08772765 
#> iter: 1  value: -95599.24 mgc: 8.748418 ustep: 1 
#> iter: 2  mgc: 8.482104e-14 
#> iter: 1  mgc: 8.482104e-14 
#> outer mgc:  157.5064 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.65214460    -2.89522374    -3.09149148    -2.50525704    -1.67349706 
#>           beta      logit_rho 
#>     0.04742715     1.44690257 
#> iter: 1  value: -95904.97 mgc: 5.94437 ustep: 1 
#> iter: 2  mgc: 3.019807e-14 
#> iter: 1  mgc: 3.019807e-14 
#> outer mgc:  640.4947 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.72023807    -2.81024606    -3.13396070    -2.52346108    -1.66248136 
#>           beta      logit_rho 
#>     0.05291554     1.82164940 
#> iter: 1  value: -97932.95 mgc: 5.976793 ustep: 1 
#> iter: 2  mgc: 5.595524e-14 
#> iter: 1  mgc: 5.595524e-14 
#> outer mgc:  155.7584 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.49373357    -2.99154457    -3.12038729    -2.57161041    -1.68761690 
#>           beta      logit_rho 
#>     0.05617952     2.08099689 
#> iter: 1  value: -94376.62 mgc: 14.05716 ustep: 1 
#> iter: 2  mgc: 3.064216e-14 
#> iter: 1  mgc: 3.064216e-14 
#> outer mgc:  592.6163 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56006236    -3.12445245    -3.05737219    -2.53159266    -1.69544225 
#>           beta      logit_rho 
#>     0.04043264     2.18443462 
#> iter: 1  value: -94158.1 mgc: 9.075887 ustep: 1 
#> iter: 2  mgc: 7.11653e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.51750719    -3.03693271    -3.09939938    -2.55896052    -1.69151695 
#>           beta      logit_rho 
#>     0.03795974     2.11604140 
#> iter: 1  value: -94337.3 mgc: 11.08029 ustep: 1 
#> iter: 2  mgc: 5.284662e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.49571325    -2.99215286    -3.12089309    -2.57294493    -1.68948623 
#>           beta      logit_rho 
#>     0.03692794     2.08105280 
#> iter: 1  value: -94433.83 mgc: 12.03298 ustep: 1 
#> iter: 2  mgc: 1.678657e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.49421830    -2.99169351    -3.12051113    -2.57193717    -1.68807461 
#>           beta      logit_rho 
#>     0.05146574     2.08101058 
#> iter: 1  value: -94395.05 mgc: 2.939343 ustep: 1 
#> iter: 2  mgc: 6.084022e-14 
#> iter: 1  mgc: 6.084022e-14 
#> outer mgc:  59.07071 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.49738022    -2.99277381    -3.12131307    -2.57404548    -1.69053100 
#>           beta      logit_rho 
#>     0.05218549     2.08107303 
#> iter: 1  value: -94502.78 mgc: 0.4254985 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  79.39652 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.50218080    -2.99956120    -3.11789400    -2.57127966    -1.69017782 
#>           beta      logit_rho 
#>     0.05084846     2.08014104 
#> iter: 1  value: -94499.76 mgc: 0.8638189 ustep: 1 
#> iter: 2  mgc: 4.840572e-14 
#> iter: 1  mgc: 4.840572e-14 
#> outer mgc:  93.93297 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.51093427    -3.01344532    -3.11066163    -2.56485327    -1.69006736 
#>           beta      logit_rho 
#>     0.05217349     2.07972393 
#> iter: 1  value: -94465.57 mgc: 0.7300091 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  79.89112 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.52362966    -3.03924538    -3.10035978    -2.55568840    -1.68997259 
#>           beta      logit_rho 
#>     0.05042233     2.10069401 
#> iter: 1  value: -94420.53 mgc: 1.263975 ustep: 1 
#> iter: 2  mgc: 2.442491e-14 
#> iter: 1  mgc: 2.442491e-14 
#> outer mgc:  154.4061 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.52320104    -3.04620627    -3.10315978    -2.55877133    -1.68929120 
#>           beta      logit_rho 
#>     0.05193347     2.13799050 
#> iter: 1  value: -94486.18 mgc: 0.8647523 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  49.53952 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.5221081     -3.0515256     -3.1075616     -2.5631621     -1.6883873 
#>           beta      logit_rho 
#>      0.0510666      2.1752738 
#> iter: 1  value: -94578.56 mgc: 0.5219072 ustep: 1 
#> iter: 2  mgc: 2.575717e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.5182925     -3.0627547     -3.1167723     -2.5716433     -1.6871283 
#>           beta      logit_rho 
#>      0.0510291      2.2572931 
#> iter: 1  value: -94739.77 mgc: 0.8612386 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.50403428    -3.10471487    -3.15119081    -2.60333551    -1.68242369 
#>           beta      logit_rho 
#>     0.05088898     2.56377911 
#> iter: 1  value: -95346.83 mgc: 3.252661 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  73.60651 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.57746759    -3.29898574    -3.24283398    -2.61022279    -1.66576322 
#>           beta      logit_rho 
#>     0.05042685     2.93468402 
#> iter: 1  value: -97923.02 mgc: 11.40281 ustep: 1 
#> iter: 2  mgc: 1.509903e-14 
#> iter: 1  mgc: 1.509903e-14 
#> outer mgc:  128.8888 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60388689    -3.41796055    -3.48225477    -2.63812109    -1.64577543 
#>           beta      logit_rho 
#>     0.05139178     3.14705700 
#> iter: 1  value: -102115 mgc: 25.5705 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  145.2116 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60973694    -3.52990592    -3.48333958    -2.62615235    -1.65133626 
#>           beta      logit_rho 
#>     0.05176069     3.47219523 
#> iter: 1  value: -102085.7 mgc: 3.771077 ustep: 1 
#> iter: 2  mgc: 1.776357e-14 
#> iter: 1  mgc: 1.776357e-14 
#> outer mgc:  38.96852 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.6134000     -3.6199175     -3.5605052     -2.6226146     -1.6496566 
#>           beta      logit_rho 
#>      0.0515544      3.6795251 
#> iter: 1  value: -103213.6 mgc: 6.998313 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  16.29983 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60995506    -3.68132304    -3.64698347    -2.62760494    -1.64799274 
#>           beta      logit_rho 
#>     0.05139861     3.77891222 
#> iter: 1  value: -104471.5 mgc: 7.588749 ustep: 1 
#> iter: 2  mgc: 1.776357e-14 
#> iter: 1  mgc: 1.776357e-14 
#> outer mgc:  1.311536 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61028659    -3.72242309    -3.71770785    -2.62069512    -1.64671929 
#>           beta      logit_rho 
#>     0.05134957     3.89901995 
#> iter: 1  value: -105418 mgc: 6.119271 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  18.979 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60851778    -3.72037199    -3.74203028    -2.62385316    -1.64671989 
#>           beta      logit_rho 
#>     0.05138244     3.88534274 
#> iter: 1  value: -105786.7 mgc: 1.972343 ustep: 1 
#> iter: 2  mgc: 2.930989e-14 
#> iter: 1  mgc: 2.930989e-14 
#> outer mgc:  5.806871 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60849096    -3.71360681    -3.76326986    -2.62469059    -1.64648811 
#>           beta      logit_rho 
#>     0.05139855     3.86807099 
#> iter: 1  value: -106110.8 mgc: 1.735876 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  1.839621 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61046329    -3.67854321    -3.83586666    -2.62557843    -1.64546009 
#>           beta      logit_rho 
#>     0.05142956     3.78914182 
#> iter: 1  value: -107231.9 mgc: 6.266385 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  3.130112 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61333195    -3.64929501    -3.89247807    -2.62442892    -1.64457166 
#>           beta      logit_rho 
#>     0.05142689     3.72880815 
#> iter: 1  value: -108107.4 mgc: 4.823385 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  1.711554 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61453202    -3.64396498    -3.90706192    -2.62335664    -1.64437390 
#>           beta      logit_rho 
#>     0.05141565     3.71995485 
#> iter: 1  value: -108331.2 mgc: 1.190504 ustep: 1 
#> iter: 2  mgc: 1.776357e-14 
#> iter: 1  mgc: 1.776357e-14 
#> outer mgc:  0.2728936 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61475255    -3.64428036    -3.90774272    -2.62316269    -1.64436261 
#>           beta      logit_rho 
#>     0.05141343     3.72087208 
#> iter: 1  value: -108342.6 mgc: 0.05508782 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  0.04405277 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61480964    -3.64439401    -3.90772786    -2.62313433    -1.64436010 
#>           beta      logit_rho 
#>     0.05141321     3.72105254 
#> iter: 1  value: -108343.1 mgc: 0.005194576 ustep: 1 
#> iter: 2  mgc: 3.419487e-14 
#> iter: 1  mgc: 3.419487e-14 
#> outer mgc:  0.01688612 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61481346    -3.64440028    -3.90773211    -2.62313580    -1.64435983 
#>           beta      logit_rho 
#>     0.05141326     3.72104382 
#> iter: 1  value: -108343.2 mgc: 0.0004241431 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  0.006399224 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61481334    -3.64440343    -3.90775810    -2.62313707    -1.64435956 
#>           beta      logit_rho 
#>     0.05141328     3.72104420 
#> iter: 1  value: -108343.6 mgc: 0.002091049 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61481346    -3.64440028    -3.90773211    -2.62313580    -1.64435983 
#>           beta      logit_rho 
#>     0.05141326     3.72104382 
#> iter: 1  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  0.006399224 
#> iter: 1  value: -108325.2 mgc: 0.06920589 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  2.986136 
#> iter: 1  value: -108361.2 mgc: 0.06934444 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  2.970662 
#> iter: 1  value: -108343.3 mgc: 0.03124589 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.09180077 
#> iter: 1  value: -108343.1 mgc: 0.03121466 ustep: 1 
#> iter: 2  mgc: 1.421085e-14 
#> outer mgc:  0.09114515 
#> iter: 1  value: -108328.4 mgc: 0.07948187 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.27622 
#> iter: 1  value: -108358 mgc: 0.07964099 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.2629101 
#> iter: 1  value: -108330.2 mgc: 0.06888544 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  2.970607 
#> iter: 1  value: -108356.2 mgc: 0.06902334 ustep: 1 
#> iter: 2  mgc: 2.575717e-14 
#> outer mgc:  2.955331 
#> iter: 1  value: -108338.9 mgc: 0.07948187 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  25.03011 
#> iter: 1  value: -108347.5 mgc: 0.07964099 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  25.0299 
#> iter: 1  value: -108343.1 mgc: 0.2970953 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  120.378 
#> iter: 1  value: -108343.1 mgc: 0.2970953 ustep: 1 
#> iter: 2  mgc: 1.865175e-14 
#> outer mgc:  120.3839 
#> iter: 1  value: -108343.2 mgc: 0.0005140475 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.03773046 
#> iter: 1  value: -108343.2 mgc: 0.0005145375 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.03712745 
#> outer mgc:  0.6132293 
#> [1] "calculating liability scale variance components"
```

`vc_a_lia`, `vc_c_par_lia`, `vc_c_sib_lia`, `gvc_s_lia`, represent the
G, P, C, and S variance compoenents respectively. The model estimates
are close to the true simulated variances.

``` r
model_output$vc_rpt
#> $outrprt
#>            param     Estimate    Std.Error
#>  1:     vc_a_lia 0.3397640449 0.0526277779
#>  2: vc_c_par_lia 0.0300942362 0.0504178327
#>  3: vc_c_sib_lia 0.2242955814 0.0296958570
#>  4:    gvc_s_lia 0.0528618148 0.0359010714
#>  5:         vc_a 0.0053555230 0.0007838456
#>  6:     vc_c_par 0.0004034475 0.0006945118
#>  7:     vc_c_sib 0.0052671198 0.0007322112
#>  8:        gvc_s 0.0006401133 0.0004641103
#>  9:    vc_pp_lia 0.0829560510 0.0560301806
#> 10:    vc_ps_lia 0.2227438373 0.0274964686
#> 11:    vc_ss_lia 0.4470394187 0.0264299955
#> 12:    gvc_s_lia 0.0528618148 0.0359010714
#> 13:        vc_pp 0.0010435608 0.0007782341
#> 14:        vc_ps 0.0033178748 0.0005161710
#> 15:        vc_ss 0.0085849946 0.0007611440
#> 16:        gvc_s 0.0006401133 0.0004641103
#> 
#> $vclia.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      0.0027696830  5.487727e-04 -4.488029e-04 -0.0012252519
#> vc_c_par  0.0005487727  2.541958e-03 -3.297723e-05 -0.0003457318
#> vc_c_sib -0.0004488029 -3.297723e-05  8.818439e-04 -0.0002452761
#> gvc_s    -0.0012252519 -3.457318e-04 -2.452761e-04  0.0012888869
#> 
#> $vcobs.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      6.144140e-07  7.931873e-08 -2.159024e-07 -1.025693e-07
#> vc_c_par  7.931873e-08  4.823467e-07 -2.339888e-08 -4.604839e-08
#> vc_c_sib -2.159024e-07 -2.339888e-08  5.361332e-07 -3.661507e-09
#> gvc_s    -1.025693e-07 -4.604839e-08 -3.661507e-09  2.153984e-07
#> 
#> $rlia.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 0.0031393811 0.0006049155 0.0003266622 0.0009431551
#> vc_ps 0.0006049155 0.0007560558 0.0002863783 0.0006762610
#> vc_ss 0.0003266622 0.0002863783 0.0006985447 0.0004309849
#> gvc_s 0.0009431551 0.0006762610 0.0004309849 0.0012888869
#> 
#> $robs.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 6.056483e-07 1.577247e-07 1.306643e-07 1.693500e-07
#> vc_ps 1.577247e-07 2.664325e-07 1.548198e-07 1.641137e-07
#> vc_ss 1.306643e-07 1.548198e-07 5.793403e-07 1.604522e-07
#> gvc_s 1.693500e-07 1.641137e-07 1.604522e-07 2.153984e-07
```

The full model report from TMB with all fitted parameters are found in
the `tmb_rpt` object.

``` r
model_output$tmb_rpt
#>                 param      Estimate   Std. Error            ll           ul
#>     1:     log_sdvc_a -2.6148134606 0.0731810540 -2.7582456909 -2.471381230
#>     2:     log_sdvc_s -3.6444002830 0.3625220897 -4.3549305225 -2.933870043
#>     3: log_sdvc_c_par -3.9077321119 0.8607214743 -5.5947152023 -2.220749021
#>     4: log_sdvc_c_sib -2.6231358009 0.0695077372 -2.7593684624 -2.486903139
#>     5:   log_sdvc_res -1.6443598274 0.0150180318 -1.6737946288 -1.614925026
#>    ---                                                                     
#> 51230:       vc_c_par  0.0004034475 0.0006945118 -0.0009577707  0.001764666
#> 51231:       vc_c_sib  0.0052671198 0.0007322112  0.0038320123  0.006702227
#> 51232:         vc_res  0.0373015777 0.0011203926  0.0351056487  0.039497507
#> 51233:           vc_s  0.0006831470 0.0004953117 -0.0002876462  0.001653940
#> 51234:            rho  0.9763635229 0.0192996032  0.9385369956  1.014190050
#>        objective                     mess fit_time sd_time modeltype
#>     1: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#>     2: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#>     3: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#>     4: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#>     5: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#>    ---                                                              
#> 51230: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#> 51231: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#> 51232: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#> 51233: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
#> 51234: -1942.728 relative convergence (4)   10.206   2.073  car_gmrf
```

The `indlevel_dt` object is a data-table showing all fitted effects and
the response variable.

``` r
model_output$indlevel_dt
#>        V1 y         fe         us           ua          upar         usib
#>     1:  1 0 0.05141326  1.1774930 -0.015295343 -0.0012443991 -0.008122989
#>     2:  1 0 0.05141326  1.1774930 -0.015295343 -0.0012443991 -0.008122989
#>     3:  1 0 0.05141326  1.1774930 -0.004975216 -0.0005300442 -0.013839750
#>     4:  1 0 0.05141326  1.1774930 -0.004975216 -0.0005300442 -0.013839750
#>     5:  1 0 0.05141326 -0.3596792 -0.007818426 -0.0006360918 -0.004152178
#>    ---                                                                   
#> 19996:  1 0 0.05141326 -0.2232779 -0.002758960 -0.0002939311 -0.007674706
#> 19997:  1 0 0.05141326 -0.1330108 -0.008920958 -0.0007257916 -0.004737706
#> 19998:  1 0 0.05141326 -0.1330108 -0.008920958 -0.0007257916 -0.004737706
#> 19999:  1 0 0.05141326 -0.1330108 -0.002901778 -0.0003091465 -0.008071988
#> 20000:  1 0 0.05141326 -0.1330108 -0.002901778 -0.0003091465 -0.008071988
#>               yhat
#>     1:  1.20424354
#>     2:  1.20424354
#>     3:  1.20956126
#>     4:  1.20956126
#>     5: -0.32087260
#>    ---            
#> 19996: -0.18259222
#> 19997: -0.09598198
#> 19998: -0.09598198
#> 19999: -0.09288044
#> 20000: -0.09288044
```
