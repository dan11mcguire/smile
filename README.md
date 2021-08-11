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
#> [1] 0.3402646
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
#>   2.000   3.000   5.000   5.665   7.000  27.000
summary(apply(sigma_s, 2, function(x) sum(x!=0)))   
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     2.0    28.0   190.0   184.3   315.0   315.0
chsigma_s<-t(chol(sigma_s))

ustar<-rnorm(nrow(wmatd), sd=1)*sqrt(sets2)
var(ustar)
#> [1] 0.06281245
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
#>    set emprel ttype pid empreltype        u_a       u_sib       u_par u_fam
#> 1:   1      1  quad   1     parent -1.7457738  0.60513039 -0.18589922     0
#> 2:   1      2  quad   2     parent  0.5603879  0.05521889 -0.18589922     0
#> 3:   1      3  quad   3        sib -0.7609064 -0.08282689 -0.35510582     0
#> 4:   1      3  quad   4        sib -0.6953841 -0.08282689 -0.17000472     0
#> 5:   2      1  quad   5     parent  0.2177877  0.51326238 -0.10783901     0
#> 6:   2      2  quad   6     parent -0.1536963 -0.86170493 -0.10783901     0
#> 7:   2      3  quad   7        sib  0.4404679 -0.68368034  0.32105597     0
#> 8:   2      3  quad   8        sib  0.2863435 -0.68368034  0.06470661     0
#>            eps fixed   kid  gid  rn          u_s          y
#> 1:  0.08461294     0 17097 1794 256  0.003037576 -1.2388921
#> 2:  0.02878910     0 17097 1794 256  0.003037576  0.4615342
#> 3: -0.76316042     0 17097 1794 256  0.003037576 -1.9589620
#> 4:  0.26704486     0 17097 1794 256  0.003037576 -0.6781333
#> 5: -0.03468727     0 25005 3223 412 -0.110243714  0.4782800
#> 6: -0.12342849     0 25005 3223 412 -0.110243714 -1.3569125
#> 7:  0.84360682     0 25005 3223 412 -0.110243714  0.8112066
#> 8: -0.28679206     0 25005 3223 412 -0.110243714 -0.7296660

ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps",names(ss$d),value=T)]
#>          u_a     u_sib      u_par u_fam       eps        u_s
#> 1: 0.4058005 0.1992955 0.04990978     0 0.3223414 0.03073576


ss$d[,lapply(.SD, var),.SDcols=c(grep("^[yu]", names(ss$d),value=T))] ## Actual variance of the outcome and the random effects
#>          u_a     u_sib      u_par u_fam        u_s        y
#> 1: 0.4058005 0.1992955 0.04990978     0 0.03073576 1.014808
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
#> -0.04815 -0.04815 -0.04815 -0.04815  0.95185 
#> 
#> Coefficients:
#>   Estimate Std. Error t value Pr(>|t|)    
#> X 0.048150   0.001514   31.81   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2141 on 19999 degrees of freedom
#> Multiple R-squared:  0.04815,    Adjusted R-squared:  0.0481 
#> F-statistic:  1012 on 1 and 19999 DF,  p-value: < 2.2e-16
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
#>      -2.581087      -2.581087      -2.581087      -2.581087      -1.887940 
#>           beta      logit_rho 
#>       0.048150       0.000000 
#> Constructing atomic qnorm1
#> Optimizing tape... Done
#> iter: 1  value: -90721.54 mgc: 83.06957 ustep: 1 
#> iter: 2  mgc: 3.947953e-13 
#> iter: 1  mgc: 3.947953e-13 
#> Matching hessian patterns... Done
#> outer mgc:  2007.034 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.3987990776  -2.6156219394  -2.4860489991  -2.4193236973  -0.9281179770 
#>           beta      logit_rho 
#>   0.1436535619  -0.0008032584 
#> iter: 1  value: -72589.18 mgc: 30.33957 ustep: 1 
#> iter: 2  mgc: 1.163514e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.558816e+00  -2.585306e+00  -2.569475e+00  -2.561323e+00  -1.770672e+00 
#>           beta      logit_rho 
#>   5.981825e-02  -9.813892e-05 
#> iter: 1  value: -89132.92 mgc: 10.93085 ustep: 1 
#> iter: 2  mgc: 1.973977e-13 
#> iter: 1  mgc: 1.973977e-13 
#> outer mgc:  1750.375 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5605790842  -2.5907793351  -2.5824375428  -2.5651746675  -1.7703064737 
#>           beta      logit_rho 
#>  -0.0614703639  -0.0001502095 
#> iter: 1  value: -88367.03 mgc: 114.0406 ustep: 1 
#> iter: 2  mgc: 7.140954e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5591314765  -2.5858218693  -2.5708337405  -2.5618301544  -1.7713791555 
#>           beta      logit_rho 
#>   0.0477226922  -0.0001026781 
#> iter: 1  value: -89176.63 mgc: 11.45208 ustep: 1 
#> iter: 2  mgc: 1.360023e-13 
#> iter: 1  mgc: 1.360023e-13 
#> outer mgc:  244.6957 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5612144720  -2.5900440538  -2.5804886076  -2.5653048708  -1.7749656855 
#>           beta      logit_rho 
#>   0.0507123890  -0.0002144348 
#> iter: 1  value: -89406.95 mgc: 2.885287 ustep: 1 
#> iter: 2  mgc: 5.406786e-14 
#> iter: 1  mgc: 5.406786e-14 
#> outer mgc:  250.6786 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5637153237  -2.6002053078  -2.6003600466  -2.5710662445  -1.7681153548 
#>           beta      logit_rho 
#>   0.0471797803  -0.0005306465 
#> iter: 1  value: -89727.86 mgc: 3.179967 ustep: 1 
#> iter: 2  mgc: 5.961898e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.566775749   -2.613297353   -2.625695890   -2.578312173   -1.758325820 
#>           beta      logit_rho 
#>    0.047662539   -0.000944995 
#> iter: 1  value: -90130.05 mgc: 2.728109 ustep: 1 
#> iter: 2  mgc: 2.4869e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.577575116   -2.659495439   -2.715098841   -2.603880984   -1.723781350 
#>           beta      logit_rho 
#>    0.049366056   -0.002407112 
#> iter: 1  value: -91581.91 mgc: 10.03395 ustep: 1 
#> iter: 2  mgc: 2.253753e-14 
#> iter: 1  mgc: 2.253753e-14 
#> outer mgc:  450.6305 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.609152839   -2.726229462   -2.836505247   -2.652885904   -1.790779981 
#>           beta      logit_rho 
#>    0.033842949   -0.004478451 
#> iter: 1  value: -94748.69 mgc: 12.59638 ustep: 1 
#> iter: 2  mgc: 2.021994e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.586072999   -2.665752971   -2.729399548   -2.613510147   -1.765327587 
#>           beta      logit_rho 
#>    0.045940024   -0.002549579 
#> iter: 1  value: -92315.69 mgc: 2.827984 ustep: 1 
#> iter: 2  mgc: 3.119727e-14 
#> iter: 1  mgc: 3.119727e-14 
#> outer mgc:  593.7781 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.588603549   -2.685644204   -2.755740561   -2.621095450   -1.741988178 
#>           beta      logit_rho 
#>    0.067105927   -0.003074405 
#> iter: 1  value: -92629.28 mgc: 17.21171 ustep: 1 
#> iter: 2  mgc: 2.260414e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.585304928   -2.666815081   -2.730263845   -2.613102152   -1.757946165 
#>           beta      logit_rho 
#>    0.056250984   -0.002567339 
#> iter: 1  value: -92258.36 mgc: 8.770278 ustep: 1 
#> iter: 2  mgc: 3.108624e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.585720844   -2.666239941   -2.729795822   -2.613323084   -1.761943258 
#>           beta      logit_rho 
#>    0.050667524   -0.002557722 
#> iter: 1  value: -92291.89 mgc: 4.055737 ustep: 1 
#> iter: 2  mgc: 3.241851e-14 
#> iter: 1  mgc: 3.241851e-14 
#> outer mgc:  358.6964 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.585296495   -2.667020576   -2.730536137   -2.613150439   -1.757366689 
#>           beta      logit_rho 
#>    0.047197126   -0.002572014 
#> iter: 1  value: -92262.95 mgc: 2.950237 ustep: 1 
#> iter: 2  mgc: 2.964295e-14 
#> iter: 1  mgc: 2.964295e-14 
#> outer mgc:  360.9547 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.586137707   -2.672098981   -2.737969110   -2.615337285   -1.751977466 
#>           beta      logit_rho 
#>    0.051869480   -0.002716574 
#> iter: 1  value: -92369.2 mgc: 3.936536 ustep: 1 
#> iter: 2  mgc: 4.307665e-14 
#> iter: 1  mgc: 4.307665e-14 
#> outer mgc:  493.6064 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.589053242   -2.683935607   -2.755625915   -2.621258864   -1.745743958 
#>           beta      logit_rho 
#>    0.047958590   -0.003053221 
#> iter: 1  value: -92691.28 mgc: 3.08793 ustep: 1 
#> iter: 2  mgc: 6.927792e-14 
#> iter: 1  mgc: 6.927792e-14 
#> outer mgc:  221.6474 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.592185786   -2.696294306   -2.773148486   -2.627353165   -1.739726731 
#>           beta      logit_rho 
#>    0.050538248   -0.003378815 
#> iter: 1  value: -93020.16 mgc: 2.135097 ustep: 1 
#> iter: 2  mgc: 7.61613e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.596638718   -2.713029580   -2.797000892   -2.635791675   -1.733198167 
#>           beta      logit_rho 
#>    0.050926853   -0.003827459 
#> iter: 1  value: -93485.97 mgc: 2.475908 ustep: 1 
#> iter: 2  mgc: 3.641532e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.608950952   -2.759302154   -2.862952141   -2.659123924   -1.715146867 
#>           beta      logit_rho 
#>    0.052001336   -0.005067947 
#> iter: 1  value: -94788.59 mgc: 7.068902 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  571.4007 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.6006785197  -2.8996769478  -2.8652898354  -2.6529054509  -1.7188533274 
#>           beta      logit_rho 
#>   0.0406829774   0.0000122619 
#> iter: 1  value: -94582.21 mgc: 6.119444 ustep: 1 
#> iter: 2  mgc: 7.66609e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.605586010   -2.821114381   -2.864605589   -2.656728568   -1.717705188 
#>           beta      logit_rho 
#>    0.041866191   -0.002847119 
#> iter: 1  value: -94719.05 mgc: 6.348955 ustep: 1 
#> iter: 2  mgc: 4.135581e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.607939773   -2.783433690   -2.864277407   -2.658562237   -1.717154509 
#>           beta      logit_rho 
#>    0.042433692   -0.004218554 
#> iter: 1  value: -94784.97 mgc: 6.479717 ustep: 1 
#> iter: 2  mgc: 3.885781e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.609310331   -2.761492822   -2.864086311   -2.659629953   -1.716833858 
#>           beta      logit_rho 
#>    0.042764138   -0.005017119 
#> iter: 1  value: -94823.45 mgc: 6.56303 ustep: 1 
#> iter: 2  mgc: 4.03011e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60910635    -2.75954785    -2.86333797    -2.65932277    -1.71571992 
#>           beta      logit_rho 
#>     0.04883055    -0.00506895 
#> iter: 1  value: -94803.37 mgc: 2.256405 ustep: 1 
#> iter: 2  mgc: 2.331468e-14 
#> iter: 1  mgc: 2.331468e-14 
#> outer mgc:  89.70804 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60972322    -2.76061464    -2.86494087    -2.66012805    -1.71782852 
#>           beta      logit_rho 
#>     0.05004407    -0.00507683 
#> iter: 1  value: -94858.81 mgc: 0.8787044 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  187.2947 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.609359131   -2.766867399   -2.865186387   -2.659798478   -1.718351942 
#>           beta      logit_rho 
#>    0.048343988   -0.004835656 
#> iter: 1  value: -94853.25 mgc: 1.162798 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  147.4545 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.609000929   -2.773271130   -2.865477579   -2.659482719   -1.718673439 
#>           beta      logit_rho 
#>    0.049388844   -0.004511259 
#> iter: 1  value: -94847.43 mgc: 0.7798896 ustep: 1 
#> iter: 2  mgc: 5.140333e-14 
#> iter: 1  mgc: 5.140333e-14 
#> outer mgc:  61.81888 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.608318933   -2.786178350   -2.866668727   -2.659277788   -1.719110424 
#>           beta      logit_rho 
#>    0.048390688   -0.003631459 
#> iter: 1  value: -94848.97 mgc: 0.6100183 ustep: 1 
#> iter: 2  mgc: 2.14273e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.606310343   -2.812438911   -2.867189408   -2.657993979   -1.719318495 
#>           beta      logit_rho 
#>    0.048257863   -0.001739287 
#> iter: 1  value: -94799.91 mgc: 0.4139964 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.597964117   -2.921558492   -2.869352974   -2.652659415   -1.720183085 
#>           beta      logit_rho 
#>    0.047705937    0.006123189 
#> iter: 1  value: -94597.58 mgc: 1.636943 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  271.9188 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59704783    -2.90417598    -2.98077301    -2.73617207    -1.69697335 
#>           beta      logit_rho 
#>     0.04985648     0.05146230 
#> iter: 1  value: -97093.41 mgc: 11.86105 ustep: 1 
#> iter: 2  mgc: 2.398082e-14 
#> iter: 1  mgc: 2.398082e-14 
#> outer mgc:  193.6 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61699684    -2.88346695    -3.07324178    -2.62720128    -1.70238866 
#>           beta      logit_rho 
#>     0.04585491     0.08298151 
#> iter: 1  value: -97347.07 mgc: 8.003251 ustep: 1 
#> iter: 2  mgc: 3.042011e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60705989    -2.89441235    -3.02624377    -2.68348658    -1.70116419 
#>           beta      logit_rho 
#>     0.04561368     0.06677939 
#> iter: 1  value: -97233.76 mgc: 4.324588 ustep: 1 
#> iter: 2  mgc: 4.396483e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60226079    -2.89969848    -3.00354586    -2.71066984    -1.70057283 
#>           beta      logit_rho 
#>     0.04549718     0.05895450 
#> iter: 1  value: -97184.5 mgc: 2.634212 ustep: 1 
#> iter: 2  mgc: 2.164935e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59943828    -2.90280744    -2.99019645    -2.72665722    -1.70022503 
#>           beta      logit_rho 
#>     0.04542866     0.05435241 
#> iter: 1  value: -97157.16 mgc: 2.650325 ustep: 1 
#> iter: 2  mgc: 2.242651e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59762811    -2.90470717    -2.98177521    -2.73661077    -1.69978184 
#>           beta      logit_rho 
#>     0.04571866     0.05147596 
#> iter: 1  value: -97137.33 mgc: 2.461304 ustep: 1 
#> iter: 2  mgc: 3.019807e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59723305    -2.90434553    -2.98109290    -2.73631210    -1.69786980 
#>           beta      logit_rho 
#>     0.04853572     0.05146666 
#> iter: 1  value: -97107.85 mgc: 0.7829112 ustep: 1 
#> iter: 2  mgc: 3.241851e-14 
#> iter: 1  mgc: 3.241851e-14 
#> outer mgc:  109.8401 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59748882    -2.90461001    -2.98157257    -2.73650523    -1.69903152 
#>           beta      logit_rho 
#>     0.04952600     0.05147325 
#> iter: 1  value: -97127.78 mgc: 0.594086 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  121.7751 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59806516    -2.90445232    -2.98367247    -2.73453946    -1.69985518 
#>           beta      logit_rho 
#>     0.04839033     0.05205046 
#> iter: 1  value: -97145.55 mgc: 0.6772067 ustep: 1 
#> iter: 2  mgc: 2.264855e-14 
#> iter: 1  mgc: 2.264855e-14 
#> outer mgc:  125.9478 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59914105    -2.90373985    -2.98782170    -2.72985129    -1.70040926 
#>           beta      logit_rho 
#>     0.04928629     0.05337667 
#> iter: 1  value: -97163.33 mgc: 0.5340967 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  68.7328 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60196447    -2.90221245    -2.99591243    -2.72038287    -1.70103397 
#>           beta      logit_rho 
#>     0.04851907     0.05629569 
#> iter: 1  value: -97204.45 mgc: 0.7996054 ustep: 1 
#> iter: 2  mgc: 4.263256e-14 
#> iter: 1  mgc: 4.263256e-14 
#> outer mgc:  99.39016 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61365409    -2.89924178    -3.00953300    -2.70305316    -1.70029171 
#>           beta      logit_rho 
#>     0.04949152     0.06447207 
#> iter: 1  value: -97370.75 mgc: 1.235917 ustep: 1 
#> iter: 2  mgc: 2.575717e-14 
#> iter: 1  mgc: 2.575717e-14 
#> outer mgc:  110.4106 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.62598724    -2.90026145    -3.01667387    -2.70099446    -1.69840091 
#>           beta      logit_rho 
#>     0.04847914     0.08653299 
#> iter: 1  value: -97653.68 mgc: 0.9787756 ustep: 1 
#> iter: 2  mgc: 2.287059e-14 
#> iter: 1  mgc: 2.287059e-14 
#> outer mgc:  109.9356 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60739839    -2.90154609    -3.02612049    -2.70190170    -1.70051661 
#>           beta      logit_rho 
#>     0.04900898     0.10258436 
#> iter: 1  value: -97476.28 mgc: 1.445775 ustep: 1 
#> iter: 2  mgc: 2.309264e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58770207    -2.90272257    -3.03581996    -2.70273789    -1.70240171 
#>           beta      logit_rho 
#>     0.04876722     0.11951414 
#> iter: 1  value: -97282.16 mgc: 1.546785 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  47.7856 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59195496    -2.91402112    -3.04297391    -2.72471653    -1.70692967 
#>           beta      logit_rho 
#>     0.05673868     0.16621787 
#> iter: 1  value: -97762.43 mgc: 4.858099 ustep: 1 
#> iter: 2  mgc: 2.930989e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58922354    -2.90450601    -3.03815678    -2.70405945    -1.70533802 
#>           beta      logit_rho 
#>     0.05212776     0.11955005 
#> iter: 1  value: -97373.55 mgc: 2.044252 ustep: 1 
#> iter: 2  mgc: 4.041212e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58796265    -2.90302802    -3.03622019    -2.70296423    -1.70290461 
#>           beta      logit_rho 
#>     0.04934278     0.11952029 
#> iter: 1  value: -97298 mgc: 0.3489379 ustep: 1 
#> iter: 2  mgc: 2.220446e-14 
#> iter: 1  mgc: 2.220446e-14 
#> outer mgc:  76.19337 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58811932    -2.90330366    -3.03650112    -2.70329996    -1.70310851 
#>           beta      logit_rho 
#>     0.04873563     0.12002474 
#> iter: 1  value: -97310.19 mgc: 0.3634486 ustep: 1 
#> iter: 2  mgc: 3.375078e-14 
#> iter: 1  mgc: 3.375078e-14 
#> outer mgc:  54.68186 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58830170    -2.90376936    -3.03688070    -2.70406731    -1.70323060 
#>           beta      logit_rho 
#>     0.04923451     0.12163114 
#> iter: 1  value: -97329.41 mgc: 0.3031111 ustep: 1 
#> iter: 2  mgc: 2.176037e-14 
#> iter: 1  mgc: 2.176037e-14 
#> outer mgc:  52.95131 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58850703    -2.90449477    -3.03738024    -2.70551138    -1.70319658 
#>           beta      logit_rho 
#>     0.04876147     0.12512130 
#> iter: 1  value: -97358.75 mgc: 0.2799902 ustep: 1 
#> iter: 2  mgc: 2.930989e-14 
#> iter: 1  mgc: 2.930989e-14 
#> outer mgc:  48.99427 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58883157    -2.90572292    -3.03820657    -2.70818302    -1.70292197 
#>           beta      logit_rho 
#>     0.04921354     0.13229929 
#> iter: 1  value: -97409.57 mgc: 0.2783803 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58911566    -2.90741274    -3.03911930    -2.71247259    -1.70226706 
#>           beta      logit_rho 
#>     0.04921717     0.14443336 
#> iter: 1  value: -97480.11 mgc: 0.3746282 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59047888    -2.91552121    -3.04349896    -2.73305577    -1.69912451 
#>           beta      logit_rho 
#>     0.04923459     0.20265783 
#> iter: 1  value: -97819.57 mgc: 1.823867 ustep: 1 
#> iter: 2  mgc: 4.263256e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.59566915    -2.94639308    -3.06017391    -2.81142336    -1.68715971 
#>           beta      logit_rho 
#>     0.04930092     0.42433934 
#> iter: 1  value: -99126.08 mgc: 7.322986 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  76.65505 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.65932755    -2.97267180    -3.16902308    -2.73847282    -1.67400384 
#>           beta      logit_rho 
#>     0.04895014     0.70937704 
#> iter: 1  value: -100755.9 mgc: 10.97526 ustep: 1 
#> iter: 2  mgc: 2.731149e-14 
#> iter: 1  mgc: 2.731149e-14 
#> outer mgc:  234.5492 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.6473840     -3.0542462     -3.2136359     -2.7407488     -1.6825953 
#>           beta      logit_rho 
#>      0.0491167      1.0167306 
#> iter: 1  value: -101226.8 mgc: 3.14008 ustep: 1 
#> iter: 2  mgc: 2.753353e-14 
#> iter: 1  mgc: 2.753353e-14 
#> outer mgc:  33.60164 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61486076    -3.20553767    -3.39458418    -2.75630691    -1.67612488 
#>           beta      logit_rho 
#>     0.04886509     0.80141824 
#> iter: 1  value: -103416.1 mgc: 18.64659 ustep: 1 
#> iter: 2  mgc: 2.153833e-14 
#> iter: 1  mgc: 2.153833e-14 
#> outer mgc:  35.4598 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61134433    -3.24578211    -3.43906324    -2.75699339    -1.67459768 
#>           beta      logit_rho 
#>     0.04883063     1.11720203 
#> iter: 1  value: -103986.7 mgc: 4.022304 ustep: 1 
#> iter: 2  mgc: 1.776357e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.60243845    -3.34770730    -3.55171323    -2.75873199    -1.67072979 
#>           beta      logit_rho 
#>     0.04874338     1.91697264 
#> iter: 1  value: -105437.3 mgc: 10.9133 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> iter: 1  mgc: 3.552714e-14 
#> outer mgc:  98.14859 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.62459647    -3.36699967    -3.61573296    -2.75243971    -1.66987910 
#>           beta      logit_rho 
#>     0.04883404     1.56046493 
#> iter: 1  value: -106694.8 mgc: 5.804832 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  24.39304 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.61965268    -3.37971656    -3.70646597    -2.75671657    -1.66878847 
#>           beta      logit_rho 
#>     0.04886916     1.59663636 
#> iter: 1  value: -107990.8 mgc: 8.475588 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  14.33749 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.6202294     -3.3786477     -3.9023487     -2.7557845     -1.6659496 
#>           beta      logit_rho 
#>      0.0488563      1.6505032 
#> iter: 1  value: -110870.3 mgc: 20.49945 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  19.81594 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.62363227    -3.35093695    -4.08286152    -2.75352810    -1.66439666 
#>           beta      logit_rho 
#>     0.04886242     1.61773070 
#> iter: 1  value: -113579.7 mgc: 18.56808 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  14.04646 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.62869142    -3.32230837    -4.25519301    -2.75154587    -1.66295165 
#>           beta      logit_rho 
#>     0.04886617     1.55530843 
#> iter: 1  value: -116210.4 mgc: 17.62952 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  13.89466 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63329854    -3.30311385    -4.43538945    -2.74919287    -1.66219334 
#>           beta      logit_rho 
#>     0.04886331     1.51556719 
#> iter: 1  value: -118951.9 mgc: 18.56443 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  5.366624 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63457911    -3.29808821    -4.61883404    -2.74899651    -1.66128464 
#>           beta      logit_rho 
#>     0.04886948     1.49778122 
#> iter: 1  value: -121712 mgc: 19.01907 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  10.53458 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.6344104     -3.3017127     -4.8029138     -2.7488724     -1.6612522 
#>           beta      logit_rho 
#>      0.0488676      1.5076548 
#> iter: 1  value: -124460.9 mgc: 19.03789 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  2.372546 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63347094    -3.30688431    -4.98697393    -2.74954035    -1.66080912 
#>           beta      logit_rho 
#>     0.04887319     1.51712106 
#> iter: 1  value: -127206.6 mgc: 19.09496 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  5.830818 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63166534    -3.31512082    -5.27178205    -2.75042139    -1.66074231 
#>           beta      logit_rho 
#>     0.04887545     1.53579716 
#> iter: 1  value: -131451.6 mgc: 32.8802 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  2.80216 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63154712    -3.31735521    -5.50441080    -2.75049709    -1.66058837 
#>           beta      logit_rho 
#>     0.04887669     1.54127392 
#> iter: 1  value: -134936.7 mgc: 25.39946 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  3.245582 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63159806    -3.31647860    -5.80363288    -2.75047770    -1.66054760 
#>           beta      logit_rho 
#>     0.04887655     1.53894246 
#> iter: 1  value: -139424 mgc: 35.11566 ustep: 1 
#> iter: 2  mgc: 4.263256e-14 
#> iter: 1  mgc: 4.263256e-14 
#> outer mgc:  1.580592 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63223890    -3.31426049    -6.10282394    -2.75017604    -1.66045469 
#>           beta      logit_rho 
#>     0.04887659     1.53454489 
#> iter: 1  value: -143918.5 mgc: 35.12453 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  1.533942 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63252294    -3.31222360    -6.40201275    -2.75003999    -1.66045496 
#>           beta      logit_rho 
#>     0.04887606     1.52986842 
#> iter: 1  value: -148409.5 mgc: 35.11979 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  0.371158 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63277846    -3.31143593    -6.70124026    -2.74991886    -1.66041530 
#>           beta      logit_rho 
#>     0.04887611     1.52836647 
#> iter: 1  value: -152900.6 mgc: 35.13451 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  0.559867 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63272284    -3.31147299    -7.00047270    -2.74994463    -1.66043091 
#>           beta      logit_rho 
#>     0.04887602     1.52833914 
#> iter: 1  value: -157388.4 mgc: 35.13146 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> iter: 1  mgc: 3.552714e-14 
#> outer mgc:  0.1341539 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63265950    -3.31192755    -7.29970303    -2.74997678    -1.66042049 
#>           beta      logit_rho 
#>     0.04887617     1.52936742 
#> iter: 1  value: -161875.9 mgc: 35.13531 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> iter: 1  mgc: 3.552714e-14 
#> outer mgc:  0.2112105 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63253501    -3.31249622    -7.66192285    -2.75003551    -1.66043323 
#>           beta      logit_rho 
#>     0.04887625     1.53059720 
#> iter: 1  value: -167307.6 mgc: 45.6074 ustep: 1 
#> iter: 2  mgc: 2.4869e-14 
#> iter: 1  mgc: 2.4869e-14 
#> outer mgc:  0.0292461 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.6324945     -3.3127166     -7.9844774     -2.7500565     -1.6604282 
#>           beta      logit_rho 
#>      0.0488763      1.5310691 
#> iter: 1  value: -172145.4 mgc: 38.86072 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  0.1150755 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63249512    -3.31272618    -8.30703232    -2.75005520    -1.66043162 
#>           beta      logit_rho 
#>     0.04887631     1.53110649 
#> iter: 1  value: -176983.7 mgc: 38.85968 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  0.03277101 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63251584    -3.31260223    -8.62958712    -2.75004591    -1.66042887 
#>           beta      logit_rho 
#>     0.04887628     1.53081951 
#> iter: 1  value: -181822.3 mgc: 38.86044 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  0.03904056 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63254170    -3.31249160    -8.95214197    -2.75003297    -1.66042886 
#>           beta      logit_rho 
#>     0.04887626     1.53059273 
#> iter: 1  value: -186661 mgc: 38.86019 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> iter: 1  mgc: 3.552714e-14 
#> outer mgc:  0.002561467 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63255228    -3.31242180    -9.29533701    -2.75002817    -1.66042795 
#>           beta      logit_rho 
#>     0.04887624     1.53042975 
#> iter: 1  value: -191809 mgc: 42.30543 ustep: 1 
#> iter: 2  mgc: 2.4869e-14 
#> iter: 1  mgc: 2.4869e-14 
#> outer mgc:  0.00748656 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63255511    -3.31242103    -9.62888898    -2.75002659    -1.66042812 
#>           beta      logit_rho 
#>     0.04887625     1.53043665 
#> iter: 1  value: -196812.3 mgc: 40.67815 ustep: 1 
#> iter: 2  mgc: 2.275957e-14 
#> iter: 1  mgc: 2.275957e-14 
#> outer mgc:  0.006993192 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63254830    -3.31244538    -9.96244094    -2.75003002    -1.66042824 
#>           beta      logit_rho 
#>     0.04887625     1.53048349 
#> iter: 1  value: -201815.5 mgc: 40.67819 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  0.004353737 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63254361    -3.31247651   -10.31479365    -2.75003218    -1.66042848 
#>           beta      logit_rho 
#>     0.04887626     1.53055585 
#> iter: 1  value: -207100.8 mgc: 43.88001 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  0.000807983 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63253990    -3.31248952   -10.64851665    -2.75003405    -1.66042859 
#>           beta      logit_rho 
#>     0.04887626     1.53058078 
#> iter: 1  value: -212106.6 mgc: 40.70676 ustep: 1 
#> iter: 2  mgc: 2.4869e-14 
#> iter: 1  mgc: 2.4869e-14 
#> outer mgc:  0.001783807 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.63253990    -3.31248952   -10.64851665    -2.75003405    -1.66042859 
#>           beta      logit_rho 
#>     0.04887626     1.53058078 
#> iter: 1  mgc: 2.4869e-14 
#> iter: 1  mgc: 2.4869e-14 
#> outer mgc:  0.001783807 
#> iter: 1  value: -212088.7 mgc: 0.07284539 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  3.138794 
#> iter: 1  value: -212124.5 mgc: 0.07299122 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  3.132504 
#> iter: 1  value: -212106.7 mgc: 0.007310772 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  0.1605147 
#> iter: 1  value: -212106.4 mgc: 0.007303464 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.1635888 
#> iter: 1  value: -212091.6 mgc: 0.0856802 ustep: 1 
#> mgc: 2.131628e-14 
#> outer mgc:  0.001784218 
#> iter: 1  value: -212121.6 mgc: 0.08585174 ustep: 1 
#> mgc: 2.131628e-14 
#> outer mgc:  0.001783397 
#> iter: 1  value: -212093.2 mgc: 0.07555746 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  2.569838 
#> iter: 1  value: -212119.9 mgc: 0.07570872 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  2.563599 
#> iter: 1  value: -212102.6 mgc: 0.0856802 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  26.346 
#> iter: 1  value: -212110.5 mgc: 0.08585174 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  26.35956 
#> iter: 1  value: -212106.4 mgc: 0.3670586 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  221.7124 
#> iter: 1  value: -212106.4 mgc: 0.3670586 ustep: 1 
#> iter: 2  mgc: 3.130829e-14 
#> outer mgc:  221.7094 
#> iter: 1  value: -212106.5 mgc: 0.0007226057 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  0.03299283 
#> iter: 1  value: -212106.6 mgc: 0.0007230713 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.03595518 
#> outer mgc:  0.639256 
#> [1] "calculating liability scale variance components"
```

`vc_a_lia`, `vc_c_par_lia`, `vc_c_sib_lia`, `gvc_s_lia`, represent the
G, P, C, and S variance compoenents respectively. The model estimates
are close to the true simulated variances.

``` r
model_output$vc_rpt
#> $outrprt
#>            param     Estimate    Std.Error
#>  1:     vc_a_lia 3.651485e-01 4.529518e-02
#>  2: vc_c_par_lia 4.979619e-08 9.412354e-05
#>  3: vc_c_sib_lia 1.944060e-01 3.099008e-02
#>  4:    gvc_s_lia 4.573282e-02 1.636946e-02
#>  5:         vc_a 5.168981e-03 7.163557e-04
#>  6:     vc_c_par 5.633988e-10 1.064923e-06
#>  7:     vc_c_sib 4.086493e-03 7.053142e-04
#>  8:        gvc_s 4.874786e-04 1.852056e-04
#>  9:    vc_pp_lia 4.573287e-02 1.636971e-02
#> 10:    vc_ps_lia 2.283071e-01 2.139486e-02
#> 11:    vc_ss_lia 4.227131e-01 2.656977e-02
#> 12:    gvc_s_lia 4.573282e-02 1.636946e-02
#> 13:        vc_pp 4.874791e-04 1.852084e-04
#> 14:        vc_ps 3.071969e-03 3.679613e-04
#> 15:        vc_ss 7.158462e-03 6.711122e-04
#> 16:        gvc_s 4.874786e-04 1.852056e-04
#> 
#> $vclia.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      2.051653e-03  1.311427e-09 -6.473260e-04 -3.231326e-04
#> vc_c_par  1.311427e-09  8.859241e-09 -3.002625e-10 -4.161724e-10
#> vc_c_sib -6.473260e-04 -3.002625e-10  9.603849e-04 -3.242300e-05
#> gvc_s    -3.231326e-04 -4.161724e-10 -3.242300e-05  2.679593e-04
#> 
#> $vcobs.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      5.131656e-07  1.962860e-13 -1.833446e-07 -2.719700e-08
#> vc_c_par  1.962860e-13  1.134060e-12 -6.390832e-14 -4.887917e-14
#> vc_c_sib -1.833446e-07 -6.390832e-14  4.974681e-07  4.362805e-10
#> gvc_s    -2.719700e-08 -4.887918e-14  4.362805e-10  3.430113e-08
#> 
#> $rlia.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 0.0002679673 0.0001063932 7.396990e-05 2.679588e-04
#> vc_ps 0.0001063932 0.0004577400 1.016540e-04 1.063930e-04
#> vc_ss 0.0000739699 0.0001016540 7.059529e-04 7.396996e-05
#> gvc_s 0.0002679588 0.0001063930 7.396996e-05 2.679593e-04
#> 
#> $robs.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 3.430216e-08 2.070268e-08 2.113889e-08 3.430108e-08
#> vc_ps 2.070268e-08 1.353955e-07 4.415952e-08 2.070263e-08
#> vc_ss 2.113889e-08 4.415952e-08 4.503916e-07 2.113891e-08
#> gvc_s 3.430108e-08 2.070263e-08 2.113891e-08 3.430113e-08
```

The full model report from TMB with all fitted parameters are found in
the `tmb_rpt` object.

``` r
model_output$tmb_rpt
#>                 param      Estimate   Std. Error            ll            ul
#>     1:     log_sdvc_a -2.632540e+00 6.929372e-02 -2.768353e+00 -2.496727e+00
#>     2:     log_sdvc_s -3.312490e+00 1.899630e-01 -3.684810e+00 -2.940169e+00
#>     3: log_sdvc_c_par -1.064852e+01 9.450878e+02 -1.862987e+03  1.841690e+03
#>     4: log_sdvc_c_sib -2.750034e+00 8.629822e-02 -2.919175e+00 -2.580893e+00
#>     5:   log_sdvc_res -1.660429e+00 1.114126e-02 -1.682265e+00 -1.638592e+00
#>    ---                                                                      
#> 51240:       vc_c_par  5.633988e-10 1.064923e-06 -2.086647e-06  2.087774e-06
#> 51241:       vc_c_sib  4.086493e-03 7.053142e-04  2.704103e-03  5.468884e-03
#> 51242:         vc_res  3.612186e-02 8.048858e-04  3.454431e-02  3.769940e-02
#> 51243:           vc_s  1.326808e-03 5.040889e-04  3.388122e-04  2.314804e-03
#> 51244:            rho  8.220913e-01 8.939701e-02  6.468763e-01  9.973062e-01
#>        objective                     mess fit_time sd_time modeltype
#>     1: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#>     2: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#>     3: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#>     4: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#>     5: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#>    ---                                                              
#> 51240: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#> 51241: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#> 51242: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#> 51243: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
#> 51244: -2557.409 relative convergence (4)   13.598    2.02  car_gmrf
```

The `indlevel_dt` object is a data-table showing all fitted effects and
the response variable.

``` r
model_output$indlevel_dt
#>        V1 y         fe          us           ua          upar         usib
#>     1:  1 0 0.04887626 -0.36204332 -0.006879988 -8.073344e-10 -0.002927914
#>     2:  1 0 0.04887626 -0.36204332 -0.006879988 -8.073344e-10 -0.002927914
#>     3:  1 0 0.04887626 -0.36204332 -0.002246116 -3.462249e-10 -0.005022536
#>     4:  1 0 0.04887626 -0.36204332 -0.002246116 -3.462249e-10 -0.005022536
#>     5:  1 0 0.04887626 -0.05514145 -0.009035055 -1.060221e-09 -0.003845044
#>    ---                                                                    
#> 19996:  1 0 0.04887626  0.11979101 -0.003350712 -5.164915e-10 -0.007492521
#> 19997:  1 0 0.04887626  0.28864290 -0.016518465  1.121443e-08 -0.004357181
#> 19998:  1 1 0.04887626  0.28864290  0.097392614  1.121443e-08  0.085698651
#> 19999:  1 0 0.04887626  0.28864290 -0.007783190 -1.199730e-09 -0.017403977
#> 20000:  1 0 0.04887626  0.28864290 -0.007783190 -1.199730e-09 -0.017403977
#>               yhat
#>     1: -0.32297497
#>     2: -0.32297497
#>     3: -0.32043572
#>     4: -0.32043572
#>     5: -0.01914529
#>    ---            
#> 19996:  0.15782404
#> 19997:  0.31664353
#> 19998:  0.52061044
#> 19999:  0.31233199
#> 20000:  0.31233199
```
