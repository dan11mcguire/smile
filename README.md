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
#> [1] 0.3176212
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
#>   2.000   4.000   5.000   5.942   7.000  30.000
summary(apply(sigma_s, 2, function(x) sum(x!=0)))   
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     2.0    46.0   205.0   179.9   289.0   298.0
chsigma_s<-t(chol(sigma_s))

ustar<-rnorm(nrow(wmatd), sd=1)*sqrt(sets2)
var(ustar)
#> [1] 0.07381008
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
#>    set emprel ttype pid empreltype         u_a       u_sib       u_par u_fam
#> 1:   1      1  quad   1     parent -0.29958885 -0.14564158 -0.18709948     0
#> 2:   1      2  quad   2     parent  0.21218133 -0.17229776 -0.18709948     0
#> 3:   1      3  quad   3        sib  0.71505496 -0.33746791  0.18989103     0
#> 4:   1      3  quad   4        sib  0.59110056 -0.33746791  0.04438050     0
#> 5:   2      1  quad   5     parent  0.76648259  0.03070334  0.05643138     0
#> 6:   2      2  quad   6     parent  0.05321685  0.21498273  0.05643138     0
#> 7:   2      3  quad   7        sib -0.22273021 -0.34555542  0.05779187     0
#> 8:   2      3  quad   8        sib  0.03553109 -0.34555542 -0.02098226     0
#>            eps fixed   kid  gid  rn         u_s          y
#> 1: -0.69001873     0 17115 1683 269 -0.04372492 -1.3660736
#> 2: -0.27053930     0 17115 1683 269 -0.04372492 -0.4614801
#> 3:  0.74017054     0 17115 1683 269 -0.04372492  1.2639237
#> 4: -0.55231558     0 17115 1683 269 -0.04372492 -0.2980274
#> 5:  1.00178214     0 06067   72  71 -0.06495228  1.7904472
#> 6: -0.08430955     0 06067   72  71 -0.06495228  0.1753691
#> 7:  0.42950432     0 06067   72  71 -0.06495228 -0.1459417
#> 8: -0.19476727     0 06067   72  71 -0.06495228 -0.5907261

ss$d[,lapply(.SD, var),.SDcols=grep("u_|eps",names(ss$d),value=T)]
#>         u_a     u_sib     u_par u_fam       eps        u_s
#> 1: 0.391497 0.1996407 0.0494618     0 0.3177995 0.02896509


ss$d[,lapply(.SD, var),.SDcols=c(grep("^[yu]", names(ss$d),value=T))] ## Actual variance of the outcome and the random effects
#>         u_a     u_sib     u_par u_fam        u_s         y
#> 1: 0.391497 0.1996407 0.0494618     0 0.02896509 0.9887033
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
#> -0.04675 -0.04675 -0.04675 -0.04675  0.95325 
#> 
#> Coefficients:
#>   Estimate Std. Error t value Pr(>|t|)    
#> X 0.046750   0.001493   31.32   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.2111 on 19999 degrees of freedom
#> Multiple R-squared:  0.04675,    Adjusted R-squared:  0.0467 
#> F-statistic: 980.8 on 1 and 19999 DF,  p-value: < 2.2e-16
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
#>      -2.595105      -2.595105      -2.595105      -2.595105      -1.901958 
#>           beta      logit_rho 
#>       0.046750       0.000000 
#> Constructing atomic qnorm1
#> Optimizing tape... Done
#> iter: 1  value: -91791.67 mgc: 85.55722 ustep: 1 
#> iter: 2  mgc: 1.607603e-13 
#> iter: 1  mgc: 1.607603e-13 
#> Matching hessian patterns... Done
#> outer mgc:  1877.476 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.406006646   -2.626987809   -2.471441576   -2.439698167   -0.947145208 
#>           beta      logit_rho 
#>   -0.063297671    0.001692176 
#> iter: 1  value: -73083.37 mgc: 25.97058 ustep: 1 
#> iter: 2  mgc: 1.865175e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5734423439  -2.5987577516  -2.5809385120  -2.5773020149  -1.7925756116 
#>           beta      logit_rho 
#>   0.0341430356   0.0001938542 
#> iter: 1  value: -90249.7 mgc: 9.997553 ustep: 1 
#> iter: 2  mgc: 9.089951e-14 
#> iter: 1  mgc: 9.089951e-14 
#> outer mgc:  2047.292 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5709590677  -2.6029260606  -2.5865355279  -2.5781002938  -1.7773018126 
#>           beta      logit_rho 
#>   0.1474343486   0.0004699439 
#> iter: 1  value: -89267.92 mgc: 88.02176 ustep: 1 
#> iter: 2  mgc: 6.972201e-13 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5735854823  -2.5991108206  -2.5817585324  -2.5777052579  -1.7930232206 
#>           beta      logit_rho 
#>   0.0457685319   0.0002181118 
#> iter: 1  value: -90282.86 mgc: 9.351349 ustep: 1 
#> iter: 2  mgc: 5.4734e-14 
#> iter: 1  mgc: 5.4734e-14 
#> outer mgc:  143.7105 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.575157679   -2.603543177   -2.590425831   -2.581988648   -1.797446415 
#>           beta      logit_rho 
#>    0.044689537    0.000349748 
#> iter: 1  value: -90507.1 mgc: 0.9059518 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5770015087  -2.6087413238  -2.6005906029  -2.5870120905  -1.8026338164 
#>           beta      logit_rho 
#>   0.0434241209   0.0005041275 
#> iter: 1  value: -90768.9 mgc: 1.072144 ustep: 1 
#> iter: 2  mgc: 3.907985e-14 
#> iter: 1  mgc: 3.907985e-14 
#> outer mgc:  406.339 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>  -2.5770963849  -2.6195575144  -2.6179018255  -2.5938970503  -1.7921664273 
#>           beta      logit_rho 
#>   0.0517825126   0.0008394454 
#> iter: 1  value: -90998.34 mgc: 6.524069 ustep: 1 
#> iter: 2  mgc: 9.126033e-14 
#> iter: 1  mgc: 9.126033e-14 
#> outer mgc:  1114.159 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.577532948   -2.630677975   -2.636222463   -2.601421144   -1.783446573 
#>           beta      logit_rho 
#>    0.044619621    0.001246403 
#> iter: 1  value: -91276.34 mgc: 5.492267 ustep: 1 
#> iter: 2  mgc: 4.823919e-14 
#> iter: 1  mgc: 4.823919e-14 
#> outer mgc:  196.2114 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.578210596   -2.642898160   -2.655170055   -2.609520249   -1.775399178 
#>           beta      logit_rho 
#>    0.046680066    0.001711803 
#> iter: 1  value: -91578.07 mgc: 2.250939 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.579251672   -2.664583561   -2.688697248   -2.623754241   -1.760182929 
#>           beta      logit_rho 
#>    0.047422093    0.002535138 
#> iter: 1  value: -92107.47 mgc: 3.887881 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  321.6957 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>   -2.579155028   -2.723703614   -2.720486931   -2.642966685   -1.765403528 
#>           beta      logit_rho 
#>    0.040930327    0.005901951 
#> iter: 1  value: -92794.71 mgc: 4.587297 ustep: 1 
#> iter: 2  mgc: 5.084821e-14 
#> iter: 1  mgc: 5.084821e-14 
#> outer mgc:  963.1166 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.57266603    -2.78386858    -2.74939233    -2.65942380    -1.76949752 
#>           beta      logit_rho 
#>     0.05089300     0.01434707 
#> iter: 1  value: -93290.15 mgc: 5.990867 ustep: 1 
#> iter: 2  mgc: 4.296563e-14 
#> iter: 1  mgc: 4.296563e-14 
#> outer mgc:  1065.065 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.55577059    -2.83259490    -2.78485081    -2.67914080    -1.76378924 
#>           beta      logit_rho 
#>     0.04481227     0.03843631 
#> iter: 1  value: -93678.95 mgc: 3.721967 ustep: 1 
#> iter: 2  mgc: 2.331468e-14 
#> iter: 1  mgc: 2.331468e-14 
#> outer mgc:  225.275 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.57678672    -2.86103812    -2.81724239    -2.71008573    -1.74530603 
#>           beta      logit_rho 
#>     0.04691276     0.07520704 
#> iter: 1  value: -94768.46 mgc: 4.306571 ustep: 1 
#> iter: 2  mgc: 2.264855e-14 
#> iter: 1  mgc: 2.264855e-14 
#> outer mgc:  226.2297 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54667429    -2.87636345    -2.86331272    -2.68122845    -1.74997107 
#>           beta      logit_rho 
#>     0.04513519     0.10408368 
#> iter: 1  value: -94510.46 mgc: 3.394228 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  166.0171 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.57133143    -2.90536566    -2.88809324    -2.69886033    -1.74723309 
#>           beta      logit_rho 
#>     0.04693205     0.15478712 
#> iter: 1  value: -95488.37 mgc: 2.395108 ustep: 1 
#> iter: 2  mgc: 2.575717e-14 
#> iter: 1  mgc: 2.575717e-14 
#> outer mgc:  230.2086 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56594739    -2.90803554    -2.88624045    -2.74115652    -1.74368282 
#>           beta      logit_rho 
#>     0.04318197     0.21046872 
#> iter: 1  value: -95902.58 mgc: 3.413381 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  611.3947 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54717456    -2.90730267    -2.89042974    -2.74957179    -1.74395058 
#>           beta      logit_rho 
#>     0.04622749     0.23854426 
#> iter: 1  value: -95744.91 mgc: 1.634204 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  97.90537 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54813793    -2.92099794    -2.90331015    -2.72526767    -1.74244473 
#>           beta      logit_rho 
#>     0.04505224     0.30184428 
#> iter: 1  value: -95604.1 mgc: 1.782492 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  191.9624 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.58790454    -2.92878414    -2.92405437    -2.69993249    -1.73731747 
#>           beta      logit_rho 
#>     0.04665578     0.43253261 
#> iter: 1  value: -96220.59 mgc: 3.593395 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  163.005 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56766000    -2.92891737    -2.93523293    -2.72791720    -1.73620411 
#>           beta      logit_rho 
#>     0.04620387     0.56855868 
#> iter: 1  value: -96376.95 mgc: 2.257663 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.50704975    -2.92895321    -2.96823569    -2.81166315    -1.73363651 
#>           beta      logit_rho 
#>     0.04853986     0.97674731 
#> iter: 1  value: -96868.89 mgc: 7.349536 ustep: 1 
#> iter: 2  mgc: 2.531308e-14 
#> iter: 1  mgc: 2.531308e-14 
#> outer mgc:  516.6218 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.49852461    -3.03563726    -2.95530966    -2.86538766    -1.73202684 
#>           beta      logit_rho 
#>     0.04279096     1.52685346 
#> iter: 1  value: -97230.26 mgc: 4.238586 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> iter: 1  mgc: 3.552714e-14 
#> outer mgc:  571.6484 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.48561455    -3.14057747    -2.93833600    -2.90299894    -1.73300507 
#>           beta      logit_rho 
#>     0.04648896     2.07845435 
#> iter: 1  value: -97264.28 mgc: 3.273299 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  82.32277 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54419664    -3.00974903    -3.04947312    -2.92766475    -1.71207509 
#>           beta      logit_rho 
#>     0.04536853     2.61063183 
#> iter: 1  value: -100057.3 mgc: 12.05796 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.50825922    -3.09026005    -2.98127020    -2.91243974    -1.72470205 
#>           beta      logit_rho 
#>     0.04557197     2.28349816 
#> iter: 1  value: -98334.12 mgc: 4.382367 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  47.95426 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.48261367    -3.20924322    -2.99357637    -2.90872592    -1.72922791 
#>           beta      logit_rho 
#>     0.04488704     2.69986160 
#> iter: 1  value: -98009.69 mgc: 4.418873 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>     -2.4332821     -3.4335064     -3.0156661     -2.9025000     -1.7398971 
#>           beta      logit_rho 
#>      0.0400687      3.4910708 
#> iter: 1  value: -97404.24 mgc: 8.670212 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  669.0283 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.66074787    -3.97329725    -3.47946377    -3.09234558    -1.64824293 
#>           beta      logit_rho 
#>     0.04466177     4.48128701 
#> iter: 1  value: -110156.6 mgc: 71.84831 ustep: 1 
#> iter: 2  mgc: 3.552714e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54550628    -3.69982163    -3.24448877    -2.99616364    -1.69467794 
#>           beta      logit_rho 
#>     0.04233477     3.97961119 
#> iter: 1  value: -103592.4 mgc: 28.23196 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  458.3039 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.53856512    -3.78960086    -3.98646679    -2.40698837    -1.71187451 
#>           beta      logit_rho 
#>     0.06178296     4.44610280 
#> iter: 1  value: -106670.7 mgc: 145.6122 ustep: 1 
#> iter: 2  mgc: 9.947598e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54428775    -3.71454391    -3.36585370    -2.89968075    -1.69751532 
#>           beta      logit_rho 
#>     0.04900837     4.05590182 
#> iter: 1  value: -103996.7 mgc: 11.90505 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  325.8707 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.55424585    -3.76184337    -3.46619391    -2.76662219    -1.69996514 
#>           beta      logit_rho 
#>     0.04744583     4.05934977 
#> iter: 1  value: -103835.3 mgc: 9.458287 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  151.4512 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.55796390    -3.81241967    -3.45603302    -2.75948228    -1.70317278 
#>           beta      logit_rho 
#>     0.04594593     4.22485182 
#> iter: 1  value: -103665.7 mgc: 2.110034 ustep: 1 
#> iter: 2  mgc: 1.421085e-14 
#> iter: 1  mgc: 1.421085e-14 
#> outer mgc:  41.47345 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.55857190    -3.86979133    -3.44982049    -2.76699860    -1.70258220 
#>           beta      logit_rho 
#>     0.04584815     4.38838208 
#> iter: 1  value: -103675.8 mgc: 1.927144 ustep: 1 
#> iter: 2  mgc: 1.776357e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56039588    -4.04190634    -3.43118290    -2.78954758    -1.70081047 
#>           beta      logit_rho 
#>     0.04555481     4.87897285 
#> iter: 1  value: -103709.1 mgc: 5.532834 ustep: 1 
#> iter: 2  mgc: 1.776357e-14 
#> iter: 1  mgc: 1.776357e-14 
#> outer mgc:  71.31105 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56993755    -4.38708777    -3.45016949    -2.77366589    -1.70533500 
#>           beta      logit_rho 
#>     0.04593614     5.48080088 
#> iter: 1  value: -103930.2 mgc: 10.04561 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  99.49835 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.56048205    -4.66003336    -3.44553606    -2.76486640    -1.70492853 
#>           beta      logit_rho 
#>     0.04596355     6.11906712 
#> iter: 1  value: -103555.8 mgc: 7.81572 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  35.51785 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54020402    -5.03969700    -3.42475101    -2.78891112    -1.70347975 
#>           beta      logit_rho 
#>     0.04590795     6.87591655 
#> iter: 1  value: -103187.8 mgc: 10.84132 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  31.13334 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.53960089    -5.21387443    -3.41134564    -2.77496289    -1.70604205 
#>           beta      logit_rho 
#>     0.04636695     7.25251760 
#> iter: 1  value: -102797.2 mgc: 5.710803 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  19.73521 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54150605    -5.31550814    -3.40260210    -2.77577122    -1.70624488 
#>           beta      logit_rho 
#>     0.04627133     7.42797903 
#> iter: 1  value: -102709.4 mgc: 3.466 ustep: 1 
#> iter: 2  mgc: 1.776357e-14 
#> iter: 1  mgc: 1.776357e-14 
#> outer mgc:  10.1258 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54177610    -5.40244581    -3.39089586    -2.77759815    -1.70648035 
#>           beta      logit_rho 
#>     0.04618126     7.59836293 
#> iter: 1  value: -102564.8 mgc: 3.135313 ustep: 1 
#> iter: 2  mgc: 2.253753e-14 
#> iter: 1  mgc: 2.253753e-14 
#> outer mgc:  0.9363197 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54181250    -5.43338053    -3.38646627    -2.77752597    -1.70661135 
#>           beta      logit_rho 
#>     0.04616764     7.66172781 
#> iter: 1  value: -102499.1 mgc: 1.154107 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> iter: 1  mgc: 2.131628e-14 
#> outer mgc:  0.3276893 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54181528    -5.43894060    -3.38574619    -2.77742881    -1.70663557 
#>           beta      logit_rho 
#>     0.04616963     7.67315425 
#> iter: 1  value: -102487.2 mgc: 0.206619 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  0.07786416 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54181455    -5.43929821    -3.38575221    -2.77741809    -1.70663522 
#>           beta      logit_rho 
#>     0.04617024     7.67386882 
#> iter: 1  value: -102487.1 mgc: 0.01178017 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> par:
#>     log_sdvc_a     log_sdvc_s log_sdvc_c_par log_sdvc_c_sib   log_sdvc_res 
#>    -2.54181528    -5.43894060    -3.38574619    -2.77742881    -1.70663557 
#>           beta      logit_rho 
#>     0.04616963     7.67315425 
#> iter: 1  mgc: 2.842171e-14 
#> iter: 1  mgc: 2.842171e-14 
#> outer mgc:  0.07786416 
#> iter: 1  value: -102469.8 mgc: 0.07408165 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  3.596237 
#> iter: 1  value: -102504.7 mgc: 0.07422996 ustep: 1 
#> iter: 2  mgc: 2.137179e-14 
#> outer mgc:  3.658041 
#> iter: 1  value: -102487.3 mgc: 0.03206214 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.0580572 
#> iter: 1  value: -102487.2 mgc: 0.03203009 ustep: 1 
#> iter: 2  mgc: 2.250457e-14 
#> outer mgc:  0.1087673 
#> iter: 1  value: -102472.8 mgc: 0.08642074 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.7783568 
#> iter: 1  value: -102501.7 mgc: 0.08659375 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.8410665 
#> iter: 1  value: -102473.9 mgc: 0.07710613 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  2.321743 
#> iter: 1  value: -102500.6 mgc: 0.07726049 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  2.383389 
#> iter: 1  value: -102482.4 mgc: 0.08642074 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  23.56982 
#> iter: 1  value: -102492 mgc: 0.08659375 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  23.64129 
#> iter: 1  value: -102487.2 mgc: 0.06072882 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  107.3356 
#> iter: 1  value: -102487.2 mgc: 0.06072882 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  107.4913 
#> iter: 1  value: -102487.2 mgc: 5.642258e-05 ustep: 1 
#> iter: 2  mgc: 2.842171e-14 
#> outer mgc:  0.05949866 
#> iter: 1  value: -102487.2 mgc: 5.647898e-05 ustep: 1 
#> iter: 2  mgc: 2.131628e-14 
#> outer mgc:  0.09626116 
#> outer mgc:  0.7346391 
#> [1] "calculating liability scale variance components"
```

`vc_a_lia`, `vc_c_par_lia`, `vc_c_sib_lia`, `gvc_s_lia`, represent the
G, P, C, and S variance compoenents respectively. The model estimates
are close to the true simulated variances.

``` r
model_output$vc_rpt
#> $outrprt
#>            param     Estimate    Std.Error
#>  1:     vc_a_lia 0.4422652241 0.0922226871
#>  2: vc_c_par_lia 0.0953335344 0.0530611478
#>  3: vc_c_sib_lia 0.1813581547 0.0334935281
#>  4:    gvc_s_lia 0.0393341731 0.1004538744
#>  5:         vc_a 0.0061973682 0.0006901193
#>  6:     vc_c_par 0.0011459831 0.0006432917
#>  7:     vc_c_sib 0.0038686192 0.0006483096
#>  8:        gvc_s 0.0003968471 0.0010681808
#>  9:    vc_pp_lia 0.1346677075 0.0908880532
#> 10:    vc_ps_lia 0.2604667852 0.0617653978
#> 11:    vc_ss_lia 0.4418249398 0.0472668689
#> 12:    gvc_s_lia 0.0393341731 0.1004538744
#> 13:        vc_pp 0.0015428303 0.0012231453
#> 14:        vc_ps 0.0034955312 0.0010930277
#> 15:        vc_ss 0.0073641505 0.0011987380
#> 16:        gvc_s 0.0003968471 0.0010681808
#> 
#> $vclia.cov
#>                  vc_a      vc_c_par      vc_c_sib        gvc_s
#> vc_a      0.008505024  0.0021834325  0.0008739610 -0.008402273
#> vc_c_par  0.002183433  0.0028154854  0.0003059209 -0.002322914
#> vc_c_sib  0.000873961  0.0003059209  0.0011218164 -0.001788292
#> gvc_s    -0.008402273 -0.0023229140 -0.0017882925  0.010090981
#> 
#> $vcobs.cov
#>                   vc_a      vc_c_par      vc_c_sib         gvc_s
#> vc_a      4.762646e-07  6.104364e-08 -1.735852e-07 -6.536681e-08
#> vc_c_par  6.104364e-08  4.138243e-07 -2.317641e-08 -2.937506e-08
#> vc_c_sib -1.735852e-07 -2.317641e-08  4.203053e-07 -2.228426e-09
#> gvc_s    -6.536681e-08 -2.937506e-08 -2.228426e-09  1.141010e-06
#> 
#> $rlia.cov
#>             vc_pp       vc_ps       vc_ss       gvc_s
#> vc_pp 0.008260638 0.004658647 0.003176275 0.007768067
#> vc_ps 0.004658647 0.003814964 0.002463652 0.005889845
#> vc_ss 0.003176275 0.002463652 0.002234157 0.004101552
#> gvc_s 0.007768067 0.005889845 0.004101552 0.010090981
#> 
#> $robs.cov
#>              vc_pp        vc_ps        vc_ss        gvc_s
#> vc_pp 1.496084e-06 1.109474e-06 1.084069e-06 1.111635e-06
#> vc_ps 1.109474e-06 1.194710e-06 1.105688e-06 1.108327e-06
#> vc_ss 1.084069e-06 1.105688e-06 1.436973e-06 1.106098e-06
#> gvc_s 1.111635e-06 1.108327e-06 1.106098e-06 1.141010e-06
```

The full model report from TMB with all fitted parameters are found in
the `tmb_rpt` object.

``` r
model_output$tmb_rpt
#>                 param      Estimate   Std. Error            ll            ul
#>     1:     log_sdvc_a -2.541815e+00 5.567841e-02 -2.650943e+00 -2.4326875966
#>     2:     log_sdvc_s -5.438941e+00 1.345834e+00 -8.076727e+00 -2.8011542532
#>     3: log_sdvc_c_par -3.385746e+00 2.806724e-01 -3.935854e+00 -2.8356383324
#>     4: log_sdvc_c_sib -2.777429e+00 8.379082e-02 -2.941656e+00 -2.6132018165
#>     5:   log_sdvc_res -1.706636e+00 1.535023e-02 -1.736721e+00 -1.6765496654
#>    ---                                                                      
#> 51294:       vc_c_par  1.145983e-03 6.432917e-04 -1.148455e-04  0.0024068118
#> 51295:       vc_c_sib  3.868619e-03 6.483096e-04  2.597956e-03  0.0051392827
#> 51296:         vc_res  3.293329e-02 1.011067e-03  3.095164e-02  0.0349149505
#> 51297:           vc_s  1.887106e-05 5.079462e-05 -8.068457e-05  0.0001184267
#> 51298:            rho  9.995351e-01 1.265331e-03  9.970551e-01  1.0020150703
#>        objective                     mess fit_time sd_time modeltype
#>     1:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#>     2:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#>     3:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#>     4:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#>     5:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#>    ---                                                              
#> 51294:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#> 51295:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#> 51296:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#> 51297:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
#> 51298:  -2894.29 relative convergence (4)    7.659   1.996  car_gmrf
```

The `indlevel_dt` object is a data-table showing all fitted effects and
the response variable.

``` r
model_output$indlevel_dt
#>        V1 y         fe         us            ua          upar         usib
#>     1:  1 0 0.04616963  0.9464858 -0.0114990864 -0.0022736906 -0.003837772
#>     2:  1 0 0.04616963  0.9464858 -0.0114990864 -0.0022736906 -0.003837772
#>     3:  1 0 0.04616963  0.9464858 -0.0037838236 -0.0009895023 -0.006680740
#>     4:  1 0 0.04616963  0.9464858 -0.0037838236 -0.0009895023 -0.006680740
#>     5:  1 1 0.04616963  0.2659179  0.1223763028  0.0235894039  0.084801347
#>    ---                                                                    
#> 19996:  1 0 0.04616963 -8.5490602 -0.0006796718 -0.0001777400 -0.001200032
#> 19997:  1 0 0.04616963 -0.3286268  0.0429797412 -0.0051428862 -0.008680699
#> 19998:  1 0 0.04616963 -0.3286268  0.0429797412 -0.0051428862 -0.008680699
#> 19999:  1 0 0.04616963 -0.3286268 -0.0187110880 -0.0048931098  0.071020386
#> 20000:  1 1 0.04616963 -0.3286268  0.0991598450  0.0259311492  0.071020386
#>               yhat
#>     1:  0.97504492
#>     2:  0.97504492
#>     3:  0.98120140
#>     4:  0.98120140
#>     5:  0.54285457
#>    ---            
#> 19996: -8.50494803
#> 19997: -0.25330100
#> 19998: -0.25330100
#> 19999: -0.23504097
#> 20000: -0.08634578
```
