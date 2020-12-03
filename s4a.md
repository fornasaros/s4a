Evaluation of SERS spectra from gingival crevicular fluid (GCF)
================

## First Glance at Data

``` r
#load("gcf.RData")
spc
```

    ## hyperSpec object
    ##    67 spectra
    ##    43 data columns
    ##    2047 data points / spectrum
    ## wavelength: paste("Raman shift (", cm^-1, ")", sep = "") [numeric] -307.40 -304.81 ... 3270.66 
    ## data:  (67 rows x 43 columns)
    ##    1. spc: Intensity (a.u.) [matrix, array2047] 32.72432 13.40092 ... 17.48172 
    ##    2. laser: laser [integer] 785 785 ... 785 
    ##    3. power: power [integer] 10 10 ... 10 
    ##    4. patient: patient [character] X10 X10 ... X7 
    ##    5. tooth: tooth [character] 11 13 ... 36 
    ##    6. ID: ID [character] X10_11 X10_13 ... X7_36 
    ##    7. Batch:  [numeric] 1 1 ... 1 
    ##    8. Tipo dente:  [numeric] 3 3 ... 1 
    ##    9. Lato:  [numeric] 1 1 ... 2 
    ##    10. Sup/Inf:  [numeric] 1 1 ... 2 
    ##    11. Genere:  [numeric] 1 1 ... 0 
    ##    12. Età:  [numeric] 57 57 ... 44 
    ##    13. Gravidanza:  [numeric] NA NA ... 0 + NA
    ##    14. Fase mestruale:  [character] NA NA ... 21 + NA
    ##    15. Fumo:  [numeric] 0 0 ... 0 
    ##    16. Patologie sistemiche:  [numeric] 0 0 ... 2 
    ##    17. Assunzione farmaci:  [numeric] 1 1 ... 2 
    ##    18. Sondaggio DV:  [numeric] 2 2 ... 6 
    ##    19. Placca DV:  [numeric] 1 1 ... 1 
    ##    20. BoP DV:  [numeric] 0 0 ... 1 
    ##    21. Sondaggio V:  [numeric] 2 2 ... 3 
    ##    22. Placca V:  [numeric] 0 1 ... 0 
    ##    23. BoP V:  [numeric] 0 0 ... 0 
    ##    24. Sondaggio MV:  [numeric] 3 3 ... 4 
    ##    25. Placca MV:  [numeric] 1 1 ... 1 
    ##    26. BoP MV:  [numeric] 0 0 ... 1 
    ##    27. Sondaggio DL:  [numeric] 2 3 ... 3 
    ##    28. Placca DL:  [numeric] 1 1 ... 1 
    ##    29. BoP DL:  [numeric] 1 1 ... 1 
    ##    30. Sondaggio L:  [numeric] 2 3 ... 3 
    ##    31. Placca L:  [numeric] 0 0 ... 1 
    ##    32. BoP L:  [numeric] 0 0 ... 1 
    ##    33. Sondaggio ML:  [numeric] 2 3 ... 4 
    ##    34. Placca ML:  [numeric] 1 1 ... 1 
    ##    35. BoP ML:  [numeric] 0 0 ... 1 
    ##    36. Diagnosi:  [character] Gingivitis    Periodontitis ... Periodontitis 
    ##    37. PSR 1:  [character] 4 4 ... 4 
    ##    38. PSR 2:  [character] 3 3 ... 3 
    ##    39. PSR 3:  [numeric] 4 4 ... 4 
    ##    40. PSR 4:  [character] 3 3 ... 4 
    ##    41. PSR 5:  [numeric] 3 3 ... 3 
    ##    42. PSR 6:  [numeric] 4 4 ... 4 
    ##    43. class: class [factor] 1 2 ... 2

### First Glance at Spectra

``` r
plot (spc) #raw spectra
```

![](s4a_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#View(spc$..) #metadata
```

## Preprocessing

### Selection of spectral range

``` r
spc.low <- spc [,, 400 ~ 1800]
```

### Baseline correction

``` r
spc.bc <- baseline(unclass(spc.low$spc), lambda=4, p= 0.001,method="als")
spc.low[[]] <- getCorrected(spc.bc)
```

### Total intensity normalization

``` r
factors <- apply (spc.low, 1, function (x){sqrt (sum (x^2))})
spc.low <- sweep (spc.low, 1, factors, "/")
plot(spc.low[,,400~1703], "spcprctile")
```

![](s4a_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Spectral fitting

``` r
#load("metabolites.RData")

# pre-processing spcmet

spcmet <- spcmet[,,400~1800]
bl <- baseline (spcmet [[]],  method='als',lambda = 4, p = 0.001, maxit = 100) # asymmetric least squares method
spcmet [[]] <- getCorrected (bl)
rm(bl)

# normalization
factor <- apply (spcmet, 1,function(x){sqrt(sum(x^2))})
spcmet <- sweep (spcmet, 1, factor, "/")
plot(spcmet, stacked=T)
```

![](s4a_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
medspc <- apply(spc.low$spc, 2, median) #obtain median GFC spectrum

# FITTING ######
y <- as.vector(medspc)
myfunc <- function(c1, c2, c3, c4,c5){
  c(spcmet[1]$spc) * c1 + c(spcmet[2]$spc) * c2 + c(spcmet[3]$spc)* c3 +
  c(spcmet[4]$spc)* c4 + c(spcmet[5]$spc)* c5
}
myfit <- nls(y ~ myfunc(c1, c2, c3, c4,c5), 
               start=list(c1=0.9,c2=0.5, c3=0.1, c4=0.9, c5=0.5))
coef(myfit)
```

    ##        c1        c2        c3        c4        c5 
    ## 0.2884725 0.2313288 0.4144316 0.0625502 0.2173057

``` r
spc.fit <- myfunc(coef(myfit)[1], coef(myfit)[2], coef(myfit)[3], 
                  coef(myfit)[4], coef(myfit)[5])

medres<- medspc - spc.fit


#calculate fitting statistics
norm_vec <- function(x) sqrt(sum(x^2))
predvals<-myfit$m$fitted()
lof<-sqrt((sum(medspc-predvals)^2)/(sum(medspc)^2))*100
rfe<- norm_vec(medspc-medres)/norm_vec(medspc)*100
r2<-(1-((sum(medspc-predvals)^2)/(sum(medspc)^2)))*100

data.frame(lof,rfe,r2)
```

    ##       lof      rfe       r2
    ## 1 5.40242 95.99903 99.70814

### Plot pure components, median gfc spectrum, best-fitted spectrum and residuals

``` r
off <- 0.2 #0.12
labcex <- 0.8 #0.8

plot(spc.low@wavelength, medspc+off*5, type="l", col="darkgrey", lwd=3,xaxs="i",
     ylim=c(0,off*5 + 0.25), xlim=c(400,1800),
     xlab = "", ylab = "", yaxt="n")
mtext(side=2, "Normalized Intensity (a.u.)",font=1, padj=-0.85)
mtext(side=1, expression(paste("Raman shift (", cm^-1, ")", sep = 
                                 "")),font=1, padj=2.25)
lines(spcmet@wavelength, spc.fit+off*5, col="black", lty=2) # fitted
lines(spc.low@wavelength, medres+off*5+0.19, type="l") # residuals
abline(h=off*5+0.19, lty=2)
lines(spcmet@wavelength, spcmet[2]$spc * coef(myfit)[2] + off*4, 
      col=rgb(0,0,0,1))
lines(spcmet@wavelength, spcmet[4]$spc * coef(myfit)[4] + off*3, 
      col=rgb(0,0,0,1))
lines(spcmet@wavelength, spcmet[3]$spc * coef(myfit)[3] + off*2, 
      col=rgb(0,0,0,1))
lines(spcmet@wavelength, spcmet[5]$spc * coef(myfit)[5] + off*1, 
      col=rgb(0,0,0,1))
lines(spcmet@wavelength, spcmet[1]$spc * coef(myfit)[1] + off*0, 
      col=rgb(0,0,0,1))# pure components
```

![](s4a_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Difference spectra

``` r
h<-spc.low[spc.low$Diagnosi=="Healthy"]
p<-spc.low[spc.low$Diagnosi=="Periodontitis"]
g<-spc.low[spc.low$Diagnosi=="Gingivitis"]

plot(p, "spcprctile", col=alpha(2,0.5)) #"#353F6B"
plot(h, "spcprctile", add=T, col=alpha("#353F6B", 0.5)) #"#d11141"
```

![](s4a_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
plot(g, "spcprctile", col=alpha("orange",0.5)) #"#353F6B"
plot(h, "spcprctile", add=T, col=alpha("#353F6B", 0.5)) #"#d11141"
```

![](s4a_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
plot(p, "spcprctile", col=alpha(2,0.5)) #"#353F6B"
plot(g, "spcprctile", add=T, col=alpha("orange", 0.5)) #"#d11141"
```

![](s4a_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
ph <- list()
for (i in 1:nrow(h)) {
  ph[[i]] <-  p - h[i]
}

ph<-hyperSpec::collapse(ph)

gh <- list()
for (i in 1:nrow(h)) {
  gh[[i]] <-  g - h[i]
}

gh<-hyperSpec::collapse(gh)

pg <- list()
for (i in 1:nrow(g)) {
  pg[[i]] <-  p - g[i]
}

pg<-hyperSpec::collapse(pg)


plot(ph, "spcprctile", col=alpha(1, 0.5))
abline(h=0, lty="dashed")
```

![](s4a_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

``` r
plot(gh, "spcprctile", col=alpha(1, 0.5))
abline(h=0, lty="dashed")
```

![](s4a_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

``` r
plot(pg, "spcprctile", col=alpha(1, 0.5))
abline(h=0, lty="dashed")
```

![](s4a_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->

# PCA

``` r
myPCA <- PCA(scale(spc.low[[]], scale=F))
row.names(myPCA$scores)<-(spc.low$ID)
scoreplot(myPCA, pc= c(1,5),show.names = T,pch=19)
```

![](s4a_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
myloadings <- decomposition (spc.low[], t(myPCA$loadings), scores = FALSE,
                             label.wavelength = "PC", label.spc = "loading I / a.u.")
myscores <- decomposition (spc.low[], myPCA$scores, label.wavelength = "PC",
                           label.spc = "score / a.u.")

#screeplot(myPCA, type = "percentage", 15, main = "Explained variance")

plot(myloadings[c(1:6)], stacked=T)
```

![](s4a_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
response<-data.frame("PC.1"=myscores[,,1]$spc, class=spc.low$class)
colnames(response)<-c("PC1", "class")
#ggboxplot(response, x = "class", y = "PC1")
res.kruskal <- response %>% kruskal_test(PC1 ~ class)
pwc <- response %>% 
  dunn_test(PC1 ~ class, p.adjust.method = "bonferroni") 
pwc
```

    ## # A tibble: 3 x 9
    ##   .y.   group1 group2    n1    n2 statistic        p   p.adj p.adj.signif
    ## * <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl>   <dbl> <chr>       
    ## 1 PC1   0      1         17    23    -2.73  0.00642  0.0192  *           
    ## 2 PC1   0      2         17    27    -3.41  0.000654 0.00196 **          
    ## 3 PC1   1      2         23    27    -0.647 0.518    1       ns

``` r
pwc <- pwc %>% add_xy_position(x = "class")
p1<-ggboxplot(response, x = "class", y = "PC1", outlier.shape = NA) +
  geom_point(aes(x=class,y=PC1,color=class,fill=class), 
             position=position_jitterdodge(jitter.width=0.1,
                                           dodge.width=0.7),
             size=4, alpha=0.7)+
  scale_colour_manual(name="",values = c("0"="#353F6B","1"="#f37735","2"="#d11141"))+
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  theme_classic()+theme(legend.position="none")

response<-data.frame("PC.5"=myscores[,,5]$spc, class=spc.low$class)
colnames(response)<-c("PC5", "class")
#ggboxplot(response, x = "class", y = "PC5")
res.kruskal <- response %>% kruskal_test(PC5 ~ class)
pwc <- response %>% 
  dunn_test(PC5 ~ class, p.adjust.method = "bonferroni") 
pwc
```

    ## # A tibble: 3 x 9
    ##   .y.   group1 group2    n1    n2 statistic       p   p.adj p.adj.signif
    ## * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
    ## 1 PC5   0      1         17    23     -1.47 0.141   0.422   ns          
    ## 2 PC5   0      2         17    27     -2.96 0.00305 0.00916 **          
    ## 3 PC5   1      2         23    27     -1.57 0.116   0.348   ns

``` r
pwc <- pwc %>% add_xy_position(x = "class")
p5<-ggboxplot(response, x = "class", y = "PC5", outlier.shape = NA) +
  geom_point(aes(x=class,y=PC5,color=class,fill=class), 
             position=position_jitterdodge(jitter.width=0.1,
                                           dodge.width=0.7),
             size=4, alpha=0.7)+
  scale_colour_manual(name="",values = c("0"="#353F6B","1"="#f37735","2"="#d11141"))+
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  theme_classic()+theme(legend.position="none")
p1
```

![](s4a_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->

``` r
p5
```

![](s4a_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->
