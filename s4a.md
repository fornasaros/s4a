Scripts used in evaluating SERS spectra from gingival crevicular fluid
(GCF)
================

## Overview

These scripts were used to process spectral data using the
[hyperSpec](http://cbeleites.github.io/hyperSpec/) package.

Please see the vignettes for a walkthrough and more complete usage and
output information.

## First Glance at Data

``` r
load("~/Dropbox/Raman/In preparation/Drafts/denti/fitting/gcf.RData")
spc
```

    ## hyperSpec object
    ##    70 spectra
    ##    49 data columns
    ##    1999 data points / spectrum
    ## wavelength: paste("Raman shift (", cm^-1, ")", sep = "") [numeric] -184.37 -181.83 ... 3270.66 
    ## data:  (70 rows x 49 columns)
    ##    1. spc: Intensity (a.u.) [matrix, array1999] 17.59814 32.06528 ... 17.48172 
    ##    2. laser: laser [integer] 785 785 ... 785 
    ##    3. power: power [integer] 10 10 ... 10 
    ##    4. patient: patient [character] X10 X10 ... X7 
    ##    5. tooth: tooth [character] 11 13 ... 36 
    ##    6. sample: sample [character] X10_11 X10_13 ... X7_36 
    ##    7. .aggregate:  [factor] X10_11 X10_13 ... X7_36 
    ##    8. ID: ID [character] X10_11 X10_13 ... X7_36 
    ##    9. Batch:  [numeric] 1 1 ... 1 
    ##    10. Estrazione:  [numeric] 0 0 ... 0 
    ##    11. Conservazione:  [numeric] 0 0 ... 0 
    ##    12. Dente:  [numeric] 11 13 ... 36 
    ##    13. Tipo dente:  [numeric] 3 3 ... 1 
    ##    14. Lato:  [numeric] 1 1 ... 2 
    ##    15. Sup/Inf:  [numeric] 1 1 ... 2 
    ##    16. Genere:  [numeric] 1 1 ... 0 
    ##    17. Età:  [numeric] 57 57 ... 44 
    ##    18. Gravidanza:  [character] / / ... 0 
    ##    19. Fase mestruale:  [character] / / ... 21 
    ##    20. Fumo:  [numeric] 0 0 ... 0 
    ##    21. Patologie sistemiche:  [numeric] 0 0 ... 2 
    ##    22. Assunzione farmaci:  [numeric] 1 1 ... 2 
    ##    23. Sondaggio DV:  [numeric] 2 2 ... 6 
    ##    24. Placca DV:  [numeric] 1 1 ... 1 
    ##    25. BoP DV:  [numeric] 0 0 ... 1 
    ##    26. Sondaggio V:  [numeric] 2 2 ... 3 
    ##    27. Placca V:  [numeric] 0 1 ... 0 
    ##    28. BoP V:  [numeric] 0 0 ... 0 
    ##    29. Sondaggio MV:  [numeric] 3 3 ... 4 
    ##    30. Placca MV:  [numeric] 1 1 ... 1 
    ##    31. BoP MV:  [numeric] 0 0 ... 1 
    ##    32. Sondaggio DL:  [numeric] 2 3 ... 3 
    ##    33. Placca DL:  [numeric] 1 1 ... 1 
    ##    34. BoP DL:  [numeric] 1 1 ... 1 
    ##    35. Sondaggio L:  [numeric] 2 3 ... 3 
    ##    36. Placca L:  [numeric] 0 0 ... 1 
    ##    37. BoP L:  [numeric] 0 0 ... 1 
    ##    38. Sondaggio ML:  [numeric] 2 3 ... 4 
    ##    39. Placca ML:  [numeric] 1 1 ... 1 
    ##    40. BoP ML:  [numeric] 0 0 ... 1 
    ##    41. Diagnosi:  [character] Gingivitis    Periodontitis ... Periodontitis 
    ##    42. PSR 1:  [character] 4 4 ... 4 
    ##    43. PSR 2:  [character] 3 3 ... 3 
    ##    44. PSR 3:  [numeric] 4 4 ... 4 
    ##    45. PSR 4:  [character] 3 3 ... 4 
    ##    46. PSR 5:  [numeric] 3 3 ... 3 
    ##    47. PSR 6:  [numeric] 4 4 ... 4 
    ##    48. age: age [character] over 40 over 40 ... over 40 
    ##    49. class: class [numeric] 1 2 ... 2

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
plot(spc.low, col=as.factor(spc.low$Diagnosi))
```

![](s4a_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Spectral fitting

``` r
load("~/Dropbox/Raman/In preparation/Drafts/denti/fitting/metabolites.RData")

# preprocessing spcmet

spcmet <- spcmet[,,400~1800]
bl <- baseline (spcmet [[]],  method='als',lambda = 4, p = 0.001, maxit = 100) #method 2
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

    ##         c1         c2         c3         c4         c5 
    ## 0.27652654 0.23379942 0.41711803 0.06533604 0.22305921

``` r
spc.fit <- myfunc(coef(myfit)[1], coef(myfit)[2], coef(myfit)[3], 
                  coef(myfit)[4], coef(myfit)[5])

medres<- medspc - spc.fit


#calculate fitting statistics
norm_vec <- function(x) sqrt(sum(x^2))
predvals<-myfit$m$fitted()
lof<-sqrt((sum(medspc-predvals)^2)/(sum(medspc)^2))*100
rfe<- norm_vec(medspc-medres)/norm_vec(medspc)
r2<-1-((sum(medspc-predvals)^2)/(sum(medspc)^2))

data.frame(lof,rfe,r2)
```

    ##      lof       rfe        r2
    ## 1 5.3428 0.9597772 0.9971454

### Plot single components, sum and residuals

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

# PCA

``` r
myPCA <- PCA(scale(spc.low[[]], scale=F))
row.names(myPCA$scores)<-(spc.low$ID)
scoreplot(myPCA, pc= c(1,4),show.names = T,pch=19)
```

![](s4a_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
myloadings <- decomposition (spc.low[], t(myPCA$loadings), scores = FALSE,
                             label.wavelength = "PC", label.spc = "loading I / a.u.")
myscores <- decomposition (spc.low[], myPCA$scores, label.wavelength = "PC",
                           label.spc = "score / a.u.")

#screeplot(myPCA, type = "percentage", 15, main = "Explained variance")

plot(myloadings[c(1:6)], stacked=T)
```

![](s4a_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
response<-data.frame("PC.1"=myscores[,,1]$spc, class=spc.low$class)
colnames(response)<-c("PC", "class")
#ggboxplot(response, x = "class", y = "PC")
res.kruskal <- response %>% kruskal_test(PC ~ class)
pwc <- response %>% 
  dunn_test(PC ~ class, p.adjust.method = "bonferroni") 
pwc
```

    ## # A tibble: 3 x 9
    ##   .y.   group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
    ## * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
    ## 1 PC    0      1         19    24  -2.75    0.00602 0.0180 *           
    ## 2 PC    0      2         19    27  -2.81    0.00489 0.0147 *           
    ## 3 PC    1      2         24    27   0.00324 0.997   1      ns

``` r
pwc <- pwc %>% add_xy_position(x = "class")
p1<-ggboxplot(response, x = "class", y = "PC") +
  geom_point(aes(x=class,y=PC,color=class,fill=class), 
             position=position_jitterdodge(jitter.width=0.1,
                                           dodge.width=0.7),
             size=4, alpha=0.7)+
  scale_colour_manual(name="",values = c("0"="#353F6B","1"="#f37735","2"="#d11141"))+
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )+
  theme_classic()+theme(legend.position="none")

response<-data.frame("PC.5"=myscores[,,5]$spc, class=spc.low$class)
colnames(response)<-c("PC", "class")
#ggboxplot(response, x = "class", y = "PC")
res.kruskal <- response %>% kruskal_test(PC ~ class)
pwc <- response %>% 
  dunn_test(PC ~ class, p.adjust.method = "bonferroni") 
pwc
```

    ## # A tibble: 3 x 9
    ##   .y.   group1 group2    n1    n2 statistic      p  p.adj p.adj.signif
    ## * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl>  <dbl> <chr>       
    ## 1 PC    0      1         19    24     0.631 0.528  1      ns          
    ## 2 PC    0      2         19    27     2.32  0.0201 0.0602 ns          
    ## 3 PC    1      2         24    27     1.79  0.0733 0.220  ns

``` r
pwc <- pwc %>% add_xy_position(x = "class")
p5<-ggboxplot(response, x = "class", y = "PC") +
  geom_point(aes(x=class,y=PC,color=class,fill=class), 
             position=position_jitterdodge(jitter.width=0.1,
                                           dodge.width=0.7),
             size=4, alpha=0.7)+
  scale_colour_manual(name="",values = c("0"="#353F6B","1"="#f37735","2"="#d11141"))+
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )+
  theme_classic()+theme(legend.position="none")
p1
```

![](s4a_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
p5
```

![](s4a_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->
