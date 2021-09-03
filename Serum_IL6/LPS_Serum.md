---
output:
  html_document: 
    df_print: default
    highlight: haddock
    keep_md: yes
    theme: paper
editor_options: 
  chunk_output_type: console
---
# **ELISA data** for serum IL-6 after LPS or SAL treatment
Data analysis of serum protein expression of IL-6 after treatment with saline (SAL) or 4500 EU/g lipopolysaccharide (LPS) for 18 hours. Data are pg/mL.

## Samples

```r
library(ggplot2)
library(cowplot)
library(lme4)
library(lmerTest)
library(multcomp)
library(tidyverse)
```

```r
data <- read.csv("210309.csv", header=TRUE, stringsAsFactors = TRUE)
data$TRT <- relevel(data$TRT, ref="SAL")
```

## Plot data
![](Fig_1_Serum_files/figure-html/unnamed-chunk-3-1.png)<!-- -->![](Fig_1_Serum_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

## Linear Mixed Effect Model MIA vs CON

```
##           CON:SAL CON:LPS MIA:SAL MIA:LPS
## # mice    37      52      45      68     
## # litters 26      34      28      43
```

```
##         Shapiro Test for Normality
## CON:SAL                      0.001
## CON:LPS                      0.000
## MIA:SAL                      0.002
## MIA:LPS                      0.000
```

```r
lme.fit <- lmer(log2(Serum.IL.6+1) ~ Group*TRT+(1|Litter), data=MvC)
summary(lme.fit)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: log2(Serum.IL.6 + 1) ~ Group * TRT + (1 | Litter)
##    Data: MvC
## 
## REML criterion at convergence: 846.7
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.5590 -0.7317 -0.1818  0.5063  2.5996 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  Litter   (Intercept) 0.7935   0.8908  
##  Residual             3.4080   1.8461  
## Number of obs: 200, groups:  Litter, 86
## 
## Fixed effects:
##                 Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)       2.2030     0.3620 156.2463   6.086 8.59e-09 ***
## GroupMIA          0.3506     0.4858 150.8619   0.722    0.472    
## TRTLPS            7.2278     0.4272 170.7738  16.919  < 2e-16 ***
## GroupMIA:TRTLPS  -0.5081     0.5726 176.4439  -0.887    0.376    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) GrpMIA TRTLPS
## GroupMIA    -0.745              
## TRTLPS      -0.717  0.535       
## GMIA:TRTLPS  0.535 -0.721 -0.746
```

## Linear Mixed Effect Model CON:CON, CON:PLX, MIA:CON, MIA:PLX all after LPS

```
##         SAL LPS # litters
## CON:CON 0   21  9        
## CON:PLX 0   30  15       
## MIA:CON 0   15  7        
## MIA:PLX 0   34  15
```

```
##         Shapiro Test for Normality
## CON:CON                      0.327
## CON:PLX                      0.034
## MIA:CON                      0.193
## MIA:PLX                      0.051
```

```r
lme.fit <- lmer(log2(Serum.IL.6+1) ~ Group+(1|Litter), data=Rescue)
summary(lme.fit)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: log2(Serum.IL.6 + 1) ~ Group + (1 | Litter)
##    Data: Rescue
## 
## REML criterion at convergence: 422.4
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.2193 -0.5886 -0.1442  0.5647  2.7495 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  Litter   (Intercept) 2.760    1.661   
##  Residual             2.609    1.615   
## Number of obs: 100, groups:  Litter, 44
## 
## Fixed effects:
##              Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  11.05676    0.67894 50.92061  16.285   <2e-16 ***
## GroupCON:PLX -0.63497    0.82633 60.38291  -0.768    0.445    
## GroupMIA:CON  0.01167    1.04422 47.63427   0.011    0.991    
## GroupMIA:PLX -0.60200    0.85880 49.57053  -0.701    0.487    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) GCON:P GMIA:C
## GropCON:PLX -0.779              
## GropMIA:CON -0.650  0.506       
## GropMIA:PLX -0.789  0.628  0.513
```

### Reproducibility info

```
## [1] "2021-03-19 11:16:47 EDT"
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 4.0.3 (2020-10-10)
##  os       macOS Big Sur 10.16         
##  system   x86_64, darwin17.0          
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       America/New_York            
##  date     2021-03-19                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package     * version    date       lib source        
##  assertthat    0.2.1      2019-03-21 [1] CRAN (R 4.0.2)
##  backports     1.2.1      2020-12-09 [1] CRAN (R 4.0.2)
##  boot          1.3-27     2021-02-12 [1] CRAN (R 4.0.2)
##  broom         0.7.5      2021-02-19 [1] CRAN (R 4.0.3)
##  cachem        1.0.4      2021-02-13 [1] CRAN (R 4.0.2)
##  callr         3.5.1      2020-10-13 [1] CRAN (R 4.0.2)
##  cellranger    1.1.0      2016-07-27 [1] CRAN (R 4.0.2)
##  cli           2.3.1      2021-02-23 [1] CRAN (R 4.0.3)
##  codetools     0.2-18     2020-11-04 [1] CRAN (R 4.0.2)
##  colorspace    2.0-0      2020-11-11 [1] CRAN (R 4.0.2)
##  cowplot     * 1.1.1      2020-12-30 [1] CRAN (R 4.0.2)
##  crayon        1.4.1      2021-02-08 [1] CRAN (R 4.0.2)
##  DBI           1.1.1      2021-01-15 [1] CRAN (R 4.0.2)
##  dbplyr        2.1.0      2021-02-03 [1] CRAN (R 4.0.3)
##  desc          1.3.0      2021-03-05 [1] CRAN (R 4.0.2)
##  devtools      2.3.2      2020-09-18 [1] CRAN (R 4.0.2)
##  digest        0.6.27     2020-10-24 [1] CRAN (R 4.0.2)
##  dplyr       * 1.0.5      2021-03-05 [1] CRAN (R 4.0.2)
##  ellipsis      0.3.1      2020-05-15 [1] CRAN (R 4.0.2)
##  evaluate      0.14       2019-05-28 [1] CRAN (R 4.0.1)
##  fansi         0.4.2      2021-01-15 [1] CRAN (R 4.0.2)
##  farver        2.1.0      2021-02-28 [1] CRAN (R 4.0.2)
##  fastmap       1.1.0      2021-01-25 [1] CRAN (R 4.0.2)
##  forcats     * 0.5.1      2021-01-27 [1] CRAN (R 4.0.3)
##  fs            1.5.0      2020-07-31 [1] CRAN (R 4.0.2)
##  generics      0.1.0      2020-10-31 [1] CRAN (R 4.0.2)
##  ggplot2     * 3.3.3      2020-12-30 [1] CRAN (R 4.0.2)
##  glue          1.4.2      2020-08-27 [1] CRAN (R 4.0.2)
##  gtable        0.3.0      2019-03-25 [1] CRAN (R 4.0.2)
##  haven         2.3.1      2020-06-01 [1] CRAN (R 4.0.2)
##  highr         0.8        2019-03-20 [1] CRAN (R 4.0.2)
##  hms           1.0.0      2021-01-13 [1] CRAN (R 4.0.2)
##  htmltools     0.5.1.1    2021-01-22 [1] CRAN (R 4.0.2)
##  httr          1.4.2      2020-07-20 [1] CRAN (R 4.0.2)
##  jsonlite      1.7.2      2020-12-09 [1] CRAN (R 4.0.2)
##  knitr         1.31       2021-01-27 [1] CRAN (R 4.0.3)
##  labeling      0.4.2      2020-10-20 [1] CRAN (R 4.0.2)
##  lattice       0.20-41    2020-04-02 [1] CRAN (R 4.0.3)
##  lifecycle     1.0.0      2021-02-15 [1] CRAN (R 4.0.2)
##  lme4        * 1.1-26     2020-12-01 [1] CRAN (R 4.0.2)
##  lmerTest    * 3.1-3      2020-10-23 [1] CRAN (R 4.0.2)
##  lubridate     1.7.10     2021-02-26 [1] CRAN (R 4.0.2)
##  magrittr      2.0.1      2020-11-17 [1] CRAN (R 4.0.2)
##  MASS        * 7.3-53.1   2021-02-12 [1] CRAN (R 4.0.2)
##  Matrix      * 1.3-2      2021-01-06 [1] CRAN (R 4.0.2)
##  memoise       2.0.0      2021-01-26 [1] CRAN (R 4.0.2)
##  minqa         1.2.4      2014-10-09 [1] CRAN (R 4.0.2)
##  modelr        0.1.8      2020-05-19 [1] CRAN (R 4.0.2)
##  multcomp    * 1.4-16     2021-02-08 [1] CRAN (R 4.0.2)
##  munsell       0.5.0      2018-06-12 [1] CRAN (R 4.0.2)
##  mvtnorm     * 1.1-1      2020-06-09 [1] CRAN (R 4.0.2)
##  nlme          3.1-152    2021-02-04 [1] CRAN (R 4.0.2)
##  nloptr        1.2.2.2    2020-07-02 [1] CRAN (R 4.0.2)
##  numDeriv      2016.8-1.1 2019-06-06 [1] CRAN (R 4.0.2)
##  pillar        1.5.1      2021-03-05 [1] CRAN (R 4.0.2)
##  pkgbuild      1.2.0      2020-12-15 [1] CRAN (R 4.0.3)
##  pkgconfig     2.0.3      2019-09-22 [1] CRAN (R 4.0.2)
##  pkgload       1.2.0      2021-02-23 [1] CRAN (R 4.0.3)
##  prettyunits   1.1.1      2020-01-24 [1] CRAN (R 4.0.2)
##  processx      3.4.5      2020-11-30 [1] CRAN (R 4.0.2)
##  ps            1.6.0      2021-02-28 [1] CRAN (R 4.0.2)
##  purrr       * 0.3.4      2020-04-17 [1] CRAN (R 4.0.2)
##  R6            2.5.0      2020-10-28 [1] CRAN (R 4.0.2)
##  Rcpp          1.0.6      2021-01-15 [1] CRAN (R 4.0.2)
##  readr       * 1.4.0      2020-10-05 [1] CRAN (R 4.0.2)
##  readxl        1.3.1      2019-03-13 [1] CRAN (R 4.0.2)
##  remotes       2.2.0      2020-07-21 [1] CRAN (R 4.0.2)
##  reprex        1.0.0      2021-01-27 [1] CRAN (R 4.0.2)
##  rlang         0.4.10     2020-12-30 [1] CRAN (R 4.0.2)
##  rmarkdown     2.7        2021-02-19 [1] CRAN (R 4.0.3)
##  rprojroot     2.0.2      2020-11-15 [1] CRAN (R 4.0.2)
##  rstudioapi    0.13       2020-11-12 [1] CRAN (R 4.0.2)
##  rvest         1.0.0      2021-03-09 [1] CRAN (R 4.0.3)
##  sandwich      3.0-0      2020-10-02 [1] CRAN (R 4.0.2)
##  scales        1.1.1      2020-05-11 [1] CRAN (R 4.0.2)
##  sessioninfo   1.1.1      2018-11-05 [1] CRAN (R 4.0.2)
##  statmod       1.4.35     2020-10-19 [1] CRAN (R 4.0.2)
##  stringi       1.5.3      2020-09-09 [1] CRAN (R 4.0.2)
##  stringr     * 1.4.0      2019-02-10 [1] CRAN (R 4.0.2)
##  survival    * 3.2-7      2020-09-28 [1] CRAN (R 4.0.3)
##  testthat      3.0.2      2021-02-14 [1] CRAN (R 4.0.2)
##  TH.data     * 1.0-10     2019-01-21 [1] CRAN (R 4.0.2)
##  tibble      * 3.1.0      2021-02-25 [1] CRAN (R 4.0.3)
##  tidyr       * 1.1.3      2021-03-03 [1] CRAN (R 4.0.2)
##  tidyselect    1.1.0      2020-05-11 [1] CRAN (R 4.0.2)
##  tidyverse   * 1.3.0      2019-11-21 [1] CRAN (R 4.0.2)
##  usethis       2.0.1      2021-02-10 [1] CRAN (R 4.0.2)
##  utf8          1.2.1      2021-03-12 [1] CRAN (R 4.0.3)
##  vctrs         0.3.6      2020-12-17 [1] CRAN (R 4.0.2)
##  withr         2.4.1      2021-01-26 [1] CRAN (R 4.0.2)
##  xfun          0.22       2021-03-11 [1] CRAN (R 4.0.2)
##  xml2          1.3.2      2020-04-23 [1] CRAN (R 4.0.2)
##  yaml          2.2.1      2020-02-01 [1] CRAN (R 4.0.2)
##  zoo           1.8-9      2021-03-09 [1] CRAN (R 4.0.3)
## 
## [1] /Library/Frameworks/R.framework/Versions/4.0/Resources/library
```
