---
title: Step by step guide
output: 
  rmarkdown::html_document:
    fig_width: 7
    fig_height: 6
    fig_caption: yes
    keep_md: yes
    self_contained: yes
    theme: readable
---


This workflow should show full strength of *RRatepol package* and serve as step by 
step guidance starting from downloading dataset from Neotoma, building age-depth
models, to estimating rate-of-change using age uncertainty.

:warning: **This workflow is only meant as example**: There are several additional steps for data reparation which should be done to really use the data from Neotoma!

## Install packages


Make a list of packages needed to from CRAN


```r
package_list <-
  c(
    "tidyverse", # general data wrangling and visualisation
    "pander", # nice tables
    "Bchron", # age-depth modeling
    "janitor", # string cleaning
    "remotes" # installing packages from GitHub
  )
```

Install all packages from CRAN using `{renv}` package


```r
lapply(
  package_list, renv::use
)
```

Install packages from GitHub


```r
# Install R-Ratepol
remotes::install_github("HOPE-UIB-BIO/R-Ratepol-package")

# Install neotoma2
remotes::install_github("NeotomaDB/neotoma2")
```

## Attach packages


```r
library(tidyverse) # general data wrangling and visualisation
library(pander) # nice tables
library(RRatepol) # rate-of-vegetation change
library(neotoma2) # obtain data from Neotoma database
library(Bchron) # age-depth modeling
library(janitor) # string cleaning
```


