
[![CRAN
status](https://www.r-pkg.org/badges/version/predicts)](https://cran.r-project.org/package=predicts)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/predicts)](http://www.r-pkg.org/pkg/predicts)

`predicts` is an *R* package for spatial predictive modeling, especially species distribution modeling. There are tutorials at [rspatial.org/sdm](https://rspatial.org/sdm/index.html). 

`predicts` uses the "terra" package for spatial data handling replaces the [dismo](https://github.com/rspatial/dismo) package that uses the old packages "raster" and "sp". 


## Installation

`predicts` is available from CRAN, so you can use `install.packages("predicts")` to get the current *released version*.

The easiest way to use the *development version* on Windows or MacOS, is to install it from the [R-universe](https://r-universe.dev/organizations/), like this:


```
install.packages('predicts', repos='https://rspatial.r-universe.dev')
```


### From source-code

To install from source-code you can use: 
```
remotes::install_github("rspatial/predicts")
```
