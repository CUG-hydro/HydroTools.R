# HydroTools.R

<!-- badges: start -->
[![R-CMD-check](https://github.com/CUG-hydro/HydroTools.R/workflows/R-CMD-check/badge.svg)](https://github.com/CUG-hydro/HydroTools.R/actions)
[![codecov](https://codecov.io/gh/CUG-hydro/HydroTools.R/branch/master/graph/badge.svg)](https://codecov.io/gh/CUG-hydro/HydroTools.R)
[![CRAN](http://www.r-pkg.org/badges/version/HydroTools.R)](https://cran.r-project.org/package=HydroTools.R)
<!-- badges: end -->

The goal of HydroTools.R is to ...

## Installation

You can install the released version of HydroTools.R from github:

``` r
# remotes::install_gitlab("r-pkgs/HydroTools.R") # github may not work, use gitlab instead
remotes::install_github("CUG-hydro/HydroTools.R") 
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(hydroTools)
## basic example code
```

## Finished

- [x] 蒸散发的基本原理，多种潜在蒸散发的计算方法
- [x] 实际蒸发：budyko理论
- [x] 实际蒸发：蒸发互补理论(CR)
- [x] 辐射与地表能量平衡
- [x] XAJ水文模型
- [x] 河网汇流：马斯京根
- [x] 土壤函数库，主要用于制作水文模型的输入数据
- [x] 气象变量处理函数，RH, q, w的转化, 湿球温度、干球温度、饱和水气压、空气密度
- [ ] 水汽通量、水汽通量散度
- [ ] 位势高度制图

> 欢迎同道中人贡献新的函数或方法
