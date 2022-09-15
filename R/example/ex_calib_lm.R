# Copyright (c) 2022 Dongdong Kong and Yuxuan Xie. All rights reserved.
# 利用SCEUA优化参数
library(ggplot2)
library(hydroTools)


## 1. write your simulation function and prepare your inputs -------------------
# simulation function
# ' @author [Yuxuan Xie](xieyx@cug.edu.cn)
fun_predict <- function(x, par) {
  x * par[1] + x * par[2]
}

#' Goal function
#' 
#' Find the parameter corresponding to the minimum value of the objective function
#' 
#' @note we need the index, the smaller, the better.
goal <- function(par, x, yobs, ...) {
  ysim = fun_predict(x, par) # ysim with the par
  hydroTools::GOF(yobs, ysim)$RMSE # one of NSE, KGE, RMSE
}

set.seed(1) # Ensure the random numbers generated each time are consistent
n = 50
x = 1:n
par = c(0.5, 2) # real parameter
yobs = par[1]*x + par[2] + rnorm(n, 1) # 实际应用过程中，yobs是已经提供的

fun_predict(x, par) # test above function
goal(par, x, yobs)

## 2. parameter optimization ---------------------------------------------------
par_u = c(6, 6)   # upper boundary of par
par_l = c(-1, -1) # lower boundary of par
par0 = c(1, 1)    # initial parameter

# 这里以sceua为例展示如何进行参数优化，其他优化函数可以自行触类旁通
opt = rtop::sceua(goal, par = par0, lower = par_l, upper = par_u, 
  x = x, yobs = yobs, # other parameters passed to `goal`
  maxn = 1e3) 
par_opt = opt$par # SCEUA生成的最优参数

## 3. visualization ------------------------------------------------------------
ysim = fun_predict(x, par_opt)
dat = data.frame(x, yobs, ysim)

xx = seq(min(x), max(x), 0.1) # xsim
yy = fun_predict(xx, par = par_opt)

plot(x, yobs)
lines(xx, yy) # 两图叠加进行验证

GOF(dat$yobs, dat$ysim)

if (require(gg.layers)) {
  # remotes::install_github("rpkgs/gg.layers")
  p = ggplot(dat, aes(yobs, ysim)) + 
    geom_point() + 
    stat_gof(show.line = TRUE)
  Ipaper::write_fig(p, "Figure1_calib_lm.pdf", 6, 5)
}
