# ' @references 包为民，水文预报（第4版）
#' routing_muskingum
#' @author Dongdong Kong
#' @export 
routing_muskingum <- function(I1, I2, Q1, deltaT, par) {
  X = par[1] #$X
  K = par[2] #$K
  # deltaT = par$deltaT
  
  c0 = (-2*X + deltaT/K) / ( 2*(1-X) + deltaT/K)
  # c1 = (2*X + deltaT/K) / ( 2*(1-X) + deltaT/K)
  c2 = (2*(1-X) - deltaT/K) / (2*(1-X) + deltaT/K)
  # c0 + c1 + c2 = 1
  # c1 = 1 - (c0 + c2)
  Q2 = c0*I2 + c1*I1 + c2*Q1
  Q2
}

#' @param deltaT in hours
#' @param par 
#' - `X`  : 0-1/2
#' - `K`  : in the same unit as `deltaT`
#' - `nt` : divide river into `nt` chunks (default 1)
#' 
#' @rdname routing_muskingum
#' @export 
routing_muskingum_nL <- function(I, Q0, deltaT, par) {
  nt = par[3]
  N = length(I)

  # 需要一段稳定流输入，warming up period
  Q = rep(Q0, N)
  for (i in 1:(N - nt - 1)) {
    Q1 = Q[i]
    I1 = I[i]
    I2 = I[i + 1]

    for (j in 1:nt) {
      Q2 = routing_muskingum(I1, I2, Q1, deltaT, par)
      # for the next loop
      I1 = I[i + j]
      I2 = I[i + 1 + j]
      Q1 = Q2
    }
    Q[i + nt] = Q2
  }
  Q
}

#' @inheritParams XAJ_calib
#' @rdname routing_muskingum
#' @export
calib_routing_muskingum_nL <- function(I, Q0, yobs, maxn = 1e3) {
  goal <- function(par, I, Q0, deltaT, yobs) {
    ysim = routing_muskingum_nL(I, Q0, deltaT, par)
    -hydroTools::GOF(yobs, ysim)$KGE # one of -NSE, -KGE, RMSE
  }

  par_l <- c(0, 6, 1) # lower boundary of par
  par_u <- c(0.5, 48, 10) # upper boundary of par
  par0 <- c(0.5, 12, 1) # initial parameter

  # 这里以sceua为例展示如何进行参数优化，其他优化函数可以自行触类旁通
  opt <- rtop::sceua(goal,
    par = par0, lower = par_l, upper = par_u,
    I = I, Q0 = Q0, deltaT = deltaT, yobs = yobs, 
    maxn = maxn
  )
  par <- opt$par # SCEUA生成的最优参数
  ysim = routing_muskingum_nL(I, Q0, deltaT, par)
  predict = data.table(yobs, ysim)
  gof = hydroTools::GOF(yobs, ysim)

  listk(par, gof, predict)
}

# TODO: 
# - 1. test those functions
# - 2. 与`randomForest`进行对比，`routing_muskingum`不见得会比`randomForest`强
