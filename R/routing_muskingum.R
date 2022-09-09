# 包为民，水文预报（第4版）
routing_muskingum <- function(I1, I2, Q1, Q2, par) {
  Δt = par$Δt
  X = par$X
  K = par$K
  
  c0 = (-2*X + Δt/K) / ( 2*(1-X) + Δt/K)
  # c1 = (2*X + Δt/K) / ( 2*(1-X) + Δt/K)
  c2 = (2*(1-X) - Δt/K) / (2*(1-X) + Δt/K)
  # c0 + c1 + c2 = 1
  # c1 = 1 - (c0 + c2)
  Q2 = c0*I2 + c1*I1 + c2*Q1
  Q2
}

# 首先通过上下游最大滞后相关系数，确定滞后时间；
# 之后采用试算法（或优化方法），确定模型系数
