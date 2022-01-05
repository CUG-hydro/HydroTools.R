library(dplyr)
# The relationship of Rn, Ts
Rn = 0:200
dat = cal_Ts(Rn, 25, D = 1, wind = 2, z.wind = 2)

plot(Rn, dat$Ts, type = "l", main = "(a) Ts ~ Rn")
plot(Rn, dat$Eeq, type = "l", main = "(b) Eeq ~ Rn")
# plot(Rn, dat$Evp, type = "l") # a constant value
dat %<>% mutate(Rn = Rn, bowen = ET0 / (Rn * 0.086400 / lambda - ET0))
plot(bowen ~ Rn, dat, type = "l", main = "(b) Eeq ~ Rn")
