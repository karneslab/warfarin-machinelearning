library("HardyWeinberg")
table(latinosdf$vkor)
x <- c(GG = 729, AG = 778, AA = 217)
HW.test <- HWChisq(x, verbose = TRUE)

table(latinosdf$cyp)
