install.packages("quid")
library(quid)
data(stroop)
resStroop <- constraintBF(formula = rtS ~ ID*cond,
                          data = stroop,
                          whichRandom = "ID",
                          ID = "ID",
                          whichConstraint = c("cond" = "2 > 1"),
                          rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10))
bfs<-resStroop@generalTestObj
bfs[4]/bfs[3]
bfs/max(bfs)

plotEffects(resStroop)
plotEffects(resStroop, .raw = TRUE)


data(ld5)

resLD5 <- constraintBF(formula = rt ~ sub * distance + side, #rt=output, Modell ist Funktion aus HE sub, HE distance, Interk sub X distance und HE side
                       data = ld5,
                       whichRandom = c("sub"), #sub 0random Factor
                       ID = "sub", #VP ID ist speichert in sub
                       whichConstraint = c("distance" = "1 > 2", "distance" = "2 > 3"), #RZ in distance 2 größer als in 1 UND RT in Distance 2 größer als in 3
                       rscaleEffects = c("sub" = 1,
                                         "side" = 1/6,
                                         "distance" = 1/6,
                                         "sub:distance" = 1/10)) #setzt spezifische priors für fixed, random und interaction effects

plotEffects(resLD5)
