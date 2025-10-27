install.packages("quid")
library(quid)
install.packages("microbenchmark")
library(microbenchmark)

data(stroop)
resStroop <- quid::constraintBF(formula = rtS ~ ID*cond,
                          data = stroop,
                          whichRandom = "ID",
                          ID = "ID",
                          whichConstraint = c("cond" = "2 > 1"),
                          rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10))
resStroop

bfs<-resStroop@generalTestObj
bfs[4]/bfs[3]
bfs/max(bfs)

plotEffects(resStroop)
plotEffects(resStroop, .raw = TRUE)

# call the function (use package namespace if needed)
prior_pass_vec <- estimatePriorProbability(
  iTheta = resStroop@designIndeces, #sollte jetzt richtig sein, wird darin gespeichert
  rscaleEffects = resStroop@generalTestObj@numerator$`ID + cond + ID:cond`@prior$rscale$effects,#verstehe die Art wie es gespeichert wird nicht
  iterationsPrior = iterationsPrior, #wird definiert aber nicht gespeichert, weiß nicht wie ich es ergänzen kann
  cleanConstraints = resStroop@constraints@cleanConstraints,
  IDorg = resStroop@generalTestObj@data$ID, #wird in observedEffects gespeichert
  effectNameOrg = resStroop@generalTestObj@data$cond) #wird in observedEffects gespeichert


##trying around

isS4(print)
showMethods(show)
# choose number of prior draws
iterationsPrior <- 2000

resStroop@totalThetas
str(resStroop)
??estimatePriorProbability
print(length(resStroop@generalTestObj@numerator$ID@prior$rscale$effects))
print(resStroop@generalTestObj@data$cond)

# Prüfe iTheta / totalThetas
str(resStroop@totalThetas)
print(length(resStroop@totalThetas[["effectLevels"]]))
print(names(resStroop@totalThetas))        # um zu sehen, welche Felder vorhanden sind

# Prüfe effectNameOrg
str(resStroop@designIndeces$indEffect)
print(length(resStroop@designIndeces$indEffect))
print(resStroop@designIndeces$indEffect)   # ist das "cond" oder ein Vektor?

# Prüfe rscaleEffects und seine Namen
rsc <- resStroop@generalTestObj@numerator$`ID + cond + ID:cond`@prior$rscale
str(rsc)


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

