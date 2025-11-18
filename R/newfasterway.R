#new estimatePrioProbability function

install.packages("quid")
library(quid)

data(stroop)
#run script with getiTheta, estimatePriorProbability and fixedFromRandomProjection function, run script utils


iTheta <- get_iTheta(formula = rtS ~ ID*cond,
                       data = stroop,
                       whichRandom = "ID",
                       ID = "ID",
                       whichConstraint = c("cond" = "2 > 1"),
                       rscaleEffects= c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                       iterationsPosterior = 3,
                       iterationsPrior = 1000,
                       burnin = 1)

totalThetas <- addThetas(thetas = iTheta$thetas, iTheta = iTheta, keep =  (1 + 1) : 1000)

prior_pass_vec <- estimatePriorProbability(iTheta = iTheta,
                                           rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                                           iterationsPrior = 1000,
                                           cleanConstraints = iTheta$cleanConstraints,
                                           totalThetas = iTheta$totalThetas,
                                           IDorg = "ID",
                                           effectNameOrg = "cond")

#(
 # iTheta = iTheta,
  #rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
  #iterationsPrior = 1000,
  #cleanConstraints = cleanConstraints,
  #IDorg = stroop$ID,
  #effectNameOrg = effectNameOrg)

# generalTestObj <- BayesFactor::generalTestBF(formula = rtS ~ ID*cond, data = stroop,
#                                             whichRandom = NULL,
#                                             rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),)
# get index of full model
# indexFullModel <- quid:::extractIndexFullModel(generalTestObj)
# sample from full model posterior
# iterationsPosterior = 10000
# thetas <- BayesFactor::posterior(generalTestObj, index = indexFullModel, iterations = iterationsPosterior)
# clean names
# colnames(thetas) <- cleanName(colnames(thetas))
# IDorg <- ID
# effectNameOrg <- unique(names(whichConstraint))
# ID <- cleanName(ID)

# iTheta <- quid:::extractIndeces(constraints = constraints,
#                                thetas = thetas,
#                                ID = thetas@data$ID,
#                                data = stroop,
#                                formula = rtS ~ ID*cond,
#                                IDorg = IDorg)
