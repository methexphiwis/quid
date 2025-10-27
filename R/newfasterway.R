#new estimatePrioProbability function

install.packages("quid")
library(quid)

data(stroop)


# get constraints
whichConstraint <- c("cond" = "2 > 1")
constraints <- quid:::createConstraints(whichConstraint = whichConstraint)
cleanConstraints <- quid:::createCleanConstraints(constraints = constraints)
effectNameOrg <- unique(names(whichConstraint))

generalTestObj <- BayesFactor::generalTestBF(formula = rtS ~ ID*cond, data = stroop,
                                             whichRandom = NULL,
                                             rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),)
# get index of full model
indexFullModel <- quid:::extractIndexFullModel(generalTestObj)
# sample from full model posterior
iterationsPosterior = 10000
thetas <- BayesFactor::posterior(generalTestObj, index = indexFullModel, iterations = iterationsPosterior)
# clean names
colnames(thetas) <- cleanName(colnames(thetas))
IDorg <- ID
effectNameOrg <- unique(names(whichConstraint))
ID <- cleanName(ID)

iTheta <- quid:::extractIndeces(constraints = constraints,
                               thetas = thetas,
                               ID = thetas@data$ID,
                               data = stroop,
                               formula = rtS ~ ID*cond,
                               IDorg = IDorg)


prior_pass_vec <- quid:::estimatePriorProbability(iTheta = iTheta,
                                                 rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                                                 iterationsPrior = 1000,
                                                 cleanConstraints = cleanConstraints,
                                                 IDorg = IDorg,
                                                 effectNameOrg = effectNameOrg)(
  iTheta = iTheta,
  rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
  iterationsPrior = 1000,
  cleanConstraints = cleanConstraints,
  IDorg = stroop$ID,
  effectNameOrg = effectNameOrg)
