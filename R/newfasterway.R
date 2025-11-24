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

totalThetas <- addThetas(thetas = iTheta$thetas, iTheta = iTheta, keep =  2 : 3)

prior_pass_vec <- estimatePriorProbability(iTheta = iTheta,
                                           rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                                           iterationsPrior = 1000,
                                           cleanConstraints = iTheta$cleanConstraints,
                                           IDorg = "ID",
                                           effectNameOrg = "cond")

priorProbability <- mean(prior_pass_vec)

# benchmarks
install.packages("remotes")
library(remotes)
remotes::install_github("joshuaulrich/microbenchmark")

pipeline <- function(x) {
  iTheta <- get_iTheta(formula = rtS ~ ID*cond,
                       data = stroop,
                       whichRandom = "ID",
                       ID = "ID",
                       whichConstraint = c("cond" = "2 > 1"),
                       rscaleEffects= c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                       iterationsPosterior = 3,
                       iterationsPrior = 1000,
                       burnin = 1)
  totalThetas <- addThetas(thetas = iTheta$thetas, iTheta = iTheta, keep =  2 : 3)
  prior_pass_vec <- estimatePriorProbability(iTheta = iTheta,
                                             rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                                             iterationsPrior = 1000,
                                             cleanConstraints = iTheta$cleanConstraints,
                                             IDorg = "ID",
                                             effectNameOrg = "cond")
  priorProbability <- mean(prior_pass_vec)

}

oldway <- function(x){
  quid::constraintBF(formula = rtS ~ ID*cond,
                    data = stroop,
                    whichRandom = "ID",
                    ID = "ID",
                    whichConstraint = c("cond" = "2 > 1"),
                    rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10))
}

res <- microbenchmark::microbenchmark(oldway(x),
                                      pipeline(x),
                                      times = 5,
                                      unit = "s")

print(res,
      unit = "s")

## Plot results:
boxplot(
  res,
  unit = "s",
  log = FALSE,
  xlab = "functions",
  horizontal = TRUE,
  main = "microbenchmark timings"
)
## Pretty plot:
if (requireNamespace("ggplot2")) {
  ggplot2::autoplot(res)
}

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
