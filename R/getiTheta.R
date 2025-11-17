#' Function to compute Bayes factors for ordinal constraints
get_iTheta <- function(formula, data, whichRandom, ID,
                        whichConstraint, rscaleEffects,
                        iterationsPosterior = 3, iterationsPrior = 1000,
                        burnin = 1, ...) {

  # Assertions
  # data <- forceDataFrame(data)
  # data <- checkUsedLevels(formula = formula, data = data)
  # ellipsis::check_dots_used()
  # checkFormulaData(formula = formula, data = data)
  # checkID(ID = ID, data = data)
  # checkConstraints(whichConstraint = whichConstraint, data = data)
  # checkPriors(rscaleEffects = rscaleEffects, formula = formula, data = data, ID = ID, whichConstraint = whichConstraint)
  # checkIterations(iterationsPosterior = iterationsPosterior, burnin = burnin)

  #run only one nessecary model
  bf_full <- BayesFactor::lmBF(
    formula = formula,
    data = data,
    whichRandom = whichRandom,
    rscaleEffects = rscaleEffects,
    iterations = iterationsPosterior,
    progress = FALSE
  )

  # no need to run all models
  # generalTestObj <- BayesFactor::generalTestBF(formula = formula, data = data,
  #                                             whichRandom = whichRandom,
  #                                             rscaleEffects = rscaleEffects,
  #                                             ...)

  # get index of full model
  # indexFullModel <- extractIndexFullModel(generalTestObj)
  # sample from full model posterior
  thetas <- BayesFactor::posterior(bf_full, iterations = iterationsPosterior)

  # clean names
  # colnames(thetas) <- cleanName(colnames(thetas))
  IDorg <- ID
  effectNameOrg <- unique(names(whichConstraint))
  # ID <- cleanName(ID)

  # get constraints with quid
  constraints <- quid:::createConstraints(whichConstraint = whichConstraint)
  cleanConstraints <- quid:::createCleanConstraints(constraints = constraints)

  # more simple way
  # cleanConstraints <- data.frame(
  #  upper = "cond_2",
  #  lower = "cond_1",
  #  stringsAsFactors = FALSE
  #)

  # get indeces for posterior (start of extractIndeces function)
  effectName <- unique(constraints$constraintEffect)
  #effectName <- cleanName(effectName)

  # get all unique values of relevant factors
  #colnames(data) <- cleanName(colnames(data))

  effectLevels <- as.factor(sort(unique(c(constraints$constraintUpper, constraints$constraintLower))))
  IDLevels <- unique(do.call(`$`, args = list(x = data, name = ID)))

  # common effect
  regexTheta0 <- paste0("^", effectName, "_", effectLevels, "$")
  iTheta0 <- sapply(regexTheta0, function(pat) grep(pattern = pat, x = colnames(thetas))) #wird nicht richtig gespeichert
  names(iTheta0) <- paste0(effectName, "_", effectLevels)

  # definitions from crossRegex function
  trms <- attr(terms(formula), "term.labels")
  #trms <- cleanName(trms)
  print(trms)
  idFirst <- paste0(ID, ":", effectName)
  print(idFirst)
  effectFirst <- paste0(effectName, ":", ID)
  print(effectFirst)

  if (idFirst %in% trms) {
    suffix <- outer(IDLevels, effectLevels, FUN = paste, sep = "_")
    regexThetaID <- apply(suffix, MARGIN = c(1, 2), function(x) paste0("^", ID, "_", effectName, "_", x, "$"))
    regexThetaID <- as.data.frame(regexThetaID)
    colnames(regexThetaID) <- paste0(effectName, "_", effectLevels)
  } else if (effectFirst %in% trms) {
    suffix <- outer(effectLevels, IDLevels, FUN = paste, sep = "_")
    suffix <- t(suffix)
    regexThetaID <- apply(suffix, MARGIN = c(1, 2), function(x) paste0("^", effectName, "_", ID, "_", x, "$"))
    regexThetaID <- as.data.frame(regexThetaID)
    colnames(regexThetaID) <- paste0(effectName, "_", effectLevels)
  } else {
    stop("Unable to match formula elements with column names from posterior samples.")
  }

  # individual effects
  iThetaID <- apply(regexThetaID, MARGIN = c(1, 2), function(pat) grep(pattern = pat, x = colnames(thetas)))


  # create iTheta (end of extractIndeces function)
  iTheta <- list(commonEffect = iTheta0,
                 indEffect = iThetaID,
                 IDLevels = IDLevels,
                 effectLevels = effectLevels,
                 ID = IDorg,
                 cleanConstraints = cleanConstraints,
                 thetas = thetas)

 }

# quid::addThetas
addThetas <- function(thetas = thetas, iTheta = iTheta, keep = keep) {

  totalTheta <- vector(mode = "list", length = length(iTheta$commonEffect))
  names(totalTheta) <- colnames(iTheta$indEffect)

  for (i in seq_along(totalTheta)) {
    totalTheta[[i]] <- thetas[keep, iTheta$commonEffect[i]] + thetas[keep, iTheta$indEffect[, i]]
  }

  return(totalTheta)
}

# quid:::estimatePriorProbability
estimatePriorProbability <- function(iTheta = iTheta,
                                     rscaleEffects = rscaleEffects,
                                     iterationsPrior = iterationsPrior,
                                     cleanConstraints = cleanConstraints,
                                     IDorg = IDorg,
                                     effectNameOrg = effectNameOrg) {
  # get parameterization
  nLevels <- length(iTheta[["effectLevels"]]) # num conditions
  I <- length(iTheta[["IDLevels"]]) # num participants
  regexEffect <- paste0("^", effectNameOrg, ":", IDorg, "$", "|", "^", IDorg, ":", effectNameOrg, "$") # creats single string containing ID:effectname or effectname:ID
  indEffect <- grep(regexEffect, names(rscaleEffects)) #vector showing the position of the values in rscalesEffect which perfectly match with values in regexEffect -> used to extract relevant elements of rscaleEffets, e.g. names(rscaleEffects) = c("age:participant1", "gender:participant2", "participant1:age") and regexEffect = "^age:participant1$|^participant1:age$" then indEffect will be c(1, 3)

  # nLevels is projected onto n-1 levels by use of the following projection matrix
  params <- fixedFromRandomProjection(nLevels, sparse = FALSE) #done to meet statistical constraints

  # sample n main effects m times
  mus <- matrix(nrow = iterationsPrior, ncol = ncol(params)) #creats empty matrix with each row = iteration and each column = parameter from params
  for(i in 1:ncol(mus)) {
    mus[, i] <- rcauchy(iterationsPrior, 0, rscaleEffects[effectNameOrg])
  } #for each column i in mus draw sample (iterationsPrior) from cauchy distribution (location =0, r from effect name) -> use it later to simulate main effects in prior

  # sample SDs for individual effects
  gID <- MCMCpack::rinvgamma(iterationsPrior, .5, .5 * rscaleEffects[indEffect]) # MCMCpack::rinvgamma(n, shape, scale): number of random samples to draw = iterationsPrios, shape parameter is .5, scale parameter is 0.5*rscaleEffects[indEffects]

  # set up for loop to estimate how often effects are in the expected direction
  pass <- 1:iterationsPrior #creates vector from 1 to iterationsPrior

  # sample n individual effects m times
  for (m in 1:iterationsPrior){
    rEffects <- matrix(nrow = I, ncol = ncol(params)) #matrix with I = ID rows and and as many columns as params has

    for (n in 1:ncol(rEffects)) {
      rEffects[, n] <- rnorm(I, 0, sqrt(gID[m]))
    } # Simulate individual-level random effects for this prior sample

    # make empty list of lenght n
    totalPrior <- vector(mode = "list", length = nLevels)
    names(totalPrior) <- colnames(iTheta$indEffect)

    # multiply effects with parameterization
    for (i in seq_along(totalPrior)) { #Iterates over each element of totalPrior, which corresponds to each effect level
      totalPrior[[i]] <- sum(mus[m, ] * params[i, ]) +
        rEffects %*% params[i, ]
    } #total prior = get the total prior effect for each individual at this effect level

    # estimate how often effects are in the expected direction
    X <- vector(mode = "list", length = nrow(cleanConstraints))
    # build up expression and then evaluate
    for (i in 1:nrow(cleanConstraints)) {
      X[[i]] <- rlang::expr(totalThetas[[!!cleanConstraints[i, "upper"]]] > totalThetas[[!!cleanConstraints[i, "lower"]]])
    }
    Y <- purrr:::reduce(X, ~ rlang::expr(!!.x & !!.y))
    constrainedPrior <- eval(Y)
    # logical vector (e.g., TRUE/FALSE for each individual) indicating whether the simulated effects meet the constraints
    pass[m] <- sum(constrainedPrior) == I #counts true for each individual
  }
  return(pass)
}

# BayesFactor:::fixedFromRandomProjection
fixedFromRandomProjection <- function (nlevRandom, sparse = FALSE) {
  centering = diag(nlevRandom) - (1/nlevRandom)
  S = as.vector((eigen(centering)$vectors)[, 1:(nlevRandom - 1)])
  return(Matrix::Matrix(S, nrow = nlevRandom, sparse = sparse))
}

