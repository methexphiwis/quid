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
    constrainedPrior <- estimateConstrainedThetas(totalThetas = totalPrior, cleanConstraints = cleanConstraints) # logical vector (e.g., TRUE/FALSE for each individual) indicating whether the simulated effects meet the constraints
    pass[m] <- sum(constrainedPrior) == I #counts true for each individual
  }
  return(pass)
}

# BayesFactor:::fixedFromRandomProjection
#' @importFrom Matrix Matrix
#'
fixedFromRandomProjection <- function (nlevRandom, sparse = FALSE) {
  centering = diag(nlevRandom) - (1/nlevRandom)
  S = as.vector((eigen(centering)$vectors)[, 1:(nlevRandom - 1)])
  return(Matrix::Matrix(S, nrow = nlevRandom, sparse = sparse))
}

