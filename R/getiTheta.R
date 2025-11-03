#' Function to compute Bayes factors for ordinal constraints
#'
#' This function uses Bayesian mixed models to estimate individual effect
#' sizes and to test theoretical order constraints.
#'
#' This function provides a way of testing whether theoretical constraints on
#' certain effects hold for all subjects. The backend is provided by the
#' \code{\link[BayesFactor]{generalTestBF}} function from the
#' \code{\link[BayesFactor]{BayesFactor-package}}. The input formula is the
#' full model to be tested. It usually contains an interaction term between
#' the subject ID and the effect for which constraints are tested (e.g.
#' \code{ID:condition}). The ID variable is to be specified in \code{ID} and is
#' usually a random factor to be specified in \code{whichRandom}.
#'
#' Order constraints on effects should be specified in \code{whichConstraint},
#' as a named character vector. Each constraint in the vector can take 2 levels
#' of the effect. They are of the form:
#' \code{"effect name" = "condition A" < "condition B"}. In order to impute more
#' than 2 levels, the same effect name has to be entered with different conditions
#' as the value. For instance, for testing whether conditions A < B < C, the
#' input should be: \code{"effect name" = "condition A" < "condition B", "effect name" = "condition B" < "condition C"}.
#' At this point, constraints can only be tested for the same effect.
#'
#' Priors have to be specified for all factors in \code{whichConstraint},
#' for \code{ID}, and for the interaction between the two. A Detailed description
#' of the models, priors and methods is given in the documentation of
#' \code{\link[BayesFactor]{anovaBF}} and more extensively in Rouder et al. (2012).
#'
#' @param formula a formula containing the full model.
#' @param data a \code{data.frame} containing the data with all variables
#'   defined in the formula.
#' @param whichRandom a character vector specifying which factors are random.
#' @param ID a character vector of length one specifying which variable
#'   holds the subject ID.
#' @param whichConstraint a named character vector specifying the constraints
#'   placed on certain factors; see Details.
#' @param rscaleEffects a named vector of prior settings for individual factors.
#'   Values are scales, names are factor names; see Details.
#' @param iterationsPosterior the number of iterations to sample from the
#'  posterior of the full model.
#' @param iterationsPrior the number of iterations to sample from the
#'  prior of the full model.
#' @param burnin the number of initial iterations to discard from posterior
#'   sampling.
#' @param ... further arguments to be passed to
#'   \code{\link[BayesFactor]{generalTestBF}}.
#'
#' @return An object of class \code{\link{BFBayesFactorConstraint-class}}.
#'
#' @references Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012)
##'   Default Bayes Factors for ANOVA Designs. Journal of Mathematical
##'   Psychology.  56.  p. 356-374.
#'
#' @importFrom stats as.formula rcauchy rnorm sd terms
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(stroop)
#'
#' resStroop <- constraintBF(rtS ~ ID*cond,
#'                           data = stroop,
#'                           whichRandom = "ID",
#'                           ID = "ID",
#'                           whichConstraint = c(cond = "2 > 1"),
#'                           rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10))
#' }
#'
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
  effectName <- cleanName(effectName)

  # get all unique values of relevant factors
  colnames(data) <- cleanName(colnames(data))

  effectLevels <- as.factor(sort(unique(c(constraints$constraintUpper, constraints$constraintLower))))
  IDLevels <- unique(do.call(`$`, args = list(x = data, name = ID)))

  # common effect
  regexTheta0 <- paste0("^", effectName, "_", effectLevels, "$")
  iTheta0 <- sapply(regexTheta0, function(pat) grep(pattern = pat, x = colnames(thetas)))
  names(iTheta0) <- paste0(effectName, "_", effectLevels)

  # definitions from crossRegex function
  trms <- attr(terms(formula), "term.labels")
  trms <- cleanName(trms)
  print(trms)
  idFirst <- paste0(ID, "_", effectName)
  print(idFirst)
  effectFirst <- paste0(effectName, "_", ID)
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
                 ID = IDorg)

  # add overall effect to individuals' deviations
  #keep <- (burnin + 1) : iterationsPosterior

  #otalThetas <- addThetas(thetas = thetas, iTheta = iTheta, keep = keep)

  # evaluate posterior probability of all thetas being positive
  # constrainedThetas <- quid:::estimateConstrainedThetas(totalThetas = totalThetas, cleanConstraints = cleanConstraints)
  # passThetas <- apply(constrainedThetas, 1, mean) == 1
  # posteriorProbability <- mean(passThetas)

  # # get prior probability of all thetas being positive
  # passPrior <- estimatePriorProbability(iTheta = iTheta,
  #                                       rscaleEffects = rscaleEffects,
  #                                       iterationsPrior = iterationsPrior,
  #                                       cleanConstraints = cleanConstraints,
  #                                       IDorg = IDorg,
  #                                       effectNameOrg = effectNameOrg)
  #
  # priorProbability <- mean(passPrior)
  #
  # # prepare return values
  # bfCU <- posteriorProbability / priorProbability
  # individualEffects <- lapply(totalThetas, colMeans)
  # posteriorSD <- sapply(individualEffects, sd)
  # posteriorMean <- colMeans(thetas[keep, iTheta$commonEffect])
  # observedEffects <- calculateObservedEffects(constraints = constraints,
  #                                             data = data,
  #                                             IDorg = IDorg,
  #                                             iTheta = iTheta,
  #                                             formula = formula,
  #                                             effectNameOrg = effectNameOrg)
  #
  # # make S4 objects
  # newConstraint <- BFConstraint(priorProbability = priorProbability,
  #                               posteriorProbability = posteriorProbability,
  #                               bayesFactor = bfCU,
  #                               constraints = constraints,
  #                               cleanConstraints = cleanConstraints)
  #
  #
  # newBFConstraint <- BFBayesFactorConstraint(generalTestObj = generalTestObj,
  #                                            constraints = newConstraint,
  #                                            individualEffects = individualEffects,
  #                                            posteriorMean = posteriorMean,
  #                                            posteriorSD = posteriorSD,
  #                                            totalThetas = totalThetas,
  #                                            mcmcFull = thetas,
  #                                            designIndeces = iTheta,
  #                                            observedEffects = observedEffects)
  #
  # return(newBFConstraint)
}
