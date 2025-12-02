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
                       iterationsPrior = 100000,
                       burnin = 1)

totalThetas <- addThetas(thetas = iTheta$thetas, iTheta = iTheta, keep =  2 : 3)

prior_pass_vec <- estimatePriorProbability(iTheta = iTheta,
                                           rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                                           iterationsPrior = 100000,
                                           cleanConstraints = iTheta$cleanConstraints,
                                           IDorg = "ID",
                                           effectNameOrg = "cond")

priorProbability <- mean(prior_pass_vec)

# microbenchmarks
install.packages("remotes")
library(remotes)
remotes::install_github("joshuaulrich/microbenchmark")

newWay <- function(x) {
  iTheta <- get_iTheta(formula = rtS ~ ID*cond,
                       data = stroop,
                       whichRandom = "ID",
                       ID = "ID",
                       whichConstraint = c("cond" = "2 > 1"),
                       rscaleEffects= c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                       iterationsPosterior = 3,
                       iterationsPrior = 100000,
                       burnin = 1)
  totalThetas <- addThetas(thetas = iTheta$thetas, iTheta = iTheta, keep =  2 : 3)
  prior_pass_vec <- estimatePriorProbability(iTheta = iTheta,
                                             rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10),
                                             iterationsPrior = 100000,
                                             cleanConstraints = iTheta$cleanConstraints,
                                             IDorg = "ID",
                                             effectNameOrg = "cond")
  priorProbability <- mean(prior_pass_vec)

}

constraintBF_function <- function(x){
  quid::constraintBF(formula = rtS ~ ID*cond,
                    data = stroop,
                    whichRandom = "ID",
                    ID = "ID",
                    whichConstraint = c("cond" = "2 > 1"),
                    rscaleEffects = c("ID" = 1, "cond" = 1/6, "ID:cond" = 1/10))
}

res <- microbenchmark::microbenchmark(constraintBF_function(x),
                                      newWay(x),
                                      times = 10,
                                      unit = "s")

print(res,
      unit = "s")

## Plot results:
boxplot(
  res,
  unit = "s",
  log = TRUE,
  xlab = "functions",
  horizontal = TRUE,
  main = "microbenchmark timings"
)
## Pretty plots:
if (requireNamespace("ggplot2")) {
  ggplot2::autoplot(res) +
  ggplot2::labs(
      title = "microbenchmark timings",
      subtitle = "constraintBF function vs. newWay"
    )
}

##boxplot with ggplot2
df <- as.data.frame(res)
df$time <- df$time / 1e9

### expr als Faktor setzen
df$expr <- factor(df$expr)

library(ggplot2)

ggplot(df, aes(x = expr, y = time, fill = expr, color = expr)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 2
  ) +
  scale_fill_manual(values = c("lightcoral", "lightgreen")) +
  scale_color_manual(values = c("firebrick3", "darkgreen")) +
  scale_y_log10() +
  coord_flip() +
  labs(
    x = "functions",
    y = "log(time) [s]",
    title = "microbenchmark timings"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

### Normales Boxplot
ggplot(df, aes(x = expr, y = time, fill = expr, color = expr)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 2,
    width = 0.6,
    fatten = 1.2
  ) +
  # Whisker mit Caps explizit zeichnen
  geom_errorbar(
    aes(ymin = ..ymin.., ymax = ..ymax..),
    stat = "boxplot",
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = c("lightcoral", "lightgreen")) +
  scale_color_manual(values = c("firebrick3", "darkgreen")) +
  scale_y_log10() +
  coord_flip() +
  labs(
    x = "functions",
    y = "log(time) [s]",
    title = "microbenchmark timings"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
