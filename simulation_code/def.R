library(partykit)

testonly <- FALSE
seed <- 1003

repl <- 100
cores <- 30
probnam <- 'hte_simulation'
n_var <- 5

n <- c(1500L)
test_patient <- 1000
rmvar <- c(NA, 'X1', 'X2', 'X3')
models <- c('normal', 'binomial', 'polr', 'weibull', 'coxph')
designs <- c('rct', 'obs')
#xmodel <- c('unif', 'normal')
xmodel <- c('normal')
progeff <- c("linear" , "nonlinear")

current <- as.integer(Sys.time())
registry <- paste('registry', current, sep = "_")

# coef_x2_p <- c(2, 5)
coef_x2_p <- 2

make_pE <- function(coef_x2_p) {
  function (x) {
    return(1/(1 + exp(coef_x2_p * x[, "X2"])))
  }
}

# coef_x2_m <- c(2, 5)
coef_x2_m <- 2
coef_x3_m <- c(2, 5)
# coef_x3_m <- 2

make_mE <- function(coef_x2_m, coef_x3_m, progeff) {
  function(x) {
    if (progeff == "linear") {
      return(2 * coef_x2_m * x[,2] + coef_x3_m * x[,3])
    } else if (progeff == "nonlinear") {
      return(2 * log(1 + exp(coef_x2_m * x[, 2] + coef_x3_m * x[, 3])))
    } 
  }
}

# coef_x1_t <- c(2, 5)
coef_x1_t <- 2
# coef_x2_t <- c(2, 5)
coef_x2_t <- 2

make_tE <- function(coef_x1_t, coef_x2_t) {
  function(x) {
    return((coef_x1_t * x[, 'X1'] + coef_x2_t * x[, 'X2']) / 2)
  }
}

evaluate_performance <- function(estimated_effect, true_effect) {
  bias_vector <- estimated_effect - true_effect
  bias_value <- mean(bias_vector) 
  variance_value <- var(estimated_effect)
  mse_value <- mean(bias_vector^2)
  est_ate <- mean(estimated_effect)
  return(list(bias_value = bias_value, variance = variance_value, mse = mse_value, est_ate = est_ate))
}

num_trees <- 500
# num_trees <- 50 use it if needs fast results for last minute
min_size_group <- 7L
min_node_size <- min_size_group * 2L

ctrl <- ctree_control(testtype = 'Univ', minsplit = 2,
                      minbucket = min_node_size,
                      mincriterion = 0, saveinfo = FALSE)

ctrl$converged <- function(mod, data, subset) {
  if (!is.factor(data$trt)) {
    trtL <- all(table(data$trt[subset] > 0) > min_size_group)
    } else {
    trtL <- all(table(data$trt[subset]) > min_size_group)
    }

  if (is.factor(data$y)) {
    yL <- all(table(data$y[subset]) > min_size_group)
  } else {
    yL <- TRUE
  }

  return(all(trtL, yL))
}

min_update <- 20L
prt <- list(replace = FALSE, fraction = 0.5)
mtry <- 20L
