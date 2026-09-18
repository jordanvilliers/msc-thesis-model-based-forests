
############# Experiments 
#----
# 0) Setup
#----

source('def.R')
set.seed(seed)
setups <- expand.grid(repl = seq(1, repl), n = n, rmvar = rmvar, models = models, designs = designs,
                      xmodel = xmodel, coef_x2_p = coef_x2_p,
                      coef_x2_m = coef_x2_m, coef_x3_m = coef_x3_m, 
                      coef_x1_t = coef_x1_t, coef_x2_t = coef_x2_t, progeff = progeff, stringsAsFactors = FALSE)
setups$seed <- round(runif(nrow(setups)) * 1e8)
pkgs <- c('batchtools', 'model4you', 'htesim', 'survival', 'MASS',
          'partykit', 'trtf', 'marginaleffects', 'gbm', 'mboost', 'grf')

#----
# 1) Load helper functions & libraries
#----

rpkgs <- sapply(pkgs, require, character.only = TRUE)
if (!all(rpkgs))
  sapply(pkgs[!rpkgs], install.packages)
rpkgs <- sapply(pkgs, require, character.only = TRUE)
if (!all(rpkgs))
  stop('could not attach all required packages')

generate_data <- function(n = 100L, rmvar = NA, models, designs, seed = NULL, xmodel = 'unif', coef_x2_p = 1,
                          coef_x2_m = 1, coef_x3_m = 1, coef_x1_t = 1, coef_x2_t = 1, progeff = "linear") {
  pE_only_x2 <- make_pE(coef_x2_p)
  mE_custom_muC <- make_mE(coef_x2_m, coef_x3_m, progeff)
  tE_strong_effect <- make_tE(coef_x1_t, coef_x2_t)
  
  if (!exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  } 
  if (is.null(seed)) {
    RNGstate <- get('.Random.seed', envir = .GlobalEnv)
  } else {
    R.seed <- get('.Random.seed', envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign('.Random.seed', R.seed, envir = .GlobalEnv))
  }
  
  if (is.na(rmvar)) {
    rmvar <- NULL
  }
  
  dgp_model_type <- ifelse(models == 'coxph', 'weibull', as.character(models))
  
  if (designs == 'rct') {
    dgp_obj <- htesim::dgp(p = 0.5, m = mE_custom_muC, t = tE_strong_effect, model = dgp_model_type, xmodel = xmodel, rmvar = rmvar)
  } else if (designs == 'obs') {
    dgp_obj <- htesim::dgp(p = pE_only_x2, m = mE_custom_muC, t = tE_strong_effect, model = dgp_model_type, xmodel = xmodel, rmvar = rmvar)
  }
  
  sim_data <- simulate(dgp_obj, nsim = n, nsimtest = test_patient, dim = n_var, seed = seed)
  
  test_data <- attr(sim_data, 'testxdf')
  true_indiv_effect <- predict(dgp_obj, newdata = test_data)[, 'tfct']
  
  attr(sim_data, 'test_data_split') <- test_data
  attr(sim_data, 'ground_truth') <- true_indiv_effect
  attr(sim_data, 'model_type') <- models
  attr(sim_data, 'design_type') <- designs
  
  return(sim_data)
}

fun_mob <- function(instance, ...) {
  test_data <- attr(instance, 'test_data_split')
  true_indiv_effect <- attr(instance, 'ground_truth')
  model_type <- attr(instance, 'model_type')
  
  if (model_type == 'normal') {
    base_model <- lm(y ~ trt, data = instance)
  } else if (model_type == 'binomial') {
    base_model <- glm(y ~ trt, data = instance, family = 'binomial')
  } else if (model_type == 'polr') {
    base_model <- polr(y ~ trt, data = instance)
  } else if (model_type == 'weibull') {
    offset <- NULL
    base_model <- as.mlt(Survreg(y ~ trt, data = instance, offset = offset))
  } else if (model_type == 'coxph') {
    base_model = survival::coxph(y ~ trt, data = instance)
  }
  
  if (model_type == 'weibull') {
    nm <- c('(Intercept)', 'trt1')
    instance$offset <- NULL
    rf_model <- trtf::traforest(base_model, formula = y | trt ~ .,
                                data = instance, ntree = num_trees,
                                perturb = prt, mtry = mtry, control = ctrl,
                                parm = nm, min_update = min_update,
                                mltargs = list(offset = offset), trace = FALSE)
    cf <- predict(rf_model, newdata = test_data, type = 'coef')
    cf <- do.call('rbind', cf)
    est_cate <- cf[, 'trt1']
  } else {
    rf_model <- pmforest(base_model, data = instance, ntree = num_trees, perturb = prt, 
                         mtry = mtry, control = ctrl, trace = FALSE)
    rf_pred <- suppressWarnings(
      pmodel(rf_model, newdata = test_data)
    )
    if (model_type == 'polr') {
      est_cate <- c(rf_pred)
    } else if (model_type == 'coxph') {
      est_cate <- -c(rf_pred)
    } else if ('trt1' %in% colnames(rf_pred)) {
      est_cate <- rf_pred[, 'trt1']
    } else {
      est_cate <- rf_pred[, 'trt']
    }
  }
  
  results <- evaluate_performance(estimated_effect = est_cate, true_effect = true_indiv_effect)
  
  return(results)
}

fun_base <- function(instance, ...) {
  test_data <- attr(instance, 'test_data_split')
  true_indiv_effect <- attr(instance, 'ground_truth')
  model_type <- attr(instance, 'model_type')
  
  if (model_type == 'normal') {
    base_model <- lm(y ~ trt, data = instance)
    est_ate <- rep(coef(base_model)['trt1'], nrow(test_data))
  } else if (model_type == 'binomial') {
    base_model <- glm(y ~ trt, data = instance, family = 'binomial')
    est_ate <- rep(coef(base_model)['trt1'], nrow(test_data))
  } else if (model_type == 'polr') {
    base_model <- polr(y ~ trt, data = instance)
    est_ate <- rep(-coef(base_model)['trt1'], nrow(test_data))
  } else if (model_type == 'weibull') {
    offset <- NULL
    base_model <- as.mlt(Survreg(y ~ trt, data = instance, offset = offset))
    est_ate <- rep(coef(base_model)['trt1'], nrow(test_data))
  } else if (model_type == 'coxph') {
    base_model = survival::coxph(y ~ trt, data = instance)
    est_ate <- rep(-coef(base_model)['trt1'], nrow(test_data))
  }
  results <- evaluate_performance(estimated_effect = est_ate, true_effect = true_indiv_effect)
  
  return(results)
}

fun_glm <- function(instance, ...) {
  test_data <- attr(instance, 'test_data_split')
  
  true_indiv_effect <- attr(instance, 'ground_truth')
  model_type <- attr(instance, 'model_type')
  
  xnam <- names(instance)[grep("X", names(instance))]
  intterms <- paste(paste0(xnam, ":trt"), collapse = " + ")
  form <- as.formula(paste("y ~ . +", intterms))
  
  if (model_type != "weibull") {
    if (model_type == 'normal') {
      base_model <- lm(form, data = instance)
      compmeas <- "difference"
      type <- NULL
    } else if (model_type == 'binomial') {
      base_model <- glm(form, data = instance, family = 'binomial')
      compmeas <- "lnor"
      type <- NULL
    } else if (model_type == 'polr') {
      base_model <- MASS::polr(form, data = instance)
      compmeas <- "lnor"
      type <- NULL
    } else if (model_type == 'coxph') {
      base_model <- survival::coxph(form, data = instance)
      compmeas <- 'difference'
      type <- 'lp'
    } 
    
    est_cate <- marginaleffects::comparisons(base_model,
                                             variables = "trt",
                                             type = type,
                                             comparison = compmeas,
                                             newdata = test_data
    )$estimate
    
    if (model_type %in% c("polr", "coxph")) {
      est_cate <- -est_cate[1:nrow(test_data)]
    }
  } else if (model_type == 'weibull') {
    base_model <- as.mlt(tram::Survreg(form, data = instance))
    rhs_form <- as.formula(paste("~", deparse(form[[3]])))
    cols_to_keep <- names(test_data)[names(test_data) %in% names(instance)]
    test_data <- test_data[, cols_to_keep, drop = FALSE]
    test_data1 <- test_data
    test_data1$trt <- factor(1, levels = levels(test_data$trt))
    test_data0 <- test_data
    test_data0$trt <- factor(0, levels = levels(test_data$trt))
    X1 <- model.matrix(rhs_form, data = test_data1)[, -1, drop = FALSE]
    X0 <- model.matrix(rhs_form, data = test_data0)[, -1, drop = FALSE]
    shift_coef <- base_model$coef[base_model$shiftcoef]
    est_cate <- as.vector((X1 - X0) %*% shift_coef)
  } 
  
  results <- evaluate_performance(estimated_effect = est_cate, true_effect = true_indiv_effect)
  return(results)
}

fun_equalized <- function(instance, ...) {
  test_data <- attr(instance, 'test_data_split')
  true_indiv_effect <- attr(instance, 'ground_truth')
  model_type <- attr(instance, 'model_type')
  design_type <- attr(instance, 'design_type')
  
  offset <- NULL
  nm <- nmt <- 'trt1'
  
  W.hat <- if (design_type == 'rct') 0.5 else NULL
  
  d0 <- instance[instance$trt == 0, ]
  d0$trt <- NULL
  d1 <- instance[instance$trt == 1, ]
  d1$trt <- NULL
  
  if (model_type == 'binomial') {
    d0$y <- c(0, 1)[d0$y]
    d1$y <- c(0, 1)[d1$y]
    boost.Y0 <- gbm::gbm(y ~ ., data = d0, distribution = 'bernoulli', interaction.depth = 2L)
    boost.Y1 <- gbm::gbm(y ~ ., data = d1, distribution = 'bernoulli', interaction.depth = 2L)
  } else if (model_type %in% c('weibull', 'coxph')) {
    boost.Y0 <- gbm::gbm(y ~ ., data = d0, distribution = 'coxph', interaction.depth = 2L)
    boost.Y1 <- gbm::gbm(y ~ ., data = d1, distribution = 'coxph', interaction.depth = 2L)
  } else if (model_type == 'polr') {
    boost.Y0 <- mboost::blackboost(y ~ ., data = d0, family = mboost::PropOdds())
    boost.Y1 <- mboost::blackboost(y ~ ., data = d1, family = mboost::PropOdds())
  }
  
  X <- as.matrix(instance[, grep('^X', colnames(instance))])
  W <- (0:1)[instance$trt]
  if (model_type != 'normal') {
    Y <- rep(0, nrow(X))
  } else {
    Y <- as.numeric(instance$y)
  }
  cf <- grf::causal_forest(X = X, Y = Y, W = W, W.hat = W.hat,
                           stabilize.splits = TRUE, min.node.size = min_size_group,
                           sample.fraction = prt$fraction, mtry = mtry, ci.group.size = 1,
                           num.trees = num_trees, honesty = FALSE)
  myW.hat <- cf$W.hat
  a <- myW.hat
  
  instance$trt <- (0:1)[instance$trt] - a
  nm <- nmt <- 'trt'
  
  if (model_type == 'normal') {
    Y.hat <- cf$Y.hat
    instance$y <- instance$y - Y.hat
  } else {
    if (model_type == 'binomial') {
      lp0 <- as.numeric(predict(boost.Y0, newdata = instance, type = 'link', n.trees = boost.Y0$n.trees))
      lp1 <- as.numeric(predict(boost.Y1, newdata = instance, type = 'link', n.trees = boost.Y1$n.trees))
    } else if (model_type %in% c('weibull', 'coxph')) {
      lp0 <- predict(boost.Y0, newdata = instance, type = 'link', n.trees = boost.Y0$n.trees)
      lp1 <- predict(boost.Y1, newdata = instance, type = 'link', n.trees = boost.Y1$n.trees)
    } else if (model_type == 'polr') {
      lp0 <- as.numeric(predict(boost.Y0, newdata = instance, type = 'link'))
      lp1 <- as.numeric(predict(boost.Y1, newdata = instance, type = 'link'))
    }
    nu <- myW.hat * lp1 + (1 - myW.hat) * lp0
    instance$offset <- offset <- nu
  }
  
  if (model_type == 'normal') {
    instance$offset <- NULL
    base_model <- lm(y ~ trt, data = instance)
  } else if (model_type == 'weibull') {
    base_model <- as.mlt(Survreg(y ~ trt, data = instance, offset = offset))
  } else if (model_type == 'polr') {
    base_model <- polr(y ~ trt + offset(offset), data = instance)
  } else if (model_type == 'binomial') {
    base_model <- glm(y ~ trt + offset(offset), data = instance, family = binomial)
  } else if (model_type == 'coxph') {
    base_model <- survival::coxph(y ~ trt + offset(offset), data = instance)
  }
  
  if (model_type == 'weibull') {
    nm <- c('(Intercept)', nmt)
    instance$offset <- NULL
    rf_model <- trtf::traforest(base_model, formula = y | trt ~ ., data = instance, ntree = num_trees,
                                perturb = prt, mtry = mtry, control = ctrl, parm = nm, min_update = min_update,
                                mltargs = list(offset = offset), trace = FALSE)
    cf <- predict(rf_model, newdata = test_data, type = 'coef')
    cf <- do.call('rbind', cf)
    est_cate <- cf[, nmt]
  } else {
    rf_model <- pmforest(base_model, data = instance, ntree = num_trees, perturb = prt,
                         mtry = mtry, control = ctrl, trace = FALSE)
    rf_pred <- suppressWarnings(pmodel(rf_model, newdata = test_data))
    if (model_type == 'polr') {
      est_cate <- c(rf_pred)
    } else if (model_type == "coxph") {
      est_cate <- -c(rf_pred)
    } else {
      est_cate <- rf_pred[, nmt]
    }
  }
  
  results <- evaluate_performance(estimated_effect = est_cate, true_effect = true_indiv_effect)
  return(results)
}


methods <- list('mob' = fun_mob,
                'glm' = fun_glm,
                'base' = fun_base,
                'equalized' = fun_equalized)

#-----
# 2) Create study environment (TEST/NO TEST)
# Result: experimental registry
#-----

if (!dir.exists("registry"))
  dir.create("registry")
if (!dir.exists("results"))
  dir.create("results")
reg <- makeExperimentRegistry(file.dir = file.path("registry", registry), 
                              packages = pkgs, seed = seed)

reg$default.resources <- list(
  ntasks = 1L,
  ncpus = 1L,
  nodes = 1L,
  clusters = "serial")
reg$cluster.functions <- makeClusterFunctionsMulticore(cores)

#----
# 3) Add problem = dataset
#----

fun <- function(job, n = 100L, rmvar = NA, models, designs, seed, xmodel = 'unif', coef_x2_p = 1,
                coef_x2_m = 1, coef_x3_m = 1, coef_x1_t = 1, coef_x2_t = 1, progeff, ...) {
  d <- generate_data(n = n, rmvar = rmvar, models = models, designs = designs, seed = seed,
                     xmodel = xmodel, coef_x2_p = coef_x2_p, coef_x2_m = coef_x2_m,
                     coef_x3_m = coef_x3_m, coef_x1_t = coef_x1_t, coef_x2_t = coef_x2_t, progeff)
  return(d)
}

addProblem(probnam, fun = fun, reg = reg)
prob.designs <- list()
prob.designs[[probnam]] <- setups

#-----
# 4) Add algorithms
#-----

algo.designs <- list()
for (method in names(methods)) {
  addAlgorithm(method, fun = methods[[method]], reg = reg)
  algo.designs[[method]] <- data.frame()
}

#----
# 5) add experiments
#----

addExperiments(prob.designs, algo.designs, reg = reg)

#----
# 6) check setup
#----

summarizeExperiments(reg = reg)
jobpars <- unwrap(getJobPars(reg = reg))

if (testonly) {
  testJob(1L, reg = reg)
  testJob(nrow(jobpars), reg = reg)
} else {
  
  #-----
  # 7) submit jobs
  #----
  
  submitJobs(reg = reg)
  waitForJobs(reg = reg)
  
  #-----
  # 8) save results
  #----
  
  res <- ijoin(
    getJobPars(reg = reg),
    reduceResultsDataTable(reg = reg, fun = function(x) as.list(x))
  )
  res <- unwrap(res, sep = ".")
  
  resname <- paste0('../results/res_', current, '.rds')
  saveRDS(res, file = resname)
}
