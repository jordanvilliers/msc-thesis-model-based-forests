library(ggplot2)
library(tidyr)
library(dplyr)

res_files <- list.files('../results', pattern = 'res_.*\\.rds', full.names = TRUE)
latest_file <- res_files[which.max(file.mtime(res_files))]
res <- readRDS(latest_file)

res <- res %>% 
  rename(
    n = prob.pars.n,
    rmvar = prob.pars.rmvar,
    models = prob.pars.models,
    designs = prob.pars.designs,
    xmodel = prob.pars.xmodel,
    mse = result.mse,
    bias_value = result.bias_value,
    variance = result.variance,
    est_ate = result.est_ate,
    coef_x2_m = prob.pars.coef_x2_m,
    coef_x3_m = prob.pars.coef_x3_m,
    coef_x1_t = prob.pars.coef_x1_t,
    coef_x2_t = prob.pars.coef_x2_t,
    progeff = prob.pars.progeff
  ) %>% 
  mutate(
    rmvar = ifelse(is.na(rmvar), 'None', as.character(rmvar)) %>% 
      factor(levels = c('None', 'X1', 'X2', 'X3'),
             labels = c('None', 'X1 (Predictive)', 'X2 (Overlays)', 'X3 (Prognostic)')),
    models = factor(models,
                    levels = c('normal', 'binomial', 'polr', 'weibull', 'coxph'),
                    labels = c('Normal', 'Binomial', 'Polr', 'Weibull', 'Coxph')),
    n = factor(n,
               levels = c(1500L),
               labels = c('1500')),
    designs = factor(designs,
                     levels = c('rct', 'obs'),
                     labels = c('RCT', 'Observational Data')),
    xmodel = factor(xmodel,
                    levels = c('unif', 'normal'),
                    labels = c('Uniform', 'Normal')),
    algorithm = ifelse(as.character(algorithm) == 'equalizedcox', 'equalized', as.character(algorithm)),
    algorithm = factor(algorithm,
                       levels = c('base', 'glm', 'mob', 'equalized'),
                       labels = c('Base', 'GLM', 'MOB', 'Equalized'))
  )

plot_function_compare <- function(data, metric, title, design, coef_filter = 2) {
  data_filtered <- data %>%
    filter(designs == design,
           coef_x3_m == coef_filter,
           !is.na(.data[[metric]]))
  
  p <- ggplot(data_filtered, aes(x = rmvar, y = .data[[metric]], fill = algorithm)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
    facet_grid(models ~ progeff, scales = 'free_y') +
    labs(title = paste(title, '-', design),
         x = 'Missing Variable', y = title, fill = 'Method') +
    scale_fill_manual(values = c('GLM' = 'darkgoldenrod2', 'MOB' = 'cyan4', 'Base' = '#CC79A7',
                                 'Equalized' = 'mediumpurple3')) +
    theme_bw() +
    theme(legend.position = 'top',
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  return(p)
}

plot_function_glm <- function(data, metric, title, design, coef_filter = 2) {
  data_filtered <- data %>%
    filter(designs == design,
           coef_x3_m == coef_filter)
  
  ggplot(data_filtered, aes(x = rmvar, y = .data[[metric]], fill = rmvar)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
    facet_grid(models ~ progeff, scales = 'free_y') +
    labs(title = paste('Method: GLM —', title),
         x = 'Missing Variable', y = title, fill = 'Omitted Variable') +
    scale_fill_manual(values = c(
      'None'            = 'darkgoldenrod2',
      'X1 (Predictive)' = 'cyan4',
      'X2 (Overlays)'   = 'green',
      'X3 (Prognostic)' = '#CC79A7'
    )) +
    theme_bw(base_size = 11) +
    theme(legend.position = 'top',
          plot.title = element_text(face = 'bold'),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_function_weibull_zoom <- function(data, metric, title, design, coef_filter = 2, ylim_range = c(0, 3)) {
  data_filtered <- data %>%
    filter(designs == design,
           coef_x3_m == coef_filter,
           models == 'Weibull',
           !is.na(.data[[metric]]))
  
  p <- ggplot(data_filtered, aes(x = rmvar, y = .data[[metric]], fill = algorithm)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
    facet_wrap(~ progeff) +
    coord_cartesian(ylim = ylim_range) +
    labs(title = paste(title, '- Weibull -', design),
         x = 'Missing Variable', y = title, fill = 'Method') +
    scale_fill_manual(values = c('GLM' = 'darkgoldenrod2', 'MOB' = 'cyan4', 'Base' = '#CC79A7',
                                 'Equalized' = 'mediumpurple3')) +
    theme_bw() +
    theme(legend.position = 'top',
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  return(p)
}

plot_specs <- expand.grid(
  metric = c('bias_value', 'variance', 'mse'),
  design = c('RCT', 'Observational Data'),
  coef_filter = c(2, 5),
  stringsAsFactors = FALSE
)
plot_specs$title <- c(bias_value = 'Bias', variance = 'Variance', mse = 'MSE')[plot_specs$metric]

for (i in seq_len(nrow(plot_specs))) {
  s <- plot_specs[i, ]
  print(plot_function_compare(res, s$metric, s$title, s$design, coef_filter = s$coef_filter))
}