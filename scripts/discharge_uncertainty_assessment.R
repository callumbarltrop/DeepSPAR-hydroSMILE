# Getting SLURM batch ID  -----------------------------------

rm(list = ls())
setwd(this.path::here())

# It is highly recommended to run this on a SLURM submission system. Set the array ids on the system to run from 1 - 400

array_id = commandArgs(trailingOnly = T) #Assuming job array submission with SLURM
array_id = as.numeric(array_id)

# If SLURM + HPC is not available, the code can be run in a for loop via
# for(array_id in 1:400){
#
#   #The rest of the script goes here
#
# }
# However be warned this will take a while

# NN and GPD tuning parameters --------------------------------------------

# Flag for fitting models. Can be set to F once models weights are saved
fit_models = T

nn.nunits = c(32, 32, 32)

quant.level = 0.85

init_shape = 0.1

l.r.vec = c(0.001, 5e-04, 0.00025, 1e-04, 5e-05)

gpd.shape.custom_activation <- function(x) {
  tf$atan(x) / pi
}

init_bias = tan(pi * (init_shape))

activation.func = "relu"

# Number of bootstraps
n_boot = 100

n_y = 52

# Packages and VE ---------------------------------------------------------

source("preamble.R")

reticulate::use_virtualenv("DeepSPAR_env", required = T)
packages = c("keras", "tensorflow", "reticulate")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

tensorflow::set_random_seed(array_id + 1)

# Data --------------------------------------------------------------------

river_names = c("ammer", "iller", "lech", "loisach")

time_periods = list(1980:2009, 2010:2039, 2040:2069, 2070:2099)

std_sigmas = readRDS(file = "../data/preprocessed_data/std_sigmas.rds")

std_mus = readRDS(file = "../data/preprocessed_data/std_mus.rds")

std_tw_dfs = readRDS(file = "../data/preprocessed_data/std_tw_dfs.rds")

os_tw_dfs = readRDS(file = "../data/preprocessed_data/os_tw_dfs.rds")

# Estimating threshold function ------------------------------------------------

indices_grid = expand.grid(1:n_boot, c(1:length(time_periods)))

print(paste("Processing task", array_id))

# Index of bootstrap iteration
boot_ind = indices_grid$Var1[array_id]

# Index of time window
tw_ind = indices_grid$Var2[array_id]

X.norm = std_tw_dfs[[tw_ind]]

# Resample rows in pre-processed data randomly
set.seed(boot_ind + 1)
X.norm = X.norm[sample(1:nrow(X.norm), nrow(X.norm), replace = T), ]

n = dim(X.norm)[1]
d = dim(X.norm)[2]

R <- apply(X.norm, 1, l2_norm)
W = X.norm / R
W = as.matrix(W)

valid.inds = sample(1:n, round(n / 10))

R.train <- R[-valid.inds]
W.train <- W[-valid.inds, ]
R.valid <- R[valid.inds]
W.valid <- W[valid.inds, ]

input_pseudo_angles <- layer_input(shape = d, name = "input_pseudo_angles")

qBranch <- input_pseudo_angles %>%
  layer_dense(
    units = nn.nunits[1],
    activation = activation.func,
    name = "q_dense1"
  )

for (i in 2:length(nn.nunits)) {
  qBranch <- qBranch %>%
    layer_dense(
      units = nn.nunits[i],
      activation = activation.func,
      name = paste0("q_dense", i)
    )
  
}

qBranch <- qBranch %>%
  layer_dense(
    units = 1,
    activation = "exponential",
    name = "q_final"
  )

quant.model <- keras_model(inputs = c(input_pseudo_angles),
                           outputs = c(qBranch))


quant.model %>%
  compile(optimizer = "adam",
          loss = tilted_loss,
          run_eagerly = T)

if (fit_models) {
  checkpoint_dir <- tempdir()
  checkpoint_path <- file.path(checkpoint_dir,
                               paste0("QR_checkpoint_bootind_", boot_ind, "_tw_", tw_ind))
  
  checkpoint <- callback_model_checkpoint(
    filepath = checkpoint_path,
    monitor = "val_loss",
    verbose = 0,
    save_best_only = TRUE,
    save_weights_only = TRUE,
    mode = "min",
    save_freq = "epoch"
  )
  n.epochs <- 500
  batch.size <- 1024
  
  dim(R) = c(length(R), 1)
  dim(R.train) = c(length(R.train), 1)
  dim(R.valid) = c(length(R.valid), 1)
  
  history <- quant.model %>%
    fit(
      list(W.train),
      R.train,
      epochs = n.epochs,
      batch_size = batch.size,
      callback = list(
        checkpoint,
        callback_early_stopping(
          monitor = "val_loss",
          min_delta = 0,
          patience = 5
        )
      ),
      validation_data = list(list(input_pseudo_angles = W.valid), R.valid)
    )
  
  quant.model <- load_model_weights_tf(quant.model, filepath = checkpoint_path)
  
  save_model_tf(
    quant.model,
    paste0(
      "../nn_parameters/QR_est/qr_fit_bootind_",
      boot_ind,
      "_tw_",
      tw_ind
    )
  )
  
  unlink(checkpoint_dir, recursive = TRUE)
  
}

quant.model <- tf$saved_model$load(paste0(
  "../nn_parameters/QR_est/qr_fit_bootind_",
  boot_ind,
  "_tw_",
  tw_ind
))

pred.quant <- k_get_value(quant.model(k_constant(W)))

# Estimating GPD parameter functions --------------------------------------

u <- pred.quant

R.train <- (R - u)[-valid.inds]
W.train <- W[-valid.inds, ]
u.train <- u[-valid.inds]
R.valid <- (R - u)[valid.inds]
W.valid <- W[valid.inds, ]
u.valid <- u[valid.inds]

W.train <- W.train[R.train > 0, ]
u.train <- u.train[R.train > 0]
R.train <- R.train[R.train > 0]
W.valid <- W.valid[R.valid > 0, ]
u.valid <- u.valid[R.valid > 0]
R.valid <- R.valid[R.valid > 0]

dim(R.train) = c(length(R.train), 1)
dim(R.valid) = c(length(R.valid), 1)

dim(u.train) = c(length(u.train), 1)
dim(u.valid) = c(length(u.valid), 1)

input_pseudo_angles <- layer_input(shape = d, name = "input_pseudo_angles")

input_u <- layer_input(shape = d, name = "input_u")

gpd.Branch <- input_pseudo_angles %>%
  layer_dense(
    units = nn.nunits[1],
    activation = activation.func,
    name = "gpd_dense1"
  )

if (length(nn.nunits) > 2) {
  for (i in 2:(length(nn.nunits) - 1)) {
    gpd.Branch <- gpd.Branch %>%
      layer_dense(
        units = nn.nunits[i],
        activation = activation.func,
        name = paste0("gpd_dense", i)
      )
    
  }
  
}

scale.Branch <- gpd.Branch %>%
  layer_dense(
    units = nn.nunits[length(nn.nunits)],
    activation = activation.func,
    name = paste0("scale_pre")
  ) %>%
  layer_dense(
    units = 1,
    activation = "exponential",
    name = "scale_final",
    weights = list(matrix(0, nrow = nn.nunits[length(nn.nunits)], ncol = 1), array(0, dim = c(1)))
  )


shape.Branch <- gpd.Branch %>%
  layer_dense(
    units = nn.nunits[length(nn.nunits)],
    activation = activation.func,
    name = paste0("shape_pre")
  ) %>%
  layer_dense(
    units = 1,
    activation = gpd.shape.custom_activation,
    name = "shape_final",
    weights = list(matrix(0, nrow = nn.nunits[length(nn.nunits)], ncol = 1), array(init_bias, dim = c(1)))
  )

output <- layer_concatenate(c(scale.Branch, shape.Branch, input_u, input_pseudo_angles))

GPD.model <- keras_model(inputs = c(input_pseudo_angles, input_u),
                         outputs = output)
summary(GPD.model)

GPD.model %>%
  compile(
    optimizer = optimizer_adam(learning_rate = l.r.vec[1]),
    loss = GPD_loss,
    run_eagerly = T
  )

if (fit_models) {
  n.epochs <- 750
  batch.size = length(R.train)
  
  checkpoint <- callback_model_checkpoint(
    filepath = paste0(
      "../nn_parameters/GPD_est/tmp_bootind_",
      boot_ind,
      "_tw_",
      tw_ind,
      "/gpd_fit_{epoch:02d}"
    ),
    monitor = "val_loss",
    verbose = 0,
    save_best_only = FALSE,
    save_weights_only = TRUE,
    mode = "min",
    save_freq = "epoch"
  )
  
  history <- GPD.model %>%
    fit(
      list(W.train, u.train),
      R.train,
      epochs = n.epochs,
      batch_size = batch.size,
      callbacks = list(
        checkpoint,
        callback_early_stopping(
          monitor = "val_loss",
          min_delta = 0,
          patience = 5
        )
      ),
      validation_data = list(list(W.valid, u.valid), R.valid)
    )
  
  optimal_epoch = which.min(history$metrics$val_loss)
  
  GPD.model <- load_model_weights_tf(
    GPD.model,
    filepath = paste0(
      "../nn_parameters/GPD_est/tmp_bootind_",
      boot_ind,
      "_tw_",
      tw_ind,
      "/gpd_fit_",
      sprintf("%02d", optimal_epoch - 1)
    )
  )
  
  save_model_weights_tf(
    GPD.model,
    paste0(
      "../nn_parameters/GPD_est/tmp_bootind_",
      boot_ind,
      "_tw_",
      tw_ind,
      "/gpd_fit_best"
    )
  )
  
  save_model_tf(
    GPD.model,
    paste0(
      "../nn_parameters/GPD_est/gpd_fit_bootind_",
      boot_ind,
      "_tw_",
      tw_ind
    )
  )
  
  best.valid.loss <- history$metrics$val_loss[which.min(history$metrics$val_loss) - 1]
  
  n.epochs <- 200
  
  for (k in 2:length(l.r.vec)) {
    print(paste0("Decreasing learning rate size to: ", l.r.vec[k]))
    
    GPD.model %>%
      compile(
        optimizer = optimizer_adam(learning_rate = l.r.vec[k]),
        loss = GPD_loss,
        run_eagerly = T
      )
    
    history <- GPD.model %>%
      fit(
        list(W.train, u.train),
        R.train,
        epochs = n.epochs,
        batch_size = batch.size,
        callback = list(
          checkpoint,
          callback_early_stopping(
            monitor = "val_loss",
            min_delta = 0,
            patience = 10
          )
        ),
        validation_data = list(
          list(input_pseudo_angles = W.valid, input_u = u.valid),
          R.valid
        )
      )
    
    if (length(which.min(history$metrics$val_loss)) > 0) {
      if (min(history$metrics$val_loss, na.rm = T) > best.valid.loss) {
        optimal_epoch = NA
        
      } else {
        optimal_epoch = which.min(history$metrics$val_loss)
        
        best.valid.loss <- history$metrics$val_loss[optimal_epoch]
        
      }
      
    } else {
      optimal_epoch = NA
      
    }
    
    if (!is.na(optimal_epoch)) {
      GPD.model <- load_model_weights_tf(
        GPD.model,
        filepath = paste0(
          "../nn_parameters/GPD_est/tmp_bootind_",
          boot_ind,
          "_tw_",
          tw_ind,
          "/gpd_fit_",
          sprintf("%02d", optimal_epoch)
        )
      )
      
      save_model_weights_tf(
        GPD.model,
        paste0(
          "../nn_parameters/GPD_est/tmp_bootind_",
          boot_ind,
          "_tw_",
          tw_ind,
          "/gpd_fit_best"
        )
      )
      
      save_model_tf(
        GPD.model,
        paste0(
          "../nn_parameters/GPD_est/gpd_fit_bootind_",
          boot_ind,
          "_tw_",
          tw_ind
        )
      )
      
    } else{
      load_model_weights_tf(
        GPD.model,
        paste0(
          "../nn_parameters/GPD_est/tmp_bootind_",
          boot_ind,
          "_tw_",
          tw_ind,
          "/gpd_fit_best"
        )
      )
      
    }
    
  }
  
  
  unlink(
    paste0(
      "../nn_parameters/GPD_est/tmp_bootind_",
      boot_ind,
      "_tw_",
      tw_ind,
      ""
    ),
    recursive = T,
    force = T
  )
  
}

GPD.model <- load_model_tf(
  filepath = paste0(
    "../nn_parameters/GPD_est/gpd_fit_bootind_",
    boot_ind,
    "_tw_",
    tw_ind
  ),
  custom_objects = list(GPD_loss = GPD_loss)
)

# Computing diagnostic quantities  ----------------------------------------

pred.GPD.exc <- k_get_value(GPD.model(list(k_constant(W[R - u > 0, ]), 
                                           k_constant(matrix(u[R - u > 0], ncol = 1)
))))

obs_quants = qexp(apply(cbind(u[R - u > 0], pred.GPD.exc[, 1:2], R[R - u >  0]),
                        1, 
                        function(x) pgpd(x[4], loc = x[1], scale = x[2] / (1 + x[3]), x[3])))

ps = readRDS(file = paste0("../computed_model_vals/gpd_qq_prob_tw_", tw_ind, ".rds"))

boot_quants = quantile(obs_quants, ps)

saveRDS(
  boot_quants,
  file = paste0(
    "../computed_model_vals/gpd_qq_quants_",
    "bootind_",
    boot_ind,
    "_tw_",
    tw_ind,
    ".rds"
  )
)

n.sim.tail = 2e6

n.sim = round((n.sim.tail / (1 - quant.level)))

X.norm.body = X.norm[which(pred.quant >= R), ]

X.sim.body = X.norm.body[sample(1:nrow(X.norm.body), n.sim - n.sim.tail, replace =  T), ]

X.norm.excess = X.norm[which(pred.quant < R), ]

W.excess = W[which(pred.quant < R), ]

W.sim = W.excess[sample(1:nrow(W.excess), n.sim.tail, replace = T), ]
u.sim <- k_get_value(quant.model(k_constant(W.sim)))

dim(u.sim) = c(length(u.sim), 1)

pred.GPD.sim <- k_get_value(GPD.model(list(k_constant(W.sim), k_constant(u.sim))))
nu.sim = pred.GPD.sim[, 1]
xi.sim = pred.GPD.sim[, 2]

unif_sample = runif(n.sim.tail, 0, 1)
R.sim = apply(cbind(unif_sample, u.sim, nu.sim, xi.sim), 1, function(x)
  qgpd(
    x[1],
    loc = x[2],
    scale = x[3] / (1 + x[4]),
    shape = x[4]
  ))

X.sim.excess = R.sim * W.sim

X.sim = rbind(X.sim.body, X.sim.excess)

rm(u.sim, pred.GPD.sim, unif_sample, nu.sim, xi.sim, R.sim, W.sim)

mus = std_mus[[tw_ind]]

sigmas = std_sigmas[[tw_ind]]

X.excess.obs.os = apply(rbind(mus, X.norm.excess), 2, function(x) {
  return(x[-1] + x[1])
})

X.excess.obs.os = apply(X.excess.obs.os, 2, function(x) {
  return(log(1 + exp(x)))
})

X.excess.obs.os = apply(rbind(sigmas, X.excess.obs.os), 2, function(x) {
  return(x[-1] * x[1])
})

X.excess.sim.os = apply(rbind(mus, X.sim.excess), 2, function(x) {
  return(x[-1] + x[1])
})

X.excess.sim.os = apply(X.excess.sim.os, 2, function(x) {
  return(log(1 + exp(x)))
})

X.excess.sim.os = apply(rbind(sigmas, X.excess.sim.os), 2, function(x) {
  return(x[-1] * x[1])
})

labels <- river_names

u_vec = seq(0.8, 0.99, by = 0.01)

p <- ncol(X.excess.sim.os)

for (j in 1:p) {
  for (i in 1:p) {
    if (i == j) {
      ps = readRDS(
        file = paste0(
          "../computed_model_vals/marg_qq_prob_",
          labels[j],
          "_tw_",
          tw_ind,
          ".rds"
        )
      )
      
      sim_quants = quantile(X.excess.sim.os[, j], ps)
      
      saveRDS(
        sim_quants,
        file = paste0(
          "../computed_model_vals/marg_qq_quants_",
          labels[j],
          "_bootind_",
          boot_ind,
          "_tw_",
          tw_ind,
          ".rds"
        )
      )
      
    } else if (j > i) {
      NULL
      
    } else {
      xs <- X.excess.sim.os[, i]
      ys <- X.excess.sim.os[, j]
      
      lam_sim <- sapply(u_vec, tail_dep, x = xs, y = ys)
      
      saveRDS(
        lam_sim,
        file = paste0(
          "../computed_model_vals/chi_plot_",
          labels[i],
          "_",
          labels[j],
          "_bootind_",
          boot_ind,
          "_tw_",
          tw_ind,
          ".rds"
        )
      )
      
    }
  }
}

for (i in 1:d) {
  title = river_names[i]
  
  fix_RPs = exp(seq(log(10), log(10000), length.out = 3000))
  
  fix_probs = 1 / (fix_RPs * n_y)
  
  model_RLs = quantile(X.sim[, i], 1 - fix_probs)
  
  model_RLs = log(1 + exp(model_RLs + mus[i])) * sigmas[i]
  
  saveRDS(
    model_RLs,
    file = paste0(
      "../computed_model_vals/rl_ests_upper_",
      title,
      "_bootind_",
      boot_ind,
      "_tw_",
      tw_ind,
      ".rds"
    )
  )
  
}

for (i in 1:d) {
  title = river_names[i]
  
  fix_RPs = exp(seq(log(10), log(10000), length.out = 3000))
  
  fix_probs = 1 / (fix_RPs * n_y)
  
  model_RLs = quantile(X.sim[, i], fix_probs)
  
  model_RLs = log(1 + exp(model_RLs + mus[i])) * sigmas[i]
  
  saveRDS(
    model_RLs,
    file = paste0(
      "../computed_model_vals/rl_ests_lower_",
      title,
      "_bootind_",
      boot_ind,
      "_tw_",
      tw_ind,
      ".rds"
    )
  )
  
}

# Estimating tail probabilities from SPAR model ----------------------------------------------------------------------

l1_seq_upper = 1000

l1_seq_lower = 30

prob_RPs = c(10)

est_RLs_upper = matrix(NA, ncol = d, nrow = length(prob_RPs))

est_RLs_lower = matrix(NA, ncol = d, nrow = length(prob_RPs))

for (i in 1:d) {
  est_RLs_upper[, i] = quantile(X.norm[, i], 1 - (1 / (n_y * prob_RPs)))
  
  est_RLs_upper[, i] = log(1 + exp(est_RLs_upper[, i] + mus[i])) * sigmas[i]
  
  est_RLs_lower[, i] = quantile(X.norm[, i], (1 / (n_y * prob_RPs)))
  
  est_RLs_lower[, i] = log(1 + exp(est_RLs_lower[, i] + mus[i])) * sigmas[i]
  
}

X.sim = apply(rbind(mus, X.sim), 2, function(x) {
  return(x[-1] + x[1])
})

X.sim = apply(X.sim, 2, function(x) {
  return(log(1 + exp(x)))
})

X.sim = apply(rbind(sigmas, X.sim), 2, function(x) {
  return(x[-1] * x[1])
})

l1_upper_ests = sapply(l1_seq_upper, function(l1_q) {
  return(l1_prob_upper(u = l1_q, sim_data = X.sim))
  
})

l1_lower_ests = sapply(l1_seq_lower, function(l1_q) {
  return(l1_prob_lower(u = l1_q, sim_data = X.sim))
  
})

upper_surv_ests = apply(est_RLs_upper, 1, function(ret_levels) {
  return(survivor_prob(u_vec = ret_levels, sim_data = X.sim))
  
})

lower_cdf_ests = apply(est_RLs_lower, 1, function(ret_levels) {
  return(cdf_prob(u_vec = ret_levels, sim_data = X.sim))
  
})

probs_list = list(
  l1_seq_upper = l1_seq_upper,
  l1_seq_lower = l1_seq_lower,
  est_RLs_upper = est_RLs_upper,
  est_RLs_lower = est_RLs_lower,
  l1_upper_ests = l1_upper_ests,
  l1_lower_ests = l1_lower_ests,
  upper_surv_ests = upper_surv_ests,
  lower_cdf_ests = lower_cdf_ests
)

saveRDS(
  probs_list,
  file = paste0(
    "../computed_model_vals/model_prob_ests_bootind_",
    boot_ind,
    "_tw_",
    tw_ind,
    ".rds"
  )
)
