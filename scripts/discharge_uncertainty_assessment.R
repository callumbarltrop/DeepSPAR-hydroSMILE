# Getting SLURM batch ID  -----------------------------------

rm(list=ls())

# It is highly recommended to run this on a SLURM submission system. Set the array_ids on the system to run from 1 - 4000

array_id=commandArgs(trailingOnly=T) #Assuming job array submission with bash/slurm
array_id=as.numeric(array_id)

# NN and GPD tuning parameters --------------------------------------------

fit_models = T

#NN layers/neurons 
nn.nunits = c(32,32,32)

#GPD quantile level 
quant.level = 0.85

#Set initial shape parameter
init_shape = 0.1

#Vector of decreasing learning rates
l.r.vec = c(0.001, 5e-04, 0.00025, 1e-04, 5e-05)  

#Link function for shape parameter
gpd.shape.custom_activation <- function(x) { #Bound to [-.5,.5]
  tf$atan(x)/pi
}

init_bias = tan(pi*(init_shape)) #Compute initial biases (to ensure valid starting parameters)

#NN activation function 
activation.func = "relu"

#Number of bootstraps 
n_boot = 100

n_y = 52 #(varies slightly depending on year)

# Packages and VE ---------------------------------------------------------

# Load some useful functions and packages
source("preamble.R")

# Load tensorflow/Keras and Virtual Environment 
reticulate::use_virtualenv("DeepSPAR_env", required = T)
packages = c("keras", "tensorflow","reticulate")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Set tensorflow seed
tensorflow::set_random_seed(array_id+1)

# Data --------------------------------------------------------------------

# Locations of stations 
river_names = c("ammer","iller","lech","loisach")

#Split into four 30 year windows
time_periods = list(1980:2009,2010:2039,2040:2069,2070:2099)

# Load in data and standardisations constants 
std_sigmas = readRDS(file = "RunOffData/std_sigmas.rds")

std_mus = readRDS(file = "RunOffData/std_mus.rds")

std_tw_dfs = readRDS(file = "RunOffData/std_tw_dfs.rds")

os_tw_dfs = readRDS(file = "RunOffData/os_tw_dfs.rds")

l1_seq_upper = 1000

l1_seq_lower = 30

prob_RPs = c(10)

# Fitting across all pairs ------------------------------------------------

indices_grid = expand.grid(1:n_boot,c(1:length(time_periods)))

print(paste("Processing task", array_id))

# Index of bootstrap iteration 
boot_ind = indices_grid$Var1[array_id]

# Index of time window. We use middle of the period for optimisation
time_window_ind = indices_grid$Var2[array_id]

# if(file.exists(file=paste0("BootVals/model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
#   stop("Script already finished")
# }

# Extract normalised data
X.norm = std_tw_dfs[[time_window_ind]]

# Resample rows in data randomly
set.seed(boot_ind+1)
X.norm = X.norm[sample(1:nrow(X.norm),nrow(X.norm),replace=T),]

polar = rect2polar(t(X.norm))  #Get polar coordinates
n = dim(X.norm)[1]
d = dim(X.norm)[2]

# Create AR decomp 
R = polar$r
W = X.norm/R
W = as.matrix(W)

# Make 10% validation data
valid.inds = sample(1:n, round(n/10))

# Build quantile regression model

R.train <- R[-valid.inds]
W.train <- W[-valid.inds,]
R.valid <- R[valid.inds]
W.valid <- W[valid.inds,]

input_pseudo_angles <- layer_input(shape = d, name = "input_pseudo_angles")

qBranch <- input_pseudo_angles %>%
  layer_dense(units = nn.nunits[1], 
              activation = activation.func, 
              name = "q_dense1", 
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))

for (i in 2:length(nn.nunits)) {
  
  qBranch <- qBranch %>%
    layer_dense(units = nn.nunits[i], 
                activation = activation.func, 
                name = paste0("q_dense", i), 
                kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))
  
}

qBranch <- qBranch %>%
  layer_dense(units = 1, 
              activation = "exponential", 
              name = "q_final", 
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))  

# Construct Keras model
quant.model <- keras_model(inputs = c(input_pseudo_angles), 
                           outputs = c(qBranch))


# Compile the model with the tilted loss and the adam optimiser
quant.model %>%
  compile(optimizer = "adam", 
          loss = tilted_loss, 
          run_eagerly = T)

if(fit_models){
  
  checkpoint <- callback_model_checkpoint(filepath = paste0("runs/QR_est/qr_fit_4d_","bootind_",boot_ind,"_tw_",time_window_ind), 
                                          monitor = "val_loss", 
                                          verbose = 0,
                                          save_best_only = TRUE, 
                                          save_weights_only = TRUE, 
                                          mode = "min", 
                                          save_freq = "epoch")
  n.epochs <- 500
  batch.size <- 1024
  
  dim(R)=c(length(R),1); dim(R.train)=c(length(R.train),1); dim(R.valid)=c(length(R.valid),1)
  
  history <- quant.model %>%
    fit(list(W.train), R.train, 
        epochs = n.epochs, 
        batch_size = batch.size, 
        callback = list(checkpoint,callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 5)), 
        validation_data = list(list(input_pseudo_angles = W.valid),R.valid))
  
  # Load the best fitting model from the checkpoint save
  quant.model <- load_model_weights_tf(quant.model, filepath = paste0("runs/QR_est/qr_fit_4d_","bootind_",boot_ind,"_tw_",time_window_ind))
  
  # Then save the best model.
  save_model_tf(quant.model, paste0("runs/QR_est/qr_fit_4d_","bootind_",boot_ind,"_tw_",time_window_ind))
  
}

quant.model <- tf$saved_model$load(paste0("runs/QR_est/qr_fit_4d_","bootind_",boot_ind,"_tw_",time_window_ind))

pred.quant <- k_get_value(quant.model(k_constant(W)))

# Compute threshold exceedances 
u <- pred.quant

# Create training and validation data
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
  layer_dense(units = nn.nunits[1], 
              activation = activation.func, 
              name = "gpd_dense1", 
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))

if (length(nn.nunits) > 2) {
  
  for (i in 2:(length(nn.nunits)-1) ) {
    
    gpd.Branch <- gpd.Branch %>%
      layer_dense(units = nn.nunits[i], 
                  activation = activation.func, 
                  name = paste0("gpd_dense", i), 
                  kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))
    
  }
  
}

scale.Branch <- gpd.Branch %>% 
  layer_dense(units = nn.nunits[length(nn.nunits)], 
              activation = activation.func, 
              name = paste0("scale_pre"), 
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04)) %>%
  layer_dense(units = 1, 
              activation = "exponential", 
              name = "scale_final", 
              weights = list(matrix(0,nrow = nn.nunits[length(nn.nunits)], ncol = 1), array(0, dim = c(1))), 
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))


shape.Branch <- gpd.Branch %>% 
  layer_dense(units = nn.nunits[length(nn.nunits)], 
              activation = activation.func, 
              name = paste0("shape_pre"), 
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04)) %>%
  layer_dense(units = 1, 
              activation = gpd.shape.custom_activation, 
              name = "shape_final", 
              weights = list(matrix(0,nrow = nn.nunits[length(nn.nunits)], ncol = 1), array(init_bias, dim = c(1))), 
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))

output <- layer_concatenate(c(scale.Branch, shape.Branch, input_u, input_pseudo_angles))

# Construct Keras model
GPD.model <- keras_model(inputs = c(input_pseudo_angles, input_u), outputs = output)
summary(GPD.model)

GPD.model %>%
  compile(optimizer = optimizer_adam(learning_rate = l.r.vec[1]), 
          loss = GPD_loss, 
          run_eagerly = T)

if(fit_models){
  
  n.epochs <- 750
  batch.size = length(R.train)
  
  checkpoint <- callback_model_checkpoint(filepath = paste0("runs/GPD_est_4d/tmp_",array_id,"/gpd_fit_{epoch:02d}"), 
                                          monitor = "val_loss", 
                                          verbose = 0, 
                                          save_best_only = FALSE, 
                                          save_weights_only = TRUE, 
                                          mode = "min",
                                          save_freq = "epoch")
  
  history <- GPD.model %>%
    fit(
      list(W.train, u.train), R.train,
      epochs = n.epochs,
      batch_size = batch.size,
      callbacks = list(
        checkpoint,
        callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 5)),
      validation_data = list(list(W.valid, u.valid), R.valid)
    )
  
  optimal_epoch = which.min(history$metrics$val_loss)
  
  GPD.model <- load_model_weights_tf(GPD.model,
                                     filepath = paste0("runs/GPD_est_4d/tmp_",array_id,"/gpd_fit_",sprintf("%02d", optimal_epoch-1)))
  
  save_model_weights_tf(GPD.model,  
                        paste0("runs/GPD_est_4d/tmp_",array_id,"/gpd_fit_best"))
  
  save_model_tf(GPD.model, 
                paste0("runs/GPD_est_4d/gpd_fit_","bootind_",boot_ind,"_tw_",time_window_ind))
  
  # Load best validation loss value 
  best.valid.loss <- history$metrics$val_loss[which.min(history$metrics$val_loss)-1]
  
  # Adjust number of epochs for further optimisation 
  n.epochs <- 200
  
  # Decreasing learning rate training scheme
  for (k in 2:length(l.r.vec)) {
    
    print(paste0("Decreasing learning rate size to: ", l.r.vec[k]))
    
    GPD.model %>%
      compile(optimizer = optimizer_adam(learning_rate = l.r.vec[k]), 
              loss = GPD_loss,
              run_eagerly = T)
    
    history <- GPD.model %>%
      fit(list(W.train, u.train), R.train,
          epochs = n.epochs,
          batch_size = batch.size,
          callback = list(checkpoint,callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 10)),
          validation_data = list(list(input_pseudo_angles = W.valid,input_u = u.valid), R.valid))
    
    # Check to see if validation loss has decreased
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
      
      GPD.model <- load_model_weights_tf(GPD.model,
                                         filepath = paste0("runs/GPD_est_4d/tmp_",array_id,"/gpd_fit_",sprintf("%02d", optimal_epoch)))
      
      save_model_weights_tf(GPD.model,  
                            paste0("runs/GPD_est_4d/tmp_",array_id,"/gpd_fit_best"))
      
      save_model_tf(GPD.model, 
                    paste0("runs/GPD_est_4d/gpd_fit_","bootind_",boot_ind,"_tw_",time_window_ind))
      
    } else{
      
      load_model_weights_tf(GPD.model,  
                            paste0("runs/GPD_est_4d/tmp_",array_id,"/gpd_fit_best"))
      
    }
    
  }
  
  
  unlink(paste0("runs/GPD_est_4d/tmp_",array_id,""),recursive = T, force = T)
  
}

GPD.model<- load_model_tf(filepath = paste0("runs/GPD_est_4d/gpd_fit_","bootind_",boot_ind,"_tw_",time_window_ind),
                          custom_objects = list(GPD_loss=GPD_loss))

#Validation 

#Evaluating model fit on thresh exc
pred.GPD.exc <- k_get_value(GPD.model(list(k_constant(W[R-u>0,]), k_constant(matrix(u[R-u>0],ncol=1)))))

obs_quants = qexp(apply(cbind(u[R-u>0], pred.GPD.exc[,1:2], R[R-u>0]), 1, function(x) pgpd(x[4], loc = x[1], scale = x[2]/(1 +x[3]), x[3])))

ps = readRDS(file = paste0("BootVals/gpd_qq_prob_tw_",time_window_ind,".rds"))

boot_quants = quantile(obs_quants,ps)

saveRDS(boot_quants,file = paste0("BootVals/gpd_qq_quants_","bootind_",boot_ind,"_tw_",time_window_ind,".rds"))

n.sim.tail = 2e6

n.sim = round((n.sim.tail/(100*(1-quant.level)))*100)

X.norm.body = X.norm[which(pred.quant >= R),]

X.sim.body = X.norm.body[sample(1:nrow(X.norm.body),n.sim-n.sim.tail,replace=T),]

X.norm.excess = X.norm[which(pred.quant < R),]

W.excess = W[which(pred.quant < R),]

W.sim = W.excess[sample(1:nrow(W.excess), n.sim.tail,replace=T),]
u.sim <- k_get_value(quant.model(k_constant(W.sim)))

dim(u.sim) = c(length(u.sim), 1)

pred.GPD.sim <- k_get_value(GPD.model(list(k_constant(W.sim), k_constant(u.sim))))
nu.sim = pred.GPD.sim[, 1]
xi.sim = pred.GPD.sim[, 2]

unif_sample = runif(n.sim.tail,0,1)
R.sim = apply(cbind(unif_sample,u.sim,nu.sim,xi.sim),1,function(x)
  qgpd(x[1], loc=x[2], scale = x[3]/(1+ x[4]), shape=x[4])
)

X.sim.excess = R.sim*W.sim

X.sim = rbind(X.sim.body,X.sim.excess)

rm(u.sim,pred.GPD.sim,unif_sample,nu.sim,xi.sim,R.sim,W.sim)

mus = std_mus[[time_window_ind]]

sigmas = std_sigmas[[time_window_ind]]

X.excess.obs.os = apply(rbind(mus,X.norm.excess),2,function(x){return( x[-1] + x[1] )})

X.excess.obs.os = apply(X.excess.obs.os,2,function(x){return( log(1 + exp(x)) )})

X.excess.obs.os = apply(rbind(sigmas,X.excess.obs.os),2,function(x){return( x[-1]*x[1] )})

X.excess.sim.os = apply(rbind(mus,X.sim.excess),2,function(x){return( x[-1] + x[1] )})

X.excess.sim.os = apply(X.excess.sim.os,2,function(x){return( log(1 + exp(x)) )})

X.excess.sim.os = apply(rbind(sigmas,X.excess.sim.os),2,function(x){return( x[-1]*x[1] )})

labels <- river_names

tail_dep <- function(u,x, y) {
  qx <- quantile(x, u, na.rm = TRUE)
  qy <- quantile(y, u, na.rm = TRUE)
  mean(x > qx & y > qy, na.rm = TRUE) / (1 - u)
}

u_vec = seq(0.8,0.99,by = 0.01)

p <- ncol(X.excess.sim.os)


for (j in 1:p) {
  for (i in 1:p) {

    if (i == j) {

      ps = readRDS(file = paste0("BootVals/marg_qq_prob_",labels[j],"_tw_",time_window_ind,".rds"))

      sim_quants = quantile(X.excess.sim.os[, j],ps)

      saveRDS(sim_quants,file = paste0("BootVals/marg_qq_quants_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))

    } else if (j > i) {

      NULL

    } else {

      xs <- X.excess.sim.os[, i]; ys <- X.excess.sim.os[, j]

      lam_sim <- sapply(u_vec,tail_dep,x=xs,y=ys)

      saveRDS(lam_sim,file = paste0("BootVals/chi_plot_",labels[i],"_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))

    }
  }
}

for(i in 1:d){

  title = river_names[i]

  fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))

  fix_probs = 1/(fix_RPs*n_y)

  model_RLs = quantile(X.sim[,i],1-fix_probs)

  model_RLs = log(1+exp(model_RLs+mus[i]))*sigmas[i]

  saveRDS(model_RLs,file = paste0("BootVals/rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))

}

for(i in 1:d){

  title = river_names[i]

  fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))

  fix_probs = 1/(fix_RPs*n_y)

  model_RLs = quantile(X.sim[,i],fix_probs)

  model_RLs = log(1+exp(model_RLs+mus[i]))*sigmas[i]

  saveRDS(model_RLs,file = paste0("BootVals/rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))

}

est_RLs_upper = matrix(NA,ncol=d,nrow=length(prob_RPs))

est_RLs_lower = matrix(NA,ncol=d,nrow=length(prob_RPs))

for(i in 1:d){

  #Upper tail return levels

  est_RLs_upper[,i] = quantile(X.norm[,i],1-(1/(n_y*prob_RPs)))

  est_RLs_upper[,i] = log(1+exp(est_RLs_upper[,i]+mus[i]))*sigmas[i]

  #Lower tail return levels

  est_RLs_lower[,i] = quantile(X.norm[,i],(1/(n_y*prob_RPs)))

  est_RLs_lower[,i] = log(1+exp(est_RLs_lower[,i]+mus[i]))*sigmas[i]

}

X.sim = apply(rbind(mus,X.sim),2,function(x){return( x[-1] + x[1] )})

X.sim = apply(X.sim,2,function(x){return( log(1 + exp(x)) )})

X.sim = apply(rbind(sigmas,X.sim),2,function(x){return( x[-1]*x[1] )})

l1_upper_ests = sapply(l1_seq_upper,function(l1_q){

  return(l1_prob_upper(u=l1_q,sim_data=X.sim))

})

l1_lower_ests = sapply(l1_seq_lower,function(l1_q){

  return(l1_prob_lower(u=l1_q,sim_data=X.sim))

})

upper_surv_ests = apply(est_RLs_upper,1,function(ret_levels){

  return(survivor_prob(u_vec=ret_levels,sim_data=X.sim))

})

lower_cdf_ests = apply(est_RLs_lower,1,function(ret_levels){

  return(cdf_prob(u_vec=ret_levels,sim_data=X.sim))

})

probs_list = list(l1_seq_upper = l1_seq_upper, l1_seq_lower = l1_seq_lower, est_RLs_upper=est_RLs_upper, est_RLs_lower=est_RLs_lower, l1_upper_ests=l1_upper_ests, l1_lower_ests=l1_lower_ests, upper_surv_ests=upper_surv_ests, lower_cdf_ests=lower_cdf_ests)

saveRDS(probs_list,file=paste0("BootVals/model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))

