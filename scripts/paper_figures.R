# Clear R environment 

rm(list=ls())

# NN and GPD tuning parameters --------------------------------------------

# NN hidden layers + neurons 

nn.nunits <- c(16,16,16)

# GPD quantile level 

quant.level <- 0.8
alpha <- 1 - quant.level

# Initial shape parameter

init_shape <- 0.01

# Vector of decreasing learning rates for SGD

l.r.vec <- c(0.001, 5e-04, 0.00025, 1e-04, 5e-05)  

# Link function for shape parameter

gpd.shape.custom_activation <- function(x) { 
  tf$atan(x)/pi
}

# Compute bias for initialising GPD MLP (to ensure valid starting parameters)

init_bias <- tan(pi*(init_shape)) 

# Specify activation function 

activation.func <- "relu"

# Flags for specifying whether or not to fit MLPs associated with quantile regression and GPD. 
# Can be set to 'F' if NN parameters already computed 

fit_QR <- T
fit_GPD <- T

# Packages and VE ---------------------------------------------------------

# Load some useful functions and packages
source("preamble.R")

# Load tensorflow/Keras and Virtual Environment 
reticulate::use_virtualenv("DeepSPAR_env", required = T)
packages <- c("keras", "tensorflow","reticulate")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Set tensorflow seed
tensorflow::set_random_seed(1)

# Bivariate Gaussian data ---------------------------------------------------

# Simulate Gaussian copula data on standard Laplace margins 

library(copula)

Laplace_inverse <- function(u){ 
  x = c()
  x[u<=0.5] = log(2*u[u<=0.5])
  x[u>0.5] = -log(2*(1-u[u>0.5]))
  return(x)
}

set.seed(1)

n <- 30000

rho <- 0.5

normc <- normalCopula(param = rho, dim = 2)

sim_data_normal <- apply(rCopula(n, copula = normc),2,Laplace_inverse)

# Decomposing the data into angular and radial components ------------------------------------------------

# Specify data matrix centered at (0,...,0)

X.sim <- sim_data_normal

# Specify sample size and dimension of data 

n <- dim(X.sim)[1]
d <- dim(X.sim)[2]

# Decompose data into radial and angular components 

R <- apply(X.sim,1,l2_norm)

W <- X.sim/R
W <- as.matrix(W)

# Fit quantile regression model  ------------------------------------------

# Compute indices of 10% validation set

valid.inds <- sample(1:n, round(n/10))

# Compute subsets of training and validation data sets

R.train <- R[-valid.inds]
W.train <- W[-valid.inds,]
R.valid <- R[valid.inds]
W.valid <- W[valid.inds,]

# Build quantile regression model via Keras 

# Specify input layer 

input_pseudo_angles <- layer_input(shape = d, name = "input_pseudo_angles")

# Specify hidden layers

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

# Specify final layer 

qBranch <- qBranch %>%
  layer_dense(units = 1,
              activation = "exponential",
              name = "q_final",
              kernel_regularizer = regularizer_l1_l2(l1 = 1e-04,l2 = 1e-04))

# Construct Keras model

quant.model <- keras_model(inputs = c(input_pseudo_angles),
                           outputs = c(qBranch))


# Compile the model with the tilted loss and the ADAM optimiser

quant.model %>%
  compile(optimizer = "adam",
          loss = tilted_loss,
          run_eagerly = T)

# Fit the quantile regression model unless fit_QR set to 'F'

if(fit_QR){
  
  # Create a temporary directory to store the checkpoint. This directory is later removed after we have saved the optimal model weights
  
  checkpoint_dir <- tempdir()
  checkpoint_path <- file.path(checkpoint_dir, "QR_checkpoint")
  
  # Create a checkpoint for fitting the model 

  checkpoint <- callback_model_checkpoint(filepath = checkpoint_path,
                                          monitor = "val_loss",
                                          verbose = 0,
                                          save_best_only = TRUE,
                                          save_weights_only = TRUE,
                                          mode = "min",
                                          save_freq = "epoch")
    
  # Specify number of epochs and batch size 
  
  n.epochs <- 500
  batch.size <- 1024

  # Convert radial vectors to matrices - this is needed for Keras to work
  
  dim(R) <- c(length(R),1)
  dim(R.train) <- c(length(R.train),1)
  dim(R.valid) <- c(length(R.valid),1)

  # Run ADAM to estimate parameters with early stopping 
  
  history <- quant.model %>%
    fit(list(W.train), R.train,
        epochs = n.epochs,
        batch_size = batch.size,
        callback = list(checkpoint,callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 5)),
        validation_data = list(list(input_pseudo_angles = W.valid),R.valid))

  # Load the best model weights from the saved checkpoint 
  
  quant.model <- load_model_weights_tf(quant.model, filepath = checkpoint_path)
  
  # Save the best tf model
  
  save_model_tf(quant.model, "../nn_parameters/QR_est/paper_figures")
  
  unlink(checkpoint_dir, recursive = TRUE)
  
}

# Load in best fitting tf model 

quant.model <- tf$saved_model$load("../nn_parameters/QR_est/paper_figures")

# Compute estimated quantile over all observed angles 

u <- k_get_value(quant.model(k_constant(W)))

# As a quick check, we compute the empirical quantile probability. This should be approximately equal to the corresponding theoretical value 

print(mean(R>u))
print(alpha)
 
# Fit GPD to threshold exceeding observations ------------------------------------------

# Create training and validation data sets for GPD model 
# Note that we only compute positive threshold exceedances 

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

# Again, convert vectors to matrices to ensure it works with Keras 

dim(R.train) <- c(length(R.train), 1)
dim(R.valid) <- c(length(R.valid), 1)

dim(u.train) <- c(length(u.train), 1)
dim(u.valid) <- c(length(u.valid), 1)

# Build GPD model via Keras 

# Specify input layers - we have two input vectors in this case

input_pseudo_angles <- layer_input(shape = d, name = "input_pseudo_angles")

input_u <- layer_input(shape = d, name = "input_u")

# Specify hidden layers

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

# Specify final layers - note we have seperate final branches for the scale and shape parameters, owing to the fact these parameters are independent 
# This is hard parameter sharing for multi-task learning 
# Also note that we initialise the weights/biases in these layers to ensure a valid starting solution

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

# Construct Keras model

output <- layer_concatenate(c(scale.Branch, shape.Branch, input_u, input_pseudo_angles))
GPD.model <- keras_model(inputs = c(input_pseudo_angles, input_u), outputs = output)

# Compile the model with the GPD log-likelihood loss and the ADAM optimiser

GPD.model %>%
  compile(optimizer = optimizer_adam(learning_rate = l.r.vec[1]),
          loss = GPD_loss,
          run_eagerly = T)

# Fit the GPD model unless fit_GPD set to 'F'

if(fit_GPD){
  
  # Create a checkpoint for fitting the model 
  
  checkpoint <- callback_model_checkpoint(filepath = paste0("../nn_parameters/GPD_est/tmp_paper_figures/gpd_fit_{epoch:02d}"),
                                          monitor = "val_loss",
                                          verbose = 0,
                                          save_best_only = FALSE,
                                          save_weights_only = TRUE,
                                          mode = "min",
                                          save_freq = "epoch")

  # Specify number of epochs and batch size 
  # Note the batch size is simply the whole training data set since this is relatively small 
  
  n.epochs <- 750
  batch.size <- length(R.train)
  
  # Run ADAM to estimate parameters with early stopping
  
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
  
  # Check which checkpoint/training step gave the best model loss 

  optimal.epoch <- which.min(history$metrics$val_loss)

  # Load the best model weights from the model fit just before the optimal epoch 
  
  GPD.model <- load_model_weights_tf(GPD.model,
                                     filepath = paste0("../nn_parameters/GPD_est/tmp_paper_figures/gpd_fit_",sprintf("%02d", optimal.epoch-1)))

  # Save the best weights and model at the (near) optimal loss value
  
  save_model_weights_tf(GPD.model,
                        paste0("../nn_parameters/GPD_est/tmp_paper_figures/gpd_fit_best"))

  save_model_tf(GPD.model,
                paste0("../nn_parameters/GPD_est/paper_figures"))

  # Get the loss value at the epoch before the optimal epoch 
  
  best.valid.loss <- history$metrics$val_loss[which.min(history$metrics$val_loss)-1]

  # Reset the number of epochs; we reduce here since we are already close to an 'optimal' solution 
  
  n.epochs <- 200

  # In the next step, we reduce the learning rate iteratively to see if we can find optimise our solution even more 
  # For each new learning rate, we load in the best current weights, re-optimise the parameters and check if we improve on loss
  # If the loss decreases, we update the model weights
  # If there is no improvement, we keep the existing weights 
  
  for (k in 2:length(l.r.vec)) {

    # Adjusting optimiser to use smaller learning rate
    
    GPD.model %>%
      compile(optimizer = optimizer_adam(learning_rate = l.r.vec[k]),
              loss = GPD_loss,
              run_eagerly = T)

    # Re-applying algorithm with lower learning rate, with existing weights as starting solution 
    
    history <- GPD.model %>%
      fit(
        list(W.train, u.train), R.train,
          epochs = n.epochs,
          batch_size = batch.size,
          callback = list(
            checkpoint,
            callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 10)),
          validation_data = list(list(input_pseudo_angles = W.valid,input_u = u.valid), R.valid)
        )

    # We check to see if validation loss has decreased
    # The first condition check if the weights give an infeasible solution, resulting in an infinite loss value 
    
    if (length(which.min(history$metrics$val_loss)) > 0) {

      # If the updated loss does not improve on the current loss, we do not update our solution
      
      if (min(history$metrics$val_loss, na.rm = T) > best.valid.loss) {

        optimal.epoch <- NA

      } else {
        
        # Find the epoch and the loss value at the optimal weights 

        optimal.epoch <- which.min(history$metrics$val_loss)

        best.valid.loss <- history$metrics$val_loss[optimal.epoch]

      }

    } else {

      # When the computed weights are infeasible, we do not update our solution 
      
      optimal.epoch <- NA

    }

    # If the loss improves with the new learning rate, we update the saved model weights 
    
    if (!is.na(optimal.epoch)) {

      # Load in weights of new 'optimal' solution 
      
      GPD.model <- load_model_weights_tf(GPD.model,
                                         filepath = paste0("../nn_parameters/GPD_est/tmp_paper_figures/gpd_fit_",sprintf("%02d", optimal.epoch)))

      # Update 'optimal' weights 
      
      save_model_weights_tf(GPD.model,
                            paste0("../nn_parameters/GPD_est/tmp_paper_figures/gpd_fit_best"))

      # Update saved model with new weights  
      
      save_model_tf(GPD.model,
                    paste0("../nn_parameters/GPD_est/paper_figures"))

    } else{

      # Re-load in best weights from previous learning rate  
      
      load_model_weights_tf(GPD.model,
                            paste0("../nn_parameters/GPD_est/tmp_paper_figures/gpd_fit_best"))

    }

  }

  # As a final step to keep the directory neat, we delete all temporary saved model weights over different epochs. 
  # We only save our 'optimal' model fit

  unlink(paste0("../nn_parameters/GPD_est/tmp_paper_figures"),recursive = T, force = T)
  
}

# Load in our best model fit 

GPD.model <- load_model_tf(filepath = paste0("../nn_parameters/GPD_est/paper_figures"),
                          custom_objects = list(GPD_loss=GPD_loss))


# Computing diagnostics using the fitted model  ---------------------------

# First we compute the GPD QQ plot 

# Estimate the GPD parameters at all threshold exceeding observations 

pred.GPD.exc <- k_get_value(GPD.model(list(k_constant(W[R-u>0,]), k_constant(matrix(u[R-u>0],ncol=1)))))

# Transform threshold exceeding radial values onto the standard exponential scale with the probability integral transform 

obs.quants <- qexp(apply(cbind(u[R-u>0], pred.GPD.exc[,1:2], R[R-u>0]), 1, function(x) pgpd(x[4], loc = x[1], scale = x[2]/(1 +x[3]), x[3])))

# Compute the corresponding theoretical quantiles, with theoretical uncertainty bounds computed from the Beta distribution (i.e., the distribution of rank statistics) 

n_p <- length(obs.quants)
ps <- (1:n_p)/(n_p + 1)
theor.quants <- qexp(ps)
lower.bound <- qexp(qbeta(0.025, 1:n_p, rev(1:n_p)))
upper.bound <- qexp(qbeta(0.975, 1:n_p, rev(1:n_p)))

# We save the GPD QQ plot in a png file 

png(
  filename = paste0("../figures/example_gpd_qq_plot.png"),
  width = 2500, height = 2500, res = 400
)

par(mgp = c(2.8, 0.7, 0))

plot(
  theor.quants, sort(obs.quants),
  pch = 16, col = "grey",
  cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.8,
  xlim = range(obs.quants, theor.quants,lower.bound,upper.bound),
  ylim = range(obs.quants, theor.quants,lower.bound,upper.bound),
  ylab = "SPAR", xlab = "Theoretical",
  main = paste0("GPD QQ plot")
)
polygon(c(rev(lower.bound), upper.bound),c(rev(sort(obs.quants)), sort(obs.quants)), col = 'grey80', border = NA)
abline(a = 0, b = 1, lwd = 4, col = 2)
points(theor.quants, sort(obs.quants), pch = 16, col = "black", cex = 1.4)

dev.off()

# Next, we compute a scatterplot using data simulated from our model 

# Specify how many observations to simulate from SPAR in the joint tail 

n.sim.tail <- 1e6

# Specify how many observations to resample from the body of data 

n.sim <- round(n.sim.tail/alpha)

# Specify non-threshold exceeding observations 

X.body <- X.sim[which(u >= R),]

# Specify threshold exceeding observations 

X.tail <- X.sim[which(u < R),]

# Resample with replacement from joint body 

X.body.resam <- X.body[sample(1:nrow(X.body),n.sim-n.sim.tail,replace=T),]

# Compute angles at threshold exceeding values 

W.excess <- W[which(u < R),]

# Resample with replacements from threshold exceeding angles 

W.sim <- W.excess[sample(1:nrow(W.excess), n.sim.tail,replace=T),]

# Compute quantile function at threshold exceeding angles 

u.sim <- k_get_value(quant.model(k_constant(W.sim)))

# Specify dimensions (for Keras)

dim(u.sim) <- c(length(u.sim), 1)

# Compute GPD parameters at threshold exceeding angles 

pred.GPD.sim <- k_get_value(GPD.model(list(k_constant(W.sim), k_constant(u.sim))))

# Simulate uniform observations 

unif_sample <- runif(n.sim.tail,0,1)

# Invert GPD to compute GPD quantiles 

R.sim <- apply(cbind(unif_sample,u.sim,pred.GPD.sim[,1:2]),1,function(x)
  qgpd(x[1], loc=x[2], scale = x[3]/(1+ x[4]), shape=x[4])
)

# Compute simulated observations on Cartesian scale 

X.tail.sim <- R.sim*W.sim

# Combine joint body and tail simulated observations 

X.spar.sim <- rbind(X.body.resam,X.tail.sim)

# Remove any objects we no longer need to save memory

rm(u.sim,pred.GPD.sim,unif_sample,R.sim,W.sim,X.body,W.excess)

# Vector of threshold for evaluating chi 

u_vec <- seq(0.8,0.99,by = 0.01)

# Specify number of variables 

p <- ncol(X.sim)

png(
  file = paste0("../figures/example_scatterplot.png"),
  width = 3500, height = 3500, res = 400
)

# 4x4 layout, all panels square, equal size
layout(matrix(1:4, nrow = 2, ncol = 2), widths = rep(1,2), heights = rep(1,2))
par(mar = c(4, 4.5, 2, 1))  # inner margin (bottom, left, top, right)
par(oma = c(1, 1, 4, 1))  # outer margin for overall labels / title

for (j in 1:p) {
  for (i in 1:p) {
    
    if (i == j) {
      
      # ---------------- Diagonal: QQ plot ----------------
      m <- min(c(length(X.tail.sim[, i]),length(X.tail[,i])))
      
      ps <- (1:m)/(m+1)
      
      sim.quants <- quantile(X.tail.sim[, i],ps)
      
      obs.quants <- quantile(X.tail[,i],ps)
      
      plot(
        sim.quants,obs.quants,
        pch=16,col="grey",
        cex.lab = 1.3, cex.axis = 1, cex.main = 1.3,
        xlim=range(obs.quants,sim.quants),
        ylim=range(obs.quants,sim.quants),
        ylab = ifelse(i==1,"Simulated data","") ,xlab = ifelse(i==1,"SPAR",""),
        main = bquote(X[.(i)]))
      abline(a=0,b=1,lwd=4,col=2)
      points(sim.quants,obs.quants,pch=16,col="black",cex=1.3)
      
    } else if (j > i) {
      # ---------------- Upper diagonal: Scatter ----------------
      xs <- X.tail.sim[, j]; ys <- X.tail.sim[, i]
      xo <- X.tail[, j]; yo <- X.tail[, i]
      
      if(i == 1){xlab=bquote(X[.(j)])}else{xlab=""}

      if(i == j-1){ylab=bquote(X[.(i)])}else{ylab=""}
      
      plot(xs[1:length(xo)],ys[1:length(xo)], pch = 16, col = rgb(0, 1, 0, 0.35), cex = 0.7,xlab=xlab,ylab=ylab,cex.lab=1.3)
      points(xo,yo, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7)
      
    } else {
      # ---------------- Lower diagonal: chi dependence coefficient on scale ----------------
      xs <- X.tail.sim[, i]; ys <- X.tail.sim[, j]
      xo <- X.tail[, i]; yo <- X.tail[, j]
      
      lam_sim <- sapply(u_vec,tail_dep,x=xs,y=ys)
      lam_obs <- sapply(u_vec,tail_dep,x=xo,y=yo)
      
      plot(u_vec,lam_obs,type="l",lwd=2,col="grey2",cex.lab=1.3, xlab=expression(u),ylab=expression(chi(u)),ylim=c(0,1))
      lines(u_vec,lam_sim,lwd=2,col="green3")
      lines(u_vec,lam_obs,type="l",lwd=2,col="grey")

      abline(h = 1,lwd=1,lty=2)

      abline(h = 0,lwd=1,lty=2)
      
    }
  }
}

dev.off()

# Finally, we compute return level plots for the upper and lower tails of both marginal variable  

# Specify number of observations per year. Given weekly maxima, this would be 52 (for the climate simulations)

n_y <- 52

# Specify fixed return periods for computing return levels from model 

fix_RPs <- exp(seq(log(10),log(1000),length.out = 3000))

# Compute upper tail probabilities associated with fixed return periods 

fix_probs <- 1/(fix_RPs*n_y)

# We save both return level plots as pngs

png(
  filename = paste0("../figures/example_upper_tail_rl_plot.png"),
  width = 2800, height = 1400, res = 400
)

par(mar = c(4, 4.5, 4, 2))

par(mgp = c(2.8, 0.7, 0))

par(mfrow=c(1,2))

for(i in 1:p){
  
  # Use rank transform to workout empirical return periods 
  
  emp_probs <- 1 - rank(X.sim[,i])/(nrow(X.sim)+1)
  
  emp_RPs <- 1/(n_y*emp_probs)
  
  # Just consider empirical return periods above 10 years
  
  emp_RLs <- X.sim[,i][emp_RPs>=10]
  
  emp_RPs <- emp_RPs[emp_RPs>=10]
  
  # Compute model return levels using data simulated from the SPAR model 
  
  model_RLs <- quantile(X.spar.sim[,i],1-fix_probs)
  
  plot(
    fix_RPs, model_RLs,
    log = "x",
    type = "l", lwd=3, col = "grey",
    xlab = "Return period (years)",
    ylab = "Return level",
    xaxt = "n",
    cex.lab=1.3,cex.main=1.5,
    ylim=range(emp_RLs,model_RLs),
    main=bquote(X[.(i)] ~ "upper tail")
  )
  
  lines(fix_RPs, model_RLs, lwd=3, col = "grey3")
  
  points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)
  
  ticks <- c(10, 100, 1000)
  
  axis(1, at = ticks, labels = as.character(ticks))

}

dev.off()

png(
  filename = paste0("../figures/example_lower_tail_rl_plot.png"),
  width = 2800, height = 1400, res = 400
)

par(mar = c(4, 4.5, 4, 2))

par(mgp = c(2.8, 0.7, 0))

par(mfrow=c(1,2))

for(i in 1:p){
  
  # Use rank transform to workout empirical return periods 
  
  emp_probs <- rank(X.sim[,i])/(nrow(X.sim)+1)
  
  emp_RPs <- 1/(n_y*emp_probs)
  
  # Just consider empirical return periods above 10 years
  
  emp_RLs <- X.sim[,i][emp_RPs>=10]
  
  emp_RPs <- emp_RPs[emp_RPs>=10]
  
  # Compute model return levels using data simulated from the SPAR model 
  
  model_RLs <- quantile(X.spar.sim[,i],fix_probs)
  
  plot(
    fix_RPs, model_RLs,
    log = "x",
    type = "l", lwd=3, col = "grey",
    xlab = "Return period (years)",
    ylab = "Return level",
    xaxt = "n",
    cex.lab=1.3,cex.main=1.5,
    ylim=range(emp_RLs,model_RLs),
    main=bquote(X[.(i)] ~ "lower tail")
  )
  
  lines(fix_RPs, model_RLs, lwd=3, col = "grey3")
  
  points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)
  
  ticks <- c(10, 100, 1000)
  
  axis(1, at = ticks, labels = as.character(ticks))
  
}

dev.off()

# Example SPAR illustration figures from paper ----------------------------

# We now recreate the illustrative figures from the article that help to explain the SPAR model formulation/paradigm  

# First, we show what discharge data could look like prior to pre-processing 

png(
  filename = paste0("../figures/SPAR_illustration_preprocess.png"),
  width = 2000, height = 2000, res = 400
)

par(mar = c(4, 4.2, 3, 1), mgp = c(2.5, 0.7, 0))

par(mfrow=c(1,1))

# We reverse the pre-processing step using an arbitrary geometric mean and standard deviations 

X.obs = apply(rbind(c(2,4),X.sim),2,function(x){return( x[-1] + x[1] )})

X.obs = apply(X.obs,2,function(x){return( log(1 + exp(x)) )})

X.obs = apply(rbind(c(10,20),X.obs),2,function(x){return( x[-1]*x[1] )})

plot(X.obs, 
     pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,
     xlab=expression(X[1]^o),ylab=expression(X[2]^o),cex.lab=1.3,
     xlim=range(X.obs[,1]),ylim=range(X.obs[,2]),
     frame.plot=F)

mtext("(a)", side = 3, line = 1, adj = 0, cex = 1.5)

points(log(1 + exp(2))*10,log(1 + exp(4))*20,col="blue",cex=1,pch=16)

dev.off()

# Next we plot the simulated data 'after' pre-processing (i.e., centered at (0,0)). We also illustrate the angular-radial decomposition 

png(
  filename = paste0("../figures/SPAR_illustration_1.png"),
  width = 2000, height = 2000, res = 400
)

par(mar = c(4, 4.2, 3, 1), mgp = c(2.5, 0.7, 0))

par(mfrow=c(1,1))

plot(X.sim, 
     pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,
     xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,
     xlim=c(-10,10),ylim=c(-10,10),
     asp=1,frame.plot=F)

mtext("(b)", side = 3, line = 1, adj = 0, cex = 1.5)

# Compute radius and bivaraite angle for a single point 

ind = 22246

phi = atan2(X.sim[ind,2],X.sim[ind,1])

r = sqrt(X.sim[ind,2]^2+X.sim[ind,1]^2)

# Illustrate point with corresponding radius and angle, alongside hypersphere (i.e., unit circle) 

arrows(0, 0, r*cos(phi), r*sin(phi), length = 0.1, col = "green3", lwd = 2.5)

points(X.sim[ind,1],X.sim[ind,2],col="orange",cex=1,pch=16)

points(0,0,col="blue",cex=1,pch=16)

pred_phis = seq(0,2*pi,len=201)

lines(cbind(cos(pred_phis),sin(pred_phis)),lwd=1.5,col=2)

points(cos(phi),sin(phi),col="purple",cex=1,pch=16)

legend(
  "topleft",
  legend = c(expression(paste("(", 0, ",", 0, ")")), 
             expression(paste("(", x[1], ",", x[2], ")")), 
             expression(paste("(", w[1], ",", w[2], ")")), 
             expression(r)),
  col = c("blue", "orange", "purple", "green3"),  # Colors for each legend item
  pch = c(16, 16, 16, NA),                # Points for the first three, no point for the last
  lty = c(0, 0, 0, 1),                    # No line for the first three, line for the last
  lwd = c(NA, NA, NA, 2.5)                # Line width for the last item
)

dev.off()

# We illustrate the threshold exceeding observations in the plot below (the joint tail)

# We compute the quantile set on a dense grid of angles 

W.circle <- cbind(cos(pred_phis),sin(pred_phis))

u.circle <- k_get_value(quant.model(k_constant(W.circle)))

quant.set <- c(u.circle)*W.circle

png(
  filename = paste0("../figures/SPAR_illustration_2.png"),
  width = 2000, height = 2000, res = 400
)

par(mar = c(4, 4.2, 3, 1), mgp = c(2.5, 0.7, 0))

par(mfrow=c(1,1))

# Plot observations in the joint body

plot(X.body,
     pch = 16,asp=1, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,
     xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,
     xlim=c(-10,10),ylim=c(-10,10),
     frame.plot=F)

# Plot observations in the joint tail 

points(X.tail,pch=16,col = rgb(0, 1, 0, 1),cex=0.7)

# Plot the quantile set 

lines(quant.set,lwd=3,col="red")

mtext("(c)", side = 3, line = 1, adj = 0, cex = 1.5)

dev.off()

# We illustrate example tail regions for which one can perform inference via the SPAR model 

# Compute the componentwise maxima of the quantile set 

quant.set.maxima <- apply(quant.set,2,max)

png(
  filename = paste0("../figures/SPAR_illustration_3.png"),
  width = 2000, height = 2000, res = 400
)

par(mar = c(4, 4.2, 3, 1), mgp = c(2.5, 0.7, 0))

par(mfrow=c(1,1))

# Plot observations in the joint tail 

plot(X.tail, 
     pch = 16,asp=1, col = rgb(0, 1, 0, 1), cex = 0.7,
     xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,
     xlim=c(-10,10),ylim=c(-10,10),
     frame.plot=F)

# Plot quantile set 

lines(quant.set,lwd=3,col="red")

mtext("(d)", side = 3, line = 1, adj = 0, cex = 1.5)

# Tail region 1

rect(quant.set.maxima,-15,15,15,col = rgb(0.5, 0.5, 0.5, 0.15),border = "black",lwd=2,lty=2)

text(7,-2,label = expression(bold(X[1] > q[1]^M)),cex=1.4,col="black",font=2)

# Tail region 2

rect(-15, 4, -2, 15, col = rgb(0, 0, 1, 0.15), border = "blue", lwd = 2, lty = 2)
text(-6.5,7,label = expression(bold(X[1] < -2) * "," * bold(X[2] > 4)),cex=1.4,col="blue")

# Tail region 3

x_poly <- c(-26, 20, -26)
y_poly <- c(20,-26,-26)

polygon(x_poly, y_poly, col = rgb(1, 0.647, 0, 0.15), border = "orange",lwd = 2, lty = 2)
text(-7,-8,label = expression(bold(X[1] + X[2] < -6)),cex=1.4,col="orange")

dev.off()

# Finally, we illustrate a range of example events sets that one can simulate via the SPAR model 

# We just look at the first 10000 simulated observations (to avoid plotting a very large number of points)

X.tail.sim2 = X.tail.sim[1:10000,]

png(
  filename = paste0("../figures/SPAR_event_sets.png"),
  width = 4500, height = 1500, res = 400
)

par(mgp = c(2.8, 0.7, 0))

par(mfrow=c(1,3))

# Event set 1: X_1 greater than 5

plot(X.tail.sim2,
     pch=16,col = rgb(0, 1, 0, 1),cex=0.7,
     xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,
     xlim=c(-12,12),ylim=c(-12,12))

lines(quant.set,lwd=3,col="red")

rect(5,-15,15,15,col = rgb(0, 0, 1, 0.15),border = "blue",lwd=2,lty=2)

text(8,-3,label = expression(bold(X[1] > 5)),cex=1.6,col="blue",font=2)

points(X.tail.sim2[X.tail.sim2[,1]>5,],pch=16,col="blue", cex = 0.7)

# Event set 2: X_1 + X_2 > 6

plot(X.tail.sim2,
     pch=16,col = rgb(0, 1, 0, 1),cex=0.7,
     xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,
     xlim=c(-12,12),ylim=c(-12,12))

lines(quant.set,lwd=3,col="red")

x_poly <- c(26, -20, 26)

y_poly <- c(-20,26,26)

polygon(x_poly, y_poly,col = rgb(0, 0, 1, 0.15),border = "blue",lwd=2,lty=2)

text(7,11,label = expression(bold(X[1] + X[2] > 6)),cex=1.6,col="blue",font=2)

points(X.tail.sim2[rowSums(X.tail.sim2)>6,],pch=16,col="blue", cex = 0.7)

# Event set 2: (X_1^2 + X_2^2) > 5^2

plot(X.tail.sim2,
     pch=16,col = rgb(0, 1, 0, 1),cex=0.7,
     xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,
     xlim=c(-12,12),ylim=c(-12,12))

lines(quant.set,lwd=3,col="red")

theta <- seq(0, 2*pi, length.out = 100)

outer_x <- 20 * cos(theta)

outer_y <- 20 * sin(theta)

inner_x <- 5 * cos(theta)

inner_y <- 5 * sin(theta)

all_x <- c(outer_x, NA, inner_x)

all_y <- c(outer_y, NA, inner_y)

polypath(all_x, all_y,col = rgb(0, 0, 1, 0.15),border = "blue",lwd=2,lty=2)

text(-7,10,label = expression(bold(X[1]^2 + X[2]^2 > 5^2)),cex=1.6,col="blue",font=2)

points(X.tail.sim2[rowSums(X.tail.sim2^2)>25,],pch=16,col="blue", cex = 0.7)

dev.off()