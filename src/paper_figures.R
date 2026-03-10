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

X.norm <- sim_data_normal

# Specify sample size and dimension of data 

n <- dim(X.norm)[1]
d <- dim(X.norm)[2]

# Decompose data into radial and angular components 

R <- apply(X.norm,1,l2_norm)

W <- X.norm/R
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
  
  # Create a checkpoint for fitting the model 
  
  checkpoint <- callback_model_checkpoint(filepath = "../QR_est/qr_fit",
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

  # Run ADAM to estimate parameters 
  
  history <- quant.model %>%
    fit(list(W.train), R.train,
        epochs = n.epochs,
        batch_size = batch.size,
        callback = list(checkpoint,callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 5)),
        validation_data = list(list(input_pseudo_angles = W.valid),R.valid))

  # Load the best model weights from the saved checkpoint 
  
  quant.model <- load_model_weights_tf(quant.model, filepath = "../QR_est/qr_fit")

  # Save the best tf model
  
  save_model_tf(quant.model, "../QR_est/qr_fit")
  
}

# Load in best fitting tf model 

quant.model <- tf$saved_model$load("../QR_est/qr_fit")

# Compute estimated quantile over all observed angles 

u <- k_get_value(quant.model(k_constant(W)))
 
# Fit GPD over angles  ------------------------------------------

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
#HERERHEREHREHR

if(fit_GPD){
  
  # Create a checkpoint for fitting the model 
  
  checkpoint <- callback_model_checkpoint(filepath = paste0("runs/GPD_est/tmp_paper_fig/gpd_fit_{epoch:02d}"),
                                          monitor = "val_loss",
                                          verbose = 0,
                                          save_best_only = FALSE,
                                          save_weights_only = TRUE,
                                          mode = "min",
                                          save_freq = "epoch")

  n.epochs <- 750
  batch.size <- length(R.train)
  
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

  optimal_epoch <- which.min(history$metrics$val_loss)

  GPD.model <- load_model_weights_tf(GPD.model,
                                     filepath = paste0("runs/GPD_est/tmp_paper_fig/gpd_fit_",sprintf("%02d", optimal_epoch-1)))

  save_model_weights_tf(GPD.model,
                        paste0("runs/GPD_est/tmp_paper_fig/gpd_fit_best"))

  save_model_tf(GPD.model,
                paste0("runs/GPD_est/gpd_fit_","paper_fig"))

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

        optimal_epoch <- NA

      } else {

        optimal_epoch <- which.min(history$metrics$val_loss)

        best.valid.loss <- history$metrics$val_loss[optimal_epoch]

      }

    } else {

      optimal_epoch <- NA

    }

    if (!is.na(optimal_epoch)) {

      GPD.model <- load_model_weights_tf(GPD.model,
                                         filepath = paste0("runs/GPD_est/tmp_paper_fig/gpd_fit_",sprintf("%02d", optimal_epoch)))

      save_model_weights_tf(GPD.model,
                            paste0("runs/GPD_est/tmp_paper_fig/gpd_fit_best"))

      save_model_tf(GPD.model,
                    paste0("runs/GPD_est/gpd_fit_","paper_fig"))

    } else{

      load_model_weights_tf(GPD.model,
                            paste0("runs/GPD_est/tmp_paper_fig/gpd_fit_best"))

    }

  }


  unlink(paste0("runs/GPD_est/tmp_paper_fig"),recursive = T, force = T)
  
}



GPD.model <- load_model_tf(filepath = paste0("runs/GPD_est/gpd_fit_","paper_fig"),
                          custom_objects = list(GPD_loss=GPD_loss))

#Validation 

#Evaluating model fit on test set 
pred.GPD.test <- k_get_value(GPD.model(list(k_constant(W.test), k_constant(u.test))))

obs_quants <- qexp(apply(cbind(u.test, pred.GPD.test[,1:2], R.test+u.test), 1, function(x) pgpd(x[4], loc = x[1], scale = x[2]/(1 +x[3]), x[3])))

n_p <- length(R.test) 
ps <- (1:n_p)/(n_p + 1) 
theor_quants <- qexp(ps)

png(
  filename = paste0("Plots/gpd_qq_plot_","paper_fig",".png"),
  width = 2500, height = 2500, res = 400
)

par(mgp = c(2.8, 0.7, 0))

plot(
  theor_quants, sort(obs_quants),
  pch = 16, col = "grey",
  cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.8,
  xlim = range(sort(obs_quants), theor_quants),
  ylim = range(sort(obs_quants), theor_quants),
  xlab = "Theoretical", ylab = "Observed",
  main = "GPD QQ plot"
)
abline(a = 0, b = 1, lwd = 4, col = 2)
points(theor_quants, sort(obs_quants), pch = 16, col = "black", cex = 1.4)

dev.off()

#We only want to evaluate marginal tails outside of the quantile set. Therefore, we just looks at marginal probabilities in this region

pred_phis <- seq(0,2*pi,length.out=1001)

W.circle <- cbind(cos(pred_phis),sin(pred_phis))

u.circle <- k_get_value(quant.model(k_constant(W.circle)))

quant.set <- c(u.circle)*W.circle

quant.set.maxima <- apply(quant.set,2,max)

quant.set.minima <- apply(quant.set,2,min)

n.sim <- 1e5

X.norm.excess <- X.norm[which(pred.quant < R),]

# W.sim <- W[sample(1:nrow(W), n.sim,replace=T),]

W.excess <- W[which(pred.quant < R),]

W.sim <- W.excess[sample(1:nrow(W.excess), n.sim,replace=T),]
u.sim <- k_get_value(quant.model(k_constant(W.sim)))

dim(u.sim) <- c(length(u.sim), 1)

pred.GPD.sim <- k_get_value(GPD.model(list(k_constant(W.sim), k_constant(u.sim))))
nu.sim <- pred.GPD.sim[, 1]
xi.sim <- pred.GPD.sim[, 2]

unif_sample <- runif(n.sim,0,1)
R.sim <- apply(cbind(unif_sample,u.sim,nu.sim,xi.sim),1,function(x)
  qgpd(x[1], loc=x[2], scale = x[3]/(1+ x[4]), shape=x[4])
)

X.sim.excess <- R.sim*W.sim

X.excess.obs.os <- X.norm.excess

X.excess.sim.os <- X.sim.excess

m = 500

for(i in 1:d){

  m = 500

  title = paste0("Variable ",i)

  X.obs.tail = X.excess.obs.os[,i][X.norm.excess[,i] > quant.set.maxima[i]]

  X.sim.tail = as.numeric(X.excess.sim.os[,i][X.sim.excess[,i] > quant.set.maxima[i]])

  if(m > length(X.obs.tail)){m = length(X.obs.tail)}

  ps  = (1:m)/(m+1)

  obs_quants = quantile(X.obs.tail,ps)
  sim_quants = quantile(X.sim.tail,ps)

  png(
    filename = paste0("Plots/marginal_tails_",title,"_upper_","paper_fig",".png"),
    width = 2500, height = 2500, res = 400
  )

  par(mgp = c(2.8, 0.7, 0))

  plot(
    sim_quants,obs_quants,
    pch=16,col="grey",
    cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.8,
    xlim=range(obs_quants,sim_quants),
    ylim=range(obs_quants,sim_quants),
    ylab = "Observed",xlab = "Model",
    main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," upper tail") )
  abline(a=0,b=1,lwd=4,col=2)
  points(sim_quants,obs_quants,pch=16,col="black",cex=1.3)

  dev.off()

}


for(i in 1:d){

  m = 500

  title = paste0("Variable ",i)

  X.obs.tail = X.excess.obs.os[,i][X.norm.excess[,i] < quant.set.minima[i]]

  X.sim.tail = as.numeric(X.excess.sim.os[,i][X.sim.excess[,i] < quant.set.minima[i]])

  if(m > length(X.obs.tail)){m = length(X.obs.tail)}

  ps  = (1:m)/(m+1)

  obs_quants = quantile(X.obs.tail,ps)
  sim_quants = quantile(X.sim.tail,ps)

  png(
    filename = paste0("Plots/marginal_tails_",title,"_lower_","paper_fig",".png"),
    width = 2500, height = 2500, res = 400
  )

  par(mgp = c(2.8, 0.7, 0))

  plot(
    sim_quants,obs_quants,
    pch=16,col="grey",
    cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.8,
    xlim=range(obs_quants,sim_quants),
    ylim=range(obs_quants,sim_quants),
    ylab = "Observed",xlab = "Model",
    main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," lower tail"))
  abline(a=0,b=1,lwd=4,col=2)
  points(sim_quants,obs_quants,pch=16,col="black",cex=1.3)

  dev.off()

}

labels <- paste0("Variable ",1:d)

tail_dep <- function(x, y, u = 0.95) {
  qx <- quantile(x, u, na.rm = TRUE)
  qy <- quantile(y, u, na.rm = TRUE)
  mean(x > qx & y > qy, na.rm = TRUE) / (1 - u)
}

png(
  file = paste0("Plots/scatterplot_", "paper_fig", ".png"),
  width = 3500, height = 3500, res = 400
)

p = ncol(X.excess.obs.os)

# 2x2 layout, all panels square, equal size
layout(matrix(1:4, nrow = 2, ncol = 2), widths = rep(1,2), heights = rep(1,2))
par(mar = c(4, 4, 2, 1))  # inner margin (bottom, left, top, right)
par(oma = c(1, 1, 4, 1))  # outer margin for overall labels / title


for (j in 1:p) {
  for (i in 1:p) {

    if (i == j) {
      # ---------------- Diagonal: QQ plot ----------------
      m = min(c(length(X.excess.sim.os[, i]),length(X.excess.obs.os[,j])))

      ps = (1:m)/(m+1)

      obs_quants = quantile(X.excess.obs.os[,i],ps)
      sim_quants = quantile(X.excess.sim.os[, j],ps)

      plot(
        sim_quants,obs_quants,
        pch=16,col="grey",
        cex.lab = 1.3, cex.axis = 1, cex.main = 1.3,
        xlim=range(obs_quants,sim_quants),
        ylim=range(obs_quants,sim_quants),
        ylab = ifelse(j==1,"Observed","") ,xlab =ifelse(i==1,"Model",""),
        main = paste0(toupper(substring(labels[j], 1, 1)), substring(labels[j], 2)," quantiles"))
      abline(a=0,b=1,lwd=4,col=2)
      points(sim_quants,obs_quants,pch=16,col="black",cex=1.3)

    } else if (j > i) {
      # ---------------- Upper diagonal: Scatter ----------------
      xs <- X.excess.sim.os[, j]; ys <- X.excess.sim.os[, i]
      xo <- X.excess.obs.os[, j]; yo <- X.excess.obs.os[, i]
      xr <- range(xs, xo); yr <- range(ys, yo)

      xlab = ifelse(i==1,paste0(toupper(substring(labels[j], 1, 1)), substring(labels[j], 2)),"")

      ylab = ifelse(i==j-1,paste0(toupper(substring(labels[i], 1, 1)), substring(labels[i], 2)),"")

      # plot(xs[1:1000], ys[1:1000], pch = 16, col = rgb(0, 1, 0, 0.35), cex = 0.7,xlab=xlab,ylab=ylab,cex.lab=1.3)
      # points(xo[1:1000], yo[1:1000], pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7)
      plot(xs, ys, pch = 16, col = rgb(0, 1, 0, 0.35), cex = 0.7,xlab=xlab,ylab=ylab,cex.lab=1.3)
      points(xo, yo, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7)

    } else {
      # ---------------- Lower diagonal: dependence coeff ----------------
      xs <- X.excess.sim.os[, i]; ys <- X.excess.sim.os[, j]
      xo <- X.excess.obs.os[, i]; yo <- X.excess.obs.os[, j]

      r_sim <- cor(xs, ys, method = "spearman", use = "complete.obs")
      r_obs <- cor(xo, yo, method = "spearman", use = "complete.obs")

      lam_sim <- tail_dep(xs, ys)
      lam_obs <- tail_dep(xo, yo)

      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0, 0.4, bquote(rho[m]==.(round(r_sim, 2))~", "~chi[m]==.(round(lam_sim, 2))),
           col="green3", cex=1.2)
      text(0, -0.4, bquote(rho[o]==.(round(r_obs, 2))~", "~chi[o]==.(round(lam_obs, 2))),
           col="grey2", cex=1.2)
      abline(h = 0, col = rgb(0,0,0,0.1))

    }
  }
}

mtext("Model vs observed extremes", outer = TRUE, side = 3, line = 1.2, cex = 1.2)
dev.off()

n.sim = 1e6

W.sim = W.excess[sample(1:nrow(W.excess), n.sim,replace=T),]
u.sim <- k_get_value(quant.model(k_constant(W.sim)))

dim(u.sim) = c(length(u.sim), 1)

pred.GPD.sim <- k_get_value(GPD.model(list(k_constant(W.sim), k_constant(u.sim))))
nu.sim = pred.GPD.sim[, 1]
xi.sim = pred.GPD.sim[, 2]

unif_sample = runif(n.sim,0,1)
R.sim = apply(cbind(unif_sample,u.sim,nu.sim,xi.sim),1,function(x)
  qgpd(x[1], loc=x[2], scale = x[3]/(1+ x[4]), shape=x[4])
)

X.sim.excess = R.sim*W.sim

n_y = 1300 #(varies slightly depending on year)

png(
  filename = paste0("Plots/rl_plot_upper_","paper_fig",".png"),
  width = 4000, height = 2000, res = 400
)

par(mgp = c(2.8, 0.7, 0))

par(mfrow=c(1,2))

for(i in 1:d){
  
  title = paste0("Variable ",i)
  
  emp_RLs = sort(X.norm.excess[,i][X.norm.excess[,i] > quant.set.maxima[i]])
  
  emp_probs = 1 - sort(rank(X.norm.excess[,i])[X.norm.excess[,i] > quant.set.maxima[i]]/(nrow(X.norm.excess)+1))
  
  emp_RPs = 1/(n_y*emp_probs*(1-quant.level))
  
  fix_RPs = exp(seq(log(emp_RPs[1]),log(1000),length.out = 3000))
  
  fix_probs = 1/(fix_RPs*n_y*(1-quant.level))
  
  model_RLs = quantile(X.sim.excess[,i],1-fix_probs)
  
  plot(
    fix_RPs, model_RLs,
    log = "x",
    type = "l", lwd=3, col = "grey",
    xlab = "Return period (log scale)",
    ylab = "",
    xaxt = "n",
    cex.lab=1.3,cex.main=1.5,
    ylim=range(emp_RLs,model_RLs),
    main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," upper tail")
  )
  
  points(emp_RPs,emp_RLs,pch=16,col="black",cex=1.4)
  
  ticks <- c(0.1, 1, 10, 100, 1000)
  
  axis(1, at = ticks, labels = as.character(ticks))
  
  if(i == 1){
    legend("topleft",legend=c("Model","Empirical"),
           col = c("grey", "black"),
           lty = c(1, NA),          # line type for line; NA for no line
           pch = c(NA, 16),         # point type for point; NA for no point
           lwd = c(3, NA),
           pt.cex = 1.2,
           cex=1.2,
           bg="white")
  }
  
}

dev.off()

png(
  filename = paste0("Plots/rl_plot_lower_","paper_fig",".png"),
  width = 4000, height = 2000, res = 400
)

par(mgp = c(2.8, 0.7, 0))

par(mfrow=c(1,2))

for(i in 1:d){
  
  title = paste0("Variable ",i)
  
  emp_RLs = sort(X.norm.excess[,i][X.norm.excess[,i] < quant.set.minima[i]])
  
  emp_probs = sort(rank(X.norm.excess[,i])[X.norm.excess[,i] < quant.set.minima[i]]/(nrow(X.norm.excess)+1))
  
  emp_RPs = 1/(n_y*emp_probs*(1-quant.level))
  
  fix_RPs = exp(seq(log(emp_RPs[length(emp_RPs)]),log(1000),length.out = 3000))
  
  fix_probs = 1/(fix_RPs*n_y*(1-quant.level))
  
  model_RLs = quantile(X.sim.excess[,i],fix_probs)
  
  plot(
    fix_RPs, model_RLs,
    log = "x",
    type = "l", lwd=3, col = "grey",
    xlab = "Return period (log scale)",
    ylab = "",
    xaxt = "n",
    cex.lab=1.3,cex.main=1.5,
    ylim=range(emp_RLs,model_RLs),
    main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," lower tail")
  )
  
  points(emp_RPs,emp_RLs,pch=16,col="black",cex=1.4)
  
  ticks <- c(0.1, 1, 10, 100, 1000)
  
  axis(1, at = ticks, labels = as.character(ticks))
  
  if(i == 1){
    legend("topright",legend=c("Model","Empirical"),
           col = c("grey", "black"),
           lty = c(1, NA),          # line type for line; NA for no line
           pch = c(NA, 16),         # point type for point; NA for no point
           lwd = c(3, NA),
           pt.cex = 1.2,
           cex=1.2,
           bg="white")
  }
  
}

dev.off()


png(
  filename = paste0("Plots/SPAR_illustration",".png"),
  width = 4500, height = 1500, res = 400
)

par(mgp = c(2.8, 0.7, 0))

par(mfrow=c(1,3))

plot(X.norm, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))

plot(X.norm[which(pred.quant >= R),], pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))
points(X.norm.excess,pch=16,col = rgb(0, 1, 0, 1),cex=0.7)
lines(quant.set,lwd=3,col="red")

plot(X.norm.excess, pch = 16, col = rgb(0, 1, 0, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))
lines(quant.set,lwd=3,col="red")
rect(quant.set.maxima,-15,15,15,col = rgb(0.5, 0.5, 0.5, 0.15),border = "black",lwd=2,lty=2)
text(7,-2,label = expression(bold(X[1] > q[1]^M)),cex=1.6,col="black",font=2)

rect(-15, 4, -2, 15, col = rgb(0, 0, 1, 0.15), border = "blue", lwd = 2, lty = 2)
text(-6.5,7,label = expression(bold(X[1] < -2) * "," * bold(X[2] > 4)),cex=1.6,col="blue")

# rect(-15, -15, 15, -4, col = rgb(1, 0.647, 0, 0.15), border = "orange", lwd = 2, lty = 2)
# text(-7,-8,label = expression(bold(X[2] < -4)),cex=1.6,col="orange")

x_poly <- c(-26, 20, -26)
y_poly <- c(20,-26,-26)

polygon(x_poly, y_poly, col = rgb(1, 0.647, 0, 0.15), border = "orange",lwd = 2, lty = 2)
text(-7,-8,label = expression(bold(X[1] + X[2] < -6)),cex=1.6,col="orange")

dev.off()

png(
  filename = paste0("Plots/SPAR_illustration_1",".png"),
  width = 2000, height = 2000, res = 400
)

par(mar = c(4, 4, 3, 1), mgp = c(2.5, 0.7, 0))

plot(X.norm, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10),asp=1)

mtext("(a)", side = 3, line = 1, adj = 0, cex = 1.5)

dev.off()


png(
  filename = paste0("Plots/SPAR_illustration_2",".png"),
  width = 2000, height = 2000, res = 400
)

par(mar = c(4, 4, 3, 1), mgp = c(2.5, 0.7, 0))

plot(X.norm[which(pred.quant >= R),], pch = 16,asp=1, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))
points(X.norm.excess,pch=16,col = rgb(0, 1, 0, 1),cex=0.7)
lines(quant.set,lwd=3,col="red")

mtext("(b)", side = 3, line = 1, adj = 0, cex = 1.5)

dev.off()


png(
  filename = paste0("Plots/SPAR_illustration_3",".png"),
  width = 2000, height = 2000, res = 400
)

par(mar = c(4, 4, 3, 1), mgp = c(2.5, 0.7, 0))

plot(X.norm.excess, pch = 16,asp=1, col = rgb(0, 1, 0, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))
lines(quant.set,lwd=3,col="red")
rect(quant.set.maxima,-15,15,15,col = rgb(0.5, 0.5, 0.5, 0.15),border = "black",lwd=2,lty=2)
text(7,-2,label = expression(bold(X[1] > q[1]^M)),cex=1.4,col="black",font=2)

mtext("(c)", side = 3, line = 1, adj = 0, cex = 1.5)

rect(-15, 4, -2, 15, col = rgb(0, 0, 1, 0.15), border = "blue", lwd = 2, lty = 2)
text(-6.5,7,label = expression(bold(X[1] < -2) * "," * bold(X[2] > 4)),cex=1.4,col="blue")

# rect(-15, -15, 15, -4, col = rgb(1, 0.647, 0, 0.15), border = "orange", lwd = 2, lty = 2)
# text(-7,-8,label = expression(bold(X[2] < -4)),cex=1.4,col="orange")

x_poly <- c(-26, 20, -26)
y_poly <- c(20,-26,-26)

polygon(x_poly, y_poly, col = rgb(1, 0.647, 0, 0.15), border = "orange",lwd = 2, lty = 2)
text(-7,-8,label = expression(bold(X[1] + X[2] < -6)),cex=1.4,col="orange")

dev.off()

n.sim = 1e4

X.norm.excess = X.norm[which(pred.quant < R),]

# W.sim = W[sample(1:nrow(W), n.sim,replace=T),]

W.excess = W[which(pred.quant < R),]

W.sim = W.excess[sample(1:nrow(W.excess), n.sim,replace=T),]
u.sim <- k_get_value(quant.model(k_constant(W.sim)))

dim(u.sim) = c(length(u.sim), 1)

pred.GPD.sim <- k_get_value(GPD.model(list(k_constant(W.sim), k_constant(u.sim))))
nu.sim = pred.GPD.sim[, 1]
xi.sim = pred.GPD.sim[, 2]

unif_sample = runif(n.sim,0,1)
R.sim = apply(cbind(unif_sample,u.sim,nu.sim,xi.sim),1,function(x)
  qgpd(x[1], loc=x[2], scale = x[3]/(1+ x[4]), shape=x[4])
)

X.sim.excess = R.sim*W.sim

png(
  filename = paste0("Plots/SPAR_event_sets",".png"),
  width = 4500, height = 1500, res = 400
)

par(mgp = c(2.8, 0.7, 0))

par(mfrow=c(1,3))

plot(X.sim.excess,pch=16,col = rgb(0, 1, 0, 1),cex=0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-12,12),ylim=c(-12,12))
lines(quant.set,lwd=3,col="red")
rect(5,-15,15,15,col = rgb(0, 0, 1, 0.15),border = "blue",lwd=2,lty=2)
text(8,-3,label = expression(bold(X[1] > 5)),cex=1.6,col="blue",font=2)
points(X.sim.excess[X.sim.excess[,1]>5,],pch=16,col="blue", cex = 0.7)

plot(X.sim.excess,pch=16,col = rgb(0, 1, 0, 1),cex=0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-12,12),ylim=c(-12,12))
lines(quant.set,lwd=3,col="red")
x_poly <- c(26, -20, 26)
y_poly <- c(-20,26,26)
polygon(x_poly, y_poly,col = rgb(0, 0, 1, 0.15),border = "blue",lwd=2,lty=2)
text(7,11,label = expression(bold(X[1] + X[2] > 6)),cex=1.6,col="blue",font=2)
points(X.sim.excess[rowSums(X.sim.excess)>6,],pch=16,col="blue", cex = 0.7)

plot(X.sim.excess,pch=16,col = rgb(0, 1, 0, 1),cex=0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-12,12),ylim=c(-12,12))
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
points(X.sim.excess[rowSums(X.sim.excess^2)>25,],pch=16,col="blue", cex = 0.7)

dev.off()

# png(
#   filename = paste0("Plots/SPAR_illustration_old",".png"),
#   width = 4500, height = 1500, res = 400
# )
# 
# par(mgp = c(2.8, 0.7, 0))
# 
# par(mfrow=c(1,3))
# 
# plot(X.norm, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))
# 
# plot(X.norm[which(pred.quant >= R),], pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))
# points(X.norm.excess,pch=16,col = rgb(0, 1, 0, 1),cex=0.7)
# lines(quant.set,lwd=3,col="red")
# 
# plot(X.norm.excess, pch = 16, col = rgb(0, 1, 0, 1), cex = 0.7,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.3,xlim=c(-10,10),ylim=c(-10,10))
# lines(quant.set,lwd=3,col="red")
# rect(quant.set.maxima,-15,15,15,col = rgb(0.5, 0.5, 0.5, 0.15),border = "black",lwd=2,lty=2)
# text(7,-5,label = expression(bold(X[1] > q[1]^u)),cex=1.8,col="black",font=2)
# # legend(5,-5,legend = expression(X[1] > q[1]^u),cex=1.8,bg="white",bty="n",text.col = "black")
# 
# rect(-15, -15, quant.set.minima[1], 15, col = rgb(0, 0, 1, 0.15), border = "blue", lwd = 2, lty = 2)
# text(-7,5,label = expression(bold(X[1] < q[1]^l)),cex=1.8,col="blue")
# 
# dev.off()
