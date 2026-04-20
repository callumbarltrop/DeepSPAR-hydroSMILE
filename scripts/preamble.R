# Functions and packages for fitting deepSPAR framework ------------------------------------------------------

# Loss function for quantile regression
tilted_loss <- function(y_true, y_pred) {
  K <- backend()

  error = y_true - y_pred
  return(K$mean(K$maximum(quant.level * error, (quant.level - 1) * error)))
}

# Loss function for GPD
GPD_loss <- function( y_true, y_pred) {
    
    K <- backend()
    
    nu=y_pred[all_dims(),1] #See Pasche and Engelke, 2023
    xi=y_pred[all_dims(),2]
    u=y_pred[all_dims(),3]
    y=y_true[all_dims(),1]
 
    #Evaluate log-likelihood
    ll1=(1/xi+1)*tf$math$log1p(xi*(xi+1)*y/nu)

    ll2= K$log(nu) -K$log(1+xi)
    
 
    return((K$sum(ll1+ll2)) )
 

}

# Algorithm for computing the geometric median 
geometric_median <- function(X, eps = 1e-6, maxit = 1000) {
  # X: n x d matrix of points
  n <- nrow(X)
  d <- ncol(X)
  y <- colMeans(X)   # initial guess
  
  for (i in 1:maxit) {
    dist <- sqrt(rowSums((X - matrix(y, n, d, byrow = TRUE))^2))
    if (any(dist == 0)) return(X[which.min(dist), ])  # exact match
    
    w <- 1 / dist
    y_new <- colSums(X * w) / sum(w)
    
    if (sqrt(sum((y_new - y)^2)) < eps) break
    y <- y_new
  }
  return(y)
}

# Function to estimate chi empirically 
tail_dep <- function(u,x, y) {
  qx <- quantile(x, u, na.rm = TRUE)
  qy <- quantile(y, u, na.rm = TRUE)
  mean(x > qx & y > qy, na.rm = TRUE) / (1 - u)
}

#Checking for required packages. This function will install any required packages if they are not already installed
packages = c("evd") 
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Algorithm for simulation uniformly on the hypersphere 
uniform_sample_dsphere = function(d,sampling_points,simulation_size = 1000){ #Using the technique of Marsaglia (1972) - see https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-43/issue-2/Choosing-a-Point-from-the-Surface-of-a-Sphere/10.1214/aoms/1177692644.full
  sphere_sample = matrix(nrow=0,ncol=d)
  while(dim(sphere_sample)[1]<sampling_points){
    uniform_sample = matrix(runif(d*simulation_size,-1,1),ncol=d)
    r = apply(uniform_sample,1,norm,type="2")
    sphere_sample = rbind(sphere_sample, (uniform_sample[r<=1,]/r[r<=1])) #Reject samples not inside sphere, then rescale to get points on sphere
  }
  return(sphere_sample[1:sampling_points,])
}

# Functions for empirically estimating probabilities associated with certain tail regions 
l1_prob_upper = function(u,sim_data,tau=0){
  return((mean(rowSums(sim_data) > u)*(1-tau)))
}

l1_prob_lower = function(u,sim_data,tau=0){
  return((mean(rowSums(sim_data) < u)*(1-tau)))
}

survivor_prob = function(u_vec,sim_data,tau=0){
  return((mean(apply(sim_data > matrix(u_vec, nrow = nrow(sim_data), ncol = ncol(sim_data), byrow = TRUE), 1, all))*(1-tau)))
}

cdf_prob = function(u_vec,sim_data,tau=0){
  return((mean(apply(sim_data < matrix(u_vec, nrow = nrow(sim_data), ncol = ncol(sim_data), byrow = TRUE), 1, all))*(1-tau)))
}

# L2 norm function 
l2_norm = function(x){
  return(sqrt(sum(x^2)))
}

