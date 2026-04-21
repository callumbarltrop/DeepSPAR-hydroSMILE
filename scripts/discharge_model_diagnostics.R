# Setting up grid for testing tuning parameters  -----------------------------------

rm(list=ls())

# NN and GPD tuning parameters --------------------------------------------

boot_vals_dir = "C:/Users/CMB/Documents/BootVals3/BootVals/"

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

# Packages and VE ---------------------------------------------------------

# Load some useful functions and packages
source("preamble.R")

# Data --------------------------------------------------------------------

# Locations of stations 
river_names = c("ammer","iller","lech","loisach")

#Split into four 30 year windows
time_periods = list(1980:2009,2010:2039,2040:2069,2070:2099)

std_sigmas = readRDS(file = "RunOffData/std_sigmas.rds")

std_mus = readRDS(file = "RunOffData/std_mus.rds")

std_tw_dfs = readRDS(file = "RunOffData/std_tw_dfs.rds")

os_tw_dfs = readRDS(file = "RunOffData/os_tw_dfs.rds")

l1_seq_upper = 1000

l1_seq_lower = 30

prob_RPs = c(10)

n_y = 52 #(varies slightly depending on year)

# Adding uncertainty bounds ------------------------------------------------

uncert_figs = T

indices_grid = expand.grid(c(1:length(time_periods)))

dir.create("PNGs/opt_fits_4d", showWarnings = FALSE)

#Number of bootstraps 
n_boot = 100

for(array_id in 1:nrow(indices_grid)){
  
  print(paste("Processing task", array_id))
  
  # Index of time window. We use middle of the period for optimisation
  time_window_ind = indices_grid$Var1[array_id]
  
  if(uncert_figs == T){
    
    X.obs = os_tw_dfs[[time_window_ind]]
    # 
    # polar = rect2polar(t(X.norm))  #Get polar coordinates
    # n = dim(X.norm)[1]
    # d = dim(X.norm)[2]
    # 
    # # Create AR decomp 
    # R = polar$r
    # W = X.norm/R
    # W = as.matrix(W)
    
    #Validation 

    # obs_quants = readRDS(file = paste0(boot_vals_dir,"gpd_qq_quants_","_tw_",time_window_ind,".rds"))
    # 
    # n_p = length(obs_quants)
    # ps = (1:n_p)/(n_p + 1)
    # theor_quants = qexp(ps)
    # 
    # boot_quant_mat = sapply(1:n_boot,function(boot_ind){
    #   if(file.exists(paste0(boot_vals_dir,"gpd_qq_quants_","bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
    #     return(readRDS(file = paste0(boot_vals_dir,"gpd_qq_quants_","bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
    #   } else {
    #     print(boot_ind);print(time_window_ind)
    #     return(rep(NA,length(ps)))
    #   }
    # })
    # 
    # png(
    #   filename = paste0("PNGs/opt_fits_4d/gpd_qq_plot_","tw_",time_window_ind,"_cis.png"),
    #   width = 2500, height = 2500, res = 400
    # )
    # 
    # par(mgp = c(2.8, 0.7, 0))
    # 
    # plot(
    #   theor_quants, sort(obs_quants),
    #   pch = 16, col = "grey",
    #   cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.8,
    #   xlim = range(sort(obs_quants), theor_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
    #   ylim = range(sort(obs_quants), theor_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
    #   ylab = "SPAR", xlab = "Theoretical",
    #   main = paste0("GPD QQ plot, ",min(time_periods[[array_id]]),"-",max(time_periods[[array_id]])),
    #   sub=""
    # )
    # polygon(c(rev(theor_quants), theor_quants), c(rev(apply(boot_quant_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_quant_mat,1,quantile,probs=0.975,na.rm=T)), col = 'grey80', border = NA)
    # abline(a = 0, b = 1, lwd = 4, col = 2)
    # points(theor_quants, sort(obs_quants), pch = 16, col = "black", cex = 1.4)
    # 
    # dev.off()
    # 
    # labels <- river_names
    # 
    # tail_dep <- function(u = 0.95,x, y) {
    #   qx <- quantile(x, u, na.rm = TRUE)
    #   qy <- quantile(y, u, na.rm = TRUE)
    #   mean(x > qx & y > qy, na.rm = TRUE) / (1 - u)
    # }
    # 
    # u_vec = seq(0.8,0.99,by = 0.01)
    # 
    # png(
    #   file = paste0("PNGs/opt_fits_4d/scatterplot_", "tw_", time_window_ind, "_cis.png"),
    #   width = 3500, height = 3500, res = 400
    # )
    # 
    # p <- 4
    # 
    # # 4x4 layout, all panels square, equal size
    # layout(matrix(1:16, nrow = 4, ncol = 4), widths = rep(1,4), heights = rep(1,4))
    # par(mar = c(4, 4.5, 2, 1))  # inner margin (bottom, left, top, right)
    # par(oma = c(1, 1, 4, 1))  # outer margin for overall labels / title
    # 
    # 
    # for (j in 1:p) {
    #   for (i in 1:p) {
    # 
    #     if (i == j) {
    # 
    #       obs_quants = readRDS(file = paste0(boot_vals_dir,"marg_qq_quants_obs_",labels[j],"_tw_",time_window_ind,".rds"))
    #       sim_quants = readRDS(file = paste0(boot_vals_dir,"marg_qq_quants_",labels[j],"_tw_",time_window_ind,".rds"))
    # 
    #       boot_quant_mat = sapply(1:n_boot,function(boot_ind){
    #         if(file.exists(paste0(boot_vals_dir,"marg_qq_quants_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
    #           return(readRDS(file = paste0(boot_vals_dir,"marg_qq_quants_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
    #         } else {
    #           print(boot_ind);print(time_window_ind)
    #           return(rep(NA,length(ps)))
    #         }
    #       })
    # 
    #       plot(
    #         sim_quants,obs_quants,
    #         pch=16,col="grey",
    #         cex.lab = 1.3, cex.axis = 1, cex.main = 1.3,
    #         xlim=range(obs_quants,sim_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
    #         ylim=range(obs_quants,sim_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
    #         ylab = ifelse(j==1,expression("WaSiM discharge (m"^3*"/s)"),"") ,xlab =ifelse(i==1,expression("SPAR discharge (m"^3*"/s)"),""),
    #         main = paste0(toupper(substring(labels[j], 1, 1)), substring(labels[j], 2)))
    #       polygon( c(rev(apply(boot_quant_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_quant_mat,1,quantile,probs=0.975,na.rm=T)), c(rev(obs_quants), obs_quants), col = 'grey80', border = NA)
    #       abline(a=0,b=1,lwd=4,col=2)
    #       points(sim_quants,obs_quants,pch=16,col="black",cex=1.3)
    # 
    #     } else if (j > i) {
    #       # ---------------- Upper diagonal: Scatter ----------------
    # 
    #       x_sim = readRDS(file = paste0(boot_vals_dir,"scatterplot_",labels[i],"_",labels[j],"_tw_",time_window_ind,".rds"))
    # 
    #       x_obs = readRDS(file = paste0(boot_vals_dir,"scatterplot_observed_",labels[i],"_",labels[j],"_tw_",time_window_ind,".rds"))
    # 
    #       label_text <- paste0(toupper(substring(labels[j], 1, 1)), substring(labels[j], 2))
    #       if(i == 1){xlab=bquote(.(label_text) ~ "discharge" ~ (m^3/s))}else{xlab=""}
    # 
    #       label_text <- paste0(toupper(substring(labels[i], 1, 1)), substring(labels[i], 2))
    #       if(i == j-1){ylab=bquote(.(label_text) ~ "discharge" ~ (m^3/s))}else{ylab=""}
    # 
    #       plot(x_sim, pch = 16, col = rgb(0, 1, 0, 0.35), cex = 0.7,xlab=xlab,ylab=ylab,cex.lab=1.3)
    #       points(x_obs, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7)
    # 
    #     } else {
    #       # ---------------- Lower diagonal: chi dependence coefficient on scale ----------------
    #       x_obs = readRDS(file = paste0(boot_vals_dir,"scatterplot_observed_",labels[j],"_",labels[i],"_tw_",time_window_ind,".rds"))
    #       xo <- x_obs[, 1]; yo <- x_obs[, 2]
    # 
    #       lam_sim <- readRDS(file = paste0(boot_vals_dir,"chi_plot_",labels[i],"_",labels[j],"_tw_",time_window_ind,".rds"))
    #       lam_obs <- readRDS(file = paste0(boot_vals_dir,"chi_plot_obs_",labels[i],"_",labels[j],"_tw_",time_window_ind,".rds"))
    # 
    #       lam_obs_mat <- sapply(1:n_boot,function(boot_ind){
    #         resample_inds = sample(1:length(xo),length(xo),replace=T)
    # 
    #         return(sapply(u_vec,tail_dep,x=xo[resample_inds],y=yo[resample_inds]))
    #       })
    # 
    #       lam_mat = sapply(1:n_boot,function(boot_ind){
    #         if(file.exists(paste0(boot_vals_dir,"chi_plot_",labels[i],"_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
    #           return(readRDS(file = paste0(boot_vals_dir,"chi_plot_",labels[i],"_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
    #         } else {
    #           print(boot_ind);print(time_window_ind)
    #           return(rep(NA,length(lam_sim)))
    #         }
    #       })
    # 
    #       plot(u_vec,lam_obs,type="l",lwd=2,col="grey2",cex.lab=1.3, xlab=expression(u),ylab=expression(chi(u)),ylim=c(0,1))
    #       polygon(c(rev(u_vec), u_vec), c(rev(apply(lam_mat,1,quantile,probs=0.025,na.rm=T)), apply(lam_mat,1,quantile,probs=0.975,na.rm=T)), col = rgb(0,205/255,0,alpha=0.15), border = NA)
    #       polygon(c(rev(u_vec), u_vec), c(rev(apply(lam_obs_mat,1,quantile,probs=0.025,na.rm=T)), apply(lam_obs_mat,1,quantile,probs=0.975,na.rm=T)), col = rgb(0.5,0.5,0.5,alpha=0.1), border = NA)
    #       lines(u_vec,lam_sim,lwd=2,col="green3")
    #       lines(u_vec,lam_obs,type="l",lwd=2,col="grey")
    # 
    #       abline(h = 1,lwd=1,lty=2)
    # 
    #       abline(h = 0,lwd=1,lty=2)
    # 
    #     }
    #   }
    # }
    # 
    # mtext(paste0(min(time_periods[[array_id]]),"-",max(time_periods[[array_id]])), outer = TRUE, side = 3, line = 1.2, cex = 1.2)
    # dev.off()

    d = 4
    
    png(
      filename = paste0("PNGs/opt_fits_4d/rl_plot_upper_","tw_",time_window_ind,"_cis.png"),
      width = 2800, height = 2800, res = 400
    )

    par(mar = c(4, 4.5, 4, 2))
    
    par(mgp = c(2.8, 0.7, 0))

    par(mfrow=c(2,2))

    for(i in 1:d){

      title = river_names[i]

      fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))
      
      emp_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_obs_",title,"_tw_",time_window_ind,".rds"))
      
      emp_RPs = readRDS(file=paste0(boot_vals_dir,"rp_ests_4d_upper_obs_",title,"_tw_",time_window_ind,".rds"))
        
      model_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_tw_",time_window_ind,".rds"))

      boot_rl_mat = sapply(1:n_boot,function(boot_ind){
        if(file.exists(paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
          return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
        } else {
          print(boot_ind);print(time_window_ind)
          return(rep(NA,length(fix_RPs)))
        }
      })

      plot(
        fix_RPs, model_RLs,
        log = "x",
        type = "l", lwd=3, col = "grey",
        xlab = "Return period (years)",
        ylab = expression("Discharge (m"^3*"/s)"),
        xaxt = "n",
        cex.lab=1.3,cex.main=1.5,
        ylim=range(emp_RLs,model_RLs,apply(boot_rl_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
        main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," upper tail")
      )

      polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat,1,quantile,probs=0.975,na.rm=T)), col = 'grey80', border = NA)

      lines(fix_RPs, model_RLs, lwd=3, col = "grey3")

      points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)

      ticks <- c(10, 100, 1000,10000)

      axis(1, at = ticks, labels = as.character(ticks))

    }

    dev.off()

    d = 4
    
    png(
      filename = paste0("PNGs/opt_fits_4d/rl_plot_lower_","tw_",time_window_ind,"_cis.png"),
      width = 2800, height = 2800, res = 400
    )

    par(mar = c(4, 4.5, 4, 2))
    
    par(mgp = c(2.8, 0.7, 0))

    par(mfrow=c(2,2))

    for(i in 1:d){

      title = river_names[i]

      fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))
      
      emp_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_obs_",title,"_tw_",time_window_ind,".rds"))

      emp_RPs = readRDS(file=paste0(boot_vals_dir,"rp_ests_4d_lower_obs_",title,"_tw_",time_window_ind,".rds"))
      
      model_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_tw_",time_window_ind,".rds"))

      boot_rl_mat = sapply(1:n_boot,function(boot_ind){
        if(file.exists(paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
          return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
        } else {
          print(boot_ind);print(time_window_ind)
          return(rep(NA,length(fix_RPs)))
        }
      })

      plot(
        fix_RPs, model_RLs,
        log = "x",
        type = "l", lwd=3, col = "grey",
        xlab = "Return period (years)",
        ylab = expression("Discharge (m"^3*"/s)"),
        xaxt = "n",
        cex.lab=1.3,cex.main=1.5,
        ylim=range(emp_RLs,model_RLs,apply(boot_rl_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
        main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," lower tail")
      )

      polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat,1,quantile,probs=0.975,na.rm=T)), col = 'grey80', border = NA)

      lines(fix_RPs, model_RLs, lwd=3, col = "grey3")

      points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)

      ticks <- c(10, 100, 1000, 10000)

      axis(1, at = ticks, labels = as.character(ticks))

      if(i == 1){
        legend("topright",legend=c("SPAR","Empirical"),
               col = c("grey3", "darkblue"),
               lty = c(1, NA),          # line type for line; NA for no line
               pch = c(NA, 16),         # point type for point; NA for no point
               lwd = c(3, NA),
               pt.cex = 1.2,
               cex=1.2,
               bg="white")
      }

    }

    dev.off()
    
  }
  
}




# Adding uncertainty bounds with different sample sizes ------------------------------------------------

n_models = c(1,2,3,5,10)

percents = 100*(n_models/50)

indices_grid = expand.grid(c(1:length(time_periods)),c(1:length(n_models)))

dir.create("PNGs/opt_fits_4d", showWarnings = FALSE)

#Number of bootstraps 
n_boot = 100

for(array_id in 1:nrow(indices_grid)){
  
  print(paste("Processing task", array_id))
  
  # Index of time window. We use middle of the period for optimisation
  time_window_ind = indices_grid$Var1[array_id]
  
  # Number of models we use to fit our NN
  n_model = n_models[indices_grid$Var2[array_id]]
  
  obs_quants = readRDS(file = paste0(boot_vals_dir,"gpd_qq_quants_tw_",time_window_ind,"_nmodel_",n_model,".rds"))
  
  n_p = length(obs_quants)
  ps = (1:n_p)/(n_p + 1)
  theor_quants = qexp(ps)
  
  boot_quant_mat = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"gpd_qq_quants_","bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"gpd_qq_quants_","bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_model,".rds")))
    } else {
      print(boot_ind);print(time_window_ind);print(n_model)
      return(rep(NA,length(ps)))
    }
  })
  
  png(
    filename = paste0("PNGs/opt_fits_4d/gpd_qq_plot_","tw_",time_window_ind,"_nmodel_",n_model,"_cis.png"),
    width = 2500, height = 2500, res = 400
  )
  
  par(mgp = c(2.8, 0.7, 0))
  
  plot(
    theor_quants, sort(obs_quants),
    pch = 16, col = "grey",
    cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.8,
    xlim = range(sort(obs_quants), theor_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
    ylim = range(sort(obs_quants), theor_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
    ylab = "SPAR", xlab = "Theoretical",
    main = paste0("GPD QQ plot, ",min(time_periods[[time_window_ind]]),"-",max(time_periods[[time_window_ind]])),
    sub=paste0(percents[indices_grid$Var2[array_id]],"% subsample")
  )
  polygon(c(rev(theor_quants), theor_quants), c(rev(apply(boot_quant_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_quant_mat,1,quantile,probs=0.975,na.rm=T)), col = 'grey80', border = NA)
  abline(a = 0, b = 1, lwd = 4, col = 2)
  points(theor_quants, sort(obs_quants), pch = 16, col = "black", cex = 1.4)
  
  dev.off()
  
  labels <- river_names

  tail_dep <- function(u = 0.95,x, y) {
    qx <- quantile(x, u, na.rm = TRUE)
    qy <- quantile(y, u, na.rm = TRUE)
    mean(x > qx & y > qy, na.rm = TRUE) / (1 - u)
  }

  u_vec = seq(0.8,0.99,by = 0.01)

  png(
    file = paste0("PNGs/opt_fits_4d/scatterplot_", "tw_", time_window_ind,"_nmodel_",n_model, "_cis.png"),
    width = 3500, height = 3500, res = 400
  )

  p <- 4

  # 4x4 layout, all panels square, equal size
  layout(matrix(1:16, nrow = 4, ncol = 4), widths = rep(1,4), heights = rep(1,4))
  par(mar = c(4, 4.5, 2, 1))  # inner margin (bottom, left, top, right)
  par(oma = c(1, 1, 4, 1))  # outer margin for overall labels / title


  for (j in 1:p) {
    for (i in 1:p) {

      if (i == j) {
        # ---------------- Diagonal: QQ plot ----------------

        obs_quants = readRDS(file = paste0(boot_vals_dir,"marg_qq_quants_obs_",labels[j],"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))
        sim_quants = readRDS(file = paste0(boot_vals_dir,"marg_qq_quants_",labels[j],"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))

        boot_quant_mat = sapply(1:n_boot,function(boot_ind){
          if(file.exists(paste0(boot_vals_dir,"marg_qq_quants_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))){
            return(readRDS(file = paste0(boot_vals_dir,"marg_qq_quants_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_model,".rds")))
          } else {
            print(boot_ind);print(time_window_ind);print(n_model)
            return(rep(NA,length(ps)))
          }
        })

        plot(
          sim_quants,obs_quants,
          pch=16,col="grey",
          cex.lab = 1.3, cex.axis = 1, cex.main = 1.3,
          xlim=range(obs_quants,sim_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
          ylim=range(obs_quants,sim_quants,apply(boot_quant_mat,1,quantile,probs=c(0.025,0.975),na.rm=T)),
          ylab = ifelse(j==1,expression("WaSiM discharge (m"^3*"/s)"),"") ,xlab =ifelse(i==1,expression("SPAR discharge (m"^3*"/s)"),""),
          main = paste0(toupper(substring(labels[j], 1, 1)), substring(labels[j], 2)))
        polygon( c(rev(apply(boot_quant_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_quant_mat,1,quantile,probs=0.975,na.rm=T)), c(rev(obs_quants), obs_quants), col = 'grey80', border = NA)
        abline(a=0,b=1,lwd=4,col=2)
        points(sim_quants,obs_quants,pch=16,col="black",cex=1.3)

      } else if (j > i) {
        # ---------------- Upper diagonal: Scatter ----------------
        x_sim = readRDS(file = paste0(boot_vals_dir,"scatterplot_",labels[i],"_",labels[j],"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))
        x_obs = readRDS(file = paste0(boot_vals_dir,"scatterplot_observed_",labels[i],"_",labels[j],"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))
        
        label_text <- paste0(toupper(substring(labels[j], 1, 1)), substring(labels[j], 2))
        if(i == 1){xlab=bquote(.(label_text) ~ "discharge" ~ (m^3/s))}else{xlab=""}
        
        label_text <- paste0(toupper(substring(labels[i], 1, 1)), substring(labels[i], 2))
        if(i == j-1){ylab=bquote(.(label_text) ~ "discharge" ~ (m^3/s))}else{ylab=""}

        plot(x_sim[1:max(nrow(x_obs),1000),], pch = 16, col = rgb(0, 1, 0, 0.35), cex = 0.7,xlab=xlab,ylab=ylab,cex.lab=1.3)
        points(x_obs, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7)

      } else {
        # ---------------- Lower diagonal: chi dependence coefficient on scale ----------------
        x_obs = readRDS(file = paste0(boot_vals_dir,"scatterplot_observed_",labels[j],"_",labels[i],"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))
        xo <- x_obs[, 1]; yo <- x_obs[, 2]

        lam_sim <- readRDS(file = paste0(boot_vals_dir,"chi_plot_",labels[i],"_",labels[j],"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))
        lam_obs <- readRDS(file = paste0(boot_vals_dir,"chi_plot_obs_",labels[i],"_",labels[j],"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))

        lam_obs_mat <- sapply(1:n_boot,function(boot_ind){
          resample_inds = sample(1:length(xo),length(xo),replace=T)
          
          return(sapply(u_vec,tail_dep,x=xo[resample_inds],y=yo[resample_inds]))
        })
        
        lam_mat = sapply(1:n_boot,function(boot_ind){
          if(file.exists(paste0(boot_vals_dir,"chi_plot_",labels[i],"_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_model,".rds"))){
            return(readRDS(file = paste0(boot_vals_dir,"chi_plot_",labels[i],"_",labels[j],"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_model,".rds")))
          } else {
            print(boot_ind);print(time_window_ind);print(n_model)
            return(rep(NA,length(lam_sim)))
          }
        })

        plot(u_vec,pmin(lam_obs,1),type="l",lwd=2,col="grey2",cex.lab=1.3, xlab=expression(u),ylab=expression(chi(u)),ylim=c(0,1))
        polygon(c(rev(u_vec), u_vec), c(rev(apply(lam_mat,1,quantile,probs=0.025,na.rm=T)), apply(lam_mat,1,quantile,probs=0.975,na.rm=T)), col = rgb(0,205/255,0,alpha=0.15), border = NA)
        polygon(c(rev(u_vec), u_vec), c(rev(apply(lam_obs_mat,1,quantile,probs=0.025,na.rm=T)), apply(lam_obs_mat,1,quantile,probs=0.975,na.rm=T)), col = rgb(0.5,0.5,0.5,alpha=0.1), border = NA)
        lines(u_vec,lam_sim,lwd=2,col="green3")
        lines(u_vec,pmin(lam_obs,1),type="l",lwd=2,col="grey2")

        abline(h = 1,lwd=1,lty=2)

        abline(h = 0,lwd=1,lty=2)

      }
    }
  }

  mtext(paste0(percents[indices_grid$Var2[array_id]],"% subsample, ",min(time_periods[[time_window_ind]]),"-",max(time_periods[[time_window_ind]])), outer = TRUE, side = 3, line = 1.2, cex = 1.2)
  dev.off()
  
}

# Comparing different sample sizes ------------------------------------------------

indices_grid = expand.grid(c(1:length(time_periods)))

for(array_id in 1:nrow(indices_grid)){
  
  print(paste("Processing task", array_id))
  
  # Index of time window. We use middle of the period for optimisation
  time_window_ind = indices_grid$Var1[array_id]
  
  d = 4
  
  png(
    filename = paste0("PNGs/opt_fits_4d/rl_plot_upper_","tw_",time_window_ind,"_comp_ss.png"),
    width = 2800, height = 2800, res = 400
  )

  par(mar = c(4, 4.5, 4, 2))

  par(mgp = c(2.8, 0.7, 0))

  par(mfrow=c(2,2))

  for(i in 1:d){

    title = river_names[i]

    fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))

    emp_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_obs_",title,"_tw_",time_window_ind,".rds"))

    emp_RPs = readRDS(file=paste0(boot_vals_dir,"rp_ests_4d_upper_obs_",title,"_tw_",time_window_ind,".rds"))

    model_RLs_1 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))

    boot_rl_mat_1 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
      }
    })

    model_RLs_2 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))

    boot_rl_mat_2 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
      }
    })

    model_RLs_3 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))

    boot_rl_mat_3 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
      }
    })

    model_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_tw_",time_window_ind,".rds"))

    boot_rl_mat = sapply(1:n_boot,function(boot_ind){
      if(file.exists(paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
      } else {
        print(boot_ind);print(time_window_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })



    cols <- c(
      "#CC79A7",   # purple
      "#D55E00",  # vermillion
      "#009E73",  # green
      "grey" #grey
    )

    plot(
      fix_RPs, model_RLs,
      log = "x",
      type = "l", lwd=3, col = cols[4],
      xlab = "Return period (years)",
      ylab = expression("Discharge (m"^3*"/s)"),
      xaxt = "n",
      cex.lab=1.3,cex.main=1.5,
      ylim=range(emp_RLs,model_RLs,apply(boot_rl_mat,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_1,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_2,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_3,1,quantile,probs=c(0.025,0.975),na.rm=T)),
      main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," upper tail")
    )

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_1,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_1,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[1], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_2,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_2,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[2], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_3,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_3,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[3], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[4], alpha.f = 0.4), border = NA)

    lines(fix_RPs, model_RLs_1, lwd=3, col =  adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs_2, lwd=3, col = adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs_3, lwd=3, col = adjustcolor(cols[3], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs, lwd=3, col = adjustcolor(cols[4], alpha.f = 1, offset = c(-0.3, -0.3, -0.3, 0)))

    points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)

    ticks <- c(10, 100, 1000,10000)

    axis(1, at = ticks, labels = as.character(ticks))

    if(i == 1){
      legend("topleft",legend=c(paste0(percents[1:3],"% subsample"),"Full data set","Empirical"),
             col = c(adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[3], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[4], alpha.f = 1, offset = c(-0.3, -0.3, -0.3, 0)), "darkblue"),
             lty = c(rep(1,4), NA),          # line type for line; NA for no line
             pch = c(rep(NA,4), 16),         # point type for point; NA for no point
             lwd = c(rep(3,4), NA),
             pt.cex = 1,
             cex=1,
             bg="white")
    }

  }

  dev.off()

  png(
    filename = paste0("PNGs/opt_fits_4d/rl_plot_lower_","tw_",time_window_ind,"_comp_ss.png"),
    width = 2800, height = 2800, res = 400
  )

  par(mar = c(4, 4.5, 4, 2))

  par(mgp = c(2.8, 0.7, 0))

  par(mfrow=c(2,2))

  for(i in 1:d){

    title = river_names[i]

    fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))

    emp_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_obs_",title,"_tw_",time_window_ind,".rds"))

    emp_RPs = readRDS(file=paste0(boot_vals_dir,"rp_ests_4d_lower_obs_",title,"_tw_",time_window_ind,".rds"))

    model_RLs_1 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))

    boot_rl_mat_1 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
      }
    })

    model_RLs_2 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))

    boot_rl_mat_2 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
      }
    })

    model_RLs_3 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))

    boot_rl_mat_3 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
      }
    })

    model_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_tw_",time_window_ind,".rds"))

    boot_rl_mat = sapply(1:n_boot,function(boot_ind){
      if(file.exists(paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
      } else {
        print(boot_ind);print(time_window_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })



    cols <- c(
      "#CC79A7",   # purple
      "#D55E00",  # vermillion
      "#009E73",  # green
      "grey" #grey
    )

    plot(
      fix_RPs, model_RLs,
      log = "x",
      type = "l", lwd=3, col = cols[4],
      xlab = "Return period (years)",
      ylab = expression("Discharge (m"^3*"/s)"),
      xaxt = "n",
      cex.lab=1.3,cex.main=1.5,
      ylim=range(emp_RLs,model_RLs,apply(boot_rl_mat,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_1,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_2,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_3,1,quantile,probs=c(0.025,0.975),na.rm=T)),
      main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," lower tail")
    )

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_1,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_1,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[1], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_2,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_2,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[2], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_3,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_3,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[3], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[4], alpha.f = 0.4), border = NA)

    lines(fix_RPs, model_RLs_1, lwd=3, col =  adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs_2, lwd=3, col = adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs_3, lwd=3, col = adjustcolor(cols[3], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs, lwd=3, col = adjustcolor(cols[4], alpha.f = 1, offset = c(-0.3, -0.3, -0.3, 0)))

    points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)

    ticks <- c(10, 100, 1000,10000)

    axis(1, at = ticks, labels = as.character(ticks))

    # if(i == 1){
    #   legend("topright",legend=c(paste0(percents[1:3],"% subsample"),"Full data set","Empirical"),
    #          col = c(adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[3], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[4], alpha.f = 1, offset = c(-0.3, -0.3, -0.3, 0)), "darkblue"),
    #          lty = c(rep(1,4), NA),          # line type for line; NA for no line
    #          pch = c(rep(NA,4), 16),         # point type for point; NA for no point
    #          lwd = c(rep(3,4), NA),
    #          pt.cex = 1,
    #          cex=1,
    #          bg="white")
    # }

  }

  dev.off()
  
  png(
    filename = paste0("PNGs/opt_fits_4d/rl_plot_upper_","tw_",time_window_ind,"_comp_bigger_ss.png"),
    width = 2800, height = 2800, res = 400
  )
  
  par(mar = c(4, 4.5, 4, 2))
  
  par(mgp = c(2.8, 0.7, 0))
  
  par(mfrow=c(2,2))
  
  for(i in 1:d){
    
    title = river_names[i]
    
    fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))
    
    emp_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_obs_",title,"_tw_",time_window_ind,".rds"))
    
    emp_RPs = readRDS(file=paste0(boot_vals_dir,"rp_ests_4d_upper_obs_",title,"_tw_",time_window_ind,".rds"))
    
    model_RLs_1 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))
    
    boot_rl_mat_1 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds")))
      } else {
        print(n_models[4]);print(boot_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })
    
    model_RLs_2 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))
    
    boot_rl_mat_2 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds")))
      } else {
        print(n_models[5]);print(boot_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })
    
    model_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_tw_",time_window_ind,".rds"))
    
    boot_rl_mat = sapply(1:n_boot,function(boot_ind){
      if(file.exists(paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_upper_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
      } else {
        print(boot_ind);print(time_window_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })
    
    cols <- c(
      "#800080",   
      "orange", 
      "grey" 
    )
    
    plot(
      fix_RPs, model_RLs,
      log = "x",
      type = "l", lwd=3, col = cols[3],
      xlab = "Return period (years)",
      ylab = expression("Discharge (m"^3*"/s)"),
      xaxt = "n",
      cex.lab=1.3,cex.main=1.5,
      ylim=range(emp_RLs,model_RLs,apply(boot_rl_mat,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_1,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_2,1,quantile,probs=c(0.025,0.975),na.rm=T)),
      main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," upper tail")
    )
    
    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_1,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_1,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[1], alpha.f = 0.4), border = NA)
    
    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_2,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_2,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[2], alpha.f = 0.4), border = NA)
    
    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[3], alpha.f = 0.4), border = NA)
    
    lines(fix_RPs, model_RLs_1, lwd=3, col =  adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))
    
    lines(fix_RPs, model_RLs_2, lwd=3, col = adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))
    
    lines(fix_RPs, model_RLs, lwd=3, col = adjustcolor(cols[3], alpha.f = 1, offset = c(-0.3, -0.3, -0.3, 0)))
    
    points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)
    
    ticks <- c(10, 100, 1000,10000)
    
    axis(1, at = ticks, labels = as.character(ticks))
    
    if(i == 1){
      legend("topleft",legend=c(paste0(percents[4:5],"% subsample"),"Full data set","Empirical"),
             col = c(adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[3], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)), "darkblue"),
             lty = c(rep(1,3), NA),          # line type for line; NA for no line
             pch = c(rep(NA,3), 16),         # point type for point; NA for no point
             lwd = c(rep(3,3), NA),
             pt.cex = 1,
             cex=1,
             bg="white")
    }
    
  }
  
  dev.off()
  
  png(
    filename = paste0("PNGs/opt_fits_4d/rl_plot_lower_","tw_",time_window_ind,"_comp_bigger_ss.png"),
    width = 2800, height = 2800, res = 400
  )

  par(mar = c(4, 4.5, 4, 2))

  par(mgp = c(2.8, 0.7, 0))

  par(mfrow=c(2,2))

  for(i in 1:d){

    title = river_names[i]

    fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))

    emp_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_obs_",title,"_tw_",time_window_ind,".rds"))

    emp_RPs = readRDS(file=paste0(boot_vals_dir,"rp_ests_4d_lower_obs_",title,"_tw_",time_window_ind,".rds"))

    model_RLs_1 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))

    boot_rl_mat_1 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })

    model_RLs_2 = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,NULL,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))

    boot_rl_mat_2 = sapply(1:n_boot,function(boot_ind){
      if(file.exists( file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds")))
      } else {
        print(n_models[1]);print(boot_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })

    model_RLs = readRDS(file=paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_tw_",time_window_ind,".rds"))

    boot_rl_mat = sapply(1:n_boot,function(boot_ind){
      if(file.exists(paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
        return(readRDS(file = paste0(boot_vals_dir,"rl_ests_4d_lower_",title,"_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
      } else {
        print(boot_ind);print(time_window_ind)
        return(rep(NA,length(fix_RPs)))
      }
    })

    cols <- c(
      "#800080",   # purple
      "orange",  # vermillion
      "grey" #grey
    )

    plot(
      fix_RPs, model_RLs,
      log = "x",
      type = "l", lwd=3, col = cols[3],
      xlab = "Return period (years)",
      ylab = expression("Discharge (m"^3*"/s)"),
      xaxt = "n",
      cex.lab=1.3,cex.main=1.5,
      ylim=range(emp_RLs,model_RLs,apply(boot_rl_mat,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_1,1,quantile,probs=c(0.025,0.975),na.rm=T),apply(boot_rl_mat_2,1,quantile,probs=c(0.025,0.975),na.rm=T)),
      main=paste0(toupper(substring(title, 1, 1)), substring(title, 2)," lower tail")
    )

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_1,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_1,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[1], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat_2,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat_2,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[2], alpha.f = 0.4), border = NA)

    polygon(c(rev(fix_RPs), fix_RPs), c(rev(apply(boot_rl_mat,1,quantile,probs=0.025,na.rm=T)), apply(boot_rl_mat,1,quantile,probs=0.975,na.rm=T)), col = adjustcolor(cols[3], alpha.f = 0.4), border = NA)

    lines(fix_RPs, model_RLs_1, lwd=3, col =  adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs_2, lwd=3, col = adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)))

    lines(fix_RPs, model_RLs, lwd=3, col = adjustcolor(cols[3], alpha.f = 1, offset = c(-0.3, -0.3, -0.3, 0)))

    points(emp_RPs,emp_RLs,pch=16,col="darkblue",cex=1.4)

    ticks <- c(10, 100, 1000,10000)

    axis(1, at = ticks, labels = as.character(ticks))

    # if(i == 1){
    #   legend("topright",legend=c(paste0(percents[4:5],"% subsample"),"Full data set","Empirical"),
    #          col = c(adjustcolor(cols[1], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[2], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)),adjustcolor(cols[3], alpha.f = 1, offset = c(-0.1, -0.1, -0.1, 0)), "darkblue"),
    #          lty = c(rep(1,4), NA),          # line type for line; NA for no line
    #          pch = c(rep(NA,4), 16),         # point type for point; NA for no point
    #          lwd = c(rep(3,4), NA),
    #          pt.cex = 1,
    #          cex=1,
    #          bg="white")
    # }

  }

  dev.off()
  
}

# Analysing probabilities 4d with uncertainty ---------------------------------------------

l1_upper_ests = matrix(NA,nrow=length(l1_seq_upper),ncol=length(time_periods))

l1_upper_97.5 = matrix(NA,nrow=length(l1_seq_upper),ncol=length(time_periods))

l1_upper_2.5 = matrix(NA,nrow=length(l1_seq_upper),ncol=length(time_periods))

l1_lower_ests = matrix(NA,nrow=length(l1_seq_lower),ncol=length(time_periods))

l1_lower_97.5 = matrix(NA,nrow=length(l1_seq_lower),ncol=length(time_periods))

l1_lower_2.5 = matrix(NA,nrow=length(l1_seq_lower),ncol=length(time_periods))

upper_surv_ests = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

upper_surv_97.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

upper_surv_2.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

lower_cdf_ests = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

lower_cdf_97.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

lower_cdf_2.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

for(t in 1:length(time_periods)){
  
  time_window_ind = t
  
  probs_list = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,".rds"))
  
  boot_prob_mat = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
    } else {
      print(boot_ind);print(time_window_ind)
      return(NULL)
    }
  },simplify=F)
  
  l1_upper_ests[,t] = probs_list$l1_upper_ests
  
  l1_upper_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_upper_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  l1_upper_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_upper_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
  l1_lower_ests[,t] = probs_list$l1_lower_ests
  
  l1_lower_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_lower_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  l1_lower_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_lower_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
  upper_surv_ests[,t] = probs_list$upper_surv_ests
  
  upper_surv_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$upper_surv_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  upper_surv_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$upper_surv_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
  lower_cdf_ests[,t] = probs_list$lower_cdf_ests
  
  lower_cdf_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$lower_cdf_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  lower_cdf_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$lower_cdf_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
}

png(file=paste0("PNGs/opt_fits_4d/rp_comparison_cis_l1_upper.png"),
    width = 3500, height = 3500, res = 400
)

par(mfrow=c(1,1))

par(mar = c(5, 5, 4, 5))

options(scipen=999)

plot(1:4,1/(n_y*l1_upper_ests[1,]),
     type="l",col="grey3",
     cex.lab=1.7, cex.axis=1.6,cex.main=1.8,
     ylim=range(0,1/(n_y*l1_upper_2.5[1,]),1/(n_y*l1_upper_97.5[1,]),na.rm = T),
     xlab = "Time period",
     xaxt = "n",
     main = expression(Pr(S > 1000 ~ (m^3/s))),
     ylab="Return period (years)",frame.plot=F)

mtext("(a)", side = 3, line = 1, adj = 0, cex = 2)

polygon(c(rev(1:4), 1:4), c(rev(1/(n_y*l1_upper_2.5[1,])), 1/(n_y*l1_upper_97.5[1,])), col = 'grey80', border = NA)

lines(1:4,1/(n_y*l1_upper_ests[1,]), lwd=3, col = "grey3")

points(1:4,1/(n_y*l1_upper_ests[1,]),
       pch=16,col="grey3",cex=1.6)

axis(1, at = 1:4,cex.axis=1.6, labels = c(paste0(min(time_periods[[1]]),"-",max(time_periods[[1]])), paste0(min(time_periods[[2]]),"-",max(time_periods[[2]])), paste0(min(time_periods[[3]]),"-",max(time_periods[[3]])), paste0(min(time_periods[[4]]),"-",max(time_periods[[4]]))))

dev.off()

png(file=paste0("PNGs/opt_fits_4d/rp_comparison_cis_l1_lower.png"),
    width = 3500, height = 3500, res = 400
)

par(mfrow=c(1,1))

par(mar = c(5, 5, 4, 5))

options(scipen=999)

plot(1:4,1/(n_y*l1_lower_ests[1,]),
     type="l",col="grey3",
     cex.lab=1.7, cex.axis=1.6,cex.main=1.8,
     ylim=range(0,1/(n_y*l1_lower_2.5[1,]),1/(n_y*l1_lower_97.5[1,]),na.rm = T),
     xlab = "Time period",
     xaxt = "n",
     main = expression(Pr(S < 30 ~ (m^3/s))),
     ylab="Return period (years)",frame.plot=F)

mtext("(b)", side = 3, line = 1, adj = 0, cex = 2)

polygon(c(rev(1:4), 1:4), c(rev(1/(n_y*l1_lower_2.5[1,])), 1/(n_y*l1_lower_97.5[1,])), col = 'grey80', border = NA)

lines(1:4,1/(n_y*l1_lower_ests[1,]), lwd=3, col = "grey3")

points(1:4,1/(n_y*l1_lower_ests[1,]),
       pch=16,col="grey3",,cex=1.6)

axis(1, at = 1:4,cex.axis=1.6, labels = c(paste0(min(time_periods[[1]]),"-",max(time_periods[[1]])), paste0(min(time_periods[[2]]),"-",max(time_periods[[2]])), paste0(min(time_periods[[3]]),"-",max(time_periods[[3]])), paste0(min(time_periods[[4]]),"-",max(time_periods[[4]]))))

par(fig = c(0.3, 0.9, 0.4, 1), new = TRUE)  # Define inset position (x1, x2, y1, y2)
par(xpd = FALSE)
zoom_x <- 2:4  # Define x range for zoom
zoom_y <- range(0,1/(n_y*l1_lower_2.5[1, 3:4]), 1/(n_y*l1_lower_97.5[1, 3:4]), na.rm = TRUE)  # Define y range for zoom

plot(zoom_x, 1/(n_y*l1_lower_ests[1, zoom_x]),
     type = "l", col = "grey3",
     cex.lab = 1.4, cex.axis = 1.4, cex.main = 1,
     ylim = zoom_y, xlim=range(4,3),
     xlab = "",xaxt="n", ylab = "", main = "", frame.plot = TRUE)

polygon(c(rev(zoom_x), zoom_x), 
        c(rev(1/(n_y*l1_lower_2.5[1, zoom_x])), 1/(n_y*l1_lower_97.5[1, zoom_x])), 
        col = 'grey80', border = NA)
lines(zoom_x, 1/(n_y*l1_lower_ests[1, zoom_x]), lwd = 2, col = "grey3")
points(zoom_x, 1/(n_y*l1_lower_ests[1, zoom_x]), pch = 16, col = "grey3", cex = 1.2)

axis(1, at = 3:4,cex.axis=1.4, labels = c(paste0(min(time_periods[[3]]),"-",max(time_periods[[3]])), paste0(min(time_periods[[4]]),"-",max(time_periods[[4]]))))

dev.off()

png(file=paste0("PNGs/opt_fits_4d/rp_comparison_cis_upper_surv.png"),
    width = 3500, height = 3500, res = 400
)

par(mfrow=c(1,1))

par(mar = c(5, 5, 4, 5))

options(scipen=999)

plot(1:4,1/(n_y*upper_surv_ests[1,]),
     type="l",col="grey3",
     cex.lab=1.7, cex.axis=1.6,cex.main=1.8,
     ylim=range(0,1/(n_y*upper_surv_2.5[1,]),1/(n_y*upper_surv_97.5[1,]),na.rm = T),
     xlab = "Time period",
     xaxt = "n",
     main = expression(Pr(X[i] > x[i]^{1-1/N}*", i = 1, ..., 4")*", N = 10 (years)"),
     ylab="Return period (years)",frame.plot=F)

mtext("(c)", side = 3, line = 1, adj = 0, cex = 2)

polygon(c(rev(1:4), 1:4), c(rev(1/(n_y*upper_surv_2.5[1,])), 1/(n_y*upper_surv_97.5[1,])), col = 'grey80', border = NA)

lines(1:4,1/(n_y*upper_surv_ests[1,]), lwd=3, col = "grey3")

points(1:4,1/(n_y*upper_surv_ests[1,]),
       pch=16,col="grey3",cex=1.6)

axis(1, at = 1:4,cex.axis=1.6, labels = c(paste0(min(time_periods[[1]]),"-",max(time_periods[[1]])), paste0(min(time_periods[[2]]),"-",max(time_periods[[2]])), paste0(min(time_periods[[3]]),"-",max(time_periods[[3]])), paste0(min(time_periods[[4]]),"-",max(time_periods[[4]]))))

dev.off()

png(file=paste0("PNGs/opt_fits_4d/rp_comparison_cis_lower_cdf.png"),
    width = 3500, height = 3500, res = 400
)

par(mfrow=c(1,1))

par(mar = c(5, 5, 4, 5))

options(scipen=999)

plot(1:4,1/(n_y*lower_cdf_ests[1,]),
     type="l",col="grey3",
     cex.lab=1.7, cex.axis=1.6,cex.main=1.8,
     ylim=range(0,1/(n_y*lower_cdf_2.5[1,]),1/(n_y*lower_cdf_97.5[1,]),na.rm = T),
     xlab = "Time period",
     xaxt = "n",
     main = expression(Pr(X[i] < x[i]^{1/N}*", i = 1, ..., 4")*", N = 10 (years)"),
     ylab="Return period (years)",frame.plot=F)

mtext("(d)", side = 3, line = 1, adj = 0, cex = 2)

polygon(c(rev(1:4), 1:4), c(rev(1/(n_y*lower_cdf_2.5[1,])), 1/(n_y*lower_cdf_97.5[1,])), col = 'grey80', border = NA)

lines(1:4,1/(n_y*lower_cdf_ests[1,]), lwd=3, col = "grey3")

points(1:4,1/(n_y*lower_cdf_ests[1,]),
       pch=16,col="grey3",cex=1.6)

axis(1, at = 1:4,cex.axis=1.6, labels = c(paste0(min(time_periods[[1]]),"-",max(time_periods[[1]])), paste0(min(time_periods[[2]]),"-",max(time_periods[[2]])), paste0(min(time_periods[[3]]),"-",max(time_periods[[3]])), paste0(min(time_periods[[4]]),"-",max(time_periods[[4]]))))

dev.off()

# Analysing probabilities 4d with uncertainty + lower sample size ---------------------------------------------

l1_upper_ests = matrix(NA,nrow=length(l1_seq_upper),ncol=length(time_periods))

l1_upper_97.5 = matrix(NA,nrow=length(l1_seq_upper),ncol=length(time_periods))

l1_upper_2.5 = matrix(NA,nrow=length(l1_seq_upper),ncol=length(time_periods))

l1_lower_ests = matrix(NA,nrow=length(l1_seq_lower),ncol=length(time_periods))

l1_lower_97.5 = matrix(NA,nrow=length(l1_seq_lower),ncol=length(time_periods))

l1_lower_2.5 = matrix(NA,nrow=length(l1_seq_lower),ncol=length(time_periods))

upper_surv_ests = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

upper_surv_97.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

upper_surv_2.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

lower_cdf_ests = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

lower_cdf_97.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

lower_cdf_2.5 = matrix(NA,nrow=length(prob_RPs),ncol=length(time_periods))

for(t in 1:length(time_periods)){
  
  time_window_ind = t
  
  probs_list = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,".rds"))
  
  boot_prob_mat = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,".rds")))
    } else {
      print(boot_ind);print(time_window_ind)
      return(NULL)
    }
  },simplify=F)
  
  l1_upper_ests[,t] = probs_list$l1_upper_ests
  
  l1_upper_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_upper_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  l1_upper_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_upper_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
  l1_lower_ests[,t] = probs_list$l1_lower_ests
  
  l1_lower_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_lower_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  l1_lower_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$l1_lower_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
  upper_surv_ests[,t] = probs_list$upper_surv_ests
  
  upper_surv_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$upper_surv_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  upper_surv_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$upper_surv_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
  lower_cdf_ests[,t] = probs_list$lower_cdf_ests
  
  lower_cdf_97.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$lower_cdf_ests)
  },simplify=F)),2,quantile,probs=0.975)
  
  lower_cdf_2.5[,t] = apply(do.call(rbind,sapply(1:n_boot,function(boot_ind){
    return(boot_prob_mat[[boot_ind]]$lower_cdf_ests)
  },simplify=F)),2,quantile,probs=0.025)
  
}

png(file=paste0("PNGs/opt_fits_4d/prob_comparison_ss_l1_upper.png"),
    width = 6000, height = 1500, res = 400
)

options(scipen=999)

par(mfrow=c(1,4))

for(t in 1:length(time_periods)){

  time_window_ind = t

  probs_1 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$l1_upper_ests

  boot_prob_mat_1 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$l1_upper_ests)
    } else {
      print(boot_ind);print(t);print(n_models[1])
      return(NULL)
    }
  },simplify=F)

  probs_2 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$l1_upper_ests

  boot_prob_mat_2 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$l1_upper_ests)
    } else {
      print(boot_ind);print(t);print(n_models[2])
      return(NULL)
    }
  },simplify=F)

  probs_3 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$l1_upper_ests

  boot_prob_mat_3 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$l1_upper_ests)
    } else {
      print(boot_ind);print(t);print(n_models[3])
      return(NULL)
    }
  },simplify=F)
  
  probs_4 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$l1_upper_ests
  
  boot_prob_mat_4 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$l1_upper_ests)
    } else {
      print(t);print(n_models[4]);print(boot_ind)
      return(NULL)
    }
  },simplify=F)
  
  probs_5 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$l1_upper_ests
  
  boot_prob_mat_5 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$l1_upper_ests)
    } else {
      print(boot_ind);print(t);print(n_models[5])
      return(NULL)
    }
  },simplify=F)

  boot_prob_mat_1 = t(do.call(rbind,boot_prob_mat_1))

  boot_prob_mat_2 = t(do.call(rbind,boot_prob_mat_2))

  boot_prob_mat_3 = t(do.call(rbind,boot_prob_mat_3))
  
  boot_prob_mat_4 = t(do.call(rbind,boot_prob_mat_4))
  
  boot_prob_mat_5 = t(do.call(rbind,boot_prob_mat_5))

  par(mar = c(5, 5, 4, 5))
  
  plot(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,l1_upper_ests[1,t]),
       type="l",col="grey3",
       cex.lab=2, cex.axis=1.4,cex.main=2,
       ylim=range(0,l1_upper_2.5[1,],l1_upper_97.5[1,],boot_prob_mat_1,boot_prob_mat_2,boot_prob_mat_3,boot_prob_mat_4,boot_prob_mat_5),
       xlab = "",#xlab = ifelse(t==3,expression(bold("Subsample")),""),
       main = paste0("Time period: ",min(time_periods[[t]]), "-", max(time_periods[[t]])),
       xaxt = "n",
       # sub = expression(Pr(S > 3000 ~ (m^3/s))),
       ylab=expression(p),frame.plot=F)
  
  
  polygon(c(rev(1:6), 1:6), c(rev(c(apply(boot_prob_mat_1,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.025,na.rm=T),l1_upper_2.5[1,t])),c(apply(boot_prob_mat_1,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.975,na.rm=T),l1_upper_97.5[1,t])), col = 'grey80', border = NA)
  
  lines(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,l1_upper_ests[1,t]), lwd=3, col = "grey3")
  
  points(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,l1_upper_ests[1,t]),
         pch=16,col="grey3",cex=1.6)
  
  axis(1, at = 1:6,cex.axis=1.4, labels = c(paste0(percents,"%"),paste0("Full")))

}

dev.off()

png(file=paste0("PNGs/opt_fits_4d/prob_comparison_ss_l1_lower.png"),
    width = 6000, height = 1500, res = 400
)

options(scipen=999)

par(mfrow=c(1,4))

for(t in 1:length(time_periods)){
  
  time_window_ind = t
  
  probs_1 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$l1_lower_ests
  
  boot_prob_mat_1 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$l1_lower_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_2 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$l1_lower_ests
  
  boot_prob_mat_2 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$l1_lower_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_3 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$l1_lower_ests
  
  boot_prob_mat_3 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$l1_lower_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_4 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$l1_lower_ests
  
  boot_prob_mat_4 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$l1_lower_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_5 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$l1_lower_ests
  
  boot_prob_mat_5 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$l1_lower_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  boot_prob_mat_1 = t(do.call(rbind,boot_prob_mat_1))
  
  boot_prob_mat_2 = t(do.call(rbind,boot_prob_mat_2))
  
  boot_prob_mat_3 = t(do.call(rbind,boot_prob_mat_3))
  
  boot_prob_mat_4 = t(do.call(rbind,boot_prob_mat_4))
  
  boot_prob_mat_5 = t(do.call(rbind,boot_prob_mat_5))
  
  par(mar = c(5, 5, 4, 5))
  
  plot(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,l1_lower_ests[1,t]),
       type="l",col="grey3",
      cex.lab=2, cex.axis=1.4,cex.main=2,
       ylim=range(0,l1_lower_2.5[1,],l1_lower_97.5[1,],l1_lower_ests[1,t],boot_prob_mat_1,boot_prob_mat_2,boot_prob_mat_3,boot_prob_mat_4,boot_prob_mat_5),
       xlab = "",#xlab = ifelse(t==3,expression(bold("Subsample")),""),
       main = paste0("Time period: ",min(time_periods[[t]]), "-", max(time_periods[[t]])),
       xaxt = "n",
       # sub = expression(Pr(S > 3000 ~ (m^3/s))),
       ylab=expression(p),frame.plot=F)
  
  
  polygon(c(rev(1:6), 1:6), c(rev(c(apply(boot_prob_mat_1,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.025,na.rm=T),l1_lower_2.5[1,t])),c(apply(boot_prob_mat_1,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.975,na.rm=T),l1_lower_97.5[1,t])), col = 'grey80', border = NA)
  
  lines(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,l1_lower_ests[1,t]), lwd=3, col = "grey3")
  
  points(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,l1_lower_ests[1,t]),
         pch=16,col="grey3",cex=1.6)
  
  axis(1, at = 1:6,cex.axis=1.4, labels = c(paste0(percents,"%"),paste0("Full")))
  
}

dev.off()

png(file=paste0("PNGs/opt_fits_4d/prob_comparison_ss_upper_surv.png"),
    width = 6000, height = 1500, res = 400
)

options(scipen=999)

par(mfrow=c(1,4))

for(t in 1:length(time_periods)){
  
  time_window_ind = t
  
  probs_1 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$upper_surv_ests
  
  boot_prob_mat_1 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$upper_surv_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_2 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$upper_surv_ests
  
  boot_prob_mat_2 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$upper_surv_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_3 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$upper_surv_ests
  
  boot_prob_mat_3 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$upper_surv_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_4 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$upper_surv_ests
  
  boot_prob_mat_4 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$upper_surv_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_5 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$upper_surv_ests
  
  boot_prob_mat_5 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$upper_surv_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  boot_prob_mat_1 = t(do.call(rbind,boot_prob_mat_1))
  
  boot_prob_mat_2 = t(do.call(rbind,boot_prob_mat_2))
  
  boot_prob_mat_3 = t(do.call(rbind,boot_prob_mat_3))
  
  boot_prob_mat_4 = t(do.call(rbind,boot_prob_mat_4))
  
  boot_prob_mat_5 = t(do.call(rbind,boot_prob_mat_5))
  
  par(mar = c(5, 5, 4, 5))
  
  plot(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,upper_surv_ests[1,t]),
       type="l",col="grey3",
      cex.lab=2, cex.axis=1.4,cex.main=2,
       ylim=range(0,upper_surv_2.5[1,],upper_surv_97.5[1,],boot_prob_mat_1,boot_prob_mat_2,boot_prob_mat_3,boot_prob_mat_4,boot_prob_mat_5),
       xlab = "",#xlab = ifelse(t==3,expression(bold("Subsample")),""),
       main = paste0("Time period: ",min(time_periods[[t]]), "-", max(time_periods[[t]])),
       xaxt = "n",
       # sub = expression(Pr(S > 3000 ~ (m^3/s))),
       ylab=expression(p),frame.plot=F)
  
  
  polygon(c(rev(1:6), 1:6), c(rev(c(apply(boot_prob_mat_1,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.025,na.rm=T),upper_surv_2.5[1,t])),c(apply(boot_prob_mat_1,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.975,na.rm=T),upper_surv_97.5[1,t])), col = 'grey80', border = NA)
  
  lines(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,upper_surv_ests[1,t]), lwd=3, col = "grey3")
  
  points(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,upper_surv_ests[1,t]),
         pch=16,col="grey3",cex=1.6)
  
  axis(1, at = 1:6,cex.axis=1.4, labels = c(paste0(percents,"%"),paste0("Full")))
  
}

dev.off()

png(file=paste0("PNGs/opt_fits_4d/prob_comparison_ss_lower_cdf.png"),
    width = 6000, height = 1500, res = 400
)

options(scipen=999)

par(mfrow=c(1,4))

for(t in 1:length(time_periods)){
  
  time_window_ind = t
  
  probs_1 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$lower_cdf_ests
  
  boot_prob_mat_1 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[1],".rds"))$lower_cdf_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_2 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$lower_cdf_ests
  
  boot_prob_mat_2 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[2],".rds"))$lower_cdf_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_3 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$lower_cdf_ests
  
  boot_prob_mat_3 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[3],".rds"))$lower_cdf_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_4 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$lower_cdf_ests
  
  boot_prob_mat_4 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[4],".rds"))$lower_cdf_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  probs_5 = readRDS(file=paste0(boot_vals_dir,"model_prob_ests_4d","_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$lower_cdf_ests
  
  boot_prob_mat_5 = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))){
      return(readRDS(file = paste0(boot_vals_dir,"model_prob_ests_4d","_bootind_",boot_ind,"_tw_",time_window_ind,"_nmodel_",n_models[5],".rds"))$lower_cdf_ests)
    } else {
      print(boot_ind);print(t)
      return(NULL)
    }
  },simplify=F)
  
  boot_prob_mat_1 = t(do.call(rbind,boot_prob_mat_1))
  
  boot_prob_mat_2 = t(do.call(rbind,boot_prob_mat_2))
  
  boot_prob_mat_3 = t(do.call(rbind,boot_prob_mat_3))
  
  boot_prob_mat_4 = t(do.call(rbind,boot_prob_mat_4))
  
  boot_prob_mat_5 = t(do.call(rbind,boot_prob_mat_5))
  
  par(mar = c(5, 5, 4, 5))
  
  plot(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,lower_cdf_ests[1,t]),
       type="l",col="grey3",
      cex.lab=2, cex.axis=1.4,cex.main=2,
       ylim=range(0,lower_cdf_2.5[1,],lower_cdf_97.5[1,],boot_prob_mat_1,boot_prob_mat_2,boot_prob_mat_3,boot_prob_mat_4,boot_prob_mat_5),
       xlab = "",#xlab = ifelse(t==3,expression(bold("Subsample")),""),
       main = paste0("Time period: ",min(time_periods[[t]]), "-", max(time_periods[[t]])),
       xaxt = "n",
       # sub = expression(Pr(S > 3000 ~ (m^3/s))),
       ylab=expression(p),frame.plot=F)
  
  
  polygon(c(rev(1:6), 1:6), c(rev(c(apply(boot_prob_mat_1,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.025,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.025,na.rm=T),lower_cdf_2.5[1,t])),c(apply(boot_prob_mat_1,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_2,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_3,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_4,1,quantile,probs=0.975,na.rm=T),apply(boot_prob_mat_5,1,quantile,probs=0.975,na.rm=T),lower_cdf_97.5[1,t])), col = 'grey80', border = NA)
  
  lines(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,lower_cdf_ests[1,t]), lwd=3, col = "grey3")
  
  points(1:6,c(probs_1,probs_2,probs_3,probs_4,probs_5,lower_cdf_ests[1,t]),
         pch=16,col="grey3",cex=1.6)
  
  axis(1, at = 1:6,cex.axis=1.4, labels = c(paste0(percents,"%"),paste0("Full")))
  
}

dev.off()
