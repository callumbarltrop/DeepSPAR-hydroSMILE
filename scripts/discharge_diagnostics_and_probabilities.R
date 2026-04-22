# Clear workspace 

rm(list=ls())

# Data --------------------------------------------------------------------

river_names = c("ammer","iller","lech","loisach")

time_periods = list(1980:2009,2010:2039,2040:2069,2070:2099)

n_boot = 100

d = 4

# Computing model diagnostics with uncertainty bounds ------------------------------------------------

for(tw_ind in 1:4){
  
  print(paste("Time window ", tw_ind))

  obs_quants = readRDS(file = paste0("../computed_model_vals/gpd_qq_quants_tw_",tw_ind,".rds"))

  n_p = length(obs_quants)
  ps = (1:n_p)/(n_p + 1)
  theor_quants = qexp(ps)

  boot_quant_mat = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0("../computed_model_vals/gpd_qq_quants_","bootind_",boot_ind,"_tw_",tw_ind,".rds"))){
      return(readRDS(file = paste0("../computed_model_vals/gpd_qq_quants_","bootind_",boot_ind,"_tw_",tw_ind,".rds")))
    } else {
      print(boot_ind);print(tw_ind)
      return(rep(NA,length(ps)))
    }
  })

  png(
    filename = paste0("../figures/gpd_qq_plot_tw_",tw_ind,".png"),
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
    main = paste0("GPD QQ plot, ",min(time_periods[[tw_ind]]),"-",max(time_periods[[tw_ind]])),
    sub=""
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
    file = paste0("../figures/scatterplot_tw_", tw_ind,".png"),
    width = 3500, height = 3500, res = 400
  )

  # 4x4 layout, all panels square, equal size
  layout(matrix(1:16, nrow = 4, ncol = 4), widths = rep(1,4), heights = rep(1,4))
  par(mar = c(4, 4.5, 2, 1))  # inner margin (bottom, left, top, right)
  par(oma = c(1, 1, 4, 1))  # outer margin for overall labels / title

  for (j in 1:d) {
    for (i in 1:d) {

      if (i == j) {
        # ---------------- Diagonal: marginal QQ plots ----------------
        obs_quants = readRDS(file = paste0("../computed_model_vals/marg_qq_quants_obs_",labels[j],"_tw_",tw_ind,".rds"))
        sim_quants = readRDS(file = paste0("../computed_model_vals/marg_qq_quants_",labels[j],"_tw_",tw_ind,".rds"))

        boot_quant_mat = sapply(1:n_boot,function(boot_ind){
          if(file.exists(paste0("../computed_model_vals/marg_qq_quants_",labels[j],"_bootind_",boot_ind,"_tw_",tw_ind,".rds"))){
            return(readRDS(file = paste0("../computed_model_vals/marg_qq_quants_",labels[j],"_bootind_",boot_ind,"_tw_",tw_ind,".rds")))
          } else {
            print(boot_ind);print(tw_ind)
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

        x_sim = readRDS(file = paste0("../computed_model_vals/scatterplot_",labels[i],"_",labels[j],"_tw_",tw_ind,".rds"))

        x_obs = readRDS(file = paste0("../computed_model_vals/scatterplot_obs_",labels[i],"_",labels[j],"_tw_",tw_ind,".rds"))

        label_text <- paste0(toupper(substring(labels[j], 1, 1)), substring(labels[j], 2))
        if(i == 1){xlab=bquote(.(label_text) ~ "discharge" ~ (m^3/s))}else{xlab=""}

        label_text <- paste0(toupper(substring(labels[i], 1, 1)), substring(labels[i], 2))
        if(i == j-1){ylab=bquote(.(label_text) ~ "discharge" ~ (m^3/s))}else{ylab=""}

        plot(x_sim, pch = 16, col = rgb(0, 1, 0, 0.35), cex = 0.7,xlab=xlab,ylab=ylab,cex.lab=1.3)
        points(x_obs, pch = 16, col = rgb(0.5, 0.5, 0.5, 1), cex = 0.7)

      } else {
        # ---------------- Lower diagonal: chi dependence coefficient on scale ----------------
        x_obs = readRDS(file = paste0("../computed_model_vals/scatterplot_obs_",labels[j],"_",labels[i],"_tw_",tw_ind,".rds"))
        xo <- x_obs[, 1]; yo <- x_obs[, 2]

        lam_sim <- readRDS(file = paste0("../computed_model_vals/chi_plot_",labels[i],"_",labels[j],"_tw_",tw_ind,".rds"))
        lam_obs <- readRDS(file = paste0("../computed_model_vals/chi_plot_obs_",labels[i],"_",labels[j],"_tw_",tw_ind,".rds"))

        lam_obs_mat <- sapply(1:n_boot,function(boot_ind){
          resample_inds = sample(1:length(xo),length(xo),replace=T)

          return(sapply(u_vec,tail_dep,x=xo[resample_inds],y=yo[resample_inds]))
        })

        lam_mat = sapply(1:n_boot,function(boot_ind){
          if(file.exists(paste0("../computed_model_vals/chi_plot_",labels[i],"_",labels[j],"_bootind_",boot_ind,"_tw_",tw_ind,".rds"))){
            return(readRDS(file = paste0("../computed_model_vals/chi_plot_",labels[i],"_",labels[j],"_bootind_",boot_ind,"_tw_",tw_ind,".rds")))
          } else {
            print(boot_ind);print(tw_ind)
            return(rep(NA,length(lam_sim)))
          }
        })

        plot(u_vec,lam_obs,type="l",lwd=2,col="grey2",cex.lab=1.3, xlab=expression(u),ylab=expression(chi(u)),ylim=c(0,1))
        polygon(c(rev(u_vec), u_vec), c(rev(apply(lam_mat,1,quantile,probs=0.025,na.rm=T)), apply(lam_mat,1,quantile,probs=0.975,na.rm=T)), col = rgb(0,205/255,0,alpha=0.15), border = NA)
        polygon(c(rev(u_vec), u_vec), c(rev(apply(lam_obs_mat,1,quantile,probs=0.025,na.rm=T)), apply(lam_obs_mat,1,quantile,probs=0.975,na.rm=T)), col = rgb(0.5,0.5,0.5,alpha=0.1), border = NA)
        lines(u_vec,lam_sim,lwd=2,col="green3")
        lines(u_vec,lam_obs,type="l",lwd=2,col="grey")

        abline(h = 1,lwd=1,lty=2)

        abline(h = 0,lwd=1,lty=2)

      }
    }
  }

  mtext(paste0(min(time_periods[[tw_ind]]),"-",max(time_periods[[tw_ind]])), outer = TRUE, side = 3, line = 1.2, cex = 1.2)
  dev.off()
  
  png(
    filename = paste0("../figures/rl_plot_upper_tw_",tw_ind,".png"),
    width = 2800, height = 2800, res = 400
  )

  par(mar = c(4, 4.5, 4, 2))
  
  par(mgp = c(2.8, 0.7, 0))

  par(mfrow=c(2,2))

  for(i in 1:d){

    title = river_names[i]

    fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))
    
    emp_RLs = readRDS(file=paste0("../computed_model_vals/rl_ests_upper_obs_",title,"_tw_",tw_ind,".rds"))
    
    emp_RPs = readRDS(file=paste0("../computed_model_vals/rp_ests_upper_obs_",title,"_tw_",tw_ind,".rds"))
      
    model_RLs = readRDS(file=paste0("../computed_model_vals/rl_ests_upper_",title,"_tw_",tw_ind,".rds"))

    boot_rl_mat = sapply(1:n_boot,function(boot_ind){
      if(file.exists(paste0("../computed_model_vals/rl_ests_upper_",title,"_bootind_",boot_ind,"_tw_",tw_ind,".rds"))){
        return(readRDS(file = paste0("../computed_model_vals/rl_ests_upper_",title,"_bootind_",boot_ind,"_tw_",tw_ind,".rds")))
      } else {
        print(boot_ind);print(tw_ind)
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
  
  png(
    filename = paste0("../figures/rl_plot_lower_tw_",tw_ind,".png"),
    width = 2800, height = 2800, res = 400
  )

  par(mar = c(4, 4.5, 4, 2))
  
  par(mgp = c(2.8, 0.7, 0))

  par(mfrow=c(2,2))

  for(i in 1:d){

    title = river_names[i]

    fix_RPs = exp(seq(log(10),log(10000),length.out = 3000))
    
    emp_RLs = readRDS(file=paste0("../computed_model_vals/rl_ests_lower_obs_",title,"_tw_",tw_ind,".rds"))

    emp_RPs = readRDS(file=paste0("../computed_model_vals/rp_ests_lower_obs_",title,"_tw_",tw_ind,".rds"))
    
    model_RLs = readRDS(file=paste0("../computed_model_vals/rl_ests_lower_",title,"_tw_",tw_ind,".rds"))

    boot_rl_mat = sapply(1:n_boot,function(boot_ind){
      if(file.exists(paste0("../computed_model_vals/rl_ests_lower_",title,"_bootind_",boot_ind,"_tw_",tw_ind,".rds"))){
        return(readRDS(file = paste0("../computed_model_vals/rl_ests_lower_",title,"_bootind_",boot_ind,"_tw_",tw_ind,".rds")))
      } else {
        print(boot_ind);print(tw_ind)
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

# Comparing probabilities over time periods with uncertainty ---------------------------------------------

l1_seq_upper = 1000

l1_seq_lower = 30

prob_RPs = 10

n_y = 52 

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
  
  tw_ind = t
  
  probs_list = readRDS(file=paste0("../computed_model_vals/model_prob_ests_tw_",tw_ind,".rds"))
  
  boot_prob_mat = sapply(1:n_boot,function(boot_ind){
    if(file.exists(paste0("../computed_model_vals/model_prob_ests_bootind_",boot_ind,"_tw_",tw_ind,".rds"))){
      return(readRDS(file = paste0("../computed_model_vals/model_prob_ests_bootind_",boot_ind,"_tw_",tw_ind,".rds")))
    } else {
      print(boot_ind);print(tw_ind)
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

png(file=paste0("../figures/rp_comparison_l1_upper.png"),
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

png(file=paste0("../figures/rp_comparison_l1_lower.png"),
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

png(file=paste0("../figures/rp_comparison_upper_surv.png"),
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

png(file=paste0("../figures/rp_comparison_lower_cdf.png"),
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
