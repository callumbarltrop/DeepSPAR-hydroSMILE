# Clear R environment 

rm(list=ls())

# Load preamble and GGally package for plotting  

source("preamble.R")

if (!require("GGally", character.only = TRUE)) {
  install.packages("GGally", dependencies = TRUE)
  library("GGally", character.only = TRUE)
}

# Combining all weekly maxima into one data frame --------------------------

# Vector of station names

river_names = c("ammer","iller","lech","loisach")

# Load in each data frame and use cbind to combine them into one

discharge_df = readRDS(paste0("../data/weekly_maxima/",river_names[1],"_cmb_weekly_max.rds"))

for(i in 2:length(river_names)){

  more_discharge_df = readRDS(paste0("../data/weekly_maxima/",river_names[i],"_cmb_weekly_max.rds"))

  discharge_df = cbind(discharge_df,more_discharge_df$weekly_max)

}

# Specify names for each river 

names(discharge_df)[3:6] = paste0(river_names,"_q_wm")

# Save combined data frame

saveRDS(discharge_df,file="../data/weekly_maxima/all_wm_df.rds")

# Plotting data on original scale -----------------------------------------------------------

# Load in all weekly maxima 

discharge_df = readRDS(file="../data/weekly_maxima/all_wm_df.rds")

# Specify time periods of study 

time_periods = list(1980:2009,2010:2039,2040:2069,2070:2099)

# Vector of station names

river_names = c("ammer","iller","lech","loisach")

# Extract the year of each 'observational' date

years = as.numeric(format(discharge_df$week_start, "%Y"))

# Plotting data on original scale for each time period 

for(i in 1:length(time_periods)){

  # Discharge data for time period i 
  
  tp_discharge_df = discharge_df[years %in% time_periods[[i]],]
  
  # ggplot of each pairwise combination 
  
  p <- ggpairs(tp_discharge_df[,3:6],
          lower = list(continuous = wrap("points", alpha = 0.4,color="blue4", size = 0.5)),
          diag = list(continuous = wrap("densityDiag")),
          upper = list(continuous = wrap("cor")),
          columnLabels = river_names) +
      labs(title = paste0("Pairwise relationships of river discharges, time period ",time_periods[[i]][1],"-",time_periods[[i]][length(time_periods[[i]])]),
         subtitle = "Lower: Scatterplots | Upper: Correlations | Diagonal: Densities") +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12))

  ggsave(paste0("../figures/pairwise_plots_tw_",i,".png"), plot = p, width = 10, height = 10, dpi = 300)

}

# Pre-processing procedure ----------------------------------------------------------

# The code below moves the data from the positive real line to the whole real line. It also centers the data at (0,0,0,0). 
# We save the data sets on the original and transformed scales, alongside the values associated with the transformation(s). 

std_tw_dfs = list()

std_sigmas = list()

std_mus = list()

os_tw_dfs = list()

for(i in 1:length(time_periods)){
  
  # Discharge data for time period i 
  
  tp_discharge_df = discharge_df[years %in% time_periods[[i]],]
  
  os_tw_dfs[[i]] = tp_discharge_df
  
  # Transforming to the whole real line 
  
  std_sigmas[[i]] = apply(tp_discharge_df[,3:6],2,sd)
  
  tp_discharge_df[,3:6] = apply(tp_discharge_df[,3:6],2,function(x){return( x/sd(x) )})
  
  tp_discharge_df[,3:6] = apply(tp_discharge_df[,3:6],2,function(x){return( log(exp(x) - 1) )})
  
  # Computing the geometric median 
  
  std_mus[[i]] = geometric_median(tp_discharge_df[,3:6])
  
  # Centering the data
  
  tp_discharge_df[,3:6] = apply(rbind(std_mus[[i]],tp_discharge_df[,3:6]),2,function(x){return( x[-1] - x[1] )})
  
  std_tw_dfs[[i]] <- tp_discharge_df[,3:6]

  # ggplot of each pairwise combination 
  
  p <- ggpairs(std_tw_dfs[[i]],
               lower = list(continuous = wrap("points", alpha = 0.4,color="blue4", size = 0.5)),
               diag = list(continuous = wrap("densityDiag")),
               upper = list(continuous = wrap("cor")),
               columnLabels = river_names) +
    labs(title = paste0("Pairwise relationships of pre-processed river discharges, time period ",time_periods[[i]][1],"-",time_periods[[i]][length(time_periods[[i]])]),
         subtitle = "Lower: Scatterplots | Upper: Correlations | Diagonal: Densities") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12))

  ggsave(paste0("../figures/pairwise_plots_preprocessed_tw_",i,".png"), plot = p, width = 10, height = 10, dpi = 300)

}

saveRDS(std_sigmas,file= "../data/preprocessed_data/std_sigmas.rds")

saveRDS(std_mus,file= "../data/preprocessed_data/std_mus.rds")

saveRDS(std_tw_dfs,file="../data/preprocessed_data/std_tw_dfs.rds")

saveRDS(os_tw_dfs,file="../data/preprocessed_data/os_tw_dfs.rds")

