## ---------------------------Title---------------------------------------------
# CIDA Shelf Life Model

## ---------------------------Description---------------------------------------
# Project: CIDA Spinach & Bacterial growth modeling (spinach)

# Description: Monte Carlo simulation of bacterial growth on spinach

## ---------------------------Packages------------------------------------------

#   Loading packages
library(dplyr); library(MASS); library(VGAM)

## ---------------------------Seed----------------------------------------------

#Set seed
set.seed(1)

## ---------------------------Inputs--------------------------------------------
#Initial_count
min_max <- read.csv("model_inputs/apc_day0_conc_dist.csv", header = TRUE)

#Load in data

apc_param <- read.csv("model_inputs/primary_model_parameters_apc_averaged.csv", header = TRUE)

mumax_secon <- read.csv("model_inputs/secondary_model_mumax_apc_averaged.csv", header = TRUE)

#Initial Count Distribution: Fitting a uniform distribution to the data

initial_count_dist <- function(n){
  ans <- runif(n, 
               min = min_max$min[min_max$lot == 8], 
               max = min_max$max[min_max$lot == 8])
  return(ans)
}

#Storage Temperature Distribution
#temp_location <- 4.06
#temp_shape <- 2.31

#temp_dist <- function(n) {
#  ans <- rlaplace(n, temp_location, temp_shape)
#  for(i in 1:n)
#  while (ans[i] < 2 | ans[i] > 14) {
#    ans[i] <- rlaplace(1, temp_location, temp_shape)
#  }
#  return(ans)
#}

#Number of strains per package


#Primary Growth Models 
#Buchanan_no_lag
buchanan_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10n0 + (day <= ((log10nmax - log10n0) * log(10) / mumax)) * mumax * day / log(10) + (day > ((log10nmax - log10n0) * log(10) / mumax)) * (log10nmax - log10n0)
}

#Secondary model 

#Mumax
NewMu <- function(temp, b, temp_min){
  if(temp > temp_min){
    sqrt_mu <- b*(temp - temp_min)
    ans <- sqrt_mu^2
    return(ans)
  } else {
    ans <- 0 
    return(ans)
  }
}


## ---------------------------Simulation Overview-------------------------------

n_lots <- 1
n_packages <- 100
n_sim <- n_packages*n_lots
n_day <- 30

## ---------------Assigning Model Inputs----------------------------------------

#Assigning initial count, and temperature to each package
packages_sim <- data.frame(matrix(nrow = n_sim, ncol = 4))
colnames(packages_sim) <- c("package", "initial_count", "storage_temp", "lot")

packages_sim$package <- 1:n_packages
packages_sim$initial_count <- initial_count_dist(n_sim)
packages_sim$storage_temp <- 10 #temp_dist(n_sim)
packages_sim$lot <- "overall"


#Expanding the dataframe to store observations
packages_sim_exp <- packages_sim %>%
  group_by(package) %>%
  mutate(lag = NA, mumax = NA, nmax = NA) 

#Assign lag to each package
for(i in 1:nrow(packages_sim_exp)){
  packages_sim_exp$lag[i] <- 0
}

#Assign mumax to each package
for(i in 1:nrow(packages_sim_exp)){
  index <- which(mumax_secon$val_data_lot == packages_sim_exp$lot[i])
  lot_mu <- NewMu(packages_sim_exp$storage_temp[i], mumax_secon$b_mu[index], mumax_secon$temp_min_mu[index])
  packages_sim_exp$mumax[i] <- lot_mu
}

mean_nmax <- mean(filter(apc_param)$nmax)
  
#Assign nmax to each package
for(i in 1:nrow(packages_sim_exp)){
  packages_sim_exp$nmax[i] <- mean_nmax
}

## ----------------Assessing Growth Over Shelf Life-----------------------------

#Creating a dataframe, to assess microbial growth over shelf life 
shelf_life_sim <- packages_sim_exp %>%
  group_by(package) %>%
  slice(rep(1, times = 30)) %>%
  mutate(day = 1:30)
  
shelf_life_sim <- shelf_life_sim %>%
  mutate(count = buchanan_no_lag(day, initial_count, nmax, mumax))

## --------------------Plotting & Analyzing Outcome-----------------------------

count_summary_by_day <- shelf_life_sim %>%
  dplyr::select(c("package", "initial_count", "storage_temp", "day", "count")) %>%
  group_by(day, package) %>%
  mutate(growth_model = "mm-apc",
         loc = "us",
         sampling_code = 8)
  
## --------------------Exporting data-------------------------------------------

file_name <- paste0("outputs/count_by_day_mm_apc_us.csv")

write.table(count_summary_by_day, file_name, sep = ",", 
            col.names = !file.exists(file_name), append = T, row.names = FALSE)

file_name <- paste0("outputs/overall_count.csv")

write.table(count_summary_by_day, file_name, sep = ",", 
            col.names = !file.exists(file_name), append = T, row.names = FALSE)


