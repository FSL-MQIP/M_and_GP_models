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
min_max <- read.csv("model_inputs/china_starting_conc.csv", header = TRUE)

#Load in data

strain_param <- read.csv("model_inputs/primary_model_parameters_strain_averaged.csv", header = TRUE)

mumax_secon <- read.csv("model_inputs/secondary_model_mumax_strain_averaged.csv", header = TRUE)

st_frequency <- read.csv("model_inputs/walmart_s1-s3_apc_st_frequency.csv", header = TRUE)

#Change strain parameters to lowercase 
strain_param$isolate <- tolower(strain_param$isolate)

mumax_secon$isolate <- tolower(mumax_secon$isolate)

#Change the formatting of a variable in the ST frequency 
st_frequency$cida_isolates <- substr(st_frequency$cida_isolates, 1, 8)
st_frequency$cida_isolates <- tolower(st_frequency$cida_isolates)

#Initial Count Distribution: Fitting a uniform distribution to the data

initial_count_dist <- function(n){
  ans <- runif(n, 
               min = 5.5, #min_max$min[min_max$sampling == "R7"], 
               max = 6.5)#min_max$max[min_max$sampling == "R7"])
  return(ans)
}

#ST Distribution
st_frequency <- st_frequency %>%
  filter(cida_isolates != "s12-0184")
isolates <- st_frequency$cida_isolates
st_prev <- st_frequency$frequency

st_dist <- function(n) {
ans <- sample(isolates, size = n, prob = st_prev, replace = F)
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

pop_dist <- function(n){
  ans <- runif(n, min = 1, max = 5)
  ans_rounded <- round(ans, digits = 0)
  return(ans_rounded)
}

#Primary Growth Models 
#Gompertz
gompertz <- function(day, log10n0, log10nmax, mumax, lag){
  log10n <- log10n0 + (log10n0 < log10nmax)*(log10nmax - log10n0) * exp(-exp(mumax * exp(1) * (lag - day) / ((log10nmax - log10n0) * log(10)) + 1))
}


#Secondary model 

#Mumax
NewMu_general <- function(temp, b, temp_min){
    if(temp > temp_min){
      sqrt_mu <- b*(temp - temp_min)
      ans <- sqrt_mu^2
      return(ans)
    } else {
    ans <- 0 
    return(ans)
  }
}

NewMu_s12_0141 <- function(temp, b, temp_min){
  if(temp < temp_min){
    sqrt_mu <- b*(temp - temp_min)
    ans <- sqrt_mu^2
    return(ans)
  } else {
    ans <- 0 
    return(ans)
  }
}

NewMu <- function(temp, b, temp_min, isolate){
  if(isolate != "s12-0141"){
    return(NewMu_general(temp, b, temp_min))
  } else {
    return(NewMu_s12_0141(temp, b, temp_min))
  }
}

lag_tmin <- mumax_secon$temp_min_mu[mumax_secon$isolate == "s12-0116"]

NewLag <- function(temp, lag) {
  ans <- lag*((6 - lag_tmin)/(temp - lag_tmin))^2
  return(ans)
}

## ---------------------------Simulation Overview-------------------------------

n_lots <- 1
n_packages <- 100
n_sim <- n_packages*n_lots
n_day <- 10

## ---------------Assigning Model Inputs----------------------------------------

#Assigning initial count, number of strains, STs and temperature to each package
packages_sim <- data.frame(matrix(nrow = n_sim, ncol = 4))
colnames(packages_sim) <- c("package", "initial_count_apc", "storage_temp", "tot_strains")

packages_sim$package <- 1:n_packages
packages_sim$initial_count_apc <- initial_count_dist(n_sim)
packages_sim$storage_temp <- 4 #temp_dist(n_sim)
packages_sim$tot_strains <- pop_dist(n_sim)

#Expanding the dataframe to store observations for each strain and day of observation
packages_sim_exp <- packages_sim %>%
  group_by(package) %>%
  slice(rep(1, times = tot_strains)) %>%
  mutate(isolate = NA, initial_count = NA, lag = NA, mumax = NA, nmax = NA) 

#Assigning isolates to each package, with replacement
for(i in 1:n_packages){
  index <- which(packages_sim_exp$package == i)
  strain_per_package <- unique(packages_sim_exp$tot_strains[index])
  packages_sim_exp$isolate[index] <- st_dist(strain_per_package) 
  packages_sim_exp$initial_count[index] <- rep(log10((10^(packages_sim$initial_count_apc[i]))/strain_per_package), times = strain_per_package)
}

#Assign lag to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == strain_param$isolate & 6 == strain_param$temp)
  strain_lag <- NewLag(packages_sim_exp$storage_temp[i], strain_param$lag[strain_index])
  packages_sim_exp$lag[i] <- strain_lag
}

#Assign mumax to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == mumax_secon$isolate)
  strain_mu <- NewMu(packages_sim_exp$storage_temp[i], mumax_secon$b_mu[strain_index], mumax_secon$temp_min_mu[strain_index], packages_sim_exp$isolate[i])
  packages_sim_exp$mumax[i] <- strain_mu
}

#Assign nmax to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == strain_param$isolate & strain_param$temp == 6)
  strain_nmax <- strain_param$nmax[strain_index]
  packages_sim_exp$nmax[i] <- strain_nmax
}

#Calculate the total nmax of a package 
packages_sim_exp$final_nmax <- vector(mode = "integer", length = nrow(packages_sim_exp))

for(i in 1:n_packages){
  index <- which(packages_sim_exp$package == i)
  packages_nmax <- packages_sim_exp$nmax[index]
  final_nmax <- log10(sum(10^packages_nmax))
  packages_sim_exp$final_nmax[index] <- final_nmax
}

## ----------------Assessing Growth Over Shelf Life-----------------------------

#Creating a dataframe, to assess microbial growth over shelf life 
shelf_life_sim <- packages_sim_exp %>%
  group_by(package, isolate) %>%
  slice(rep(1, times = 10)) %>%
  mutate(day = 1:10)
  
shelf_life_sim <- shelf_life_sim %>%
  mutate(count = gompertz(day, initial_count, nmax, mumax, lag))

## --------------------Plotting & Analyzing Outcome-----------------------------

count_summary_by_day <- shelf_life_sim %>%
  group_by(day, package) %>%
  mutate(total_count = log10(sum(10^(count)))) %>%
  distinct(package, initial_count_apc, storage_temp, day, total_count) %>%
  mutate(growth_model = "mm-strain",
         loc = "china",
         sampling_code = "overall")

colnames(count_summary_by_day)[2] <- "initial_count"
colnames(count_summary_by_day)[5] <- "count"

## --------------------Exporting data-------------------------------------------

file_name <- paste0("outputs/count_by_day_mm_strain_china.csv")

write.table(count_summary_by_day, file_name, sep = ",", 
            col.names = !file.exists(file_name), append = T, row.names = FALSE)

file_name <- paste0("outputs/overall_count.csv")

write.table(count_summary_by_day, file_name, sep = ",", 
            col.names = !file.exists(file_name), append = T, row.names = FALSE)

