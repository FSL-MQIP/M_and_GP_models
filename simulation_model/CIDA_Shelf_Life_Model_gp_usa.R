## ---------------------------Title---------------------------------------------
# CIDA Shelf Life Model

## ---------------------------Description---------------------------------------
# Project: CIDA Spinach & Bacterial growth modeling (spinach)

# Description: Monte Carlo simulation of bacterial growth on spinach

## ---------------------------Packages------------------------------------------

#   Loading packages
library(dplyr); library(MASS); library(VGAM); library(nlme)

## ---------------------------Seed----------------------------------------------

#Set seed
set.seed(1)

## ---------------------------Inputs--------------------------------------------
#Initial_count
min_max <- read.csv("model_inputs/apc_day0_conc_dist.csv", header = TRUE)

#Load in data

mod <- readRDS("model_inputs/ml_models/gp_l1.RDS")

#Initial Count Distribution: Fitting a uniform distribution to the data

initial_count_dist <- function(n){
  ans <- runif(n, 
               min = min_max$min[min_max$lot == 1], 
               max = min_max$max[min_max$lot == 1])
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

## ---------------------------Simulation Overview-------------------------------

n_lots <- 1
n_packages <- 100
n_sim <- n_packages*n_lots
n_day <- 30

## ---------------Assigning Model Inputs----------------------------------------

#Assigning initial count, and temperature to each package
packages_sim <- data.frame(matrix(nrow = n_sim, ncol = 3))
colnames(packages_sim) <- c("package", "initial_count", "storage_temp")

packages_sim$package <- 1:n_packages
packages_sim$initial_count<- initial_count_dist(n_sim)
packages_sim$storage_temp <- 6 #temp_dist(n_sim)

#Expanding the dataframe to store observations 
packages_sim_exp <- packages_sim 

## ----------------Assessing Growth Over Shelf Life-----------------------------

#Creating a dataframe, to assess microbial growth over shelf life 
shelf_life_sim <- packages_sim_exp %>%
  group_by(package) %>%
  dplyr::slice(rep(1, times = 30)) %>%
  mutate(day = 1:30) %>%
  mutate(day2 = I(day^2), day.log_d0_conc = day*initial_count, 
         day2.log_d0_conc = day2*initial_count) %>%
  mutate(count = NA)

colnames(shelf_life_sim)[1] <- "lot"
colnames(shelf_life_sim)[2] <- "lot_d0_conc"
colnames(shelf_life_sim)[3] <- "temperature"
colnames(shelf_life_sim)[6] <- "day:lot_d0_conc"
colnames(shelf_life_sim)[7] <- "lot_d0_conc:day2"

summary(mod)

shelf_life_sim$count <- predict(mod, shelf_life_sim[, c(4, 2, 5, 3, 6, 7)], level = 0)

colnames(shelf_life_sim)[1] <- "package"
colnames(shelf_life_sim)[2] <- "initial_count"
colnames(shelf_life_sim)[3] <- "storage_temp"

## --------------------Plotting & Analyzing Outcome-----------------------------

count_summary_by_day <- shelf_life_sim %>%
  dplyr::select(-c("day2", "day:lot_d0_conc", "lot_d0_conc:day2")) %>%
  mutate(growth_model = "gp-apc",
         loc = "us",
         sampling_code = 1)
  
## --------------------Exporting data-------------------------------------------

file_name <- paste0("outputs/count_by_day_gp_apc_us.csv")

write.table(count_summary_by_day, file_name, sep = ",", 
            col.names = !file.exists(file_name), append = T, row.names = FALSE)

file_name <- paste0("outputs/overall_count.csv")

write.table(count_summary_by_day, file_name, sep = ",", 
            col.names = !file.exists(file_name), append = T, row.names = FALSE)

