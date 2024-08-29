##  ----------------------------Title-------------------------------------------
#   Fitting the Ratkowsky square root model to the strain growth data (averaged growth parameters) 

##  --------------------------Description---------------------------------------
#   Project: Bacterial growth modeling (spinach)

#   Script description: Fitting a secondary model, to assess change in mu, based on primary growth parameters

##  --------------------------Packages------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(ggplot2)

##  --------------------------Data----------------------------------------------
#Read in data 
primary_param <- read.csv("data/raw/primary_model_parameters_apc_averaged.csv", 
                          header = TRUE)

## ------------------Defining Secondary Model-----------------------------------

#Ratkowsky
ratkowsky <- function(temp, b, temp_min){
  sqrt_new_mu <- b*(temp - temp_min)
  new_mu <- sqrt_new_mu^2
  return(new_mu)
}

#End of defining secondary model

## ------------------Data Wrangling---------------------------------------------

#Making a dataframe to store the fitted model parameter estimates
no_of_models <- 9
secondary_param <- matrix(nrow = no_of_models, ncol = 4)
secondary_param <- data.frame(secondary_param)
colnames(secondary_param) <- c("pop", "b_mu", "temp_min_mu", "val_data_lot")

primary_param <- primary_param %>%
  mutate(sqrt_mu = sqrt(mumax))

#End of data wrangling

## ------------------Model Fitting----------------------------------------------
#Lot 1
mod_mu_1 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 1))
summary(mod_mu_1)

secondary_param[1, 1] <- "apc"
secondary_param[1, 2] <- mod_mu_1$coefficients[2]
secondary_param[1, 3] <- (-mod_mu_1$coefficients[1])/mod_mu_1$coefficients[2]
secondary_param[1, 4] <- 1

#Lot 2
mod_mu_2 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 2))
summary(mod_mu_2)

secondary_param[2, 1] <- "apc"
secondary_param[2, 2] <- mod_mu_2$coefficients[2]
secondary_param[2, 3] <- (-mod_mu_2$coefficients[1])/mod_mu_2$coefficients[2]
secondary_param[2, 4] <- 2

#Lot 3
mod_mu_3 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 3))
summary(mod_mu_3)

secondary_param[3, 1] <- "apc"
secondary_param[3, 2] <- mod_mu_3$coefficients[2]
secondary_param[3, 3] <- (-mod_mu_3$coefficients[1])/mod_mu_3$coefficients[2]
secondary_param[3, 4] <- 3

#Lot 4
mod_mu_4 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 4))
summary(mod_mu_4)

secondary_param[4, 1] <- "apc"
secondary_param[4, 2] <- mod_mu_4$coefficients[2]
secondary_param[4, 3] <- (-mod_mu_4$coefficients[1])/mod_mu_4$coefficients[2]
secondary_param[4, 4] <- 4

#Lot 5
mod_mu_5 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 5))
summary(mod_mu_5)

secondary_param[5, 1] <- "apc"
secondary_param[5, 2] <- mod_mu_5$coefficients[2]
secondary_param[5, 3] <- (-mod_mu_5$coefficients[1])/mod_mu_5$coefficients[2]
secondary_param[5, 4] <- 5

#Lot 6
mod_mu_6 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 6))
summary(mod_mu_6)

secondary_param[6, 1] <- "apc"
secondary_param[6, 2] <- mod_mu_6$coefficients[2]
secondary_param[6, 3] <- (-mod_mu_6$coefficients[1])/mod_mu_6$coefficients[2]
secondary_param[6, 4] <- 6

#Lot 7
mod_mu_7 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 7))
summary(mod_mu_7)

secondary_param[7, 1] <- "apc"
secondary_param[7, 2] <- mod_mu_7$coefficients[2]
secondary_param[7, 3] <- (-mod_mu_7$coefficients[1])/mod_mu_7$coefficients[2]
secondary_param[7, 4] <- 7

#Lot 8
mod_mu_8 <- lm(sqrt_mu ~ temp, 
               data = filter(primary_param, lot != 8))
summary(mod_mu_8)

secondary_param[8, 1] <- "apc"
secondary_param[8, 2] <- mod_mu_8$coefficients[2]
secondary_param[8, 3] <- (-mod_mu_8$coefficients[1])/mod_mu_8$coefficients[2]
secondary_param[8, 4] <- 8

#Overall
mod_mu_overall <- lm(sqrt_mu ~ temp, 
                     data = primary_param)
summary(mod_mu_overall)

secondary_param[9, 1] <- "apc"
secondary_param[9, 2] <- mod_mu_overall$coefficients[2]
secondary_param[9, 3] <- (-mod_mu_overall$coefficients[1])/mod_mu_overall$coefficients[2]
secondary_param[9, 4] <- "overall"

## ------------------------Predicting mu at new temperatures--------------------

apc <- rep(secondary_param$pop, each = 28)
temperature <- rep(3:30, times = length(secondary_param$pop))
predicted_mu <- vector(mode = "logical", length = length(apc))
lot <- rep(secondary_param$val_data_lot, each = 28)

pred_new_mu <- bind_cols(apc, temperature, predicted_mu, lot)
colnames(pred_new_mu) <- c("apc", "temp", "predicted_mu", "left_out_lot")

for(i in 1:nrow(pred_new_mu)){
  index <- which(pred_new_mu$left_out_lot[i] == secondary_param$val_data_lot)
  pred_new_mu$predicted_mu[i] <- ratkowsky(pred_new_mu$temp[i], secondary_param$b_mu[index], secondary_param$temp_min_mu[index])
}

#Lot1
l1 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 1)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 1)) + 
  labs(title = "Lot left out: 1")

#Lot2
l2 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 2)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 2)) + 
  labs(title = "Lot left out: 2")

#Lot3
l3 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 3)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 3)) + 
  labs(title = "Lot left out: 3")

#Lot4
l4 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 4)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 4)) + 
  labs(title = "Lot left out: 4")

#Lot5
l5 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 5)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 5)) + 
  labs(title = "Lot left out: 5")

#Lot6
l6 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 6)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 6)) + 
  labs(title = "Lot left out: 6")

#Lot7
l7 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 7)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 7)) + 
  labs(title = "Lot left out: 7")

#Lot8
l8 <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 8)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot != 8)) + 
  labs(title = "Lot left out: 8")

#overall
overall_plot <- ggplot() + geom_line(mapping = aes(x = temp, y = predicted_mu), data = filter(pred_new_mu, lot == 8)) + 
  geom_point(mapping = aes(x = temp, y = mumax), data = filter(primary_param, lot == "overall")) + 
  labs(title = "Lot left out: none")

## -------------------------------Push data back onto server--------------------
#Push the data back to the R project 

#Parameters
write.csv(secondary_param, paste0("output/secondary_model_mumax_apc_averaged.csv"), row.names = FALSE)

# End of saving data into R Project folder 









