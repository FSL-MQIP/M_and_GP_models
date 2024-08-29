##  ----------------------------Title-------------------------------------------
#   Getting the initial conc for the U.S-based supply chain

##  --------------------------Description---------------------------------------
#   Project: Bacterial growth modeling (spinach)

#   Script description: See title; 

##  --------------------------Packages------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(minpack.lm)

##  --------------------------Data----------------------------------------------
#Read in data 
apc_data <- read.csv("data/raw/apc_strain_averaged_wrangled_01_openrefine_deathphase.csv", 
                     header = TRUE)

##  --------------------------Data Wrangling------------------------------------

apc_data <- apc_data %>%
  mutate(lot = NA)

for(i in 1:nrow(apc_data)){
  if(apc_data$biorep[i] == 2 & 
     apc_data$isolate[i] %in% c("S12-0116", "S12-0132", "S12-0141") &
     apc_data$temperature[i] == 6){
    
    apc_data$lot[i] <- 1
    
  } else if(apc_data$biorep[i] == 3 & 
            apc_data$isolate[i] %in% c("S12-0116", "S12-0132", "S12-0141") &
            apc_data$temperature[i] == 6){
    
    apc_data$lot[i] <- 2
    
  } else if(apc_data$biorep[i] == 4 & 
            apc_data$isolate[i] %in% c("S12-0116", "S12-0132", "S12-0141") &
            apc_data$temperature[i] == 6){
    
    apc_data$lot[i] <- 3
    
  } else if(apc_data$biorep[i] == 1 & 
            apc_data$isolate[i] %in% c("S12-0166", "S12-0180", "S12-0184") & 
            apc_data$temperature[i] == 6){
    
    apc_data$lot[i] <- 4
    
  } else if(apc_data$biorep[i] == 2 & 
            apc_data$isolate[i] %in% c("S12-0166", "S12-0180", "S12-0184") & 
            apc_data$temperature[i] == 6){
    
    apc_data$lot[i] <- 5
    
  }else if(apc_data$biorep[i] == 3 & 
           apc_data$isolate[i] %in% c("S12-0166", "S12-0180", "S12-0184") & 
           apc_data$temperature[i] == 6){
    
    apc_data$lot[i] <- 6
    
  } else if(apc_data$biorep[i] == 1 & 
           apc_data$isolate[i] %in% c("S12-0116", "S12-0132", "S12-0141") &
           apc_data$temperature[i] == 10){
    
    apc_data$lot[i] <- 7
    
  } else if(apc_data$biorep[i] == 1 & 
            apc_data$isolate[i] %in% c("S12-0166", "S12-0180", "S12-0184") &
            apc_data$temperature[i] == 10){
    
    apc_data$lot[i] <- 8
  }}
  

apc_data_2 <- apc_data %>%
  filter((isolate %in% c("S12-0166", "S12-0180", "S12-0184") & time == "36") |
           (isolate %in% c("S12-0116", "S12-0132", "S12-0141")))

#Removing the BHIrif observation that did not have a corresponding 
#APC observation
apc_strain_data_bhi_adj <- apc_data_2 %>%
  filter(!(isolate == "S12-0184" & batch == "B7" & day == "2" & temperature == "6")) %>%
  group_by(temperature, batch, day, isolate) %>%
  mutate(average_wrangled_conc_adj = average_wrangled_conc[media == "BHI"] - average_wrangled_conc[media == "BHIrif"])

#Obtaining the min-max day 0 concentration for each lot
d0_bhi_conc <- apc_strain_data_bhi_adj %>% 
  filter(media == "BHI") %>%
  filter(day == 0 )

lot_min_max <- d0_bhi_conc %>%
  filter(!(average_wrangled_conc_adj <= 0)) %>% #Leaving out the observation where the day 0 concentration was 0 (S12-0132, lot 1)
  mutate(log_adj_conc = log10(average_wrangled_conc_adj)) %>% 
  group_by(lot) %>%
  summarize(min = min(log_adj_conc), max = max(log_adj_conc))

#Mean and range of day 0 concentration (geometric mean)
summary_lot <- d0_bhi_conc %>%
  filter(!(average_wrangled_conc_adj <= 0)) %>%
  group_by(lot) %>%
  mutate(log_adj_conc = log10(average_wrangled_conc_adj)) %>%
  ungroup() %>%
  summarize(min = min(log_adj_conc), max = max(log_adj_conc), mean = mean(log_adj_conc))


##---------------------------Export data----------------------------------------

write.csv(lot_min_max, paste0("output/apc_day0_conc_dist.csv"), row.names = FALSE)


