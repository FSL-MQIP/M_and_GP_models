## -----------------------------Title-------------------------------------------
# Wrangling Raw Data to fit ML models in R

## -----------------------------Description-------------------------------------
# Project: Bacterial growth modeling (spinach)

# Script description: ensuring all data frames (i.e., growth curve data,
# walmart data) have the same columns: lot, temperature, day, log_conc_geom, lot_d0_conc

## -----------------------------Packages----------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); 

## -----------------------------Reading in Data---------------------------------
growth_curve <- read.csv("data/raw/apc_strain_averaged_wrangled_01_openrefine_deathphase.csv", 
                         header = TRUE)


walmart_sampling_data <- read.csv("data/raw/concentration_by_day_modeling_data.csv",
                                  header = TRUE)

## -------------------Growth Curves: Wrangling----------------------------------

growth_curve_wrang <- growth_curve %>%
  mutate(lot = case_when(
    (isolate == "S12-0116" & biorep == "2" & temperature == 6) ~ "1",
    (isolate == "S12-0132" & biorep == "2" & temperature == 6) ~ "1",
    (isolate == "S12-0141" & biorep == "2" & temperature == 6) ~ "1",
    
    (isolate == "S12-0116" & biorep == "3" & temperature == 6) ~ "2",
    (isolate == "S12-0132" & biorep == "3" & temperature == 6) ~ "2",
    (isolate == "S12-0141" & biorep == "3" & temperature == 6) ~ "2",
    
    (isolate == "S12-0116" & biorep == "4" & temperature == 6) ~ "3",
    (isolate == "S12-0132" & biorep == "4" & temperature == 6) ~ "3",
    (isolate == "S12-0141" & biorep == "4" & temperature == 6) ~ "3",
    
    (isolate == "S12-0166" & biorep == "1" & temperature == 6) ~ "4",
    (isolate == "S12-0180" & biorep == "1" & temperature == 6) ~ "4",
    (isolate == "S12-0184" & biorep == "1" & temperature == 6) ~ "4",
    
    (isolate == "S12-0166" & biorep == "2" & temperature == 6) ~ "5",
    (isolate == "S12-0180" & biorep == "2" & temperature == 6) ~ "5",
    (isolate == "S12-0184" & biorep == "2" & temperature == 6) ~ "5",
    
    (isolate == "S12-0166" & biorep == "3" & temperature == 6) ~ "6",
    (isolate == "S12-0180" & biorep == "3" & temperature == 6) ~ "6",
    (isolate == "S12-0184" & biorep == "3" & temperature == 6) ~ "6",
    
    (isolate == "S12-0116" & biorep == "1" & temperature == 10) ~ "7",
    (isolate == "S12-0132" & biorep == "1" & temperature == 10) ~ "7",
    (isolate == "S12-0141" & biorep == "1" & temperature == 10) ~ "7",
    
    (isolate == "S12-0166" & biorep == "1" & temperature == 10) ~ "8",
    (isolate == "S12-0180" & biorep == "1" & temperature == 10) ~ "8",
    (isolate == "S12-0184" & biorep == "1" & temperature == 10) ~ "8",
  )) %>%
  mutate(tech_rep = case_when(
    (isolate == "S12-0116" & biorep == "2" & temperature == 6) ~ "1",
    (isolate == "S12-0132" & biorep == "2" & temperature == 6) ~ "2",
    (isolate == "S12-0141" & biorep == "2" & temperature == 6) ~ "3",
    
    (isolate == "S12-0116" & biorep == "3" & temperature == 6) ~ "1",
    (isolate == "S12-0132" & biorep == "3" & temperature == 6) ~ "2",
    (isolate == "S12-0141" & biorep == "3" & temperature == 6) ~ "3",
    
    (isolate == "S12-0116" & biorep == "4" & temperature == 6) ~ "1",
    (isolate == "S12-0132" & biorep == "4" & temperature == 6) ~ "2",
    (isolate == "S12-0141" & biorep == "4" & temperature == 6) ~ "3",
    
    (isolate == "S12-0166" & biorep == "1" & temperature == 6) ~ "1",
    (isolate == "S12-0180" & biorep == "1" & temperature == 6) ~ "2",
    (isolate == "S12-0184" & biorep == "1" & temperature == 6) ~ "3",
    
    (isolate == "S12-0166" & biorep == "2" & temperature == 6) ~ "1",
    (isolate == "S12-0180" & biorep == "2" & temperature == 6) ~ "2",
    (isolate == "S12-0184" & biorep == "2" & temperature == 6) ~ "3",
    
    (isolate == "S12-0166" & biorep == "3" & temperature == 6) ~ "1",
    (isolate == "S12-0180" & biorep == "3" & temperature == 6) ~ "2",
    (isolate == "S12-0184" & biorep == "3" & temperature == 6) ~ "3",
    
    (isolate == "S12-0116" & biorep == "1" & temperature == 10) ~ "1",
    (isolate == "S12-0132" & biorep == "1" & temperature == 10) ~ "2",
    (isolate == "S12-0141" & biorep == "1" & temperature == 10) ~ "3",
    
    (isolate == "S12-0166" & biorep == "1" & temperature == 10) ~ "1",
    (isolate == "S12-0180" & biorep == "1" & temperature == 10) ~ "2",
    (isolate == "S12-0184" & biorep == "1" & temperature == 10) ~ "3",
  ))

#Subtracting strain-specific counts from the APC counts

growth_curve_wrang_2 <- growth_curve_wrang %>%
  filter((isolate %in% c("S12-0166", "S12-0180", "S12-0184") & time == "36") |
         (isolate %in% c("S12-0116", "S12-0132", "S12-0141")))

growth_curve_wrang_3 <- growth_curve_wrang_2 %>%
  filter(!(isolate == "S12-0184" & batch == "B7" & day == "2" & temperature == "6")) %>%
  group_by(temperature, batch, day, isolate) %>%
  mutate(average_wrangled_conc_adj = average_wrangled_conc[media == "BHI"] - average_wrangled_conc[media == "BHIrif"])

growth_curve_wrang_4 <- growth_curve_wrang_3 %>%
  filter(media == "BHI") %>%
  select(c("batch", "temperature", "day", "isolate", "media", "time", "death_phase", 
           "lot", "tech_rep", "average_wrangled_conc_adj"))

colnames(growth_curve_wrang_4)[10] <- "average_wrangled_conc"

growth_curve_wrang_5 <- growth_curve_wrang_4 %>%
  filter(average_wrangled_conc > 0) %>%
  mutate(log_average_wrangled_conc = log10(average_wrangled_conc))

growth_curve_wrang_6 <- growth_curve_wrang_5 %>%
  group_by(lot, day) %>%
  mutate(conc_geom = mean(log_average_wrangled_conc)) %>%
  distinct(lot, day, temperature, conc_geom) %>%
  ungroup() %>%
  group_by(lot) %>%
  mutate(lot_d0_conc = conc_geom[day == 0]) %>%
  ungroup()

min_max_d0_us <- growth_curve_wrang_6 %>%
  filter(day == 0) %>%
  summarize(min = round(min(conc_geom), 2), 
            max = round(max(conc_geom), 2))

## -------------------Walmart Sampling: Wrangling-------------------------------

walmart_dat <- walmart_sampling_data %>%
  filter(grepl("D", stage)) %>%
  filter(test == "APC") %>%
  dplyr::select(-c("season")) %>%
  mutate(temperature = 4) %>%
  filter(!(sampling %in% c("S4", "S7", "S9")))

walmart_dat$day <- as.numeric(substr(walmart_dat$stage, 2, as.character(length(walmart_dat$stage))))
  
walmart_dat <- walmart_dat %>%
  dplyr::select(-c("stage", "test")) %>%
  group_by(sampling) %>%
  mutate(lot_d0_conc = conc_geom[day == 0])

colnames(walmart_dat)[1:2] <- c("lot", "log_conc_geom")

min_max_d0_china <- walmart_dat %>%
  ungroup() %>%
  filter(day == 0) %>%
  summarize(min = round(min(log_conc_geom), 2), 
            max = round(max(log_conc_geom), 2))

## -----------------------------------------------------------------------------
#Save wrangled data to R Project folder 

#Push the wrangled data back to the R project 
write.csv(growth_curve_wrang_6, paste("data/wrangled/growth_curve_biorep.csv", sep = ""), row.names = FALSE)
write.csv(walmart_dat, paste("data/wrangled/walmart_biorep_01.csv", sep = ""), row.names = FALSE)

# End of saving data into R Project folder 

