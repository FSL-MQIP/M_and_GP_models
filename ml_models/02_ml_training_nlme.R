## -----------------------------Title-------------------------------------------
# Training GP Models

## -----------------------------Description-------------------------------------
# Project: Bacterial growth modeling (spinach)

# Script description: Training GP models to data from the US and China based supply chains 

## -----------------------------Packages----------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(nlme); library(gpboost); 
library(iml); library(Metrics)

## -----------------------------Reading in Data---------------------------------
growth_curves <- read.csv("data/wrangled/growth_curve_biorep.csv", header = TRUE)

growth_curves$lot <- as.factor(growth_curves$lot)
growth_curves$temperature <- as.numeric(growth_curves$temperature)
growth_curves$day <- as.numeric(growth_curves$day)
growth_curves$conc_geom <- as.numeric(growth_curves$conc_geom)
growth_curves$lot_d0_conc <- as.numeric(growth_curves$lot_d0_conc)

walmart <- read.csv("data/wrangled/walmart_biorep.csv", header = TRUE)

walmart$lot <- as.factor(walmart$lot)
walmart$temperature <- as.numeric(walmart$temperature)
walmart$day <- as.numeric(walmart$day)
walmart$log_conc_geom <- as.numeric(walmart$log_conc_geom)
walmart$lot_d0_conc <- as.numeric(walmart$lot_d0_conc)

## -----------------------------Wrangling---------------------------------------

gc_train_data <- growth_curves %>%
  filter(day != 0) %>%
  mutate(day2 = day^2) %>%
  mutate(day3 = day^3)

gc_train_data <- gc_train_data[, c(4, 2, 6, 3, 5, 7, 1)]

w_train_data <- walmart %>%
  filter(day != 0) %>%
  mutate(day2 = I(day^2))

w_train_data <- w_train_data[, c(2, 4, 6, 5, 1)]

## -----------------------------Gaussian Process--------------------------------

#US Supply Chain
#Overall
gc_mod_o1 <- lme(conc_geom ~ day + day2  +  temperature + lot_d0_conc, 
                    data = gc_train_data,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_o1)

gc_mod_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc  +  temperature + lot_d0_conc, 
                    data = gc_train_data,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_o2)
anova(gc_mod_o1, gc_mod_o2) 

rm(gc_mod_o1) #Use gc_mod_o2

#W/o l1
gc_train_data_l1 <- filter(gc_train_data, lot != "1")

gc_mod_l1_o1 <- lme(conc_geom ~ day + day2  +  temperature + lot_d0_conc, 
                 data = gc_train_data_l1,
                 correlation = corExp(form = ~day|lot), random = ~1|lot,
                 method = "ML")
summary(gc_mod_l1_o1)

gc_mod_l1_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc  +  temperature + lot_d0_conc, 
                    data = gc_train_data_l1,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l1_o2)
anova(gc_mod_l1_o1, gc_mod_l1_o2) 

rm(gc_mod_l1_o1) #Use gc_mod_l1_o2

#W/o l2
gc_train_data_l2 <- filter(gc_train_data, lot != "2")

gc_mod_l2_o1 <- lme(conc_geom ~ day + day2 + temperature + lot_d0_conc, 
                 data = gc_train_data_l2,
                 correlation = corExp(form = ~day|lot), random = ~1|lot,
                 method = "ML")
summary(gc_mod_l2_o1)

gc_mod_l2_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + temperature + lot_d0_conc, 
                    data = gc_train_data_l2,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l2_o2)
anova(gc_mod_l2_o1, gc_mod_l2_o2)

rm(gc_mod_l2_o1) #Use gc_mod_l2_o2

#W/o l3
gc_train_data_l3 <- filter(gc_train_data, lot != "3")

gc_mod_l3_o1 <- lme(conc_geom ~ day + day2 + temperature + lot_d0_conc, 
                    data = gc_train_data_l3,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l3_o1)

gc_mod_l3_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + temperature + lot_d0_conc, 
                    data = gc_train_data_l3,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l3_o2)
anova(gc_mod_l3_o1, gc_mod_l3_o2)

rm(gc_mod_l3_o1) #Use gc_mod_l3_o2

#W/o l4
gc_train_data_l4 <- filter(gc_train_data, lot != "4")

gc_mod_l4_o1 <- lme(conc_geom ~ day + day2 + temperature + lot_d0_conc, 
                    data = gc_train_data_l4,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")

gc_mod_l4_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + temperature + lot_d0_conc, 
                    data = gc_train_data_l4,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l4_o2)
anova(gc_mod_l4_o1, gc_mod_l4_o2)

rm(gc_mod_l4_o1) #Use gc_mod_l4_o2

#W/o l5
gc_train_data_l5 <- filter(gc_train_data, lot != "5")

gc_mod_l5_o1 <- lme(conc_geom ~ day + day2 + temperature + lot_d0_conc, 
                    data = gc_train_data_l5,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")

gc_mod_l5_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + temperature + lot_d0_conc, 
                    data = gc_train_data_l5,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l5_o2)
anova(gc_mod_l5_o1, gc_mod_l5_o2)

rm(gc_mod_l5_o1) #Use gc_mod_l5_o1

#W/o l6
gc_train_data_l6 <- filter(gc_train_data, lot != "6")

gc_mod_l6_o1 <- lme(conc_geom ~ day + day2 + temperature + lot_d0_conc, 
                    data = gc_train_data_l6,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l6_o1)

gc_mod_l6_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + temperature + lot_d0_conc, 
                    data = gc_train_data_l6,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l6_o2)
anova(gc_mod_l6_o1, gc_mod_l6_o2)

rm(gc_mod_l6_o1)  #Use gc_mod_l6_o2

#W/o l7
gc_train_data_l7 <- filter(gc_train_data, lot != "7")

gc_mod_l7_o1 <- lme(conc_geom ~ day + day2 + temperature + lot_d0_conc, 
                    data = gc_train_data_l7,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l7_o1)

gc_mod_l7_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + temperature + lot_d0_conc, 
                    data = gc_train_data_l7,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l7_o2)
anova(gc_mod_l7_o1, gc_mod_l7_o2)

rm(gc_mod_l7_o1) #Use gc_mod_l7_o2

#W/o l8
gc_train_data_l8 <- filter(gc_train_data, lot != "8")

gc_mod_l8_o1 <- lme(conc_geom ~ day + day2 + temperature + lot_d0_conc, 
                    data = gc_train_data_l8,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l8_o1)

gc_mod_l8_o2 <- lme(conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + temperature + lot_d0_conc, 
                    data = gc_train_data_l8,
                    correlation = corExp(form = ~day|lot), random = ~1|lot,
                    method = "ML")
summary(gc_mod_l8_o2)
anova(gc_mod_l8_o1, gc_mod_l8_o2)

rm(gc_mod_l8_o1) #Use gc_mod_l8_o2

#China supply chain
#Overall
w_mod_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                   data = w_train_data,
                   correlation = corExp(form = ~day|lot), random = ~1|lot,
                   method = "ML")
summary(w_mod_o1)

w_mod_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                   data = w_train_data,
                   correlation = corExp(form = ~day|lot), random = ~1|lot,
                   method = "ML")
summary(w_mod_o2)
anova(w_mod_o1, w_mod_o2)

rm(w_mod_o2) #Use w_mod_s1_o1

#S1 
w_train_data_s1 <- filter(w_train_data, lot != "S1")

w_mod_s1_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_s1,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_s1_o1)

w_mod_s1_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_s1,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_s1_o2)
anova(w_mod_s1_o1, w_mod_s1_o2)

rm(w_mod_s1_o2) #Use w_mod_s1_o1

#S2
w_train_data_s2 <- filter(w_train_data, lot != "S2")

w_mod_s2_o1 <- lme(log_conc_geom ~ day  + day2 + lot_d0_conc, 
                data = w_train_data_s2,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_s2_o1)

w_mod_s2_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_s2,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_s2_o2)
anova(w_mod_s2_o1, w_mod_s2_o2)

rm(w_mod_s2_o2) #Use w_mod_s2_o1

#S3
w_train_data_s3 <- filter(w_train_data, lot != "S3")

w_mod_s3_o1 <- lme(log_conc_geom ~ day + lot_d0_conc, 
                data = w_train_data_s3,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_s3_o1)

w_mod_s3_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_s3,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_s3_o2)
anova(w_mod_s3_o1, w_mod_s3_o2)

rm(w_mod_s3_o2) #Use w_mod_s3_o1

#F1
w_train_data_f1 <- filter(w_train_data, lot != "F1")

w_mod_f1_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_f1,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_f1_o1)

w_mod_f1_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_f1,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_f1_o2)
anova(w_mod_f1_o1, w_mod_f1_o2)

rm(w_mod_f1_o2) #Use w_mod_f1_o1

#F2
w_train_data_f2 <- filter(w_train_data, lot != "F2")

w_mod_f2_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_f2,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_f2_o1)

w_mod_f2_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                   data = w_train_data_f2,
                   correlation = corExp(form = ~day|lot), random = ~1|lot,
                   method = "ML")
summary(w_mod_f2_o2)
anova(w_mod_f2_o1, w_mod_f2_o2)

rm(w_mod_f2_o2) #Use w_mod_f2_o1

#F3
w_train_data_f3 <- filter(w_train_data, lot != "F3")

w_mod_f3_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_f3,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_f3_o1)

w_mod_f3_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_f3,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_f3_o2)
anova(w_mod_f3_o1, w_mod_f3_o2)

rm(w_mod_f3_o2) #Use w_mod_f3_o1

#F4
w_train_data_f4 <- filter(w_train_data, lot != "F4")

w_mod_f4_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_f4,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_f4_o1)

w_mod_f4_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_f4,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_f4_o2)

anova(w_mod_f4_o1, w_mod_f4_o2)

rm(w_mod_f4_o2) #Use w_mod_f4_o1

#R1
w_train_data_r1 <- filter(w_train_data, lot != "R1")

w_mod_r1_o1 <- lme(log_conc_geom ~ day  + day2 + lot_d0_conc, 
                data = w_train_data_r1,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r1_o1)

w_mod_r1_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_r1,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r1_o2)
anova(w_mod_r1_o1, w_mod_r1_o2)

rm(w_mod_r1_o2) #Use w_mod_r1_o1

#R2
w_train_data_r2 <- filter(w_train_data, lot != "R2")

w_mod_r2_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_r2,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")

w_mod_r2_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                   data = w_train_data_r2,
                   correlation = corExp(form = ~day|lot), random = ~1|lot,
                   method = "ML")
summary(w_mod_r2_o2)
anova(w_mod_r2_o1, w_mod_r2_o2)

rm(w_mod_r2_o2) #Use w_mod_r2_o1

#R3
w_train_data_r3 <- filter(w_train_data, lot != "R3")

w_mod_r3_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_r3,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r3_o1)

w_mod_r3_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_r3,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r3_o2)
anova(w_mod_r3_o1, w_mod_r3_o2)

rm(w_mod_r3_o2) #Use w_mod_r3_o1

#R4
w_train_data_r4 <- filter(w_train_data, lot != "R4")

w_mod_r4_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_r4,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r4_o1)

w_mod_r4_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_r4,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r4_o2)
anova(w_mod_r4_o1, w_mod_r4_o2)

rm(w_mod_r4_o2) #Use w_mod_r4_o1

#R5
w_train_data_r5 <- filter(w_train_data, lot != "R5")

w_mod_r5_o1 <- lme(log_conc_geom ~ day  + day2 + lot_d0_conc, 
                data = w_train_data_r5,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r5_o1)

w_mod_r5_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_r5,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r5_o2)
anova(w_mod_r5_o1, w_mod_r5_o2)

rm(w_mod_r5_o2) #Use w_mod_r5_o1

#R6
w_train_data_r6 <- filter(w_train_data, lot != "R6")

w_mod_r6_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_r6,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r6_o1)

w_mod_r6_o2 <- lme(log_conc_geom ~ day*lot_d0_conc  + day2*lot_d0_conc + lot_d0_conc, 
                   data = w_train_data_r6,
                   correlation = corExp(form = ~day|lot), random = ~1|lot,
                   method = "ML")
summary(w_mod_r6_o2)

anova(w_mod_r6_o1, w_mod_r6_o2)

rm(w_mod_r6_o2) #Use w_mod_r6_o1

#R7
w_train_data_r7 <- filter(w_train_data, lot != "R7")

w_mod_r7_o1 <- lme(log_conc_geom ~ day + day2 + lot_d0_conc, 
                data = w_train_data_r7,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r7_o1)

w_mod_r7_o2 <- lme(log_conc_geom ~ day*lot_d0_conc + day2*lot_d0_conc + lot_d0_conc, 
                data = w_train_data_r7,
                correlation = corExp(form = ~day|lot), random = ~1|lot,
                method = "ML")
summary(w_mod_r7_o2)
anova(w_mod_r7_o1, w_mod_r7_o2)

rm(w_mod_r7_o1) #Use w_mod_r7_o2

## -----------------------------------------------------------------------------
#Save wrangled data to R Project folder 

saveRDS(gc_mod_o2, paste("output/fit_ml/gp_us_overall.RDS", sep = ""))
saveRDS(gc_mod_l1_o2, paste("output/fit_ml/gp_l1.RDS", sep = ""))
saveRDS(gc_mod_l2_o2, paste("output/fit_ml/gp_l2.RDS", sep = ""))
saveRDS(gc_mod_l3_o2, paste("output/fit_ml/gp_l3.RDS", sep = ""))
saveRDS(gc_mod_l4_o2, paste("output/fit_ml/gp_l4.RDS", sep = ""))
saveRDS(gc_mod_l5_o2, paste("output/fit_ml/gp_l5.RDS", sep = ""))
saveRDS(gc_mod_l6_o2, paste("output/fit_ml/gp_l6.RDS", sep = ""))
saveRDS(gc_mod_l7_o2, paste("output/fit_ml/gp_l7.RDS", sep = ""))
saveRDS(gc_mod_l8_o2, paste("output/fit_ml/gp_l8.RDS", sep = ""))

saveRDS(w_mod_o1, paste("output/fit_ml/gp_china_overall.RDS", sep = ""))
saveRDS(w_mod_s1_o1, paste("output/fit_ml/gp_s1.RDS", sep = ""))
saveRDS(w_mod_s2_o1, paste("output/fit_ml/gp_s2.RDS", sep = ""))
saveRDS(w_mod_s3_o1, paste("output/fit_ml/gp_s3.RDS", sep = ""))
saveRDS(w_mod_f1_o1, paste("output/fit_ml/gp_f1.RDS", sep = ""))
saveRDS(w_mod_f2_o1, paste("output/fit_ml/gp_f2.RDS", sep = ""))
saveRDS(w_mod_f3_o1, paste("output/fit_ml/gp_f3.RDS", sep = ""))
saveRDS(w_mod_f4_o1, paste("output/fit_ml/gp_f4.RDS", sep = ""))
saveRDS(w_mod_r1_o1, paste("output/fit_ml/gp_r1.RDS", sep = ""))
saveRDS(w_mod_r2_o1, paste("output/fit_ml/gp_r2.RDS", sep = ""))
saveRDS(w_mod_r3_o1, paste("output/fit_ml/gp_r3.RDS", sep = ""))
saveRDS(w_mod_r4_o1, paste("output/fit_ml/gp_r4.RDS", sep = ""))
saveRDS(w_mod_r5_o1, paste("output/fit_ml/gp_r5.RDS", sep = ""))
saveRDS(w_mod_r6_o1, paste("output/fit_ml/gp_r6.RDS", sep = ""))
saveRDS(w_mod_r7_o2, paste("output/fit_ml/gp_r7.RDS", sep = ""))

# End of saving data into R Project folder 

