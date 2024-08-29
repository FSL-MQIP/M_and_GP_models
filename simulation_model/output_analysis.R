## ------------------------------Title------------------------------------------
# Analyzing the modeling data

## ------------------------------Description------------------------------------
# Project: Bacterial growth modeling (spinach)

# Script description: See title

## ------------------------------Packages---------------------------------------

library(dplyr); library(lmerTest); library(ggplot2); library(Metrics)
library(reshape2); library(tidyverse); library(RColorBrewer)
library(ggpubr)

## ------------------------------Read in raw data-------------------------------

data_geom_china <- read.csv("model_inputs/concentration_by_day_modeling_data.csv")
data_geom_us <- read.csv("model_inputs/growth_curve_biorep.csv")

mm_strain_us <- read.csv("outputs/count_by_day_mm_strain_us.csv", header = TRUE)
mm_strain_china <- read.csv("outputs/count_by_day_mm_strain_china.csv", header = TRUE)

mm_apc_us <- read.csv("outputs/count_by_day_mm_apc_us.csv", header = TRUE)

gp_apc_us <- read.csv("outputs/count_by_day_gp_apc_us.csv", header = TRUE)
gp_apc_china <- read.csv("outputs/count_by_day_gp_apc_china.csv", header = TRUE)

overall <- read.csv("outputs/overall_count.csv", header = TRUE)

## ------------------------------Wrangling - Exp--------------------------------

data_obs_china <- data_geom_china %>%
  filter(test == "APC") %>%
  filter(!(stage %in% c("H", "T", "D0"))) %>%
  ungroup() %>%
  dplyr::select(sampling, stage, conc_geom) %>%
  mutate(source = "obs")

colnames(data_obs_china) <- c("sampling_code", "day", "obs_conc", "source")

data_obs_china$day <- as.numeric(gsub("D", "", data_obs_china$day))
data_obs_china$sampday <- paste0(data_obs_china$sampling_code, "_", data_obs_china$day)

rm(data_geom_china)

data_obs_us <- data_geom_us %>%
  filter(day != 0) %>%
  ungroup() %>%
  dplyr::select(lot, day, conc_geom) %>%
  mutate(source = "obs")

colnames(data_obs_us) <- c("sampling_code", "day", "obs_conc", "source")

data_obs_us$day <- as.numeric(gsub("D", "", data_obs_us$day))
data_obs_us$sampday <- paste0(data_obs_us$sampling_code, "_", data_obs_us$day)

rm(data_geom_us)

## ------------------------------Wrangling - Model------------------------------
#Overall

overall_summary_stats <- overall %>%
  group_by(growth_model, loc, day) %>%
  summarize(median = round(median(count), 2), lq = round(quantile(count)[2], 2),
            uq = round(quantile(count)[4], 2)) %>%
  filter(day %in% c(5, 10))

#Model Validation

mm_strain_china_sum <- mm_strain_china %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "mm_strain_china") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(mm_strain_china_sum)[3] <- "pred_conc"

mm_strain_us_sum <- mm_strain_us %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "mm_strain_us") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(mm_strain_us_sum)[3] <- "pred_conc"

mm_apc_us_sum <- mm_apc_us %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "mm_apc_us") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(mm_apc_us_sum)[3] <- "pred_conc"

gp_apc_us_sum <- gp_apc_us %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "gp_apc_us") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(gp_apc_us_sum)[3] <- "pred_conc"

gp_apc_china_sum <- gp_apc_china %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "gp_apc_china") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(gp_apc_china_sum)[3] <- "pred_conc"

us_lots <- c(1:8, "overall")
us_rmse <- vector(mode = "logical", length = length(us_lots))
us_af <- vector(mode = "logical", length = length(us_lots))
us_bf <- vector(mode = "logical", length = length(us_lots))

us_perf_met <- bind_cols(us_lots, us_rmse, us_af, us_bf)
colnames(us_perf_met) <- c("lots", "rmse", "af", "bf")

china_lots <- c("S1", "S2", "S3", "F1", "F2", "F3", "F4", "R1", "R2", "R3", 
                "R4", "R5", "R6", "R7", "overall")
china_rmse <- vector(mode = "logical", length = length(china_lots))
china_af <- vector(mode = "logical", length = length(china_lots))
china_bf <- vector(mode = "logical", length = length(china_lots))

china_perf_met <- bind_cols(china_lots, china_rmse, china_af, china_bf)
colnames(china_perf_met) <- c("lots", "rmse", "af", "bf")

met_mm_strain_us <-  us_perf_met
met_mm_apc_us <-  us_perf_met
met_gp_apc_us <-  us_perf_met

met_mm_strain_china <-  china_perf_met
met_gp_apc_china <-  china_perf_met

## --------------------Combining Model and Experimental Data--------------------

#Assigning predicted conc to the observed conc

#MM Strain US

mm_strain_us_sum$obs_conc <- NA
for(i in 1:nrow(mm_strain_us_sum)){
  if(mm_strain_us_sum$sampday[i] %in% data_obs_us$sampday){
    index <- which(mm_strain_us_sum$sampday[i] == data_obs_us$sampday)
    mm_strain_us_sum$obs_conc[i] <- data_obs_us$obs_conc[index]
  }
}

mm_strain_us_sum <- mm_strain_us_sum %>%
  filter(!(is.na(obs_conc)))

for(i in 1:nrow(met_mm_strain_us)){
  if(met_mm_strain_us$lots[i] != "overall"){
    index <- met_mm_strain_us$lots[i]
    dat <- filter(mm_strain_us_sum, sampling_code == index)
    
    met_mm_strain_us$rmse[i] <- rmse(dat$obs_conc, dat$pred_conc)
    
    af_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(dat)))
    met_mm_strain_us$af[i] <- af
    
    bf_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc)))
    bf <- exp(bf_sum/nrow(dat))
    met_mm_strain_us$bf[i] <- bf
    
  } else {
    
    met_mm_strain_us$rmse[i] <- rmse(mm_strain_us_sum$obs_conc, mm_strain_us_sum$pred_conc)
    
    af_sum <- sum((log(mm_strain_us_sum$pred_conc) - log(mm_strain_us_sum$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(mm_strain_us_sum)))
    met_mm_strain_us$af[i] <- af
    
    bf_sum <- sum((log(mm_strain_us_sum$pred_conc) - log(mm_strain_us_sum$obs_conc)))
    bf <- exp(bf_sum/nrow(mm_strain_us_sum))
    met_mm_strain_us$bf[i] <- bf
    
  }
}

#MM Strain China

mm_strain_china_sum$obs_conc <- NA
for(i in 1:nrow(mm_strain_china_sum)){
  if(mm_strain_china_sum$sampday[i] %in% data_obs_china$sampday){
    index <- which(mm_strain_china_sum$sampday[i] == data_obs_china$sampday)
    mm_strain_china_sum$obs_conc[i] <- data_obs_china$obs_conc[index]
  }
}

mm_strain_china_sum <- mm_strain_china_sum %>%
  filter(!(is.na(obs_conc)))

for(i in 1:nrow(met_mm_strain_china)){
  if(met_mm_strain_china$lots[i] != "overall"){
    index <- met_mm_strain_china$lots[i]
    dat <- filter(mm_strain_china_sum, sampling_code == index)
    
    met_mm_strain_china$rmse[i] <- rmse(dat$obs_conc, dat$pred_conc)
    
    af_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(dat)))
    met_mm_strain_china$af[i] <- af
    
    bf_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc)))
    bf <- exp(bf_sum/nrow(dat))
    met_mm_strain_china$bf[i] <- bf
    
  } else {
    
    met_mm_strain_china$rmse[i] <- rmse(mm_strain_china_sum$obs_conc, mm_strain_china_sum$pred_conc)
    
    af_sum <- sum((log(mm_strain_china_sum$pred_conc) - log(mm_strain_china_sum$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(mm_strain_china_sum)))
    met_mm_strain_china$af[i] <- af
    
    bf_sum <- sum((log(mm_strain_china_sum$pred_conc) - log(mm_strain_china_sum$obs_conc)))
    bf <- exp(bf_sum/nrow(mm_strain_china_sum))
    met_mm_strain_china$bf[i] <- bf
    
  }
}

#MM APC US

mm_apc_us_sum$obs_conc <- NA
for(i in 1:nrow(mm_apc_us_sum)){
  if(mm_apc_us_sum$sampday[i] %in% data_obs_us$sampday){
    index <- which(mm_apc_us_sum$sampday[i] == data_obs_us$sampday)
    mm_apc_us_sum$obs_conc[i] <- data_obs_us$obs_conc[index]
  }
}

mm_apc_us_sum <- mm_apc_us_sum %>%
  filter(!(is.na(obs_conc)))

for(i in 1:nrow(met_mm_apc_us)){
  if(met_mm_apc_us$lots[i] != "overall"){
    index <- met_mm_apc_us$lots[i]
    dat <- filter(mm_apc_us_sum, sampling_code == index)
    
    met_mm_apc_us$rmse[i] <- rmse(dat$obs_conc, dat$pred_conc)
    
    af_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(dat)))
    met_mm_apc_us$af[i] <- af
    
    bf_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc)))
    bf <- exp(bf_sum/nrow(dat))
    met_mm_apc_us$bf[i] <- bf
    
  } else {
    
    met_mm_apc_us$rmse[i] <- rmse(mm_apc_us_sum$obs_conc, mm_apc_us_sum$pred_conc)
    
    af_sum <- sum((log(mm_apc_us_sum$pred_conc) - log(mm_apc_us_sum$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(mm_apc_us_sum)))
    met_mm_apc_us$af[i] <- af
    
    bf_sum <- sum((log(mm_apc_us_sum$pred_conc) - log(mm_apc_us_sum$obs_conc)))
    bf <- exp(bf_sum/nrow(mm_apc_us_sum))
    met_mm_apc_us$bf[i] <- bf
    
  }
}

#GP APC US

gp_apc_us_sum$obs_conc <- NA
for(i in 1:nrow(gp_apc_us_sum)){
  if(gp_apc_us_sum$sampday[i] %in% data_obs_us$sampday){
    index <- which(gp_apc_us_sum$sampday[i] == data_obs_us$sampday)
    gp_apc_us_sum$obs_conc[i] <- data_obs_us$obs_conc[index]
  }
}

gp_apc_us_sum <- gp_apc_us_sum %>%
  filter(!(is.na(obs_conc)))

for(i in 1:nrow(met_gp_apc_us)){
  if(met_gp_apc_us$lots[i] != "overall"){
    index <- met_gp_apc_us$lots[i]
    dat <- filter(gp_apc_us_sum, sampling_code == index)
    
    met_gp_apc_us$rmse[i] <- rmse(dat$obs_conc, dat$pred_conc)
    
    af_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(dat)))
    met_gp_apc_us$af[i] <- af
    
    bf_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc)))
    bf <- exp(bf_sum/nrow(dat))
    met_gp_apc_us$bf[i] <- bf
    
  } else {
    
    met_gp_apc_us$rmse[i] <- rmse(gp_apc_us_sum$obs_conc, gp_apc_us_sum$pred_conc)
    
    af_sum <- sum((log(gp_apc_us_sum$pred_conc) - log(gp_apc_us_sum$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(gp_apc_us_sum)))
    met_gp_apc_us$af[i] <- af
    
    bf_sum <- sum((log(gp_apc_us_sum$pred_conc) - log(gp_apc_us_sum$obs_conc)))
    bf <- exp(bf_sum/nrow(gp_apc_us_sum))
    met_gp_apc_us$bf[i] <- bf
    
  }
}

#GP APC China

gp_apc_china_sum$obs_conc <- NA
for(i in 1:nrow(gp_apc_china_sum)){
  if(gp_apc_china_sum$sampday[i] %in% data_obs_china$sampday){
    index <- which(gp_apc_china_sum$sampday[i] == data_obs_china$sampday)
    gp_apc_china_sum$obs_conc[i] <- data_obs_china$obs_conc[index]
  }
}

gp_apc_china_sum <- gp_apc_china_sum %>%
  filter(!(is.na(obs_conc)))

for(i in 1:nrow(met_gp_apc_china)){
  if(met_gp_apc_china$lots[i] != "overall"){
    index <- met_gp_apc_china$lots[i]
    dat <- filter(gp_apc_china_sum, sampling_code == index)
    
    met_gp_apc_china$rmse[i] <- rmse(dat$obs_conc, dat$pred_conc)
    
    af_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(dat)))
    met_gp_apc_china$af[i] <- af
    
    bf_sum <- sum((log(dat$pred_conc) - log(dat$obs_conc)))
    bf <- exp(bf_sum/nrow(dat))
    met_gp_apc_china$bf[i] <- bf
    
  } else {
    
    met_gp_apc_china$rmse[i] <- rmse(gp_apc_china_sum$obs_conc, gp_apc_china_sum$pred_conc)
    
    af_sum <- sum((log(gp_apc_china_sum$pred_conc) - log(gp_apc_china_sum$obs_conc))^2)
    af <- exp(sqrt(af_sum/nrow(gp_apc_china_sum)))
    met_gp_apc_china$af[i] <- af
    
    bf_sum <- sum((log(gp_apc_china_sum$pred_conc) - log(gp_apc_china_sum$obs_conc)))
    bf <- exp(bf_sum/nrow(gp_apc_china_sum))
    met_gp_apc_china$bf[i] <- bf
    
  }
}

## ------------------------------Plotting---------------------------------------

#Summarizing data for plotting
mm_strain_china_plot <- mm_strain_china %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "mm_strain_china") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(mm_strain_china_plot)[3] <- "pred_conc"

mm_strain_us_plot <- mm_strain_us %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "mm_strain_us") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(mm_strain_us_plot)[3] <- "pred_conc"

mm_apc_us_plot <- mm_apc_us %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "mm_apc_us") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(mm_apc_us_plot)[3] <- "pred_conc"

gp_apc_us_plot <- gp_apc_us %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "gp_apc_us") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(gp_apc_us_plot)[3] <- "pred_conc"

gp_apc_china_plot <- gp_apc_china %>%
  group_by(sampling_code, day) %>%
  summarize(med_conc = median(count),
            lq = quantile(count)[2],
            uq = quantile(count)[4]) %>%
  mutate(source = "gp_apc_china") %>%
  mutate(sampday = paste0(sampling_code, "_", day))
colnames(gp_apc_china_plot)[3] <- "pred_conc"

#Plotting
us_pred_plot <- bind_rows(mm_strain_us_plot,
                          mm_apc_us_plot,
                          gp_apc_us_plot)

us_pred_plot <- us_pred_plot %>%
  mutate(source = case_when(
    source %in% c("mm_strain_us") ~ "Approach 1: MM, Strain (US)",
    source %in% c("mm_apc_us") ~ "Approach 2: MM, APC (US)",
    source %in% c("gp_apc_us") ~ "Approach 3: GP, APC (US)",
    source %in% c("obs") ~ "Observed"
  ))

colnames(us_pred_plot)[6] <- "Source"

ggplot() + 
  geom_line(data = us_pred_plot, aes(x = day, y = pred_conc, group = Source, 
                                     linetype = Source)) + 
  facet_grid(rows = "sampling_code") + 
  geom_point(data = data_obs_us, aes(x = day, y = obs_conc), size = 0.5) +
  ylab("Concentration, log10 CFU/g") + xlab("Day") +
  scale_linetype_discrete(name = "Source") +
  scale_linetype_manual(values=c("solid","dashed", "dotted", "dotdash"))+
  theme_bw()


china_pred_plot <- bind_rows(mm_strain_china_plot,
                          gp_apc_china_plot)

china_pred_plot <- china_pred_plot %>%
  mutate(source = case_when(
    source %in% c("mm_strain_china") ~ "Approach 1: MM, Strain (China)",
    source %in% c("gp_apc_china") ~ "Approach 4: GP, APC (China)",
    source %in% c("obs") ~ "Observed"
  ))

colnames(china_pred_plot)[6] <- "Source"

plota <- ggplot() + 
  geom_line(data = filter(china_pred_plot,
                          sampling_code %in% c("S1", "S2", "S3",
                                               "F1", "F2", "F3", "F4")), 
                          aes(x = day, y = pred_conc, group = Source, 
                                     linetype = Source)) + 
  facet_grid(rows = "sampling_code") + 
  geom_point(data = filter(data_obs_china,
                           sampling_code %in% c("S1", "S2", "S3",
                                                "F1", "F2", "F3", "F4")), 
             aes(x = day, y = obs_conc), size = 0.5) +
  ylab("Concentration, log10 CFU/g") + xlab("Day") +
  scale_linetype_discrete(name = "source") +
  scale_linetype_manual(values=c("solid","dashed", "dotted", "dotdash"))+
  theme_bw()

plotb <-ggplot() + 
  geom_line(data = filter(china_pred_plot, 
                          sampling_code %in% c("R2", "R3",
                                               "R4", "R5", "R6",
                                               "R7")),
            aes(x = day, y = pred_conc, group = Source, 
                                        linetype = Source)) + 
  facet_grid(rows = "sampling_code") + 
  geom_point(data = filter(data_obs_china, 
                           sampling_code %in% c("R2", "R3",
                                                "R4", "R5", "R6",
                                                "R7")), 
             aes(x = day, y = obs_conc), size = 0.5) +
  ylab("Concentration, log10 CFU/g") + xlab("Day") +
  scale_linetype_discrete(name = "source") +
  scale_linetype_manual(values=c("solid","dashed", "dotted", "dotdash"))+
  theme_bw()

ggarrange(plota, plotb + rremove("ylab"), common.legend =  TRUE)
