## ------------------------------Title------------------------------------------
# Analyzing the China-supply chain data

## ------------------------------Description------------------------------------
# Project: Bacterial growth modeling (spinach)

# Script description: Analyzing data from the following samplings:
#S4, S5, S6, S7, S8, S9, F1, F2, F3, F4, R1, R2, R3, R4, R5, R6, R7

## ------------------------------Packages---------------------------------------

library(dplyr); library(lmerTest); library(ggplot2); library(tidyverse)
library(reshape2); library(forcats)

## ------------------------------Read in raw data-------------------------------

counts <- read.csv("sampling_data/raw/microdata_fullV2.csv", header = TRUE, stringsAsFactors = FALSE)
sample_mass <- read.csv("sampling_data/raw/testweights_full.csv", header = TRUE, stringsAsFactors = FALSE)

#S1 - S3 data, from Murphy et al. 2022 (https://github.com/IvanekLab/SpinachSpoilage)
s1_s3 <- read.csv("sampling_data/raw/ShelfLifeMicroData.csv", header = TRUE, stringsAsFactors = FALSE)

## ------------------------------Wrangling--------------------------------------

#Ensure consistent variable-level names
#Change L3, SL5, Sl7 and SL10 to D3, D5, D7 and D10, respectively 

sample_mass <- sample_mass %>%
  mutate(stage = case_when(
    stage == "SL3" ~ "D3",
    stage == "SL5" ~ "D5",
    stage == "SL7" ~ "D7",
    stage == "SL10" ~ "D10",
    .default = as.character(stage)
  ))

counts <- counts %>%
  mutate(stage = case_when(
    stage == "SL3" ~ "D3",
    stage == "SL5" ~ "D5",
    stage == "SL7" ~ "D7",
    stage == "SL10" ~ "D10",
    .default = as.character(stage)
  ))

#Treating P as D0 for lots: F1, F2, F3, F4, R1, R2, R3, R4, R5, R6, R7

sample_mass$stage <- ifelse(
  (sample_mass$sampling %in% c("F1", "F2", "F3", "F4", "R1",
                              "R2", "R3", "R4", "R5", "R6", "R7") &
    sample_mass$stage == "P"),
  "D0",
  sample_mass$stage)

counts$stage <- ifelse(
  (counts$sampling %in% c("F1", "F2", "F3", "F4", "R1",
                              "R2", "R3", "R4", "R5", "R6", "R7") &
    counts$stage == "P"),
  "D0",
  counts$stage)

#Removing samples that were not packed by the company, thus we were unable to 
#collect shelf life data. These are the D0 samples for: S5, S6 and S8

#Removed samples that were not tested from S2 and S3. 

sample_mass_removed <- sample_mass %>%
  filter(grepl("NotPacked", weight_before_dilution)) %>%
  filter(grepl("NA", weight_before_dilution))

sample_mass_v2 <-  sample_mass %>%
  filter(!grepl("NotPacked", weight_before_dilution)) %>%
  filter(!grepl("NA", weight_before_dilution))

counts_removed <- counts %>%
  filter(grepl("Not Packed", notes)) %>%
  filter(grepl("NA", raw_count))

counts_v2 <-  counts %>%
  filter(!grepl("Not Packed", notes)) %>%
  filter(!grepl("NA", raw_count))

#Checking whether each plate replicate was counted from the same dilution
#The plate replicates are from the same dilution, so the plate count data can be
#averaged

counts_v2 %>%
  group_by(location, sampling, stage, sample_rep, test) %>%
  mutate(same_dilution = dilution[1] == dilution[2]) %>%
  with(table(same_dilution))

#Four samples were TNTC. These were substituted with a 25% increase of the upper countable limit
#The upper countable limit was 300 for AC Petrifilm (based on 3M's Interpretation Guide)
#The upper countable limit was 150 for CC Petrifilm (based on 3M's Interpretation Guide)
#The observations that were TNTC were:
#(i) location L, sampling F4, stage D10, sample_rep 1, APC and GN

counts_v2$raw_count[counts_v2$raw_count == "TNTC" & counts_v2$test == "APC"] <- 300*1.25
counts_v2$raw_count[counts_v2$raw_count == "TNTC" & counts_v2$test == "GN"] <- 150*1.25

counts_v2$raw_count <- as.numeric(counts_v2$raw_count)

counts_v3 <- counts_v2 %>%
  group_by(location, sampling, stage, sample_rep, test, dilution) %>%
  mutate(averaged_count = mean(raw_count)) %>%
  distinct(location, sampling, stage, sample_rep, test, dilution, averaged_count)

counts_v3$merge_col <- paste0(counts_v3$location, "_",
                              counts_v3$sampling, "_",
                              counts_v3$stage, "_",
                              counts_v3$sample_rep)


sample_mass_v2$merge_col <- paste0(sample_mass_v2$location, "_",
                                   sample_mass_v2$sampling, "_",
                                   sample_mass_v2$stage, "_",
                                   sample_mass_v2$sample_rep)

data <- merge(counts_v3, sample_mass_v2, by = "merge_col",) %>%
  dplyr::select(location.x, sampling.x, stage.x, sample_rep.x, test, weight_before_dilution,
         weight_after_dilution, dilution, averaged_count)

colnames(data)[1:4] <- c("location", "sampling", "stage", "sample_rep")

data[, c("weight_before_dilution", "weight_after_dilution",
         "dilution", "averaged_count")] <- sapply(data[, c("weight_before_dilution", "weight_after_dilution",
                                                           "dilution", "averaged_count")],
                                                  as.numeric)

#The final mass, after dilution, was not recorded for some samples; for these samples,
#I have added the mean buffer mass, of the remaining samples, to sample mass

buffer_mass <- data %>%
  filter(!is.na(weight_after_dilution)) %>%
  mutate(buffer_mass = weight_after_dilution - weight_before_dilution)

mean_buffer_mass <- mean(buffer_mass$buffer_mass)

data$weight_after_dilution <- ifelse(is.na(data$weight_after_dilution),
                                     (data$weight_before_dilution + mean_buffer_mass),
                                     data$weight_after_dilution)

#Calculating the concentration of the samples replicates

data_arith <- data %>%
  mutate(conc_arith = log10(averaged_count*
           (weight_after_dilution/weight_before_dilution)*
           (10^dilution))) %>%
  dplyr::select(-c("weight_before_dilution", "weight_after_dilution", "dilution",
                   "averaged_count"))

data_geom <- data_arith %>%
  group_by(location, sampling, stage, test) %>%
  mutate(conc_geom = mean(conc_arith)) %>%
  distinct(location, sampling, stage, test, conc_geom)


data_geom$season <- ifelse(grepl("S", data_geom$sampling), 
                           "winter",
                           "summer")

## ------------------------------Analysis---------------------------------------

#These data (i.e., samplings S1, S2, S3) were published in Murphy et al. 2022.
#Thus, they will be left out of analysis here
data_geom_analysis <- data_geom %>%
  filter(!(sampling %in% c("S1", "S2", "S3")))
#Test whether bacterial concentration of the eCommerce samples was significantly 
#different at the end of shelf life

e_com_test_data <- data_geom_analysis %>%
  filter(stage == "D10") %>%
  filter(!sampling %in% c("R7", "F3", "R1")) 

e_com_test_data %>%
  group_by(test, location) %>%
  summarize(median = median(conc_geom), 
            percentile25 = summary(conc_geom)[2],
            percentile75 = summary(conc_geom)[5])

e_com_test_data_subset <- e_com_test_data %>%
  filter(test == "APC", location == "E")

wilcox.test(conc_geom ~ location, data = filter(e_com_test_data, test == "APC"), paired = TRUE)
wilcox.test(conc_geom ~ location, data = filter(e_com_test_data, test == "GN"), paired = TRUE)
wilcox.test(conc_geom ~ location, data = filter(e_com_test_data, test == "YM"), paired = TRUE)


## ------------------------------Summary Stats----------------------------------

summary_stats <- data_geom_analysis %>%
  group_by(location, stage, test) %>%
  summarize(median = quantile(conc_geom, type = 7)[3],
            lq = quantile(conc_geom, type = 7)[2],
            uq = quantile(conc_geom, type = 7)[4],
            n = n())

#Comparing net growth by distribution channel

net_growth_h_d0 <- data_geom_analysis %>%
  filter(!(sampling %in% c("S5", "S6", "S8"))) %>%
  filter(stage %in% c("H", "D0")) %>%
  group_by(sampling, test) %>%
  mutate(h_conc = conc_geom[stage == "H"]) %>%
  ungroup() %>%
  filter(stage == "D0") %>%
  mutate(net_growth = conc_geom - h_conc)

net_growth_summary_h_d0 <- net_growth_h_d0 %>%
  group_by(location, test) %>%
  summarize(mean = mean(net_growth), 
            min = min(net_growth),
            max = max(net_growth),
            n = n())

net_growth_d0_d10 <- data_geom_analysis %>%
  filter(!(sampling %in% c("S4", "S5", "S6", "S7", "S8", "S9"))) %>%
  filter(stage %in% c("D0", "D10")) %>%
  group_by(sampling, test) %>%
  mutate(day0_conc = conc_geom[stage == "D0"]) %>%
  ungroup() %>%
  filter(stage == "D10") %>%
  mutate(net_growth = conc_geom - day0_conc)

net_growth_summary_d0_d10 <- net_growth_d0_d10 %>%
  group_by(location, test) %>%
  summarize(mean = mean(net_growth), 
            min = min(net_growth),
            max = max(net_growth),
            n = n())

paired_net_growth_d0_d10 <- net_growth_d0_d10 %>%
  filter(sampling %in% c("F1", "F2", "F4", "R2", 
                         "R3", "R4", "R5", "R6")) %>%
  dplyr::select("sampling", "test", "net_growth", "location") %>%
  pivot_wider(names_from = location, values_from = net_growth)


#Net growth by distribution channel 
paired_apc_net_growth <- filter(paired_net_growth_d0_d10, test == "APC")
wilcox.test(paired_apc_net_growth$E, paired_apc_net_growth$L, paired = TRUE)

paired_gn_net_growth <- filter(paired_net_growth_d0_d10, test == "GN")
wilcox.test(paired_gn_net_growth$E, paired_gn_net_growth$L, paired = TRUE)

paired_ym_net_growth <- filter(paired_net_growth_d0_d10, test == "YM")
wilcox.test(paired_ym_net_growth$E, paired_ym_net_growth$L, paired = TRUE)


## ------------------------------Plotting---------------------------------------

data_geom$stage <- factor(data_geom$stage, levels = c("H", "T", "R", "D0", 
                                                      "D2", "D3", "D4", 
                                                      "D5", "D6", "D7",
                                                      "D8", "D10"))

data_geom_plot <- filter(data_geom, !(sampling %in% c("S1", "S2", "S3")))
data_geom_plot$location <- ifelse(data_geom_plot$location == "E",
                                  "eCommerce", 
                                  "Local")
data_geom_plot$location <- factor(data_geom_plot$location, 
                                  levels = c("Local", "eCommerce"))

data_geom_boxplot <- filter(data_geom_plot, 
                            stage %in% c("H", "D0", "D10"))


ggplot() + 
  geom_point(data = data_geom_plot, 
             aes(x = factor(stage), y = conc_geom, color = location, group = location),
             alpha = 0, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6)) + 
  facet_wrap( ~ test, nrow = 3) + 
  scale_color_brewer(palette = "Dark2", name = "Distribution channel") + 
  theme(strip.text = element_text(size = 7), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) +
  geom_boxplot(data = data_geom_boxplot, 
               aes(x = factor(stage), y = conc_geom, color = location),
               outlier.size=2, outlier.shape = NA, linewidth = 0.3, width = 0.4,
               position = position_dodge2(preserve = "single")) + 
  geom_point(data = data_geom_plot, 
             aes(x = factor(stage), y = conc_geom, color = location, group = location),
             size = 0.3,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) +
  ylab(expression("Concentration, log"[10]*"CFU/g")) + xlab("Sample") + 
  theme_bw()

## ---------------------------Modeling Data Wrangling---------------------------

s1_s3_geom_by_location <- s1_s3 %>%
  mutate(lot = case_when(
    lot == "L1" ~ "S1",
    lot == "L2" ~ "S2",
    lot == "L3" ~ "S3"
  )) %>%
  group_by(location, lot, stage, test) %>%
  summarise(conc_geom = mean(log10count)) %>%
  mutate(season = "winter")

#The counts from the s1_s3_geom_by_location df match those in
#the data_geom dataframe for local s1, s2, s3;

#For modeling, will average the shelf life data over distribution channel (except for s1, s2, s3)
#Will use the data from the s1_s3 df for: s1, s2, s3. For these, will use the local distribution data
#For the other samplings, will use the information in the data_geom sheet

s1_s3_geom_mod <- s1_s3 %>%
  mutate(lot = case_when(
    lot == "L1" ~ "S1",
    lot == "L2" ~ "S2",
    lot == "L3" ~ "S3"
  )) %>%
  filter(location != "E") %>%
  group_by(lot, stage, test) %>%
  summarise(conc_geom = mean(log10count)) %>%
  mutate(season = "winter")

s1_s3_geom_mod$stage <- ifelse(s1_s3_geom_mod$stage == "DI", 
                           "D3",
                           s1_s3_geom_mod$stage)

colnames(s1_s3_geom_mod)[1] <- "sampling"

data_geom_mod <- data_arith %>%
  filter(!(sampling %in% c("S1", "S2", "S3"))) %>%
  group_by(sampling, stage, test) %>%
  mutate(conc_geom = mean(conc_arith)) %>%
  distinct(sampling, stage, test, conc_geom)

data_geom_mod$season <- ifelse(grepl("S", data_geom_mod$sampling), 
                           "winter",
                           "summer")

mod_data <- bind_rows(s1_s3_geom_mod, data_geom_mod)

mod_data_summary <- mod_data %>%
  group_by(stage, test) %>%
  summarize(mean = mean(conc_geom),
            min = min(conc_geom),
            max = max(conc_geom),
            n = n())

#Sampling_specific starting concentration
starting_conc <- data_arith %>%
  filter(stage == "D0" & test == "APC") %>%
  group_by(sampling) %>%
  summarize(min = min(conc_arith), max = max(conc_arith))

## ------------------------------Export-----------------------------------------


#Export data

write.csv(mod_data, paste0("sampling_data/wrangled/concentration_by_day_modeling_data.csv"), 
          row.names = FALSE)
write.csv(starting_conc, paste0("sampling_data/wrangled/china_starting_conc.csv"), 
          row.names = FALSE)
