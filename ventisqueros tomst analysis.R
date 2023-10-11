## Ventisqueros TOMST analyses

library(tidyverse) #need lubridate (is part of tidyverse)
library(plotrix)
library(lme4)
library(lmerTest)
tomstdata <- read.csv("TOMST_master_oct20_june22.csv")

#pull out data since Oct 2021 - March 2022

# 2020-08-10 09:45:00

tomstdata$tidydate_CR <- ymd_hms(tomstdata$datetime_CR)

tomstdata <- tomstdata %>% 
  mutate(tidydate_CR = ymd_hms(tomstdata$datetime_CR),
          Day = date(tomstdata$tidydate_CR),
         Hour = hour(tomstdata$tidydate_CR),
         Month = month(tomstdata$tidydate_CR),
         Year = year(tomstdata$tidydate_CR)) %>% 
  filter(Day >= '2021-10-01', Day <='2022-02-02') %>% 
  rename(Collar_ID = Collar)
  
  
tomst <- tomstdata %>%  select(sensorID, subsoilT, surfaceT, airT, soilm_count, filename, Collar_ID, Block, Plot, Treatment, tidydate_CR, Day, Hour, Month, Year)


### soil moisture calculations
#from TOMST_calibrations_calcs_w_curves.xlsx

# manual VWC curve for both  soil types: y = 0.0002x - 0.1795. R2 = 0.93
# manual gravimetric soil m (GSM) curve for both soil types: y = -2E-07x2 + 0.001x - 0.7437, R2 = 0.93
# these equations are from excel; they don't actually map that well onto the trendline in excel even though
  # the equation and trendline go together.

tomst <- tomst %>% 
  mutate(VWC_cm3cm3 = ((0.0002*soilm_count)-0.1795), 
         GSM_gg = ((-0.0000002*(soilm_count^2)) + (0.001*soilm_count) - 0.7437),
         timeofday = if_else(Hour >=7 & Hour <=20, "day", "night"),
         season = if_else(Month %in% c("10","11","12"), "Wet", "Dry"))

#now we are going to summarize data by Hour, Day, Month, and Year so that we are ready to merge wiht CO2 data

## within-hour measurements averaged for each collar

tidy_temp <- tomst %>% 
  select(Hour,Day, Month, Year, Treatment, Collar_ID, subsoilT, surfaceT, airT) %>% 
  mutate(season = if_else(Month %in% c("10","11","12"), "Wet", "Dry")) %>% 
  pivot_longer(cols = c( subsoilT, surfaceT, airT), # The columns we want to collapse into two new columns
               names_to = "sensor", # The name we want to assign to the column that will contain our current column names
               values_to = "temp_C")

tidy_mois <- tomst %>% 
  select(Hour,Day, Month, Year, Treatment, Collar_ID, VWC_cm3cm3, GSM_gg) %>% 
  mutate(season = if_else(Month %in% c("10","11","12"), "Wet", "Dry")) %>% 
  pivot_longer(cols = c(VWC_cm3cm3, GSM_gg), # The columns we want to collapse into two new columns
               names_to = "sensor", # The name we want to assign to the column that will contain our current column names
               values_to = "soilm")

tidy_temp %>% 
  group_by(Hour, Treatment, sensor) %>% 
  ggplot(aes(x=Hour, y=temp_C, fill=Treatment))+
  geom_bar(stat = "summary", fun = "mean") +
  facet_wrap(~sensor) +
  ggtitle("Avg hourly temperatures by treatment")


tidy_temp %>% 
  filter(Hour == 13) %>% 
  group_by(Treatment, sensor) %>% 
  summarise(
    sd = sd(temp_C, na.rm = TRUE),
    temp_C = mean(temp_C)) %>%
  ggplot(aes(x=Treatment, y=temp_C, fill=Treatment))+
  geom_bar(stat = "summary", fun = "mean", position =position_dodge(width=0.9), width=0.7) +
  geom_errorbar(aes(ymin=temp_C-sd, ymax=temp_C+sd), width=0.2, position=position_dodge(width=0.9), width=0.7) +
          facet_wrap(~sensor) +
  labs(title = "Afternoon temperatures by treatment", caption = "error bars represent sd")
  
tidy_temp %>% 
  filter(Hour == 4) %>% 
  group_by(Treatment, sensor) %>% 
  summarise(
    sd = sd(temp_C, na.rm = TRUE),
    temp_C = mean(temp_C)) %>%
  ggplot(aes(x=Treatment, y=temp_C, fill=Treatment))+
  geom_bar(stat = "summary", fun = "mean", position =position_dodge(width=0.9), width=0.7) +
  geom_errorbar(aes(ymin=temp_C-sd, ymax=temp_C+sd), width=0.2, position=position_dodge(width=0.9), width=0.7) +
  facet_wrap(~sensor) +
  labs(title = "Nighttime temperatures by treatment", caption = "error bars represent sd")

tidy_mois %>% 
  ggplot(aes(x=as.integer(Month), y=soilm, fill=Treatment))+
  geom_bar(stat = "summary", fun = "mean") +
  facet_wrap(~sensor) +
  ylab("soil moisture (cm3/cm3 or g H20/gsoil)") +
  ggtitle("Soil moisture by month")

tidy_mois %>% 
  group_by(season, Treatment, sensor) %>% 
  summarise(
    sd = sd(soilm, na.rm = TRUE),
    soilm = mean(soilm)) %>% 
  ggplot(aes(x=sensor, y=soilm, fill=Treatment))+
  geom_bar(stat = "summary", fun = "mean", position =position_dodge(width=0.9), width=0.7) +
  geom_errorbar(aes(ymin=soilm-sd, ymax=soilm+sd), width=0.2, position=position_dodge(width=0.9), width=0.7) +
  facet_wrap(~season) + 
  ylab("soil moisture (cm3/cm3 or g H20/gsoil)") +
  ggtitle("Soil moisture by treatment") +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

### calculate the difference between control and treatment and plot difference along time of day

head(tomst)
head(tidy_temp)

# convert to wide format with columns for OTC temps and columns for contorl temps to then be able to subtract
widetomst <- tomst %>% 
  select(subsoilT, surfaceT, airT,VWC_cm3cm3, GSM_gg, Block, Plot, Treatment, Day, Hour, Month, Year) %>% 
  group_by(Block, Plot, Treatment, Day, Hour, Month, Year) %>% 
  summarise(
    mean_soilT = mean(subsoilT),
    mean_surfaceT = mean(surfaceT),
    mean_airT = mean(airT),
    mean_VWC = mean(VWC_cm3cm3),
    mean_GSM = mean(GSM_gg)
              ) %>% 
  group_by(Block, Plot, Day) %>% 
  pivot_wider(names_from = Treatment, values_from = c(mean_soilT, mean_surfaceT, mean_airT, mean_VWC, mean_GSM))

# make new columns for the difference between warmed and control for each temperature
# diff = warm - control. if positive, warm is hotter, if negative, warm is cooler
widetomst <- widetomst %>% 
  mutate(soilTdiff = mean_soilT_warming - mean_soilT_control,
         surfaceTdiff = mean_surfaceT_warming - mean_surfaceT_control,
         airTdiff = mean_airT_warming - mean_airT_control,
         VWCdiff = mean_VWC_warming - mean_VWC_control,
         GSMdiff = mean_GSM_warming - mean_GSM_control)

## plot difference for each temp by hour; need to put into long again


longtemp <- widetomst %>% 
  select(Block, Plot, Day, Hour, Month, Year, soilTdiff, surfaceTdiff, airTdiff) %>% 
  pivot_longer(cols=c(soilTdiff, surfaceTdiff, airTdiff), names_to = "sensor", values_to ="temp_C_diff")

longmois <- widetomst %>% 
  select(Block, Plot, Day, Hour, Month, Year, VWCdiff, GSMdiff) %>% 
  pivot_longer(cols=c(VWCdiff, GSMdiff), names_to = "soilm_meas", values_to ="soilm_diff")

# there are some NA rows I think because those plots just don't have sensors either at all or in both OTCs and control

longtemp %>% 
  group_by(Hour, sensor) %>% 
  summarise(
    sd = sd(temp_C_diff, na.rm = TRUE),
    mean_temp_diff = mean(temp_C_diff, na.rm=TRUE)) %>% 
  ggplot(aes(x=Hour, y= mean_temp_diff)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_temp_diff-sd, ymax=mean_temp_diff+sd), width=0.2) +
  facet_wrap(~sensor, ncol=1) +
  geom_hline(yintercept=0, linetype="dashed", color="blue", size=1) +
  labs(title = "Temperature difference warmed - control by hour", caption = "error bars represent sd")
## want to change order to air, surface, then soil top to bottom
  
longmois %>% 
  group_by(Month, soilm_meas) %>% 
  summarise(
    sd = sd(soilm_diff, na.rm = TRUE),
    mean_mois_diff = mean(soilm_diff, na.rm=TRUE)) %>% 
  ggplot(aes(x=Month, y= mean_mois_diff)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_mois_diff-sd, ymax=mean_mois_diff+sd), width=0.2) +
  facet_wrap(~soilm_meas, ncol=1) +
  geom_hline(yintercept=0, linetype="dashed", color="blue", size=1) +
  labs(title = "Soil moisture difference warmed - control by season and measure", caption = "error bars represent sd")

## want to change month numbers to names and order Oct, Nov, Dec, Jan, Feb


##### summarize differences between control and warming temps and soil m
# group by trt, block, hour for temp and trt block, month for soil m


temp_summary <- 
  tomst %>% 
  group_by(Treatment, timeofday) %>% 
  summarise(
    mean_soilT = mean(subsoilT, na.rm =TRUE),
    sd_soilT = sd(subsoilT, na.rm = TRUE),
    mean_surfaceT = mean(surfaceT, na.rm = TRUE),
    sd_surfaceT = sd(surfaceT, na.rm = TRUE),
    mean_airT = mean(airT, na.rm=TRUE),
    sd_airT =sd(airT, na.rm =TRUE)
  )

write.csv(x=temp_summary, "temp_summary.csv")  


airTmod<-lmer(airT~Treatment*timeofday + (1|Collar_ID), data=tomst)
summary(airTmod)

surfaceTmod<-lmer(surfaceT~Treatment*timeofday + (1|Collar_ID), data=tomst)
summary(surfaceTmod)

soilTmod<-lmer(subsoilT~Treatment*timeofday + (1|Collar_ID), data=tomst)
summary(soilTmod)


mois_summary <- 
  tomst %>% 
  group_by(season, Treatment) %>% 
  summarise(
    mean_VWC = mean(VWC_cm3cm3, na.rm =TRUE),
    sd_VWC = sd(VWC_cm3cm3, na.rm = TRUE),
    mean_GSM = mean(GSM_gg, na.rm = TRUE),
    sd_GSM = sd(GSM_gg, na.rm = TRUE)
  )

write.csv(x=mois_summary, "mois_summary.csv")  

VWCmod<-lmer(VWC_cm3cm3~Treatment*season + (1|Collar_ID), data=tomst)
summary(VWCmod)

GSMmod<-lmer(GSM_gg~Treatment*season + (1|Collar_ID), data=tomst)
summary(GSMmod)
