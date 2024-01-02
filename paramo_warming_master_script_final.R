## cleaned up, final paramo master script
# 12/31/23
library(tidyverse)
library(lme4)
library(lmerTest)
library(rstatix)
library(car)
library(patchwork)
library(flux)
library(performance)
library(interactions)
library(corrplot)
library(MuMIn)
library(ggplot2)

tomst <- read.csv("TOMST_master_oct20_june22.csv")
CO2 <- read.csv("CO2_master_nov2023.csv")
MBC <- read.csv("MBC.csv")
CN <- read.csv("CN_inorgNmaster.csv")
veg <- read.csv("collar_veg_areas.csv")
pH <- read.csv("chirripo-soil-pH.csv")

#### clean up CO2 data

CO2$tidydate = ymd_hms(CO2$Date) #format date

CO2 <- CO2 %>% #extract date objects
  mutate(Day = date(CO2$tidydate),
         Hour = hour(CO2$tidydate),
         Month = month(CO2$tidydate),
         Year = year(CO2$tidydate)) 

cleanCO2 <- CO2 %>% 
  filter(Exp_R2 > 0.8) %>% # subset out fluxes with R2 of less than 0.8
  mutate(CO2flux_umols_m2_s = if_else(Exp_Flux <0, NA_integer_, Exp_Flux), ## want positive fluxes only; negative fluxes are instrument error
         season = if_else(Day < "2021-12-01", "wet", "dry"), ## assign seasons
         month_ch = case_when( #assign months
           Month == "1" ~ "January",
           Month == "2" ~ "February",
           Month == "10" ~ "October",
           Month == "11" ~ "November",
           Month == "12" ~ "December"),
         CO2_mgC_m2_h = CO2flux_umols_m2_s/1000000*12.001*1000*3600 # convert umols/m2/s to mgC/m2/h
  )

cleanCO2 <- cleanCO2 %>% #select columns we will use
  select(Date, Block, Plot, Treatment, sensorID, tidydate, Day, Hour, Month, Year, Collar, CO2flux_umols_m2_s, CO2_mgC_m2_h, season, month_ch)

cleanCO2 <- cleanCO2 %>% # combine days of measurements into single timepoint- since some CO2 measurements happened over 2-3 days
  mutate(timepoint = case_when(Day == "2021-10-04" | Day == "2021-10-03" ~ "10/3/21",
                               Day== "2021-10-25" | Day =="2021-10-26" | Day =="2021-10-24" ~ "10/25/21",
                               Day == "2021-11-05" | Day == "2021-11-06"  ~ "11/05/21",
                               Day =="2021-11-16" | Day == "2021-11-17" ~ "11/16/21",
                               Day =="2021-12-01" | Day == "2021-12-03" ~ "12/1/21",
                               Day =="2021-12-17" ~ "12/17/21",
                               Day =="2022-01-04" ~ "1/04/22",
                               Day =="2022-01-21" | Day == "2022-01-23" ~ "1/21/22",
                               Day =="2022-02-02" | Day == "2022-02-03" ~ "2/02/22")) 

cleanCO2 <- cleanCO2 %>% 
  mutate(timepoint_tidy = mdy(cleanCO2$timepoint),
         tp_day = date(timepoint_tidy),#make timepoint a Date object
         plotid = paste0(Block, "-", Plot)) #plotid object in case we want it later

cleanCO2$Block <- as.factor(cleanCO2$Block)
cleanCO2$Plot <- as.factor(cleanCO2$Plot)


######## Clean up tomst data

tomst$tidydate <- ymd_hms(tomst$datetime_CR)

## Make new columns for Day, Hour, Month, Year
tomst <- tomst %>% 
  mutate(tidydate = ymd_hms(tomst$datetime_CR),
         Day = date(tomst$tidydate),
         Hour = hour(tomst$tidydate),
         Month = month(tomst$tidydate),
         Year = year(tomst$tidydate)) %>% 
  filter(Day >= '2021-10-03', Day <='2022-02-02')

tomst <- tomst %>% 
  mutate(VWC_cm3cm3 = ((0.00000000011338*(soilm_count^2))+(0.0002821*soilm_count))-0.22013, 
         GSM_gg = ((-0.0000002*(soilm_count^2)) + (0.001*soilm_count) - 0.7437),
         timeofday = if_else(Hour >=7 & Hour <=20, "day", "night"),
         season = if_else(Day < "2021-12-01", "wet", "dry"))

tomst <- tomst %>% 
  mutate(VWC_cm3cm3 = if_else(VWC_cm3cm3 <0, NA_integer_, VWC_cm3cm3))

### How does microclimate differ by treatment?

airTmod<-lmer(airT~Treatment*season + Treatment*timeofday + (1|Block/Collar), data=tomst)
summary(airTmod)

surfaceTmod<-lmer(surfaceT~Treatment*season + Treatment*timeofday + (1|Block/Collar), data=tomst)
summary(surfaceTmod)

soilTmod<-lmer(subsoilT~Treatment*season + Treatment*timeofday + (1|Block/Collar), data=tomst)
summary(soilTmod)

VWCmod<-lmer(VWC_cm3cm3~Treatment*season + (1|Block/Collar), data=tomst)
summary(VWCmod)

### OTC effects on temp and moisture figures

      #### air temperature ####
p1 <- tomst %>% 
  group_by(timeofday, Treatment) %>% 
  summarize(meanairT = mean(airT, na.rm=TRUE),
            se = sd(airT, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=timeofday, y=meanairT, color=Treatment)) + 
  geom_point(size=4) +
  ylim(0,13) +
  geom_errorbar(aes(ymin=meanairT-se, ymax=meanairT+se), width=0.2) +
  ggtitle("air temperature") +
  ylab("air temperature (C)") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

p2 <- tomst %>% 
  group_by(season, Treatment) %>% 
  summarize(meanairT = mean(airT, na.rm=TRUE),
            se = sd(airT, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=season, y=meanairT, color=Treatment)) + 
  geom_point(size=4) +
  ylim(0,13) +
  geom_errorbar(aes(ymin=meanairT-se, ymax=meanairT+se), width=0.2) +
  ggtitle("air temperature") +
  ylab("air temperature (C)") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

  ### soil temp ###
p3 <- tomst %>% 
  group_by(timeofday, Treatment) %>% 
  summarize(meansoilT = mean(subsoilT, na.rm=TRUE),
            se = sd(subsoilT, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=timeofday, y=meansoilT, color=Treatment)) + 
  geom_point(size=4) +
  ylim(0,13) +
  geom_errorbar(aes(ymin=meansoilT-se, ymax=meansoilT+se), width=0.2) +
  ggtitle("soil temperature") +
  ylab("soil temperature (C)") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

p4 <- tomst %>% 
  group_by(season, Treatment) %>% 
  summarize(meansoilT = mean(subsoilT, na.rm=TRUE),
            se = sd(airT, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=season, y=meansoilT, color=Treatment)) + 
  geom_point(size=4) +
  ylim(0,13) +
  geom_errorbar(aes(ymin=meansoilT-se, ymax=meansoilT+se), width=0.2) +
  ggtitle("soil temperature") +
  ylab("soil temperature (C)") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

  #### soil moisture ###
p5 <- tomst %>% 
  group_by(season, Treatment) %>% 
  summarize(meanVWC = mean(VWC_cm3cm3, na.rm=TRUE),
            se = sd(VWC_cm3cm3, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=season, y=meanVWC, color=Treatment)) + 
  geom_point(size=4) +
  geom_errorbar(aes(ymin=meanVWC-se, ymax=meanVWC+se), width=0.2) +
  ggtitle("soil moisture") +
  ylab("volumetric soil moisture (cm3/cm3)") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

microclim_panel  <- p1 + p3 + p2 + p4 + p5 & theme(legend.position = "bottom")

microclim_panel + plot_layout(guides = "collect", nrow=3)

##### alternative way (more concise) to present temperature and moisture data

tidy_temp <- tomst %>% 
  select(Hour,Day, Month, Year, Treatment, timeofday, Collar, subsoilT, surfaceT, airT) %>% 
  mutate(season = if_else(Month %in% c("10","11","12"), "Wet", "Dry")) %>% 
  pivot_longer(cols = c( subsoilT, surfaceT, airT), # The columns we want to collapse into two new columns
               names_to = "sensor", # The name we want to assign to the column that will contain our current column names
               values_to = "temp_C")

tidy_temp$Treatment <- fct_relevel(tidy_temp$Treatment, "warming") # to get warm cols behind con

##put both air and soil temp in same plots

p1 <- tidy_temp %>% 
  filter(sensor != "surfaceT") %>% 
  group_by(season, Treatment,sensor) %>% 
  summarize(meanT = mean(temp_C, na.rm=TRUE),
            se = sd(temp_C, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=season, y=meanT, color=Treatment, shape=sensor)) + 
  geom_point(size=8) +
  ylim(0,15) +
  geom_errorbar(aes(ymin=meanT-se, ymax=meanT+se), width=0.2) +
  #ggtitle("Seasonal effects of OTCs \non temperature") +
  ylab("Temperature (C)") +
  xlab("Season") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

p2 <- tidy_temp %>% 
  filter(sensor != "surfaceT") %>% 
  group_by(Treatment,timeofday,sensor) %>% 
  summarize(meanT = mean(temp_C, na.rm=TRUE),
            se = sd(temp_C, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=timeofday, y=meanT, color=Treatment, shape=sensor)) + 
  geom_point(size=8) +
  ylim(0,15) +
  geom_errorbar(aes(ymin=meanT-se, ymax=meanT+se), width=0.2) +
  #ggtitle("Diurnal effects of OTCs \non temperature") +
  ylab("Temperature (C)") +
  xlab("Time of Day") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))


tidy_mois <- tomst %>% 
  select(Hour,Day, Month, Year, Treatment, Collar, VWC_cm3cm3, GSM_gg) %>% 
  mutate(season = if_else(Month %in% c("10","11","12"), "Wet", "Dry")) %>% 
  pivot_longer(cols = c(VWC_cm3cm3, GSM_gg), # The columns we want to collapse into two new columns
               names_to = "sensor", # The name we want to assign to the column that will contain our current column names
               values_to = "soilm")

p3 <- tidy_mois %>% 
  filter(sensor == "VWC_cm3cm3") %>% 
  group_by(season, Treatment, sensor) %>% 
  summarise(se = sd(soilm, na.rm=TRUE)/sqrt(n()),
            soilm = mean(soilm)) %>% 
  ggplot(aes(x=Treatment, y=soilm, fill=Treatment))+
  geom_bar(stat = "summary", fun = "mean", position =position_dodge(width=0.9), width=0.7) +
  geom_errorbar(aes(ymin=soilm-se, ymax=soilm+se), width=0.2, position=position_dodge(width=0.9), width=1) +
  facet_wrap(~season) + 
  ylab("Soil moisture (cm3/cm3)") +
  labs(caption = "error bars represent se") +
 #ggtitle("Seasonal effects of OTCs \non soil moisture") +
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))


microclim_panel  <- p1 + p2 + p3 

microclim_panel + plot_layout(guides = "collect", nrow=2) & theme(legend.position = "right")

## bar plot temperatures by hour. for appendix
tidy_temp %>% 
  filter(sensor!= "surfaceT") %>% 
  group_by(Hour, Treatment, sensor) %>% 
  ggplot(aes(x=Hour, y=temp_C, color=Treatment))+
  geom_bar(stat = "summary", fun = "mean", position = "identity", fill = "white", alpha=0.2) +
  facet_wrap(~sensor) +
  ggtitle("Avg hourly temps by treatment") +
  ylab("Degrees Celsius") +
  theme(legend.position = "none")+
  theme_bw()+
  theme(axis.text = element_text(size=20), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))


### graph of soil moisture by timepoint that shows when dry vs wet season is. 
    ### this can go in the appendix.
tomst %>% 
  group_by(Day) %>% 
  summarize(mean_VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            se = sd(VWC_cm3cm3, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=Day, y=mean_VWC))+
  geom_point() +
  geom_errorbar(aes(ymin=mean_VWC-se, ymax=mean_VWC+se), width=1) +
  ylab("VWC cm3/cm3") +
  xlab("Month") +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3")) +
  labs(title = "VWC ~ month")+
  scale_x_date(date_breaks = "2 weeks") +
  theme_bw()+
  theme(axis.text = element_text(size=20), axis.text.x = element_text(angle=45, vjust=1, hjust=1), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24))


### Summary tables

tomst_summary_seasonoverall <- tomst %>% 
  group_by(season) %>% 
  summarize(mean_soilT = mean(subsoilT, na.rm=TRUE),
            sd_soilT = sd(subsoilT, na.rm=TRUE),
            mean_airT = mean(airT, na.rm=TRUE),
            sd_airT =sd(airT, na.rm=TRUE),
            mean_VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            sd_VWC = sd(VWC_cm3cm3, na.rm=TRUE))

tomst_summary_todoverall <- tomst %>% 
  group_by(timeofday) %>% 
  summarize(mean_soilT = mean(subsoilT, na.rm=TRUE),
            sd_soilT = sd(subsoilT, na.rm=TRUE),
            mean_airT = mean(airT, na.rm=TRUE),
            sd_airT =sd(airT, na.rm=TRUE),
            mean_VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            sd_VWC = sd(VWC_cm3cm3, na.rm=TRUE))

tomst_summary_trt <- tomst %>% 
  group_by(Treatment) %>% 
  summarize(mean_soilT = mean(subsoilT, na.rm=TRUE),
            sd_soilT = sd(subsoilT, na.rm=TRUE),
            mean_airT = mean(airT, na.rm=TRUE),
            sd_airT =sd(airT, na.rm=TRUE),
            mean_VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            sd_VWC = sd(VWC_cm3cm3, na.rm=TRUE))

tomst_summary_tod <- tomst %>% 
  group_by(timeofday, Treatment) %>% 
  summarize(mean_soilT = mean(subsoilT, na.rm=TRUE),
            sd_soilT = sd(subsoilT, na.rm=TRUE),
            mean_airT = mean(airT, na.rm=TRUE),
            sd_airT =sd(airT, na.rm=TRUE),
            mean_VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            sd_VWC = sd(VWC_cm3cm3, na.rm=TRUE))

tomst_summary_season <- tomst %>% 
  group_by(season, Treatment) %>% 
  summarize(mean_soilT = mean(subsoilT, na.rm=TRUE),
            sd_soilT = sd(subsoilT, na.rm=TRUE),
            mean_airT = mean(airT, na.rm=TRUE),
            sd_airT =sd(airT, na.rm=TRUE),
            mean_VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            sd_VWC = sd(VWC_cm3cm3, na.rm=TRUE))


############ MERGE CO2 & TOMST ############ 


## averaging microclim data by the hour so it can be merged with CO2.

hourlytomst <- tomst %>%  
  group_by(Collar, Day, Hour, Month, Year, season, timeofday) %>% 
  summarise(soilT = mean(subsoilT, na.rm=TRUE),
            surfaceT = mean(surfaceT, na.rm=TRUE),
            airT = mean(airT, na.rm = TRUE),
            VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            sd_VWC = sd(VWC_cm3cm3, na.rm=TRUE),
            GSM = mean(GSM_gg, na.rm = TRUE))

## merge

tomstCO2 <- left_join(cleanCO2, hourlytomst, by=c("Collar", "Day", "Hour", "Month", "Year", "season"))


#### average observations by timepoint (at least 3 observations were taken each timepoint, these are meant to be averaged)
### 9 timepoints * 35 collars = 315 (actual N= 292 due to collar issues, missed observations, etc)

continuous_CO2_dat <- tomstCO2 %>% 
  group_by(tp_day, Block, Plot, Treatment, season, Collar) %>% 
  summarize(mean_CO2 = mean(CO2_mgC_m2_h, na.rm = TRUE),
            mean_flux_uumolsm2h = mean(CO2flux_umols_m2_s, na.rm=TRUE),
            se_CO2 = sd(CO2_mgC_m2_h, na.rm=TRUE)/sqrt(n()),
            mean_airT = mean(airT, na.rm=TRUE),
            mean_soilT = mean(soilT, na.rm=TRUE),
            mean_surfaceT = mean(surfaceT, na.rm=TRUE),
            mean_VWC = mean(VWC, na.rm=TRUE),
            mean_VWC_sd = mean(sd_VWC, na.rm=TRUE))

continuous_CO2_dat <- ungroup(continuous_CO2_dat) ## scale gives NaN if df is grouped

continuous_CO2_dat <-continuous_CO2_dat %>% 
  mutate(sairT = scale(mean_airT, center=TRUE, scale=TRUE),
         ssoilT = scale(mean_soilT, center=TRUE, scale=TRUE),
         ssurfaceT = scale(mean_surfaceT, center=TRUE, scale=TRUE),
         sVWC = scale(mean_VWC, center=TRUE, scale = TRUE),
         sVWC_sd = scale(mean_VWC_sd, center=TRUE, scale=TRUE),
         ssoilTday = scale(soilT_day, center=TRUE, scale=TRUE),
         ssoilTnight = scale(soilT_night, center=TRUE, scale=TRUE),
         ssurfaceTnight = scale(surfaceT_night, center=TRUE, scale=TRUE),
         ssurfaceTday = scale(surfaceT_day, center=TRUE, scale=TRUE),
         sairtTday = scale(airT_day, center=TRUE, scale=TRUE),
         sairTnight = scale(airT_night, center=TRUE, scale=TRUE))


## get the average temp and moisture values for later merge with soil vars
tomst_tod_wide <- tomst_tod_summary %>% 
  pivot_wider(id_cols=c(Collar, season), values_from = c(soilT, surfaceT, airT, VWC), names_from=timeofday)

## merge
continuous_CO2_dat <- continuous_CO2_dat %>% 
  left_join(tomst_tod_wide, by= c("Collar", "season"))

## look at distributions

continuous_CO2_dat %>%
  ggplot(aes(x=mean_CO2))+
  geom_histogram() # positively skewed, with one apparent outlier above 225 and maybe a group of them after 160

continuous_CO2_dat <-continuous_CO2_dat %>% 
  mutate(sqrt_CO2= sqrt(mean_CO2)) 

continuous_CO2_dat %>%
  ggplot(aes(x=sqrt_CO2))+
  geom_histogram() #this fixes skew, only one seeming outlier remaining

out <- boxplot.stats(continuous_CO2_dat$sqrt_CO2)$out
out_ind <- which(continuous_CO2_dat$sqrt_CO2 %in% c(out))
out_ind ## this confirms that row # 40 (2-17-con on 10-25-2021) is the outlier

continuous_CO2_dat <- continuous_CO2_dat %>% 
  mutate(mean_CO2 = if_else(mean_CO2 > 225, NA_integer_, mean_CO2)) ## can't just do "NA" because this is a numeric column. https://stackoverflow.com/questions/53636644/r-if-else-assign-na-value

#############################################
#### Rs models #############################


## timepoint*treatment
continuous_CO2_dat$tp_factor = as.factor(continuous_CO2_dat$tp_day)

summary(lmer(sqrt_CO2 ~ tp_day*Treatment + (1|Block/Collar), data=continuous_CO2_dat))
r.squaredGLMM(lmer(sqrt_CO2 ~ tp_day*Treatment + (1|Block/Collar), data=continuous_CO2_dat))

summary(lmer(sqrt_CO2 ~ season*Treatment + (1|Block/Collar), data=continuous_CO2_dat))
r.squaredGLMM(lmer(sqrt_CO2 ~ season*Treatment + (1|Block/Collar), data=continuous_CO2_dat))

##### ln (log quadratic) model of rs ~ soilT, sensu Carey et al paper, table S3
continuous_CO2_dat <- continuous_CO2_dat %>% 
  mutate(ln_Rs = log10(mean_CO2),
         meansoilT2 = mean_soilT^2,
         ssoilT2= ssoilT^2)

summary(lmer(ln_Rs ~ ssoilT + ssoilT2 + (1|Block/Collar), data=continuous_CO2_dat))
r.squaredGLMM(lmer(ln_Rs ~ ssoilT + ssoilT2 + (1|Block/Collar), data=continuous_CO2_dat))

summary(lmer(ln_Rs ~ (ssoilT + ssoilT2) * Treatment + (1|Block/Collar), data=continuous_CO2_dat))
r.squaredGLMM(lmer(ln_Rs ~ (ssoilT + ssoilT2) * Treatment + (1|Block/Collar), data=continuous_CO2_dat))

summary(lmer(ln_Rs ~ (ssoilT + ssoilT2) * sVWC * Treatment + (1|Block/Collar), data=continuous_CO2_dat))
r.squaredGLMM(lmer(ln_Rs ~ (ssoilT + ssoilT2) * sVWC + (1|Block/Collar), data=continuous_CO2_dat))



##### graphs to accompany Rs models

#### Rs ~ timepoint * treatment
continuous_CO2_dat %>% 
  group_by(tp_day, Treatment) %>% 
  summarize(CO2_mean = mean(mean_CO2, na.rm=TRUE),
            se = sd(mean_CO2, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=tp_day, y = CO2_mean, color = Treatment, group=Treatment)) +
  geom_point(size=3) +
  geom_line(linewidth=2) +
  geom_vline(xintercept=as.numeric(as.Date("2021-12-01")), linetype="dashed") +
  scale_x_date(date_breaks = "1 month") + ## posix error happens if you have date on x and don't use scale_x_date()
  geom_errorbar(aes(ymin=CO2_mean-se, ymax=CO2_mean+se), width=.5) +
  labs(caption = "error bars represent se", y="CO2-C mg/m2/h", x = "measurement date") + 
  ggtitle("CO2 respiration across time") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

### Rs summary

continuous_CO2_dat %>% 
  group_by(season) %>% 
  summarize(meanCO2 = mean(mean_CO2, na.rm=TRUE),
            sdCO2 = sd(mean_CO2, na.rm=TRUE))


### spatial variation in Rs

continuous_CO2_dat %>% 
  filter(Plot != 13) %>% 
  group_by(Plot, Treatment) %>% 
  summarize(mean = mean(mean_CO2, na.rm=TRUE),
            se = sd(mean_CO2, na.rm=TRUE)/sqrt(n())) %>% 
  ggplot(aes(x=Plot, y=mean, fill=Treatment)) +
  geom_bar(stat = "summary", fun = "mean", position =position_dodge(width=0.9), width=0.7) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(width=0.9), width=1) +
  xlab("Experimental Plot") +
  ylab("CO2 respiration (mg CO2-C / m2/h)") +
  theme_bw()+
  theme(axis.text = element_text(size=24), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

continuous_CO2_dat %>% 
  filter(Plot != 13) %>% 
  group_by(Plot, Treatment) %>% 
  ggplot(aes(x=Plot, y=mean_CO2, fill=Treatment)) +
  geom_boxplot() +
  xlab("Experimental Plot") +
  ylab("CO2 respiration (mg CO2-C / m2/h)") +
  theme_bw()+
  theme(axis.text = element_text(size=24), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))
############################

### cumulative CO2

##############################

continuous_CO2_dat$num_date <- as.numeric(continuous_CO2_dat$tp_day)*24 # hours since 1/1/1970
auc <- continuous_CO2_dat %>%
  filter(Collar !="3-13-con") %>% 
  select(Collar, Treatment, Block, Plot, num_date, mean_CO2) %>% 
  group_by(Collar, Treatment, Block, Plot) %>% 
  summarise(cumulativeCO2_mgC_m2 = auc(num_date, mean_CO2))

aucmod <- aov(cumulativeCO2_mgC_m2 ~ Treatment, data=auc)
summary(aucmod)


##################### SOIL VARIABLES #######################
## pare down each dataframe to the necessary variables, add any calculated variables, 
  # remove outliers, and then run the ANOVA (using lm to be consistent with reporting t test statistics)
########## MBC ##############
MBC$Block <- as.factor(MBC$Block)

## hadn't previously considered scaling soil vars but I think it is a good idea since they are orders of magnitude different!

MBC <- MBC %>% 
  select(Collar, Block, Plot, Treatment, mgC_gsoil_uf, mgN_gsoil_uf, MBC, MBN) %>% 
  mutate(plotid = paste0(Block, "-", Plot),
         DCN = mgC_gsoil_uf/mgN_gsoil_uf,
         MBCN = MBC/MBN)

## DOC
MBC %>%
  ggplot(aes(x=mgC_gsoil_uf))+
  geom_histogram() 

out <- boxplot.stats(MBC$mgC_gsoil_uf)$out
out_ind <- which(MBC$mgC_gsoil_uf %in% c(out))
out_ind ## outlier >0.3, 1 observation removed (this removed the negative MBC value as well)



MBC <- MBC %>% 
  mutate(mgC_gsoil_uf = if_else(mgC_gsoil_uf >0.3, NA_integer_, mgC_gsoil_uf))

MBC <-MBC %>% 
  mutate(sqrt_DOC= sqrt(mgC_gsoil_uf))

MBC %>% 
  ggplot(aes(x=Treatment, y = mgC_gsoil_uf, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "DOC~Treatment by block") +
  ylab("DOC mgC/gsoil") +
  labs(title = "DOC by treatment")+
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_DOC~Treatment, data=MBC)) # no trt diff in DOC

## TDN

MBC %>%
  ggplot(aes(x=mgN_gsoil_uf))+
  geom_histogram()

out <- boxplot.stats(MBC$mgN_gsoil_uf)$out
out_ind <- which(MBC$mgN_gsoil_uf %in% c(out))
out_ind ## no outliers

MBC <-MBC %>% 
  mutate(sqrt_TDN= sqrt(mgN_gsoil_uf))

MBC %>% 
  ggplot(aes(x=Treatment, y = mgN_gsoil_uf, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "TDN~Treatment by block") +
  ylab("TDN mgN/gsoil") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_TDN~Treatment, data=MBC)) # no trt diff in TDN

### DCN
MBC %>%
  ggplot(aes(x=DCN))+
  geom_histogram()

out <- boxplot.stats(MBC$DCN)$out
out_ind <- which(MBC$DCN %in% c(out))
out_ind ## four outliers but I only see one really

MBC <- MBC %>% 
  mutate(DCN = if_else(DCN >25, NA_integer_, DCN))

MBC <-MBC %>% 
  mutate(sqrt_DCN= sqrt(DCN))


MBC %>% 
  ggplot(aes(x=Treatment, y = DCN, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "DCN~Treatment by block") +
  ylab("dissolved C/N ratio") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_DCN~Treatment, data=MBC)) # significant trt difference in DCN

### MBC

MBC %>%
  ggplot(aes(x=MBC))+
  geom_histogram()

out <- boxplot.stats(MBC$MBC)$out
out_ind <- which(MBC$MBC %in% c(out))
out_ind ## 31 is an outlier 

MBC <- MBC %>% 
  mutate(MBC = if_else(MBC < 0, NA_integer_, MBC))

MBC <-MBC %>% 
  mutate(sqrt_MBC= sqrt(MBC))

MBC %>% 
  ggplot(aes(x=Treatment, y = MBC, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "Microbial biomass carbon") +
  ylab("Microbial biomass C mgC/gsoil") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_MBC~Treatment, data=MBC)) # no trt diff in MBC

## MBN

MBC %>%
  ggplot(aes(x=MBN))+
  geom_histogram() #apparent outlier over 0.08

out <- boxplot.stats(MBC$MBN)$out
out_ind <- which(MBC$MBN %in% c(out))
out_ind ## no outliers according to this test


MBC <-MBC %>% 
  mutate(sqrt_MBN= sqrt(MBN))

MBC %>% 
  ggplot(aes(x=Treatment, y = MBN, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "Microbial biomass nitrogen") +
  ylab("Microbial biomass N mgC/gsoil") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_MBN~Treatment, data=MBC)) # no trt diff in MBN

### MBCN


MBC %>%
  ggplot(aes(x=MBCN))+
  geom_histogram() 

out <- boxplot.stats(MBC$MBCN)$out
out_ind <- which(MBC$MBCN %in% c(out))
out_ind ## three outliers, one less than 0 and 1 over 20


MBC <- MBC %>% 
  mutate(MBCN = if_else(MBCN > 20, NA_integer_, MBCN))

MBC <- MBC %>% 
  mutate(MBCN = if_else(MBCN < 0, NA_integer_, MBCN))

MBC <-MBC %>% 
  mutate(sqrt_MBCN= sqrt(MBCN)) 

MBC %>% 
  ggplot(aes(x=Treatment, y = MBCN, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "Microbial biomass C/N") +
  ylab("Microbial biomass C/N") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_MBCN~Treatment, data=MBC)) # no trt diff in MBN


######### vegetation cover ###########
veg <- veg %>% 
  select(Collar, grid_1, grid_2, grid_3) %>% 
  mutate_if(is.integer, as.numeric)
veg<- veg%>% 
  mutate(avg_cov = rowMeans(select(.,c(grid_1, grid_2, grid_3))))

veg %>%
  ggplot(aes(x=avg_cov))+
  geom_histogram() # not exactly normal but no apparent outliers


veg <- veg %>% 
  select(Collar, avg_cov) %>% 
  mutate(sqrt_cov = sqrt(avg_cov))

############# pH data ############

pH <- pH %>% 
  select(block, unit,treat, soil.pH) %>% 
  mutate(Collar = paste(block, unit, treat, sep="-"))

pH <- pH %>% 
  select(Collar, soil.pH)

## veg and pH aren't combined with metadata yet so will look at trt differences after it is merged into big dataset

################# CN ################

CN$Block <- as.factor(CN$Block)
CN <- CN %>% select(Collar, Block, Plot, Treatment, mgNH4_gsoil, mgNO3_gsoil, mgNH4_kgsoil, mgNO3_kgsoil, percC, percN, C.N)

## percC

CN %>%
  ggplot(aes(x=percC))+
  geom_histogram() 

out <- boxplot.stats(CN$percC)$out
out_ind <- which(CN$percC %in% c(out))
out_ind ## no outliers according to this test


CN <-CN %>% 
  mutate(sqrt_percC= sqrt(percC))

CN %>% 
  ggplot(aes(x=Treatment, y = percC, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "soil % C") +
  ylab("soil % C") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_percC~Treatment, data=CN)) # no trt diff in %C

## percN
CN %>%
  ggplot(aes(x=percN))+
  geom_histogram() 

out <- boxplot.stats(CN$percN)$out
out_ind <- which(CN$percN %in% c(out))
out_ind ## no outliers according to this test


CN <-CN %>% 
  mutate(sqrt_percN= sqrt(percN))

CN %>% 
  ggplot(aes(x=Treatment, y = percN, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "soil % N") +
  ylab("soil % N") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_percN~Treatment, data=CN)) # no trt diff in %N
## NO3

CN %>%
  ggplot(aes(x=mgNO3_gsoil))+
  geom_histogram() 

out <- boxplot.stats(CN$mgNO3_gsoil)$out
out_ind <- which(CN$mgNO3_gsoil %in% c(out))
out_ind ## no outliers according to this test


CN <-CN %>% 
  mutate(sqrt_NO3= sqrt(mgNO3_gsoil))

CN %>% 
  ggplot(aes(x=Treatment, y = mgNO3_gsoil, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "NO3") +
  ylab("NO3, mg N / g soil") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_NO3~Treatment, data=CN)) # no trt diff in %N

## C/N

CN %>%
  ggplot(aes(x=C.N))+
  geom_histogram() 

out <- boxplot.stats(CN$C.N)$out
out_ind <- which(CN$C.N %in% c(out))
out_ind ## no outliers according to this test


CN <-CN %>% 
  mutate(sqrt_CN= sqrt(C.N))

CN %>% 
  ggplot(aes(x=Treatment, y = C.N, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "soil C/N") +
  ylab("soil C/N ratio") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_CN~Treatment, data=CN)) # no trt diff in %N

### NH4
CN %>% 
  ggplot(aes(x=mgNH4_kgsoil)) +
  geom_histogram()

out <- boxplot.stats(CN$mgNH4_kgsoil)$out
out_ind <- which(CN$mgNH4_kgsoil %in% c(out))
out_ind ## two outliers according to this test, 3-3-war and 4-11-war

CN <- CN %>% 
  mutate(mgNH4_kgsoil = if_else(mgNH4_kgsoil > 14, NA_integer_, mgNH4_kgsoil))

CN <-CN %>% 
  mutate(sqrt_NH4= sqrt(mgNH4_kgsoil)) 


CN %>% 
  ggplot(aes(x=Treatment, y = mgNH4_kgsoil, fill=Treatment)) +
  geom_boxplot() +
  labs(title = "NH4") +
  ylab("NH4, mg N / g soil") +
  theme_bw()+
  theme(axis.text = element_text(size=16), strip.text = element_text(size=20), 
        axis.title= element_text(size=20), plot.title = element_text(size=24)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))

summary(lm(sqrt_NH4~Treatment, data=CN)) # Nh4 higher in warmed plots, even with outliers removed (which were both high vals for warmed plots)


################## make our master dataframe #########################

## get average CO2 rates
CO2_avg_rates <- continuous_CO2_dat %>%  
  group_by(Collar, Treatment) %>% 
  summarise(mean_CO2_mg_m2h = mean(mean_CO2, na.rm=TRUE),
            mean_airtemp = mean(mean_airT, na.rm=TRUE),
            mean_soiltemp = mean(mean_soilT, na.rm=TRUE),
            mean_soilm = mean(mean_VWC, na.rm=TRUE))

### start joining
# cumulative CO2 with average rates
alldat <- left_join(CO2_avg_rates, auc, by = c("Collar", "Treatment"))

alldat <- left_join(alldat, daily_temps_wide, by = c("Collar", "Treatment"))

MBC$Plot <- as.factor(MBC$Plot)
#MBC
alldat <- left_join(alldat, MBC, by= c("Collar", "Block", "Plot", "Treatment") )

CN$Plot <- as.factor(CN$Plot)
#CN
alldat <- left_join(alldat, CN, by=c("Collar", "Block", "Plot", "Treatment"))

#veg
alldat <- left_join(alldat, veg, by="Collar")

#pH
alldat <- left_join(alldat, pH, by="Collar")

### add BD

alldat <-alldat %>% 
  mutate(soil_bd = case_when(Block == "1" ~ 1.00,
                             Block == "2" ~ 1.00,
                             Block == "3" ~ 0.66,
                             Block == "4" ~ 0.87))

## roughly calculate mass-specific respiration rates. bulk density data isn't super reliable though
alldat <- alldat %>% 
  mutate(cumulativeCO2_mgC_cm3 = cumulativeCO2_mgC_m2/10000/5,
         mean_CO2_mg_cm3 =mean_CO2_mg_m2h/10000/5)

alldat <- alldat %>% 
  mutate(cumulativeCO2_mgC_gsoil = cumulativeCO2_mgC_cm3 / soil_bd,
         mean_CO2_mg_gsoil = mean_CO2_mg_cm3/ soil_bd)

alldat <- alldat %>% 
  mutate(cumulativeCO2_mgC02C_mgMBC = cumulativeCO2_mgC_gsoil / MBC,
         mean_CO2_mgCO2C_mgMBC = mean_CO2_mg_gsoil/ MBC)

### roughly calculate SOC stocks. same bulk density problem

alldat <- alldat %>% 
  mutate(SOC_stock_kgha = percC/100 * soil_bd* 5 * 100000,
         SOC_stock_Mgha = percC/100 * soil_bd* 5 * 100,
         SOC_stock_kgm2 = percC/100 * soil_bd * 5 * 10,
         resp_kgC_ha_y = cumulativeCO2_mgC_m2/1000000/0.33699) #perc C (g/g) * bd (g/cm3) * sampling depth of 5 cm 

alldat %>% 
  summarise(sum_respkgC_ha_y = sum(resp_kgC_ha_y, na.rm=TRUE))

## do the models for pH and veg cover
summary(lm(soil.pH~Treatment, data=alldat))
summary(lm(sqrt_cov~Treatment, data=alldat))

### Does SOC determine C losses?

summary(lm(cumulativeCO2_mgC_m2~ percC*Treatment, data=alldat)) 

#### Summary table for variables

alldat <- ungroup(alldat)
trt_means <- alldat %>% 
  group_by(Treatment) %>% 
  summarise(DOC = mean(mgC_gsoil_uf, na.rm=TRUE),
            sdDOC =sd(mgC_gsoil_uf, na.rm=TRUE),
            TDN = mean(mgN_gsoil_uf, na.rm=TRUE),
            sdTDN =sd(mgN_gsoil_uf, na.rm=TRUE),
            MBC_mean = mean(MBC, na.rm=TRUE),
            sdMBC = sd(MBC, na.rm=TRUE),
            MBN_mean = mean(MBN, na.rm=TRUE),
            sdMBN =sd(MBN, na.rm=TRUE),
            meanMBCN = mean(MBCN, na.rm=TRUE),
            sdMBCN = sd(MBCN, na.rm=TRUE),
            meanDCN = mean(DCN, na.rm=TRUE),
            sdDCN = sd(DCN, na.rm=TRUE),
            NH4 = mean(mgNH4_kgsoil, na.rm=TRUE),
            sdNH4 = sd(mgNH4_kgsoil, na.rm=TRUE),
            NO3 = mean(mgNO3_kgsoil, na.rm=TRUE),
            sdNO3 = sd(mgNO3_kgsoil, na.rm=TRUE),
            percC_mean = mean(percC, na.rm=TRUE),
            sdpercC = sd(percC, na.rm=TRUE),
            percN_mean = mean(percN, na.rm=TRUE),
            sdpercN = sd(percN, na.rm=TRUE),
            CN = mean(C.N, na.rm=TRUE),
            sdCN = sd(C.N, na.rm=TRUE),
            mean_CO2_flux = mean(mean_CO2_mg_m2h, na.rm=TRUE),
            sd_CO2 = sd(mean_CO2_mg_m2h, na.rm=TRUE),
            mean_cumCO2 = mean(cumulativeCO2_mgC_m2, na.rm=TRUE),
            sd_cumCO2 = sd(cumulativeCO2_mgC_m2, na.rm=TRUE),
            mean_cov = mean(avg_cov, na.rm=TRUE),
            sd_cov = sd(avg_cov, na.rm=TRUE),
            soil_pH = mean(soil.pH, na.rm=TRUE),
            sd_pH = sd(soil.pH, na.rm=TRUE))


overall_means <- alldat %>% 
  summarise(DOC = mean(mgC_gsoil_uf, na.rm=TRUE),
            sdDOC =sd(mgC_gsoil_uf, na.rm=TRUE),
            TDN = mean(mgN_gsoil_uf, na.rm=TRUE),
            sdTDN =sd(mgN_gsoil_uf, na.rm=TRUE),
            MBC_mean = mean(MBC, na.rm=TRUE),
            sdMBC = sd(MBC, na.rm=TRUE),
            MBN_mean = mean(MBN, na.rm=TRUE),
            sdMBN =sd(MBN, na.rm=TRUE),
            NH4 = mean(mgNH4_kgsoil, na.rm=TRUE),
            sdNH4 = sd(mgNH4_kgsoil, na.rm=TRUE),
            NO3 = mean(mgNO3_kgsoil, na.rm=TRUE),
            sdNO3 = sd(mgNO3_kgsoil, na.rm=TRUE),
            percC_mean = mean(percC, na.rm=TRUE),
            sdpercC = sd(percC, na.rm=TRUE),
            percN_mean = mean(percN, na.rm=TRUE),
            sdpercN = sd(percN, na.rm=TRUE),
            CN = mean(C.N, na.rm=TRUE),
            sdCN = sd(C.N, na.rm=TRUE),
            dry_CO2_flux = mean(dry_mean_CO2_mg_m2h, na.rm=TRUE),
            sd_dryCO2 = sd(dry_mean_CO2_mg_m2h, na.rm=TRUE),
            wet_CO2_flux = mean(wet_mean_CO2_mg_m2h, na.rm=TRUE),
            sd_wetCO2 =sd(wet_mean_CO2_mg_m2h, na.rm=TRUE),
            mean_cumCO2 = mean(cumulativeCO2_mgC_m2, na.rm=TRUE),
            sd_cumCO2 = sd(cumulativeCO2_mgC_m2, na.rm=TRUE),
            mean_cov = mean(avg_cov, na.rm=TRUE),
            sd_cov = sd(avg_cov, na.rm=TRUE),
            soil_pH = mean(soil.pH, na.rm=TRUE),
            sd_pH = sd(soil.pH, na.rm=TRUE))


################### correlations ##########################
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software#google_vignette

# spearman's correlation is unaffected by transformation; no need to use sqrt versions of data


corrdat <- alldat %>% 
  select(c(cumulativeCO2_mgC_m2,mgC_gsoil_uf,mgN_gsoil_uf,
           MBC,MBN,DCN,MBCN,mgNH4_gsoil,mgNO3_gsoil,percC,percN,
           C.N,avg_cov, soil.pH))
corrdat <- ungroup(corrdat) #otherwise there is an error

scor <- corrdat %>% 
  select(where(is.numeric)) %>%  
  cor(method = "spearman", use="complete.obs")
plot.new()

corrplot(scor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 1)
## visualizing correlation matrix https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

testRes <- corrdat %>%
  select(where(is.numeric)) %>%  
  cor.mtest(.)

plot.new()
corrplot(scor, p.mat=testRes$p, method='color', insig='label_sig', sig.level=c(0.001, 0.05), pch.cex=0.9, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, col = COL2('BrBG'), tl.cex = 0.5)

corrplot(scor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, col = COL2('BrBG'), tl.cex = 1.2)
