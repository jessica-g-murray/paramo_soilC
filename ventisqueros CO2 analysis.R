# Ventisqueros CO2 analysis Fall 2023

library(tidyverse)
library(plotrix)
library(lme4)
library(lmerTest)
library(rstatix)

data <- read.csv("CO2_master_oct2023.csv") ### flux was calculated along with offset using SoilFluxPro software.
head(data)

data$tidydate = ymd_hms(data$Date)

data <- data %>% 
  mutate(Day = date(data$tidydate),
         Hour = hour(data$tidydate),
         Month = month(data$tidydate),
         Year = year(data$tidydate)) 

### inspect CO2 data

min(data$Exp_Flux, na.rm=TRUE)
max(data$Exp_Flux, na.rm=TRUE)

min(data$Exp_R2, na.rm=TRUE) 
max(data$Exp_R2, na.rm=TRUE)

hist(data$Exp_Flux, na.rm=TRUE)

#remove negative fluxes and fluxes with R2 less than 0.8

cleandat <- data %>% 
  filter(Exp_R2 > 0.8) %>% 
  mutate(CO2flux = if_else(Exp_Flux <0, NA, Exp_Flux),
         season = if_else(Month %in% c("10","11","12"), "Wet", "Dry"),
         month_ch = case_when(
           Month == "1" ~ "January",
           Month == "2" ~ "February",
           Month == "10" ~ "October",
           Month == "11" ~ "November",
           Month == "12" ~ "December"
           )
         )


# Plot flux by treatment, month

cleandat %>% 
  mutate(month_ch = fct_relevel(month_ch, "October", "November", "December", "January", "February")) %>% 
  group_by(Treatment, month_ch, Block) %>% 
  summarise(mean_CO2 = mean(CO2flux, na.rm = TRUE),
           sd_CO2 = sd(CO2flux, na.rm=TRUE)) %>% 
  ggplot(aes(x=month_ch, y=mean_CO2, color=Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_CO2-sd_CO2, ymax=mean_CO2+sd_CO2), width=.2)+
  facet_wrap(~Block)


CO2aov <- aov(CO2flux~Treatment * month_ch*Block + Error(Collar), data=cleandat)
summary(CO2aov)
## tells me I need to control for block and month in all analyses (to investigate mechanisms behind different C cycling responses)

## https://landscape-agroecology.com/mixed-effects-models-in-r/
CO2lme <-lmer(CO2flux~Treatment *season * Block + (1 | Collar), data=cleandat)
summary(CO2lme)

CO2summary <-
  cleandat %>% 
  group_by(season, Treatment) %>% 
  summarise(
    mean_CO2 = mean(Exp_Flux, na.rm =TRUE),
    sd_CO2 = sd(Exp_Flux, na.rm = TRUE)
  )
write.csv(x=CO2summary, "CO2_summary.csv")  
