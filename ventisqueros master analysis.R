### ventisqueros master analyses
library(tidyverse)
library(lme4)
library(lmerTest)
library(rstatix)
library(car)
#library(broom)

tomst <- read.csv("TOMST_master_oct20_june22.csv")
CO2 <- read.csv("CO2_master_oct2023.csv")
MBC <- read.csv("MBC.csv")
CN <- read.csv("CN_inorgNmaster.csv")

CO2$tidydate = ymd_hms(CO2$Date)

CO2 <- CO2 %>% 
  mutate(Day = date(CO2$tidydate),
         Hour = hour(CO2$tidydate),
         Month = month(CO2$tidydate),
         Year = year(CO2$tidydate)) 

cleanCO2 <- CO2 %>% 
  filter(Exp_R2 > 0.8) %>% 
  mutate(CO2flux_umols_m2_s = if_else(Exp_Flux <0, NA, Exp_Flux),
         season = if_else(Month %in% c("10","11","12"), "Wet", "Dry"),
         month_ch = case_when(
           Month == "1" ~ "January",
           Month == "2" ~ "February",
           Month == "10" ~ "October",
           Month == "11" ~ "November",
           Month == "12" ~ "December"),
         CO2_mgC_m2_h = CO2flux_umols_m2_s/1000000*12.001*1000*3600
           )

cleanCO2 <- cleanCO2 %>% 
  select(Date, Block, Plot, Treatment, sensorID, tidydate, Day, Hour, Month, Year, Collar, CO2flux_umols_m2_s, CO2_mgC_m2_h, season, month_ch)

tomst$tidydate <- ymd_hms(tomst$datetime_CR)
tomst <- tomst %>% 
  mutate(tidydate = ymd_hms(tomst$datetime_CR),
         Day = date(tomst$tidydate),
         Hour = hour(tomst$tidydate),
         Month = month(tomst$tidydate),
         Year = year(tomst$tidydate)) %>% 
  filter(Day >= '2021-10-01', Day <='2022-02-02')

# manual VWC curve for both  soil types: y = 0.0002x - 0.1795. R2 = 0.93
# manual gravimetric soil m (GSM) curve for both soil types: y = -2E-07x2 + 0.001x - 0.7437, R2 = 0.93
# these equations are from excel; they don't actually map that well onto the trendline in excel even though
# the equation and trendline go together.

tomst <- tomst %>% 
  mutate(VWC_cm3cm3 = ((0.0002*soilm_count)-0.1795), 
         GSM_gg = ((-0.0000002*(soilm_count^2)) + (0.001*soilm_count) - 0.7437),
         timeofday = if_else(Hour >=7 & Hour <=20, "day", "night"),
         season = if_else(Month %in% c("10","11","12"), "Wet", "Dry"))

hourlytomst <- tomst %>%  
  group_by(Collar, Day, Hour, Month, Year) %>% 
  summarise(soilT = mean(subsoilT, na.rm=TRUE),
            surfaceT = mean(surfaceT, na.rm=TRUE),
            airT = mean(airT, na.rm = TRUE),
            VWC = mean(VWC_cm3cm3, na.rm=TRUE),
            GSM = mean(GSM_gg, na.rm = TRUE))

tomstCO2 <- left_join(cleanCO2, hourlytomst, relationship = "many-to-many")

tomstCO2 <- tomstCO2 %>% 
  mutate(timepoint = case_when(Day== "2021-10-25" | Day =="2021-10-26" |Day =="2021-10-24" ~ "10/25/21",
                               Day == "2021-11-05" | Day == "2021-11-06"  ~ "11/05/21",
                               Day =="2021-11-16" ~ "11/16/21",
                               Day =="2021-12-01" | Day == "2021-12-03" ~ "12/1/21",
                               Day =="2021-12-17" ~ "12/17/21",
                               Day =="2022-01-04" ~ "1/04/22",
                               Day =="2022-01-21" | Day == "2022-01-23" ~ "1/21/22",
                               Day =="2022-02-02" ~ "2/02/22")) 
tomstCO2$timepoint.tidy <- mdy(tomstCO2$timepoint)
### can relate CO2 and soilm, soilT here

tomstCO2 %>% 
  group_by(Treatment, timepoint.tidy) %>% 
  summarise(mean_CO2 = mean(CO2_mgC_m2_h, na.rm=TRUE),
            sd_CO2 = sd(CO2_mgC_m2_h, na.rm=TRUE)) %>% 
  ggplot(aes(x=timepoint.tidy, y=mean_CO2, color=Treatment))+
  geom_point() +
  ylab("mg CO2-C / m2/ h") +
  xlab("Timepoint") +
  geom_errorbar(aes(ymin=mean_CO2-sd_CO2, ymax=mean_CO2+sd_CO2), width=1) +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3")) +
  labs(title = "CO2 by treatment and timepoint", caption = "error bars represent sd")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))
  

tomstCO2 %>% 
  ggplot(aes(x=soilT, y=CO2_mgC_m2_h, color=Treatment))+
  geom_point() +
  ylab("mg CO2-C / m2/ h") +
  xlab("soil temperature C") +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3")) +
  labs(title = "CO2 ~ soil temp", caption = "error bars represent sd")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

tomstCO2 %>% 
  ggplot(aes(x=VWC, y=CO2_mgC_m2_h, color=Treatment))+
  geom_point() +
  ylab("mg CO2-C / m2/ h") +
  xlab("Volumetric water content (cm3/cm3") +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3")) +
  labs(title = "CO2 ~ soil m", caption = "error bars represent sd")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

tomstCO2 %>% # summarising because I joined by hour, which means there are repeat values of CO2 and soilm counts for a given timepoint for each collar
  group_by(Day, Hour, Treatment) %>% 
  summarise(mean_CO2 = mean(CO2_mgC_m2_h),
            mean_GSM = mean(GSM)) %>% 
  ggplot(aes(x=mean_GSM, y=mean_CO2, color=Treatment))+
  geom_point() +
  ylab("mg CO2-C / m2/ h") +
  xlab("Gravimetric soil moisture (g/g)") +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3")) +
  labs(title = "CO2 ~ soil m", caption = "error bars represent sd")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))


tomstCO2 %>% 
  ggplot(aes(x=airT, y=CO2_mgC_m2_h, color=Treatment))+
  geom_point() +
  ylab("mg CO2-C / m2/ h") +
  xlab("air temperature C") +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3")) +
  labs(title = "CO2 ~ air temp", caption = "error bars represent sd")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))



tomstCO2day <- tomstCO2 %>% 
  group_by(Treatment, Block, Plot, Collar, Day, Month, Year, season) %>% 
  summarise(mean_CO2_mgCm2h = mean(CO2_mgC_m2_h),
            mean_soilT = mean(soilT),
            mean_surfaceT = mean(surfaceT),
            mean_airT = mean(airT),
            mean_gsm = mean(GSM),
            mean_VWC = mean(VWC))

tempCO2mod <- lmer(mean_CO2_mgCm2h ~ mean_soilT + mean_surfaceT + mean_airT + mean_gsm + (1|Collar), 
                   data=tomstCO2day)
vif(tempCO2mod)

## surfaceT and airT are multicolinear so just one should be included, if either.

tempCO2mod <- (lmer(mean_CO2_mgCm2h ~ mean_soilT + mean_airT + mean_gsm + (1|Collar), 
             data=tomstCO2day))

vif(tempCO2mod)
summary(tempCO2mod)

### look at TOC and TDN
MBC$Block <- as.factor(MBC$Block)
MBC <- MBC %>% 
  mutate(plotid = paste0(Block, "-", Plot))
MBC %>% 
  ggplot(aes(x=Treatment, y = mgC_gsoil_uf, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Block) +
  labs(title = "TOC~Treatment by block") +
  ylab("TOC mgC/gsoil") +
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))


MBC %>% 
  ggplot(aes(x=Treatment, y = mgN_gsoil_uf, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Block) +
  labs(title = "TDN~Treatment by block") +
  ylab("TDN mgN/gsoil") +
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12)) +
  scale_fill_manual(values=c("control" = "skyblue3", "warming" = "brown3"))
  


### joining CO2 and MBC data

# maybe do dry season CO2 and wet season CO2?
  
  #this is all CO2 averaged across time within each month
CO2sum <- cleanCO2 %>%  
  mutate(plotid = paste0(Block, "-", Plot)) %>% 
  group_by(season, plotid, Treatment) %>% 
  summarise(mean_CO2_umolsm2s = mean(CO2flux_umols_m2_s, na.rm=TRUE), 
            mean_CO2_mg_m2h = mean(CO2_mgC_m2_h, na.rm=TRUE))

write.csv(CO2sum, "CO2summary.csv")


CO2MBC <- left_join(CO2sum, MBC) # this has repeat observations for MBC and MBN for each season; should be subsetted

CO2MBC %>% 
  ggplot(aes(x=MBC, y=mean_CO2_mg_m2h, color=Treatment)) +
  facet_wrap(~season) +
  geom_point() +
  ggtitle("CO2 ~ MBC")

CO2MBC %>% 
  filter(season=="Wet") %>% 
  ggplot(aes(x=MBC, y=mean_CO2_mg_m2h, color=Treatment)) +
  geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Wet season CO2 ~ MBC") +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3"))



CO2MBC %>% 
  ggplot(aes(x=mgC_gsoil_uf, y=mean_CO2_mg_m2h, color=Treatment)) +
  facet_wrap(~season) +
  geom_point() +
  ggtitle("CO2 ~ TOC")


CO2MBC %>% 
  ggplot(aes(x=mgN_gsoil_uf, y=mean_CO2_mg_m2h, color=Treatment)) +
  facet_wrap(~season) +
  geom_point() +
  ggtitle("CO2 ~ TDN")

  CO2MBC %>% 
  filter(season == "Wet") %>% 
  lm(mean_CO2_mg_m2h~Treatment *MBC, data=.) %>% 
  summary()

  ## bringing in CN data
CN$Block <- as.factor(CN$Block)
CN <- CN %>% select(-X)
CO2_MBC_CN <- left_join(CO2MBC, CN)


CO2_MBC_CN %>% 
  ggplot(aes(x=percC, y = mean_CO2_mg_m2h, color = Treatment))+
  geom_point() +
  geom_smooth(method="lm")

wet_data <- CO2_MBC_CN %>% 
  filter(season == "Wet",
         plotid!= "2-17") 

dry_data <-CO2_MBC_CN %>% 
  filter(season == "Dry",
         plotid!= "2-17")
 
summary(lm(mean_CO2_mg_m2h ~ Treatment *percC + , data=CO2_MBC_CN))


summary(lm(mean_CO2_mg_m2h ~ Treatment * mgC_gsoil_uf, data=dry_data))
 
wet_data %>% 
  ggplot(aes(x=MBC, y=mean_CO2_mg_m2h, color = Treatment)) +
  geom_point() +
  geom_smooth(method="lm") +
  scale_color_manual(values=c("control" = "skyblue3", "warming" = "brown3")) 
 ### CO2 analyses
  
CO2day <- tomstCO2 %>% 
  group_by(Block, Plot, Treatment, Day, Month, Year, Collar, season) %>% 
  summarise(mean_CO2 = mean(CO2_mgC_m2_h, na.rm=TRUE))

mod <- lmer(mean_CO2 ~ Treatment + season + Block + (1|Collar), data=CO2day)  


#https://yury-zablotski.netlify.app/post/mixed-effects-models-2/
#Collar can be nested within Block if Block is a random effect, but if we are interested in the effect of Block it should be a fixed effect
summary(mod)

summary(CO2day$mean_CO2)
hist(CO2day$mean_CO2,
     xlab = "CO2",
     main = "CO2 averaged within each observation",
     breaks = sqrt(nrow(CO2day))
)
boxplot(CO2day$mean_CO2)
out <- boxplot.stats(CO2day$mean_CO2)$out
out_ind <-which(CO2day$mean_CO2 %in% c(out))





### calculate percent change 
CO2_MBC_CN <- CO2_MBC_CN %>% 
  mutate(plotid = paste0(Block, "-", Plot),
         inorgN_mgN_gsoil = mgNO3_gsoil+mgNH4_gsoil)

widedat <- CO2_MBC_CN %>% 
  select(season, plotid, Block, Treatment, mean_CO2_mg_m2h, mgC_gsoil_uf, 
         mgN_gsoil_uf, MBC, MBN, mgNH4_gsoil, mgNO3_gsoil, percC, percN, C.N) %>% 
  pivot_wider(id_cols = c(plotid, season, Block),
                names_from = Treatment,
              values_from = c(mean_CO2_mg_m2h, mgC_gsoil_uf, 
                              mgN_gsoil_uf, MBC, MBN, mgNH4_gsoil, mgNO3_gsoil,
                              percC, percN, C.N))


## look at how soil properties influence the response to warming

# what drives magnitude and direction of CO2 response?

d_data <- widedat %>% 
  mutate(dCO2 = ((mean_CO2_mg_m2h_warming - mean_CO2_mg_m2h_control)/mean_CO2_mg_m2h_control)*100,
         dDOC = ((mgC_gsoil_uf_warming - mgC_gsoil_uf_control) / mgC_gsoil_uf_control)*100,
         dTDN = ((mgN_gsoil_uf_warming - mgN_gsoil_uf_control) / mgN_gsoil_uf_control)*100,
         dMBC = ((MBC_warming - MBC_control)/MBC_control)*100,
         dMBN = ((MBN_warming - MBN_control)/MBN_control)*100,
         dNH4 = ((mgNH4_gsoil_warming - mgNH4_gsoil_control)/mgNH4_gsoil_control)*100,
         dNO3 = ((mgNO3_gsoil_warming - mgNO3_gsoil_control) / mgNO3_gsoil_control)*100,
         inorgn_warm = mgNH4_gsoil_warming + mgNO3_gsoil_warming,
         inorgn_control = mgNH4_gsoil_control + mgNO3_gsoil_control,
         dinorgN = ((inorgn_warm - inorgn_control) / inorgn_control)*100,
         dpercC = ((percC_warming - percC_control) / percC_control)*100,
         dpercN = ((percN_warming - percN_control) / percN_control)*100,
         dCN = ((C.N_warming - C.N_control)/C.N_control)*100)
         
d_data %>% 
  ggplot(aes(x=percC_control, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=1) +
  ylab("% change CO2") +
  xlab("% soil C, control") +
  labs(title = "%change CO2 from warmed to control in relation to soil C", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

d_data %>% 
  ggplot(aes(x=dpercC, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=1) +
  geom_vline(xintercept=0, linetype = "dashed", color="black", size=0.5) +
  ylab("% change CO2") +
  xlab("% change soil C") +
  labs(title = "%change CO2 from warmed to control in relation to % change soil C", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

d_data %>% 
  ggplot(aes(x=MBC_control, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=1) +
  labs(title = "%change CO2 from warmed to control in relation to microbial biomass C control", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))


d_data %>% 
  ggplot(aes(x=dMBC, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=1) +
  labs(title = "%change CO2 from warmed to control with %change MBC", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

d_data %>% 
  ggplot(aes(x=dMBN, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=1) +
  labs(title = "%change CO2 from warmed to control with %change MBN", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

d_data %>% 
  ggplot(aes(x=dDOC, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  ylab("% change CO2") +
  xlab("% change dissolved org C") +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype = "dashed", color="black", size=0.5) +
  labs(title = "%change CO2 from warmed to control with %change dissolved org C", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

d_data %>% 
  ggplot(aes(x=mgC_gsoil_uf_control, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  ylab("% change CO2") +
  xlab("dissolved org C mg C/ gsoil, control") +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  labs(title = "%change CO2 from warmed to control with dissolved org C in control", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

d_data %>% 
  ggplot(aes(x=dTDN, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  ylab("% change CO2") +
  xlab("% change total dissolved N") +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype = "dashed", color="black", size=0.5) +
  labs(title = "%change CO2 from warmed to control with %change total dissolved N", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))


d_data %>% 
  ggplot(aes(x=inorgn_control, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  ylab("% change CO2") +
  xlab("total inorg N, mg/g soil") +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  labs(title = "%change CO2 from warmed to control with total inorganic N (control)", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))


d_data %>% 
  ggplot(aes(x=dCN, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  ylab("% change CO2") +
  xlab("% change soil C/N") +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  labs(title = "%change CO2 from warmed to control with % change soil C:N", caption = "positive y = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

### linear regressions for these data

# don't need to do mixed effects because I am looking at one value of CO2 for each season and one value of each soil var

# percent change data is heavily right skewed with negative values, so this creates a problem for transforming it
# https://stats.stackexchange.com/questions/296990/transforming-positively-skewed-data-with-positive-and-negative-values

# could also just do warm/con, with values greater than 1 indicating warm is greater than control

p_data <- widedat %>% 
  mutate(dCO2 = mean_CO2_mg_m2h_warming / mean_CO2_mg_m2h_control,
         dDOC = mgC_gsoil_uf_warming / mgC_gsoil_uf_control,
         dTDN = mgN_gsoil_uf_warming / mgN_gsoil_uf_control,
         dMBC = MBC_warming /MBC_control,
         dMBN = MBN_warming / MBN_control,
         dNH4 = mgNH4_gsoil_warming / mgNH4_gsoil_control,
         dNO3 = mgNO3_gsoil_warming / mgNO3_gsoil_control,
         inorgn_warm = mgNH4_gsoil_warming + mgNO3_gsoil_warming,
         inorgn_control = mgNH4_gsoil_control + mgNO3_gsoil_control,
         dinorgN = inorgn_warm / inorgn_control,
         dpercC = percC_warming / percC_control,
         dpercN = percN_warming / percN_control,
         dCN = C.N_warming /C.N_control
  )

p_data %>% 
  ggplot(aes(x=percC_control, y = dCO2)) +
  geom_point(aes(shape = Block), size=3) +
  facet_wrap(~season) +
  geom_hline(yintercept=1, linetype="dashed", color="black", size=1) +
  ylab("ratio warm/control CO2") +
  xlab("% soil C, control") +
  labs(title = "warm/con CO2 in relation to soil C", caption = "y>1 = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

d_data %>% 
  filter(season == "Wet") %>% 
  ggplot(aes(x=percC_control, y = dpercC)) +
  geom_point(aes(shape = Block), size=3) +
  geom_hline(yintercept=1, linetype="dashed", color="black", size=1) +
  ylab("% change in soil C") +
  xlab("% soil C, control") +
  labs(title = "changes in soil C in relation to soil C in the control", caption = "y>1 = warm > con")+
  theme(axis.text = element_text(size=12), strip.text = element_text(size=12))

summary(lm(dpercC ~ percC_control, data=p_data))

### proportions (warm/control) right skewed and can easily be fixed with log transformation because no neg. valeus

# models

## CO2 ~ microclimate

tempCO2mod <- lmer(mean_CO2_mgCm2h ~ mean_soilT + mean_gsm + (1|Collar), 
                   data=tomstCO2day)
summary(tempCO2mod)

tomstCO2day %>% 
  ggplot(aes(x=mean_gsm, y=mean_CO2_mgCm2h, color=as.factor(Month))) +
  geom_point()

## CO2 ~ C pools (MBC, DOC)

summary(lm(mean_CO2_mg_m2h ~ MBC  +mgC_gsoil_uf, data=CO2_MBC_CN))
# DOC best predictor of CO2
## CO2 ~ N pools (MBN, TDN, inorgN)

summary(lm(mean_CO2_mg_m2h ~ percN +inorgN_mgN_gsoil + MBN +mgN_gsoil_uf + Treatment ,data=CO2_MBC_CN))
# TDN best predictor of CO2

# but when all vars in the model, nothing emerges as significant

# hypo-driven model:

# CO2 production should be directly related to how much microbial biomass there is, how much
  # C substrate (and perhaps N) there is for the microbial biomass to utilize, and the abiotic conditions that
  # may mediate microbial metabolism.

# so perhaps a new calculated value could represent substrate per unit MBC? 

wet_data <- CO2_MBC_CN %>% 
  filter(plotid != "2-17") %>% #because there is the weird DOC outlier
  filter(season == "Wet") %>% 
  mutate(DOC_per_MBC = mgC_gsoil_uf/MBC, #substrate C per unit MBC gives unit-less index
         DCN = mgC_gsoil_uf/mgN_gsoil_uf, #CN of K2SO4-extractable substrate
         MBCN = MBC/MBN, #CN of microbial biomass
         CN_sub_by_MB = DCN/MBCN) #an index of stoichiometry of substrate divided by microbial biomass

summary(lm(mean_CO2_mg_m2h~DOC_per_MBC * Treatment, data=wet_data)) 


summary(lm(mean_CO2_mg_m2h ~  MBC + mgC_gsoil_uf * Treatment, data=wet_data)) 

# trying to capture size of microbial biomass, the CN of the dissolved substrate ("quality")

# MBC + DCN + Treatment, only MBC significant
# MBC and DOC are fairly correlated (vif = 2.7)
# maybe Co2 per unit biomass ~ soil vars? i.e. how much each microbe respires perhaps depends on substrate?

## should I do a PCA with all soil vars?

CO2_MBC_CN %>% 
  ggplot(aes(x=mgC_gsoil_uf, y=mean_CO2_mg_m2h, color=Treatment))+
  geom_point() +
  theme_bw()


# 
# ### hypo-driven model:
# 
# * CO2 production should be directly related to how much microbial biomass there is, how much C substrate (and perhaps N) there is for the microbial biomass to utilize, and the abiotic conditions that may mediate microbial metabolism. 
# 
# * How to structure the models?
#   
#   + CO2 as the response
# + dCO2 as the response (warm/control, e.g.)