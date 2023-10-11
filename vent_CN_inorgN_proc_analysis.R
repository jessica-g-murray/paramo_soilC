## ventisqueros CN and inorg N data processing and analysis

data <- read.csv("CN_raw_all.csv")
head(data)


metadata <- read.csv("ventisqueros cn metadata.csv")
head(metadata)

CN <- left_join(data, metadata)
head(CN)

CN <- CN %>% 
  select(-ID.LAB) %>% 
  rename(Collar = ID.USUARIO,
         percC = C,
         percN = N)

dat <- read.csv("inorgN_raw_all.csv")
dat <- dat %>% 
  rename(NH4 = N.NH4.,
         NO3 = N.NO3.) %>% 
  select(Collar, NH4, NO3)

metadat <- read.csv("ventisqueros inorgN metadata.csv")
inorg <- full_join(dat, metadat)
head(inorg)

inorg <- inorg %>% 
  mutate(soilm_dry_g = soilm_dry_tin - soilm_tin, #dry soil from soil m
         soilm_grav = ((soilm_wet_g - soilm_dry_g) / soilm_dry_g), #gravimetric soil m
         drysoil_sample_g = (freshsoil_g/(1+soilm_grav)), #dry soil that was extracted
         liquid_sample_L = (((freshsoil_g-drysoil_sample_g)*0.001) + KCL_L), #total solution associated with sample: K2SO4 and also soil water
         mgNH4_gsoil = (NH4 * liquid_sample_L / drysoil_sample_g), #calc mg/g soil
         mgNO3_gsoil = (NO3 * liquid_sample_L / drysoil_sample_g),
         mgNH4_kgsoil = mgNH4_gsoil*1000,
         mgNO3_kgsoil = mgNO3_gsoil*1000
  )

CN_inorg <- left_join(inorg, CN)
write.csv(CN_inorg, "CN_inorgNmaster.csv")

