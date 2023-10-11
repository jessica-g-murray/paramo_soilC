### ventisqueros MBC and TOC analysis

library(ggplot2)

### going to read in raw data and do calculations here

# Drift: I'm notconvinced that the variation in the check standards is "drift" - it looks more like instrument error
  # to me because the variation isn't directional across the run and it changes the values a good bit. I'm going to
  # do analyses on the original values (i.e. not drift corrected) but can change that later if needed.

folder <- list.files(path = "/Users/yessicamurray/Box/J Murray/Projects/Paramo 2021/Analyses/Paramo C cycling warming experiment/all TOC data",  
                     # Now follows a regular expression that matches:
                     pattern = ".csv",
                     #          |            |        the standard file extension of comma-separated values
                     #          |            the variable parts (two digits, each between 0 and 9)
                     #          the static part of the filenames
                     full.names = TRUE)
TOClist = list()

for (i in 1:length(folder)) {
  d <- read.csv(file=folder[i], header=FALSE, stringsAsFactors = FALSE)
  tmp <- d[-c(1:9),]
  tmp = tmp %>%
    mutate(row_num = as.numeric(row.names(.)))
    colnames(tmp) = tmp[1,]
  tmp <- tmp[-1,]

  data <- tmp %>% 
    select(`Sample Name`, `Result(NPOC)`, `Result(TN)`) %>% 
    rename(sample_id = `Sample Name`) %>% 
    filter(!grepl("CHECK|BKL|STD|Untitled", sample_id))
  
  TOClist[[i]] <- data
}

tocdat = do.call(rbind, TOClist)  

metadata <- read.csv("vent_toc_metadata.csv")

toc_meta <- left_join(tocdat, metadata)

## skipping blank correction for now



#ventisqueros TOC analyses, MBC calcs

names(toc_meta)
toc_meta <- toc_meta %>% 
  rename(NPOC_mgL = `Result(NPOC)`,
         TDN_mgL = `Result(TN)`) %>% 
  mutate_at(c('NPOC_mgL', 'TDN_mgL'), as.numeric)

TOC <- toc_meta %>% 
  mutate(soilm_dry_g = soilm_dry_tin - soilm_tin, #dry soil from soil m
         soilm_grav = ((soilm_wet_g - soilm_dry_g) / soilm_dry_g), #gravimetric soil m
         drysoil_sample_g = (freshsoil_g/(1+soilm_grav)), #dry soil that was extracted
         liquid_sample_L = (((freshsoil_g-drysoil_sample_g)*0.001) + k2so4_L), #total solution associated with sample: K2SO4 and also soil water
         mgC_gsoil = (NPOC_mgL * liquid_sample_L / drysoil_sample_g), #calc mg/g soil
         mgN_gsoil = (TDN_mgL * liquid_sample_L / drysoil_sample_g)
  )

TOC <- TOC %>%
  na.omit()
## change to wide format to calculate MBC (want separate columns for FUM and UF)
wideTOC <- TOC %>% 
  select(Collar, Block, Plot, Treatment,FUM, mgC_gsoil, mgN_gsoil) %>% 
  pivot_wider(names_from = FUM, values_from = c(mgC_gsoil, mgN_gsoil))

wideTOC <- wideTOC %>% 
  mutate(MBC = (mgC_gsoil_fum/0.45) - mgC_gsoil_uf, #0.45 correction
         MBN = (mgN_gsoil_fum/0.45) - mgN_gsoil_uf)

wideTOC %>% 
ggplot(aes(x=Treatment, y = MBC, fill=Treatment)) +
  geom_boxplot()

wideTOC %>% 
  ggplot(aes(x=Treatment, y = MBN, fill=Treatment)) +
  geom_boxplot()

write.csv(wideTOC, "MBC.csv")

