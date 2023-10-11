## ventisqeros CO2 loop

library(lubridate)
library(dplyr)
library(readr)


d <- read_tsv("~./all CO2 data/LI-COR2-241021_offset.txt")

folder <- list.files(path = "/Users/yessicamurray/Box/J Murray/Projects/Paramo 2021/datos!/all CO2 data",  
                        # Now follows a regular expression that matches:
                        pattern = "_offset.txt",
                        #          |            |        the standard file extension of comma-separated values
                        #          |            the variable parts (two digits, each between 0 and 9)
                        #          the static part of the filenames
                        full.names = TRUE)


  CO2list = list()
  
  for (i in 1:length(folder)) {
    co2dat <- read_tsv(file=folder[i])
    co2dat$path = folder[i]
    co2dat <- co2dat %>% 
      select(Label, `File Name`, Date_IV, `Obs#`, Offset, Area, `Instrument Name`, `Pre-purge`, `Dead Band`, `Observation Length`, `Post-purge`, Exp_Flux, Exp_R2, Tcham_IV, Pressure_IV, RH_IV, path)
    co2dat$Date <- as.character(co2dat$Date_IV)
    co2dat$Date <- ymd_hms(co2dat$Date)
    CO2list[[i]] <- co2dat
        }
  CO2_data = do.call(rbind, CO2list)  

# there are some NAs in the Date column because not all were formatted the same (some ymd, some mdy)
  # pull out the NAs
  nalist <- CO2_data %>% 
    filter(is.na(Date))
  #format them
  nalist$Date <- mdy_hm(nalist$Date_IV)
  
  #pull out the ones without NAs
  datelist <- CO2_data %>% 
    filter(!is.na(Date))
  
#bind the two dataframes with correctly formatted dates together
  CO2 <- bind_rows(datelist, nalist)
  
  CO2 <- CO2 %>%
    rename(Collar =Label)
  
  
  
  metadata <- read.csv("ventisqueros CO2 metadata.csv")
  
  CO2 <- inner_join(CO2, metadata, by="Collar")
  
write.csv(x=CO2, "CO2_master_oct2023.csv")
  