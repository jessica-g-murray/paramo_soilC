## analysing TOMST Data

#Semicolon separated
#0;31.10.2013 11:45;0;21.5625;22.0625;23.125;148;1;0

# 0	index of the measure
# 31.10.2013 11:45	date and time in UTC
# 0	time zone
# 21.5625	T1
# 22.0625	T2
# 23.125	T3
# 148	soil moisture count (raw moisture data)
# 1	shake – values 1 (shake) or 0. Only for old model TMS-3.
# 0	errFlag (if =1 the device couldn’t convert time from PCF chip)

library(lubridate)
library(dplyr)


setwd("~/Box/J Murray/projects/Paramo 2021/datos!/TOMST data/june 2022/")

## sensor name is in the file name, so I'm going to need to 
    # employ a loop similar to that used for Idgie's files
# https://swcarpentry.github.io/r-novice-inflammation/03-loops-R/

##create an object containing a list of all the TOMST files you want to extract data from
filenames <- list.files(path = "~/Box/J Murray/projects/Paramo 2021/datos!/TOMST data/june 2022",  
                        # Now follows a regular expression that matches:
                        pattern = "data_942[0-1]6[1-6][0-9]{2}_[0-4].csv",
                        #          |            |        the standard file extension of comma-separated values
                        #          |            the variable parts (two digits, each between 0 and 9)
                        #          the static part of the filenames
                        full.names = TRUE)




datalist = list()

for (i in 1:length(filenames)) {
  # ... make some data
  dat <- read.csv2(file = filenames[i], sep=';', header=FALSE)
  header <- c("index", "datetime",'timezone', 'T1', 'T2', 'T3', 'soilm_count')
  names(dat) <- header
  dat <- dat[-c(8:10)]
  dat$path = filenames[i]
  dat$filename = sapply(as.character(dat$path), function(y) strsplit(y, "/")[[1]][11]) #makes a column with the filename
  dat$sensorID = sapply(as.character(dat$filename), function(y) strsplit(y, "_")[[1]][2]) #makes a column with the sensor ID
  datalist[[i]] <- dat # add it to your list
}



big_data = do.call(rbind, datalist)


big_data$date <- strptime(big_data$datetime, "%Y.%m.%d %H:%M") 


          
setwd("~/Box/J Murray/Projects/Paramo 2021/Analyses/Paramo C cycling warming experiment")
metadata <- read.csv("ventisqueros CO2 metadata.csv")

tomstdata <- merge(big_data, metadata, by="sensorID")

tomstdata <- tomstdata %>%
  rename(subsoilT =T1,
         surfaceT =T2,
         airT =T3)



## Let's fix TZ
# time zone
tomstdata$datetime_px_tz<-strptime(tomstdata$datetime, "%Y.%m.%d %H:%M", tz="GMT")
      # tell it the time zone

tomstdata$datetime_CR <-with_tz(tomstdata$datetime_px_tz, tzone = "Etc/GMT+6")
      #confusingly you need to put +6 to get GMT -6 because it is inversed. see https://en.wikipedia.org/wiki/List_of_tz_database_time_zones

write.csv(x=tomstdata, "~/Box/J Murray/Projects/Paramo 2021/Analyses/Paramo C cycling warming experiment/TOMST_master_oct20_june22.csv")







