FileName <- "IllinoisClimateNetworkDaily_cmi_EC"
Variables <- c("max_wind__gust","avg_wind_speed","avg_wind_dir","sol_rad","max_air_temp","min_air_temp","avg_air_temp",
                  "max_rel_hum","min_rel_hum","avg_rel_hum","avg_dewpt_temp","precip","pot_evapot","max_soiltemp_4in","min_soiltemp_4in",
                  "avg_soiltemp_4in","max_soiltemp_8in","min_soiltemp_8in","avg_soiltemp_8in","DayOfYear")


#Read in the data
   GetData <- function(FileName){
     DATA <- read.csv(paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\",FileName,".csv",sep=""))
     #Get rid of rows without a date
     DATA <- DATA[!is.na(DATA[,"Date"]),]
     DATA
   }

DayOfYearSorter <- function(DATA,Variables){   
   ###This code will parse the data into mean values for each variable for each day of the year
   DayOfYear <- DATA[,"Date"] %% 10000
   DATA <- cbind(DATA,DayOfYear)
   ListOfDays <- unique(DayOfYear)
   DailyAverages <- matrix(NA,length(ListOfDays),length(Variables))
   for (i in 1:length(ListOfDays)){
      DailySubSet <- DATA[DATA[,"DayOfYear"]==ListOfDays[i],]
      for j in 1:length(Variables)){
         DailyAverages[i,j] <- mean(DailySubSet[,Variables[j]])
      }
   }
   colnames(DailyAverages) <- Variables
   DailyAverages
}

SerialCorrChecker <- function(DATA,NumDays,Variables){
  CorrArray <- matrix(NA,NumDays,length(Variables))

  for (i in 1:(NumDays)){
     for (j in 1:length(Variables)){
        Current <- DATA[1:(nrow(DATA)-i),Variables[j]]
        Past <- DATA[(i+1):nrow(DATA),Variables[j]]
        CorrArray[i,j] <- cor(Current,Past)
     }
  }
  colnames(CorrArray) <- Variables
  CorrArray
}

OutputName <- "SerialCorr"
write.csv(CorrArray,paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\",OutputName,".csv",sep=""))
 
 
 PredictiveCorrChecker <- function(DATA,NumDays,Variables){
   for(i in 1:NumDays){
      print(paste("i =",i))
      PredDays <- matrix(NA,length(Variables),length(Variables))
      colnames(PredDays) <- Variables
      rownames(PredDays) <- Variables
      for(j in 1:length(Variables)){
         print(paste("j =",j))
         for(k in 1:length(Variables)){
            print(paste("k =",k))
            Current <- DATA[1:(nrow(DATA)-i),Variables[j]]
            Past    <- DATA[(i+1):nrow(DATA),Variables[k]]
            PredDays[j,k] <- cor(Current,Past)
         }
      }
   OutputName <- paste("PredCorr",i,sep="_")
   write.csv(PredDays,paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\",OutputName,".csv",sep="")) 
   }
 }
 
    
 
  
