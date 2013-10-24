#Functions for soil moisture probe data file reformatting
#Evan Coopersmith, JDTIC, July 2012

#Read in a file, given a directory
GetFile <- function(Name,Path){
   filename <- paste(Name,".csv",sep="")
   FullPath <- paste(Path,filename,sep="")
   Data <- read.csv(FullPath)
   Data
}

#Takes in a column that divides "soil moisture" rows by "sensor 10" or
#"sensor 20" etc, then processes them into appropriate row names
GetDepths <- function(DATA,Default){
  #Format to be addressed as the proper sensor depth
  Depths <- as.character(DATA[,"Type"])
  ValidName <- Depths[1] #This should be something like "Sensor 10"
  for(i in 1:length(Depths)){
     if(Depths[i]==Default){
       Depths[i] <- ValidName
     }else{
       ValidName <- Depths[i]
     }
  }
  Depths
}

#Reads in the time column " 06-03-2012 13:22:00 " and returns:
#Year = 2012, "Day of Year = 154.51(e.g.)
ProcessTimes <- function(DATA,Leap){
   L <- nrow(DATA)
   DaysInMonths <- c(31,28,31,30,31,30,31,31,30,31,30,31)
   if(Leap){ c(31,29,31,30,31,30,31,31,30,31,30,31) }
   CumSumMonths <- cumsum(DaysInMonths)
   Times <- as.character(DATA[,"Time"])
   MyDatesAndTimes <- matrix(NA,L,3)
   for(i in 1:L){
        SplitData <- strsplit(Times[i]," ")
     if(length(SplitData[[1]])>1){
        SplitData <- c(SplitData[[1]][1],SplitData[[1]][2])
        MDY   <- strsplit(SplitData[1],"-")
        Month <- as.numeric(MDY[[1]][1])
        Day   <- as.numeric(MDY[[1]][2])
        Year  <- as.numeric(MDY[[1]][3])
        if(Month > 1){
           DOY <- CumSumMonths[Month-1] + Day
        }else{
           DOY <- Day
        }
        HourMinSec <- strsplit(SplitData[2],":")
        Hour <- as.numeric(HourMinSec[[1]][1])
        Min  <- as.numeric(HourMinSec[[1]][2])
        Sec  <- as.numeric(HourMinSec[[1]][3])
        TOD <- Hour/24 + Min/24/60 + Sec/24/60/60
        DecimalDate <- DOY + TOD
        MyDatesAndTimes[i,] <- c(Year,Month,DecimalDate)
      }
   }
   colnames(MyDatesAndTimes) <- c("Year","Month","DayOfYear")
   MyDatesAndTimes
}

#Reads in the time column " 06-03-2012 13:22:00 " and returns:
#Year = 2012, "Day of Year = 154.51(e.g.)
ProcessTimes_v2 <- function(DATA,Leap){
   L <- nrow(DATA)
   DaysInMonths <- c(31,28,31,30,31,30,31,31,30,31,30,31)
   if(Leap){ c(31,29,31,30,31,30,31,31,30,31,30,31) }
   CumSumMonths <- cumsum(DaysInMonths)
   Times <- as.character(DATA[,"Time"])
   MyDatesAndTimes <- matrix(NA,L,3)
   for(i in 1:L){
        SplitData <- strsplit(Times[i]," ")
     if(length(SplitData[[1]])>1){
        SplitData <- c(SplitData[[1]][1],SplitData[[1]][2])
        if(substring(SplitData[1],first=3,last=3)=="-"){
          SplitChar = "-"
        }else{
          SplitChar = "/"
        }
        MDY   <- strsplit(SplitData[1],SplitChar)
        Month <- as.numeric(MDY[[1]][1])
        Day   <- as.numeric(MDY[[1]][2])
        Year  <- as.numeric(MDY[[1]][3])
        if(Month > 1){
           DOY <- CumSumMonths[Month-1] + Day
        }else{
           DOY <- Day
        }
        HourMinSec <- strsplit(SplitData[2],":")
        Hour <- as.numeric(HourMinSec[[1]][1])
        Min  <- as.numeric(HourMinSec[[1]][2])
        TOD <- Hour/24 + Min/24/60
        DecimalDate <- DOY + TOD
        MyDatesAndTimes[i,] <- c(Year,Month,DecimalDate)
      }
   }
   colnames(MyDatesAndTimes) <- c("Year","Month","DayOfYear")
   MyDatesAndTimes
}


#Given a 0-364.9999 date, return the month (1-12)
ConvertDateToMonth <- function(Date, CMonths){
   M <- min(which(CMonths >= Date))
   M
}

#Given a date of of the format "MM/DD/YYYY HH:MM", return MM/DD HH:MM AM
FormatDay <- function(DATA,Leap){

}