rm(list=ls())

#Put necessary functions in memory (You'll need to change the directory to
#wherever these soil probe functions are going to be stored on the host
#computer)
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SoilProbeFunctions.r")

#Define the higher directory
UpperDir <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\CropSense Files 10-3\\"

#Put list of sites into memory - this will tell us in what folders to find the relevant files
SiteListName <- "ListOfSites"
SiteList <- GetFile(SiteListName,UpperDir)

#Put list of probes by site into memory - this will tell us in what folders to find the relevant files
ListBySiteName <- "ListBySite"
ListBySite <- GetFile(ListBySiteName,UpperDir)

NumSites <- length(as.character(as.list(SiteList)[[1]]))  #How many sites?

Year      <- 2012  #Cosmetic in our output file
StartDate <- 150   #Starting day of year (0-365)
Interval  <- 30    #How many minutes between samples?
EndDate   <- 280
Leap      <- T     #Leap year? (T/F)
#Build template, to be filled and data-merged as we read individual probe files
LTemp <- (EndDate-StartDate)*24*60/Interval
MonthCuts <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
   if(Leap){ MonthCuts[2] <- MonthCuts[2] + 1 }
   CMonths <- cumsum(MonthCuts) #Cumulative summation

for(i in 1:NumSites){  #Iterate over all possible university sites

  SiteArray <- matrix(NA,LTemp,3)
  colnames(SiteArray) <- c("Year", "Month", "DayOfYear")
  #Fill the template
  for(x in 1:LTemp){
     CurrentTime <- StartDate + (x-1)*Interval/60/24
     SiteArray[x,1] <- Year
     SiteArray[x,2] <- ConvertDateToMonth(CurrentTime, CMonths)
     SiteArray[x,3] <- CurrentTime
  }

  SiteName <- as.character(as.list(SiteList)[[1]])[i] 
  print(paste("Currently analyzing",SiteName))
  ProbeList <- as.character(ListBySite[,SiteName])
  ProbeList <- ProbeList[ProbeList != "XXXX"]
  for(j in 1:length(ProbeList)){ #Iterate over all probes at that site 
      Name <- as.character(ProbeList[j]) 
      print(Name)    
      Path <- paste(UpperDir,SiteName,"\\",sep="")
      DATA <- GetFile(Name,Path) #GetFile is a function in SoilProbeFunctions
      Default <- "SoilMoisture" #What is the label for MOST of the rows?
      #Grab the sensor information, then hang onto it...we'll need it later
      SensorInfo <- DATA[4:13,(ncol(DATA)-1):(ncol(DATA))]
        #Store the sensor code, and we'll assume the rest is unnecessary for now
        ProbeCode <- as.character(SensorInfo[,2])[3]
      #Chop down to size
      DATA <- DATA[,c("Type","Time","Calibration")]
    
      #Grab the appropriate "depths" for each row
      Depths <- GetDepths(DATA,Default)  #From SoilProbeFunctions
    
      #Get year and day
      YearMonthDays <- ProcessTimes_v2(DATA,Leap) #From SoilProbeFunctions
    
      #Grab the readings column
      SMReadings <- as.numeric(DATA[,"Calibration"])
      
      #Build the final database
      U <- length(unique(Depths)) #How many depths does this sensor track?
      L <- length(Depths)         #How many total rows?
      RepeatedChunk <- L/U        #Assumes each sensor has the same number of readings
      NewFile <- matrix(NA,RepeatedChunk,U+3)
      colnames(NewFile) <- c(colnames(YearMonthDays),unique(Depths))
      NewFile[,1:3] <- YearMonthDays[1:RepeatedChunk,]
      for(m in 1:(RepeatedChunk)){
         for(n in 1:U){
            NewFile[m,3+n] <- SMReadings[(n-1)*RepeatedChunk+m]
         }
      }
      NewFile <- NewFile[-1,]
      
      #Add probe code to the sensor columns...
      colnames(NewFile)[4:length(colnames(NewFile))] <- 
         paste(colnames(NewFile)[4:length(colnames(NewFile))],ProbeCode)
         
      #Attach the "NewFile" to the existing template
      NewColumns <- ncol(NewFile)-3
      SiteArray <- cbind(SiteArray, matrix(NA,nrow(SiteArray),NewColumns))
      colnames(SiteArray)[(ncol(SiteArray)-NewColumns+1):ncol(SiteArray)] <-
         colnames(NewFile)[4:ncol(NewFile)]
         
         #Line up new file with the existing time stamps (missing stamps will
         #(receive the last valid reading)  
         print("Merging") 
         for(k in 1:nrow(SiteArray)){
             Now <- SiteArray[k,"DayOfYear"]
            if(Now > max(NewFile[,"DayOfYear"])){ #If we're after the end of the file
               Index <- nrow(NewFile)
            }else{
               Index <- min(which(NewFile[,"DayOfYear"]>=Now)) #First index after
            }
            SiteArray[k,(ncol(SiteArray)-NewColumns+1):ncol(SiteArray)] <-
               NewFile[Index,4:ncol(NewFile)]
         }

  }
  write.csv(SiteArray,paste(UpperDir,SiteName,"Combined",".csv",sep=""))
}

