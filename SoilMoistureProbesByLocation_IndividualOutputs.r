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

Leap      <- T     #Leap year? (T/F)
OutputTag <- "Out" #To name our output file

for(i in 1:NumSites){  #Iterate over all possible university sites
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
      NewFile[,1:2] <- YearMonthDays[1:RepeatedChunk,1:2]
      NewFile[,"DayOfYear"] <- as.character(DATA[1:RepeatedChunk,"Time"])
      for(i in 1:(RepeatedChunk)){
         for(j in 1:U){
            NewFile[i,3+j] <- SMReadings[(j-1)*RepeatedChunk+i]
         }
      }
      NewFile <- NewFile[-1,]
      
      #Stitch the sensor info
      NewColumns <- matrix(NA,nrow(NewFile),3)
      for(i in 1:nrow(SensorInfo)){
        for(j in 1:ncol(SensorInfo)){
           NewColumns[i,j+1] <- as.character(SensorInfo[i,j])
        }
      }
      NewColumns[,1] <- " "
      NewFile <- cbind(NewFile,NewColumns)
      NewFile[is.na(NewFile)] <- " "
      #Write the re-formatted file to the same directory
      write.csv(NewFile,paste(Path,OutputTag,Name,".csv",sep=""))
  
    }
}

