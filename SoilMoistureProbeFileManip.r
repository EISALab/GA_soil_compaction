rm(list=ls())

#Put necessary functions in memory (You'll need to change the directory to
#wherever these soil probe functions are going to be stored on the host
#computer)
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SoilProbeFunctions.r")

#Read soil sensor file
#Name <- "SensorsEC"  #this is the input file name (must be saved as .csv)

OutputTag <- "Out" #how shall we tag the output file?
Path <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SoilMoistureProbeFiles\\"

#Read in the list of files that require conversion
FNameDir <- "SoilProbeNames"
FileNameList <- GetFile(FNameDir,Path)

    
for(i in 1:nrow(FileNameList)){
    Name <- as.character(FileNameList[i,"FileName"]) 
    print(Name)    
    DATA <- GetFile(Name,Path) #GetFile is a function in SoilProbeFunctions
    Default <- "SoilMoisture" #What is the label for MOST of the rows?
      
    #Grab the sensor information, then hang onto it...we'll need it later
    SensorInfo <- DATA[4:13,(ncol(DATA)-1):(ncol(DATA))]
    
    #Chop down to size
    DATA <- DATA[,c("Type","Time","Calibration")]
    
    #Grab the appropriate "depths" for each row
    Depths <- GetDepths(DATA,Default)  #From SoilProbeFunctions
    
    #Get year and day
    YearMonthDays <- ProcessTimes(DATA) #From SoilProbeFunctions
    
    #Grab the readings column
    SMReadings <- as.numeric(DATA[,"Calibration"])
    
    #Build the final database
    U <- length(unique(Depths)) #How many depths does this sensor track?
    L <- length(Depths)         #How many total rows?
    RepeatedChunk <- L/U        #Assumes each sensor has the same number of readings
    NewFile <- matrix(NA,RepeatedChunk,U+3)
    colnames(NewFile) <- c(colnames(YearMonthDays),unique(Depths))
    NewFile[,1:3] <- YearMonthDays[1:RepeatedChunk,]
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
