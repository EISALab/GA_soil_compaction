####Given the name of a soil moisture label, return its longitude and latitude
GetSoilCoordinates <- function(MoistureLabel,SoilDetails){
   ParsedLabel <- strsplit(MoistureLabel,"_")
   if(ParsedLabel[[1]][2]=="Avg"){
      FileName <- paste("Complete",ParsedLabel[[1]][3],sep="_")
   }else{
      for(i in 1:nrow(SoilDetails)){
         SplitCheckedName <- strsplit(as.character(SoilDetails[i,"FileName"]),"_")
         if(length(SplitCheckedName[[1]])>2){
             if(SplitCheckedName[[1]][3]==ParsedLabel[[1]][2]){
                FileName <- paste(SplitCheckedName[[1]][1],"SoilMoisture",SplitCheckedName[[1]][3],sep="_")
             }
         }
      }
   }
   Index     <- which(SoilDetails[,"FileName"]==FileName)
   LatLong   <- SoilDetails[Index,c("Lat","Long")]
   LatLong
}

####Given an AllSoilData matrix, containing all of the soil moisture data to which
####we have access, create a matrix which lists the latitude and longitudes of the soil
####locations in the appropriate order
SoilLatLongs <- function(AllSoilData,SoilDetails){
  LatLongs  <- matrix(NA,(ncol(AllSoilData)-1),2)
  for(i in 1:nrow(LatLongs)){
     MoistureLabel      <- colnames(AllSoilData)[i+1]
     SoilLatLong        <- GetSoilCoordinates(MoistureLabel,SoilDetails)
     LatLongs[i,1]      <- as.numeric(SoilLatLong["Lat"])
     LatLongs[i,2]      <- as.numeric(SoilLatLong["Long"])
  }
  colnames(LatLongs) <- c("Lat","Long")
  LatLongs
}

###Add soil moisture estimates to the soil compaction data
###Note: this will add the profile at the depth selected when AllSoilDATA is called
###Each time this function is called, a new column is added to our dataset
###for the profile desired along with two other columns for a 2-day history.
AddSoilMoisture <- function(SoilCompactionDATA, AllSoilDATA, SoilDetails, Name){
  SoilMoistureData   <- array(NA,nrow(SoilCompactionDATA))
  SoilMoistureData_1 <- array(NA,nrow(SoilCompactionDATA))
  SoilMoistureData_2 <- array(NA,nrow(SoilCompactionDATA))
  LatLongs <- SoilLatLongs(AllSoilDATA,SoilDetails)
  for(i in 1:nrow(SoilCompactionDATA)){
     Lat        <- SoilCompactionDATA[i,"Lat"]
     Long       <- SoilCompactionDATA[i,"Long"]
     SQ_Dists   <- array(9999,nrow(LatLongs))
     Date       <- SoilCompactionDATA[i,"Date"]
     DateIndex  <- which(AllSoilDATA[,"Date"]==Date)
     SensorVals <- as.numeric(AllSoilDATA[DateIndex,2:ncol(AllSoilDATA)])
     for(j in 1:nrow(LatLongs)){
        if(!(is.na(SensorVals[j]))){
          SQ_Dists[j] <- (Lat-as.numeric(LatLongs[j,"Lat"]))^2 + (Long-as.numeric(LatLongs[j,"Long"]))^2
        }
     }
     DistanceWeights     <- 1/SQ_Dists
     NormWeights         <- DistanceWeights/sum(DistanceWeights)
     if(length(DateIndex)==1){                                     
        SoilMoistureEst     <- sum(NormWeights*SensorVals,na.rm=T)
     }else{
        SoilMoistureEst     <- NA 
     }
     SoilMoistureData[i] <- SoilMoistureEst
     if(i < nrow(SoilCompactionDATA)){
        SoilMoistureData_1[i+1] <- SoilMoistureEst
     }
     if(i+1 < nrow(SoilCompactionDATA)){
        SoilMoistureData_2[i+2] <- SoilMoistureEst
     }
  }
  SMDATA             <- cbind(SoilMoistureData,SoilMoistureData_1,SoilMoistureData_2)
  colnames(SMDATA)   <- c(Name,paste(Name,"1",sep="_"),paste(Name,"2",sep="_"))
  SoilCompactionDATA <- cbind(SoilCompactionDATA,SMDATA)
  ###Cut the data after which there are no more soil moisture readings!
  SoilCompactionDATA <- SoilCompactionDATA[!is.na(SoilCompactionDATA[,Name]),]
  SoilCompactionDATA
}


####Create a covariance matrix between soil moisture readings, being careful to eliminate points
####in which one or more of the pairs are missing
####ERROR HANDLING IS REQUIRED FOR LARGE VALUES...CORRELATIONS CAN APPROACH 1 OTHERWISE!
####This should be done in the GetAllSoilData function
SoilMoistureCorMat <- function(AllSoilDATA){
   CorMat      <- matrix(NA,ncol(AllSoilDATA)-1,ncol(AllSoilDATA)-1)
   CorCountMat <- matrix(NA,nrow(CorMat),nrow(CorMat))
   for(i in 1:nrow(CorMat)){
      for(j in 1:nrow(CorMat)){
         Vector1 <- AllSoilDATA[,i+1]
         Vector2 <- AllSoilDATA[,j+1]
         CorMat[i,j] = round(cor(Vector1,Vector2,use="complete.obs"),3)
      }
   }
   colnames(CorMat) <- colnames(AllSoilDATA)[2:ncol(AllSoilDATA)]
   rownames(CorMat) <- colnames(AllSoilDATA)[2:ncol(AllSoilDATA)]
   CorMat
}

####Read in all available soil moisture data, line up all present daily data.
####Cut at specific times of day
GetAllSoilData <- function(TimeOfDay,MoistureReading,BeginDate,EndDate,SoilCompactionDATA,SoilDetails){

    SoilDATA              <- matrix(NA,(EndDate-BeginDate+1),2)
    AllSoilDATA           <- matrix(NA,(EndDate-BeginDate+1),2)
    colnames(SoilDATA)    <- c("Date","_")
    colnames(AllSoilDATA) <- c("Date","_")
    SoilDATA[,"Date"]     <- seq(BeginDate,EndDate)
    AllSoilDATA[,"Date"]  <- seq(BeginDate,EndDate)


    for(i in 1:nrow(SoilDetails)){
       FileName      <- as.character(SoilDetails[i,"FileName"])
         if(FileName != "NO_FILE"){          ###Make sure there is a file to open
             SoilFile      <- GetEBIData(FileName)
             MinDate       <- 40298 #floor(min(SoilFile[,"TIMESTAMP"]))
             MaxDate       <- floor(max(SoilFile[,"TIMESTAMP"]))
             if(MaxDate > MinDate){
               Code          <- strsplit(FileName,"_")[[1]][3]
               #Now, slice the dataset so as only to include datapoints at the appropriate time
               #We will need to iterate over the dates of the soil file to determine Jordan's key times
               #Check if the relevant date was one where Jordan took data...
               TimeStamps <- array(0,(MaxDate-MinDate+1))
               for(j in 1:length(TimeStamps)){
                 DateIndex <- which(SoilCompactionDATA[,"Date"]==(MinDate+j-1))
                 if(length(DateIndex)!=0){
                    if(SoilCompactionDATA[DateIndex,"Time"]==0){
                       TimeStamps[j] <- TimeOfDay/2400
                    }else{
                       TimeStamps[j] <- SoilCompactionDATA[DateIndex,"Time"]
                    }
                 }else{
                    TimeStamps[j] <- TimeOfDay/2400
                 }
               }
               ##Now, we have the timestamps for all the relevant times of day at which Jordan measured...
               
               ##So now, we grab the soil moisture data at the approrpriate timestamps...and ONLY those.
               SnippedSoilFile <- matrix(NA,(MaxDate-MinDate+1),ncol(SoilFile))
               colnames(SnippedSoilFile) <- colnames(SoilFile)
               for(k in 1:nrow(SnippedSoilFile)){
                 Time   <- TimeStamps[k] + MinDate + k -1
                 KeyRow <- SoilFile[SoilFile[,"TIMESTAMP"]<(Time+.005) & SoilFile[,"TIMESTAMP"]>(Time-.005),]
                 if(nrow(KeyRow)>1){
                   KeyRow <- KeyRow[1,]
                 }
                 SnippedSoilFile[k,] <- as.numeric(KeyRow) 
               }
               
               SoilFile       <- SnippedSoilFile
               #SoilFile      <- SoilFile[(SoilFile[,"TIMESTAMP"]-floor(SoilFile[,"TIMESTAMP"]))<(TimeOfDay/2400+.005) & (SoilFile[,"TIMESTAMP"]-floor(SoilFile[,"TIMESTAMP"]))>(TimeOfDay/2400-.005),]
               if(is.na(Code)){
                 ColumnCode    <- paste("Moisture",MoistureReading,sep="")
               }else{
                 ColumnCode    <- paste("M",Code,MoistureReading,sep="_")
               }
               print(ColumnCode)
               DailySoilData <- cbind(floor(SoilFile[,"TIMESTAMP"]),SoilFile[,ColumnCode])
               if(is.na(Code)){
                 ColumnCode <- paste(ColumnCode,strsplit(FileName,"_")[[1]][2],sep="_")
                 colnames(DailySoilData) <- c("Date",ColumnCode)
               }else{
                 colnames(DailySoilData) <- c("Date",ColumnCode)
               }
               DailySoilData[is.na(DailySoilData[,ColumnCode]),ColumnCode] <- 0
               MissingFrac   <- sum(DailySoilData[,ColumnCode]==0)/nrow(DailySoilData)
               ####If substantial fractions of data are missing, do NOT add this to the dataset.
               if(MissingFrac < 0.5){
                  for(j in 1:nrow(DailySoilData)){
                     DateIndex <- which(DailySoilData[j,"Date"]==SoilDATA[,"Date"])
                     if(j>2){
                         ###If data values are not "missing," but simply never updated...do not include
                         if(!is.na(as.numeric(DailySoilData[j,ColumnCode])) && !is.na(as.numeric(DailySoilData[j-1,ColumnCode])) && !is.na(as.numeric(DailySoilData[j-2,ColumnCode]))){
                               ###Also, if the value exceeds 1, the reading is not valid.
                               if(!(as.numeric(DailySoilData[j,ColumnCode])==as.numeric(DailySoilData[j-1,ColumnCode]) && as.numeric(DailySoilData[j-1,ColumnCode])==as.numeric(DailySoilData[j-2,ColumnCode])) && as.numeric(DailySoilData[j,ColumnCode]<1)){
                                  SoilDATA[DateIndex,"_"] <- DailySoilData[j,ColumnCode]
                               }
                         }
                     }
                  }
                  if(i == 1){
                     AllSoilDATA[,"_"] <- SoilDATA[,"_"]
                     colnames(AllSoilDATA)[which(colnames(AllSoilDATA)=="_")] <- ColumnCode
                  }else{
                     AllSoilDATA <- cbind(AllSoilDATA,SoilDATA[,"_"])
                     colnames(AllSoilDATA)[which(colnames(AllSoilDATA)=="")] <- ColumnCode
                  }
               }
          }
        }
    }
    AllSoilDATA
}

####Add precipitation data (using spatial interpolation)
####Add potential evaporation data
####Add Jordan's scale data
####Add a two-day history for all of these variables (more if needed)
AddPrecip_PotEvap_Scale_2DayHistory <- function(SoilCompactionDATA, PrecipDATA, ClimateDATA, COCORAHS=T, CCDATA=0, PublicData=T, ICN=F){
    ######Latitude and longitude for precipitation data...
    PrecipCoordinates <- matrix(NA,4,2)
    rownames(PrecipCoordinates) <- colnames(PrecipDATA)[2:ncol(PrecipDATA)]
    colnames(PrecipCoordinates) <- c("Lat","Long")
    ###Hard code the locations of the sensors, so we can query them by name!
    PRECIPS <- matrix(NA,2,4)
    colnames(PRECIPS) <- c("ICN","NE","SE","CEN")
    PRECIPS[,"ICN"]   <- c(40.08,-88.23) 
    PRECIPS[,"NE"]    <- c(40.06731,-88.19346)
    PRECIPS[,"SE"]    <- c(40.06377,-88.19342)
    PRECIPS[,"CEN"]   <- c(40.06363,-88.19737)
    PrecipCoordinates["ICN",] <- PRECIPS[,"ICN"]
    PrecipCoordinates["NE",] <- PRECIPS[,"NE"]
    PrecipCoordinates["SE",] <- PRECIPS[,"SE"]
    PrecipCoordinates["CEN",] <- PRECIPS[,"CEN"]
    if(COCORAHS==T){
       PrecipCoordinates <- AddCOCOCoordinates(CCDATA,PrecipCoordinates)
    }

    #For each spatial location, determine the estimated precipitation (inverse squared distance)
    PrecipEstimates   <- array(NA,nrow(SoilCompactionDATA))
    PrecipEstimates_1 <- array(NA,nrow(SoilCompactionDATA))
    PrecipEstimates_2 <- array(NA,nrow(SoilCompactionDATA))
    for(i in 1:nrow(SoilCompactionDATA)){
       Date            <- SoilCompactionDATA[i,"Date"]
       PrecipIndex     <- which(PrecipDATA[,"Date"]==Date) ###Index of the ICN and EBI sensors
       COCOIndex       <- which(CCDATA[,"X"] == Date) ###Index of COCORAHS dataset.
       Lat             <- SoilCompactionDATA[i,"Lat"]
       Long            <- SoilCompactionDATA[i,"Long"]
       DistanceWeights <- array(0,nrow(PrecipCoordinates))
       ###Get the inverse distance weights 
       for(j in 1:length(DistanceWeights)){
          DistanceWeights[j] <- 1/((Lat-PrecipCoordinates[j,"Lat"])^2 + (Long-PrecipCoordinates[j,"Long"])^2)
       }
       ###If we are using only public data, we need to remove the EBI sensors...
       if(PublicData == T){
         EBINames <- c("CEN","NE","SE")
         ZeroOut  <- which(rownames(PrecipCoordinates)==EBINames)
         DistanceWeights[ZeroOut] <- 0
       }
       ###Removes ICN data, so that only real-time, public sources can be used
       if(PublicData == T){
         EBINames <- c("CEN","NE","SE")
         ZeroOut  <- which(rownames(PrecipCoordinates)==EBINames)
         DistanceWeights[ZeroOut] <- 0
       }
       if(ICN == F){
         ZeroOut  <- which(rownames(PrecipCoordinates)=="ICN")
         DistanceWeights[ZeroOut] <- 0
       }
       
       
       #Find out which sensors have real data on this data
       UsefulICN_EBIVals <- as.numeric(!is.na(PrecipDATA[PrecipIndex,2:ncol(PrecipDATA)]))
       if(COCORAHS==T){
          UsefulCOCOVals  <- as.numeric(!is.na(CCDATA[COCOIndex,2:ncol(CCDATA)]))
          UsefulVals <- c(UsefulICN_EBIVals,UsefulCOCOVals)
       }else{
          UsefulVals <- UsefulICN_EBIVals 
       }
       TotalWeight  <- sum(DistanceWeights*UsefulVals)
       NormWeights  <- DistanceWeights/TotalWeight
       
       #Estimates on the day in question (or earlier) for ICN and EBI
       ICN_EBI_Rain   <- as.numeric(PrecipDATA[PrecipIndex,colnames(PrecipDATA)[2:ncol(PrecipDATA)]])
       ICN_EBI_Rain_1 <- as.numeric(PrecipDATA[PrecipIndex-1,colnames(PrecipDATA)[2:ncol(PrecipDATA)]])
       ICN_EBI_Rain_2 <- as.numeric(PrecipDATA[PrecipIndex-2,colnames(PrecipDATA)[2:ncol(PrecipDATA)]])
       #Estimates on the day in question for COCORAHS
       if(COCORAHS==T){
         COCO_Rain      <- as.numeric(CCDATA[COCOIndex,colnames(CCDATA)[2:ncol(CCDATA)]])
         COCO_Rain_1    <- as.numeric(CCDATA[COCOIndex-1,colnames(CCDATA)[2:ncol(CCDATA)]])
         COCO_Rain_2    <- as.numeric(CCDATA[COCOIndex-2,colnames(CCDATA)[2:ncol(CCDATA)]])
         DailyRain      <- c(ICN_EBI_Rain,COCO_Rain)
         DailyRain_1    <- c(ICN_EBI_Rain_1,COCO_Rain_1)
         DailyRain_2    <- c(ICN_EBI_Rain_2,COCO_Rain_2)
       }else{
         DailyRain      <- ICN_EBI_Rain
         DailyRain_1    <- ICN_EBI_Rain_1
         DailyRain_2    <- ICN_EBI_Rain_2
       }
       
       PrecipEstimates[i]   <- sum(NormWeights*DailyRain*UsefulVals,na.rm=T)
       PrecipEstimates_1[i] <- sum(NormWeights*DailyRain_1*UsefulVals,na.rm=T)
       PrecipEstimates_2[i] <- sum(NormWeights*DailyRain_2*UsefulVals,na.rm=T)
       #####FOR STUDYING SYSTEM MEMORY, NOT FOR PREDICTION...#####
    #PrecipEstimates_3[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-3),c("ICN","SE","NE","CEN")])
    #PrecipEstimates_4[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-4),c("ICN","SE","NE","CEN")])
    #PrecipEstimates_5[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-5),c("ICN","SE","NE","CEN")])
    #PrecipEstimates_6[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-6),c("ICN","SE","NE","CEN")])
    #PrecipEstimates_7[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-7),c("ICN","SE","NE","CEN")])
    #PrecipEstimates_8[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-8),c("ICN","SE","NE","CEN")])
    #PrecipEstimates_9[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-9),c("ICN","SE","NE","CEN")])
    #PrecipEstimates_10[i] <- sum(NormWeights*PrecipDATA[(PrecipIndex-10),c("ICN","SE","NE","CEN")])
    #########################################################################
    }

    #List of all dates from which soil compaction readings were taken
    UniqueDates        <- unique(SoilCompactionDATA[,"Date"])
    UniqueClimateDates <- unique(ClimateDATA[,"Date"])
    #UniqueRatingDates  <- unique(SoilRatingDATA[,"Date"])    ##Not needed without compaction

    #Integrate the data into one larger array
    #Ratings   <- matrix(NA,nrow(SoilCompactionDATA),2)
    Precip    <- array(NA,nrow(SoilCompactionDATA))
    Precip_1  <- array(NA,nrow(SoilCompactionDATA))
    Precip_2  <- array(NA,nrow(SoilCompactionDATA))

    ###################################################
    ##############FOR MEMORY CHECK#####################
    #Precip_3  <- array(NA,nrow(SoilCompactionDATA))
    #Precip_4  <- array(NA,nrow(SoilCompactionDATA))
    #Precip_5  <- array(NA,nrow(SoilCompactionDATA))
    #Precip_6  <- array(NA,nrow(SoilCompactionDATA))
    #Precip_7  <- array(NA,nrow(SoilCompactionDATA))
    #Precip_8  <- array(NA,nrow(SoilCompactionDATA))
    #Precip_9  <- array(NA,nrow(SoilCompactionDATA))
    #Precip_10  <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_3 <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_4 <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_5 <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_6 <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_7 <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_8 <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_9 <- array(NA,nrow(SoilCompactionDATA))
    #PotEvap_10 <- array(NA,nrow(SoilCompactionDATA))
    ###################################################
    ###################################################

    PotEvap <- array(NA,nrow(SoilCompactionDATA))
    PotEvap_1 <- array(NA,nrow(SoilCompactionDATA))
    PotEvap_2 <- array(NA,nrow(SoilCompactionDATA))
    for(i in 1:nrow(SoilCompactionDATA)){
       DateIndex       <- which(UniqueClimateDates==SoilCompactionDATA[i,"Date"])
       Date <- SoilCompactionDATA[i,"Date"]  ###Get Jordan's Rating for each day...
       #Ratings[i,1] <- SoilRatingDATA[which(SoilRatingDATA[,"Date"]==Date),"Scale.1"]
       #Ratings[i,2] <- SoilRatingDATA[which(SoilRatingDATA[,"Date"]==Date),"Scale.2"]
       Precip[i]    <- PrecipEstimates[i]
       Precip_1[i]  <- PrecipEstimates_1[i]   ###PRECIP DATA REPRESENTS THE SPATIAL ESTIMATES
       Precip_2[i]  <- PrecipEstimates_2[i]
          ###########AGAIN, FOR SYSTEM MEMORY CHECK################
          #########################################################
       #Precip_3[i]   <- PrecipEstimates_3[i]
       #Precip_4[i]   <- PrecipEstimates_4[i]
       #Precip_5[i]   <- PrecipEstimates_5[i]
       #Precip_6[i]   <- PrecipEstimates_6[i]
       #Precip_7[i]   <- PrecipEstimates_7[i]
       #Precip_8[i]   <- PrecipEstimates_8[i]
       #Precip_9[i]   <- PrecipEstimates_9[i]
       #Precip_10[i]  <- PrecipEstimates_10[i]
       #ExtraPrecip   <- cbind(Precip_3,Precip_4,Precip_5,Precip_6,Precip_7,Precip_8,Precip_9,Precip_10)
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       #PotEvap_3[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-3],"pot_evapot"]
       ###############################################################
       ###############################################################
       PotEvap[i]   <- ClimateDATA[which(ClimateDATA[,"Date"]==Date),"pot_evapot"]
       if(DateIndex > 2){ ##If we're after the second date, we can look back two days and know we have values
          PotEvap_1[i] <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-1],"pot_evapot"]
          PotEvap_2[i] <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-2],"pot_evapot"]
       }else if(DateIndex > 1){
          PotEvap_1[i] <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-1],"pot_evapot"]
          PotEvap_2[i] <- mean(ClimateDATA[,"pot_evapot"])
       }else{
          PotEvap_1[i] <- mean(ClimateDATA[,"pot_evapot"])
          PotEvap_2[i] <- mean(ClimateDATA[,"pot_evapot"])
       }

    }
    #colnames(Ratings) <-  c("Scale1","Scale2")
    DATA <- cbind(SoilCompactionDATA,Precip,PotEvap,Precip_1,Precip_2,PotEvap_1,PotEvap_2)#Ratings,)

    #Because we have multiple samples on one type of terrain in a given day, we must assess
    #the previous day's reading by averaging all readings on that terrain type on the most
    #recent day for which such values are available.

#    Depths       <- c("D0","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")
#      #Our estimate for yesterday's soil compaction
#    T_Minus_1    <- matrix(NA,nrow(SoilCompactionDATA),11)
#      #Our estimate for two days ago, soil compaction
#    T_Minus_2    <- matrix(NA,nrow(SoilCompactionDATA),11)
#    for(i in 1:nrow(SoilCompactionDATA)){
#       DateIndex       <- which(UniqueDates==SoilCompactionDATA[i,"Date"])
#       for(j in 1:11){
#          if(DateIndex==1){
#             T_Minus_1[i,j] <- mean(SoilCompactionDATA[SoilCompactionDATA[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],Depths[j]])
#             T_Minus_2[i,j] <- mean(SoilCompactionDATA[SoilCompactionDATA[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],Depths[j]])
#          }else if(DateIndex==2){
#             T_Minus_2[i,j] <- mean(SoilCompactionDATA[SoilCompactionDATA[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],Depths[j]])
#             DateSubset      <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]==UniqueDates[DateIndex-1],]
#             TerrainSubset   <- DateSubset[DateSubset[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],]
#             T_Minus_1[i,j]  <- mean(TerrainSubset[,Depths[j]])
#          }else{
#             DateSubset      <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]==UniqueDates[DateIndex-1],]
#             TerrainSubset   <- DateSubset[DateSubset[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],]
#             T_Minus_1[i,j]  <- mean(TerrainSubset[,Depths[j]])
#
#             DateSubset2     <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]==UniqueDates[DateIndex-2],]
#             TerrainSubset2  <- DateSubset2[DateSubset2[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],]
#             T_Minus_2[i,j]  <- mean(TerrainSubset2[,Depths[j]])
#          }
#       }
#    }
#    colnames(T_Minus_1) <- paste(Depths,"1",sep="_")
#    colnames(T_Minus_2) <- paste(Depths,"2",sep="_")
#    DATA <- cbind(DATA,T_Minus_1,T_Minus_2)

    #####Let us add soil ratings from previous days into the dataset...
#    Scale1_1 <- array(NA,nrow(SoilCompactionDATA))
#    Scale1_2 <- array(NA,nrow(SoilCompactionDATA))
#    Scale2_1 <- array(NA,nrow(SoilCompactionDATA))
#    Scale2_2 <- array(NA,nrow(SoilCompactionDATA))
#
#    for(i in 1:nrow(SoilCompactionDATA)){
#       DateIndex       <- which(UniqueRatingDates==SoilCompactionDATA[i,"Date"])
#          if(DateIndex==1){
#             Scale1_1[i] <- mean(SoilRatingDATA[,"Scale.1"])
#             Scale1_2[i] <- mean(SoilRatingDATA[,"Scale.1"])
#             Scale2_1[i] <- mean(SoilRatingDATA[,"Scale.2"])
#             Scale2_2[i] <- mean(SoilRatingDATA[,"Scale.2"])
#          }else if(DateIndex==2){
#             Scale1_1[i] <- SoilRatingDATA[SoilRatingDATA[,"Date"]==UniqueRatingDates[DateIndex-1],"Scale.1"]
#             Scale1_2[i] <- mean(SoilRatingDATA[,"Scale.1"])
#             Scale2_1[i] <- SoilRatingDATA[SoilRatingDATA[,"Date"]==UniqueRatingDates[DateIndex-1],"Scale.2"]
#             Scale2_2[i] <- mean(SoilRatingDATA[,"Scale.2"])
#          }else{
#             Scale1_1[i] <- SoilRatingDATA[SoilRatingDATA[,"Date"]==UniqueRatingDates[DateIndex-1],"Scale.1"]
#             Scale1_2[i] <- SoilRatingDATA[SoilRatingDATA[,"Date"]==UniqueRatingDates[DateIndex-2],"Scale.1"]
#             Scale2_1[i] <- SoilRatingDATA[SoilRatingDATA[,"Date"]==UniqueRatingDates[DateIndex-1],"Scale.2"]
#             Scale2_2[i] <- SoilRatingDATA[SoilRatingDATA[,"Date"]==UniqueRatingDates[DateIndex-2],"Scale.2"]
#          }
#    }
#    DATA <- cbind(DATA,Scale1_1,Scale1_2,Scale2_1,Scale2_2)
#    DATA
}

####Define readiness characteristics
###Create an array of "READINESS" (2-Class-Only)
GetReadiness <- function(DATA,DepVar,Threshhold){
    Readiness  <- array(NA,nrow(DATA))
    Readiness  <- DATA[,DepVar]>Threshhold
    Readiness[Readiness]          <- "READY"
    Readiness[Readiness=="FALSE"] <- "NOT READY"
     #DATA <- cbind(DATA,Readiness)
    DATA[,"Readiness"] <- Readiness
    DATA
}

###Split readiness into three classes, with two division points
###for example c(2.5,4.5) will define 1,2 = "GREEN", 3,4 = "YELLOW", and 5 = "RED"
GetReadiness3 <- function(DATA,DepVar,Threshholds){
   Readiness <- array(NA,nrow(DATA))
   Greens      <- DATA[,DepVar]>Threshholds[2]
   Reds    <- DATA[,DepVar]<Threshholds[1]
   Yellows   <- DATA[,DepVar]>=Threshholds[1]&DATA[,DepVar]<=Threshholds[2]
   Readiness[Reds]    <- "RED"
   Readiness[Yellows] <- "YELLOW"
   Readiness[Greens]  <- "GREEN"
   DATA <- cbind(DATA,Readiness)
}

###Grab the precipitation data
###If daily averaging is allowable, use all sensors.  
###Otherwise, use only sensors for which 15-minute data is available.
GetPrecipData <- function(SoilCompactionDATA,DefaultTimeOfDay,BeginDate,EndDate,DailyAveraged=F){
  if(DailyAveraged==T){
    PrecipDATA         <- GetData("2010Precip_Spatial")         #From ICN & 3 EBI Sensors
  }else{
     CEN <- GetEBIData("WeatherCEN")
     NE  <- GetEBIData("WeatherNE")
     SE  <- GetEBIData("WeatherSE")
     PrecipDATA <- matrix(NA,(EndDate-BeginDate+1),5)
     colnames(PrecipDATA) <- c("Date","CEN","NE","SE","ICN")
     PrecipDATA[,"ICN"]   <- NA 
     for(i in 1:nrow(PrecipDATA)){
        Date <- BeginDate + i - 1
        #Find the time of day at which Jordan takes his measurements
        if(length(SoilCompactionDATA[SoilCompactionDATA[,"Date"]==Date,"Time"])==0){
           TimeOfDay <- DefaultTimeOfDay/2400
        } else if (SoilCompactionDATA[SoilCompactionDATA[,"Date"]==Date,"Time"]==0){
           TimeOfDay <- DefaultTimeOfDay/2400
        } else {
           TimeOfDay <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]==Date,"Time"]
        }
        #print(TimeOfDay)
        PrecipDATA[i,"Date"] <- Date
        if(min(CEN[,"TIMESTAMP"]) < Date-1){
           DailyChunkCEN <- CEN[CEN[,"TIMESTAMP"]>=Date-1+TimeOfDay/2400 & CEN[,"TIMESTAMP"]<=Date+TimeOfDay/2400,] 
           PrecipDATA[i,"CEN"] <- sum(DailyChunkCEN[,"RAIN_Tot"])/100
        }
        if(min(NE[,"TIMESTAMP"]) < Date-1){
           DailyChunkNE <- NE[NE[,"TIMESTAMP"]>=Date-1+TimeOfDay/2400 & NE[,"TIMESTAMP"]<=Date+TimeOfDay/2400,] 
           PrecipDATA[i,"NE"] <- sum(DailyChunkNE[,"RAIN_Tot"])/100
        }
        if(min(SE[,"TIMESTAMP"]) < Date-1){
           DailyChunkSE <- SE[SE[,"TIMESTAMP"]>=Date-1+TimeOfDay/2400 & SE[,"TIMESTAMP"]<=Date+TimeOfDay/2400,] 
           PrecipDATA[i,"SE"] <- sum(DailyChunkSE[,"RAIN_Tot"])/100
        }
        
     }
  }
  PrecipDATA
}

MakePrecipICN_Only <- function(DATA,PrecipDATA,ClimateDATA){
  Precip   <- array(0,nrow(DATA))
  Precip_1 <- array(0,nrow(DATA)) 
  for(i in 1:nrow(DATA)){
     Date        <- DATA[i,"Date"]
     DateIndex   <- which(ClimateDATA[,"Date"]==Date)
     Precip[i]   <- ClimateDATA[(DateIndex-1),"precip"]
     Precip_1[i] <- ClimateDATA[(DateIndex-2),"precip"] 
  }
  DATA[,"Precip"]   <- Precip
  DATA[,"Precip_1"] <- Precip_1
  DATA
}

AddCOCOCoordinates <- function(CCDATA, PrecipCoordinates){
  NewStations <- colnames(CCDATA[2:ncol(CCDATA)])
  NewLats     <- array(0,ncol(CCDATA)-1)
  NewLongs     <- array(0,ncol(CCDATA)-1)
  for(i in 1:(ncol(CCDATA)-1)){
    NewLats[i]     <- as.numeric(coerce(CCDATA[1,(i+1)],to="numeric"))   
    NewLongs[i]    <- as.numeric(coerce(CCDATA[2,(i+1)],to="numeric"))
  }
  AllLats     <- c(as.numeric(PrecipCoordinates[,"Lat"]),NewLats)
  AllLongs    <- c(as.numeric(PrecipCoordinates[,"Long"]),NewLongs)
  AllStations <- c(rownames(PrecipCoordinates),NewStations) 
  NewPrecipCoordinates <- cbind(AllLats, AllLongs)
  colnames(NewPrecipCoordinates) <- c("Lat","Long")
  rownames(NewPrecipCoordinates) <- AllStations
  NewPrecipCoordinates
}

#read.delim("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\COCORAHS_DATA.txt", header = TRUE, sep = "\t", quote="\"", dec=".",fill = TRUE, comment.char="")

HighestFracOf24 <- function(Num){
###Takes in a decimal between 0 and 1, returns the highest fraction of 24 for 
#which data could be known.  For instance, at 12:30PM, the fraction is 25/48, 
#but we could only know hourly data until 12:00PM...thus this function returns
#the value of 24/48 = 0.5.
   Out = 0
   for(i in 1:24){
      if(i/24 <= Num){
        Out = i/24
      }
   }
   Out
}

Get48HoursChunk <- function(NexradDATA, Date, Frac){
###Takes in a date, a fractional time of day, and Nexrad data...
#returns the appropriate 48 rows.
   MaxTime <- Date+Frac
   Chunk   <- NexradDATA[NexradDATA[,"Date"] <= MaxTime,]
   L       <- nrow(Chunk)
   Chunk   <- Chunk[(L-47):L,]
   Chunk
}

GetNexradCoordinates <- function(NexradDATA){
###Takes a .csv file of Nexrad radar data, and returns the latitudes
#and longitudes
   LatLongs <- NexradDATA[1:2,2:ncol(NexradDATA)]
   LatLongs
}

GetDistWeights <- function(Lat,Long,LatLongs){
###Given a latitude and longitude of a specific location, return the 
#inverse square distance weights based upon the pythagorean distance to 
#each known precip location
   SqDiffs <- array(NA,ncol(LatLongs))
   for(i in 1:length(SqDiffs)){
       SqDiffs[i] <- 1/sqrt((LatLongs[1,i]-Lat)^2 + (LatLongs[2,i]-Long)^2)
   }
   TotalDiff   <- sum(SqDiffs)
   DistWeights <- SqDiffs/TotalDiff
}

ApplyWeights <- function(Row,Weights){
#Simply multiplies an array by a vector of weights and sums
  Output <- sum(Row*Weights)
  Output
}

WeightAndAggregatePrecip <- function(PrecipEsts, Day1Decay, Day2Decay){
###Takes in the spatially-weighted precips, then aggregates them into 
   Day1Weights <- (seq(1,24)/24)^Day1Decay
   Day1Weights <- Day1Weights*24/sum(Day1Weights)
   Day2Weights <- (seq(1,24)/24)^Day2Decay
   Day2Weights <- Day2Weights*24/sum(Day2Weights)
   Precip_1    <- sum(PrecipEsts[25:48] * Day1Weights)
   #Use ALL 48 Hours - we do not want hours 20-24 devalued due to decay.
   Precip_2    <- sum(PrecipEsts[1:48] * Day2Weights)
   PrecipHist  <- c(Precip_1,Precip_2)
   PrecipHist
}

GetNEXRADHourlyPrecipData <- function(SoilCompactionDATA,NexradDATA,DefaultTimeOfDay,BeginDate,EndDate,Day1Decay,Day2Decay){
#For each time and location within our validation set:
  #1.  Estimate the inverse distance weighted average of the rainfall at that
  #    location over the most recent 48 hours.
  #2.  Use an inverse temporal distance summation function to overweight
  #    rain which has occurred more recently at the location.
  
  #Store Precipitation Data
     Precips           <- matrix(NA,nrow(SoilCompactionDATA),2)
     colnames(Precips) <- c("Precip_1","Precip_2")

  #Clip SoilCompactionDATA
     SoilCompactionDATA <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]>= BeginDate,]
     SoilCompactionDATA <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]<= EndDate,]
  
  for(i in 1:nrow(SoilCompactionDATA)){
     if(SoilCompactionDATA[i,"Time"] == 0){
       LastHourDecimal = HighestFracOf24(DefaultTimeOfDay/2400)
     }else{
       LastHourDecimal = HighestFracOf24(SoilCompactionDATA[i,"Time"])
     }  
     
     #Grab the relevant Chunk   
     Chunk <- Get48HoursChunk(NexradDATA,SoilCompactionDATA[i,"Date"],LastHourDecimal)   
     #Clean it - set missing data to 0.
     Chunk[Chunk==-3200] <- 0
     
     #Perform the spatial inverse-distance weighting using the needed lat/long
     LatLongs        <- GetNexradCoordinates(NexradDATA)
     Lat             <- SoilCompactionDATA[i,"Lat"]
     Long            <- SoilCompactionDATA[i,"Long"]
     DistWeights     <- GetDistWeights(Lat,Long,LatLongs)
     PrecipEsts      <- apply(Chunk[,2:ncol(Chunk)],1,ApplyWeights,Weights=DistWeights)
     WeightedPrecip  <- WeightAndAggregatePrecip(PrecipEsts, Day1Decay, Day2Decay)
     Precips[i,]     <- WeightedPrecip
  } 
  SoilCompactionDATA <- cbind(SoilCompactionDATA,Precips)
  SoilCompactionDATA
}

AddPotEvap <- function(SoilCompactionDATA,ClimateDATA,Hours){
###Add potential evaporation data to the existing dataset.  Use the most recent
###user-chosen number of hours
   PotEvaps <- matrix(NA,nrow(SoilCompactionDATA),2)
   colnames(PotEvaps) <- c("PotEvap_1","PotEvap_2")
   for(i in 1:nrow(SoilCompactionDATA)){
       Date          <- SoilCompactionDATA[i,"Date"]
       Time          <- SoilCompactionDATA[i,"Time"]
       AdjTime       <- floor(Time*24)*100
       if(AdjTime==0){
          AdjTime <- 900
       }
       KeyIndex      <- which(ClimateDATA[,"Date"]==Date & ClimateDATA[,"Hour"]==AdjTime)
       PotEvaps[i,1] <- sum(ClimateDATA[(KeyIndex-Hours+1):KeyIndex,"PotEvap"])
       PotEvaps[i,2] <- sum(ClimateDATA[(KeyIndex-Hours*2+1):(KeyIndex-Hours),"PotEvap"])
   }
   SoilCompactionDATA  <- cbind(SoilCompactionDATA,PotEvaps)
   SoilCompactionDATA
}  


##Pre-Processing of DATA, add the precip data based on hourly time-stamps of
##validation measurements, and pot.evapotranspiration in the same fashion
PreProcess <- function(SoilCompactionDATA,NexradDATA,DefaultTimeOfDay,BeginDate,EndDate,Day1Decay,Day2Decay,ClimateDATA,Hrs){
  SoilCompactionDATA <- SoilCompactionDATA[SoilCompactionDATA[,"Long"]>-88.2,] ##Hard-Coded, ensures we just grab location 1.
  SoilCompactionDATA <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]<=EndDate,] ##This way we will not run beyond the end of ICN

  #Use Nexrad data to gather precipitation information for each revelant timestep
  SoilCompactionDATA <- GetNEXRADHourlyPrecipData(SoilCompactionDATA,NexradDATA,DefaultTimeOfDay,BeginDate,EndDate,Day1Decay,Day2Decay)
  #Add the potential evaporation data from ICN
  SoilCompactionDATA <- AddPotEvap(SoilCompactionDATA,ClimateDATA,Hrs)
  SoilCompactionDATA
  }
  
DefineClasses <- function(DATA,DepVar,Threshhold,Threshholds,Readiness){
  if(NumClasses==2){
     DATA <- GetReadiness(DATA,DepVar,Threshhold)  
  }else if(NumClasses==3){
     DATA <- GetReadiness3(DATA,DepVar,Threshholds)  
  }else if(NumClasses==5){
     Readiness <- DATA[,"Scale1"]
    DATA <- cbind(DATA,Readiness)
  }
  DATA
}

RemoveEasyDays <- function(MDWR,Remove,DATA){
  for(i in (MDWR+1):nrow(DATA)){
     DryCounter = 0
     for(j in 1:MDWR){
       if(DATA[(i-j),"Precip_1"]==0){
         DryCounter = DryCounter + 1
        }
     }
     if(DryCounter==MDWR & DATA[i,"Precip_1"]==0){                  
       Remove <- c(Remove,i)
     }
  } 
  Remove <- Remove[2:length(Remove)]
  DATA <- DATA[-Remove,]
  DATA
}

#####THIS IS A CODE DUMP FROM THE MAIN WHICH WILL LIKELY NEVER BE USED AGAIN!
###Add yesterday's estimate of readiness, obscuring where necessary...
#Obscure<-F
#Scale1_1 <- array(NA,nrow(DATA))
#for(i in 1:nrow(DATA)){
#  if(i==1){
#    Scale1_1[i] <- 3
#  }else{
#    if(DATA[i-1,"Scale1"]==1 || DATA[i-1,"Scale1"]==5 || Obscure==F){
#      Scale1_1[i] <- DATA[i-1,"Scale1"]
#    }else{
#      Scale1_1[i] <- 3
#    }
#  }                                                                      
#}
#
#DATA <- cbind(DATA,Scale1_1)