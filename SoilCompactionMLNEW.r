rm(list=ls())

source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\MLCompendium_EC.r")

#Read in our data.
SoilCompactionDATA <- GetData("SoilCompaction")      #From Jordan Pitcher (Location 1, begins 4/20/10)
ClimateDATA <- GetData("2010Precip_PotEvap")         #From ICN, Champaign Sensor
SoilRatingDATA <- GetData("SoilRatings")             #From Jordan Pitcher (begins 5/11/10)

#Note:  Terrain=1, HillTop / Terrain=2, SlightSlope / Terrain=3, BottomSlope

#List of all dates from which soil compaction readings were taken
UniqueDates        <- unique(SoilCompactionDATA[,"Date"])
UniqueClimateDates <- unique(ClimateDATA[,"Date"])

#Integrate the data into one larger array
Ratings   <- matrix(NA,nrow(SoilCompactionDATA),2)
Precip    <- array(NA,nrow(SoilCompactionDATA))
Precip_1  <- array(NA,nrow(SoilCompactionDATA))
Precip_2  <- array(NA,nrow(SoilCompactionDATA))
PotEvap <- array(NA,nrow(SoilCompactionDATA))
PotEvap_1 <- array(NA,nrow(SoilCompactionDATA))
PotEvap_2 <- array(NA,nrow(SoilCompactionDATA))
for(i in 1:nrow(SoilCompactionDATA)){
   DateIndex       <- which(UniqueClimateDates==SoilCompactionDATA[i,"Date"])
   Date <- SoilCompactionDATA[i,"Date"]  ###Get Jordan's Rating for each day...
   Ratings[i,1] <- SoilRatingDATA[which(SoilRatingDATA[,"Date"]==Date),"Scale.1"]
   Ratings[i,2] <- SoilRatingDATA[which(SoilRatingDATA[,"Date"]==Date),"Scale.2"]
   Precip[i]    <- ClimateDATA[which(ClimateDATA[,"Date"]==Date),"precip"]
   PotEvap[i]   <- ClimateDATA[which(ClimateDATA[,"Date"]==Date),"pot_evapot"]
   if(DateIndex > 2){ ##If we're after the second date, we can look back two days and know we have values
      Precip_1[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-1],"precip"]
      Precip_2[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-2],"precip"]
      PotEvap_1[i] <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-1],"pot_evapot"]
      PotEvap_2[i] <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-2],"pot_evapot"]
   }else if(DateIndex > 1){
      Precip_1[i]  <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-1],"precip"]
      Precip_2[i]  <- mean(ClimateDATA[,"precip"])
      PotEvap_1[i] <- ClimateDATA[ClimateDATA[,"Date"]==UniqueClimateDates[DateIndex-1],"pot_evapot"]
      PotEvap_2[i] <- mean(ClimateDATA[,"pot_evapot"])      
   }else{
      Precip_1[i]  <- mean(ClimateDATA[,"precip"])
      Precip_2[i]  <- mean(ClimateDATA[,"precip"])
      PotEvap_1[i] <- mean(ClimateDATA[,"pot_evapot"])  
      PotEvap_2[i] <- mean(ClimateDATA[,"pot_evapot"])            
   }
   
}
colnames(Ratings) <-  c("Scale1","Scale2")
DATA <- cbind(SoilCompactionDATA,Ratings,Precip,PotEvap,Precip_1,Precip_2,PotEvap_1,PotEvap_2)

#Because we have multiple samples on one type of terrain in a given day, we must assess
#the previous day's reading by averaging all readings on that terrain type on the most
#recent day for which such values are available.

Depths       <- c("D0","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")
  #Our estimate for yesterday's soil compaction
T_Minus_1    <- matrix(NA,nrow(SoilCompactionDATA),11)
  #Our estimate for two days ago, soil compaction
T_Minus_2    <- matrix(NA,nrow(SoilCompactionDATA),11)
for(i in 1:nrow(SoilCompactionDATA)){
   DateIndex       <- which(UniqueDates==SoilCompactionDATA[i,"Date"])
   for(j in 1:11){ 
      if(DateIndex==1){
         T_Minus_1[i,j] <- mean(SoilCompactionDATA[SoilCompactionDATA[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],Depths[j]])
         T_Minus_2[i,j] <- mean(SoilCompactionDATA[SoilCompactionDATA[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],Depths[j]])  
      }else if(DateIndex==2){
         T_Minus_2[i,j] <- mean(SoilCompactionDATA[SoilCompactionDATA[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],Depths[j]])  
         DateSubset      <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]==UniqueDates[DateIndex-1],]
         TerrainSubset   <- DateSubset[DateSubset[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],] 
         T_Minus_1[i,j]  <- mean(TerrainSubset[,Depths[j]]) 
      }else{
         DateSubset      <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]==UniqueDates[DateIndex-1],]
         TerrainSubset   <- DateSubset[DateSubset[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],] 
         T_Minus_1[i,j]  <- mean(TerrainSubset[,Depths[j]])
         
         DateSubset2     <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]==UniqueDates[DateIndex-2],]
         TerrainSubset2  <- DateSubset2[DateSubset2[,"Terrain"]==SoilCompactionDATA[i,"Terrain"],]
         T_Minus_2[i,j]  <- mean(TerrainSubset2[,Depths[j]])
      }
   }
}
colnames(T_Minus_1) <- paste(Depths,"1",sep="_")
colnames(T_Minus_2) <- paste(Depths,"2",sep="_")
DATA <- cbind(DATA,T_Minus_1,T_Minus_2)

#####THIS IS WHERE WE DEFINE READINESS CHARACTERISTICS#####
DepVar     <- "D3"  #We will try to predict readiness based on 3in depth readings
Threshhold <- 120   #This determines whether or not the field is "READY"

###Create an array of "READINESS"
Readiness  <- array(NA,nrow(DATA))
Readiness  <- DATA[,DepVar]>Threshhold
Readiness[Readiness]          <- "READY"
Readiness[Readiness=="FALSE"] <- "NOT READY"
 #DATA <- cbind(DATA,Readiness)
DATA[,"Readiness"] <- Readiness

###Define Independent Variables
Yesterday   <- paste(DepVar,"1",sep="_")
Yesterday_2 <- paste(DepVar,"2",sep="_")
IndVars <- c(Yesterday,Yesterday_2,"Precip_1","PotEvap_1","Terrain")

MinNodeSize <- 2
MinSplitSize <- 8
HistEndDate <- 40336
EndTestDate <- 40365

write.csv(DATA,"C:\\Users\\Evan\\Desktop\\Dump.csv")
DATA <- read.csv("C:\\Users\\Evan\\Desktop\\Dump.csv")        

###Actual Results
ActualResults <- DATA[DATA[,"Date"]>HistEndDate,"Readiness"]

###############RUN THE REGRESSION TREE#####################
MyTree <- TrainAndTestRegTrees(IndVars, c("Readiness"), DATA, MinNodeSize, MinSplitSize, HistEndDate, EndTestDate)

SuccessRateRTREE <- sum(MyTree==ActualResults)/length(ActualResults)
KindSR_RTREE <- array(NA,length(ActualResults))
for(i in 1:length(ActualResults)){
   if((MyTree[i] == ActualResults[i]) | (MyTree[i] == "NOT READY" & ActualResults[i] == "READY IN TWO DAYS") | (MyTree[i] == "READY IN TWO DAYS" & ActualResults[i] == "NOT READY")){
      KindSR_RTREE[i] = 1
   }else{
      KindSR_RTREE[i] = 0
   }
}
OutputTREE <- cbind(MyTree,ActualResults)
                                 
                                 
#For KNN - set classify to true for most likely outcome.
#RUN THE KNN Algorithm

k <- 125
RegType <- "KNN"
EndForecastDate <- 40365
IndVars <- c(Yesterday, "Precip_1","PotEvap_1","Terrain")
###Model Estimates
MyKNN <- KNN(IndVars,"Readiness",k,DATA,HistEndDate,EndForecastDate,RegType,Unif=F,Classify=T)

SuccessRateKNN   <- sum(MyKNN==ActualResults)/length(ActualResults)
KindSR_KNN <- array(NA,length(ActualResults))
for(i in 1:length(ActualResults)){
   if((MyKNN[i] == ActualResults[i]) | (MyKNN[i] == "NOT READY" & ActualResults[i] == "READY IN TWO DAYS") | (MyKNN[i] == "READY IN TWO DAYS" & ActualResults[i] == "NOT READY")){
      KindSR_KNN[i] = 1
   }else{
      KindSR_KNN[i] = 0
   }
}
OutputKNN <- cbind(MyKNN,ActualResults)


                                                               