rm(list=ls())

source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\MLCompendium_EC.r")
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\AssortedSoilFunctions.r")

BeginDate              <- 40300
EndDate                <- 40445  ##Currently limited by ICN data...
DefaultTimeOfDay       <- 900 ###Measured as 0 to 2400
COCORAHS               <- 1   ###If 1, use all available COCORAHS sensors...
UsePublicDataOnly      <- T
UseICN                 <- F
ICNOnly                <- F

#Read in our data.
SoilCompactionDATA <- GetData("SoilReadiness")              #From Jordan Pitcher (Location 1, begins 4/20/10)
    #Now no longer compaction, just readiness               #From Jordan Pitcher (Location 2, begins 5/20/10, we begin 6/29)
    #Now Used!#                                             #From Jordan Pitcher (Location 3, begins 5/20/10, we begin 6/29) 
                                                                
ClimateDATA        <- GetData("2010Precip_PotEvap")         #From ICN, Champaign Sensor
#SoilRatingDATA     <- GetData("SoilRatings")               #From Jordan Pitcher (begins 5/11/10)

DailyAveraged = F
##
SoilCompactionDATA <- SoilCompactionDATA[SoilCompactionDATA[,"Long"]>-88.2,] ##Hard-Coded, ensures we just grab location 1.
SoilCompactionDATA <- SoilCompactionDATA[SoilCompactionDATA[,"Date"]<=EndDate,] ##This way we will not run beyond the end of ICN
##
PrecipDATA <- GetPrecipData(SoilCompactionDATA,DefaultTimeOfDay,BeginDate,EndDate,DailyAveraged)

if(COCORAHS == 1){
   COCORAHS_DATA <- GetData("COCORAHS_DATA")
}else{
   COCORAHS_DATA = 0 ###Just so we have some dummy argument to pass to the AddPrecip function
}


#Note:  Terrain=1, HillTop / Terrain=2, SlightSlope / Terrain=3, BottomSlope
 
#This accesses the data file with relevant information regarding the small plot files
SoilDetails        <- GetEBIData("SmallPlotInfo") 

#Iterate over all sources of soil data and construct an appropriate soil moisture dataset.
MoistureReading <- "A5_Avg" ###Determines which set of readings we must analyze
AllSoilDATA     <- GetAllSoilData(DefaultTimeOfDay,MoistureReading,BeginDate,EndDate,SoilCompactionDATA,SoilDetails)
Name            <- paste("SM",MoistureReading,sep="_")

####Add precipitation data (using spatial interpolation)
####Add potential evaporation data
####Add Jordan's scale data
####Add a two-day history for all of these variables (more if needed)
DATA            <- AddPrecip_PotEvap_Scale_2DayHistory(SoilCompactionDATA, PrecipDATA, ClimateDATA, COCORAHS, COCORAHS_DATA, UsePublicDataOnly, UseICN)
DATA[DATA[,"Time"]==0,"Time"] <- DefaultTimeOfDay/2400

####If we insist upon only using publically available data, replace "Precip" and "Precip_1"
####with only the values used for the ICN
if(ICNOnly==T){
  DATA <- MakePrecipICN_Only(DATA,PrecipDATA,ClimateDATA)
}

####Determining the Memory of the System
#DepVar   <- DATA[,"Scale1"]
#IndVars  <- cbind(DATA[,"Precip_1"],DATA[,"PotEvap_1"],DATA[,"Precip_2"],DATA[,"PotEvap_2"])
#M_Reg    <- lm(DepVar ~ IndVars)
#PredVals <- as.numeric(M_Reg[["fitted.values"]])

#SSE       <- sum((DepVar-PredVals)^2)
#SSR       <- sum((DepVar-mean(DepVar))^2)
#R_Squared <- 1-SSE/SSR


####FINDING A PSI/SoilMoisture THRESHOLD######
#MaxAccuracy <- 0

#for (i in seq(0.2,0.4,0.001)){
#   i=0.205
#   Correct   <- sum(DATA[,Name]>i & DATA[,"Scale1"]>2.5) + sum(DATA[,Name]<i & DATA[,"Scale1"]<2.5) 
#   Incorrect <- sum(DATA[,Name]>i & DATA[,"Scale1"]<2.5) + sum(DATA[,Name]<i & DATA[,"Scale1"]>2.5)
#   Accuracy  <- Correct/(Correct + Incorrect)
#   if(Accuracy > MaxAccuracy){
#      MaxAccuracy <- Accuracy
#      print(i)
#      print(paste("Accuracy = ",Accuracy))
#   }
#}
######


#####THIS IS WHERE WE DEFINE READINESS CHARACTERISTICS#####
#DepVar     <- "D3"  #We will try to predict readiness based on 3in depth readings
#Threshhold <- 97
DepVar      <- "Scale1"
Threshhold  <- 2.5
Threshholds <- c(2.5,4.5)
NumClasses  <- 2
if(NumClasses==2){
   DATA <- GetReadiness(DATA,DepVar,Threshhold)  
}else if(NumClasses==3){
   DATA <- GetReadiness3(DATA,DepVar,Threshholds)  
}else if(NumClasses==5){
   Readiness <- DATA[,"Scale1"]
   DATA <- cbind(DATA,Readiness)
}

###Add soil moisture readings, using spatial estimates for every data point
DATA <- AddSoilMoisture(DATA, AllSoilDATA, SoilDetails, Name)

###If we are eliminating DATA beyond N consecutive dry days.
MDWR = 3 ##MaximumDaysWithoutRain
Remove = 1
for(i in (MDWR+1):nrow(DATA)){
print(i)
   DryCounter = 0
   for(j in 1:MDWR){
      if(DATA[(i-j),"Precip"]==0){
         DryCounter = DryCounter + 1
      }
   }
   if(DryCounter==MDWR & DATA[i,"Precip"]==0){                  
      Remove <- c(Remove,i)
   }
}
Remove <- Remove[2:length(Remove)]
DATA <- DATA[-Remove,]

###Add yesterday's estimate of readiness, obscuring where necessary...
Obscure<-F
Scale1_1 <- array(NA,nrow(DATA))
for(i in 1:nrow(DATA)){
  if(i==1){
    Scale1_1[i] <- 3
  }else{
    if(DATA[i-1,"Scale1"]==1 || DATA[i-1,"Scale1"]==5 || Obscure==F){
      Scale1_1[i] <- DATA[i-1,"Scale1"]
    }else{
      Scale1_1[i] <- 3
    }
  }                                                                      
}

DATA <- cbind(DATA,Scale1_1)

###Define Independent Variables
Yesterday   <- paste(DepVar,"1",sep="_")
Yesterday_2 <- paste(DepVar,"2",sep="_")
SM_1        <- paste(Name,"1",sep="_")
SM_2        <- paste(Name,"2",sep="_")
IndVars <- c("Precip","Precip_1","PotEvap_1")

MinNodeSize <- 2
MinSplitSize <- 8
HistEndDate <- 40357
EndTestDate <- 40445

write.csv(DATA,"C:\\Users\\Evan\\Desktop\\Dump.csv")
DATA <- read.csv("C:\\Users\\Evan\\Desktop\\Dump.csv")      

SmallDATA <- DATA[,c("Date","Precip","PotEvap_1","Precip_1","SM_A5_Avg","Scale1","Readiness")] 
Values <- array(0,nrow(SmallDATA))                                                       
for(i in 1:nrow(SmallDATA)){
   if(SmallDATA[i,"Readiness"]=="READY"){
      Values[i] <- 1
   } 
}
SmallDATA <- cbind(SmallDATA,Values)

 

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
###Current Best, 85.5%###

k <- 5
RegType <- "KNN"
HistEndDate <- 40353
EndForecastDate <- 40445
IndVars <- c("Precip","PotEvap_1","Precip_1")
###Model Estimates
MyKNN <- KNN(IndVars,"Readiness",k,DATA,HistEndDate,EndForecastDate,RegType,Unif=F,Classify=F)
                                              
#This is the data we want to analyze to assess possible sources of error
KeyVars <- c("Time","Date","Precip","Precip_1","PotEvap_1","Scale1","Scale2")
ShowData <- DATA[DATA[,"Date"]<=EndForecastDate,KeyVars]
                                                                                                                
###Actual Results
ActualResults <- DATA[DATA[,"Date"]>HistEndDate,"Readiness"]

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


                                                               