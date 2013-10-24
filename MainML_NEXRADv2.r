rm(list=ls())

source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\MLCompendium_EC2.r")
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\AssortedSoilFunctions.r")

BeginDate              <- 40300
EndDate                <- 40445  ##Currently limited by ICN data...
DefaultTimeOfDay       <- 900 ###Measured as 0 to 2400
UsePublicDataOnly      <- T
UseICN                 <- F
ICNOnly                <- F
Day1Decay              <- 3   #Describe how strongly we overweight more recent rain
Day2Decay              <- 0   #If set to 0 - all rain data is weighted equally
Hrs                    <- 48  #How far to look back for PotEvap?
LongHist               <- 48  #How far back to look
iterations             <- 10  #For boosted perceptron only
Cycles                 <- 20 #For SGD boosted perceptron only
Thresh                 <- 2.5 #Also BP
SGD                    <- F   #Type of fitting in BP

#Read in our data.
SoilCompactionDATA <- GetData("SoilReadiness")              #From Jordan Pitcher (Location 1, begins 4/20/10)
    #Now no longer compaction, just readiness               #From Jordan Pitcher (Location 2, begins 5/20/10, we begin 6/29)
    #Now Used!#                                             #From Jordan Pitcher (Location 3, begins 5/20/10, we begin 6/29) 
                                                                
ClimateDATA        <- GetData("HourlyPotEvap")              #From ICN, Champaign Sensor (HOURLY)
NexradDATA         <- GetData("UrbanaHourlyNEXRAD")         #From Q2 Data Query
  ##########AT THIS POINT, WE CAN READ IN ADDITIONAL NEXRAD FILES I'VE GATHERED
  ##########AND MERGE THEM WITH THIS URBANA FILE.  HOWEVER, AT THE MOMENT, WE HAVE
  ##########ONLY ONE SET OF VALIDATION DATA - SO THIS IS ALL THAT IS NEEDEDs

DailyAveraged = F
DATA <- PreProcess(SoilCompactionDATA,NexradDATA,DefaultTimeOfDay,BeginDate,EndDate,Day1Decay,Day2Decay,ClimateDATA,Hrs)

#####THIS IS WHERE WE DEFINE READINESS CHARACTERISTICS#####
#DepVar     <- "D3"  #We will try to predict readiness based on 3in depth readings
#Threshhold <- 97
DepVar      <- "Scale1"
Threshhold  <- 2.5
Threshholds <- c(2.5,4.5)
NumClasses  <- 2
DATA <- DefineClasses(DATA,DepVar,Threshhold,Threshholds,Readiness)


###If we are eliminating DATA beyond N consecutive dry days.
MDWR = 3 ##MaximumDaysWithoutRain
Remove = 1
DATA <- RemoveEasyDays(MDWR,Remove,DATA)


###Define Independent Variables
Yesterday   <- paste(DepVar,"1",sep="_")
Yesterday_2 <- paste(DepVar,"2",sep="_")
IndVars <- c("Precip_1","Precip_2","PotEvap_1")

MinNodeSize <- 2
MinSplitSize <- 8
HistEndDate <- 40353    ###Bring this number forward if not running a backtest.
EndTestDate <- 40445

write.csv(DATA,"C:\\Users\\Evan\\Desktop\\Dump.csv")
DATA <- read.csv("C:\\Users\\Evan\\Desktop\\Dump.csv")      

###Actual Results
ActualResults <- DATA[DATA[,"Date"]>HistEndDate,"Readiness"]
###################
RunBacktest <- 0###  DO WE WANT TO RUN TESTS OR INDIVIDUAL CONDITIONS?
###################



################## ALWAYS RUN WHAT IS ABOVE #########################
#####################################################################
################# CHOOSE THE ALGORITHM BELOW ########################

if(RunBacktest == 1){
###############RUN THE REGRESSION TREE#####################
####USE THIS IF WE WISH TO BACKTEST A CERTAIN TIMEFRAME############
MyTree <- TrainAndTestRegTrees(IndVars, c("Readiness"), DATA, MinNodeSize, MinSplitSize, HistEndDate, EndTestDate)
OutputTREE <- cbind(DATA[DATA[,"Date"]>HistEndDate,],MyTree)
colnames(OutputTREE)[which(colnames(OutputTREE)=="MyTree")] <- "ModelPred"
}else{
  RTree <- PutTreeInMemory(IndVars, c("Readiness"), DATA, MinNodeSize, MinSplitSize, HistEndDate, EndTestDate)
}

####TEST AN INDIVIDUAL SET OF CONDITIONS ON THE REGRESSION TREE
#The "c(0,0,0)" represents a sample set of input conditions for the variables
#specified in the variables IndVars.  Essentially, we define variables, then
#input a given set of values for those variables, then test on our tree
Class <- TestConditionsRTree(RTree,IndVars,c(0,0,0),Classify=F)
                                 
                                 
#For KNN - set classify to true for most likely outcome.
#RUN THE KNN Algorithm
###Current Best, 85.5%###

k <- 11
RegType <- "KNN"
HistEndDate <- 40353
EndForecastDate <- 40445
IndVars <- c("Precip_1","PotEvap_1","Precip_2")

if(RunBacktest == 1){
###############RUN THE KNN#####################
####USE THIS IF WE WISH TO BACKTEST A CERTAIN TIMEFRAME############
  ###Model Estimates
  MyKNN <- KNN(IndVars,"Readiness",k,DATA,HistEndDate,EndForecastDate,RegType,Unif=F,Classify=T)
  K_Counts <- as.numeric(MyKNN[,"K_Counts"])
  K_Prob = (K_Counts + 1)/(k+2)
  
  ###Actual Results
  ActualResults <- DATA[DATA[,"Date"]>HistEndDate,"Readiness"]
  
  OutputKNN <- cbind(DATA[DATA[,"Date"]>HistEndDate,],MyKNN,K_Prob)
  colnames(OutputKNN)[which(colnames(OutputKNN)=="DepEstimates")] <- "ModelPred"

}else{
   HistoricalSet <- PrepKNN(IndVars,"Readiness",k,DATA,HistEndDate,EndForecastDate,RegType,Unif,Classify)
}

####TEST AN INDIVIDUAL SET OF CONDITIONS ON THE KNN
#The "c(0,0,0)" represents a sample set of input conditions for the variables
#specified in the variables IndVars.  Essentially, we define variables, then
#input a given set of values for those variables, then test on KNN
Class <- TestConditionsKNN(HistoricalSet,IndVars,"Readiness",k,RegType,Unif=F,c(0,0,0))


if(RunBacktest==1){
###############RUN THE BOOSTED PERCEPTRON#####################
####USE THIS IF WE WISH TO BACKTEST A CERTAIN TIMEFRAME############                                                                   BP <- BoostedPerceptron(DATA, IndVars, DepVar, Thresh, HistEndDate, iterations, Cycles,SGD)
  ModelPred <- array(NA,length(BP))
  ModelPred[BP >= .5] <- "READY" ; ModelPred[BP < .5] <- "NOT READY"
  BP_Prob <- BP
  OutputBP <- cbind(DATA[DATA[,"Date"]>HistEndDate,],ModelPred,BP_Prob) 

}else{
  BPReg <- PrepBoostedPercep(DATA, IndVars, DepVar, Thresh, HistEndDate, iterations, Cycles, SGD)
}
     
     
####TEST AN INDIVIDUAL SET OF CONDITIONS ON THE BOOSTED PERCEPTRON
#The "c(0,0,0)" represents a sample set of input conditions for the variables
#specified in the variables IndVars.  Essentially, we define variables, then
#input a given set of values for those variables, then test on the BP algorithm
Class <- TestConditionsBP(BPReg,c(0,0,0))
                                                               