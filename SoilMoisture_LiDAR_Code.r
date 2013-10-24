#Read Nexrad Files, produce appropriate precipitation series values at every
#lat-lon location.

rm(list=ls())
#Place functions into memory that enable soil moisture ML algorithms
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\NationalSoilMoistureFunctions.r")

############ INPUTS #####################

#Acceptable NEXRAD names, ("Illinois", "IowaState", "Purdue")
PrecipSite <- "IowaState"
#Acceptable folder names ("University.of.Illinois", Purdue", "Iowa.State")
SMFolder <- "Iowa.State"
#Path to NEXRAD data and soil moisture data
TestSitePath <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\CropSense Files 10-3\\"
#If we are missing precipitation values, to what shall we set that time stamp?
Default <- 0
#Name for the site within the Illinois directory
SoilMoistureSite <- "201A"
#Depth of soil profile in which we are interested
zTrain <- 4
zTest  <- 4
#Name of lat/lon list
CoordinatesBySitePath <- "SiteList_w_LatLon"
#Path for the storage location of the clean datafiles and parameter list
EndPath <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SCAN_Data_Long_LinConnect_4\\"
#Which site shall we use for our calibrated parameters?
ID <- 2111
LongTime <- 2000 #Long enough for ANY beta series
#If we adjust based on topography, we can input a new value here...0 is a default
newC4 <- 0
#Independent and dependent variables
IndVars  <- c("TimeStamps","TimeOfDays","BetaDiffs","SMSeries","BetaSeries") 
DepVar   <- "ErrorSeq"  #What we are attempting to model in the ML overlay phase

##(In mm per hour)...if we see rainfall rates above this, we do not "count it"
MaxInfiltration <- 1000
MaxInfiltration <- MaxInfiltration/10/2.54
##########################################
                                                 
SMDepthVar <- paste("SMS.I.1..",zTrain,sep="")                

#Read in the parameter list
ParameterList <- read.csv(paste(EndPath,"ParamList.csv",sep=""))

#Get the latitude and longitude of the soil moisture site
LatLong <- GetSMLatLon(CoordinatesBySitePath,TestSitePath, SoilMoistureSite)
SiteLat <- LatLong[1]
SiteLon <- LatLong[2]

#Get the precipitation time series(inverse-distance-weighted), dates, & times
PrecipsDatesAndTimes <- GetPrecipSeries(PrecipSite,TestSitePath,SiteLat,SiteLon,Default)
Precips     <- PrecipsDatesAndTimes[,"SitePrecip"]
#####HANDLING MAX INFILTRATION######
Precips[Precips>MaxInfiltration] <- MaxInfiltration

PrecipDates <- floor(PrecipsDatesAndTimes[,"RainDates"])
PrecipDates  <- FinanceDateToDOYFormat(PrecipDates) #Convert to YYYYDOY format
PrecipTimes <- PrecipsDatesAndTimes[,"RainTimes"]

#Get the soil moisture file
SMData <- GetSMFile(TestSitePath,SMFolder,SoilMoistureSite)

#Obtain date and time stamps from SMDA / Convert the soil moisture file dates
#to date/time stamps of YYYYMMDD,.58333 e.g.
DatesAndTimes <- GetDatesAndTimes(SMData)
DatesAndTimes[,"Date"] <- FinanceDateToDOYFormat(DatesAndTimes[,"Date"])

#Grab the soil moisture series
SMSensorLabel <- paste("Sensor.",zTest,sep="")
SMSeries <- SMData[,SMSensorLabel] #In (in / 4*in)
SMSeries <- SMSeries/4 * 100 #Converts to 0 to 100 moisture levels

#Merge the precipitation series (which is hourly) with the soil moisture series
#which is half-hourly, and somewhat skewed
minDate <- DatesAndTimes[1,"Date"] #first day of soil moisture
minTime <- DatesAndTimes[1,"Time"] #earliest soil moisture time stamp
L <- nrow(SMData)
maxDate <- DatesAndTimes[L,"Date"]
maxTime <- DatesAndTimes[L,"Time"]
SMMerged <- MergeSMToPrecip(Precips,PrecipDates,PrecipTimes,DatesAndTimes,SMSeries,
                            minDate, maxDate, minTime, maxTime)

#Build the "SiteData" matrix used in soil moisture functions
SiteData <- BuildSiteData(ID, SMMerged, PrecipDates, PrecipTimes, Precips, SMDepthVar)

#Running the SM Model from Pan et al with parameters from previous calibration
i <- which(ParameterList[,"SiteID"]==ID) #relevant index in parameter list
n <- ParameterList[i,"n"] #length of window (in hours)
PrecipSeries <- SiteData[1:nrow(SiteData),"PRCP.H.1..in."]
TimeStamps   <- SiteData[1:nrow(SiteData),"Date"] %% 1000 + 
                SiteData[1:nrow(SiteData),"Time"] - 1
TimeOfDays   <- SiteData[1:nrow(SiteData),"Time"]               
TrueSeq      <- SiteData[1:nrow(SiteData),SMDepthVar]
InvalidSeq   <- SiteData[1:nrow(SiteData),"Invalid"]
EtaSeries    <- AddEtaSeries(ParameterList[i,"C1"], ParameterList[i,"C2"], 
                ParameterList[i,"C3"], SiteData)
BetaSeries   <- GetBeta(PrecipSeries, EtaSeries, zTest, n) 
BetaDiffs    <- GetBetaDiffs(PrecipSeries,EtaSeries,zTest,LongTime,BetaSeries) 
  #Are we doing an adjustment based on elevation?
  if(newC4==0){ C4 = ParameterList[i,"C4"] } else {C4 = newC4}                   
SMSeries     <- GetSMSeries(BetaSeries, ParameterList[i,"Porosity"], ParameterList[i,"ResSM"], C4)
ErrorSeq     <- TrueSeq - SMSeries
#when does the soil moisture data begin and end?
CleanSeq     <- which(TrueSeq!=9999)
#StartInd <- min(which(TrueSeq!=9999))               
#EndInd   <- max(which(TrueSeq!=9999))

#How'd we do?#         
    if(PrecipSite=="IowaState"){
       cor(TrueSeq[CleanSeq][1:408],SMSeries[CleanSeq][1:408])
    }else{
       cor(TrueSeq[CleanSeq],SMSeries[CleanSeq])
    }

#Grab a "training set" from the calibration index
TrainingSet <- GetTrainingSetFromOther(ID, EndPath,ParameterList,zTrain,TSH=0,newC4)
#Reconfigure the SiteData array into a "testing set"
TestingSet <- cbind(TimeStamps, BetaSeries, BetaDiffs, TimeOfDays, TrueSeq, 
                    SMSeries, InvalidSeq, ErrorSeq)  
TestingSet <- TestingSet[CleanSeq,]

k = floor(nrow(TrainingSet)/100)
  
KNNOut <- KNN_SM(IndVars,DepVar,k,TrainingSet,TestingSet,"KNN",Unif=F,Classify=F)
if(DepVar == "ErrorSeq"){
  ES <- KNNOut 
  KNNOut <- ES + TestingSet[,"SMSeries"]
} 
if(PrecipSite=="IowaState"){
  print(cor(KNNOut[1:408],TestingSet[1:408,"TrueSeq"]))
  Output <- cbind(KNNOut[1:408],TestingSet[1:408,"TrueSeq"])
}else{
  print(cor(KNNOut,TestingSet[,"TrueSeq"]))
}
Output <- cbind(KNNOut,TestingSet)
write.csv(Output,paste("C:\\Users\\Evan\\Desktop\\",SoilMoistureSite,ID,MaxInfiltration,".csv",sep=""))

##########################################################################
##########################################################################                                                                          
###This adjusts the results based upon the elevation class some of a brief 
###calibration period ####################################################

#Where did we store all the sensors?
AllSensorsName <- "IowaSensors_LatLonLearning.csv"

#Read in list of all the sensors containing KNN_errors and elevation info
AllSensors <- read.csv(paste(TestSitePath,AllSensorsName,sep=""))

#Training list examples
NumTrainingHours <- 1200 #1200 for Illinois, less for Iowa 
#Last valid day of sensor data (tillage, sensor failure, etc)
LastDay <- 185

#Elevation class of the point to adjust (0 - Pit, 1 - Slope, 2 - Peak)
ElevClass <- 0 

#New independent and dependent variables
IndVarsElev <- c("BetaSeries","KNNOut","SMSeries")
DepVar <- "KNN_Error"

#What is the earliest time stamp and when do we stop training?
FirstTime <- min(AllSensors[,"TimeStamps"])
EndTrain <- FirstTime + NumTrainingHours/24

#Remove ALL data entries after a given date (useful for IowaState, whose data
#becomes invalid after tillage occurs...)
AllSensors <- AllSensors[AllSensors[,"TimeStamps"] < LastDay,]

#Dice up all sensory data to obtain the TrainingSet
#TrainingSet <- AllSensors[AllSensors[,"TimeStamps"]<EndTrain & AllSensors[,"ElevClass"]==ElevClass,]
   #As an experiment try 1/2, odd and even number days (even days for training)
   EvenDays <- floor(AllSensors[,"TimeStamps"])/2 == floor(floor(AllSensors[,"TimeStamps"])/2)
   TrainingSet <- AllSensors[EvenDays & AllSensors[,"ElevClass"]==ElevClass,]

#The existing testing set needs to be paired down...
KNN_Error <- Output[,"TrueSeq"] - Output[,"KNNOut"] 
Output <- cbind(Output,KNN_Error)
Output <- Output[Output[,"TimeStamps"] < LastDay,]

#TestingSet <- Output[Output[,"TimeStamps"]>=EndTrain,]
   #Odd days for testing
   OddDays <- floor(Output[,"TimeStamps"])/2 != floor(floor(Output[,"TimeStamps"])/2)
   TestingSet <- Output[OddDays,]

#NumSimilar matches
k <- floor(nrow(TrainingSet)/1000)
ElevKNN <- KNN_SM(IndVarsElev,DepVar,k,TrainingSet,TestingSet,"KNN",Unif=F,Classify=F)
ElevAdj <- ElevKNN + TestingSet[,"KNNOut"]

#Before and after the elevation adjustment correlations
cor(TestingSet[,"KNNOut"],TestingSet[,"TrueSeq"])
cor(ElevAdj,TestingSet[,"TrueSeq"])
write.csv(ElevAdj,paste("C:\\Users\\Evan\\Desktop\\",SoilMoistureSite,ID,"ElevAdj",".csv",sep=""))
