#Read in all soil moisture series at a given location, merge them by date, then
#integrate with existing LiDAR classifications.

rm(list=ls())
#Place functions into memory that enable soil moisture ML algorithms
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\NationalSoilMoistureFunctions.r")

#Acceptable NEXRAD names, ("Illinois", "IowaState", "Purdue")
PrecipSite <- "Illinois"
#Acceptable folder names ("University.of.Illinois", Purdue", "Iowa.State")
SMFolder <- "University.of.Illinois"
#Path to NEXRAD data and soil moisture data
TestSitePath <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\CropSense Files 10-3\\"
CoordinatesBySitePath <- "SiteList_w_LatLon"
ID <- 2036
zTrain <- 4
zTest <- 4
Default <- 0
LongTime <- 2000 #Long enough for ANY beta series
EndPath <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SCAN_Data_Long_LinConnect_4\\"

SensorsAtSite <- GetFile("ListBySite",TestSitePath)
Sensors <- SensorsAtSite[,SMFolder]
Sensors <- Sensors[which(Sensors!="XXXX")]

for(i in 1:length(Sensors)){
    SensorName <- as.character(Sensors[[i]])
    print(SensorName)
    #Get the soil moisture file
    SMData <- GetSMFile(TestSitePath,SMFolder,SensorName)
    #Get the latitude and longitude of the soil moisture site
    LatLong <- GetSMLatLon(CoordinatesBySitePath,TestSitePath, SensorName)
    SiteLat <- LatLong[1]
    SiteLon <- LatLong[2]
    ProcessedData <- SensorPrecipAndSM(ID, zTrain, zTest, Default, EndPath,
                                       LongTime, PrecipSite, SiteLat, SiteLon)
    colnames(ProcessedData) <- paste(colnames(ProcessedData),SensorName,sep="_")
    if(i == 1){
      AllSensors <- ProcessedData
    }else{
      AllSensors <- cbind(AllSensors, ProcessedData[,3:ncol(ProcessedData)])
    }
}

write.csv(AllSensors,paste(EndPath,"IllinoisSensors",sep="_"))

