rm(list=ls())

source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\MLCompendium_EC2.r")
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\AssortedSoilFunctions.r")

##Finding the datafile for the Champaign NOAA sensor, located at Willard Airport
FileName <- "UNIVERSI OF IL WILLARD APT (94870)"

GetNOAAData <- function(FileName){
  ##Reads in a NOAA data file
  DATA <- read.csv(paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\NOAA_DATA\\",FileName,".csv",sep=""))
  DATA
}

##Read in the data file
DATA <- GetNOAAData(FileName)

FirstDate <- 40179
ConvertTime <- function(DATA,FirstDate){

}

