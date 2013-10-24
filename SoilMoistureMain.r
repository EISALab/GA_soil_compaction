#This is the simplified "main module" for predicting and testing soil moisture

#This set of scripts will ultimately model soil moisture at ~150 locations
#across the USA.
rm(list=ls())

#Place functions into memory that enable soil moisture ML algorithms
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\NationalSoilMoistureFunctions.r")

#Place GA functions into memory
source("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\GeneticAlgorithmCode.r")

#Where to get an uncleaned file
Path <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SCAN_Data_Long\\"

#Path for the storage location of the clean datafiles and parameter list
EndPath <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SCAN_Data_Long_LinConnect_4\\"
 
#Read in the parameter list
ParameterList <- read.csv(paste(EndPath,"ParamList.csv",sep=""))

#################################
###PRE-SET, USER INPUT VALUES####
################################# 

  ####WHAT MODULES NEED TO BE RUN?
  PreCleaned <- T #Are we going to load a cleaned file with missing data handled
  #already, or are we going to need to go through this process from scratch?
  FittingEta <- F #Do we need to fit the sinusoidal parameters for Eta, or 
  #do we already have them?
  FittingDiag <- F #Is the diagnostic soil equation's parameters already fit,
  #or do we need to fit that?
  Validation <- T #Shall we validate the findings after the parameters are fit?
  CrossValidation <- T #Validate one site with information from another
  

  TestingSetHours <- 8760*2 #Leave two years on which we will not calibrate
  SMDepthVar      <- "SMS.I.1..4" 
  
  #Necessary Variables for SM analysis using the method from Pan et al.
  KeyVars <- c("Site.Id","Date","Time","PRCP.H.1..in.",SMDepthVar)
  #,"SMS.I.1..2", "SMS.I.1..4", "SMS.I.1..8", "SMS.I.1..20", "SMS.I.1..40")
                                                   
  #Don't judge fitness if we're lacking these variables                                                                               
  InvalidIfMissing <- c("PRCP.H.1..in.",SMDepthVar)
  
  #"KNN", "ClassTree" (Used for validating results)
  WhatAlgo <- "KNN"
  IndVars  <- c("TimeStamps","TimeOfDays","BetaDiffs","SMSeries","BetaSeries")
  DepVar   <- "ErrorSeq"  #"SoilClass" should be done with Classify=T
  Classify <- F            #"TrueSeq" should be done with Classify=F
  RegType  <- "KNN" #Can be "KNN" or "Kernel"
  k        <- 501 #For convenience, make this odd
  
  InitialPop  <- 20
  Generations <- 1000
  ParamRangesE <- rbind(c(0.1, 0.1, 365),c(0.00001,0.00001,0.01))  #NO ZEROS
  ParamRangesD <- rbind(c(80,50,10),c(0,0,0))
  GenType     <- "LatinHypercube"  #"LatinHypercube" or "MonteCarlo"
  RepType     <- "Two-Parent"      #"Two-Parent" or "Multi-Parent"
  RGen        <- 4
  NumMates    <- 12
  MutProb     <- 0.5
  MutMag      <- 0.03
  DeathRate   <- 0.2
  z           <- 4
  c2constr    <- T #For the Beta modeling, the amplitude of the sinusoid must be
                   #less than the the vertical shift to avoid producing negative
                   #beta values.  This parameter insists upon this  
  
  #To determine an optimal system memory, we need an EtaSeries and a BetaSeries
  #These needn't be optimal, just approximate.  So here are parameters that are,
  #at a minimum, in the ballpark and should allow system memory to be modeled
  DefaultTriad <- c(0.015, 0.01, 235)
  Interval <- 50  #Minimum memory (I set this to 50 hours...)
  LongTime <- 2000 #Long-enough system memory for ANY beta series.
  Thresh   <- 0.97 #How correlated must our shorter beta series be with the 
  #"long" version before we consider the short one to approximate the long one?
  
  #Where should we draw our dividing lines in terms of percentiles?
  QCuts <- c(0.25,0.75)
  
  #Cross valid indices (CalibratedSite, UnCalibratedSite)
  CrossValid <- c(25,25)
  
##############
##############                   
if(CrossValidation == F){
  for(i in 1:nrow(ParameterList)){
    SiteID <- ParameterList[i,"SiteID"]
  
    if(PreCleaned == T){                                                   
       #Read in the site's data file
       SiteData <- GetCleanSCANFile(SiteID, EndPath)
    }else{
       SiteData <- PreProcessor(SiteID,Path,KeyVars,EndPath,InvalidIfMissing)
    }
    
  #Basic site info 
  FirstValid <- ParameterList[i,"FirstValid"]
  #Stick to growing season...
  GrowStart <- 100#ParameterList[i,"Grow_Start"]
  GrowEnd <- 300#ParameterList[i,"Grow_End"]
  PrecipSeries <- SiteData[FirstValid:(nrow(SiteData)-TestingSetHours),"PRCP.H.1..in."]
  TimeStamps   <- SiteData[FirstValid:(nrow(SiteData)-TestingSetHours),"Date"] %% 1000 + 
                       SiteData[FirstValid:(nrow(SiteData)-TestingSetHours),"Time"] - 1
  TimeOfDays   <- SiteData[FirstValid:(nrow(SiteData)-TestingSetHours),"Time"]               
  TrueSeq      <- SiteData[FirstValid:(nrow(SiteData)-TestingSetHours),SMDepthVar]
  InvalidSeq   <- SiteData[FirstValid:(nrow(SiteData)-TestingSetHours),"Invalid"]
  
      
    if(ParameterList[i,"n"] == 0){
      print("Calculating System Memory")
      n <- GetSystemMemory(DefaultTriad,SiteData,FirstValid,TestingSetHours,
                                PrecipSeries,z,LongTime,Interval,Thresh)
      ParameterList[i,"n"] <- n
    }
    #Which process are we accomplishing?
    if(FittingEta){
       #Using GA to fit the three sinusoidal parameters for eta...
       n <- ParameterList[i,"n"]
       BEST <- FullGA(InitialPop, Generations, ParamRangesE, GenType, RepType, 
                      RGen, NumMates, TrueSeq, TimeStamps, MutProb, MutMag, 
                      DeathRate, PrecipSeries, z, n, c2constr, BetaSeries=0, 
                      GrowStart, GrowEnd, InvalidSeq, "FittingEta")
       ParameterList[i,c("C1","C2","C3")] <- BEST[nrow(BEST),1:3]
    }
    
    if(FittingDiag){
       #Given that we've already fit the three eta parameters, fit the three
       #parameters used in the diagnostic soil equation  
       n <- ParameterList[i,"n"]
       EtaSeries    <- AddEtaSeries(ParameterList[i,"C1"], ParameterList[i,"C2"], 
                       ParameterList[i,"C3"], SiteData[FirstValid:(nrow(SiteData)
                       -TestingSetHours),])
       BetaSeries <- GetBeta(PrecipSeries, EtaSeries, z, n)                
       BEST <- FullGA(InitialPop, Generations, ParamRangesD, GenType, RepType, 
                      RGen, NumMates, TrueSeq, TimeStamps, MutProb, MutMag, 
                      DeathRate, PrecipSeries, z, n, c2constr, BetaSeries, 
                      GrowStart, GrowEnd, InvalidSeq, "FittingDiag")
       ParameterList[i,c("Porosity","ResSM","C4")] <- BEST[nrow(BEST),1:3]              
    }
    
    if(Validation){
       #With all six parameters fit, test on the validation data for the purposes
       #of practical decision-support
       n <- ParameterList[i,"n"]
       #The FULL LENGTH versions
       PrecipSeries <- SiteData[FirstValid:nrow(SiteData),"PRCP.H.1..in."]
       TimeStamps   <- SiteData[FirstValid:nrow(SiteData),"Date"] %% 1000 + 
                       SiteData[FirstValid:nrow(SiteData),"Time"] - 1
       TimeOfDays   <- SiteData[FirstValid:nrow(SiteData),"Time"]               
       TrueSeq      <- SiteData[FirstValid:nrow(SiteData),SMDepthVar]
       InvalidSeq   <- SiteData[FirstValid:nrow(SiteData),"Invalid"]
       #Redundant if we've already done in the previous module
       if(FittingDiag == F){
          EtaSeries    <- AddEtaSeries(ParameterList[i,"C1"], ParameterList[i,"C2"], 
                      ParameterList[i,"C3"], SiteData[FirstValid:nrow(SiteData),])
          BetaSeries <- GetBeta(PrecipSeries, EtaSeries, z, n)                   
       }
       
       BetaDiffs <- GetBetaDiffs(PrecipSeries,EtaSeries,z,LongTime,BetaSeries)         
  
       SMSeries <- GetSMSeries(BetaSeries, ParameterList[i,"Porosity"], 
                               ParameterList[i,"ResSM"], ParameterList[i,"C4"])
       ErrorSeq <- TrueSeq - SMSeries
       
       StartInd <- n+1
       EndInd   <- length(TrueSeq)
       EndTrain <- EndInd - TestingSetHours
       #Cutoff everything before the start of growing season and after the end
       GrowingInd  <- which(TimeStamps[StartInd:EndTrain] >= GrowStart &
                           TimeStamps[StartInd:EndTrain] <= GrowEnd)
       MetaSiteData <- cbind(TimeStamps,BetaSeries,BetaDiffs,
                             TimeOfDays,TrueSeq,SMSeries,InvalidSeq,ErrorSeq)
       SnippedInvalid <- InvalidSeq[StartInd:EndTrain][GrowingInd]
       
       #This is where the "question" is asked...where the "classes" are defined.
       #Users will specify the boundaries between the classes they would prefer.                      
       MetaSiteData <- AddSoilClass(TrueSeq,StartInd,EndTrain,GrowingInd,
                                    QCuts,SnippedInvalid,MetaSiteData)
       
       #Break into training and testing datasets
       TrainIndices <- seq(from=StartInd,to=EndTrain)
       TestIndices  <- seq(from=(EndTrain+1),to=EndInd)
       TrainData    <- MetaSiteData[TrainIndices,]                
       TestData     <- MetaSiteData[TestIndices,]
       #Filter out any examples that are not within the appropriate growing days
       TestGInd     <-  which(TestData[,"TimeStamps"] >= GrowStart &
                        TestData[,"TimeStamps"] <= GrowEnd)   
       TestInvalid  <- InvalidSeq[TestIndices][TestGInd]                 
       TrainingSet  <- TrainData[GrowingInd,][!SnippedInvalid,]
       TestingSet   <- TestData[TestGInd,][!TestInvalid,]             
                                           
       if(WhatAlgo == "KNN"){
          KNNOut <- KNN_SM(IndVars,DepVar,k,TrainingSet,TestingSet,RegType,Unif=F,Classify)
          if(DepVar == "ErrorSeq"){
             ES <- KNNOut 
             KNNOut <- ES + TestingSet[,"SMSeries"]
          } 
         #Evaluation Battery
          if(Classify){
            KN <- as.numeric(KNNOut[,"DepEstimates"])
            CLASSES <- unique(TrainingSet[,DepVar])
            AccGrid <- matrix(NA,length(CLASSES),length(CLASSES))        
            for(m in 1:max(CLASSES)){
              for(n in 1:max(CLASSES)){
                 AccGrid[m,n] <- sum(KN==m & TestingSet[,DepVar]==n)
              }
            }        
            #Allows further scrutiny of results in excel
            Output <- cbind(KN,as.numeric(KNNOut[,"K_Counts"]),TestingSet)
            write.csv(Output,"C:\\Users\\Evan\\Desktop\\Out.csv")
            print(cor(TestingSet[,"SMSeries"],TestingSet[,"TrueSeq"]))
          }else{
            print(cor(KNNOut,TestingSet[,"TrueSeq"]))
            Output <- cbind(KNNOut,TestingSet)
            write.csv(Output,"C:\\Users\\Evan\\Desktop\\Out.csv")
          }
       }
       if(WhatAlgo == "ClassTree"){
           #Not implemented...
       }                                                                      
    }
  }                
}
if(CrossValidation==T){
  TrainI <- CrossValid[1]
  TestI  <- CrossValid[2]
  SiteIDTrain <- ParameterList[TrainI,"SiteID"]
  #Read in the site's data file
  SiteDataTrain <- GetCleanSCANFile(SiteIDTrain, EndPath)
  
  SiteIDTest <- ParameterList[TestI,"SiteID"]
  #Read in the site's data file
  SiteDataTest <- GetCleanSCANFile(SiteIDTest, EndPath)
  ########################################################### 
  TrainingSet <- GetTrainingOrTestingSet(ParameterList,TrainI,SiteDataTrain,"Train")
  TestingSet  <- GetTrainingOrTestingSet(ParameterList,TrainI,SiteDataTest,"Test")
  
  #We won't know the optimal k, so we'll just set it based on the size of the 
  #training set... 
  k = floor(nrow(TrainingSet)/100)
  
  KNNOut <- KNN_SM(IndVars,DepVar,k,TrainingSet,TestingSet,RegType,Unif=F,Classify)
  if(DepVar == "ErrorSeq"){
    ES <- KNNOut 
    KNNOut <- ES + TestingSet[,"SMSeries"]
  } 
  print(cor(KNNOut,TestingSet[,"TrueSeq"]))
  Output <- cbind(KNNOut,TestingSet)
  write.csv(Output,"C:\\Users\\Evan\\Desktop\\Out.csv")
}  
