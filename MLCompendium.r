#This file is a compendium of a variety of machine learning algorithms.
#Each function will begin with a matrix of one dependent variable and eight
#independent variables.  The user will select which of the eight independent
#variables should be included and the model which will be used.

#Read in the data
   GetData <- function(FileName){
     DATA <- read.csv(paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\",FileName,".csv",sep=""))
     DATA
   }

 ########################   KNN   #####################################################
IndVars <- c("X1","X2","X3")
DepVar <- "Y"
k <- 5
HistEndDate <- 100
EndForecastDate <- 150
RegType <- "KNN"
Unif <- F

 #The generic version of KNN. For historical record included, this function
 #should return an estimate of the dependent variables
    KNN <- function(IndVars,DepVar,k,DATA,HistEndDate,EndForecastDate,RegType,Unif=F){
           #This counts from the end of the historical set to the end of the predicted set
           #NumRows will equal the number of data points for which we are making a prediction
           NumRows                 <- sum(DATA[,"Date"]<EndForecastDate)-sum(DATA[,"Date"]<HistEndDate)
           #The predictions made
           DepEstimates            <- array(NA,NumRows)
           #A relic from old baseball code
           #MusAndSigmas            <- matrix(NA,NumRows,(2*length(IndVars)+1))
           HistoricalSet <- DATA[DATA[,"Date"]<=HistEndDate,c("Date",DepVar,IndVars)]

           for(i in 1:NumRows){
            print(i)
            #The independent variables for the data point to be predicted
            X_i <- DATA[i+nrow(HistoricalSet),IndVars]
            #Call GetDistances, locate the nearest neighbors
            DistanceSet   <- GetDistances(X_i,HistoricalSet[,c(DepVar,IndVars)],IndVars,k,RegType,Unif)

            #A relic from old baseball code to record properties of data
            #MusAndSigmas[i,1]      <-  DATA[i,"Date"]
            #MusAndSigmas[i,2:(length(IndVars)+1)]    <-  colMeans(HistoricalSet[,IndVars])
            #MusAndSigmas[i,(length(IndVars)+2):ncol(MusAndSigmas)]   <-  apply(HistoricalSet[,IndVars],2,sd)

            #Average our nearest neighbors - obtain our estimate, then output
            DepEstimates[i]           <- mean(DistanceSet[,DepVar])        }
            DepEstimates
    }

#Given a matrix of independent variables (and a dependent variable which will be snipped),
#we compute a distance value (S_ij similarity) for every point in the historical database.
GetDistances <- function(X_i,HistoricalSet,Vars,k,RegType,Unif=F){

      DistanceCol    <- array(NA,nrow(HistoricalSet))
      Sigmas         <- apply(HistoricalSet[,Vars],2,sd)
      Mus            <- colMeans(HistoricalSet[,Vars])
      MeanZeroMatrix <- t(t(HistoricalSet[,Vars])-Mus)
      NormMatrix     <- t(t(MeanZeroMatrix[,Vars])/Sigmas)

      #If we would prefer to work with a uniform dist rather than a multivariate Gaussian...
      if(Unif){
        NormMatrix <- pnorm(NormMatrix)
      }

      NormX_i     <- as.numeric((X_i-Mus)/Sigmas)
      DiffMatrix  <- t(t(NormMatrix[,Vars])-NormX_i)
      DiffMatrix  <- DiffMatrix^2

      DistanceCol   <- sqrt(rowSums(DiffMatrix))
      HistoricalSet <- cbind(HistoricalSet,DistanceCol)
      HistoricalSet <- HistoricalSet[order(DistanceCol),]
      if(RegType=="KNN"){
        Matches       <- HistoricalSet[1:k,]
      }
      if(RegType=="Kernel"){
        Matches       <- HistoricalSet[HistoricalSet[,"DistanceCol"]<k,]
      }
      Matches
    }

########################   NEURAL NETS   #####################################################

IndVars <- c("X1","X2","X3")
DepVar <- "Y"
OutputNodes <- 2
HiddenNodes <- 5
HistEndDate <- 100
EndForecastDate <- 150
Gamma <- 0.1
Alpha <- 0

#HONN models require combinatoric terms...this function will enumerate all needed combinations
  GetCombinations <- function(N){
    CombinationMat <- matrix(NA,(factorial(N)/(factorial(N-2)*factorial(2))),2)
    #This will be our counter over our rows of CombinationMat
    CurrentRow <- 1
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        CombinationMat[CurrentRow,] <- c(i,j)
        CurrentRow <- CurrentRow + 1
      }
    }
    CombinationMat
  }

#This function is the main for the Standard Neural Netword
  NeuralNet <- function(IndVars, DepVar, DATA, OutputNodes, HiddenNodes, HistEndDate, EndForecastDate, Gamma, Alpha, InitializeWeights=T, ITHWs, HTOWs){
    #Either all weights begin with a value of one, or they are given as inputs
    if(InitializeWeights==T){
      #Weights are ordered as followed: 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3, 1.4, 2.4, 3.4, 1.5, 2.5, 3.5
      #a.b represents the link from node a in the input (hidden) layer to node
      #b in the hidden (output) layer.
      InputToHiddenWeights  <- array(1,(length(IndVars)*HiddenNodes))
      HiddenToOutputWeights <- array(1,(HiddenNodes*OutputNodes))
    }else{
      InputToHiddenWeights  <- ITHWs
      HiddenToOutputWeights <- HTOWs
    }
  }
  
  #Using the inputs from the input nodes and the relevant weights, use sigmoid
  #functions to determine the activation of each hidden node.
  GetHiddenNodeValues <- function(InputToHiddenWeights, Inputs, HiddenNodes){
    HiddenNodeValues <- array(NA,HiddenNodes)
    LinksPerHiddenNode <- length(Inputs)
    for(i in 1:HiddenNodes){
       SumOfAllWeights <- sum(InputToHiddenWeights[(LinksPerHiddenNode*(i-1)+1):(LinksPerHiddenNode*i)]*Inputs)
       HiddenNodeValues[i] <- Sigmoid(SumOfAllWeights)
    }
    HiddenNodeValues
  }
  
  ##Using the values from the hidden nodes and the relevant weights, use sigmoid
  ##functions to determine the activation of each output node.
  GetOutputNodeValues <- function(HiddenToOutputWeights, HiddenNodeValues, OutputNodes){
    OutputNodeValues <- array(NA,OutputNodes)
    LinksPerOutputNode <- length(HiddenNodeValues)
    for(i in 1:OutputNodes){
       SumOfAllWeights <- sum(HiddenToOutputWeights[(LinksPerOutputNode*(i-1)+1):(LinksPerOutputNode*i)]*HiddenNodeValues)
       OutputNodeValues[i] <- Sigmoid(SumOfAllWeights)
    }
    OutputNodeValues
  }
  
  #Using the output node values along with the 'true' values, determine the
  #values for delta_o
  GetDeltaOutputValues <- function(OutputNodeValues, ActualValues){
    DeltaOutputValues <- (ActualValues - OutputNodeValues)*OutputNodeValues*(1-OutputNodeValues)
  }
  
  #Using the hidden node values, and the previously calculated delta output
  #values, determine the delta values for the hidden layer
  GetDeltaHiddenValues <- function(HiddenNodeValues,DeltaOutputValues,HiddenToOutputWeights){
    DeltaHiddenValues <- array(NA,length(HiddenNodeValues))
    for(i in 1:length(HiddenNodeValues)){
      SumProduct <- 0
      for(j in 1:length(DeltaOutputValues)){
        SumProduct <- SumProduct + DeltaOutputValues[j]*HiddenToOutputWeights[(1+length(HiddenNodeValues)*(j-1))]
      }
      DeltaHiddenValues[i] <- HiddenNodeValues[i]*(1-HiddenNodeValues[i])*SumProduct
    }
    DeltaHiddenValues
  }
  
  #With DeltaOutputValues, Update the InputToHiddenWeights
  UpdateInputToHiddenWeights <- function(InputToHiddenWeights,DeltaHiddenValues,Gamma,HiddenNodeValues){
    for(i in 1:length(InputToHiddenWeights){
      HiddenNodeInQuestion <- floor(i/length(DeltaHiddenValues))+1
      InputToHiddenWeights[i] <- InputToHiddenWeights[i] + Gamma*DeltaHiddenValues[HiddenNodeInQuestion]*HiddenNodeValues[HiddenNodeInQuestion]
    }
    InputToHiddenWeights
  }
  
  #Compute a simple sigmoid function
  Sigmoid <- function(Sp){
    Sig <- 1/(1+exp(Sp))
    Sig
  }