#This file is a compendium of a variety of machine learning algorithms.
#Each function will begin with a matrix of one dependent variable and eight
#independent variables.  The user will select which of the eight independent
#variables should be included and the model which will be used.

#Read in the data
   GetData <- function(FileName){
     DATA <- read.csv(paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\",FileName,".csv",sep=""))
     DATA
   }

#Analog to previous function, reads EBI soil data
    GetEBIData <- function(FileName){
      DATA <- read.csv(paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\EBI-data\\",FileName,".csv",sep=""))
      DATA
    }
   
#Normalize a matrix, set every column to mean=0, sd=1
   Normalize <- function(DATA,Vars){
      Sigmas         <- apply(DATA[,Vars],2,sd)
      Mus            <- colMeans(DATA[,Vars])
      MeanZeroMatrix <- t(t(DATA[,Vars])-Mus)
      NormMatrix     <- t(t(MeanZeroMatrix[,Vars])/Sigmas)
      NormMatrix
   }

######################################################################################
########################   KNN   #####################################################
######################################################################################
IndVars <- c("X1","X2","X3")
DepVar <- "Y"
k <- 5
HistEndDate <- 100
EndForecastDate <- 150
RegType <- "KNN"
Unif <- F

 #The generic version of KNN. For historical record included, this function
 #should return an estimate of the dependent variables
    KNN <- function(IndVars,DepVar,k,DATA,HistEndDate,EndForecastDate,RegType,Unif=F,Classify=F){
           #This counts from the end of the historical set to the end of the predicted set
           #NumRows will equal the number of data points for which we are making a prediction
           NumRows                 <- sum(as.numeric(DATA[,"Date"])<EndForecastDate)-sum(as.numeric(DATA[,"Date"])<HistEndDate)
           #The predictions made
           DepEstimates            <- array(NA,NumRows)
           if(Classify){
              K_Counts             <- array(NA,NumRows)
           }
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
            if(Classify==F){
               DepEstimates[i]           <- mean(DistanceSet[,DepVar])        
            }else if(Classify==T){
               DepEstimates[i]           <- names(sort(-table(DistanceSet[,DepVar])))[1]
               K_Counts[i]               <- sum(DistanceSet[,DepVar]=="READY")
            }
     
            }
            if(Classify){
               DepEstimates <- cbind(DepEstimates,K_Counts)
            }
            DepEstimates
    }

PrepKNN <- function(IndVars,DepVar,k,DATA,HistEndDate,EndForecastDate,RegType,Unif=F,Classify=F){
#Place historical set into memory appropriately for KNN use
           HistoricalSet <- DATA[DATA[,"Date"]<=HistEndDate,c("Date",DepVar,IndVars)]
           HistoricalSet  
}

TestConditionsKNN <- function(HistoricalSet,IndVars,DepVar,k,RegType,Unif=F,X_i,Classify=T){
            #Call GetDistances, locate the nearest neighbors
            DistanceSet   <- GetDistances(X_i,HistoricalSet[,c(DepVar,IndVars)],IndVars,k,RegType,Unif)

            #Average our nearest neighbors - obtain our estimate, then output
            if(Classify==F){
               DepEstimate           <- mean(DistanceSet[,DepVar])        
            }else if(Classify==T){
               DepEstimate           <- names(sort(-table(DistanceSet[,DepVar])))[1]
               K_Count               <- sum(DistanceSet[,DepVar]=="READY")
            }
            if(Classify){
               DepEstimate <- cbind(DepEstimate,K_Count)
            }
            DepEstimate
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
    
##############################################################################################
########################   NEURAL NETS   #####################################################
##############################################################################################

IndVars <- c("X1","X2","X3")
DepVar <- "Y"
OutputNodes <- 2
HiddenNodes <- 5
DATA <- GetData("MyRawData")
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
  
#This function trains a higher order neural network
#Note - it is important to be cautious with the relative scales of the variables
#especially with respect to the combinatoric terms.  If there multiplicative terms
#grow too large, our sigmoid function is ALWAYS going to set values to 0 or 1 
#and prevent our weights from changing appropriately
  TrainHONN <- function(IndVars,DepVar,Threshholds,DATA,OutputNodes,HistEndDate,LearningRate,InitializeWeights=T,Iterations,ITOWs=NULL,Alpha=0){
    #Either all weights begin with a value of one, or they are given as inputs
    CombinationMat <- GetCombinations(length(IndVars))
    if(InitializeWeights==T){
        InputToOutputWeights  <- matrix(0,(length(IndVars)+nrow(CombinationMat)),OutputNodes)
     }else{
      InputToOutputWeights  <- ITOWs
    }
    
    HistoricalSet <- DATA[DATA[,"Date"]<=HistEndDate,c("Date",DepVar,IndVars)]
    #Add the combinatoric columns to the historical set
    CHistoricalSet <- BuildCombinatoricHistoricalSet(CombinationMat,HistoricalSet)

    #Iterate over the training data             
    for(i in 1:Iterations){
      print(paste("Iteration #",i))                                                                                 
      for(j in 1:nrow(HistoricalSet)){                                                    
        print(paste("Data Point #",j))
        #This is our independent variables + the combinatoric terms
        AllIndVars <- as.numeric(CHistoricalSet[j,(3:(ncol(CHistoricalSet)))])
        #Determine values for output nodes                                    
        OutputNodeValues <- GetOutputNodeValues(InputToOutputWeights,AllIndVars)
        #Get actual values 
        ActualValues <- GetActualValues(CHistoricalSet[j,DepVar],Threshholds)
        #Get delta_output values
        DeltaOutputValues <- GetDeltaOutputValues(OutputNodeValues,ActualValues)
        #Update the weights between hidden nodes and output nodes
        InputToOutputWeights <- UpdateHiddenToOutputWeights(InputToOutputWeights,DeltaOutputValues,LearningRate,OutputNodeValues,AllIndVars)
      }
    }
    InputToOutputWeights
  } 
  
#This function is the main for the standard HONN testing.
 TestHONN <- function(DATA, StartDate, EndDate, IndVars, DepVar, Threshholds, OutputNodes, InputToOutputWeights, Train=FALSE, LearningRate=NULL, Alpha=0){
        TestingSet <- DATA[as.numeric(DATA[,"Date"])<=EndDate & as.numeric(DATA[,"Date"])>=StartDate,c("Date",DepVar,IndVars)]
        #Add the combinatoric columns to the testing set
        CTestingSet <- BuildCombinatoricHistoricalSet(CombinationMat,TestingSet)
    
        NumPredictions <- nrow(TestingSet)
        HONN_Forecasts <- array(NA,NumPredictions) 
        for(j in 1:NumPredictions){
            print(paste("DataPoint #: ",j))
            #This is our independent variables + the combinatoric terms
            AllIndVars <- as.numeric(CTestingSet[j,(3:(ncol(CTestingSet)))])
            #Determine values for output nodes                                    
            OutputNodeValues <- GetOutputNodeValues(InputToOutputWeights,AllIndVars)
            #Each forecast is an integer between 1 and the number of output nodes
            HONN_Forecasts[j] <- max(which(OutputNodeValues==max(OutputNodeValues))) 
            if(Train==T){
                #Get actual values 
                ActualValues <- GetActualValues(CTestingSet[j,DepVar],Threshholds)
                #Get delta_output values
                DeltaOutputValues <- GetDeltaOutputValues(OutputNodeValues,ActualValues)
                #Update the weights between input nodes and output nodes
                InputToOutputWeights <- UpdateHiddenToOutputWeights(InputToOutputWeights,DeltaOutputValues,LearningRate,OutputNodeValues,AllIndVars)
            }
        }
        NN_Forecasts
 }   
  
  
#This function takes a historical set and a combination matrix and assembles
#the full-size historical set.
  BuildCombinatoricHistoricalSet <- function(CombinationMat,HistoricalSet){
    for(i in 1:nrow(CombinationMat)){
      NewColumn <- HistoricalSet[,(CombinationMat[i,1]+2)]*HistoricalSet[,(CombinationMat[i,2]+2)]
      HistoricalSet <- cbind(HistoricalSet,NewColumn)
    }
    HistoricalSet 
  }
   

### TEST <- TrainNeuralNet(IndVars,DepVar,100,DATA,2,5,200,0.1,T,2)
#This function is the main for the standard Artificial Neural Network training.
  TrainNeuralNet <- function(IndVars, DepVar, Threshholds, DATA, OutputNodes, HiddenNodes, HistEndDate, LearningRate, InitializeWeights=T, Iterations, ITHWs=NULL, HTOWs=NULL, Alpha=0){
    #Either all weights begin with a value of one, or they are given as inputs
    if(InitializeWeights==T){
      InputToHiddenWeights  <- matrix(0,length(IndVars),HiddenNodes)
      HiddenToOutputWeights <- matrix(0,HiddenNodes,OutputNodes)                      
    }else{
      InputToHiddenWeights  <- ITHWs
      HiddenToOutputWeights <- HTOWs
    }
    
    #Grab the segment of our dataset used to train the network
    HistoricalSet <- DATA[DATA[,"Date"]<=HistEndDate,c("Date",DepVar,IndVars)]
    #Iterate over the training data             
    for(i in 1:Iterations){
      print(paste("Iteration #",i))                                                                                 
      for(j in 1:nrow(HistoricalSet)){                                                    
        print(paste("Data Point #",j))
        #Determine values for hidden nodes
        HiddenNodeValues <- GetHiddenNodeValues(InputToHiddenWeights,HistoricalSet[j,IndVars])
        #Determine values for output nodes                                    
        OutputNodeValues <- GetOutputNodeValues(HiddenToOutputWeights,HiddenNodeValues)
        #Get actual values 
        ActualValues <- GetActualValues(HistoricalSet[j,DepVar],Threshholds)
        #Get delta_output values
        DeltaOutputValues <- GetDeltaOutputValues(OutputNodeValues,ActualValues)
        #Get delta_hidden values
        DeltaHiddenValues <- GetDeltaHiddenValues(HiddenNodeValues,DeltaOutputValues,HiddenToOutputWeights)
        #Update the weights between hidden nodes and output nodes
        HiddenToOutputWeights <- UpdateHiddenToOutputWeights(HiddenToOutputWeights,DeltaOutputValues,LearningRate,OutputNodeValues,HiddenNodeValues)
        #Update the weights between input nodes and hidden nodes
        InputToHiddenWeights <- UpdateInputToHiddenWeights(InputToHiddenWeights,DeltaHiddenValues,LearningRate,HiddenNodeValues, HistoricalSet[j,IndVars])
      }
    }
    FinalWeights <- c(as.numeric(InputToHiddenWeights),as.numeric(HiddenToOutputWeights))
    FinalWeights
  }
  
 #This function is the main for the standard Artificial Neural Network testing.
 TestNeuralNetwork <- function(DATA, StartDate, EndDate, IndVars, DepVar, Threshholds, HiddenNodes, OutputNodes, AllWeights, Train=FALSE, LearningRate=NULL, Alpha=0){
        #Grab the relevant pieces from AllWeights
        InputNodes <- length(IndVars)
        InputToHiddenWeights <- FormInputToHiddenWeights(AllWeights,InputNodes,HiddenNodes)
        HiddenToOutputWeights <- FormHiddenToOutputWeights(AllWeights,InputNodes,HiddenNodes,OutputNodes)
        
        TestingSet <- DATA[as.numeric(DATA[,"Date"])<=EndDate & as.numeric(DATA[,"Date"])>=StartDate,c("Date",DepVar,IndVars)]
        NumPredictions <- nrow(TestingSet)
        NN_Forecasts <- array(NA,NumPredictions) 
        for(j in 1:NumPredictions){
            print(paste("DataPoint #: ",j))
            #Determine values for hidden nodes
            HiddenNodeValues <- GetHiddenNodeValues(InputToHiddenWeights,TestingSet[j,IndVars])
            #Determine values for output nodes                                    
            OutputNodeValues <- GetOutputNodeValues(HiddenToOutputWeights,HiddenNodeValues)
            #Each forecast is an integer between 1 and the number of output nodes
            NN_Forecasts[j] <- max(which(OutputNodeValues==max(OutputNodeValues))) 
            if(Train==T){
                #Get actual values 
                ActualValues <- GetActualValues(TestingSet[j,DepVar],Threshholds)
                #Get delta_output values
                DeltaOutputValues <- GetDeltaOutputValues(OutputNodeValues,ActualValues)
                #Get delta_hidden values
                DeltaHiddenValues <- GetDeltaHiddenValues(HiddenNodeValues,DeltaOutputValues,HiddenToOutputWeights)
                #Update the weights between hidden nodes and output nodes
                HiddenToOutputWeights <- UpdateHiddenToOutputWeights(HiddenToOutputWeights,DeltaOutputValues,LearningRate,OutputNodeValues,HiddenNodeValues)
                #Update the weights between input nodes and hidden nodes
                InputToHiddenWeights <- UpdateInputToHiddenWeights(InputToHiddenWeights,DeltaHiddenValues,LearningRate,HiddenNodeValues, HistoricalSet[j,IndVars])
            }
        }
        NN_Forecasts
 } 
 
 #Using the string of all network weights, return InputToHiddenWeights
 FormInputToHiddenWeights <- function(AllWeights,InputNodes,HiddenNodes){
   InputToHiddenWeights <- matrix(NA,InputNodes,HiddenNodes)
   for(i in 1:InputNodes){
     for(j in 1:HiddenNodes){
       InputToHiddenWeights[i,j] <- AllWeights[i+(j-1)*InputNodes] 
     }
   } 
   InputToHiddenWeights
 }
 
 #Using the string of all network weights, return InputToHiddenWeights
 FormHiddenToOutputWeights <- function(AllWeights,InputNodes,HiddenNodes,OutputNodes){
   ExtraNodes <- InputNodes * HiddenNodes 
   HiddenToOutputWeights <- matrix(NA,HiddenNodes,OutputNodes)
   for(i in 1:HiddenNodes){
     for(j in 1:OutputNodes){
       HiddenToOutputWeights[i,j] <- AllWeights[ExtraNodes+i+(j-1)*HiddenNodes] 
     }
   } 
   HiddenToOutputWeights
 }  
  #Using the inputs from the input nodes and the relevant weights, use sigmoid
  #functions to determine the activation of each hidden node.
  GetHiddenNodeValues <- function(InputToHiddenWeights, Inputs){
    HiddenNodeValues <- array(NA,ncol(InputToHiddenWeights))
    for(j in 1:length(HiddenNodeValues)){
       Z_j <- sum(InputToHiddenWeights[,j]*Inputs)
       HiddenNodeValues[j] <- Sigmoid(Z_j)
    }
    HiddenNodeValues
  }
  
  ##Using the values from the hidden nodes and the relevant weights, use sigmoid
  ##functions to determine the activation of each output node.
  GetOutputNodeValues <- function(HiddenToOutputWeights, HiddenNodeValues){
    OutputNodeValues <- array(NA,ncol(HiddenToOutputWeights))
    for(j in 1:length(OutputNodeValues)){
       Z_j <- sum(HiddenToOutputWeights[,j]*HiddenNodeValues)
       OutputNodeValues[j] <- Sigmoid(Z_j)
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
       DeltaHiddenValues[i] <- HiddenNodeValues[i]*(1-HiddenNodeValues[i])*sum(HiddenToOutputWeights[i,]*DeltaOutputValues)  
    }
    DeltaHiddenValues
  }
       
  #With DeltaHiddenValues, Update the InputToHiddenWeights
  UpdateInputToHiddenWeights <- function(InputToHiddenWeights,DeltaHiddenValues,LearningRate,HiddenNodeValues,Inputs){
    LR <- as.numeric(LearningRate)
    for(i in 1:nrow(InputToHiddenWeights)){
      for(j in 1:ncol(InputToHiddenWeights)){
         InputToHiddenWeights[i,j] <- InputToHiddenWeights[i,j] + LR*DeltaHiddenValues[j]*as.numeric(Inputs[i])#*InputToHiddenWeights[i,j]
      }
    }
    InputToHiddenWeights
  }
  #With DeltaOutputValues, Update the HiddenToOutputWeights
  UpdateHiddenToOutputWeights <- function(HiddenToOutputWeights,DeltaOutputValues,LearningRate,OutputNodeValues,HiddenNodeValues){
    LR <- as.numeric(LearningRate)
    for(i in 1:nrow(HiddenToOutputWeights)){
       for(j in 1:ncol(HiddenToOutputWeights)){
          HiddenToOutputWeights[i,j] <- HiddenToOutputWeights[i,j] + LR*DeltaOutputValues[j]*HiddenNodeValues[i]#*HiddenToOutputWeights[i,j] 
       }
    }
    HiddenToOutputWeights
  }
  
  #Compute a simple sigmoid function
  Sigmoid <- function(Sp){
    Sig <- 1/(1+exp(-Sp))
    Sig
  }
  
  #Depending on the number of classifications in the neural network, different
  #functions will be required to determine the 'correct' answer to which the
  #network must train itself.  For instance, if there are only two possible 
  #outputs, "ready"/"not ready," then there will be only two output nodes, and
  #the desired result will be (0,1) or (1,0).  However, if an additional 
  #classification is possible, i.e. "ready"/"somewhat ready"/"not ready," then
  #desired output would be (0,0,1) or (0,1,0) or (1,0,0).  This function will
  #construe the dependent variable and determine this appropriate output
  GetActualValues <- function(Result, Threshholds){
     #If we have n output nodes, n-1 threshholds will be required
     #We will sort them highest to lowest
     if(length(Threshholds)>1){
       Threshholds <- Threshholds[order(Threshholds,decreasing=TRUE)]
     }
     ActualValues <- array(NA,(1+length(Threshholds)))
     Found <- 0  #Determines if a threshhold has previous been breached
     for(i in 1:length(Threshholds)){
        if(i < length(Threshholds)){
            if(Result > Threshholds[i] & Found == 0){
              ActualValues[i] <- 1
              Found = 1
            }else{
              ActualValues[i] <- 0  
            }
        }else{
            if(Result > Threshholds[i] & Found == 0){
              ActualValues[i] <- 1
              ActualValues[i+1] <- 0
            }else if(Result < Threshholds[i] & Found == 0){
              ActualValues[i] <- 0
              ActualValues[i+1] <- 1
            }else{
              ActualValues[i] <- 0
              ActualValues[i+1] <- 0          
            }
            
        }    
     }
     ActualValues
  }
  
  ##########################################################################################
  #################### REGRESSION TREES ####################################################
  ##########################################################################################
  
  TrainAndTestRegTrees <- function(IndVars, DepVar, DATA, MinNodeSize, MinSplitSize, HistEndDate, EndTestDate, Classify=T){
      #Break data into the set used to build the tree and the test set
      HistoricalSet <- DATA[as.numeric(DATA[,"Date"])<=HistEndDate,c("Date",DepVar,IndVars)] 
      TestingSet <- DATA[as.numeric(DATA[,"Date"])<=EndTestDate & as.numeric(DATA[,"Date"])> HistEndDate,c("Date",DepVar,IndVars)]

      ###Using the rpart library, construct and manipulate regression trees
      library(rpart)             
      #minbucket determines the minimum size of a terminal node, 
      #minsplit determines the minimum size for which to attempt a split
      #Once we make the tree, the object can be saved and stored
      RTree <- rpart(HistoricalSet[,DepVar] ~ .,HistoricalSet[,IndVars], minbucket = MinNodeSize, minsplit = MinSplitSize)
      plot(RTree,uniform=T,branch=1, margin=0.1, cex=0.9)
      text(RTree,cex=0.75)
      
      #Our forecasts
      RTreeForecasts <- array(NA,nrow(TestingSet))
      #Send each set of inputs into the tree to generate a prediction
      for(i in 1:nrow(TestingSet)){
         #print(paste("Data point #",i))
         NewData = list()
         for(j in 1:length(IndVars)){                                                          
           NewData[[IndVars[j]]] = TestingSet[i,IndVars[j]]
         }
         if(Classify==F){
            RTreeForecasts[i] <- predict(RTree,NewData) 
         }else{
            Predictions       <- predict(RTree,NewData)
            Classification    <- which(predict(RTree,NewData)==max(predict(RTree,NewData)))
            RTreeForecasts[i] <- colnames(Predictions)[Classification]
         }
      }
  RTreeForecasts
      #Prune specific nodes:
      #RTreeSnipped <- snip.rpart(RTree,c(4,7))
  }

PutTreeInMemory <- function(IndVars, DepVar, DATA, MinNodeSize, MinSplitSize, HistEndDate, EndTestDate, Classify=T){
#Returns a regression tree into memory for later use
      #Break data into the set used to build the tree and the test set
      HistoricalSet <- DATA[as.numeric(DATA[,"Date"])<=HistEndDate,c("Date",DepVar,IndVars)] 
      TestingSet <- DATA[as.numeric(DATA[,"Date"])<=EndTestDate & as.numeric(DATA[,"Date"])> HistEndDate,c("Date",DepVar,IndVars)]

      ###Using the rpart library, construct and manipulate regression trees
      library(rpart)             
      #minbucket determines the minimum size of a terminal node, 
      #minsplit determines the minimum size for which to attempt a split
      #Once we make the tree, the object can be saved and stored
      RTree <- rpart(HistoricalSet[,DepVar] ~ .,HistoricalSet[,IndVars], minbucket = MinNodeSize, minsplit = MinSplitSize)
      plot(RTree,uniform=T,branch=1, margin=0.1, cex=0.9)
      text(RTree,cex=0.75)
      RTree
}

TestConditionsRTree <- function(RTree,IndVars,Vals,Classify=T){
#Tests one condition using a regression tree stored in memory
    NewData = list()
    for(j in 1:length(IndVars)){                                                          
        NewData[[IndVars[j]]] = Vals[j]
    }
    if(Classify==F){
        RTreeForecast <- predict(RTree,NewData) 
    }else{
      Predictions       <- predict(RTree,NewData)
      Classification    <- which(predict(RTree,NewData)==max(predict(RTree,NewData)))
      RTreeForecast <- colnames(Predictions)[Classification]
      RTreeForecast
    }
}   

#########BOOSTED PERCEPTRON############
#######################################

BoostedPerceptron <- function(DATA, IndVars, DepVar, Thresh, LastTrainDate, iterations, Cycles,SGD){

  SeparatorsAndAlphas <- PrepBoostedPercep(DATA, IndVars, DepVar, Thresh, LastTrainDate, iterations, Cycles, SGD)
  
  ###Test these separators on new data.          `
  TestDATA   <- DATA[DATA[,"Date"] > LastTrainDate,]
  BoostPreds <- array(NA,nrow(TestDATA))
  for(i in 1:nrow(TestDATA)){
     #BoostPreds[i] <- BoostedPred(SeparatorsAndAlphas,cbind(1,TestDATA[i,IndVars]))
     BoostPreds[i] <- TestConditionsBP(SeparatorsAndAlphas, as.numeric(TestDATA[i,IndVars]))
  }
  
  BoostPreds
}

PrepBoostedPercep <- function(DATA, IndVars, DepVar, Thresh, LastTrainDate, iterations, Cycles, SGD){
  TrainDATA <- DATA[DATA[,"Date"]<=LastTrainDate,]
  DepClasses = TrainDATA[,DepVar]>Thresh 
  
  ##Define the space for grid search
  Lows  <- c(0, -2.5, 0, -2.5)
  Highs <- c(1, 0, 30, 0)
  Cuts  <- 15 
  ###                                                                                  
  
  DayWeights <- array(1/length(DepClasses),length(DepClasses))
  Separators <- matrix(NA,iterations,(length(IndVars)+1))
  Alphas     <- array(NA,iterations)
  
  ##Create a decided number of weighted linear separators
  for(i in 1:iterations){
    print(paste("Iteration #", i))
    if(SGD==F){
      ClassesAndW <- GetPercepClasses(DepClasses,Lows,Highs,Cuts,TrainDATA,DayWeights)
    } else {
      ClassesAndW <- BoostedPerceptronSGD(TrainDATA, IndVars, DepVar, DepClasses, LastTrainDate, Cycles, DayWeights)
    }
    Separators[i,] <- ClassesAndW[(nrow(TrainDATA)+1):length(ClassesAndW)]
    Correct = ClassesAndW[1:(nrow(TrainDATA))]==DepClasses
    Error = sum(DayWeights*(!Correct))
    Alpha = 0.5*log((1-Error)/Error)
    Alphas[i] = Alpha
    DayWeights[Correct]  <- DayWeights[Correct]*exp(-Alpha)
    DayWeights[!Correct] <- DayWeights[!Correct]*exp(Alpha)
    NormSum              <- sum(DayWeights)
    DayWeights           <- DayWeights/NormSum
  }
  ###Obtain the optimal separators and their respective weights
  SeparatorsAndAlphas <- cbind(Separators,Alphas)
  colnames(SeparatorsAndAlphas)[1:ncol(Separators)] <- c("Intercept",IndVars)
  SeparatorsAndAlphas
}

BoostedPred <- function(SeparatorsAndAlphas,Vector){
   #Given a day's input from the test set, return the boosted prediction
   WeightedTotal <- 0
   Dim           <- ncol(SeparatorsAndAlphas)-1
   SumAlphas     <- sum(SeparatorsAndAlphas[,(Dim+1)])
   
   for(i in 1:nrow(SeparatorsAndAlphas)){
      Prediction <- sum(Vector*SeparatorsAndAlphas[i,1:Dim])
      if(Prediction > 0){
        WeightedTotal <- WeightedTotal + SeparatorsAndAlphas[i,(Dim+1)] 
      }
   }
   WeightedTotal <- WeightedTotal/SumAlphas
   WeightedTotal
}
###Weighted Prediction, uses all separators and Alphas

TestConditionsBP <- function(SeparatorsAndAlphas,Vals){
  #Run the boosted perceptron on an individual set of conditions
  BoostPred <- BoostedPred(SeparatorsAndAlphas,c(1,Vals))
  BoostPred
}

BoostedPerceptronSGD <- function(TrainDATA, IndVars, DepVar, DepClasses, LastTrainDate, Cycles, DayWeights){
  #Implements the boosted perceptron algorithm, using stochastic gradient descent
  #to fit the relevant weights for each hypothesis rather than gridsearch
  WeightStore <- matrix(NA,Cycles,length(IndVars)+1)
  ErrorStore  <- array(NA,Cycles)
  Weights     <- rep(0,length(IndVars)+1)
  for(n in 1:Cycles){
      Predictions <- array(NA,nrow(TrainDATA))
      for(i in 1:nrow(TrainDATA)){
         Predictions[i] <- max(sign(sum(Weights*as.numeric(c(1,TrainDATA[i,IndVars])))),0)
      }
         Errors = DayWeights*(DepClasses - Predictions)
      for(j in 1:length(Weights)){
         if(j == 1){        
           Weights[j] <- Weights[j] + 1/n*(sum(Errors))
         }else{
           Weights[j] <- Weights[j] + 1/n*(sum(TrainDATA[,IndVars[j-1]]*Errors))
         }
      }
      WeightStore[n,]    <- Weights
      ErrorStore[n]      <- sum(abs(Errors))
    }
                                
  BestHypothesis     <- which(ErrorStore==min(ErrorStore))[1]
  if(BestHypothesis==1){  ###Handle Non-Improvement Case
      if(length(which(ErrorStore==min(ErrorStore)))>1){
         BestHypothesis <- which(ErrorStore==min(ErrorStore))[2]
         Weights            <- WeightStore[BestHypothesis-1,]
      }else{
         Weights <- c(1,1,1,1)
      }
  }else{
      Weights            <- WeightStore[BestHypothesis-1,]
  }  
  for(i in 1:nrow(TrainDATA)){
      Predictions[i] <- max(sign(sum(Weights*as.numeric(c(1,TrainDATA[i,IndVars])))),0)
  }

  ClassesAndW <- c(Predictions,Weights)
  ClassesAndW
}

GetPercepClasses <- function(DepClasses, Lows, Highs, Cuts, DATA, DayWeights){
   Range <- matrix(NA,length(Lows),Cuts)
   #Define the grid search space
   for(i in 1:length(Lows)){
      Step <- (Highs[i]-Lows[i])/(Cuts-1)
      Range[i,] <- seq(Lows[i],Highs[i],Step)
   }
   LowestError <- 99999
   FinalPredClasses <- array(NA,length(DayWeights))
   #Test each combination, see which minimizes errors
   for(a in 1:ncol(Range)){
      print(a)
      for(b in 1:ncol(Range)){
         for(cc in 1:ncol(Range)){
            for(d in 1:ncol(Range)){
                Intercept <- Range[1,a]
                Precip    <- Range[2,b]
                Precip_1  <- Range[3,cc]
                PotEvap   <- Range[4,d]
                Weights = c(Intercept,Precip,Precip_1,PotEvap)
                Examples <- cbind(1,DATA[,IndVars])
                Predictions <- apply(Examples,1,PerformReg,Weights)
                PredClasses <- Predictions>0.5
             
                
                Errors <- (DepClasses-PredClasses)^2
                TotalError = sum(Errors*DayWeights)
                if(TotalError < LowestError){
                   print(paste("The Best Total Error is now:", TotalError))
                   LowestError <- TotalError
                   FinalPredClasses <- PredClasses
                   W <- Weights 
                }
            }
         }
      }
   }  
   c(FinalPredClasses,W)
}

###No, they don't line up, but this returns all necessary info
PerformReg <- function(Vec,W){
  ###Just multiplies inputs by the perceptron weight vector and sums
  Output <- sum(Vec*W)
  Output
}
