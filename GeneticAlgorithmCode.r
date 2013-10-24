#This set of scripts will contains functions required for genetic algorithm
#implementation for parameter fitting.  Hopefully, these functions will be
#written with sufficient generalization to be used for any parameter-fitting
#task that will be required.

FullGA <- function(InitialPop, Generations, ParamRanges, GenType, RepType, RGen,
                   NumMates, TrueSeq, TimeStamps, MutProb, MutMag, DeathRate,
                   PrecipSeries, z, n, c2constr, BetaSeries=0, GrowStart, GrowEnd,
                   InvalidSeq, WhatProcess){
                   
  #Create a store for the best creature
  BEST <- array(NA, ncol(ParamRanges)+2)
  BEST[length(BEST)-1] <- -9999999 #Starts off without any fitness whatsoever
  BEST[length(BEST)] <- 0
  names(BEST) <- c(rep("Param",ncol(ParamRanges)),"Fitness","Generation")

  #Initialize the Population
  Population <- InitializePopulation(InitialPop,ParamRanges,GenType,c2constr)

  #Determine the fitness of each member of the initial population
  #This arguments may need to be varied from use to use...in this case we need
  #times, but others may need precip, or other parameters
  
  if(WhatProcess == "FittingEta"){
     Population_w_Fitness <- GetFitnessE(Population, TrueSeq, TimeStamps, z, 
                                         GrowStart, GrowEnd, n, InvalidSeq)
  }
  if(WhatProcess == "FittingDiag"){
     Population_w_Fitness <- GetFitnessD(Population, TrueSeq, TimeStamps,
                                         GrowStart, GrowEnd, BetaSeries, n,
                                         InvalidSeq)
  }

  for(i in 1:Generations){
    print(paste("Generation", i))
    #Get fitnesses, post-mutation
    if(i > 1){
      if(WhatProcess == "FittingEta"){
         Population_w_Fitness <- GetFitnessE(Population_w_Fitness, TrueSeq, TimeStamps, z, 
                                             GrowStart, GrowEnd, n, InvalidSeq)
      }
      if(WhatProcess == "FittingDiag"){
         Population_w_Fitness <- GetFitnessD(Population_w_Fitness, TrueSeq, TimeStamps,
                                             GrowStart, GrowEnd, BetaSeries, n,
                                             InvalidSeq)
      }
    }
    #Store the "fittest" creature we've found thus far
    TopFit <- max(Population_w_Fitness[,"Fitnesses"],na.rm=T)
    if(length(BEST)>length(ParamRanges)+2){  #If it's 2+ rows...
      if(TopFit > BEST[nrow(BEST),"Fitness"]){  #If we found a better creature
        BestIndex <- which(Population_w_Fitness[,"Fitnesses"] == TopFit)
        NewBest <- c(Population_w_Fitness[BestIndex,],i)
        BEST <- rbind(BEST,NewBest)
        print(BEST) #SHOW THE CURRENT BEST SOLUTION
      }
    }else{ #If "BEST" is only 1 row
        BestIndex <- which(Population_w_Fitness[,"Fitnesses"] == TopFit)
        NewBest <- c(Population_w_Fitness[BestIndex,],i)
        BEST <- rbind(BEST,NewBest)
        print(BEST) #SHOW THE CURRENT BEST SOLUTION
    }

    #Perform the selection process, in this case, fitness proportional
    MatingPopulation <- Selection(Population_w_Fitness, NumMates)

    #Create new creatures, based on those organisms that have made it into the
    #selection set.  Two techniques are implemented.  The first involves the
    #standard, two-parent method, with a random crossover point.  The second
    #involves a slightly more progressive multi-parent technique, where a given
    #"chromosome" can be chosen from any of the available "parents"
    NewCreatures     <- Reproduction(MatingPopulation, RepType)
    
    #Generate a couple random creatures to avoid sticking at a local max
    if(RGen > 0){
      RandomNew        <- InitializePopulation(RGen,ParamRanges,GenType,c2constr)
      NewCreatures     <- rbind(NewCreatures, RandomNew)
    }
    
    #Get the fitnesses of the new creatures
      if(WhatProcess == "FittingEta"){
         NewCreatures <- GetFitnessE(NewCreatures, TrueSeq, TimeStamps, z, 
                                             GrowStart, GrowEnd, n, InvalidSeq)
      }
      if(WhatProcess == "FittingDiag"){
         NewCreatures <- GetFitnessD(NewCreatures, TrueSeq, TimeStamps,
                                             GrowStart, GrowEnd, BetaSeries, n, 
                                             InvalidSeq)
      }

    #Build the bigger, full population
    Population_w_Fitness <- rbind(Population_w_Fitness,NewCreatures)

    #Kill off weaker creatures
    Population_w_Fitness <- Death(Population_w_Fitness, DeathRate)

    #Allow for mutation over all creatures
    Population_w_Fitness <- Mutation(Population_w_Fitness, MutProb, MutMag, ParamRanges, c2constr)

  }
  BEST
}

#Generate initial creatures within a given parameter space, using either
#Monte Carlo generation or Latin Hypercube methodology (user-specified).
#ParamRanges is a (2 x dim) matrix with the max of the ranges on the 1st row
#and the min of the ranges on the 2nd row
InitializePopulation <- function(InitialPop, ParamRanges, GenType, c2constr){
    Dim <- ncol(ParamRanges)
    Ranges <- ParamRanges[1,] - ParamRanges[2,]
    Diffs <- Ranges/InitialPop
    Population <- matrix(NA, InitialPop, Dim)
    RandVals <- matrix(runif(InitialPop*Dim), InitialPop, Dim)
    if(GenType == "MonteCarlo"){
        #just multiply by range
       RandVals <- sweep(RandVals, 2, Ranges, '*')
        #Then add the min
       RandVals  <- sweep(RandVals, 2, ParamRanges[2,], "+")
       Population <- RandVals
    }
    if(GenType == "LatinHypercube"){  #We don't "fill the cube" but rather, we
    #ensure that no 'box' is filled with more than one creature.
       for(j in 1:Dim){
          Order <- order(RandVals[,j])-1 #Sort the random values
          NewRand <- runif(InitialPop) #Generate new random values
          Population[,j] <- as.numeric(ParamRanges[2,j]) + Order*Diffs[j] + NewRand*Diffs[j]
       }
    }
    if(c2constr == T){ #If c2 < c1
      for(i in 1:InitialPop){
         if(Population[i,1] <= Population[i,2]){
            Population[i,2] <- runif(1)*Population[i,1] #between 0 and c1
         }
      }
    }
    Population
}

#Obtain a fitness metric for every possible solution within the population
################################################################################
#  This is the function that is likely to vary with each new implementation of
#  the genetic algorithm.  For modeling a time series, we need to use each
#  of the parameters to generate an equivalent time series, then use least
#  square errors to determine a fitness value.  For soil moisture modeling, for
#  example, we'll need to use parameters to generate years of estimates, which
#  will require a precipitation time series as well, and thus, a new function.
#
#  Unfortunately, this means this function cannot be standarized.  It will need
#  to be versioned for each type of use.
################################################################################

#Version #1, fitting a curve
#GetFitness <- function(Population, TrueSeq, TimeStamps){
#   Fitnesses <- array(NA,nrow(Population))
#   for(i in 1:length(Fitnesses)){
#      #Generate the sequence using the parameters from one creature
#      SimSeq     <- sapply(TimeStamps,EtaSin,Population[i,1],Population[i,2],Population[i,3],Leap=0)
#      Errors     <- SimSeq - TrueSeq
#      Fitnesses[i] <- sum(Errors^2)*-1 #Simple least-square errors
#   }
#   PopAndFit <- cbind(Population, Fitnesses)
#   PopAndFit #Negative so that fewer errors are "better"
#}
#####################################
#Version #2, maximizing correlation between two series
GetFitnessE <- function(Population, TrueSeq, TimeStamps, z, GrowStart, GrowEnd, 
                        n, InvalidSeq){
   Fitnesses <- array(NA,nrow(Population))
   StartInd <- n+1
   EndInd   <- length(TrueSeq)
   #Cutoff everything before the start of growing season and after the end
   GrowingInd  <- which(TimeStamps[StartInd:EndInd] >= GrowStart &
                        TimeStamps[StartInd:EndInd] <= GrowEnd)
   #This can be done once...it's the same for all creatures
   SnippedInvalidSeq <- InvalidSeq[StartInd:EndInd][GrowingInd]
   SnippedTrueSeq    <- TrueSeq[StartInd:EndInd][GrowingInd]

   for(i in 1:length(Fitnesses)){
      #Generate the sequence using the parameters from one creature
      SimSeq     <- sapply(TimeStamps,EtaSin,Population[i,1],Population[i,2],Population[i,3],Leap=0)
      BetaSeries <- GetBeta(PrecipSeries, SimSeq, z, n)
      SnippedBetaSeries <- BetaSeries[StartInd:EndInd][GrowingInd]

      Fitnesses[i] <- cor(SnippedBetaSeries[!SnippedInvalidSeq], 
                          SnippedTrueSeq[!SnippedInvalidSeq])-1
      #correlation between SM & Beta
      #The "-1" at the end is simply so fitnesses are negative, with a score of
      #zero indicating a "perfect" score.
   }
   PopAndFit <- cbind(Population, Fitnesses)
   PopAndFit #Negative so that fewer errors are "better"
}
######################################
#Version #3, minimizing the square errors between an predicted and empirical
#soil moisture series.
GetFitnessD <- function(Population, TrueSeq, TimeStamps, GrowStart, GrowEnd, 
                        BetaSeries, n, InvalidSeq){
   Fitnesses <- array(NA,nrow(Population))
   StartInd <- n+1
   EndInd   <- length(TrueSeq)
   #Cutoff everything before the start of growing season and after the end
   GrowingInd  <- which(TimeStamps[StartInd:EndInd] >= GrowStart &
                       TimeStamps[StartInd:EndInd] <= GrowEnd)
   #This can be done once...it's the same for all creatures
   SnippedInvalidSeq <- InvalidSeq[StartInd:EndInd][GrowingInd]


   for(i in 1:length(Fitnesses)){
      #Generate the sequence using the parameters from one creature
      SMSeries <- GetSMSeries(BetaSeries, Population[i,1], Population[i,2], Population[i,3])
      Errors   <- SMSeries[StartInd:EndInd][GrowingInd][!SnippedInvalidSeq] - TrueSeq[StartInd:EndInd][GrowingInd][!SnippedInvalidSeq]
      SSE      <- sum(Errors^2)
      Fitnesses[i] <- SSE * -1
      #sum of square errors
      #The "-1" at the end is simply so fitnesses are negative, with a score of
      #zero indicating a "perfect" score.
   }
   PopAndFit <- cbind(Population, Fitnesses)
   PopAndFit #Negative so that fewer errors are "better"
}


#Given a population of possible solutions along with their fitness metrics,
#select the group that will reproduce.  This uses fitness proportionate selection
Selection <- function(PopAndFit,NumMates){
   Fitnesses <- PopAndFit[,"Fitnesses"]
   MinFit <- min(Fitnesses)
   #Because fitnesses are negative...this makes all of them positive AND non-zero
   Fitnesses <- Fitnesses - 1.01*MinFit
   PopAndFit[,"Fitnesses"] <- Fitnesses
   PopAndFit <- PopAndFit[order(PopAndFit[,"Fitnesses"]),]
   #We will fill this up with the members of the population that are chosen
   SelectedSet <- matrix(NA, NumMates, ncol(PopAndFit))

   #Iterate, randomly choosing a member of the population, then eliminating it
   #so we can choose from those which remain
   for(i in 1:NumMates){
      SumFit <- sum(PopAndFit[,"Fitnesses"])
      Ps     <- PopAndFit[,"Fitnesses"]/SumFit  #Normalized, sum of all Ps = 1
      RouletteWheel <- cumsum(Ps) #Thus we can pick with a random 0-1
      MatePick <- runif(1) #Which population gets picked?
      Index <- min(which(RouletteWheel > MatePick))
      SelectedSet[i,] <- PopAndFit[Index,]  #Grab creature
      PopAndFit <- PopAndFit[-Index,] #Remove that creature from the next selection...
   }
   SelectedSet
}

#Takes the "selected set," those creatures chosen by fitness to be worthy
#of reproduction, and produces 1/2 that number in new offspring (rounded down)
Reproduction <- function(SelectionSet, Type="Two-Parent"){
  NumCreatures <- nrow(SelectionSet)
  NumOffspring <- floor(NumCreatures/2)
  NewCreatures <- matrix(NA,NumOffspring,(ncol(SelectionSet)-1))
     #Now generate new creatures
     if(Type=="Two-Parent"){
       i = 1
       while(nrow(SelectionSet)>1){
         #Choose the two parents
         ParentPick <- order(runif(nrow(SelectionSet)))[1:2]
         #Choose the cross-over point
         CrossoverPoint <- floor(runif(1)*(ncol(SelectionSet)-2))+1
         #Create the child
         NewChild <- c(SelectionSet[ParentPick[1],(1:CrossoverPoint)],
                       SelectionSet[ParentPick[2],((CrossoverPoint+1):
                       (ncol(SelectionSet)-1))])
         #Store the child and raise the counter
         NewCreatures[i,] <- NewChild
         i = i + 1
         #Remove the parents from the SelectionSet
         SelectionSet <- SelectionSet[-ParentPick,]
       }

     }
     if(Type=="Multi-Parent"){
        for(i in 1:NumOffspring){
           for(j in 1:(ncol(SelectionSet)-1)){
               #For each chromosome, grab from any random parent
               WhichParent <- floor(runif(1)*nrow(SelectionSet)+1)
               NewCreatures[i,j] <- SelectionSet[WhichParent,j]
           }
        }
     }
     NewCreatures
}

#Take a population, and allow for random mutations for every chromosome of
#every creature.  MutProb is the probability that any given chromosome will be
#allowed to mutate.  All mutations are Gaussian, centered upon the current
#value of the chromosome.  Mutations will not be permitted to exit the allowable
#ranges.  MutMag is the standard deviation of the Gaussian, expressed as a
#fraction of that variable's allowable range.  For instance, for date (0-365), a
#MutMag of 0.3 implies a mutation of N ~ (0, [365-0]*.3), i.e. N ~ (0,109.5)
Mutation <- function(Population_w_Fitness, MutProb, MutMag, ParamRanges, c2constr){
  Chromosomes <- ncol(Population_w_Fitness) - 1
  Creatures   <- nrow(Population_w_Fitness)
  for(i in 1:Creatures){
     for(j in 1:Chromosomes){
        MutationCheck <- runif(1)
        if(MutationCheck < MutProb){ #does this chromosome mutate?
           Population_w_Fitness[i,j] <- min(max((Population_w_Fitness[i,j] +
           rnorm(1)*MutMag*(ParamRanges[2,j]-ParamRanges[1,j])),
           ParamRanges[2,j]),ParamRanges[1,j])
        }
     }
  }
  #Delete existing fitness metric - they need to be recalculated after mutation
  Population_w_Fitness <- Population_w_Fitness[,-ncol(Population_w_Fitness)]

  #Ensure that c1 is still > c2
  if(c2constr == T){
    for(i in 1:Creatures){
       if(Population_w_Fitness[i,1] <= Population_w_Fitness[i,2]){
          Population_w_Fitness[i,2] <- Population_w_Fitness[i,1] - 0.0001
       }
    }
  }

  Population_w_Fitness
}


#This function kills off a fixed proportion of the existing creatures, as the
#user defines by the parameter "DeathRate" (rounded DOWN).  Be warned, if the
#DeathRate value exceeds 0.5*NumMates/(0.5*NumMates + InitialPop), numbers will
#eventually dwindle to zero.  This function is essentially a means of preventing
#the population from growing so large as to create issues with memory
#requirements
Death <- function(Population_w_Fitness, DeathRate){
   #Error catch...remove rows with NA for fitness
   Population_w_Fitness <- Population_w_Fitness[!is.na(Population_w_Fitness[,"Fitnesses"]),]
   NumDeaths <- floor(DeathRate*nrow(Population_w_Fitness))
   if(NumDeaths > 0){
     for(i in 1:NumDeaths){
        TotFit <- sum(Population_w_Fitness[,"Fitnesses"],na.rm=T)
        P_Deaths <- Population_w_Fitness[,"Fitnesses"]/TotFit #Likelihood of dying
        RouletteWheel <- cumsum(P_Deaths)
        Scythe <- min(which(RouletteWheel >= runif(1)),na.rm=T) #Which creature dies
        Population_w_Fitness <- Population_w_Fitness[-Scythe,]
     }
   }
   Population_w_Fitness
}




#TimeStamps <- seq(from=0,to=730,by=(1/24))