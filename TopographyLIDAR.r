rm(list=ls())

Path <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\Topographies\\"

GetLIDAR <- function(Path,FileName){
#Reads in a comma separated text file (x, y, elev)
  Direction = paste(Path,FileName,sep="")
  DATA = read.table(Direction,header=F,sep=",")
  DATA
}

#One of the largest challenges is in processing a dataset of this size in a
#manner that is not excruciatingly inefficient.  Functionality should be enabled
#to pick an x,y location, several user-defined distance measures, then create a
#dense grid in close proximity to that point, a semi-sparse grid at a second level,
#then a sparse grid farther away...

GetDistance <- function(KL,Lat,Long, LatScaleFac){
   #Given a two locations (one user-chosen), and the other presumably
   #from a LIDAR dataset, return the euclidian distance between them.
   #"LatScaleFac" addresses the different distances implied by 1 deg.
   #of latitude vs. 1 deg. of longitude
   LatDist = Lat - KL[1]
   LongDist = (Long - KL[2])*LatScaleFac
   Dist = sqrt(LatDist^2 + LongDist^2)
   Dist
}

DefineDistances <- function(KL,LIDARData,LatScaleFac){
#Given a location, and a LIDAR dataset (x,y,elev), return the distance
#from every point to the key location
  Distances = array(NA,L)
  for (i in 1:L){
    Distances[i] <- GetDistance(KL,LIDARData[i,"Lat"],LIDARData[i,"Long"]
                    ,LatScaleFac)
  }
  Distances
}

MakeSparse <- function(LData,Distances,SparseFac,DistanceCutoffs,L){
#Given a LIDAR dataset, a distance array, sparse factors and their distance
#boundaries, remove every nth point when needed to create a more managable array
  SparseCounters <- array(1,length(DistanceCutoffs)+1)
  DistanceCutoffs <- c(DistanceCutoffs,9999) #Everything is closer than 9999
  SubLData <- LData[1,]
  for(i in 1:L){
    DistClass <- min(which(DistanceCutoffs>Distances[i]))
    if(SparseCounters[DistClass]==SparseFac[DistClass]){
       print(paste("DistClass",DistClass,"row", i))
       SubLData <- rbind(SubLData,LData[i,]) #Append
    }
    SparseCounters[DistClass] <- SparseCounters[DistClass] + 1 #Cycle
    if(SparseCounters[DistClass] > SparseFac[DistClass]){
      SparseCounters[DistClass] <- 1
    }
  }
  SubLData
}

#Now, let's create a heuristic to take a point, and classify it's topography
GetTopoClass <- function(Lat,Long,Ranges,SData,BandThick,k,LatLongScale){
   #Lat/Long is the location of the point we are examining
   #Range is the distance, in deg lat/lon, used to scale our characterization
   #SData is our sparse dataset, BandThick is used for N,S,W,E calculations
   
   #Compare to current elevation
   #MiniData is a VERY small dataset, only immediately surrounding the key point
   while(k > 0){
     MiniData <- SData[SData[,"Long"]<=Long+BandThick &
                 SData[,"Long"]>=Long-BandThick &
                 SData[,"Lat"]>=Lat-BandThick*LatLongScale &
                 SData[,"Lat"]<=Lat+BandThick*LatLongScale,]
     k = min(k,nrow(MiniData))
     if(k == 0){ BandThick <- BandThick*2
       k=1
     } else {
       k = 0 #If we can't find a close point, search wider
     }
   }

   KeyElev <- GetElevationEst(Lat,Long,MiniData,k,LatLongScale)
   
   
   ElevClasses <- array(NA,length(Ranges))
   for(i in 1:length(Ranges)){
     #Make sure we don't exit the boundaries of the area we are examining
     N <- SData[SData[,"Long"]<=Long+BandThick & SData[,"Long"]>=Long-BandThick
        & SData[,"Lat"]>=Lat & SData[,"Lat"]<=Lat+Ranges[i]*LatLongScale,]
     S <- SData[SData[,"Long"]<=Long+BandThick & SData[,"Long"]>=Long-BandThick
        & SData[,"Lat"]<=Lat & SData[,"Lat"]>=Lat-Ranges[i]*LatLongScale,]
     W <- SData[SData[,"Lat"]<=Lat+BandThick*LatLongScale
        & SData[,"Lat"]>=Lat-BandThick*LatLongScale
        & SData[,"Long"]<=Long & SData[,"Long"]>= Long-Ranges[i],]
     E <- SData[SData[,"Lat"]<=Lat+BandThick*LatLongScale
        & SData[,"Lat"]>=Lat-BandThick*LatLongScale
        & SData[,"Long"]>=Long & SData[,"Long"]<= Long+Ranges[i],]
     #If we're so sparse that we're lacking...
     if(nrow(N)==0){
       N <- c(Lat,Long,KeyElev)
       names(N) <- c("Lat","Long","Elev")
       N <- rbind(N,N)
     }
     if(nrow(S)==0){
       S <- c(Lat,Long,KeyElev)
       names(S) <- c("Lat","Long","Elev")
       S <- rbind(N,N)
     }
     if(nrow(W)==0){
       W <- c(Lat,Long,KeyElev)
       names(W) <- c("Lat","Long","Elev")
       W <- rbind(N,N)
     }
     if(nrow(E)==0){
       E <- c(Lat,Long,KeyElev)
       names(E) <- c("Lat","Long","Elev")
       E <- rbind(N,N)
     }

     #Count how many directions are higher
     HigherCount <- 0
     if(mean(N[,"Elev"])>KeyElev){ HigherCount <- HigherCount + 1 }
     if(mean(S[,"Elev"])>KeyElev){ HigherCount <- HigherCount + 1 }
     if(mean(E[,"Elev"])>KeyElev){ HigherCount <- HigherCount + 1 }
     if(mean(W[,"Elev"])>KeyElev){ HigherCount <- HigherCount + 1 }

     if(HigherCount <= 1){ Class <- 2 #Peak
     }else if(HigherCount <= 2){ Class <- 1 #Intermediate
     }else{ Class <- 0 }#Valley
     ElevClasses[i] <- Class
   }

   OUT <- c(ElevClasses,KeyElev)
   OUT #return both the class AND the elevation
}

GetElevationEst <- function(Lat,Long,SData,k,LatScaleFac){
  #Given a lat,lon coordinate, determine the elevation of that point by grabbing
  #the k-nearest-neighbors (k should be small to ensure a strong estimate)
  N = nrow(SData)
  Distance_to_site <- array(NA,N)
  for (i in 1:N){
     Distance_to_site[i] <- GetDistance(c(Lat,Long),SData[i,"Lat"],
                            SData[i,"Long"],LatScaleFac)
  }
  SData <- SData[order(Distance_to_site),]
  SimSet <- SData[1:k,]
  Elev <- mean(SimSet[,"Elev"])
  Elev
}

#Take each constituent class-layer from a matrix of layer-wise classifications,
#then create a number ... i.e classes (1,0,2) becomes 102.
GetOverallClass <- function(LayerMat,Layers){
  ClassMat <- matrix(0,nrow(LayerMat),ncol(LayerMat))
  L <- length(Layers)
  for(i in 1:L){
     ClassMat <- ClassMat + LayerMat[,,(L-i+1)]*10^(L-i)
  }
  ClassMat
}

GetElev_and_Class <- function(Ranges,SData,BandThick,LatScaleFac,MaxLat, MinLat,
                           MaxLong, MinLong, LatGrid, LongGrid, Path, Filename){
  #For all of the points in a grid of user-specified granularity, determine
  #an estimated elevation and lat/lon class
  LatCut = (MaxLat - MinLat)/(LatGrid+1)
  LongCut = (MaxLong - MinLong)/(LongGrid+1)
  ClassLayers <- array(NA,dim=c(LatGrid,LongGrid, length(Ranges)))
  Elevations <- matrix(NA,LatGrid,LongGrid)
  for(j in 1:LatGrid){
    print(j)
    for(k in 1:LongGrid){
      Lat = MaxLat - LatCut*j
      Long = MinLong + LongCut*k
      C_and_E <- GetTopoClass(Lat,Long,Ranges,SData,BandThick,k,LatLongScale)
      Elevations[j,k] <- C_and_E[length(C_and_E)]
      ClassLayers[j,k,] <- C_and_E[1:(length(C_and_E)-1)]
    }
  }
  #Write the elevation to file
  write.csv(Elevations, paste(Path,FileName,"Elev",LatGrid,LongGrid,".csv",sep="_"))
  for(i in 1:length(Ranges)){
    ClassLayer <- ClassLayers[,,i]
    write.csv(ClassLayer,paste(Path,strsplit(FileName,split="_")[[1]][1],
              "Range",Ranges[i],".csv",sep="_"))
  }
  ClassLayers
}

#Converts every base 10 value of an array to a user-defined base
BaseConversion <- function(ClassMat,NumClasses){
  for(i in 1:nrow(ClassMat)){
    for(j in 1:ncol(ClassMat)){
       ClassMat[i,j] <- GetBase(ClassMat[i,j],NumClasses)
    }
  }
  ClassMat
}

#Change a number base X to base 10
GetBase <- function(Number,Base){
  if(Number==0){ NumDigits  <- 1 } else {
    NumDigits <- floor(log10(Number)+1) #How many digits in this number?
  }
  NewNum    <- 0
  for(i in 1:NumDigits){
     Digit <- floor((Number %% 10^i)/10^(i-1))
     NewNum <- NewNum + Digit*Base^(i-1)
  }
  NewNum
}

#Change a number from base 10 to base X
FromBaseTen <- function(Number,Base,AsNum=F){
  if(Number == 0){
    NewNum <- array(0,1)
    NumDigits <- 1
  }else{
    NumDigits <- floor(log(Number,Base))+1 #How many digits in base X
    NewNum <- array(NA,NumDigits)
    for(i in 1:NumDigits){
      VDigit <- Base^(NumDigits-i)  #Value of the current digit
      NewNum[i] <- floor(Number/VDigit)
      Number = Number - NewNum[i]*VDigit
    }
  }
  #Now return the "number" R will think it is base 10,
  #but we can treat it otherwise
  if(AsNum==T){
    OutputNum = 0
    for(i in 1:NumDigits){
       OutputNum = OutputNum + 10^(NumDigits-i) * NewNum[i]
    }
    NewNum = OutputNum
  }
  NewNum
}


#Read data from file
FileName = "MahometElevations_v2.txt"
LIDARData <- GetLIDAR(Path,FileName)
colnames(LIDARData) <- c("Lat","Long","Elev")

#Distances are measured in degrees lat/long
DistanceCutoffs <- c(.001,.01)
#Between 0(grab no points beyond this distance) to N (grab every Nth point)
SparseFac <- c(7,41,211)
            #c(1,11,53)
Ranges    <- c(0.001, 0.005, 0.01) #Used to create elevation classes
L = nrow(LIDARData)
KeyLocation <- c(40.205,88.45)
BandThick <- 0.0001
Distances <- DefineDistances(KeyLocation, LIDARData, LatScaleFac)
#Now make the dataset "sparse", so as to be more managably sorted and searched
SparseData <- MakeSparse(LData,Distances,SparseFac,DistanceCutoffs,L)

LatGrid = 100
LongGrid = 100
MaxLat = max(SparseData[,"Lat"])
MinLat = min(SparseData[,"Lat"])
MaxLong = max(SparseData[,"Long"])
MinLong = min(SparseData[,"Long"])
k=4

ClassLayers <- GetElev_and_Class(Ranges,SData,BandThick,LatScaleFac,MaxLat,
                     MinLat,MaxLong, MinLong, LatGrid, LongGrid, Path, Filename)
                     
ClassMat <- GetOverallClass(ClassLayers,Ranges) #Convert these layers to classes
ClassMat <- BaseConversion(ClassMat,length(Ranges)) #Converse to base X where
write.csv(ClassMat,paste(Path,FileName,"AllClasses.csv",sep="_"))
Elevations <- read.csv(paste(Path,FileName,"Elev",LatGrid,LongGrid,".csv",sep="_"))


SimpleFlatFile <- function(MaxLat,MinLat,MaxLong,MinLong,Elevations,
                           ClassMat,LatGrid,LongGrid,Ranges){
  #Takes various n x n array, and creates a much longer flat file...
   LatCut = (MaxLat - MinLat)/(LatGrid+1)
   LongCut = (MaxLong - MinLong)/(LongGrid+1)
   OutputMat <- matrix(NA,LatGrid*LongGrid,(length(Ranges)+4))
   colnames(OutputMat) <- c("Lat","Long","Elevations",paste("Ranges",Ranges),"AllClasses")
   for(i in 1:LatGrid){
      print(i)
      for(j in 1:LongGrid){
         Row = j + (i-1)*LongGrid
         OutputMat[Row,"Lat"] <- MaxLat - LatCut*i
         OutputMat[Row,"Long"] <- MinLong + LongCut*j
         OutputMat[Row,"Elevations"] <- Elevations[i,j+1] #The first column is 1,2,3...
         DiscreteClasses <- FromBaseTen(ClassMat[i,j], length(Ranges), F)
         OutputMat[Row,paste("Ranges",Ranges)] <-
         c(rep(0,length(Ranges)-length(DiscreteClasses)),DiscreteClasses)
         OutputMat[Row,"AllClasses"] <- ClassMat[i,j]
      }
   }
   write.csv(OutputMat,paste(Path,FileName,"LongFlat.csv",sep="_"))
}


