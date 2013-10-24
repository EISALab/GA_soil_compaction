#Generate a hypothetical random topography given a simple 2D psuedo-random-walk
#algorithm, then classify the various surfaces...
Path <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\Topographies\\"

XD <- 100 #Number of points in x-direction
YD <- 100 #Number of points in y-direction
k1 <- 0.2 #Slope parameters for the smooth plane put in place below the random- 
k2 <- 0.5 #walk generated surface
NumClasses <- 3  #3 Classes (Valley, Intermediate, and Peak)

OutName     <- "10K2"
Layers      <- c(1,10,50)
SmoothFac   <- 1 #How much shall we smooth the surface (0 to 1) 

Elevation <- matrix(NA,YD,XD) #matrix to be filled by randomly-generated elevs
Elevation <- FillMatrix(XD,YD) #fill the matrix
Elevation <- SmoothMatrix(XD,YD,Elevation,SmoothFac) #Smooth the matrix
Elevation <- AddPlane(XD,YD,k1,k2) #Add a smooth plane beneath the surface

LayerMat  <- GetLayerMat(XD,YD,Elevation,Layers) #classify every point, by layer
WriteLayerMat(LayerMat,Layers,Path,OutName) #Write each layer to file

ClassMat <- GetOverallClass(LayerMat,Layers) #Convert these layers to classes
ClassMat <- BaseConversion(ClassMat,NumClasses) #Converse to base X where
#X is the number of possible classifications 

write.csv(Elevation,paste(Path,OutName,"Elevation",".csv",sep=""))
write.csv(ClassMat,paste(Path,OutName,"Classes",".csv",sep=""))


##########################################################################
################################END MAIN##################################
##########################START FUNCTION LIST#############################
##########################################################################

#Build a classification matrix with a prescribed range
GetClassMat <- function(Elevation,Range){
  ClassMat <- matrix(NA,nrow(Elevation),ncol(Elevation))
  for(i in 1:YD){
    for(j in 1:XD){
       ClassMat[i,j] <- GetClass(i,j,Elevation,Range)
    }
  }
  ClassMat
}

#Classify a point topographically based on a certain range around that point
GetClass <- function(i,j, Elevation, Range){
   #Make sure we don't exit the boundaries of the area we are examining
   N <- max(1,i-Range)
   S <- min(YD,i+Range)
   E <- min(XD,j+Range)
   W <- max(1,j-Range)
   #What series are we averaging to determing relative position
   NSeries <- Elevation[N:i,j]
   SSeries <- Elevation[i:S,j]
   ESeries <- Elevation[i,j:E]
   WSeries <- Elevation[i,W:j]
   #Compare to current elevation
   KeyElev <- Elevation[i,j]
   #Count how many directions are higher
   HigherCount <- 0
   if(mean(NSeries)>KeyElev){ HigherCount <- HigherCount + 1 }
   if(mean(SSeries)>KeyElev){ HigherCount <- HigherCount + 1 }
   if(mean(ESeries)>KeyElev){ HigherCount <- HigherCount + 1 }
   if(mean(WSeries)>KeyElev){ HigherCount <- HigherCount + 1 }

   if(HigherCount <= 1){ Class <- 2 #Peak
   }else if(HigherCount <= 2){ Class <- 1 #Intermediate
   }else{ Class <- 0 }#Valley
   Class
}

#Fill the matrix using a 2D random walk
FillMatrix <- function(XD,YD){
    Elevation[1,1] <- 0 #Begin random walk in NW corner, at elev = 0

    #Fill the matrix
    #Fill the western edge, with single dependency (North)
    for(i in 2:YD){
      j=1
      Elevation[i,j] <- Elevation[i-1,j] + rnorm(1)
    }
    #Fill the northern edge, with single dependency (West)
    for(j in 2:XD){
       i=1
       Elevation[i,j] <- Elevation[i,j-1] + rnorm(1)
    }
    #Everything else, with two dependencies, North & West
    for(i in 2:YD){
      for(j in 2:XD){
         Elevation[i,j] <- 0.5*(Elevation[i-1,j]+Elevation[i,j-1])+rnorm(1)
      }
    }
  Elevation
}

#Smooth the matrix (optional)...SFac is a smoothing factor between 0 and 1...
#a value of of 1 will render the space fully smooth, 0 will do nothing
SmoothMatrix <- function(XD,YD,Elevation,SFac){
   #Smooth the western edge
   for(i in 2:(YD-1)){
     j=1
     SmoothVal <- 0.5*(Elevation[i-1,j]+Elevation[i+1,j]) 
     Elevation[i,j] <- SFac*SmoothVal + (1-SFac)*Elevation[i,j]
   }
   #Smooth the northern edge
   for(j in 2:(XD-1)){
     i=1 
     SmoothVal <- 0.5*(Elevation[i,j-1]+Elevation[i,j+1]) 
     Elevation[i,j] <- SFac*SmoothVal + (1-SFac)*Elevation[i,j]
   }
   #Smooth the "middle"
   for(i in 2:(YD-1)){
     for(j in 2:(XD-1)){
       SmoothVal <- 0.25*(Elevation[i-1,j-1]+Elevation[i-1,j+1]+
                          Elevation[i+1,j-1]+Elevation[i+1,j+1])
       Elevation[i,j] <- SFac*SmoothVal + (1-SFac)*Elevation[i,j]                 
     }
   }
   #Smooth the eastern edge (SmoothFac doesn't matter)
   for(i in 2:(YD-1)){
      j=XD
      Elevation[i,j] <- (Elevation[i,1]-Elevation[i,j])/YD + Elevation[i,j] 
   }
   #Smooth the southern edge (SmoothFac doesn't matter)
   for(j in 2:(XD-1)){
      i=YD
      Elevation[i,j] <- (Elevation[i,j]-Elevation[i,1])/XD + Elevation[i,j]
   }
   #Smooth the southeastern corner (SmoothFac doesn't matter)
   Elevation[YD,XD] <- 0.5*(Elevation[YD,XD-1]+Elevation[YD-1,XD])
   Elevation
}

#Build a 3D matrix of topographical classifications based on the layers listed
GetLayerMat <- function(XD,YD,Elevation,Layers){
  LayerMat <- array(NA,dim=c(XD,YD,length(Layers)))
  for(i in 1:length(Layers)){
    LayerMat[,,i] <- GetClassMat(Elevation,Layers[i])
  }
  LayerMat
}

#Write each individual layer class matrix to file, named appropriately
WriteLayerMat <- function(LayerMat,Layers,Dir,Out){
   for(i in 1:length(Layers)){
      Range <- Layers[i]
      FileName <- paste(Dir,Out,"_",Range,".csv",sep="")
      write.csv(LayerMat[,,i],FileName)
   }
}


#This adds a smooth function to Elevation, a plane which makes the jagged random
#walk appear less pronounced.
AddPlane <- function(XD,YD,k1,k2){
  Plane <- matrix(NA,XD,YD)
  for(i in 1:YD){
    for(j in 1:XD){
       Plane[i,j] <- k1*i + k2*j
    }
  }
  Elevation <- Elevation + Plane
  Elevation
}

#Change a number (base 10) to base X
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

#Take each constituent class-layer from a matrix of layer-wise classifications, 
#then create a number ... i.e classes (1,0,2) becomes 102.
GetOverallClass <- function(LayerMat,Layers){
  ClassMat <- matrix(0,YD,XD)
  L <- length(Layers)
  for(i in 1:L){
     ClassMat <- ClassMat + LayerMat[,,(L-i+1)]*10^(L-i)  
  }
  ClassMat
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