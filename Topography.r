#Generate a hypothetical random topography given a simple 2D psuedo-random-walk
#algorithm, then classify the various surfaces...
Path <- "C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\Topographies\\"

XD <- 100 #Number of points in x-direction
YD <- 100 #Number of points in y-direction

OutName <- "10K"
Layers <- c(1,5,10,50)

Elevation <- matrix(NA,XD,YD) #matrix to be filled by randomly-generated elevs
Elevation <- FillMatrix(XD,YD) #fill the matrix
LayerMat  <- GetLayerMat(XD,YD,Elevation,Layers) #classify every point, by layer
WriteLayerMat(LayerMat,Layers,Path,OutName)
write.csv(Elevation,paste(Path,OutName,"Elevation",".csv",sep=""))

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


