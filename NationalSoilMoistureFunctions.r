############################################################################
#############END MAIN#############START FUNCTIONS###########################
############################################################################

#Read in an "un-cleaned file", perform the cleaning and processing tasks,
#then write the cleaned file to disk for later use
PreProcessor <- function(SiteID,Path,KeyVars,EndPath,InvalidIfMissing){     
     SiteData <- GetCleanSCANFile(SiteID,Path)
     ##Process the file to standardize its column names...
     print("Standardizing Columns")
     SiteData <- ColumnNameProcess(SiteData)
     
     #Convert dates and times to appropriate numbers
     print("Converting Dates")
     SiteData[,"Date"] <- GetNumericalDates(SiteData[,"Date"])
     print("Converting Times")
     SiteData[,"Time"] <- GetNumericalTimes(SiteData[,"Time"])
     
     #Check for redundant time stamps (remove) and add filler times as needed...
     print("Cleaning Time Stamps")
     SiteData <- TimeStampCleaner(SiteData)
     
     #Pare this down the relevant variables...the next step is slow enough!
     print("Selecting Necessary Variables")
     SiteData <- SiteData[,KeyVars]
     
     ##Handle missing data
     print("Handling Missing Data")
     SiteData <- CleanMissingData(SiteData,"LinearConnect",InvalidIfMissing)
     
     ##Handle the "Frozen Sensor" case
     print("Frozen Precip Sensor Check")
     SiteData <- FrozenSensorFilter(SiteData,"PRCP.H.1..in.")
     
     ##The "LastValid" process is rather slow...thus, let's get these files
     ##onto the hard disk pre-formatted.  Thus, they can be opened quickly for 
     ##future use:
     WriteSCANFile(SiteData,SiteID,EndPath)
     SiteData
}

#Opens a cleaned SCAN file (this is a better function because it can take any
#path as an input, where the next function "GetSCANFile" cannot
GetCleanSCANFile <- function(SiteID,Path){
   filename <- paste("Site",SiteID,".csv",sep="")
   FullPath <- paste(Path,filename,sep="")
   SiteData <- read.csv(FullPath)
   SiteData
}

#Opens a given file from one of the SCAN sites (crappy...hard-coded path)
GetSCANFile <- function(SiteID){
   filename = paste("Site",SiteID,".csv",sep="")
   File <-read.csv(paste("C:\\Users\\Evan\\Desktop\\U. ILLINOIS WORK\\PhDResearch\\JohnDeere\\SCAN_DataFiles\\",filename,sep=""))
   File
}

#Writes a cleaned SCAN file to disk
WriteSCANFile <- function(SiteData,SiteID,Path){
   filename <- paste("Site",SiteID,".csv",sep="")
   FullPath <- paste(Path,filename,sep="")
   write.csv(SiteData,FullPath)
}

#Unfortunately, the column names for these site files from the SCAN network are
#not perfectly standarized...for instance, the soil moisture column might
#contain information about the soil type.  This is useful, but in terms of
#writing an algorithm that requests a specific variable, this is a pain.  To
#this end, this function standardizes the name given to each variable
ColumnNameProcess <- function(SiteData){
   ColNames <- colnames(SiteData)
   for(i in 1:length(ColNames)){
      Name <- ColNames[i]    #Grab the name for the given variable
      SplitName <- unlist(strsplit(Name, "\\."))  #Split on "."
      if(length(SplitName)>5){ #If it's non-standarized
        NewName <- paste(SplitName[1],SplitName[2],SplitName[3],SplitName[4],SplitName[5],sep=".")
        ColNames[i] <- NewName #Make the new name only what was needed
      }
   }
   colnames(SiteData) <- ColNames
   SiteData
}

#This site uses a value of -99.9 to denote a missing value, there are a couple
#of relevant options here.  The first is simple, take all the missing values,
#and set them equal to zero.  This could, in certain cases, make sense.  I.e.
#we could assume no precipitation or no solar radiation, and if this only
#occured for an occasional hour or two, this might be minor.  However, for
#certain variables, this could be problematic (i.e. temp of 0 during summer)
CleanMissingData <- function(SiteData,Procedure="ZeroOut",InvalidIfMissing){
     #Grab the indices where invalidating missing data is present BEFORE 
     #fixing them and filling in gaps.
   N <- nrow(SiteData)
   Invalid <- array(F,N)
   for(i in 1:N){  
     if(min(SiteData[i,InvalidIfMissing]) < 0){ 
        Invalid[i] <- T
     }                       
   }
     #This takes < 1 second
   if(Procedure=="ZeroOut"){    #If we're just zeroing things out
      SiteData[SiteData<-90]<-0
   }
     #This takes a file ~10-15 minutes...
   if(Procedure=="LastValid"){  #Replace -99.9 with the last non -99.9 value
      for(i in 1:ncol(SiteData)){
      print(colnames(SiteData)[i])
         LastValid <- 0;
         for(j in 1:nrow(SiteData)){
            if(as.numeric(SiteData[j,i]) > -90 || is.na(SiteData[j,i])){
               LastValid <- SiteData[j,i]
            }else{
               SiteData[j,i] <- LastValid  #Replace with most recent valid #
            }
         }
      }
   }
   if(Procedure=="LinearConnect"){  #This won't be fast either...
      for(i in 1:ncol(SiteData)){
      print(colnames(SiteData)[i])
          for(j in 1:nrow(SiteData)){
            if(as.numeric(SiteData[j,i]) < -90 || is.na(SiteData[j,i])){
                GapInfo <- GetNextValid(SiteData[,i],j)
               if(j > 1){
                  NewVals <- LinearReplace(j,GapInfo[2],SiteData[j-i,i],GapInfo[1]) #LinConnect
               }else{  #If it's the first time stamp and j-1 is NULL...
                  NewVals <- LinearReplace(j,GapInfo[2],0,GapInfo[1])
               }
               SiteData[(j):(GapInfo[2]),i] <- NewVals
               j <- GapInfo[2] #Jump ahead in the loop to the next real value
            }
         }
      }
   }
   #Add in the column of indices that denote invalid data
   if(nrow(SiteData)>N){
      SiteData <- SiteData[1:N,]
   }
   SiteData <- cbind(SiteData, array(NA,N))
   SiteData[,ncol(SiteData)] <- Invalid
   colnames(SiteData)[ncol(SiteData)] <- "Invalid"
   SiteData
}

#This function takes, as arguments, the index at which we first have an invalid
#data point, the index at which we regain validity, the last valid value of the
#variable, and the first valid value thereafter and creates a linear sequence
#"connecting the dots," so to speak"
LinearReplace <- function(StartIndex, EndIndex, LastValid, NextValid){
  Diff <- NextValid - LastValid #Distance over which a linear change must occur
  NumSteps <- EndIndex - StartIndex + 1
  NewVals <- array(NA,NumSteps)
  Step <- Diff/(NumSteps+1) #How big each step must be
  for(i in 1:NumSteps){
     NewVals[i] <- LastValid + i*Step
  }
  NewVals
}

#Given a string of "-99.9" or another such string of missing values, find the
#first value AFTER this string that is valid and its index, in other words
#"how many hours forward must we search?"
GetNextValid <- function(Array, Index){
  NextValid <- 0
  while(Index <= length(Array) && (Array[Index] < -99 || is.na(Array[Index]))){
     Index <- Index + 1
  }
  if(Index <= length(Array)){
     NextValid <- Array[Index] #This becomes the next valid
  }else{
     NextValid <- 0 #If there is no "next valid" value, just use 0.
  }
  Output <- c(NextValid, Index)  #What's the next valid point, and where is it?
  Output
}

#Unfortunately, at some of these SCAN sites, the hourly precip sensor becomes
#"frozen" and records the same hourly rainfall for 50-100+ consecutive hours.
#This will ruin the results, as every hour in that string after the first is 
#incorrect and will result in the model observing "rainfall" that does not make
#the soil wet.  These must be converted to 0 values.  Sadly, this introduces 
#another error...as if the exact same amount of rainfall (in hundreths of inches)
#does legitimately occur in consecutive hours, the latter will be zeroed.
FrozenSensorFilter <- function(SiteData,FVar){
   L = nrow(SiteData)
   for(i in 2:L){
      if(SiteData[i-1,FVar]==SiteData[i,FVar] && SiteData[i,FVar] > 0){
         SiteData[i,FVar] <- 0
      }
   }
   SiteData   
}


#Take a SCAN data file, and create new columns containing data aggregated over
#the past N hours.  For instance, create a column that contains total rainfall
#over the past 24 hours.  A weighting algorithm is offered to weight more recent
#values more heavily in the aggregation process.  For instance, rain that fell
#in the past 2 hours may play a more substantive role than rain from last night.
#This function can also provide an instaneous value for N-hours previously or an
#average value over the past N-hours.
AddNewVars <- function(SiteData, OldVars, HrsBack, Flavors, Weights){
   MaxBack <- max(HrsBack)  #We must shrink our dataset such that every point
                            #does have access to the key values N hrs ago...
   L = nrow(SiteData)-MaxBack
   NewColumns <- matrix(NA, L, length(OldVars))
   for(i in 1:length(OldVars)){
      KeyVar <- OldVars[i]
      print(KeyVar)
      for(j in (1+MaxBack):nrow(SiteData)){
          if(Flavors[i]=="Inst"){  #Grab N hours ago...
             NewColumns[j-MaxBack,i] <- SiteData[j-HrsBack[i],KeyVar]
          }
          if(Flavors[i]=="Ave"){ #Average the last N hours.
             KeyChunk <- SiteData[(j-HrsBack[i]):(j-1),KeyVar]
             NewColumns[j-MaxBack,i] <- mean(KeyChunk)
          }
          if(Flavors[i]=="Agg"){ #Aggregate last N hours, polynomial weighting
             KeyChunk <- SiteData[(j-HrsBack[i]):(j-1),KeyVar]
             Seq <- seq(1:length(KeyChunk))
             PolySeq <- Seq^Weights[i]
             NormPolySeq <- PolySeq/mean(PolySeq)
             NewColumns[j-MaxBack,i] <- sum(KeyChunk * NormPolySeq)
          }
      }
   }
   #Next, create descriptive column names for the "new variables"
   colnames(NewColumns) <- rep("NULL",length(OldVars))
   for(i in 1:length(OldVars)){
      NewName <- paste(OldVars[i],HrsBack[i],Flavors[i],sep="_")
      colnames(NewColumns)[i] <- NewName
   }

   #We need to chop off the rows of SiteData for which we would be missing data
   #then affix these new variables to build a workable dataset
   LiveRows <- SiteData[(1+MaxBack):nrow(SiteData),]
   SiteData <- cbind(LiveRows,NewColumns)
   SiteData
}

#The "date" column of each SCAN file is currently in "MM/DD/YYYY" format.  This
#is ultimately going to prove inconvenient.  Convert to a 1-365 value for the
#day of year, placing the year up front.  In otherwords, "01/22/2011" -> 2011022
#"04/23/2011" -> 2011113 (113th day).
GetNumericalDates <- function(DateCol){
   Leaps <- c(1988,1992,1996,2000,2004,2008,2012,2016)  #Leap Years
   DaysInMonths <- c(31,28,31,30,31,30,31,31,30,31,30,31)
   NewDates <- array(NA,length(DateCol))
   for(i in 1:length(DateCol)){
      #Split on "/"
      SplitDate = unlist(strsplit(as.character(DateCol[i]), "\\/"))
      Year <- as.numeric(SplitDate[3])
      Month <- as.numeric(SplitDate[1])
      Day <- as.numeric(SplitDate[2])
      if(Year %in% Leaps && Month > 2){  #Leap year, after Feb.
         DOY <- sum(DaysInMonths[1:(Month-1)])+Day+1
      }else{ #Not a leap-year-influenced date
         if(Month==1){
            DOY <- Day  #For January dates, we needn't do the summation
         }else{
            DOY <- sum(DaysInMonths[1:(Month-1)])+Day
         }
      }
      NewDates[i] <- Year*1000 + DOY #This is the new "date" value
   }
   NewDates
}

#The "time" column of each SCAN file is currently in "H:MM" format.  This
#function will convert it to a decimal between 0 (12AM) and 0.9999 (11:59:59...)
GetNumericalTimes <- function(TimeCol){
   NewTimes <- array(NA,length(TimeCol))
   for(i in 1:length(TimeCol)){
      SplitTime <- unlist(strsplit(as.character(TimeCol[i]), "\\:"))
      Hour <- as.numeric(SplitTime[1])
      Minute <- as.numeric(SplitTime[2])
      NewTimes[i] <- Hour/24+Minute/60/24
   }
   NewTimes
}

#Some of these files contain some odd time-stamps (23:59 e.g.), redundant
#time stamps (23:59 & 0:00, which are separated by only a minute), or missing
#time stamps.  This function will make a pass over the data and address each
#issue as it arises.  This may need editing as the process develops.
#Note: The input to this function must have ALREADY been date & time adjusted
#to yield numerical values, but not -99.9 cleaned, because this creates filler
TimeStampCleaner <- function(SiteData){
  #One time array, from 2011000.000 to 2011364.999
  CurrentTimes <- SiteData[,"Date"]+SiteData[,"Time"]-1 #so we end at 365...
  #This is for debug purposes, it lets me know where the small and large gaps
  #ultimately fall...feel free to remove once this works!
  DIFF <- CurrentTimes[2:length(CurrentTimes)] - CurrentTimes[1:(length(CurrentTimes)-1)]
  Smalls <- which(DIFF < 0.04)
  Larges <- which(DIFF > 0.05)
  OneHour <- 1/24
  TotLargeGap <- 0  #This prevents our i-notation from becoming off-line when
                   #when extra rows are introduced and the row numbers are "off"
  SmallGap <- 0  #Avoids the problem of finding a "gap" immediately
                 #after removing a row.  For instance:
                 # t = 0.958333 --> 0.999 --> 0.000 is 1 gap, not two...we only
                 # need to remove the 0.999 and we'll be all set

  #Get out all the large gaps, leaving the small gaps for later
  for(i in 2:nrow(SiteData)){
    if(floor(CurrentTimes[i]/1000) != floor(CurrentTimes[i-1]/1000)){
         FullMissingYears <- floor(CurrentTimes[i]/1000)-floor(CurrentTimes[i-1]/1000)-1
         #If these two time stamps straddle the end of one year...
         LastDay <- SiteData[i-1+TotLargeGap,"Date"] %% 1000
         if(LastDay == 366){
             #Leap year check
             TimeDiff <- CurrentTimes[i] - CurrentTimes[i-1] - 634 - FullMissingYears*635
         }else{                           ##"FullMissingYears kicks in if missing a full year or more
             TimeDiff <- CurrentTimes[i] - CurrentTimes[i-1] - 635 - FullMissingYears*635 
         }
    }else{
        TimeDiff <- CurrentTimes[i] - CurrentTimes[i-1]
    }
    if((TimeDiff - OneHour)>0.00001){ #If the gap is more than one hour
          GapHours <- floor(TimeDiff*24-0.999) #Create (Gap-1) new (dummy) rows
          SD <- SiteData[i-1+TotLargeGap,"Date"]  #Where to start building
          ST <- SiteData[i-1+TotLargeGap,"Time"]  #our new rows
          ColNames <- colnames(SiteData)
          ID <- SiteData[i,"Site.Id"]
          #This is the "filler" we must create
          NewRows <- CreateDummyRows(SD, ST, GapHours, ColNames, ID)
          #This is the chunk before the "filler" is needed
          EarlyChunk <- SiteData[1:(i-1+TotLargeGap),]
          #This is the chunk after the "filler" is needed
          LateChunk  <- SiteData[(i+TotLargeGap):nrow(SiteData),]
          #And now bind the early chunk to the late chunk with the "filler"
          SiteData <- rbind(EarlyChunk, NewRows, LateChunk)
          TotLargeGap <- TotLargeGap + GapHours #We're inserting new rows,
          #this impacts our indexing.  "TotLargeGap" accounts for this
    }
  }

  #Now, let's take a run at the small gaps (less than 1 hr between data points)
  CurrentTimes <- SiteData[,"Date"]+SiteData[,"Time"]-1 #RESET THIS
  for(i in 2:nrow(SiteData)){
     if(SmallGap == 0){ #If it's the first "too small" gap
       TimeDiff <- CurrentTimes[i] - CurrentTimes[i-1]
    }else{ #If it's the second "too small" gap and we've fixed it already
       TimeDiff <- CurrentTimes[i+1] - CurrentTimes[i]
       SmallGap <- 0 #Reset
    }
    if(TimeDiff - OneHour + 0.001 < 0){ #Probably a redundant stamp, REMOVE
      SiteData <- SiteData[-(i-1),] #remove the entire offending row
      SmallGap <- 1 #make sure not to note the second gap
    }
  }
  SiteData
}
#This function creates N filler rows into SiteData.  It will give these rows
#appropriate dates and times, but fill all the missing data values with -99.9
#the "CleanMissingData" function will then use either the ZeroOut or LastValid
#technique to fill these values appropriately.
CreateDummyRows <- function(StartDate,StartTime,N,ColNames,ID){
  if(N > 1){  #If we need to fill many rows, use a matrix
     DummyRows <- matrix(-99.9,N,length(ColNames))
     colnames(DummyRows) <- ColNames
  }else{     #Otherwise, use an array (R gets finnicky about arras vs. matrices)
     DummyRows <- array(-99.9,length(ColNames))
     names(DummyRows) <- ColNames
  }
  SD <- StartDate #We will update these as we move through to determine what
  ST <- StartTime #values must be introduced into the dummy rows
  if(N > 1){  #Matrix operations
    for(i in 1:N){
       if(floor(ST + 1/24 + .000001) > floor(ST)){ #If it's the last hour of the day
          if(floor(SD %% 1000) >= 364){
            #print(c(i,"IN HERE!"))
            ST <- 0
            SD <- (floor(SD/1000)+1)*1000
            DummyRows[i,"Date"] <- SD
            DummyRows[i,"Time"] <- ST
          }else{
              ST <- 0
              SD <- SD + 1
              DummyRows[i,"Date"] <- SD
              DummyRows[i,"Time"] <- ST
          }
       }else{
          ST <- ST + 1/24
          DummyRows[i,"Date"] <- SD
          DummyRows[i,"Time"] <- ST
       }
       DummyRows[i,"Site.Id"] <- ID
       DummyRows[i,4:length(ColNames)] <- -99.9 #Filler to be removed later
    }
  }else{
       if(floor(ST + 1/24 + .000001) > floor(ST)){ #If it's the last hour of the day
          if(floor(SD %% 1000) >= 364){ #If it's also the last hour of the year
            ST <- 0
            SD <- (floor(SD/1000)+1)*1000
            DummyRows["Date"] <- SD
            DummyRows["Time"] <- ST
          }else{
            ST <- 0
            SD <- SD + 1
            DummyRows["Date"] <- SD
            DummyRows["Time"] <- ST
          }
       }else{
          ST <- ST + 1/24
          DummyRows["Date"] <- SD
          DummyRows["Time"] <- ST
       }
       DummyRows["Site.Id"] <- ID
       DummyRows[4:length(ColNames)] <- -99.9 #Filler to be removed later
  }
  DummyRows
}

#This function receives the column names of the given cleaned data file, then
#checks the corresponding variables that must be aggregated/averaged/etc
#for each, it returns a T/F for its existence
DoVariablesExist <- function(ColNames,ProposedVars){
   N <- length(ProposedVars)
   CheckArray <- array(NA,N)
   for(i in 1:N){
      if(ProposedVars[i] %in% ColNames){ #If the variable is found in SiteData
         CheckArray[i] <- T
      }else{
         CheckArray[i] <- F
      }
   }
   CheckArray
}

#Read in EVERY cleaned SCAN file, calculate auto-correlation correlation
#coefficients at, say daily intervals, from 0 to Days.  Store these results for
#creation of visuals and later mathematical insights

AutoCorrelationCycle <- function(SensorList,Path,DepVar,Days){
    AutoCorMat <- matrix(NA,nrow(SensorList),Days+5)
    #Name the columns of this autocorrelation matrix
    ColNames <- rep(NA,Days+5)
    ColNames[1:5] <- colnames(SensorList)[1:5] #General site information
    for(i in 1:Days){
       ColNames[5+i] <- paste(i,"Days",sep="")
    }
    colnames(AutoCorMat) <- ColNames

   for(j in 1:nrow(SensorList)){
     ID <- SensorList[j,"SiteID"]
     print(c(j,ID))

     if(SensorList[j,"ReportingSince"] < 20110000){  #Not all files have 2011 data
        #Read in the cleaned file
        SiteData <- GetCleanSCANFile(ID,Path)
        N <- nrow(SiteData)
          if(N > 2400 && DepVar %in% colnames(SiteData)){
          #Record the site data
          AutoCorMat[j,1:5] <- as.numeric(SensorList[j,1:5])
          #Now perform the correlations
          for(i in 1:Days){
             Early <- SiteData[1:(N-i*24),DepVar] #Staggered back by i days
             Late <- SiteData[(i*24+1):N,DepVar]
             AutoCorMat[j,i+5] <- as.numeric(cor(Early,Late))
          }
        }
     }
   }
   AutoCorMat
}

#This function will divide a given cleaned SCAN file into a training and
#testing set.  Keep in mind, there are issues with anachronistic train/test
#splits when hourly data is involved.  For this reason the breakdown should be
#performed, at most, by month.  This function takes, as arguments, either a list
#of months to be used as training data (the rest can be used to test) or
#alternately, a split date where the training data is simply everything before
#The dataset returned will contain ONLY the independent and dependent variables
#along with a new column that is marked "TRAIN" or "TEST"
PreTestPrep <- function(SiteData, SplitType, IndVars, DepVar, TrainMonths=c(1,3,5,7,9,11), DivDate=2011200){
   #The smaller dataset we will ultimately return
   SubData <- SiteData[,c("Site.Id","Date","Time",IndVars,DepVar)]
   N <- nrow(SubData)

   #Whether a given record is entered into the training set or testing set
   #We may as well fill them up with "TEST" labels, and then assign "TRAIN"
   TrainTest <- array("TEST",N)

   #If we're going to split on a single date
   if(SplitType=="DateSplit"){
      #Set everything on or before that date to "TRAIN," else leave as "TEST"
      TrainIndices <- which(SubData[,"Date"]<=DivDate)
      TrainTest[TrainIndices] <- "TRAIN"
   }
   #If we're going to make certain months "training" and others "testing"
   if(SplitType=="MonthSplit"){
      #Take the YYYYDOY formatted numbers, and return, a month (0-12)
      DOY <- SubData[,"Date"] %% 1000
      Months <- sapply(DOY,DayOfYearToMonth)
      #Set only those months to "TRAIN," leaving the rest as "TEST"
      Trains <- which(Months %in% TrainMonths)
      TrainTest[Trains] <- "TRAIN"
   }
   #Now that every record is labeled, add this column back onto Subdata
   SubData <- cbind(SubData,TrainTest)
   SubData
}

#Take a day of the year (1 to 365) and returns the its month (1-12)
DayOfYearToMonth <- function(DOY){
   DaysInMonths <- c(31,28,31,30,31,30,31,31,30,31,30,31)
   CumulativeDays <- cumsum(DaysInMonths)  #No leap-year handling here...
   Month <- min(which(CumulativeDays>=DOY))
   Month
}

#The generic version of KNN. For historical record included, this function
#should return an estimate of the dependent variables.
#THIS IS MODIFIED FROM THE SOURCED VERSION in ML_Compendium.  For this version,
#we separate based on column with "TRAIN" and "TEST" as labels, rather than the
#date explicitly.  Subtle functional changes are made as well.
    KNN_SM <- function(IndVars,DepVar,k,HistoricalSet,TestSet,RegType,Unif=F,
                       Classify=F){
           #Construction around the test set
           NumRows                 <- nrow(TestSet)
           #The predictions made
           DepEstimates            <- array(NA,NumRows)
           if(Classify){
              K_Counts             <- array(NA,NumRows)
           }
           #What is the most common class, "guess" this without a convincing 
           #majority...if 40% of our k-votes say "LOW" and 30% say "HIGH", we
           #should aim for the middle
           ClassList <- unique(HistoricalSet[,DepVar])
           Default <- c(NA,0)
           for(i in 1:length(ClassList)){
              ClassCheck <- ClassList[i]
              TotalInClass <- sum(HistoricalSet[,DepVar]==ClassCheck)
              if(TotalInClass > Default[2]){
                 Default <- c(ClassCheck,TotalInClass)
              }
           }

           for(i in 1:NumRows){
            if(floor(i/100) == i/100){ print(i) }
            #The independent variables for the data point to be predicted
            X_i <- TestSet[i,IndVars]
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
               #How many of the similar set agreed on the answer?
               K_Counts[i]               <- sum(DistanceSet[,DepVar]==DepEstimates[i])
               if(as.numeric(K_Counts[i]) < k/2){
                   #Don't pick anything but the default without a majority
                   DepEstimates[i] <- Default[1]
               }
            }

            }
            if(Classify){
               DepEstimates <- cbind(DepEstimates,K_Counts)
            }
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



#Constructs a regression tree, for the purposes of soil moisture modeling, then
#tests the examaples from the testing set

TrainAndTestRegTrees_SM <- function(IndVars, DepVar, DATA, MinNodeSize, MinSplitSize, Classify=F){
      #Break data into the set used to build the tree and the test set
      TestSet  <- DATA[DATA[,"TrainTest"]=="TEST",c("Date",DepVar,IndVars)]
      TrainSet <- DATA[DATA[,"TrainTest"]=="TRAIN",c("Date",DepVar,IndVars)]

      ###Using the rpart library, construct and manipulate regression trees
      library(rpart)
      #minbucket determines the minimum size of a terminal node,
      #minsplit determines the minimum size for which to attempt a split
      #Once we make the tree, the object can be saved and stored
      RTree <- rpart(TrainSet[,DepVar] ~ .,TrainSet[,IndVars], minbucket = MinNodeSize, minsplit = MinSplitSize)
      plot(RTree,uniform=T,branch=1, margin=0.1, cex=0.9)
      text(RTree,cex=0.75)

      #Our forecasts
      RTreeForecasts <- array(NA,nrow(TestSet))
      #Send each set of inputs into the tree to generate a prediction
      for(i in 1:nrow(TestSet)){
         print(i)
         #print(paste("Data point #",i))
         NewData = list()
         for(j in 1:length(IndVars)){
           NewData[[IndVars[j]]] = TestSet[i,IndVars[j]]
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

#Generates a Beta value, as specified in Pan et al, 2011, equation 9:
#The first argument takes a vector of hourly precipitation data, which must be
#of length n or greater.  The second argument, the Eta series, is going to
#require the 3 parameter sinusoidal function matched up with the date/time
#stamps of the precipitation series.  Z is simply the sensor's soil depth.
#n represents the length we will look back.  This should be the minimum value
#necessary to stabilize Beta.  
GetBeta <- function(PrecipSeries, EtaSeries, z, n){
    WetDryRatioSeries  <- PrecipSeries/EtaSeries 
    OneMinusExpSeries  <- 1-exp(-1*EtaSeries/z)
    WetDry_w_ExpSeries <- WetDryRatioSeries * OneMinusExpSeries
    EtaOverZSeries     <- EtaSeries/z 
    ReverseEZSeries    <- rev(EtaOverZSeries)
    CumRevEZSeries     <- cumsum(ReverseEZSeries)#Sum backwards
    EZSeries           <- rev(CumRevEZSeries) #Reverse the summation...
    EtaRevMat          <- matrix(NA,length(EtaSeries),n-2)
    L                  <- length(EZSeries)
    for(j in 1:(n-2)){
       EtaRevMat[,j]   <- EZSeries - c(rep(0,j), EZSeries[1:(L-j)])
    }
    ExpEtaRevMat <- exp(EtaRevMat)
    #This is the big summation, eq. 9, Pan et al.
    #I know I can do this O(L)...can we do this with O(n)?
    BetaSeries <- c(rep(n,0),rep(NA,L-n))
    for(i in (n+1):L){
       BetaSeries[i] <- sum(rev(WetDry_w_ExpSeries[(i-n+2):(i-1)]) *
                        ExpEtaRevMat[i-1,]) + WetDry_w_ExpSeries[i]                     
    }  
    BetaSeries
} 


#This function will take SiteData, determine the time within the year of each 
#individual date, then calculate an Eta value for EVERY time stamp given three
#sinusoidal parameters c1,c2,c3
AddEtaSeries <- function(c1, c2, c3, SiteData){   
   Leaps <- c(1988,1992,1996,2000,2004,2008,2012,2016)  #Leap Years
      #Counts from 0.000 to 364.9583
   TimeInYearVec <- SiteData[,"Date"] %% 1000 + SiteData[,"Time"] - 1
   EtaSeries <- sapply(TimeInYearVec,EtaSin,c1,c2,c3,Leap=0)
   EtaSeries
}

#Runs the sinusoid function...if it's a leap year, the denominator becomes 366
EtaSin <- function(Time,c1,c2,c3,Leap=0){
  Value <- c1 + c2*sin(2*pi*(Time + c3)/(365+Leap))
  Value
}

#Runs the diagnostic equation.
DiagSoilEq <- function(B, e, re, c4){
   Value <- re + (e-re)*(1-exp(-c4*B))
   Value
}

#Given a beta series, and three parameters, return a soil moisture series
#for the time period in question.  This will be used to fit the final three
#parameters for this algorithm.
GetSMSeries <- function(BetaSeries, Porosity, ResidualSM, c4){
   SMSeries <- sapply(BetaSeries,DiagSoilEq,Porosity,ResidualSM,c4)
   SMSeries
}       

#Generates a "long-memory" beta series, of length 2000, then returns the 
#difference between its values and the "short-memory" version.  For cases where
#the long-memory version cannot be calcaulated (the long window precedes the 
#beginning of the dataset), the difference is assumed to be 0.    
GetBetaDiffs <- function(PrecipSeries,EtaSeries,z,LongTime,BetaSeries){
  #In case we've been in an especially rainy/dry period...this allows the
  #system to have a memory it otherwise might lack  
  Beta2000   <- GetBeta(PrecipSeries, EtaSeries, z, LongTime)
  #Since this look back is longer, we'll assume Beta2000 = BetaSeries for 
  #those cases where it is missing
  Beta2000[is.na(Beta2000)] <- BetaSeries[is.na(Beta2000)]  
  BetaDiffs <- Beta2000-BetaSeries
}                                                                                         

#Create a column of soil classifications (integers 1,2,...n) based on user-
#specified percentiles and bind this column to the meta site data.
AddSoilClass <- function(TrueSeq,StartInd,EndTrain,GrowingInd,QCuts=c(0.25,.75),
                         SnippedInvalid,MetaSiteData){
   L <- length(QCuts)
   SoilClass    <- array(1,length(TrueSeq))
   Cutoffs <- array(NA,L)
   for(i in 1:L){
      Cutoffs[i] <- quantile(TrueSeq[StartInd:EndTrain][GrowingInd][!SnippedInvalid],QCuts[i])
   }            
   for(i in 2:(L+1)){
      SoilClass[TrueSeq >= Cutoffs[i-1]] <- i
   }
   MetaSiteData <- cbind(MetaSiteData,SoilClass)
   MetaSiteData         
}

  
  ###WRITE A FUNCTION TO DETERMINE AN APPROPRIATE VALUE FOR N#######
GetSystemMemory <- function(DefaultTriad,SiteData,FirstValid,TestingSetHours,
                            PrecipSeries,z,LongTime,Interval,Thresh){
   DefaultEtas <- AddEtaSeries(DefaultTriad[1],DefaultTriad[2],DefaultTriad[3],
                         SiteData[FirstValid:(nrow(SiteData)-TestingSetHours),])
   L        <- length(DefaultEtas)             
   BetaLong <- GetBeta(PrecipSeries,DefaultEtas,z,LongTime)        
   for(i in 1:floor(LongTime/Interval)){
      n <- i*Interval
      TestBeta <- GetBeta(PrecipSeries,DefaultEtas, z, n)
      Correl   <- cor(TestBeta[(LongTime+1):L],BetaLong[(LongTime+1):L])
      if(Correl > Thresh){ break }
   }
   n
}
#begTime <- Sys.time()
#EtaSeries <- sapply(TimeInYearVec,EtaSin,c1,c2,c3,Leap=0)
#runTime <- Sys.time()-begTime
# <- 0

#Read in a file, given a directory
GetFile <- function(Name,Path){
   filename <- paste(Name,".csv",sep="")
   FullPath <- paste(Path,filename,sep="")
   Data <- read.csv(FullPath)
   Data
}

#Takes in a column that divides "soil moisture" rows by "sensor 10" or
#"sensor 20" etc, then processes them into appropriate row names
GetDepths <- function(DATA,Default){
  #Format to be addressed as the proper sensor depth
  Depths <- as.character(DATA[,"Type"])
  ValidName <- Depths[1] #This should be something like "Sensor 10"
  for(i in 1:length(Depths)){
     if(Depths[i]==Default){
       Depths[i] <- ValidName
     }else{
       ValidName <- Depths[i]
     }
  }
  Depths
}

#Reads in the time column " 06-03-2012 13:22:00 " and returns:
#Year = 2012, "Day of Year = 154.51(e.g.)
ProcessTimes <- function(DATA){
   L <- nrow(DATA)
   DaysInMonths <- c(31,28,31,30,31,30,31,31,30,31,30,31)
   CumSumMonths <- cumsum(DaysInMonths)
   Times <- as.character(DATA[,"Time"])
   MyDatesAndTimes <- matrix(NA,L,3)
   for(i in 1:L){
        SplitData <- strsplit(Times[i]," ")
     if(length(SplitData[[1]])>1){     
        SplitData <- c(SplitData[[1]][1],SplitData[[1]][2])
        MDY   <- strsplit(SplitData[1],"-")
        Month <- as.numeric(MDY[[1]][1])
        Day   <- as.numeric(MDY[[1]][2])
        Year  <- as.numeric(MDY[[1]][3]) 
        if(Month > 1){
           DOY <- CumSumMonths[Month-1] + Day
        }else{
           DOY <- Day
        }
        HourMinSec <- strsplit(SplitData[2],":") 
        Hour <- as.numeric(HourMinSec[[1]][1])
        Min  <- as.numeric(HourMinSec[[1]][2])
        Sec  <- as.numeric(HourMinSec[[1]][3]) 
        TOD <- Hour/24 + Min/24/60 + Sec/24/60/60
        DecimalDate <- DOY + TOD
        MyDatesAndTimes[i,] <- c(Year,Month,DecimalDate)
      }
   }
   colnames(MyDatesAndTimes) <- c("Year","Month","DayOfYear")
   MyDatesAndTimes
}

#Given appropriate parameters, and a "SiteData" matrix for training, return 
#a training set to be used for a site to be determined later
GetTrainingOrTestingSet <- function(ParameterList,TrainI,SiteDataTrain,Which,z,TestingSetHours=0,newC4){
  #Basic training site info 
  FirstValid <- ParameterList[TrainI,"FirstValid"]
  #Stick to growing season...
  GrowStart <- 100#ParameterList[i,"Grow_Start"]
  GrowEnd <- 300#ParameterList[i,"Grow_End"]
  LongTime <- 2000#Default
  n <- ParameterList[TrainI,"n"]
  #The FULL LENGTH versions
  PrecipSeries <- SiteDataTrain[FirstValid:nrow(SiteDataTrain),"PRCP.H.1..in."]
  TimeStamps   <- SiteDataTrain[FirstValid:nrow(SiteDataTrain),"Date"] %% 1000 + 
                  SiteDataTrain[FirstValid:nrow(SiteDataTrain),"Time"] - 1
  TimeOfDays   <- SiteDataTrain[FirstValid:nrow(SiteDataTrain),"Time"]               
  TrueSeq      <- SiteDataTrain[FirstValid:nrow(SiteDataTrain),SMDepthVar]
  InvalidSeq   <- SiteDataTrain[FirstValid:nrow(SiteDataTrain),"Invalid"]

  #Generate for the training set
  EtaSeries    <- AddEtaSeries(ParameterList[TrainI,"C1"], ParameterList[TrainI,"C2"], 
                      ParameterList[TrainI,"C3"], SiteDataTrain[FirstValid:nrow(SiteDataTrain),])
  BetaSeries <- GetBeta(PrecipSeries, EtaSeries, z, n)    
  
  BetaDiffs <- GetBetaDiffs(PrecipSeries,EtaSeries,z,LongTime,BetaSeries)
  #Are we doing an adjustment based on elevation?
  if(newC4==0){ C4 = ParameterList[TrainI,"C4"] } else {C4 = newC4}         
  SMSeries <- GetSMSeries(BetaSeries, ParameterList[TrainI,"Porosity"], 
                               ParameterList[TrainI,"ResSM"], C4)
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
  #MetaSiteData <- AddSoilClass(TrueSeq,StartInd,EndTrain,GrowingInd,
                                    #QCuts,SnippedInvalid,MetaSiteData)
       
  #Break into training and testing datasets
  TrainIndices <- seq(from=StartInd,to=EndTrain)
  TestIndices  <- seq(from=(EndTrain+1),to=EndInd)
  TrainData    <- MetaSiteData[TrainIndices,]                
  if(Which=="Test"){
    TestData     <- MetaSiteData[TestIndices,]
    #Filter out any examples that are not within the appropriate growing days
    TestGInd     <-  which(TestData[,"TimeStamps"] >= GrowStart &
                        TestData[,"TimeStamps"] <= GrowEnd)   
    TestInvalid  <- InvalidSeq[TestIndices][TestGInd]
    TestingSet   <- TestData[TestGInd,][!TestInvalid,]
  }                 
  TrainingSet  <- TrainData[GrowingInd,][!SnippedInvalid,]
  
  if(Which=="Train"){
    Output <- TrainingSet
  }else if(Which=="Test"){
    Output <- TestingSet
  }
  Output
}  

#Determine the inverse-square-distance weights for each of the 25 rain sensors
GetWeights <- function(LongVals,LatVals,Long,Lat){
   LongDiffs <- LongVals - Long
   LatDiffs  <- LatVals - Lat
   SquareDists <- LongDiffs^2 + LatDiffs^2
   NormFac <- sum(SquareDists)
   Weights <- SquareDists / NormFac
   Weights
}

###Given a lat and long and an appropriate nexrad file path, estimate the
##precipitation levels at each hour using inverse-square distance weighting
##also returns the rain dates (YYYYMMDD) and times (expressed as 0 to 1
##fraction, 12:00PM = 0.5, e.g.)
GetPrecipSeries <- function(PrecipSite,TestSitePath,SiteLat,SiteLon,Default){
  #Define the search path for the NEXRAD file
  NexradFileName <- paste(PrecipSite,"NEXRAD",sep="_")
  NexradPrecipFile <- GetFile(NexradFileName,TestSitePath)

  #Process this file into a format more conducive to coding
  Cols  <-      ncol(NexradPrecipFile)
  N     <-      nrow(NexradPrecipFile)
  Longs <-      as.numeric(NexradPrecipFile[2,2:Cols])
  Lats  <-      as.numeric(NexradPrecipFile[1,2:Cols])
  PrecipData <- NexradPrecipFile[3:N,2:Cols]  #Grab 25 precip sensors
  PrecipData[PrecipData < 0] <- Default  #Set missing data values to a default
  PrecipData <- PrecipData * 0.03937 #Convert from mm to inches
  RainDates <- NexradPrecipFile[3:N,1] #Grab dates
  RainTimes <- (RainDates-floor(RainDates))/0.24 #Convert to fractions of a day
  PrecipWeights <- GetWeights(Longs,Lats,SiteLon,SiteLat) #Inverse distance wts.
  WeightedPrecips <- sweep(PrecipData,MARGIN=2,PrecipWeights,'*')
  SitePrecip <- apply(WeightedPrecips,1,sum) #Weight the 25 sensors...and sum
  Output <- cbind(SitePrecip,RainDates,RainTimes)
  Output
}

#Get lats and longs of SM test site
GetSMLatLon <- function(CoordinatesBySitePath,TestSitePath,SMS){
  CoordinateList <- GetFile(CoordinatesBySitePath,TestSitePath)
  SoilMoistureFileList <- as.character(CoordinateList[,"FileName"]) #list of sites
  FileIndex <- which(SoilMoistureFileList == SMS) #which row?
  #Lat and Long where testing is needed
  SiteLat   <- as.numeric(CoordinateList[FileIndex,"Lat"])
  SiteLon   <- as.numeric(CoordinateList[FileIndex,"Lon"])
  LatLong <- c(SiteLat,SiteLon)
  LatLong
}

#Define path to soil moisture file
GetSMFile <- function(TestSitePath,SMFolder,SoilMoistureSite){
  SMPath <- paste(TestSitePath,SMFolder,"\\",sep="")
  SMFile <- paste("Out",SoilMoistureSite,sep="")
  SMData <- GetFile(SMFile,SMPath)
}

#Convert the soil moisture file dates to date/time stamps of YYYYMMDD,.58333 e.g.
GetDatesAndTimes <- function(SMData){
  DayOfYearSeries <- SMData[,"DayOfYear"]
  N <- nrow(SMData)
  DatesAndTimes   <- matrix(NA,N,2) #Fill this with YYYYMMDD Dates
  colnames(DatesAndTimes) <- c("Date","Time") #and 0 to 1 Times
  for(i in 1:N){
    DayOfYear <- as.character(SMData[i,"DayOfYear"])
    SplitDOY  <- strsplit(DayOfYear,"-")
    Month     <- SplitDOY[[1]][1]
    Day       <- SplitDOY[[1]][2]
    YearTime  <- SplitDOY[[1]][3]
    SplitYT   <- strsplit(YearTime, " ")
    Year      <- SplitYT[[1]][1]
    TimeOfDay <- SplitYT[[1]][2]
    Date      <- as.numeric(paste(Year,Month,Day,sep=""))
    SplitTOD  <- strsplit(TimeOfDay,":")
    Hour      <- SplitTOD[[1]][1]
    Min       <- SplitTOD[[1]][2]
    Sec       <- SplitTOD[[1]][3]
    Time      <- as.numeric(Hour)/24 +
                 as.numeric(Min)/24/60 +
                 as.numeric(Sec)/24/60/60
    DatesAndTimes[i,] <- c(Date,Time)
  }
  DatesAndTimes
}

#Take SM data, with its dates and times, and align this with existing 
#precipitation date to facilitate the subsequent analysis
MergeSMToPrecip <- function(Precips,PrecipDates,PrecipTimes,DatesAndTimes,SMSeries,
                            minDate, maxDate, minTime, maxTime){
  P <- length(Precips)
  SMMergedToPrecip <- array(9999,P)
  for(i in 1:P){
    PrecipDate <- as.numeric(PrecipDates[i])
    PrecipTime <- as.numeric(PrecipTimes[i])
    if(PrecipDate >= minDate){
      if((PrecipDate == minDate && PrecipTime >= minTime && PrecipDate <= maxDate 
          && PrecipTime <= maxTime) || (PrecipDate > minDate && PrecipDate < maxDate)){
      #which of the soil moisture dates are the same as the precip dates?
        DateInds    <- which(floor(DatesAndTimes[,"Date"])==PrecipDate)
        if(length(DateInds) > 0){ #handles missing soil moisture days
          PrecipTime  <- as.numeric(PrecipTimes[i])
          KeyChunk    <- DatesAndTimes[DateInds,"Time"]
          Diffs       <- abs(KeyChunk-PrecipTime)
          ClosestTime <- which(Diffs == min(Diffs))
          SMInd       <- DateInds[ClosestTime]
          SMVal       <- SMSeries[SMInd]
          SMMergedToPrecip[i] <- SMVal
        }else{
          SMMergedToPrecip[i] <- 9999 
        }
      }
    }
  }
  SMMergedToPrecip
}

#Build the "SiteData" matrix used in soil moisture functions
BuildSiteData <- function(ID,SMMerged, PrecipDates, PrecipTimes, Precips, SMDepthVar){
  L <- length(SMMerged) 
  SiteData <- matrix(F,L,6)
  colnames(SiteData) <- c("Site.Id", "Date","Time","PRCP.H.1..in.",SMDepthVar,"Invalid")
  SiteData[,"Site.Id"] <- rep(ID,L)
  SiteData[,"Date"] <- PrecipDates
  SiteData[,"Time"] <- PrecipTimes
  SiteData[,"PRCP.H.1..in."] <- Precips
  SiteData[,SMDepthVar] <- SMMerged
  Invalid <- which(SMMerged==9999)
  SiteData[Invalid,"Invalid"] <- T
  SiteData
}

#Converts a YYYYMMDD date into YYYYDOY format.  Ex 20120104 becomes 2012004
FinanceDateToDOYFormat <- function(DateList){
   Leaps <- c(1988,1992,1996,2000,2004,2008,2012,2016)  #Leap Years
   DaysInMonths <- c(31,28,31,30,31,30,31,31,30,31,30,31)
   L <- length(DateList)
   NewDates <- array(NA,L)
   for(i in 1:L){
      FDate <- DateList[i]
      Year  <- floor(FDate/10000)
      Month <- floor((FDate - Year*10000)/100)
      Day   <- FDate - Year*10000-Month*100
      if(Year %in% Leaps && Month > 2){  #Leap year, after Feb.
         DOY <- sum(DaysInMonths[1:(Month-1)])+Day+1
      }else{ #Not a leap-year-influenced date
         if(Month==1){
            DOY <- Day  #For January dates, we needn't do the summation
         }else{
            DOY <- sum(DaysInMonths[1:(Month-1)])+Day
         }
      }
      NewDates[i] <- Year*1000 + DOY #This is the new "date" value
   }
   NewDates
}

#Grab a "training set" from the calibration index
GetTrainingSetFromOther <- function(SiteIDTrain, EndPath, ParameterList,z,TSH=0,newC4){
  TrainI <- which(ParameterList[,"SiteID"]==SiteIDTrain)
  #Read in the site's data file
  SiteDataTrain <- GetCleanSCANFile(SiteIDTrain, EndPath)
  TrainingSet <- GetTrainingOrTestingSet(ParameterList,TrainI,SiteDataTrain,"Train",z,TSH,newC4)
  TrainingSet 
}



#Read in a comma-separated text file, either "ALL" or the number of rows desired
GetPartialText <- function(Name,Path,NumRows="ALL", ext){
   filename <- paste(Name,ext,sep="")
   FullPath <- paste(Path,filename,sep="")
   if(NumRows=="ALL"){
      Data <- tab5rows <- read.table(FullPath, header = FALSE, sep = ",")
   }else{
      Data <- tab5rows <- read.table(FullPath, header = FALSE, sep = ",", NumRows)
   }
   Data
}

#Given a lat/lon pair, determine if it is within a user-specified box defined by
#its NW(Lat,Lon) and SE(Lat,Lon) corners.  1 implies "in" 0 is "out"
is.in.box <- function(Lat,Lon,NWLat,NWLon,SELat,SELon){
  if(Lat <= NWLat && Lat >= SELat && Lon >= NWLon && Lon <= SELon){
    In <- 1
  }
  else{
    In <- 0
  }
  In
}

#Given a chunk of LiDAR data, remove every nth row if we deem that data to be
#not needed at the extremely high granularity at which it is provided
Sparsify <- function(Data, Good, SFac){
   L <- length(Good)
   GoodToKeep <- which(Good==1) #which ones are safe from random cuts
   Seq <- seq(1:L)
   Keepers <- Seq[floor(Seq/SFac)==Seq/SFac]  #which are kept
   Data <- Data[unique(c(Keepers, GoodToKeep)),] #keep every nth ("SFacth")
   Data
}

###Read in each of the partial files and append the needed rows
ParseLiDAR <- function(FileName,LiDARPath,SmallBox_NW,SmallBox,SE,SparseFac,NumFiles,ext){
  for(i in 0:(NumFiles-1)){
     FName <- paste(FileName,i,sep="_")
     print(FName)
     DataChunk <- GetPartialText(FName,LiDARPath,"ALL",ext) #Grab segment
     print("Performing Box Check")
     colnames(DataChunk) <- c("Lat","Lon","Elev")
     InBoxVec  <- mapply(is.in.box,DataChunk[,"Lat"],DataChunk[,"Lon"],SmallBox_NW[1],
                  SmallBox_NW[2],SmallBox_SE[1],SmallBox_SE[2]) #In small box?
     print("Sparsification")
     DataChunk <- Sparsify(DataChunk, InBoxVec, SparseFac)
     if(i == 0){
       LiDARSubset <- DataChunk
     }else{
       LiDARSubset <- rbind(LiDARSubset, DataChunk)
     }
  }
  write.csv(LiDARSubset, paste(LiDARPath,"Output.csv",sep=""))
  LiDARSubset
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

#Given an individual sensor, merge with NEXRAD precipitation data, and using
#the chosen soil moisture sensor calibration from SCAN, predict the SMSeries
#and merge...output a matrix containing the relevant information
SensorPrecipAndSM <- function(ID, zTrain, zTest, Default, EndPath,
                              LongTime, PrecipSite, SiteLat, SiteLon){

      SMDepthVar <- paste("SMS.I.1..",zTrain,sep="")
      ParameterList <- read.csv(paste(EndPath,"ParamList.csv",sep=""))

      #Get the precipitation time series(inverse-distance-weighted), dates, & times
      PrecipsDatesAndTimes <- GetPrecipSeries(PrecipSite,TestSitePath,SiteLat,SiteLon,Default)
      Precips     <- PrecipsDatesAndTimes[,"SitePrecip"]
      PrecipDates <- floor(PrecipsDatesAndTimes[,"RainDates"])
      PrecipDates  <- FinanceDateToDOYFormat(PrecipDates) #Convert to YYYYDOY format
      PrecipTimes <- PrecipsDatesAndTimes[,"RainTimes"]

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
      SMSeries     <- GetSMSeries(BetaSeries, ParameterList[i,"Porosity"],
                      ParameterList[i,"ResSM"], ParameterList[i,"C4"])
      ErrorSeq     <- TrueSeq - SMSeries
      #when does the soil moisture data begin and end?
      CleanSeq     <- which(TrueSeq!=9999)

      ProcessedData <- cbind(PrecipDates,PrecipTimes,PrecipSeries,SMSeries,
                             BetaSeries,EtaSeries,BetaDiffs,InvalidSeq,TrueSeq)
}