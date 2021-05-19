### SummaryStats - Calculates overall statistics for Cubist outputs
#
# Version 3.2 - 19/05/2021 - Fix bug that occurs when model runs only have one soil attribute and/or no depths
# Version 3.1 - 16/05/2017 - Implemented a 'user adjustable project directory' & Fixed ModelData.csv bug
# Version 3.0 - 23/02/2017 - Added model run number to filename ModelRun.csv
#
### Script description 
####### Part 1 - combines the median stats from individual tables,
####### Part 2 - averages the variable (covariate) importances
####### Part 3 - creates a raster for each predicted attribute that represents the difference between the values for the 95th percentile and 5th percentile (range),
####### Part 3 - creates rasters for the average range for each attribute aswell the overall range
####### Script was initially designed to expect multiple soil attributes and six soil depths
####### It has since been updated to handle case where there is only one soil attribute run aswell as
####### cases where there is less than six soil depths
####### Make sure only model runs that have completed successfully are in the /TrainingDataOutputs/ directory
###
### Part 1 & 2 can be run after Module 1
### Part 3 can be run after Module 2
###

### USER INPUT 
## If your covariates are not GTiff .tif then lines 37 and 38 need to be remarked out and lines 39 and 40 need to be made active

## Adjust model run number and project directory below (Lines 25 & 26)
ModelRun = "BMSE1"

ProjectDir = "D://"
##ProjectDir = "//lands//data//DSITI//LandSciences//NAS//slr_soils//Projects//PMap//Modelling//Stage2//"
##ProjectDir = "//athenasmb/scratchDSITI/zundp/"

## Script processing starts here
library(plyr)
library(matrixStats) #for std deviation function rowSds()

# The following code is for a preliminary step in generating an excel summary file of 'Variable Importance'(Part 2)
# Creates a blank dataframe with a list of covariates used in the model run to which data will be added later.
CovariateDirectory = paste(ProjectDir, ModelRun, "//Covariates", sep="") 
setwd(CovariateDirectory)
Covariate <- list.files(getwd(), full.names=FALSE, pattern=paste0(gsub("[.]","",".tif"),"$")) #Obtain list of covariates
Covariate <- gsub('.{4}$', '', Covariate) #delete last five characters (file extension)
#Covariate <- list.files(getwd(), full.names=FALSE, pattern=paste0(gsub("[.]","",".sdat"),"$")) #Obtain list of covariates
#Covariate <- gsub('.{5}$', '', Covariate) #delete last five characters (file extension)

# Part 1 - Create a summary csv file showing all 'median value' stats for all attributes modelled
ModelRunDirectory = paste(ProjectDir, ModelRun, "//TrainingDataOutputs", sep="") 
setwd(ModelRunDirectory)

Attribute <- list.files(getwd(), full.names=FALSE) #Obtain list of attributes modelled
AllSummary = data.frame(Kfold_level = numeric(0), Calib_RMSE = numeric(0), Calib_R2 = numeric(0), Calib_Bias = numeric(0), Calib_CC = numeric(0), Valid_RMSE = numeric(0), Valid_R2 = numeric(0), Valid_Bias = numeric(0), Valid_CC = numeric(0), Perc.within.UPP.LOW.limits = numeric(0), Depth = character(0), Attribute = character(0)) # Create a blank df to bind to later
blankrow = data.frame(Kfold_level = numeric(0), Calib_RMSE = numeric(0), Calib_R2 = numeric(0), Calib_Bias = numeric(0), Calib_CC = numeric(0), Valid_RMSE = numeric(0), Valid_R2 = numeric(0), Valid_Bias = numeric(0), Valid_CC = numeric(0), Perc.within.UPP.LOW.limits = numeric(0), Depth = character(0), Attribute = character(0)) # Create a blank spacer row
blankrow <- ""

#Iterate over all soil attributes in ModelRun Directory
for (g in 1:length(Attribute)){
  
  # Part 1 combines the median stats from the individual diagnostic tables outputed by Webbs Module 1 Cubist script
  InputDirectoryDiag=paste((ModelRunDirectory),"//", (Attribute[g]), sep="") 
  setwd(InputDirectoryDiag)
  Depths <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE) #depths modelled
  if(length(Depths) == 6){ Depths <- Depths[c(1, 5, 3, 4, 6, 2)]} #reorder columns if there are six depths
  Summary = data.frame(Kfold_level = numeric(0), Calib_RMSE = numeric(0), Calib_R2 = numeric(0), Calib_Bias = numeric(0), Calib_CC = numeric(0), Valid_RMSE = numeric(0), Valid_R2 = numeric(0), Valid_Bias = numeric(0), Valid_CC = numeric(0), Perc.within.UPP.LOW.limits = numeric(0), Depth = character(0), Attribute = character(0))
  
  for (d in 1:length(Depths)){
    Filename = paste((Attribute[g]), (Depths[d]), "_Diagnostics.csv", sep = "")
    diagnostics <- read.table(Filename, sep=",", header=TRUE)[12, ]
    diagnostics$Depth <- Depths[d]
    Summary <- rbind(Summary, diagnostics )
  }
  Summary$Attribute <- (Attribute[g]) # Add an attribute column 
  AllSummary <- rbind(AllSummary, Summary) #Create a df with all attribute stats
  AllSummary <- rbind(AllSummary, blankrow) #Add a blank row between attributes
  #reorder columns and export csv
  Summary <- Summary[c(12, 11, 1:10)]
  Filename = paste((Attribute[g]), "overall_diagnostic_stats.csv", sep = "")
  write.table(Summary, file = Filename, sep = ",", row.names = FALSE, col.names = TRUE)
}

#Export all attribute stats to csv
ModelRunDirectoryUpper = paste(ProjectDir, ModelRun, "//", sep="")
setwd(ModelRunDirectoryUpper)
#reorder columns and export csv
AllSummary <- AllSummary[c(12, 11, 1:10)]
Filename = paste((ModelRun), "_attribute_diagnostic_stats.csv", sep = "")
write.table(AllSummary, file = Filename, sep = ",", row.names = FALSE, col.names = TRUE)

#Add model run code to modeldata.csv and remove old file (i.e. modeldata.csv to modeldata_FSE11.csv)
ModelData <- read.table("modeldata.csv", sep=",", header=TRUE)
MDFilename <- paste("modeldata_", ModelRun, ".csv", sep = "")
write.csv(ModelData,file = MDFilename)
file.remove(paste(ModelRunDirectoryUpper, "modeldata.csv", sep=""))

### End of Part 1

### Part 2 #### 
### Combine the covariate importance stats from the individual variable importance tables outputed by Webbs Module 2 Cubist script as modified by pz
### Calculate mean importance, standard deviation and rank and order by rank
###

# Kfold list (default 10)
Kfold <- c("Kfold_1", "Kfold_2", "Kfold_3", "Kfold_4", "Kfold_5", "Kfold_6", "Kfold_7", "Kfold_8", "Kfold_9", "Kfold_10")

for (a in 1:length(Attribute)){
  AttributeDirectory=paste(ModelRunDirectory,"//", (Attribute[a]), sep="") 
  setwd(AttributeDirectory)
  Depths <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE) #depths modelled
  if(length(Depths) == 6){ Depths <- Depths[c(1, 5, 3, 4, 6, 2)]}
  Allvi <- as.data.frame(Covariate)
  
  for (b in 1:length(Depths)){
    DepthDirectory=paste(AttributeDirectory,"//", (Depths[b]), sep="")
    setwd(DepthDirectory)
    
    #Import variable importsance from each kfold
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_1_CubistVariableImp.csv", sep = "")    
    Kfold1 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_2_CubistVariableImp.csv", sep = "")    
    Kfold2 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_3_CubistVariableImp.csv", sep = "")    
    Kfold3 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_4_CubistVariableImp.csv", sep = "")    
    Kfold4 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_5_CubistVariableImp.csv", sep = "")    
    Kfold5 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_6_CubistVariableImp.csv", sep = "")    
    Kfold6 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_7_CubistVariableImp.csv", sep = "")    
    Kfold7 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_8_CubistVariableImp.csv", sep = "")    
    Kfold8 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_9_CubistVariableImp.csv", sep = "")    
    Kfold9 <- read.table(Filename, header=TRUE, sep=",")
    Filename = paste((Attribute[a]), (Depths[b]), "_Kfold_10_CubistVariableImp.csv", sep = "")    
    Kfold10 <- read.table(Filename, header=TRUE, sep=",")
    
    # aggregate each kfold variable importance into one table, calculate average, reorder columns and export csv
    m1 <- merge(Kfold1, Kfold2, by = "row.names")
    m2 <- merge(Kfold3, Kfold4, by = "row.names")
    m3 <- merge(Kfold5, Kfold6, by = "row.names")
    m4 <- merge(Kfold7, Kfold8, by = "row.names")
    m5 <- merge(Kfold9, Kfold10, by = "row.names")
    vi <- join_all(list(m1,m2, m3, m4, m5), by = "Row.names")  
    vi$Average <- rowMeans(vi[,2:11]) #Calculate mean Importance
    vi$SD <- rowSds(as.matrix(vi[,2:11])) #Calculate standard deviation across all 10 Kfolds
    vi$Rank <- (vi[,12]/vi[,13]) #Calculate Rank based on mean/sd
    vi <- vi[order(-vi$Rank),] #Order according to rank
    #vi <- na.omit(vi)
    names(vi) <- c("Covariate", "Kfold1", "Kfold2", "Kfold3", "Kfold4", "Kfold5", "Kfold6", "Kfold7", "Kfold8", "Kfold9", "Kfold10", "Average", "SD", "Rank") 
    Filename = paste((Attribute[a]), "average_variable_importance.csv", sep = "")    
    write.table(vi, file = Filename, sep = ",", row.names = FALSE, col.names = TRUE)
    vii <- vi[,c(1,12)]
    Allvi <- merge(Allvi, vii, by = "Covariate")
  }
  
  #Create depth summary table for variable importances 
  if(length(Depths) > 1){
    Allvi$Average <- rowMeans(Allvi[,2:(length(Depths)+1)]) #Calculate mean accross mean depth importances
    Allvi$SD <- rowSds(as.matrix(Allvi[,2:(length(Depths)+1)])) #Calculate standard deviation accross mean depth importances
    Allvi$Rank <- (Allvi[,(length(Depths)+2)]/Allvi[,(length(Depths)+3)]) #Calulate Rank based on mean/sd
    Allvi <- Allvi[order(-Allvi$Rank),] #Order according to rank
    if(length(Depths)==6){names(Allvi) <- c("Covariate", "0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm", "Average", "SD", "Rank") }
    setwd(AttributeDirectory)
    Filename = paste((Attribute[a]), "average_acorss_depth_importance.csv", sep = "")
    write.table(Allvi, file = Filename, sep = ",", row.names = FALSE, col.names = TRUE)
  }
}
### END OF PART 2

### PART 3 - Create range rasters - ####

library(rgdal)
library(sp)
library(raster)

# Set temp directory for raster package
TempDirectory = paste(ProjectDir, ModelRun, "//TrainingDatatemp", sep="") 
rasterOptions(tmpdir=TempDirectory)

# Create individual range (uncertainity) rasters
for (a in 1:length(Attribute)){
  AttributeDirectory=paste(ModelRunDirectory,"//", (Attribute[a]), sep="") 
  setwd(AttributeDirectory)
  Depths <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE) #depths modelled
  
  for (b in 1:length(Depths)){
    Filename = paste((Attribute[a]), (Depths[b]), "_upper_mosaic.tif", sep = "")    
    ur <- raster(Filename)
    Filename = paste((Attribute[a]), (Depths[b]), "_lower_mosaic.tif", sep = "")    
    lr <- raster(Filename)
    uncertainity <- overlay(ur, lr, fun=function(r1, r2){return(r1-r2)})
    stduncert <- uncertainity/maxValue(uncertainity)
    writeRaster(stduncert, filename = paste(Attribute[a], Depths[b], "_uncertainity", ".tif", sep=""), format = "GTiff", overwrite = TRUE)
    }
  }
  
# Create a "mean of uncertainty raster" across depths 
for (a in 1:length(Attribute)){
  AttributeDirectory=paste(ModelRunDirectory,"//", (Attribute[a]), sep="") 
  setwd(AttributeDirectory)
  Depths <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE) #depths modelled
  rangeStack <- stack()
  for (b in 1:length(Depths)){
    Filename = paste((Attribute[a]), (Depths[b]), "_uncertainity.tif", sep = "")    
    tempraster <- raster(Filename)
    rangeStack <- stack(rangeStack, tempraster)
    }
  meanrange <- calc(rangeStack, mean)
  writeRaster(meanrange, filename = paste(Attribute[a],"average_uncertainity.tif", sep=""), format = "GTiff", overwrite = TRUE)
  }

# Create a "mean of uncertainty raster" across soil attributes
meanStack <- stack()
for (a in 1:length(Attribute)){
  AttributeDirectory=paste(ModelRunDirectory,"//", (Attribute[a]), sep="") 
  setwd(AttributeDirectory)
  Filename = paste((Attribute[a]),"average_uncertainity.tif", sep = "")
  tempraster <- raster(Filename)
  meanStack <- stack(meanStack, tempraster)
  }
meanrange <- calc(meanStack, mean)
setwd(ModelRunDirectoryUpper)
writeRaster(meanrange, filename = paste(ModelRun, "_average_uncertainity.tif", sep=""), format = "GTiff", overwrite = TRUE)

#Remove temp raster files 
rasfiles<- list.files(TempDirectory, pattern='$')
for(u in 1: length(rasfiles)){
  delete<-paste(TempDirectory,"//",rasfiles[u],sep="")
  file.remove(delete)
}

### END OF PART 3
###
### END OF THIS CODE

