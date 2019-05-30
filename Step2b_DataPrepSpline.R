####################################################
## Step 3. - Version 1 - 13/3/2017 
##
## This script is Step 3. in a multi-stepped process to extract SALI lab data used in Cubist modelling of soil attributes (SAMFTAGA_RT_Kfolds_windows_Modlue1.R) for the project 'Mapping soil erodibility in the Fitzroy'
##
## Summary of all steps
## Scripts are located in \\lands\data\DSITI\LandSciences\NAS\slr_soils\Projects\Fitzroy_ErodSoil\Modelling\SiteData\scripts
## For each new data extraction run
## 1. Run Inserts.sql - Outputs Anew_depths.csv & Bnew_depths.csv
## For each soil attribute (Clay, Silt, FS, CS, ESP, Salinity, Cat_Ca, Cat_Mg, Clay_Act)
## 2. Run SiteData.sql - Output labdata.csv
## 3. Run DataPrepSpline.R - Output soilattribute_with_inserts.csv 
## 4. Run SplineTool.exe - Output soilattribute_stdout.txt
## 5. Run DataPrepCubist.R - Output soilattribute.csv
##
## Summary of what this script does
##    Merges the fictous depths created in Step 1 with the real lab data extracted in Step 2
##    Attributes each fictous depth with a value for the particular soil attribute taken from the nearest real sample within the same design master horizon
##    Excludes sites with just one sample
##    Removes samples that are duplicated or overlap
##
## Required user input
##    Line 71 - adjust filename of output csv with attribute name e.g. file = "soilattribute_with_inserts.csv"
##
######################################################

## CODE STARTS HERE

setwd("//lands/data/DSITI/LandSciences/NAS/slr_soils/Projects/Fitzroy_ErodSoil/Modelling/SiteData/Harmonized_data") #set working directory
lab.data <- read.table("labdata.csv",header=TRUE,sep=",") #read all the data in
library(reshape)
library(plyr)
labdata1 <- rename(lab.data, c(UD = "SUD", LD = "SLD")) #rename sample upper depth and sample lower depth

##Remove sites with only one analysed depth
count1 <- count(labdata1, vars = "ID")
morethan1 <- subset(count1, freq > 1)
labdata <- merge(labdata1, morethan1, by = c("ID"))

## A inserts
Inserts <- read.table("Anew_depths.csv",header=TRUE,sep=",") #Import Anew_depths.csv data generated in SQLDev for each duplex soil 
Allvalue <- merge(Inserts, labdata, by = c("PROJECT_CODE", "SITE_ID")) #Add lab data to new_depths (one to one)
Allvalue$diff <- Allvalue$UD - Allvalue$SUD #Add a 'diff' field to Allvalue df which is the 'difference between sample depth and A/B horizon depth'
Value <- subset(Allvalue, diff >= 0 & UD >= 10, select = c( ID, UD, LD, VALUE, diff)) #Limit records to lab samples in A horizon using 'diff' field and only where A horizon >= 10cm depth
Nearsamp<-aggregate(x = Value$diff, by = list(ID = Value$ID), min) #Limit records to sample nearest to A/B change
Nearsamp1 <- rename(Nearsamp, c(x = "diff")) #Rename column 'x' to 'diff' in 'Nearsamp' df to conform to 'Value' df
Aresult <- join(Value, Nearsamp1, by = c("ID", "diff"), type = "right", match = "first") #Combine df with nearest sample 'nearsamp' with df with lab data 'Value'
Aresult <- subset(Aresult, select = c(ID, UD, LD, VALUE)) #Remove 'diff' in prep for Spline program

## B inserts
Inserts <- read.table("Bnew_depths.csv",header=TRUE,sep=",") #Import Bnew_depths.csv data generated in SQLDev for each duplex soil 
Allvalue <- merge(Inserts, labdata, by = c("PROJECT_CODE", "SITE_ID")) #Add lab data to new_depths (one to one)
Allvalue$diff <- Allvalue$LD - Allvalue$SLD #Add a 'diff' field to Allvalue df which is the 'difference between sample depth and A/B horizon depth'
Value <- subset(Allvalue, diff <= 0 & UD >= 10, select = c( ID, UD, LD, VALUE, diff)) #Limit records to lab samples in B horizon using 'diff' field and only where B horizon Upper Depth >= 10cm
Nearsamp <- aggregate(x = Value$diff, by = list(ID = Value$ID), max) #Limit records to sample nearest to A/B change
Nearsamp1 <- rename(Nearsamp, c(x = "diff")) #Rename column 'x' to 'diff' in 'Nearsamp' df to conform to 'Value' df
Bresult <- join(Value, Nearsamp1, by = c("ID", "diff"), type = "right", match = "first") #Combine df with nearest sample 'nearsamp' with df with lab data 'Value'
Bresult <- subset(Bresult, select = c(ID, UD, LD, VALUE)) #Remove 'diff' in prep for Spline program

## Real sample depths
Realdata <- na.omit(subset(labdata, select = c(ID, SUD, SLD, VALUE))) ##select only fields required for spline and omit nulls 
DupsRemoved <- unique(Realdata[,1:2]) ## delete duplicates or overlaps
Realdata <- join(Realdata, DupsRemoved, by = c("ID", "SUD"), type = "right", match = "first")
Realdata <- rename(Realdata, c(SUD = "UD", SLD = "LD")) #rename depth columns to conform with Aresult and Bresult df

##merge Real samples with A & B inserts
Result1 <- rbind(Realdata,Aresult,Bresult) #Append data from each result above into the one df
Result2 <- aggregate(x = Result1$VALUE, by = list(ID = Result1$ID, UD = Result1$UD, LD = Result1$LD), mean) #For identical sample depths with more than one measured value, take the average of these values
Result2 <- rename(Result2, c( x = "VALUE"))
Result <- Result2[order(Result2$ID, Result2$UD, Result2$LD),] #Order df records on ID, upper depth
write.table(Result, file = "TN_with_inserts.csv", sep = ",", col.names = TRUE, row.names = FALSE) #Export result to project directory

## CODE ENDS HERE