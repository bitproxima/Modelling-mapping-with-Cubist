####################################################
## Step 5. Version 1 - 13/3/2017
##
## This script is Step 5. in a multi-stepped process to extract SALI lab data used in Cubist modelling of soil attributes (SAMFTAGA_RT_Kfolds_windows_Modlue1.R) for the project 'Mapping soil erodibility in the Fitzroy'
##
## Summary of all steps
## Scripts are located in \\lands\data\DSITI\LandSciences\NAS\slr_soils\Projects\Fitzroy_ErodSoil\Modelling\SiteData\scripts
## For each new run
## 1. Run Inserts.sql - Outputs Anew_depths.csv & Bnew_depths.csv
## For each soil attribute (Clay, Silt, FS, CS, ESP, Salinity, Cat_Ca, Cat_Mg, Clay_Act)
## 2. Run SiteData.sql - Output labdata.csv
## 3. Run DataPrepSpline.R - Output soilattribute_with_inserts.csv 
## 4. Run SplineTool.exe - Output soilattribute_stdout.txt
## 5. Run DataPrepCubist.R - Output soilattribute.csv
##
## Summary of what this script does
##    Pivots Spline output
##    Adds spatial references lost when the data was put through the spline
##    Change no data value from -9999 to NA    
##    
## Required user input
##    Line 32 - adjust filename of input txt with attribute name e.g. file = "soilattribute_stdout.txt"
##    Line 59 - adjust filename of output csv with attribute name e.g. file = "soilattribute.csv"
##
######################################################

## CODE STARTS HERE

library(plyr)
library(reshape)
setwd("//lands/data/DSITI/LandSciences/NAS/slr_soils/Projects/Fitzroy_ErodSoil/Modelling/SiteData/Harmonized_data") #set working directory
splineoutput<-read.table("CS_stdout.txt",header=TRUE,sep=",") #read all the data in
a <- subset(splineoutput, UpperDepth == 0 & LowerDepth == 5, select = c(Id, Value))
a <- rename(a, c(Value = "x0to5cm"))
b <- subset(splineoutput, UpperDepth == 5 & LowerDepth == 15, select = c(Id, Value))
b <- rename(b, c(Value = "x5to15cm"))
c <- subset(splineoutput, UpperDepth == 15 & LowerDepth == 30, select = c(Id, Value))
c <- rename(c, c(Value = "x15to30cm"))
d <- subset(splineoutput, UpperDepth == 30 & LowerDepth == 60, select = c(Id, Value))
d <- rename(d, c(Value = "x30to60cm"))
e <- subset(splineoutput, UpperDepth == 60 & LowerDepth == 100, select = c(Id, Value))
e <- rename(e, c(Value = "x60to100cm"))
f <- subset(splineoutput, UpperDepth == 100 & LowerDepth == 200, select = c(Id, Value))
f <- rename(f, c(Value = "x100to200cm"))
ff <- cbind(a,b[2], c[2], d[2], e[2], f[2])
ff <- rename(ff, c(Id = "ID"))
cords<- read.table("labdata.csv", header = TRUE, sep = ",")
cords <- subset(cords, select = c(ID, LATITUDE, LONGITUDE))
fff <- join(ff, cords, by = c("ID"), type = "left", match = "first")
fff <- fff[c(9,8,1,2,3,4,5,6,7)]
fff <- rename(fff, c(LATITUDE = "Y",LONGITUDE = "X" ))
fff[fff==-9999] <- NA
fff$x0to5cm[fff$x0to5cm <0] <- NA
fff$x5to15cm[fff$x5to15cm <0] <- NA
fff$x15to30cm[fff$x15to30cm <0] <- NA
fff$x30to60cm[fff$x30to60cm <0] <- NA
fff$x60to100cm[fff$x60to100cm <0] <- NA
fff$x100to200cm[fff$x100to200cm <0] <- NA
write.table(fff, "CS.csv", sep = ",", row.names = FALSE)

### END OF CODE
###
### END OF DATA PREP FOR MODELLING
