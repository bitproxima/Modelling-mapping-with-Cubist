#-------------------------------------------------------------------------------------------------------------------------------------------#
###SOIL ATTRIBUTE MODELLING USING CUBIST & K-FOLD CROSS VALIDATION INCLUDING UNCERTAINTY CALCULATION (ESTIMATION OF UPPER & LOWER LIMITS)####
#-------------------------------------------------------------------------------------------------------------------------------------------#

##SUMMARY: 
#Module 2 - Will map out the cubist rules -tile by tile- at each k-fold iteration as defined from module 1. The upper/lower predictions will also be mapped. The final surfaces (including the upper/lower prediction rasters) 
#are derived from  averaging or taking the median of all K-fold surfaces. The module will then mosaic all the final predictions and upper/lower prediction raster tiles to derive the final product. 

##SCRIPT SETTINGS##
#NB: The settings for CovariateDirectory, DataDirectory, rastext, unbi, extr and ptile will need to be the same as stipulated for module 1.
CovariateDirectory="C:\\Temp\\K3\\Covariates" #Set the Covariate directory
DataDirectory="C:\\Temp\\K3\\TrainingData"        #Set the TrainingData directory
DEM="relief_dem3s.tif.sdat"   #Insert name of the DEM or similar raster that best represent the covariate space extent. This will be used in configuring the tile extents. 
rastext=".sdat"  #Set correct extension name of the covariate datasets 
splitdim=5      #Set appropriate tile dimension number for mapping outputs by tiles (minimum is 2), e.g. 5 will produce a fishnet of 5x5 equal parts (25 tiles) and outputs will be mapped sequentially (tile by tile) until all 25 tiles are mapped. 
#This is required to split the modelling when mapping, for example, large areas at high resolutions. After each tile has been mapped, the outputs will be merged to produce a final output at the original mapping extent (i.e. the covariate space extent) 

##Cubist Control parameters:
unbi<-F  # logical (T/F); should unbiased rules be used? 
extr<-0  # a number between 0 and 100: since Cubist uses linear models, predictions can be outside of the outside of the range as seen in the training set. This parameter controls how much rule predictions are adjusted to be consistent with the training set.
ptile<-90 # Set the percentile to which the upper/lower limits are estimated. GSM specifies 90%.

output<- "median" #"mean" or "median" raster outputs 

####################################################################################################################################################################################################
###CODE BEGINS HERE#################################################################################################################################################################################
####################################################################################################################################################################################################

#Bring in Modules
library(raster)
library(sp)
library(epiR)
library(gstat)
library(Cubist)
library(snow)
library(cvTools)
library(rgdal)
library(doParallel)
library(foreach)
library(plyr)
library(data.table)
library(miscTools)
#library(bigmemory)
#library(biganalytics)

RootDirectory<-gsub(gsub("^.*/","",DataDirectory),"",DataDirectory)
dir.create(paste(RootDirectory,"temp",sep=""))
TempDirectory=paste(RootDirectory,"temp",sep="")
rasterOptions(tmpdir=TempDirectory)
write(paste0("TMPDIR = ","'",TempDirectory,"'"), file=file.path(Sys.getenv('R_USER'), '.Renviron'))
dir.create(paste(RootDirectory,"Outputs",sep=""))
OutputDirectory=paste(RootDirectory,"Outputs",sep="")

#concordance Function
ccc<- function(observed, predicted){
  mx=mean(observed)
  my=mean(predicted)
  s2x=var(observed)
  s2y=var(predicted)
  sxy=mean((observed-mx)*(predicted-my))
  ccc=2*sxy/(s2x+s2y+(mx-my)^2 )
  return(ccc)}

setwd(DataDirectory) 

#Obtain list of soil data files
DataFiles<- list.files(getwd(),  pattern=".csv", full.names=FALSE)

##############################################################################################################################################################################
##############################################################################################################################################################################
##MAPPING OF THE CUBIST RULES TILE BY TILE####################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

#Establish covariate space prediction tiles
#Derive extent of the covariates
setwd(CovariateDirectory)
extent_rast<-raster(DEM)
rastres<-res(extent_rast)[1] #Obtain raster solution for raster modelling below
x_diff<-xmax(extent_rast)-xmin(extent_rast)
y_diff<-ymax(extent_rast)-ymin(extent_rast)
x_crop<-x_diff/splitdim
y_crop<-y_diff/splitdim
x_min<-xmin(extent_rast)
x_max<-x_min+x_crop

#Determine tiles that can be discarded, i.e. NA tiles
#Tiles that are relevant will have their extent parameters written into a matrix to be used later on when determining upper/lower limits with multiple rules.
crop_mat<-data.frame(matrix(NA,ncol=4,nrow=splitdim*splitdim))
cnt<-0
for(c in 1:splitdim){
  y_min<-ymin(extent_rast)
  y_max<-y_min+y_crop
  for(d in 1:splitdim){
    cnt<-cnt+1
    #Import covariates as a covariate stack
    setwd(CovariateDirectory) 
    extent_rast_extent<-extent(x_min,x_max,y_min,y_max)
    extent_rast_crop<-crop(extent_rast,extent_rast_extent)
    extent_rast_crop_max<-cellStats(extent_rast_crop, stat="max", na.rm=TRUE)
    if(is.infinite(extent_rast_crop_max)){
    }else{
      crop_mat[cnt,1]<-x_min
      crop_mat[cnt,2]<-x_max
      crop_mat[cnt,3]<-y_min
      crop_mat[cnt,4]<-y_max
    }
    y_min<-y_max
    y_max<-y_max+y_crop
  }
  x_min<-x_max
  x_max<-x_max+x_crop
}
#Drop NA tiles
crop_mat<-na.omit(crop_mat)

#For every tile, map out the cubist rules per data file, then by k fold, followed by variable name. 
for(y in 1:nrow(crop_mat)){
  #Determine cropping extent and crop the covariate stack to this extent
  covstack_extent<-extent(crop_mat[y,1],crop_mat[y,2],crop_mat[y,3],crop_mat[y,4])
  print(paste0("Cropping extent... tile=",y))
  
  setwd(CovariateDirectory) 
  list.files(getwd(),  pattern=paste0(gsub("[.]","",rastext),"$"), full.names=FALSE)
  files<- list.files(getwd(), pattern=paste0(gsub("[.]","",rastext),"$"))
  r1<- raster(files[1])
  for(i in 2:length(files)){
    r1<- stack(r1,files[i])}
  covstack_crop<-crop(r1,covstack_extent)
  
  #Convert the cropped covariate stack to pnts so we can determing rules governing each raster cells
  covstack_crop_pnts<-data.frame(na.omit(rasterToPoints(covstack_crop)))
  covstack_crop_pnts$ID<-row.names(covstack_crop_pnts)
  
  #Iterate over all soil data files in Data Directory
  for (g in 1:length(DataFiles)){
    #Read in data from soil data file to obtain variable names
    setwd(DataDirectory)
    SoilVar<- gsub(substr(DataFiles[g], nchar(DataFiles[g])-3, nchar(DataFiles[g])),"_",DataFiles[g])
    data<- read.table(DataFiles[g], sep=",",header=TRUE) # all data
    diagnam<-names(data)[4:ncol(data)]
    
    #Determine directories
    if(file.exists(paste((OutputDirectory),"/",SoilVar,sep=""))){
      OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
    }else{
      dir.create(paste((OutputDirectory),"/",SoilVar,sep=""))
      OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
    }
    
    #Iterate over each K fold rasters
    for(q in 1:10){
      #Iterate over each soil depth/variable
      for(k in 1:(ncol(data)-3)){
        Tvar<- names(data)[k+3]
        
        #Determine directorie
        OutputDirectory2=paste((OutputDirectory),"/",SoilVar,"/",Tvar,sep="")
        setwd(OutputDirectory2)
        #Import training data, upper/lower limits, cubist rules and submodel equations
        mod.dat<- read.table(paste(SoilVar,Tvar,"_Kfold_",q,"_Training.csv",sep=""), sep=",",header=TRUE) # all data
        r.ulPI_mat<- read.table(paste(SoilVar,Tvar,"_Kfold_",q,"_UpperLowerLimits.csv",sep=""), sep=",",header=TRUE) # all data
        cubrule_mat<- read.table(paste(SoilVar,Tvar,"_Kfold_",q,"_CubistRules.csv",sep=""), sep=",",header=TRUE) # all data
        cubcoef_mat<-read.table(paste(SoilVar,Tvar,"_Kfold_",q,"_Submodels.csv",sep=""), sep=",",header=TRUE) # all data
        
        #Produce cubist model from training data
        #Set target vaiable calls
        TvarC<-paste("mod.dat$",Tvar,sep="")
        target.C<- eval(parse(text=TvarC)) 
        cubistPred_drain<-cubist(x= mod.dat[,5:ncol(mod.dat)], y=target.C,cubistControl(unbiased = unbi,extrapolation = extr, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
        
        #Depending on the number of cubist paritions, map out the cubist submodel equations.
        if(cubrule_mat[1,1]==-9999){ #If there is no Cubist partitons, map out the single submodel equation to the entire cropped covariate space 
          print("Mapping single rules...")
          
          # Make submodel rule into function and determine predicted values
          tbl<-data.frame(covstack_crop_pnts)
          cubcoef_mat[] <- lapply(cubcoef_mat, as.character)
          cubeq<-cubcoef_mat
          
          CubistSubRule<-gsub("[j,]","",toString(cubeq))
          CubistSubRule<-gsub("\\[|\\]","",CubistSubRule)
          
          predval<-eval(parse(text=CubistSubRule))
          #Remove extreme values, i.e. ensure min and max prediction values are constrained to the min/max of the training values
          extr_val<-(extr/100)*(max(mod.dat[,4]) - min(mod.dat[,4]))
          maxval<-max(mod.dat[,4])+extr_val
          minval<-min(mod.dat[,4])-extr_val
          predval<-ifelse(predval>maxval,maxval,ifelse(predval<minval,minval,predval))
          
          #Produce xyz data frame for mapping  
          tbl$FP<-as.numeric(predval)
          pred_mat<-cbind(tbl[1:2],tbl$FP)   
          
          ##Create sub folders to place predicted mapping output
          if(file.exists(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))){
            OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
          }else{
            dir.create(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))
            OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
          }
          if(file.exists(paste((OutputDirectory2),"/predicted/tile",y,sep=""))){
            PredictedDirectory=paste((OutputDirectory2),"/predicted/tile",y,sep="")
          }else{
            dir.create(paste((OutputDirectory2),"/predicted",sep=""))
            dir.create(paste((OutputDirectory2),"/predicted/tile",y,sep=""))
            PredictedDirectory=paste((OutputDirectory2),"/predicted/tile",y,sep="")
          }
          #Cubist Mapping#################################################################################################
          setwd(PredictedDirectory)
          ##Map out predicted values
          cubistMapName<-paste(SoilVar,Tvar,"_Kfold_",q,"_tile_",y,"_predicted",sep="")
          pred_mat_rast<-rasterFromXYZ(pred_mat, res=rastres, crs=NA, digits=2)
          cubistMap<-writeRaster(pred_mat_rast, filename =paste0(cubistMapName,".tif"), res=rastres, format = "GTiff", overwrite = TRUE)
          
          ##Place lower predictions in separate subfolder
          if(file.exists(paste((OutputDirectory2),"/lower/tile",y,sep=""))){
            LowerPredDirectory=paste((OutputDirectory2),"/lower/tile",y,sep="")
          }else{
            dir.create(paste((OutputDirectory2),"/lower",sep=""))
            dir.create(paste((OutputDirectory2),"/lower/tile",y,sep=""))
            LowerPredDirectory=paste((OutputDirectory2),"/lower/tile",y,sep="")
          }
          setwd(LowerPredDirectory)
          ##Map out lower predicted values using upper/lower matrix and predicted raster
          writeRaster(cubistMap+r.ulPI_mat[1,1], filename =paste(SoilVar,Tvar,"_Kfold_",q,"_tile_",y,"_lowerpredtest.tif",sep=""), res=rastres, format = "GTiff", overwrite = TRUE)
          
          ##Place upper predictions in separate subfolder
          if(file.exists(paste((OutputDirectory2),"/upper/tile",y,sep=""))){
            UpperPredDirectory=paste((OutputDirectory2),"/upper/tile",y,sep="")
          }else{
            dir.create(paste((OutputDirectory2),"/upper",sep=""))
            dir.create(paste((OutputDirectory2),"/upper/tile",y,sep=""))
            UpperPredDirectory=paste((OutputDirectory2),"/upper/tile",y,sep="")
          }
          setwd(UpperPredDirectory)
          ##Map out upper predicted values using upper/lower matrix and predicted raster
          writeRaster(cubistMap+r.ulPI_mat[1,2], filename =paste(SoilVar,Tvar,"_Kfold_",q,"_tile_",y,"_upperpredtest.tif",sep=""), res=rastres, format = "GTiff", overwrite = TRUE)
          #Remove temp raster files 
          rasfiles<- list.files(TempDirectory, pattern='$')
          for(u in 1: length(rasfiles)){
            delete<-paste(TempDirectory,"/",rasfiles[u],sep="")
            file.remove(delete)
          }
          gc()  
        }else{ #Each cubist rule to be segmented and the corresponding submodel equations to be initiated separately.  
          print("Processing multiple rule tiles...")
          
          #Set up initial covariate stack matrix. After each rule is processed, data will be rbinded to this.
          covstack_crop_pnts2<-covstack_crop_pnts[1,]
          covstack_crop_pnts2$FP<- -9999
          covstack_crop_pnts2$rule<- -9999
          covstack_crop_pnts2$FP_lowerPI=-9999
          covstack_crop_pnts2$FP_upperPI=-9999
          covstack_crop_pnts2$FP_lower<--9999
          covstack_crop_pnts2$FP_upper<--9999
          
          #For each rule, run corresponding submodel equation
          for(i in 1:nrow(cubrule_mat)){
            
            #Obtain relevant cubist rule partition and convert it to a function
            CubistRule<- function(tbl){
              cr <- eval(parse(text=toString(cubrule_mat[i,])))
              return(cr)
            }
            
            #Obtain relevant submodel equation 
            cubeq<-cubcoef_mat[i,1]
            
            #Set up tbl so it can be used by the cubist rule partition function to determine relevant raster cells to undergo prediction
            tbl<-data.frame(covstack_crop_pnts)
            tbl$FP<-0
            #Run cubist rule partition function, determine which cells/rows are relvant
            tbl$rule<-CubistRule(tbl)
            tbl<-tbl[tbl$rule==TRUE,]
            if(nrow(tbl)==0){
            } else{
              
              tbl$rule[which(tbl$rule==TRUE)]<-i
              
              #convert submodel equation to a function
              CubistSubRule<-gsub("[j,]","",toString(cubeq))
              CubistSubRule<-gsub("\\[|\\]","",CubistSubRule)
              
              #Run submodel function, determine predicted values
              predval<-eval(parse(text=CubistSubRule))
              #Remove extreme values, i.e. ensure min and max prediction values are constrained to the min/max of the training values
              extr_val<-(extr/100)*(max(mod.dat[,4]) - min(mod.dat[,4]))
              maxval<-max(mod.dat[,4])+extr_val
              minval<-min(mod.dat[,4])-extr_val
              predval<-ifelse(predval>maxval,maxval,ifelse(predval<minval,minval,predval))
              
              #Attach predicted values back to tbl
              tbl$FP<-as.numeric(predval)
              
              #Attach upper/lower limits
              tbl$FP_lowerPI=0
              tbl$FP_upperPI=0 
              tbl$FP_lowerPI[which(tbl$rule==i)]<- r.ulPI_mat[i,1]
              tbl$FP_upperPI[which(tbl$rule==i)]<- r.ulPI_mat[i,2]
              
              #Determine upper/lower limts prediction intervals
              tbl$FP_lower<-tbl$FP + tbl$FP_lowerPI
              tbl$FP_upper<-tbl$FP + tbl$FP_upperPI
              
              #Attach/rbind to the initial covariate stack matrix
              covstack_crop_pnts2<-rbind(covstack_crop_pnts2,tbl)
            }
          }
          #Clean up the newly amalgamated covariate stack matrix
          covstack_crop_pnts2<-covstack_crop_pnts2[-1,] 
          covstack_crop_pnts2<-na.omit(covstack_crop_pnts2)
          covstack_crop_pnts2$ID<-as.numeric(covstack_crop_pnts2$ID)
          
          #Make new data frame by removing the covariate values retaining only x, y, predictions and upper/lower predictions
          covstack_crop_pnts1<-covstack_crop_pnts2[order(covstack_crop_pnts2$ID),]
          covstack_crop_pnts1<-data.frame(cbind(covstack_crop_pnts1$x,covstack_crop_pnts1$y,covstack_crop_pnts1$ID,covstack_crop_pnts1$FP_lower,covstack_crop_pnts1$FP,covstack_crop_pnts1$FP_upper))
          colnames(covstack_crop_pnts1)<-c("x","y","ID","FP_lower","FP","FP_upper")
          
          #We need to determine cells/rows that are predicted by more than one rule/submodel equation. we do this by separating these cells/rows from the data frame.
          #Once determined we will average the predicted values and reamalgamate back to the data frame.
          covstack_crop_pnts1_id<-data.frame(covstack_crop_pnts1$ID)
          #Determine cell/rows that are duplicated (this will create a single column giving T/F)
          covstack_crop_pnts1i<-duplicated(covstack_crop_pnts1_id) | duplicated(covstack_crop_pnts1_id[nrow(covstack_crop_pnts1_id):1, ])[nrow(covstack_crop_pnts1_id):1]
          covstack_crop_pnts1<-cbind(covstack_crop_pnts1,covstack_crop_pnts1i) #reattach to master data frame
          #Split this data frame based on T/F
          covstack_crop_pnts1_sp<-split(covstack_crop_pnts1,covstack_crop_pnts1[,ncol(covstack_crop_pnts1)])
          #If FALSE means we do not need to do averaging
          covstack_crop_pnts4<-covstack_crop_pnts1_sp$'FALSE'  
          #In case all cells need to be averaged, we will derive a sinle row from the master data frame to attach the averaged cells after averaging
          if(is.null(covstack_crop_pnts4)){   
            covstack_crop_pnts4<-covstack_crop_pnts1[1,]
            covstack_crop_pnts4[1,]<-NA
          } else{}
          
          covstack_crop_pnts4<-covstack_crop_pnts4[,-ncol(covstack_crop_pnts4)] #remove T/F column
          #If TRUE means we need to do averaging
          covstack_crop_pnts3<-covstack_crop_pnts1_sp$'TRUE'
          if(is.null(covstack_crop_pnts3)){ #Incase there are no 'TRUE' cells, go straight to creating raster data frames
            avg_mat<-na.omit(covstack_crop_pnts4)
            pred_mat<-cbind(avg_mat[,1:2],avg_mat$FP)
            lower_mat<-cbind(avg_mat[,1:2],avg_mat$FP_lower)
            upper_mat<-cbind(avg_mat[,1:2],avg_mat$FP_upper)
          } else{ #If TRUE, cases with 2 or more rules are averaged
            covstack_crop_pnts3<-covstack_crop_pnts3[order(covstack_crop_pnts3$ID),-ncol(covstack_crop_pnts3)]            
            uniq<-unique(covstack_crop_pnts3[,1:3]) #derive unique ID's with X and Ys
            #Average out duplicated ID's (goverend by more than one rule). The function below was vary handy and speedy and found at:
            #(http://stats.stackexchange.com/questions/8225/how-to-summarize-data-by-group-in-r)
            covstack_crop_pnts3<-data.table(covstack_crop_pnts3)
            setkey(covstack_crop_pnts3,ID)
            FP_lower<- covstack_crop_pnts3[,list(mean=mean(FP_lower)),by=ID]
            FP<- covstack_crop_pnts3[,list(mean=mean(FP)),by=ID]
            FP_upper<- covstack_crop_pnts3[,list(mean=mean(FP_upper)),by=ID]
            
            #Below is a slower version that does not need the data.table package
            #FP_lower<-ddply(covstack_crop_pnts3,~ID,summarise,mean=mean(FP_lower))
            #FP<-ddply(covstack_crop_pnts3,~ID,summarise,mean=mean(FP))
            #FP_upper<-ddply(covstack_crop_pnts3,~ID,summarise,mean=mean(FP_upper))
            
            #Combine the unique ID's, X, Y and avaerged values (prediction, upper/lower limit predictions)
            uniq_c<-cbind(uniq,FP_lower[,2],FP[,2],FP_upper[,2])
            colnames(uniq_c)<-colnames(covstack_crop_pnts4)
            
            #Create raster data frames
            avg_mat<-rbind(covstack_crop_pnts4,uniq_c)
            avg_mat<-avg_mat[order(avg_mat$ID),]
            pred_mat<-cbind(avg_mat[,1:2],avg_mat$FP)
            lower_mat<-cbind(avg_mat[,1:2],avg_mat$FP_lower)
            upper_mat<-cbind(avg_mat[,1:2],avg_mat$FP_upper)
          }
          
          ##Create sub folders and map out predictions
          if(file.exists(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))){
            OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
          }else{
            dir.create(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))
            OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
          }
          
          #Cubist Mapping#################################################################################################
          print("Mapping multiple rules...")
          #Map out predicted cubist values
          pred_mat_rast<-rasterFromXYZ(pred_mat, res=rastres, crs=NA, digits=2)
          lower_mat_rast<-rasterFromXYZ(lower_mat, res=rastres, crs=NA, digits=2)
          upper_mat_rast<-rasterFromXYZ(upper_mat, res=rastres, crs=NA, digits=2)
          
          if(file.exists(paste((OutputDirectory2),"/predicted/tile",y,sep=""))){
            PredictedDirectory=paste((OutputDirectory2),"/predicted/tile",y,sep="")
          }else{
            dir.create(paste((OutputDirectory2),"/predicted",sep=""))
            dir.create(paste((OutputDirectory2),"/predicted/tile",y,sep=""))
            PredictedDirectory=paste((OutputDirectory2),"/predicted/tile",y,sep="")
          }
          setwd(PredictedDirectory)
          cubistMapName<-paste(SoilVar,Tvar,"_Kfold_",q,"_tile_",y,"_predicted",sep="")
          ##Map out predicted values
          writeRaster(pred_mat_rast, filename =paste0(cubistMapName,".tif"), res=rastres, format = "GTiff", overwrite = TRUE)
          
          ##Place lower predictions in separate subfolder
          if(file.exists(paste((OutputDirectory2),"/lower/tile",y,sep=""))){
            LowerPredDirectory=paste((OutputDirectory2),"/lower/tile",y,sep="")
          }else{
            dir.create(paste((OutputDirectory2),"/lower",sep=""))
            dir.create(paste((OutputDirectory2),"/lower/tile",y,sep=""))
            LowerPredDirectory=paste((OutputDirectory2),"/lower/tile",y,sep="")
          }
          setwd(LowerPredDirectory)
          ##Map out lower predicted values using upper/lower matrix and predicted raster
          writeRaster(lower_mat_rast, filename =paste(SoilVar,Tvar,"_Kfold_",q,"_tile_",y,"_lowerpred.tif",sep=""), res=rastres, format = "GTiff", overwrite = TRUE)
          
          ##Place upper predictions in separate subfolder
          if(file.exists(paste((OutputDirectory2),"/upper/tile",y,sep=""))){
            UpperPredDirectory=paste((OutputDirectory2),"/upper/tile",y,sep="")
          }else{
            dir.create(paste((OutputDirectory2),"/upper",sep=""))
            dir.create(paste((OutputDirectory2),"/upper/tile",y,sep=""))
            UpperPredDirectory=paste((OutputDirectory2),"/upper/tile",y,sep="")
          }
          setwd(UpperPredDirectory)
          ##Map out upper predicted values using upper/lower matrix and predicted raster
          writeRaster(upper_mat_rast, filename =paste(SoilVar,Tvar,"_Kfold_",q,"_tile_",y,"_upperpred.tif",sep=""), res=rastres, format = "GTiff", overwrite = TRUE)
          #Remove temp raster files 
          rasfiles<- list.files(TempDirectory, pattern='$')
          for(u in 1: length(rasfiles)){
            delete<-paste(TempDirectory,"/",rasfiles[u],sep="")
            file.remove(delete)
          }
          gc()
        }     #else bracket
      }       #variable name loop bracket
    }         #K fold loop bracket
    gc()
    
    if(output=="mean"){ #20150731
      ########################################
      #####Average out K-fold predictions#####
      ########################################
      #Predicted rasters
      for(p in diagnam){ #For each variable name within the Data file, average out the kfolds
        #Determine correct tile folder
        OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/predicted/tile",y,sep="")
        setwd(OutputDirectory3)
        #Place k fold rasters in stack
        files<- list.files(getwd(), pattern='.tif$')
        predstack<- raster(files[1])
        for(o in 2:length(files)){
          predstack<- stack(predstack,crop(raster(files[o]),covstack_extent))}
        #Obtain mean values
        f1 <- function(x) calc(x, mean)
        #Determine output folder
        PredictedDirectory=paste((OutputDirectorySoilVar),"/",p,"/predicted",sep="")
        setwd(PredictedDirectory)
        #Map out averaged values
        beginCluster()
        pred <- clusterR(predstack, fun=f1, filename =paste(SoilVar,p,"_tile_",y,"_predicted_mean.tif",sep=""),format="GTiff",progress="text",overwrite=T)
        endCluster()
        #Delete tile directory
        DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,"/predicted",sep="")
        setwd(DeleteDirectory)
        unlink(paste0("tile",y), recursive = TRUE, force = TRUE)
      } 
      #Lower predicted rasters
      for(p in diagnam){#For each variable name within the Data file, average out the kfolds
        OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/lower/tile",y,sep="")
        setwd(OutputDirectory3)
        #Place k fold rasters in stack
        files<- list.files(getwd(), pattern='.tif$')
        lowerstack<- raster(files[1])
        for(o in 2:length(files)){
          lowerstack<- stack(lowerstack,files[o])}
        #Obtain mean values
        f1 <- function(x) calc(x, mean)
        #Determine output folder
        lowerDirectory=paste((OutputDirectorySoilVar),"/",p,"/lower",sep="")
        setwd(lowerDirectory)
        #Map out averaged values
        beginCluster()
        lower_pred <- clusterR(lowerstack, fun=f1, filename =paste(SoilVar,p,"_tile_",y,"_lowerpred_mean.tif",sep=""),format="GTiff",progress="text",overwrite=T)
        endCluster()
        #Delete tile directory
        DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,"/lower",sep="")
        setwd(DeleteDirectory)
        unlink(paste0("tile",y), recursive = TRUE, force = TRUE)
      }  
      #Upper predicted rasters
      for(p in diagnam){#For each variable name within the Data file, average out the kfolds
        OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/upper/tile",y,sep="")
        setwd(OutputDirectory3)
        #Place k fold rasters in stack
        files<- list.files(getwd(), pattern='.tif$')
        upperstack<- raster(files[1])
        for(o in 2:length(files)){
          upperstack<- stack(upperstack,files[o])}
        #Obtain mean values
        f1 <- function(x) calc(x, mean)
        #Determine output folder
        upperDirectory=paste((OutputDirectorySoilVar),"/",p,"/upper",sep="")
        setwd(upperDirectory)
        #Map out averaged values
        beginCluster()
        upper_pred <- clusterR(upperstack, fun=f1, filename =paste(SoilVar,p,"_tile_",y,"_upperpred_mean.tif",sep=""),format="GTiff",progress="text",overwrite=T)
        endCluster()
        #Delete tile directory
        DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,"/upper",sep="")
        setwd(DeleteDirectory)
        unlink(paste0("tile",y), recursive = TRUE, force = TRUE)
      }  
    } else{ #20150731
      ########################################
      #####K-fold prediction median###########
      ########################################
      #Predicted rasters
      for(p in diagnam){ #For each variable name within the Data file, average out the kfolds
        #Determine correct tile folder
        OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/predicted/tile",y,sep="")
        setwd(OutputDirectory3)
        #Place k fold rasters in stack
        files<- list.files(getwd(), pattern='.tif$')
        predstack<- raster(files[1])
        for(o in 2:length(files)){
          predstack<- stack(predstack,crop(raster(files[o]),covstack_extent))}
        #Obtain mean values
        f1 <- function(x) calc(x, median)
        #Determine output folder
        PredictedDirectory=paste((OutputDirectorySoilVar),"/",p,"/predicted",sep="")
        setwd(PredictedDirectory)
        #Map out averaged values
        beginCluster()
        pred <- clusterR(predstack, fun=f1, filename =paste(SoilVar,p,"_tile_",y,"_predicted_median.tif",sep=""),format="GTiff",progress="text",overwrite=T)
        endCluster()
        #Delete tile directory
        DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,"/predicted",sep="")
        setwd(DeleteDirectory)
        unlink(paste0("tile",y), recursive = TRUE, force = TRUE)
      } 
      #Lower predicted rasters
      for(p in diagnam){#For each variable name within the Data file, average out the kfolds
        OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/lower/tile",y,sep="")
        setwd(OutputDirectory3)
        #Place k fold rasters in stack
        files<- list.files(getwd(), pattern='.tif$')
        lowerstack<- raster(files[1])
        for(o in 2:length(files)){
          lowerstack<- stack(lowerstack,files[o])}
        #Obtain mean values
        f1 <- function(x) calc(x, median)
        #Determine output folder
        lowerDirectory=paste((OutputDirectorySoilVar),"/",p,"/lower",sep="")
        setwd(lowerDirectory)
        #Map out averaged values
        beginCluster()
        lower_pred <- clusterR(lowerstack, fun=f1, filename =paste(SoilVar,p,"_tile_",y,"_lowerpred_median.tif",sep=""),format="GTiff",progress="text",overwrite=T)
        endCluster()
        #Delete tile directory
        DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,"/lower",sep="")
        setwd(DeleteDirectory)
        unlink(paste0("tile",y), recursive = TRUE, force = TRUE)
      }  
      #Upper predicted rasters
      for(p in diagnam){#For each variable name within the Data file, average out the kfolds
        OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/upper/tile",y,sep="")
        setwd(OutputDirectory3)
        #Place k fold rasters in stack
        files<- list.files(getwd(), pattern='.tif$')
        upperstack<- raster(files[1])
        for(o in 2:length(files)){
          upperstack<- stack(upperstack,files[o])}
        #Obtain mean values
        f1 <- function(x) calc(x, median)
        #Determine output folder
        upperDirectory=paste((OutputDirectorySoilVar),"/",p,"/upper",sep="")
        setwd(upperDirectory)
        #Map out averaged values
        beginCluster()
        upper_pred <- clusterR(upperstack, fun=f1, filename =paste(SoilVar,p,"_tile_",y,"_upperpred_median.tif",sep=""),format="GTiff",progress="text",overwrite=T)
        endCluster()
        #Delete tile directory
        DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,"/upper",sep="")
        setwd(DeleteDirectory)
        unlink(paste0("tile",y), recursive = TRUE, force = TRUE)
      } 
    } #mean/median output
  } #Data file loop
}   #Crop matrix loop

##############################################################################################################################################################################
##############################################################################################################################################################################
##MOSAIC RASTER OUTPUT TILES##################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

for (g in 1:length(DataFiles)){
  #Read in one raster covariate layer to which we will extend the resulting mosaic output to maintain consistency
  setwd(CovariateDirectory) 
  list.files(getwd(),  pattern=paste0(gsub("[.]","",rastext),"$"), full.names=FALSE)
  files<- list.files(getwd(), pattern=paste0(gsub("[.]","",rastext),"$"))
  r1<- raster(files[1])
  rastres<-res(r1)[1] #Obtain raster solution for raster modelling below
  #Read in data file to obtain names
  SoilVar<- gsub(substr(DataFiles[g], nchar(DataFiles[g])-3, nchar(DataFiles[g])),"_",DataFiles[g])
  setwd(DataDirectory)
  data<- read.table(DataFiles[g], sep=",",header=TRUE) # all data
  diagnam<-names(data)[4:ncol(data)]
  
  #Define correct directories
  if(file.exists(paste((OutputDirectory),"/",SoilVar,sep=""))){
    OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
  }else{
    dir.create(paste((OutputDirectory),"/",SoilVar,sep=""))
    OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
  }
  for(p in diagnam){ #For every variable name within a data file, mosaic the tiles within the predicted, upper, lower directories
    setwd(CovariateDirectory) 
    list.files(getwd(),  pattern=paste0(gsub("[.]","",rastext),"$"), full.names=FALSE)
    files<- list.files(getwd(), pattern=paste0(gsub("[.]","",rastext),"$"))
    r1<- raster(files[1])
    #Determine prediction directory
    PredDirectory=paste((OutputDirectorySoilVar),"/",p,"/predicted",sep="")
    setwd(PredDirectory)
    files<- list.files(getwd(), pattern='.tif$')
    pred_merge<-raster(files[1])
    for(i in 2:length(files)){
      pred_merge<-merge(pred_merge,raster(files[i]))
    }
    #Place mosaic in data file directory
    OutputDirectorySoilVar<-paste0(OutputDirectory,"/",SoilVar)
    setwd(OutputDirectorySoilVar)
    #Ensure mosaic is extended to the same extent as the covariates
    writeRaster(extend(pred_merge,r1), filename =paste(SoilVar,p,"_predicted_mosaic.tif",sep=""), res=rastres, format = "GTiff", overwrite = TRUE)
    #Remove temp raster files 
    rasfiles<- list.files(TempDirectory, pattern='$')
    for(u in 1: length(rasfiles)){
      delete<-paste(TempDirectory,"/",rasfiles[u],sep="")
      file.remove(delete)
    }
    gc()
    
    setwd(CovariateDirectory) 
    list.files(getwd(),  pattern=paste0(gsub("[.]","",rastext),"$"), full.names=FALSE)
    files<- list.files(getwd(), pattern=paste0(gsub("[.]","",rastext),"$"))
    r1<- raster(files[1])
    #Determine lower directory
    LowerDirectory=paste((OutputDirectorySoilVar),"/",p,"/lower",sep="")
    setwd(LowerDirectory)
    files<- list.files(getwd(), pattern='.tif$')
    pred_merge<-raster(files[1])
    for(i in 2:length(files)){
      pred_merge<-merge(pred_merge,raster(files[i]))
    }
    #Place mosaic in data file directory
    OutputDirectorySoilVar<-paste0(OutputDirectory,"/",SoilVar)
    setwd(OutputDirectorySoilVar)
    #Ensure mosaic is extended to the same extent as the covariates
    writeRaster(extend(pred_merge,r1), filename =paste(SoilVar,p,"_lower_mosaic.tif",sep=""), res=rastres, format = "GTiff", overwrite = TRUE)
    #Remove temp raster files 
    rasfiles<- list.files(TempDirectory, pattern='$')
    for(u in 1: length(rasfiles)){
      delete<-paste(TempDirectory,"/",rasfiles[u],sep="")
      file.remove(delete)
    }
    gc()
    
    setwd(CovariateDirectory) 
    list.files(getwd(),  pattern=paste0(gsub("[.]","",rastext),"$"), full.names=FALSE)
    files<- list.files(getwd(), pattern=paste0(gsub("[.]","",rastext),"$"))
    r1<- raster(files[1])
    #Determine upper directory
    UpperDirectory=paste((OutputDirectorySoilVar),"/",p,"/upper",sep="")
    setwd(UpperDirectory)
    files<- list.files(getwd(), pattern='.tif$')
    pred_merge<-raster(files[1])
    for(i in 2:length(files)){
      pred_merge<-merge(extend(pred_merge,r1),raster(files[i]))
    }
    #Place mosaic in data file directory
    OutputDirectorySoilVar<-paste0(OutputDirectory,"/",SoilVar)
    setwd(OutputDirectorySoilVar)
    #Ensure mosaic is extended to the same extent as the covariates
    writeRaster(pred_merge, filename =paste(SoilVar,p,"_upper_mosaic.tif",sep=""), res=rastres, format = "GTiff", overwrite = TRUE)
    #Remove temp raster files 
    rasfiles<- list.files(TempDirectory, pattern='$')
    for(u in 1: length(rasfiles)){
      delete<-paste(TempDirectory,"/",rasfiles[u],sep="")
      file.remove(delete)
    }
    gc()
  }      
}
setwd(RootDirectory)
unlink("temp", recursive = TRUE, force = TRUE)
setwd(Sys.getenv('R_USER'))
unlink(".Renviron", recursive = TRUE, force = TRUE)

####END OF SCRIPT####