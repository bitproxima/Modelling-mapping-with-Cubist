#   Module One - Version 2.3 - 3/11/2017
#   What changed? Changed write.csv for Modeldata to not add a row number column
#
#   Previous versions
#   Version 2.2 17/05/2017 - Changed directory where ModelData.csv is exported to .../Cubist/Job/
#   Version 2.1 23/02/2017 - Changed directory where ModelData.csv is exported to .../Cubist/Job/Covariates/
#   Version 2.0 ??/??/2017 - Changed directory where ModelData.csv is exported to .../Cubist/Job/TrainingDataOutputs/
#
#-------------------------------------------------------------------------------------------------------------------------------------------#
###SOIL ATTRIBUTE MODELLING USING CUBIST & K-FOLD CROSS VALIDATION INCLUDING UNCERTAINTY CALCULATION (ESTIMATION OF UPPER & LOWER LIMITS)####
#-------------------------------------------------------------------------------------------------------------------------------------------#

##SUMMARY: 
#Module 1 -  Will fit a CUBIST model, perform k-fold cross validation, and calculate the upper/lower prediction limits (using leave one-out cross validation) to the specified CI interval (normally set at 90%). 
#The outputs are written to file which include cross validation results as well as each cubist partition rule and linear sub-models. Directories to where outputs are written are created automatically within the script.
#The module and overall script can handle multiple training files, i.e. more than one .csv file can be placed in the training directory (This is considered batch processing). 
#Subsequently, each new input training file will have their corresponding output files placed in separate directories (refer to folder structure diagram the instructions document). 
#In addition, each training file may include additional target variables such as, for example, multiple pH depths:
# 	 X       	Y          ID     x0to5cm        x5to15cm      x15to30cm
#             540512  5377883   1     44.38902517  50.79748922  NA
#             540712  5376883   2     18.38972517  27.39758926  74.51045857
#             531512  5375983   3     8.389025179  29.19448962  NA
#             535112  5371183   4     51.98907517  65.294489234 84.91145859 
#As such, each target variable outputs files are also place in separate subfolder directories under the corresponding training folder directory.

##SCRIPT SETTINGS##
CovariateDirectory="C:\\Temp\\K3\\Covariates" #Set the Covariate directory
DataDirectory="C:\\Temp\\K3\\TrainingData"        #Set the TrainingData directory
rastext=".sdat"  #Set correct extension name of the covariate datasets 

##Cubist Control parameters:
unbi<-F  # logical (T/F); should unbiased rules be used? 
extr<-100  # a number between 0 and 100: since Cubist uses linear models, predictions can be outside of the outside of the range as seen in the training set. This parameter controls how much rule predictions are adjusted to be consistent with the training set.
ptile<-90 # Set the percentile to which the upper/lower limits are estimated. GSM specifies 90%.

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
library(caret)

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
##MODEL TESTING - DIAGNOSTICS AND UPPER/LOWER LIMIT DETERMINATION#############################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

#Iterate over all soil data files in Data Directory
for (g in 1:length(DataFiles)){
  
  #Read in data from soil data file 
  setwd(DataDirectory)
  data<- read.table(DataFiles[g], sep=",",header=TRUE) # all data
  
  #Define soil data name
  SoilVar<- gsub(substr(DataFiles[g], nchar(DataFiles[g])-3, nchar(DataFiles[g])),"_",DataFiles[g])
  
  #Import covariates as a covariate stack
  setwd(CovariateDirectory) 
  list.files(getwd(),  pattern=paste0(gsub("[.]","",rastext),"$"), full.names=FALSE)
  files<- list.files(getwd(), pattern=paste0(gsub("[.]","",rastext),"$"))
  r1<- raster(files[1])
  for(i in 2:length(files)){
    r1<- stack(r1,files[i])}
  rastres<-res(r1)[1] #Obtain raster solution for raster modelling below
  #Extract covariate values to the training data
  dat.df<-data.frame(data)
  dat.sdf<-SpatialPointsDataFrame(coords=dat.df[,1:2], data=dat.df)
  dat.sdf.ext<-extract(r1,dat.sdf)
  datacov<-cbind(data,dat.sdf.ext)
  
  if(file.exists(paste((OutputDirectory),"/",SoilVar,sep=""))){
    OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
  }else{
    dir.create(paste((OutputDirectory),"/",SoilVar,sep=""))
    OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
  }
  
  ##Create K-fold datasets to assess model uncertainty
  
  k<-10
  folds<- cvFolds(n=nrow(data), K=k, R = 1, type = "random") 
  ksize<- trunc(nrow(data)/10)
  kmat<- matrix(NA,nrow=ksize,ncol=10)
  for (i in 0:(k-1)){
    st= (ksize*i)+1
    en= (ksize*i)+ksize
    kmat[,i+1]<- folds$subsets[st:en]
  }
  if(ksize==nrow(data)/10){
  } else{
    kmat_ext<- matrix(folds$subsets[(en+1):(nrow(data))])
    kmat2<- matrix(NA,nrow=1,ncol=10)
    for (i in 1:nrow(kmat_ext)){
      kmat2[1,i]<-kmat_ext[i,1]
    }
    kmat<- rbind(kmat,kmat2)
  }
  
  #Create diagnostic matrix (to be written)
  diagnostic_mat<-matrix(NA,nrow=((ncol(data)-3)*k),ncol=11)
  colnames(diagnostic_mat)<-c("Kfold_level","VarName","Calib_RMSE","Calib_R2","Calib_Bias","Calib_CC","Valid_RMSE","Valid_R2","Valid_Bias","Valid_CC","Perc within UPP/LOW limits")
  cnt<-0
  for(q in 1:10){
    print(q)
    valdat<- kmat[,q]
    valdat<-na.omit(valdat)
    
    tradat_mat<- kmat[,-q]
    tradat<- c(tradat_mat[,1],tradat_mat[,2],tradat_mat[,3],tradat_mat[,4],tradat_mat[,5],tradat_mat[,6],tradat_mat[,7],tradat_mat[,8],tradat_mat[,9])
    tradat<-na.omit(tradat)
    #Iterate over each soil depth in data file
    for(k in 1:(ncol(data)-3)){
      cnt<-cnt+1
      diagnostic_mat[cnt,1]<-q
      #Define soil depth name
      Tvar<- names(data)[k+3]
      diagnostic_mat[cnt,2]<-Tvar
      
      #Define soil data depth training table
      soil.dat <- datacov[,c(1:3,k+3,(ncol(data)+1):ncol(datacov))]
      soil.dat <- na.omit(soil.dat)
      soil.dat <- soil.dat[which(soil.dat[,4]>(-9999) | soil.dat[,4]<(9999)), ]   
      soil.dat <- round(soil.dat,2)
      
      #Randomly sample training data/validation data
      mod.dat<- soil.dat[tradat,] #calibration set
      mod.dat <- na.omit(mod.dat)
      val.dat<- soil.dat[valdat,] #validation set
      val.dat <- na.omit(val.dat)
        
      #Set target vaiable calls
      TvarC<-paste("mod.dat$",Tvar,sep="")
      target.C<- eval(parse(text=TvarC))  
      TvarV<-paste("val.dat$",Tvar,sep="")
      target.V<-eval(parse(text=TvarV))
      
      #Produce cubist model and predict onto training and validation datasets
      cubistPred_drain<-cubist(x= mod.dat[,5:ncol(mod.dat)], y=target.C,cubistControl(unbiased = unbi,extrapolation = extr, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
      setwd(OutputDirectorySoilVar)
      mod.pred<- predict(cubistPred_drain, newdata = mod.dat,neighbors = 0, )
      mod.pred.V<- predict(cubistPred_drain, newdata = val.dat,neighbors = 0, )
      
      #DIOGNOSTICS (calibration)
      diagnostic_mat[cnt,3]<- sqrt(mean((target.C -mod.pred)^2))
      R2.cal<- lm(mod.pred ~target.C)
      diagnostic_mat[cnt,4]<-as.matrix(summary(R2.cal)$adj.r.squared)
      diagnostic_mat[cnt,5]<-mean(mod.pred)- mean(target.C)
      diagnostic_mat[cnt,6]<- ccc(target.C,mod.pred)
           
      #Add residuals to predictions
      val.dat$FP<- -9999 #mod.pred.V
      
      #sep up cubist rule iterations in matrix
      rules<-cubistPred_drain$splits[2:5]
      
      #If 'rules' result in NULL value indicates only one cubist rule established and the cubist model can be run for the entire dataset during LOCV
      ##otherwise each rule will be sequentially run in the 'else' argument 
      if(is.null(rules)){ 
        
        #Derive submodel rule and place into matrix
        tbl<-mod.dat
        cubcoef<-cubistPred_drain$coefficients
        cubcoef<-cubcoef[1,1:(ncol(cubcoef)-2)]
        cubcoef<-cubcoef[, colSums(is.na(cubcoef)) != nrow(cubcoef)]
        cubcoefcubcoef_int<-cubcoef[1,1]#######################################changed on 20150630#######################
        cubeq<-toString(cubcoef[1,1])
        for(f in 2:ncol(cubcoef)){
          cubeq<-paste(cubeq,"+(",cubcoef[1,f],"*tbl","[j,]","$",colnames(cubcoef)[f],")",sep="")
        }
        #Set up coeficients matrix
        cubcoef_mat<-data.frame(matrix(NA,ncol=1,nrow=1))
        colnames(cubcoef_mat)<-"SubModels"
        cubcoef_mat[1,1]<-cubeq
        
        #resid_mat<-data.frame(matrix(NA,nrow=nrow(tbl),ncol=1))
        #for(j in 1:nrow(tbl1)){
        #  resid_mat[j,1]<-tbl1[j,4]-eval(parse(text=cubeq))
        #}
        
        #RULE-BASED LOCV - this will estimate the Upper/Lower limits
        ########## Leave-one-out cross validation #####################################################################################################
        tbl=mod.dat
        looresfunc<- function(j){
          loocubistPred<-cubist(x= tbl[-j,5:(ncol(tbl)-1)], y=tbl[-j,4],cubistControl(unbiased = unbi,extrapolation = extr, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
          looPred<-predict( loocubistPred, newdata = tbl[j,],neighbors = 0)
          looRes<- tbl[j,4] - (looPred)#+int.resids1.r1) 
          return(looRes)
        }
        print("Single Rule. Running parallel processing function on leave one out cross validation loop...")
        cl<-makeCluster(detectCores())
        registerDoParallel(cl)
        looResiduals<-foreach(j=1:nrow(tbl), .packages=library(Cubist), .combine='rbind') %dopar%{looresfunc(j)}
        stopCluster(cl)
        ################################################################################################################################################
        looResiduals<-as.numeric(looResiduals)
        ptile_val<-((100-ptile)/2)/100
        lower_ptile<-ptile_val
        upper_ptile<-1-ptile_val
        r.ulPI<- quantile(looResiduals, probs = c(lower_ptile,upper_ptile), na.rm = FALSE,names = F, type = 7) # rule lower and upper PI
        r.ulPI_mat<-data.frame(matrix(NA,nrow=1,ncol=2))
        colnames(r.ulPI_mat)<-c("LowerLimitRange","UpperLimitRange")
        r.ulPI_mat[1,1]<-r.ulPI[1]
        r.ulPI_mat[1,2]<-r.ulPI[2]
        
        #Use the upper/lower limits on the validation matrix to assess the upper/lower limit accuracy accuracy and to write validation stats
        tbl=val.dat
        if(is.null(dim(cubcoef))){ #estimate predicted values using the submodel equation on all validation data
          tbl[,ncol(tbl)]<-cubcoef_int
        } else{
          extr_val<-(extr/100)*(max(mod.dat[,4]) - min(mod.dat[,4]))
          maxval<-max(mod.dat[,4])+extr_val
          minval<-min(mod.dat[,4])-extr_val
          for(j in 1:nrow(tbl)){
            predval<-eval(parse(text=cubeq))
            #Remove extreme values, i.e. ensure min and max prediction values are constrained to the min/max of the training values
            tbl[j,ncol(tbl)]<-ifelse(predval>maxval,maxval,ifelse(predval<minval,minval,predval))
          }
        }
        
        #calculate upper/lower limits on validation data and the upper/lower limits on the validation matrix to assess its accuracy          
        tbl$FP_lowerPI=r.ulPI[1]
        tbl$FP_upperPI=r.ulPI[2]
        tbl$FP_lower<-tbl$FP + tbl$FP_lowerPI
        tbl$FP_upper<-tbl$FP + tbl$FP_upperPI
        val.dat2<-tbl
        
        #Determine the proportion of predicted values that reside within the upper/lower limits within the validation dataset and write to diagnostic matrix.
        TvarV<-paste("val.dat2$",Tvar,sep="")
        target.V<-eval(parse(text=TvarV))  
        target.V.FP_lower<-eval(parse(text="val.dat2$FP_lower"))
        target.V.FP_upper<-eval(parse(text="val.dat2$FP_upper"))
        val.dat2$PICP<- as.numeric(target.V >= target.V.FP_lower & target.V <= target.V.FP_upper)
        #summary(as.factor(val.dat2$PICP))
        #diagnostic_mat[cnt,11]<-summary(as.factor(val.dat2$PICP))[2]/nrow(val.dat2)
        diagnostic_mat[cnt,11]<-sum(val.dat2$PICP)/nrow(val.dat2) #######################################changed on 20150630#######################
        #Derive diagnostics stats (validation)
        diagnostic_mat[cnt,7]<- sqrt(mean((target.V -val.dat2$FP)^2))
        R2.val<- lm(val.dat2$FP~target.V)
        diagnostic_mat[cnt,8]<-as.matrix(summary(R2.val)$adj.r.squared)
        diagnostic_mat[cnt,9]<-mean(val.dat2$FP)- mean(target.V)
        diagnostic_mat[cnt,10]<- ccc(target.V,val.dat2$FP)
        
        cubrule_mat<-data.frame(matrix(NA,ncol=1,nrow=1))
        cubrule_mat[,1]<--9999
        
      }else{ #Each rule to be segmented and run sequentially
        
        #Form the rules matrix
        rulesdf<- data.frame(rules)
        
        cubsum<-unlist(strsplit(as.character(summary(cubistPred_drain)),"\t"))
        
        cubsum1<-cubsum[1]
        for (p in 2:length(cubsum)){
          cubsum1<-paste(cubsum1,cubsum[p])
        }
        
        cubsum2<-gsub(" ","",cubsum1)   
        
        rulesdf_sp<-split(rulesdf,rulesdf$rule)
        rulesnum<-data.frame(unique(rulesdf[,1]))
        cnt1<-0
        for(o in 1:nrow(rulesnum)){
          rulesdf_sp1<-data.frame(rulesdf_sp[o])
          for (p in 1:nrow(rulesdf_sp1)){
            cnt1<-cnt1+1
            rule<- paste0(as.character(rulesdf_sp1[p,2]),as.character(rulesdf_sp1[p,3]))
            cubsum2i<-unlist(strsplit(cubsum2,paste0("Rule",o)))
            cubsum2ii<-cubsum2i[2]
            cubsum3<-unlist(strsplit(cubsum2ii,rule))
            cubsum4<-unlist(strsplit(cubsum3[2],"\n"))
            rulesdf[cnt1,4]<-as.numeric(cubsum4[1])
          }
        }
        
        #Set up validation matrix to where the upper and lower limit attributes will be written
        val.dat2<-val.dat[1,]
        val.dat2$rule<- -9999
        val.dat2$FP_lowerPI=-9999
        val.dat2$FP_upperPI=-9999
        val.dat2$FP_lower<--9999
        val.dat2$FP_upper<--9999
        
        #set up Upper/Lower value matrix
        numrules<-max(rules[,1])
        r.ulPI_mat<-matrix(NA,nrow=numrules,ncol=2)
        colnames(r.ulPI_mat)<-c("LowerLimitRange","UpperLimitRange")
        
        #Set up cubist rule matrix
        cubrule_mat<-data.frame(matrix(NA,ncol=1,nrow=numrules))
        colnames(cubrule_mat)<-"CubistRules"
        
        #Set up coeficients matrix
        cubcoef_mat<-data.frame(matrix(NA,ncol=1,nrow=numrules))
        colnames(cubcoef_mat)<-"SubModels"
        
        #For each rule derive partition rules and submodel rules aswell as the upper/lower limits (via LOOCV) and place in matrix.
        for(i in 1:numrules){
          
          #Derive a single cubist rule from the cubist rule matrix and manipulate it into a function
          cubrule<-rulesdf[rulesdf$rule==i,]
          cubrule<-cubrule[order(cubrule$variable),]
          cubrule<-data.frame(lapply(cubrule[,-1], as.character), stringsAsFactors=FALSE)
          
          #Ensure same name variables are together (eg. DEM<220 & DEM>230)
          cubrule_spl<-split(cubrule,cubrule$variable)
          cubrule_uni<-unique(cubrule$variable)
          cubrule_matr<-matrix(NA,nrow=length(cubrule_uni),ncol=1)
          for(j in 1:length(cubrule_uni)){
            uvar<- paste("cubrule_spl$",cubrule_uni[j],sep="")
            uvar1<-eval(parse(text=uvar))
            if(nrow(uvar1)>1){
              uvar1_mat<- matrix(NA,nrow=nrow(uvar1),ncol=1)
              for (f in 1:nrow(uvar1_mat)){
                uvar1_mat[f,1]<- paste("tbl$",uvar1[f,1],uvar1[f,2],uvar1[f,3], sep="")
              }
              uvar2<- ""
              for (f in 1:nrow(uvar1_mat)){
                uvar2<- paste(uvar2,"&",uvar1_mat[f,1])
              }
              uvar2<-substr(uvar2, 4, nchar(uvar2)) 
              cubrule_matr[j,1]<-paste("(",uvar2,")",sep="")
            } else{cubrule_matr[j,1]<-paste("tbl$",uvar1[1,1],uvar1[1,2],uvar1[1,3], sep="")}
          }
          
          #Make cubist rule into 1 line of query to form the function
          cubrule1<- ""
          for (j in 1:nrow(cubrule_matr)){
            cubrule1<- paste(cubrule1,"&",cubrule_matr[j,1])
          }
          
          #Write rule to matrix
          cubrule_mat[i,]<-substr(cubrule1, 4, nchar(cubrule1))
          
          #Produce the function
          CubistRule<- function(tbl){
            cr <- eval(parse(text=cubrule_mat[i,]))
            return(cr)
          }
          
          #Derive submodel rule and place into matrix
          tbl=mod.dat
          tbl$rule<-CubistRule(tbl)
          tbl<-tbl[tbl$rule==TRUE,]
          #mod.dat_ii<-tbl[tbl$rule==FALSE,-ncol(tbl)]
          cubcoef<-cubistPred_drain$coefficients
          cubcoef<-cubcoef[i,1:(ncol(cubcoef)-2)]
          cubcoef_int<-cubcoef[1,1]
          cubcoef<-cubcoef[, colSums(is.na(cubcoef)) != nrow(cubcoef)]
            if(is.null(dim(cubcoef))){
              #resid_mat<-data.frame(matrix(NA,nrow=nrow(tbl1),ncol=1))
              #resid_mat[,1]<-tbl1[,4]-cubcoef_int
              cubcoef_mat[i,1]<-cubcoef_int
            } else{
              cubeq<-toString(cubcoef[1,1])
              for(f in 2:ncol(cubcoef)){
                cubeq<-paste(cubeq,"+(",cubcoef[1,f],"*tbl","[j,]","$",colnames(cubcoef)[f],")",sep="")
              }
              cubcoef_mat[i,1]<-cubeq
              #resid_mat<-data.frame(matrix(NA,nrow=nrow(tbl1),ncol=1))
              #for(j in 1:nrow(tbl1)){
              #  resid_mat[j,1]<-tbl1[j,4]-eval(parse(text=cubeq))
              #}
            }
            
          #RULE-BASED LOCV - this will estimate the Upper/Lower limits
          ########## Leave-one-out cross validation #####################################################################################################
          looresfunc<- function(j){
            loocubistPred<-cubist(x= tbl[-j,5:(ncol(tbl)-1)], y=tbl[-j,4],cubistControl(unbiased = unbi,extrapolation = extr, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
            looPred<-predict( loocubistPred, newdata = tbl[j,],neighbors = 0)
            looRes<- tbl[j,4] - (looPred)#+int.resids1.r1) 
            return(looRes)
          }
          print("Multiple Rules. Running parallel processing function on leave one out cross validation loop...")
          cl<-makeCluster(detectCores())
          registerDoParallel(cl)
          #strt<-Sys.time()
          looResiduals<-foreach(j=1:nrow(tbl), .packages=library(Cubist), .combine='rbind') %dopar%{looresfunc(j)}
          #print(Sys.time()-strt)
          stopCluster(cl)
          #################################################################################################################################################
          looResiduals<-as.numeric(looResiduals)
          ptile_val<-((100-ptile)/2)/100
          lower_ptile<-ptile_val
          upper_ptile<-1-ptile_val
          r.ulPI<- quantile(looResiduals, probs = c(lower_ptile,upper_ptile), na.rm = FALSE,names = F, type = 7) # rule lower and upper PI
          r.ulPI_mat[i,1]<-r.ulPI[1]
          r.ulPI_mat[i,2]<-r.ulPI[2]
            
          #Use the upper/lower limits on the validation matrix to assess the upper/lower limit accuracy accuracy
          tbl=val.dat
          if(is.null(dim(cubcoef))){ #estimate predicted values using the submodel equation on all validation data
              tbl[,ncol(tbl)]<-cubcoef_int
            } else{
              extr_val<-(extr/100)*(max(mod.dat[,4]) - min(mod.dat[,4]))
              maxval<-max(mod.dat[,4])+extr_val
              minval<-min(mod.dat[,4])-extr_val
              for(j in 1:nrow(tbl)){
                predval<-eval(parse(text=cubeq))
                #Remove extreme values, i.e. ensure min and max prediction values are constrained to the min/max of the training values
                tbl[j,ncol(tbl)]<-ifelse(predval>maxval,maxval,ifelse(predval<minval,minval,predval))
              }
            }
            tbl<-tbl
            tbl$rule<-CubistRule(tbl) #derive correct validation values according to the cubist parition rule
            tbl1<-tbl[tbl$rule==TRUE,]
            tbl1$rule[which(tbl1$rule==TRUE)]<-i
            
            if(nrow(tbl1)==0){ #calculate upper/lower limits on validation data
            } else {
              
              tbl1$FP_lowerPI=0
              tbl1$FP_upperPI=0 
              
              tbl1$FP_lowerPI[which(tbl1$rule==i)]<- r.ulPI[1]
              tbl1$FP_upperPI[which(tbl1$rule==i)]<- r.ulPI[2]
              
              tbl1$FP_lower<-tbl1$FP + tbl1$FP_lowerPI
              tbl1$FP_upper<-tbl1$FP + tbl1$FP_upperPI
              
              val.dat2<-rbind(val.dat2,tbl1)
            }
            if (nrow(mod.dat)==0) break
          }
        
        r.ulPI_mat<-na.omit(r.ulPI_mat)
        cubrule_mat<-na.omit(cubrule_mat)
          
        #Determine the proportion of predicted values that reside within the upper/lower limits within the validation dataset and write to diagnostic matrix.
        val.dat2<-val.dat2[-1,] 
        TvarV<-paste("val.dat2$",Tvar,sep="")
        target.V<-eval(parse(text=TvarV))  
        target.V.FP_lower<-eval(parse(text="val.dat2$FP_lower"))
        target.V.FP_upper<-eval(parse(text="val.dat2$FP_upper"))
        val.dat2$PICP<- as.numeric(target.V >= target.V.FP_lower & target.V <= target.V.FP_upper)
          #summary(as.factor(val.dat2$PICP))
        #Where more than two rules dictate a case; average out the values including the predicted and upper/lower limit values
        val.dat1<-val.dat2[order(val.dat2[,3]),]
        val.dat3<-val.dat1[1,]
        val.dat4<-val.dat1
        for(h in 1:(nrow(val.dat1)-1)){ #Determine case where more than two rules are applied
          #ifelse(val.dat1[h,3]==val.dat1[h+1,3],val.dat3<-rbind(val.dat3,val.dat1[h,],val.dat1[h+1,]),val.dat4<-rbind(val.dat4,val.dat1[h,]))
          if(val.dat1[h,3]==val.dat1[h+1,3]){
            val.dat3<-rbind(val.dat3,val.dat1[h,],val.dat1[h+1,])
            val.dat4[h,3]<-NA
            val.dat4[h+1,3]<-NA
          } else{}
        }
        val.dat4<-na.omit(val.dat4)
        val.dat3<-val.dat3[-1,]
        if(nrow(val.dat3)==0){ 
          #diagnostic_mat[cnt,11]<-summary(as.factor(val.dat2$PICP))[2]/nrow(val.dat2) #if no cases with 2 or more rules are found
          diagnostic_mat[cnt,11]<-sum(val.dat2$PICP)/nrow(val.dat2) #######################################changed on 20150630#######################
        } else{ #if cases with 2 or more rules are found
          uid<-unique(val.dat3$ID)
          val.dat3_sp<-split(val.dat3,val.dat3$ID)
          for(h in uid){ #Average out the values and combine with new val.dat matrix
            pas<- paste("val.dat3_sp$'",h,"'",sep="")
            valdatspl<-eval(parse(text=pas))
            valdatspl1<-valdatspl[1,]
            valdatspl1[1,(ncol(valdatspl)-7):ncol(valdatspl)]<-colMeans(valdatspl[,(ncol(valdatspl)-7):ncol(valdatspl)])
            val.dat4<-rbind(val.dat4,valdatspl1[1,])
          }
          #Determine % within up/lo
          TvarV<-paste("val.dat4$",Tvar,sep="")
          target.V<-eval(parse(text=TvarV))  
          target.V.FP_lower<-eval(parse(text="val.dat4$FP_lower"))
          target.V.FP_upper<-eval(parse(text="val.dat4$FP_upper"))
          val.dat4[,ncol(val.dat4)]<- as.numeric(target.V >= target.V.FP_lower & target.V <= target.V.FP_upper)
          #diagnostic_mat[cnt,11]<-summary(as.factor(val.dat4$PICP))[2]/nrow(val.dat4)
          diagnostic_mat[cnt,11]<-sum(val.dat4$PICP)/nrow(val.dat4) #######################################changed on 20150630#######################
        }
        
        #Derive diagnostics stats (validation)
        TvarV<-paste("val.dat4$",Tvar,sep="")
        target.V<-eval(parse(text=TvarV))
        diagnostic_mat[cnt,7]<- sqrt(mean((target.V -val.dat4$FP)^2))
        R2.val<- lm(val.dat4$FP~target.V)
        diagnostic_mat[cnt,8]<-as.matrix(summary(R2.val)$adj.r.squared)
        diagnostic_mat[cnt,9]<-mean(val.dat4$FP)- mean(target.V)
        diagnostic_mat[cnt,10]<- ccc(target.V,val.dat4$FP)
        
        setwd(OutputDirectorySoilVar)
      }
      
      ##Create sub folders and write out model results and diagnostics
      if(file.exists(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))){
        OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
      }else{
        dir.create(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))
        OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
      }
      
      setwd(OutputDirectory2)
      if(is.null(rules)){}else{write.table(rulesdf,file=paste(SoilVar,Tvar,"_Kfold_",q,"_CubistPartitions.csv",sep=""),sep=",", col.names=T,row.names=F)}
      write.table(cubcoef_mat,file=paste(SoilVar,Tvar,"_Kfold_",q,"_Submodels.csv",sep=""),sep=",", col.names=T,row.names=F)
      write.table(cubistPred_drain$usage,file=paste(SoilVar,Tvar,"_Kfold_",q,"_CubistVariableUsage.csv",sep=""),sep=",", col.names=T,row.names=F)
      diagnostic_mat2<-matrix(NA,nrow=1,ncol=11)
      write.table((varImp(cubistPred_drain)), file=paste(SoilVar,Tvar,"_Kfold_",q,"_CubistVariableImp.csv",sep=""),sep=",", col.names=T,row.names=T)
      colnames(diagnostic_mat2)<-c("Kfold_level","VarName","Calib_RMSE","Calib_R2","Calib_Bias","Calib_CC","Valid_RMSE","Valid_R2","Valid_Bias","Valid_CC","Perc within UPP/LOW limits")
      diagnostic_mat2[1,]<-diagnostic_mat[cnt,]
      write.table(diagnostic_mat2,file=paste(SoilVar,Tvar,"_Kfold_",q,"_Diagnostics.csv",sep=""),sep=",", col.names=T,row.names=F)
      write.table(mod.dat,file=paste(SoilVar,Tvar,"_Kfold_",q,"_Training.csv",sep=""),sep=",", col.names=T,row.names=F)
      write.table(r.ulPI_mat,file=paste(SoilVar,Tvar,"_Kfold_",q,"_UpperLowerLimits.csv",sep=""),sep=",", col.names=T,row.names=F)
      write.table(cubrule_mat,file=paste(SoilVar,Tvar,"_Kfold_",q,"_CubistRules.csv",sep=""),sep=",", col.names=T,row.names=F)
      gc()
    }
  }
  
  #Write out final Soil Diagnostics and obtain mean k-fold stats
  setwd(OutputDirectorySoilVar) #Place in Data file directory
  diagnostic_mat1<-data.frame(diagnostic_mat[order(diagnostic_mat[,2], diagnostic_mat[,1]),]) 
  diagnostic_mat2<- split(diagnostic_mat1,diagnostic_mat1[,2])
  diagnam<-names(diagnostic_mat2)
  for(p in diagnam){ #For each variable in the data file
    #Subset diagnostic matrix relating to diagnam
    diagvar<-paste("diagnostic_mat2$",p,sep="")
    diagmat<-eval(parse(text=diagvar))
    #Derive individual stats from the variable (diagnam) in question 
    diagmat1<-cbind(diagmat$Kfold_level,diagmat[,3:ncol(diagmat)])
    diagmat2<-matrix(NA,nrow=nrow(diagmat1),ncol=ncol(diagmat1))
    for(s in 1:ncol(diagmat1)){
      for(t in 1:nrow(diagmat1)){
        diagmat2[t,s]<-as.numeric(as.character(diagmat1[t,s]))
      }
    }
    #Derive K-fold means
    diagmat2<-rbind(diagmat2,NA)
    diagmat2[(nrow(diagmat2)),2:ncol(diagmat2)]<-colMeans(diagmat2[1:(nrow(diagmat2)-1),2:ncol(diagmat2)]) 
    diagmat2<-diagmat2[order(diagmat2[,1]),]
    #Derive K-fold median
    diagmat2<-rbind(diagmat2,NA)
    diagmat2[(nrow(diagmat2)),2:ncol(diagmat2)]<-colMedians(diagmat2[1:(nrow(diagmat2)-2),2:ncol(diagmat2)]) 
    diagmat2<-diagmat2[order(diagmat2[,1]),]
    #Derive K-fold standard deviations
    diagmat2<-rbind(diagmat2,NA)
    for(o in 2:ncol(diagmat2)){
      diagmat2[nrow(diagmat2),o]<-sd(diagmat2[1:(nrow(diagmat2)-3),o]) #Derive K-fold standard dviations
    }
    #Write final diagnostics table
    diagmat2<-cbind(c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","MEAN","MEDIAN","STDEV"),diagmat2[,2:ncol(diagmat2)])
    colnames(diagmat2)<- c(names(diagmat[1]),names(diagmat[3:length(names(diagmat))]))
    write.table(diagmat2,file=paste(SoilVar,p,"_Diagnostics.csv",sep=""),sep=",", col.names=T,row.names=F)
  }
}

#Added by PZ - Export model data - 17/5/17
JobDir <- gsub('.{12}$', '', DataDirectory)
setwd(JobDir)
write.csv(datacov, file = "modeldata.csv", row.names = FALSE) 

#End of Module 1


