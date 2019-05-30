#########################################################################
# Code for Soil Erodability Index (for Fitzroy, only slightly modified from code developed Aug 2013 for Burdekin)
# Revisions:
# 15-Feb-2017,  Keith: adapted from Burdekin_Landsuit_20150310.py
#########################################################################

import time
import os, sys, platform
import numpy
import qgis.core
from osgeo import gdal
from osgeo.gdalconst import *
print 'Running QGis: ',qgis.core.QGis.QGIS_VERSION
print 'using Python ',platform.python_version(),"   Numpy ", numpy.version.version,"   GDAL ", gdal.VersionInfo(),'\n'

#The numpy data types most useful for numerical work are int8, int16, int32, int64, float32, float64
#for GDAL the names are gdal.GDT_Int16  GDT_Byte  GDT_Int16  GDT_Float32 ...

# folder structures are different on starship Vs Keith's Vs Peter Zunds machines just recording them all here
# root directory for input files please dont change these
KMRootDirectory = r'D:\__Fitzroy_SoilErod\20170131'     #Keith
StRootDirectory = r'\\starship\era\projects\Fitzroy_ErodSoil\Modelling\InputModelData' # starship, NB these are zipped
PZRootDirectory = r'C:\DATA\Temp'                          #Peter

# root directory for outputs, dont change these
KMOutputDirectory = KMRootDirectory + r'/outputs'  # Keith
PZOutputDirectory = r'C:\DATA\Temp'                          # Peter & Jim

# set which one you want to use here
RootDirectory = StRootDirectory
OutputDirectory = PZOutputDirectory

# *** runSteps can be used to only run the code to a particular step so that outputs from that step can be viewed
# NB if set to 5 or higher then outputs from steps 3,4,5 will be produced.
#runSteps = 1   # run to step 1 and stop, out tiffs will be from step 1             outputs named like out1A...
#runSteps = 2   # run step 1 and 2 then stop, output tiffs will be from step 2,     outputs named like out2A...
#runSteps = 3   # run step 1..3 then stop, output tiffs will be from step 3,        outputs named like out3A...
#runSteps = 4   # run step 1..4 then stop, output tiffs will be from steps 3,4      outputs out3A..., out4A...
#runSteps = 5   # run step 1..5 then stop, output tiffs will be from steps 3,4,5    outputs out3A*, out4A*, outOverall

# to get all outputs run with runSteps= 1, then with runSteps = 2, then with runSteps = 99
runSteps = 6

print "Running steps up to step ", runSteps,'\n'
print 'Reading from: ' + RootDirectory + ' outputs going to: ' + OutputDirectory

# tileSize constant
tileSize = 64

##### Set file names ####
# Now set up file names for all the input and output files using above 2 'paths'
# inputs. name indicates depth of upper boundary of the layer. e.g. Clay15F is layer from 0.15 to 0.30m

CaMg0F         = RootDirectory + r'/Ca_Mg_0to5cm.tif'                  # 0.00 to 0.05 m
CaMg5F         = RootDirectory + r'/Ca_Mg_5to15cm.tif'
CaMg15F        = RootDirectory + r'/Ca_Mg_15to30cm.tif'
CaMg30F        = RootDirectory + r'/Ca_Mg_30to60cm.tif'
CaMg60F        = RootDirectory + r'/Ca_Mg_60to100cm.tif'               #0.60m to 1.0m
CaMg100F       = RootDirectory + r'/Ca_Mg_100to200cm.tif'

Clay0F         = RootDirectory + r'/Clay_x0to5cm_predicted_mosaic.tif'
Clay5F         = RootDirectory + r'/Clay_x5to15cm_predicted_mosaic.tif'
Clay15F        = RootDirectory + r'/Clay_x15to30cm_predicted_mosaic.tif'
Clay30F        = RootDirectory + r'/Clay_x30to60cm_predicted_mosaic.tif'
Clay60F        = RootDirectory + r'/Clay_x60to100cm_predicted_mosaic.tif'
Clay100F       = RootDirectory + r'/Clay_x100to200cm_predicted_mosaic.tif'

ClayAct30F     = RootDirectory + r'/ClayActivity_30_60cm.tif'

Depth2BF       = RootDirectory + r'/D2B_DEPTH_predicted_mosaic.tif'
Depth2RockF    = RootDirectory + r'/D2R_DEPTH_predicted_mosaic.tif'

EC0F          = RootDirectory + r'/EC_x0to5cm_predicted_mosaic.tif'
EC5F          = RootDirectory + r'/EC_x5to15cm_predicted_mosaic.tif'
EC15F         = RootDirectory + r'/EC_x15to30cm_predicted_mosaic.tif'
EC30F         = RootDirectory + r'/EC_x30to60cm_predicted_mosaic.tif'
EC60F         = RootDirectory + r'/EC_x60to100cm_predicted_mosaic.tif'
EC100F        = RootDirectory + r'/EC_x100to200cm_predicted_mosaic.tif'

ESP0F         = RootDirectory + r'/ESP_x0to5cm_predicted_mosaic.tif'
ESP5F         = RootDirectory + r'/ESP_x5to15cm_predicted_mosaic.tif'
ESP15F        = RootDirectory + r'/ESP_x15to30cm_predicted_mosaic.tif'
ESP30F        = RootDirectory + r'/ESP_x30to60cm_predicted_mosaic.tif'
ESP60F        = RootDirectory + r'/ESP_x60to100cm_predicted_mosaic.tif'
ESP100F       = RootDirectory + r'/ESP_x100to200cm_predicted_mosaic.tif'


# Outputs
if ( runSteps <= 1):
    N = '1'  # N used in output filenames
elif (runSteps == 2):
    N = '2'
else:
    N = '3'    
outA_0         = OutputDirectory + r'/out' +N+ r'A_0.tif'
outA_5         = OutputDirectory + r'/out' +N+ r'A_5.tif'
outA_15        = OutputDirectory + r'/out' +N+ r'A_15.tif'
outA_30        = OutputDirectory + r'/out' +N+ r'A_30.tif'
outA_60        = OutputDirectory + r'/out' +N+ r'A_60.tif'
outA_100       = OutputDirectory + r'/out' +N+ r'A_100.tif'

outB_0         = OutputDirectory + r'/out' +N+ r'B_0.tif'
outB_5         = OutputDirectory + r'/out' +N+ r'B_5.tif'
outB_15        = OutputDirectory + r'/out' +N+ r'B_15.tif'
outB_30        = OutputDirectory + r'/out' +N+ r'B_30.tif'
outB_60        = OutputDirectory + r'/out' +N+ r'B_60.tif'
outB_100       = OutputDirectory + r'/out' +N+ r'B_100.tif'

outA4overallF   = OutputDirectory + r'/outA4overall.tif'
outB4overallF   = OutputDirectory + r'/outB4overall.tif'
outOverallF     = OutputDirectory + r'/outOverall.tif'
outOverall6DAFF = OutputDirectory + r'/outOverallDAFF.tif'
outOverall6FORG = OutputDirectory + r'/outOverallForage.tif'


#--------- Functions v ----------
# Some of these functions may not be in use in this code, just kept in in case they are needed in future
#
# encode NB everything between 1st and 2nd values is assigned class 1, etc. So zero never assigned.
# used by codeRaster  iff lvlType=1
def encode( datArray, lvls):
    """
    returns a new array of same dimensions as datArray with values encoded according to a passed array of levels,
    everything between the 1st and 2nd values in 'lvls' is assigned 1, everything between the 2nd and 3rd values in lvls is assigned 2, etc NB zero is never assigned.
    """
    code = numpy.zeros(datArray.shape, dtype=numpy.int8)
    for ii in range(1, len(lvls)):
       code[ (datArray >= lvls[ii-1] ) & (datArray < lvls[ii]) ] = ii
    return code

def encodeRangeInt( datArray, lvls):
    """
    returns a new int8 array of same dimensions as datArray with values encoded using 'lvls' (string),
    lvls consists of triplets e.g. "3, 1, 3; 5, 6, 8"  NB each triplet = valueToAssign, lowerLimit, upperLimit
    Note: lowerLimit and upperLimit expected to both be integers
    """
    code = numpy.zeros(datArray.shape, dtype=numpy.int8)
    tmp = lvls.split(';')
    for lvl in tmp:
        xx = lvl.split(',')
        v = int(xx[0])
        a = int(xx[1])
        b = int(xx[2])
        if (a > b):
            t = a; a = b; b = t
        code[ (datArray >= a ) & (datArray <= b) ] = v
    return code

#encodeRangeFlt not in use atm
def encodeRangeFlt( datArray, lvls):
    """
    returns a new int8 array of same dimensions as datArray with values encoded using 'lvls' (string),
    lvls consists of triplets e.g. "3, 1.5, 3.3; 5, 5.6, 7.7"  NB each triplet = valueToAssign, lowerLimit, upperLimit
    """
    code = numpy.zeros(datArray.shape, dtype=numpy.int8)
    tmp = lvls.split(';')
    for lvl in tmp:
        xx = lvl.split(',')
        v = int(xx[0])
        a = float(xx[1])
        b = float(xx[2])
        if (a > b):
            t = a; a = b; b = t
        code[ (datArray >= a ) & (datArray <= b) ] = v
    return code

#encodeFromOne not in use atm. 
def encodeFromOne( datArray, lvls):
    code = numpy.zeros(datArray.shape, dtype=numpy.int8)
    tmp = lvls.split(';')
    for lvl in tmp:
        xx = lvl.split(',')
        v = int(xx[0])
        a = float(xx[1])
        code[ (datArray == a ) ] = v
    return code

#encodeFromTwo not in use atm
def encodeFromTwo( datArray1, datArray2, lvls):
    """
    returns a new array of same dimensions as datArray with values encoded using array lvls
    lvls consists of triplets e.g. "3, 1.5, 3.3; 5, 5.6, 7.7"  NB each triplet = valueToAssign, lowerLimit, upperLimit
    """
    code = numpy.zeros(datArray1.shape, dtype=numpy.int8)
    tmp = lvls.split(';')
    for lvl in tmp:
        xx = lvl.split(',')
        v = int(xx[0])
        a = float(xx[1])
        b = float(xx[2])
        code[ (datArray1 == a ) & (datArray2 == b) ] = v
    return code

def checkSame( dataSet, masterDS ):
    """
    checks transform and dimensions of datasets are compatable, used by openRead
    """
    failMsg = ''
    if dataSet.RasterXSize != masterDS.RasterXSize:
        failMsg = failMsg + 'Columns '
    if dataSet.RasterYSize != masterDS.RasterYSize:
        failMsg = failMsg + 'Rows '
    geo1 = masterDS.GetGeoTransform()
    geo2 = dataSet.GetGeoTransform()
    if abs(geo1[0] - geo2[0]) > 0.00001:
        failMsg = failMsg + 'X Origin '
    if abs(geo1[3] - geo2[3]) > 0.00001:
        failMsg = failMsg + 'Y Origin '
    if abs(geo1[1] - geo2[1]) > 0.00001:
        failMsg = failMsg + 'pixel width '
    if abs(geo1[5] - geo2[5]) > 0.00001:
        failMsg = failMsg + 'pixel height '
    if abs(geo1[2] - geo2[2]) > 0.00001:
        failMsg = failMsg + 'rotation '
    if abs(geo1[4] - geo2[4]) > 0.00001:
        failMsg = failMsg + 'rotation '
    if failMsg != '':
        raise Exception( 'Datasets had critical differences ' + failMsg )
    return True

def openRead( fileName, master ):
    """
    opens fileName for ReadOnly. returns handle to opened file
    errors if geotransform and size of raster are different to 'master' dataset (NB skips this is master is None)
    """
    h = gdal.Open( fileName, GA_ReadOnly)
    if h is None:
        print 'Could not open ' + fileName
        return None
    
    #Next 3 lines usful when testing
    #b = h.GetRasterBand(1)   #1st band 1 indexed            
    #print "Opened: ", fileName, "(", h.RasterYSize, h.RasterXSize, h.RasterCount, gdal.GetDataTypeName(b.DataType),")"
    #b = None
    
    if (master is None):
        return h
    
    OK = checkSame( h, master )
    if OK == False:
        raise Exception( fileName + " dimensions or location are different to master")
    
    return h

def codeRaster( inDataSet, inBandNum, levels, lvlType, outFile ):
    """
    Returns the open result dataset.
    Encode inBandNum of inDataSet with results as Int8 in band 1 of outFile.
    geotransform, rows, columns etc preserved.
    
    levels is an array used to encode the data as follows
    if lvlType is 1 then [-1.0, 0.1, 0.2, 0.5, 999.9] will result in values from -1.0 to 0.1 being encoded as 1, 0.2 to 0.5 being encoded as 3, etc
    if lvlType is 3 then it expects a series of triplets seperated by ';' with each triplet internally seperated with commas
    e.g. '0, -999999, 999999;  1,  101, 103;   2, 111, 113' will result in values from 101 to 103 being encodes as 1, etc. NB triplets processed from left to right
    """
    cols = inDataSet.RasterXSize
    rows = inDataSet.RasterYSize
    inpBand = inDataSet.GetRasterBand(inBandNum)   #band 1 indexed
    #inpNoData = inpBand.GetNoDataValue()
    
    #os.remove(outFile)
    dsout = driver.Create( outFile, cols, rows, 1, gdal.GDT_Byte )  
    dsout.SetProjection( inDataSet.GetProjection() )
    dsout.SetGeoTransform( inDataSet.GetGeoTransform() )
    outband = dsout.GetRasterBand(1)
    
    ray = numpy.zeros( (tileSize,cols), dtype=numpy.int8 )
    
    start = time.time()
    twoPC = rows/50
    for rcs in range(0, rows, tileSize):
        # read tileSize rows from each raster and classify
        size = min(tileSize, rows-rcs) 
        data1 = inpBand.ReadAsArray(0, rcs, cols, size)    # to read chunks of 'size' rows starting from rs
        
        if lvlType == 1:
            ray = encode(data1, levels )
        if lvlType == 3:
            ray = encodeRangeInt(data1, levels )
        
        outband.WriteArray( ray[:size,:], 0, rcs)      #( array, xoff, yoff )
        if rcs%twoPC <tileSize:
            sys.stdout.write('.')
    
    #Finished
    end = time.time()
    print 'codeRaster completed, {0:.2f}'.format(end-start)
    
    outband.FlushCache()
    outband = None
    inpBand = None
    return dsout

def code_A_layers( clayP, clayActivP, ESPP, DepthToBP, DepthToRockP, depth, outFile):
    cols = clayP.RasterXSize
    rows = clayP.RasterYSize
    clayB = clayP.GetRasterBand(1)   #band 1 indexed
    clayActivB = clayActivP.GetRasterBand(1)
    ESPB = ESPP.GetRasterBand(1)
    depth2BB = DepthToBP.GetRasterBand(1)
    depth2RockB = DepthToRockP.GetRasterBand(1)
    
    # from the passed in depth value, determine the layer number
    #if (depth >= 0.0 and depth < 0.0499):      # 0 to 0.05
    #    layerNum=1
    #elif ( depth < 0.1499 ):  # 0.05 to 0.15
    #    layerNum = 2
    #elif ( depth < 0.2999 ):  # 0.15 to 0.30
    #    layerNum = 3
    #elif ( depth < 0.5999 ):  # 0.30 to 0.60
    #    layerNum = 4
    #elif ( depth < 0.9999 ):  # 0.60 to 1.0
    #    layerNum = 5
    #elif ( depth < 1.9999 ):  # 1.0 to 2.0
    #    layerNum = 6
    #else:
    #    layerNum = 9
    
    #os.remove(outFile)
    dsout = driver.Create( outFile, cols, rows, 1, gdal.GDT_Byte )  # gdal.GDT_Int16  GDT_Byte  GDT_Int16  GDT_Float32 ...
    dsout.SetGeoTransform( clayP.GetGeoTransform() )
    dsout.SetProjection( clayP.GetProjection() )
    outband = dsout.GetRasterBand(1)
    
    start = time.time()
    twoPC = rows/50
    for rcs in range(0, rows, tileSize):
        # read bSize rows from each raster and classify
        size = min(tileSize, rows-rcs) 
        clay = clayB.ReadAsArray(0, rcs, cols, size)    # to read chunks of 'size' rows starting from rs
        clayActiv = clayActivB.ReadAsArray(0, rcs, cols, size) 
        ESP  = ESPB.ReadAsArray(0, rcs, cols, size)
        depthToB  = depth2BB.ReadAsArray(0, rcs, cols, size)
        depthToRock  = depth2RockB.ReadAsArray(0, rcs, cols, size)
        
        # First create tile and take care of NoData outside of the catchment. Using clay as a mask
        result = numpy.zeros( (size,cols), dtype=numpy.int8 )   # everything is 0 (including areas outside mask)
        #result[ (result==0) ] = 0 # make everything -1  (includes nulls)
        result[ (clay >= 0) & (clay <= 100) ] = 99              # make non-nulls 99 (use clay to apply a mask)
        
        # Now apply the rules
        result[ ~(clay > 20.0) ] = 2                            # 2 done
        result[ ((result==99) & ~(ESP>=6.0) )] = 1              # 2,1 done
        result[ ((result==99) & ~(clayActiv>0.6) )] = 3         # 3,2,1 done
        result[ (result==99) ] = 4                              # all done
        
        if ( runSteps >=2 ):
            #Step 2: if this layer (depth) is deeper than depthToB in this cell
            result[ ((result>0) &(depth>depthToB)) ] = 8            # B layer
        
        if ( runSteps >=3):
            #Step 3: if this layer (depth) is deeper than depthToRock in this cell
            result[ ((result>0) & (depth > depthToRock )) ] = 9  # rock
        
        # Now save
        outband.WriteArray( result[:size,:], 0, rcs)      #( array, xoff, yoff )
        if rcs%twoPC <tileSize:
            sys.stdout.write('.')
    
    #Finished
    end = time.time()
    print '  A bands completed, {0:.2f}'.format(end-start)
    outband.FlushCache()
    
    outband = None    
    clayB = None
    clayActivB = None
    ESPB = None
    depth2BB = None
    depth2RockB = None
    return dsout

def code_B_layers( clayP, ESPP, ECP, CaMgP, DepthtoBP, DepthToRockP, depth, outFile ):
    cols = ESPP.RasterXSize
    rows = ESPP.RasterYSize
    #if ( clayP != None):
    clayB = clayP.GetRasterBand(1)   #band 1 indexed
    ESPB = ESPP.GetRasterBand(1)
    ECB = ECP.GetRasterBand(1)
    #if ( CaMgP != None):
    CaMgB = CaMgP.GetRasterBand(1)
    depth2BB = DepthtoBP.GetRasterBand(1)
    depth2RockB = DepthToRockP.GetRasterBand(1)
    
    
    #os.remove(outFile)
    dsout = driver.Create( outFile, cols, rows, 1, gdal.GDT_Byte )  # gdal.GDT_Int16  GDT_Byte  GDT_Int16  GDT_Float32 ...
    dsout.SetGeoTransform( ESPP.GetGeoTransform() )
    dsout.SetProjection( ESPP.GetProjection() )
    outband = dsout.GetRasterBand(1)
    
    start = time.time()
    twoPC = rows/50
    
    for rcs in range(0, rows, tileSize):
        # read bSize rows from each raster and classify
        size = min(tileSize, rows-rcs)
        #if ( clayP != None):
        clay = clayB.ReadAsArray(0, rcs, cols, size)    # to read chunks of 'size' rows starting from rs
        ESP  = ESPB.ReadAsArray(0, rcs, cols, size)
        EC = ECB.ReadAsArray(0, rcs, cols, size)
        #if ( CaMgP != None):
        CaMg = CaMgB.ReadAsArray(0, rcs, cols, size)
        depthToB  = depth2BB.ReadAsArray(0, rcs, cols, size)
        depthToRock  = depth2RockB.ReadAsArray(0, rcs, cols, size)
        
        # First create tile and take care of NoData outside of the catchment. Using clay as a mask
        result = numpy.zeros( (size,cols), dtype=numpy.int8 )   # everything is zero
        result[ (ESP >= 0) & (ESP <= 100) ] = 99                # make non-nulls 99, null value is 0
        
        # now apply the rules. NB coding style here is to keep it as similar as possible (readable/checkable) to way it is in the .doc
        #Step 1:
        #if ( clayP != None):
        result[ ~(clay > 10.0) ] = 1                        # 1 part done, rest 0
        result[ ((result==99) & ~(ESP > 15.0) )] = 1            # 1 done, rest 0
        result[ ((result==99) & ~( EC < 0.5 ) )] = 2            # 1,2 done, rest 0
        #if ( CaMgP != None):
        result[ ((result==99) & ~(CaMg<1) )] = 3            # 1,2,3 done, just 4 left
        result[  (result==99) ] = 4                         # anything still 99 will be a 4
        #else:
        #    result[ (result==99) ] = 3
        
        if ( runSteps >=2 ):
            #Step 2: exclude if part of A layers
            result[ ((result>0) & (depth<=depthToB)) ] = 8          # A layer
        if ( runSteps >=3 ):
            #Step 3: if this layer (depth) is deeper than depthToRock in this cell
            result[ ((result>0) & (depth > depthToRock )) ] = 9  # rock
            result[ ((result>0) & (depthToB > (depthToRock-0.01))) ] = 9
            #NB if this changed then reconsider effect of  if ( filterTopLayer==1 ):  in step 4 B wrt rock
        
        #Save result
        outband.WriteArray( result[:size,:], 0, rcs)      #( array, xoff, yoff )
        if rcs%twoPC <tileSize:
            sys.stdout.write('.')
        
    #Finished
    end = time.time()
    
    print '  B bands completed, {0:.2f}'.format(end-start)
    outband.FlushCache()
    outband = None; clayB = None; ESPB = None; ECB = None; depth2BB = None; depth2RockB = None
    return dsout

#codeDepth not in use atm
def codeDepth6( soilProbab, outFile):
    """
    turn the 6 bands(layers) of SoilProbability into a single band 1..6 that indicates the shallowest band that has rock
    according to the "if cell value in a layer is less than supplied threshold then layer is rock"
    """
    cols = soilProbab.RasterXSize
    rows = soilProbab.RasterYSize
    B1 = soilProbab.GetRasterBand(1) # 0-0.05m depth. band 1 indexed
    B2 = soilProbab.GetRasterBand(2)
    B3 = soilProbab.GetRasterBand(3)
    B4 = soilProbab.GetRasterBand(4)
    B5 = soilProbab.GetRasterBand(5)
    B6 = soilProbab.GetRasterBand(6) # 1-2m depth
    
    #os.remove(outFile)
    dsout = driver.Create( outFile, cols, rows, 1, gdal.GDT_Byte )  # gdal.GDT_Int16  GDT_Byte  GDT_Int16  GDT_Float32 ...
    dsout.SetGeoTransform( soilProbab.GetGeoTransform() )
    dsout.SetProjection( soilProbab.GetProjection() )
    outband = dsout.GetRasterBand(1)
    
    start = time.time()
    twoPC = rows/50
    for rcs in range(0, rows, tileSize):
        # read tileSize rows from each band and classify
        size = min(tileSize, rows-rcs) 
        aa = B1.ReadAsArray(0, rcs, cols, size)    # to read chunks of 'size' rows starting from rs
        bb = B2.ReadAsArray(0, rcs, cols, size)
        cc = B3.ReadAsArray(0, rcs, cols, size)
        dd = B4.ReadAsArray(0, rcs, cols, size)
        ee = B5.ReadAsArray(0, rcs, cols, size)
        ff = B6.ReadAsArray(0, rcs, cols, size)
        
        # First create tile and take care of NoData outside of the catchment. Using clay as a mask
        result = numpy.zeros( (size,cols), dtype=numpy.int8 )   # everything is zero
        result[ (result==0) ] = 9                               # make everything 9  (includes nulls)
        result[ ((aa>= -1.00) & (aa <= 2.0)) ] = 7              # call everything inside catchment a 7 for now
        
        result[ (aa<0.34) ] = 1                             # layer 1 is rock
        result[ (result==7) & (bb<0.44) ] = 2               # layer 2 is rock
        result[ (result==7) & (cc<0.66) ] = 3               # layer 3 is rock
        result[ (result==7) & (dd<0.64) ] = 4               # layer 4 is rock
        result[ (result==7) & (ee<0.67) ] = 5               # layer 5 is rock
        result[ (result==7) & (ff<0.63) ] = 6               # layer 6 is rock
        
        outband.WriteArray( result[:size,:], 0, rcs)      #( array, xoff, yoff )
        if rcs%twoPC <tileSize:
            sys.stdout.write('.')
        
    #Finished
    end = time.time()
    print 'code depths completed, {0:.2f}'.format(end-start)
    
    outband.FlushCache()
    outband = None
    #dsout = None    # close band but not dataset
    B1 = None
    B2 = None
    B3 = None
    B4 = None
    B5 = None
    B6 = None
    return dsout

def newMergeLayers( D0, D5, D15, D30, D60, D100, filterTopLayer, outFile):
    """
    for each cell in raster the result is the highest value for that cell in any of the 5 input datasets D0,D5,D15,D30,D60,D100
    """
    cols = D0.RasterXSize
    rows = D0.RasterYSize
    RB0 = D0.GetRasterBand(1) # 0-0.05m depth. band 1 indexed
    RB5 = D5.GetRasterBand(1)
    RB15 = D15.GetRasterBand(1)
    RB30 = D30.GetRasterBand(1)
    RB60 = D60.GetRasterBand(1)
    if ( D100 != None ):
        RB100 = D100.GetRasterBand(1)
    
    #os.remove(outFile)
    dsout = driver.Create( outFile, cols, rows, 1, gdal.GDT_Int16 )  # gdal.GDT_Int16  GDT_Byte  GDT_Int16  GDT_Float32 ...
    dsout.SetGeoTransform( D0.GetGeoTransform() )
    dsout.SetProjection( D0.GetProjection() )
    outband = dsout.GetRasterBand(1)
    
    start = time.time()
    twoPC = rows/50
    for rcs in range(0, rows, tileSize):
        # read bSize rows from each raster and classify
        size = min(tileSize, rows-rcs)
        #L0 refers to layer0-5cm, L5 is layer5-15cm,..., L100 is layer 100-200cm if present
        L0 = RB0.ReadAsArray(0, rcs, cols, size)    # to read chunks of 'size' rows starting from rs
        L5 = RB5.ReadAsArray(0, rcs, cols, size)
        L15 = RB15.ReadAsArray(0, rcs, cols, size)
        L30 = RB30.ReadAsArray(0, rcs, cols, size)
        L60 = RB60.ReadAsArray(0, rcs, cols, size)
        if ( D100 != None ):
            L100 = RB100.ReadAsArray(0, rcs, cols, size)
        
        # filter layers just read to exclude values that mean rock etc
        # have to do this or these become maximums
        if ( filterTopLayer==1 ):            #Dont always want to apply to surface layer
            L0[ (L0>6) ] = 0 
        if ( filterTopLayer==2 ):
            L0[ (L0==8) ] = 0
        L5[ (L5>6) ] = 0
        L15[ (L15>6) ] = 0
        L30[ (L30>6) ] = 0
        L60[ (L60>6) ] = 0
        
        t2 = numpy.maximum(L0, L5)
        t3 = numpy.maximum(t2, L15)
        t4 = numpy.maximum(t3, L30)
        t5 = numpy.maximum(t4, L60)
        if ( D100 != None ):
            L100[ (L100>6) ] = 0
            t6 = numpy.maximum(t5, L100)
        else:
            t6 = t5
        
        outband.WriteArray( t5[:size,:], 0, rcs)      #( array, xoff, yoff )
        if rcs%twoPC <tileSize:
            sys.stdout.write('.')
        
    #Finished
    end = time.time()
    print 'merge bands completed, {0:.2f}'.format(end-start)
    
    #Close bands but not the result dataset
    outband.FlushCache()
    outband = None
    RB0   = None
    RB5   = None
    RB15  = None
    RB30  = None
    RB60  = None
    RB100 = None
    return dsout

# when dealing with the Aoverall output an area may be
# 0 - outside area or a water body
# 1-4 A value
# 8   B layer   ( if depthToB == 0)  but this never happens at present
# 9   Rock      ( if depthToB >= (depthToRock-0.001) )

#when dealing with B overall outout an area may be
# when dealing with the Aoverall output an area may be
# 0 - outside area or a water body
# 1-4 B value
# 8   A layer   ( if depthToB > 1 metre)  
# 9   Rock      ( if depthToB >= (depthToRock-0.001) )  BUT different thresholds for this in each layer!
def mergeLayers( D1, D2, D3, D4, D5, D6, DepthToRock, DepthToB, depthFlag, outFile):
    """
    for each cell in raster the result is the highest value for that cell in any of the 5 input datasets D1..D5
    """
    cols = D1.RasterXSize
    rows = D1.RasterYSize
    B1 = D1.GetRasterBand(1) # 0-0.05m depth. band 1 indexed
    B2 = D2.GetRasterBand(1)
    B3 = D3.GetRasterBand(1)
    B4 = D4.GetRasterBand(1)
    depthThreshold = 1.0      # for Layer 5 of 60 to 100 cm
    B5 = D5.GetRasterBand(1)
    if ( D6 != None ):
        depthThreshold = 2.0  # for layer 6 of 1 to 2m
        B6 = D6.GetRasterBand(1)
    depth2Rock = DepthToRock.GetRasterBand(1)
    depth2B  = DepthToB.GetRasterBand(1)
    
    #os.remove(outFile)
    dsout = driver.Create( outFile, cols, rows, 1, gdal.GDT_Int16 )  # gdal.GDT_Int16  GDT_Byte  GDT_Int16  GDT_Float32 ...
    dsout.SetGeoTransform( D1.GetGeoTransform() )
    dsout.SetProjection( D1.GetProjection() )
    outband = dsout.GetRasterBand(1)
    d2RockNoData = depth2Rock.GetNoDataValue()
    
    start = time.time()
    twoPC = rows/50
    for rcs in range(0, rows, tileSize):
        # read bSize rows from each raster and classify
        size = min(tileSize, rows-rcs)
        
        aa = B1.ReadAsArray(0, rcs, cols, size)    # to read chunks of 'size' rows starting from rs
        bb = B2.ReadAsArray(0, rcs, cols, size)
        cc = B3.ReadAsArray(0, rcs, cols, size)
        dd = B4.ReadAsArray(0, rcs, cols, size)
        ee = B5.ReadAsArray(0, rcs, cols, size)
        if ( D6 != None ):
            ff = B6.ReadAsArray(0, rcs, cols, size)
        d2rock = depth2Rock.ReadAsArray(0, rcs, cols, size)
        d2b    = depth2B.ReadAsArray(0, rcs, cols, size)
        
        aa[ (aa>6) ] = 0 
        bb[ (bb>6) ] = 0
        cc[ (cc>6) ] = 0
        dd[ (dd>6) ] = 0
        ee[ (ee>6) ] = 0
        
        t2 = numpy.maximum(aa, bb)
        t3 = numpy.maximum(t2, cc)
        t4 = numpy.maximum(t3, dd)
        t5 = numpy.maximum(t4, ee)
        if ( D6 != None ):
            ff[ (ff>6) ] = 0
            t6 = numpy.maximum(t5, ff)
        else:
            t6 = t5
        
        if ( depthFlag == "A" ):
            t6[ (d2b < 0.01)] = 8      # B layer  (no or almost no A layer at all)
            t6[ (d2rock <= 0.025)] = 9  # needs to match code in _AB_layers for A0-5 layer
        else:
            #t6[ (d2b > d2rock)] = 9
            t6[ (d2b > 0.8) ] = 8              # A layer  (A layer more than 1m thick)
            # Now test for rock in each of the 6 standard depths 'layers' 
            # NB these tests have to be done using same depth as was passed to code_A_layers() and code_B_layers() e.g. 0.225
            t6[ (t6==0) & (d2b <= 0.025) & (d2rock < 0.025) ] = 9
            t6[ (t6==0) & (d2b >= 0.025) & (d2b <= 0.10 ) & (d2rock < 0.10 ) ] = 9
            t6[ (t6==0) & (d2b >= 0.100) & (d2b <= 0.225) & (d2rock < 0.225) ] = 9
            t6[ (t6==0) & (d2b >= 0.225) & (d2b <= 0.45 ) & (d2rock < 0.45 ) ] = 9
            t6[ (t6==0) & (d2b >= 0.450) & (d2b <= 0.80 ) & (d2rock < 0.80 ) ] = 9
            
            t6[ (t6==0) & (d2b >= 0.80) & (d2b <= 1.50 ) & (d2rock < 1.50 ) ] = 9  # added this as we didn't have deepest layer for Burdekin
        #tried isnan that didnt work
        #tried this t6[ (d2rock is not None) ] = 0  but it gave 'lines'
        
        t6[ (d2rock == d2RockNoData) ] = 0
        outband.WriteArray( t6[:size,:], 0, rcs )      #( array, xoff, yoff )
        if rcs%twoPC <tileSize:
            sys.stdout.write('.')
        
    #Finished
    end = time.time()
    print 'merge bands completed, {0:.2f}'.format(end-start)
    
    #Close bands but not the result dataset
    outband.FlushCache()
    outband = None
    B1 = None
    B2 = None
    B3 = None
    B4 = None
    B5 = None
    B6 = None
    return dsout

def merge_ABdepth( AD, BD, depthToB, outFile):
    """
    merge A and B and depthToB to give final result
    """
    cols = AD.RasterXSize
    rows = AD.RasterYSize
    A = AD.GetRasterBand(1)
    B = BD.GetRasterBand(1)
    D = depthToB.GetRasterBand(1)
    
    #os.remove(outFile)
    dsout = driver.Create( outFile, cols, rows, 1, gdal.GDT_Int16 )  # gdal.GDT_Int16  GDT_Byte  GDT_Int16  GDT_Float32 ...
    dsout.SetGeoTransform( AD.GetGeoTransform() )
    dsout.SetProjection( AD.GetProjection() )
    outband = dsout.GetRasterBand(1)
    
    start = time.time()
    twoPC = rows/50
    for rcs in range(0, rows, tileSize):
        size = min(tileSize, rows-rcs)
        
        #aa is a horizon overall result for A horizon, bb is overall result for B horizon dd is depth to rock
        aa = A.ReadAsArray(0, rcs, cols, size)    # to read chunks of 'size' rows starting from rs
        bb = B.ReadAsArray(0, rcs, cols, size)
        dd = D.ReadAsArray(0, rcs, cols, size)
        
        # Make a temporary array of the depthToB as a 1..3 value
        tmp = numpy.zeros( (size,cols), dtype=numpy.int16 )      # everything is zero
        # need to keep next 3 lines in right order
        tmp[(dd>=0.0)] = 3  
        tmp[(dd>=0.1)] = 2
        tmp[(dd>=0.3)] = 1
        
        result = tmp + aa*100 + bb*10
        
        outband.WriteArray( result[:size,:], 0, rcs)      #( array, xoff, yoff )
        if rcs%twoPC <tileSize:
            sys.stdout.write('.')
        
    #Finished
    end = time.time()
    print 'final result completed, {0:.2f}'.format(end-start)
    outband.FlushCache()
    outband = None
    A = None
    B = None
    D = None
    return dsout

#--------- End of Functions ----------


#
# All rasters need to be the same shape and cell sizes so that all cells match exactly
#
print "Starting",'\n'
driver = gdal.GetDriverByName('GTiff')   # HFA for imagine
driver.Register()

#
# Stage 1: Generate Limitation Layers
#
# inputs: (from \\starship\era\projects\Burdekin_DSM\GIS\Raster\Attribute_Outputs )
#

#===========
# run A & B layer models
#===========
depthB = openRead( Depth2BF, None )
if depthB is None:
    raise Exception( 'Could not open ' + Depth2BF )
    
depthRockB = openRead( Depth2RockF, depthB )
if depthRockB is None:
    raise Exception( 'Could not open ' + Depth2RockF )
    
if ( runSteps>=1 ):
    #single Clay Activity used for all layers
    clayActiv = openRead( ClayAct30F, depthB )
    if clayActiv is None:
        raise Exception( 'Could not open ' + ClayAct30F)
    
    # layer
    print '\n\n','L1: Starting 0.0 to 0.05m bands','\n'
    clay = openRead( Clay0F, depthB )
    if clay is None:
        raise Exception( 'Could not open ' + Clay0F)
    
    esp = openRead( ESP0F, depthB )
    if esp is None:
        raise Exception( 'Could not open ' + ESP0F)
    
    ec = openRead( EC0F, depthB )
    if ec is None:
        raise Exception( 'Could not open ' + EC0F)
    
    caMg = openRead( CaMg0F, depthB )
    if caMg is None:
        raise Exception( 'Could not open ' + CaMg0F)
    
    # using inputs of clay, esp, depthToB, rock,  and the depth of top of the layer ( in this case 0.0m)
    # generate output result which is stored in file out3A_0  and also as an open dataset A0
    # using average depth  i.e. (0.0+0.05)/2
    A0 = code_A_layers( clay, clayActiv, esp, depthB, depthRockB, 0.025, outA_0)
    
    #similarly for B def code_B_layers( clayP, ESPP, ECP, CaMgP, DepthtoBP, DepthToRockP, depth, outFile ):
    B0 = code_B_layers( clay, esp, ec, caMg, depthB, depthRockB, 0.025, outB_0)    # keep this as needed by merge
    
    # layer
    print '\n\n','L2: Starting 0.05 to 0.15m bands','\n'
    clay = openRead( Clay5F, depthB )
    if clay is None:
        raise Exception( 'Could not open ' + Clay5F )
    
    esp = openRead( ESP5F, depthB )
    if esp is None:
        raise Exception( 'Could not open ' + ESP5F)
    
    ec = openRead( EC5F, depthB )
    if ec is None:
        raise Exception( 'Could not open ' + EC5F)
    
    caMg = openRead( CaMg5F, depthB )
    if caMg is None:
        raise Exception( 'Could not open ' + CaMg5F)
    
    # using average depth i.e. (0.15+0.3)/2.0        
    A5 = code_A_layers( clay, clayActiv, esp, depthB, depthRockB, 0.1, outA_5)    # NB A5 is output for 5-15 layer i.e. named according to depth of top of layer
    B5 = code_B_layers( clay, esp, ec, caMg, depthB, depthRockB,  0.1, outB_5)  
    
    # layer
    print '\n\n','L3: Starting 0.15 to 0.30m bands','\n'
    clay = openRead( Clay15F, depthB )
    if clay is None:
        raise Exception( 'Could not open ' + Clay15F)
    
    esp = openRead( ESP15F, depthB )
    if esp is None:
        raise Exception( 'Could not open ' + ESP15F)
    
    ec = openRead( EC15F, depthB )
    if ec is None:
        raise Exception( 'Could not open ' + EC15F)
    
    caMg = openRead( CaMg15F, depthB )
    if caMg is None:
        raise Exception( 'Could not open ' + CaMg15F)
    
    # using average depth i.e. (0.15+0.3)/2.0
    A15 = code_A_layers( clay, clayActiv,  esp, depthB, depthRockB, 0.225, outA_15)
    B15 = code_B_layers( clay, esp,   ec, caMg, depthB, depthRockB, 0.225, outB_15)
    
    # layer
    print '\n\n','L4: Starting 0.30 to 0.60m bands','\n'
    clay = openRead( Clay30F, depthB )
    if clay is None:
        raise Exception( 'Could not open ' + Clay30F)
    
    esp = openRead( ESP30F, depthB )
    if esp is None:
        raise Exception( 'Could not open ' + ESP30F)
    
    ec = openRead( EC30F, depthB )
    if ec is None:
        raise Exception( 'Could not open ' + EC30F)
    
    caMg = openRead( CaMg30F, depthB )
    if caMg is None:
        raise Exception( 'Could not open ' + CaMg30F)
    
    A30 = code_A_layers( clay, clayActiv,  esp, depthB, depthRockB, 0.45, outA_30)
    B30 = code_B_layers( clay, esp,   ec, caMg, depthB, depthRockB, 0.45, outB_30)
    
    # layer
    print '\n\n','L5: Starting 0.60 to 1.0m bands','\n'
    clay = openRead( Clay60F, depthB )
    if clay is None:
        raise Exception( 'Could not open ' + Clay60F)
    
    esp = openRead( ESP60F, depthB )
    if esp is None:
        raise Exception( 'Could not open ' + ESP60F)
    
    ec = openRead( EC60F, depthB )
    if ec is None:
        raise Exception( 'Could not open ' + EC60F)
    
    caMg = openRead( CaMg60F, depthB )
    if caMg is None:
        raise Exception( 'Could not open ' + CaMg60F)
    
    A60 = code_A_layers( clay, clayActiv,  esp, depthB, depthRockB, 0.80, outA_60)
    B60 = code_B_layers( clay, esp,   ec, caMg, depthB, depthRockB, 0.80, outB_60)
    
    # layer
    print '\n\n','L6: Starting 1.0 to 2.0m bands','\n'
    # Note: code in mergeLayers assumes that L5 is 60-100 cm needs to be changed there s well if you change that
    # it also assumes that if a L6 is passed to it that that is 100-200cm,  see depthThreshold variable in mergeLayers( )
    clay = openRead( Clay100F, depthB )
    if clay is None:
        raise Exception( 'Could not open ' + Clay100F)
    
    esp = openRead( ESP100F, depthB )
    if esp is None:
        raise Exception( 'Could not open ' + ESP100F)
    
    ec = openRead( EC100F, depthB )
    if ec is None:
        raise Exception( 'Could not open ' + EC100F)
    
    caMg = openRead( CaMg100F, depthB )
    if caMg is None:
        raise Exception( 'Could not open ' + CaMg100F)
    
    A100 = code_A_layers( clay, clayActiv,  esp, depthB, depthRockB, 1.50, outA_100)
    B100 = code_B_layers( clay, esp,   ec, caMg, depthB, depthRockB, 1.50, outB_100)

clay = None
clayActiv = None
esp = None
ec = None
caMg = None

# Now step 4 merge the 5 different layers for A horizon into 1 output, similarly for B
if ( runSteps>=4 ):
    print '\n','merging all A bands and all B bands','\n'

    #mergeLayers( D1, D2, D3, D4, D5, D6, DepthToRock, DepthToB, depthFlag, outFile
    Aoverall = mergeLayers( A0, A5, A15, A30, A60, A100, depthRockB, depthB, "A", outA4overallF)
    Boverall = mergeLayers( B0, B5, B15, B30, B60, B100, depthRockB, depthB, "B", outB4overallF)
    
A0 = None; B0 = None
A5 = None; B5 = None
A15 = None; B15 = None
A30 = None; B30 = None
A60 = None; B60 = None
A100 = None; B100 = None


if ( runSteps>=5 ):
    print '\n','merging A and B to give final'
    final = merge_ABdepth( Aoverall, Boverall, depthB, outOverallF)

Aoverall = None; Boverall = None; depthB = None
final = None # close as no further procesing at this time

# Step 6
if runSteps >= 6:
    outOverall = openRead( outOverallF, None )
    if outOverall is None:
        raise Exception('Could not open ' + outOverallF)

    # simpler output with 17 categories
    levels  = '0, -999999, 999999;  1,  191, 193;   2, 111, 113;  3, 121, 123;  4, 211, 213;  5,  291, 293;   6, 131, 133;  7, 221, 223;  8, 491, 493;  9,  141, 143;  10, 231, 233;  11, 311, 313;  11, 321, 323;  11, 391, 393;  12, 241, 243;  13, 331, 333;  14, 411, 413;  14, 421, 423;  15, 341, 343;  16, 431, 433; 17, 441, 443'
    #levels = '0, -999999, 999999;  1,  101, 103;   2, 111, 113;  3, 121, 123;  4, 211, 213;  5,  201, 203;   6, 131, 133;  7, 221, 223;  8, 401, 403;  9,  141, 143;  10, 231, 233;  11, 301, 303;  11, 311, 313;  11, 321, 323;  12, 241, 243;  13, 331, 333;  14, 411, 413;  14, 421, 423;  15, 341, 343;  16, 431, 433; 17, 441, 443; 18, 901, 999'    
    dsout = codeRaster( outOverall, 1, levels, 3, outOverall6DAFF )  # currently this expects outOverall to be Int
    dsout = None;

    # and even fewer again
    levels  = '0, -999999, 999999;  21,  191, 193;   21, 111, 113;  21, 121, 123;  22, 211, 213;  22,  291, 293;   22, 131, 133;  22, 221, 223;  23, 491, 493;  23,  141, 143;  23, 231, 233;  23, 311, 313;  23, 321, 323;  23, 391, 393;  24, 241, 243;  24, 331, 333;  24, 411, 413;  24, 421, 423;  25, 341, 343;  25, 431, 433; 25, 441, 443'
    #levels = '0, -999999, 999999;  21,  101, 103;   21, 111, 113;  21, 121, 123;  22, 211, 213;  22,  201, 203;   22, 131, 133;  22, 221, 223;  23, 401, 403;  23,  141, 143;  23, 231, 233;  23, 301, 303;  23, 311, 313;  23, 321, 323;  24, 241, 243;  24, 331, 333;  24, 411, 413;  24, 421, 423;  25, 341, 343;  25, 431, 433; 25, 441, 443; 26, 901, 999'    
    dsout = codeRaster( outOverall, 1, levels, 3, outOverall6FORG )  # currently this expects outOverall to be Int
dsout = None;
outOverall = None;

print '\nFinished :)'
#End
