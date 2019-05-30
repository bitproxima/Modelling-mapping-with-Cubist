# Modelling-mapping-with-Cubist
This is the SLR implementation of cubist 

## Summary
This workflow was orginally developed to digitally map vulnerability to soil erosion in the Fitzroy NRM region. It includes a variety of scripts written in SQL, R and python. The main R script was provided through TERN by Kidd and Webb in 2014. Vulnerability to soil erosion is mapped according to a decision tree framework developed for a similar project in the Burdekin Dry Tropics region. It is based on soil attributes that were mapped individually using a method based on a cubist data mining algorthium. The soil attribute maps use soil data from SALI at specfic points in the region. The workflow and tools provide in this repository can be and have been use in other regions of Queensland. For example the Upper Brisbane Valley and Logan-Albert catchments.

To assist with the task, a workflow has been defined, necessary detail provided to complete each step including background info and scripts to automate some of the tasks.

## Workflow
1.  Extract soil attribute data from SALI (Clay, ESP, Salinity, Cat_Ca, Cat_Mg, Clay_Act).
1.  Prepare soil attribute data for the spline by adjusting it if soil profile has a texture-contrast.
1.  Standardise soil data to global soil map standard depths using a spline.
1.  Format standardised data for modelling in cubist.
1.  Fit a cubist model to the data, perform k-fold cross validation and calculate the upper/lower prediction limits to the specfied confidence intervals.
1.  Map the cubist rules to produce the lower, upper and mediam predictions of each soil attribute.
1.  Summaries model statistics and determine uncertainity of final maps.
1.  Produce vulnerability to soil erosion maps according the framework.

**Installation:** To download all scripts and files in this repository, click green **Clone or Download** button, select **Download ZIP**, save zip file to C://Temp, extract files. Always source files from this respository so that you have the latest version. You will need SQL Developer (linked to SALI), RStudio, R and Pyhon installed.

**Usage:** This project is wriiten to assist with digital soil mapping.

**Contributing:** If you can see improvements that can be made, please do a *pull request* or raise an *issue*.

**Credits:** The main cubist modelling, validation and uncertainity calculations were scripted by Darren Kidd, Mathew Webb and Brenden Malone (2014) for TERN. Keith Moodie converted the vulnerability to soil erosion framework into a python script. Peter Zund conceived the framework and wrote the pre and post processing scripts to implement the whole process and is the author of this repository.

**License:** Free to use with acknowledgement.

### Further detail
1.  
