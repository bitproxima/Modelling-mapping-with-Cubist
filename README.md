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
1.  Extracts lab data from SALI with sample id, sample depths and site location (Datum 3). Calulate cation ratios (observing minimum thereholds for individual calulations). For sites with results from more than one method, query selects most appropriate method as per cation SSA guidelines (2014) and other methods. Save results as \\...\Modelling\SiteData\Harmonised_data\labdata.csv

1.  There are two separate sql queries in this script, one for A horizon and one for B horizon. The queries create a list of fictious sample depths at the change between A and B horizons in duplex soils to influenece the ASRIS Spline v2.0 tool. Duplex soils are identified by their assigned soil classification in either ASC, PPF or GSG (see below for included classifications). Sites without a soil classification are not considered. Buried horizons are not considered. Save results from each query in \\...\Modelling\SiteData\Harmonised_data\Anew_depths.csv or as Bnew_depths.csv

1.  
