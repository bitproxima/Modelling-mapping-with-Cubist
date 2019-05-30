/*
 Step 1.
 
 These sql queries are Step 1. in a multi-stepped process to extract SALI lab data used in Cubist modelling of soil attributes (SAMFTAGA_RT_Kfolds_windows_Modlue1.R) for the project 'Mapping soil erodibility in the Fitzroy'
 
 Summary of all steps
  Scripts are located in \\lands\data\DSITI\LandSciences\NAS\slr_soils\Projects\Fitzroy_ErodSoil\Modelling\SiteData\scripts
  For each new dat extraction run
  1. Run Inserts.sql - Outputs Anew_depths.csv & Bnew_depths.csv
  For each soil attribute (Clay, Silt, FS, CS, ESP, Salinity, Cat_Ca, Cat_Mg, Clay_Act)
  2. Run SiteData.sql - Output labdata.csv
  3. Run DataPrepSpline.R - Output soilattribute_with_inserts.csv 
  4. Run SplineTool.exe - Output soilattribute_stdout.txt
  5. Run DataPrepCubist.R - Output soilattribute.csv
 
 Summary of what this sql query does
 There are two separate sql queries in this script, one for A horizon and one for B horizon
 The queries create a list of fictious sample depths at the change between A and B horizons in duplex soils to influenece the ASRIS Spline v2.0 tool
 Duplex soils are identified by their assigned soil classification in either ASC, PPF or GSG (see below for included classifications). Sites without a soil classification are not considered
 Buried horizons are not considered and results are restricted to the Fitzroy modelling area
 Save results from each query in \\lands\data\DSITI\LandSciences\NAS\slr_soils\Projects\Fitzroy_ErodSoil\Modelling\SiteData\Harmonised_data\Anew_depths.csv or as Bnew_depths.csv
  
*/

--Create one fictious sample immediately above the A/B change in the A horizon
select
  h.project_code,
  h.site_ID,
  MAX(h.lower_depth-0.01)*100 UD,
  MAX(h.lower_depth) *100 LD
  
from
  SIT_OBSERVATIONS o,
  SIT_LOCATIONS c,
  SIT_SOIL_CLASSIFICATIONS sc,
  SIT_HORIZONS h
  
where
  o.project_code = h.project_code
  and o.site_id = h.site_id
  and o.obs_no = h.obs_no
  and o.project_code = c.project_code
  and o.site_id = c.site_id
  and o.obs_no = c.obs_no
  and o.project_code = sc.project_code
  and o.site_id = sc.site_id
  and o.obs_no = sc.obs_no
  and DESIGN_MASTER like 'A'
  and h.DESIGN_NUM_PREFIX is null
  and c.latitude between -26.476401 and -21.205394 and c.LONGITUDE between 146.557398 and 151.586142 and c.DATUM = 3 --Fitzroy modelling area
  and (sc.asc_ord in ('CH', 'KU', 'SO') or sc.PPF like 'D%' or sc.exhibit_GSG_CODE in ('SDS', 'BP', 'LP', 'SH', 'YP', 'GP', 'GBP', 'SB', 'SC', 'RP', 'SK', 'SZ' ))--Duplex profiles only
    
GROUP BY
  h.PROJECT_CODE,
  h.SITE_ID  
  
ORDER BY
  h.PROJECT_CODE,
  h.SITE_ID
;

--Create one fictious sample immediately below the A/B change in the B horizon
select
  h.project_code,
  h.site_ID,
  MIN(h.UPPER_DEPTH) *100 UD,
  MIN(h.UPPER_DEPTH+0.01)*100 LD
  
from
  SIT_OBSERVATIONS o,
  SIT_LOCATIONS c,
  SIT_SOIL_CLASSIFICATIONS sc,
  SIT_HORIZONS h
  
where
  o.project_code = h.project_code
  and o.site_id = h.site_id
  and o.obs_no = h.obs_no
  and o.project_code = c.project_code
  and o.site_id = c.site_id
  and o.obs_no = c.obs_no
  and o.project_code = sc.project_code
  and o.site_id = sc.site_id
  and o.obs_no = sc.obs_no
  and DESIGN_MASTER like 'B'
  and h.DESIGN_NUM_PREFIX is null
  and c.latitude between -26.476401 and -21.205394 and c.LONGITUDE between 146.557398 and 151.586142 and c.DATUM = 3 --Fitzroy modelling area
  and (sc.asc_ord in ('CH', 'KU', 'SO') or sc.PPF like 'D%' or sc.exhibit_GSG_CODE in ('SDS', 'BP', 'LP', 'SH', 'YP', 'GP', 'GBP', 'SB', 'SC', 'RP', 'SK', 'SZ' ))--Duplex profiles only
    
GROUP BY
  h.PROJECT_CODE,
  h.SITE_ID  
  
ORDER BY
  h.PROJECT_CODE,
  h.SITE_ID
;