This readme.txt file was generated on 03Jan2019 by Michael Verhoeven


-------------------
GENERAL INFORMATION
-------------------


1. COMPLETE DATA AND ANALYSIS FOR: Verhoeven MR, Larkin DJ, Newman RM. Constraining invader dominance: Effects of repeated herbicidal management and environmental factors on curlyleaf pondweed dynamics in 50 Minnesota lakes.

2. Author Information

  Principal Investigator Contact Information
        Name: Michael R Verhoeven
           Institution: University of Minnesota
           Address: 135 Skok Hall; 2003 Upper Buford Circle; St Paul MN 55108
           Email: michael.verhoeven.mrv@gmail.com
	   ORCID: https://orcid.org/0000-0002-6340-9490

  Associate or Co-investigator Contact Information
        Name: Daniel J Larkin
           Institution: University of Minnesota
           Address: 135 Skok Hall; 2003 Upper Buford Circle; St Paul MN 55108
           Email: djlarkin@umn.edu
	   ORCID: https://orcid.org/0000-0001-6378-0495

  Associate or Co-investigator Contact Information
        Name: Raymond M Newman
	   Institution: University of Minnesota
           Address: 135 Skok Hall; 2003 Upper Buford Circle; St Paul MN 55108
           Email: newma004@umn.edu
	   ORCID: https://orcid.org/0000-0002-1170-3217

3. Date of data collection: 
	Plant management and surveys from 2006 through 2015, compiled in 2017.

4. Geographic location of data collection: Minnesota, United States


5. Information about funding sources that supported the collection of the data: Funding for this research was provided through the Minnesota Aquatic Invasive Species Research Center from the Minnesota Environment and Natural Resources Trust Fund with additional support from the USDA National Institute of Food and Agriculture, Hatch grant MIN-41-081.


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 


1. Licenses/restrictions placed on the data: https://creativecommons.org/licenses/by-nc-sa/3.0/us/


2. Links to publications that cite or use the data: https://doi.org/10.1111/fwb.13468


3. Links to other publicly accessible locations of the data: none
	

4. Links/relationships to ancillary data sets: none
	

5. Was data derived from another source? no
           If yes, list source(s): n/a


6. Recommended citation for the data:
	Verhoeven MR, Larkin DJ, Newman RM. (2020) Complete data and analysis for: Constraining invader dominance: Effects of
	repeated herbicidal management and environmental factors on curlyleaf pondweed dynamics in 50 Minnesota
	lakes. Retrieved from the Data Repository for the University of Minnesota. https://doi.org/10.13020/aw92-e606



---------------------
DATA & FILE OVERVIEW
---------------------


1. File List
   A. Filename:   P_crispus_dataset_prep.R     
      Short description:   Script to combine raw data
      Description: Code used to bring environmental data, management data, and plant survey data together for analysis.      
 
   B. Filename:  h_p_crispus_analysis_manuscript.R      
      Short description:   Code used to pull in complete datasets, print summary statistics, analyze the data, and create visualizations.    
      
   C. Filename:  h_p_crispus_analysis_manuscript.html      
      Short description: Easily readable version of code used to pull in complete datasets, print summary statistics, analyze the data, and create visualizations.    

   D. Filename:  pcrispus_wDNRtrtdat.csv      
      Short description: Lake treatment/management and environmental data
      Description: Dataset containing management information, ice cover, snow cover, and Secchi clarity for lakes and years included in the analysis.

   E. Filename:  surveylktaxastat.csv      
      Short description: Plant Survey Data 
      Description: Dataset containing plant abundance data from species observed in contributed surveys that were collated for this analysis.

   F. Filename:  pcri.e_akDNRmv_combined.csv  
      Short description: Complete dataset for analysis of APRIL/MAY surveys, only includes Potamogeton crispus abundance data

   G. Filename:  pcri.p_akDNRmv_combined.csv  
      Short description: Complete dataset for analysis of JUNE surveys, only includes Potamogeton crispus abundance data

2. Relationship between files:   The raw data ("pcrispus_wDNRtrtdat.csv" & "surveylktaxastat.csv") are read into "P_crispus_dataset_prep.R", joined, cleaned, prepped, and saved (as "pcri.e_akDNRmv_combined.csv" & "pcri.p_akDNRmv_combined.csv" ). These resulting files ("pcri.e_akDNRmv_combined.csv" & "pcri.p_akDNRmv_combined.csv") are read into "h_p_crispus_analysis_manuscript.R" where statistical summaries are printed and saved, analyses are developed, tested, and applied to the data, and finally	visualizations are produced from those analyses.  All code was written in .R, analysis code was subsequently converted to .Rmd file then to the .html file in this repository. 

3. Additional related data collected that was not included in the current data package:
	Surveys submitted by data contributors in their raw formats (various) are not included here. 

4. Are there multiple versions of the dataset? NO
 
--------------------------
METHODOLOGICAL INFORMATION
--------------------------


1. Description of methods used for collection/generation of data: 
	See Methods in: Verhoeven MR, Larkin DJ, Newman RM. Constraining invader dominance: Effects of
	repeated herbicidal management and environmental factors on curlyleaf pondweed dynamics in 50 Minnesota
	lakes. Freshwater Biology. 2020; 00:1-14. https://doi.org/10.1111/fwb.13468

2. Methods for processing the data: 
	Raw point-intercept survey data were compiled into a single database. Submitted names of observed taxa were corrected based on keys submitted by the 
	data contributors, our expertise, and in consultation with accepted commonly used ID guides (e.g. Crow and Hellquist 2000, Skawinski 2014, Chadde 2012).
	From this taxonomy-corrected dataset, we calcuated various statistics for each taxon within each survey (Verhoeven, Larkin & Newman 2020). 
	These metrics are contained in "surveylktaxastat.csv". All other data development steps are described in Verhoeven, Larkin & Newman 2020. 

3. Instrument- or software-specific information needed to interpret the data:
	Program R (Full details of R session information including the packages used are available at the bottom of each output html file)

4. Standards and calibration information, if appropriate: NA
	
5. Environmental/experimental conditions: NA

6. Describe any quality-assurance procedures performed on the data:
	See "P_crispus_dataset_prep.R" for quality assurance and control steps. 

7. People involved with sample collection, processing, analysis and/or submission:
	R.M. Newman gathered existing plant, weather, and secchi data (see Acknowledgement section of Verhoeven, Larkin & Newman 2020 for data contribution details).
	M.R. Verhoeven collected existing herbicide treatment data, collated datasets, processed data to prep for analysis. M.R. Verhoeven conducted analyses under
	supervision of D.J. Larkin and R.M. Newman, and organized project files for submission to DRUM.


-----------------------------------------
DATA-SPECIFIC INFORMATION FOR: surveylktaxastat.csv
-----------------------------------------

1. Number of variables: 13

2. Number of cases/rows: 5482

3. Missing data codes:
        Where "datemv" is January 1 (e.g., 1/1/2012) survey contributed had only year completed.
	Where taxa in "taxon" have the prefix "unk." the surveyors intended taxon identity was not interpretable or was unknown. 

4. Variable List
        
    A. Name: datemv
       Description: date of survey (yyyy/mm/dd)
                    	A unique combination of datemv, datasourcemv, and Lake.ID serve as a key for unique surveys.

    B. Name: datasourcemv
       Description: abbreviated name of person/org contributing data
			Contact michael.verhoeven.mrv@gmail.com for full information and contact information for these contributors. 

    C. Name: lknamemv
       Description: lake name used in analysis
                    	These names are unique names for each waterbody included.

    D. Name: Lake.ID
       Description: Minnesota Department of Waters waterbody identifier number (DOW ID)

    E. Name: taxon
       Description: scientific name of taxon for which remaining variables were calculated
                    	Where taxa in "taxon" have the prefix "unk.", the surveyors intended taxon identity was not interpretable or was unknown.
			Where taxa in "taxon" have the suffix ".spp", the surveyors was not able to identify the observed organism to a lower taxonomic level.

    F. Name: litfoc
       Description: frequency of occurrence in sampled sites < 4.6 meter depth
                     
    G. Name: rd
       Description: mean rake density assigned to taxon at sites where it was observed in this survey
                    
    H. Name: maxdep
       Description: maximum depth at which this species was observed in this survey (meters)
                   
    I. Name: totfoc
       Description: frequency of occurrence across all sites sampled in survey
                    
    J. Name: tnsamp
       Description: total number of points sampled in survey
    
    K. Name: lnsamp
       Description: number of sites sampled at < 4.6 meter depth
                    
    L. Name: lslsamp
       Description: number of sites sampled at depths < 2*maximum observed Secchi depth for each lake
                    
    M. Name: lslfoc
       Description: frequency of occurrence in sites sampled with depths < 2*(maximum observed Secchi depth) for each lake

-----------------------------------------
DATA-SPECIFIC INFORMATION FOR: pcrispus_wDNRtrtdat.csv
-----------------------------------------

1. Number of variables:35

2. Number of cases/rows: 238

3. Missing data codes:
        NA = Data not submitted, not available, or not calculable. 
	
4. Variable List
               
    A. Name: Lake.ID
       Description: Minnesota Department of Waters waterbody identifier number (DOW ID)
                    	
    B. Name: Year
       Description: Year of record

    C. Name: Lake.x
       Description: lake name used in analysis
                    	These names are unique names for each waterbody included.

    D. Name: Trt
       Description: Indicates if lake manager reported treatment of lake in that year with herbicides targeting curlyleaf pondweed.
			T = treated
			U = untreated

    E. Name: YrsTrt
       Description: Number of years the lake has previously been treated with herbicides targeting curlyleaf pondweed.
                    	
    F. Name: PrevAUGSecchi
       Description: Previous year's August Secchi depth measurement (meters)
                     
    G. Name: IceIn
       Description: First date of full ice cover on last in previous fall (dd-MMM)
                    
    H. Name: IceOut
       Description: Last date of full ice cover in spring (dd-MMM)
                   
    I. Name: IceCover
       Description: number of days lake was covered in ice in the winter leading into datapoint year
                    
    J. Name: Snow
       Description: mean depth of snow cover across dates of IceIn and IceOut (inches)
    
    K. Name: APRMAYPcri
       Description: P. crispus frequency of occurrence at sites w/ depth < 4.6m as reported by lake manager/data contributor for April or May survey
                    
    L. Name: APRMAY3.7
       Description: P. crispus frequency of occurrence at sites w/ depth < 3.7m as reported by lake manager/data contributor for April or May survey
                    
    M. Name: APRMAYDIFF
       Description: difference between APRMAYPcri and APRMAY3.7

    N. Name: APRMAYRelDens
       Description: Mean Rake density of P. crispus at sites where it was found (excludes null occurrences) for April or May survey

    O. Name: APRMAYRelDensall
       Description: Mean Rake density of P. crispus across all sites sampled in the survey (includes Null occurrences) for April or May survey

    P. Name: JUNPcri
       Description: P. crispus frequency of occurrence at sites w/ depth < 4.6m as reported by lake manager/data contributor for June survey

    Q. Name: JUN3.7
       Description: P. crispus frequency of occurrence at sites w/ depth < 3.7m as reported by lake manager/data contributor for June survey

    R. Name: JUNDIFF
       Description: difference between JUNPcri and JUN3.7

    S. Name: JUNRelDens
       Description: Mean Rake density of P. crispus at sites where it was found (excludes null occurrences) for June survey

    T. Name: JUNRelDensall
       Description: Mean Rake density of P. crispus across all sites sampled in the survey (includes Null occurrences) for June survey

    U. Name: DIFF.4.6
       Description: difference between APRMAYPcri and JUNPcri

    V. Name: DIFF.3.7
       Description: difference between APRMAY3.7 and JUN3.7

    W. Name: DIFFRelDens
       Description: difference between APRMAYRelDens and JUNRelDens

    X. Name: County
       Description: County where lake is located

    Y. Name: Latitude.y
       Description: Latitude of lake centroid (decimal degrees)

    Z. Name: Longitude
       Description: Longitude of lake centroid (decimal degrees)

    AA. Name: Ecoregion
       Description: Minnesota Pollution Control Agency ecoregion classification (https://www.pca.state.mn.us/quick-links/eda-guide-typical-minnesota-water-quality-conditions)
			NCHF = North Central Hardwood Forests
			NGP  = Northern Glaciated Plains
			NLF  = Northern Lakes and Forests
			WCBP = Western Corn Belt Plains

    AB. Name: Lake.area..ha.
       Description: lake surface area (hectares)

    AC. Name: Littoral.area..ha.
       Description: Area of lake surface with depth < 4.6m (hectares)
    
    AD. Name: Proportion.littoral
       Description: Proportion of lake total area that is littoral =(Littoral.area..ha./Lake.area..ha.)

    AE. Name: Max.depth..m.
       Description: Maximum depth of waterbody as reported in MN Department of Natural Resources Hydrography GIS Layer (meters)

    AF. Name: Acres_Chem
       Description: Acres of lake treated with herbicides targeting P. crispus. Data from MN Department of Natural Resources permitting records.

    AG. Name: Acres_Harvested
       Description: Acres of lake harvested by commercial harvester targeting P. crispus. Data from MN Department of Natural Resources permitting records.

    AC. Name: litt.ac
       Description: Area of lake surface with depth < 4.6m (acres)

    AC. Name: perc.lit.chem
       Description: Proportion of lake littoral area (<4.6m) treated with herbicides targeting P. crispus.



-----------------------------------------
DATA-SPECIFIC INFORMATION FOR: pcri.e_akDNRmv_combined.csv
-----------------------------------------

1. Number of variables: 42

2. Number of cases/rows: 238

3. Missing data codes:
        NA = Data not submitted, not available, or not calculable. 
	
4. Variable List
               
    A. Name: Lake.ID
       Description: Minnesota Department of Waters waterbody identifier number (DOW ID)
                    	
    B. Name: Year
       Description: Year of record

    C. Name: Lake
       Description: lake name used in analysis
                    	These names are unique names for each waterbody included.

    D. Name: Latitude.x
       Description: Latitude of lake centroid (decimal degrees)

    E. Name: Trt
       Description: Indicates if lake manager reported treatment of lake in that year with herbicides targeting curlyleaf pondweed.
			1 = treated
			0 = untreated

    F. Name: YrsTrt
       Description: Number of years the lake has previously been treated with herbicides targeting curlyleaf pondweed (cumulative, reported by lake manager).
                    	
    G. Name: PrevAUGSecchi
       Description: Previous year's August Secchi depth measurement (meters)
                     
    H. Name: IceIn
       Description: First date of full ice cover on last in previous fall (DD-MMM)
                    
    I. Name: IceOut
       Description: Last date of full ice cover in spring (DD-MMM)
                   
    J. Name: IceCover
       Description: number of days lake was covered in ice in the winter leading into datapoint year
                    
    K. Name: Snow
       Description: mean depth of snow cover across dates of IceIn and IceOut (inches)
    
    L. Name: APRMAYPcri
       Description: P. crispus frequency of occurrence at sites w/ depth < 4.6m as reported by lake manager/data contributor for April or May survey
                    
    M. Name: APRMAY3.7
       Description: P. crispus frequency of occurrence at sites w/ depth < 3.7m as reported by lake manager/data contributor for April or May survey
                    
    N. Name: APRMAYDIFF
       Description: difference between APRMAYPcri and APRMAY3.7

    O. Name: APRMAYRelDens
       Description: Mean Rake density of P. crispus at sites where it was found (excludes null occurrences) for April or May survey as reported by lake manager/data contributor

    P. Name: APRMAYRelDensall
       Description: Mean Rake density of P. crispus across all sites sampled in the survey (includes Null occurrences) for April or May survey as reported by lake manager/data contributor

    Q. Name: County
       Description: County where lake is located

    R. Name: Latitude.y
       Description: Latitude of lake centroid (decimal degrees)

    S. Name: Longitude
       Description: Longitude of lake centroid (decimal degrees)

    T. Name: Ecoregion
       Description: Minnesota Pollution Control Agency ecoregion classification (https://www.pca.state.mn.us/quick-links/eda-guide-typical-minnesota-water-quality-conditions)
			NCHF = North Central Hardwood Forests
			NGP  = Northern Glaciated Plains
			NLF  = Northern Lakes and Forests
			WCBP = Western Corn Belt Plains

    U. Name: Lake.area..ha.
       Description: lake surface area (hectares)

    V. Name: Littoral.area..ha.
       Description: Area of lake surface with depth < 4.6m (hectares)
    
    W. Name: Proportion.littoral
       Description: Proportion of lake total area that is littoral =("Littoral.area..ha."/"Lake.area..ha.")

    X. Name: Max.depth..m.
       Description: Maximum depth of waterbody as reported in MN Department of Natural Resources Hydrography GIS Layer (meters)

    Y. Name: Lake.y
       Description: lake name used in analysis
                    	These names are unique names for each waterbody included.

    Z. Name: Acres_Chem
       Description: Acres of lake treated with herbicides targeting P. crispus. Data from MN Department of Natural Resources permitting records.

    AA. Name: Acres_Harvested
       Description: Acres of lake harvested by commercial harvester targeting P. crispus. Data from MN Department of Natural Resources permitting records.

    AB. Name: litt.ac
       Description: Area of lake surface with depth < 4.6m (acres)

    AC. Name: perc.lit.chem
       Description: Proportion of lake littoral area (<4.6m) treated with herbicides targeting P. crispus.

    AD. Name: Snow.sc
       Description: z-transformed (centered relative to mean and scaled relative to standard deviation) variable "Snow"

    AE. Name: PrevAUGSecchi.sc
       Description: z-transformed (centered relative to mean and scaled relative to standard deviation) variable "PrevAUGSecchi"
    
    AF. Name: IceCover.sc
       Description: z-transformed (centered relative to mean and scaled relative to standard deviation) variable "IceCover"

    AG. Name: lake_year
       Description: lake from "Lake" and year of observation from "Year." Serves as unique key for this dataset.

    AH. Name: nsamp
       Description: number of sites sampled in depths < 4.6 meters calculated from submitted point-intercept survey data

    AI. Name: nocc_pcri
       Description: number of sites sampled in depths < 4.6 meters with P. crispus present calculated from submitted point-intercept survey data

    AJ. Name: rd_pcri
       Description: mean rake density of P. crispus at survey sited where present calculated from submitted point-intercept survey data

    AK. Name: datasource
       Description: abbreviated name of person/org contributing data
			Contact michael.verhoeven.mrv@gmail.com for full information and contact information for these contributors.

    AL. Name: foc_clp
       Description: Frequency of occurrence of P. crispus across sampled sites < 4.6 m depth  = ("nocc_pcri"/"rd_pcri")

    AM. Name: cytrt
       Description: Number of consecutive years that lake was treated with herbicides targeting P. crispus prior to previous year.
			e.g., For "lakeA" in 2008, which was treated in 2006 & 2007, cytrt = 1

    AN. Name: yutrt
       Description: Years untreated. Number of years since last treatment.

    AO. Name: Trte
       Description: treated with herbicides targeting P. crispus in the previous year

    AP. Name: perc.cheme
       Description: Proportion of lake littoral area (<4.6m) treated with herbicides targeting P. crispus in the previous year

    AQ. Name:  [blank]
       Description: duplicate of field "X", below

    AR. Name: X
       Description: row identification number from write.csv function

-----------------------------------------
DATA-SPECIFIC INFORMATION FOR: pcri.p_akDNRmv_combined.csv
-----------------------------------------

1. Number of variables: 48

2. Number of cases/rows: 238

3. Missing data codes:
        NA = Data not submitted, not available, or not calculable. 
	
4. Variable List
               
    A. Name: Lake.ID
       Description: Minnesota Department of Waters waterbody identifier number (DOW ID)
                    	
    B. Name: Year
       Description: Year of record

    C. Name: Lake
       Description: lake name used in analysis
                    	These names are unique names for each waterbody included.

    D. Name: Latitude.x
       Description: Latitude of lake centroid (decimal degrees)

    E. Name: Trt
       Description: Indicates if lake manager reported treatment of lake in that year with herbicides targeting curlyleaf pondweed.
			1 = treated
			0 = untreated

    F. Name: YrsTrt
       Description: Number of years the lake has previously been treated with herbicides targeting curlyleaf pondweed (cumulative, reported by lake manager).
                    	
    G. Name: PrevAUGSecchi
       Description: Previous year's August Secchi depth measurement (meters)
                     
    H. Name: IceIn
       Description: First date of full ice cover on last in previous fall (DD-MMM)
                    
    I. Name: IceOut
       Description: Last date of full ice cover in spring (DD-MMM)
                   
    J. Name: IceCover
       Description: number of days lake was covered in ice in the winter leading into datapoint year
                    
    K. Name: Snow
       Description: mean depth of snow cover across dates of IceIn and IceOut (inches)
    
    L. Name: APRMAYPcri
       Description: P. crispus frequency of occurrence at sites w/ depth < 4.6m as reported by lake manager/data contributor for April or May survey
                    
    M. Name: APRMAY3.7
       Description: P. crispus frequency of occurrence at sites w/ depth < 3.7m as reported by lake manager/data contributor for April or May survey
                    
    N. Name: APRMAYDIFF
       Description: difference between APRMAYPcri and APRMAY3.7

    O. Name: APRMAYRelDens
       Description: Mean Rake density of P. crispus at sites where it was found (excludes null occurrences) for April or May survey as reported by lake manager/data contributor

    P. Name: APRMAYRelDensall
       Description: Mean Rake density of P. crispus across all sites sampled in the survey (includes Null occurrences) for April or May survey as reported by lake manager/data contributor

    Q. Name: JUNPcri
       Description: P. crispus frequency of occurrence at sites w/ depth < 4.6m as reported by lake manager/data contributor for June survey

    R. Name: JUN3.7
       Description: P. crispus frequency of occurrence at sites w/ depth < 3.7m as reported by lake manager/data contributor for June survey

    S. Name: JUNDIFF
       Description: difference between JUNPcri and JUN3.7

    T. Name: JUNRelDens
       Description: Mean Rake density of P. crispus at sites where it was found (excludes null occurrences) for June survey as reported by lake manager/data contributor

    U. Name: JUNRelDensall
       Description: Mean Rake density of P. crispus across all sites sampled in the survey (includes Null occurrences) for June survey as reported by lake manager/data contributor

    V. Name: DIFF.4.6
       Description: difference between APRMAYPcri and JUNPcri

    W. Name: DIFF.3.7
       Description: difference between APRMAY3.7 and JUN3.7

    X. Name: DIFFRelDens
       Description: difference between APRMAYRelDens and JUNRelDens

    Y. Name: County
       Description: County where lake is located

    Z. Name: Latitude.y
       Description: Latitude of lake centroid (decimal degrees)

    AA. Name: Longitude
       Description: Longitude of lake centroid (decimal degrees)

    AB. Name: Ecoregion
       Description: Minnesota Pollution Control Agency ecoregion classification (https://www.pca.state.mn.us/quick-links/eda-guide-typical-minnesota-water-quality-conditions)
			NCHF = North Central Hardwood Forests
			NGP  = Northern Glaciated Plains
			NLF  = Northern Lakes and Forests
			WCBP = Western Corn Belt Plains

    AC. Name: Lake.area..ha.
       Description: lake surface area (hectares)

    AD. Name: Littoral.area..ha.
       Description: Area of lake surface with depth < 4.6m (hectares)
    
    AE. Name: Proportion.littoral
       Description: Proportion of lake total area that is littoral =("Littoral.area..ha."/"Lake.area..ha.")

    AF. Name: Max.depth..m.
       Description: Maximum depth of waterbody as reported in MN Department of Natural Resources Hydrography GIS Layer (meters)

    AG. Name: Lake.y
       Description: lake name used in analysis
                    	These names are unique names for each waterbody included.

    AH. Name: Acres_Chem
       Description: Acres of lake treated with herbicides targeting P. crispus. Data from MN Department of Natural Resources permitting records.

    AI. Name: Acres_Harvested
       Description: Acres of lake harvested by commercial harvester targeting P. crispus. Data from MN Department of Natural Resources permitting records.

    AJ. Name: litt.ac
       Description: Area of lake surface with depth < 4.6m (acres)

    AK. Name: perc.lit.chem
       Description: Proportion of lake littoral area (<4.6m) treated with herbicides targeting P. crispus.

    AL. Name: Snow.sc
       Description: z-transformed (centered relative to mean and scaled relative to standard deviation) variable "Snow"

    AM. Name: PrevAUGSecchi.sc
       Description: z-transformed (centered relative to mean and scaled relative to standard deviation) variable "PrevAUGSecchi"
    
    AN. Name: IceCover.sc
       Description: z-transformed (centered relative to mean and scaled relative to standard deviation) variable "IceCover"

    AO. Name: lake_year
       Description: lake from "Lake" and year of observation from "Year." Serves as unique key for this dataset.

    AP. Name: nsamp
       Description: number of sites sampled in depths < 4.6 meters calculated from submitted JUNE point-intercept survey data

    AQ. Name: nocc_pcri
       Description: number of sites sampled in depths < 4.6 meters with P. crispus present calculated from JUNE submitted point-intercept survey data

    AR. Name: rd_pcri
       Description: mean rake density of P. crispus at survey sited where present calculated from submitted JUNE point-intercept survey data

    AS. Name: datasource
       Description: abbreviated name of person/org contributing data
			Contact michael.verhoeven.mrv@gmail.com for full information and contact information for these contributors.

    AT. Name: foc_clp
       Description: Frequency of occurrence of P. crispus across sampled sites < 4.6 m depth = "nocc_pcri"/"nsamp"

    AU. Name: cytrt
       Description: Number of consecutive years that lake was treated with herbicides targeting P. crispus prior to current year.
			e.g., For "lakeA" in 2008, which was treated in 2006 & 2007, cytrt = 2 (regardless of current treatment status)

    AV. Name: yutrt
       Description: Years untreated. Number of years since last treatment.

    AW. Name:  [blank]
       Description: duplicate of field "X", below

    AX. Name: X
       Description: row identification number from write.csv function





