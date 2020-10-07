# EZTr

The pipeline automates conversion of gravimetric sensors or load cells (LC) time series collected at HTP platforms into 
Evapotranspiration (ETr) and Transpiration (Tr) time series.

The EZTr_main function iterates over the different time points (at 15 min interval) of
each sector or pot to generate a matrix of Tr values. Fifteen biologically
relevant features are extracted from the Tr time series for each day,
followed by computation of feature heritability.

The entire process comprises four steps: 
Step1 - Conversion of LC data to ETr.
Step2 - Calculation of reference i.e. Penman Monteith ET for the given weather condition, and filtering ETr based on reference ET.
Step3 - Generating smooth ETr time series. 
Step4 - Extraction of Tr from smooth ETr, Tr features and each feature's broad-sense heritability estimate.

Functions used for each step:
Step1 - extractRawLCmatrix(), curateRawLCgetET(), filterETrExtremeCols().
Step2 - extractWthrVar(), prepcsWthr(), calculateETref(),generateETrRatio(), genThreshVal(), dataPART(), threshETr(), ordFiltETr().
Step3 - smoothETr().
Step4 - calculateTr(), getFeatures(), getFeatureHe().

# Input Files required:
1. Loadcells
2. Loadcells metadata
3. Sensor climate data
4. Sensor unit map
5. Plant eye 
(Merge all the above data files into R data by executing the script createExpData.R to generate corresponding RData)

# Results folder
In the results folder, create 2 sub-folders: rawFeaturesTimeSeries and smthFeaturesTimeSeries 

# EZTr_main.R execution
Run the pipeline script EZTr_main.R using the following experiment-specific inputs:
L 35. load("./data/HTP_data_Exp41_PM.RData"); 
L 36. allData <- HTP_data_Exp41_PM; 
L 38. LastDate="2020-01-28"; 
L 39. irrg.dts <- c("2019-12-28", "2019-12-30", "2020-01-02", "2020-01-04", "2020-01-06",
"2020-01-08", "2020-01-10", "2020-01-11", "2020-01-13", "2020-01-14",
"2020-01-16", "2020-01-17", "2020-01-18", "2020-01-20", "2020-01-22",
"2020-01-23", "2020-01-25", "2020-01-27"); 
L 43. Date1="2019-12-27 23:46:00"; 
L 44. Date2="2020-01-28 23:45:00"; 
L 45. opPATH <- "./results/Exp41-PM-NOirrg/"; 
L 46. opPATH.smth="D:/EZTr-master/results/Exp41-PM-NOirrg/smthFeaturesTimeSeries/"; 
L 47. opPATH.raw="D:/EZTr-master/results/Exp41-PM-NOirrg/rawFeaturesTimeSeries/";
In function calculateTr: L 38. LAI.mat[i, ] <- ((((sec.lai.tmp/100))/0.36)*(1/0.26)/10000)
