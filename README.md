# EZTr

The pipeline automates conversion of gravimetric sensors or load cells (LC) time series collected at HTP platforms into 
Evapotranspiration (ETr) and Transpiration (Tr) time series.

The EZTr_main function iterates over the different time points (at 15 min interval) of
each sector or pot to generate a matrix of Tr values. Fifteen biologically
relevant features are extracted from the Tr time series for each day,
followed by computation of feature heritability.

The entire process comprises four steps: 
Step1 - Conversion of LC data to ETr.
Step2 - Calculation of reference i.e. Penman Monteith ET for the given weather
condition, and filtering ETr based on reference ET.
Step3 - Generating smooth ETr time series. 
Step4 - Extraction of Tr from smooth ETr, Tr features and each feature's broad-sense heritability estimate.

Functions used for each step:
Step1 - extractRawLCmatrix(), curateRawLCgetET(), filterETrExtremeCols()
Step2 - extractWthrVar(), prepcsWthr(), calculateETref(),generateETrRatio(), genThreshVal(), dataPART(), threshETr(), ordFiltETr()
Step3 - smoothETr()
Step4 - calculateTr(), getFeatures(), getFeatureHe()
