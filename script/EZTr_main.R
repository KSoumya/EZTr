##################
# EZTr_main #
# Author: Soumyashree Kar

# This framework is developed as part of the research article:
# Automated Discretization of ‘Transpiration Restriction to Increasing VPD’ Features from Outdoors High-Throughput Phenotyping Data
#
# Soumyashree Kar, Ryokei Tanaka, Lijalem Balcha Korbu, Jana Kholová, Hiroyoshi Iwata, Surya S. Durbha, J. Adinarayana, Vincent Vadez
##################

#' HTP data processing
#'
#' Convert  gravimetric sensors or load cells (LC) time series into
#' Evapotranspiration (ETr) and Transpiration (Tr) time series.
#'
#' The function iterates over the different time points (at 15 min interval) of
#' each sector or pot to generate a matrix of Tr values. Fifteen biologically
#' relevant features are extracted from the Tr time series for each day,
#' followed by computation of feature heritability.
#'
#' The entire process comprises four steps: 1st - Conversion of LC data to ETr;
#' 2nd - Calculation of reference i.e. Penman Monteith ET for the given weather
#' condition, and filtering ETr based on reference ET; 3rd - Generating smooth
#' ETr time series; 4th - Extraction of Tr from smooth ETr, Tr features and
#' each feature's broad-sense heritability estimate.
#' 
#' Input Files required:
#' 1.Loadcells
#' 2.Loadcells metadata
#' 3.Sensor climate data
#' 4.Sensor unit map
#' 5.Plant eye
#'
#' #' Functions used for each step:
#' Step1 - extractRawLCmatrix(), curateRawLCgetET(), filterETrExtremeCols();
#' Step2 - extractWthrVar(), prepcsWthr(), calculateETref(),generateETrRatio(),
#'         genThreshVal(), dataPART(), threshETr(), ordFiltETr();
#' Step3 - smoothETr();
#' Step4 - calculateTr(), getFeatures(), getFeatureHe().
#'
#'
#' The function includes seven inputs, 'data' includes LC data,
#' experimental design and genotype-replicate information in metadata, weather
#' data and genotype-replicate specific leaf area data; first (Date1) and
#' last (Date2) dates of the experiment; 'irrg.dts' is a vector of dates when
#' plants were irrigated; 'opPATH' and opPATH.smth to store intermediate
#' feature-specific files.
#'
#' @param HTP_data \code{list} of 5 dataframes composed of input files
#' a) loadcell data with unit, genotype, g_alias, treatment, timestamp
#' and Mass(g) columns
#' b) metadata with unit, old_unit, Experiment, Treatment, Species
#' Genotype, G. Alias and Replicates columns
#' c) weather data with sensor, variable, timestamp, value columns
#' d) solar radiation data is stored in separate dataframe in format c)
#' e) leaf area data with Sector, Experiment, Treatment, Species,
#' Genotype, G..Alias, Replicates, timestamp and leaf area columns
#'
#' @param lastDate \code{string} containing last date of experiment
#' in YYYY-MM-DD format.
#'
#' @param irrg.dts \code{vector} of irrigation dates, each date mentioned as
#' string in YYYY-MM-DD format. If no such date is needed, same as the last date.
#'
#' @param Date1 \code{string} value specifying first date of time series in
#' YYYY-MM-DD hh:mm:ss format.
#'
#' @param Date2 \code{string} value specifying last date of time series in
#' YYYY-MM-DD hh:mm:ss format.
#'
#' @param opPATH \code{string} value specifying the directory path for
#' intermediate results to be stored.
#'
#' @param opPATH.smth \code{string} value specifying the directory path for
#' feature-specific results to be stored.
#'
#'
#' @return Return:
#'
#' The function returns a list of fifteen outcomes:
#' LCraw_tsmeta = Matrix of time stamp values present in raw load cell data,
#' LCraw = Matrix of raw load cell data,
#' ETrmm_Obs = Matrix of raw or observed ETr in mm
#' ETr_Obs_ERR.SEC.nms = Matrix of sector names with very high proportion of extreme values,
#' ETmm_Obs_FINAL = Matrix of ETr in mm after removing the erroneous sectors,
#' Wthr_agg15min = Matrix of all weather variables aggregated for 15 minutes interval,
#' Wthr.ETref.ETobs = Matrix of combined time series of Reference ET, ETr, weather values,
#' ETrRatio_TW1_ThreshVALS = Matrix of filtered values of Day-time ref ET/ETr ratios,
#' ETrRatio_TW2_ThreshVALS = Matrix of filtered values of Night-time ref ET/ETr ratios,
#' IrrgFilt_ETr = Matrix of Irrigation filtered ETr matrix,
#' IrrgFilt_ETr_Imptd = Matrix of irrigation filtered ETr matrix imputed,
#' ETr_smth = Matrix of smooth ETr time series,
#' Tr = Matrix of Tr time series,
#' featureH2 = Matrix of each feature heritability estimate on each day,
#' eachFeature_TS = List of each feature's time series for all genotypes.
#'
#' @author Soumyashree Kar email<ksoumya2301@gmail.com>
#'
#' @references
#'
#' Vadez, V., Kholová, J., Hummel, G., Zhokhavets, U., Gupta, S. K., &
#' Hash, C. T. (2015). LeasyScan: a novel concept combining 3D imaging and
#' lysimetry for high-throughput phenotyping of traits controlling plant water
#' budget. Journal of Experimental Botany, 66(18), 5581-5593.

library(easypackages)
libraries("readxl", "hms", "xts", "dplyr","mgcv","PerformanceAnalytics",   
          "signal", "tidyverse", "zoo", "h2o", "sqldf", "ggplot2", "plyr",
          "lubridate", "BioFTF", "plantecophys", "highfrequency", "stringr",
          "nonlinearTseries", "tsfeatures", "splitstackshape")

setwd("D:/EZTr_DEMO/")

# Load all the functions #
# Functions needed for processing Stage-I: LC to ETr extraction
source('./functions/extractRawLCmatrix.R')
source('./functions/curateRawLC.R')
source('./functions/filterLCExtremeCols.R')
source('./functions/getETr.R')
source('./functions/filterETrExtremeCols.R')

# Functions needed for processing Stage-II: ETref extraction and ETr thresholding
source('./functions/extractWthrVar.R')
source('./functions/prepcsWthr.R')
source('./functions/calculateETref.R')
source('./functions/generateETrRatio.R')
source('./functions/genThreshVal.R')
source('./functions/dataPART.R')
source('./functions/threshETr.R')
source('./functions/ordFiltETr.R')

# Functions needed for processing Stage-III: raw and smooth Tr, Tr-features and feature-H2 extraction
source('./functions/calculateTr.R')
source('./functions/getFeatures.R')
source('./functions/getFeatureHe.R')
source('./functions/smoothETr.R')

# Load data
load("./data/HTP_data_Exp41_PM.RData")
allData <- HTP_data_Exp41_PM

lastDate="2020-01-28"
irrg.dts <- c("2019-12-28", "2019-12-30", "2020-01-02", "2020-01-04", "2020-01-06",
                      "2020-01-08", "2020-01-10", "2020-01-11", "2020-01-13", "2020-01-14",
                      "2020-01-16", "2020-01-17", "2020-01-18", "2020-01-20", "2020-01-22",
                      "2020-01-23", "2020-01-25", "2020-01-27")
Date1="2019-12-27 23:46:00"
Date2="2020-01-28 23:45:00"
opPATH <- "./results/Exp41-PM/"
opPATH.smth="D:/EZTr_DEMO/results/Exp41-PM/smthFeaturesTimeSeries/"
opPATH.raw="D:/EZTr_DEMO/results/Exp41-PM/rawFeaturesTimeSeries/"

# Get load cells data i.e. weights of sector
m.lc <- allData$m.lc

# Get Genotype and Exp design metadata
meta.d <- allData$meta.d

# Get 3D leaf area values
pe.df <- allData$pe.df

# Get the list of species from the metadata
species.nm <- unique(meta.d$Species)

# Include the species ID e.g. 1 for 'Pearl Millet' as per the data
meta.d.sp <- meta.d[meta.d$Species==species.nm[1], ]

# Set the last date for the LC data from the exp duration
lastDate <- lastDate

# Get the dates of irrigation in Date format
irrg.dts <- as.Date(irrg.dts)

# Find sectors with missing metadata
noEntrySecs <- which(!unique(m.lc$unit) %in% unique(meta.d.sp$unit))

noEntrySecNms <- unique(m.lc$unit)[noEntrySecs]

# Remove sectors-without-metadata from original loadcells file
m.lc <- m.lc[! m.lc$unit %in% noEntrySecNms, ]


####### Stage-I: Process LC data and generate ETr matrix #######
st.time <- Sys.time()

# Extract matrix of loadcell data #
LC.MAT.OP <- extractRawLCmatrix(x = m.lc, y = meta.d.sp, z = lastDate)

LC.MAT.f <- LC.MAT.OP$LC.MAT.f

LC.MAT.TSinfo <- LC.MAT.OP$LC_tsmeta


# Add metadata to LC matrix #
# select from Metadata: "unit", "old_unit", "Genotype, "G..Alias", "Replicates"
meta.d.LCmat <- meta.d.sp[meta.d.sp$unit %in% colnames(LC.MAT.f)[-1], 
                          c(1,2,6,7,8)] 

LC.MAT.f.t <- as.data.frame(t(LC.MAT.f))

colnames(LC.MAT.f.t) <- LC.MAT.f$TS

LC.MAT.f.t <- LC.MAT.f.t[-1,]

# reorder rows of 'meta.d.sp' according to rownames/unit of LC.MAT.f.t
meta.LCDF <- meta.d.LCmat[order(match(meta.d.LCmat$unit, rownames(LC.MAT.f.t))), ]

LC.MAT.raw <- as.data.frame(cbind(meta.LCDF, LC.MAT.f.t))

write.csv(LC.MAT.raw, paste0(opPATH, "LCraw_wNA.csv"))


# Start outlier detection, removal and imputation of LC Matric to generate ETr profiles #
# Pre-process raw LC data: outliers removal and imputation #
imputed.DF.final <- curateRawLC(x = LC.MAT.f, y = meta.LCDF)

write.csv(imputed.DF.final, paste0(opPATH, "LC_olrm_imputed.csv"))


# Identify the highly extreme valued sectors #
err.sec.info <- filterLCExtremeCols(x = imputed.DF.final, y = meta.LCDF)

err.sec.nm <- err.sec.info$err.sec.NM

err.sec.meta <- err.sec.info$err.sec.META

write.csv(err.sec.meta, paste0(opPATH, "LCimp_errorUnits.csv"))


# Remove the err.cols i.e. sectors with extreme values
impData.errSEC.rmvd <- imputed.DF.final[!imputed.DF.final$unit %in% err.sec.nm, ]

write.csv(impData.errSEC.rmvd, paste0(opPATH, "LCimp_ETr_IP.csv"))


# Generate ETr profiles from "impData.errSEC.rmvd" dataframe #
et.vals <- getETr(x = impData.errSEC.rmvd)

et.obs <- et.vals$obsETr_core

ETr_Meta <- et.vals$obsETr_meta

write.csv(ETr_Meta, paste0(opPATH, "ET_Obs.csv"))


# Convert ETr in grams to mm #
ETr_Meta[ ,6:ncol(ETr_Meta)] <- et.obs*4/1000

write.csv(ETr_Meta, paste0(opPATH, "ETmm_Obs.csv"))


# Identify error plots from ETr values using the similar method as above #
ETr_err.sec.info <- filterETrExtremeCols(x = ETr_Meta, y = meta.LCDF)

err.sec.nm <- ETr_err.sec.info$ETr_err.sec.NM

err.sec.meta <- ETr_err.sec.info$ETr_err.sec.META

write.csv(err.sec.meta,  paste0(opPATH, "ET_Obs_ERR.SEC.nms.csv"))


# Remove the err.cols i.e. sectors with extreme values #
ETr_Meta_ERRsec.rmvd <- ETr_Meta

ETr_Meta_ERRsec.rmvd <- ETr_Meta_ERRsec.rmvd[!ETr_Meta_ERRsec.rmvd$unit %in% err.sec.nm, ]

write.csv(ETr_Meta_ERRsec.rmvd, paste0(opPATH, "ETmm_Obs_FINAL.csv"))
##### END ######


####### Stage-II: Process Weather data to obtain ETref and ETr ratio matrix #######
clm.df <- allData$clm.df
sensor.unit.df <- allData$sensor.unit.df

clm.df.mapped <- clm.df[clm.df$sensor %in% sensor.unit.df$sensor, ]
unq.clm.var <- unique(clm.df.mapped$variable)


# Extract weather variables individually; Assign numbers appropriately #
temperature.DF <- extractWthrVar(x = 1, y = clm.df.mapped);
temperature.DF$ts <- temperature.DF$ts + 5.5*60*60

relHUM.DF <- extractWthrVar(x = 2, y = clm.df.mapped)
relHUM.DF$ts <- relHUM.DF$ts + 5.5*60*60

windS.DF <- extractWthrVar(x = 3, y = clm.df.mapped)
windS.DF$ts <- windS.DF$ts + 5.5*60*60

solarRad.DF <- extractWthrVar(x = 5, y = clm.df.mapped)
solarRad.DF$ts <- solarRad.DF$ts + 5.5*60*60


# Read the ET obs data and get the first-last time stamps #
et.obs <- ETr_Meta_ERRsec.rmvd
et.f.ts <- colnames(et.obs[6]); 
et.l.ts <- colnames(et.obs[ncol(et.obs)])

Date1<-ymd_hms(Date1)  ## Make sure it's a Complete Cycle
Date2<-ymd_hms(Date2)

dates <- seq(Date1, Date2, by="min")


# Subset weather data #
temperature.DF.filt <- subset(temperature.DF, ts %in% dates)
relHUM.DF.filt <- subset(relHUM.DF, ts %in% dates)
solarRad.DF.filt <- subset(solarRad.DF, ts %in% dates)
windS.DF.filt <- subset(windS.DF, ts %in% dates)


# Create dataframe for Subsetted weather data #
weather.DF <- as.data.frame(matrix(nrow = length(dates), ncol = 6))

names(weather.DF) <- c("TS", "Temp", "RH", "VPD", "SR", "WS")

weather.DF[,1]<-dates

i<-nrow(weather.DF)
pbar <- create_progress_bar('text')
pbar$init(i)
for(i in 1:nrow(weather.DF))
  {
  if(i == 1)print("Weather matrix timestamp mapping status")
  ifelse(weather.DF$TS[i] %in% temperature.DF.filt$ts,
         weather.DF[i, 2] <- temperature.DF.filt$value[temperature.DF.filt$ts==weather.DF$TS[i]],
         weather.DF[i, 2] <- NA)
  
  ifelse(weather.DF$TS[i] %in% relHUM.DF.filt$ts,
         weather.DF[i, 3] <- relHUM.DF.filt$value[relHUM.DF.filt$ts==weather.DF$TS[i]],
         weather.DF[i, 3] <- NA)
  
  ifelse(weather.DF$TS[i] %in% solarRad.DF.filt$ts,
         weather.DF[i, 5] <- solarRad.DF.filt$value[solarRad.DF.filt$ts==weather.DF$TS[i]],
         weather.DF[i, 5] <- NA)
  
  ifelse(weather.DF$TS[i] %in% windS.DF.filt$ts,
         weather.DF[i, 6] <- windS.DF.filt$value[windS.DF.filt$ts==weather.DF$TS[i]],
         weather.DF[i, 6] <- NA)
  pbar$step()
}

write.csv(weather.DF, paste0(opPATH, "RAW_weather.TS_5.5.csv"))


# Preprocess each weather variable except VPD #
weather.DF[ ,2] <- prepcsWthr(x = weather.DF, y = 2) # temperature
weather.DF[ ,3] <- prepcsWthr(x = weather.DF, y = 3) # relative humidity
weather.DF[ ,5] <- prepcsWthr(x = weather.DF, y = 5) # solar radiation
weather.DF[ ,6] <- prepcsWthr(x = weather.DF, y = 6) # wind speed


# Comput VPD and insert into the weather DF #
SVP <- 610.7*(10^(7.5*weather.DF[ ,2]/(237.3+weather.DF[ ,2])))
VPD <- ((1 - (weather.DF[ ,3]/100))*SVP)/1000
weather.DF[ ,4] <- VPD

write.csv(weather.DF, paste0(opPATH, "prepcs_weather.TS_5.5.csv"))

wthr.DFxts.TS <- xts(weather.DF, order.by = as.POSIXct(weather.DF$TS, format="%Y-%m-%d %H:%M"))
wthr.DFagg15min = highfrequency::aggregatets(wthr.DFxts.TS, on="minutes",k=15, dropna=TRUE)

### Run below if TS is like 00:14:00, 00:29:00, etc. Else save directly ###
ts.x <- ymd_hms(wthr.DFagg15min$TS)
ts.x <- ts.x + hms("00:01:00")
wthr.DFagg15min$TS <- as.character(ts.x)

# Check and Filter duration of weather data as per observed ET
et.obs.TS <- colnames(et.obs)[6:ncol(et.obs)]

# IMP: Check and set the date format of et.obs.TS #
et.obs.TS <- ymd_hms(et.obs.TS)                                 
et.obs.TS.chr <- as.character(et.obs.TS)
wthr.DFagg15min <- wthr.DFagg15min[wthr.DFagg15min$TS %in% et.obs.TS.chr, ]

write.csv(wthr.DFagg15min, paste0(opPATH, "weather.DF.agg15min_5.5.csv"))


# Calculate Penman Monteith ET #
wthr.df1 <- calculateETref(wthr.DFagg15min)
wthr.ETref.df <- as.data.frame(wthr.df1)
empty.MAT <- matrix(nrow = 8, ncol = (ncol(et.obs)-nrow(wthr.df1)))

# select columns "Temp"  "RH"    "VPD"   "SR"    "WS"    "Tmax"  "Tmin"  "ETref"
empty.MAT.wthr.ETref <- as.data.frame(cbind(empty.MAT, t(wthr.ETref.df[,c(2:6, 9:11)])))
colnames(empty.MAT.wthr.ETref) <- colnames(et.obs)
wthr.ETref.ETobs <- as.data.frame(rbind(empty.MAT.wthr.ETref, et.obs))
write.csv(wthr.ETref.ETobs, paste0(opPATH, "wthr.ETref.ETobs.csv"))


### Start ET0 Ratio-based Thresholding ###
ET_ratio_mat <- generateETrRatio(x = wthr.ETref.ETobs)
write.csv(ET_ratio_mat, paste0(opPATH, "wthr.ETref.ETobs.Ratio.csv"))

ETr_baseFILE <- wthr.ETref.ETobs
baseFILE <- ET_ratio_mat

# Make the original date sequence in the LC time series from which irrg. dates will be filtered
act.dts <- seq(date(wthr.DFagg15min$TS[1]), date(wthr.DFagg15min$TS[nrow(wthr.DFagg15min)]), 1)

act.dts.ts <- rep(act.dts, each = 96)

act.dts.ts.df <- data.frame(org.dts.ts = c(rep(NA, 
                 (ncol(ETr_baseFILE) - length(act.dts.ts))),
                                           as.character(act.dts.ts)))
act.dts.ts.df$org.dts.ts <- as.Date(act.dts.ts.df$org.dts.ts)

rownames(act.dts.ts.df) <- colnames(ETr_baseFILE)

# Identify columns from ETr raw data which need to be filtered
irrg.cols <- which(act.dts.ts.df$org.dts.ts %in% irrg.dts)

# Prepare a copy of raw data as the Input File for the Smoothing process
ETr_smth_IP <- ETr_baseFILE
Ratio_smth_IP<- baseFILE

ETr_smth_IP <- ETr_smth_IP[ ,-irrg.cols]
Ratio_smth_IP <- Ratio_smth_IP[ ,-irrg.cols]


# Thresholding using 2 time-windows 06:30 to 18:30 and remaining #
by_Genotype <- Ratio_smth_IP[-c(1:8) ,c(3, 6:ncol(Ratio_smth_IP))] %>% group_by(Genotype)

by_Genotype_Mean <- by_Genotype  %>% summarise_all(~mean(.))

by_Genotype_Mean$Genotype <- factor(by_Genotype_Mean$Genotype)

geno.ETr <- as.data.frame(t(by_Genotype_Mean))

geno.ETr <- geno.ETr[-1, ]

geno.ETr <- as.data.frame(apply(geno.ETr, 2, function(x) {as.numeric(as.character(x))}))

colnames(geno.ETr) <- paste0("G_", by_Genotype_Mean$Genotype)

ETref <- as.data.frame(Ratio_smth_IP[8, 6:ncol(Ratio_smth_IP)])

baseDF <- as.data.frame(cbind(t(ETref), geno.ETr))
colnames(baseDF)[1] <- c("ETref")


# Check date format in rownames and run accordingly
baseDF$TS <- ymd_hms(rownames(baseDF)); 
baseDF$date <- date(baseDF$TS); 
baseDF$time <- strftime(baseDF$TS, format="%H:%M:%S", tz="UTC")
baseDF$solarRAD <- as.numeric(as.character(Ratio_smth_IP[4, 6:ncol(Ratio_smth_IP)])) ## Need to convert row to vector using c()

baseDF <- baseDF[ ,c((ncol(baseDF)-3), (ncol(baseDF)-2), (ncol(baseDF)-1), ncol(baseDF), 1:(ncol(baseDF)-4))]


### Steps to find the Thresholds and filter raw ETr ###

# 1. First find the set of unique dates.
# 2. Divide the whole data into 2 time-windows, TW (06:30 - 18:30, Rest).
# 3. For each TW, get the threshold values.
# 4. Partition ETr as per TWs and apply thresholds.
# 5. Filter ETr, impute and merge the TWs.


# 1. First find the set of unique dates
unq.dts <- unique(baseDF$date)

# 2. Divide the whole data i.e. baseDF into 2 time-windows, TW (06:30 - 18:30, Rest)
base_TW1 <- baseDF[baseDF$time >= "06:30:00" & baseDF$time < "18:30:00", ]
base_TW2 <- baseDF[!baseDF$time %in% base_TW1$time, ]

# 3. For each TW, get the threshold values. 
baseTW1_ThreshVALS <- genThreshVal(x = unq.dts, y = base_TW1)
baseTW2_ThreshVALS <- genThreshVal(x = unq.dts, y = base_TW2)

TW1.thresh <- floor(median(baseTW1_ThreshVALS$Q_75, na.rm = TRUE))
TW2.thresh <- ceiling(median(baseTW2_ThreshVALS$Q_75, na.rm = TRUE))

boxplot(baseTW1_ThreshVALS[,-ncol(baseTW1_ThreshVALS)], 
        main = paste0("TW1: Threshold = ", TW1.thresh))$out

boxplot(baseTW2_ThreshVALS[,-ncol(baseTW2_ThreshVALS)], 
        main = paste0("TW2: Threshold = ", TW2.thresh))$out

write.csv(baseTW1_ThreshVALS, paste0(opPATH, "ETrRatio_TW1_ThreshVALS.csv"))
write.csv(baseTW2_ThreshVALS, paste0(opPATH, "ETrRatio_TW2_ThreshVALS.csv"))


# 4. Make ETr partitions
ETr_TWs <- dataPART(x = ETr_smth_IP)

ETr_P1 <- ETr_TWs$P1 
ETr_P2 <- ETr_TWs$P2


# Make ETr Ratio partitions
ETr_Ratio_TWs <- dataPART(x = Ratio_smth_IP)

ETr_Ratio_P1 <- ETr_Ratio_TWs$P1
ETr_Ratio_P2 <- ETr_Ratio_TWs$P2


# Apply TW-specific thresholds
ETr_filt_TW1<- threshETr(x = TW1.thresh, y = ETr_Ratio_P1, z = ETr_P1)
ETr_filt_TW2<- threshETr(x = TW2.thresh, y = ETr_Ratio_P2, z = ETr_P2)

# Save filtered ETr data
ETr_filt_TW1_ord <- ordFiltETr(x = ETr_filt_TW1, y = ETr_smth_IP)
ETr_filt_TW2_ord <- ordFiltETr(x = ETr_filt_TW2, y = ETr_smth_IP)

ETr_filt_ord <- as.data.frame(rbind(ETr_filt_TW1_ord, ETr_filt_TW2_ord))

ETr_filt_ord$TS <- ymd_hms(rownames(ETr_filt_ord))

ETr_filt_ord <- ETr_filt_ord[order(ETr_filt_ord$TS), ]

ETr_filt_FILE <- ETr_smth_IP
subs.d <- as.matrix(t(ETr_filt_ord[,-1]))

ETr_filt_FILE[9:nrow(ETr_filt_FILE), 6:ncol(ETr_filt_FILE)] <- subs.d

write.csv(ETr_filt_FILE, paste0(opPATH, "ETr_Thresh_Filt.csv"))


# 5. Thresholded/Filtered ETr interpolation
ETr_filt_imputed_FILE <- ETr_filt_FILE
subs.d.imp <- subs.d
na.sum <- c()

for(i in 1:nrow(subs.d.imp))
{
  subs.d.imp[i, ] <- na.aggregate.default(subs.d.imp[i, ])
  na.sum[i] <- sum(is.na(subs.d.imp[i, ]))
}
print(paste0("#Sectors still with NA: ", (length(na.sum) - length(!na.sum == 0))))

ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE), 6:ncol(ETr_filt_imputed_FILE)] <- subs.d.imp

write.csv(ETr_filt_imputed_FILE, paste0(opPATH, "ETr_Thresh_Impt.csv"))



### Start Tr extraction to get raw Tr profiles ###
# Check format and set TS
pe.df$TS <- ymd_hms(pe.df$timestamp)

pe.df$date <- date(pe.df$TS)

sel.secs <- unique(na.omit(ETr_filt_imputed_FILE$old_unit))

pe.df.ETr <- pe.df[pe.df$date %in% unq.dts & pe.df$Sector %in% sel.secs, ]

# Select the following columns: "Sector", "Genotype", "Replicates", "LeafArea3D", "TS", "date"
pe.df.ETr <- pe.df.ETr[ ,c(1, 5, 7, 12, 17, 18)]

colnames(pe.df.ETr) <- c("old_unit", "Genotype", "Replicates", "LeafArea3D", "TS", "date")
pe.df.ETr$Genotype <- factor(pe.df.ETr$Genotype)
pe.df.ETr$Replicates <- factor(pe.df.ETr$Replicates)

pe.ETr.grpDT <- pe.df.ETr %>% group_by(old_unit, date, Genotype, Replicates) %>% 
                dplyr::summarise(Max = max(LeafArea3D, na.rm=TRUE))

names(pe.ETr.grpDT)[5] <- "LeafArea3D"

# make a matrix of LA3D of dim = ETr_core data, then calculate TR #
## LAI=((((3DLA/100) - 29.9)/0.36)*(1/0.26))/10000 ##
## T = (1-exp(-0.463*LAI))*ETr ##

LAI.mat <- matrix(NA, nrow = length(sel.secs), ncol = ncol(subs.d))
rownames(LAI.mat) <- sel.secs

# sort LAI.mat rownames as per the ETR_smooth file 
LAI.mat <- LAI.mat[order(match(rownames(LAI.mat), 
                  ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE), "old_unit"])),]

LAI.all.dates <- LAI.mat; unq.dts.copy <- unq.dts

LAI.all.dates <- LAI.all.dates[ , c(1:length(unq.dts))]
colnames(LAI.all.dates) <- c(as.character(unq.dts.copy))


# Calculate raw Transpiration, Tr
Tr_OP <- calculateTr(x = ETr_filt_imputed_FILE, y = pe.df.ETr, z = LAI.mat, d = unq.dts)

raw.trans.mat <- Tr_OP$Trans.mat

LA3D.all.dates <- Tr_OP$LA3D_TS

LAI.mat <- Tr_OP$LAI.mat


raw.trans <- ETr_filt_imputed_FILE

raw.trans[9:nrow(raw.trans), 6:ncol(raw.trans)] <- raw.trans.mat

write.csv(raw.trans, paste0(opPATH, "raw_Tr.csv"))

write.csv(LA3D.all.dates, paste0(opPATH, "3DLA_TS.csv"))

write.csv(LAI.mat, paste0(opPATH, "LAI_TS.csv"))


### Feature Extraction of RAW Transpiration Data ###
featuresRES <- getFeatures(x = raw.trans)

allFeatures <- featuresRES$allFeatures

# create H2 dataframe to store H2 est. of each feature for each day
F.He <- as.data.frame(matrix(NA, nrow = length(unq.dts), ncol = 15)) # Date-ROW, feature-COL
colnames(F.He) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                    "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                    "auc.prop.07-19", "auc.night", "cos.sim.index")
rownames(F.He) <- unq.dts

featureHeRES <- getFeatureHe(x = allFeatures, y = raw.trans, d = unq.dts, p = opPATH.raw)
write.csv(featureHeRES, paste0(opPATH, "rawTr_featureH2.csv"))

### save all features as feature Time Series ###
### Each feature set: dim(length(unq.dts) x (nrow(raw.trans)-8)) ###

## Prepare data for 'each feature'
maxET <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.maxET.6 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.07maxET <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.00.07 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.19.2345 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
curvmaxET <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
total.auc <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.10.15 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
sd.10.15 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.prop.10.15 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.07.19 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
sd.07.19 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.prop.07.19 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.night <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
cos.sim.index <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))

for (j in 1:(nrow(raw.trans)-8)){
  
  # "maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
  # "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
  # "auc.prop.07-19", "auc.night", "cos.sim.index"
  
  for(d in 1:length(unq.dts))
  {maxET[j, d] <- data.frame(allFeatures[[j]][d, 1])
  slope.maxET.6[j, d] <- data.frame(allFeatures[[j]][d, 2])
  slope.07maxET[j, d] <- data.frame(allFeatures[[j]][d, 3])
  slope.00.07 [j, d] <- data.frame(allFeatures[[j]][d, 4])
  slope.19.2345[j, d] <- data.frame(allFeatures[[j]][d, 5])
  
  curvmaxET[j, d] <- data.frame(allFeatures[[j]][d, 6])
  total.auc[j, d] <- data.frame(allFeatures[[j]][d, 7])
  auc.10.15[j, d] <- data.frame(allFeatures[[j]][d, 8])
  sd.10.15 [j, d] <- data.frame(allFeatures[[j]][d, 9])
  auc.prop.10.15[j, d] <- data.frame(allFeatures[[j]][d, 10])
  
  auc.07.19[j, d] <- data.frame(allFeatures[[j]][d, 11])
  sd.07.19[j, d] <- data.frame(allFeatures[[j]][d, 12])
  auc.prop.07.19[j, d] <- data.frame(allFeatures[[j]][d, 13])
  auc.night [j, d] <- data.frame(allFeatures[[j]][d, 14])
  cos.sim.index[j, d] <- data.frame(allFeatures[[j]][d, 15])
  }
  
} 

names(maxET)=names(slope.maxET.6)=names(slope.07maxET)=names(slope.00.07)=names(slope.19.2345)<-unq.dts
names(curvmaxET)=names(total.auc)=names(auc.10.15)=names(sd.10.15)=names(auc.prop.10.15)<-unq.dts
names(auc.07.19)=names(sd.07.19)=names(auc.prop.07.19)=names(auc.night)=names(cos.sim.index)<-unq.dts


write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], maxET)), 
          paste0(opPATH.raw, "maxET.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], slope.maxET.6)), 
          paste0(opPATH.raw, "slope.maxET.6.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], slope.00.07)), 
          paste0(opPATH.raw, "slope.00.07.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], slope.07maxET)), 
          paste0(opPATH.raw, "slope.07maxET.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], slope.19.2345)), 
          paste0(opPATH.raw, "slope.19.2345.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], curvmaxET)), 
          paste0(opPATH.raw, "curvmaxET.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], total.auc)), 
          paste0(opPATH.raw, "total.auc.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], auc.10.15)), 
          paste0(opPATH.raw, "auc.10.15.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], sd.10.15)), 
          paste0(opPATH.raw, "sd.10.15.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], auc.prop.10.15)), 
          paste0(opPATH.raw, "auc.prop.10.15.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], auc.07.19)), 
          paste0(opPATH.raw, "auc.07.19.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], sd.07.19)), 
          paste0(opPATH.raw, "sd.07.19.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], auc.prop.07.19)), 
          paste0(opPATH.raw, "auc.prop.07.19.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], auc.night)), 
          paste0(opPATH.raw, "auc.night.csv"))

write.csv(as.data.frame(cbind(raw.trans[9:nrow(raw.trans), 1:5], cos.sim.index)), 
          paste0(opPATH.raw, "cos.sim.index.csv"))

##### END #####


###### Stage III: Generate smooth ETr and obtain smooth Tr features ######
# ETr smoothing
smoothETrMAT <- smoothETr(x = subs.d.imp)

ETr_smoothFILE <- ETr_filt_imputed_FILE
ETr_smoothFILE[9:nrow(ETr_smoothFILE) , 6:ncol(ETr_smoothFILE)] <- smoothETrMAT/1000

write.csv(ETr_smoothFILE, paste0(opPATH, "ETr_smth.csv"))


# Calculate Tr from smooth ETr
Tr_OP <- calculateTr(x = ETr_smoothFILE, y = pe.df.ETr, z = LAI.mat, d = unq.dts)

smth.trans.mat <- Tr_OP$Trans.mat

smth.trans <- ETr_smoothFILE

smth.trans[9:nrow(smth.trans), 6:ncol(smth.trans)] <- smth.trans.mat

write.csv(smth.trans, paste0(opPATH, "smth_Tr.csv"))


### Feature Extraction of SMOOTH Transpiration Data ###
featuresRES <- getFeatures(x = smth.trans)

allFeatures <- featuresRES$allFeatures

# create H2 dataframe to store H2 est. of each feature for each day
F.He <- as.data.frame(matrix(NA, nrow = length(unq.dts), ncol = 15)) # Date-ROW, feature-COL
colnames(F.He) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                    "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                    "auc.prop.07-19", "auc.night", "cos.sim.index")
rownames(F.He) <- unq.dts

featureHeRES <- getFeatureHe(x = allFeatures, y = smth.trans, d = unq.dts, p = opPATH.smth)
write.csv(featureHeRES, paste0(opPATH, "smthTr_featureH2.csv"))


### save all features as feature Time Series ###
### Each feature set: dim(length(unq.dts) x (nrow(raw.trans)-8)) ###

## Prepare data for 'each feature'
maxET <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.maxET.6 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.07maxET <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.00.07 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
slope.19.2345 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
curvmaxET <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
total.auc <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.10.15 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
sd.10.15 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.prop.10.15 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.07.19 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
sd.07.19 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.prop.07.19 <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
auc.night <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))
cos.sim.index <- as.data.frame(matrix(nr = (nrow(raw.trans)-8), nc = length(unq.dts)))

for (j in 1:(nrow(raw.trans)-8)){
  
  # "maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
  # "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
  # "auc.prop.07-19", "auc.night", "cos.sim.index"
  
  for(d in 1:length(unq.dts))
  {maxET[j, d] <- data.frame(allFeatures[[j]][d, 1])
  slope.maxET.6[j, d] <- data.frame(allFeatures[[j]][d, 2])
  slope.07maxET[j, d] <- data.frame(allFeatures[[j]][d, 3])
  slope.00.07 [j, d] <- data.frame(allFeatures[[j]][d, 4])
  slope.19.2345[j, d] <- data.frame(allFeatures[[j]][d, 5])
  
  curvmaxET[j, d] <- data.frame(allFeatures[[j]][d, 6])
  total.auc[j, d] <- data.frame(allFeatures[[j]][d, 7])
  auc.10.15[j, d] <- data.frame(allFeatures[[j]][d, 8])
  sd.10.15 [j, d] <- data.frame(allFeatures[[j]][d, 9])
  auc.prop.10.15[j, d] <- data.frame(allFeatures[[j]][d, 10])
  
  auc.07.19[j, d] <- data.frame(allFeatures[[j]][d, 11])
  sd.07.19[j, d] <- data.frame(allFeatures[[j]][d, 12])
  auc.prop.07.19[j, d] <- data.frame(allFeatures[[j]][d, 13])
  auc.night [j, d] <- data.frame(allFeatures[[j]][d, 14])
  cos.sim.index[j, d] <- data.frame(allFeatures[[j]][d, 15])
  }
  
} 

names(maxET)=names(slope.maxET.6)=names(slope.07maxET)=names(slope.00.07)=names(slope.19.2345)<-unq.dts
names(curvmaxET)=names(total.auc)=names(auc.10.15)=names(sd.10.15)=names(auc.prop.10.15)<-unq.dts
names(auc.07.19)=names(sd.07.19)=names(auc.prop.07.19)=names(auc.night)=names(cos.sim.index)<-unq.dts


write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], maxET)), 
          paste0(opPATH.smth, "maxET.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], slope.maxET.6)), 
          paste0(opPATH.smth, "slope.maxET.6.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], slope.00.07)), 
          paste0(opPATH.smth, "slope.00.07.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], slope.07maxET)), 
          paste0(opPATH.smth, "slope.07maxET.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], slope.19.2345)), 
          paste0(opPATH.smth, "slope.19.2345.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], curvmaxET)), 
          paste0(opPATH.smth, "curvmaxET.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], total.auc)), 
          paste0(opPATH.smth, "total.auc.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], auc.10.15)), 
          paste0(opPATH.smth, "auc.10.15.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], sd.10.15)), 
          paste0(opPATH.smth, "sd.10.15.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], auc.prop.10.15)), 
          paste0(opPATH.smth, "auc.prop.10.15.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], auc.07.19)), 
          paste0(opPATH.smth, "auc.07.19.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], sd.07.19)), 
          paste0(opPATH.smth, "sd.07.19.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], auc.prop.07.19)), 
          paste0(opPATH.smth, "auc.prop.07.19.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], auc.night)), 
          paste0(opPATH.smth, "auc.night.csv"))

write.csv(as.data.frame(cbind(smth.trans[9:nrow(smth.trans), 1:5], cos.sim.index)), 
          paste0(opPATH.smth, "cos.sim.index.csv"))

end.time <- Sys.time()

print(paste0("Complete processing executed in: ", round((end.time-st.time), 2), "minutes"))
##### END #####

