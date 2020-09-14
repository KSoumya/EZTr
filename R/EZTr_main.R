##################
# EZTr_main #
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
#' Functions used for each step:
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
#' @param HTP_data \code{list} of 5 dataframes
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
#' @author Soumyashree Kar
#'
#' @references
#'
#' Vadez, V., Kholová, J., Hummel, G., Zhokhavets, U., Gupta, S. K., &
#' Hash, C. T. (2015). LeasyScan: a novel concept combining 3D imaging and
#' lysimetry for high-throughput phenotyping of traits controlling plant water
#' budget. Journal of Experimental Botany, 66(18), 5581-5593.
#'
#'
#' @examples
#'
#' load("HTP_data.RData")
#'
#' m.lc <- HTP_data$m.lc
#' meta.d <- HTP_data$meta.d
#' t.rh.ws <- HTP_data$t.rh.ws
#' solRAD <- HTP_data$solRAD
#' pe.df <- HTP_data$pe.df
#'
#'
#' \dontrun{
#'
#' EZTr_OP <- EZTr_main(HTP_data, lastDate="2017-03-07",
#'                      irrg.dts=c("2017-03-07"),
#'                      Date1="2017-02-18 23:46:00",
#'                      Date2="2017-03-07 23:45:00",
#'                      opPATH="D:/EZTr/RES/",
#'                      opPATH.smth="D:/EZTr/RES/FeaturesTimeSeries/")
#'
#' }
#'
#' @export
#'



EZTr_main <- function(HTP_data, lastDate="2017-03-07",
                      irrg.dts=c("2017-03-07"),
                      Date1="2017-02-18 23:46:00",
                      Date2="2017-03-07 23:45:00",
                      opPATH="D:/EZTr/RES/",
                      opPATH.smth="D:/EZTr/RES/FeaturesTimeSeries/"){
  # Get load cells data i.e. weights of sector
  m.lc <- HTP_data$m.lc

  # Get Genotype and Exp design metadata
  meta.d <- HTP_data$meta.d

  # Get Temperature, Relative humidity% and Wind Speed
  t.rh.ws <- HTP_data$t.rh.ws

  # get Solar Radiation
  solRAD <- HTP_data$solRAD

  # Get 3D leaf area values
  pe.df <- HTP_data$pe.df

  # Get the list of species from the metadata
  species.nm <- unique(meta.d$Species)

  # Include the species ID e.g. 1 for 'Chickpea' as per the data
  meta.d.sp <- meta.d[meta.d$Species==species.nm[1], ]

  # Set the last date for the LC data from the exp duration
  lastDate <- lastDate

  # Input the dates of irrigation
  irrg.dts <- as.Date(irrg.dts)


  ####### Stage-I: Start processing LC data to generate ETr matrix #######
  st.time <- Sys.time()

  # Extract matrix of load cells data #
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
  meta.LCDF <- meta.d.LCmat[order(match(meta.d.LCmat$unit,
                                        rownames(LC.MAT.f.t))), ]

  LC.MAT.raw <- as.data.frame(cbind(meta.LCDF, LC.MAT.f.t))

  # Convert load cell time series to ETr profiles
  ETr_Meta <- curateRawLCgetET(x = LC.MAT.f, y = meta.LCDF)

  meta.LCDF <- ETr_Meta[, 1:5]

  # Convert ETr in grams to mm #
  ETr.mm_Meta <- ETr_Meta
  et.obs <- ETr.mm_Meta[ ,6:ncol(ETr.mm_Meta)]
  ETr.mm_Meta[ ,6:ncol(ETr.mm_Meta)] <- et.obs*4/1000

  # Identify error plots from ETr values using the similar method as above #
  ETr_err.sec.info <- filterETrExtremeCols(x = ETr.mm_Meta, y = meta.LCDF)

  err.sec.nm <- ETr_err.sec.info$ETr_err.sec.NM

  err.sec.meta <- ETr_err.sec.info$ETr_err.sec.META

  # Remove the err.cols i.e. sectors with extreme values #
  ETr_Meta_ERRsec.rmvd <- ETr.mm_Meta

  ETr_Meta_ERRsec.rmvd <- ETr_Meta_ERRsec.rmvd[!ETr_Meta_ERRsec.rmvd$unit %in% err.sec.nm, ]
  #####################################################################


  ####### Stage-II: Start processing Weather data to obtain ETref and ETr ratio matrix #######
  clm.df <- as.data.frame(rbind(t.rh.ws, solRAD))
  unq.clm.var <- unique(clm.df$variable)

  # Set time adjustment factor e.g."- 5.5*60*60", based on apt Time Zone
  #Extract weather variables individually #
  temperature.DF <- extractWthrVar(x = 1, y = clm.df)
  temperature.DF$ts <- temperature.DF$ts - 5.5*60*60

  relHUM.DF <- extractWthrVar(x = 2, y = clm.df)
  relHUM.DF$ts <- relHUM.DF$ts - 5.5*60*60

  windS.DF <- extractWthrVar(x = 3, y = clm.df)
  windS.DF$ts <- windS.DF$ts - 5.5*60*60

  solarRad.DF <- extractWthrVar(x = 4, y = clm.df)
  solarRad.DF$ts <- solarRad.DF$ts - 5.5*60*60


  # Read the ET obs data and get the first-last time stamps #
  et.obs <- ETr_Meta_ERRsec.rmvd
  et.f.ts <- colnames(et.obs[6]); et.l.ts <- colnames(et.obs[ncol(et.obs)])

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

  for(i in 1:nrow(weather.DF))
  {
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
  }

  # Preprocess each weather variable except VPD #
  weather.DF[ ,2] <- prepcsWthr(x = weather.DF, y = 2) # temperature
  weather.DF[ ,3] <- prepcsWthr(x = weather.DF, y = 3) # relative humidity
  weather.DF[ ,5] <- prepcsWthr(x = weather.DF, y = 5) # solar radiation
  weather.DF[ ,6] <- prepcsWthr(x = weather.DF, y = 6) # wind speed


  # Compute VPD and insert into the weather DF #
  SVP <- 610.7*(10^(7.5*weather.DF[ ,2]/(237.3+weather.DF[ ,2])))
  VPD <- ((1 - (weather.DF[ ,3]/100))*SVP)/1000
  weather.DF[ ,4] <- VPD

  wthr.DFxts.TS <- xts(weather.DF, order.by = as.POSIXct(weather.DF$TS, format="%Y-%m-%d %H:%M"))
  wthr.DFagg15min = highfrequency::aggregatets(wthr.DFxts.TS, on="minutes",k=15, dropna=TRUE)

  # Need these 3 lines if TS is like 00:14:00 - 00:29:00, etc.
  # Else save directly
  ts.x <- ymd_hms(wthr.DFagg15min$TS)
  ts.x <- ts.x + hms("00:01:00")
  wthr.DFagg15min$TS <- as.character(ts.x)

  # Check and Filter duration of weather data as per observed ET
  et.obs.TS <- colnames(et.obs)[6:ncol(et.obs)]

  # IMP: Check and set the date format of et.obs.TS #
  et.obs.TS <- ymd_hms(et.obs.TS)
  et.obs.TS.chr <- as.character(et.obs.TS)
  wthr.DFagg15min <- wthr.DFagg15min[wthr.DFagg15min$TS %in% et.obs.TS.chr, ]

  # Calculate Penman Monteith ET #
  wthr.df1 <- calculateETref(wthr.DFagg15min)

  wthr.ETref.df <- as.data.frame(wthr.df1)

  empty.MAT <- matrix(nrow = 8, ncol = (ncol(et.obs)-nrow(wthr.df1)))

  # select columns "Temp"  "RH"    "VPD"   "SR"    "WS"    "Tmax"  "Tmin"  "ETref"
  empty.MAT.wthr.ETref <- as.data.frame(cbind(empty.MAT, t(wthr.ETref.df[,c(2:6, 9:11)])))

  colnames(empty.MAT.wthr.ETref) <- colnames(et.obs)

  wthr.ETref.ETobs <- as.data.frame(rbind(empty.MAT.wthr.ETref, et.obs))


  ##### Start ET0 Ratio-based Thresholding #####
  # Generate the ratio #
  ET_ratio_mat <- generateETrRatio(x = wthr.ETref.ETobs)

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


  ######### Thresholding using 2 time-windows 06:30 to 18:30 and remaining ##########

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


  ##### Steps to find the Thresholds #####

  # 1. First find the set of unique dates.
  # 2. Divide the whole data into 2 time-windows, TW (06:30 - 18:30, Rest).
  # 3. For each TW, need to compute OLS regression between ETref and the genotypes, get slope and intercept,
  #    before and after boxplot-outlier removal, the 5 quantile values and mean of solarRAD. Hence, 10 variables.
  #    For each TW create a matrix of size: #unique_dates X #variables.
  # 4. Obtain 2 matrices of time series of 10 variables.
  # 5. Extract the thresholds for each TW and perform filtering.


  # 1. First find the set of unique dates
  unq.dts <- unique(baseDF$date)

  # 2. Divide the whole data i.e. baseDF into 2 time-windows, TW (06:30 - 18:30, Rest)
  base_TW1 <- baseDF[baseDF$time >= "06:30:00" & baseDF$time < "18:30:00", ]
  base_TW2 <- baseDF[!baseDF$time %in% base_TW1$time, ]

  # 3. Find Thresholds
  # For each TW, need to compute OLS regression between ETref and the genotypes, get slope and intercept,
  # before and after boxplot-outlier removal, the 5 quantile values and mean of solarRAD. Hence, 10 variables.
  baseTW1_ThreshVALS <- genThreshVal(x = unq.dts, y = base_TW1)
  baseTW2_ThreshVALS <- genThreshVal(x = unq.dts, y = base_TW2)

  TW1.thresh <- ceiling(median(baseTW1_ThreshVALS$Q_75, na.rm = TRUE))
  TW2.thresh <- ceiling(median(baseTW2_ThreshVALS$Q_75, na.rm = TRUE))

  boxplot(baseTW1_ThreshVALS[,-ncol(baseTW1_ThreshVALS)],
          main = paste0("TW1: Threshold = ", TW1.thresh))$out

  boxplot(baseTW2_ThreshVALS[,-ncol(baseTW2_ThreshVALS)],
          main = paste0("TW2: Threshold = ", TW2.thresh))$out

  ##### Start ETr data partioning, thresholding and interpolation #####
  # Make ETr partitions
  ETr_TWs <- dataPART(x = ETr_smth_IP)

  ETr_P1 <- ETr_TWs$P1
  ETr_P2 <- ETr_TWs$P2

  # Make ETr Ratio partitions
  ETr_Ratio_TWs <- dataPART(x = Ratio_smth_IP)

  ETr_Ratio_P1 <- ETr_Ratio_TWs$P1
  ETr_Ratio_P2 <- ETr_Ratio_TWs$P2

  # Step3: Applying TW-specific thresholds
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

  # Step4: Interpolation of Thresholded/Filtered ETr
  ETr_filt_imputed_FILE <- ETr_filt_FILE
  subs.d.imp <- subs.d
  na.sum <- c()

  for(i in 1:nrow(subs.d.imp))
  {
    subs.d.imp[i, ] <- na.spline(subs.d.imp[i, ])

    na.sum[i] <- sum(is.na(subs.d.imp[i, ]))
  }
  print(paste0("#Sectors still with NA: ", (length(na.sum) - length(!na.sum == 0))))

  ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE),
                        6:ncol(ETr_filt_imputed_FILE)] <- subs.d.imp

  #####################################################################


  ### Stage III ###
  ####### Start Smoothing and smooth Tr Extraction #######

  #Step5: ETr smoothing #
  smoothETrMAT <- smoothETr(x = subs.d.imp)

  ETr_smoothFILE <- ETr_filt_imputed_FILE

  ETr_smoothFILE[9:nrow(ETr_smoothFILE), 6:ncol(ETr_smoothFILE)] <- smoothETrMAT/1000
  ###########################################################


  ##### Start Tr Calculation #####
  pe.df$TS <- dmy_hm(pe.df$timestamp)

  pe.df$date <- date(pe.df$TS)

  sel.secs <- unique(na.omit(ETr_filt_imputed_FILE$unit))

  pe.df.ETr <- pe.df[pe.df$date %in% unq.dts & pe.df$Sector %in% sel.secs, ]

  pe.df.ETr <- pe.df.ETr[ ,c(1, 5, 7, 12, 17, 18)]

  colnames(pe.df.ETr) <- c("old_unit", "Genotype", "Replicates", "LeafArea3D", "TS", "date")
  pe.df.ETr$Genotype <- factor(pe.df.ETr$Genotype)
  pe.df.ETr$Replicates <- factor(pe.df.ETr$Replicates)

  pe.ETr.grpDT <- pe.df.ETr %>% group_by(old_unit, date, Genotype, Replicates) %>%
    dplyr::summarise(Max = max(LeafArea3D, na.rm=TRUE))

  names(pe.ETr.grpDT)[5] <- "LeafArea3D"

  LAI.mat <- matrix(NA, nrow = length(sel.secs), ncol = ncol(subs.d))

  rownames(LAI.mat) <- sel.secs

  # sort LAI.mat rownames as per the ETR_smooth file
  LAI.mat <- LAI.mat[order(match(rownames(LAI.mat),
                                 ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE), "old_unit"])),]

  LAI.all.dates <- LAI.mat; unq.dts.copy <- unq.dts

  LAI.all.dates <- LAI.all.dates[ , c(1:length(unq.dts))]
  colnames(LAI.all.dates) <- c(as.character(unq.dts.copy))


  # Calculate TR from smooth ETr
  Tr_OP <- calculateTr(x = ETr_smoothFILE, y = pe.df.ETr, z = LAI.mat, d = unq.dts)

  smth.trans.mat <- Tr_OP$Trans.mat

  smth.TR.mat <- Tr_OP$TR.mat

  LAI.OP <- Tr_OP$LAI.mat

  LA_OP <- Tr_OP$LA3D_TS

  smth.trans.rate <- smth.trans <- ETr_smoothFILE

  smth.trans[9:nrow(smth.trans), 6:ncol(smth.trans)] <- smth.trans.mat

  smth.trans.rate[9:nrow(smth.trans.rate), 6:ncol(smth.trans.rate)] <- smth.TR.mat


  ##### Feature Extraction of SMOOTH Transpiration Data #####
  featuresRES <- getFeatures(x = smth.trans)

  allFeatures <- featuresRES$allFeatures

  # create H2 dataframe to store H2 est. of each feature for each day
  F.He <- as.data.frame(matrix(NA, nrow = length(unq.dts), ncol = 15)) # Date-ROW, feature-COL
  colnames(F.He) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET",
                      "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",
                      "auc.prop.07-19", "auc.night", "cos.sim.index")
  rownames(F.He) <- unq.dts

  # Get Feature Heritability
  featureHeRES <- getFeatureHe(x = allFeatures, y = smth.trans,
                               d = unq.dts, p = opPATH.smth)

  ##### save all features as feature Time Series #####
  # Each feature set: dim(length(unq.dts) x (nrow(smth.trans)-8)) #

  ## Prepare data for 'each feature'
  maxET <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  slope.maxET.6 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  slope.07maxET <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  slope.00.07 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  slope.19.2345 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  curvmaxET <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  total.auc <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  auc.10.15 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  sd.10.15 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  auc.prop.10.15 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  auc.07.19 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  sd.07.19 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  auc.prop.07.19 <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  auc.night <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))
  cos.sim.index <- as.data.frame(matrix(nr = (nrow(smth.trans)-8), nc = length(unq.dts)))

  names(maxET)=names(slope.maxET.6)=names(slope.07maxET)=names(slope.00.07)=names(slope.19.2345)<-unq.dts
  names(curvmaxET)=names(total.auc)=names(auc.10.15)=names(sd.10.15)=names(auc.prop.10.15)<-unq.dts
  names(auc.07.19)=names(sd.07.19)=names(auc.prop.07.19)=names(auc.night)=names(cos.sim.index)<-unq.dts

  for (j in 1:(nrow(smth.trans)-8)){

    # Features:
    # "maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET",
    # "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",
    # "auc.prop.07-19", "auc.night", "cos.sim.index"

    for(d in 1:length(unq.dts)){
    maxET[j, d] <- data.frame(allFeatures[[j]][d, 1])
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

  featureList <- list(maxET, slope.maxET.6, slope.07maxET, slope.00.07,
                      slope.19.2345, curvmaxET, total.auc, sd.10.15,
                      auc.prop.10.15, auc.07.19, auc.prop.07.19, auc.night,
                      cos.sim.index)

  end.time <- Sys.time()

  print(paste0("Complete processing executed in: ", round((end.time-st.time), 2), "minutes"))
  ##### END #####

  list(LCraw_tsmeta = LC.MAT.TSinfo,
       LCraw = LC.MAT.raw,
       ETrmm_Obs = ETr.mm_Meta,
       ETr_Obs_ERR.SEC.nms = err.sec.meta,
       ETmm_Obs_FINAL = ETr_Meta_ERRsec.rmvd,
       Wthr_agg15min = wthr.DFagg15min,
       Wthr.ETref.ETobs = wthr.ETref.ETobs,
       ETrRatio_TW1_ThreshVALS = baseTW1_ThreshVALS,
       ETrRatio_TW2_ThreshVALS = baseTW2_ThreshVALS,
       IrrgFilt_ETr = ETr_filt_FILE,
       IrrgFilt_ETr_Imptd = ETr_filt_imputed_FILE,
       ETr_smth = ETr_smoothFILE,
       Tr = smth.trans,
       featureH2 = featureHeRES,
       eachFeature_TS=featureList)
}

