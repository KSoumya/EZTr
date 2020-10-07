extractRawLCmatrix <- function(x, y, z) {

  a <- Sys.time()

  unq.secs <- unique(x$unit)

  m.lc.df <- x

  # check the date/time format from head(m.lc) and run ymd_hms/dmy_hm
  m.lc.df$ts <- dmy_hm(m.lc.df$timestamp)
  # m.lc.df <- m.lc.df[,-3] # required only if timestamp == ts

  m.lc.df$date <- date(m.lc.df$ts)

  m.lc.df$time <- strftime(m.lc.df$ts, format="%H:%M:%S", tz="UTC")

  head(m.lc.df)

  unq.sec.dates <- list() # create list to store the sorted unique dates of each sector

  unq.sec.dates.ln <- c() # create vector to store the number of unique dates of each sector

  # store date values
  for (i in 1:length(unq.secs))
  {
    unq.sec.dates[[i]] <- sort(unique(m.lc.df$date[m.lc.df$unit==unq.secs[i]]))
    unq.sec.dates.ln <- c(unq.sec.dates.ln, length(unq.sec.dates[[i]]))
  }

  # find mode of the length vector: unq.sec.dates.ln and remove sectors less than mode
  sec.dtLen <- data.frame(sec = unq.secs, dt.len = unq.sec.dates.ln)

  # create function to calculate mode:
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  max.dt <- as.numeric(getmode(sec.dtLen$dt.len))

  # select from Metadata: "unit", "old_unit", "Genotype, "G..Alias", "Replicates"
  sec.dtLenMETA <- y[y$unit %in%  sec.dtLen$sec, c(1,2,6,7,8)]

  # reorder rows of 'meta.d.sp' according to rownames/unit of LC.MAT.f.t
  sec.dtLenMETA <- sec.dtLenMETA[order(match(sec.dtLenMETA$unit, sec.dtLen$sec)), ]

  # Append no. of days data available to Metadata of each unit in the LC file
  sec.dtLenMETA <- as.data.frame(cbind(sec.dtLenMETA, sec.dtLen$dt.len))

  colnames(sec.dtLenMETA)[ncol(sec.dtLenMETA)] <- "#DaysData"

  # Remove sectors with very few days data
  err.sec <- sec.dtLen$sec[sec.dtLen$dt.len < floor(0.25*(max.dt))] # find the error sectors

  err.sec.loc <- which(sec.dtLen$dt.len < floor(0.25*(max.dt))) # find location of the error sectors

  m.lc.df <- m.lc.df[!m.lc.df$unit %in% err.sec, ] # remove from data

  unq.secs.good <- unq.secs[!unq.secs %in% err.sec] # create good unique sectors vector


  # get the first date and last date of each of the sector

  for (i in 1:length(unq.secs))
  {
    if(!(i %in% err.sec.loc))

    {

      sec.dtLen$F.dt[i] <- as.character(unq.sec.dates[[i]][1])

      sec.dtLen$L.dt[i] <- as.character(unq.sec.dates[[i]][length(unq.sec.dates[[i]])])
    }

    else
    {sec.dtLen$F.dt[i] <- NA; sec.dtLen$L.dt[i] <- NA}
  }

  # update sec.dtLen dataframe
  sec.dtLen <- sec.dtLen[is.na(sec.dtLen$F.dt)==FALSE, ]

  # Append metadata to updated sec.dtLen dataframe
  sec.dtLenMETA.FullUnits <- as.data.frame(cbind(sec.dtLenMETA[sec.dtLenMETA$unit %in% sec.dtLen$sec, ],
                                                 sec.dtLen[ ,-c(1:2)])) # Remove 'unit', 'dt.len' cols

  # Change time-stamps minutes != 00:00 or 00:15 or 00:30 or 00:45

  ## Create a sequence of time-stamps with 15min interval for 1 day
  TS_base<-as.data.frame(as.character(seq(ymd_hm(paste0(sec.dtLen$F.dt[1]," ",'00:00')),
                                          ymd_hm(paste0(sec.dtLen$F.dt[1]," ",'23:45')), by = '15 mins')))
  names(TS_base)[1]<-c("int.val")

  TS_base$time <- strftime(TS_base$int.val, format="%H:%M:%S", tz="UTC")

  # Since, the hms values slightly differ in the orginal dataset than the ideal 15min interval values,
  # replace them with the TS_base strftime (hms) values

  hms.ts.base<-unique(TS_base$time)
  names(hms.ts.base)<-c("time")

  print("Loadcell matrix timestamp mapping status")

  i<-nrow(m.lc.df)
  pbar <- create_progress_bar('text')
  pbar$init(i)

  for(i in 1:nrow(m.lc.df))
  {
    if(! m.lc.df$time[i] %in% hms.ts.base)
    {
      j <- which.min(abs(strptime(TS_base$time, "%H:%M") - strptime(m.lc.df$time[i], "%H:%M"))) # find which value in ts-base-vector is nearest to each-DP

      m.lc.df$time[i] <- TS_base$time[j]} # assign the nearest ts-base-vector value to that DP

    pbar$step()
  }

  # combine date and new time and convert to time object
  m.lc.df$TS.n <- ymd_hms(paste(m.lc.df$date, m.lc.df$time))

  # Remove rows with the same TS as the previous, so that T.diff.new != 0
  m.lc.df<-m.lc.df[!duplicated(m.lc.df[,c('TS.n', 'unit')]),]

  # create an ideal TS for ENTIRE DURATION of the LC dataset
  if(length(unique(sec.dtLen$F.dt))==1 & length(unique(sec.dtLen$L.dt))==1)
  {

    F.dt <- unique(sec.dtLen$F.dt); L.dt <- unique(sec.dtLen$L.dt)

  }else

  {F.dt <- min(sec.dtLen$F.dt); L.dt <- max(sec.dtLen$L.dt)}


  # Set Last date as per experiment details: For Exp-41, last date = 2020-01-28

  L.dt <- z

  TS_ALL<-as.data.frame(as.character(seq(ymd_hm(paste0(F.dt," ",'00:00')),
                                         ymd_hm(paste0(L.dt," ",'23:45')), by = '15 mins')))

  names(TS_ALL)[1]<-c("TS.n") # MUST be the same as in original data set m.lc.df

  TS_ALL$TS.n <- ymd_hms(TS_ALL$TS.n)

  # create empty dataframe to store all values

  LC.MAT <- as.data.frame(matrix(nrow = length(TS_ALL$TS.n), ncol = nrow(sec.dtLen)+1))

  LC.MAT[ ,1] <- TS_ALL$TS.n; names(LC.MAT)[1] <- "TS"


  # extract TS.n and mass data for each sector and merge with TS_ALL
  # First ensure column indices of 'TS.n', 'unit', 'mass' from "m.lc.df" dataframe and then execute the loop
  for (n in 1:(nrow(sec.dtLen)))
  {
    tmp.lc <- m.lc.df[(m.lc.df$unit == sec.dtLen$sec[n] &

                         (m.lc.df$date < L.dt && m.lc.df$TS.n > m.lc.df$TS.n[1] - hms("00:15:00"))),

                      c(10, 1, 6)] #('TS.n', 'unit', 'mass')

    df <- merge(x = tmp.lc, y = TS_ALL, by = "TS.n", all = TRUE) # perform outer join to merge by id=TS.n

    df <- df[df$TS.n %in% TS_ALL$TS.n, ] # subset after outer join to ensure that NAs don't add extra rows

    LC.MAT[ ,n+1] <- df[df$unit == sec.dtLen$sec[n], 3] #'Mass'

    names(LC.MAT)[n+1] <- as.character(sec.dtLen$sec[n])

  }

  LC.MAT.f <- LC.MAT[LC.MAT$TS > (ymd_hms(paste0(F.dt," ",'23:45:00'))), ]

  b <- Sys.time()

  print(paste0("Loadcell matrix created in: ", round((b-a), 2), " minutes"))

  list(LC.MAT.f = LC.MAT.f, LC_tsmeta <- sec.dtLenMETA.FullUnits)
}
