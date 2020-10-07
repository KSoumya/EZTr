library(readxl)

### Exp41-PM ###
m.lc <- read.csv("./Exp41-PM/loadcells.csv")
meta.d <- read.csv("./Exp41-PM/Metadata_BC_Exp41.csv")
clm.df <- read.csv("./Exp41-PM/sensor_climate_data.csv")
sensor.unit.df <- read.csv("./Exp41-PM/sensor_unit_map.csv")
pe.df <- read.csv("./Exp41-PM/Exp41 PMiGAP Theo VV Dec 2019 - B_pe.csv")

HTP_data_Exp41_PM <- list(m.lc = m.lc, meta.d = meta.d, clm.df = clm.df,
                          sensor.unit.df = sensor.unit.df, pe.df = pe.df)
save(HTP_data_Exp41_PM, file = "HTP_data_Exp41_PM.RData")
