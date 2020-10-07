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


### Exp41-SG ###
m.lc <- read.csv("./Exp41-SG/loadcells.csv")
meta.d <- read.csv("./Exp41-SG/Metadata_BC_Exp41.csv")
clm.df <- read.csv("./Exp41-SG/sensor_climate_data.csv")
sensor.unit.df <- read.csv("./Exp41-SG/sensor_unit_map.csv")
pe.df <- read.csv("./Exp41-SG/Exp41 Sorghum High-Low-TE Theno Dec 2019 - C_pe.csv")

HTP_data_Exp41_SG <- list(m.lc = m.lc, meta.d = meta.d, clm.df = clm.df,
                          sensor.unit.df = sensor.unit.df, pe.df = pe.df)
save(HTP_data_Exp41_SG, file = "HTP_data_Exp41_SG.RData")


### Exp42-PM ###
m.lc <- as.data.frame(read_xlsx("./Exp42-PM/Exp42 Loadcell Sorghum and millet Theo expt.xlsx", 2))
meta.d <- read.csv("./Exp42-PM/Exp42 Millet metdata.csv")
t.rh.ws <- as.data.frame(read_xlsx("./Exp42-PM/Exp42 climate TempRHWindSpeed.xlsx", 1))
solRAD <- as.data.frame(read_xlsx("./Exp42-PM/Exp42 climate solar radiation.xlsx", 1))
clm.df <- as.data.frame(rbind(t.rh.ws, solRAD))
sensor.unit.df <- read.csv("./Exp42-PM/Exp42_sensor_unit_map.csv")
pe.df <- as.data.frame(read_xlsx("./Exp42-PM/Exp42_PearlMilletPE.xlsx", 2))

HTP_data_Exp42_PM <- list(m.lc = m.lc, meta.d = meta.d, clm.df = clm.df,
                          sensor.unit.df = sensor.unit.df, pe.df = pe.df)
save(HTP_data_Exp42_PM, file = "HTP_data_Exp42_PM.RData")


### Exp42-SG ###
m.lc <- as.data.frame(read_xlsx("./Exp42-SG/Exp42 Loadcell Sorghum and millet Theo expt.xlsx", 1))
meta.d <- read.csv("./Exp42-SG/Exp42 sorghum Metdata.csv")
t.rh.ws <- as.data.frame(read_xlsx("./Exp42-SG/Exp42 climate TempRHWindSpeed.xlsx", 1))
solRAD <- as.data.frame(read_xlsx("./Exp42-SG/Exp42 climate solar radiation.xlsx", 1))
clm.df <- as.data.frame(rbind(t.rh.ws, solRAD))
sensor.unit.df <- read.csv("./Exp42-SG/Exp42_sensor_unit_map.csv")
pe.df <- as.data.frame(read_xlsx("./Exp42-SG/Exp42_SorghumPE.xlsx", 1))

HTP_data_Exp42_SG <- list(m.lc = m.lc, meta.d = meta.d, clm.df = clm.df,
                          sensor.unit.df = sensor.unit.df, pe.df = pe.df)
save(HTP_data_Exp42_SG, file = "HTP_data_Exp42_SG.RData")
