library(dplyr)
library(ggplot2)
library(mgcv)
library(tictoc)
library(lubridate)
library(tidyr)
library(purrr)
library(gratia) 
library(stringr)
#remotes::install_github("christyray/sciscales")
library(sciscales)

# NOTE - SEPARATE rda FILE SUPPLIED - DON'T NEED TO QUERY API
# ITS UNDERGONE AN UPDATE AND ddspwq PACKAGE NOT CURRENTLY WORKING
# SKIP TO ZZZ BELOW

# using ddspWQ package to access data on Open WIMS (WQ data explorer)
# more details on the patform here: https://environment.data.gov.uk/water-quality/view/landing
#remotes::install_github('a-jone5/ddspWQ')
#library(ddspWQ)


# If can't download ddspWQ package from github then can used the saved version of the data (.rda file) 
# skip to ZZZ comment below

# obviously change this to a suitable location for yourself
#my.save.location <- 'O:/NCES Team/Evidence Synthesis/R/Trend modelling examples/HMS from open data.rda'


# Code by Mike Dunbar to illustrate fitting a gamm with censored normal family to 
# water quality time series from multiple sites (Harmonised Monitoring Scheme sites for England)
# to look at smoothed trends over time and the effects of seasonality
# extracted from the EA WQ data explorer
# separate models for each determinand
# site (sampling point) is considered as a random intercept
# we fit one simple model to each det

# we log-transform the raw data, then fit, then transform back as there is no
# censored lognormal family available this does mean that we are effectively
# modelling the geometric mean of the data

# models can get more complex, e.g. 
# factor smooth models (separate trajectory for each site) 
# seasonality changing over time, separate models for RBDs, including flow as a
# covariate etc but we keep it simple here we do however look at a couple of
# different predictions for each model: with and without the seasonality term
# we extract data for oprthophosphate, dissolved zinc and PFOS

# NOTE - PFOS data not officially part of this network so there's alot less
# data as it's collected at a subset of sites for various other monitoring
# programmes. It does however illustrate the effect of sample size in the final
# predictions

# TO DO:
# look into neighbourhood cross-validation (NCV) and randomised quantile residuals (RQR) for model evaluation
# Investigate bcg() family for censored box-cox (in mgcv 1.9-4 onwards)

# sites we want to extact - these are the HMS sites for England and data goes
# back a long time
# although data in the WQ archive not available before yr 2000
my.sites.df <-
  structure(list(
    WIMSCode =
    c("26M31", "51M01", "53M14", "ANCOC", "BL01", "BUR120", "CH01", "CL01",
      "NENE550W", "NENE640D", "ST03", "ST11", "WELL280T", "WELL420C", "WEN250",
      "WITHM", "00025085", "50022", "00055140", "04778460", "13598380",
      "23314180", "26944690", "36701570", "36741880", "36768280", "38473020",
      "46247300", "49690300", "54509300", "59000500", "70256300", "41000079",
      "42300080", "42500030", "43100155", "43200040", "43300001", "43400033",
      "44100119", "44500006", "45200010", "45400013", "49000137", "49100479",
      "49100488", "49200090", "49301422", "49301589", "49301624", "49301842",
      "49400409", "49400424", "49500343", "49600142", "49700156", "3",
      "88000879", "88002001", "88002065", "88002348", "88002634", "88002884",
      "88003147", "88003442", "88003521", "88003532", "88003561", "88003872",
      "88004024", "88004397", "88004563", "88005740", "88006220", "88006264",
      "88006451", "88006479", "88020376", "88021071", "E0000362", "E0000807",
      "E0001255", "E0001545", "F0002075", "F0002151", "F0002783", "F0002886",
      "G0003786", "G0003885", "G0003989", "PCHR0016", "PCNR0025", "PKER0025",
      "PLDR0029", "PLER0053", "PLER0057", "PLER0067", "PLER0076", "PMLR0022",
      "PRGR0038", "PTAR0022", "PTHR0075", "PTHR0107", "PTHR0113", "PWER0030",
      "50280271", "50370169", "50450129", "50590127", "60250424", "70220159",
      "70420116", "70540110", "70540224", "70620154", "70720104", "70826005",
      "72920121", "73020127", "73030120", "73080442", "81120120", "81231133",
      "81250144", "81520205", "81930120", "81950522", "82310103", "82528005",
      "91251605", "A1260103", "A3190103", "E1008100", "Z1010706", "Z1012401"),
   Region =
    c("AN", "AN", "AN", "AN", "AN", "AN", "AN", "AN", "AN", "AN", "AN", "AN",
      "AN", "AN", "AN", "AN", "MI", "MI", "MI", "MI", "MI", "MI", "MI", "MI",
      "MI", "MI", "MI", "MI", "MI", "MI", "MI", "MI", "NE", "NE", "NE", "NE",
      "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NE",
      "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NW", "NW", "NW", "NW",
      "NW", "NW", "NW", "NW", "NW", "NW", "NW", "NW", "NW", "NW", "NW", "NW",
      "NW", "NW", "NW", "NW", "NW", "NW", "NW", "SO", "SO", "SO", "SO", "SO",
      "SO", "SO", "SO", "SO", "SO", "SO", "TH", "TH", "TH", "TH", "TH", "TH",
      "TH", "TH", "TH", "TH", "TH", "TH", "TH", "TH", "TH", "SW", "SW", "SW",
      "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW",
      "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW", "SW",
      "SW", "SW", "SW")),
  row.names = c(NA, -135L), class = "data.frame") |>
  mutate(site.notation = paste(Region, WIMSCode, sep = '-'))

# dets to extract - later we use the combination of the code and the det name
det.list <- c('0117', # nitrate - should be minimal censoring
              '0180', # orthophosphate
              '3408', # zinc, dissolved
              '9276') # PFOS (not routinely monitored at these sites but could be some available)

# tic()
# samples <- 
#   fetch_sample_res(site_notation = my.sites.df$site.notation, dets = det.list)
# toc() # can take a little while

save(samples, det.list, my.sites.df, file = 'HMS from open data.rda')

