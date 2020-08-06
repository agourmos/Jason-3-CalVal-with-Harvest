# Jason-3-CalVal-with-Harvest
This code set pulls Jason-3 SLA data from the RADS database by remko and calibrate and validates it using SLA data aquired from a lidar located on the Harvest platform off the coast of California.

# NOTES

# CalVal_Main.m

- The main code takes multipel sets of data from different sources and outputs them into a figure to compare the data sets. The data comes from NOAA bubbler, lidar, radar, Jason-3 (from RADs), and lidar data from CU's very own lidar located on the Harvest oil platform.

- SECTION A: In this section of the code, majority of variables are defined to use throughout the code. This includes empty vectors and matrices as well as the longitude and       latitude of HArvest and the latitude of the bins (which was hardcoded)

- SECTION B: This section extracts data from the lidar data at the Harvest platform from a .csv file. The .csv file is generated from code from located on the CODs. The code takes in overflight times, referenced to Jason-3 flying over Harvest, and pulls the data from the CODs server and inputs it into a .csv file form those exact times. The .csv file also contains bubbler and radar data from NOAA on the Harvest platform. This section also removes outliars from the these variables and replaces them with NaN values.

- SECTION C: This section pulls data from the netCDF file which comes from RADs. The first step in this section is to reformat the time given in all data sets so they can be compared. The next step compares those dates in order to retrieve comparable data.

- SECTION D: In this section, the data is is split into cycles through the process of binning. This section also further corrects the data by averaging the 4, 9, and 16 closest points to harvets into a single value. However one of the graphs created shows that the more points you average over, the less precise the data becomes. This section also interpolates the data to show the trend of rms of the data. Fianlly it also fills in amy dates that are missing according to the lidar data when compared Jason-3 data.

- SECTION E: These are the figure numbers and what they represent.
  -- figure 1 shows the flight path, Jason3 vs. Lidar data, and the difference of the 2 set "HarvestCalVal.png"
  -- figure 2 shows the difference between the closest point to Harvest against 4, 9, and 16 closest points to Harvest "HarvestCalValAnalysis.png"
  -- figure 3 shows strictly the closest SLA values at Harvest compared to the closest 4 points to Harvest "CalVal1ptv4pt.png"
  -- figure 4 shows strictly the closest SLA values at Harvest compared to the closest 9 points to Harvest "CalVal1ptv9pt.png"
  -- figure 5 shows strictly the closest SLA values at Harvest compared to the closest 16 points to Harvest "CalVal1ptv16pt.png"
  -- figure 6 shows the compilation of all point averages when compared to the closest point to Harvest "CalValAnalysisAll.png"
  -- figure 7 shows the precision of n-point means by plotting n-mean vs. rms value "PrecisionOfAnalysis.png"
  -- figure 8 shows the exact features as figure 1 but instead of lidar data it shows NOAA lidar data "HarvestCalVal_NOAA.png"
  -- figure 9 shows Jason-3, lidar at Harvest, and NOAA lidar data comapred to one another "HarvestCalValAllDaa.png"
  -- figure 10 shows the attempt to extrapolate cdf data and compare it to the original dataset "HarvestCalValExtrapProof.png"
  -- figure 11 shows the flight path of Jason-3 along with its bins seperated by latitude "Jason3OverflightRoute.png"
  
  # minDist.m
  
- This function finds the distance between two points given the latitude and longitude of each points.
