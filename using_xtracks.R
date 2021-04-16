# Guide to using xtracks (TODO: make this an RMarkdown document)


# source the code file; change file path as needed.
# please note, that are a number of libraries used in the source file that you probably will need to install.
source("R/xtracks.R")

## creating an xtrack object

# load the raw data ; change file path as needed.
exd_1 <- read.csv("~/Dropbox/Hadza Data Dropbox/GIS Data and Maps/code/xtracks/exd_1.csv", stringsAsFactors=F)
exd_2 <- read.csv("~/Dropbox/Hadza Data Dropbox/GIS Data and Maps/code/xtracks/exd_2.csv", stringsAsFactors=F)
usethis::use_data(exd_1, exd_1)
usethis::use_data(exd_2, exd_2)

load("Data/d1.rda")
load("Data/d2.rda")

# lat, lon, elevation_m, in_camp, unix_time, distance_from_camp_m, utm_epsg

To construct an xtrack object, one must specify:
the lat, lon, elevation, in-camp status, time, distance from camp centroid, whether each
trackpoint is "in camp" or not, and the utm_epsg code.
lat and lon are expected to be in decimal-degree, WGS 84 format,
which is the default in most GPS devices mobile devices.
elevation is expected to be in meters above sealevel.
the "in_camp" parameter refers to whether each trackpoint is within or outside the boundaries
of a residential area, which in Wood et al. 2021 refers to the spatial boundaries of a Hadza camp;
but could more generally be considered the boundaries of a residential or habitation area, something
like a village or a camp, as appropriate in a given field setting. This is useful for indicating
travel for the purpose of aquiring resources -- AKA foraging travel, and needed for sinuosity measures.
Distance from camp centroid is expected to be the as-the-crow-flies distance
from the center of a residential area / 'camp' in meters (also needed for sinuosity measures).
The epsg code identifies the UTM zone of your study location.
This is needed for projecting lat / lon coordinates into UTM space.
To find the epsg code for your study location region of your track,
check out https://spatialreference.org/ref/epsg/

xt_1 <- xtrack(lat=d1$lat, lon=d1$lon, elevation_m=d1$elevation_m, in_camp=d1$in_camp, unix_time=d1$unix_time, distance_from_camp_m=d1$distance_from_camp_m, utm_epsg=32736)
xt_2 <- xtrack(lat=d2$lat, lon=d2$lon, elevation_m=d2$elevation_m, in_camp=d2$in_camp, unix_time=d2$unix_time, distance_from_camp_m=d2$distance_from_camp_m, utm_epsg=32736)

# *******************
# Analysis Functions
# *******************


#
# Segmentation of travel into bouts of out of camp travel.
#
# An out of camp bout is when an individual leaves camp, travels some distance,
# and then returns to camp.

out_of_camp_bout_records_1 <- xt_1$get_out_of_camp_bout_records()
out_of_camp_bout_records_2 <- xt_2$get_out_of_camp_bout_records()

head(out_of_camp_bout_records_1)


# How to get the total length of the xtrack in kilometers

xt_1$track_length_km

# how to get the total duration of the xtrack in hours
xt_1$track_duration_hr

# how to get the trackpoints of the longest duration bout
trackpoints_of_longest_duration_bout <- xt_1$get_longest_bout_trackpoints()
head(trackpoints_of_longest_duration_bout)

# how to test if the xtrack has data sufficient to enable sinuosity calculations
# of the manner carried out in Wood et al. 2021.

xt_1$has_data_for_sinuosity_measures()
xt_2$has_data_for_sinuosity_measures()

# how to get inbound and outbound sinuosity following the methods of Wood et al. 2021
xt_2$get_inbound_sinuosity()
xt_2$get_outbound_sinuosity()

# How to get more measures related to inbound and outbound sinuosity. These include:
# - The length (km) of the outbound and inbound segments as traveled,
# - The length (km) of the 'as the crow flies' distance from the
#   point of leaving camp to the most distant point (sp_distance_outbound_km).
# - The length of the 'as the crow flies' distance from the
#   the most distant point to the point of returning to camp (sp_distance_inbound_km).
# - The mean sinuosity of the inbound and outbound segments,
# - The trackpoint_id of the most distant (from camp centroid) trackpoint.

xt_2$get_sinuosity_measures()


#
# Raster analysis to categorize places on the landscape as visited or not visited.
# This is a binary raster representation of the xtrack,
# where cells that are visited are given value 1,
# and those not visited given value 0.
# The length and width of the raster cells in meters is determined by the parameter
# cell_size_m and is by default 10.
# This function accepts a parameter called selected_trackpoints which determines
# which of the trackpoints are rasterized. The acceptable values are
# all, in_camp, or out_of_camp. The default is all, as used in Wood et al. 2021.

bin_ras_1_10m <- xt_1$as_raster_of_habitat_visited_binary()
bin_ras_1_20m <- xt_1$as_raster_of_habitat_visited_binary(cell_size_m=20)
class(bin_ras_1_10m)
plot(bin_ras_1_10m)

# Raster analysis to categorize variable visitation intensity of places on the landscape
# This is a integer raster representation of the xtrack,
# where the count for each cell represents the number of trackpoints that fell within
# that cell's boundaries. Assuming that trackpoints are logged at regular time intervals,
# this raster provides a measure of the amount of time spent within each cell.
# un-visited cells are given value 0.
# As with the binary raster representation,
# The length and width of the raster cells in meters is determined by the parameter
# cell_size_m and is by default 10.

bin_ras_2 <- xt_2$as_raster_of_habitat_visited_counts()
class(bin_ras_2)

## Analysis of rates of habitat exploration

# Following Wood et al. 2021, this analysis computes the habitat visited / explored
# each day, the marginal 'new' habitat visited for each day, and the cumulative habitat explored across all days.
# In a real research application modeled on our paper, this analysis should be done with a temporally-sorted
# list of xtracks, with each xtrack representing one day of travel of the same person.
# In the list, the first day of data should be in position 1.
# For demonstration purposes, I construct below some dummy data that does not actually represent 8 days
# of real travel, but instead, segments 1 day into 8 partially overlapping (temporally and spatially) segments.
# Imagine however that xt_3 through xt_10 represent 8 days of travel; the code works the same.

dx <- read.csv("Data/d1.rda", stringsAsFactors=F)
d1 <- dx[1:1000,]
d2 <- dx[500:2000,]
d3 <- dx[1500:3000,]
d4 <- dx[2500:4000,]
d5 <- dx[3500:5000,]
d6 <- dx[4500:6000,]
d7 <- dx[5500:7000,]
d8 <- dx[6500:8260,]


## for building documentation
usethis::use_data(dx, dx)
usethis::use_data(d1, d1)
usethis::use_data(d2, d2)
usethis::use_data(d3, d3)
usethis::use_data(d4, d4)
usethis::use_data(d5, d5)
usethis::use_data(d6, d6)
usethis::use_data(d7, d7)
usethis::use_data(d8, d8)

xt_1 <- xtrack(lat=d1$lat, lon=d1$lon, elevation_m=d1$elevation_m, in_camp=d1$in_camp, unix_time=d1$unix_time, distance_from_camp_m=d1$distance_from_camp_m, utm_epsg=32736)
xt_2 <- xtrack(lat=d2$lat, lon=d2$lon, elevation_m=d2$elevation_m, in_camp=d2$in_camp, unix_time=d2$unix_time, distance_from_camp_m=d2$distance_from_camp_m, utm_epsg=32736)
xt_3 <- xtrack(lat=d3$lat, lon=d3$lon, elevation_m=d3$elevation_m, in_camp=d3$in_camp, unix_time=d3$unix_time, distance_from_camp_m=d3$distance_from_camp_m, utm_epsg=32736)
xt_4 <- xtrack(lat=d4$lat, lon=d4$lon, elevation_m=d4$elevation_m, in_camp=d4$in_camp, unix_time=d4$unix_time, distance_from_camp_m=d4$distance_from_camp_m, utm_epsg=32736)
xt_5 <- xtrack(lat=d5$lat, lon=d5$lon, elevation_m=d5$elevation_m, in_camp=d5$in_camp, unix_time=d5$unix_time, distance_from_camp_m=d5$distance_from_camp_m, utm_epsg=32736)
xt_6 <- xtrack(lat=d6$lat, lon=d6$lon, elevation_m=d6$elevation_m, in_camp=d6$in_camp, unix_time=d6$unix_time, distance_from_camp_m=d6$distance_from_camp_m, utm_epsg=32736)
xt_7 <- xtrack(lat=d7$lat, lon=d7$lon, elevation_m=d7$elevation_m, in_camp=d7$in_camp, unix_time=d7$unix_time, distance_from_camp_m=d7$distance_from_camp_m, utm_epsg=32736)
xt_8 <- xtrack(lat=d8$lat, lon=d8$lon, elevation_m=d8$elevation_m, in_camp=d8$in_camp, unix_time=d8$unix_time, distance_from_camp_m=d8$distance_from_camp_m, utm_epsg=32736)


list_of_xtracks <- list(xt_3, xt_4, xt_5, xt_6, xt_7, xt_8, xt_9, xt_10)

hab_exp_results <- get_hab_exp_across_days(list_of_xtracks, cell_size_m = 10)

#results are handed back in a hopefully-easy to understand data frame format.
hab_exp_results

# The plot below is similar to Figure 3 in Wood et al. 2021, though not as fancy.
# What is plotted is just one individual's travel across 'days'
# The y-axis here is also in units of square meters, not square kilometers.
plot(x=hab_exp_results$day, y=hab_exp_results$cum_sum_square_meters_visited_across_days, type="l", xlab="Day", ylab="Cummulative land explored (meters squared)", main="Habitat Explored Across Days")


## ************************************************
## Exporting xtracks to other file / object formats
## ************************************************

# data frame
# how to get all the trackpoints from an xtrack in dataframe format.
# This representation has more columns / more information and annotations
# that the 'raw' data used to construct an xtrack.
# This includes columns for the time between each trackpoint, meters traveled between trackpoints,
# speed of travel between trackpoints, and the utm coordinates of each trackpoint.
trackpoints_1 <- xt_1$get_trackpoints()
head(trackpoints_1)

# KML
# (KML files are used in Google Earth and elsewhere)
xt_1$write_kml_file(kml_file_name = "xt_1_c.kml")

# GPX
xt_1$write_gpx_file(gpx_file_name="xt_1.gpx", gpx_track_name="XT1")

# SpatialLinesDataFrame
# This is an object type in the sp package, an important R package for spatial analysis
xt_1_sldf <- xt_1$as_spatial_lines_dataframe()
class(xt_1_sldf)

## **********************
## Plotting options ##
## **********************

# Plot options provided by xtracks

# A call to plot_nice_map_of_track creates a 'nice map'
# that harnesses ggplot2 functions.
# It is a clean plot of the xtrack's travel path with
# a simple ggplot2 black and white theme,
# a customized scale bar, and some metadata displayed in the subtitle area.
# This function accepts parameters for a title (the_title),
# and the color of the line representing the xtrack (line_color).
xt_1$plot_nice_map_of_track()
xt_1$plot_nice_map_of_track(the_title="Day 249")
xt_1$plot_nice_map_of_track(the_title="Day 249", line_color="red")
xt_2$plot_nice_map_of_track(the_title="Day 250", line_color="blue")

# A call to plot_sinuosity_map creates a visual representation of the
# outbound travel segment (red), the inbound segment (blue),
# travel in camp or during shorter out of camp segments (green),
# the 'as the crow flies' shortest path segments used to
# calculate sinuosity measures (grey dashed line).
# The trackpoint that is maximally distant from the camp centroid is plotted
# in yellow.
# The sinuosity measures themselves are plotted in the subtitle and the
# plot accepts a parameter for the map title.

xt_2$plot_sinuosity_map(the_title="Sinuosity map, Day 1324")


# A call to as_mapview is a thin interface to the package mapview,
# producing a dynamic plot that harnesses the power of package mapview.
# this function accepts all parameters that can be fed to the
# function mapview in the package mapview, such as layer.name, color, etc.
xt_1$as_mapview(format="line", layer.name="Hot track", color="red")

# Using the power of mapview, we can display multiple xtracks on a dynamic map background
# as follows:
a <- xt_1$as_mapview(format="line", layer.name="Cool track", color="blue")
b <- xt_2$as_mapview(format="line", layer.name="Hot track", color="red")
a + b

# Plotting options harnessing the power of 'plot' functions defined in other packages.
# Below are simple examples of using the plot functions defined for
# SP and Raster; these plots can be sweetened by consulting those packages'
# documentation.

# A plot of the spatialLinesDataFrame representation of the xtrack
plot(xt_1$as_spatial_lines_dataframe())

# A plot of the binary raster representation of the xtrack
plot(xt_1$as_raster_of_habitat_visited_binary())

# A plot of the integer count raster representation of the xtrack
plot(xt_1$as_raster_of_habitat_visited_counts())


