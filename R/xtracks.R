#library(devtools)
#library(testthat)
#library(usethis)
#library(roxygen2)
#library(knitr)
#library(rmarkdown)
# if("sp" %in% rownames(installed.packages()) == FALSE) {install.packages("sp")}
# if("sf" %in% rownames(installed.packages()) == FALSE) {install.packages("sf")}
# if("geosphere" %in% rownames(installed.packages()) == FALSE) {install.packages("geosphere")}
# if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
# if("raster" %in% rownames(installed.packages()) == FALSE) {install.packages("raster")}
# if("ggsn" %in% rownames(installed.packages()) == FALSE) {install.packages("ggsn")}
# if("maptools" %in% rownames(installed.packages()) == FALSE) {install.packages("maptools")}
# if("mapview" %in% rownames(installed.packages()) == FALSE) {install.packages("mapview")}
#
 # library(sp)
 # library(sf)
 # library(geosphere)
 # library(ggplot2)
 # library(raster)
 # library(ggsn)
 # library(maptools)
 # library(mapview)


#' An xtrack object represents the movement of one individual on one day.
#'
#' @field trackpoints A dataframe of trackpoints
#' @field track_length_km The length of the trajectory in km
#' @field track_duration_hr the total duration of the xtrack in hours
#' @export xtrack
#' @exportClass xtrack
xtrack <- setRefClass("xtrack",
                     fields = list(pk_track_id="numeric", int_track_id="numeric", int_res_sec="numeric", person_id="character", age="numeric", sex="character",
                                   trackpoints="data.frame", out_of_camp_bout_records="data.frame",
                                   utm_epsg="numeric", utm_proj_args="character", sfc_linestring_object ="data.frame",
                                   camp="character", date="character", camp_lat="numeric", camp_lon="numeric",
                                   track_length_km="numeric", track_duration_hr="numeric",
                                   min_x_utm="numeric", max_x_utm="numeric", min_y_utm="numeric", max_y_utm="numeric",
                                   furthest_trackpoint_id="numeric", outbound_section_starting_trackpoint_id="numeric", inbound_section_ending_trackpoint_id="numeric",
                                   length_outbound_section_km="numeric",length_inbound_section_km="numeric",  sp_distance_outbound_km="numeric",
                                   sp_distance_inbound_km="numeric", outbound_sinuosity="numeric", inbound_sinuosity="numeric",mean_sinuosity="numeric", has_bout_appropriate_for_sinuosity_measures="numeric"
                     ),methods=list(
                       initialize=function(lat, lon, elevation_m, in_camp, unix_time, distance_from_camp_m, utm_epsg)
                       {
                         "Creates an xtrack object. To construct an xtrack object, one must specify: the lat, lon, elevation, in-camp status, time, distance from camp centroid, whether each trackpoint is \'in camp\' or not, and the utm_epsg code. lat and lon are expected to be in decimal-degree, WGS 84 format, which is the default in most GPS devices mobile devices. Elevation is expected to be in meters above sealevel. the \"in_camp\" parameter refers to whether each trackpoint is within or outside the boundaries of a residential area, which in Wood et al. 2021 refers to the spatial boundaries of a Hadza camp; but could more generally be considered the boundaries of a residential or habitation area, something like a village or a camp, as appropriate in a given field setting. This is useful for indicating travel for the purpose of aquiring resources -- AKA foraging travel, and needed for sinuosity measures.Distance from camp centroid is expected to be the as-the-crow-flies distance from the center of a residential area  or 'camp' in meters (also needed for sinuosity measures). The epsg code identifies the UTM zone of your study location. This is needed for projecting lat / lon coordinates into UTM space. To find the epsg code for your study location region of your track, check out <https://spatialreference.org/ref/epsg/>"
                         trackpoints <<- data.frame(lat=lat, lon=lon, elevation_m=elevation_m, in_camp=in_camp, unix_time=unix_time, distance_from_camp_m=distance_from_camp_m)
                         trackpoints <<- trackpoints[order(trackpoints$unix_time),]
                         trackpoints$pk_trackpoint_id <<- 1:nrow(trackpoints)

                         pts = matrix(0, length(lon), 2)
                         pts[,1] = lon
                         pts[,2] = lat
                         ls = sf::st_sfc(sf::st_linestring(pts), crs = 4326)
                         sfc_linestring_object <<- data.frame(ls)

                         track_length_km <<- round(sp::LineLength(pts, longlat=TRUE, sum=TRUE), 3)

                         #print(paste0("the length from LineLength: ", track_length_km))
                         trackpoints$seconds_since_prior_trackpoint <<- c(0,diff(trackpoints$unix_time))
                         trackpoints$meters_from_prior_trackpoint <<- 0
                         distances <- geoshere::distHaversine(p1=cbind(trackpoints$lon[1:(nrow(trackpoints)-1)], trackpoints$lat[1:(nrow(trackpoints)-1)]), p2=cbind(trackpoints$lon[2:nrow(trackpoints)], trackpoints$lat[2:nrow(trackpoints)]))
                         trackpoints$meters_from_prior_trackpoint <<- c(0,distances)
                         trackpoints$speed_m_s_from_prior_trackpoint <<- trackpoints$meters_from_prior_trackpoint/trackpoints$seconds_since_prior_trackpoint
                         trackpoints$speed_m_s_from_prior_trackpoint[1] <<- 0
                         #person_and_camp <- dbGetQuery(con_ti, paste("SELECT `fk_camp`, `fk_person` FROM `track` WHERE `pk_track_id`='",pk_track_id, "'", sep=""))
                         track_duration_hr <<- (trackpoints$unix_time[nrow(trackpoints)]-trackpoints$unix_time[1])/60/60
                         utm_epsg <<- utm_epsg
                         # this code is specific to where I work in Tanzania. Need to add this as a general function to translate data to UTM


                         #y = lat
                         #x = lon



                         tps <- data.frame(lon=trackpoints$lon, lat=trackpoints$lat)
                         lat_lon_args = "+proj=longlat +datum=WGS84 +ellps=WGS84"
                         utm_proj_args <<- paste0("+init=epsg:",utm_epsg)

                         #utm <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
                         coordinates(tps) = cbind("lon", "lat")
                         proj4string(tps) <- rgdal::CRS(lat_lon_args)
                         #lines <- SpatialLines(list(Lines(list(Line(tps)), "id")))
                         #proj4string(lines) <- CRS(lat_lon_args)
                         trackpoints_utm <- sp::spTransform(tps,  rgdal::CRS(utm_proj_args))
                          #head(coordinates(trackpoints_utm))
                         trackpoints$utm_x <<- sp::coordinates(trackpoints_utm)[,1]
                         trackpoints$utm_y <<- sp::coordinates(trackpoints_utm)[,2]
                         range_utm_x <- range(trackpoints$utm_x)
                         range_utm_y <- range(trackpoints$utm_y)
                         min_x_utm <<- range_utm_x[1]
                         max_x_utm <<- range_utm_x[2]
                         min_y_utm <<- range_utm_y[1]
                         max_y_utm <<- range_utm_y[2]


                         set_bout_records()
                         set_sinuosity_indices()

                         #print("finished set_sinuosity_indices")
                         #print("calling setBearingsAlongTrack")

                         #to do
                         #set_bearings_along_track()
                         #print("finished setBearingsAlongTrack")

                       },
                       as_mapview = function(format=c("line"), ...)
                       {
                         "# A call to as_mapview is a thin interface to the package mapview, producing a dynamic plot that harnesses the power of package mapview. This function accepts all parameters that can be fed to the function mapview in the package mapview, such as layer.name, color, etc."
                         if(format=="line")
                         {
                           mapview::mapview(.self$as_sfc_linestring(), ...=...)

                         } else if (format=="points")
                         {
                           return(mapview::mapview(x=data.frame(lon=trackpoints$lon, lat=trackpoints$lat), xcol=c("lon"), ycol=c("lat"), crs=c("WGS84"), ...=...))
                         }

                       },
                       # as_spatial_lines <- function()
                       # {
                       #     #t <- .self$as_sfc_linestring()
                       #     sfc_ls <- .self$as_sfc_linestring()
                       #     t2 <- sf::as_Spatial(from=sfc_ls, cast=FALSE)
                       #     return(t2)
                       # },
                       as_spatial_lines_dataframe=function()
                       {"Provides a SpatialLinesDataFrame representation of the track. This is an object type in the sp package, an important R package for spatial analysis."
                        the_spatial_lines <- sp::SpatialLines(list(Lines(Line(cbind(trackpoints$lon,trackpoints$lat)), ID="a")))
                        emptyData <- data.frame(matrix(0, ncol = 2, nrow = length(the_spatial_lines)))
                        the_spatialLinesDataFrame <- sp::SpatialLinesDataFrame(sl=the_spatial_lines, data=emptyData, match.ID=FALSE)
                        return(the_spatialLinesDataFrame)
                       },
                       as_raster_of_habitat_visited_binary=function(cell_size_m=10, xmin=NULL, xmax=NULL,ymin=NULL,ymax=NULL,selected_trackpoints="all")
                       {

                         "Raster analysis to categorize places on the landscape as visited or not visited. This is a binary raster representation of the xtrack, where cells that are visited / intersected are given value 1, and those not visited given value 0. The length and width of the raster cells in meters is determined by the parameter cell_size_m and is by default 10. This function accepts a parameter called selected_trackpoints which determines which of the trackpoints are rasterized. The acceptable values are all, in_camp, or out_of_camp. The default is all, as used in Wood et al. 2021."

                         the_extent <- NA


                         if(selected_trackpoints=="all")
                         {
                           selected_trackpoints <- trackpoints
                         } else if(selected_trackpoints=="in_camp")
                         {
                           selected_trackpoints <- trackpoints[trackpoints$in_camp==1,]
                         } else if(selected_trackpoints=="out_of_camp")
                         {
                           selected_trackpoints <- trackpoints[trackpoints$in_camp==0,]
                         } else {
                           print("value supplied for selected_trackpoints is invalid (should be 'all', 'in_camp', or 'out_of_camp')")
                           break
                         }

                         if(is.null(xmin)|is.null(xmax)|is.null(ymin)|is.null(ymax))
                         {
                           xmin <- min(selected_trackpoints$utm_x)
                           xmax <- max(selected_trackpoints$utm_x)
                           ymin <- min(selected_trackpoints$utm_y)
                           ymax <- max(selected_trackpoints$utm_y)
                           the_extent <- raster::extent(xmin, xmax, ymin, ymax)
                         } else {
                           the_extent <- raster::extent(xmin, xmax, ymin, ymax)
                         }

                          tp_r <- raster::raster(the_extent)
                          raster::res(tp_r) <- cell_size_m
                          rasterized_trackpoints <- raster::rasterize(x=selected_trackpoints[c("utm_x", "utm_y")], y=tp_r, fun='count')
                          cells_visited_today <- rasterized_trackpoints>0
                          cells_visited_today[is.na(cells_visited_today[])] <- 0

                          return(cells_visited_today)
                       },
                       as_raster_of_habitat_visited_counts=function(cell_size_m=10, xmin=NULL, xmax=NULL,ymin=NULL,ymax=NULL, selected_trackpoints="all")
                       {
                         "Raster analysis to categorize variable visitation intensity of places on the landscape. This is a integer raster representation of the xtrack, where the count for each cell represents the number of trackpoints that fell within that cell's boundaries. Assuming that trackpoints are logged at regular time intervals, this raster provides a measure of the amount of time spent within each cell. un-visited cells are given value 0. As with the binary raster representation, the length and width of the raster cells in meters is determined by the parameter cell_size_m and is by default 10."

                         the_extent <- NA

                         if(selected_trackpoints=="all")
                         {
                           selected_trackpoints <- trackpoints
                         } else if(selected_trackpoints=="in_camp")
                         {
                           selected_trackpoints <- trackpoints[trackpoints$in_camp==1,]
                         } else if(selected_trackpoints=="out_of_camp")
                         {
                           selected_trackpoints <- trackpoints[trackpoints$in_camp==0,]
                         } else {
                           print("value supplied for selected_trackpoints is invalid (should be 'all', 'in_camp', or 'out_of_camp')")
                           break
                         }

                         if(is.null(xmin)|is.null(xmax)|is.null(ymin)|is.null(ymax))
                         {
                           xmin <- min(selected_trackpoints$utm_x)
                           xmax <- max(selected_trackpoints$utm_x)
                           ymin <- min(selected_trackpoints$utm_y)
                           ymax <- max(selected_trackpoints$utm_y)
                           the_extent <- raster::extent(xmin, xmax, ymin, ymax)
                         } else {
                           the_extent <- raster::extent(xmin, xmax, ymin, ymax)
                         }

                         tp_r <- raster::raster(the_extent)
                         raster::res(tp_r) <- cell_size_m

                         rasterized_trackpoints <- raster::rasterize(x=selected_trackpoints[c("utm_x", "utm_y")], y=tp_r, fun='count')
                         rasterized_trackpoints[is.na(rasterized_trackpoints[])] <- 0
                         return(rasterized_trackpoints)
                       },
                       as_sfc_linestring=function()
                       {
                         # #print("as_sfc_linestring being called")
                         # y = trackpoints$lat
                         # x = trackpoints$lon
                         # pts = matrix(0, length(y), 2)
                         # pts[,1] = x
                         # pts[,2] = y
                         # ls = st_sfc(st_linestring(pts), crs = 4326)
                         #sfc_linestring_object[1,1]
                         return(.self$sfc_linestring_object[1,1])
                       },
                       get_trackpoints = function()
                       {"Returns a data frame of trackpoints. This representation has more columns / more information and annotations that the 'raw' data used to construct an xtrack. This includes columns for the time between each trackpoint, meters traveled between trackpoints, speed of travel between trackpoints, and the utm coordinates of each trackpoint."
                         return(trackpoints)
                       },
                       get_out_of_camp_bout_records = function()
                       {"Segmentation of travel into bouts of out of camp travel. An out of camp bout is when an individual leaves camp, travels some distance, and then returns to camp."
                         return(out_of_camp_bout_records)
                       },
                       get_longest_bout_trackpoints = function()
                       { "Get the trackpoints of the longest duration out of camp bout"


                         if(nrow(out_of_camp_bout_records)==0)
                         {
                           #print("there are no out of camp bouts")
                           return(out_of_camp_bout_records)
                         }

                         duration_longest_bout <- max(out_of_camp_bout_records$total_time_of_bout_seconds)
                         id_of_longest_bout <- out_of_camp_bout_records$bout_number[out_of_camp_bout_records$total_time_of_bout_seconds==duration_longest_bout]

                         start_time_lb <- out_of_camp_bout_records$time_leaving_camp[out_of_camp_bout_records$bout_number==id_of_longest_bout]
                         end_time_lb <- out_of_camp_bout_records$time_returning_to_camp[out_of_camp_bout_records$bout_number==id_of_longest_bout]

                         after_leave_time <- trackpoints$unix_time >= start_time_lb
                         before_return_time <- trackpoints$unix_time <= end_time_lb
                         during_longest_bout <- after_leave_time & before_return_time
                         trackpoints_of_lb <- trackpoints[during_longest_bout,]
                         return(trackpoints_of_lb)
                       },
                       get_outbound_sinuosity = function()
                       {"Get outbound sinuosity following the methods of Wood et al. 2021"
                         return(outbound_sinuosity)
                       },
                       get_inbound_sinuosity = function()
                       {"Get inbound sinuosity following the methods of Wood et al. 2021"
                         return(inbound_sinuosity)
                       },
                       get_sinuosity_measures=function()
                       {
                         "Get more measures related to inbound and outbound sinuosity. These include: The length (km) of the outbound and inbound segments as traveled; the length (km) of the 'as the crow flies' distance from the point of leaving camp to the most distant point (sp_distance_outbound_km).
; The length of the 'as the crow flies' distance from the most distant point to the point of returning to camp (sp_distance_inbound_km); The mean sinuosity of the inbound and outbound segments; The trackpoint_id of the most distant (from camp centroid) trackpoint."
                         return_value <- data.frame(cbind(length_outbound_section_km,length_inbound_section_km, sp_distance_outbound_km,sp_distance_inbound_km,outbound_sinuosity,inbound_sinuosity,mean_sinuosity, furthest_trackpoint_id))
                         return_value
                       },
                       plot_nice_map_of_track=function(the_title="", line_color="black")
                       {"A call to plot_nice_map_of_track creates a 'nice map' that harnesses ggplot2 functions. It is a clean plot of the xtrack's travel path with a simple ggplot2 black and white theme, a customized scale bar, and some metadata displayed in the subtitle area. This function accepts parameters for a title (the_title), and the color of the line representing the xtrack (line_color)."


                         dist_m <- track_length_km*1000#track_length_km <- ti$track_length_km
                         time_seconds <- track_duration_hr*60*60#track_duration_hr <- ti$track_duration_hr
                         average_speed_m_second <- round(dist_m/time_seconds,2)
                         #the_label <- paste("track_id:",  pk_track_id)
                         n_trackpoints <- nrow(trackpoints)
                         the_subtitle <- paste("Km traveled: ", round(track_length_km,2), ", Hours worn: ",round(track_duration_hr,1), "\nn trackpoints: ", n_trackpoints, sep="")
                         range_lon <- range(trackpoints$lon)
                         range_lat <- range(trackpoints$lat)
                         min_x <- range_lon[1]
                         max_x <- range_lon[2]
                         min_y <- range_lat[1]
                         max_y <- range_lat[2]

                         x_lim<-c(min_x, max_x)
                         y_lim<-c(min_y, max_y)

                         res <- equalizeAxesLimits(x_lim, y_lim)
                         x_lim <- res$x_lim
                         y_lim <- res$y_lim
                         new_dist_y_axis_km <- res$new_y_axis_distance_m/1000

                         #this creates a 10% buffer at the bottom of the plot area
                         #to permit a scale bar from being placed and not interfering with any
                         #points or lines on the map.

                         footer_fraction <- .1
                         footer_km <- new_dist_y_axis_km*footer_fraction
                         footer_m <- footer_km*1000
                         ll_with_footer <- geosphere::destPoint(p=c(x_lim[1], y_lim[1]), b=180, d=footer_m)
                         y_min_with_footer <- ll_with_footer[2]
                         y_lim[1] <- y_min_with_footer
                         res <- equalizeAxesLimits(x_lim, y_lim)
                         x_lim <- res$x_lim
                         y_lim <- res$y_lim

                         dist_y_axis_m_with_footer <- res$new_x_axis_distance_m
                         map_width_m <- res$new_x_axis_distance_m
                         map_width_km <- map_width_m/1000

                         right_nudge_distance_m <- map_width_m*.05
                         scale_lon_nudged_right <- geosphere::destPoint(p=c(x_lim[1], y_lim[1]), b=90, d=right_nudge_distance_m)[1]
                         scale_lon <- scale_lon_nudged_right
                         scale_lat <- y_lim[1]

                         new_dist_y_axis_km <- res$new_y_axis_distance_m/1000
                         one_third_map_width_km <- map_width_km/3
                         scale_width_km <- 0
                         scale_units <- "km"

                         if(one_third_map_width_km>1)
                         {
                           scale_width_km<-1
                         } else if(one_third_map_width_km>0.6) {
                           scale_width_km<-.25
                         } else if(one_third_map_width_km>0.2) {
                           scale_width_km<-.1
                         } else if(one_third_map_width_km>0.15) {
                           scale_width_km<-.15
                         } else if(one_third_map_width_km>0.10) {
                           scale_width_km<-.10
                         } else if(one_third_map_width_km>0.05) {
                           scale_width_km<-.05
                         } else if(one_third_map_width_km>0.01) {
                           scale_width_km<-.01
                         }


                         the_white_map <- ggplot(trackpoints) +
                           geom_path(aes(x=lon, y=lat), color=line_color, show.legend = TRUE, data=trackpoints) +
                           coord_equal() + xlim(x_lim) + ylim(y_lim) +
                           scale_bar(lon = scale_lon, lat = scale_lat,
                                     distance_lon = scale_width_km, distance_lat = new_dist_y_axis_km*.02, distance_legend = new_dist_y_axis_km*.05,
                                     dist_unit = "km", orientation = FALSE) +
                           labs(title=the_title, subtitle=the_subtitle, x="Longitude", y="Latitude", fill="") +
                           theme_bw(base_size = 12) +
                           theme(legend.title = element_blank()) +
                           scale_color_manual(values = c("red", "blue"), guide = guide_legend(override.aes = list(
                             linetype = c("solid", "blank"), shape = c(NA, 16))))


                           return(the_white_map)

                       },
                       plot_sinuosity_map=function(the_title="Sinuosity map")
                       {
                         "A call to plot_sinuosity_map creates a visual representation of the outbound travel segment (red), the inbound segment (blue),  travel in camp or during shorter out of camp segments (green), the 'as the crow flies' shortest path segments used to  calculate sinuosity measures (grey dashed line). The trackpoint that is maximally distant from the camp centroid is plotted  in yellow. The sinuosity measures themselves are plotted in the subtitle and the plot accepts a parameter for the map title."
                         #t<-track(3149)
                         #trackpoints <- t$trackpoints
                         #out_of_camp_bout_records <- t$out_of_camp_bout_records
                         #furthest_trackpoint_id <- t$furthest_trackpoint_id
                         #out_of_camp_bouts <- t$
                         #con_gsm <- opencon()

                         sin_rec <- get_sinuosity_measures()

                         min_x <- min(trackpoints$lon)
                         max_x <- max(trackpoints$lon)
                         min_y <- min(trackpoints$lat)
                         max_y <- max(trackpoints$lat)

                         x_lim<-c(min_x, max_x)
                         y_lim<-c(min_y, max_y)

                         res <- equalizeAxesLimits(x_lim, y_lim)
                         x_lim <- res$x_lim
                         y_lim <- res$y_lim

                         if(sin_rec$inbound_sinuosity==0 | sin_rec$outbound_sinuosity==0)
                         {
                           p <- ggplot() + geom_path(show.legend = TRUE, aes(x=trackpoints$lon, y=trackpoints$lat), color='yellow') +
                             coord_equal() + xlim(x_lim) + ylim(y_lim) +
                             theme(legend.title = element_blank()) +
                             labs(title=the_title, subtitle="There are no bouts out of camp to display", x="Longitude", y="Latitude", fill="") + theme_bw(base_size = 12)

                           return(p)
                         } else { #if there are sinuosity calculations to work with

                           trackpoints$path_segment <<- "In camp / other"
                           #trackpoints$path_segment <- "In camp / other"
                           trackpoints$path_segment[trackpoints$is_in_outbound_segment]<<-"Outbound segment"
                           #trackpoints$path_segment[trackpoints$is_in_outbound_segment]<-"Outbound segment"
                           trackpoints$path_segment[trackpoints$is_in_inbound_segment]<<-"Inbound segment"
                           #trackpoints$path_segment[trackpoints$is_in_inbound_segment]<-"Inbound segment"


                           sorted_bouts <- out_of_camp_bout_records[order(-out_of_camp_bout_records$total_length_of_bout_m),]
                           most_distant_bout <- sorted_bouts[1,]

                           most_distant_lon <- most_distant_bout$max_dist_during_bout_lon
                           most_distant_lat <- most_distant_bout$max_dist_during_bout_lat

                           start_outbound_lon <- trackpoints$lon[trackpoints$unix_time==most_distant_bout$time_leaving_camp]
                           start_outbound_lat <- trackpoints$lat[trackpoints$unix_time==most_distant_bout$time_leaving_camp]
                           end_inbound_lon <- trackpoints$lon[trackpoints$unix_time==most_distant_bout$time_returning_to_camp]
                           end_inbound_lat <- trackpoints$lat[trackpoints$unix_time==most_distant_bout$time_returning_to_camp]

                           longest_bout <- get_longest_bout_trackpoints()
                           outbound_sp <- data.frame(x=c(start_outbound_lon, most_distant_lon), y=c(start_outbound_lat, most_distant_lat))
                           inbound_sp <- data.frame(x=c(most_distant_lon, end_inbound_lon), y=c(most_distant_lat, end_inbound_lat))
                           outbound_section <- longest_bout[longest_bout$is_in_outbound_segment,]
                           inbound_section <- longest_bout[longest_bout$is_in_inbound_segment,]
                           outbound_sp$spath_segment<-"outbound shortest path"
                           inbound_sp$spath_segment<-"inbound shortest path"

                           the_subtitle<- paste("Outbound sinuosity:", round(outbound_sinuosity,1), "\ninbound sinuosity:", round(inbound_sinuosity,1))
                           #trackpoints <- t4235$trackpoints

                           trackpoints$path_segment <<- factor(trackpoints$path_segment, levels=c("Outbound segment", "Inbound segment", "In camp / other"), labels=c("Outbound segment", "Inbound segment", "In camp / other"))
                           trackpoints$`Track segment` <<- trackpoints$path_segment
                           p <- ggplot(trackpoints, aes(lon, lat, color=`Track segment`))+ geom_path() +
                             scale_colour_manual(values=c("red","blue","darkgreen")) +
                             coord_equal() + xlim(x_lim) + ylim(y_lim) +
                             theme(legend.title = element_blank()) +
                             geom_path(data=outbound_sp, aes(x=x, y=y, color=spath_segment), linetype=2, color='grey') +
                             geom_path(data=inbound_sp, aes(x=x, y=y, color=spath_segment), linetype=2, color='grey') +
                             geom_point(aes(x=most_distant_lon, y=most_distant_lat), color='yellow') +
                             labs(title=paste(the_title), subtitle=the_subtitle, x="Longitude", y="Latitude", fill="") + theme_bw(base_size = 12)
                           return(p)
                         }

                       },
                       scale_bar = function(lon, lat, distance_lon, distance_lat, distance_legend, dist_unit = "km", rec_fill = "white", rec_colour = "black", rec2_fill = "black", rec2_colour = "black", legend_colour = "black", legend_size = 3, orientation = TRUE, arrow_length = 500, arrow_distance = 300, arrow_north_size = 6)
                       {
                         the_scale_bar <- create_scale_bar(lon = lon, lat = lat, distance_lon = distance_lon, distance_lat = distance_lat, distance_legend = distance_legend, dist_unit = dist_unit)
                         #the_scale_bar$legend
                         # First rectangle
                         rectangle1 <- geom_polygon(data = the_scale_bar$rectangle, aes(x = lon, y = lat), fill = rec_fill, colour = rec_colour)

                         # Second rectangle
                         rectangle2 <- geom_polygon(data = the_scale_bar$rectangle2, aes(x = lon, y = lat), fill = rec2_fill, colour = rec2_colour)

                         # Legend
                         scale_bar_legend <- annotate("text", label = paste(the_scale_bar$legend[,"text"], dist_unit, sep=""), x = the_scale_bar$legend[,"long"], y = the_scale_bar$legend[,"lat"], size = legend_size, colour = legend_colour)

                         res <- list(rectangle1, rectangle2, scale_bar_legend)

                         if(orientation){# Add an arrow pointing North
                           coords_arrow <- create_orientation_arrow(scale_bar = the_scale_bar, length = arrow_length, distance = arrow_distance, dist_unit = dist_unit)
                           arrow <- list(geom_segment(data = coords_arrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coords_arrow$coords_n[1,"x"], y = coords_arrow$coords_n[1,"y"], size = arrow_north_size, colour = "black"))
                           res <- c(res, arrow)
                         }
                         return(res)
                       },
                       create_scale_bar = function(lon,lat,distance_lon,distance_lat,distance_legend, dist_units = "km")
                       {
                         # First rectangle
                         bottom_right <- maptools::gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon, dist.units = dist_units, model = "WGS84")

                         topLeft <- maptools::gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_lat, dist.units = dist_units, model = "WGS84")
                         rectangle <- cbind(lon=c(lon, lon, bottom_right[1,"long"], bottom_right[1,"long"], lon),
                                            lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
                         rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)

                         # Second rectangle t right of the first rectangle
                         bottom_right2 <- maptools::gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon*2, dist.units = dist_units, model = "WGS84")
                         rectangle2 <- cbind(lon = c(bottom_right[1,"long"], bottom_right[1,"long"], bottom_right2[1,"long"], bottom_right2[1,"long"], bottom_right[1,"long"]),
                                             lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
                         rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)

                         # Now let's deal with the text
                         on_top <- maptools::gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_legend, dist.units = dist_units, model = "WGS84")
                         on_top2 <- on_top3 <- on_top
                         on_top2[1,"long"] <- bottom_right[1,"long"]
                         on_top3[1,"long"] <- bottom_right2[1,"long"]

                         legend <- rbind(on_top, on_top2, on_top3)
                         legend <- data.frame(cbind(legend, text = c(0, distance_lon, distance_lon*2)), stringsAsFactors = FALSE, row.names = NULL)
                         return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
                       },
                       calc_and_save_raster_of_habitat_visited_binary=function(cell_size_m=10, interpolated=TRUE, int_time_res=5)
                       {
                         #this function calculates and saves a raster of habvis binary.

                         #cell_size_m=10
                         #interpolated=TRUE
                         #int_time_res=5
                         #t <- track(4416, interpolated=TRUE, int_time_res=5)
                         #camp <- t$camp
                         #pk_track_id <- t$pk_track_id
                         #trackpoints <- t$trackpoints
                         #constructs a SpatialLines object from the trackpoints,
                         #then overlays this on a raster background, and each cell
                         #that is intersected > 0 times receives a value of 1, otherwise NA

                         file_name <- NA
                         if(interpolated)
                         {
                           file_name <- paste("track_", pk_track_id, "_", cell_size_m, "_m_raster_of_visitation_binary_interpolated_", int_time_res, ".rds", sep="")

                         } else {
                           file_name <- paste("track_", pk_track_id, "_", cell_size_m, "_m_raster_of_visitation_binary.rds", sep="")
                         }

                         setwd("~/Dropbox/Manuscripts/Age and sex differences in Hadza hunter-gatherer space use/data/rasterized_tracks/")


                         camp_raster <- getRasterByCamp(camp, cell_size_m)
                         #camp_raster <- getRasterByCamp("hukumako_", cell_size_m)
                         extent_of_camp_raster <- extent(camp_raster)



                         if(file.exists(file_name))
                         {
                           r <- readRDS(file=file_name)
                           extent_of_saved_raster <- extent(r)
                           extent_is_equal_to_that_of_camp <- extent_of_camp_raster==extent_of_saved_raster

                           if(extent_is_equal_to_that_of_camp)
                           {
                             print("this habvis has already been calculated, but making another now")
                             camp_raster <- getRasterByCamp(camp, cell_size_m)
                             #camp_raster <- getRasterByCamp("hukumako_", cell_size_m)
                             cr_e <- extent(camp_raster)
                             tps <- data.frame(lon=trackpoints$lon, lat=trackpoints$lat)

                             latlong = "+proj=longlat +datum=WGS84 +ellps=WGS84"
                             utm <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
                             coordinates(tps) = cbind("lon", "lat")
                             proj4string(tps) <- CRS(latlong)
                             lines <- SpatialLines(list(Lines(list(Line(tps)), "id")))
                             proj4string(lines) <- CRS(latlong)
                             trackpoints_utm <- spTransform(lines,  CRS(utm))
                             #class(trackpoints_utm)

                             tp_e <- extent(trackpoints_utm)
                             #plot(trackpoints_utm, add=TRUE)
                             tp_r <- raster(cr_e)
                             res(tp_r) <- cell_size_m
                             #plot(trackpoints_utm)
                             rasterized_trackpoints <- rasterize(x=trackpoints_utm, y=tp_r, fun='count')
                             cells_visited_today <- rasterized_trackpoints>0
                             saveRDS(cells_visited_today, file=file_name)
                           } else {
                             print("saved raster version of that track found, but wrong extent, so making a new one")
                             camp_raster <- getRasterByCamp(camp, cell_size_m)
                             #camp_raster <- getRasterByCamp("hukumako_", cell_size_m)
                             cr_e <- extent(camp_raster)
                             #trackpoints <- t4424$trackpoints
                             tps <- data.frame(lon=trackpoints$lon, lat=trackpoints$lat)

                             latlong = "+proj=longlat +datum=WGS84 +ellps=WGS84"
                             utm <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
                             coordinates(tps) = cbind("lon", "lat")
                             proj4string(tps) <- CRS(latlong)
                             lines <- SpatialLines(list(Lines(list(Line(tps)), "id")))
                             proj4string(lines) <- CRS(latlong)
                             trackpoints_utm <- spTransform(lines,  CRS(utm))
                             tp_e <- extent(trackpoints_utm)
                             #plot(trackpoints_utm, add=TRUE)
                             tp_r <- raster(cr_e)


                             res(tp_r) <- cell_size_m

                             rasterized_trackpoints <- rasterize(x=trackpoints_utm, y=tp_r, fun='count')
                             extent(rasterized_trackpoints)
                             cells_visited_today <- rasterized_trackpoints>0
                             saveRDS(cells_visited_today, file=file_name)
                             #return(cells_visited_today)
                           }
                         } else {
                           print("no saved raster version of that track found, so creating one")

                           camp_raster <- getRasterByCamp(camp, cell_size_m)
                           #camp_raster <- getRasterByCamp("hukumako_", cell_size_m)
                           cr_e <- extent(camp_raster)
                           tps <- data.frame(lon=trackpoints$lon, lat=trackpoints$lat)

                           latlong = "+proj=longlat +datum=WGS84 +ellps=WGS84"
                           utm <- "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
                           coordinates(tps) = cbind("lon", "lat")
                           proj4string(tps) <- CRS(latlong)
                           lines <- SpatialLines(list(Lines(list(Line(tps)), "id")))
                           proj4string(lines) <- CRS(latlong)
                           trackpoints_utm <- spTransform(lines,  CRS(utm))
                           #class(trackpoints_utm)

                           tp_e <- extent(trackpoints_utm)
                           #plot(trackpoints_utm, add=TRUE)
                           tp_r <- raster(cr_e)
                           res(tp_r) <- cell_size_m
                           #plot(trackpoints_utm)
                           rasterized_trackpoints <- rasterize(x=trackpoints_utm, y=tp_r, fun='count')
                           cells_visited_today <- rasterized_trackpoints>0
                           saveRDS(cells_visited_today, file=file_name)
                           #return(cells_visited_today)

                         }
                       },
                       has_data_for_sinuosity_measures = function()
                       {"Test if the xtrack has data sufficient to enable sinuosity calculations of the manner carried out in Wood et al. 2021."
                         return(has_bout_appropriate_for_sinuosity_measures==1)
                       },
                       set_bout_records = function()
                       {
                         #for testing
                         # d <- read.csv("data/track_data.csv", stringsAsFactors=F)
                         # xt_1 <- xtrack(lat=d$lat, lon=d$lon, elevation_m=d$elevation_m, in_camp=d$in_camp, unix_time=d$unix_time, distance_from_camp_m=d$distance_from_camp_m)
                         # trackpoints <- xt_1$trackpoints
                         # trackpoints$bout_number <- "in camp"

                         trackpoints$bout_number <<- "in camp"


                         bout_records <- data.frame(time_leaving_camp=numeric(100),
                                                    time_returning_to_camp=numeric(100),
                                                    total_time_of_bout_seconds=numeric(100),
                                                    bout_identified=numeric(100),
                                                    trackpoint_id_max_dist_from_camp_during_bout=numeric(100),
                                                    max_dist_from_camp_during_bout_m=numeric(100),
                                                    mean_dist_from_camp_during_bout_m=numeric(100),
                                                    bout_number=numeric(100),
                                                    max_dist_during_bout_lon=numeric(100),
                                                    max_dist_during_bout_lat=numeric(100),
                                                    max_dist_unix_time=numeric(100),
                                                    total_length_of_bout_m=numeric(100),
                                                    stringsAsFactors=FALSE)
#
#                          bout_records$track_id <- pk_track_id
#                          #bout_records$track_id <- t$pk_track_id
#
#                          bout_records$person_id <- person_id
#                          #bout_records$person_id <- t$person_id
#
#                          bout_records$age <- age
#                          #bout_records$age <- t$age
#
#                          bout_records$sex <- sex
#                          #bout_records$sex <- t$sex
#
#                          bout_records$camp <- camp
#                          #bout_records$camp <- t$camp
#
#                          bout_records$date <- date
#                          #bout_records$date <- t$date
#                          #trackpoints <- suppressWarnings(dbGetQuery(con, paste("SELECT lat, lon, distance_from_camp_m, unix_time, in_camp, pk_trackpoint_id FROM trackpoint WHERE fk_track_id='", track_id, "' ORDER BY unix_time", sep="")))
#                          #head(trackpoints)
#
#                          #table(trackpoints$in_camp)

                         bout_record_index <- 1

                         #this code identifies each time the track 'leaves camp'
                         for(x in 2:nrow(trackpoints))
                         {
                           #print(x)
                           in_camp_time_n_minus_1 <- trackpoints$in_camp[x-1]
                           in_camp_time_n <- trackpoints$in_camp[x]


                           if(in_camp_time_n_minus_1==1 & in_camp_time_n==0)
                           {
                             #print(paste("left camp at time", trackpoints$unix_time[x], "inserted in row: ", bout_record_index, "of bout records"))
                             bout_records$time_leaving_camp[bout_record_index] <- trackpoints$unix_time[x]
                             bout_record_index <- bout_record_index+1
                           }

                         }

                         bout_records <- bout_records[bout_records$time_leaving_camp!=0,]
                         #print(paste("there are", nrow(bout_records), "bout records to scan"))

                         if(nrow(bout_records)==0)
                         {
                           #out_of_camp_bout_records <<- NULL
                           print("Warning: In this track, there were no bouts of out of camp travel by which to calculate sinuosity of out-of-camp travel")
                           return()
                         }

                         #this code figures out, for each time the track leaves camp,
                         # the next time in which the track re-enters camp.
                         for(i in 1:nrow(bout_records))
                         {

                           trackpoints_to_scan <- trackpoints[(trackpoints$unix_time > bout_records$time_leaving_camp[i]),]
                           #if they left camp on the very last trackpoint of a track file

                           if(nrow(trackpoints_to_scan)==0)
                           {
                             bout_records$time_returning_to_camp[i]<- bout_records$time_leaving_camp[i]
                             bout_records$total_time_of_bout_seconds[i] <- 0
                             bout_records$bout_identified[i] <- 0
                             break
                           }

                           for(j in 1:nrow(trackpoints_to_scan))
                           {
                             #print(j)
                             in_camp_status_now <- trackpoints_to_scan$in_camp[j]

                             #if(is.na(in_camp_status_now)) {in_camp_status_now<-0}

                             if(in_camp_status_now==1)
                             {
                               #print(paste("returned to camp at time:", trackpoints_to_scan$unix_time[j], "inserted in row", i, "of bout records" ))
                               bout_records$time_returning_to_camp[i]<- trackpoints_to_scan$unix_time[j]
                               bout_records$total_time_of_bout_seconds[i] <- bout_records$time_returning_to_camp[i] - bout_records$time_leaving_camp[i]
                               bout_records$bout_identified[i] <- 1
                               break
                             }
                           }
                         }


                         bout_records <- bout_records[bout_records$bout_identified==1,]

                         if(nrow(bout_records)==0)
                         {
                           #out_of_camp_bout_records <<- NULL
                           return()
                         }

                         bout_records$bout_number <- 1:nrow(bout_records)

                         for(b in 1:nrow(bout_records))
                         {

                           time_left <- bout_records$time_leaving_camp[b]
                           time_ret <- bout_records$time_returning_to_camp[b]

                           trackpoints_of_bout <- trackpoints[(trackpoints$unix_time>=time_left&trackpoints$unix_time<=time_ret),]
                           length_of_bout_m <- sum(trackpoints_of_bout$meters_from_prior_trackpoint)-trackpoints_of_bout$meters_from_prior_trackpoint[1]
                           max_dist_from_camp_during_bout_m <- max(trackpoints_of_bout$distance_from_camp_m)
                           trackpoint_id_max_dist <- trackpoints_of_bout$pk_trackpoint_id[trackpoints_of_bout$distance_from_camp_m==max_dist_from_camp_during_bout_m][1]
                           max_dist_lon_during_bout <- trackpoints_of_bout$lon[trackpoints_of_bout$pk_trackpoint_id==trackpoint_id_max_dist]
                           max_dist_lat_during_bout <- trackpoints_of_bout$lat[trackpoints_of_bout$pk_trackpoint_id==trackpoint_id_max_dist]
                           max_dist_unix_time <- trackpoints_of_bout$unix_time[trackpoints_of_bout$pk_trackpoint_id==trackpoint_id_max_dist]
                           mean_dist_during_bout_m <- mean(trackpoints_of_bout$distance_from_camp_m)

                           bout_records$trackpoint_id_max_dist_from_camp_during_bout[b]<-trackpoint_id_max_dist
                           bout_records$max_dist_from_camp_during_bout_m[b]<-max_dist_from_camp_during_bout_m
                           bout_records$max_dist_during_bout_lon[b]<-max_dist_lon_during_bout
                           bout_records$max_dist_during_bout_lat[b]<-max_dist_lat_during_bout
                           bout_records$max_dist_unix_time[b]<-max_dist_unix_time
                           bout_records$mean_dist_from_camp_during_bout_m[b]<-mean_dist_during_bout_m
                           bout_records$total_length_of_bout_m[b]<-length_of_bout_m
                         }

                         out_of_camp_bout_records <<- bout_records

                         n_bouts <- nrow(out_of_camp_bout_records)

                         for(i in 1:n_bouts)
                         {
                           s_time <- out_of_camp_bout_records$time_leaving_camp[i]
                           e_time <- out_of_camp_bout_records$time_returning_to_camp[i]
                           distance_m <- round(out_of_camp_bout_records$total_length_of_bout_m[i],0)
                           trackpoints$bout_number[trackpoints$unix_time>=s_time&trackpoints$unix_time<=e_time]<<-as.character(i)
                         }
                       },
                      set_sinuosity_indices=function()
                      {
                        #print("in set sinuosity indices")
                        # true at initialization, may be set false if this track does not have a bout long enough to qualify for sinuosity measures
                        has_bout_appropriate_for_sinuosity_measures <<- 1
                        furthest_trackpoint_id <<- trackpoints$pk_trackpoint_id[trackpoints$distance_from_camp_m==max(trackpoints$distance_from_camp_m)][1]




                        #pp("setting sinuosity indices for track", pk_track_id)
                        #con<-getHadzaspatialRemoteCon()
                        #test_track<-track(4090)
                        #trackpoints <- test_track$trackpoints

                        trackpoints_of_longest_bout <- get_longest_bout_trackpoints()
                        #pk_track_id <- test_track$pk_track_id


                        # print(paste("setting sinuousity indices for track", pk_track_id))

                        if(nrow(trackpoints_of_longest_bout)==0)
                        {
                          has_bout_appropriate_for_sinuosity_measures <<- 0
                          length_outbound_section_km<<-0
                          length_inbound_section_km<<-0
                          sp_distance_outbound_km<<-0
                          sp_distance_inbound_km<<-0
                          outbound_sinuosity<<-0
                          inbound_sinuosity<<-0
                          mean_sinuosity<<-0
                          outbound_section_starting_trackpoint_id<<- -99
                          inbound_section_ending_trackpoint_id<<- -99
                          return()
                        }



                        #out_of_camp_bout_records <- xt_1$get_out_of_camp_bout_records()
                        if(nrow(out_of_camp_bout_records)>0)
                        {
                          #this means there is at least one bout record
                          #trackpoints_of_longest_bout <- t1$getLongestDistanceBoutOnDay()


                          sorted_bout_records <- out_of_camp_bout_records[order(-out_of_camp_bout_records$total_length_of_bout_m),]
                          stats_on_longest_distance_bout <- sorted_bout_records[1,]
                          max_dist_on_longest_dist_bout <- stats_on_longest_distance_bout$max_dist_from_camp_during_bout_m
                          total_length_of_bout <- stats_on_longest_distance_bout$total_length_of_bout_m

                          if(max_dist_on_longest_dist_bout < 500 | total_length_of_bout < 500)
                          {
                            #this means the longest distance bout is too short to consider for calculations,
                            #following Raichlen et al. 2014, we only consider bouts where the
                            #subject has traveled at least 500 m and has gone 500 m from camp.


                            has_bout_appropriate_for_sinuosity_measures <<- 0
                            length_outbound_section_km<<-0
                            length_inbound_section_km<<-0
                            sp_distance_outbound_km<<-0
                            sp_distance_inbound_km<<-0
                            outbound_sinuosity<<-0
                            inbound_sinuosity<<-0
                            mean_sinuosity<<-0
                            outbound_section_starting_trackpoint_id<<- -99
                            inbound_section_ending_trackpoint_id<<- -99
                          } else {
                            #this means there is a bout, and that it is long enough to
                            #merit calculation of sinuosity measures


                            start_trackpoint_of_bout <- trackpoints_of_longest_bout[1,]
                            end_trackpoint_of_bout <- trackpoints_of_longest_bout[nrow(trackpoints_of_longest_bout),]
                            #here
                            outbound_section_starting_trackpoint_id<<-start_trackpoint_of_bout$pk_trackpoint_id
                            inbound_section_ending_trackpoint_id<<-end_trackpoint_of_bout$pk_trackpoint_id
                            #outbound_section_starting_trackpoint_id<-start_trackpoint_of_bout$pk_trackpoint_id
                            #inbound_section_ending_trackpoint_id<-end_trackpoint_of_bout$pk_trackpoint_id


                            furthest_trackpoint <- trackpoints_of_longest_bout[trackpoints_of_longest_bout$unix_time==stats_on_longest_distance_bout$max_dist_unix_time,]


                            furthest_trackpoint_id <<- furthest_trackpoint$pk_trackpoint_id
                            furthest_trackpoint_unix_time <- furthest_trackpoint$unix_time
                            #furthest_trackpoint_id <- furthest_trackpoint$pk_trackpoint_id

                            outbound_section <- trackpoints_of_longest_bout[(trackpoints_of_longest_bout$unix_time < furthest_trackpoint$unix_time),]
                            inbound_section <- trackpoints_of_longest_bout[(trackpoints_of_longest_bout$unix_time >= furthest_trackpoint$unix_time),]

                            #length(outbound_section)

                            trackpoints$is_in_outbound_segment <<- trackpoints$pk_trackpoint_id %in% outbound_section$pk_trackpoint_id
                            #trackpoints$is_in_outbound_segment <- trackpoints$pk_trackpoint_id %in% outbound_section$pk_trackpoint_id
                            #table(trackpoints$is_in_outbound_segment)

                            trackpoints$is_in_inbound_segment <<- trackpoints$pk_trackpoint_id %in% inbound_section$pk_trackpoint_id
                            #trackpoints$is_in_inbound_segment <- trackpoints$pk_trackpoint_id %in% inbound_section$pk_trackpoint_id

                            n_trackpoints_outbound_segment <- sum(trackpoints$is_in_outbound_segment)
                            n_trackpoints_inbound_segment <- sum(trackpoints$is_in_inbound_segment)


                            #if there is no segment to be considered inbound or outbound --
                            #this can happen with very short trips out of camp.
                            if(n_trackpoints_outbound_segment<=1 | n_trackpoints_inbound_segment<=1)
                            {

                              #pp("outbound or inbound segment less than 1 trackpoint long, so no sinuosity calculated")
                              has_bout_appropriate_for_sinuosity_measures <<- 0
                              length_outbound_section_km<<-0
                              length_inbound_section_km<<-0
                              sp_distance_outbound_km <<- 0
                              sp_distance_inbound_km <<- 0
                              outbound_sinuosity <<- 0
                              inbound_sinuosity <<- 0
                              mean_sinuosity <<- 0
                              outbound_section_starting_trackpoint_id<<- -99
                              inbound_section_ending_trackpoint_id<<- -99

                            } else{
                              #pp("sinuosity calculated successfully")
                              length_outbound_section_km <<- round(LineLength(as.matrix(cbind(outbound_section$lon, outbound_section$lat)), longlat=TRUE, sum=TRUE), 2)
                              length_inbound_section_km <<- round(LineLength(as.matrix(cbind(inbound_section$lon, inbound_section$lat)), longlat=TRUE, sum=TRUE), 2)
                              start_trackpoint_lat <- start_trackpoint_of_bout$lat
                              start_trackpoint_lon <- start_trackpoint_of_bout$lon
                              furthest_trackpoint_lat <- furthest_trackpoint$lat
                              furthest_trackpoint_lon <- furthest_trackpoint$lon
                              end_trackpoint_lat <- end_trackpoint_of_bout$lat
                              end_trackpoint_lon <- end_trackpoint_of_bout$lon
                              start_coords <- as.matrix(cbind(lon=start_trackpoint_of_bout$lon, lat=start_trackpoint_of_bout$lat))
                              furthest_coords <- as.matrix(cbind(furthest_trackpoint$lon, lat=furthest_trackpoint$lat))
                              end_coords <- as.matrix(cbind(lon=end_trackpoint_of_bout$lon, lat=end_trackpoint_of_bout$lat))
                              sp_distance_outbound_km <<- spDistsN1(start_coords,furthest_coords, longlat=TRUE)
                              sp_distance_inbound_km <<- spDistsN1(furthest_coords,end_coords, longlat=TRUE)
                              outbound_sinuosity <<- length_outbound_section_km/sp_distance_outbound_km
                              inbound_sinuosity <<- length_inbound_section_km/sp_distance_inbound_km
                              mean_sinuosity <<- mean(c(outbound_sinuosity,inbound_sinuosity))
                            }
                          }




                        } else {
                          #this means there are no out of camp bout records ...
                          has_bout_appropriate_for_sinuosity_measures <<- 0
                          length_outbound_section_km<<-0
                          length_inbound_section_km<<-0
                          sp_distance_outbound_km<<-0
                          sp_distance_inbound_km<<-0
                          outbound_sinuosity<<-0
                          inbound_sinuosity<<-0
                          mean_sinuosity<<-0
                          outbound_section_starting_trackpoint_id<<- -99
                          inbound_section_ending_trackpoint_id<<- -99
                        }

                        #print(paste("in the longest distance bout there are:", nrow(trackpoints_of_longest_bout), "trackpoints"))





                      },
                      set_habvis_binary_raster=function(habvis_raster)
                      {

                      },
                      write_gpx_file=function(gpx_file_name="test.gpx", gpx_track_name="test")
                      { "writes a GPX file representation of the track. GPX files are widely used in GIS and GPS applications."
                        #trackpoints <- trackpoints_2
                        elevation = trackpoints$elevation_m
                        lat = trackpoints$lat
                        lon = trackpoints$lon
                        utc_datetime = as.POSIXlt(trackpoints$unix_time, origin="1970-01-01", tz="UTC")
                        utc_datetime_format="%Y-%m-%d %H:%M:%S"
                        utc_datetime_formatted <- format(utc_datetime, format=utc_datetime_format)

                        # gpx_track_name="test_track"
                        # gpx_file_name="noname.gpx"

                        iso_8601_date_time <- format(strptime(utc_datetime, format=utc_datetime_format, tz = "UTC"), format="%Y-%m-%dT%H:%M:%SZ")
                        gpx_open <- "<gpx><metadata></metadata>"#if the metadata node is not present, then this gpx file
                        #will not be able to be read using read_xml -- so that is why it is there (after lots of headaches testing)
                        track_open <- paste("<trk><name>", gpx_track_name, "</name>", sep="")
                        trackseg_open <- "<trkseg>"
                        the_trackpoints<-paste("<trkpt lat=\"",lat, "\" lon=\"", lon, "\"><ele>", elevation, "</ele><time>",iso_8601_date_time,"</time></trkpt>",sep="")
                        trackpoints_collapsed <- paste(the_trackpoints, collapse="", sep="")
                        #cat(trackpoints_collapsed)
                        trackseg_close <- "</trkseg>"
                        track_close <- "</trk>"
                        gpx_close <- "</gpx>"
                        all_gpx <- paste(gpx_open, track_open, trackseg_open, trackpoints_collapsed,trackseg_close,track_close,gpx_close, collapse="", sep="")

                        #write_file(x = all_gpx, path = "test_lat_lon_with_quotes.gpx")

                        # all_gpx_with_metadata<-paste("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?><gpx xmlns=\"http://www.topografix.com/GPX/1/1\" xmlns:gpxx=\"http://www.garmin.com/xmlschemas/WaypointExtension/v1\" xmlns:gpxtrx=\"http://www.garmin.com/xmlschemas/GpxExtensions/v3\" xmlns:gpxtpx=\"http://www.garmin.com/xmlschemas/TrackPointExtension/v1\" creator=\"Oregon 550t\" version=\"1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.topografix.com/GPX/1/1 http://www.topografix.com/GPX/1/1/gpx.xsd http://www.garmin.com/xmlschemas/WaypointExtension/v1 http://www8.garmin.com/xmlschemas/WaypointExtensionv1.xsd http://www.garmin.com/xmlschemas/TrackPointExtension/v1 http://www.garmin.com/xmlschemas/TrackPointExtensionv1.xsd\"><metadata><link href=\"http://www.garmin.com\"><text>Garmin International</text></link><time>2017-07-13T16:53:53Z</time></metadata>
                        # ", all_gpx, collapse="", sep="")

                        write(x = all_gpx, file = gpx_file_name)
                      },
                      write_kml_file=function(kml_file_name="test.kml", kml_track_name="test", color="red", lwd=3, kml_name="", kml_description="")
                      {
                        "Writes a KML file representation of the xtrack. KML files are used in Google Earth and elsewhere."
                        # kml_file_name="test.kml"
                        # kml_track_name="test"
                        # color="red"
                        # lwd=3
                        # kml_name=""
                        # kml_description=""
                        # lon <- xt_1$trackpoints$lon
                        # lat <- xt_1$trackpoints$lat
                        the_spatial_lines <- SpatialLines(list(Lines(Line(cbind(trackpoints$lon,trackpoints$lat)), ID="a")))
                        emptyData <- data.frame(matrix(0, ncol = 2, nrow = length(the_spatial_lines)))
                        the_spatialLinesDataFrame <- SpatialLinesDataFrame(sl=the_spatial_lines, data=emptyData, match.ID=FALSE)
                        kmlLines(obj=the_spatialLinesDataFrame, kmlfile=kml_file_name, name=kml_name, col=color,lwd=lwd,kmlname=kml_name, kmldescription=kml_description)

                      }




))



#' Analysis of rates of habitat exploration.
#'
#' Following Wood et al. 2021, this analysis computes the habitat visited / explored
#' each day, the marginal 'new' habitat visited for each day, and the cumulative habitat explored across all days.
#'
#' @param list_of_xtracks This analysis should be done with a temporally-sorted list of xtrack objects, with each xtrack representing
#' one day of travel of the same person. In the list, the first day of data should be in position 1.
#' @param cell_size_m The x and y size of each raster cell in the raster representation of the landscape, in meters.
#' @return A dataframe (format described below) in which rows represent days, and columns provide analysis outcomes
#' \describe{
#'   \item{day}{day number}
#'   \item{sum_cells_visited_this_day}{the number of cells intersected on that day}
#'   \item{cum_sum_cells_visited_across_days}{cumulative number of cells intersected from day 1 to current day}
#'   \item{n_new_cells_visited_this_day}{number of cells visited that day and not on prior days}
#'   \item{square_meters_per_cell}{the area of each cell in square meters}
#'   \item{sum_square_meters_visited_this_day}{area in square meters of all cells visited that day}
#'   \item{cum_sum_square_meters_visited_across_days}{cummulative area of cells visited from day 1 to the current day}
#'   \item{square_meters_new_habitat_visited_this_day}{the square meters of cells visited on that day and not on prior days.}
#' }
#' @examples
#' xt_1 <- xtrack(lat=d1$lat, lon=d1$lon, elevation_m=d1$elevation_m, in_camp=d1$in_camp, unix_time=d1$unix_time, distance_from_camp_m=d1$distance_from_camp_m, utm_epsg=32736)
#' xt_2 <- xtrack(lat=d2$lat, lon=d2$lon, elevation_m=d2$elevation_m, in_camp=d2$in_camp, unix_time=d2$unix_time, distance_from_camp_m=d2$distance_from_camp_m, utm_epsg=32736)
#' xt_3 <- xtrack(lat=d3$lat, lon=d3$lon, elevation_m=d3$elevation_m, in_camp=d3$in_camp, unix_time=d3$unix_time, distance_from_camp_m=d3$distance_from_camp_m, utm_epsg=32736)
#' xt_4 <- xtrack(lat=d4$lat, lon=d4$lon, elevation_m=d4$elevation_m, in_camp=d4$in_camp, unix_time=d4$unix_time, distance_from_camp_m=d4$distance_from_camp_m, utm_epsg=32736)
#' xt_5 <- xtrack(lat=d5$lat, lon=d5$lon, elevation_m=d5$elevation_m, in_camp=d5$in_camp, unix_time=d5$unix_time, distance_from_camp_m=d5$distance_from_camp_m, utm_epsg=32736)
#' xt_6 <- xtrack(lat=d6$lat, lon=d6$lon, elevation_m=d6$elevation_m, in_camp=d6$in_camp, unix_time=d6$unix_time, distance_from_camp_m=d6$distance_from_camp_m, utm_epsg=32736)
#' xt_7 <- xtrack(lat=d7$lat, lon=d7$lon, elevation_m=d7$elevation_m, in_camp=d7$in_camp, unix_time=d7$unix_time, distance_from_camp_m=d7$distance_from_camp_m, utm_epsg=32736)
#' xt_8 <- xtrack(lat=d8$lat, lon=d8$lon, elevation_m=d8$elevation_m, in_camp=d8$in_camp, unix_time=d8$unix_time, distance_from_camp_m=d8$distance_from_camp_m, utm_epsg=32736)
#' list_of_xtracks <- list(xt_1, xt_2, xt_3, xt_4, xt_5, xt_6, xt_7, xt_8)
#' hab_exp_results <- get_hab_exp_across_days(list_of_xtracks, cell_size_m = 10)
#' @seealso For Land exploration analysis described in Wood et al. 2021 publication see \url{https://www.nature.com/articles/s41562-020-01002-7#Sec15}
#' @export

get_hab_exp_across_days<-function(list_of_xtracks, cell_size_m=10)
{

  n_xtracks <- length(list_of_xtracks)
  # calculate the extent needed to contain all the supplied tracks

  min_y_all_xtracks <- NA
  max_y_all_xtracks <- NA
  min_x_all_xtracks <- NA
  max_x_all_xtracks <- NA

  for(i in 1:n_xtracks)
  {
      min_y_all_xtracks <- min(min_y_all_xtracks, list_of_xtracks[[i]]$min_y_utm, na.rm=TRUE)
      max_y_all_xtracks <- max(max_y_all_xtracks, list_of_xtracks[[i]]$max_y_utm, na.rm=TRUE)
      min_x_all_xtracks <- min(min_x_all_xtracks, list_of_xtracks[[i]]$min_x_utm, na.rm=TRUE)
      max_x_all_xtracks <- max(max_x_all_xtracks, list_of_xtracks[[i]]$max_x_utm, na.rm=TRUE)
  }

  #this is just for the aesthetics when the raster's are plotted it is
  #good to have a little buffer / margin for sanity's sake
  min_y_all_xtracks <-min_y_all_xtracks-20
  max_y_all_xtracks <-max_y_all_xtracks+20
  min_x_all_xtracks <-min_x_all_xtracks-20
  max_x_all_xtracks <-max_x_all_xtracks+20

  #create list to hold raster representations of xtracks
  list_of_rasters <- vector(length=n_xtracks, mode = "list")

  for(i in 1:n_xtracks)
  {
    list_of_rasters[[i]] <- list_of_xtracks[[i]]$as_raster_of_habitat_visited_binary(cell_size_m = cell_size_m, xmin=min_x_all_xtracks, xmax=max_x_all_xtracks, ymin=min_y_all_xtracks, ymax=max_y_all_xtracks)
  }

  #create data frame to hold results of habitat exploration analysis

  res <- data.frame(day=1:n_xtracks)
  res$sum_cells_visited_this_day <- NA
  res$cum_sum_cells_visited_across_days <- NA
  res$n_new_cells_visited_this_day <- NA

  cum_sum_raster <- list_of_rasters[[1]]

  res$sum_cells_visited_this_day[1] <- sum(values(list_of_rasters[[1]]))
  res$cum_sum_cells_visited_across_days[1] <- sum(values(cum_sum_raster))
  res$n_new_cells_visited_this_day[1] <- sum(values(list_of_rasters[[1]]))

  #the heart of the matter is here
  for(i in 2:n_xtracks)
  {
    cum_sum_raster <- cum_sum_raster + list_of_rasters[[i]]
    cum_sum_raster[cum_sum_raster[]>1] <- 1
    res$sum_cells_visited_this_day[i] <- sum(values(list_of_rasters[[i]]))
    res$cum_sum_cells_visited_across_days[i] <- sum(values(cum_sum_raster))
    res$n_new_cells_visited_this_day[i] <- res$cum_sum_cells_visited_across_days[i] - res$cum_sum_cells_visited_across_days[i-1]
  }

  square_meters_per_cell <- cell_size_m * cell_size_m

  res$square_meters_per_cell <- square_meters_per_cell

  res$sum_square_meters_visited_this_day <- res$sum_cells_visited_this_day * square_meters_per_cell

  res$cum_sum_square_meters_visited_across_days <- res$cum_sum_cells_visited_across_days * square_meters_per_cell

  res$square_meters_new_habitat_visited_this_day <- res$n_new_cells_visited_this_day * square_meters_per_cell

  #the concept of 'new habitat visited today' doesn't make sense on the first day, hence NA
  res$n_new_cells_visited_this_day[1] <- NA
  res$square_meters_new_habitat_visited_this_day[1] <- NA

  return(res)
}


equalizeAxesLimits<- function(x_lim, y_lim)
{

  x_axis_meters <- distm(x = c(x_lim[1], y_lim[1]), y=c(x_lim[2], y_lim[1]))[1,1]
  y_axis_meters <- distm(x= c(x_lim[1], y_lim[1]), y=c(x_lim[1], y_lim[2]))[1,1]

  new_min_x <- 0
  new_max_x <- 0
  new_min_y <- 0
  new_max_y <- 0
  new_x_axis_distance_m <- 0
  new_y_axis_distance_m <- 0
  new_x_lim <- c(0,0)
  new_y_lim <- c(0,0)

  if(x_axis_meters<=y_axis_meters)
  {
    mid_x <- (x_lim[1]+x_lim[2])/2
    new_min_x <- destPoint(p = c(mid_x, y_lim[1]), b = 270, d = y_axis_meters/2)[1]
    new_max_x <- destPoint(p = c(mid_x, y_lim[1]), b = 90, d = y_axis_meters/2)[1]
    new_x_axis_distance_m <- distm(c(new_min_x, y_lim[1]), c(new_max_x, y_lim[1]))[1,1]
    new_y_axis_distance_m <- y_axis_meters
    new_x_lim <- c(new_min_x, new_max_x)
    new_y_lim <- y_lim
  } else if(x_axis_meters>y_axis_meters) {

    mid_y <- (y_lim[1]+y_lim[2])/2
    new_min_y <- destPoint(p=c(x_lim[1], mid_y), b=180, d=x_axis_meters/2)[2]
    new_max_y <- destPoint(p=c(x_lim[1], mid_y), b=0, d=x_axis_meters/2)[2]
    new_x_axis_distance_m <- x_axis_meters
    new_y_axis_distance_m <- distm(c(x_lim[1], new_min_y), c(x_lim[1], new_max_y))[1,1]
    new_x_lim <- x_lim
    new_y_lim <- c(new_min_y, new_max_y)
  }
  return(list(x_lim=new_x_lim, y_lim=new_y_lim, new_x_axis_distance_m=new_x_axis_distance_m, new_y_axis_distance_m=new_y_axis_distance_m))
}

