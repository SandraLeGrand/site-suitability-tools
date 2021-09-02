#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Geomorphic Oscillation Assessment Tool
# Author: Mike Ekegren ERDC 04 FEB 2020
# Modified by S. LeGrand ERDC 18 MAR 2020;
#        11 JAN 2021 - added code to print plots and 
#                      set directory paths
#      
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# comment this out if you do not want the script to
# clear your R global environment
rm(list=ls()) # remove old variables from memory

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# set the input data directory path (required)
iDir <- "."

# set the output data directory path (required)
oDir <- "."

# set file names (required)
dsm_file    <- "dsm.tif"            # DSM GeoTiff file
dtm_file    <- "dtm.tif"            # DSM GeoTiff file
nlcd_file   <- "NLCD.tif"           # Landcover GeoTiff file
output_file <- "GOATout.csv"        # output csv file

# set domain attributes (use decimal degrees; required)
BASE_LON <- -116.7928  # domain origin-point longitude (deg) 
BASE_LAT <- 35.3870    # domain origin-point latitude (deg)
subDist <- 400         # cardinal direction distance (m)

# tunable parameters
sampRes <- 25  # domain Matrix sampling resolution (m)
sampRad <- 35  # circular sampling radius (m)
slopeAcc <- 7  # maximum macroslope threshold (m)

####### end of user input  

# additional prescribed values
sampInt <- 2.5             # circular sampling interval (deg)
OOC <- 2000                # nonsuitable site indicator (deg)
vsn <- 0.001               # very small number
nLCF <- 95                 # total number of land cover classes
LCF_hazard <- c(11,90,95)  # water-related NLCD land cover classes

# set binary land cover filter vector
LCF <- c(rep(1,nLCF)) 
LCF[LCF_hazard] <- 0 

# set vector of circle sampling direction angles
sampDir <- c(seq(0,359, sampInt))

# load required libraries
library(raster)       # performs distance/direction for sampling
library(rgdal)        # geospatial library
library(geosphere)    # performs distance/direction for sampling
library(dplyr)        # used for lag function 
library(ggplot2)      # only needed for R plotting
library(RColorBrewer) # only needed for R plotting

# set file paths 
# may need to reverse slash depending on machine OS
dtm_path <- paste0(iDir,"/",dtm_file)        # DTM GeoTiff path 
dsm_path <- paste0(iDir,"/",dsm_file)        # DSM GeoTiff path
landcover_path <- paste0(iDir,"/",nlcd_file) # NLCD GeoTiff path
output_path <- paste0(oDir,"/",output_file)  # output path

# read in raster files
dtm <- raster(dtm_path)
dsm <- raster(dsm_path)
landcover <- raster(landcover_path)

# establish key domain placement attributes
# includes central point, NE corner, and NW corner of the Matrix
basePoint <- matrix(c( BASE_LON, BASE_LAT ), nrow=1)
northeast <- destPoint(basePoint, 45, sqrt((subDist^2)+(subDist^2)))
northwest <- destPoint(basePoint, 315, sqrt((subDist^2)+(subDist^2)))

# set frontage lon/lat ends and depth of zone
fZf <- matrix(c(northwest, northeast), nrow=2)
fZd <- subDist*2

# determine lon/lat pairs along frontage line
fZfBearing <- bearing(fZf[,1], fZf[,2])
fZfDist    <- pointDistance(fZf[,1], fZf[,2], lonlat=TRUE) 
fZSamp     <- floor(fZfDist/sampRes)
fZfPoint   <- gcIntermediate(fZf[,1], fZf[,2], fZSamp,
                             addStartEnd=TRUE, sp=FALSE, sepNA=FALSE)

# determine zone depth bearing
fZdBearing <- fZfBearing + 90 
fZdBearing <- ifelse(fZdBearing > 359, fZdBearing - 360, fZdBearing)

# build points along rear (south) boundary, same number as fZfPoint
fZdRear <- destPoint(fZfPoint, fZdBearing, fZd)

# fill in the rest of the domain matrix lon/lat values
sampPoints <- matrix(ncol=2) 
for(i in 1:nrow(fZfPoint)){ 
  newPoints <- gcIntermediate(fZfPoint[i,], fZdRear[i,], fZSamp, 
                              addStartEnd=TRUE, sp=FALSE, sepNA=FALSE)
  sampPoints <- rbind(sampPoints, newPoints) 
}

#remove the first row which only has NA values
sampPoints <- sampPoints[-1,]

# crop the input datasets
calcExt <- c(xmin=min(sampPoints[,1])-vsn, xmax=max(sampPoints[,1])+vsn, 
             ymin=min(sampPoints[,2])-vsn, ymax=max(sampPoints[,2])+vsn) 
dtm <- crop(dtm, calcExt) 
dsm <- crop(dsm, calcExt)
landcover <- crop(landcover, calcExt, snap='out')

# create master data frame
spdf <- as.data.frame(sampPoints) 
colnames(spdf) <- c("Lon", "Lat") 
spdf$DTM <- dtm[cellFromXY(dtm,spdf)] 
spdf$DSM <- dsm[cellFromXY(dsm,spdf)] 
spdf$LC <- LCF[landcover[cellFromXY(landcover,spdf)]] 

# establish analytical Point Function
Topo <- function(evalPoint){ 
  
  # sample points around circle
  localPts <- destPoint(evalPoint[1,1:2], sampDir, sampRad)
  pointDis <- pointDistance(localPts[1,], localPts[2,], lonlat=TRUE) 
  
  # establish temporary data frame
  local_df <- as.data.frame(localPts)
  colnames(local_df) <- c("Lon", "Lat")
  
  # set up dtm values
  local_df$DTM <- dtm[cellFromXY(dtm,local_df)]
  local_df$DTML <- lag(local_df$DTM) # lag by 1 row
  local_df$DTML[1] <- local_df$DTM[nrow(local_df)] # apply bottom # row to top
  
  # calculate adjacent circle point elevation difference
  local_df$ZDelta <- local_df$DTM - local_df$DTML
  
  # estimate initial roughness
  vertChange <- sum(abs(local_df$ZDelta))
  
  # calculate the slope in degrees to each circle sample point
  local_df$DTM_SLP <- atan((local_df$DTM - evalPoint[1,3]) / sampRad) * (180/pi)
  
  # determine steepest slope and location
  steepestSlope <- max(abs(local_df$DTM_SLP))
  steepestSlopeNum <- which.max(abs(local_df$DTM_SLP))
  
  # determine influential slope
  influenceSlope <- c(steepestSlopeNum, 
                      steepestSlopeNum-(nrow(local_df)/2))
  influenceSlope <- ifelse(influenceSlope < 1,
                           influenceSlope + nrow(local_df), 
                           influenceSlope)
  slopeInf <- atan((abs(local_df$DTM[influenceSlope[1]] - 
                        local_df$DTM[influenceSlope[2]])) / 
                   (sampRad*2)) * (180/pi)
  
  # estimate deviation from planar
  planeDev <- mean(local_df$DTM_SLP)
  
  # sample dsm at each radial point around the circle
  local_df$DSM <- dsm[cellFromXY(dsm, local_df)] 
  
  # determine degree of climb for each radial point
  local_df$CLIMB <- atan((local_df$DSM - evalPoint[1,3]) / sampRad) * (180/pi)
  maxClimb <- max(local_df$CLIMB) # steepest climb
  
  # identify any tree line hazards
  treeLine <- ifelse(any(local_df$DSM - local_df$DTM > 3), 0, 1)
  
  # establish output matrix
  outPut <- matrix(c(vertChange, steepestSlope, slopeInf, 
                     planeDev, treeLine, maxClimb), ncol=6)
  return(outPut)
}

# establish data passing function
pointValue <- function(spdfRows){ 
  obsMat <- matrix(nrow=nrow(spdfRows), ncol=6)       
  for(i in 1:nrow(spdfRows)){
    obsMat[i,1:6] <- Topo(spdfRows[i,])}
  return(obsMat)
}

# run calculation
localObs <- pointValue(spdf)

# add variables from return matrix to master data frame
spdf$VertDelta <- localObs[,1] # initial roughness
spdf$SteepSlp  <- localObs[,2] # steepest slope
spdf$InflSlp   <- localObs[,3] # influential slope
spdf$Plane     <- localObs[,4] # planer deviation
spdf$Obstacle  <- localObs[,5] # tree line obstacle indicator
spdf$MaxClimb  <- localObs[,6] # maximum climb

# calculate roughness diagnostic
spdf$Rough <- ifelse((spdf$SteepSlp <= slopeAcc), 
                     spdf$VertDelta - 
                     ((spdf$InflSlp -abs(spdf$Plane))* 2.447438), OOC)

# incorporate land cover and vertical obstructions 
spdf$Suitability <- ifelse(((spdf$LC==1) & 
                            abs(spdf$DSM-spdf$DTM < 1) & 
                            (spdf$Obstacle==1)), 
                           spdf$Rough, OOC)

####### plot results
# show macro-slope
plot_SteepSlp <- ggplot(spdf, aes(x=Lon, y=Lat, colour=SteepSlp))+ 
  geom_point(shape=16, size=1.5)+ 
  theme_bw()+
  theme(legend.position="bottom", 
        legend.direction = "horizontal")+
  guides(colour = guide_colourbar(title.position = "top",
                                barwidth = 15, 
                                barheight = 0.75))+
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_colour_gradient2(name = "Steepest Slope (decimal degress)",
    low="forestgreen", mid="yellow", high="red2", 
    midpoint=(slopeAcc/2), limits=c(0,slopeAcc))

# show initial terrain roughness estimate
plot_VertDelta <- ggplot(spdf, aes(x=Lon, y=Lat, colour=VertDelta))+ 
  geom_point(shape=16, size=1.5)+ 
  theme_bw()+
  theme(legend.position="bottom", 
        legend.direction = "horizontal")+
  guides(colour = guide_colourbar(title.position = "top",
                                  barwidth = 15, 
                                  barheight = 0.75))+
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_colour_gradient2(name = "Initial Roughness (m)",
    low="forestgreen", mid="yellow", high="red2", 
    midpoint=10, limits=c(0,20)) 

# show adjusted terrain roughness
plot_Rough <- ggplot(spdf, aes(x=Lon, y=Lat, colour=Rough))+ 
  geom_point(shape=16, size=1.5)+ 
  theme_bw()+
  theme(legend.position="bottom", 
        legend.direction = "horizontal")+
  guides(colour = guide_colourbar(title.position = "top",
                                  barwidth = 15, 
                                  barheight = 0.75))+
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_colour_gradient2(name = "Terrain Roughness (m)",
    low="forestgreen", mid="yellow", high="red2", 
    midpoint=2, limits=c(0,4)) 

# show integrated terrain suitability
plot_Suitability <- ggplot(spdf, aes(x=Lon, y=Lat, colour=Suitability))+ 
  geom_point(shape=16, size=1.5)+ 
  theme_bw()+
  theme(legend.position="bottom", 
        legend.direction = "horizontal")+
  guides(colour = guide_colourbar(title.position = "top",
                                  barwidth = 15, 
                                  barheight = 0.75))+
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_colour_gradient2(name = "HLZ Suitability",
    low="forestgreen", mid="yellow", high="red2", 
    midpoint=2, limits=c(0,4))

# show maximum climb
plot_MaxClimb <- ggplot(spdf, aes(x=Lon, y=Lat, colour=MaxClimb))+ 
  geom_point(shape=16, size=1.5)+ 
  theme_bw()+
  theme(legend.position="bottom", 
        legend.direction = "horizontal")+
  guides(colour = guide_colourbar(title.position = "top",
                                  barwidth = 15, 
                                  barheight = 0.75))+
  labs(x = "Longitude", y = "Latitude")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_colour_gradient2(name = "Maximum climb (decimal degrees)",
    low="forestgreen", mid="yellow", high="red2", 
    midpoint=7.5, limits=c(0,15))

####### write results

# announce progress
print("Writing output csv file.")
write.csv(spdf, output_path)

# announce progress
print("Saving plots.")

# create output file paths for plots
out_SteepSlp <- paste0(oDir,"/plot_SteepSlp.jpg")
out_VertDelta <- paste0(oDir,"/plot_VertDelta.jpg")
out_Rough <- paste0(oDir,"/plot_Rough.jpg")
out_Suitability <- paste0(oDir,"/plot_Suitability.jpg")
out_MaxClimb <- paste0(oDir,"/plot_MaxClimb.jpg")
  
# save plots
ggsave(out_SteepSlp, plot = plot_SteepSlp, width = 6.5, height = 7)
ggsave(out_VertDelta, plot = plot_VertDelta, width = 6.5, height = 7)
ggsave(out_Rough, plot = plot_Rough, width = 6.5, height = 7)
ggsave(out_Suitability, plot = plot_Suitability, width = 6.5, height = 7)
ggsave(out_MaxClimb, plot = plot_MaxClimb, width = 6.5, height = 7)

# announce progress
print("Processing complete.")
