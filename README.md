# site-suitability-tools
Terrain characterization tools for site suitability assessment.

## Authors
* Michael T. Ekegren (ERDC)
* Sandra L. LeGrand (ERDC)

## Created 
Summer 2021

## Purpose
This repository includes functional area terrain evalation tools to support site suitability assessment applications.

## Code
**Geomorphic Oscillation Assessment Tool (GOAT.R)**<br/>
The Geomorphic Oscillation Assessment Tool (GOAT) quantifies terrain roughness as a mechanism to better explain helicopter landing zone (HLZ) suitability for aviation applications. Surface roughness is a critical discriminator for site utility in complex terrain. GOAT uses a spatial sampling of high-resolution elevation and land cover data to construct data frames, which enable a relational analysis of component and aggregate site suitability. By incorporating multiple criteria from various doctrinal sources, GOAT produces a composite, first-order quality assessment of the areal options available for use. 

This script uses a circular sampling technique to assess the general roughness of a location's surroundings. A user can modify the circle’s radius to suit their particular application needs, as long as the radius is greater than or equal to the required physical footprint of the helicopter (i.e., the area encompassed by the rotor disk and the tail). By assessing the variability in elevation from adjacent points that make up the circle, the tool can draw inferences about the general roughness of the landing zone. This circular sampling approach increases the likelihood that linear features, such as drainages, in the proximity of a potential area of interest are considered.

Please see the report by Ekegren and LeGrand (2021) for a detailed overview of the script equations and code. 

Ekegren, M. T. and LeGrand S. L.: Incorporating terrain roughness into helicopter landing zone site selection by using the Geomorphic Oscillation Assessment Tool (GOAT) v1.0. ERDC/CRREL SR-20-DRAFT, U.S. Army Engineer Research and Development Center, Hanover, New Hampshire, USA, (in press).

