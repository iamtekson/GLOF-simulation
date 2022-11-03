#!/bin/sh

# Lake map (m)
r.in.gdal -o --overwrite input=hrelease3.tif output=hrelease3


# Elevation map (m)
r.in.gdal -o --overwrite input=elevation11.tif output=elevation11
r.in.gdal -o --overwrite input=elevation12.tif output=elevation12
r.in.gdal -o --overwrite input=elevation13.tif output=elevation13
r.in.gdal -o --overwrite input=elevation14.tif output=elevation14

r.in.gdal -o --overwrite input=elevation21.tif output=elevation21
r.in.gdal -o --overwrite input=elevation22.tif output=elevation22
r.in.gdal -o --overwrite input=elevation23.tif output=elevation23
r.in.gdal -o --overwrite input=elevation24.tif output=elevation24

r.in.gdal -o --overwrite input=elevation31.tif output=elevation31
r.in.gdal -o --overwrite input=elevation32.tif output=elevation32
r.in.gdal -o --overwrite input=elevation33.tif output=elevation33
r.in.gdal -o --overwrite input=elevation34.tif output=elevation34


#Landslide maps (m)
r.in.gdal -o --overwrite input=ls11.tif output=ls11
r.in.gdal -o --overwrite input=ls12.tif output=ls12
r.in.gdal -o --overwrite input=ls13.tif output=ls13
r.in.gdal -o --overwrite input=ls14.tif output=ls14

r.in.gdal -o --overwrite input=ls21.tif output=ls21
r.in.gdal -o --overwrite input=ls22.tif output=ls22
r.in.gdal -o --overwrite input=ls23.tif output=ls23
r.in.gdal -o --overwrite input=ls24.tif output=ls24

r.in.gdal -o --overwrite input=ls31.tif output=ls31
r.in.gdal -o --overwrite input=ls32.tif output=ls32
r.in.gdal -o --overwrite input=ls33.tif output=ls33
r.in.gdal -o --overwrite input=ls34.tif output=ls34


# Set the region of the grass database to the elevation map
g.region -s rast=elevation11


# Printout the region
g.region -p


#########################################################
############# Hydrocoordinates ##########################
# 500m # 500977,5581969,1000,78.69
# 1000m # 501242,5581614,1000,315
# 2000m # 502128,5581708,2000,180

# the 250, 500, 750, 1000, 2000 hydrocoordinates are used for the hydrological model
#########################################################
hydrocoords=500977,5581969,1000,78.69,501242,5581614,1000,315,501141,5581815,100,180,502128,5581708,2000,180


# Set of simulation to run the model
r.avaflow prefix=31 phases=m elevation=elevation31 hrelease1=ls31 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=32 phases=m elevation=elevation32 hrelease1=ls32 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=33 phases=m elevation=elevation33 hrelease1=ls33 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=34 phases=m elevation=elevation34 hrelease1=ls34 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords

r.avaflow prefix=21 phases=m elevation=elevation21 hrelease1=ls21 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=22 phases=m elevation=elevation22 hrelease1=ls22 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=23 phases=m elevation=elevation23 hrelease1=ls23 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=24 phases=m elevation=elevation24 hrelease1=ls24 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords

r.avaflow prefix=11 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=12 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords 
r.avaflow prefix=13 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords
r.avaflow prefix=14 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=hrelease3 controls=0,0,0,1,0,0 hydrocoords=$hydrocoords





