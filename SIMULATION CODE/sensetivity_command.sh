#!/bin/sh

# Lake map (m)
r.in.gdal -o --overwrite input=lake/lake.tif output=lake

# Elevation map (m)
r.in.gdal -o --overwrite input=elevation/elevation11.tif output=elevation11
r.in.gdal -o --overwrite input=elevation/elevation12.tif output=elevation12
r.in.gdal -o --overwrite input=elevation/elevation13.tif output=elevation13
r.in.gdal -o --overwrite input=elevation/elevation14.tif output=elevation14

r.in.gdal -o --overwrite input=elevation/elevation21.tif output=elevation21
r.in.gdal -o --overwrite input=elevation/elevation22.tif output=elevation22
r.in.gdal -o --overwrite input=elevation/elevation23.tif output=elevation23

r.in.gdal -o --overwrite input=elevation/elevation31.tif output=elevation31
r.in.gdal -o --overwrite input=elevation/elevation32.tif output=elevation32
r.in.gdal -o --overwrite input=elevation/elevation33.tif output=elevation33

r.in.gdal -o --overwrite input=elevation/elevation41.tif output=elevation41



#Landslide maps (m)
r.in.gdal -o --overwrite input=ls/ls11.tif output=ls11
r.in.gdal -o --overwrite input=ls/ls12.tif output=ls12
r.in.gdal -o --overwrite input=ls/ls13.tif output=ls13
r.in.gdal -o --overwrite input=ls/ls14.tif output=ls14

r.in.gdal -o --overwrite input=ls/ls21.tif output=ls21
r.in.gdal -o --overwrite input=ls/ls22.tif output=ls22
r.in.gdal -o --overwrite input=ls/ls23.tif output=ls23

r.in.gdal -o --overwrite input=ls/ls31.tif output=ls31
r.in.gdal -o --overwrite input=ls/ls32.tif output=ls32
r.in.gdal -o --overwrite input=ls/ls33.tif output=ls33

r.in.gdal -o --overwrite input=ls/ls41.tif output=ls41


# Set the region of the grass database to the elevation map
g.region -s rast=elevation11


# Printout the region
g.region -p

#########################################################
############# Hydrocoordinates ##########################
# 0m   # 500556,5581835,1000,180 # updated coordinate # 500559,5582923,700,180
# 500m # 500977,5581969,1000,78.69
# 1000m # 501242,5581614,1000,315
# 2000m # 502128,5581708,2000,180
# 4000m # 503788,5581056,2000,315
#########################################################
hydrocoords=500559,5581930,700,180,501242,5581614,1000,135,502128,5581708,2000,180

friction=35,20,3,0,0,0,0.05


# Set of simulation to run the model
r.avaflow prefix=35_20 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=40,20,3,0,0,0,0.05

r.avaflow prefix=40_20 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=40_20 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=40_20 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=40_20 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=30,20,3,0,0,0,0.05

r.avaflow prefix=30_20 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=30_20 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=30_20 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=30_20 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=35,25,3,0,0,0,0.05

r.avaflow prefix=35_25 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_25 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_25 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_25 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=35,15,3,0,0,0,0.05

r.avaflow prefix=35_15 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_15 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_15 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_15 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=35,20,3,0,0,0,0.03

r.avaflow prefix=35_20_03 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_03 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_03 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_03 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=35,20,3,0,0,0,0.07

r.avaflow prefix=35_20_07 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_07 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_07 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_07 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=35,20,3,0,0,0,0.04

r.avaflow prefix=35_20_04 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_04 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_04 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_04 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction

friction=35,20,3,0,0,0,0.06

r.avaflow prefix=35_20_06 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_06 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_06 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction
r.avaflow prefix=35_20_06 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords friction=$friction