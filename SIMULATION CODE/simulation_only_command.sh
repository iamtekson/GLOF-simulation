#!/bin/sh

#########################################################
############# Hydrocoordinates ##########################
# 0m   # 500556,5581835,1000,180 # updated coordinate # 500559,5582923,700,180
# 500m # 500977,5581969,1000,78.69
# 1000m # 501242,5581614,1000,315
# 2000m # 502128,5581708,2000,180
# 4000m # 503788,5581056,2000,315
#########################################################
hydrocoords=500559,5581930,700,180,501242,5581614,1000,135,502128,5581708,2000,180


# Set of simulation to run the model
r.avaflow prefix=11 phases=m elevation=elevation11 hrelease1=ls11 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
r.avaflow prefix=12 phases=m elevation=elevation12 hrelease1=ls12 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
r.avaflow prefix=13 phases=m elevation=elevation13 hrelease1=ls13 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
r.avaflow prefix=14 phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
#r.avaflow prefix=14b phases=m elevation=elevation14 hrelease1=ls14 hrelease3=lake time=50,3600 hydrocoords=$hydrocoords

r.avaflow prefix=21 phases=m elevation=elevation21 hrelease1=ls21 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
r.avaflow prefix=22 phases=m elevation=elevation22 hrelease1=ls22 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
r.avaflow prefix=23 phases=m elevation=elevation23 hrelease1=ls23 hrelease3=lake time=5,300 hydrocoords=$hydrocoords

r.avaflow prefix=31 phases=m elevation=elevation31 hrelease1=ls31 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
r.avaflow prefix=32 phases=m elevation=elevation32 hrelease1=ls32 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
r.avaflow prefix=33 phases=m elevation=elevation33 hrelease1=ls33 hrelease3=lake time=5,300 hydrocoords=$hydrocoords

r.avaflow prefix=41 phases=m elevation=elevation41 hrelease1=ls41 hrelease3=lake time=5,300 hydrocoords=$hydrocoords
