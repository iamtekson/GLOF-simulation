#!/usr/bin/R

##############################################################################
#
# MODULE:       r.avaflow.map.R
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The simulation model for avalanche and debris flows
#               Script for the creation of map plots of simulation results
#
# COPYRIGHT:    (c) 2013 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 1993 - 2021 by the R Development Core Team
#
# VERSION:      20210525 (25 May 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2).
#
##############################################################################

# Loading libraries
library(maptools)
library(sp)
library(rgdal)
library(raster)

# Importing arguments defined in r.avaflow.py
args <- commandArgs(trailingOnly = TRUE)
temppath <- args[1]  #temporary directory
prefix <- args[2]  #prefix for file names
model <- as.integer(args[3])  #model
ntimesteps <- as.numeric(args[4])  #number of time step
fill <- args[5]  #complete string for time step
contstep <- args[6]  #interval of flow height contours
thrs5 <- as.numeric(args[7])  #legend thresholds
thrs4 <- as.numeric(args[8])
thrs3 <- as.numeric(args[9])
thrs2 <- as.numeric(args[10])
thrs1 <- as.numeric(args[11])
thrsm <- as.numeric(args[12])
thrs5n <- as.numeric(args[13])
thrs4n <- as.numeric(args[14])
thrs3n <- as.numeric(args[15])
thrs2n <- as.numeric(args[16])
thrs1n <- as.numeric(args[17])
thrsmn <- as.numeric(args[18])
ymax <- as.numeric(args[19])  #coordinates defining map boundaries
ymin <- as.numeric(args[20])
xmin <- as.numeric(args[21])
xmax <- as.numeric(args[22])
cellsize <- as.numeric(args[23])  #pixel size
tint <- as.numeric(args[24])  #length of time step
tstop <- as.numeric(args[25])  #length of total simulation
impdef <- as.integer(args[26])  #control for observed impact area map
depdef <- as.integer(args[27])  #control for observed deposit map
maxvflow <- as.numeric(args[28])  #maximum flow velocity
mapprof <- as.integer(args[29])  #availability of profile data
hydrograph <- as.integer(args[30])  #hydrograph
releasemass <- as.integer(args[31])  #presence of release mass
ninhyd <- as.integer(args[32])  #number of input hydrographs
nouthyd <- as.integer(args[33])  #number of output hydrographs
mstring <- args[34]  #parameter to display
mconv <- args[35]  #unit conversion factor for parameter
mdig <- as.integer(args[36])  #number of digits after comma for colour bar labels
ntimemax <- as.integer(args[37])  #total number of time steps
ctrlpts <- as.integer(args[38])  #control for control points
ortho <- as.character(args[39])  #path to orthophoto
phase1 <- as.character(args[40])  #first phase
phase2 <- as.character(args[41])  #second phase
phase3 <- as.character(args[42])  #third phase

# Reading flow velocity and direction data from file defining control variable:
if (ntimesteps > 0 && ntimesteps <= ntimemax) {
    velocity <- 1
} else {
    velocity <- 0
}
if (velocity == 1) {

    dirline <- ntimesteps - 1

    intable = paste(prefix, "_results/", prefix, "_files/", prefix, "_directions1.txt", 
        sep = "")
    dirx <- scan(intable, nlines = 1)  #(PHASE 1) x coordinate
    diry <- scan(intable, skip = 1, nlines = 1)  #(PHASE 1) y coordinate

    if (model <= 3) {

        dirh <- scan(intable, skip = 2 + 3 * dirline, nlines = 1)  #flow height
        dirvx <- scan(intable, skip = 3 + 3 * dirline, nlines = 1)  #flow velocity in x direction
        dirvy <- scan(intable, skip = 4 + 3 * dirline, nlines = 1)  #flow velocity in y direction

    } else if (model > 3 && model < 7) {

        intable2 = paste(prefix, "_results/", prefix, "_files/", prefix, "_directions2.txt", 
            sep = "")
        dirxf <- scan(intable2, nlines = 1)  #PHASE 2 x coordinate
        diryf <- scan(intable2, skip = 1, nlines = 1)  #PHASE 2 y coordinate

        dirh <- scan(intable, skip = 2 + 3 * dirline, nlines = 1)  #PHASE 1 flow height
        dirvx <- scan(intable, skip = 3 + 3 * dirline, nlines = 1)  #PHASE 1 flow velocity in x direction
        dirvy <- scan(intable, skip = 4 + 3 * dirline, nlines = 1)  #PHASE 1 flow velocity in y direction
        dirhf <- scan(intable2, skip = 2 + 3 * dirline, nlines = 1)  #PHASE 2 flow height
        dirvxf <- scan(intable2, skip = 3 + 3 * dirline, nlines = 1)  #PHASE 2 flow velocity in x direction
        dirvyf <- scan(intable2, skip = 4 + 3 * dirline, nlines = 1)  #PHASE 2 flow velocity in y direction

        dirhf[dirhf > 0] <- 1

    } else if (model == 7) {

        intable2 = paste(prefix, "_results/", prefix, "_files/", prefix, "_directions2.txt", 
            sep = "")
        dirxf <- scan(intable2, nlines = 1)  #PHASE 2 x coordinate
        diryf <- scan(intable2, skip = 1, nlines = 1)  #PHASE 2 y coordinate

        dirh <- scan(intable, skip = 2 + 3 * dirline, nlines = 1)  #PHASE 1 flow height
        dirvx <- scan(intable, skip = 3 + 3 * dirline, nlines = 1)  #PHASE 1 flow velocity in x direction
        dirvy <- scan(intable, skip = 4 + 3 * dirline, nlines = 1)  #PHASE 1 flow velocity in y direction
        dirhf <- scan(intable2, skip = 2 + 3 * dirline, nlines = 1)  #PHASE 2 flow height
        dirvxf <- scan(intable2, skip = 3 + 3 * dirline, nlines = 1)  #PHASE 2 flow velocity in x direction
        dirvyf <- scan(intable2, skip = 4 + 3 * dirline, nlines = 1)  #PHASE 2 flow velocity in y direction

        dirhf[dirhf > 0] <- 1

        intable3 = paste(prefix, "_results/", prefix, "_files/", prefix, "_directions3.txt", 
            sep = "")
        dirxw <- scan(intable3, nlines = 1)  #PHASE 3 x coordinate
        diryw <- scan(intable3, skip = 1, nlines = 1)  #PHASE 3 y coordinate

        dirhw <- scan(intable3, skip = 2 + 3 * dirline, nlines = 1)  #PHASE 3 flow height
        dirvxw <- scan(intable3, skip = 3 + 3 * dirline, nlines = 1)  #PHASE 3 flow velocity in x direction
        dirvyw <- scan(intable3, skip = 4 + 3 * dirline, nlines = 1)  #PHASE 3 flow velocity in y direction

        dirhw[dirhw > 0] <- 1
    }

    dirh[dirh > 0] <- 1
}

# Computing map extent
xdiff <- xmax - xmin  #extent in x direction
ydiff <- ymax - ymin  #extent in y direction

# if ( ydiff > xdiff ) { #ensuring minimum extent in x direction

# xmin<-xmin-(ydiff-xdiff)/2 xmax<-xmax+(ydiff-xdiff)/2 xdiff<-xmax-xmin

# } else

if (xdiff > 1.6 * ydiff) {
    # ensuring minimum extent in y direction

    ymin <- ymin - (xdiff/1.6 - ydiff)/2
    ymax <- ymax + (xdiff/1.6 - ydiff)/2
    ydiff <- ymax - ymin
}

vunit <- 0.05 * min(xdiff, ydiff)/maxvflow

# Building map geometry
asprat <- xdiff/(1.4 * ydiff)  #x/y ratio
margx <- 6  #margin in x direction
margy <- 2.5  #margin in y direction
dispx <- 15  #total width of output image
dispy <- (dispx - margx)/asprat + margy  #total height of output image

# Adapting boundaries for vector layers
xsmin <- xmin + 0.035 * xdiff
xsmax <- xmax - 0.035 * xdiff
ysmin <- ymin + 0.035 * ydiff
ysmax <- ymax - 0.035 * ydiff

# Defining graphic parameters for raster plot
if (model <= 3 || fill == "iscore" || mstring == "_tsun" || mstring == "_treach") {

    if (mstring == "_hflow" || mstring == "_treach" || mstring == "_iii_hflow" || 
        mstring == "_dii") {

        if (fill == "iscore") {
            munit <- ""
            if (mstring == "_iii_hflow") 
                mlabel <- "Impact indicator index"
            if (mstring == "_dii") 
                mlabel <- "Deposition indicator index"
            rcol1 <- rgb(0, 0, 0.8, 0.4)
            rcol7 <- rgb(0.15, 0.1, 0.6, 0.5)
            rcol13 <- rgb(0.3, 0.2, 0.4, 0.6)
            rcol19 <- rgb(0.45, 0.3, 0.2, 0.7)
            rcol25 <- rgb(0.6, 0.4, 0, 0.8)
        } else if (mstring == "_hflow") {
            if (mconv == "1") 
                munit <- "m"
            if (mconv == "0.001") 
                munit <- "km"
            mlabel <- "Flow height"
            rcol1 <- rgb(0.5, 0.5, 0, 0.2)
            rcol7 <- rgb(0.75, 0.25, 0, 0.4)
            rcol13 <- rgb(1, 0, 0, 0.6)
            rcol19 <- rgb(0.8, 0, 0.2, 0.8)
            rcol25 <- rgb(0.6, 0, 0.4, 1)
        } else if (mstring == "_treach") {
            munit <- "s"
            mlabel <- "Time of reach"
            rcol1 <- rgb(0.7, 0.3, 0.1, 1)
            rcol7 <- rgb(0.55, 0.275, 0.2, 0.8)
            rcol13 <- rgb(0.4, 0.25, 0.3, 0.6)
            rcol19 <- rgb(0.25, 0.225, 0.4, 0.4)
            rcol25 <- rgb(0.1, 0.2, 0.5, 0.2)
        }
    } else if (mstring == "_tflow" || mstring == "_iii_tflow") {

        if (fill != "iscore") {
            if (mconv == "1") 
                munit <- "J"
            if (mconv == "0.001") 
                munit <- "kJ"
            if (mconv == "0.000001") 
                munit <- "MJ"
            if (mconv == "0.000000001") 
                munit <- "GJ"
            if (mconv == "0.000000000001") 
                munit <- "TJ"
            mlabel <- "Flow kinetic energy"
            rcol1 <- rgb(0, 0.5, 0.5, 0.2)
            rcol7 <- rgb(0.15, 0.45, 0.4, 0.4)
            rcol13 <- rgb(0.3, 0.4, 0.3, 0.6)
            rcol19 <- rgb(0.45, 0.35, 0.2, 0.8)
            rcol25 <- rgb(0.6, 0.3, 0.1, 1)
        } else {
            munit <- ""
            mlabel <- "Impact indicator index"
            rcol1 <- rgb(0, 0.5, 0.5, 0.4)
            rcol7 <- rgb(0.15, 0.45, 0.4, 0.5)
            rcol13 <- rgb(0.3, 0.4, 0.3, 0.6)
            rcol19 <- rgb(0.45, 0.35, 0.2, 0.7)
            rcol25 <- rgb(0.6, 0.3, 0.1, 0.8)
        }
    } else if (mstring == "_pflow" || mstring == "_iii_pflow") {

        if (fill != "iscore") {
            if (mconv == "1") 
                munit <- "Pa"
            if (mconv == "0.001") 
                munit <- "kPa"
            if (mconv == "0.000001") 
                munit <- "MPa"
            if (mconv == "0.000000001") 
                munit <- "GPa"
            if (mconv == "0.000000000001") 
                munit <- "TPa"
            mlabel <- "Flow pressure"
            rcol1 <- rgb(0.2, 0, 0.6, 0.2)
            rcol7 <- rgb(0.4, 0.3, 0.35, 0.4)
            rcol13 <- rgb(0.4, 0.3, 0.3, 0.6)
            rcol19 <- rgb(0.5, 0.45, 0.15, 0.8)
            rcol25 <- rgb(0.6, 0.6, 0, 1)
        } else {
            munit <- ""
            mlabel <- "Impact indicator index"
            rcol1 <- rgb(0.2, 0, 0.6, 0.4)
            rcol7 <- rgb(0.4, 0.3, 0.35, 0.5)
            rcol13 <- rgb(0.4, 0.3, 0.3, 0.6)
            rcol19 <- rgb(0.5, 0.45, 0.15, 0.7)
            rcol25 <- rgb(0.6, 0.6, 0, 0.8)
        }
    } else if (mstring == "_basechange") {

        if (mconv == "1") 
            munit <- "m"
        if (mconv == "0.001") 
            munit <- "km"
        mlabel <- "Change of basal topography"

        rcol3 <- rgb(0.1, 0.2, 0.6, 1)
        rcol7 <- rgb(0.1, 0.2, 0.6, 0.75)
        rcol11 <- rgb(0.1, 0.2, 0.6, 0.5)
        rcol15 <- rgb(0.1, 0.2, 0.6, 0.25)
        rcol17 <- rgb(1, 1, 1, 0)
        rcol20 <- rgb(0.2, 0.6, 0.1, 0.25)
        rcol24 <- rgb(0.2, 0.6, 0.1, 0.5)
        rcol28 <- rgb(0.2, 0.6, 0.1, 0.75)
        rcol32 <- rgb(0.2, 0.6, 0.1, 1)

    } else if (mstring == "_tsun") {

        if (mconv == "1") 
            munit <- "m"
        if (mconv == "0.001") 
            munit <- "km"
        mlabel <- "Tsunami height"

        rcol3 <- rgb(0, 0.6, 0.1, 1)
        rcol7 <- rgb(0, 0.6, 0.1, 0.75)
        rcol11 <- rgb(0, 0.6, 0.1, 0.5)
        rcol15 <- rgb(0, 0.6, 0.1, 0.25)
        rcol17 <- rgb(1, 1, 1, 0)
        rcol20 <- rgb(0, 0.1, 0.6, 0.25)
        rcol24 <- rgb(0, 0.1, 0.6, 0.5)
        rcol28 <- rgb(0, 0.1, 0.6, 0.75)
        rcol32 <- rgb(0, 0.1, 0.6, 1)
    }

    if (mstring == "_basechange" || mstring == "_tsun") {

        ctext = vector("expression", 10)  #vector for colour bar labels

        if (mstring == "_basechange") {
            ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
        } else {
            ctext[1] <- ""
        }
        ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsmn)), 
            fmunit = format(munit)))[2]
        ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
            fmunit = format(munit)))[2]
        ctext[7] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[8] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[9] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[10] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]

        rthresholds <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)
        rlabels <- ctext
        rcolours <- c(rcol3, rcol7, rcol11, rcol15, rcol17, rcol20, rcol24, rcol28, 
            rcol32)

    } else {

        ctext = vector("expression", 6)  #vector for colour bar labels

        ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
            fmunit = format(munit)))[2]
        ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]

        rthresholds <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
        rlabels <- ctext
        rcolours <- c(rcol1, rcol7, rcol13, rcol19, rcol25)
    }

} else if (model > 3 && model < 7) {

    if (mstring == "_hflow") {

        if (mconv == "1") 
            munit <- "m"
        if (mconv == "0.001") 
            munit <- "km"
        mlabel <- paste("Flow height, Volumetric ratio of", phase1, sep = " ")

        rcol1 <- rgb(0, 0, 0.8, 0.2)
        rcol2 <- rgb(0.15, 0.1, 0.6, 0.2)
        rcol3 <- rgb(0.3, 0.2, 0.4, 0.2)
        rcol4 <- rgb(0.45, 0.3, 0.2, 0.2)
        rcol5 <- rgb(0.6, 0.4, 0, 0.2)
        rcol6 <- rgb(0, 0, 0.8, 0.4)
        rcol7 <- rgb(0.15, 0.1, 0.6, 0.4)
        rcol8 <- rgb(0.3, 0.2, 0.4, 0.4)
        rcol9 <- rgb(0.45, 0.3, 0.2, 0.4)
        rcol10 <- rgb(0.6, 0.4, 0, 0.4)
        rcol11 <- rgb(0, 0, 0.8, 0.6)
        rcol12 <- rgb(0.15, 0.1, 0.6, 0.6)
        rcol13 <- rgb(0.3, 0.2, 0.4, 0.6)
        rcol14 <- rgb(0.45, 0.3, 0.2, 0.6)
        rcol15 <- rgb(0.6, 0.4, 0, 0.6)
        rcol16 <- rgb(0, 0, 0.8, 0.8)
        rcol17 <- rgb(0.15, 0.1, 0.6, 0.8)
        rcol18 <- rgb(0.3, 0.2, 0.4, 0.8)
        rcol19 <- rgb(0.45, 0.3, 0.2, 0.8)
        rcol20 <- rgb(0.6, 0.4, 0, 0.8)
        rcol21 <- rgb(0, 0, 0.8, 1)
        rcol22 <- rgb(0.15, 0.1, 0.6, 1)
        rcol23 <- rgb(0.3, 0.2, 0.4, 1)
        rcol24 <- rgb(0.45, 0.3, 0.2, 1)
        rcol25 <- rgb(0.6, 0.4, 0, 1)

    } else if (mstring == "_tflow") {

        if (mconv == "1") 
            munit <- "J"
        if (mconv == "0.001") 
            munit <- "kJ"
        if (mconv == "0.000001") 
            munit <- "MJ"
        if (mconv == "0.000000001") 
            munit <- "GJ"
        if (mconv == "0.000000000001") 
            munit <- "TJ"
        mlabel <- paste("Flow kinetic energy, Ratio of", phase1, "component", sep = " ")

        rcol1 <- rgb(0, 0.5, 0.5, 0.2)
        rcol2 <- rgb(0.15, 0.45, 0.4, 0.2)
        rcol3 <- rgb(0.3, 0.4, 0.3, 0.2)
        rcol4 <- rgb(0.45, 0.35, 0.2, 0.2)
        rcol5 <- rgb(0.6, 0.3, 0.1, 0.2)
        rcol6 <- rgb(0, 0.5, 0.5, 0.4)
        rcol7 <- rgb(0.15, 0.45, 0.4, 0.4)
        rcol8 <- rgb(0.3, 0.4, 0.3, 0.4)
        rcol9 <- rgb(0.45, 0.35, 0.2, 0.4)
        rcol10 <- rgb(0.6, 0.3, 0.1, 0.4)
        rcol11 <- rgb(0, 0.5, 0.5, 0.6)
        rcol12 <- rgb(0.15, 0.45, 0.4, 0.6)
        rcol13 <- rgb(0.3, 0.4, 0.3, 0.6)
        rcol14 <- rgb(0.45, 0.35, 0.2, 0.6)
        rcol15 <- rgb(0.6, 0.3, 0.1, 0.6)
        rcol16 <- rgb(0, 0.5, 0.5, 0.8)
        rcol17 <- rgb(0.15, 0.45, 0.4, 0.8)
        rcol18 <- rgb(0.3, 0.4, 0.3, 0.8)
        rcol19 <- rgb(0.45, 0.35, 0.2, 0.8)
        rcol20 <- rgb(0.6, 0.3, 0.1, 0.8)
        rcol21 <- rgb(0, 0.5, 0.5, 1)
        rcol22 <- rgb(0.15, 0.45, 0.4, 1)
        rcol23 <- rgb(0.3, 0.4, 0.3, 1)
        rcol24 <- rgb(0.45, 0.35, 0.2, 1)
        rcol25 <- rgb(0.6, 0.3, 0.1, 1)

    } else if (mstring == "_pflow") {

        if (mconv == "1") 
            munit <- "Pa"
        if (mconv == "0.001") 
            munit <- "kPa"
        if (mconv == "0.000001") 
            munit <- "MPa"
        if (mconv == "0.000000001") 
            munit <- "GPa"
        if (mconv == "0.000000000001") 
            munit <- "TPa"
        mlabel <- paste("Flow pressure, Ratio of", phase1, "component", sep = " ")

        rcol1 <- rgb(0.2, 0, 0.6, 0.2)
        rcol2 <- rgb(0.3, 0.15, 0.45, 0.2)
        rcol3 <- rgb(0.4, 0.3, 0.3, 0.2)
        rcol4 <- rgb(0.5, 0.45, 0.15, 0.2)
        rcol5 <- rgb(0.6, 0.6, 0, 0.2)
        rcol6 <- rgb(0.2, 0, 0.6, 0.4)
        rcol7 <- rgb(0.4, 0.3, 0.35, 0.4)
        rcol8 <- rgb(0.3, 0.2, 0.4, 0.4)
        rcol9 <- rgb(0.5, 0.45, 0.15, 0.4)
        rcol10 <- rgb(0.6, 0.6, 0, 0.4)
        rcol11 <- rgb(0.2, 0, 0.6, 0.6)
        rcol12 <- rgb(0.3, 0.15, 0.45, 0.6)
        rcol13 <- rgb(0.4, 0.3, 0.3, 0.6)
        rcol14 <- rgb(0.5, 0.45, 0.15, 0.6)
        rcol15 <- rgb(0.6, 0.6, 0, 0.6)
        rcol16 <- rgb(0.2, 0, 0.6, 0.8)
        rcol17 <- rgb(0.3, 0.15, 0.45, 0.8)
        rcol18 <- rgb(0.4, 0.3, 0.3, 0.8)
        rcol19 <- rgb(0.5, 0.45, 0.15, 0.8)
        rcol20 <- rgb(0.6, 0.6, 0, 0.8)
        rcol21 <- rgb(0.2, 0, 0.6, 1)
        rcol22 <- rgb(0.3, 0.15, 0.45, 1)
        rcol23 <- rgb(0.4, 0.3, 0.3, 1)
        rcol24 <- rgb(0.5, 0.45, 0.15, 1)
        rcol25 <- rgb(0.6, 0.6, 0, 1)

    } else if (mstring == "_basechange") {

        if (mconv == "1") 
            munit <- "m"
        if (mconv == "0.001") 
            munit <- "km"
        mlabel <- paste("Change of basal topography, Volumetric ratio of", phase1, 
            sep = " ")

        rcol1 <- rgb(0, 0.4, 0.4, 1)
        rcol2 <- rgb(0.2, 0.3, 0.3, 1)
        rcol3 <- rgb(0.6, 0.2, 0.1, 1)
        rcol4 <- rgb(0.8, 0.1, 0, 1)
        rcol5 <- rgb(0, 0.4, 0.4, 0.75)
        rcol6 <- rgb(0.2, 0.3, 0.3, 0.75)
        rcol7 <- rgb(0.6, 0.2, 0.1, 0.75)
        rcol8 <- rgb(0.8, 0.1, 0, 0.75)
        rcol9 <- rgb(0, 0.4, 0.4, 0.5)
        rcol10 <- rgb(0.2, 0.3, 0.3, 0.5)
        rcol11 <- rgb(0.6, 0.2, 0.1, 0.5)
        rcol12 <- rgb(0.8, 0.1, 0, 0.5)
        rcol13 <- rgb(0, 0.4, 0.4, 0.25)
        rcol14 <- rgb(0.2, 0.3, 0.3, 0.25)
        rcol15 <- rgb(0.6, 0.2, 0.1, 0.25)
        rcol16 <- rgb(0.8, 0.1, 0, 0.25)
        rcol17 <- rgb(1, 1, 1, 0)
        rcol18 <- rgb(0, 0.2, 0.9, 0.25)
        rcol19 <- rgb(0.3, 0.3, 0.6, 0.25)
        rcol20 <- rgb(0.6, 0.3, 0.3, 0.25)
        rcol21 <- rgb(0.9, 0.4, 0, 0.25)
        rcol22 <- rgb(0, 0.2, 0.9, 0.5)
        rcol23 <- rgb(0.3, 0.3, 0.6, 0.5)
        rcol24 <- rgb(0.6, 0.3, 0.3, 0.5)
        rcol25 <- rgb(0.9, 0.4, 0, 0.5)
        rcol26 <- rgb(0, 0.2, 0.9, 0.75)
        rcol27 <- rgb(0.3, 0.3, 0.6, 0.75)
        rcol28 <- rgb(0.6, 0.3, 0.3, 0.75)
        rcol29 <- rgb(0.9, 0.4, 0, 0.75)
        rcol30 <- rgb(0, 0.2, 0.9, 1)
        rcol31 <- rgb(0.3, 0.3, 0.6, 1)
        rcol32 <- rgb(0.6, 0.3, 0.3, 1)
        rcol33 <- rgb(0.9, 0.4, 0, 1)
    }

    if (mstring == "_basechange") {

        ctext = vector("expression", 10)  #vector for colour bar labels

        ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsmn)), 
            fmunit = format(munit)))[2]
        ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
            fmunit = format(munit)))[2]
        ctext[7] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[8] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[9] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[10] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]

        rthresholds <- c(-16.5, -15.5, -14.5, -13.5, -12.5, -11.5, -10.5, -9.5, -8.5, 
            -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 
            5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5)
        rlabels = c(ctext[1], " 0.25", " 0.50", " 0.75", ctext[2], " 0.25", " 0.50", 
            " 0.75", ctext[3], " 0.25", " 0.50", " 0.75", ctext[4], " 0.25", " 0.50", 
            " 0.75", ctext[5], ctext[6], " 0.25", " 0.50", " 0.75", ctext[7], " 0.25", 
            " 0.50", " 0.75", ctext[8], " 0.25", " 0.50", " 0.75", ctext[9], " 0.25", 
            " 0.50", " 0.75", ctext[10])
        rcolours <- c(rcol1, rcol2, rcol3, rcol4, rcol5, rcol6, rcol7, rcol8, rcol9, 
            rcol10, rcol11, rcol12, rcol13, rcol14, rcol15, rcol16, rcol17, rcol18, 
            rcol19, rcol20, rcol21, rcol22, rcol23, rcol24, rcol25, rcol26, rcol27, 
            rcol28, rcol29, rcol30, rcol31, rcol32, rcol33)

    } else {

        ctext = vector("expression", 6)  #vector for colour bar labels

        ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
            fmunit = format(munit)))[2]
        ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]

        rthresholds <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 
            11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 
            23.5, 24.5, 25.5)
        rlabels = c(ctext[1], " 0.2", " 0.4", " 0.6", " 0.8", ctext[2], " 0.2", " 0.4", 
            " 0.6", " 0.8", ctext[3], " 0.2", " 0.4", " 0.6", " 0.8", ctext[4], " 0.2", 
            " 0.4", " 0.6", " 0.8", ctext[5], " 0.2", " 0.4", " 0.6", " 0.8", ctext[6])
        rcolours <- c(rcol1, rcol2, rcol3, rcol4, rcol5, rcol6, rcol7, rcol8, rcol9, 
            rcol10, rcol11, rcol12, rcol13, rcol14, rcol15, rcol16, rcol17, rcol18, 
            rcol19, rcol20, rcol21, rcol22, rcol23, rcol24, rcol25)
    }
} else if (model == 7) {

    if (mstring == "_hflow") {

        if (mconv == "1") 
            munit <- "m"
        if (mconv == "0.001") 
            munit <- "km"
        mlabel <- "Flow height"

        rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
            0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
        rthresholds <- rep(0:625) + 0.5

    } else if (mstring == "_tflow") {

        if (mconv == "1") 
            munit <- "J"
        if (mconv == "0.001") 
            munit <- "kJ"
        if (mconv == "0.000001") 
            munit <- "MJ"
        if (mconv == "0.000000001") 
            munit <- "GJ"
        if (mconv == "0.000000000001") 
            munit <- "TJ"
        mlabel <- "Flow kinetic\nenergy"

        rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
            0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
        rthresholds <- rep(0:625) + 0.5

    } else if (mstring == "_pflow") {

        if (mconv == "1") 
            munit <- "Pa"
        if (mconv == "0.001") 
            munit <- "kPa"
        if (mconv == "0.000001") 
            munit <- "MPa"
        if (mconv == "0.000000001") 
            munit <- "GPa"
        if (mconv == "0.000000000001") 
            munit <- "TPa"
        mlabel <- "Flow pressure"

        rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
            0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
        rthresholds <- rep(0:625) + 0.5

    } else if (mstring == "_basechange") {

        if (mconv == "1") 
            munit <- "m"
        if (mconv == "0.001") 
            munit <- "km"
        mlabel <- "Change of\ntopography"

        rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
            0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
        rthresholds <- rep(0:625) + 0.5

    } else if (mstring == "_tsun") {

        if (mconv == "1") 
            munit <- "m"
        if (mconv == "0.001") 
            munit <- "km"
        mlabel <- "Height of tsunami"

        rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
            0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
        rthresholds <- rep(0:625) + 0.5
    }

    if (mstring == "_basechange") {

        ctext = vector("expression", 6)  #vector for contour line legend labels

        ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
            fmunit = format(munit)))[2]
        ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]

        ctext_neg = vector("expression", 6)  #vector for negative contour line legend labels

        ctext_neg[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrs5n)), 
            fmunit = format(munit)))[2]
        ctext_neg[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext_neg[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext_neg[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext_neg[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1n), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext_neg[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrsmn), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]

    } else {

        ctext = vector("expression", 6)  #vector for contour line legend labels

        ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
            fmunit = format(munit)))[2]
        ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
        ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), 
            mdig), nsmall = mdig), fmunit = format(munit)))[2]
    }
}


# Importing raster layers
hflowmap <- paste(temppath, "/a2.asc", sep = "")  #reclassified parameter map
rastplot <- raster(hflowmap)  #raster

if (ortho != "0") orthophoto <- stack(ortho)  #orthophoto
hillshmap = paste(temppath, "/hillshade.asc", sep = "")  #hillshade
hillshade <- raster(hillshmap)  #hillshade

# Creating plot file
if (fill == "iscore") {
    mapplot = paste(prefix, "_results/", prefix, "_plots/", prefix, mstring, ".png", 
        sep = "")
} else if (ntimesteps <= ntimemax) {
    mapplot = paste(prefix, "_results/", prefix, "_plots/", prefix, "_maps_timesteps/", 
        prefix, mstring, fill, ".png", sep = "")
} else if (mstring == "_basechange") {
    mapplot = paste(prefix, "_results/", prefix, "_plots/", prefix, mstring, "_fin.png", 
        sep = "")
} else if (mstring == "_tsun") {
    mapplot = paste(prefix, "_results/", prefix, "_plots/", prefix, mstring, "_max.png", 
        sep = "")
} else if (mstring == "_treach") {
    mapplot = paste(prefix, "_results/", prefix, "_plots/", prefix, mstring, ".png", 
        sep = "")
} else {
    mapplot = paste(prefix, "_results/", prefix, "_plots/", prefix, mstring, "_max.png", 
        sep = "")
}
png(filename = mapplot, width = dispx, height = dispy, units = "cm", res = 300)

# Defining margins
par(mar = c(2, 1.3, 0.5, 2.5))
par(oma = c(0, 0, 0, 1.5))
clip(xmin, xmax, ymin, ymax)  #constraining drawing area

if (ydiff > xdiff * 1.2) {
    # shrink factor and scaling of labels for colour bar legend
    lshrink <- 0.7
    lcex <- 1
} else if (ydiff > xdiff/1.2) {
    lshrink <- 0.85
    lcex <- 1
} else {
    lshrink <- 1
    lcex <- 0.8
}

# Plotting raster layers
if (ortho != 0) {
    plot(hillshade, legend.width = 1, legend = FALSE, col = rgb(0.75, 0.75, 0.75, 
        1), axes = FALSE, box = FALSE, xlab = NA, ylab = NA, xlim = c(xmin, xmax), 
        ylim = c(ymin, ymax), useRaster = TRUE)  #hillshade
    plotRGB(orthophoto, add = T, r = 3, g = 2, b = 1, stretch = "hist", alpha = 150)  #orthophoto
} else {
    plot(hillshade, legend.width = 1, col = gray((max(0, cellStats(hillshade, "min")):min(255, 
        cellStats(hillshade, "max")))/255), legend = FALSE, axes = FALSE, box = FALSE, 
        xlab = NA, ylab = NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), useRaster = TRUE)  #hillshade
}

par(new = TRUE)
par(cex = lcex)

if (model < 7 || fill == "iscore" || mstring == "_tsun" || mstring == "_treach") {
    plot(rastplot, legend.width = 1, useRaster = TRUE, col = rcolours, breaks = rthresholds, 
        xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE, xlab = NA, ylab = NA, 
        box = FALSE, legend.args = list(text = mlabel, side = 4, line = -2.1, cex = lcex), 
        axis.args = list(labels = rlabels), legend.shrink = lshrink)  #flow parameter
} else {
    plot(rastplot, legend.width = 1, useRaster = TRUE, col = rcolours, breaks = rthresholds, 
        xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE, xlab = NA, ylab = NA, 
        box = FALSE, legend = FALSE)  #flow parameter
}

par(cex = 1)

# Plotting vector layers
if (model <= 3 && releasemass > 0) {
    startshape = readOGR(temppath, layer = "d_lrelease")
    par(new = TRUE)
    plot(startshape, bg = "transparent", col = "red", lty = 3, lwd = 1.5, xlim = c(xsmin, 
        xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, ylab = NA)
    # release area
} else if (model > 3 && model < 7) {
    if (releasemass > 1) {
        startshape2 = readOGR(temppath, layer = "d_lrelease2")
        par(new = TRUE)
        plot(startshape2, bg = "transparent", col = "blue", lty = 3, lwd = 1.5, xlim = c(xsmin, 
            xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, ylab = NA)
        # PHASE 2 release area
    }
    if (releasemass == 1 || releasemass == 3) {
        startshape = readOGR(temppath, layer = "d_lrelease")
        par(new = TRUE)
        plot(startshape, bg = "transparent", col = "brown", lty = 3, lwd = 1.5, xlim = c(xsmin, 
            xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, ylab = NA)
        # PHASE 1 release area
    }
} else if (model == 7) {
    if (releasemass > 3) {
        startshape3 = readOGR(temppath, layer = "d_lrelease3")
        par(new = TRUE)
        plot(startshape3, bg = "transparent", col = rgb(0, 0, 1), lty = 3, lwd = 1.5, 
            xlim = c(xsmin, xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, 
            ylab = NA)
        # PHASE 3 release area
    }
    if (releasemass == 2 || releasemass == 3 || releasemass == 7) {
        startshape2 = readOGR(temppath, layer = "d_lrelease2")
        par(new = TRUE)
        plot(startshape2, bg = "transparent", col = rgb(0, 1, 0), lty = 3, lwd = 1.5, 
            xlim = c(xsmin, xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, 
            ylab = NA)
        # PHASE 2 release area
    }
    if (releasemass == 1 || releasemass == 3 || releasemass == 5 || releasemass == 
        7) {
        startshape = readOGR(temppath, layer = "d_lrelease")
        par(new = TRUE)
        plot(startshape, bg = "transparent", col = rgb(1, 0, 0), lty = 3, lwd = 1.5, 
            xlim = c(xsmin, xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, 
            ylab = NA)
        # PHASE 1 release area
    }
}

if (impdef == 1 && (fill != "iscore" || mlabel == "Impact indicator index")) {

    par(new = TRUE)
    impactareashape = readOGR(temppath, layer = "d_limpactarea")
    plot(impactareashape, bg = "transparent", col = "red", lty = 1, lwd = 1.5, xlim = c(xsmin, 
        xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, ylab = NA)
    # observed impact area
}

if (depdef == 1 && (fill != "iscore" || mlabel == "Deposition indicator index")) {

    par(new = TRUE)
    depositshape = readOGR(temppath, layer = "d_ldeposit")
    plot(depositshape, bg = "transparent", col = "orange", lty = 1, lwd = 1.5, xlim = c(xsmin, 
        xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, ylab = NA)
    # observed deposit
}

if (mapprof == 1) {

    par(new = TRUE)
    intable <- paste(temppath, "/fmapprof.txt", sep = "")
    profilex <- scan(intable, nlines = 1)  #profile x coordinate
    profiley <- scan(intable, skip = 1, nlines = 1)  #profile y coordinate
    lines(x = profilex, y = profiley, col = "yellow", lwd = 1.5, lty = 2)  #profile
}

if (cellStats(rastplot, "max") != cellStats(rastplot, "min")) {

    par(new = TRUE)
    contourshape = try(readOGR(temppath, layer = "d_contours"))
    try(plot(contourshape, bg = "transparent", col = "black", lty = 1, lwd = 0.5, 
        xlim = c(xsmin, xsmax), ylim = c(ysmin, ysmax), axes = F, xlab = NA, ylab = NA))  #contour lines of flow height

    if (model == 7 && mstring == "_basechange") {

        par(new = TRUE)
        contourshape_neg = try(readOGR(temppath, layer = "d_contours_neg"))
        try(plot(contourshape_neg, bg = "transparent", col = "gainsboro", lty = 1, 
            lwd = 0.5, xlim = c(xsmin, xsmax), ylim = c(ysmin, ysmax), axes = F, 
            xlab = NA, ylab = NA))  #contour lines of flow height

    }
}

if (velocity == 1) {

    if (model <= 3) {
        par(new = TRUE)
        arrows(x0 = dirx, y0 = diry, x1 = dirx + vunit * dirvy, y1 = diry - vunit * 
            dirvx, length = 0.025, code = 2, angle = 45, col = "red", lwd = 0.5 * 
            dirh)  #directions
    } else if (model > 3 && model < 7) {
        par(new = TRUE)
        arrows(x0 = dirxf, y0 = diryf, x1 = dirxf + vunit * dirvyf, y1 = diryf - 
            vunit * dirvxf, length = 0.025, code = 2, angle = 45, col = "blue", lwd = 0.5 * 
            dirhf)  #PHASE 2 directions
        par(new = TRUE)
        arrows(x0 = dirx, y0 = diry, x1 = dirx + vunit * dirvy, y1 = diry - vunit * 
            dirvx, length = 0.025, code = 2, angle = 45, col = "brown", lwd = 0.5 * 
            dirh)  #PHASE 1 directions
    } else if (model == 7) {
        par(new = TRUE)
        arrows(x0 = dirxw, y0 = diryw, x1 = dirxw + vunit * dirvyw, y1 = diryw - 
            vunit * dirvxw, length = 0.025, code = 2, angle = 45, col = rgb(0, 0, 
            1), lwd = 0.5 * dirhw)  #PHASE 3 directions
        par(new = TRUE)
        arrows(x0 = dirxf, y0 = diryf, x1 = dirxf + vunit * dirvyf, y1 = diryf - 
            vunit * dirvxf, length = 0.025, code = 2, angle = 45, col = rgb(0, 1, 
            0), lwd = 0.5 * dirhf)  #PHASE 2 directions
        par(new = TRUE)
        arrows(x0 = dirx, y0 = diry, x1 = dirx + vunit * dirvy, y1 = diry - vunit * 
            dirvx, length = 0.025, code = 2, angle = 45, col = rgb(1, 0, 0), lwd = 0.5 * 
            dirh)  #PHASE 1 directions
    }
}

# Plotting control points
if (ctrlpts > 0) {

    par(new = TRUE)
    inctrl <- paste(prefix, "_results/", prefix, "_files/", prefix, "_ctrlpoints.txt", 
        sep = "")

    j <- 1  #initializing counter for loop over all control points
    repeat {
        # starting loop over all control points
        if (j > ctrlpts) {
            break  #break if last control point was reached
        } else {
            ctrlx <- read.table(inctrl, header = FALSE)[(j - 1) * (ntimemax + 2) + 
                2, 3]  #control point x coordinate
            ctrly <- read.table(inctrl, header = FALSE)[(j - 1) * (ntimemax + 2) + 
                2, 4]  #control point y coordinate
            points(x = as.vector(t(ctrlx)), y = as.vector(t(ctrly)), col = "yellow", 
                pch = 3)  #plotting control point
            ctrltext <- as.character(read.table(inctrl, header = FALSE)[(j - 1) * 
                (ntimemax + 2) + 2, 1])

            text(x = as.numeric(as.vector(t(ctrlx))) + 0.035 * (xmax - xmin), y = as.numeric(as.vector(t(ctrly))) + 
                0.035 * (xmax - xmin), labels = ctrltext, cex = 1, col = "yellow", 
                font = 1)  #control point label

            j <- j + 1
        }
    }
}

# Plotting hydrograph profiles
if (hydrograph > 0) {

    par(new = TRUE)
    inhyd <- paste(prefix, "_results/", prefix, "_files/", prefix, "_hydprofiles.txt", 
        sep = "")

    j <- 1  #initializing counter for loop over all hydrographs
    repeat {
        # starting loop over all hydrographs
        if (j > (ninhyd + nouthyd)) {
            break  #break if last hydrograph was reached
        } else {
            if (j <= ninhyd) {
                # color of hydrograph profile depending on type:
                hydcol <- "green"
            } else {
                hydcol <- "purple"
            }

            hydx <- read.table(inhyd, header = FALSE)[j + 1, 2:4]  #hydrograph x coordinates
            hydy <- read.table(inhyd, header = FALSE)[j + 1, 5:7]  #hydrograph y coordinates
            lines(x = as.vector(t(hydx)), y = as.vector(t(hydy)), col = hydcol, lwd = 1.5, 
                lty = 3)  #profile
            points(x = as.vector(t(hydx))[2], y = as.vector(t(hydy))[2], col = hydcol, 
                pch = 19, cex = 0.5)  #centre point
            hydtext <- as.character(read.table(inhyd, header = FALSE)[j + 1, 1])

            hyddx1 <- as.numeric(as.vector(t(hydx))[1])  #coordinates for hydrograph label:
            hyddx3 <- as.numeric(as.vector(t(hydx))[3])
            if (hyddx1 > hyddx3) {
                hyddx <- hyddx1
                hyddy <- as.numeric(as.vector(t(hydy))[1])
            } else {
                hyddx <- hyddx3
                hyddy <- as.numeric(as.vector(t(hydy))[3])
            }

            text(x = hyddx, y = hyddy, labels = hydtext, cex = 1, col = hydcol, font = 1, 
                pos = 4, offset = 0.3)  #hydrograph label

            j <- j + 1
        }
    }
}

# Plotting bounding box and axes
box()
par(cex = 0.8)
axis(side = 1, tck = -0.01, labels = NA)  #x axis
axis(side = 2, tck = -0.01, labels = NA)  #y axis
axis(side = 1, lwd = 0, line = -0.6)  #x axis
axis(side = 2, lwd = 0, line = -0.6)  #y axis
par(cex = 1)

# Plotting header text
htext = vector("expression", 1)

if (fill != "iscore") {
    # defining x position and label:
    if (ntimesteps <= ntimemax) {
        secpassed <- tint * ntimesteps  #time passed since start of simulation
        if (secpassed > tstop) 
            secpassed <- tstop
        htext[1] <- substitute(expression(italic(t) == secformat ~ s), list(secformat = format(round(secpassed, 
            1), nsmall = 1)))[2]
        posx <- xmin + xdiff/5
    } else if (mstring == "_hflow") {
        htext[1] <- "Maximum flow height"
        posx <- xmin + xdiff/2
    } else if (mstring == "_tflow") {
        htext[1] <- "Maximum flow kinetic energy"
        posx <- xmin + xdiff/2
    } else if (mstring == "_pflow") {
        htext[1] <- "Maximum flow pressure"
        posx <- xmin + xdiff/2
    } else if (mstring == "_basechange") {
        htext[1] <- "Final change of basal topography"
        posx <- xmin + xdiff/2
    } else if (mstring == "_tsun") {
        htext[1] <- "Maximum tsunami height"
        posx <- xmin + xdiff/2
    } else if (mstring == "_treach") {
        htext[1] <- "Time of reach"
        posx <- xmin + xdiff/2
    }
} else if (mstring == "_iii_hflow") {
    htext[1] <- "Impact indicator index"
    posx <- xmin + xdiff/2
} else if (mstring == "_iii_tflow") {
    htext[1] <- "Impact indicator index"
    posx <- xmin + xdiff/2
} else if (mstring == "_iii_pflow") {
    htext[1] <- "Impact indicator index"
    posx <- xmin + xdiff/2
} else if (mstring == "_dii") {
    htext[1] <- "Deposition indicator index"
    posx <- xmin + xdiff/2
}

if (ydiff > xdiff * 1.2) {
    # y positions and scaling of header text and header legend
    posy1 <- ymax + ydiff/14
    posy2 <- ymax + ydiff/10
    posy3 <- ymax + ydiff/14
    posy4 <- ymax + ydiff/26
    lcex <- 1
} else if (ydiff > xdiff/1.2) {
    posy1 <- ymax + ydiff/12
    posy2 <- ymax + ydiff/8
    posy3 <- ymax + ydiff/13
    posy4 <- ymax + ydiff/34
    lcex <- 1
} else {
    posy1 <- ymax + ydiff/11
    posy2 <- ymax + ydiff/7
    posy3 <- ymax + ydiff/11.5
    posy4 <- ymax + ydiff/27
    lcex <- 0.8
}

text(x = posx, y = posy1, labels = htext[1], cex = 1.4 * lcex, col = "black")  #printing text

# Plotting header legend (flow velocity)
if (fill != "iscore" && ntimesteps <= ntimemax) {

    par(new = TRUE)
    if (model <= 3) {
        arrows(x0 = xmin + 0.62 * xdiff, y0 = posy3, x1 = xmin + 0.62 * xdiff + 0.05 * 
            min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, 
            col = "red", lwd = 0.5)  #legend for velocity
    } else if (model > 3 && model < 7) {
        arrows(x0 = xmin + 0.62 * xdiff, y0 = posy3, x1 = xmin + 0.62 * xdiff + 0.05 * 
            min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, 
            col = "brown", lwd = 0.5)  #legend for PHASE 1 velocity
        arrows(x0 = xmin + 0.62 * xdiff + 0.05 * min(xdiff, ydiff), y0 = posy4, x1 = xmin + 
            0.62 * xdiff, y1 = posy4, length = 0.025, code = 2, angle = 45, col = "blue", 
            lwd = 0.5)  #legend for PHASE 2 velocity
    } else if (model == 7) {
        arrows(x0 = xmin + 0.58 * xdiff, y0 = posy3, x1 = xmin + 0.58 * xdiff + 0.05 * 
            min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, 
            col = rgb(1, 0, 0), lwd = 0.5)  #legend for PHASE 1 velocity
        arrows(x0 = xmin + 0.58 * xdiff, y0 = posy4, x1 = xmin + 0.58 * xdiff + 0.05 * 
            min(xdiff, ydiff), y1 = posy4, length = 0.025, code = 2, angle = 45, 
            col = rgb(0, 1, 0), lwd = 0.5)  #legend for PHASE 2 velocity
        arrows(x0 = xmin + 0.68 * xdiff, y0 = (posy3 + posy4)/2, x1 = xmin + 0.68 * 
            xdiff + 0.05 * min(xdiff, ydiff), y1 = (posy3 + posy4)/2, length = 0.025, 
            code = 2, angle = 45, col = rgb(0, 0, 1), lwd = 0.5)  #legend for PHASE 3 velocity
    }

    textarrow = vector("expression", 1)
    textarrow[1] = substitute(expression(italic(v)[max] == vformat ~ m/s), list(vformat = format(round(maxvflow, 
        1), nsmall = 1)))[2]  #defining label for flow velocity legend
    text(x = xmax - xdiff/2.15, y = posy2, labels = textarrow[1], col = "black", 
        pos = 4, cex = lcex)

    if (model > 3 && model < 7) {
        text(x = xmax - xdiff/2.95, y = posy3, labels = "Solid", col = "brown", pos = 4, 
            cex = lcex)  #printing legend text for PHASE 1 velocity
        text(x = xmax - xdiff/2.95, y = posy4, labels = "Fluid", col = "blue", pos = 4, 
            cex = lcex)  #printing legend text for PHASE 2 velocity
    } else if (model == 7) {
        text(x = xmin + 0.61 * xdiff, y = posy3, labels = "P1", col = rgb(1, 0, 0), 
            pos = 4, cex = lcex)  #printing legend text for PHASE 1 velocity
        text(x = xmin + 0.61 * xdiff, y = posy4, labels = "P2", col = rgb(0, 1, 0), 
            pos = 4, cex = lcex)  #printing legend text for PHASE 2 velocity
        text(x = xmin + 0.71 * xdiff, y = (posy3 + posy4)/2, labels = "P3", col = rgb(0, 
            0, 1), pos = 4, cex = lcex)  #printing legend text for PHASE 3 velocity
    }
}

# Plotting footer legend (release heights, hydrographs and observed impact area
# and deposit) x positions:
if (releasemass > 0) {
    if (model < 7) {
        pos1 <- 0.035
    } else {
        pos1 <- 0.01
        pos1w <- 0.125
    }
    pos2 <- 0.335
} else {
    pos2 <- 0.035
}
pos3 <- 0.635

if (ydiff > xdiff * 1.2) {
    # y positions and scaling of labels
    posy1 <- 0.05
    posy2 <- 0.025
    posy3 <- 0
    lcex <- 1
} else if (ydiff > xdiff/1.2) {
    posy1 <- 0.066
    posy2 <- 0.026
    posy3 <- -0.014
    lcex <- 1
} else {
    posy1 <- 0.071
    posy2 <- 0.029
    posy3 <- -0.013
    lcey <- 0.8
}

if (releasemass > 0) {
    if (model <= 3) {
        legend("bottomleft", lwd = NA, legend = "Release area", text.col = "black", 
            bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1 - 0.1, posy1), 
            cex = lcex)
        legend("bottomleft", lty = 3, lwd = 1.5, col = "red", legend = "", text.col = "black", 
            bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1, posy2), cex = lcex)  #legend for release area
    } else if (model > 3 && model < 7) {
        legend("bottomleft", lwd = NA, legend = "Release areas", text.col = "black", 
            bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1 - 0.1, posy1), 
            cex = lcex)
        legend("bottomleft", lty = 3, lwd = 1.5, col = "brown", legend = "Solid", 
            text.col = "brown", bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1, 
                posy2), cex = lcex)  #legend for PHASE 1 release area
        legend("bottomleft", lty = 3, lwd = 1.5, col = "blue", legend = "Fluid", 
            text.col = "blue", bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1, 
                posy3), cex = lcex)  #legend for PHASE 2 release area
    } else if (model == 7) {
        legend("bottomleft", lwd = NA, legend = "Release areas", text.col = "black", 
            bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1 - 0.075, posy1), 
            cex = lcex)
        legend("bottomleft", lty = 3, lwd = 1.5, col = rgb(1, 0, 0), legend = "P1", 
            text.col = "red", bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1, 
                posy2), cex = lcex)  #legend for PHASE 1 release area
        legend("bottomleft", lty = 3, lwd = 1.5, col = rgb(0, 1, 0), legend = "P2", 
            text.col = "green", bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1w, 
                (posy2 + posy3)/2), cex = lcex)  #legend for PHASE 2 release area
        legend("bottomleft", lty = 3, lwd = 1.5, col = rgb(0, 0, 1), legend = "P3", 
            text.col = "blue", bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos1, 
                posy3), cex = lcex)  #legend for PHASE 3 release area
    }
}

if (model == 7 && mstring == "_basechange" && fill != "iscore") {
    par(xpd = TRUE)

    legend(x = xmax + xdiff/50, y = ymax - ydiff/25, xjust = 0, title = mlabel, legend = NA, 
        bty = "n", horiz = FALSE, cex = lcex)  #legend title

    legend(x = xmax, y = ymax - ydiff/10, xjust = 0, legend = c(ctext), lty = c(1, 
        1, 1, 1, 1, 1), lwd = 0.5, col = c("white", "black", "black", "black", "black", 
        "black"), text.col = "black", bty = "n", horiz = FALSE, cex = lcex)  #main legend

    legend(x = xmax, y = ymax - ydiff/2, xjust = 0, legend = c(ctext_neg), lty = c(1, 
        1, 1, 1, 1, 1), lwd = 0.5, col = c("gainsboro", "gainsboro", "gainsboro", 
        "gainsboro", "gainsboro", "white"), text.col = "black", bty = "n", horiz = FALSE, 
        cex = lcex)  #main legend for negative values

    legend(x = xmax, y = ymax - ydiff/1.1, xjust = 0, legend = c(paste("P1", phase1, 
        sep = " "), paste("P2", phase2, sep = " "), paste("P3", phase3, sep = " ")), 
        fill = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), border = "black", 
        text.col = "black", bty = "n", horiz = FALSE, cex = lcex)  #legend for fraction

    # legend(x=xmax+xdiff/50, y=ymax-ydiff/9, xjust=0, title='Max', legend=NA,
    # bty='n', horiz=FALSE, cex=lcex) #maximum legend(x=xmax+xdiff/50,
    # y=ymax-ydiff/1.3275, xjust=0, title='Min', legend=NA, bty='n', horiz=FALSE,
    # cex=lcex) #minimum

} else if (model == 7 && mstring != "_tsun" && mstring != "_treach" && fill != "iscore") {

    par(xpd = TRUE)

    legend(x = xmax + xdiff/50, y = ymax - ydiff/25, xjust = 0, title = mlabel, legend = NA, 
        bty = "n", horiz = FALSE, cex = lcex)  #legend title

    legend(x = xmax, y = ymax - ydiff/10, xjust = 0, legend = c(ctext), lty = c(1, 
        1, 1, 1, 1, 1), lwd = 0.5, col = c("white", "black", "black", "black", "black", 
        "black"), text.col = "black", bty = "n", horiz = FALSE, cex = lcex)  #main legend

    legend(x = xmax, y = ymax - ydiff/1.75, xjust = 0, legend = c(paste("P1", phase1, 
        sep = " "), paste("P2", phase2, sep = " "), paste("P3", phase3, sep = " ")), 
        fill = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), border = "black", 
        text.col = "black", bty = "n", horiz = FALSE, cex = lcex)  #legend for fraction

    legend(x = xmax + xdiff/50, y = ymax - ydiff/9, xjust = 0, title = "Max", legend = NA, 
        bty = "n", horiz = FALSE, cex = lcex)  #maximum
}

if (hydrograph > 0) {
    legend("bottomleft", lwd = NA, legend = "Hydrographs", text.col = "black", bty = "n", 
        horiz = TRUE, x.intersp = 0.25, inset = c(pos2 - 0.1, posy1), cex = lcex)
    legend("bottomleft", lty = 3, lwd = 1.5, pch = 19, pt.cex = 0.5, col = "green", 
        legend = "Input", text.col = "green", bty = "n", horiz = TRUE, x.intersp = 0.25, 
        inset = c(pos2, posy2), cex = lcex)  #legend for input hydrograph
    legend("bottomleft", lty = 3, lwd = 1.5, pch = 19, pt.cex = 0.5, col = "purple", 
        legend = "Output", text.col = "purple", bty = "n", horiz = TRUE, x.intersp = 0.25, 
        inset = c(pos2, posy3), cex = lcex)  #legend for output hydrograph
}
if (depdef == 1 || impdef == 1) {
    legend("bottomleft", lwd = NA, legend = "Observation", text.col = "black", bty = "n", 
        horiz = TRUE, x.intersp = 0.25, inset = c(pos3 - 0.1, posy1), cex = lcex)
    legend("bottomleft", lty = 1, lwd = 1.5, col = "red", legend = "Impact area", 
        text.col = "red", bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos3, 
            posy2), cex = lcex)  #legend for observed impact area
    legend("bottomleft", lty = 1, lwd = 1.5, col = "orange", legend = "Deposit", 
        text.col = "orange", bty = "n", horiz = TRUE, x.intersp = 0.25, inset = c(pos3, 
            posy3), cex = lcex)  #legend for observed deposit
}

# Closing plot file
dev.off()
