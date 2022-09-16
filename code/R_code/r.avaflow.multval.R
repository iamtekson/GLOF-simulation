#!/usr/bin/R

##############################################################################
#
# MODULE:       r.avaflow.multval.R
# AUTHOR:       Martin Mergili
#
# PURPOSE:      Script for systematic analysis
#               of multiple evaluation parameters
#
# COPYRIGHT:    (c) 2016 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2016 - 2020 by the BOKU University, Vienna
#               (c) 2016 - 2020 by the University of Vienna
#               (c) 1993 - 2021 by the R Development Core Team
#
# VERSION:      20210525 (25 May 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2).
#
##############################################################################

# Loading library
library(stats)

# Importing arguments (defined in r.ranger.py)
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
ipar1 <- as.integer(args[2])
ipar2 <- as.integer(args[3])
obstype <- as.character(args[4])
ictrlpoint <- as.character(args[5])

# Defining input parameter file name
inparam = paste(prefix, "_results/", prefix, "_files/", prefix, "_params.txt", sep = "")

# Defining evaluation parameter file name
invalid = paste(prefix, "_results/", prefix, "_files/", prefix, "_evaluation.txt", 
    sep = "")

# Creating vectors from input parameters
par1 <- read.table(inparam, skip = 1)[, ipar1 + 1]  #1st parameter
par2 <- read.table(inparam, skip = 1)[, ipar2 + 1]  #2nd parameter

# Building geometry of profile plot
maxx <- max(par1, na.rm = TRUE)  #maximum of horizontal coordinate
minx <- min(par1, na.rm = TRUE)  #minimum of horizontal coordinate
maxy <- max(par2, na.rm = TRUE)  #maximum of vertical coordinate
miny <- min(par2, na.rm = TRUE)  #minimum of vertical coordinate

# Defining vectors characterizing each evaluation parameter
if (obstype == "i") {
    coldata <- c(16, 23, 24, 26, 27, 28)
} else if (obstype == "d") {
    coldata <- c(18, 33, 34, 36, 37, 38)
} else if (obstype == "t") {
    coldata <- c(38 + 3 * as.numeric(ictrlpoint))
}

if (obstype == "t") {

    colname <- c(paste("treachrat", ictrlpoint, sep = ""))
    coltext <- c("Time of reach ratio")
    mdig <- c(2)
    unit <- c("")
    ngraphs <- 1

} else {

    colname <- c("Lratio", "CSI", "HSS", "D2PC", "FoC", "SPI")
    coltext <- c("Excess travel distance ratio", "Critical success index", "Heidke skill score", 
        "Distance to perfect classification", "Factor of conservativeness", "Synthetic Performance Index")
    mdig <- c(2, 2, 2, 2, 2, 2)
    unit <- c("", "", "", "", "", "")
    ngraphs <- 6
}

# Starting loop over all evaluation parameters
for (i in 1:ngraphs) {

    # Reading evaluation parameter-specific data from file
    data <- read.table(invalid, skip = 1)[, coldata[i]]

    # Creating plot file
    profileplot = paste(prefix, "_results/", prefix, "_plots/", prefix, "_multval_", 
        obstype, "_", colname[i], ".png", sep = "")
    png(filename = profileplot, width = 10, height = 10, units = "cm", res = 300)

    # Defining margins
    par(mar = c(3.2, 3.2, 3.2, 3.2))

    # Defining thresholds and colours
    if (obstype == "t") {

        data1 <- 0.4
        data2 <- 0.7
        data3 <- 1
        data4 <- 1.4
        data5 <- 2.5

        plotblue <- 2 * log10(data)
        plotred <- -2 * log10(data)
        plotblue[plotblue < 0] <- 0
        plotblue[plotblue > 1] <- 1
        plotred[plotred < 0] <- 0
        plotred[plotred > 1] <- 1
        plotgreen <- 1 - plotred - plotblue

        colleg <- c(rgb(0.796, 0.204, 0), rgb(0.31, 0.69, 0), rgb(0, 1, 0), rgb(0, 
            0.708, 0.292), rgb(0, 0.204, 0.796))

    } else if (i == 1) {

        maxdata <- 0.5
        mindata <- -0.5

        data1 <- mindata
        data2 <- mindata * 0.5
        data3 <- 0
        data4 <- maxdata * 0.5
        data5 <- maxdata

        plotblue <- data/maxdata
        plotred <- data/mindata
        plotblue[plotblue > 1] <- 1
        plotblue[plotblue < 0] <- 0
        plotred[plotred > 1] <- 1
        plotred[plotred < 0] <- 0
        plotgreen <- 1 - plotblue - plotred

        colleg <- c(rgb(1, 0, 0), rgb(0.5, 0.5, 0), rgb(0, 1, 0), rgb(0, 0.5, 0.5), 
            rgb(0, 0, 1))

    } else if (i < 4 || i == 6) {

        maxdata = 1
        mindata = 0

        data1 <- 0
        data2 <- 0.25
        data3 <- 0.5
        data4 <- 0.75
        data5 <- 1

        plotred <- (maxdata - data)/(maxdata - mindata)
        plotred[plotred > 1] <- 1
        plotred[plotred < 0] <- 0
        plotgreen <- (data - mindata)/(maxdata - mindata)
        plotgreen[plotgreen > 1] <- 1
        plotgreen[plotgreen < 0] <- 0
        plotblue <- 0 * data

        colleg <- c(rgb(1, 0, 0), rgb(0.75, 0.25, 0), rgb(0.5, 0.5, 0), rgb(0.25, 
            0.75, 0), rgb(0, 1, 0))

    } else if (i < 5) {

        maxdata = 1
        mindata = 0

        data1 <- 0
        data2 <- 0.25
        data3 <- 0.5
        data4 <- 0.75
        data5 <- 1

        plotred <- (data - mindata)/(maxdata - mindata)
        plotred[plotred > 1] <- 1
        plotred[plotred < 0] <- 0
        plotgreen <- (maxdata - data)/(maxdata - mindata)
        plotgreen[plotgreen > 1] <- 1
        plotgreen[plotgreen < 0] <- 0
        plotblue <- 0 * data

        colleg <- c(rgb(0, 1, 0), rgb(0.25, 0.75, 0), rgb(0.5, 0.5, 0), rgb(0.75, 
            0.25, 0), rgb(1, 0, 0))

    } else {

        data1 <- 0.1
        data2 <- 0.32
        data3 <- 1
        data4 <- 3.16
        data5 <- 10

        plotblue <- log10(data)
        plotred <- -log10(data)
        plotblue[plotblue < 0] <- 0
        plotblue[plotblue > 1] <- 1
        plotred[plotred < 0] <- 0
        plotred[plotred > 1] <- 1
        plotgreen <- 1 - plotred - plotblue

        colleg <- c(rgb(1, 0, 0), rgb(0.5, 0.5, 0), rgb(0, 1, 0), rgb(0, 0.5, 0.5), 
            rgb(0, 0, 1))
    }

    plotblue[data == -9999] <- 0.9
    plotred[data == -9999] <- 0.9
    plotgreen[data == -9999] <- 0.9
    plotcol <- rgb(plotred, plotgreen, plotblue)

    # Plotting points
    plot(x = par1, y = par2, axes = F, xlab = NA, ylab = NA, type = "p", xlim = c(minx, 
        maxx), ylim = c(miny, maxy), pch = 15, cex = 3, col = plotcol)

    # Plotting header
    par(xpd = TRUE)
    text(x = minx + (maxx - minx)/2, y = maxy + (maxy - miny)/4, labels = coltext[i], 
        cex = 1.2, col = "black")

    # Plotting legend
    ltext = vector("expression", 5)
    ltext[1] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data1), 
        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]
    ltext[2] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data2), 
        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]
    ltext[3] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data3), 
        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]
    ltext[4] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data4), 
        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]
    ltext[5] <- substitute(expression(text ~ funit), list(text = format(round(as.numeric(data5), 
        mdig[i]), nsmall = mdig[i]), funit = format(unit[i])))[2]

    legend("topleft", legend = ltext, pch = 15, pt.cex = 3, col = colleg, text.col = "black", 
        bty = "n", inset = c(-0.2, -0.15), horiz = TRUE, x.intersp = 1)

    # Plotting bounding box and axes
    box()
    axis(side = 1, tck = -0.02, labels = NA)  #x axis
    axis(side = 1, lwd = 0, line = -0.4)
    axis(side = 2, tck = -0.02, labels = NA)  #y axis
    axis(side = 2, lwd = 0, line = -0.4)

    # Plotting axis labels
    mtext(side = 1, paste("Parameter", as.character(ipar1), sep = " "), line = 1.8)  #x axis
    mtext(side = 2, paste("Parameter", as.character(ipar2), sep = " "), line = 1.8)  #y axis

    # Closing plot file
    dev.off()
}
