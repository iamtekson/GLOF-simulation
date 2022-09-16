#!/usr/bin/R

##############################################################################
#
# MODULE:       r.avaflow.profile.R
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The simulation model for avalanche and debris flows
#               Script for the creation of profile plots of model results
#
# COPYRIGHT:    (c) 2013 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2020 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
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

# Importing arguments defined in r.avaflow.py
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]  #prefix for file names
model <- as.integer(args[2])  #prefix for file names
ntimesteps <- as.integer(args[3])  #number of time step
cellsize <- as.integer(args[4])  #pixel size
hmax <- as.numeric(args[5])  #maximum flow height for display
smax <- as.numeric(args[6])  #maximum (PHASE 1) flow velocity for display
smaxf <- as.numeric(args[7])  #maximum PHASE 2 flow velocity for display
smaxw <- as.numeric(args[8])  #maximum PHASE 3 flow velocity for display
hexagg <- as.integer(args[9])  #exaggeration of flow height for display
fill <- as.character(args[10])  #complete string for time step
tint <- as.numeric(args[11])  #length of time step
tstop <- as.numeric(args[12])  #total time of simulation
depdef <- as.integer(args[13])  #control for observed deposit map
sparam <- as.integer(args[14])  #control for parameter to be displayed as bar plot
mconv <- args[15]  #conversion factor for bar plot

# Defining data file name
intable = paste(prefix, "_results/", prefix, "_files/", prefix, "_profile.txt", sep = "")

# Creating vectors from data
horpos <- read.table(intable, skip = 1)[, 1]  #horizontal position
elev <- read.table(intable, skip = 1)[, 2]  #elevation
if (depdef == 1) {
    hdeposit <- read.table(intable, skip = 1)[, 3]  #height of observed deposit
    hdeposit[hdeposit < 0] <- 0
}

stext = vector("expression", 7)  #initializing right y axis labels
ltext = vector("expression", 8)  #initializing legend labels
if (model <= 3) {
    # parameter to be displayed as bar diagram:
    ltext[2] = substitute(expression(italic(H)))[2]
    ltext[4] = substitute(expression(italic(H)[d]))[2]
    if (sparam == 1) {
        scol <- 6
        sstring <- "_tflow"
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
        stext[1] = substitute(expression(Flow ~ kinetic ~ energy ~ italic(T) ~ (funit)), 
            list(funit = format(munit)))[2]
        ltext[5] = substitute(expression(italic(T)))[2]
        col <- rgb(0.6, 0.3, 0.1, 0.5)
    } else if (sparam == 2) {
        scol <- 7
        sstring <- "_pflow"
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
        stext[1] = substitute(expression(Flow ~ pressure ~ italic(P) ~ (funit)), 
            list(funit = format(munit)))[2]
        ltext[5] = substitute(expression(italic(P)))[2]
        col <- rgb(0.6, 0.6, 0, 0.5)
    } else {
        scol <- 5
        sstring <- "_vflow"
        stext[1] = substitute(expression(Flow ~ velocity ~ italic(V) ~ (m/s)))[2]
        ltext[5] = substitute(expression(italic(T)))[2]
        col <- rgb(0, 0, 0.8, 0.5)
    }
} else if (model > 3 && model < 7) {
    ltext[1] = substitute(expression(italic(H)[P1]))[2]
    ltext[2] = substitute(expression(italic(H)[P2]))[2]
    ltext[5] = substitute(expression(italic(H)[P3]))[2]
    if (sparam == 1) {
        scol <- 8
        sstring <- "_tflow"
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
        stext[1] = substitute(expression(Flow ~ kinetic ~ energy ~ italic(T) ~ (funit)), 
            list(funit = format(munit)))[2]
        ltext[6] = substitute(expression(italic(T)[s]))[2]
        ltext[7] = substitute(expression(-italic(T)[f]))[2]
        cols <- rgb(0.6, 0.3, 0.1, 0.5)
        colf <- rgb(0, 0.5, 0.5, 0.5)
    } else if (sparam == 2) {
        scol <- 10
        sstring <- "_pflow"
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
        stext[1] = substitute(expression(Flow ~ pressure ~ italic(P) ~ (funit)), 
            list(funit = format(munit)))[2]
        ltext[6] = substitute(expression(italic(P)[s]))[2]
        ltext[7] = substitute(expression(-italic(P)[f]))[2]
        cols <- rgb(0.6, 0.6, 0, 0.5)
        colf <- rgb(0.2, 0, 0.6, 0.5)
    } else {
        scol <- 6
        sstring <- "_vflow"
        stext[1] = substitute(expression(Flow ~ velocity ~ italic(V) ~ (m/s)))[2]
        ltext[6] = substitute(expression(italic(V)[s]))[2]
        ltext[7] = substitute(expression(-italic(V)[f]))[2]
        colf <- rgb(0, 0, 0.8, 0.5)
        cols <- rgb(0.6, 0.4, 0, 0.5)
    }
} else if (model == 7) {
    ltext[1] = substitute(expression(italic(H)[P1]))[2]
    ltext[2] = substitute(expression(italic(H)[P2]))[2]
    ltext[3] = substitute(expression(-italic(H)[P3]))[2]
    ltext[5] = substitute(expression(italic(H)[d]))[2]
    if (sparam == 1) {
        scol <- 10
        sstring <- "_tflow"
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
        stext[1] = substitute(expression(Flow ~ kinetic ~ energy ~ italic(T) ~ (funit)), 
            list(funit = format(munit)))[2]
        ltext[6] = substitute(expression(italic(T)[P1]))[2]
        ltext[7] = substitute(expression(italic(T)[P2]))[2]
        ltext[8] = substitute(expression(-italic(T)[P3]))[2]
        cols <- rgb(0.8, 0, 0, 0.5)
        colf <- rgb(0, 0.8, 0, 0.5)
        colw <- rgb(0, 0, 0.8, 0.5)
    } else if (sparam == 2) {
        scol <- 13
        sstring <- "_pflow"
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
        stext[1] = substitute(expression(Flow ~ pressure ~ italic(P) ~ (funit)), 
            list(funit = format(munit)))[2]
        ltext[6] = substitute(expression(italic(P)[P1]))[2]
        ltext[7] = substitute(expression(italic(P)[P2]))[2]
        ltext[8] = substitute(expression(-italic(P)[P3]))[2]
        cols <- rgb(0.8, 0, 0, 0.5)
        colf <- rgb(0, 0.8, 0, 0.5)
        colw <- rgb(0, 0, 0.8, 0.5)
    } else {
        scol <- 7
        sstring <- "_vflow"
        stext[1] = substitute(expression(Flow ~ velocity ~ italic(V) ~ (m/s)))[2]
        ltext[6] = substitute(expression(italic(V)[P1]))[2]
        ltext[7] = substitute(expression(italic(V)[P2]))[2]
        ltext[8] = substitute(expression(-italic(V)[P3]))[2]
        cols <- rgb(0.8, 0, 0, 0.5)
        colf <- rgb(0, 0.8, 0, 0.5)
        colw <- rgb(0, 0, 0.8, 0.5)
    }
}
if (model <= 3) {

    hrelease <- read.table(intable, skip = 1)[, 3 + depdef]  #release height
    hflow <- read.table(intable, skip = 1)[, 5 * ntimesteps + 3 + depdef]  #flow height
    sflow <- read.table(intable, skip = 1)[, 5 * ntimesteps + scol - 1 + depdef]  #flow parameter
    hentr <- read.table(intable, skip = 1)[, 5 * ntimesteps + 7 + depdef]  #entrained/deposited depth

} else if (model > 3 && model < 7) {

    hreleases <- read.table(intable, skip = 1)[, 3 + depdef]  #PHASE 1 release height
    hreleasef <- read.table(intable, skip = 1)[, 4 + depdef]  #fluid release height
    hflows <- read.table(intable, skip = 1)[, 10 * ntimesteps + 3 + depdef]  #PHASE 1 flow height
    hflowf <- read.table(intable, skip = 1)[, 10 * ntimesteps + 4 + depdef]  #fluid flow height
    sflows <- read.table(intable, skip = 1)[, 10 * ntimesteps + scol - 1 + depdef]  #PHASE 1 flow parameter
    sflowf <- read.table(intable, skip = 1)[, 10 * ntimesteps + scol + depdef]  #fluid flow parameter
    hentrs <- read.table(intable, skip = 1)[, 10 * ntimesteps + 11 + depdef]  #PHASE 1 entrained/deposited depth
    hentrf <- read.table(intable, skip = 1)[, 10 * ntimesteps + 12 + depdef]  #fluid entrained/deposited depth
} else if (model == 7) {

    hreleases <- read.table(intable, skip = 1)[, 3 + depdef]  #PHASE 1 release height
    hreleasef <- read.table(intable, skip = 1)[, 4 + depdef]  #PHASE 2 release height
    hreleasew <- read.table(intable, skip = 1)[, 5 + depdef]  #PHASE 3 release height
    hflows <- read.table(intable, skip = 1)[, 15 * ntimesteps + 3 + depdef]  #PHASE 1 flow height
    hflowf <- read.table(intable, skip = 1)[, 15 * ntimesteps + 4 + depdef]  #PHASE 2 flow height
    hfloww <- read.table(intable, skip = 1)[, 15 * ntimesteps + 5 + depdef]  #PHASE 3 flow height
    sflows <- read.table(intable, skip = 1)[, 15 * ntimesteps + scol - 1 + depdef]  #PHASE 1 flow parameter
    sflowf <- read.table(intable, skip = 1)[, 15 * ntimesteps + scol + depdef]  #PHASE 2 flow parameter
    sfloww <- read.table(intable, skip = 1)[, 15 * ntimesteps + scol + 1 + depdef]  #PHASE 3 flow parameter
    hentrs <- read.table(intable, skip = 1)[, 15 * ntimesteps + 15 + depdef]  #PHASE 1 entrained/deposited depth
    hentrf <- read.table(intable, skip = 1)[, 15 * ntimesteps + 16 + depdef]  #PHASE 2 entrained/deposited depth
    hentrw <- read.table(intable, skip = 1)[, 15 * ntimesteps + 17 + depdef]  #PHASE 3 entrained/deposited depth
}

# Building geometry of profile plot
maxhorpos <- max(horpos, na.rm = TRUE)  #maximum of horizontal coordinate
minhorpos <- min(horpos, na.rm = TRUE)  #minimum of horizontal coordinate
topelev <- elev + hexagg * hmax  #maximum elevation of each pixel
maxelev <- max(topelev, na.rm = TRUE)  #maximum elevation over entire area
minelev <- min(elev, na.rm = TRUE)  #minimum elevation over entire area
maxelev <- maxelev + (maxelev - minelev) * 0.1  #updating maximum elevation over entire area (for legend)
minelev <- minelev - (maxelev - minelev) * 0.1  #updating maximum elevation over entire area (for texts)

# Creating plot file
profileplot = paste(prefix, "_results/", prefix, "_plots/", prefix, "_profiles_timesteps/", 
    prefix, sstring, fill, ".png", sep = "")
png(filename = profileplot, width = 15, height = 10, units = "cm", res = 300)

# Defining margins
par(mar = c(3.2, 3.2, 1, 3.2))

# Plotting bars for all time steps except initial condition:
if (ntimesteps > 0) {
    if (model <= 3) {
        barplot(sflow * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-smax * 
            1.4, smax * 1.4), border = NA, col = col)
        # flow parameter
    } else if (model > 3 && model < 7) {
        barplot(sflows * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 
            1.4, -smaxf * 1.4), max(smax * 1.4, smaxf * 1.4)), border = NA, col = cols)  #PHASE 1 flow parameter
        par(new = TRUE)
        barplot(-sflowf * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 
            1.4, -smaxf * 1.4), max(smax * 1.4, smaxf * 1.4)), border = NA, col = colf)  #PHASE 2 flow parameter
    } else if (model == 7) {
        barplot((sflows + sflowf) * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, 
            ylim = c(min(-smax * 1.4 - smaxf * 1.4, -smaxw * 1.4), max(smax * 1.4 + 
                smaxf * 1.4, smaxw * 1.4)), border = NA, col = cols)  #PHASE 1 flow parameter
        par(new = TRUE)
        barplot(sflowf * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 
            1.4 - smaxf * 1.4, -smaxw * 1.4), max(smax * 1.4 + smaxf * 1.4, smaxw * 
            1.4)), border = NA, col = colf)  #PHASE 2 flow parameter
        par(new = TRUE)
        barplot(-sfloww * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(min(-smax * 
            1.4 - smaxf * 1.4, -smaxw * 1.4), max(smax * 1.4 + smaxf * 1.4, smaxw * 
            1.4)), border = NA, col = colw)  #PHASE 3 flow parameter
    }
    axis(side = 4, tck = -0.02, labels = NA)  #y axis for velocity
    axis(side = 4, lwd = 0, line = -0.4)
    par(new = TRUE)
}

# Plotting lines
plot(x = horpos, y = elev, axes = F, xlab = NA, ylab = NA, type = "l", ylim = c(minelev, 
    maxelev), lty = 1, lwd = 1, col = "gray")  #elevation

if (depdef == 1) lines(x = horpos, y = elev + hexagg * hdeposit, lty = 3, lwd = 1, 
    col = "black")  #height of deposit

if (model <= 3) {
    lines(x = horpos, y = elev + hexagg * hrelease, lty = 3, lwd = 1, col = "red")  #initial flow height
    lines(x = horpos, y = elev + hexagg * hentr, lty = 1, lwd = 1, col = "black")  #entrained/deposited depth
    lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] + 
        hexagg * (hflow[2:(length(horpos) - 1)] + hentr[2:(length(horpos) - 1)]) + 
        0 * (1/(hflow[2:(length(horpos) - 1)] + hflow[1:(length(horpos) - 2)] + hflow[3:length(horpos)])), 
        lty = 1, lwd = 2, col = "red")  #flow height
} else if (model > 3 && model < 7) {
    lines(x = horpos, y = elev + hexagg * (hreleases + hreleasef), lty = 3, lwd = 1, 
        col = rgb(0, 0, 0.8))  #initial PHASE 2 flow height
    lines(x = horpos, y = elev + hexagg * hreleases, lty = 3, lwd = 1, col = rgb(0.45, 
        0.3, 0))  #initial PHASE 1 flow height
    lines(x = horpos, y = elev + hexagg * (hentrs + hentrf), lty = 1, lwd = 1, col = "black")  #entrained/deposited PHASE 1 depth
    lines(x = horpos, y = elev + hexagg * hentrf, lty = 1, lwd = 1, col = "black")  #entrained/deposited PHASE 2 depth
    lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] + 
        hexagg * (hflows[2:(length(horpos) - 1)] + hflowf[2:(length(horpos) - 1)] + 
            hentrs[2:(length(horpos) - 1)] + hentrf[2:(length(horpos) - 1)]) + 0 * 
        (1/(hflowf[2:(length(horpos) - 1)] + hflowf[1:(length(horpos) - 2)] + hflowf[3:length(horpos)])), 
        lty = 1, lwd = 2, col = rgb(0, 0, 0.8))  #PHASE 2 flow height
    lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] + 
        hexagg * (hflows[2:(length(horpos) - 1)] + hentrs[2:(length(horpos) - 1)] + 
            hentrf[2:(length(horpos) - 1)]) + 0 * (1/(hflows[2:(length(horpos) - 
        1)] + hflows[1:(length(horpos) - 2)] + hflows[3:length(horpos)])), lty = 1, 
        lwd = 2, col = rgb(0.45, 0.3, 0))  #PHASE 1 flow height
} else if (model == 7) {
    lines(x = horpos, y = elev + hexagg * (hreleases + hreleasef + hreleasew), lty = 3, 
        lwd = 1, col = rgb(0, 0, 0.8))  #initial PHASE 3 flow height
    lines(x = horpos, y = elev + hexagg * (hreleases + hreleasef), lty = 3, lwd = 1, 
        col = rgb(0, 0.8, 0))  #initial PHASE 2 flow height
    lines(x = horpos, y = elev + hexagg * hreleases, lty = 3, lwd = 1, col = rgb(0.8, 
        0, 0))  #initial PHASE 1 flow height
    lines(x = horpos, y = elev + hexagg * (hentrs + hentrf + hentrw), lty = 1, lwd = 1, 
        col = "black")  #entrained/deposited PHASE 3 depth
    lines(x = horpos, y = elev + hexagg * (hentrs + hentrf), lty = 1, lwd = 1, col = "black")  #entrained/deposited PHASE 1 depth
    lines(x = horpos, y = elev + hexagg * hentrf, lty = 1, lwd = 1, col = "black")  #entrained/deposited PHASE 2 depth
    lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] + 
        hexagg * (hflows[2:(length(horpos) - 1)] + hflowf[2:(length(horpos) - 1)] + 
            hfloww[2:(length(horpos) - 1)] + hentrs[2:(length(horpos) - 1)] + hentrf[2:(length(horpos) - 
            1)] + hentrw[2:(length(horpos) - 1)]) + 0 * (1/(hfloww[2:(length(horpos) - 
        1)] + hfloww[1:(length(horpos) - 2)] + hfloww[3:length(horpos)])), lty = 1, 
        lwd = 2, col = rgb(0, 0, 0.8))  #PHASE 3 flow height
    lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] + 
        hexagg * (hflows[2:(length(horpos) - 1)] + hflowf[2:(length(horpos) - 1)] + 
            hentrs[2:(length(horpos) - 1)] + hentrf[2:(length(horpos) - 1)] + hentrw[2:(length(horpos) - 
            1)]) + 0 * (1/(hflowf[2:(length(horpos) - 1)] + hflowf[1:(length(horpos) - 
        2)] + hflowf[3:length(horpos)])), lty = 1, lwd = 2, col = rgb(0, 0.8, 0))  #PHASE 2 flow height
    lines(x = horpos[2:(length(horpos) - 1)], y = elev[2:(length(horpos) - 1)] + 
        hexagg * (hflows[2:(length(horpos) - 1)] + hentrs[2:(length(horpos) - 1)] + 
            hentrf[2:(length(horpos) - 1)] + hentrw[2:(length(horpos) - 1)]) + 0 * 
        (1/(hflows[2:(length(horpos) - 1)] + hflows[1:(length(horpos) - 2)] + hflows[3:length(horpos)])), 
        lty = 1, lwd = 2, col = rgb(0.8, 0, 0))  #PHASE 1 flow height
}

# Legends

if (model <= 3) {

    legend("topright", legend = ltext[2], lty = 1, lwd = 2, col = "red", text.col = "red", 
        bty = "n", inset = c(0.825, -0.025), horiz = TRUE, x.intersp = 0.5)  #flow height

    if (depdef == 1) 
        legend("topright", legend = ltext[4], lty = 3, lwd = 2, col = "black", text.col = "black", 
            bty = "n", inset = c(0.425, -0.025), horiz = TRUE, x.intersp = 0.5)  #height of deposit

    if (ntimesteps > 0) {
        legend("topright", legend = ltext[5], border = NA, fill = col, text.col = col, 
            bty = "n", inset = c(0.225, -0.025), horiz = TRUE)  #bar plot parameter
    }

} else if (model > 3 && model < 7) {

    legend("topright", legend = ltext[1], lty = 1, lwd = 2, col = rgb(0.45, 0.3, 
        0), text.col = rgb(0.45, 0.3, 0), bty = "n", inset = c(0.825, -0.025), horiz = TRUE, 
        x.intersp = 0.5)  #PHASE 1 flow height
    legend("topright", legend = ltext[2], lty = 1, lwd = 2, col = rgb(0, 0, 0.8), 
        text.col = rgb(0, 0, 0.8), bty = "n", inset = c(0.625, -0.025), horiz = TRUE, 
        x.intersp = 0.5)  #PHASE 2 flow height

    if (depdef == 1) 
        legend("topright", legend = ltext[5], lty = 3, lwd = 2, col = "black", text.col = "black", 
            bty = "n", inset = c(0.425, -0.025), horiz = TRUE, x.intersp = 0.5)  #height of deposit

    if (ntimesteps > 0) {
        legend("topright", legend = ltext[6], border = NA, fill = cols, text.col = cols, 
            bty = "n", inset = c(0.225, -0.025), horiz = TRUE)  #PHASE 1 bar plot parameter
        legend("topright", legend = ltext[7], border = NA, fill = colf, text.col = colf, 
            bty = "n", inset = c(0.025, -0.025), horiz = TRUE)  #PHASE 2 bar plot parameter
    }
} else if (model == 7) {

    par(cex = 0.75)
    legend("topright", legend = ltext[1], lty = 1, lwd = 2, col = rgb(0.8, 0, 0), 
        text.col = rgb(0.8, 0, 0), bty = "n", inset = c(0.85, -0.025), horiz = TRUE, 
        x.intersp = 0.5)  #PHASE 1 flow height
    legend("topright", legend = ltext[2], lty = 1, lwd = 2, col = rgb(0, 0.8, 0), 
        text.col = rgb(0, 0.8, 0), bty = "n", inset = c(0.72, -0.025), horiz = TRUE, 
        x.intersp = 0.5)  #PHASE 2 flow height
    legend("topright", legend = ltext[3], lty = 1, lwd = 2, col = rgb(0, 0, 0.8), 
        text.col = rgb(0, 0, 0.8), bty = "n", inset = c(0.575, -0.025), horiz = TRUE, 
        x.intersp = 0.5)  #PHASE 3 flow height

    if (depdef == 1) 
        legend("topright", legend = ltext[5], lty = 3, lwd = 2, col = "black", text.col = "black", 
            bty = "n", inset = c(0.425, -0.025), horiz = TRUE, x.intersp = 0.5)  #height of deposit

    if (ntimesteps > 0) {
        legend("topright", legend = ltext[6], border = NA, fill = cols, text.col = cols, 
            bty = "n", inset = c(0.275, -0.025), horiz = TRUE)  #PHASE 1 bar plot parameter
        legend("topright", legend = ltext[7], border = NA, fill = colf, text.col = colf, 
            bty = "n", inset = c(0.165, -0.025), horiz = TRUE)  #PHASE 2 bar plot parameter
        legend("topright", legend = ltext[8], border = NA, fill = colw, text.col = colw, 
            bty = "n", inset = c(0.025, -0.025), horiz = TRUE)  #PHASE 3 bar plot parameter
    }
    par(cex = 1)
}

# Plotting bounding box and axes
box()
axis(side = 1, tck = -0.02, labels = NA)  #x axis
axis(side = 1, lwd = 0, line = -0.4)
axis(side = 2, tck = -0.02, labels = NA)  #y axis for lines
axis(side = 2, lwd = 0, line = -0.4)

# Plotting axis labels
mtext(side = 1, "Horizontal distance (m)", line = 2)  #x axis
mtext(side = 2, "Elevation (m)", line = 2)  #y axis for lines
if (ntimesteps > 0) mtext(side = 4, stext[1], line = 2)  #y axis for bars


# Texts showing time passed since start of flow and exaggeration of flow height
secpassed <- tint * ntimesteps  #time passed since start of flow
if (secpassed > tstop) secpassed <- tstop
textpassed = vector("expression", 1)  #initializing text label
textpassed[1] = substitute(expression(italic(t) == secformat ~ s), list(secformat = format(round(secpassed, 
    1), nsmall = 1)))[2]  #defining text label
text(x = minhorpos + (maxhorpos - minhorpos)/6, y = minelev + (maxelev - minelev)/25, 
    labels = textpassed[1], cex = 1.4, col = "black", font = 2)  #printing text for time passed since start of flow

if (hexagg != 1) {
    texthexagg <- paste(hexagg, "-fold exaggeration of flow height", sep = "")  #text label for exaggeration of flow height
    text(x = minhorpos + (maxhorpos - minhorpos) * 4/6, y = minelev + (maxelev - 
        minelev)/35, labels = texthexagg, col = "black")  #printing text for exaggeration of flow height
}

# Closing plot file
dev.off()
